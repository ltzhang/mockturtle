/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2024  EPFL
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/*!
  \file emap.hpp
  \brief An extended technology mapper

  \author Alessandro Tempia Calvino
*/

#pragma once

#include <chrono>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <cstdint>
#include <limits>
#include <memory>
#include <cmath>
#include <stdexcept>
#include <type_traits>
#include <string>
#include <unordered_map>
#include <vector>

#include <kitty/constructors.hpp>
#include <kitty/dynamic_truth_table.hpp>
#include <kitty/hash.hpp>
#include <kitty/static_truth_table.hpp>

#include <fmt/format.h>
#include <parallel_hashmap/phmap.h>

#include "../networks/aig.hpp"
#include "../networks/block.hpp"
#include "../networks/klut.hpp"
#include "../utils/cuts.hpp"
#include "../utils/load_curve.hpp"
#include "../utils/node_map.hpp"
#include "../utils/stopwatch.hpp"
#include "../utils/tech_library.hpp"
#include "../views/binding_view.hpp"
#include "../views/cell_view.hpp"
#include "../views/choice_view.hpp"
#include "../views/topo_view.hpp"
#include "cleanup.hpp"
#include "cut_enumeration.hpp"
#include "detail/mffc_utils.hpp"
#include "detail/switching_activity.hpp"

namespace mockturtle
{

/*! \brief Parameters for emap.
 *
 * The data structure `emap_params` holds configurable parameters
 * with default arguments for `emap`.
 */
struct emap_params
{
  emap_params()
  {
    cut_enumeration_ps.cut_limit = 16;
    cut_enumeration_ps.minimize_truth_table = true;
  }

  /*! \brief Parameters for cut enumeration
   *
   * The default cut limit is 16.
   * The maximum cut limit is 19.
   * By default, truth table minimization
   * is performed.
   */
  cut_enumeration_params cut_enumeration_ps{};

  /*! \brief Do area-oriented mapping. */
  bool area_oriented_mapping{ false };

  /*! \brief Maps using multi-output gates */
  bool map_multioutput{ false };

  /*! \brief Matching mode
   *
   * Boolean uses Boolean matching (up to 6-input cells),
   * Structural uses pattern matching for fully-DSD cells,
   * Hybrid combines the two.
   */
  enum matching_mode_t
  {
    boolean,
    structural,
    hybrid
  } matching_mode = hybrid;

  /*! \brief Target required time (for each PO). */
  double required_time{ 0.0f };

  /*! \brief Required time relaxation in percentage (10 = 10%). */
  double relax_required{ 0.0f };

  /*! \brief Custom input arrival times. */
  std::vector<double> arrival_times{};

  /*! \brief Custom output required times. */
  std::vector<double> required_times{};

  /*! \brief Number of rounds for area flow optimization. */
  uint32_t area_flow_rounds{ 3u };

  /*! \brief Number of rounds for exact area optimization. */
  uint32_t ela_rounds{ 2u };

  /*! \brief Number of rounds for exact switching power optimization. */
  uint32_t eswp_rounds{ 0u };

  /*! \brief Number of patterns for switching activity computation. */
  uint32_t switching_activity_patterns{ 2048u };

  /*! \brief Electrical (load-aware) cost: price a candidate's pin delay as
   * block + slope × the capacitive load the CURRENT cover imposes on the node's output, prefer
   * candidates whose output limit (max_load) covers that load, and refresh the per-node load
   * estimate from the committed choices after every mapping round. Requires a tech_library built
   * with tech_library_params::electrical_model (otherwise every slope/cap is 0 and this flag has
   * no effect). Off: byte-identical mapping. */
  bool electrical_model{ false };

  /*! \brief Wire capacitance charged per fanout edge in the electrical load estimate. */
  double wire_cap_per_fanout{ 0.0 };

  /*! Optional total wire capacitance C(F), in the library's capacitance units.
   * Driver excluded, output terminals included. Replaces the scalar when set.
   * Total must be finite/nonnegative; a table need not be monotone. */
  std::function<double(std::uint32_t)> wire_cap_for_fanout;

  /*! \brief Bounded per-node area-delay Pareto frontier (electrical model only). Each matched
   * node keeps up to this many non-dominated (arrival, area) points per phase, and at each
   * non-ELA round end every covered node is re-bound to the cheapest same-cut point that still
   * meets its required time — now known from the COMPLETE cover, which is what per-match
   * selection cannot see. 0 (default) = one best match per phase, byte-identical. Hard-capped
   * at 8 points (see max_curve_points). */
  uint32_t curve_points{ 0 };

  /*! \brief Load-indexed arrival curves (electrical model only): store each node's arrival as a
   * FUNCTION of its output load, sampled at this many library-anchored load points, and price a
   * candidate's leaves at the capacitance that candidate itself presents rather than at the one
   * scalar estimate the leaf happens to carry.
   *
   * Without this, a candidate's input pin capacitance — the one electrical term matching could
   * know exactly — is invisible: a gate with large input pins slows its fanins and is never
   * charged for it. 0 (default) or 1 = one arrival per node, byte-identical (a single load point
   * IS the scalar model). Hard-capped at `max_load_points`; a library declaring no drive limit
   * anywhere has no ceiling to span and leaves the curves inactive. */
  uint32_t load_points{ 0 };

  /*! \brief Optional single-input covering stage (electrical model only): after each load
   * refresh, every covered node may CHOOSE a library buffer on its output — the driver then sees
   * only the buffer's input capacitance while the buffer carries the sink load — priced by the
   * same block + slope × load model as every other candidate, with "no buffer" (the zero-cost
   * identity) always competing. Off (default): byte-identical mapping. */
  bool cover_buffer{ false };

  /*! \brief Frontier through the covering DP (delay-oriented mapping only). After area recovery,
   * every node keeps a bounded Pareto frontier of (arrival, area-flow) points per output phase,
   * built by merging its leaves' frontiers through EVERY cut and match -- not only the chosen cut,
   * which is all `curve_points` re-binds within. A backward pass then resolves each node against
   * its required time (the cheapest point that meets it), which fixes the node's cut and match and
   * pushes a required time onto each leaf, so the cover can trade a slower, cheaper leaf for a
   * faster match above it. One exact-area round cleans up after. 0 (default) = off, byte-identical.
   * Hard-capped at `max_frontier_points`. */
  uint32_t frontier_points{ 0 };

  /*! \brief How a frontier point charges a leaf's cost: `flow` divides it by the blended reference
   * estimate (area flow), `cover` by the entering cover's actual reference count -- a leaf the cover
   * does not use pays its whole cost, which keeps the resolution from buying new, unshared logic
   * on the strength of fanout that is not there. `cover` measured better on real libraries. */
  enum class frontier_share_t
  {
    flow,
    cover
  };
  frontier_share_t frontier_share{ frontier_share_t::cover };

  /*! \brief Resolution rounds: each re-derives the sharing from the cover the previous round kept,
   * and the loop stops at the first round that keeps nothing. Clamped to 1..8. */
  uint32_t frontier_rounds{ 4 };

  /*! \brief Memory ceiling for total incremental frontier storage, in MiB. Checked BEFORE allocation; a network
   * over it keeps the ordinary cover and reports the decline. */
  double frontier_mem_budget_mb{ 1024.0 };

  /*! \brief Compute area-oriented alternative matches */
  bool use_match_alternatives{ true };

  /*! \brief Remove the cuts that are contained in others */
  bool remove_dominated_cuts{ false };

  /*! \brief Remove overlapping multi-output cuts */
  bool remove_overlapping_multicuts{ false };

  /*! \brief Be verbose. */
  bool verbose{ false };
};

/*! \brief Statistics for emap.
 *
 * The data structure `emap_stats` provides data collected by running
 * `emap`.
 */
struct emap_stats
{
  /*! \brief Area result. */
  double area{ 0 };
  /*! \brief Worst delay result. */
  double delay{ 0 };
  /*! \brief Power result. */
  /*! \brief Bound cells whose ESTIMATED load exceeds their own drive limit (electrical model
   * only). Reported, never a rejection: if no candidate for a node is drive-legal, the strongest
   * is bound and counted here. A cell that declared no max_load is excluded (0 = unconstrained). */
  uint32_t max_load_violations{ 0 };
  /*! \brief Worst estimated load / drive limit ratio over constrained drivers (electrical model
   * only). */
  double worst_load_ratio{ 0 };
  /*! \brief Output buffers the optional single-input covering stage chose (cover_buffer only). */
  uint32_t cover_buffers{ 0 };
  /*! \brief Load points the arrival curves were actually sampled at (`load_points` only). 0 means
   * the curves stayed inactive — the request was off, a single point, or the library declared no
   * drive limit to anchor a ceiling on — so a caller can tell "asked for" from "got". */
  uint32_t load_points{ 0 };
  /*! \brief Floor and ceiling of the load ladder, in the library's capacitance unit (0 when the
   * curves are inactive). */
  double load_curve_lo{ 0 };
  double load_curve_hi{ 0 };
  /*! \brief Bytes the curve store occupies (0 when inactive). */
  uint64_t load_curve_bytes{ 0 };

  /*! \brief Frontier through the covering DP (`frontier_points`). `frontier_ran` says the round
   * committed a cover; otherwise `frontier_declined` names why it did not (empty when the round was
   * never requested or not applicable, e.g. an area-oriented run). */
  bool frontier_ran{ false };
  /*! \brief The resolved cover was kept: it met the entering cover's delay target at a smaller
   * exact area. When false after a run, the entering cover was restored unchanged. */
  bool frontier_kept{ false };
  uint32_t frontier_rounds_kept{ 0 }; /* resolution rounds that were kept */
  std::string frontier_declined{};
  uint32_t frontier_nodes{ 0 };        /* nodes that carried a frontier */
  uint32_t frontier_max_size{ 0 };     /* largest per-phase frontier kept */
  uint32_t frontier_trimmed{ 0 };      /* per-phase frontiers the bound cut down */
  uint64_t frontier_bytes{ 0 };        /* total conservative incremental admission estimate */
  uint64_t frontier_bytes_points{ 0 };
  uint64_t frontier_bytes_headers{ 0 };
  uint64_t frontier_bytes_saved_matches{ 0 };
  uint64_t frontier_bytes_sharing{ 0 };
  uint64_t frontier_bytes_resolution{ 0 };
  uint64_t frontier_bytes_scratch{ 0 };
  uint64_t frontier_bytes_allocator{ 0 }; /* page rounding + allocator allowance, not baseline RSS */
  uint32_t frontier_cut_changes{ 0 };  /* covered (node, phase)s resolved onto a different cut */
  uint32_t frontier_gate_changes{ 0 }; /* covered (node, phase)s resolved onto a different match */
  uint32_t frontier_unmet{ 0 };        /* demanded (node, phase)s no point could meet (fastest used) */
  double frontier_area_before{ 0 };    /* cover area entering the round */
  double frontier_area_after{ 0 };     /* after the resolution and its exact-area clean-up */
  double frontier_delay_before{ 0 };
  double frontier_delay_after{ 0 };

  double power{ 0 };
  /*! \brief Power result. */
  uint32_t inverters{ 0 };

  /*! \brief Mapped multi-output gates. */
  uint32_t multioutput_gates{ 0 };

  /*! \brief Runtime for multi-output matching. */
  stopwatch<>::duration time_multioutput{ 0 };
  /*! \brief Total runtime. */
  stopwatch<>::duration time_total{ 0 };

  /*! \brief Cut enumeration stats. */
  cut_enumeration_stats cut_enumeration_st{};

  /*! \brief Delay and area stats for each round. */
  std::vector<std::string> round_stats{};

  /*! \brief Mapping error. */
  bool mapping_error{ false };

  void report() const
  {
    for ( auto const& stat : round_stats )
    {
      std::cout << stat;
    }
    std::cout << fmt::format( "[i] Area = {:>5.2f}; Delay = {:>5.2f};", area, delay );
    if ( power != 0 )
      std::cout << fmt::format( " Power = {:>5.2f};\n", power );
    else
      std::cout << "\n";
    if ( multioutput_gates )
    {
      std::cout << fmt::format( "[i] Multi-output gates   = {:>5}\n", multioutput_gates );
      std::cout << fmt::format( "[i] Multi-output runtime = {:>5.2f} secs\n", to_seconds( time_multioutput ) );
    }
    std::cout << fmt::format( "[i] Total runtime        = {:>5.2f} secs\n", to_seconds( time_total ) );
  }
};

namespace detail
{

/* Conservative per-array live allocation bound, qualified for glibc and mimalloc.
 * Include up to 64 bytes of array cookie/alignment padding BEFORE size-class rounding,
 * then one page for allocator metadata/rounding. This is not an arena or RSS bound. */
inline bool frontier_allocation_bytes( uint64_t request, uint64_t& bytes )
{
  constexpr uint64_t maximum = std::numeric_limits<uint64_t>::max();
  if ( request > maximum - 64 )
    return false;
  uint64_t const padded = request + 64;
  uint64_t rounded = 1;
  while ( rounded < padded )
  {
    if ( rounded > maximum / 2 )
      return false;
    rounded *= 2;
  }
  if ( rounded > maximum - 4096 )
    return false;
  bytes = rounded + 4096;
  return true;
}

#pragma region cut set
template<unsigned NInputs>
struct cut_enumeration_emap_cut
{
  /* stats */
  uint32_t delay;
  float flow;
  bool ignore;

  /* pattern index for structural matching*/
  uint32_t pattern_index;

  /* function */
  kitty::static_truth_table<6> function;

  /* list of supergates matching the cut for positive and negative output phases */
  std::array<std::vector<supergate<NInputs>> const*, 2> supergates;
  /* input negations, 0: pos, 1: neg */
  std::array<uint16_t, 2> negations;
};

struct cut_enumeration_emap_multi_cut
{
  /* stats */
  uint64_t id{ 0 };
};

enum class emap_cut_sort_type
{
  DELAY = 0,
  DELAY2 = 1,
  AREA = 2,
  AREA2 = 3,
  NONE = 4
};

template<typename CutType, uint32_t MaxCuts>
class emap_cut_set
{
public:
  /*! \brief Standard constructor.
   */
  emap_cut_set()
  {
    clear();
  }

  /*! \brief Assignment operator.
   */
  emap_cut_set& operator=( emap_cut_set const& other )
  {
    if ( this != &other )
    {
      _pcend = _pend = _pcuts.begin();
      _set_limit = other._set_limit;

      auto it = other.begin();
      while ( it != other.end() )
      {
        **_pend++ = **it++;
        ++_pcend;
      }
    }

    return *this;
  }

  /*! \brief Clears a cut set.
   */
  void clear()
  {
    _pcend = _pend = _pcuts.begin();
    auto pit = _pcuts.begin();
    for ( auto& c : _cuts )
    {
      *pit++ = &c;
    }
  }

  /*! \brief Sets the cut limit.
   */
  void set_cut_limit( uint32_t limit )
  {
    _set_limit = std::min( MaxCuts, limit );
  }

  /*! \brief Adds a cut to the end of the set.
   *
   * This function should only be called to create a set of cuts which is known
   * to be sorted and irredundant (i.e., no cut in the set dominates another
   * cut).
   *
   * \param begin Begin iterator to leaf indexes
   * \param end End iterator (exclusive) to leaf indexes
   * \return Reference to the added cut
   */
  template<typename Iterator>
  CutType& add_cut( Iterator begin, Iterator end )
  {
    assert( _pend != _pcuts.end() );

    auto& cut = **_pend++;
    cut.set_leaves( begin, end );

    ++_pcend;
    return cut;
  }

  /*! \brief Appends a cut to the end of the set.
   *
   * This function should only be called to create a set of cuts which is known
   * to be sorted and irredundant (i.e., no cut in the set dominates another
   * cut).
   *
   * \param cut Cut to insert
   */
  void append_cut( CutType const& cut )
  {
    assert( _pend != _pcuts.end() );

    **_pend++ = cut;
    ++_pcend;
  }

  /*! \brief Checks whether cut is dominates by any cut in the set.
   *
   * \param cut Cut outside of the set
   */
  bool is_dominated( CutType const& cut ) const
  {
    return std::find_if( _pcuts.begin(), _pcend, [&cut]( auto const* other ) { return other->dominates( cut ); } ) != _pcend;
  }

  static bool sort_delay( CutType const& c1, CutType const& c2 )
  {
    constexpr auto eps{ 0.005f };
    if ( !c1->ignore && c2->ignore )
      return true;
    if ( c1->ignore && !c2->ignore )
      return false;
    if ( c1->delay < c2->delay - eps )
      return true;
    if ( c1->delay > c2->delay + eps )
      return false;
    if ( c1->flow < c2->flow - eps )
      return true;
    if ( c1->flow > c2->flow + eps )
      return false;
    return c1.size() < c2.size();
  }

  static bool sort_delay2( CutType const& c1, CutType const& c2 )
  {
    constexpr auto eps{ 0.005f };
    if ( !c1->ignore && c2->ignore )
      return true;
    if ( c1->ignore && !c2->ignore )
      return false;
    if ( c1.size() < c2.size() )
      return true;
    if ( c1.size() > c2.size() )
      return false;
    if ( c1->delay < c2->delay - eps )
      return true;
    if ( c1->delay > c2->delay + eps )
      return false;
    return c1->flow < c2->flow - eps;
  }

  static bool sort_area( CutType const& c1, CutType const& c2 )
  {
    constexpr auto eps{ 0.005f };
    if ( !c1->ignore && c2->ignore )
      return true;
    if ( c1->ignore && !c2->ignore )
      return false;
    if ( c1->flow < c2->flow - eps )
      return true;
    if ( c1->flow > c2->flow + eps )
      return false;
    if ( c1.size() < c2.size() )
      return true;
    if ( c1.size() > c2.size() )
      return false;
    return c1->delay < c2->delay - eps;
  }

  static bool sort_area2( CutType const& c1, CutType const& c2 )
  {
    constexpr auto eps{ 0.005f };
    if ( !c1->ignore && c2->ignore )
      return true;
    if ( c1->ignore && !c2->ignore )
      return false;
    if ( c1->flow < c2->flow - eps )
      return true;
    if ( c1->flow > c2->flow + eps )
      return false;
    if ( c1->delay < c2->delay - eps )
      return true;
    if ( c1->delay > c2->delay + eps )
      return false;
    return c1.size() < c2.size();
  }

  /*! \brief Compare two cuts using sorting functions.
   *
   * This method compares two cuts using a sorting function.
   *
   * \param cut1 first cut.
   * \param cut2 second cut.
   * \param sort sorting function.
   */
  static bool compare( CutType const& cut1, CutType const& cut2, emap_cut_sort_type sort = emap_cut_sort_type::NONE )
  {
    if ( sort == emap_cut_sort_type::DELAY )
    {
      return sort_delay( cut1, cut2 );
    }
    else if ( sort == emap_cut_sort_type::DELAY2 )
    {
      return sort_delay2( cut1, cut2 );
    }
    else if ( sort == emap_cut_sort_type::AREA )
    {
      return sort_area( cut1, cut2 );
    }
    else if ( sort == emap_cut_sort_type::AREA2 )
    {
      return sort_area2( cut1, cut2 );
    }
    else
    {
      return false;
    }
  }

  /*! \brief Inserts a cut into a set without checking dominance.
   *
   * This method will insert a cut into a set and maintain an order.  This
   * method doesn't remove the cuts that are dominated by `cut`.
   *
   * If `cut` is dominated by any of the cuts in the set, it will still be
   * inserted.  The caller is responsible to check whether `cut` is dominated
   * before inserting it into the set.
   *
   * \param cut Cut to insert.
   * \param sort Cut prioritization function.
   */
  void simple_insert( CutType const& cut, emap_cut_sort_type sort = emap_cut_sort_type::NONE )
  {
    /* insert cut in a sorted way */
    typename std::array<CutType*, MaxCuts>::iterator ipos = _pcuts.begin();

    bool limit_reached = std::distance( _pcuts.begin(), _pend ) >= _set_limit;

    /* do not insert if worst than set_limit */
    if ( limit_reached )
    {
      if ( sort == emap_cut_sort_type::AREA && !sort_area( cut, **( ipos + _set_limit - 1 ) ) )
      {
        return;
      }
      else if ( sort != emap_cut_sort_type::AREA )
      {
        return;
      }
    }

    if ( sort == emap_cut_sort_type::NONE )
    {
      ipos = _pend;
    }
    else /* AREA */
    {
      ipos = std::upper_bound( _pcuts.begin(), _pend, &cut, []( auto a, auto b ) { return sort_area( *a, *b ); } );
    }

    /* check for redundant cut */
    typename std::array<CutType*, MaxCuts>::iterator jpos = ipos;
    if ( cut->ignore )
    {
      while ( jpos != _pcuts.begin() )
      {
        --jpos;
        if ( ( *jpos )->size() < cut.size() )
          break;
        if ( ( *jpos )->signature() == cut.signature() && std::equal( cut.begin(), cut.end(), ( *jpos )->begin() ) )
          return;
      }
    }
    else if ( ipos != _pcuts.begin() )
    {
      if ( ( *( ipos - 1 ) )->signature() == cut.signature() && std::equal( cut.begin(), cut.end(), ( *( ipos - 1 ) )->begin() ) )
      {
        return;
      }
    }

    /* too many cuts, we need to remove one */
    if ( _pend == _pcuts.end() || limit_reached )
    {
      /* cut to be inserted is worse than all the others, return */
      if ( ipos == _pend )
      {
        return;
      }
      else
      {
        /* remove last cut */
        --_pend;
        --_pcend;
      }
    }

    /* copy cut */
    auto& icut = *_pend;
    icut->set_leaves( cut.begin(), cut.end() );
    icut->data() = cut.data();

    if ( ipos != _pend )
    {
      auto it = _pend;
      while ( it > ipos )
      {
        std::swap( *it, *( it - 1 ) );
        --it;
      }
    }

    /* update iterators */
    _pcend++;
    _pend++;
  }

  /*! \brief Inserts a cut into a set.
   *
   * This method will insert a cut into a set and maintain an order.  Before the
   * cut is inserted into the correct position, it will remove all cuts that are
   * dominated by `cut`. Variable `skip0` tell to skip the dominance check on
   * cut zero.
   *
   * If `cut` is dominated by any of the cuts in the set, it will still be
   * inserted.  The caller is responsible to check whether `cut` is dominated
   * before inserting it into the set.
   *
   * \param cut Cut to insert.
   * \param skip0 Skip dominance check on cut zero.
   * \param sort Cut prioritization function.
   */
  void insert( CutType const& cut, bool skip0 = false, emap_cut_sort_type sort = emap_cut_sort_type::NONE )
  {
    auto begin = _pcuts.begin();

    if ( skip0 && _pend != _pcuts.begin() )
      ++begin;

    /* remove elements that are dominated by new cut */
    _pcend = _pend = std::stable_partition( begin, _pend, [&cut]( auto const* other ) { return !cut.dominates( *other ); } );

    /* insert cut in a sorted way */
    simple_insert( cut, sort );
  }

  /*! \brief Replaces a cut of the set.
   *
   * This method replaces the cut at position `index` in the set by `cut`
   * and maintains the cuts order. The function does not check whether
   * index is in the valid range.
   *
   * \param index Index of the cut to replace.
   * \param cut Cut to insert.
   */
  void replace( uint32_t index, CutType const& cut )
  {
    *_pcuts[index] = cut;
  }

  /*! \brief Begin iterator (constant).
   *
   * The iterator will point to a cut pointer.
   */
  auto begin() const { return _pcuts.begin(); }

  /*! \brief End iterator (constant). */
  auto end() const { return _pcend; }

  /*! \brief Begin iterator (mutable).
   *
   * The iterator will point to a cut pointer.
   */
  auto begin() { return _pcuts.begin(); }

  /*! \brief End iterator (mutable). */
  auto end() { return _pend; }

  /*! \brief Number of cuts in the set. */
  auto size() const { return _pcend - _pcuts.begin(); }

  /*! \brief Returns reference to cut at index.
   *
   * This function does not return the cut pointer but dereferences it and
   * returns a reference.  The function does not check whether index is in the
   * valid range.
   *
   * \param index Index
   */
  auto const& operator[]( uint32_t index ) const { return *_pcuts[index]; }

  /*! \brief Returns the best cut, i.e., the first cut.
   */
  auto const& best() const { return *_pcuts[0]; }

  /*! \brief Updates the best cut.
   *
   * This method will set the cut at index `index` to be the best cut.  All
   * cuts before `index` will be moved one position higher.
   *
   * \param index Index of new best cut
   */
  void update_best( uint32_t index )
  {
    auto* best = _pcuts[index];
    for ( auto i = index; i > 0; --i )
    {
      _pcuts[i] = _pcuts[i - 1];
    }
    _pcuts[0] = best;
  }

  /*! \brief Resize the cut set, if it is too large.
   *
   * This method will resize the cut set to `size` only if the cut set has more
   * than `size` elements.  Otherwise, the size will remain the same.
   */
  void limit( uint32_t size )
  {
    if ( std::distance( _pcuts.begin(), _pend ) > static_cast<long>( size ) )
    {
      _pcend = _pend = _pcuts.begin() + size;
    }
  }

  /*! \brief Prints a cut set. */
  friend std::ostream& operator<<( std::ostream& os, emap_cut_set const& set )
  {
    for ( auto const& c : set )
    {
      os << *c << "\n";
    }
    return os;
  }

  /*! \brief Returns if the cut set contains already `cut`. */
  bool is_contained( CutType const& cut )
  {
    typename std::array<CutType*, MaxCuts>::iterator ipos = _pcuts.begin();

    while ( ipos != _pend )
    {
      if ( ( *ipos )->signature() == cut.signature() && std::equal( cut.begin(), cut.end(), ( *ipos )->begin() ) )
        return true;
      ++ipos;
    }

    return false;
  }

private:
  std::array<CutType, MaxCuts> _cuts;
  std::array<CutType*, MaxCuts> _pcuts;
  typename std::array<CutType*, MaxCuts>::const_iterator _pcend{ _pcuts.begin() };
  typename std::array<CutType*, MaxCuts>::iterator _pend{ _pcuts.begin() };
  uint32_t _set_limit{ MaxCuts };
};
#pragma endregion

#pragma region Hashing
template<uint32_t max_multioutput_cut_size>
struct emap_triple_hash
{
  inline uint64_t operator()( const std::array<uint32_t, max_multioutput_cut_size>& p ) const
  {
    uint64_t seed = hash_block( p[0] );

    for ( uint32_t i = 1; i < max_multioutput_cut_size; ++i )
    {
      hash_combine( seed, hash_block( p[i] ) );
    }

    return seed;
  }
};
#pragma endregion

template<unsigned NInputs>
struct best_gate_emap
{
  supergate<NInputs> const* gate;
  double arrival;
  float area;
  float flow;
  unsigned phase : 16;
  unsigned cut : 12;
  unsigned size : 4;
};

template<unsigned NInputs>
struct node_match_emap
{
  /* best gate match for positive and negative output phases */
  supergate<NInputs> const* best_gate[2];
  /* alternative best gate for positibe and negative output phase */
  best_gate_emap<NInputs> best_alternative[2];
  /* fanin pin phases for both output phases */
  uint16_t phase[2];
  /* best cut index for both phases */
  uint16_t best_cut[2];
  /* node is mapped using only one phase */
  bool same_match;
  /* node is mapped to a multi-output gate */
  bool multioutput_match[2];

  /* arrival time at node output */
  double arrival[2];
  /* required time at node output */
  double required[2];
  /* area of the best matches */
  float area[2];

  /* number of references in the cover 0: pos, 1: neg */
  uint32_t map_refs[2];
  /* references estimation */
  float est_refs[2];
  /* area flow */
  float flows[2];
};

template<class Ntk, unsigned CutSize, unsigned NInputs, classification_type Configuration>
class emap_impl
{
private:
  union multi_match_data
  {
    uint64_t data{ 0 };
    struct
    {
      uint64_t in_tfi : 1;
      uint64_t cut_index : 31;
      uint64_t node_index : 32;
    };
  };
  union multioutput_info
  {
    uint32_t data;
    struct
    {
      unsigned index : 29;
      unsigned lowest_index : 1;
      unsigned highest_index : 1;
      unsigned has_info : 1;
    };
  };

public:
  static constexpr float epsilon = 0.0005;
  static constexpr uint32_t max_cut_num = 20;
  using cut_t = cut<CutSize, cut_enumeration_emap_cut<NInputs>>;
  using cut_set_t = emap_cut_set<cut_t, max_cut_num>;
  using cut_merge_t = typename std::array<cut_set_t*, Ntk::max_fanin_size + 1>;
  using fanin_cut_t = typename std::array<cut_t const*, Ntk::max_fanin_size>;
  using support_t = typename std::array<uint8_t, CutSize>;
  using TT = kitty::static_truth_table<6>;
  using truth_compute_t = typename std::array<TT, CutSize>;
  using node_match_t = std::vector<node_match_emap<NInputs>>;
  using klut_map = std::unordered_map<uint32_t, std::array<signal<klut_network>, 2>>;
  using block_map = std::unordered_map<uint32_t, std::array<signal<block_network>, 2>>;

  static constexpr uint32_t max_multioutput_cut_size = 3;
  static constexpr uint32_t max_multioutput_output_size = 2;
  using multi_cuts_t = fast_network_cuts<Ntk, max_multioutput_cut_size, true, cut_enumeration_emap_multi_cut>;
  using multi_cut_t = typename multi_cuts_t::cut_t;
  using multi_leaves_set_t = std::array<uint32_t, max_multioutput_cut_size>;
  using multi_output_set_t = std::vector<multi_match_data>;
  using multi_hash_t = phmap::flat_hash_map<multi_leaves_set_t, multi_output_set_t, emap_triple_hash<max_multioutput_cut_size>>;
  using multi_match_t = std::array<multi_match_data, max_multioutput_output_size>;
  using multi_cut_set_t = std::vector<std::array<cut_t, max_multioutput_output_size>>;
  using multi_single_matches_t = std::vector<multi_match_t>;
  using multi_matches_t = std::vector<std::vector<multi_match_t>>;

  // Expose the actual fixed per-node tuple payload for external preallocation admission.
  // Variable multi-output matches and optional stores are deliberately separate.
  static constexpr std::size_t tuple_storage_bytes_per_node()
  {
    return sizeof( multioutput_info );
  }

  using clock = typename std::chrono::steady_clock;
  using time_point = typename clock::time_point;

public:
  explicit emap_impl( Ntk const& ntk, tech_library<NInputs, Configuration> const& library, emap_params const& ps, emap_stats& st )
      : ntk( ntk ),
        library( library ),
        ps( ps ),
        st( st ),
        node_match( ntk.size() ),
        node_tuple_match( ntk.size() ),
        switch_activity( ps.eswp_rounds ? switching_activity( ntk, ps.switching_activity_patterns ) : std::vector<float>( 0 ) ),
        cuts( ntk.size() )
  {
    if ( ps.wire_cap_for_fanout && ( !ps.electrical_model || ps.map_multioutput || ps.cover_buffer ) )
      throw std::invalid_argument( "wire capacitance callback requires electrical single-output mapping without cover-buffer" );
    std::memset( node_tuple_match.data(), 0, sizeof( multioutput_info ) * ntk.size() );
    std::tie( lib_inv_area, lib_inv_delay, lib_inv_id ) = library.get_inverter_info();
    std::tie( lib_buf_area, lib_buf_delay, lib_buf_id ) = library.get_buffer_info();
    std::tie( lib_inv_cap, lib_inv_slope ) = library.get_inverter_electrical();
    curve_budget = std::min( ps.curve_points, max_curve_points );
    if ( curve_budget > 0 )
      curves.resize( ntk.size() );
    init_load_ladder();
    buffer_stage = ps.cover_buffer && ps.electrical_model && library.has_buffer_drives();
    if ( buffer_stage )
      node_buffers.resize( ntk.size() );
    tmp_visited.reserve( 100 );
  }

  explicit emap_impl( Ntk const& ntk, tech_library<NInputs, Configuration> const& library, std::vector<float> const& switch_activity, emap_params const& ps, emap_stats& st )
      : ntk( ntk ),
        library( library ),
        ps( ps ),
        st( st ),
        node_match( ntk.size() ),
        node_tuple_match( ntk.size() ),
        switch_activity( switch_activity ),
        cuts( ntk.size() )
  {
    if ( ps.wire_cap_for_fanout && ( !ps.electrical_model || ps.map_multioutput || ps.cover_buffer ) )
      throw std::invalid_argument( "wire capacitance callback requires electrical single-output mapping without cover-buffer" );
    std::memset( node_tuple_match.data(), 0, sizeof( multioutput_info ) * ntk.size() );
    std::tie( lib_inv_area, lib_inv_delay, lib_inv_id ) = library.get_inverter_info();
    std::tie( lib_buf_area, lib_buf_delay, lib_buf_id ) = library.get_buffer_info();
    std::tie( lib_inv_cap, lib_inv_slope ) = library.get_inverter_electrical();
    curve_budget = std::min( ps.curve_points, max_curve_points );
    if ( curve_budget > 0 )
      curves.resize( ntk.size() );
    init_load_ladder();
    buffer_stage = ps.cover_buffer && ps.electrical_model && library.has_buffer_drives();
    if ( buffer_stage )
      node_buffers.resize( ntk.size() );
    tmp_visited.reserve( 100 );
  }

  cell_view<block_network> run_block()
  {
    time_begin = clock::now();

    auto [res, old2new] = initialize_block_network();

    /* multi-output initialization */
    if ( ps.map_multioutput && ps.matching_mode != emap_params::structural )
    {
      compute_multioutput_match();
    }

    /* compute and save topological order */
    init_topo_order();

    /* init arrival time */
    if ( !init_arrivals() )
      return res;

    /* search for large matches */
    if ( ps.matching_mode == emap_params::structural || CutSize > 6 )
    {
      if ( !compute_struct_match() )
      {
        return res;
      }
    }

    /* compute cuts, matches, and initial mapping */
    if ( !ps.area_oriented_mapping )
    {
      if ( !compute_mapping_match<false>() )
      {
        return res;
      }
    }
    else
    {
      if ( !compute_mapping_match<true>() )
      {
        return res;
      }
    }

    /* run area recovery */
    if ( !improve_mapping() )
      return res;

    /* insert buffers for POs driven by PIs */
    insert_buffers();

    /* generate the output network */
    finalize_cover_block( res, old2new );
    st.time_total = ( clock::now() - time_begin );

    return res;
  }

  binding_view<klut_network> run_klut()
  {
    time_begin = clock::now();

    auto [res, old2new] = initialize_map_network();

    /* multi-output initialization */
    if ( ps.map_multioutput && ps.matching_mode != emap_params::structural )
    {
      compute_multioutput_match();
    }

    /* compute and save topological order */
    init_topo_order();

    /* init arrival time */
    if ( !init_arrivals() )
      return res;

    /* search for large matches */
    if ( ps.matching_mode == emap_params::structural || CutSize > 6 )
    {
      if ( !compute_struct_match() )
      {
        return res;
      }
    }

    /* compute cuts, matches, and initial mapping */
    if ( !ps.area_oriented_mapping )
    {
      if ( !compute_mapping_match<false>() )
      {
        return res;
      }
    }
    else
    {
      if ( !compute_mapping_match<true>() )
      {
        return res;
      }
    }

    /* run area recovery */
    if ( !improve_mapping() )
      return res;

    /* insert buffers for POs driven by PIs */
    insert_buffers();

    /* generate the output network */
    finalize_cover( res, old2new );
    st.time_total = ( clock::now() - time_begin );

    return res;
  }

  binding_view<klut_network> run_node_map()
  {
    time_begin = clock::now();

    auto [res, old2new] = initialize_map_network();

    /* [i] multi-output support is currently not implemented */

    /* compute and save topological order */
    init_topo_order();

    /* init arrival time */
    if ( !init_arrivals() )
      return res;

    /* compute cuts, matches, and initial mapping */
    if ( !ps.area_oriented_mapping )
    {
      if ( !compute_mapping_match_node<false>() )
      {
        return res;
      }
    }
    else
    {
      if ( !compute_mapping_match_node<true>() )
      {
        return res;
      }
    }

    /* run area recovery */
    if ( !improve_mapping() )
      return res;

    /* insert buffers for POs driven by PIs */
    insert_buffers();

    /* generate the output network */
    finalize_cover( res, old2new );
    st.time_total = ( clock::now() - time_begin );

    return res;
  }

private:
  bool improve_mapping()
  {
    /* compute mapping using global area flow */
    uint32_t i = 0;
    while ( i++ < ps.area_flow_rounds )
    {
      if ( !compute_mapping<true>() )
      {
        return false;
      }
    }

    /* compute mapping using exact area */
    i = 0;
    compute_required_time( true );
    while ( i++ < ps.ela_rounds )
    {
      if ( !compute_mapping_exact_reversed<false>() )
      {
        return false;
      }
    }

    /* compute mapping using exact switching activity estimation */
    i = 0;
    while ( i++ < ps.eswp_rounds )
    {
      if ( !compute_mapping_exact_reversed<true>() )
      {
        return false;
      }
    }

    /* frontier through the covering DP (opt-in; see emap_params::frontier_points) */
    if ( ps.frontier_points > 0 && !ps.area_oriented_mapping )
    {
      uint32_t const rounds = std::max( 1u, std::min( ps.frontier_rounds, 8u ) );
      double first_area = 0.0, first_delay = 0.0;
      uint32_t cut_changes = 0, gate_changes = 0;
      for ( uint32_t r = 0; r < rounds; ++r )
      {
        if ( !compute_mapping_frontier() )
          return false;
        if ( r == 0 )
        {
          first_area = st.frontier_area_before;
          first_delay = st.frontier_delay_before;
        }
        if ( !st.frontier_ran || !st.frontier_kept )
          break;
        ++st.frontier_rounds_kept;
        cut_changes += st.frontier_cut_changes;
        gate_changes += st.frontier_gate_changes;
      }
      /* report the whole loop: from the cover it entered with to the one it leaves */
      if ( st.frontier_rounds_kept > 0 )
      {
        st.frontier_kept = true;
        st.frontier_area_after = area;
        st.frontier_delay_after = frontier_current_delay();
        st.frontier_cut_changes = cut_changes;
        st.frontier_gate_changes = gate_changes;
      }
      st.frontier_area_before = first_area;
      st.frontier_delay_before = first_delay;
    }

    return true;
  }

#pragma region Core
  template<bool DO_AREA>
  bool compute_mapping_match()
  {
    bool warning_box = false;

    for ( auto const& n : topo_order )
    {
      auto const index = ntk.node_to_index( n );

      if ( !compute_matches_node<DO_AREA>( n, warning_box ) )
      {
        continue;
      }

      /* load multi-output cuts and data */
      if ( ps.map_multioutput && node_tuple_match[index].has_info )
      {
        match_multi_add_cuts( n );
      }

      /* match positive phase */
      match_phase<DO_AREA>( n, 0u );

      /* match negative phase */
      match_phase<DO_AREA>( n, 1u );

      /* try to drop one phase */
      match_drop_phase<DO_AREA, false>( n );

      /* select alternative matches to use */
      select_alternatives<DO_AREA>( n );

      /* try multi-output matches */
      if constexpr ( DO_AREA )
      {
        if ( ps.map_multioutput && node_tuple_match[index].highest_index )
        {
          if ( match_multioutput<DO_AREA>( n ) )
            multi_node_update<DO_AREA>( n );
        }
      }

      /* this node's phases have settled: refresh the curves its consumers read (see the same
       * call in compute_mapping) */
      update_node_curves( index );
    }

    double area_old = area;
    bool success = set_mapping_refs_and_req<DO_AREA, false>();

    if ( warning_box )
    {
      std::cerr << "[i] MAP WARNING: not mapped don't touch gates are treated as sequential black boxes\n";
    }

    /* round stats */
    if ( ps.verbose )
    {
      std::stringstream stats{};
      float area_gain = 0.0f;

      if ( iteration != 1 )
        area_gain = float( ( area_old - area ) / area_old * 100 );

      if constexpr ( DO_AREA )
      {
        stats << fmt::format( "[i] AreaFlow : Delay = {:>12.2f}  Area = {:>12.2f}  Gain = {:>5.2f} %  Inverters = {:>5}  Time = {:>5.2f}\n", delay, area, area_gain, inv, to_seconds( clock::now() - time_begin ) );
      }
      else
      {
        stats << fmt::format( "[i] Delay    : Delay = {:>12.2f}  Area = {:>12.2f}  Gain = {:>5.2f} %  Inverters = {:>5}  Time = {:>5.2f}\n", delay, area, area_gain, inv, to_seconds( clock::now() - time_begin ) );
      }
      st.round_stats.push_back( stats.str() );
    }

    return success;
  }

  template<bool DO_AREA>
  inline bool compute_matches_node( node<Ntk> const& n, bool& warning_box )
  {
    auto const index = ntk.node_to_index( n );
    auto& node_data = node_match[index];

    node_data.est_refs[0] = node_data.est_refs[1] = static_cast<double>( ntk.fanout_size( n ) );
    node_data.map_refs[0] = node_data.map_refs[1] = 0;
    node_data.required[0] = node_data.required[1] = std::numeric_limits<float>::max();

    if ( ntk.is_constant( n ) )
    {
      /* all terminals have flow 0.0 */
      node_data.flows[0] = node_data.flows[1] = 0.0f;
      node_data.best_alternative[0].flow = node_data.best_alternative[1].flow = 0.0f;
      node_data.arrival[0] = node_data.arrival[1] = 0.0f;
      node_data.best_alternative[0].arrival = node_data.best_alternative[1].arrival = 0.0f;
      /* skip if cuts have been computed before */
      if ( cuts[index].size() == 0 )
      {
        add_zero_cut( index );
        match_constants( index );
      }
      return false;
    }
    else if ( ntk.is_pi( n ) )
    {
      node_data.flows[0] = 0.0f;
      node_data.best_alternative[0].flow = 0.0f;
      /* PIs have the negative phase implemented with an inverter */
      node_data.flows[1] = lib_inv_area / node_data.est_refs[1];
      node_data.best_alternative[1].flow = lib_inv_area / node_data.est_refs[1];
      /* skip if cuts have been computed before */
      if ( cuts[index].size() == 0 )
      {
        add_unit_cut( index );
      }
      return false;
    }

    if ( ps.matching_mode == emap_params::structural )
      return true;

    /* don't touch box */
    if constexpr ( has_is_dont_touch_v<Ntk> )
    {
      if ( ntk.is_dont_touch( n ) )
      {
        warning_box |= initialize_box( n );
        return false;
      }
    }

    /* compute cuts for node */
    if constexpr ( Ntk::min_fanin_size == 2 && Ntk::max_fanin_size == 2 )
    {
      merge_cuts2<DO_AREA>( n );
    }
    else
    {
      merge_cuts<DO_AREA>( n );
    }

    return true;
  }

  template<bool DO_AREA>
  void merge_cuts2( node<Ntk> const& n )
  {
    static constexpr uint32_t max_cut_size = CutSize > 6 ? 6 : CutSize;

    auto index = ntk.node_to_index( n );
    auto& node_data = node_match[index];
    emap_cut_sort_type sort = emap_cut_sort_type::AREA;

    /* compute cuts */
    const auto fanin = 2;
    ntk.foreach_fanin( ntk.index_to_node( index ), [this]( auto child, auto i ) {
      lcuts[i] = &cuts[ntk.node_to_index( ntk.get_node( child ) )];
    } );
    lcuts[2] = &cuts[index];
    auto& rcuts = *lcuts[fanin];

    /* move pre-computed structural cuts to a temporary cutset */
    bool reinsert_cuts = false;
    if ( rcuts.size() )
    {
      temp_cuts.clear();
      for ( auto& cut : rcuts )
      {
        if ( ( *cut )->ignore )
          continue;
        recompute_cut_data( *cut, n );
        temp_cuts.simple_insert( *cut );
        reinsert_cuts = true;
      }
      rcuts.clear();
    }

    /* set cut limit for run-time optimization*/
    rcuts.set_cut_limit( ps.cut_enumeration_ps.cut_limit );

    cut_t new_cut;
    new_cut->pattern_index = 0;
    fanin_cut_t vcuts;

    for ( auto const& c1 : *lcuts[0] )
    {
      /* skip cuts of pattern matching */
      if ( ( *c1 )->pattern_index > 1 )
        continue;
      vcuts[0] = c1;

      for ( auto const& c2 : *lcuts[1] )
      {
        /* skip cuts of pattern matching */
        if ( ( *c2 )->pattern_index > 1 )
          continue;

        if ( !c1->merge( *c2, new_cut, max_cut_size ) )
        {
          continue;
        }

        if ( ps.remove_dominated_cuts && rcuts.is_dominated( new_cut ) )
        {
          continue;
        }

        /* compute function */
        vcuts[1] = c2;
        compute_truth_table( index, vcuts, fanin, new_cut );

        /* match cut and compute data */
        compute_cut_data( new_cut, n );

        if ( ps.remove_dominated_cuts )
          rcuts.insert( new_cut, false, sort );
        else
          rcuts.simple_insert( new_cut, sort );
      }
    }

    if ( reinsert_cuts )
    {
      for ( auto const& cut : temp_cuts )
      {
        rcuts.simple_insert( *cut, sort );
      }
    }

    cuts_total += rcuts.size();

    /* limit the maximum number of cuts */
    rcuts.limit( ps.cut_enumeration_ps.cut_limit );

    /* add trivial cut */
    if ( rcuts.size() > 1 || ( *rcuts.begin() )->size() > 1 )
    {
      add_unit_cut( index );
    }
  }

  template<bool DO_AREA>
  void merge_cuts( node<Ntk> const& n )
  {
    static constexpr uint32_t max_cut_size = CutSize > 6 ? 6 : CutSize;

    auto index = ntk.node_to_index( n );
    auto& node_data = node_match[index];
    emap_cut_sort_type sort = emap_cut_sort_type::AREA;
    cut_t best_cut;

    /* compute cuts */
    std::vector<uint32_t> cut_sizes;
    ntk.foreach_fanin( ntk.index_to_node( index ), [this, &cut_sizes]( auto child, auto i ) {
      lcuts[i] = &cuts[ntk.node_to_index( ntk.get_node( child ) )];
      cut_sizes.push_back( static_cast<uint32_t>( lcuts[i]->size() ) );
    } );
    const auto fanin = cut_sizes.size();
    lcuts[fanin] = &cuts[index];
    auto& rcuts = *lcuts[fanin];

    /* set cut limit for run-time optimization*/
    rcuts.set_cut_limit( ps.cut_enumeration_ps.cut_limit );
    fanin_cut_t vcuts;

    if ( fanin > 1 && fanin <= ps.cut_enumeration_ps.fanin_limit )
    {
      cut_t new_cut, tmp_cut;

      foreach_mixed_radix_tuple( cut_sizes.begin(), cut_sizes.end(), [&]( auto begin, auto end ) {
        auto it = vcuts.begin();
        auto i = 0u;
        while ( begin != end )
        {
          *it++ = &( ( *lcuts[i++] )[*begin++] );
        }

        if ( !vcuts[0]->merge( *vcuts[1], new_cut, max_cut_size ) )
        {
          return true; /* continue */
        }

        for ( i = 2; i < fanin; ++i )
        {
          tmp_cut = new_cut;
          if ( !vcuts[i]->merge( tmp_cut, new_cut, max_cut_size ) )
          {
            return true; /* continue */
          }
        }

        if ( ps.remove_dominated_cuts && rcuts.is_dominated( new_cut ) )
        {
          return true; /* continue */
        }

        compute_truth_table( index, vcuts, fanin, new_cut );

        /* match cut and compute data */
        compute_cut_data( new_cut, n );

        if ( ps.remove_dominated_cuts )
          rcuts.insert( new_cut, false, sort );
        else
          rcuts.simple_insert( new_cut, sort );

        return true;
      } );

      /* limit the maximum number of cuts */
      rcuts.limit( ps.cut_enumeration_ps.cut_limit );
    }
    else if ( fanin == 1 )
    {
      for ( auto const& cut : *lcuts[0] )
      {
        cut_t new_cut = *cut;
        vcuts[0] = cut;

        compute_truth_table( index, vcuts, fanin, new_cut );

        /* match cut and compute data */
        compute_cut_data( new_cut, n );

        if ( ps.remove_dominated_cuts )
          rcuts.insert( new_cut, false, sort );
        else
          rcuts.simple_insert( new_cut, sort );
      }

      /* limit the maximum number of cuts */
      rcuts.limit( ps.cut_enumeration_ps.cut_limit );
    }

    cuts_total += rcuts.size();

    add_unit_cut( index );
  }

  bool compute_struct_match()
  {
    if ( ps.matching_mode == emap_params::boolean )
      return true;

    /* compatible only with AIGs */
    if constexpr ( !is_aig_network_type_v<Ntk> )
    {
      if ( ps.matching_mode == emap_params::structural )
      {
        std::cerr << "[e] MAP ERROR: structural library works only with AIGs\n";
        return false;
      }
      return true;
    }

    /* no large gates identified */
    if ( library.num_structural_gates() == 0 )
    {
      if ( ps.matching_mode == emap_params::structural )
      {
        std::cerr << "[e] MAP ERROR: structural library is empty\n";
        return false;
      }
      return true;
    }

    bool warning_box = false;
    for ( auto const& n : topo_order )
    {
      auto const index = ntk.node_to_index( n );
      auto& node_data = node_match[index];

      if ( ntk.is_constant( n ) )
      {
        add_zero_cut( index );
        match_constants( index );
        continue;
      }
      else if ( ntk.is_pi( n ) )
      {
        add_unit_cut( index );
        continue;
      }

      /* don't touch box */
      if constexpr ( has_is_dont_touch_v<Ntk> )
      {
        if ( ntk.is_dont_touch( n ) )
        {
          add_unit_cut( index );
          continue;
        }
      }

      /* compute cuts for node */
      merge_cuts_structural( n );
    }

    if ( warning_box )
    {
      std::cerr << "[i] MAP WARNING: not mapped don't touch gates are treated as sequential black boxes\n";
    }

    /* round stats */
    if ( ps.verbose )
    {
      st.round_stats.push_back( fmt::format( "[i] SCuts    : Cuts  = {:>12d}  Time = {:>12.2f}\n", cuts_total, to_seconds( clock::now() - time_begin ) ) );
    }

    return true;
  }

  void merge_cuts_structural( node<Ntk> const& n )
  {
    auto index = ntk.node_to_index( n );
    auto& node_data = node_match[index];
    emap_cut_sort_type sort = emap_cut_sort_type::AREA;

    /* compute cuts */
    const auto fanin = 2;
    std::array<uint32_t, 2> children_phase;
    ntk.foreach_fanin( ntk.index_to_node( index ), [&]( auto child, auto i ) {
      lcuts[i] = &cuts[ntk.node_to_index( ntk.get_node( child ) )];
      children_phase[i] = ntk.is_complemented( child ) ? 1 : 0;
    } );
    lcuts[2] = &cuts[index];
    auto& rcuts = *lcuts[fanin];

    /* set cut limit for run-time optimization*/
    rcuts.set_cut_limit( ps.cut_enumeration_ps.cut_limit );

    cut_t new_cut;
    std::vector<cut_t const*> vcuts( fanin );

    for ( auto const& c1 : *lcuts[0] )
    {
      for ( auto const& c2 : *lcuts[1] )
      {
        /* filter large cuts */
        if ( c1->size() + c2->size() > CutSize || c1->size() + c2->size() > NInputs )
          continue;
        /* filter cuts involving constants */
        if ( ( *c1 )->pattern_index == 0 || ( *c2 )->pattern_index == 0 )
          continue;

        vcuts[0] = c1;
        vcuts[1] = c2;
        uint32_t pattern_id1 = ( ( *c1 )->pattern_index << 1 ) | children_phase[0];
        uint32_t pattern_id2 = ( ( *c2 )->pattern_index << 1 ) | children_phase[1];
        if ( pattern_id1 > pattern_id2 )
        {
          std::swap( vcuts[0], vcuts[1] );
          std::swap( pattern_id1, pattern_id2 );
        }

        uint32_t new_pattern = library.get_pattern_id( pattern_id1, pattern_id2 );

        /* pattern not matched */
        if ( new_pattern == UINT32_MAX )
          continue;

        create_structural_cut( new_cut, vcuts, new_pattern, pattern_id1, pattern_id2 );

        if ( ps.remove_dominated_cuts && rcuts.is_dominated( new_cut ) )
          continue;

        /* match cut and compute data */
        compute_cut_data_structural( new_cut, n );

        if ( ps.remove_dominated_cuts )
          rcuts.insert( new_cut, false, sort );
        else
          rcuts.simple_insert( new_cut, sort );
      }
    }

    cuts_total += rcuts.size();

    /* limit the maximum number of cuts */
    rcuts.limit( ps.cut_enumeration_ps.cut_limit );

    /* add trivial cut */
    if ( rcuts.size() > 1 || ( *rcuts.begin() )->size() > 1 )
    {
      add_unit_cut( index );
    }
  }

  template<bool DO_AREA>
  bool compute_mapping_match_node()
  {
    for ( auto const& n : topo_order )
    {
      auto const index = ntk.node_to_index( n );
      auto& node_data = node_match[index];

      node_data.best_gate[0] = node_data.best_gate[1] = nullptr;
      node_data.same_match = 0;
      node_data.multioutput_match[0] = node_data.multioutput_match[1] = false;
      node_data.required[0] = node_data.required[1] = std::numeric_limits<float>::max();
      node_data.map_refs[0] = node_data.map_refs[1] = 0;
      node_data.est_refs[0] = node_data.est_refs[1] = static_cast<float>( ntk.fanout_size( n ) );

      if ( ntk.is_constant( n ) )
      {
        /* all terminals have flow 0 */
        node_data.flows[0] = node_data.flows[1] = 0.0f;
        node_data.arrival[0] = node_data.arrival[1] = 0.0f;
        add_zero_cut( index );
        match_constants( index );
        continue;
      }
      else if ( ntk.is_pi( n ) )
      {
        /* all terminals have flow 0 */
        node_data.flows[0] = 0.0f;
        /* PIs have the negative phase implemented with an inverter */
        node_data.flows[1] = lib_inv_area / node_data.est_refs[1];
        add_unit_cut( index );
        continue;
      }

      /* compute the node mapping */
      add_node_cut<DO_AREA>( n );

      /* match positive phase */
      match_phase<DO_AREA>( n, 0u );

      /* match negative phase */
      match_phase<DO_AREA>( n, 1u );

      /* try to drop one phase */
      match_drop_phase<DO_AREA, false>( n );

      /* select alternative matches to use */
      select_alternatives<DO_AREA>( n );
    }
    double area_old = area;
    bool success = set_mapping_refs_and_req<DO_AREA, false>();

    /* round stats */
    if ( ps.verbose )
    {
      std::stringstream stats{};
      float area_gain = 0.0f;

      if ( iteration != 1 )
        area_gain = float( ( area_old - area ) / area_old * 100 );

      if constexpr ( DO_AREA )
      {
        stats << fmt::format( "[i] AreaFlow : Delay = {:>12.2f}  Area = {:>12.2f}  Gain = {:>5.2f} %  Inverters = {:>5}  Time = {:>5.2f}\n", delay, area, area_gain, inv, to_seconds( clock::now() - time_begin ) );
      }
      else
      {
        stats << fmt::format( "[i] Delay    : Delay = {:>12.2f}  Area = {:>12.2f}  Gain = {:>5.2f} %  Inverters = {:>5}  Time = {:>5.2f}\n", delay, area, area_gain, inv, to_seconds( clock::now() - time_begin ) );
      }
      st.round_stats.push_back( stats.str() );
    }

    return success;
  }

  template<bool DO_AREA>
  void add_node_cut( node<Ntk> const& n )
  {
    auto index = ntk.node_to_index( n );
    auto& node_data = node_match[index];
    auto& rcuts = &cuts[index];

    std::vector<uint32_t> fanin_indexes;
    fanin_indexes.reserve( Ntk::max_fanin_size );

    ntk.foreach_fanin( n, [&]( auto const& f ) {
      fanin_indexes.push_back( ntk.node_to_index( ntk.get_node( f ) ) );
    } );

    assert( fanin_indexes.size() <= CutSize );

    cut_t new_cut = rcuts.add_cut( fanin_indexes.begin(), fanin_indexes.end() );
    new_cut->function = kitty::extend_to<6>( ntk.node_function( n ) );

    /* match cut and compute data */
    compute_cut_data( new_cut, n );

    ++cuts_total;
  }

  template<bool DO_AREA>
  bool compute_mapping()
  {
    for ( auto const& n : topo_order )
    {
      uint32_t index = ntk.node_to_index( n );

      /* reset mapping */
      node_match[index].map_refs[0] = node_match[index].map_refs[1] = 0u;

      if ( ntk.is_constant( n ) )
        continue;
      if ( ntk.is_pi( n ) )
      {
        node_match[index].flows[1] = lib_inv_area / node_match[index].est_refs[1];
        node_match[index].best_alternative[1].flow = lib_inv_area / node_match[index].est_refs[1];
        continue;
      }

      /* don't touch box */
      if constexpr ( has_is_dont_touch_v<Ntk> )
      {
        if ( ntk.is_dont_touch( n ) )
        {
          if constexpr ( has_has_binding_v<Ntk> )
          {
            propagate_data_forward_white_box( n );
          }
          continue;
        }
      }

      /* match positive phase */
      match_phase<DO_AREA>( n, 0u );

      /* match negative phase */
      match_phase<DO_AREA>( n, 1u );

      /* try to drop one phase */
      match_drop_phase<DO_AREA, false>( n );

      /* try a multi-output match */
      if constexpr ( DO_AREA )
      {
        if ( ps.map_multioutput && node_tuple_match[index].highest_index )
        {
          bool multi_success = match_multioutput<DO_AREA>( n );
          if ( multi_success )
            multi_node_update<DO_AREA>( n );
        }
      }

      /* this node's phases have settled: refresh the arrival curves its consumers read, in the
       * same topological pass, so a parent later in this round prices itself against a leaf
       * curve from THIS round exactly as it already does for scalar arrivals */
      update_node_curves( index );

      assert( node_match[index].arrival[0] < node_match[index].required[0] + epsilon );
      assert( node_match[index].arrival[1] < node_match[index].required[1] + epsilon );
    }

    double area_old = area;
    bool success = set_mapping_refs_and_req<DO_AREA, false>();

    /* round stats */
    if ( ps.verbose )
    {
      std::stringstream stats{};
      float area_gain = 0.0f;

      if ( iteration != 1 )
        area_gain = float( ( area_old - area ) / area_old * 100 );

      if constexpr ( DO_AREA )
      {
        stats << fmt::format( "[i] AreaFlow : Delay = {:>12.2f}  Area = {:>12.2f}  Gain = {:>5.2f} %  Inverters = {:>5}  Time = {:>5.2f}\n", delay, area, area_gain, inv, to_seconds( clock::now() - time_begin ) );
      }
      else
      {
        stats << fmt::format( "[i] Delay    : Delay = {:>12.2f}  Area = {:>12.2f}  Gain = {:>5.2f} %  Inverters = {:>5}  Time = {:>5.2f}\n", delay, area, area_gain, inv, to_seconds( clock::now() - time_begin ) );
      }
      st.round_stats.push_back( stats.str() );
    }

    return success;
  }

  template<bool SwitchActivity>
  bool compute_mapping_exact_reversed( bool record_stats = true )
  {
    for ( auto it = topo_order.rbegin(); it != topo_order.rend(); ++it )
    {
      if ( ntk.is_constant( *it ) || ntk.is_pi( *it ) )
        continue;

      const auto index = ntk.node_to_index( *it );
      auto& node_data = node_match[index];

      /* skip not mapped nodes */
      if ( !node_data.map_refs[0] && !node_data.map_refs[1] )
        continue;

      /* don't touch box */
      if constexpr ( has_is_dont_touch_v<Ntk> )
      {
        node<Ntk> n = ntk.index_to_node( index );
        if ( ntk.is_dont_touch( n ) )
        {
          if constexpr ( has_has_binding_v<Ntk> )
          {
            propagate_data_backward_white_box( n );
          }
          continue;
        }
      }

      /* recursively deselect the best cut shared between
       * the two phases if in use in the cover */
      uint8_t use_phase = node_data.best_gate[0] != nullptr ? 0 : 1;
      double old_required = -1;
      if ( node_data.same_match )
      {
        auto const& best_cut = cuts[index][node_data.best_cut[use_phase]];
        cut_deref<SwitchActivity>( best_cut, *it, use_phase );

        /* propagate required time over the output inverter if present */
        if ( node_data.map_refs[use_phase ^ 1] > 0 )
        {
          old_required = node_data.required[use_phase];
          node_data.required[use_phase] = std::min( node_data.required[use_phase], node_data.required[use_phase ^ 1] - inv_delay_at( index, use_phase ^ 1 ) );
        }
      }
      else if ( !node_data.map_refs[0] || !node_data.map_refs[1] )
      {
        use_phase = node_data.map_refs[0] ? 0 : 1;
        auto const& best_cut = cuts[index][node_data.best_cut[use_phase]];
        cut_deref<SwitchActivity>( best_cut, *it, use_phase );
        node_data.same_match = true;
      }

      /* match positive phase */
      match_phase_exact<SwitchActivity>( *it, 0u );

      /* match negative phase */
      match_phase_exact<SwitchActivity>( *it, 1u );

      /* restore required time */
      if ( old_required > 0 )
      {
        node_data.required[use_phase] = old_required;
      }

      /* try to drop one phase */
      match_drop_phase<true, true, SwitchActivity>( *it );

      /* try a multi-output match */ /* TODO: fix the required time*/
      if ( ps.map_multioutput && node_tuple_match[index].lowest_index )
      {
        bool mapped = match_multioutput_exact<SwitchActivity>( *it, true );

        /* propagate required time for the selected gates */
        if ( mapped )
        {
          match_multioutput_propagate_required( *it );
        }
        else
        {
          match_propagate_required( index );
        }
      }
      else
      {
        match_propagate_required( index );
      }
    }

    double area_old = area;

    propagate_arrival_times();

    /* round stats */
    if ( ps.verbose && record_stats )
    {
      float area_gain = float( ( area_old - area ) / area_old * 100 );
      std::stringstream stats{};
      if constexpr ( SwitchActivity )
        stats << fmt::format( "[i] Switching: Delay = {:>12.2f}  Area = {:>12.2f}  Gain = {:>5.2f} %  Inverters = {:>5}  Time = {:>5.2f}\n", delay, area, area_gain, inv, to_seconds( clock::now() - time_begin ) );
      else
        stats << fmt::format( "[i] Area Rev : Delay = {:>12.2f}  Area = {:>12.2f}  Gain = {:>5.2f} %  Inverters = {:>5}  Time = {:>5.2f}\n", delay, area, area_gain, inv, to_seconds( clock::now() - time_begin ) );
      st.round_stats.push_back( stats.str() );
    }

    return true;
  }

  /* --- electrical (load-aware) model ------------------------------------------------------ */

  /* Load-aware pin delay of supergate @p g implementing (index, phase): block + slope × the
   * load the current cover imposes on that output signal. Exactly g->tdelay[ctr] with the
   * model off. */
  inline double gate_delay( supergate<NInputs> const* g, uint32_t ctr, uint32_t index, uint8_t phase ) const
  {
    if ( !ps.electrical_model )
      return g->tdelay[ctr];
    return g->tdelay[ctr] + g->slope[ctr] * node_loads[index][phase];
  }

  /* Delay of a polarity inverter whose OUTPUT is (index, inv_phase)'s signal. */
  inline double inv_delay_at( uint32_t index, uint8_t inv_phase ) const
  {
    if ( !ps.electrical_model )
      return lib_inv_delay;
    return lib_inv_delay + lib_inv_slope * node_loads[index][inv_phase];
  }

  /* Can @p g legally drive (index, phase)'s current load? A gate without a declared limit
   * (max_load 0: no data / unconstrained) is always legal. Trivially true with the model off. */
  inline bool gate_load_legal( supergate<NInputs> const* g, uint32_t index, uint8_t phase ) const
  {
    if ( !ps.electrical_model || g->max_load <= 0.0f )
      return true;
    return node_loads[index][phase] <= g->max_load + epsilon;
  }

  /* Extra delay of the optional output buffer on (index, phase)'s signal (cover_buffer only):
   * block + slope × the SINK load the buffer carries. 0 when the stage is off or the identity
   * ("no buffer") was chosen — the zero-cost wire branch of the two-level covering decision. */
  inline double buffer_delay_at( uint32_t index, uint8_t phase ) const
  {
    if ( !buffer_stage || !node_buffers[index][phase].active )
      return 0.0;
    auto const& b = node_buffers[index][phase];
    return b.drive.block + b.drive.slope * b.sink_load;
  }

  /* Area of the optional output buffer on (index, phase)'s signal; 0 when none. */
  inline float buffer_area_at( uint32_t index, uint8_t phase ) const
  {
    if ( !buffer_stage || !node_buffers[index][phase].active )
      return 0.0f;
    return node_buffers[index][phase].drive.area;
  }

  /* The binding for a MATERIALIZED polarity inverter on (index, inv_phase)'s signal: under the
   * electrical model, the smallest inverter drive that legally carries that signal's real load
   * (tech_library::select_inverter) — a hardcoded smallest inverter on a high-fanout net is the
   * classic min-size drive violation. Load-blind: the plain smallest inverter. */
  inline uint32_t select_polarity_inverter( uint32_t index, uint8_t inv_phase ) const
  {
    if ( !ps.electrical_model )
      return lib_inv_id;
    return library.select_inverter( node_loads[index][inv_phase] );
  }

  double wire_cap_at( uint32_t fanout ) const
  {
    double const cap = ps.wire_cap_for_fanout( fanout );
    if ( !std::isfinite( cap ) || cap < 0.0 || cap > std::numeric_limits<float>::max() )
      throw std::invalid_argument( "wire capacitance callback must return finite nonnegative mapper-range total capacitance" );
    return cap;
  }

  float incremental_wire_cap( uint32_t index, uint8_t phase ) const
  {
    if ( !ps.wire_cap_for_fanout )
      return static_cast<float>( ps.wire_cap_per_fanout );
    uint32_t const count = node_fanouts[index][phase];
    if ( count == std::numeric_limits<uint32_t>::max() )
      throw std::overflow_error( "wire capacitance fanout count overflow" );
    // A new sink REPLACES C(F) with C(F+1). The delta can be negative for
    // legitimate nonmonotone Liberty tables; never clamp it to zero.
    return static_cast<float>( wire_cap_at( count + 1 ) - wire_cap_at( count ) );
  }

  /* --- load-indexed arrival curves (ADR-0050) ---------------------------------------------- */

  /* Build the ladder a node's arrival is sampled on and size the store. The anchors are library
   * facts: the smallest load any sink can present (plus the wire term, which every real edge
   * carries) and the largest a driver declares it may legally carry.
   *
   * The curves stay INACTIVE — every query then answers with the node's scalar arrival, exactly
   * as before this feature — in three cases the caller must be able to tell apart from a working
   * ladder, which is why `st.load_points` reports what was built rather than what was asked for:
   * the model is off, a single point was requested (one load point IS the scalar model), or the
   * library declares no drive limit anywhere and so has no ceiling to span. */
  void init_load_ladder()
  {
    if ( !ps.electrical_model || ps.load_points < 2 )
      return;
    float const wire = ps.wire_cap_for_fanout ? static_cast<float>( wire_cap_at( 1 ) )
                                               : static_cast<float>( ps.wire_cap_per_fanout );
    ladder = build_load_ladder( library.get_min_pin_cap() + wire, library.get_max_drive_load(),
                                ps.load_points );
    if ( ladder.size < 2 )
    {
      ladder = load_ladder{};
      return;
    }
    node_curves.resize( ntk.size() );
    st.load_points = ladder.size;
    st.load_curve_lo = ladder.lo();
    st.load_curve_hi = ladder.hi();
    st.load_curve_bytes = load_curve_store_bytes( ntk.size(), ladder.size, sizeof( double ) );
  }

  /* Arrival at (@p leaf, @p leaf_phase) once @p extra_cap — the input capacitance of the
   * candidate pin asking the question, plus its wire — is added to that signal's load.
   *
   * This is the whole point of the curve. The scalar model answers with the number the leaf
   * stored under its own load estimate, which is the same answer for every candidate, so a
   * gate with heavy input pins is never charged for the slowdown it causes upstream.
   *
   * The query load is the leaf's currently estimated TOTAL load plus this pin, which
   * over-states it by whatever this same edge contributed to the estimate last round — bounded
   * by one pin capacitance, always in the pessimistic direction, and identical in sign for
   * every candidate, so the ORDER the cover cares about is preserved. Subtracting the exact
   * previous share would need a per-edge charge record; that refinement is recorded in
   * ADR-0050 rather than guessed at here. */
  inline double leaf_arrival( uint32_t leaf, uint8_t leaf_phase, float extra_cap ) const
  {
    double const scalar = node_match[leaf].arrival[leaf_phase];
    if ( node_curves.empty() )
      return scalar;
    float const load = node_loads[leaf][leaf_phase] + extra_cap;
    if ( ps.wire_cap_for_fanout && ( !std::isfinite( load ) || load < 0.0f ) )
      throw std::invalid_argument( "nonlinear wire candidate load must remain finite nonnegative" );
    double const curved = interpolate_curve( ladder, node_curves[leaf][leaf_phase], load, scalar );
    /* Under the scalar wire estimate, query load is the leaf's own estimate PLUS a
     * pin, and arrival rises with load, so a
     * curve can only ever answer at or above the scalar. Taking the max enforces that as an
     * invariant rather than trusting it: a curve that has not been refreshed yet (a node the
     * current round has not reached, a phase no match settled on) then answers exactly what the
     * scalar model would, which is conservative. Answering BELOW the scalar would be optimism
     * the rest of the mapper reconciles against and cannot absorb — required times derive from
     * these arrivals. */
    if ( ps.wire_cap_for_fanout && extra_cap < 0.0f )
      return curved; // a table decrease can outweigh the new sink's pin capacitance
    return curved > scalar ? curved : scalar;
  }

  /* Same query for the pin @p ctr of candidate @p g. */
  inline double leaf_arrival_pin( uint32_t leaf, uint8_t leaf_phase, supergate<NInputs> const* g,
                                  uint32_t ctr ) const
  {
    if ( node_curves.empty() )
      return node_match[leaf].arrival[leaf_phase];
    return leaf_arrival( leaf, leaf_phase, g->cap[ctr] + incremental_wire_cap( leaf, leaf_phase ) );
  }

  /* How much later this pin makes its leaf arrive: the load penalty the candidate imposes,
   * 0 with the curves off.
   *
   * Required-time propagation must subtract this beside the gate delay. Arrival says
   * `A_parent = A_leaf + penalty + gate_delay`, so the leaf's deadline is
   * `required_parent - gate_delay - penalty` — anything looser and the two directions of the
   * same inequality would use different delay models, which is how a cover ends up asserting
   * that an arrival it just computed exceeds a required time it derived from the old one. */
  inline double leaf_load_penalty( uint32_t leaf, uint8_t leaf_phase, supergate<NInputs> const* g,
                                   uint32_t ctr ) const
  {
    if ( node_curves.empty() )
      return 0.0;
    return leaf_arrival_pin( leaf, leaf_phase, g, ctr ) - node_match[leaf].arrival[leaf_phase];
  }

  /* Seed every curve with the node's scalar arrival: a node that never gets a match (a PI, a
   * constant, a don't-touch box, an unmatched phase) then answers every query exactly as the
   * scalar model would, so activating the ladder cannot perturb a node the DP does not price. */
  void init_node_curves()
  {
    if ( node_curves.empty() )
      return;
    for ( uint32_t i = 0; i < node_curves.size(); ++i )
    {
      for ( uint8_t phase = 0; phase < 2; ++phase )
        node_curves[i][phase].fill( node_match[i].arrival[phase] );
    }

    /* a primary input's own arrival is fixed, but its NEGATED phase is a materialized inverter,
     * and an inverter's delay depends on what it drives like any other cell's — the one place a
     * load curve is not the mapper's to build later, since the DP never matches a PI */
    ntk.foreach_pi( [&]( auto const& n ) {
      uint32_t const index = ntk.node_to_index( n );
      for ( uint32_t k = 0; k < ladder.size; ++k )
        node_curves[index][1][k] = node_match[index].arrival[0] + lib_inv_delay + lib_inv_slope * ladder.points[k];
    } );
  }

  /* Rebuild (index)'s curves from the match it just committed: what this node delivers at each
   * sampled output load, given the cell it bound and the leaves it reads.
   *
   * The curve describes the BOUND match, not the best match available at each load. That keeps
   * it consistent with `node_data.arrival` — the number the rest of the mapper reconciles
   * against — instead of handing parents an optimistic envelope no committed cover realizes.
   * A leaf's stronger alternative still reaches this node, one round later, through the load
   * refresh the parent's own choice feeds (ADR-0050).
   *
   * Called in topological order right after the node's phases settle, so a parent always reads
   * a leaf curve from the current round, matching how scalar arrivals already flow. */
  void update_node_curves( uint32_t index )
  {
    if ( node_curves.empty() )
      return;

    auto const& node_data = node_match[index];
    uint8_t const use_phase = node_data.best_gate[0] != nullptr ? 0u : 1u;

    for ( uint8_t p = 0; p < 2; ++p )
    {
      uint8_t const phase = static_cast<uint8_t>( use_phase ^ p ); /* the implemented phase first */
      auto& curve = node_curves[index][phase];
      supergate<NInputs> const* g = node_data.best_gate[phase];

      if ( g == nullptr )
      {
        if ( p == 1 && node_data.best_gate[phase ^ 1] != nullptr )
        {
          /* One cell serving both polarities: this phase is a materialized inverter reading the
           * implemented signal, so its curve is that signal's arrival plus the inverter's own
           * block and slope at each sampled load — the same shape as `inv_delay_at`, with the
           * load axis exposed.
           *
           * The base is the implemented phase's SCALAR arrival, deliberately, not its curve at
           * the inverter's input capacitance: `update_node_loads` already charges that
           * capacitance onto the implemented signal, so it is inside that arrival. Querying the
           * curve here would bill the same pin twice. */
          double const in = node_match[index].arrival[phase ^ 1];
          for ( uint32_t k = 0; k < ladder.size; ++k )
            curve[k] = in + lib_inv_delay + lib_inv_slope * ladder.points[k];
        }
        else
        {
          curve.fill( node_data.arrival[phase] );
        }
        continue;
      }

      /* per-pin line: (leaf arrival at this pin's capacitance + block) + slope x load */
      std::array<double, NInputs> base{};
      std::array<float, NInputs> slope{};
      uint32_t pins = 0;
      for ( auto l : cuts[index][node_data.best_cut[phase]] )
      {
        if ( pins >= NInputs )
          break;
        uint8_t const leaf_phase = ( node_data.phase[phase] >> pins ) & 1;
        base[pins] = node_match[l].arrival[leaf_phase] + g->tdelay[pins];
        slope[pins] = g->slope[pins];
        ++pins;
      }

      double const buffered = buffer_delay_at( index, phase );
      for ( uint32_t k = 0; k < ladder.size; ++k )
      {
        double worst = 0.0;
        for ( uint32_t i = 0; i < pins; ++i )
          worst = std::max( worst, base[i] + slope[i] * ladder.points[k] );
        curve[k] = worst + buffered;
      }
    }
  }

  void init_node_loads()
  {
    if ( !ps.electrical_model )
      return;
    /* round-0 seed: fanout count × the library's median input capacitance — replaced by the
     * cover's actual pin capacitances after the first round (update_node_loads) */
    float const seed_cap = library.get_median_pin_cap() + static_cast<float>( ps.wire_cap_per_fanout );
    node_loads.assign( ntk.size(), { { 0.0f, 0.0f } } );
    if ( ps.wire_cap_for_fanout )
      node_fanouts.assign( ntk.size(), { { 0u, 0u } } );
    ntk.foreach_node( [&]( auto const& n ) {
      uint32_t const index = ntk.node_to_index( n );
      uint32_t const count = ntk.fanout_size( n ); // network includes virtual output refs
      float const load = ps.wire_cap_for_fanout
                           ? library.get_median_pin_cap() * count + static_cast<float>( wire_cap_at( count ) )
                           : seed_cap * count;
      node_loads[index][0] = node_loads[index][1] = load;
      if ( ps.wire_cap_for_fanout )
        node_fanouts[index][0] = node_fanouts[index][1] = count;
    } );
  }

  /* Recompute each signal's load from the cover's committed choices. The traversal mirrors
   * set_mapping_refs_and_req's counting walk EXACTLY — a cell is present when
   * `same_match || map_refs[use_phase] > 0` (with same_match a node's references can sit
   * entirely in the OTHER phase while the cell lives at use_phase, so pairing map_refs[p] with
   * best_gate[p] would charge nothing on such nodes and silently decay the whole estimate);
   * with same_match the other phase is served by a materialized inverter whose input
   * capacitance loads the implemented signal; a negated PI materializes an inverter reading
   * the PI's positive signal. */
  void update_node_loads()
  {
    if ( !ps.electrical_model )
      return;
    float const wire = ps.wire_cap_for_fanout ? 0.0f : static_cast<float>( ps.wire_cap_per_fanout );
    for ( auto& l : node_loads )
      l[0] = l[1] = 0.0f;
    if ( ps.wire_cap_for_fanout )
      for ( auto& count : node_fanouts )
        count[0] = count[1] = 0u;

    auto charge_cut = [&]( uint32_t index, uint8_t phase ) {
      auto const& node_data = node_match[index];
      supergate<NInputs> const* sg = node_data.best_gate[phase];
      auto const& best_cut = cuts[index][node_data.best_cut[phase]];
      uint32_t ctr = 0u;
      for ( auto leaf : best_cut )
      {
        if ( ctr >= NInputs )
          break;
        uint8_t const leaf_phase = ( node_data.phase[phase] >> ctr ) & 1;
        node_loads[leaf][leaf_phase] += sg->cap[ctr] + wire;
        if ( ps.wire_cap_for_fanout )
          ++node_fanouts[leaf][leaf_phase];
        ++ctr;
      }
    };

    for ( auto const& n : topo_order )
    {
      uint32_t const index = ntk.node_to_index( n );
      auto const& node_data = node_match[index];

      if ( ntk.is_constant( n ) )
        continue;
      if ( ntk.is_pi( n ) )
      {
        if ( node_data.map_refs[1] > 0 )
        {
          node_loads[index][0] += lib_inv_cap + wire; /* negated PI: an inverter reads the PI */
          if ( ps.wire_cap_for_fanout )
            ++node_fanouts[index][0];
        }
        continue;
      }
      if constexpr ( has_is_dont_touch_v<Ntk> )
      {
        if ( ntk.is_dont_touch( n ) )
          continue; /* box internals unknown: charge nothing (conservative) */
      }
      if ( !node_data.map_refs[0] && !node_data.map_refs[1] )
        continue;

      uint8_t const use_phase = node_data.best_gate[0] == nullptr ? 1u : 0u;
      if ( node_data.best_gate[use_phase] == nullptr )
        continue;

      if ( node_data.same_match || node_data.map_refs[use_phase] > 0 )
      {
        charge_cut( index, use_phase );
        /* one cell serving both polarities: the materialized inverter's input loads the
         * implemented signal */
        if ( node_data.same_match && node_data.map_refs[use_phase ^ 1] > 0 )
        {
          node_loads[index][use_phase] += lib_inv_cap + wire;
          if ( ps.wire_cap_for_fanout )
            ++node_fanouts[index][use_phase];
        }
      }
      if ( !node_data.same_match && node_data.map_refs[use_phase ^ 1] > 0 &&
           node_data.best_gate[use_phase ^ 1] != nullptr )
      {
        charge_cut( index, use_phase ^ 1 );
      }
    }

    if ( ps.wire_cap_for_fanout )
    {
      ntk.foreach_po( [&]( auto const& signal ) {
        uint32_t const index = ntk.node_to_index( ntk.get_node( signal ) );
        uint8_t const phase = ntk.is_complemented( signal ) ? 1u : 0u;
        ++node_fanouts[index][phase];
      } );
      for ( uint32_t index = 0; index < node_loads.size(); ++index )
        for ( uint8_t phase = 0; phase < 2; ++phase )
          node_loads[index][phase] += static_cast<float>( wire_cap_at( node_fanouts[index][phase] ) );
    }

    select_output_buffers();
  }

  /* Optional single-input covering stage (cover_buffer): once the raw sink sums are known, every
   * covered signal decides between the zero-cost identity ("no buffer" — the wire branch) and a
   * library buffer drive, priced by the same block + slope × load model as every other candidate.
   * An active buffer REWRITES node_loads[index][phase] to the load the DRIVER actually sees (the
   * buffer's input capacitance + wire), so every existing reader — gate_delay, gate_load_legal,
   * the curve re-bind — prices the buffered driver without knowing buffers exist; the sink load
   * the buffer carries is kept in the choice for legality reporting and buffer_delay_at.
   *
   * Decision rule, mirroring the match pick's legality-first shape: a buffer is chosen when it
   * turns an illegal driver legal, or when both branches are legal and it strictly reduces the
   * driver-path delay. Hysteresis: an active buffer is kept unless strictly worse (and never
   * while it is the only legal branch) — the load refresh is exact and undamped, and a
   * discontinuous Σsinks → cap flip is exactly what an undamped fixed point can oscillate on.
   *
   * v1 scope (recorded residual): gate outputs serving ONE polarity only — PI signals, constant
   * signals, and same_match nodes whose complement is also referenced (the materialized polarity
   * inverter would have to re-read the buffered signal and its own load model) keep bare nets. */
  void select_output_buffers()
  {
    if ( !buffer_stage )
      return;
    float const wire = static_cast<float>( ps.wire_cap_per_fanout );
    uint32_t count = 0;
    for ( auto const& n : topo_order )
    {
      uint32_t const index = ntk.node_to_index( n );
      if ( ntk.is_constant( n ) || ntk.is_pi( n ) )
      {
        node_buffers[index][0].active = node_buffers[index][1].active = false;
        continue;
      }
      auto const& node_data = node_match[index];
      bool const dual_service = node_data.same_match && node_data.map_refs[0] > 0 && node_data.map_refs[1] > 0;
      for ( uint8_t phase = 0; phase < 2; ++phase )
      {
        auto& choice = node_buffers[index][phase];
        supergate<NInputs> const* g = node_data.best_gate[phase];
        float const raw_load = node_loads[index][phase];
        if ( g == nullptr || dual_service || node_data.map_refs[phase] == 0 || raw_load <= 0.0f )
        {
          choice.active = false;
          continue;
        }
        auto const drive = library.select_buffer( raw_load );
        if ( drive.id == UINT32_MAX )
        {
          choice.active = false;
          continue;
        }
        float const driver_load_buffered = drive.cap + wire;
        /* the driver's worst output-load sensitivity over the bound cut's pins */
        float slope_g = 0.0f;
        {
          auto const& best_cut = cuts[index][node_data.best_cut[phase]];
          uint32_t const pins = std::min<uint32_t>( static_cast<uint32_t>( best_cut.size() ), NInputs );
          for ( uint32_t ctr = 0; ctr < pins; ++ctr )
            slope_g = std::max( slope_g, g->slope[ctr] );
        }
        bool const legal_plain = g->max_load <= 0.0f || raw_load <= g->max_load + epsilon;
        bool const legal_buf = ( g->max_load <= 0.0f || driver_load_buffered <= g->max_load + epsilon ) &&
                               ( drive.max_load <= 0.0f || raw_load <= drive.max_load + epsilon );
        double const d_plain = static_cast<double>( slope_g ) * raw_load;
        double const d_buf = static_cast<double>( slope_g ) * driver_load_buffered + drive.block + drive.slope * raw_load;
        /* a delay-driven buffer is offered only where timing actually needs help: on a met signal
         * the identity (wire) branch wins even when the buffer's flat slope models faster —
         * otherwise a steep-slope library would buy buffers on every met multi-sink net, paying
         * area for slack nobody asked for (the required time IS the ask) */
        bool const critical = node_data.arrival[phase] + epsilon > node_data.required[phase];
        bool take;
        if ( legal_buf != legal_plain )
        {
          take = legal_buf; /* legality first, the model's own rule */
        }
        else if ( choice.active )
        {
          /* hysteresis: keep unless strictly worse — criticality gates only FRESH insertions, or
           * a buffer that just met the node's timing would be removed for having met it */
          take = d_buf <= d_plain + epsilon;
        }
        else
        {
          take = critical && d_buf + epsilon < d_plain; /* insert only on strict improvement */
        }
        if ( take )
        {
          choice.active = true;
          choice.drive = drive;
          choice.sink_load = raw_load;
          node_loads[index][phase] = driver_load_buffered;
          ++count;
        }
        else
        {
          choice.active = false;
        }
      }
    }
    st.cover_buffers = count;
  }

  /*! \brief Re-derive arrivals and required times of the committed cover under the current load
   * estimate. Uses the existing per-node arrival propagation, so it introduces no second copy of
   * the arrival walk and has none of `propagate_arrival_times`'s area/iteration side effects. */
  void recompute_arrivals_and_required()
  {
    if ( !ps.electrical_model )
      return;
    for ( auto const& n : topo_order )
    {
      if ( ntk.is_constant( n ) || ntk.is_pi( n ) )
        continue;
      uint32_t const index = ntk.node_to_index( n );
      auto const& node_data = node_match[index];

      /* only nodes actually IN the cover, mirroring compute_required_time's own guard */
      if ( !node_data.map_refs[0] && !node_data.map_refs[1] )
        continue;

      /* propagate_arrival_node's precondition: it derives the second phase either through the
       * shared-match inverter or from that phase's own gate, so a node matched in ONE phase with
       * same_match false has no gate to read there. Honoring this here is not optional — an
       * asserts-enabled build traps, and a release build dereferences null. */
      if ( !node_data.same_match &&
           ( node_data.best_gate[0] == nullptr || node_data.best_gate[1] == nullptr ) )
        continue;
      if ( node_data.best_gate[0] == nullptr && node_data.best_gate[1] == nullptr )
        continue;

      propagate_arrival_node( n );
    }

    delay = 0.0;
    ntk.foreach_po( [this]( auto s ) {
      const auto index = ntk.node_to_index( ntk.get_node( s ) );
      delay = std::max( delay, node_match[index].arrival[ntk.is_complemented( s ) ? 1 : 0] );
    } );

    compute_required_time();
  }

  /*! \brief Measure the committed cover's drive legality into the stats (electrical model only).
   *
   * The presence test mirrors the counting walk (`same_match || map_refs[use_phase]`, then the
   * other phase): with `same_match` a node's references can sit entirely in the phase whose
   * best_gate is null (the inverter-served one), so the intuitive `map_refs[p] && best_gate[p]`
   * pairing would skip such a node in BOTH phases and silently report a clean 0. */
  void report_drive_legality()
  {
    if ( !ps.electrical_model )
      return;
    st.max_load_violations = 0;
    st.worst_load_ratio = 0.0;
    auto check = [&]( uint32_t index, uint8_t phase ) {
      auto const& node_data = node_match[index];
      if ( node_data.best_gate[phase] == nullptr )
        return;
      /* two-segment check under an active output buffer: the driver sees the buffer's input
       * (already stored in node_loads by select_output_buffers), the buffer carries the sinks */
      if ( buffer_stage && node_buffers[index][phase].active )
      {
        auto const& b = node_buffers[index][phase];
        if ( b.drive.max_load > 0.0f )
        {
          st.worst_load_ratio =
              std::max( st.worst_load_ratio, static_cast<double>( b.sink_load ) / b.drive.max_load );
          if ( b.sink_load > b.drive.max_load )
            ++st.max_load_violations;
        }
      }
      float const limit = node_data.best_gate[phase]->max_load;
      if ( limit <= 0.0f )
        return;
      st.worst_load_ratio =
          std::max( st.worst_load_ratio, static_cast<double>( node_loads[index][phase] ) / limit );
      if ( node_loads[index][phase] > limit )
        ++st.max_load_violations;
    };
    for ( auto const& n : topo_order )
    {
      if ( ntk.is_constant( n ) || ntk.is_pi( n ) )
        continue;
      uint32_t const index = ntk.node_to_index( n );
      auto const& node_data = node_match[index];
      if ( !node_data.map_refs[0] && !node_data.map_refs[1] )
        continue;
      uint8_t const use_phase = node_data.best_gate[0] == nullptr ? 1u : 0u;
      if ( node_data.best_gate[use_phase] == nullptr )
        continue;
      if ( node_data.same_match || node_data.map_refs[use_phase] > 0 )
        check( index, use_phase );
      if ( !node_data.same_match && node_data.map_refs[use_phase ^ 1] > 0 )
        check( index, use_phase ^ 1 );
    }
  }

  /* --- area-delay frontier (electrical model; ps.curve_points > 0) ------------------------- */

  /*! \brief Offer a matched candidate to a node's bounded (arrival, area) frontier.
   *
   * Dominance-pruned: a point no better in BOTH arrival and area than an existing one is
   * dropped, and a point that dominates existing ones evicts them. When the budget is full the
   * worst-area point is evicted only if the newcomer is cheaper, so the frontier keeps its
   * extremes. */
  void curve_offer( uint32_t index, uint8_t phase, supergate<NInputs> const* gate, double arrival,
                    float area, float flow, uint16_t polarity, uint32_t cut_index, uint32_t size )
  {
    if ( curve_budget == 0 || index >= curves.size() )
      return;

    curve_set& cs = curves[index][phase];

    for ( uint8_t i = 0; i < cs.size; ++i )
    {
      auto const& p = cs.points[i];
      if ( p.cut != cut_index )
        continue;
      /* an existing point on this cut is at least as good in both dimensions */
      if ( p.arrival <= arrival + epsilon && p.area <= area + epsilon )
        return;
    }

    /* drop the points this candidate dominates */
    uint8_t w = 0;
    for ( uint8_t i = 0; i < cs.size; ++i )
    {
      auto const& p = cs.points[i];
      const bool dominated =
          ( p.cut == cut_index ) && ( arrival <= p.arrival + epsilon ) && ( area <= p.area + epsilon );
      if ( !dominated )
        cs.points[w++] = cs.points[i];
    }
    cs.size = w;

    if ( cs.size >= max_curve_points )
    {
      /* full: replace the largest-area point, and only if this one is cheaper */
      uint8_t worst = 0;
      for ( uint8_t i = 1; i < cs.size; ++i )
        if ( cs.points[i].area > cs.points[worst].area )
          worst = i;
      if ( cs.points[worst].area <= area )
        return;
      cs.points[worst] = best_gate_emap<NInputs>{ gate, arrival, area, flow, polarity, cut_index, size };
      return;
    }

    cs.points[cs.size++] =
        best_gate_emap<NInputs>{ gate, arrival, area, flow, polarity, cut_index, size };
  }

  /*! \brief Re-bind each covered node to the cheapest frontier point that still meets its
   * required time, now that the required times come from the COMPLETE cover.
   *
   * Only points on the node's chosen cut are eligible, so a re-bind never changes which leaves
   * the node reads: references, area accounting and the cover structure stay valid and only the
   * cell changes. Each candidate's arrival is recomputed from the current leaf arrivals rather
   * than reusing the value stored at match time, which was measured under a different load
   * estimate.
   *
   * This is what `use_match_alternatives` cannot do: it also keeps a second point per phase, but
   * picks it during matching, against the PREVIOUS round's required times. */
  void reselect_from_curves()
  {
    if ( curve_budget == 0 )
      return;
    uint32_t rebinds = 0;

    for ( auto const& n : topo_order )
    {
      if ( ntk.is_constant( n ) || ntk.is_pi( n ) )
        continue;

      uint32_t const index = ntk.node_to_index( n );
      if ( index >= curves.size() )
        continue;
      auto& node_data = node_match[index];

      for ( uint8_t phase = 0; phase < 2; ++phase )
      {
        if ( node_data.map_refs[phase] == 0 || node_data.best_gate[phase] == nullptr )
          continue;

        curve_set const& cs = curves[index][phase];
        if ( cs.size < 2 )
          continue;

        double const required = node_data.required[phase];
        auto const& cut = cuts[index][node_data.best_cut[phase]];

        supergate<NInputs> const* best = node_data.best_gate[phase];
        float best_area = node_data.area[phase];
        double best_arrival = node_data.arrival[phase];
        uint16_t best_polarity = node_data.phase[phase];

        for ( uint8_t i = 0; i < cs.size; ++i )
        {
          auto const& p = cs.points[i];
          if ( p.gate == nullptr || p.cut != node_data.best_cut[phase] )
            continue;

          /* recompute this candidate's arrival under the CURRENT load estimate */
          double arrival = 0.0;
          uint32_t ctr = 0u;
          for ( auto l : cut )
          {
            uint8_t const leaf_phase = ( p.phase >> ctr ) & 1;
            arrival = std::max( arrival, node_match[l].arrival[leaf_phase] +
                                             gate_delay( p.gate, ctr, index, phase ) );
            ++ctr;
          }

          if ( arrival > required + epsilon )
            continue; /* does not meet the node's required time */
          if ( ps.electrical_model && !gate_load_legal( p.gate, index, phase ) &&
               gate_load_legal( best, index, phase ) )
            continue; /* never trade a drive-legal binding for an illegal one */
          if ( p.area + epsilon >= best_area )
            continue; /* not cheaper */

          best = p.gate;
          best_area = p.area;
          best_arrival = arrival;
          best_polarity = p.phase;
        }

        if ( best != node_data.best_gate[phase] )
        {
          node_data.best_gate[phase] = best;
          node_data.area[phase] = best_area;
          node_data.arrival[phase] = best_arrival;
          node_data.phase[phase] = best_polarity;
          ++rebinds;
        }
      }
    }

    if ( ps.verbose && rebinds > 0 )
    {
      std::stringstream stats{};
      stats << fmt::format( "[i] Frontier : re-bound {} node(s) against the complete cover's required times\n", rebinds );
      st.round_stats.push_back( stats.str() );
    }
  }

  /* --- end electrical model --------------------------------------------------------------- */

  /* --- frontier through the covering DP (emap_params::frontier_points) --------------------- */

  /* One point of a node-phase frontier: implementing the phase this way can arrive at `arrival`
   * for an estimated area-flow `cost` (the gate plus each leaf's cost shared over its estimated
   * references). A `via_inv` point is the opposite phase's DIRECT point `cut`/`gate`/`polarity`
   * plus the library inverter; it never refers to another via_inv point, so resolution cannot
   * loop between phases. */
  struct frontier_point
  {
    double arrival;
    float cost;
    supergate<NInputs> const* gate;
    uint16_t cut;
    uint16_t polarity;
    bool via_inv;
    uint64_t order{ 0 }; /* stable tie order without stable_sort's unbounded temporary allocation */
  };

  /* Non-owning bounded view. One exact-size point slab owns every stored frontier; candidate
   * slabs are also sized at admission. No insertion can trigger hidden capacity growth. */
  struct frontier_list
  {
    frontier_point* data{ nullptr };
    std::size_t count{ 0 }, capacity{ 0 };
    frontier_point* begin() const { return data; }
    frontier_point* end() const { return count ? data + count : data; }
    std::size_t size() const { return count; }
    bool empty() const { return count == 0; }
    void clear() { count = 0; }
    frontier_point& operator[]( std::size_t i ) const { return data[i]; }
    void push_back( frontier_point const& p )
    {
      if ( count == capacity )
        throw std::length_error( "frontier candidate capacity exhausted" );
      data[count++] = p;
    }
  };
  struct frontier_decision
  {
    bool used[2]{ false, false };
    frontier_point pick[2]{};
    int inverted{ -1 };
  };

  /* Exact-size iterative DFS frames preserve recursive child-return float accumulation and
   * the same-match increment AFTER the child returns. Only the opt-in round uses this path. */
  struct frontier_walk_frame
  {
    cut_t const* cut{ nullptr };
    uint32_t next{ 0 };
    uint16_t polarity{ 0 };
    float count{ 0 };
    uint64_t pending{ std::numeric_limits<uint64_t>::max() };
  };
  enum class frontier_walk_mode
  {
    reference,
    dereference,
    visit
  };
  frontier_walk_frame* frontier_frames{ nullptr };
  uint64_t* frontier_visits{ nullptr };
  std::size_t frontier_frame_capacity{ 0 }, frontier_visit_capacity{ 0 }, frontier_visit_count{ 0 };

  float frontier_cut_walk( cut_t const& root_cut, node<Ntk> const& root, uint8_t phase, frontier_walk_mode mode )
  {
    std::size_t depth = 0;
    auto push = [&]( cut_t const& cut, uint32_t index, uint8_t ph ) {
      if ( depth == frontier_frame_capacity )
        throw std::length_error( "frontier traversal frame capacity exhausted" );
      frontier_frames[depth++] = { &cut, 0, node_match[index].phase[ph], node_match[index].area[ph], std::numeric_limits<uint64_t>::max() };
    };
    auto add_reference = [&]( frontier_walk_frame& frame, uint32_t leaf, uint8_t ph ) {
      auto& d = node_match[leaf];
      if ( d.map_refs[ph]++ == 0u && d.best_gate[ph] == nullptr )
        frame.count += lib_inv_area;
    };
    push( root_cut, ntk.node_to_index( root ), phase );
    while ( depth )
    {
      auto& frame = frontier_frames[depth - 1];
      if ( frame.pending != std::numeric_limits<uint64_t>::max() )
      {
        add_reference( frame, static_cast<uint32_t>( frame.pending >> 1 ), frame.pending & 1 );
        frame.pending = std::numeric_limits<uint64_t>::max();
      }
      if ( frame.next == frame.cut->size() )
      {
        float const count = frame.count;
        if ( --depth == 0 )
          return count;
        frontier_frames[depth - 1].count += count;
        continue;
      }
      uint32_t const position = frame.next++;
      uint32_t const leaf = *( frame.cut->begin() + position );
      /* Polarity is captured from the current gate when its frame is pushed. */
      uint8_t const ph = ( frame.polarity >> position ) & 1;
      auto const n = ntk.index_to_node( leaf );
      if ( ntk.is_constant( n ) )
        continue;
      if ( mode == frontier_walk_mode::visit )
      {
        if ( frontier_visit_count == frontier_visit_capacity )
          throw std::length_error( "frontier visit log capacity exhausted" );
        frontier_visits[frontier_visit_count++] = ( uint64_t( leaf ) << 1 ) | ph;
      }
      auto& d = node_match[leaf];
      if ( ntk.is_pi( n ) )
      {
        bool const last = mode == frontier_walk_mode::dereference ? --d.map_refs[ph] == 0u : d.map_refs[ph]++ == 0u;
        if ( ph && last )
          frame.count += lib_inv_area;
        continue;
      }
      bool descend = false;
      if ( mode == frontier_walk_mode::dereference )
      {
        if ( d.same_match )
        {
          if ( --d.map_refs[ph] == 0u && d.best_gate[ph] == nullptr )
            frame.count += lib_inv_area;
          descend = !d.map_refs[0] && !d.map_refs[1];
        }
        else
          descend = --d.map_refs[ph] == 0u;
      }
      else if ( d.same_match )
      {
        descend = !d.map_refs[0] && !d.map_refs[1];
        if ( descend )
          frame.pending = ( uint64_t( leaf ) << 1 ) | ph;
        else
          add_reference( frame, leaf, ph );
      }
      else
        descend = d.map_refs[ph]++ == 0u;
      if ( descend )
        push( cuts[leaf][d.best_cut[ph]], leaf, ph );
    }
    return 0;
  }

  /* Keep the non-dominated points, fastest first (cost strictly falling), then cut the list to
   * `limit`: the fastest and the cheapest always stay, the rest are taken evenly by rank. The
   * sort is stable over insertion order (cut order, then match order), so ties resolve the same
   * way every run. */
  bool frontier_prune( frontier_list& pts, uint32_t limit ) const
  {
    if ( pts.empty() )
      return false;
    for ( std::size_t i = 0; i < pts.size(); ++i )
      pts[i].order = i;
    std::sort( pts.begin(), pts.end(), []( frontier_point const& a, frontier_point const& b ) {
      if ( a.arrival != b.arrival )
        return a.arrival < b.arrival;
      if ( a.cost != b.cost )
        return a.cost < b.cost;
      return a.order < b.order;
    } );
    std::size_t kept = 0;
    float best_cost = std::numeric_limits<float>::max();
    for ( std::size_t i = 0; i < pts.size(); ++i )
    {
      if ( pts[i].cost < best_cost - epsilon )
      {
        best_cost = pts[i].cost;
        pts[kept++] = pts[i];
      }
    }
    bool const trimmed = kept > limit;
    if ( trimmed )
    {
      /* Selected source indices increase and are >= destination indices: safe in place. */
      for ( uint32_t j = 0; j < limit; ++j )
      {
        std::size_t const at = ( uint64_t( j ) * ( kept - 1 ) + ( limit - 1 ) / 2 ) / ( limit - 1 );
        pts[j] = pts[at];
      }
      kept = limit;
    }
    pts.count = kept;
    return trimmed;
  }

  /* The frontier a consumer sees for (leaf, phase): a PI or constant is one fixed point (a negated
   * PI costs the inverter, shared like any other leaf cost); an internal node is its own list. */
  frontier_list const& frontier_of( uint32_t leaf, uint8_t phase, std::array<frontier_list, 2> const* fr,
                                    frontier_list& scratch ) const
  {
    auto const n = ntk.index_to_node( leaf );
    if ( ntk.is_pi( n ) || ntk.is_constant( n ) )
    {
      scratch.clear();
      float const cost = ( ntk.is_pi( n ) && phase == 1 ) ? static_cast<float>( lib_inv_area ) : 0.0f;
      scratch.push_back( frontier_point{ node_match[leaf].arrival[phase], cost, nullptr, 0, 0, false } );
      return scratch;
    }
    return fr[leaf][phase];
  }

  /* Merge the leaves' frontiers through one match. Every leaf starts at its cheapest point; the
   * merged arrival is the latest shifted leaf arrival, so the only way to be faster is to move
   * every leaf that attains it to its next faster point. Each step emits one merged point, and the
   * walk ends when a critical leaf has no faster point: O(sum of leaf frontier sizes). */
  void frontier_merge_match( uint32_t index, uint8_t phase, supergate<NInputs> const& gate, uint16_t polarity,
                             uint16_t cut_index, cut_t const& cut, std::array<frontier_list, 2> const* fr,
                             frontier_list& out, std::array<float, 2> const* frontier_share )
  {
    uint32_t const k = cut.size();
    std::array<frontier_list const*, NInputs> lists{};
    std::array<frontier_point, NInputs> scratch_points{};
    std::array<frontier_list, NInputs> scratch{};
    for ( uint32_t j = 0; j < NInputs; ++j )
      scratch[j] = { &scratch_points[j], 0, 1 };
    std::array<double, NInputs> shift{};
    std::array<float, NInputs> share{};
    std::array<uint32_t, NInputs> at{};
    uint32_t i = 0;
    for ( auto l : cut )
    {
      uint8_t const lph = ( polarity >> i ) & 1;
      lists[i] = &frontier_of( l, lph, fr, scratch[i] );
      if ( lists[i]->empty() )
        return; /* the leaf phase cannot be implemented: neither can this match */
      shift[i] = gate_delay( &gate, i, index, phase );
      share[i] = frontier_share[l][lph];
      at[i] = static_cast<uint32_t>( lists[i]->size() - 1 ); /* cheapest */
      ++i;
    }
    while ( true )
    {
      double t = 0.0;
      float cost = static_cast<float>( gate.area );
      for ( uint32_t j = 0; j < k; ++j )
      {
        auto const& p = ( *lists[j] )[at[j]];
        t = std::max( t, p.arrival + shift[j] );
        cost += p.cost / share[j];
      }
      if ( t >= std::numeric_limits<float>::max() )
        return;
      out.push_back( frontier_point{ t, cost, &gate, cut_index, polarity, false } );
      bool moved = false;
      for ( uint32_t j = 0; j < k; ++j )
      {
        if ( ( *lists[j] )[at[j]].arrival + shift[j] < t - epsilon )
          continue; /* not critical */
        if ( at[j] == 0 )
          return; /* a critical leaf is already at its fastest point */
        --at[j];
        moved = true;
      }
      if ( !moved )
        return;
    }
  }

  /* Reference (leaf, phase) as one more consumer -- the per-leaf body of cut_ref, for a primary
   * output, which has no cut of its own. */
  void frontier_ref_signal( uint32_t leaf, uint8_t leaf_phase )
  {
    auto const n = ntk.index_to_node( leaf );
    if ( ntk.is_constant( n ) )
      return;
    if ( ntk.is_pi( n ) )
    {
      ++node_match[leaf].map_refs[leaf_phase];
      return;
    }
    auto& d = node_match[leaf];
    if ( d.same_match )
    {
      if ( !d.map_refs[0] && !d.map_refs[1] )
        cut_ref<false>( cuts[leaf][d.best_cut[leaf_phase]], n, leaf_phase );
      ++d.map_refs[leaf_phase];
    }
    else if ( d.map_refs[leaf_phase]++ == 0u )
    {
      cut_ref<false>( cuts[leaf][d.best_cut[leaf_phase]], n, leaf_phase );
    }
  }

  double frontier_current_delay() const
  {
    double worst = 0.0;
    ntk.foreach_po( [&]( auto const& s ) {
      auto const index = ntk.node_to_index( ntk.get_node( s ) );
      worst = std::max( worst, node_match[index].arrival[ntk.is_complemented( s ) ? 1 : 0] );
    } );
    return worst;
  }

  bool compute_mapping_frontier()
  {
    /* Shapes this round does not model keep the ordinary cover, and say so. */
    auto decline = [&]( char const* why ) {
      st.frontier_declined = why;
      return true;
    };
    if ( ps.map_multioutput )
      return decline( "multi-output matching is on (a multi-output match has no per-phase frontier)" );
    if ( ps.cover_buffer )
      return decline( "cover-buffer is on (a chosen buffer is not a frontier point)" );
    if ( ps.load_points > 1 )
      return decline( "load-indexed curves are on (their arrivals are not scalar points)" );
    if constexpr ( has_is_dont_touch_v<Ntk> )
    {
      bool boxed = false;
      ntk.foreach_node( [&]( auto const& n ) { boxed = boxed || ntk.is_dont_touch( n ); } );
      if ( boxed )
        return decline( "the network holds dont-touch boxes" );
    }

    uint32_t const limit = std::max<uint32_t>( 2u, std::min( ps.frontier_points, max_frontier_points ) );
    /* Derive scratch capacities from existing cuts/matches, without allocating or mutating
     * the entering cover. A match emits at most 1 + sum(leaf_size - 1) points. */
    uint64_t merged_capacity = 0, direct_capacity = 0;
    auto checked_add = []( uint64_t a, uint64_t b, uint64_t& result ) {
      if ( b > std::numeric_limits<uint64_t>::max() - a )
        return false;
      result = a + b;
      return true;
    };
    auto checked_mul = []( uint64_t a, uint64_t b, uint64_t& result ) {
      if ( b && a > std::numeric_limits<uint64_t>::max() / b )
        return false;
      result = a * b;
      return true;
    };
    for ( auto const& n : topo_order )
    {
      if ( ntk.is_constant( n ) || ntk.is_pi( n ) )
        continue;
      auto const index = ntk.node_to_index( n );
      for ( uint8_t ph = 0; ph < 2; ++ph )
      {
        uint64_t phase_capacity = 0;
        for ( auto const& cut : cuts[index] )
        {
          if ( ( *cut )->ignore || ( *cut )->supergates[ph] == nullptr )
            continue;
          uint64_t cap = 0;
          if ( !checked_mul( ( *cut )->supergates[ph]->size(), 1 + uint64_t( cut->size() ) * ( limit - 1 ), cap ) ||
               !checked_add( phase_capacity, limit, phase_capacity ) )
            return decline( "frontier memory budget accounting overflow" );
          merged_capacity = std::max( merged_capacity, cap );
        }
        direct_capacity = std::max( direct_capacity, phase_capacity );
      }
    }
    uint64_t const nodes = ntk.size();
    uint64_t point_count = 0, scratch_count = 0, visit_capacity = 0, walk_bytes = 0, visit_bytes = 0;
    if ( !checked_mul( nodes, 2 * limit, point_count ) ||
         !checked_mul( direct_capacity, 2, scratch_count ) ||
         !checked_add( scratch_count, merged_capacity, scratch_count ) ||
         !checked_add( scratch_count, 2 * limit, scratch_count ) ||
         !checked_mul( point_count, sizeof( frontier_point ), st.frontier_bytes_points ) ||
         !checked_mul( nodes, sizeof( std::array<frontier_list, 2> ), st.frontier_bytes_headers ) ||
         !checked_mul( nodes, sizeof( node_match_emap<NInputs> ), st.frontier_bytes_saved_matches ) ||
         !checked_mul( nodes, sizeof( std::array<float, 2> ), st.frontier_bytes_sharing ) ||
         !checked_mul( nodes, 2 * sizeof( std::array<double, 2> ) + sizeof( std::array<bool, 2> ) + sizeof( frontier_decision ), st.frontier_bytes_resolution ) ||
         !checked_mul( scratch_count, sizeof( frontier_point ), st.frontier_bytes_scratch ) )
      return decline( "frontier memory budget accounting overflow" );
    /* Each phase implementation expands at most once during a visit trial; every expanded
     * cut adds at most CutSize leaf records. References prevent repeated shared expansion,
     * including same_match after its DFS child returns. An acyclic path has <= nodes frames. */
    if ( !checked_mul( nodes, 2 * uint64_t( CutSize ), visit_capacity ) ||
         !checked_mul( visit_capacity, sizeof( uint64_t ), visit_bytes ) ||
         !checked_mul( nodes, sizeof( frontier_walk_frame ), walk_bytes ) ||
         !checked_add( st.frontier_bytes_scratch, visit_bytes, st.frontier_bytes_scratch ) ||
         !checked_add( st.frontier_bytes_scratch, walk_bytes, st.frontier_bytes_scratch ) )
      return decline( "frontier memory budget accounting overflow" );
    /* Fixed allowance for local merge/resolve frames and allocation-free sort's logarithmic
     * recursion. Fourteen exact-size heap arrays are rounded individually below, including
     * array cookies/padding before size-class rounding. This is not a process-RSS limit;
     * allocator arenas, pre-existing mapper/network memory and OS residency are excluded. */
    if ( !checked_add( st.frontier_bytes_scratch, 16384 + uint64_t( NInputs ) * ( sizeof( frontier_point ) + sizeof( frontier_list ) + sizeof( frontier_list* ) + sizeof( double ) + sizeof( float ) + sizeof( uint32_t ) ), st.frontier_bytes_scratch ) )
      return decline( "frontier memory budget accounting overflow" );
    static_assert( alignof( frontier_point ) <= 64 && alignof( frontier_walk_frame ) <= 64 &&
                       alignof( frontier_decision ) <= 64 && alignof( frontier_list ) <= 64 &&
                       alignof( node_match_emap<NInputs> ) <= 64,
                   "frontier array padding bound requires alignment <= 64" );
    static_assert( std::is_trivially_destructible_v<frontier_point> &&
                       std::is_trivially_destructible_v<frontier_walk_frame> &&
                       std::is_trivially_destructible_v<frontier_decision> &&
                       std::is_trivially_destructible_v<frontier_list> &&
                       std::is_trivially_destructible_v<node_match_emap<NInputs>>,
                   "frontier accounting assumes trivial array destruction" );
    uint64_t req_bytes = 0, demand_bytes = 0, decision_bytes = 0, direct_bytes = 0, merged_bytes = 0;
    if ( !checked_mul( nodes, sizeof( std::array<double, 2> ), req_bytes ) ||
         !checked_mul( nodes, sizeof( std::array<bool, 2> ), demand_bytes ) ||
         !checked_mul( nodes, sizeof( frontier_decision ), decision_bytes ) ||
         !checked_mul( direct_capacity, sizeof( frontier_point ), direct_bytes ) ||
         !checked_mul( merged_capacity, sizeof( frontier_point ), merged_bytes ) )
      return decline( "frontier memory budget accounting overflow" );
    std::array<uint64_t, 14> const allocation_requests{
        st.frontier_bytes_sharing, st.frontier_bytes_saved_matches, req_bytes,
        st.frontier_bytes_points, st.frontier_bytes_headers, direct_bytes, direct_bytes,
        merged_bytes, uint64_t( 2 * limit ) * sizeof( frontier_point ), walk_bytes,
        visit_bytes, req_bytes, demand_bytes, decision_bytes };
    st.frontier_bytes_allocator = 0;
    for ( auto request : allocation_requests )
    {
      uint64_t accounted = 0;
      if ( !frontier_allocation_bytes( request, accounted ) ||
           !checked_add( st.frontier_bytes_allocator, accounted - request, st.frontier_bytes_allocator ) )
        return decline( "frontier memory budget accounting overflow" );
    }
    st.frontier_bytes = 0;
    for ( uint64_t component : { st.frontier_bytes_points, st.frontier_bytes_headers,
                                 st.frontier_bytes_saved_matches, st.frontier_bytes_sharing, st.frontier_bytes_resolution,
                                 st.frontier_bytes_scratch, st.frontier_bytes_allocator } )
      if ( !checked_add( st.frontier_bytes, component, st.frontier_bytes ) )
        return decline( "frontier memory budget accounting overflow" );
    long double const budget_bytes = static_cast<long double>( ps.frontier_mem_budget_mb ) * 1048576.0L;
    if ( !std::isfinite( ps.frontier_mem_budget_mb ) || ps.frontier_mem_budget_mb < 0 ||
         budget_bytes >= static_cast<long double>( std::numeric_limits<uint64_t>::max() ) )
      return decline( "invalid or overflowing frontier_mem_budget_mb budget" );
    if ( st.frontier_bytes > std::numeric_limits<std::size_t>::max() ||
         static_cast<long double>( st.frontier_bytes ) > budget_bytes )
      return decline( "total incremental frontier storage exceeds frontier_mem_budget_mb budget" );

    /* All owning arrays have exact extents. They die at round exit (including declines), so
     * repeated rounds cannot overlap retained auxiliary capacity with the next admission. */
    auto frontier_share = std::make_unique<std::array<float, 2>[]>( nodes );
    for ( uint32_t i = 0; i < ntk.size(); ++i )
      for ( uint8_t ph = 0; ph < 2; ++ph )
        frontier_share[i][ph] = ps.frontier_share == emap_params::frontier_share_t::cover
                                    ? std::max( 1.0f, static_cast<float>( node_match[i].map_refs[ph] ) )
                                    : std::max( 1.0f, node_match[i].est_refs[ph] );

    st.frontier_area_before = area;
    st.frontier_delay_before = frontier_current_delay();
    /* the entering cover, restored whole if the resolution does not improve on it */
    auto saved_match = std::make_unique<node_match_emap<NInputs>[]>( nodes );
    std::copy( node_match.begin(), node_match.end(), saved_match.get() );
    double const saved_area = area;
    uint32_t const saved_inv = inv;

    /* required times at the outputs of the CURRENT cover (the delay target the round keeps) */
    delay = st.frontier_delay_before;
    compute_required_time( true );
    auto po_req = std::make_unique<std::array<double, 2>[]>( nodes );
    std::fill_n( po_req.get(), nodes, std::array<double, 2>{ std::numeric_limits<double>::max(), std::numeric_limits<double>::max() } );
    ntk.foreach_po( [&]( auto const& s ) {
      auto const index = ntk.node_to_index( ntk.get_node( s ) );
      uint8_t const ph = ntk.is_complemented( s ) ? 1 : 0;
      po_req[index][ph] = std::min<double>( po_req[index][ph], node_match[index].required[ph] );
    } );

    /* forward: per-node frontiers in topological order */
    auto points = std::make_unique<frontier_point[]>( point_count );
    auto fr = std::make_unique<std::array<frontier_list, 2>[]>( nodes );
    for ( uint64_t i = 0; i < nodes; ++i )
      for ( uint8_t ph = 0; ph < 2; ++ph )
        fr[i][ph] = { points.get() + ( 2 * i + ph ) * limit, 0, limit };
    auto direct_points0 = std::make_unique<frontier_point[]>( direct_capacity );
    auto direct_points1 = std::make_unique<frontier_point[]>( direct_capacity );
    auto merged_points = std::make_unique<frontier_point[]>( merged_capacity );
    auto all_points = std::make_unique<frontier_point[]>( 2 * limit );
    auto walk_frames = std::make_unique<frontier_walk_frame[]>( nodes );
    auto visit_log = std::make_unique<uint64_t[]>( visit_capacity );
    std::array<frontier_list, 2> direct{ frontier_list{ direct_points0.get(), 0, direct_capacity },
                                         frontier_list{ direct_points1.get(), 0, direct_capacity } };
    frontier_list merged{ merged_points.get(), 0, merged_capacity };
    frontier_list all{ all_points.get(), 0, 2 * limit };
    for ( auto const& n : topo_order )
    {
      if ( ntk.is_constant( n ) || ntk.is_pi( n ) )
        continue;
      uint32_t const index = ntk.node_to_index( n );
      direct[0].clear();
      direct[1].clear();
      for ( uint8_t ph = 0; ph < 2; ++ph )
      {
        uint16_t cut_index = 0;
        for ( auto& cut : cuts[index] )
        {
          if ( ( *cut )->ignore || ( *cut )->supergates[ph] == nullptr )
          {
            ++cut_index;
            continue;
          }
          auto const negation = ( *cut )->negations[ph];
          merged.clear();
          for ( auto const& gate : *( *cut )->supergates[ph] )
            frontier_merge_match( index, ph, gate, static_cast<uint16_t>( gate.polarity ^ negation ), cut_index, *cut, fr.get(), merged, frontier_share.get() );
          frontier_prune( merged, limit );
          for ( auto const& p : merged )
            direct[ph].push_back( p );
          ++cut_index;
        }
        if ( frontier_prune( direct[ph], limit ) )
          ++st.frontier_trimmed;
      }
      for ( uint8_t ph = 0; ph < 2; ++ph )
      {
        all.clear();
        for ( auto const& p : direct[ph] )
          all.push_back( p );
        double const inv_delay = inv_delay_at( index, ph );
        for ( auto const& p : direct[ph ^ 1] )
          all.push_back( frontier_point{ p.arrival + inv_delay, p.cost + static_cast<float>( lib_inv_area ), p.gate, p.cut, p.polarity, true } );
        if ( frontier_prune( all, limit ) )
          ++st.frontier_trimmed;
        st.frontier_max_size = std::max<uint32_t>( st.frontier_max_size, static_cast<uint32_t>( all.size() ) );
        for ( auto const& p : all )
          fr[index][ph].push_back( p );
      }
      if ( !fr[index][0].empty() || !fr[index][1].empty() )
        ++st.frontier_nodes;
    }

    /* backward: resolve each demanded (node, phase) against its required time */
    auto req = std::make_unique<std::array<double, 2>[]>( nodes );
    std::copy_n( po_req.get(), nodes, req.get() );
    auto demanded = std::make_unique<std::array<bool, 2>[]>( nodes );
    ntk.foreach_po( [&]( auto const& s ) {
      demanded[ntk.node_to_index( ntk.get_node( s ) )][ntk.is_complemented( s ) ? 1 : 0] = true;
    } );
    auto dec = std::make_unique<frontier_decision[]>( nodes );
    auto cheapest_meeting = [&]( frontier_list const& pts, double r, bool& met ) -> frontier_point const* {
      frontier_point const* best = nullptr;
      for ( auto const& p : pts )
      {
        if ( p.via_inv )
          continue;
        if ( p.arrival <= r + epsilon && ( best == nullptr || p.cost < best->cost ) )
          best = &p;
      }
      met = best != nullptr;
      if ( best == nullptr ) /* nothing meets it: the fastest direct point */
        for ( auto const& p : pts )
          if ( !p.via_inv && ( best == nullptr || p.arrival < best->arrival ) )
            best = &p;
      return best;
    };
    uint32_t unmet = 0;
    for ( auto it = topo_order.rbegin(); it != topo_order.rend(); ++it )
    {
      if ( ntk.is_constant( *it ) || ntk.is_pi( *it ) )
        continue;
      uint32_t const index = ntk.node_to_index( *it );
      if ( !demanded[index][0] && !demanded[index][1] )
        continue;
      double const inf = std::numeric_limits<double>::max();
      double const r0 = demanded[index][0] ? req[index][0] : inf;
      double const r1 = demanded[index][1] ? req[index][1] : inf;

      /* option: each demanded phase direct */
      struct option
      {
        frontier_point const* p[2]{ nullptr, nullptr };
        int inverted{ -1 };
        bool met{ true };
        float cost{ 0 };
        bool ok{ true };
      };
      auto direct_option = [&]() {
        option o;
        for ( uint8_t ph = 0; ph < 2; ++ph )
        {
          if ( !demanded[index][ph] )
            continue;
          bool met = false;
          o.p[ph] = cheapest_meeting( fr[index][ph], ph ? r1 : r0, met );
          if ( o.p[ph] == nullptr )
          {
            o.ok = false;
            return o;
          }
          o.met = o.met && met;
          o.cost += o.p[ph]->cost;
        }
        return o;
      };
      /* option: `src` direct, the other phase from it through the inverter */
      auto inverter_option = [&]( uint8_t src ) {
        option o;
        uint8_t const dst = src ^ 1;
        double const need = std::min( src ? r1 : r0, ( dst ? r1 : r0 ) - inv_delay_at( index, dst ) );
        bool met = false;
        o.p[src] = cheapest_meeting( fr[index][src], need, met );
        if ( o.p[src] == nullptr || !demanded[index][dst] )
        {
          o.ok = false;
          return o;
        }
        o.met = met;
        o.inverted = dst;
        o.cost = o.p[src]->cost + static_cast<float>( lib_inv_area );
        return o;
      };
      std::array<option, 3> opts{ direct_option(), inverter_option( 0 ), inverter_option( 1 ) };
      option const* best = nullptr;
      for ( auto const& o : opts )
      {
        if ( !o.ok )
          continue;
        if ( best == nullptr || ( o.met && !best->met ) || ( o.met == best->met && o.cost < best->cost - epsilon ) )
          best = &o;
      }
      if ( best == nullptr )
      {
        std::copy_n( saved_match.get(), nodes, node_match.begin() );
        return decline( "a covered node has no implementable phase" );
      }
      if ( !best->met )
        ++unmet;

      frontier_decision& d = dec[index];
      d.inverted = best->inverted;
      for ( uint8_t ph = 0; ph < 2; ++ph )
      {
        if ( best->p[ph] == nullptr )
          continue;
        d.used[ph] = true;
        d.pick[ph] = *best->p[ph];
        /* the phase's own requirement, tightened by the inverted phase it also feeds */
        double r = ph ? r1 : r0;
        if ( best->inverted == ( ph ^ 1 ) )
          r = std::min( r, ( ph ? r0 : r1 ) - inv_delay_at( index, ph ^ 1 ) );
        auto const& cut = cuts[index][d.pick[ph].cut];
        uint32_t i = 0;
        for ( auto l : cut )
        {
          uint8_t const lph = ( d.pick[ph].polarity >> i ) & 1;
          req[l][lph] = std::min( req[l][lph], r - gate_delay( d.pick[ph].gate, i, index, ph ) );
          demanded[l][lph] = true;
          ++i;
        }
      }
    }

    /* Use bounded iterative references for BOTH rebuilding and every exact-cleanup trial.
     * Scoped reset also protects mapper state on errors or early exits. */
    struct workspace_reset
    {
      emap_impl& owner;
      ~workspace_reset()
      {
        owner.frontier_frames = nullptr;
        owner.frontier_visits = nullptr;
        owner.frontier_frame_capacity = owner.frontier_visit_capacity = owner.frontier_visit_count = 0;
      }
    } reset{ *this };
    frontier_frames = walk_frames.get();
    frontier_visits = visit_log.get();
    frontier_frame_capacity = nodes;
    frontier_visit_capacity = visit_capacity;
    frontier_visit_count = 0;

    /* apply: rewrite every resolved node, then rebuild references from the outputs down */
    uint32_t cut_changes = 0, gate_changes = 0;
    for ( auto const& n : topo_order )
    {
      if ( ntk.is_constant( n ) || ntk.is_pi( n ) )
        continue;
      uint32_t const index = ntk.node_to_index( n );
      frontier_decision const& d = dec[index];
      if ( !d.used[0] && !d.used[1] )
        continue;
      auto& nd = node_match[index];
      for ( uint8_t ph = 0; ph < 2; ++ph )
      {
        if ( !d.used[ph] )
          continue;
        bool const was_live = nd.map_refs[ph] > 0 && nd.best_gate[ph] != nullptr;
        if ( was_live && nd.best_cut[ph] != d.pick[ph].cut )
          ++cut_changes;
        else if ( was_live && nd.best_gate[ph] != d.pick[ph].gate )
          ++gate_changes;
        nd.best_gate[ph] = d.pick[ph].gate;
        nd.best_cut[ph] = d.pick[ph].cut;
        nd.phase[ph] = d.pick[ph].polarity;
        nd.area[ph] = static_cast<float>( d.pick[ph].gate->area );
        nd.flows[ph] = d.pick[ph].cost;
      }
      nd.multioutput_match[0] = nd.multioutput_match[1] = false;
      if ( d.used[0] && d.used[1] )
      {
        nd.same_match = false;
      }
      else
      {
        /* one direct phase: the other (if demanded) is its complement through the inverter */
        uint8_t const src = d.used[0] ? 0 : 1;
        uint8_t const dst = src ^ 1;
        nd.same_match = true;
        nd.best_gate[dst] = nullptr;
        nd.best_cut[dst] = nd.best_cut[src];
        nd.phase[dst] = nd.phase[src];
        nd.area[dst] = nd.area[src];
      }
    }
    for ( auto& nd : node_match )
      nd.map_refs[0] = nd.map_refs[1] = 0u;
    ntk.foreach_po( [&]( auto const& s ) {
      frontier_ref_signal( ntk.node_to_index( ntk.get_node( s ) ), ntk.is_complemented( s ) ? 1 : 0 );
    } );
    propagate_arrival_times();

    /* exact-area clean-up against the resolved cover's own required times */
    for ( uint32_t round = 0; round < std::max( 1u, ps.ela_rounds ); ++round )
    {
      delay = frontier_current_delay();
      compute_required_time( true );
      if ( !compute_mapping_exact_reversed<false>( false ) )
        return false;
    }

    st.frontier_ran = true;
    st.frontier_cut_changes = cut_changes;
    st.frontier_gate_changes = gate_changes;
    st.frontier_unmet = unmet;
    st.frontier_area_after = area;
    st.frontier_delay_after = frontier_current_delay();

    /* Keep the resolution only where it did what it is for: the same delay target (every output
     * still meets the required time it entered with) at a smaller exact area. Otherwise the
     * entering cover comes back unchanged and the stats say so. */
    bool meets = unmet == 0u;
    ntk.foreach_po( [&]( auto const& s ) {
      auto const index = ntk.node_to_index( ntk.get_node( s ) );
      uint8_t const ph = ntk.is_complemented( s ) ? 1 : 0;
      if ( node_match[index].arrival[ph] > po_req[index][ph] + epsilon )
        meets = false;
    } );
    st.frontier_kept = meets && area < saved_area - epsilon;
    if ( !st.frontier_kept )
    {
      std::copy_n( saved_match.get(), nodes, node_match.begin() );
      inv = saved_inv;
      propagate_arrival_times(); /* restores loads, arrivals and area from the restored choices */
      area = saved_area;
    }
    delay = frontier_current_delay();
    return true;
  }

  /* --- end frontier ------------------------------------------------------------------------ */

  inline void match_propagate_required( uint32_t index )
  {
    /* don't touch box */
    if constexpr ( has_is_dont_touch_v<Ntk> )
    {
      node<Ntk> n = ntk.index_to_node( index );
      if ( ntk.is_dont_touch( n ) )
      {
        if constexpr ( has_has_binding_v<Ntk> )
        {
          propagate_data_backward_white_box( n );
        }
        return;
      }
    }

    auto& node_data = node_match[index];

    /* propagate required time through the leaves */
    unsigned use_phase = node_data.best_gate[0] == nullptr ? 1u : 0u;
    unsigned other_phase = use_phase ^ 1;

    assert( node_data.best_gate[0] != nullptr || node_data.best_gate[1] != nullptr );
    if ( node_data.best_gate[0] == nullptr && node_data.best_gate[1] == nullptr )
    {
      /* cover-completeness invariant broken (should be unreachable: match_phase_relaxed keeps
       * every referenced node matched) — report loudly instead of dereferencing null below */
      std::cerr << fmt::format( "[e] emap: node {} has no match in either phase; skipping required-time propagation\n", index );
      return;
    }
    // assert( node_data.map_refs[0] || node_data.map_refs[1] );

    /* propagate required time over the output inverter if present */
    if ( node_data.same_match && node_data.map_refs[use_phase ^ 1] > 0 )
    {
      node_data.required[use_phase] = std::min( node_data.required[use_phase], node_data.required[other_phase] - inv_delay_at( index, other_phase ) );
    }


    if ( node_data.map_refs[0] )
      assert( node_data.arrival[0] < node_data.required[0] + epsilon );
    if ( node_data.map_refs[1] )
      assert( node_data.arrival[1] < node_data.required[1] + epsilon );

    if ( node_data.same_match || node_data.map_refs[use_phase] > 0 )
    {
      auto ctr = 0u;
      auto const& best_cut = cuts[index][node_data.best_cut[use_phase]];
      auto const& supergate = node_data.best_gate[use_phase];
      /* the stored required is at the (optionally buffered) output the sinks see; the gate must
       * finish a buffer-delay earlier, so the buffer's delay comes off before the gate's */
      double const buf_use = buffer_delay_at( index, static_cast<uint8_t>( use_phase ) );
      for ( auto leaf : best_cut )
      {
        auto phase = ( node_data.phase[use_phase] >> ctr ) & 1;
        node_match[leaf].required[phase] = std::min( node_match[leaf].required[phase], node_data.required[use_phase] - buf_use - gate_delay( supergate, ctr, index, use_phase ) );
        ++ctr;
      }
    }

    if ( !node_data.same_match && node_data.map_refs[other_phase] > 0 )
    {
      auto ctr = 0u;
      auto const& best_cut = cuts[index][node_data.best_cut[other_phase]];
      auto const& supergate = node_data.best_gate[other_phase];
      double const buf_other = buffer_delay_at( index, static_cast<uint8_t>( other_phase ) );
      for ( auto leaf : best_cut )
      {
        auto phase = ( node_data.phase[other_phase] >> ctr ) & 1;
        node_match[leaf].required[phase] = std::min( node_match[leaf].required[phase], node_data.required[other_phase] - buf_other - gate_delay( supergate, ctr, index, other_phase ) );
        ++ctr;
      }
    }
  }

  template<bool ELA>
  bool set_mapping_refs()
  {
    /* compute the current worst delay and update the mapping refs */
    delay = 0.0f;
    ntk.foreach_po( [this]( auto s ) {
      const auto index = ntk.node_to_index( ntk.get_node( s ) );

      if ( ntk.is_complemented( s ) )
        delay = std::max( delay, node_match[index].arrival[1] );
      else
        delay = std::max( delay, node_match[index].arrival[0] );

      if constexpr ( !ELA )
      {
        if ( ntk.is_complemented( s ) )
          node_match[index].map_refs[1]++;
        else
          node_match[index].map_refs[0]++;
      }
    } );

    /* compute current area and update mapping refs in top-down order */
    area = 0.0f;
    inv = 0;
    for ( auto it = topo_order.rbegin(); it != topo_order.rend(); ++it )
    {
      const auto index = ntk.node_to_index( *it );
      auto& node_data = node_match[index];

      /* skip constants and PIs */
      if ( ntk.is_constant( *it ) )
      {
        if ( node_data.map_refs[0] || node_data.map_refs[1] )
        {
          /* if used and not available in the library launch a mapping error */
          if ( node_data.best_gate[0] == nullptr && node_data.best_gate[1] == nullptr )
          {
            std::cerr << "[e] MAP ERROR: technology library does not contain constant gates, impossible to perform mapping" << std::endl;
            st.mapping_error = true;
            return false;
          }
        }
        continue;
      }
      else if ( ntk.is_pi( *it ) )
      {
        if ( node_match[index].map_refs[1] > 0u )
        {
          /* Add inverter area over the negated fanins */
          area += lib_inv_area;
          ++inv;
        }
        continue;
      }

      /* continue if not referenced in the cover */
      if ( !node_match[index].map_refs[0] && !node_match[index].map_refs[1] )
        continue;

      /* don't touch box */
      if constexpr ( has_is_dont_touch_v<Ntk> )
      {
        if ( ntk.is_dont_touch( *it ) )
        {
          set_mapping_refs_dont_touch<ELA>( *it );
          continue;
        }
      }

      unsigned use_phase = node_data.best_gate[0] == nullptr ? 1u : 0u;

      if ( node_data.best_gate[use_phase] == nullptr )
      {
        /* Library is not complete, mapping is not possible */
        std::cerr << "[e] MAP ERROR: technology library is not complete, impossible to perform mapping" << std::endl;
        st.mapping_error = true;
        return false;
      }

      if ( node_data.same_match || node_data.map_refs[use_phase] > 0 )
      {
        if constexpr ( !ELA )
        {
          auto const& best_cut = cuts[index][node_data.best_cut[use_phase]];
          auto ctr = 0u;

          for ( auto const leaf : best_cut )
          {
            if ( ( node_data.phase[use_phase] >> ctr++ ) & 1 )
              node_match[leaf].map_refs[1]++;
            else
              node_match[leaf].map_refs[0]++;
          }
        }
        area += node_data.area[use_phase] + buffer_area_at( index, static_cast<uint8_t>( use_phase ) );
        if ( node_data.same_match && node_data.map_refs[use_phase ^ 1] > 0 )
        {
          if ( iteration < ps.area_flow_rounds )
          {
            ++node_data.map_refs[use_phase];
          }
          area += lib_inv_area;
          ++inv;
        }
      }

      /* invert the phase */
      use_phase = use_phase ^ 1;

      /* if both phases are implemented and used */
      if ( !node_data.same_match && node_data.map_refs[use_phase] > 0 )
      {
        if constexpr ( !ELA )
        {
          auto const& best_cut = cuts[index][node_data.best_cut[use_phase]];

          auto ctr = 0u;
          for ( auto const leaf : best_cut )
          {
            if ( ( node_data.phase[use_phase] >> ctr++ ) & 1 )
              node_match[leaf].map_refs[1]++;
            else
              node_match[leaf].map_refs[0]++;
          }
        }
        area += node_data.area[use_phase] + buffer_area_at( index, static_cast<uint8_t>( use_phase ) );
      }
    }

    ++iteration;

    if constexpr ( ELA )
    {
      /* ELA rounds refresh the load estimate too (the cover just changed under them), and must
       * restore the (arrival, required) consistency the moved delay model broke — stale required
       * times would send every heavily loaded node into the ignore-required fallback next round. */
      update_node_loads();
      recompute_arrivals_and_required();
      return true;
    }

    /* blend estimated references */
    float const coef = 1.0f / ( ( iteration + 1.0f ) * ( iteration + 1.0f ) );
    for ( auto i = 0u; i < ntk.size(); ++i )
    {
      node_match[i].est_refs[0] = std::max( 1.0f, coef * node_match[i].est_refs[0] + ( 1 - coef ) * node_match[i].map_refs[0] );
      node_match[i].est_refs[1] = std::max( 1.0f, coef * node_match[i].est_refs[1] + ( 1 - coef ) * node_match[i].map_refs[1] );
    }

    /* refresh the electrical load estimate from this round's committed cover (EXACT, not damped —
     * measured better than a damped blend on the ibex drive-legality benchmark), then restore the
     * (arrival, required) consistency the moved delay model broke. */
    update_node_loads();
    recompute_arrivals_and_required();

    /* Re-bind each covered node against the required times of the COMPLETE cover (area-delay
     * frontier), then restore consistency again: the re-selection changes arrivals, and leaving
     * the required times derived from the old ones is the same staleness trap the load
     * refinement above has to avoid. */
    if ( curve_budget > 0 )
    {
      reselect_from_curves();
      recompute_arrivals_and_required();
    }

    return true;
  }

  template<bool DO_AREA, bool ELA>
  bool set_mapping_refs_and_req()
  {
    for ( auto i = 0u; i < node_match.size(); ++i )
    {
      node_match[i].required[0] = node_match[i].required[1] = std::numeric_limits<float>::max();
    }

    /* compute the current worst delay and update the mapping refs */
    delay = 0.0f;
    ntk.foreach_po( [this]( auto s ) {
      const auto index = ntk.node_to_index( ntk.get_node( s ) );

      if ( ntk.is_complemented( s ) )
        delay = std::max( delay, node_match[index].arrival[1] );
      else
        delay = std::max( delay, node_match[index].arrival[0] );

      if constexpr ( !ELA )
      {
        if ( ntk.is_complemented( s ) )
          node_match[index].map_refs[1]++;
        else
          node_match[index].map_refs[0]++;
      }
    } );

    set_output_required_time( iteration == 0 );

    /* compute current area and update mapping refs in top-down order */
    area = 0.0f;
    inv = 0;
    for ( auto it = topo_order.rbegin(); it != topo_order.rend(); ++it )
    {
      const auto index = ntk.node_to_index( *it );
      auto& node_data = node_match[index];

      /* skip constants and PIs */
      if ( ntk.is_constant( *it ) )
      {
        if ( node_match[index].map_refs[0] || node_match[index].map_refs[1] )
        {
          /* if used and not available in the library launch a mapping error */
          if ( node_data.best_gate[0] == nullptr && node_data.best_gate[1] == nullptr )
          {
            std::cerr << "[e] MAP ERROR: technology library does not contain constant gates, impossible to perform mapping" << std::endl;
            st.mapping_error = true;
            return false;
          }
        }
        continue;
      }
      else if ( ntk.is_pi( *it ) )
      {
        if ( node_match[index].map_refs[1] > 0u )
        {
          /* Add inverter area over the negated fanins */
          area += lib_inv_area;
          ++inv;
        }
        continue;
      }

      /* continue if not referenced in the cover */
      if ( !node_match[index].map_refs[0] && !node_match[index].map_refs[1] )
        continue;

      /* don't touch box */
      if constexpr ( has_is_dont_touch_v<Ntk> )
      {
        if ( ntk.is_dont_touch( *it ) )
        {
          set_mapping_refs_dont_touch<ELA>( *it );
          continue;
        }
      }

      /* refine best matches with alternatives */
      if constexpr ( !DO_AREA )
      {
        if ( ps.use_match_alternatives )
          refine_best_matches( *it );
      }

      unsigned use_phase = node_data.best_gate[0] == nullptr ? 1u : 0u;
      if ( node_data.best_gate[use_phase] == nullptr )
      {
        /* Library is not complete, mapping is not possible */
        std::cerr << "[e] MAP ERROR: technology library is not complete, impossible to perform mapping" << std::endl;
        st.mapping_error = true;
        return false;
      }

      if ( node_data.same_match || node_data.map_refs[use_phase] > 0 )
      {
        if constexpr ( !ELA )
        {
          auto const& best_cut = cuts[index][node_data.best_cut[use_phase]];
          auto ctr = 0u;

          for ( auto const leaf : best_cut )
          {
            if ( ( node_data.phase[use_phase] >> ctr++ ) & 1 )
              node_match[leaf].map_refs[1]++;
            else
              node_match[leaf].map_refs[0]++;
          }
        }
        area += node_data.area[use_phase] + buffer_area_at( index, static_cast<uint8_t>( use_phase ) );
        if ( node_data.same_match && node_data.map_refs[use_phase ^ 1] > 0 )
        {
          if ( iteration < ps.area_flow_rounds )
          {
            ++node_data.map_refs[use_phase];
          }
          area += lib_inv_area;
          ++inv;
        }
      }

      /* invert the phase */
      use_phase = use_phase ^ 1;

      /* if both phases are implemented and used */
      if ( !node_data.same_match && node_data.map_refs[use_phase] > 0 )
      {
        if constexpr ( !ELA )
        {
          auto const& best_cut = cuts[index][node_data.best_cut[use_phase]];

          auto ctr = 0u;
          for ( auto const leaf : best_cut )
          {
            if ( ( node_data.phase[use_phase] >> ctr++ ) & 1 )
              node_match[leaf].map_refs[1]++;
            else
              node_match[leaf].map_refs[0]++;
          }
        }
        area += node_data.area[use_phase] + buffer_area_at( index, static_cast<uint8_t>( use_phase ) );
      }

      if ( !ps.area_oriented_mapping )
      {
        match_propagate_required( index );
      }
    }

    ++iteration;

    if constexpr ( ELA )
    {
      /* ELA rounds refresh the load estimate too (the cover just changed under them), and must
       * restore the (arrival, required) consistency the moved delay model broke — stale required
       * times would send every heavily loaded node into the ignore-required fallback next round. */
      update_node_loads();
      recompute_arrivals_and_required();
      return true;
    }

    /* blend estimated references */
    float const coef = 1.0f / ( ( iteration + 1.0f ) * ( iteration + 1.0f ) );
    for ( auto i = 0u; i < ntk.size(); ++i )
    {
      node_match[i].est_refs[0] = std::max( 1.0f, coef * node_match[i].est_refs[0] + ( 1 - coef ) * node_match[i].map_refs[0] );
      node_match[i].est_refs[1] = std::max( 1.0f, coef * node_match[i].est_refs[1] + ( 1 - coef ) * node_match[i].map_refs[1] );
    }

    /* refresh the electrical load estimate from this round's committed cover (EXACT, not damped —
     * measured better than a damped blend on the ibex drive-legality benchmark), then restore the
     * (arrival, required) consistency the moved delay model broke. */
    update_node_loads();
    recompute_arrivals_and_required();

    /* Re-bind each covered node against the required times of the COMPLETE cover (area-delay
     * frontier), then restore consistency again: the re-selection changes arrivals, and leaving
     * the required times derived from the old ones is the same staleness trap the load
     * refinement above has to avoid. */
    if ( curve_budget > 0 )
    {
      reselect_from_curves();
      recompute_arrivals_and_required();
    }

    return true;
  }

  template<bool ELA>
  inline void set_mapping_refs_dont_touch( node<Ntk> const& n )
  {
    if constexpr ( !ELA )
    {
      /* reference node */
      ntk.foreach_fanin( n, [&]( auto const& f ) {
        uint32_t leaf = ntk.node_to_index( ntk.get_node( f ) );
        uint8_t phase = ntk.is_complemented( f ) ? 1 : 0;
        node_match[leaf].map_refs[phase]++;
      } );
    }

    const auto index = ntk.node_to_index( n );

    if constexpr ( has_has_binding_v<Ntk> )
    {
      /* increase area */
      area += node_match[index].area[0];
      if ( node_match[index].map_refs[1] )
      {
        if ( iteration < ps.area_flow_rounds )
        {
          ++node_match[index].map_refs[0];
        }
        area += lib_inv_area;
        ++inv;
      }
    }
  }

  void set_output_required_time( bool warning )
  {
    double required = delay;
    /* relax delay constraints */
    if ( iteration == 0 && ps.required_time == 0.0f && ps.required_times.empty() && ps.relax_required > 0.0f )
    {
      required *= ( 100.0 + ps.relax_required ) / 100.0;
    }

    /* Global target time constraint */
    if ( ps.required_times.empty() )
    {
      if ( ps.required_time != 0.0f )
      {
        if ( ps.required_time < delay - epsilon )
        {
          if ( warning )
            std::cerr << fmt::format( "[i] MAP WARNING: cannot meet the target required time of {:.2f}", ps.required_time ) << std::endl;
        }
        else
        {
          required = ps.required_time;
        }
      }

      /* set the required time at POs */
      ntk.foreach_po( [&]( auto const& s ) {
        const auto index = ntk.node_to_index( ntk.get_node( s ) );
        if ( ntk.is_complemented( s ) )
          node_match[index].required[1] = required;
        else
          node_match[index].required[0] = required;
      } );

      return;
    }

    /* Output-specific target time constraint */
    ntk.foreach_po( [&]( auto const& s, uint32_t i ) {
      const auto index = ntk.node_to_index( ntk.get_node( s ) );
      uint8_t phase = ntk.is_complemented( s ) ? 1 : 0;
      if ( node_match[index].arrival[phase] > ps.required_times[i] + epsilon )
      {
        /* maintain the same delay */
        node_match[index].required[phase] = node_match[index].arrival[phase];
        if ( warning )
          std::cerr << fmt::format( "[i] MAP WARNING: cannot meet the target required time of {:.2f} at output {}", ps.required_times[i], i ) << std::endl;
      }
      else
      {
        node_match[index].required[phase] = ps.required_times[i];
      }
    } );
  }

  void compute_required_time( bool exit_early = false )
  {
    for ( auto i = 0u; i < node_match.size(); ++i )
    {
      node_match[i].required[0] = node_match[i].required[1] = std::numeric_limits<float>::max();
    }

    /* return if mapping is area oriented */
    if ( ps.area_oriented_mapping )
      return;

    set_output_required_time( iteration == 1 );

    if ( exit_early )
      return;

    /* propagate required time to the PIs */
    for ( auto it = topo_order.rbegin(); it != topo_order.rend(); ++it )
    {
      if ( ntk.is_pi( *it ) || ntk.is_constant( *it ) )
        break;

      const auto index = ntk.node_to_index( *it );

      if ( !node_match[index].map_refs[0] && !node_match[index].map_refs[1] )
        continue;

      match_propagate_required( index );
    }
  }

  void propagate_arrival_times()
  {
    /* the walk below derives arrivals, so the load estimate must be refreshed BEFORE it — a
     * refresh after the walk would leave the arrivals it just computed stale against the loads */
    update_node_loads();

    area = 0.0f;
    inv = 0;
    for ( auto const& n : topo_order )
    {
      auto index = ntk.node_to_index( n );
      auto& node_data = node_match[index];

      /* measure area */
      if ( ntk.is_constant( n ) )
      {
        continue;
      }
      else if ( ntk.is_pi( n ) )
      {
        if ( node_data.map_refs[1] > 0u )
        {
          /* Add inverter area over the negated fanins */
          area += lib_inv_area;
          ++inv;
        }
        continue;
      }

      /* reset required time */
      node_data.required[0] = std::numeric_limits<float>::max();
      node_data.required[1] = std::numeric_limits<float>::max();

      /* don't touch box */
      if constexpr ( has_is_dont_touch_v<Ntk> )
      {
        node<Ntk> n = ntk.index_to_node( index );
        if ( ntk.is_dont_touch( n ) )
        {
          if constexpr ( has_has_binding_v<Ntk> )
          {
            propagate_data_forward_white_box( n );
            if ( node_match[index].map_refs[0] || node_match[index].map_refs[1] )
              area += node_data.area[0];
            if ( node_data.map_refs[1] )
            {
              area += lib_inv_area;
              ++inv;
            }
          }
          continue;
        }
      }

      uint8_t use_phase = node_data.best_gate[0] != nullptr ? 0 : 1;

      /* compute arrival of use_phase */
      supergate<NInputs> const* best_gate = node_data.best_gate[use_phase];
      double worst_arrival = 0;
      uint16_t best_phase = node_data.phase[use_phase];
      auto ctr = 0u;
      for ( auto l : cuts[index][node_data.best_cut[use_phase]] )
      {
        double arrival_pin = node_match[l].arrival[( best_phase >> ctr ) & 1] + gate_delay( best_gate, ctr, index, use_phase );
        worst_arrival = std::max( worst_arrival, arrival_pin );
        ++ctr;
      }

      /* through the optional output buffer if one is chosen (cover_buffer) */
      worst_arrival += buffer_delay_at( index, use_phase );
      node_data.arrival[use_phase] = worst_arrival;

      /* compute area */
      if ( node_data.map_refs[use_phase] > 0 || ( node_data.same_match && ( node_match[index].map_refs[0] || node_match[index].map_refs[1] ) ) )
      {
        area += node_data.area[use_phase];
        if ( buffer_stage && node_buffers[index][use_phase].active )
          area += node_buffers[index][use_phase].drive.area;
        if ( node_data.same_match && node_data.map_refs[use_phase ^ 1] > 0 )
        {
          area += lib_inv_area;
          ++inv;
        }
      }

      /* compute arrival of the other phase */
      use_phase ^= 1;
      if ( node_data.same_match )
      {
        node_data.arrival[use_phase] = worst_arrival + inv_delay_at( index, use_phase );
        continue;
      }

      assert( node_data.best_gate[use_phase] != nullptr );

      best_gate = node_data.best_gate[use_phase];
      worst_arrival = 0;
      best_phase = node_data.phase[use_phase];
      ctr = 0u;
      for ( auto l : cuts[index][node_data.best_cut[use_phase]] )
      {
        double arrival_pin = node_match[l].arrival[( best_phase >> ctr ) & 1] + gate_delay( best_gate, ctr, index, use_phase );
        worst_arrival = std::max( worst_arrival, arrival_pin );
        ++ctr;
      }

      node_data.arrival[use_phase] = worst_arrival + buffer_delay_at( index, use_phase );

      if ( node_data.map_refs[use_phase] > 0 )
      {
        area += node_data.area[use_phase];
        if ( buffer_stage && node_buffers[index][use_phase].active )
          area += node_buffers[index][use_phase].drive.area;
      }
    }

    /* compute the current worst delay */
    delay = 0.0f;
    ntk.foreach_po( [this]( auto s ) {
      const auto index = ntk.node_to_index( ntk.get_node( s ) );

      if ( ntk.is_complemented( s ) )
        delay = std::max( delay, node_match[index].arrival[1] );
      else
        delay = std::max( delay, node_match[index].arrival[0] );
    } );

    /* return if mapping is area oriented */
    ++iteration;
    if ( ps.area_oriented_mapping )
      return;

    /* set the required time at POs */
    ntk.foreach_po( [&]( auto const& s ) {
      const auto index = ntk.node_to_index( ntk.get_node( s ) );
      if ( ntk.is_complemented( s ) )
        node_match[index].required[1] = delay;
      else
        node_match[index].required[0] = delay;
    } );
  }

  void propagate_arrival_node( node<Ntk> const& n )
  {
    uint32_t index = ntk.node_to_index( n );
    auto& node_data = node_match[index];
    uint8_t use_phase = node_data.best_gate[0] != nullptr ? 0 : 1;

    /* compute arrival of use_phase */
    supergate<NInputs> const* best_gate = node_data.best_gate[use_phase];
    double worst_arrival = 0;
    uint16_t best_phase = node_data.phase[use_phase];
    auto ctr = 0u;
    for ( auto l : cuts[index][node_data.best_cut[use_phase]] )
    {
      double arrival_pin = node_match[l].arrival[( best_phase >> ctr ) & 1] + gate_delay( best_gate, ctr, index, use_phase );
      worst_arrival = std::max( worst_arrival, arrival_pin );
      ++ctr;
    }
    /* the stored arrival is what the SINKS see: through the optional output buffer if one is
     * chosen (cover_buffer) */
    worst_arrival += buffer_delay_at( index, use_phase );
    node_data.arrival[use_phase] = worst_arrival;

    /* compute arrival of the other phase */
    use_phase ^= 1;
    if ( node_data.same_match )
    {
      node_data.arrival[use_phase] = worst_arrival + inv_delay_at( index, use_phase );
      return;
    }

    assert( node_data.best_gate[0] != nullptr );

    best_gate = node_data.best_gate[use_phase];
    worst_arrival = 0;
    best_phase = node_data.phase[use_phase];
    ctr = 0u;
    for ( auto l : cuts[index][node_data.best_cut[use_phase]] )
    {
      double arrival_pin = node_match[l].arrival[( best_phase >> ctr ) & 1] + gate_delay( best_gate, ctr, index, use_phase );
      worst_arrival = std::max( worst_arrival, arrival_pin );
      ++ctr;
    }

    node_data.arrival[use_phase] = worst_arrival + buffer_delay_at( index, use_phase );
  }

  template<bool DO_AREA>
  void match_phase( node<Ntk> const& n, uint8_t phase )
  {
    auto index = ntk.node_to_index( n );
    auto& node_data = node_match[index];
    uint32_t cut_index = 0u;

    node_data.best_gate[phase] = nullptr;
    node_data.arrival[phase] = std::numeric_limits<float>::max();
    node_data.flows[phase] = std::numeric_limits<float>::max();
    node_data.area[phase] = std::numeric_limits<float>::max();
    uint32_t best_size = UINT32_MAX;

    best_gate_emap<NInputs>& gA = node_data.best_alternative[phase];
    gA.gate = nullptr;
    gA.arrival = std::numeric_limits<float>::max();
    gA.flow = std::numeric_limits<float>::max();
    uint32_t best_sizeA = UINT32_MAX;
    /* electrical model: whether the current best can legally drive this signal's load (the
     * alternative gA keeps the pure cost comparison — it competes under required-time checks
     * downstream, and its arrival already carries the load term) */
    bool best_legal = false;
    /* the incumbent's RANKING cost. It shadows node_data.arrival[phase] in the comparison only:
     * the stored arrival stays the true one, so required times, area recovery and the cover's
     * own bookkeeping keep reconciling against a single honest number while the choice between
     * candidates accounts for the load each of them imposes upstream (ADR-0050). Identical to
     * the arrival with the curves off, which is what keeps the default cover byte-identical. */
    double best_rank = std::numeric_limits<float>::max();
    /* the frontier is rebuilt with the match: a point from a previous round would carry a cut
     * index that no longer means the same thing */
    if ( curve_budget > 0 && index < curves.size() )
      curves[index][phase].size = 0;

    /* unmap multioutput */
    node_data.multioutput_match[phase] = false;

    /* foreach cut */
    for ( auto& cut : cuts[index] )
    {
      /* trivial cuts or not matched cuts */
      if ( ( *cut )->ignore )
      {
        ++cut_index;
        continue;
      }

      auto const& supergates = ( *cut )->supergates;
      auto const negation = ( *cut )->negations[phase];

      if ( supergates[phase] == nullptr )
      {
        ++cut_index;
        continue;
      }

      /* match each gate and take the best one */
      for ( auto const& gate : *supergates[phase] )
      {
        uint16_t gate_polarity = gate.polarity ^ negation;
        double worst_arrival = 0.0f;
        /* the same arrival with each leaf charged for the capacitance THIS candidate puts on it
         * (ADR-0050) — the ranking cost, never a stored time */
        double rank_arrival = 0.0f;
        double worst_arrivalA = 0.0f;
        float area_local = gate.area;
        float area_localA = gate.area;

        auto ctr = 0u;
        for ( auto l : *cut )
        {
          uint8_t leaf_phase = ( gate_polarity >> ctr ) & 1;

          double arrival_pinA = node_match[l].best_alternative[leaf_phase].arrival + gate_delay( &gate, ctr, index, phase );
          worst_arrivalA = std::max( worst_arrivalA, arrival_pinA );

          // if constexpr ( DO_AREA )
          // {
          //   if ( worst_arrivalA > node_data.required[phase] + epsilon || worst_arrivalA >= std::numeric_limits<float>::max() )
          //     break;
          // }

          double arrival_pin = node_match[l].arrival[leaf_phase] + gate_delay( &gate, ctr, index, phase );
          worst_arrival = std::max( worst_arrival, arrival_pin );
          rank_arrival = std::max( rank_arrival, arrival_pin + leaf_load_penalty( l, leaf_phase, &gate, ctr ) );

          area_local += node_match[l].flows[leaf_phase];
          area_localA += node_match[l].best_alternative[leaf_phase].flow;
          ++ctr;
        }

        bool skip = false;
        if constexpr ( DO_AREA )
        {
          if ( ctr < cut->size() )
            continue;
          if ( worst_arrival > node_data.required[phase] + epsilon || worst_arrival >= std::numeric_limits<float>::max() )
            skip = true;
        }

        /* every feasible candidate is offered to the node's frontier, so a post-cover
         * re-selection has real alternatives to choose from (see reselect_from_curves) */
        if ( !skip )
          curve_offer( index, phase, &gate, worst_arrival, gate.area, area_local, gate_polarity,
                       cut_index, cut->size() );

        /* electrical model: a drive that can legally carry the signal's load outranks one
         * that cannot; among equals the usual cost comparison decides. Inert with the model
         * off (every candidate is trivially legal, so pick == compare_map). */
        bool pick;
        bool const cand_legal = gate_load_legal( &gate, index, phase );
        if ( ps.electrical_model && cand_legal != best_legal )
          pick = cand_legal;
        else
          pick = compare_map<DO_AREA>( rank_arrival, best_rank, area_local, node_data.flows[phase], cut->size(), best_size );
        if ( !skip && pick )
        {
          node_data.best_gate[phase] = &gate;
          node_data.arrival[phase] = worst_arrival;
          best_rank = rank_arrival;
          node_data.flows[phase] = area_local;
          node_data.best_cut[phase] = cut_index;
          node_data.area[phase] = gate.area;
          node_data.phase[phase] = gate_polarity;
          best_size = cut->size();
          best_legal = cand_legal;
        }

        /* compute the alternative */
        if ( compare_map<!DO_AREA>( worst_arrivalA, gA.arrival, area_localA, gA.flow, cut->size(), best_sizeA ) )
        {
          gA.gate = &gate;
          gA.arrival = worst_arrivalA;
          gA.area = gate.area;
          gA.flow = area_localA;
          gA.phase = gate_polarity;
          gA.cut = cut_index;
          best_sizeA = cut->size();
          gA.size = cut->size();
        }
      }

      ++cut_index;
    }

    if constexpr ( DO_AREA )
    {
      if ( node_data.best_gate[phase] == nullptr )
      {
        /* no candidate met the required time (multi-output squeeze, see match_phase_relaxed):
         * re-match ignoring required times to keep the node covered */
        match_phase_relaxed<false, false>( n, phase );
      }
    }
  }

  /* Fallback single-output match that ignores required times.
   *
   * With multi-output mapping, a node's required time can be set by a multi-output pin that is
   * faster than ANY single-output implementation of that node. This is possible when the
   * multi-output cell has no delay-matching single-output cover in the library (the library
   * loader warns about such cells). The required-filtered searches in match_phase /
   * match_phase_exact then find no candidate at all, leaving a referenced node without a match
   * in either phase — and the later required-time propagation would dereference a null gate.
   * This fallback restores the cover-completeness invariant (every referenced node has a valid
   * match) by picking the fastest implementable match, accepting the (already warned) local
   * required-time violation. The multi-output re-match runs after the phase matches, so when the
   * multi-output gate is still usable it supersedes this fallback and the violation disappears. */
  template<bool ExactArea, bool SwitchActivity>
  void match_phase_relaxed( node<Ntk> const& n, uint8_t phase )
  {
    double best_arrival = std::numeric_limits<double>::max();
    double best_rank = std::numeric_limits<double>::max(); /* ranking cost only; see match_phase */
    float best_flow = std::numeric_limits<float>::max();
    float best_area = std::numeric_limits<float>::max();
    uint32_t best_size = UINT32_MAX;
    uint32_t best_cut = 0u;
    uint16_t best_phase = 0u;
    supergate<NInputs> const* best_gate = nullptr;
    uint32_t cut_index = 0u;
    auto index = ntk.node_to_index( n );

    auto& node_data = node_match[index];

    /* foreach cut */
    for ( auto& cut : cuts[index] )
    {
      /* trivial cuts or not matched cuts */
      if ( ( *cut )->ignore )
      {
        ++cut_index;
        continue;
      }

      auto const& supergates = ( *cut )->supergates;
      auto const negation = ( *cut )->negations[phase];

      if ( supergates[phase] == nullptr )
      {
        ++cut_index;
        continue;
      }

      /* match each gate and take the fastest one (smallest required-time violation) */
      for ( auto const& gate : *supergates[phase] )
      {
        uint16_t gate_polarity = gate.polarity ^ negation;
        double worst_arrival = 0.0f;
        double rank_arrival = 0.0f; /* ranking cost only (ADR-0050); see match_phase */
        float flow = gate.area;

        auto ctr = 0u;
        for ( auto l : *cut )
        {
          uint8_t leaf_phase = ( gate_polarity >> ctr ) & 1;
          double const arrival_pin = node_match[l].arrival[leaf_phase] + gate_delay( &gate, ctr, index, phase );
          worst_arrival = std::max( worst_arrival, arrival_pin );
          rank_arrival = std::max( rank_arrival, arrival_pin + leaf_load_penalty( l, leaf_phase, &gate, ctr ) );
          if constexpr ( !ExactArea )
            flow += node_match[l].flows[leaf_phase];
          ++ctr;
        }

        if ( worst_arrival >= std::numeric_limits<float>::max() )
          continue;

        if constexpr ( ExactArea )
        {
          node_data.phase[phase] = gate_polarity;
          node_data.area[phase] = gate.area;
          flow = cut_measure_mffc<SwitchActivity>( *cut, n, phase );
        }

        if ( compare_map<false>( rank_arrival, best_rank, flow, best_flow, cut->size(), best_size ) )
        {
          best_rank = rank_arrival;
          best_arrival = worst_arrival;
          best_flow = flow;
          best_area = gate.area;
          best_size = cut->size();
          best_cut = cut_index;
          best_phase = gate_polarity;
          best_gate = &gate;
        }
      }

      ++cut_index;
    }

    if ( best_gate == nullptr )
      return; /* no single-output implementation at all: leave the (pre-existing) unmatched state */

    node_data.flows[phase] = best_flow;
    node_data.arrival[phase] = best_arrival;
    node_data.area[phase] = best_area;
    node_data.best_cut[phase] = best_cut;
    node_data.phase[phase] = best_phase;
    node_data.best_gate[phase] = best_gate;
  }

  template<bool SwitchActivity>
  void match_phase_exact( node<Ntk> const& n, uint8_t phase )
  {
    double best_arrival = std::numeric_limits<float>::max();
    double best_rank = std::numeric_limits<float>::max(); /* ranking cost only; see match_phase */
    float best_exact_area = std::numeric_limits<float>::max();
    float best_area = std::numeric_limits<float>::max();
    uint32_t best_size = UINT32_MAX;
    uint8_t best_cut = 0u;
    uint16_t best_phase = 0u;
    uint8_t cut_index = 0u;
    bool best_legal = false;
    auto index = ntk.node_to_index( n );

    auto& node_data = node_match[index];
    supergate<NInputs> const* best_gate = node_data.best_gate[phase];

    /* unmap multioutput */
    if ( node_data.multioutput_match[phase] )
    {
      /* dereference multi-output */
      if ( !node_data.same_match && best_gate != nullptr && node_data.map_refs[phase] )
      {
        auto const& cut = multi_cut_set[node_data.best_cut[phase]][0];
        cut_deref<SwitchActivity>( cut, n, phase );
      }
      best_gate = nullptr;
      node_data.multioutput_match[phase] = false;
    }

    /* recompute best match info */
    if ( best_gate != nullptr )
    {
      /* if cut is implemented, remove it from the cover */
      if ( !node_data.same_match && node_data.map_refs[phase] )
      {
        auto const& cut = cuts[index][node_data.best_cut[phase]];
        cut_deref<SwitchActivity>( cut, n, phase );
      }
    }

    /* the search below starts from scratch: a failed search must not keep the previous gate
     * paired with an unrelated best_cut (index 0) — it falls back to a relaxed re-match instead */
    best_gate = nullptr;

    /* foreach cut */
    for ( auto& cut : cuts[index] )
    {
      /* trivial cuts or not matched cuts */
      if ( ( *cut )->ignore )
      {
        ++cut_index;
        continue;
      }

      auto const& supergates = ( *cut )->supergates;
      auto const negation = ( *cut )->negations[phase];

      if ( supergates[phase] == nullptr )
      {
        ++cut_index;
        continue;
      }

      /* match each gate and take the best one */
      for ( auto const& gate : *supergates[phase] )
      {
        uint16_t gate_polarity = gate.polarity ^ negation;
        double worst_arrival = 0.0f;
        double rank_arrival = 0.0f; /* ranking cost only (ADR-0050); see match_phase */

        auto ctr = 0u;
        for ( auto l : *cut )
        {
          uint8_t const leaf_phase = ( gate_polarity >> ctr ) & 1;
          double arrival_pin = node_match[l].arrival[leaf_phase] + gate_delay( &gate, ctr, index, phase );
          worst_arrival = std::max( worst_arrival, arrival_pin );
          rank_arrival = std::max( rank_arrival, arrival_pin + leaf_load_penalty( l, leaf_phase, &gate, ctr ) );
          ++ctr;
        }

        if ( worst_arrival > node_data.required[phase] + epsilon || worst_arrival >= std::numeric_limits<float>::max() )
          continue;

        node_data.phase[phase] = gate_polarity;
        node_data.area[phase] = gate.area;
        float area_exact = cut_measure_mffc<SwitchActivity>( *cut, n, phase );

        /* electrical model: legal-drive preference (see match_phase) */
        bool pick;
        bool const cand_legal = gate_load_legal( &gate, index, phase );
        if ( ps.electrical_model && cand_legal != best_legal )
          pick = cand_legal;
        else
          pick = compare_map<true>( rank_arrival, best_rank, area_exact, best_exact_area, cut->size(), best_size );
        if ( pick )
        {
          best_arrival = worst_arrival;
          best_rank = rank_arrival;
          best_exact_area = area_exact;
          best_area = gate.area;
          best_size = cut->size();
          best_cut = cut_index;
          best_phase = gate_polarity;
          best_gate = &gate;
          best_legal = cand_legal;
        }
      }

      ++cut_index;
    }

    if ( best_gate == nullptr )
    {
      /* no candidate met the required time (multi-output squeeze, see match_phase_relaxed):
       * re-match ignoring required times to keep the node covered */
      node_data.best_gate[phase] = nullptr;
      node_data.arrival[phase] = best_arrival;
      node_data.flows[phase] = best_exact_area;
      match_phase_relaxed<true, SwitchActivity>( n, phase );
      if ( node_data.best_gate[phase] != nullptr && !node_data.same_match && node_data.map_refs[phase] )
      {
        cut_ref<SwitchActivity>( cuts[index][node_data.best_cut[phase]], n, phase );
      }
      return;
    }

    node_data.flows[phase] = best_exact_area;
    node_data.arrival[phase] = best_arrival;
    node_data.area[phase] = best_area;
    node_data.best_cut[phase] = best_cut;
    node_data.phase[phase] = best_phase;
    node_data.best_gate[phase] = best_gate;

    if ( !node_data.same_match && node_data.map_refs[phase] )
    {
      best_exact_area = cut_ref<SwitchActivity>( cuts[index][best_cut], n, phase );
    }
  }

  template<bool DO_AREA, bool ELA, bool SwitchActivity = false>
  void match_drop_phase( node<Ntk> const& n )
  {
    auto index = ntk.node_to_index( n );
    auto& node_data = node_match[index];

    /* compute arrival adding an inverter to the other match phase */
    double worst_arrival_npos = node_data.arrival[1] + inv_delay_at( index, 0 );
    double worst_arrival_nneg = node_data.arrival[0] + inv_delay_at( index, 1 );
    bool use_zero = false;
    bool use_one = false;

    /* only one phase is matched */
    if ( node_data.best_gate[0] == nullptr )
    {
      set_match_complemented_phase( index, 1, worst_arrival_npos );
      if constexpr ( ELA )
      {
        if ( node_data.map_refs[0] || node_data.map_refs[1] )
          cut_ref<false>( cuts[index][node_data.best_cut[1]], n, 1 );
      }
      return;
    }
    else if ( node_data.best_gate[1] == nullptr )
    {
      set_match_complemented_phase( index, 0, worst_arrival_nneg );
      if constexpr ( ELA )
      {
        if ( node_data.map_refs[0] || node_data.map_refs[1] )
          cut_ref<false>( cuts[index][node_data.best_cut[0]], n, 0 );
      }
      return;
    }

    /* try to use only one match to cover both phases */
    if constexpr ( !DO_AREA )
    {
      /* if arrival improves matching the other phase and inserting an inverter */
      if ( worst_arrival_npos < node_data.arrival[0] + epsilon )
      {
        use_one = true;
      }
      if ( worst_arrival_nneg < node_data.arrival[1] + epsilon )
      {
        use_zero = true;
      }
    }
    else
    {
      /* check if both phases + inverter meet the required time */
      use_zero = worst_arrival_nneg < ( node_data.required[1] + epsilon );
      use_one = worst_arrival_npos < ( node_data.required[0] + epsilon );
    }

    /* condition on not used phases, evaluate a substitution during exact area recovery */
    if constexpr ( ELA )
    {
      if ( node_data.map_refs[0] == 0 || node_data.map_refs[1] == 0 )
      {
        /* select the used match */
        auto phase = 0;
        auto nphase = 0;
        if ( node_data.map_refs[0] == 0 )
        {
          phase = 1;
          use_one = true;
          use_zero = false;
        }
        else
        {
          nphase = 1;
          use_one = false;
          use_zero = true;
        }
        /* select the not used match instead if it leads to area improvement and doesn't violate the required time */
        if ( node_data.arrival[nphase] + inv_delay_at( index, phase ) < node_data.required[phase] + epsilon )
        {
          auto size_phase = cuts[index][node_data.best_cut[phase]].size();
          auto size_nphase = cuts[index][node_data.best_cut[nphase]].size();

          if ( compare_map<DO_AREA>( node_data.arrival[nphase] + inv_delay_at( index, phase ), node_data.arrival[phase], node_data.flows[nphase] + lib_inv_area, node_data.flows[phase], size_nphase, size_phase ) )
          {
            /* invert the choice */
            use_zero = !use_zero;
            use_one = !use_one;
          }
        }
      }
    }

    if ( ( !use_zero && !use_one ) )
    {
      /* use both phases */
      node_data.flows[0] = node_data.flows[0] / node_data.est_refs[0];
      node_data.flows[1] = node_data.flows[1] / node_data.est_refs[1];
      node_data.same_match = false;
      return;
    }

    /* use area flow as a tiebreaker */
    if ( use_zero && use_one )
    {
      auto size_zero = cuts[index][node_data.best_cut[0]].size();
      auto size_one = cuts[index][node_data.best_cut[1]].size();

      if constexpr ( ELA )
      {
        if ( !node_data.same_match )
        {
          /* both phases were implemented --> evaluate substitution */
          cut_deref<false>( cuts[index][node_data.best_cut[0]], n, 0 );
          node_data.flows[1] = cut_deref<false>( cuts[index][node_data.best_cut[1]], n, 1 );
          node_data.flows[0] = cut_ref<false>( cuts[index][node_data.best_cut[0]], n, 0 );
          cut_ref<false>( cuts[index][node_data.best_cut[1]], n, 1 );
        }
        /* evaluate based on inverter cost */
        if constexpr ( !SwitchActivity )
        {
          use_zero = lib_inv_area < node_data.flows[1] + epsilon;
          use_one = lib_inv_area < node_data.flows[0] + epsilon;
        }

        if ( use_one && use_zero )
        {
          if ( compare_map<DO_AREA>( worst_arrival_nneg, worst_arrival_npos, node_data.flows[0], node_data.flows[1], size_zero, size_one ) )
            use_one = false;
          else
            use_zero = false;
        }
        else if ( !use_one && !use_zero && node_data.same_match )
        {
          node_data.same_match = false;
          cut_ref<false>( cuts[index][node_data.best_cut[0]], n, 0 );
          cut_ref<false>( cuts[index][node_data.best_cut[1]], n, 1 );
          return;
        }
      }
      else
      {
        /* compare flows by looking at the most convinient and referenced */
        if ( node_data.flows[0] / node_data.est_refs[0] + lib_inv_area < node_data.flows[1] / node_data.est_refs[1] + epsilon )
        {
          use_one = false;
        }
        else if ( node_data.flows[1] / node_data.est_refs[1] + lib_inv_area < node_data.flows[0] / node_data.est_refs[0] + epsilon )
        {
          use_zero = false;
        }
        else
        {
          /* delay the decision on what to keep --> wait for better estimations */
          node_data.flows[0] = node_data.flows[0] / node_data.est_refs[0];
          node_data.flows[1] = node_data.flows[1] / node_data.est_refs[1];
          node_data.same_match = false;
          return;
        }
      }
    }

    if ( use_zero )
    {
      if constexpr ( ELA )
      {
        /* set cut references */
        if ( !node_data.same_match )
        {
          /* dereference the negative phase cut if in use */
          if ( node_data.map_refs[1] > 0 )
            cut_deref<false>( cuts[index][node_data.best_cut[1]], n, 1 );
          /* reference the positive cut if not in use before */
          if ( node_data.map_refs[0] == 0 && node_data.map_refs[1] > 0 )
            cut_ref<false>( cuts[index][node_data.best_cut[0]], n, 0 );
        }
        else if ( node_data.map_refs[0] || node_data.map_refs[1] )
          cut_ref<false>( cuts[index][node_data.best_cut[0]], n, 0 );
      }
      set_match_complemented_phase( index, 0, worst_arrival_nneg );
    }
    else
    {
      if constexpr ( ELA )
      {
        /* set cut references */
        if ( !node_data.same_match )
        {
          /* dereference the positive phase cut if in use */
          if ( node_data.map_refs[0] > 0 )
            cut_deref<false>( cuts[index][node_data.best_cut[0]], n, 0 );
          /* reference the negative cut if not in use before */
          if ( node_data.map_refs[1] == 0 && node_data.map_refs[0] > 0 )
            cut_ref<false>( cuts[index][node_data.best_cut[1]], n, 1 );
        }
        else if ( node_data.map_refs[0] || node_data.map_refs[1] )
          cut_ref<false>( cuts[index][node_data.best_cut[1]], n, 1 );
      }
      set_match_complemented_phase( index, 1, worst_arrival_npos );
    }
  }

  inline void set_match_complemented_phase( uint32_t index, uint8_t phase, double worst_arrival_n )
  {
    auto& node_data = node_match[index];
    auto phase_n = phase ^ 1;
    node_data.same_match = true;
    node_data.best_gate[phase_n] = nullptr;
    node_data.best_cut[phase_n] = node_data.best_cut[phase];
    node_data.phase[phase_n] = node_data.phase[phase];
    node_data.arrival[phase_n] = worst_arrival_n;
    node_data.area[phase_n] = node_data.area[phase];
    node_data.flows[phase_n] = ( node_data.flows[phase] + lib_inv_area ) / node_data.est_refs[phase_n];
    node_data.flows[phase] = node_data.flows[phase] / node_data.est_refs[phase];
  }

  template<bool DO_AREA>
  inline void select_alternatives( node<Ntk> const& n )
  {
    if constexpr ( DO_AREA )
      return;

    if ( !ps.use_match_alternatives )
      return;

    auto index = ntk.node_to_index( n );
    auto& node_data = node_match[index];

    best_gate_emap<NInputs>& g0 = node_data.best_alternative[0];
    best_gate_emap<NInputs>& g1 = node_data.best_alternative[1];
    float g0flow = g0.flow / node_data.est_refs[0];
    float g1flow = g1.flow / node_data.est_refs[1];

    /* process for best area */ /* removed check on required since this is executed only during a delay pass */
    if ( g0.gate != nullptr && g0flow + lib_inv_area < g1flow + epsilon )
    {
      g1 = g0;
      g1.gate = nullptr;
      g1.arrival += inv_delay_at( index, 1 );
      g1.flow = ( g1.flow + lib_inv_area ) / node_data.est_refs[1];
      g0.flow = g0flow;
      return;
    }
    else if ( g1.gate != nullptr && g1flow + lib_inv_area < g0flow + epsilon )
    {
      g0 = g1;
      g0.gate = nullptr;
      g0.arrival += inv_delay_at( index, 0 );
      g0.flow = ( g0.flow + lib_inv_area ) / node_data.est_refs[0];
      g1.flow = g1flow;
      return;
    }

    g0.flow = g0flow;
    g1.flow = g1flow;
  }

  inline void refine_best_matches( node<Ntk> const& n )
  {
    auto index = ntk.node_to_index( n );
    auto& node_data = node_match[index];

    /* evaluate to change the best matches with the best alternative */
    best_gate_emap<NInputs>& g0 = node_data.best_alternative[0];
    best_gate_emap<NInputs>& g1 = node_data.best_alternative[1];

    if ( node_data.map_refs[0] && node_data.map_refs[1] )
    {
      if ( node_data.same_match )
      {
        /* pick best implementation between the two alternatives */
        unsigned best_match_phase = node_data.best_gate[0] == nullptr ? 1 : 0;
        unsigned use_phase = g0.gate == nullptr ? 1 : 0;
        if ( g0.gate != nullptr && g1.gate != nullptr )
        {
          if ( g0.arrival > node_data.required[0] + epsilon || g1.arrival > node_data.required[1] + epsilon )
            return;

          refine_best_matches_copy_refinement( n, 0, false );
          refine_best_matches_copy_refinement( n, 1, false );
          node_data.same_match = false;
          return;
        }
        else
        {
          best_gate_emap<NInputs>& gUse = node_data.best_alternative[use_phase];
          if ( gUse.arrival > node_data.required[use_phase] + epsilon || gUse.arrival + inv_delay_at( index, use_phase ^ 1 ) > node_data.required[use_phase ^ 1] + epsilon )
          {
            return;
          }
          refine_best_matches_copy_refinement( n, use_phase, true );
          return;
        }
      }
      else
      {
        /* not same match: evaluate both zero and one phase */
        if ( g0.gate != nullptr && g0.arrival < node_data.required[0] + epsilon )
        {
          node_data.same_match = false;
          refine_best_matches_copy_refinement( n, 0, g1.gate == nullptr && g0.arrival + inv_delay_at( index, 1 ) < node_data.required[1] + epsilon );
        }
        if ( g1.gate != nullptr && g1.arrival < node_data.required[1] + epsilon )
        {
          node_data.same_match = false;
          refine_best_matches_copy_refinement( n, 1, g0.gate == nullptr && g1.arrival + inv_delay_at( index, 0 ) < node_data.required[0] + epsilon );
        }
      }
    }
    else if ( node_data.map_refs[0] )
    {
      if ( g0.gate != nullptr && g0.arrival < node_data.required[0] + epsilon )
      {
        node_data.same_match = false;
        refine_best_matches_copy_refinement( n, 0, false );
      }
      else if ( g0.gate == nullptr && g1.arrival + inv_delay_at( index, 0 ) < node_data.required[0] + epsilon )
      {
        refine_best_matches_copy_refinement( n, 1, true );
      }
    }
    else
    {
      if ( g1.gate != nullptr && g1.arrival < node_data.required[1] + epsilon )
      {
        node_data.same_match = false;
        refine_best_matches_copy_refinement( n, 1, false );
      }
      else if ( g1.gate == nullptr && g0.arrival + inv_delay_at( index, 1 ) < node_data.required[1] + epsilon )
      {
        refine_best_matches_copy_refinement( n, 0, true );
      }
    }
  }

  inline void refine_best_matches_copy_refinement( node<Ntk> const& n, unsigned phase, bool both_phases )
  {
    auto index = ntk.node_to_index( n );
    auto& node_data = node_match[index];
    best_gate_emap<NInputs>& bg = node_data.best_alternative[phase];

    node_data.best_gate[phase] = bg.gate;
    node_data.phase[phase] = bg.phase;
    node_data.best_cut[phase] = bg.cut;
    node_data.arrival[phase] = bg.arrival;
    node_data.area[phase] = bg.area;
    node_data.flows[phase] = bg.flow;

    if ( !both_phases )
      return;

    node_data.same_match = true;
    phase ^= 1;
    node_data.best_gate[phase] = nullptr;
    node_data.phase[phase] = bg.phase;
    node_data.best_cut[phase] = bg.cut;
    node_data.arrival[phase] = bg.arrival + inv_delay_at( index, phase );
    node_data.area[phase] = bg.area;
    node_data.flows[phase] = ( bg.flow * node_data.est_refs[phase ^ 1] + lib_inv_area ) / node_data.est_refs[phase];
  }

  bool initialize_box( node<Ntk> const& n )
  {
    uint32_t index = ntk.node_to_index( n );

    if ( cuts[index].size() == 0 )
      add_unit_cut( index );

    auto& node_data = node_match[index];
    node_data.same_match = true;

    /* if it has mapping data propagate the delays and measure the data */
    if constexpr ( has_has_binding_v<Ntk> )
    {
      propagate_data_forward_white_box( n );
      return false;
    }

    /* consider as a black box */
    node_data.flows[0] = 0.0f;
    node_data.flows[1] = lib_inv_area / node_data.est_ref[1];
    node_data.arrival[0] = 0.0f;
    node_data.arrival[1] = lib_inv_delay;
    node_data.area[0] = node_data.area[1] = 0;

    return true;
  }

  void propagate_data_forward_white_box( node<Ntk> const& n )
  {
    uint32_t index = ntk.node_to_index( n );
    auto& node_data = node_match[index];
    auto const& gate = ntk.get_binding( n );

    /* propagate arrival time */
    double arrival = 0;
    ntk.foreach_fanin( n, [&]( auto const& f, auto i ) {
      uint32_t f_index = ntk.node_to_index( ntk.get_node( f ) );
      uint8_t phase = ntk.is_complemented( f ) ? 1 : 0;
      double propagation_delay = std::max( gate.pins[i].rise_block_delay, gate.pins[i].fall_block_delay );
      arrival = std::max( arrival, node_match[f_index].arrival[phase] + propagation_delay );
    } );

    /* set data */
    node_data.arrival[0] = arrival;
    node_data.arrival[1] = arrival + lib_inv_delay;
    node_data.area[0] = node_data.area[1] = gate.area;
    node_data.flows[1] = ( node_data.flows[0] + lib_inv_area ) / node_data.est_refs[1];
    node_data.flows[0] = node_data.area[0] / node_data.est_refs[0];
  }

  void propagate_data_backward_white_box( node<Ntk> const& n )
  {
    uint32_t index = ntk.node_to_index( n );
    auto& node_data = node_match[index];
    auto const& gate = ntk.get_binding( n );

    assert( node_data.map_refs[0] || node_data.map_refs[1] );

    /* propagate required time over the output inverter if present */
    if ( node_data.map_refs[1] > 0 )
    {
      node_data.required[0] = std::min( node_data.required[0], node_data.required[1] - lib_inv_delay );
    }

    if ( node_data.map_refs[0] )
      assert( node_data.arrival[0] < node_data.required[0] + epsilon );
    if ( node_data.map_refs[1] )
      assert( node_data.arrival[1] < node_data.required[1] + epsilon );

    ntk.foreach_fanin( n, [&]( auto const& f, auto i ) {
      uint32_t f_index = ntk.node_to_index( ntk.get_node( f ) );
      uint8_t phase = ntk.is_complemented( f ) ? 1 : 0;
      double propagation_delay = std::max( gate.pins[i].rise_block_delay, gate.pins[i].fall_block_delay );
      node_match[f_index].required[phase] = std::min( node_match[f_index].required[phase], node_data.required[0] - propagation_delay );
    } );
  }

  void match_constants( uint32_t index )
  {
    auto& node_data = node_match[index];

    kitty::static_truth_table<6> zero_tt;
    auto const supergates_zero = library.get_supergates( zero_tt );
    auto const supergates_one = library.get_supergates( ~zero_tt );

    /* Not available in the library */
    if ( supergates_zero == nullptr && supergates_one == nullptr )
    {
      return;
    }
    /* if only one is available, the other is obtained using an inverter */
    if ( supergates_zero != nullptr )
    {
      node_data.best_gate[0] = &( ( *supergates_zero )[0] );
      node_data.arrival[0] = node_data.best_gate[0]->tdelay[0];
      node_data.area[0] = node_data.best_gate[0]->area;
      node_data.phase[0] = 0;
    }
    if ( supergates_one != nullptr )
    {
      node_data.best_gate[1] = &( ( *supergates_one )[0] );
      node_data.arrival[1] = node_data.best_gate[1]->tdelay[0];
      node_data.area[1] = node_data.best_gate[1]->area;
      node_data.phase[1] = 0;
    }
    else
    {
      node_data.same_match = true;
      node_data.arrival[1] = node_data.arrival[0] + lib_inv_delay;
      node_data.area[1] = node_data.area[0] + lib_inv_area;
      node_data.phase[1] = 1;
    }
    if ( supergates_zero == nullptr )
    {
      node_data.same_match = true;
      node_data.arrival[0] = node_data.arrival[1] + lib_inv_delay;
      node_data.area[0] = node_data.area[1] + lib_inv_area;
      node_data.phase[0] = 1;
    }
  }

  template<bool DO_AREA>
  bool match_multioutput( node<Ntk> const& n )
  {
    /* extract outputs tuple */
    uint32_t index = ntk.node_to_index( n );
    multi_match_t const& tuple_data = multi_node_match[node_tuple_match[index].index][0];

    /* get the cut */
    auto const& cut0 = cuts[tuple_data[0].node_index][tuple_data[0].cut_index];

    /* local values storage */
    std::array<double, max_multioutput_output_size> arrival;
    std::array<float, max_multioutput_output_size> area_flow;
    std::array<float, max_multioutput_output_size> area;
    std::array<uint8_t, max_multioutput_output_size> phase;
    std::array<uint16_t, max_multioutput_output_size> pin_phase;
    std::array<double, max_multioutput_output_size> est_refs;
    std::array<uint32_t, max_multioutput_output_size> cut_index;
    bool mapped_multioutput = false;

    uint8_t iteration_phase = cut0->supergates[0] == nullptr ? 1 : 0;

    /* iterate for each possible match */
    for ( auto i = 0; i < cut0->supergates[iteration_phase]->size(); ++i )
    {
      /* store local validity and comparison info */
      bool valid = true;
      bool is_best = true;
      bool respects_required = true;
      double old_flow_sum = 0;

      /* iterate for each output of the multi-output gate */
      for ( auto j = 0; j < max_multioutput_output_size; ++j )
      {
        uint32_t node_index = tuple_data[j].node_index;
        cut_index[j] = tuple_data[j].cut_index;
        auto& node_data = node_match[node_index];
        auto const& cut = cuts[node_index][cut_index[j]];
        uint8_t phase_inverted = cut->supergates[0] == nullptr ? 1 : 0;
        supergate<NInputs> const& gate = ( *( cut->supergates[phase_inverted] ) )[i];

        /* protection on complicated duplicated nodes to remap to multioutput */
        if ( !node_data.same_match )
          return false;

        /* get the output phase */
        pin_phase[j] = gate.polarity;
        phase[j] = ( gate.polarity >> NInputs ) ^ phase_inverted;

        /* compute arrival */
        arrival[j] = 0.0;
        auto ctr = 0u;
        for ( auto l : cut )
        {
          /* multi-output matching is block-delay-only and stays outside the electrical model
           * (ADR-0047), so this leaf read is deliberately the scalar one: a load-indexed query
           * here would mix a load-aware leaf term into a load-blind candidate cost */
          double arrival_pin = node_match[l].arrival[( gate.polarity >> ctr ) & 1] + gate.tdelay[ctr];
          arrival[j] = std::max( arrival[j], arrival_pin );
          ++ctr;
        }

        /* check required time: same_match is true */
        if constexpr ( DO_AREA )
        {
          if ( arrival[j] > node_data.required[phase[j]] + epsilon )
          {
            valid = false;
            break;
          }
          if ( arrival[j] + lib_inv_delay > node_data.required[phase[j] ^ 1] + epsilon )
          {
            valid = false;
            break;
          }
        }

        /* check required time of the current solution */
        if ( node_data.arrival[phase[j]] > node_data.required[phase[j]] )
          respects_required = false;
        if ( node_data.same_match && node_data.arrival[phase[j] ^ 1] > node_data.required[phase[j] ^ 1] )
          respects_required = false;

        /* compute area flow */
        if ( j == 0 || !node_data.multioutput_match[0] )
        {
          uint8_t current_phase = node_data.best_gate[0] == nullptr ? 1 : 0;
          old_flow_sum += node_data.flows[current_phase];
        }
        uint8_t old_phase = node_data.phase[phase[j]];
        node_data.phase[phase[j]] = gate.polarity;
        area[j] = gate.area;
        area_flow[j] = gate.area + cut_leaves_flow( cut, n, phase[j] );
        node_data.phase[phase[j]] = old_phase;

        /* current version may lead to delay increase */
        est_refs[j] = node_data.est_refs[phase[j]];
      }

      /* not better than individual gates */
      if ( !valid )
        continue;

      if constexpr ( !DO_AREA )
      {
        if ( !is_best )
          continue;
      }

      /* combine evaluation for precise area flow estimantion */
      /* compute equation AF(n) = ( Area(G) + |roots| * SUM_{l in leaves} AF(l) ) / SUM_{p in roots} est_refs( p ) */
      float flow_sum_pos = 0, flow_sum_neg;
      float combined_est_refs = 0;
      for ( auto j = 0; j < max_multioutput_output_size; ++j )
      {
        flow_sum_pos += area_flow[j];
        combined_est_refs += est_refs[j];
      }
      flow_sum_neg = flow_sum_pos;
      flow_sum_pos /= combined_est_refs;

      /* not better than individual gates */
      if ( respects_required && ( flow_sum_pos > old_flow_sum + epsilon ) )
        continue;

      mapped_multioutput = true;
      flow_sum_neg = ( flow_sum_neg + lib_inv_area ) / combined_est_refs;

      /* commit multi-output gate */
      for ( uint32_t j = 0; j < max_multioutput_output_size; ++j )
      {
        uint32_t node_index = tuple_data[j].node_index;
        auto& node_data = node_match[node_index];
        auto const& cut = cuts[node_index][cut_index[j]];
        uint8_t phase_inverted = cut->supergates[0] == nullptr ? 1 : 0;
        supergate<NInputs> const& gate = ( *( cut->supergates[phase_inverted] ) )[i];

        uint8_t mapped_phase = phase[j];
        node_data.multioutput_match[mapped_phase] = true;

        node_data.best_gate[mapped_phase] = &gate;
        node_data.best_cut[mapped_phase] = cut_index[j];
        node_data.phase[mapped_phase] = pin_phase[j];
        node_data.arrival[mapped_phase] = arrival[j];
        node_data.area[mapped_phase] = area[j]; /* partial area contribution */
        node_data.flows[mapped_phase] = flow_sum_pos;

        assert( node_data.arrival[mapped_phase] < node_data.required[mapped_phase] + epsilon );

        /* select opposite phase */
        mapped_phase ^= 1;
        node_data.multioutput_match[mapped_phase] = true;
        node_data.best_gate[mapped_phase] = nullptr;
        node_data.best_cut[mapped_phase] = cut_index[j];
        node_data.phase[mapped_phase] = pin_phase[j];
        node_data.arrival[mapped_phase] = arrival[j] + lib_inv_delay;
        node_data.area[mapped_phase] = area[j]; /* partial area contribution */
        node_data.flows[mapped_phase] = flow_sum_neg;

        assert( node_data.arrival[mapped_phase] < node_data.required[mapped_phase] + epsilon );
      }
    }

    return mapped_multioutput;
  }

  template<bool SwitchActivity>
  bool match_multioutput_exact( node<Ntk> const& n, bool last_round )
  {
    /* extract outputs tuple */
    uint32_t index = ntk.node_to_index( n );
    multi_match_t const& tuple_data = multi_node_match[node_tuple_match[index].index][0];

    /* local values storage */
    std::array<float, max_multioutput_output_size> best_exact_area;

    for ( int j = max_multioutput_output_size - 1; j >= 0; --j )
    {
      /* protection on complicated duplicated nodes to remap to multioutput */
      if ( !node_match[tuple_data[j].node_index].same_match )
        return false;
    }

    /* if one of the outputs is not referenced, do not use multi-output gate */
    if ( last_round )
    {
      for ( uint32_t j = 0; j < max_multioutput_output_size; ++j )
      {
        uint32_t node_index = tuple_data[j].node_index;
        if ( !node_match[node_index].map_refs[0] && !node_match[node_index].map_refs[1] )
        {
          return false;
        }
      }
    }

    /* if "same match" and used in the cover dereference the leaves (reverse topo order) */
    for ( int j = max_multioutput_output_size - 1; j >= 0; --j )
    {
      uint32_t node_index = tuple_data[j].node_index;
      uint8_t selected_phase = node_match[node_index].best_gate[0] == nullptr ? 1 : 0;

      if ( node_match[node_index].map_refs[0] || node_match[node_index].map_refs[1] )
      {
        /* match is always single output here */
        auto const& cut = cuts[node_index][node_match[node_index].best_cut[0]];
        uint8_t use_phase = node_match[node_index].best_gate[0] != nullptr ? 0 : 1;
        best_exact_area[j] = cut_deref<SwitchActivity>( cut, ntk.index_to_node( node_index ), use_phase );

        /* mapping a non referenced phase */
        if ( node_match[node_index].map_refs[selected_phase] == 0 )
          best_exact_area[j] += lib_inv_area;
      }
    }

    /* perform mapping */
    bool mapped_multioutput = false;
    mapped_multioutput = match_multioutput_exact_core<SwitchActivity>( tuple_data, best_exact_area );

    /* if "same match" and used in the cover reference the leaves (topo order) */
    for ( auto j = 0; j < max_multioutput_output_size; ++j )
    {
      uint32_t node_index = tuple_data[j].node_index;

      if ( node_match[node_index].map_refs[0] || node_match[node_index].map_refs[1] )
      {
        uint8_t use_phase = node_match[node_index].best_gate[0] != nullptr ? 0 : 1;
        auto const& best_cut = cuts[node_index][node_match[node_index].best_cut[use_phase]];
        cut_ref<SwitchActivity>( best_cut, ntk.index_to_node( node_index ), use_phase );
      }
    }

    return mapped_multioutput;
  }

  template<bool SwitchActivity>
  inline bool match_multioutput_exact_core( multi_match_t const& tuple_data, std::array<float, max_multioutput_output_size>& best_exact_area )
  {
    /* get the cut representative */
    auto const& cut0 = cuts[tuple_data[0].node_index][tuple_data[0].cut_index];

    /* local values storage */
    std::array<double, max_multioutput_output_size> arrival;
    std::array<float, max_multioutput_output_size> area_exact;
    std::array<float, max_multioutput_output_size> area;
    std::array<uint8_t, max_multioutput_output_size> phase;
    std::array<uint16_t, max_multioutput_output_size> pin_phase;
    std::array<uint32_t, max_multioutput_output_size> cut_index;

    uint8_t iteration_phase = cut0->supergates[0] == nullptr ? 1 : 0;

    bool mapped_multioutput = false;

    /* iterate for each possible match */
    for ( auto i = 0; i < cut0->supergates[iteration_phase]->size(); ++i )
    {
      /* store local validity and comparison info */
      bool valid = true;
      bool is_best = true;
      bool respects_required = true;
      uint32_t it_counter = 0;

      /* iterate for each output of the multi-output gate (reverse topo order) */
      for ( int j = max_multioutput_output_size - 1; j >= 0; --j )
      {
        uint32_t node_index = tuple_data[j].node_index;
        cut_index[j] = tuple_data[j].cut_index;
        auto& node_data = node_match[node_index];
        auto const& cut = cuts[node_index][cut_index[j]];
        uint8_t phase_inverted = cut->supergates[0] == nullptr ? 1 : 0;
        supergate<NInputs> const& gate = ( *( cut->supergates[phase_inverted] ) )[i];
        ++it_counter;

        /* get the output phase and area */
        pin_phase[j] = gate.polarity;
        phase[j] = ( gate.polarity >> NInputs ) ^ phase_inverted;
        area[j] = gate.area;

        /* compute arrival */
        arrival[j] = 0.0;
        auto ctr = 0u;
        for ( auto l : cut )
        {
          /* multi-output matching is block-delay-only and stays outside the electrical model
           * (ADR-0047), so this leaf read is deliberately the scalar one: a load-indexed query
           * here would mix a load-aware leaf term into a load-blind candidate cost */
          double arrival_pin = node_match[l].arrival[( gate.polarity >> ctr ) & 1] + gate.tdelay[ctr];
          arrival[j] = std::max( arrival[j], arrival_pin );
          ++ctr;
        }

        /* check required time */
        if ( arrival[j] > node_data.required[phase[j]] + epsilon )
        {
          valid = false;
          break;
        }
        if ( arrival[j] + lib_inv_delay > node_data.required[phase[j] ^ 1] + epsilon )
        {
          valid = false;
          break;
        }

        /* check required time of current solution */
        if ( node_data.arrival[phase[j]] > node_data.required[phase[j]] )
          respects_required = false;
        if ( node_data.arrival[phase[j] ^ 1] > node_data.required[phase[j] ^ 1] )
          respects_required = false;

        /* compute exact area for match: needed only for the first node (leaves are shared) */
        if ( it_counter == 1 )
        {
          auto old_phase = node_data.phase[phase[j]];
          auto old_area = node_data.area[phase[j]];
          node_data.phase[phase[j]] = pin_phase[j];
          node_data.area[phase[j]] = area[j];
          area_exact[j] = cut_measure_mffc<SwitchActivity>( cut, ntk.index_to_node( node_index ), phase[j] );
          node_data.phase[phase[j]] = old_phase;
          node_data.area[phase[j]] = old_area;
        }
        else
        {
          area_exact[j] = area[j];
        }

        /* Add output inverter cost if mapping a non referenced phase */
        if ( node_data.map_refs[phase[j]] == 0 && node_data.map_refs[phase[j] ^ 1] > 0 )
        {
          area_exact[j] += lib_inv_area;
        }
      }

      /* check quality: TODO add output inverter in the cost if necessary */
      float best_exact_area_total = 0;
      float area_exact_total = 0;
      for ( auto j = 0; j < max_multioutput_output_size; ++j )
      {
        best_exact_area_total += best_exact_area[j];
        area_exact_total += area_exact[j];
      }

      /* not better than individual gates */
      if ( !valid || ( area_exact_total > best_exact_area_total - epsilon && respects_required ) )
      {
        continue;
      }

      mapped_multioutput = true;

      /* commit multi-output gate (topo order) */
      for ( uint32_t j = 0; j < max_multioutput_output_size; ++j )
      {
        uint32_t node_index = tuple_data[j].node_index;
        auto& node_data = node_match[node_index];
        auto const& cut = cuts[node_index][cut_index[j]];
        uint8_t phase_inverted = cut->supergates[0] == nullptr ? 1 : 0;
        supergate<NInputs> const& gate = ( *( cut->supergates[phase_inverted] ) )[i];

        uint8_t mapped_phase = phase[j];
        best_exact_area[j] = area_exact[j];

        if ( node_data.map_refs[phase[j]] == 0 && node_data.map_refs[phase[j] ^ 1] > 0 )
        {
          best_exact_area[j] += lib_inv_area;
        }

        /* write data */
        node_data.multioutput_match[mapped_phase] = true;
        node_data.best_gate[mapped_phase] = &gate;
        node_data.best_cut[mapped_phase] = cut_index[j];
        node_data.phase[mapped_phase] = pin_phase[j];
        node_data.arrival[mapped_phase] = arrival[j];
        node_data.area[mapped_phase] = area[j]; /* partial area contribution */

        node_data.flows[mapped_phase] = area_exact[j]; /* partial exact area contribution */
        /* select opposite phase */
        mapped_phase ^= 1;
        node_data.multioutput_match[mapped_phase] = true;
        node_data.best_gate[mapped_phase] = nullptr;
        node_data.best_cut[mapped_phase] = cut_index[j];
        node_data.phase[mapped_phase] = pin_phase[j];
        node_data.arrival[mapped_phase] = arrival[j] + lib_inv_delay;
        node_data.area[mapped_phase] = area[j]; /* partial area contribution */
        node_data.flows[mapped_phase] = area_exact[j];

        assert( node_data.arrival[mapped_phase] < node_data.required[mapped_phase] + epsilon );
      }
    }

    return mapped_multioutput;
  }

  template<bool DO_AREA>
  void multi_node_update( node<Ntk> const& n )
  {
    uint32_t check_index = ntk.node_to_index( n );
    multi_match_t const& tuple_data = multi_node_match[node_tuple_match[ntk.node_to_index( n )].index][0];
    uint64_t signature = 0;

    /* check if a node is in TFI: there is a path of length > 1 */
    bool in_tfi = false;
    node<Ntk> min_node = n;
    for ( auto j = 0; j < max_multioutput_output_size - 1; ++j )
    {
      if ( tuple_data[j].in_tfi )
      {
        min_node = ntk.index_to_node( tuple_data[j].node_index );
        in_tfi = true;
        signature |= UINT64_C( 1 ) << ( tuple_data[j].node_index & 0x3f );
      }
    }

    if ( !in_tfi )
      return;

    /* recompute data in between: should I mark the leaves? (not necessary under some assumptions) */
    ntk.incr_trav_id();
    ntk.foreach_fanin( n, [&]( auto const& f ) {
      /* TODO: this recursion works as it is for a maximum multioutput value of 2 */
      multi_node_update_rec<DO_AREA>( ntk.get_node( f ), min_node + 1, signature );
    } );
  }

  template<bool DO_AREA>
  void multi_node_update_rec( node<Ntk> const& n, uint32_t min_index, uint64_t& signature )
  {
    uint32_t index = ntk.node_to_index( n );

    if ( index < min_index )
      return;
    if ( ntk.visited( n ) == ntk.trav_id() )
      return;

    ntk.set_visited( n, ntk.trav_id() );
    ntk.foreach_fanin( n, [&]( auto const& f ) {
      multi_node_update_rec<DO_AREA>( ntk.get_node( f ), min_index, signature );
    } );

    /* update the node if uses an updated leaf */
    auto& node_data = node_match[index];
    bool leaf_used = multi_node_update_cut_check( index, signature, 0 );

    if ( !node_data.same_match )
      leaf_used |= multi_node_update_cut_check( index, signature, 1 );

    if ( !leaf_used )
      return;

    signature |= UINT64_C( 1 ) << ( index & 0x3f );

    /* avoid cycles by recomputing arrival times for multi-output gates or decomposing them */
    if ( node_data.same_match && node_data.multioutput_match[0] )
    {
      propagate_arrival_node( n );
      /* check required time */
      if ( node_data.arrival[0] < node_data.required[0] + epsilon && node_data.arrival[1] < node_data.required[1] + epsilon )
        return;
    }

    /* match positive phase */
    match_phase<DO_AREA>( n, 0u );

    /* match negative phase */
    match_phase<DO_AREA>( n, 1u );

    /* try to drop one phase */
    match_drop_phase<DO_AREA, false>( n );

    assert( node_data.arrival[0] < node_data.required[0] + epsilon );
    assert( node_data.arrival[1] < node_data.required[1] + epsilon );
  }

  template<bool SwitchActivity>
  void multi_node_update_exact( node<Ntk> const& n )
  {
    uint32_t check_index = ntk.node_to_index( n );
    multi_match_t const& tuple_data = multi_node_match[node_tuple_match[ntk.node_to_index( n )].index][0];
    uint64_t signature = 0;

    /* check if a node is in TFI: there is a path of length > 1 */
    bool in_tfi = false;
    node<Ntk> min_node = n;
    for ( auto j = 0; j < max_multioutput_output_size - 1; ++j )
    {
      if ( tuple_data[j].in_tfi )
      {
        min_node = ntk.index_to_node( tuple_data[j].node_index );
        in_tfi = true;
        signature |= UINT64_C( 1 ) << ( tuple_data[j].node_index & 0x3f );
      }
    }

    if ( !in_tfi )
      return;

    /* recompute data in between: should I mark the leaves? (not necessary under some assumptions) */
    ntk.incr_trav_id();
    ntk.foreach_fanin( n, [&]( auto const& f ) {
      /* TODO: this recursion works as it is for a maximum multioutput value of 2 */
      multi_node_update_exact_rec<SwitchActivity>( ntk.get_node( f ), min_node + 1, signature );
    } );
  }

  template<bool SwitchActivity>
  void multi_node_update_exact_rec( node<Ntk> const& n, uint32_t min_index, uint64_t& signature )
  {
    uint32_t index = ntk.node_to_index( n );

    if ( index < min_index )
      return;
    if ( ntk.visited( n ) == ntk.trav_id() )
      return;

    ntk.set_visited( n, ntk.trav_id() );
    ntk.foreach_fanin( n, [&]( auto const& f ) {
      multi_node_update_exact_rec<SwitchActivity>( ntk.get_node( f ), min_index, signature );
    } );

    /* update the node if uses an updated leaf */
    auto& node_data = node_match[index];
    bool leaf_used = multi_node_update_cut_check( index, signature, 0 );

    if ( !node_data.same_match )
      leaf_used |= multi_node_update_cut_check( index, signature, 1 );

    if ( !leaf_used )
      return;

    signature |= UINT64_C( 1 ) << ( index & 0x3f );

    assert( !node_data.multioutput_match[0] );
    assert( !node_data.multioutput_match[1] );

    if ( node_data.same_match && ( node_data.map_refs[0] || node_data.map_refs[1] ) )
    {
      uint8_t use_phase = node_data.best_gate[0] != nullptr ? 0 : 1;
      auto const& best_cut = cuts[index][node_data.best_cut[use_phase]];
      cut_deref<SwitchActivity>( best_cut, n, use_phase );
    }

    /* match positive phase */
    match_phase_exact<SwitchActivity>( n, 0u );

    /* match negative phase */
    match_phase_exact<SwitchActivity>( n, 1u );

    /* try to drop one phase */
    match_drop_phase<true, true>( n );

    assert( node_data.arrival[0] < std::numeric_limits<float>::max() );
    assert( node_data.arrival[1] < std::numeric_limits<float>::max() );
  }

  inline void match_multioutput_propagate_required( node<Ntk> const& n )
  {
    /* extract outputs tuple */
    uint32_t index = ntk.node_to_index( n );
    multi_match_t const& tuple_data = multi_node_match[node_tuple_match[index].index][0];

    for ( int j = max_multioutput_output_size - 1; j >= 0; --j )
    {
      const auto node_index = tuple_data[j].node_index;
      match_propagate_required( node_index );
    }
  }

  void match_multi_add_cuts( node<Ntk> const& n )
  {
    /* assume a single cut (current version) */
    uint32_t index = ntk.node_to_index( n );
    multi_match_t& matches = multi_node_match[node_tuple_match[index].index][0];

    /* find the corresponding cut */
    uint32_t cut_p = 0;
    while ( matches[cut_p].node_index != index )
      ++cut_p;

    assert( cut_p < matches.size() );
    uint32_t cut_index = matches[cut_p].cut_index;
    auto& cut = multi_cut_set[cut_index][cut_p];
    auto single_cut = multi_cut_set[cut_index][cut_p];
    auto& rcuts = cuts[index];

    /* not enough space in the data structure: abort */
    if ( rcuts.size() == max_cut_num )
    {
      match_multi_add_cuts_remove_entry( matches );
      return;
    }

    /* insert single cut variation if unique (for delay preservation) */
    if ( !rcuts.is_contained( single_cut ) )
    {
      single_cut->pattern_index = 0;
      compute_cut_data( single_cut, ntk.index_to_node( index ) );
      rcuts.append_cut( single_cut );

      /* not enough space in the data structure: abort */
      if ( rcuts.size() == max_cut_num )
      {
        rcuts.limit( rcuts.size() - 1 );
        match_multi_add_cuts_remove_entry( matches );
        return;
      }
    }

    /* add multi-output cut */
    uint32_t num_cuts_pre = rcuts.size();
    cut->ignore = true;
    rcuts.append_cut( cut );

    uint32_t num_cuts_after = rcuts.size();
    assert( num_cuts_after == num_cuts_pre + 1 );

    rcuts.limit( num_cuts_pre );

    /* update tuple data */
    matches[cut_p].cut_index = num_cuts_pre;
  }

  inline void match_multi_add_cuts_remove_entry( multi_match_t const& matches )
  {
    /* reset matches */
    for ( multi_match_data const& entry : matches )
    {
      node_tuple_match[entry.node_index].data = 0;
    }
  }

  inline bool multi_node_update_cut_check( uint32_t index, uint64_t signature, uint8_t phase )
  {
    auto const& cut = cuts[index][node_match[index].best_cut[phase]];

    if ( ( signature & cut.signature() ) > 0 )
      return true;

    return false;
  }
#pragma endregion

#pragma region Mapping utils
  inline double cut_leaves_flow( cut_t const& cut, node<Ntk> const& n, uint8_t phase )
  {
    double flow{ 0.0f };
    auto const& node_data = node_match[ntk.node_to_index( n )];

    uint8_t ctr = 0u;
    for ( auto leaf : cut )
    {
      uint8_t leaf_phase = ( node_data.phase[phase] >> ctr++ ) & 1;
      flow += node_match[leaf].flows[leaf_phase];
    }

    return flow;
  }

  template<bool SwitchActivity>
  float cut_ref( cut_t const& cut, node<Ntk> const& n, uint8_t phase )
  {
    if constexpr ( !SwitchActivity )
      if ( frontier_frames ) return frontier_cut_walk( cut, n, phase, frontier_walk_mode::reference );
    auto const& node_data = node_match[ntk.node_to_index( n )];
    float count;

    if constexpr ( SwitchActivity )
      count = switch_activity[ntk.node_to_index( n )];
    else
      count = node_data.area[phase];

    /* don't touch box */
    if constexpr ( has_is_dont_touch_v<Ntk> )
    {
      if ( ntk.is_dont_touch( n ) )
      {
        return count;
      }
    }

    uint8_t ctr = 0;
    for ( auto leaf : cut )
    {
      /* compute leaf phase using the current gate */
      uint8_t leaf_phase = ( node_data.phase[phase] >> ctr++ ) & 1;

      if ( ntk.is_constant( ntk.index_to_node( leaf ) ) )
      {
        continue;
      }
      else if ( ntk.is_pi( ntk.index_to_node( leaf ) ) )
      {
        /* reference PIs, add inverter cost for negative phase */
        if ( leaf_phase == 1u )
        {
          if ( node_match[leaf].map_refs[1]++ == 0u )
          {
            if constexpr ( SwitchActivity )
              count += switch_activity[leaf];
            else
              count += lib_inv_area;
          }
        }
        else
        {
          ++node_match[leaf].map_refs[0];
        }
        continue;
      }

      if ( node_match[leaf].same_match )
      {
        /* Recursive referencing if leaf was not referenced */
        if ( !node_match[leaf].map_refs[0] && !node_match[leaf].map_refs[1] )
        {
          auto const& best_cut = cuts[leaf][node_match[leaf].best_cut[leaf_phase]];
          count += cut_ref<SwitchActivity>( best_cut, ntk.index_to_node( leaf ), leaf_phase );
        }

        /* Add inverter area if not present yet and leaf node is implemented in the opposite phase */
        if ( node_match[leaf].map_refs[leaf_phase]++ == 0u && node_match[leaf].best_gate[leaf_phase] == nullptr )
        {
          if constexpr ( SwitchActivity )
            count += switch_activity[leaf];
          else
            count += lib_inv_area;
        }
      }
      else
      {
        if ( node_match[leaf].map_refs[leaf_phase]++ == 0u )
        {
          auto const& best_cut = cuts[leaf][node_match[leaf].best_cut[leaf_phase]];
          count += cut_ref<SwitchActivity>( best_cut, ntk.index_to_node( leaf ), leaf_phase );
        }
      }
    }
    return count;
  }

  template<bool SwitchActivity>
  float cut_deref( cut_t const& cut, node<Ntk> const& n, uint8_t phase )
  {
    if constexpr ( !SwitchActivity )
      if ( frontier_frames )
        return frontier_cut_walk( cut, n, phase, frontier_walk_mode::dereference );
    auto const& node_data = node_match[ntk.node_to_index( n )];
    float count;

    if constexpr ( SwitchActivity )
      count = switch_activity[ntk.node_to_index( n )];
    else
      count = node_data.area[phase];

    /* don't touch box */
    if constexpr ( has_is_dont_touch_v<Ntk> )
    {
      if ( ntk.is_dont_touch( n ) )
      {
        return count;
      }
    }

    uint8_t ctr = 0;
    for ( auto leaf : cut )
    {
      /* compute leaf phase using the current gate */
      uint8_t leaf_phase = ( node_data.phase[phase] >> ctr++ ) & 1;

      if ( ntk.is_constant( ntk.index_to_node( leaf ) ) )
      {
        continue;
      }
      else if ( ntk.is_pi( ntk.index_to_node( leaf ) ) )
      {
        /* dereference PIs, add inverter cost for negative phase */
        if ( leaf_phase == 1u )
        {
          if ( --node_match[leaf].map_refs[1] == 0u )
          {
            if constexpr ( SwitchActivity )
              count += switch_activity[leaf];
            else
              count += lib_inv_area;
          }
        }
        else
        {
          --node_match[leaf].map_refs[0];
        }
        continue;
      }

      if ( node_match[leaf].same_match )
      {
        /* Add inverter area if it is used only by the current gate and leaf node is implemented in the opposite phase */
        if ( --node_match[leaf].map_refs[leaf_phase] == 0u && node_match[leaf].best_gate[leaf_phase] == nullptr )
        {
          if constexpr ( SwitchActivity )
            count += switch_activity[leaf];
          else
            count += lib_inv_area;
        }
        /* Recursive dereferencing */
        if ( !node_match[leaf].map_refs[0] && !node_match[leaf].map_refs[1] )
        {
          auto const& best_cut = cuts[leaf][node_match[leaf].best_cut[leaf_phase]];
          count += cut_deref<SwitchActivity>( best_cut, ntk.index_to_node( leaf ), leaf_phase );
        }
      }
      else
      {
        if ( --node_match[leaf].map_refs[leaf_phase] == 0u )
        {
          auto const& best_cut = cuts[leaf][node_match[leaf].best_cut[leaf_phase]];
          count += cut_deref<SwitchActivity>( best_cut, ntk.index_to_node( leaf ), leaf_phase );
        }
      }
    }
    return count;
  }

  template<bool SwitchActivity>
  float cut_measure_mffc( cut_t const& cut, node<Ntk> const& n, uint8_t phase )
  {
    if constexpr ( !SwitchActivity )
    {
      if ( frontier_frames )
      {
        frontier_visit_count = 0;
        float const count = frontier_cut_walk( cut, n, phase, frontier_walk_mode::visit );
        for ( std::size_t i = 0; i < frontier_visit_count; ++i )
          --node_match[frontier_visits[i] >> 1].map_refs[frontier_visits[i] & 1];
        return count;
      }
    }
    tmp_visited.clear();

    float count = cut_ref_visit<SwitchActivity>( cut, n, phase );

    /* dereference visited */
    for ( auto s : tmp_visited )
    {
      uint32_t leaf = s >> 1;
      --node_match[leaf].map_refs[s & 1];
    }

    return count;
  }

  template<bool SwitchActivity>
  float cut_ref_visit( cut_t const& cut, node<Ntk> const& n, uint8_t phase )
  {
    if constexpr ( !SwitchActivity )
      if ( frontier_frames )
        return frontier_cut_walk( cut, n, phase, frontier_walk_mode::visit );
    auto const& node_data = node_match[ntk.node_to_index( n )];
    float count;

    if constexpr ( SwitchActivity )
      count = switch_activity[ntk.node_to_index( n )];
    else
      count = node_data.area[phase];

    /* don't touch box */
    if constexpr ( has_is_dont_touch_v<Ntk> )
    {
      if ( ntk.is_dont_touch( n ) )
      {
        return count;
      }
    }

    uint8_t ctr = 0;
    for ( auto leaf : cut )
    {
      /* compute leaf phase using the current gate */
      uint8_t leaf_phase = ( node_data.phase[phase] >> ctr++ ) & 1;

      if ( ntk.is_constant( ntk.index_to_node( leaf ) ) )
      {
        continue;
      }

      /* add to visited */
      tmp_visited.push_back( ( static_cast<uint64_t>( leaf ) << 1 ) | leaf_phase );

      if ( ntk.is_pi( ntk.index_to_node( leaf ) ) )
      {
        /* reference PIs, add inverter cost for negative phase */
        if ( leaf_phase == 1u )
        {
          if ( node_match[leaf].map_refs[1]++ == 0u )
          {
            if constexpr ( SwitchActivity )
              count += switch_activity[leaf];
            else
              count += lib_inv_area;
          }
        }
        else
        {
          ++node_match[leaf].map_refs[0];
        }
        continue;
      }

      if ( node_match[leaf].same_match )
      {
        /* Recursive referencing if leaf was not referenced */
        if ( !node_match[leaf].map_refs[0] && !node_match[leaf].map_refs[1] )
        {
          auto const& best_cut = cuts[leaf][node_match[leaf].best_cut[leaf_phase]];
          count += cut_ref_visit<SwitchActivity>( best_cut, ntk.index_to_node( leaf ), leaf_phase );
        }

        /* Add inverter area if not present yet and leaf node is implemented in the opposite phase */
        if ( node_match[leaf].map_refs[leaf_phase]++ == 0u && node_match[leaf].best_gate[leaf_phase] == nullptr )
        {
          if constexpr ( SwitchActivity )
            count += switch_activity[leaf];
          else
            count += lib_inv_area;
        }
      }
      else
      {
        if ( node_match[leaf].map_refs[leaf_phase]++ == 0u )
        {
          auto const& best_cut = cuts[leaf][node_match[leaf].best_cut[leaf_phase]];
          count += cut_ref_visit<SwitchActivity>( best_cut, ntk.index_to_node( leaf ), leaf_phase );
        }
      }
    }
    return count;
  }
#pragma endregion

#pragma region Initialize and dump the mapped network
  void insert_buffers()
  {
    if ( lib_buf_id != UINT32_MAX )
    {
      double area_old = area;
      bool buffers = false;

      ntk.foreach_po( [&]( auto const& f ) {
        auto const& n = ntk.get_node( f );
        if ( !ntk.is_constant( n ) && ntk.is_pi( n ) && !ntk.is_complemented( f ) )
        {
          area += lib_buf_area;
          delay = std::max( delay, node_match[ntk.node_to_index( n )].arrival[0] + lib_inv_delay );
          buffers = true;
        }
      } );

      /* round stats */
      if ( ps.verbose && buffers )
      {
        std::stringstream stats{};
        float area_gain = 0.0f;

        area_gain = float( ( area_old - area ) / area_old * 100 );

        stats << fmt::format( "[i] Buffering: Delay = {:>12.2f}  Area = {:>12.2f}  Gain = {:>5.2f} %  Inverters = {:>5}  Time = {:>5.2f}\n", delay, area, area_gain, inv, to_seconds( clock::now() - time_begin ) );
        st.round_stats.push_back( stats.str() );
      }
    }
  }

  std::pair<binding_view<klut_network>, klut_map> initialize_map_network()
  {
    binding_view<klut_network> dest( library.get_gates() );
    klut_map old2new;

    old2new[ntk.node_to_index( ntk.get_node( ntk.get_constant( false ) ) )][0] = dest.get_constant( false );
    old2new[ntk.node_to_index( ntk.get_node( ntk.get_constant( false ) ) )][1] = dest.get_constant( true );

    ntk.foreach_pi( [&]( auto const& n ) {
      old2new[ntk.node_to_index( n )][0] = dest.create_pi();
    } );
    return { dest, old2new };
  }

  std::pair<cell_view<block_network>, block_map> initialize_block_network()
  {
    cell_view<block_network> dest( library.get_cells() );
    block_map old2new;

    old2new[ntk.node_to_index( ntk.get_node( ntk.get_constant( false ) ) )][0] = dest.get_constant( false );
    old2new[ntk.node_to_index( ntk.get_node( ntk.get_constant( false ) ) )][1] = dest.get_constant( true );

    ntk.foreach_pi( [&]( auto const& n ) {
      old2new[ntk.node_to_index( n )][0] = dest.create_pi();
    } );
    return { dest, old2new };
  }

  void init_topo_order()
  {
    topo_order.reserve( ntk.size() );

    if ( multi_node_match.size() > 0 )
    {
      multi_init_topo_order();
      return;
    }

    topo_view<Ntk>( ntk ).foreach_node( [this]( auto n ) {
      topo_order.push_back( n );
    } );
  }

  bool init_arrivals()
  {
    if ( ps.required_times.size() && ps.required_times.size() != ntk.num_pos() )
    {
      std::cerr << "[e] MAP ERROR: required time vector does not match the output size of the network" << std::endl;
      st.mapping_error = true;
      return false;
    }

    /* seed the electrical load estimate before any arrival uses it (no-op with the model off) */
    init_node_loads();

    if ( ps.arrival_times.empty() )
    {
      ntk.foreach_pi( [&]( auto const& n ) {
        uint32_t const index = ntk.node_to_index( n );
        auto& node_data = node_match[index];
        node_data.arrival[0] = node_data.best_alternative[0].arrival = 0;
        node_data.arrival[1] = node_data.best_alternative[1].arrival = inv_delay_at( index, 1 );
      } );
      init_node_curves();
      return true;
    }

    if ( ps.arrival_times.size() != ntk.num_pis() )
    {
      std::cerr << "[e] MAP ERROR: arrival time vector does not match the input size of the network" << std::endl;
      st.mapping_error = true;
      return false;
    }

    ntk.foreach_pi( [&]( auto const& n, uint32_t i ) {
      uint32_t const index = ntk.node_to_index( n );
      auto& node_data = node_match[index];
      node_data.arrival[0] = node_data.best_alternative[0].arrival = ps.arrival_times[i];
      node_data.arrival[1] = node_data.best_alternative[1].arrival = ps.arrival_times[i] + inv_delay_at( index, 1 );
    } );

    init_node_curves();
    return true;
  }

  /* Materialize the covering stage's chosen output buffer on (index, phase): a real 1-input
   * identity node bound to the SELECTED drive, rebinding old2new so every later consumer (the
   * walk is topological) — and the same_match polarity inverter, created after this — reads the
   * buffered signal, exactly matching the DP's arrival model. The recipe is the PO-buffer one
   * (create_node over tt 0x2 + add_binding); klut create_buf is identity and creates no node. */
  void materialize_cover_buffer( binding_view<klut_network>& res, klut_map& old2new, uint32_t index, uint8_t phase )
  {
    if ( !buffer_stage || !node_buffers[index][phase].active )
      return;
    static uint64_t _buf_tt = 0x2;
    kitty::dynamic_truth_table tt_buf( 1 );
    kitty::create_from_words( tt_buf, &_buf_tt, &_buf_tt + 1 );
    const auto buf = res.create_node( { old2new[index][phase] }, tt_buf );
    res.add_binding( res.get_node( buf ), node_buffers[index][phase].drive.id );
    old2new[index][phase] = buf;
  }

  void finalize_cover( binding_view<klut_network>& res, klut_map& old2new )
  {
    uint32_t multioutput_count = 0;

    for ( auto const& n : topo_order )
    {
      auto index = ntk.node_to_index( n );
      auto const& node_data = node_match[index];

      /* add inverter at PI if needed */
      if ( ntk.is_constant( n ) )
      {
        if ( node_data.best_gate[0] == nullptr && node_data.best_gate[1] == nullptr )
          continue;
      }
      else if ( ntk.is_pi( n ) )
      {
        if ( node_data.map_refs[1] > 0 )
        {
          old2new[index][1] = res.create_not( old2new[n][0] );
          res.add_binding( res.get_node( old2new[index][1] ),
                           select_polarity_inverter( index, 1 ) );
        }
        continue;
      }

      /* continue if cut is not in the cover */
      if ( !node_data.map_refs[0] && !node_data.map_refs[1] )
        continue;

      /* don't touch box */
      if constexpr ( has_is_dont_touch_v<Ntk> )
      {
        if ( ntk.is_dont_touch( n ) )
        {
          clone_box( res, old2new, index );
          continue;
        }
      }

      unsigned phase = ( node_data.best_gate[0] != nullptr ) ? 0 : 1;

      /* add used cut */
      if ( node_data.same_match || node_data.map_refs[phase] > 0 )
      {
        create_lut_for_gate( res, old2new, index, phase );
        materialize_cover_buffer( res, old2new, index, static_cast<uint8_t>( phase ) );

        /* add inverted version if used */
        if ( node_data.same_match && node_data.map_refs[phase ^ 1] > 0 )
        {
          old2new[index][phase ^ 1] = res.create_not( old2new[index][phase] );
          res.add_binding( res.get_node( old2new[index][phase ^ 1] ),
                           select_polarity_inverter( index, phase ^ 1 ) );
        }

        /* count multioutput gates */
        if ( ps.map_multioutput && node_tuple_match[index].lowest_index && node_data.multioutput_match[phase] )
        {
          ++multioutput_count;
        }
      }

      phase = phase ^ 1;
      /* add the optional other match if used */
      if ( !node_data.same_match && node_data.map_refs[phase] > 0 )
      {
        create_lut_for_gate( res, old2new, index, phase );
        materialize_cover_buffer( res, old2new, index, static_cast<uint8_t>( phase ) );

        /* count multioutput gates */
        if ( ps.map_multioutput && node_tuple_match[index].lowest_index && node_data.multioutput_match[phase] )
        {
          ++multioutput_count;
        }
      }

      st.multioutput_gates = multioutput_count;
    }

    /* create POs */
    ntk.foreach_po( [&]( auto const& f ) {
      if ( ntk.is_complemented( f ) )
      {
        res.create_po( old2new[ntk.node_to_index( ntk.get_node( f ) )][1] );
      }
      else if ( !ntk.is_constant( ntk.get_node( f ) ) && ntk.is_pi( ntk.get_node( f ) ) && lib_buf_id != UINT32_MAX )
      {
        /* create buffers for POs */
        static uint64_t _buf = 0x2;
        kitty::dynamic_truth_table tt_buf( 1 );
        kitty::create_from_words( tt_buf, &_buf, &_buf + 1 );
        const auto buf = res.create_node( { old2new[ntk.node_to_index( ntk.get_node( f ) )][0] }, tt_buf );
        res.create_po( buf );
        res.add_binding( res.get_node( buf ), lib_buf_id );
      }
      else
      {
        res.create_po( old2new[ntk.node_to_index( ntk.get_node( f ) )][0] );
      }
    } );

    /* write final results */
    st.area = area;
    st.delay = delay;
    report_drive_legality();
    if ( ps.eswp_rounds )
      st.power = compute_switching_power();
  }

  void finalize_cover_block( cell_view<block_network>& res, block_map& old2new )
  {
    uint32_t multioutput_count = 0;

    /* get standard cells */
    std::vector<standard_cell> const& lib = res.get_library();

    /* get translation ID from GENLIB to STD_CELL */
    std::vector<uint32_t> genlib_to_cell( library.get_gates().size() );
    for ( standard_cell const& cell : lib )
    {
      for ( gate const& g : cell.gates )
      {
        genlib_to_cell[g.id] = cell.id;
      }
    }

    for ( auto const& n : topo_order )
    {
      auto index = ntk.node_to_index( n );
      auto const& node_data = node_match[index];

      /* add inverter at PI if needed */
      if ( ntk.is_constant( n ) )
      {
        if ( node_data.best_gate[0] == nullptr && node_data.best_gate[1] == nullptr )
          continue;
      }
      else if ( ntk.is_pi( n ) )
      {
        if ( node_data.map_refs[1] > 0 )
        {
          old2new[index][1] = res.create_not( old2new[n][0] );
          res.add_cell( res.get_node( old2new[index][1] ), genlib_to_cell[lib_inv_id] );
        }
        continue;
      }

      /* continue if cut is not in the cover */
      if ( !node_data.map_refs[0] && !node_data.map_refs[1] )
        continue;

      /* don't touch box */
      if constexpr ( has_is_dont_touch_v<Ntk> )
      {
        if ( ntk.is_dont_touch( n ) )
        {
          clone_box2( res, old2new, index, genlib_to_cell );
          continue;
        }
      }

      unsigned phase = ( node_data.best_gate[0] != nullptr ) ? 0 : 1;

      /* add used cut */
      if ( node_data.same_match || node_data.map_refs[phase] > 0 )
      {
        /* create multioutput gates */
        if ( ps.map_multioutput && node_data.multioutput_match[phase] )
        {
          assert( node_data.same_match == true );

          if ( node_tuple_match[index].has_info && node_tuple_match[index].lowest_index )
          {
            ++multioutput_count;
            create_block_for_gate( res, old2new, index, phase, genlib_to_cell );
          }
          continue;
        }

        create_lut_for_gate2( res, old2new, index, phase, genlib_to_cell );

        /* add inverted version if used */
        if ( node_data.same_match && node_data.map_refs[phase ^ 1] > 0 )
        {
          old2new[index][phase ^ 1] = res.create_not( old2new[index][phase] );
          res.add_cell( res.get_node( old2new[index][phase ^ 1] ), genlib_to_cell[lib_inv_id] );
        }
      }

      phase = phase ^ 1;
      /* add the optional other match if used */
      if ( !node_data.same_match && node_data.map_refs[phase] > 0 )
      {
        assert( !ps.map_multioutput || !node_data.multioutput_match[phase] );
        create_lut_for_gate2( res, old2new, index, phase, genlib_to_cell );
      }
    }

    /* create POs */
    ntk.foreach_po( [&]( auto const& f ) {
      if ( ntk.is_complemented( f ) )
      {
        res.create_po( old2new[ntk.node_to_index( ntk.get_node( f ) )][1] );
      }
      else if ( !ntk.is_constant( ntk.get_node( f ) ) && ntk.is_pi( ntk.get_node( f ) ) && lib_buf_id != UINT32_MAX )
      {
        /* create buffers for POs */
        static uint64_t _buf = 0x2;
        kitty::dynamic_truth_table tt_buf( 1 );
        kitty::create_from_words( tt_buf, &_buf, &_buf + 1 );
        const auto buf = res.create_node( { old2new[ntk.node_to_index( ntk.get_node( f ) )][0] }, tt_buf );
        res.create_po( buf );
        res.add_cell( res.get_node( buf ), genlib_to_cell[lib_buf_id] );
      }
      else
      {
        res.create_po( old2new[ntk.node_to_index( ntk.get_node( f ) )][0] );
      }
    } );

    /* write final results */
    st.area = area;
    st.delay = delay;
    st.multioutput_gates = multioutput_count;
    report_drive_legality();
    if ( ps.eswp_rounds )
      st.power = compute_switching_power();
  }

  void create_lut_for_gate( binding_view<klut_network>& res, klut_map& old2new, uint32_t index, unsigned phase )
  {
    auto const& node_data = node_match[index];
    auto const& best_cut = cuts[index][node_data.best_cut[phase]];
    auto const& gate = node_data.best_gate[phase]->root;

    /* permutate and negate to obtain the matched gate truth table */
    std::vector<signal<klut_network>> children( gate->num_vars );

    auto ctr = 0u;
    for ( auto l : best_cut )
    {
      if ( ctr >= gate->num_vars )
        break;
      children[node_data.best_gate[phase]->permutation[ctr]] = old2new[l][( node_data.phase[phase] >> ctr ) & 1];
      ++ctr;
    }

    if ( !gate->is_super )
    {
      /* create the node */
      auto f = res.create_node( children, gate->function );
      res.add_binding( res.get_node( f ), gate->root->id );

      /* add the node in the data structure */
      old2new[index][phase] = f;
    }
    else
    {
      /* supergate, create sub-gates */
      auto f = create_lut_for_gate_rec( res, *gate, children );

      /* add the node in the data structure */
      old2new[index][phase] = f;
    }
  }

  signal<klut_network> create_lut_for_gate_rec( binding_view<klut_network>& res, composed_gate<NInputs> const& gate, std::vector<signal<klut_network>> const& children )
  {
    std::vector<signal<klut_network>> children_local( gate.fanin.size() );

    auto i = 0u;
    for ( auto const fanin : gate.fanin )
    {
      if ( fanin->root == nullptr )
      {
        /* terminal condition */
        children_local[i] = children[fanin->id];
      }
      else
      {
        children_local[i] = create_lut_for_gate_rec( res, *fanin, children );
      }
      ++i;
    }

    auto f = res.create_node( children_local, gate.root->function );
    res.add_binding( res.get_node( f ), gate.root->id );
    return f;
  }

  void create_lut_for_gate2( cell_view<block_network>& res, block_map& old2new, uint32_t index, unsigned phase, std::vector<uint32_t> const& genlib_to_cell )
  {
    auto const& node_data = node_match[index];
    auto const& best_cut = cuts[index][node_data.best_cut[phase]];
    auto const& gate = node_data.best_gate[phase]->root;

    /* permutate and negate to obtain the matched gate truth table */
    std::vector<signal<block_network>> children( gate->num_vars );

    auto ctr = 0u;
    for ( auto l : best_cut )
    {
      if ( ctr >= gate->num_vars )
        break;
      children[node_data.best_gate[phase]->permutation[ctr]] = old2new[l][( node_data.phase[phase] >> ctr ) & 1];
      ++ctr;
    }

    if ( !gate->is_super )
    {
      /* create the node */
      auto f = res.create_node( children, gate->function );
      res.add_cell( res.get_node( f ), genlib_to_cell.at( gate->root->id ) );

      /* add the node in the data structure */
      old2new[index][phase] = f;
    }
    else
    {
      /* supergate, create sub-gates */
      auto f = create_lut_for_gate2_rec( res, *gate, children, genlib_to_cell );

      /* add the node in the data structure */
      old2new[index][phase] = f;
    }
  }

  signal<block_network> create_lut_for_gate2_rec( cell_view<block_network>& res, composed_gate<NInputs> const& gate, std::vector<signal<block_network>> const& children, std::vector<uint32_t> const& genlib_to_cell )
  {
    std::vector<signal<block_network>> children_local( gate.fanin.size() );

    auto i = 0u;
    for ( auto const fanin : gate.fanin )
    {
      if ( fanin->root == nullptr )
      {
        /* terminal condition */
        children_local[i] = children[fanin->id];
      }
      else
      {
        children_local[i] = create_lut_for_gate2_rec( res, *fanin, children, genlib_to_cell );
      }
      ++i;
    }

    auto f = res.create_node( children_local, gate.root->function );
    res.add_cell( res.get_node( f ), genlib_to_cell.at( gate.root->id ) );
    return f;
  }

  void create_block_for_gate( cell_view<block_network>& res, block_map& old2new, uint32_t index, unsigned phase, std::vector<uint32_t> const& genlib_to_cell )
  {
    std::vector<standard_cell> const& lib = res.get_library();
    composed_gate<NInputs> const* local_gate = node_match[index].best_gate[phase]->root;
    standard_cell const& cell = lib[genlib_to_cell.at( local_gate->root->id )];

    assert( !local_gate->is_super );
    auto const& best_cut = cuts[index][node_match[index].best_cut[phase]];

    /* permutate and negate to obtain the matched gate truth table */
    std::vector<signal<block_network>> children( cell.gates.front().num_vars );

    /* output negations have already been assigned by the mapper */
    auto ctr = 0u;
    for ( auto l : best_cut )
    {
      if ( ctr >= local_gate->num_vars )
        break;
      children[node_match[index].best_gate[phase]->permutation[ctr]] = old2new[l][( node_match[index].phase[phase] >> ctr ) & 1];
      ++ctr;
    }

    multi_match_t const& tuple_data = multi_node_match[node_tuple_match[index].index][0];
    std::vector<uint32_t> outputs;
    std::vector<kitty::dynamic_truth_table> functions;

    /* re-order outputs to match the ones of the cell */
    for ( gate const& g : cell.gates )
    {
      /* find the correct node */
      for ( auto j = 0; j < max_multioutput_output_size; ++j )
      {
        uint32_t node_index = tuple_data[j].node_index;
        assert( node_match[node_index].same_match );
        uint8_t node_phase = node_match[node_index].best_gate[0] != nullptr ? 0 : 1;
        assert( node_match[node_index].multioutput_match[node_phase] );

        gate const* node_gate = node_match[node_index].best_gate[node_phase]->root->root;

        /* wrong output */
        if ( node_gate->id != g.id )
          continue;

        outputs.push_back( node_index );
        functions.push_back( g.function );
      }
    }

    assert( outputs.size() == cell.gates.size() );

    /* create the block */
    auto f = res.create_node( children, functions );
    res.add_cell( res.get_node( f ), genlib_to_cell.at( local_gate->root->id ) );

    for ( uint32_t s : outputs )
    {
      /* add inverted version if used */
      uint8_t node_phase = node_match[s].best_gate[0] != nullptr ? 0 : 1;
      assert( node_match[s].same_match );

      /* add the node in the data structure */
      old2new[s][node_phase] = f;

      if ( node_match[s].map_refs[node_phase ^ 1] > 0 )
      {
        old2new[s][node_phase ^ 1] = res.create_not( f );
        res.add_cell( res.get_node( old2new[s][node_phase ^ 1] ), genlib_to_cell.at( lib_inv_id ) );
      }

      f = res.next_output_pin( f );
    }
  }

  void clone_box( binding_view<klut_network>& res, klut_map& old2new, uint32_t index )
  {
    node<Ntk> n = ntk.index_to_node( index );
    std::vector<signal<klut_network>> children;

    ntk.foreach_fanin( n, [&]( auto const& f ) {
      children.push_back( old2new[ntk.get_node( f )][ntk.is_complemented( f ) ? 1 : 0] );
    } );

    /* create the node */
    auto const& tt = ntk.node_function( n );
    auto f = res.create_node( children, tt );

    /* add the node in the data structure */
    old2new[index][0] = f;
    if ( node_match[index].map_refs[1] )
    {
      old2new[index][1] = res.create_not( f );
      res.add_binding( res.get_node( old2new[index][1] ), lib_inv_id );
    }

    if constexpr ( has_has_binding_v<Ntk> )
    {
      if ( ntk.has_binding( n ) )
        res.add_binding( res.get_node( f ), ntk.get_binding_index( n ) );
    }
  }

  void clone_box2( cell_view<block_network>& res, klut_map& old2new, uint32_t index, std::vector<uint32_t> const& genlib_to_cell )
  {
    node<Ntk> n = ntk.index_to_node( index );
    std::vector<signal<block_network>> children;

    ntk.foreach_fanin( n, [&]( auto const& f ) {
      children.push_back( old2new[ntk.get_node( f )][ntk.is_complemented( f ) ? 1 : 0] );
    } );

    /* check if multi-output */
    std::vector<standard_cell> const& lib = res.get_library();
    if constexpr ( has_has_binding_v<Ntk> )
    {
      bool is_multioutput = false;
      if ( ntk.has_binding( n ) )
      {
        uint32_t cell_id = genlib_to_cell.at( ntk.get_binding_index( n ) );
        if ( lib.at( cell_id ).gates.size() > 1 )
          is_multioutput = true;
      }

      /* create the multioutput node (partially dangling) */
      if ( is_multioutput )
      {
        standard_cell const& cell = lib.at( genlib_to_cell.at( ntk.get_binding_index( n ) ) );
        std::vector<kitty::dynamic_truth_table> functions;
        for ( auto const& g : cell.gates )
        {
          functions.push_back( g.function );
        }

        auto f = res.create_node( children, functions );

        /* find and connect the correct pin */
        for ( auto const& g : cell.gates )
        {
          if ( g.id == cell.id )
            break;
          res.next_output_pin( f );
        }

        old2new[index][0] = f;
        res.add_cell( res.get_node( f ), cell.id );
        if ( node_match[index].map_refs[1] )
        {
          old2new[index][1] = res.create_not( f );
          res.add_cell( res.get_node( old2new[index][1] ), genlib_to_cell.at( lib_inv_id ) );
        }
        return;
      }
    }

    /* create the single-output node */
    auto const& tt = ntk.node_function( n );
    auto f = res.create_node( children, tt );

    /* add the node in the data structure */
    old2new[index][0] = f;
    if ( node_match[index].map_refs[1] )
    {
      old2new[index][1] = res.create_not( f );
      res.add_cell( res.get_node( old2new[index][1] ), genlib_to_cell.at( lib_inv_id ) );
    }

    if constexpr ( has_has_binding_v<Ntk> )
    {
      if ( ntk.has_binding( n ) )
        res.add_cell( res.get_node( f ), genlib_to_cell.at( ntk.get_binding_index( n ) ) );
    }
  }
#pragma endregion

#pragma region Cuts and matching utils
  void compute_cut_data( cut_t& cut, node<Ntk> const& n )
  {
    cut->delay = std::numeric_limits<uint32_t>::max();
    cut->flow = std::numeric_limits<float>::max();
    cut->ignore = false;

    if ( cut.size() > NInputs || cut.size() > 6 )
    {
      /* Ignore cuts too big to be mapped using the library */
      cut->ignore = true;
      return;
    }

    const auto tt = cut->function;
    const kitty::static_truth_table<6> fe = kitty::extend_to<6>( tt );
    auto fe_canon = fe;

    uint16_t negations_pos = 0;
    uint16_t negations_neg = 0;

    /* match positive polarity */
    if constexpr ( Configuration == classification_type::p_configurations )
    {
      auto canon = kitty::exact_n_canonization_support( fe, cut.size() );
      fe_canon = std::get<0>( canon );
      negations_pos = std::get<1>( canon );
    }

    auto const supergates_pos = library.get_supergates( fe_canon );

    /* match negative polarity */
    if constexpr ( Configuration == classification_type::p_configurations )
    {
      auto canon = kitty::exact_n_canonization_support( ~fe, cut.size() );
      fe_canon = std::get<0>( canon );
      negations_neg = std::get<1>( canon );
    }
    else
    {
      fe_canon = ~fe;
    }

    auto const supergates_neg = library.get_supergates( fe_canon );

    if ( supergates_pos != nullptr || supergates_neg != nullptr )
    {
      cut->supergates = { supergates_pos, supergates_neg };
      cut->negations = { negations_pos, negations_neg };
    }
    else
    {
      /* Ignore not matched cuts */
      cut->ignore = true;
      return;
    }

    /* compute cut cost based on LUT area */
    recompute_cut_data( cut, n );
  }

  void compute_cut_data_structural( cut_t& cut, node<Ntk> const& n )
  {
    cut->delay = std::numeric_limits<uint32_t>::max();
    cut->flow = std::numeric_limits<float>::max();
    cut->ignore = false;

    assert( cut.size() <= NInputs );

    const auto supergates_pos = library.get_supergates_pattern( cut->pattern_index, false );
    const auto supergates_neg = library.get_supergates_pattern( cut->pattern_index, true );

    if ( supergates_pos != nullptr || supergates_neg != nullptr )
    {
      cut->supergates = { supergates_pos, supergates_neg };
    }
    else
    {
      /* Ignore not matched cuts */
      cut->ignore = true;
      return;
    }

    /* compute cut cost based on LUT area */
    recompute_cut_data( cut, n );
  }

  void recompute_cut_data( cut_t& cut, node<Ntk> const& n )
  {
    /* compute cut cost based on LUT area */
    uint32_t best_arrival = 0;
    float best_area_flow = cut.size() > 1 ? cut.size() : 0;

    for ( auto leaf : cut )
    {
      const auto& best_leaf_cut = cuts[leaf][0];
      best_arrival = std::max( best_arrival, best_leaf_cut->delay );
      best_area_flow += best_leaf_cut->flow;
    }

    cut->delay = best_arrival + ( cut.size() > 1 ) ? 1 : 0;
    cut->flow = best_area_flow / ntk.fanout_size( n );
  }

  /* compute positions of leave indices in cut `sub` (subset) with respect to
   * leaves in cut `sup` (super set).
   *
   * Example:
   *   compute_truth_table_support( {1, 3, 6}, {0, 1, 2, 3, 6, 7} ) = {1, 3, 4}
   */
  void compute_truth_table_support( cut_t const& sub, cut_t const& sup, TT& tt )
  {
    size_t j = 0;
    auto itp = sup.begin();
    for ( auto i : sub )
    {
      itp = std::find( itp, sup.end(), i );
      lsupport[j++] = static_cast<uint8_t>( std::distance( sup.begin(), itp ) );
    }

    /* swap variables in the truth table */
    for ( int i = j - 1; i >= 0; --i )
    {
      assert( i <= lsupport[i] );
      kitty::swap_inplace( tt, i, lsupport[i] );
    }
  }

  void add_zero_cut( uint32_t index )
  {
    auto& cut = cuts[index].add_cut( &index, &index ); /* fake iterator for emptyness */
    cut->ignore = true;
    cut->delay = 0;
    cut->flow = 0;
    cut->pattern_index = 0;
    cut->negations[0] = cut->negations[1] = 0;
  }

  void add_unit_cut( uint32_t index )
  {
    auto& cut = cuts[index].add_cut( &index, &index + 1 );

    kitty::create_nth_var( cut->function, 0 );
    cut->ignore = true;
    cut->delay = 0;
    cut->flow = 0;
    cut->pattern_index = 1;
    cut->negations[0] = cut->negations[1] = 0;
  }

  inline void create_structural_cut( cut_t& new_cut, std::vector<cut_t const*> const& vcuts, uint32_t new_pattern, uint32_t pattern_id1, uint32_t pattern_id2 )
  {
    new_cut.set_leaves( *vcuts[0] );
    new_cut.add_leaves( vcuts[1]->begin(), vcuts[1]->end() );
    new_cut->pattern_index = new_pattern;

    /* get the polarity of the leaves of the new cut */
    uint16_t neg_l = 0, neg_r = 0;
    if ( ( *vcuts[0] )->pattern_index == 1 )
    {
      neg_r = static_cast<uint16_t>( pattern_id1 & 1 );
    }
    else
    {
      neg_r = ( *vcuts[0] )->negations[0];
    }
    if ( ( *vcuts[1] )->pattern_index == 1 )
    {
      neg_l = static_cast<uint16_t>( pattern_id2 & 1 );
    }
    else
    {
      neg_l = ( *vcuts[1] )->negations[0];
    }

    new_cut->negations[0] = ( neg_l << vcuts[0]->size() ) | neg_r;
    new_cut->negations[1] = new_cut->negations[0];
  }

  inline bool fast_support_minimization( TT const& tt, cut_t& res )
  {
    uint32_t support = 0u;
    uint32_t support_size = 0u;
    for ( uint32_t i = 0u; i < tt.num_vars(); ++i )
    {
      if ( kitty::has_var( tt, i ) )
      {
        support |= 1u << i;
        ++support_size;
      }
    }

    /* has not minimized support? */
    if ( ( support & ( support + 1u ) ) != 0u )
    {
      return false;
    }

    /* variables not in the support are the most significative */
    if ( support_size != res.size() )
    {
      std::vector<uint32_t> leaves( res.begin(), res.begin() + support_size );
      res.set_leaves( leaves.begin(), leaves.end() );
    }

    return true;
  }

  void compute_truth_table( uint32_t index, fanin_cut_t const& vcuts, uint32_t fanin, cut_t& res )
  {
    for ( uint32_t i = 0; i < fanin; ++i )
    {
      cut_t const* cut = vcuts[i];
      ltruth[i] = ( *cut )->function;
      compute_truth_table_support( *cut, res, ltruth[i] );
    }

    auto tt_res = ntk.compute( ntk.index_to_node( index ), ltruth.begin(), ltruth.begin() + fanin );

    if ( ps.cut_enumeration_ps.minimize_truth_table && !fast_support_minimization( tt_res, res ) )
    {
      const auto support = kitty::min_base_inplace( tt_res );

      std::vector<uint32_t> leaves_before( res.begin(), res.end() );
      std::vector<uint32_t> leaves_after( support.size() );

      auto it_support = support.begin();
      auto it_leaves = leaves_after.begin();
      while ( it_support != support.end() )
      {
        *it_leaves++ = leaves_before[*it_support++];
      }
      res.set_leaves( leaves_after.begin(), leaves_after.end() );
    }

    res->function = tt_res;
  }
#pragma endregion

  template<bool DO_AREA>
  inline bool compare_map( double arrival, double best_arrival, float area_flow, float best_area_flow, uint32_t size, uint32_t best_size )
  {
    if constexpr ( DO_AREA )
    {
      if ( area_flow < best_area_flow - epsilon )
      {
        return true;
      }
      else if ( area_flow > best_area_flow + epsilon )
      {
        return false;
      }
      else if ( arrival < best_arrival - epsilon )
      {
        return true;
      }
      else if ( arrival > best_arrival + epsilon )
      {
        return false;
      }
      return size < best_size;
    }
    else
    {
      if ( arrival < best_arrival - epsilon )
      {
        return true;
      }
      else if ( arrival > best_arrival + epsilon )
      {
        return false;
      }
      else if ( area_flow < best_area_flow - epsilon )
      {
        return true;
      }
      else if ( area_flow > best_area_flow + epsilon )
      {
        return false;
      }
      return size < best_size;
    }
  }

  double compute_switching_power()
  {
    double power = 0.0f;

    for ( auto const& n : topo_order )
    {
      const auto index = ntk.node_to_index( n );
      auto& node_data = node_match[index];

      if ( ntk.is_constant( n ) )
      {
        if ( node_data.best_gate[0] == nullptr && node_data.best_gate[1] == nullptr )
          continue;
      }
      else if ( ntk.is_pi( n ) )
      {
        if ( node_data.map_refs[1] > 0 )
          power += switch_activity[ntk.node_to_index( n )];
        continue;
      }

      /* continue if cut is not in the cover */
      if ( !node_data.map_refs[0] && !node_data.map_refs[1] )
        continue;

      unsigned phase = ( node_data.best_gate[0] != nullptr ) ? 0 : 1;

      if ( node_data.same_match || node_data.map_refs[phase] > 0 )
      {
        power += switch_activity[ntk.node_to_index( n )];

        if ( node_data.same_match && node_data.map_refs[phase ^ 1] > 0 )
          power += switch_activity[ntk.node_to_index( n )];
      }

      phase = phase ^ 1;
      if ( !node_data.same_match && node_data.map_refs[phase] > 0 )
      {
        power += switch_activity[ntk.node_to_index( n )];
      }
    }

    return power;
  }

#pragma region multioutput
  /* Experimental code */
  void compute_multioutput_match()
  {
    stopwatch t( st.time_multioutput );

    if ( library.num_multioutput_gates() == 0 )
      return;

    /* compute cuts: first simple method without proper matching */
    cut_enumeration_params multi_ps;
    multi_ps.minimize_truth_table = false;
    multi_cuts_t multi_cuts = fast_cut_enumeration<Ntk, max_multioutput_cut_size, true, cut_enumeration_emap_multi_cut>( ntk, multi_ps );

    /* cuts leaves classes */
    multi_hash_t multi_cuts_classes;
    multi_cuts_classes.reserve( 2000 );

    /* Multi-output matching */
    multi_enumerate_matches( multi_cuts, multi_cuts_classes );

    multi_single_matches_t multi_node_match_local;
    multi_node_match_local.reserve( multi_cuts_classes.size() );

    multi_compute_matches( multi_cuts, multi_cuts_classes, multi_node_match_local );

    if ( ps.remove_overlapping_multicuts )
      multi_filter_and_match<true>( multi_cuts, multi_node_match_local ); /* it also adds the tuple for node mapping */
    else
      multi_filter_and_match<false>( multi_cuts, multi_node_match_local ); /* it also adds the tuple for node mapping */
  }

  void multi_init_topo_order()
  {
    /* create and initialize a choice view to store the tuples */
    choice_view<Ntk> choice_ntk{ ntk };
    multi_add_choices( choice_ntk );

    ntk.incr_trav_id();
    ntk.incr_trav_id();

    /* add constants and CIs */
    const auto c0 = ntk.get_node( ntk.get_constant( false ) );
    topo_order.push_back( c0 );
    ntk.set_visited( c0, ntk.trav_id() );

    if ( const auto c1 = ntk.get_node( ntk.get_constant( true ) ); ntk.visited( c1 ) != ntk.trav_id() )
    {
      topo_order.push_back( c1 );
      ntk.set_visited( c1, ntk.trav_id() );
    }

    ntk.foreach_ci( [&]( auto const& n ) {
      if ( ntk.visited( n ) != ntk.trav_id() )
      {
        topo_order.push_back( n );
        ntk.set_visited( n, ntk.trav_id() );
      }
    } );

    /* sort topologically */
    ntk.foreach_co( [&]( auto const& f ) {
      if ( ntk.visited( ntk.get_node( f ) ) == ntk.trav_id() )
        return;
      multi_topo_sort_rec( choice_ntk, ntk.get_node( f ) );
    } );
  }

  /* Experimental code resticted to only half adders and full adders */
  void multi_enumerate_matches( multi_cuts_t const& multi_cuts, multi_hash_t& multi_cuts_classes )
  {
    static_assert( max_multioutput_cut_size > 1 && max_multioutput_cut_size < 7 );

    uint32_t counter = 0;
    multi_leaves_set_t leaves = { 0 };

    ntk.foreach_gate( [&]( auto const& n ) {
      uint32_t cut_index = 0;
      for ( auto& cut : multi_cuts.cuts( ntk.node_to_index( n ) ) )
      {
        kitty::static_truth_table<max_multioutput_cut_size> tt = multi_cuts.truth_table( *cut );
        /* reduce support for matching ID */
        uint64_t tt_id = ( cut->size() < 3 ) ? ( tt._bits & 0xF ) : tt._bits;
        uint64_t id = library.get_multi_function_id( tt_id );

        if ( !id )
        {
          ++cut_index;
          continue;
        }

        ( *cut )->data.id = id;

        multi_match_data data;
        data.node_index = ntk.node_to_index( n );
        data.cut_index = cut_index;
        leaves[2] = 0;
        uint32_t i = 0;
        for ( auto l : *cut )
          leaves[i++] = l;

        /* add to hash table */
        multi_cuts_classes[leaves].push_back( data );

        ++cut_index;
      }
    } );
  }

  /* Experimental code */
  void multi_compute_matches( multi_cuts_t const& multi_cuts, multi_hash_t& multi_cuts_classes, multi_single_matches_t& multi_node_match_local )
  {
    ntk.clear_values();

    /* copy set and sort by gate size: improve, too slow */
    std::vector<std::pair<multi_leaves_set_t, multi_output_set_t>> class_list;
    class_list.reserve( multi_cuts_classes.size() );
    for ( auto& it : multi_cuts_classes )
    {
      /* insert multiple occurring cuts */
      if ( it.second.size() > 1 )
        class_list.push_back( it );
    }

    std::stable_sort( class_list.begin(), class_list.end(), [&]( auto const& a, auto const& b ) {
      return a.first[2] > b.first[2];
    } );

    /* combine and match: specific code for 2-output cells */
    for ( auto it : class_list )
    {
      for ( uint32_t i = 0; i < it.second.size() - 1; ++i )
      {
        multi_match_data data_i = it.second[i];
        uint32_t index_i = data_i.node_index;
        uint32_t cut_index_i = data_i.cut_index;
        auto const& cut_i = multi_cuts.cuts( index_i )[cut_index_i];

        for ( uint32_t j = i + 1; j < it.second.size(); ++j )
        {
          multi_match_data data_j = it.second[j];
          uint32_t index_j = data_j.node_index;
          uint32_t cut_index_j = data_j.cut_index;
          auto const& cut_j = multi_cuts.cuts( index_j )[cut_index_j];

          /* not compatible -> TODO: change */
          if ( cut_i->data.id == cut_j->data.id )
            continue;

          /* check compatibility */
          if ( !multi_check_partally_dangling( index_i, index_j, cut_i ) )
            continue;

          multi_node_match_local.push_back( { data_i, data_j } );
        }
      }
    }
  }

  /* Experimental code */
  template<bool OverlapFilter>
  void multi_filter_and_match( multi_cuts_t const& multi_cuts, multi_single_matches_t const& multi_node_match_local )
  {
    multi_cut_set.reserve( multi_node_match_local.size() );
    multi_node_match.reserve( multi_node_match_local.size() );

    ntk.incr_trav_id();

    for ( auto& pair : multi_node_match_local )
    {
      uint32_t index1 = pair[0].node_index;
      uint32_t index2 = pair[1].node_index;
      uint32_t cut_index1 = pair[0].cut_index;
      uint32_t cut_index2 = pair[1].cut_index;
      multi_cut_t const& cut1 = multi_cuts.cuts( index1 )[cut_index1];
      multi_cut_t const& cut2 = multi_cuts.cuts( index2 )[cut_index2];

      assert( index1 < index2 );

      /* remove incompatible multi-output cuts */
      bool is_new = true;
      uint32_t insertion_index = multi_node_match.size();
      if constexpr ( OverlapFilter )
      {
        if ( multi_gate_check_overlapping( index1, index2, cut1 ) )
          continue;
      }
      else
      {
        if ( multi_gate_check_incompatible( index1, index2, is_new, insertion_index ) )
          continue;
        // if ( is_new && multi_gate_check_overlapping( index1, index2, cut1 ) )
        //   continue;
      }

      /* copy cuts */
      cut_t new_cut1, new_cut2;
      new_cut1.set_leaves( cut1.begin(), cut1.end() );
      new_cut2.set_leaves( cut2.begin(), cut2.end() );
      new_cut1->function = kitty::extend_to<6>( multi_cuts.truth_table( cut1 ) );
      new_cut2->function = kitty::extend_to<6>( multi_cuts.truth_table( cut2 ) );

      /* Multi-output Boolean matching, continue if no match */
      std::array<cut_t, max_multioutput_output_size> cut_pair = { new_cut1, new_cut2 };
      if ( !multi_compute_cut_data( cut_pair ) )
        continue;

      /* mark multioutput gate */
      if constexpr ( OverlapFilter )
      {
        multi_gate_mark_visited( index1, index2, cut1 );
        node_tuple_match[index1].has_info = 1;
        node_tuple_match[index1].lowest_index = 1;
        node_tuple_match[index1].index = multi_node_match.size();
        node_tuple_match[index2].has_info = 1;
        node_tuple_match[index2].highest_index = 1;
        node_tuple_match[index2].index = multi_node_match.size();
      }
      else
      {
        // multi_gate_mark_visited( index1, index2, cut1 );
        multi_gate_mark_compatibility( index1, index2, insertion_index );
      }

      /* add cut */
      multi_cut_set.push_back( cut_pair );

      /* re-index data */
      multi_match_data new_data1, new_data2;
      new_data1.node_index = index1;
      new_data1.cut_index = multi_cut_set.size() - 1;
      new_data2.node_index = index2;
      new_data2.cut_index = multi_cut_set.size() - 1;
      multi_match_t p = { new_data1, new_data2 };

      /* add cuts to the correct bucket */
      if ( is_new )
      {
        multi_node_match.push_back( { p } );
      }
      else
      {
        multi_node_match[insertion_index].push_back( p );
      }
    }
  }

  bool multi_compute_cut_data( std::array<cut_t, max_multioutput_output_size>& cut_tuple )
  {
    std::array<kitty::static_truth_table<6>, max_multioutput_output_size> tts;
    std::array<kitty::static_truth_table<6>, max_multioutput_output_size> tts_order;
    std::array<size_t, max_multioutput_output_size> order = {};
    std::array<uint16_t, max_multioutput_output_size> phase = { 0 };
    std::array<uint8_t, max_multioutput_output_size> phase_order;

    std::iota( order.begin(), order.end(), 0 );

    for ( auto i = 0; i < max_multioutput_output_size; ++i )
    {
      tts[i] = kitty::extend_to<6>( cut_tuple[i]->function );
      if ( ( tts[i]._bits & 1 ) == 1 )
      {
        tts[i] = ~tts[i];
        phase[i] = 1;
      }
    }

    std::stable_sort( order.begin(), order.end(), [&]( size_t a, size_t b ) {
      return tts[a] < tts[b];
    } );

    std::transform( order.begin(), order.end(), tts_order.begin(), [&]( size_t a ) {
      return tts[a];
    } );

    std::transform( order.begin(), order.end(), phase_order.begin(), [&]( uint8_t a ) {
      return phase[a];
    } );

    auto const multigates_match = library.get_multi_supergates( tts_order );

    /* Ignore not matched cuts */
    if ( multigates_match == nullptr )
      return false;

    /* add cut matches */
    for ( auto i = 0; i < max_multioutput_output_size; ++i )
    {
      cut_tuple[order[i]]->supergates[0] = nullptr;
      cut_tuple[order[i]]->supergates[1] = nullptr;
      cut_tuple[order[i]]->ignore = false;
      std::vector<supergate<NInputs>> const* multigate = &( ( *multigates_match )[i] );
      cut_tuple[order[i]]->supergates[phase_order[i]] = multigate;
    }

    return true;
  }

  inline bool multi_check_partally_dangling( uint32_t index1, uint32_t index2, multi_cut_t const& cut1 )
  {
    bool valid = true;

    /* check containment of cut1 in cut2 and viceversa */
    if ( index1 > index2 )
    {
      std::swap( index1, index2 );
    }

    ntk.foreach_fanin( ntk.index_to_node( index2 ), [&]( auto const& f ) {
      auto g = ntk.get_node( f );
      if ( ntk.node_to_index( g ) == index1 && ntk.fanout_size( g ) == 1 )
      {
        valid = false;
      }
      return valid;
    } );

    if ( !valid )
      return false;

    if ( !is_contained_mffc( ntk.index_to_node( index2 ), ntk.index_to_node( index1 ), cut1 ) )
      return false;

    return true;
  }

  inline bool multi_gate_check_overlapping( uint32_t index1, uint32_t index2, multi_cut_t const& cut )
  {
    bool contained = false;

    /* mark leaves */
    for ( auto leaf : cut )
    {
      ntk.incr_value( ntk.index_to_node( leaf ) );
    }

    contained = multi_mark_visited_rec<false>( ntk.index_to_node( index1 ) );
    contained |= multi_mark_visited_rec<false>( ntk.index_to_node( index2 ) );

    /* unmark leaves */
    for ( auto leaf : cut )
    {
      ntk.decr_value( ntk.index_to_node( leaf ) );
    }

    return contained;
  }

  inline bool multi_gate_check_incompatible( uint32_t index1, uint32_t index2, bool& is_new, uint32_t& data_index )
  {
    /* check cut assigned cut outputs, specialized code for 2 outputs */
    if ( !node_tuple_match[index1].has_info && !node_tuple_match[index2].has_info )
      return false;

    if ( node_tuple_match[index1].has_info && node_tuple_match[index2].has_info )
    {
      uint32_t current_assignment = node_tuple_match[index1].index;
      if ( current_assignment != node_tuple_match[index2].index )
        return true;
      is_new = false;
      data_index = current_assignment;
      return false;
    }

    return true;
  }

  inline void multi_gate_mark_compatibility( uint32_t index1, uint32_t index2, uint32_t mark_value )
  {
    node_tuple_match[index1].has_info = 1;
    node_tuple_match[index1].lowest_index = 1;
    node_tuple_match[index1].index = mark_value;
    node_tuple_match[index2].has_info = 1;
    node_tuple_match[index2].highest_index = 1;
    node_tuple_match[index2].index = mark_value;
  }

  inline void multi_gate_mark_visited( uint32_t index1, uint32_t index2, multi_cut_t const& cut )
  {
    /* mark leaves */
    for ( auto leaf : cut )
    {
      ntk.incr_value( ntk.index_to_node( leaf ) );
    }

    /* mark */
    multi_mark_visited_rec<true>( ntk.index_to_node( index1 ) );
    multi_mark_visited_rec<true>( ntk.index_to_node( index2 ) );

    /* unmark leaves */
    for ( auto leaf : cut )
    {
      ntk.decr_value( ntk.index_to_node( leaf ) );
    }
  }

  template<bool MARK>
  bool multi_mark_visited_rec( node<Ntk> const& n )
  {
    /* leaf */
    if ( ntk.value( n ) )
      return false;

    /* already visited */
    if ( ntk.visited( n ) == ntk.trav_id() )
      return true;

    if constexpr ( MARK )
    {
      ntk.set_visited( n, ntk.trav_id() );
    }

    bool contained = false;
    ntk.foreach_fanin( n, [&]( auto const& f ) {
      contained |= multi_mark_visited_rec<MARK>( ntk.get_node( f ) );

      if constexpr ( !MARK )
      {
        if ( contained )
          return false;
      }

      return true;
    } );

    return contained;
  }

  bool is_contained_mffc( node<Ntk> root, node<Ntk> n, multi_cut_t const& cut )
  {
    /* reference cut leaves */
    for ( auto leaf : cut )
    {
      ntk.incr_value( ntk.index_to_node( leaf ) );
    }

    bool valid = true;
    tmp_visited.clear();
    dereference_node_rec( root );

    if ( ntk.fanout_size( n ) == 0 )
      valid = false;

    for ( uint64_t g : tmp_visited )
      ntk.incr_fanout_size( ntk.index_to_node( g ) );

    /* dereference leaves */
    for ( auto leaf : cut )
    {
      ntk.decr_value( ntk.index_to_node( leaf ) );
    }

    return valid;
  }

  void dereference_node_rec( node<Ntk> const& n )
  {
    /* leaf */
    if ( ntk.value( n ) )
      return;

    ntk.foreach_fanin( n, [&]( auto const& f ) {
      node<Ntk> g = ntk.get_node( f );
      if ( ntk.decr_fanout_size( g ) == 0 )
      {
        dereference_node_rec( g );
      }
      tmp_visited.push_back( ntk.node_to_index( g ) );
    } );
  }

  void multi_add_choices( choice_view<Ntk>& choice_ntk )
  {
    for ( auto& field : multi_node_match )
    {
      auto& pair = field.front();
      uint32_t index1 = pair[0].node_index;
      uint32_t index2 = pair[1].node_index;
      uint32_t cut_index1 = pair[0].cut_index;
      cut_t const& cut = multi_cut_set[cut_index1][0];

      /* don't add choice if in TFI, set TFI bit */
      if ( multi_is_in_tfi( ntk.index_to_node( index2 ), ntk.index_to_node( index1 ), cut ) )
      {
        /* if there is a path of length > 1 linking node 1 and 2, save as TFI node */
        uint32_t in_tfi = multi_is_in_direct_tfi( ntk.index_to_node( index2 ), ntk.index_to_node( index1 ) ) ? 0 : 1;
        for ( auto& match : field )
          match[0].in_tfi = in_tfi;
        /* Propagate the TFI dependency to the cut boundary so DFS cannot temporarily mark an
           intermediate node before it visits the tuple's upper output. */
        multi_set_tfi_dependency( ntk.index_to_node( index2 ), ntk.index_to_node( index1 ), cut );
        continue;
      }

      choice_ntk.add_choice( ntk.index_to_node( index1 ), ntk.index_to_node( index2 ) );

      assert( choice_ntk.count_choices( ntk.index_to_node( index1 ) ) == 2 );
    }
  }

  bool multi_topo_sort_rec( choice_view<Ntk>& choice_ntk, node<Ntk> const& n )
  {
    /* is permanently marked? */
    if ( ntk.visited( n ) == ntk.trav_id() )
      return true;

    /* loop detected: backtrack to remove the cause */
    if ( ntk.visited( n ) == ntk.trav_id() - 1 )
      return false;

    /* get the representative (smallest index) */
    node<Ntk> repr = choice_ntk.get_choice_representative( n );

    /* loop detected: backtrack to remove the cause */
    if ( ntk.visited( repr ) == ntk.trav_id() - 1 )
      return false;

    /* solve the TFI dependency first */
    node<Ntk> dependency_node = ntk.index_to_node( ntk.value( n ) );
    if ( dependency_node > 0 && ntk.visited( dependency_node ) != ntk.trav_id() - 1 )
    {
      if ( !multi_topo_sort_rec( choice_ntk, dependency_node ) )
        return false;
      assert( ntk.visited( n ) == ntk.trav_id() );
      return true;
    }

    /* for all the choices */
    uint32_t i = 0;
    bool check = true;
    choice_ntk.foreach_choice( repr, [&]( auto const& g ) {
      /* ensure that the node is not visited or temporarily marked */
      assert( ntk.visited( g ) != ntk.trav_id() );
      assert( ntk.visited( g ) != ntk.trav_id() - 1 );

      /* mark node temporarily */
      ntk.set_visited( g, ntk.trav_id() - 1 );

      /* mark children */
      ntk.foreach_fanin( g, [&]( auto const& f ) {
        check = multi_topo_sort_rec( choice_ntk, ntk.get_node( f ) );
        return check;
      } );

      /* cycle detected: backtrack to the last choice jump */
      if ( !check )
      {
        /* revert visited */
        ntk.set_visited( g, ntk.trav_id() - 2 );
        if ( i > 0 && n == repr )
        {
          /* fix cycle: remove multi-output match */
          choice_ntk.foreach_choice( repr, [&]( auto const& p ) {
            node_tuple_match[ntk.node_to_index( p )].data = 0;
            return true;
          } );
          choice_ntk.remove_choice( g );
          check = true;
        }
        return false;
      }

      ++i;
      return true;
    } );

    if ( !check )
    {
      return false;
    }

    choice_ntk.foreach_choice( repr, [&]( auto const& g ) {
      /* ensure that the node is not visited */
      assert( ntk.visited( g ) != ntk.trav_id() );

      /* mark node n permanently */
      ntk.set_visited( g, ntk.trav_id() );

      /* visit node */
      topo_order.push_back( g );

      return true;
    } );

    return true;
  }

  inline bool multi_is_in_tfi( node<Ntk> const& root, node<Ntk> const& n, cut_t const& cut )
  {
    /* reference cut leaves */
    for ( auto leaf : cut )
    {
      ntk.incr_value( ntk.index_to_node( leaf ) );
    }

    ntk.incr_trav_id();
    multi_mark_visited_rec<true>( root );
    bool contained = ntk.visited( n ) == ntk.trav_id();

    /* dereference leaves */
    for ( auto leaf : cut )
    {
      ntk.decr_value( ntk.index_to_node( leaf ) );
    }

    return contained;
  }

  inline bool multi_is_in_direct_tfi( node<Ntk> const& root, node<Ntk> const& n )
  {
    bool contained = false;

    ntk.foreach_fanin( root, [&]( auto const& f ) {
      if ( ntk.get_node( f ) == n )
        contained = true;
    } );

    return contained;
  }

  inline void multi_set_tfi_dependency( node<Ntk> const& root, node<Ntk> const& n, cut_t const& cut )
  {
    /* reference cut leaves */
    for ( auto leaf : cut )
    {
      ntk.incr_value( ntk.index_to_node( leaf ) );
    }

    ntk.incr_trav_id();

    /* add a TFI dependencies */
    ntk.set_value( n, ntk.node_to_index( root ) );
    ntk.set_visited( n, ntk.trav_id() );
    multi_set_tfi_dependency_rec( root, ntk.node_to_index( root ) );

    /* reset root's dependency info */
    ntk.set_value( root, 0 );

    /* dereference leaves */
    for ( auto leaf : cut )
    {
      ntk.decr_value( ntk.index_to_node( leaf ) );
    }
  }

  void multi_set_tfi_dependency_rec( node<Ntk> const& n, uint32_t const dependency_info )
  {
    /* leaf */
    if ( ntk.value( n ) )
      return;

    /* already visited */
    if ( ntk.visited( n ) == ntk.trav_id() )
      return;

    ntk.set_visited( n, ntk.trav_id() );
    ntk.set_value( n, dependency_info );

    ntk.foreach_fanin( n, [&]( auto const& f ) {
      multi_set_tfi_dependency_rec( ntk.get_node( f ), dependency_info );
    } );
  }
#pragma endregion

private:
  Ntk const& ntk;
  tech_library<NInputs, Configuration> const& library;
  emap_params const& ps;
  emap_stats& st;

  uint32_t iteration{ 0 }; /* current mapping iteration */
  double delay{ 0.0f };    /* current delay of the mapping */
  double area{ 0.0f };     /* current area of the mapping */
  uint32_t inv{ 0 };       /* current inverter count */

  /* lib inverter info */
  float lib_inv_area;
  float lib_inv_delay;
  uint32_t lib_inv_id;

  /* electrical (load-aware) model: the smallest inverter's input capacitance and delay slope,
   * and the per-(node, phase) capacitive load the current cover imposes on each signal. The
   * loads are seeded from fanout counts (init_node_loads) and refreshed from the committed
   * cover after every round (update_node_loads). All zero unless ps.electrical_model. */
  float lib_inv_cap{ 0 };
  float lib_inv_slope{ 0 };
  std::vector<std::array<float, 2>> node_loads;
  // Only allocated for the optional nonlinear model; tracks real sink pins and POs.
  std::vector<std::array<uint32_t, 2>> node_fanouts;

  /* optional single-input covering stage (`cover_buffer` only; empty otherwise). One choice per
   * (node, phase): the selected buffer drive, or the identity (id == UINT32_MAX = no buffer). The
   * buffered node's DRIVER load becomes the buffer's input cap; the buffer itself carries the
   * sink load (kept in sink_load for legality reporting and the buffer's own drive selection). */
  struct buffer_choice
  {
    typename tech_library<NInputs, Configuration>::buffer_drive drive{};
    float sink_load{ 0.0f };
    bool active{ false };
  };
  std::vector<std::array<buffer_choice, 2>> node_buffers;
  bool buffer_stage{ false };

  /* bounded per-node, per-phase area-delay frontier (`curve_points` only; empty otherwise). The
   * frontier school is the memory-hungry one, so the budget is hard-capped. */
  static constexpr uint32_t max_curve_points = 8;
  static constexpr uint32_t max_frontier_points = 16;
  struct curve_set
  {
    std::array<best_gate_emap<NInputs>, max_curve_points> points{};
    uint8_t size{ 0 };
  };
  std::vector<std::array<curve_set, 2>> curves;
  uint32_t curve_budget{ 0 };

  /* load-indexed arrival curves (`load_points` only; both empty otherwise). One arrival per
   * (node, phase, ladder point): what that signal delivers when its output carries that load.
   * Unlike the area-delay frontier above this keeps no match identity — the node binds exactly
   * one gate, and the curve exists so that node's CONSUMERS can price themselves against the
   * load they each impose (ADR-0050). Eight doubles per phase is 128 B/node, an order of
   * magnitude under a frontier of the same width. */
  load_ladder ladder{};
  std::vector<std::array<std::array<double, max_load_points>, 2>> node_curves;

  /* lib buffer info */
  float lib_buf_area;
  float lib_buf_delay;
  uint32_t lib_buf_id;

  std::vector<node<Ntk>> topo_order;
  node_match_t node_match;
  std::vector<multioutput_info> node_tuple_match;
  std::vector<float> switch_activity;
  std::vector<uint64_t> tmp_visited;

  /* cut computation */
  std::vector<cut_set_t> cuts; /* compressed representation of cuts */
  cut_merge_t lcuts;           /* cut merger container */
  cut_set_t temp_cuts;         /* temporary cut set container */
  truth_compute_t ltruth;      /* truth table merger container */
  support_t lsupport;          /* support merger container */
  uint32_t cuts_total{ 0 };    /* current computed cuts */

  /* multi-output matching */
  multi_cut_set_t multi_cut_set;    /* set of multi-output cuts */
  multi_matches_t multi_node_match; /* matched multi-output gates */

  time_point time_begin;
};

} /* namespace detail */

/*! \brief Technology mapping.
 *
 * This function implements a technology mapping algorithm.
 *
 * The function takes the size of the cuts in the template parameter `CutSize`.
 *
 * The function returns a block network that supports multi-output cells.
 *
 * The novelties of this mapper are contained in 2 publications:
 * - A. Tempia Calvino and G. De Micheli, "Technology Mapping Using Multi-Output Library Cells," ICCAD, 2023.
 * - G. Radi, A. Tempia Calvino, and G. De Micheli, "In Medio Stat Virtus: Combining Boolean and Pattern Matching," ASP-DAC, 2024.
 *
 * **Required network functions:**
 * - `size`
 * - `is_pi`
 * - `is_constant`
 * - `node_to_index`
 * - `index_to_node`
 * - `get_node`
 * - `foreach_po`
 * - `foreach_node`
 * - `fanout_size`
 *
 * \param ntk Network
 * \param library Technology library
 * \param ps Mapping params
 * \param pst Mapping statistics
 *
 */
template<unsigned CutSize = 6u, class Ntk, unsigned NInputs, classification_type Configuration>
cell_view<block_network> emap( Ntk const& ntk, tech_library<NInputs, Configuration> const& library, emap_params const& ps = {}, emap_stats* pst = nullptr )
{
  static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
  static_assert( has_size_v<Ntk>, "Ntk does not implement the size method" );
  static_assert( has_is_pi_v<Ntk>, "Ntk does not implement the is_pi method" );
  static_assert( has_is_constant_v<Ntk>, "Ntk does not implement the is_constant method" );
  static_assert( has_node_to_index_v<Ntk>, "Ntk does not implement the node_to_index method" );
  static_assert( has_index_to_node_v<Ntk>, "Ntk does not implement the index_to_node method" );
  static_assert( has_get_node_v<Ntk>, "Ntk does not implement the get_node method" );
  static_assert( has_foreach_po_v<Ntk>, "Ntk does not implement the foreach_po method" );
  static_assert( has_foreach_node_v<Ntk>, "Ntk does not implement the foreach_node method" );
  static_assert( has_fanout_size_v<Ntk>, "Ntk does not implement the fanout_size method" );

  emap_stats st;
  detail::emap_impl<Ntk, CutSize, NInputs, Configuration> p( ntk, library, ps, st );
  auto res = p.run_block();

  if ( ps.verbose && !st.mapping_error )
  {
    st.report();
  }

  if ( pst )
  {
    *pst = st;
  }
  return res;
}

/*! \brief Technology mapping.
 *
 * This function implements a technology mapping algorithm.
 *
 * The function takes the size of the cuts in the template parameter `CutSize`.
 *
 * The function returns a k-LUT network. Each LUT abstacts a gate of the technology library.
 *
 * The novelties of this mapper are contained in 2 publications:
 * - A. Tempia Calvino and G. De Micheli, "Technology Mapping Using Multi-Output Library Cells," ICCAD, 2023.
 * - G. Radi, A. Tempia Calvino, and G. De Micheli, "In Medio Stat Virtus: Combining Boolean and Pattern Matching," ASP-DAC, 2024.
 *
 * **Required network functions:**
 * - `size`
 * - `is_pi`
 * - `is_constant`
 * - `node_to_index`
 * - `index_to_node`
 * - `get_node`
 * - `foreach_po`
 * - `foreach_node`
 * - `fanout_size`
 *
 * \param ntk Network
 * \param library Technology library
 * \param ps Mapping params
 * \param pst Mapping statistics
 *
 */
template<unsigned CutSize = 6u, class Ntk, unsigned NInputs, classification_type Configuration>
binding_view<klut_network> emap_klut( Ntk const& ntk, tech_library<NInputs, Configuration> const& library, emap_params const& ps = {}, emap_stats* pst = nullptr )
{
  static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
  static_assert( has_size_v<Ntk>, "Ntk does not implement the size method" );
  static_assert( has_is_pi_v<Ntk>, "Ntk does not implement the is_pi method" );
  static_assert( has_is_constant_v<Ntk>, "Ntk does not implement the is_constant method" );
  static_assert( has_node_to_index_v<Ntk>, "Ntk does not implement the node_to_index method" );
  static_assert( has_index_to_node_v<Ntk>, "Ntk does not implement the index_to_node method" );
  static_assert( has_get_node_v<Ntk>, "Ntk does not implement the get_node method" );
  static_assert( has_foreach_po_v<Ntk>, "Ntk does not implement the foreach_po method" );
  static_assert( has_foreach_node_v<Ntk>, "Ntk does not implement the foreach_node method" );
  static_assert( has_fanout_size_v<Ntk>, "Ntk does not implement the fanout_size method" );

  emap_stats st;
  detail::emap_impl<Ntk, CutSize, NInputs, Configuration> p( ntk, library, ps, st );
  auto res = p.run_klut();

  if ( ps.verbose && !st.mapping_error )
  {
    st.report();
  }

  if ( pst )
  {
    *pst = st;
  }
  return res;
}

/*! \brief Technology node mapping.
 *
 * This function implements a simple technology mapping algorithm.
 * The algorithm maps each node to the best implementation in the technology library.
 *
 * **Required network functions:**
 * - `size`
 * - `is_pi`
 * - `is_constant`
 * - `node_to_index`
 * - `index_to_node`
 * - `get_node`
 * - `foreach_po`
 * - `foreach_node`
 * - `fanout_size`
 * - `has_binding`
 *
 * \param ntk Network
 * \param library Technology library
 * \param ps Mapping params
 * \param pst Mapping statistics
 *
 */
template<unsigned CutSize = 6u, class Ntk, unsigned NInputs, classification_type Configuration>
binding_view<klut_network> emap_node_map( Ntk const& ntk, tech_library<NInputs, Configuration> const& library, emap_params const& ps = {}, emap_stats* pst = nullptr )
{
  static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
  static_assert( has_size_v<Ntk>, "Ntk does not implement the size method" );
  static_assert( has_is_pi_v<Ntk>, "Ntk does not implement the is_pi method" );
  static_assert( has_is_constant_v<Ntk>, "Ntk does not implement the is_constant method" );
  static_assert( has_node_to_index_v<Ntk>, "Ntk does not implement the node_to_index method" );
  static_assert( has_index_to_node_v<Ntk>, "Ntk does not implement the index_to_node method" );
  static_assert( has_get_node_v<Ntk>, "Ntk does not implement the get_node method" );
  static_assert( has_foreach_po_v<Ntk>, "Ntk does not implement the foreach_po method" );
  static_assert( has_foreach_node_v<Ntk>, "Ntk does not implement the foreach_node method" );
  static_assert( has_has_binding_v<Ntk>, "Ntk does not implement the has_binding method" );

  emap_stats st;
  detail::emap_impl<Ntk, CutSize, NInputs, Configuration> p( ntk, library, ps, st );
  auto res = p.run_node_map();

  if ( ps.verbose && !st.mapping_error )
  {
    st.report();
  }

  if ( pst )
  {
    *pst = st;
  }
  return res;
}

/*! \brief Technology node mapping.
 *
 * This function implements a simple technology mapping algorithm.
 * The algorithm maps each node to the first implementation in the technology library.
 *
 * The input must be a binding_view with the gates correctly loaded.
 *
 * **Required network functions:**
 * - `size`
 * - `is_pi`
 * - `is_constant`
 * - `node_to_index`
 * - `index_to_node`
 * - `get_node`
 * - `foreach_po`
 * - `foreach_node`
 * - `fanout_size`
 * - `has_binding`
 *
 * \param ntk Network
 *
 */
template<class Ntk>
void emap_load_mapping( Ntk& ntk )
{
  static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
  static_assert( has_size_v<Ntk>, "Ntk does not implement the size method" );
  static_assert( has_is_pi_v<Ntk>, "Ntk does not implement the is_pi method" );
  static_assert( has_is_constant_v<Ntk>, "Ntk does not implement the is_constant method" );
  static_assert( has_node_to_index_v<Ntk>, "Ntk does not implement the node_to_index method" );
  static_assert( has_index_to_node_v<Ntk>, "Ntk does not implement the index_to_node method" );
  static_assert( has_get_node_v<Ntk>, "Ntk does not implement the get_node method" );
  static_assert( has_foreach_po_v<Ntk>, "Ntk does not implement the foreach_po method" );
  static_assert( has_foreach_node_v<Ntk>, "Ntk does not implement the foreach_node method" );
  static_assert( has_has_binding_v<Ntk>, "Ntk does not implement the has_binding method" );

  /* build the library map */
  using lib_t = std::unordered_map<kitty::dynamic_truth_table, uint32_t, kitty::hash<kitty::dynamic_truth_table>>;
  lib_t tt_to_gate;

  for ( auto const& g : ntk.get_library() )
  {
    tt_to_gate[g.function] = g.id;
  }

  ntk.foreach_gate( [&]( auto const& n ) {
    if ( auto it = tt_to_gate.find( ntk.node_function( n ) ); it != tt_to_gate.end() )
    {
      ntk.add_binding( n, it->second );
    }
    else
    {
      std::cout << fmt::format( "[e] node mapping for node {} failed: no match in the tech library\n", ntk.node_to_index( n ) );
    }
  } );
}

} /* namespace mockturtle */
