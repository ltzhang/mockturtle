/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2023  EPFL
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
  \file load_curve.hpp
  \brief Load ladder and arrival-curve evaluation for load-indexed technology mapping.

  A load-aware mapper that keeps ONE arrival number per node has priced that node at one
  assumed output load. Its consumers then read that number no matter what capacitance they
  themselves present, so the one electrical term the mapper could know exactly at match time —
  a candidate's own input pin capacitance — is the term it drops.

  The remedy is to store arrival as a FUNCTION of the node's output load, sampled on a small
  ladder of load points, and to evaluate a leaf at the capacitance the querying candidate
  actually presents. This file holds the two pure pieces of that: how the ladder is built from
  a library's capacitance range, and how a per-point value vector is read at an arbitrary load.
  Both are free functions over plain data so they can be reasoned about and tested without a
  mapper, a network, or a library.
*/

#pragma once

#include <array>
#include <cmath>
#include <cstdint>

namespace mockturtle
{

/*! \brief Hard cap on ladder points.
 *
 * The per-node store costs `points x 2 x sizeof(point)` bytes, so this cap is a memory bound
 * before it is anything else. Sixteen points span a decade-wide library range at better than
 * 20% per step, which is finer than the accuracy of the linearized delay model they index.
 */
static constexpr uint32_t max_load_points = 16u;

/*! \brief The load points a curve is sampled at, ascending.
 *
 * `size == 0` means the mechanism is off; every evaluation then returns the caller's scalar
 * fallback, so a disabled ladder cannot perturb a single number.
 */
struct load_ladder
{
  std::array<float, max_load_points> points{};
  uint32_t size{ 0 };

  /*! \brief Smallest sampled load (0 on an empty ladder). */
  float lo() const { return size == 0 ? 0.0f : points[0]; }
  /*! \brief Largest sampled load (0 on an empty ladder). */
  float hi() const { return size == 0 ? 0.0f : points[size - 1]; }
};

/*! \brief Build a ladder of @p points loads spanning [@p lo, @p hi] geometrically.
 *
 * Both anchors sit exactly ON the ladder: the smallest realizable load and the largest legal
 * one are the two a match most needs priced exactly, not interpolated. The interior is
 * geometric rather than uniform because a cell library's capacitance range spans decades — a
 * uniform ladder would put nearly every point near the ceiling, where the choices no longer
 * differ.
 *
 * The ladder is a function of the library alone: not of the design, not of a round's state,
 * not of iteration order. Two constructions from the same anchors are bit-identical.
 *
 * Degenerate anchors (equal, inverted, or a non-positive floor) are not an error and not a
 * guess. The ladder collapses to the single load it can defend, which makes the whole
 * mechanism reduce exactly to the scalar model rather than inventing a range. @p points == 0
 * returns an empty ladder (the feature off); anything above `max_load_points` is clamped.
 */
inline load_ladder build_load_ladder( float lo, float hi, uint32_t points )
{
  load_ladder L{};
  if ( points == 0 )
    return L;
  if ( points > max_load_points )
    points = max_load_points;

  if ( !( lo > 0.0f ) || !( hi > lo ) )
  {
    /* prefer the floor when it is realizable — it is the anchor every sink can present */
    L.points[0] = lo > 0.0f ? lo : ( hi > 0.0f ? hi : 0.0f );
    L.size = 1;
    return L;
  }

  if ( points == 1 )
  {
    /* the balanced single sample of a geometric range */
    L.points[0] = static_cast<float>( std::sqrt( static_cast<double>( lo ) * static_cast<double>( hi ) ) );
    L.size = 1;
    return L;
  }

  double const ratio = static_cast<double>( hi ) / static_cast<double>( lo );
  for ( uint32_t k = 0; k < points; ++k )
  {
    double const t = static_cast<double>( k ) / static_cast<double>( points - 1 );
    L.points[k] = static_cast<float>( static_cast<double>( lo ) * std::pow( ratio, t ) );
  }
  /* pin the anchors so no rounding can move them off the range they define */
  L.points[0] = lo;
  L.points[points - 1] = hi;
  L.size = points;
  return L;
}

/*! \brief Read the per-point values @p v of a curve on ladder @p L at output load @p load.
 *
 * Exact at the knots; linear between them. Outside the range the two directions are
 * deliberately asymmetric:
 *
 * - **Below the floor: clamp.** The floor is the smallest load a real sink can present, so a
 *   query under it is a rounding artifact. Extrapolating downward would invent a delay smaller
 *   than any the library can deliver.
 * - **Above the ceiling: extrapolate on the last segment.** Clamping there would price an
 *   over-loaded driver as if load stopped mattering — optimistic exactly where the cover must
 *   be pessimistic. Delay is linear in load, so the last segment carries the honest slope.
 *
 * An empty ladder (the mechanism off) returns @p fallback unchanged.
 */
inline double interpolate_curve( load_ladder const& L, std::array<double, max_load_points> const& v,
                                 float load, double fallback = 0.0 )
{
  if ( L.size == 0 )
    return fallback;
  if ( L.size == 1 || load <= L.points[0] )
    return v[0];

  uint32_t const n = L.size;
  if ( load >= L.points[n - 1] )
  {
    float const span = L.points[n - 1] - L.points[n - 2];
    if ( !( span > 0.0f ) )
      return v[n - 1];
    double const slope = ( v[n - 1] - v[n - 2] ) / static_cast<double>( span );
    return v[n - 1] + slope * ( static_cast<double>( load ) - static_cast<double>( L.points[n - 1] ) );
  }

  /* at most 16 points: a scan beats a search and keeps the order of evaluation fixed */
  uint32_t k = 1;
  while ( k + 1 < n && L.points[k] < load )
    ++k;

  float const span = L.points[k] - L.points[k - 1];
  if ( !( span > 0.0f ) )
    return v[k];
  double const t = ( static_cast<double>( load ) - static_cast<double>( L.points[k - 1] ) ) / static_cast<double>( span );
  return v[k - 1] + t * ( v[k] - v[k - 1] );
}

/*! \brief Does candidate (@p arrival, @p area, @p size) beat the incumbent at one ladder point?
 *
 * The store keeps ONE match per ladder point — the best at that load — rather than a Pareto set,
 * which is what makes its size exactly `points x 2` and its maintenance a comparison rather than
 * a search. @p has_current false means the point is empty and any candidate takes it.
 *
 * Arrival decides first: this is the load axis, and at a fixed load the earlier signal is the
 * better match. Ties (within @p eps, since two drives of one function differ in the last bits of
 * a linearized delay) fall to area, then to cut size — a smaller cut references fewer leaves and
 * so costs less to keep alive. Every remaining tie keeps the INCUMBENT, so the outcome depends
 * only on the candidate enumeration order the mapper already fixes, never on which of two equal
 * matches was offered last (ADR-0042).
 */
inline bool load_curve_prefer( double arrival, float area, uint32_t size, bool has_current,
                               double cur_arrival, float cur_area, uint32_t cur_size,
                               double eps = 1e-9 )
{
  if ( !has_current )
    return true;
  if ( arrival < cur_arrival - eps )
    return true;
  if ( arrival > cur_arrival + eps )
    return false;
  if ( area < cur_area )
    return true;
  if ( area > cur_area )
    return false;
  return size < cur_size;
}

/*! \brief Bytes a load-indexed match store costs: one point per (node, phase, ladder point).
 *
 * Computed BEFORE allocating (ADR-0050 clause 6) so an over-budget request is refused with a
 * number rather than discovered at the OOM killer. Returns 0 when either dimension is 0.
 */
inline uint64_t load_curve_store_bytes( uint64_t nodes, uint32_t points, uint64_t point_bytes )
{
  return nodes * 2ull * static_cast<uint64_t>( points ) * point_bytes;
}

} // namespace mockturtle
