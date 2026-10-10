/* Standalone bounded frontier allocation probe (not a Catch test).
 * Compile as documented in the qualification report; run fresh subprocesses for RSS,
 * because Linux ru_maxrss survives exec. Hooks measure live C++ allocations using glibc's
 * malloc_usable_size; standard malloc/OS allocations outside C++ new are not counted.
 */
#include <algorithm>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <malloc.h>
#include <mockturtle/algorithms/emap.hpp>
#include <mockturtle/generators/arithmetic.hpp>
#include <mockturtle/io/genlib_reader.hpp>
#include <sys/resource.h>
static size_t live_bytes, peak_bytes;
void* operator new( std::size_t n )
{
  void* p = std::malloc( n ? n : 1 );
  if ( !p )
    throw std::bad_alloc();
  live_bytes += malloc_usable_size( p );
  peak_bytes = std::max( peak_bytes, live_bytes );
  return p;
}
void operator delete( void* p ) noexcept
{
  if ( p )
  {
    live_bytes -= malloc_usable_size( p );
    std::free( p );
  }
}
void* operator new[]( std::size_t n ) { return ::operator new( n ); }
void operator delete[]( void* p ) noexcept { ::operator delete( p ); }
void operator delete( void* p, std::size_t ) noexcept { ::operator delete( p ); }
void operator delete[]( void* p, std::size_t ) noexcept { ::operator delete( p ); }
using namespace mockturtle;
int main( int argc, char** argv )
try
{
  if ( argc < 3 || argc > 5 || ( argc > 3 && std::string( argv[3] ) != "ties" && std::string( argv[3] ) != "plain" ) ||
       ( argc > 4 && std::string( argv[4] ) != "multiply" && std::string( argv[4] ) != "carry" ) )
  {
    std::cerr << "usage: probe bits budget_mib [ties|plain] [multiply|carry]\n";
    return 2;
  }
  std::size_t end_bits = 0, end_budget = 0;
  auto const parsed_bits = std::stoull( argv[1], &end_bits );
  double const budget = std::stod( argv[2], &end_budget );
  if ( parsed_bits == 0 || parsed_bits > 8192 || end_bits != std::string( argv[1] ).size() ||
       end_budget != std::string( argv[2] ).size() || !std::isfinite( budget ) ||
       ( budget < 0 && budget != -1 ) || ( argc > 4 && std::string( argv[4] ) == "multiply" && parsed_bits > 32 ) )
    return 2;
  uint32_t const bits = parsed_bits;
  bool const ties = argc > 3 && std::string( argv[3] ) == "ties";
  bool const multiply = argc > 4 && std::string( argv[4] ) == "multiply";
  bool const carry_only = argc > 4 && std::string( argv[4] ) == "carry";
  std::string genlib = "GATE inv 1 O=!a; PIN * INV 1 999 0.9 0.3 0.9 0.3\n"
                       "GATE nand 2 O=!(a*b); PIN * INV 1 999 1.0 0.2 1.0 0.2\n"
                       "GATE and 3 O=a*b; PIN * INV 1 999 1.7 0.2 1.7 0.2\n"
                       "GATE xor 4 O=a^b; PIN * UNKNOWN 2 999 1.9 0.5 1.9 0.5\n"
                       "GATE maj 3 O=a*b+a*c+b*c; PIN * INV 1 999 2.0 0.2 2.0 0.2\n"
                       "GATE buf 2 O=a; PIN * NONINV 1 999 1.0 0 1.0 0\n"
                       "GATE zero 0 O=CONST0;\nGATE one 0 O=CONST1;\n";
  std::vector<gate> gates;
  std::istringstream in( genlib );
  if ( lorina::read_genlib( in, genlib_reader( gates ) ) != lorina::return_code::success )
    return 2;
  if ( ties )
  {
    /* Equal matches plus just-under/just-over epsilon costs/delays, permuted in insertion order. */
    auto original = gates;
    for ( uint32_t variant = 0; variant < 4; ++variant )
      for ( auto g : original )
      {
        g.name += std::to_string( variant );
        g.id = gates.size();
        g.area += variant * 0.0004;
        for ( auto& p : g.pins )
        {
          p.rise_block_delay += variant * 0.0004;
          p.fall_block_delay += variant * 0.0004;
        }
        gates.push_back( g );
      }
  }
  tech_library<3> lib( gates );
  aig_network aig;
  std::vector<aig_network::signal> a( bits ), b( bits );
  std::generate( a.begin(), a.end(), [&] { return aig.create_pi(); } );
  std::generate( b.begin(), b.end(), [&] { return aig.create_pi(); } );
  if ( multiply )
  {
    for ( auto o : carry_ripple_multiplier( aig, a, b ) )
      aig.create_po( o );
  }
  else
  {
    auto carry = aig.get_constant( false );
    carry_ripple_adder_inplace( aig, a, b, carry );
    if ( !carry_only )
      for ( auto s : a )
        aig.create_po( s );
    aig.create_po( carry );
  }
  emap_params ps;
  ps.frontier_points = budget < 0 ? 0 : 8;
  ps.frontier_rounds = 1;
  if ( carry_only )
    ps.ela_rounds = 0;
  ps.frontier_mem_budget_mb = budget < 0 ? 1024 : budget;
  auto entering = live_bytes;
  peak_bytes = entering;
  emap_stats st;
  auto mapped = emap_klut( aig, lib, ps, &st );
  auto peak = peak_bytes;
  struct rusage r;
  if ( getrusage( RUSAGE_SELF, &r ) != 0 )
    return 2;
  std::cout << std::setprecision( 17 ) << "nodes=" << aig.size() << " points=" << ps.frontier_points << " budget=" << budget << " ran=" << st.frontier_ran << " kept=" << st.frontier_kept << " estimate=" << st.frontier_bytes << " points_bytes=" << st.frontier_bytes_points
            << " headers_bytes=" << st.frontier_bytes_headers
            << " saved_matches_bytes=" << st.frontier_bytes_saved_matches
            << " sharing_bytes=" << st.frontier_bytes_sharing
            << " resolution_bytes=" << st.frontier_bytes_resolution
            << " scratch_bytes=" << st.frontier_bytes_scratch
            << " allocator_bytes=" << st.frontier_bytes_allocator
            << " allocation_peak=" << peak << " mapping_increment=" << peak - entering << " rss_peak_kib=" << r.ru_maxrss << " area=" << st.area << " delay=" << st.delay << " cut_changes=" << st.frontier_cut_changes << " gate_changes=" << st.frontier_gate_changes << " attempt_area=" << st.frontier_area_after << " attempt_delay=" << st.frontier_delay_after << "\n";
  if ( ties )
    mapped.foreach_node( [&]( auto n ) {
      if ( mapped.has_binding( n ) )
        std::cout << mapped.node_to_index( n ) << ":" << mapped.get_binding_index( n ) << ";";
    } );
}
catch ( std::exception const& error )
{
  std::cerr << "probe failed: " << error.what() << "\n";
  return 2;
}
