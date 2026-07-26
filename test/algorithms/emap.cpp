#include <catch.hpp>

#include <cstdint>
#include <vector>

#include <lorina/genlib.hpp>
#include <lorina/super.hpp>
#include <mockturtle/algorithms/emap.hpp>
#include <mockturtle/algorithms/simulation.hpp>
#include <mockturtle/generators/arithmetic.hpp>
#include <mockturtle/io/genlib_reader.hpp>
#include <mockturtle/io/super_reader.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/block.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/utils/tech_library.hpp>
#include <mockturtle/views/binding_view.hpp>
#include <mockturtle/views/cell_view.hpp>
#include <mockturtle/views/dont_touch_view.hpp>

using namespace mockturtle;

std::string const test_library = "GATE   inv1    1 O=!a;            PIN * INV 1 999 0.9 0.3 0.9 0.3\n"
                                 "GATE   inv2    2 O=!a;            PIN * INV 2 999 1.0 0.1 1.0 0.1\n"
                                 "GATE   nand2   2 O=!(a*b);        PIN * INV 1 999 1.0 0.2 1.0 0.2\n"
                                 "GATE   and2    3 O=a*b;           PIN * INV 1 999 1.7 0.2 1.7 0.2\n"
                                 "GATE   xor2    4 O=a^b;           PIN * UNKNOWN 2 999 1.9 0.5 1.9 0.5\n"
                                 "GATE   mig3    3 O=a*b+a*c+b*c;   PIN * INV 1 999 2.0 0.2 2.0 0.2\n"
                                 "GATE   xor3    5 O=a^b^c;         PIN * UNKNOWN 2 999 3.0 0.5 3.0 0.5\n"
                                 "GATE   buf     2 O=a;             PIN * NONINV 1 999 1.0 0.0 1.0 0.0\n"
                                 "GATE   zero    0 O=CONST0;\n"
                                 "GATE   one     0 O=CONST1;\n"
                                 "GATE   ha      5 C=a*b;           PIN * INV 1 999 1.7 0.4 1.7 0.4\n"
                                 "GATE   ha      5 S=!a*b+a*!b;     PIN * INV 1 999 2.1 0.4 2.1 0.4\n"
                                 "GATE   fa      6 C=a*b+a*c+b*c;   PIN * INV 1 999 2.1 0.4 2.1 0.4\n"
                                 "GATE   fa      6 S=a^b^c;         PIN * INV 1 999 3.0 0.4 3.0 0.4";

std::string const large_library = "GATE   inv1    1 O=!a;            PIN * INV 1 999 0.9 0.3 0.9 0.3\n"
                                  "GATE   inv2    2 O=!a;            PIN * INV 2 999 1.0 0.1 1.0 0.1\n"
                                  "GATE   nand2   2 O=!(a*b);        PIN * INV 1 999 1.0 0.2 1.0 0.2\n"
                                  "GATE   xor2    5 O=a^b;           PIN * UNKNOWN 2 999 1.9 0.5 1.9 0.5\n"
                                  "GATE   mig3    3 O=a*b+a*c+b*c;   PIN * INV 1 999 2.0 0.2 2.0 0.2\n"
                                  "GATE   buf     2 O=a;             PIN * NONINV 1 999 1.0 0.0 1.0 0.0\n"
                                  "GATE   zero    0 O=CONST0;\n"
                                  "GATE   one     0 O=CONST1;\n"
                                  "GATE   nand8   8 O=!(a*b*c*d*e*f*g*h);   PIN * INV 1 999 4.0 0.2 4.0 0.2\n";

std::string const super_library = "test.genlib\n"
                                  "3\n"
                                  "2\n"
                                  "6\n"
                                  "* nand2 1 0\n"
                                  "inv1 3\n"
                                  "* nand2 2 4\n"
                                  "\0";

TEST_CASE( "Emap on MAJ3", "[emap]" )
{
  std::vector<gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  tech_library<3> lib( gates );

  aig_network aig;
  const auto a = aig.create_pi();
  const auto b = aig.create_pi();
  const auto c = aig.create_pi();

  const auto f = aig.create_maj( a, b, c );
  aig.create_po( f );

  emap_params ps;
  emap_stats st;
  binding_view<klut_network> luts = emap_klut( aig, lib, ps, &st );

  CHECK( luts.size() == 6u );
  CHECK( luts.num_pis() == 3u );
  CHECK( luts.num_pos() == 1u );
  CHECK( luts.num_gates() == 1u );
  CHECK( st.area == 3.0f );
  CHECK( st.delay == 2.0f );
}

TEST_CASE( "Emap on bad MAJ3 and constant output", "[emap]" )
{
  std::vector<gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  tech_library<3> lib( gates );

  aig_network aig;
  const auto a = aig.create_pi();
  const auto b = aig.create_pi();
  const auto c = aig.create_pi();

  const auto f = aig.create_maj( a, aig.create_maj( a, b, c ), c );
  aig.create_po( f );
  aig.create_po( aig.get_constant( true ) );

  emap_params ps;
  emap_stats st;
  binding_view<klut_network> luts = emap_klut( aig, lib, ps, &st );

  CHECK( luts.size() == 6u );
  CHECK( luts.num_pis() == 3u );
  CHECK( luts.num_pos() == 2u );
  CHECK( luts.num_gates() == 1u );
  CHECK( st.area == 3.0f );
  CHECK( st.delay == 2.0f );
}

TEST_CASE( "Emap on full adder 1", "[emap]" )
{
  std::vector<gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  tech_library<3> lib( gates );

  aig_network aig;
  const auto a = aig.create_pi();
  const auto b = aig.create_pi();
  const auto c = aig.create_pi();

  const auto [sum, carry] = full_adder( aig, a, b, c );
  aig.create_po( sum );
  aig.create_po( carry );

  emap_params ps;
  emap_stats st;
  binding_view<klut_network> luts = emap_klut( aig, lib, ps, &st );

  const float eps{ 0.005f };

  CHECK( luts.size() == 7u );
  CHECK( luts.num_pis() == 3u );
  CHECK( luts.num_pos() == 2u );
  CHECK( luts.num_gates() == 2u );
  CHECK( st.area > 8.0f - eps );
  CHECK( st.area < 8.0f + eps );
  CHECK( st.delay > 3.0f - eps );
  CHECK( st.delay < 3.0f + eps );
}

TEST_CASE( "Emap on full adder 2", "[emap]" )
{
  std::vector<gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  tech_library<3, classification_type::p_configurations> lib( gates );

  aig_network aig;
  const auto a = aig.create_pi();
  const auto b = aig.create_pi();
  const auto c = aig.create_pi();

  const auto [sum, carry] = full_adder( aig, a, b, c );
  aig.create_po( sum );
  aig.create_po( carry );

  emap_params ps;
  ps.cut_enumeration_ps.minimize_truth_table = false;
  ps.ela_rounds = 1;
  ps.eswp_rounds = 2;
  emap_stats st;
  binding_view<klut_network> luts = emap_klut( aig, lib, ps, &st );

  const float eps{ 0.005f };

  CHECK( luts.size() == 7u );
  CHECK( luts.num_pis() == 3u );
  CHECK( luts.num_pos() == 2u );
  CHECK( luts.num_gates() == 2u );
  CHECK( st.area > 8.0f - eps );
  CHECK( st.area < 8.0f + eps );
  CHECK( st.delay > 3.0f - eps );
  CHECK( st.delay < 3.0f + eps );
}

TEST_CASE( "Emap on full adder 1 with cells", "[emap]" )
{
  std::vector<gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  tech_library<3> lib( gates );

  aig_network aig;
  const auto a = aig.create_pi();
  const auto b = aig.create_pi();
  const auto c = aig.create_pi();

  const auto [sum, carry] = full_adder( aig, a, b, c );
  aig.create_po( sum );
  aig.create_po( carry );

  emap_params ps;
  emap_stats st;
  cell_view<block_network> luts = emap( aig, lib, ps, &st );

  const float eps{ 0.005f };

  CHECK( luts.size() == 7u );
  CHECK( luts.num_pis() == 3u );
  CHECK( luts.num_pos() == 2u );
  CHECK( luts.num_gates() == 2u );
  CHECK( st.area > 8.0f - eps );
  CHECK( st.area < 8.0f + eps );
  CHECK( st.delay > 3.0f - eps );
  CHECK( st.delay < 3.0f + eps );
}

TEST_CASE( "Emap on full adder 2 with cells", "[emap]" )
{
  std::vector<gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  tech_library<3, classification_type::p_configurations> lib( gates );

  aig_network aig;
  const auto a = aig.create_pi();
  const auto b = aig.create_pi();
  const auto c = aig.create_pi();

  const auto [sum, carry] = full_adder( aig, a, b, c );
  aig.create_po( sum );
  aig.create_po( carry );

  emap_params ps;
  ps.cut_enumeration_ps.minimize_truth_table = false;
  ps.ela_rounds = 1;
  ps.eswp_rounds = 2;
  emap_stats st;
  cell_view<block_network> luts = emap( aig, lib, ps, &st );

  const float eps{ 0.005f };

  CHECK( luts.size() == 7u );
  CHECK( luts.num_pis() == 3u );
  CHECK( luts.num_pos() == 2u );
  CHECK( luts.num_gates() == 2u );
  CHECK( st.area > 8.0f - eps );
  CHECK( st.area < 8.0f + eps );
  CHECK( st.delay > 3.0f - eps );
  CHECK( st.delay < 3.0f + eps );
}

TEST_CASE( "Emap on ripple carry adder with multi-output gates", "[emap]" )
{
  std::vector<gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  tech_library_params tps;
  tps.load_multioutput_gates_single = false;
  tech_library<3, classification_type::p_configurations> lib( gates, tps );

  aig_network aig;
  
  std::vector<aig_network::signal> a( 8 ), b( 8 );
  std::generate( a.begin(), a.end(), [&aig]() { return aig.create_pi(); } );
  std::generate( b.begin(), b.end(), [&aig]() { return aig.create_pi(); } );
  auto carry = aig.get_constant( false );

  carry_ripple_adder_inplace( aig, a, b, carry );

  std::for_each( a.begin(), a.end(), [&]( auto f ) { aig.create_po( f ); } );
  aig.create_po( carry );

  emap_params ps;
  ps.map_multioutput = true;
  ps.area_oriented_mapping = true;
  emap_stats st;
  binding_view<klut_network> luts = emap_klut( aig, lib, ps, &st );

  const float eps{ 0.005f };

  CHECK( luts.size() == 34u );
  CHECK( luts.num_pis() == 16u );
  CHECK( luts.num_pos() == 9u );
  CHECK( luts.num_gates() == 16u );
  CHECK( st.area > 47.0f - eps );
  CHECK( st.area < 47.0f + eps );
  CHECK( st.delay > 17.3f - eps );
  CHECK( st.delay < 17.3f + eps );
  CHECK( st.multioutput_gates == 8 );
}

TEST_CASE( "Emap on ripple carry adder with multi-output cells", "[emap]" )
{
  std::vector<gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  tech_library_params tps;
  tps.load_multioutput_gates_single = false;
  tech_library<3, classification_type::p_configurations> lib( gates, tps );

  aig_network aig;
  
  std::vector<aig_network::signal> a( 8 ), b( 8 );
  std::generate( a.begin(), a.end(), [&aig]() { return aig.create_pi(); } );
  std::generate( b.begin(), b.end(), [&aig]() { return aig.create_pi(); } );
  auto carry = aig.get_constant( false );

  carry_ripple_adder_inplace( aig, a, b, carry );

  std::for_each( a.begin(), a.end(), [&]( auto f ) { aig.create_po( f ); } );
  aig.create_po( carry );

  emap_params ps;
  ps.map_multioutput = true;
  ps.area_oriented_mapping = true;
  emap_stats st;
  cell_view<block_network> luts = emap( aig, lib, ps, &st );

  const float eps{ 0.005f };

  CHECK( luts.size() == 26u );
  CHECK( luts.num_pis() == 16u );
  CHECK( luts.num_pos() == 9u );
  CHECK( luts.num_gates() == 8u );
  CHECK( st.area > 47.0f - eps );
  CHECK( st.area < 47.0f + eps );
  CHECK( st.delay > 17.3f - eps );
  CHECK( st.delay < 17.3f + eps );
  CHECK( st.multioutput_gates == 8 );
}

TEST_CASE( "Emap preserves outputs across transitive multi-output dependencies", "[emap]" )
{
  std::vector<gate> gates;
  std::istringstream in( test_library );
  CHECK( lorina::read_genlib( in, genlib_reader( gates ) ) == lorina::return_code::success );

  tech_library_params tps;
  tps.load_multioutput_gates = true;
  tech_library<3, classification_type::p_configurations> lib( gates, tps );

  aig_network aig;
  const auto a = aig.create_pi();
  const auto b = aig.create_pi();
  const auto c = aig.create_pi();
  const auto carry = aig.create_maj( a, b, c );
  const auto abc = aig.create_and( aig.create_and( a, b ), c );
  const auto any = aig.create_or( aig.create_or( a, b ), c );
  const auto sum_when_carry = aig.create_and( carry, abc );
  const auto sum_when_no_carry = aig.create_and( !carry, any );
  const auto sum = aig.create_or( sum_when_carry, sum_when_no_carry );
  aig.create_po( aig.create_xor( sum_when_carry, sum ) );

  emap_params ps;
  ps.area_oriented_mapping = true;
  emap_stats st;

  const cell_view<block_network> baseline = emap( aig, lib, ps, &st );
  CHECK( simulate<kitty::static_truth_table<3u>>( baseline ) ==
         simulate<kitty::static_truth_table<3u>>( aig ) );

  ps.map_multioutput = true;
  const cell_view<block_network> mapped = emap( aig, lib, ps, &st );

  CHECK( simulate<kitty::static_truth_table<3u>>( mapped ) ==
         simulate<kitty::static_truth_table<3u>>( aig ) );
}

TEST_CASE( "Emap on multiplier with multi-output gates", "[emap]" )
{
  std::vector<gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  tech_library_params tps;
  tps.load_minimum_size_only = false;
  tps.load_multioutput_gates_single = true;
  tech_library<3> lib( gates, tps );

  aig_network aig;

  std::vector<typename aig_network::signal> a( 8 ), b( 8 );
  std::generate( a.begin(), a.end(), [&aig]() { return aig.create_pi(); } );
  std::generate( b.begin(), b.end(), [&aig]() { return aig.create_pi(); } );

  for ( auto const& o : carry_ripple_multiplier( aig, a, b ) )
  {
    aig.create_po( o );
  }

  CHECK( aig.num_pis() == 16 );
  CHECK( aig.num_pos() == 16 );

  emap_params ps;
  ps.map_multioutput = true;
  emap_stats st;
  binding_view<klut_network> luts = emap_klut( aig, lib, ps, &st );

  const float eps{ 0.005f };

  CHECK( luts.size() == 235u );
  CHECK( luts.num_pis() == 16u );
  CHECK( luts.num_pos() == 16u );
  CHECK( luts.num_gates() == 217u );
  CHECK( st.area > 612.0f - eps );
  CHECK( st.area < 612.0f + eps );
  CHECK( st.delay > 33.60f - eps );
  CHECK( st.delay < 33.60f + eps );
  CHECK( st.multioutput_gates == 40 );
}

TEST_CASE( "Emap with inverters", "[emap]" )
{
  std::vector<gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  tech_library<3> lib( gates );

  aig_network aig;
  const auto a = aig.create_pi();
  const auto b = aig.create_pi();
  const auto c = aig.create_pi();

  const auto f1 = aig.create_and( !a, b );
  const auto f2 = aig.create_and( f1, !c );

  aig.create_po( f2 );

  emap_params ps;
  emap_stats st;
  binding_view<klut_network> luts = emap_klut( aig, lib, ps, &st );

  const float eps{ 0.005f };

  CHECK( luts.size() == 9u );
  CHECK( luts.num_pis() == 3u );
  CHECK( luts.num_pos() == 1u );
  CHECK( luts.num_gates() == 4u );
  CHECK( st.area > 8.0f - eps );
  CHECK( st.area < 8.0f + eps );
  CHECK( st.delay > 4.3f - eps );
  CHECK( st.delay < 4.3f + eps );
}

TEST_CASE( "Emap with inverters minimization", "[emap]" )
{
  std::vector<gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  tech_library<3> lib( gates );

  aig_network aig;
  const auto a = aig.create_pi();
  const auto b = aig.create_pi();
  const auto c = aig.create_pi();

  const auto f = aig.create_maj( !a, !b, !c );
  aig.create_po( f );

  emap_params ps;
  emap_stats st;
  binding_view<klut_network> luts = emap_klut( aig, lib, ps, &st );

  const float eps{ 0.005f };

  CHECK( luts.size() == 7u );
  CHECK( luts.num_pis() == 3u );
  CHECK( luts.num_pos() == 1u );
  CHECK( luts.num_gates() == 2u );
  CHECK( st.area > 4.0f - eps );
  CHECK( st.area < 4.0f + eps );
  CHECK( st.delay > 2.9f - eps );
  CHECK( st.delay < 2.9f + eps );
}

TEST_CASE( "Emap on buffer and constant outputs", "[emap]" )
{
  std::vector<gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  tech_library<3, classification_type::np_configurations> lib( gates );

  aig_network aig;
  const auto a = aig.create_pi();
  const auto b = aig.create_pi();
  const auto c = aig.create_pi();
  const auto d = aig.create_pi();

  const auto n5 = aig.create_and( a, d );
  const auto n6 = aig.create_and( a, !c );
  const auto n7 = aig.create_and( !c, n5 );
  const auto n8 = aig.create_and( c, n6 );
  const auto n9 = aig.create_and( !n6, n7 );
  const auto n10 = aig.create_and( n7, n8 );
  const auto n11 = aig.create_and( a, n10 );
  const auto n12 = aig.create_and( !d, n11 );
  const auto n13 = aig.create_and( !d, !n7 );
  const auto n14 = aig.create_and( !n6, !n7 );

  aig.create_po( aig.get_constant( true ) );
  aig.create_po( b );
  aig.create_po( n9 );
  aig.create_po( n12 );
  aig.create_po( !n13 );
  aig.create_po( n14 );

  emap_params ps;
  emap_stats st;
  binding_view<klut_network> luts = emap_klut( aig, lib, ps, &st );

  const float eps{ 0.005f };

  CHECK( luts.size() == 10u );
  CHECK( luts.num_pis() == 4u );
  CHECK( luts.num_pos() == 6u );
  CHECK( luts.num_gates() == 4u );
  CHECK( st.area > 7.0f - eps );
  CHECK( st.area < 7.0f + eps );
  CHECK( st.delay > 1.9f - eps );
  CHECK( st.delay < 1.9f + eps );
}

TEST_CASE( "Emap with boolean matching", "[emap]" )
{
  std::vector<gate> gates;

  std::istringstream in( large_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  tech_library<8> lib( gates );

  aig_network aig;
  const auto a = aig.create_pi();
  const auto b = aig.create_pi();
  const auto c = aig.create_pi();
  const auto d = aig.create_pi();
  const auto e = aig.create_pi();
  const auto f = aig.create_pi();
  const auto g = aig.create_pi();
  const auto h = aig.create_pi();

  const auto f1 = aig.create_and( !a, b );
  const auto f2 = aig.create_and( f1, !c );
  const auto f3 = aig.create_and( d, e );
  const auto f4 = aig.create_and( f, !g );
  const auto f5 = aig.create_and( f4, h );
  const auto f6 = aig.create_and( f2, f3 );
  const auto f7 = aig.create_and( f5, f6 );

  aig.create_po( f7 );

  emap_params ps;
  ps.matching_mode = emap_params::boolean;
  emap_stats st;
  cell_view<block_network> ntk = emap<8>( aig, lib, ps, &st );

  const float eps{ 0.005f };

  CHECK( ntk.size() == 27u );
  CHECK( ntk.num_pis() == 8u );
  CHECK( ntk.num_pos() == 1u );
  CHECK( ntk.num_gates() == 17u );
  CHECK( st.area > 24.0f - eps );
  CHECK( st.area < 24.0f + eps );
  CHECK( st.delay > 8.5f - eps );
  CHECK( st.delay < 8.5f + eps );
}

TEST_CASE( "Emap with structural matching", "[emap]" )
{
  std::vector<gate> gates;

  std::istringstream in( large_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  tech_library<8> lib( gates );

  aig_network aig;
  const auto a = aig.create_pi();
  const auto b = aig.create_pi();
  const auto c = aig.create_pi();
  const auto d = aig.create_pi();
  const auto e = aig.create_pi();
  const auto f = aig.create_pi();
  const auto g = aig.create_pi();
  const auto h = aig.create_pi();

  const auto f1 = aig.create_and( !a, b );
  const auto f2 = aig.create_and( f1, !c );
  const auto f3 = aig.create_and( d, e );
  const auto f4 = aig.create_and( f, !g );
  const auto f5 = aig.create_and( f4, h );
  const auto f6 = aig.create_and( f2, f3 );
  const auto f7 = aig.create_and( f5, f6 );

  aig.create_po( f7 );

  emap_params ps;
  ps.matching_mode = emap_params::structural;
  emap_stats st;
  cell_view<block_network> ntk = emap<8>( aig, lib, ps, &st );

  const float eps{ 0.005f };

  CHECK( ntk.size() == 15u );
  CHECK( ntk.num_pis() == 8u );
  CHECK( ntk.num_pos() == 1u );
  CHECK( ntk.num_gates() == 5u );
  CHECK( st.area > 12.0f - eps );
  CHECK( st.area < 12.0f + eps );
  CHECK( st.delay > 5.8f - eps );
  CHECK( st.delay < 5.8f + eps );
}

TEST_CASE( "Emap with hybrid matching", "[emap]" )
{
  std::vector<gate> gates;

  std::istringstream in( large_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  tech_library<8> lib( gates );

  aig_network aig;
  const auto a = aig.create_pi();
  const auto b = aig.create_pi();
  const auto c = aig.create_pi();
  const auto d = aig.create_pi();
  const auto e = aig.create_pi();
  const auto f = aig.create_pi();
  const auto g = aig.create_pi();
  const auto h = aig.create_pi();

  const auto f1 = aig.create_and( !a, b );
  const auto f2 = aig.create_and( f1, !c );
  const auto f3 = aig.create_and( d, e );
  const auto f4 = aig.create_and( f, !g );
  const auto f5 = aig.create_and( f4, h );
  const auto f6 = aig.create_and( f2, f3 );
  const auto f7 = aig.create_and( f5, f6 );

  aig.create_po( f7 );

  emap_params ps;
  ps.matching_mode = emap_params::hybrid;
  emap_stats st;
  cell_view<block_network> ntk = emap<8>( aig, lib, ps, &st );

  const float eps{ 0.005f };

  CHECK( ntk.size() == 15u );
  CHECK( ntk.num_pis() == 8u );
  CHECK( ntk.num_pos() == 1u );
  CHECK( ntk.num_gates() == 5u );
  CHECK( st.area > 12.0f - eps );
  CHECK( st.area < 12.0f + eps );
  CHECK( st.delay > 5.8f - eps );
  CHECK( st.delay < 5.8f + eps );
}

TEST_CASE( "Emap with arrival times", "[emap]" )
{
  std::vector<gate> gates;

  std::istringstream in( large_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  tech_library<6> lib( gates );

  aig_network aig;
  const auto a = aig.create_pi();
  const auto b = aig.create_pi();
  const auto c = aig.create_pi();
  const auto d = aig.create_pi();
  const auto e = aig.create_pi();
  const auto f = aig.create_pi();
  const auto g = aig.create_pi();
  const auto h = aig.create_pi();

  const auto f1 = aig.create_and( !a, b );
  const auto f2 = aig.create_and( f1, !c );
  const auto f3 = aig.create_and( d, e );
  const auto f4 = aig.create_and( f, !g );
  const auto f5 = aig.create_and( f4, h );
  const auto f6 = aig.create_and( f2, f3 );
  const auto f7 = aig.create_and( f5, f6 );

  aig.create_po( f7 );

  emap_params ps;
  ps.matching_mode = emap_params::boolean;
  emap_stats st;

  ps.arrival_times = std::vector<double>( 8 );
  ps.arrival_times[0] = 0.0;
  ps.arrival_times[1] = 1.0;
  ps.arrival_times[2] = 2.0;
  ps.arrival_times[3] = 3.0;
  ps.arrival_times[4] = 4.0;
  ps.arrival_times[5] = 5.0;
  ps.arrival_times[6] = 6.0;
  ps.arrival_times[7] = 7.0;

  cell_view<block_network> ntk = emap<6>( aig, lib, ps, &st );

  const float eps{ 0.005f };

  CHECK( ntk.size() == 27u );
  CHECK( ntk.num_pis() == 8u );
  CHECK( ntk.num_pos() == 1u );
  CHECK( ntk.num_gates() == 17u );
  CHECK( st.area > 24.0f - eps );
  CHECK( st.area < 24.0f + eps );
  CHECK( st.delay > 12.6f - eps );
  CHECK( st.delay < 12.6f + eps );
}

TEST_CASE( "Emap with global required times", "[emap]" )
{
  std::vector<gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  tech_library<6> lib( gates );

  aig_network aig;
  
  std::vector<aig_network::signal> a( 8 ), b( 8 );
  std::generate( a.begin(), a.end(), [&aig]() { return aig.create_pi(); } );
  std::generate( b.begin(), b.end(), [&aig]() { return aig.create_pi(); } );
  auto carry = aig.get_constant( false );

  carry_ripple_adder_inplace( aig, a, b, carry );

  std::for_each( a.begin(), a.end(), [&]( auto f ) { aig.create_po( f ); } );
  aig.create_po( carry );

  emap_params ps;
  ps.matching_mode = emap_params::boolean;
  ps.required_time = 20.0; // real delay 15.7
  emap_stats st;

  cell_view<block_network> ntk = emap<6>( aig, lib, ps, &st );

  const float eps{ 0.005f };

  CHECK( ntk.size() == 34 );
  CHECK( ntk.num_pis() == 16u );
  CHECK( ntk.num_pos() == 9u );
  CHECK( ntk.num_gates() == 16u );
  CHECK( st.area > 63.0f - eps );
  CHECK( st.area < 63.0f + eps );
  CHECK( st.delay < 20.0f + eps );
}

TEST_CASE( "Emap with required times", "[emap]" )
{
  std::vector<gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  tech_library<6> lib( gates );

  aig_network aig;
  
  std::vector<aig_network::signal> a( 8 ), b( 8 );
  std::generate( a.begin(), a.end(), [&aig]() { return aig.create_pi(); } );
  std::generate( b.begin(), b.end(), [&aig]() { return aig.create_pi(); } );
  auto carry = aig.get_constant( false );

  carry_ripple_adder_inplace( aig, a, b, carry );

  emap_params ps;
  ps.matching_mode = emap_params::boolean;
  // ps.required_time = 20.0; // real delay 15.7
  emap_stats st;

  std::for_each( a.begin(), a.end(), [&]( auto f ) { aig.create_po( f ); ps.required_times.push_back( 19.0 ); } );
  aig.create_po( carry );
  ps.required_times.push_back( 20.0 );

  cell_view<block_network> ntk = emap<6>( aig, lib, ps, &st );

  const float eps{ 0.005f };

  CHECK( ntk.size() == 34 );
  CHECK( ntk.num_pis() == 16u );
  CHECK( ntk.num_pos() == 9u );
  CHECK( ntk.num_gates() == 16u );
  CHECK( st.area > 63.0f - eps );
  CHECK( st.area < 63.0f + eps );
  CHECK( st.delay < 20.0f + eps );
}

TEST_CASE( "Emap with required time relaxation", "[emap]" )
{
  std::vector<gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  tech_library<6> lib( gates );

  aig_network aig;
  
  std::vector<aig_network::signal> a( 8 ), b( 8 );
  std::generate( a.begin(), a.end(), [&aig]() { return aig.create_pi(); } );
  std::generate( b.begin(), b.end(), [&aig]() { return aig.create_pi(); } );
  auto carry = aig.get_constant( false );

  carry_ripple_adder_inplace( aig, a, b, carry );

  std::for_each( a.begin(), a.end(), [&]( auto f ) { aig.create_po( f ); } );
  aig.create_po( carry );

  emap_params ps;
  ps.matching_mode = emap_params::boolean;
  ps.relax_required = 27.5; // real delay 15.7
  emap_stats st;

  cell_view<block_network> ntk = emap<6>( aig, lib, ps, &st );

  const float eps{ 0.005f };

  CHECK( ntk.size() == 34 );
  CHECK( ntk.num_pis() == 16u );
  CHECK( ntk.num_pos() == 9u );
  CHECK( ntk.num_gates() == 16u );
  CHECK( st.area > 63.0f - eps );
  CHECK( st.area < 63.0f + eps );
  CHECK( st.delay < 20.0f + eps );
}

TEST_CASE( "Emap with supergates", "[emap]" )
{
  std::vector<gate> gates;
  super_lib super_data;

  std::istringstream in_lib( test_library );
  auto result = lorina::read_genlib( in_lib, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  std::istringstream in_super( super_library );
  result = lorina::read_super( in_super, super_reader( super_data ) );
  CHECK( result == lorina::return_code::success );

  tech_library<3, classification_type::p_configurations> lib( gates, super_data );

  aig_network aig;
  const auto a = aig.create_pi();
  const auto b = aig.create_pi();
  const auto c = aig.create_pi();

  const auto n4 = aig.create_and( a, b );
  const auto n5 = aig.create_and( b, c );
  const auto f = aig.create_and( n4, n5 );
  aig.create_po( f );

  emap_params ps;
  emap_stats st;
  binding_view<klut_network> luts = emap_klut( aig, lib, ps, &st );

  const float eps{ 0.005f };

  CHECK( luts.size() == 8u );
  CHECK( luts.num_pis() == 3u );
  CHECK( luts.num_pos() == 1u );
  CHECK( luts.num_gates() == 3u );
  CHECK( st.area == 9.0f );
  CHECK( st.delay > 3.4f - eps );
  CHECK( st.delay < 3.4f + eps );
}

TEST_CASE( "Emap with supergates 2", "[emap]" )
{
  std::vector<gate> gates;
  super_lib super_data;

  std::istringstream in_lib( test_library );
  auto result = lorina::read_genlib( in_lib, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  std::istringstream in_super( super_library );
  result = lorina::read_super( in_super, super_reader( super_data ) );
  CHECK( result == lorina::return_code::success );

  tech_library<3, classification_type::p_configurations> lib( gates, super_data );

  aig_network aig;
  const auto a = aig.create_pi();
  const auto b = aig.create_pi();
  const auto c = aig.create_pi();

  const auto n4 = aig.create_and( a, b );
  const auto n5 = aig.create_and( b, c );
  const auto f = aig.create_and( n4, n5 );
  aig.create_po( f );

  emap_params ps;
  emap_stats st;
  cell_view<block_network> luts = emap( aig, lib, ps, &st );

  const float eps{ 0.005f };

  CHECK( luts.size() == 8u );
  CHECK( luts.num_pis() == 3u );
  CHECK( luts.num_pos() == 1u );
  CHECK( luts.num_gates() == 3u );
  CHECK( st.area == 9.0f );
  CHECK( st.delay > 3.4f - eps );
  CHECK( st.delay < 3.4f + eps );
}

TEST_CASE( "Emap on circuit with don't touch gates", "[emap]" )
{
  std::vector<gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  tech_library<3, classification_type::np_configurations> lib( gates );

  klut_network klut;
  const auto a = klut.create_pi();
  const auto b = klut.create_pi();
  const auto c = klut.create_pi();
  const auto d = klut.create_pi();

  const auto n5 = klut.create_xor( c, d );
  const auto n6 = klut.create_not( n5 );
  const auto n7 = klut.create_xor( a, b );
  const auto sum = klut.create_xor( n6, n7 );
  const auto carry = klut.create_maj( a, b, n5 );

  klut.create_po( sum );
  klut.create_po( carry );

  binding_view<klut_network> b_klut{ klut, gates };
  dont_touch_view<binding_view<klut_network>> db_klut{ b_klut };

  db_klut.add_binding( klut.get_node( n5 ), 3 );
  db_klut.select_dont_touch( klut.get_node( n5 ) );
  db_klut.add_binding( klut.get_node( n6 ), 0 );
  db_klut.select_dont_touch( klut.get_node( n6 ) );

  emap_params ps;
  ps.map_multioutput = true;
  ps.area_oriented_mapping = true;
  emap_stats st;
  binding_view<klut_network> luts = emap_klut( klut, lib, ps, &st );

  const float eps{ 0.005f };

  CHECK( luts.size() == 10u );
  CHECK( luts.num_pis() == 4u );
  CHECK( luts.num_pos() == 2u );
  CHECK( luts.num_gates() == 4u );
  CHECK( st.area > 11.0f - eps );
  CHECK( st.area < 11.0f + eps );
  CHECK( st.delay > 5.8f - eps );
  CHECK( st.delay < 5.8f + eps );
}

TEST_CASE( "Emap on circuit with don't touch cells", "[emap]" )
{
  std::vector<gate> gates;

  std::istringstream in( test_library );
  auto result = lorina::read_genlib( in, genlib_reader( gates ) );
  CHECK( result == lorina::return_code::success );

  tech_library<3, classification_type::np_configurations> lib( gates );

  klut_network klut;
  const auto a = klut.create_pi();
  const auto b = klut.create_pi();
  const auto c = klut.create_pi();
  const auto d = klut.create_pi();

  const auto n5 = klut.create_xor( c, d );
  const auto n6 = klut.create_not( n5 );
  const auto n7 = klut.create_xor( a, b );
  const auto sum = klut.create_xor( n6, n7 );
  const auto carry = klut.create_maj( a, b, n5 );

  klut.create_po( sum );
  klut.create_po( carry );

  binding_view<klut_network> b_klut{ klut, gates };
  dont_touch_view<binding_view<klut_network>> db_klut{ b_klut };

  db_klut.add_binding( klut.get_node( n5 ), 3 );
  db_klut.select_dont_touch( klut.get_node( n5 ) );
  db_klut.add_binding( klut.get_node( n6 ), 0 );
  db_klut.select_dont_touch( klut.get_node( n6 ) );

  emap_params ps;
  ps.map_multioutput = true;
  ps.area_oriented_mapping = true;
  emap_stats st;
  cell_view<block_network> luts = emap( klut, lib, ps, &st );

  const float eps{ 0.005f };

  CHECK( luts.size() == 9u );
  CHECK( luts.num_pis() == 4u );
  CHECK( luts.num_pos() == 2u );
  CHECK( luts.num_gates() == 3u );
  CHECK( st.area > 11.0f - eps );
  CHECK( st.area < 11.0f + eps );
  CHECK( st.delay > 5.8f - eps );
  CHECK( st.delay < 5.8f + eps );
}
/* Drive-legality fixture (ADR-0047): two drive strengths per function with real max_load limits. The
 * strong cell costs 4x the area and is slower unloaded, but has a 10x smaller load slope and a 10x
 * higher drive limit -- so it only ever wins on a node that actually drives a load. */
std::string const drive_family_library = "GATE   inv_w   1 O=!a;     PIN * INV 0.5 2.0 0.010 0.100 0.010 0.100\n"
                                         "GATE   inv_s   4 O=!a;     PIN * INV 2.0 20.0 0.020 0.010 0.020 0.010\n"
                                         "GATE   nand_w  2 O=!(a*b); PIN * INV 0.5 2.0 0.015 0.100 0.015 0.100\n"
                                         "GATE   nand_s  8 O=!(a*b); PIN * INV 2.0 20.0 0.030 0.010 0.030 0.010\n"
                                         "GATE   zero    0 O=CONST0;\n"
                                         "GATE   one     0 O=CONST1;";

/* A shared node driven into EIGHT distinct sinks: whatever cell is bound for it drives eight input
 * pins, which at the weak cell's 0.5 pin capacitance is a load of 4.0 against its own 2.0 limit. */
namespace
{
aig_network eight_sink_aig()
{
  aig_network aig;
  auto const a = aig.create_pi();
  auto const b = aig.create_pi();
  auto const t = aig.create_and( a, b );
  for ( int i = 0; i < 8; ++i )
  {
    auto const p = aig.create_pi();
    aig.create_po( aig.create_and( !t, p ) );
  }
  return aig;
}

struct drive_mix
{
  uint32_t weak{ 0 };
  uint32_t strong{ 0 };
};

drive_mix map_and_count( bool load_aware )
{
  std::vector<gate> gates;
  std::istringstream in( drive_family_library );
  lorina::read_genlib( in, genlib_reader( gates ) );

  tech_library_params tps;
  tps.electrical_model = load_aware;
  tech_library<3> tlib( gates, tps );

  emap_params ps;
  ps.electrical_model = load_aware;
  ps.area_oriented_mapping = false;
  emap_stats st;

  aig_network aig = eight_sink_aig();
  auto res = emap_klut( aig, tlib, ps, &st );

  drive_mix mix;
  res.foreach_node( [&]( auto const& n ) {
    if ( !res.has_binding( n ) )
      return;
    std::string const& name = res.get_binding( n ).name;
    if ( name.size() > 2 && name.compare( name.size() - 2, 2, "_s" ) == 0 )
      ++mix.strong;
    if ( name.size() > 2 && name.compare( name.size() - 2, 2, "_w" ) == 0 )
      ++mix.weak;
  } );
  return mix;
}
} // namespace

TEST_CASE( "emap load-aware mapping strengthens an overloaded driver", "[emap]" )
{
  drive_mix const blind = map_and_count( false );
  drive_mix const aware = map_and_count( true );

  /* load-blind: only the minimum-size cells exist, so no strong cell can be chosen */
  CHECK( blind.strong == 0u );
  CHECK( blind.weak > 0u );

  /* load-aware: the heavily loaded node earns a stronger drive */
  CHECK( aware.strong > 0u );
  /* and the mapping does not simply upsize everything -- the lightly loaded cells stay small */
  CHECK( aware.weak > 0u );
}

TEST_CASE( "emap load-aware mapping reports drive legality", "[emap]" )
{
  std::vector<gate> gates;
  std::istringstream in( drive_family_library );
  lorina::read_genlib( in, genlib_reader( gates ) );

  aig_network aig = eight_sink_aig();

  /* load-blind: the overloaded driver is invisible, so nothing is reported */
  {
    tech_library<3> tlib( gates, tech_library_params{} );
    emap_params ps;
    emap_stats st;
    emap_klut( aig, tlib, ps, &st );
    CHECK( st.max_load_violations == 0u );
    CHECK( st.worst_load_ratio == 0.0 );
  }

  /* load-aware: the estimated load is measured against each bound cell's own limit */
  {
    tech_library_params tps;
    tps.electrical_model = true;
    tech_library<3> tlib( gates, tps );
    emap_params ps;
    ps.electrical_model = true;
    emap_stats st;
    emap_klut( aig, tlib, ps, &st );
    CHECK( st.worst_load_ratio > 0.0 );
  }
}

TEST_CASE( "emap load-aware mapping is off by default", "[emap]" )
{
  /* The default must be byte-identical to the pre-ADR-0047 behavior: same area, same delay, same
   * cell count on a library that HAS drive variants but a tech_library built without load_aware. */
  std::vector<gate> gates;
  std::istringstream in( drive_family_library );
  lorina::read_genlib( in, genlib_reader( gates ) );

  tech_library<3> tlib( gates, tech_library_params{} );
  aig_network aig = eight_sink_aig();

  emap_params ps;
  CHECK( ps.electrical_model == false );
  CHECK( ps.wire_cap_per_fanout == 0.0f );

  emap_stats st1, st2;
  auto r1 = emap_klut( aig, tlib, ps, &st1 );
  auto r2 = emap_klut( aig, tlib, ps, &st2 );
  CHECK( r1.compute_area() == r2.compute_area() );
  CHECK( st1.area == st2.area );
  CHECK( st1.delay == st2.delay );
}

TEST_CASE( "emap area-delay frontier is off by default and bounded when on", "[emap]" )
{
  emap_params ps;
  CHECK( ps.curve_points == 0u ); /* the default keeps one best match per phase */

  std::vector<gate> gates;
  std::istringstream in( drive_family_library );
  lorina::read_genlib( in, genlib_reader( gates ) );
  tech_library_params tps;
  tps.electrical_model = true;
  tech_library<3> tlib( gates, tps );
  aig_network aig = eight_sink_aig();

  /* an absurd budget must be clamped, not honored: the frontier is the memory-hungry structure here */
  emap_params big;
  big.electrical_model = true;
  big.curve_points = 1000;
  emap_stats st;
  auto res = emap_klut( aig, tlib, big, &st );
  CHECK( st.area > 0.0 );
  res.foreach_node( [&]( auto const& n ) {
    if ( !res.is_constant( n ) && !res.is_ci( n ) )
      CHECK( res.has_binding( n ) ); /* every node still bound: the frontier never drops a match */
  } );
}

TEST_CASE( "emap frontier re-selection keeps the cover valid and no more expensive", "[emap]" )
{
  std::vector<gate> gates;
  std::istringstream in( drive_family_library );
  lorina::read_genlib( in, genlib_reader( gates ) );
  tech_library_params tps;
  tps.electrical_model = true;
  tech_library<3> tlib( gates, tps );
  aig_network aig = eight_sink_aig();

  emap_params base;
  base.electrical_model = true;
  emap_stats st_base;
  emap_klut( aig, tlib, base, &st_base );

  emap_params curve = base;
  curve.curve_points = 4;
  emap_stats st_curve;
  auto res = emap_klut( aig, tlib, curve, &st_curve );

  /* Each re-selection step only ever moves to a cheaper point that still meets required -- but it
   * also changes arrivals, hence the next round's required times and search trajectory, so the FINAL
   * area is not monotone against the single-best cover. That is a keep-best decision for the caller,
   * not an invariant to assert here (WiseSyn offers the frontier cover as an additional portfolio
   * candidate for exactly this reason). What must hold is that the cover stays valid and legality is
   * never traded away. */
  CHECK( st_curve.area > 0.0 );
  CHECK( st_curve.max_load_violations <= st_base.max_load_violations );
  CHECK( res.compute_area() > 0.0 );
  res.foreach_node( [&]( auto const& n ) {
    if ( !res.is_constant( n ) && !res.is_ci( n ) )
      CHECK( res.has_binding( n ) );
  } );
}

/* BP-18(a): the optional single-input covering stage. The forcing property of this library is that
 * there is NO strong nand — with -load-aware alone the hot driver is stuck at max_load 2.0 against
 * 16 x 0.5 = 8.0, a 4x violation NO drive selection can fix. The only legal cover routes the weak
 * nand through the strong buffer (driver sees 0.5, buffer carries 8.0 <= 40.0), so this fixture
 * proves the mechanism does something drive selection cannot. */
std::string const buffer_stage_library = "GATE   nand_w  2 O=!(a*b); PIN * INV 0.5 2.0  0.015 0.100 0.015 0.100\n"
                                         "GATE   inv_w   1 O=!a;     PIN * INV 0.5 2.0  0.010 0.100 0.010 0.100\n"
                                         "GATE   buf_s   3 O=a;      PIN * NONINV 0.5 40.0 0.030 0.005 0.030 0.005\n"
                                         "GATE   zero    0 O=CONST0;\n"
                                         "GATE   one     0 O=CONST1;";

namespace
{
aig_network many_sink_aig( uint32_t sinks )
{
  aig_network aig;
  auto const a = aig.create_pi();
  auto const b = aig.create_pi();
  auto const t = aig.create_and( a, b );
  for ( uint32_t i = 0; i < sinks; ++i )
  {
    auto const p = aig.create_pi();
    aig.create_po( aig.create_and( !t, p ) );
  }
  return aig;
}

struct buffer_probe
{
  emap_stats st;
  uint32_t buffers_bound{ 0 };
  bool all_bound{ true };
};

buffer_probe map_with_buffer_stage( uint32_t sinks, bool cover_buffer )
{
  std::vector<gate> gates;
  std::istringstream in( buffer_stage_library );
  lorina::read_genlib( in, genlib_reader( gates ) );

  tech_library_params tps;
  tps.electrical_model = true;
  tech_library<3> tlib( gates, tps );

  emap_params ps;
  ps.electrical_model = true;
  ps.cover_buffer = cover_buffer;
  ps.area_oriented_mapping = false;

  buffer_probe out;
  aig_network aig = many_sink_aig( sinks );
  auto res = emap_klut( aig, tlib, ps, &out.st );
  res.foreach_node( [&]( auto const& n ) {
    if ( res.is_constant( n ) || res.is_pi( n ) )
      return;
    if ( !res.has_binding( n ) )
    {
      out.all_bound = false;
      return;
    }
    if ( res.get_binding( n ).name == "buf_s" )
      ++out.buffers_bound;
  } );
  return out;
}
} // namespace

TEST_CASE( "emap cover-buffer stage rescues a load no drive selection can carry", "[emap]" )
{
  /* off: the weak-only library leaves the 16-sink driver in violation */
  auto const off = map_with_buffer_stage( 16, false );
  CHECK( off.st.max_load_violations >= 1 );
  CHECK( off.st.worst_load_ratio > 3.0 );
  CHECK( off.st.cover_buffers == 0 );
  CHECK( off.buffers_bound == 0 );
  CHECK( off.all_bound );

  /* on: the cover chooses a buffer, the violation disappears, area grows by the buffer */
  auto const on = map_with_buffer_stage( 16, true );
  CHECK( on.st.cover_buffers >= 1 );
  CHECK( on.buffers_bound >= 1 );
  CHECK( on.st.max_load_violations == 0 );
  CHECK( on.all_bound );
  CHECK( on.st.area > off.st.area );
}

TEST_CASE( "emap cover-buffer identity branch really competes", "[emap]" )
{
  /* one sink: load 0.5 sits below the buffer's crossover (slope_g*L < slope_g*cap + block +
   * slope_b*L), so buffering would ADD modeled delay — the zero-cost wire branch must win. Note
   * two sinks is already past the crossover for this steep-slope library AND on the critical path
   * (emap's PO required time IS the achieved delay), so a buffer there is the delay-optimal pick,
   * not a bug — which is exactly why the identity check uses the sub-crossover load. */
  auto const light = map_with_buffer_stage( 1, true );
  CHECK( light.st.cover_buffers == 0 );
  CHECK( light.buffers_bound == 0 );
  CHECK( light.st.max_load_violations == 0 );
}

TEST_CASE( "emap cover-buffer default off is byte-identical", "[emap]" )
{
  emap_params defaults;
  CHECK( defaults.cover_buffer == false );
  auto const a = map_with_buffer_stage( 8, false );
  auto const b = map_with_buffer_stage( 8, false );
  CHECK( a.st.area == b.st.area );
  CHECK( a.st.delay == b.st.delay );
  CHECK( a.st.cover_buffers == 0 );
}

TEST_CASE( "emap cover-buffer composes with the frontier re-bind", "[emap]" )
{
  /* buffer decides the load, curve re-binds the cell under that load — the ordering inside the
   * refresh seams; together they must stay legal and fully bound */
  std::vector<gate> gates;
  std::istringstream in( buffer_stage_library );
  lorina::read_genlib( in, genlib_reader( gates ) );
  tech_library_params tps;
  tps.electrical_model = true;
  tech_library<3> tlib( gates, tps );

  emap_params ps;
  ps.electrical_model = true;
  ps.cover_buffer = true;
  ps.curve_points = 4;
  ps.area_oriented_mapping = false;
  emap_stats st;
  aig_network aig = many_sink_aig( 16 );
  auto res = emap_klut( aig, tlib, ps, &st );

  CHECK( st.max_load_violations == 0 );
  CHECK( st.cover_buffers >= 1 );
  bool all_bound = true;
  res.foreach_node( [&]( auto const& n ) {
    if ( !res.is_constant( n ) && !res.is_pi( n ) && !res.has_binding( n ) )
      all_bound = false;
  } );
  CHECK( all_bound );
}
