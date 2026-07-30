#include <catch.hpp>

#include <array>
#include <cmath>

#include <mockturtle/utils/load_curve.hpp>

using namespace mockturtle;

TEST_CASE( "load ladder spans the library range geometrically", "[load_curve]" )
{
  auto const L = build_load_ladder( 1.0f, 1000.0f, 4 );

  REQUIRE( L.size == 4 );
  /* both anchors are ON the ladder — the smallest realizable load and the largest legal one
   * are exactly the points a match is priced at, not interpolated between */
  CHECK( L.points[0] == Approx( 1.0f ) );
  CHECK( L.points[3] == Approx( 1000.0f ) );
  /* geometric: equal ratios, so the resolution follows the decades a cell library actually
   * spans rather than crowding every point near the top */
  CHECK( L.points[1] == Approx( 10.0f ) );
  CHECK( L.points[2] == Approx( 100.0f ) );

  for ( uint32_t i = 1; i < L.size; ++i )
    CHECK( L.points[i] > L.points[i - 1] );
}

TEST_CASE( "load ladder is bounded and deterministic", "[load_curve]" )
{
  /* the hard cap is the store's memory bound (ADR-0050): a request above it is clamped, not
   * honored, and never allocates more than the declared budget formula predicts */
  auto const big = build_load_ladder( 0.5f, 500.0f, 1000 );
  CHECK( big.size == max_load_points );

  auto const off = build_load_ladder( 0.5f, 500.0f, 0 );
  CHECK( off.size == 0 );

  /* two constructions from the same anchors are bit-identical (ADR-0042) */
  auto const a = build_load_ladder( 0.7f, 123.4f, 8 );
  auto const b = build_load_ladder( 0.7f, 123.4f, 8 );
  for ( uint32_t i = 0; i < a.size; ++i )
    CHECK( a.points[i] == b.points[i] );
}

TEST_CASE( "a degenerate library range collapses to one point", "[load_curve]" )
{
  /* no usable spread (equal anchors, an inverted range, or a non-positive floor) is not an
   * error and not a guess: the ladder degenerates to the single load it can defend, which
   * makes the mechanism reduce exactly to the scalar model instead of inventing a range */
  CHECK( build_load_ladder( 4.0f, 4.0f, 8 ).size == 1 );
  CHECK( build_load_ladder( 9.0f, 2.0f, 8 ).size == 1 );
  CHECK( build_load_ladder( 0.0f, 10.0f, 8 ).size == 1 );

  auto const one = build_load_ladder( 4.0f, 4.0f, 8 );
  CHECK( one.points[0] == Approx( 4.0f ) );

  /* one requested point over a real range samples the geometric mean — the balanced single
   * sample of a range whose points are geometric */
  auto const mean = build_load_ladder( 1.0f, 100.0f, 1 );
  REQUIRE( mean.size == 1 );
  CHECK( mean.points[0] == Approx( 10.0f ) );
}

TEST_CASE( "curve evaluation is exact at the knots", "[load_curve]" )
{
  auto const L = build_load_ladder( 1.0f, 1000.0f, 4 );
  std::array<double, max_load_points> v{ 10.0, 20.0, 40.0, 80.0 };

  for ( uint32_t i = 0; i < L.size; ++i )
    CHECK( interpolate_curve( L, v, L.points[i] ) == Approx( v[i] ) );
}

TEST_CASE( "curve evaluation interpolates linearly between knots", "[load_curve]" )
{
  load_ladder L{};
  L.size = 3;
  L.points[0] = 0.0f;
  L.points[1] = 10.0f;
  L.points[2] = 20.0f;
  std::array<double, max_load_points> v{ 100.0, 200.0, 400.0 };

  CHECK( interpolate_curve( L, v, 5.0f ) == Approx( 150.0 ) );
  CHECK( interpolate_curve( L, v, 15.0f ) == Approx( 300.0 ) );
}

TEST_CASE( "curve evaluation clamps below and extrapolates above", "[load_curve]" )
{
  load_ladder L{};
  L.size = 2;
  L.points[0] = 10.0f;
  L.points[1] = 20.0f;
  std::array<double, max_load_points> v{ 100.0, 200.0 };

  /* below the floor: clamp. The floor is the smallest load any real sink can present, so a
   * query under it is a rounding artifact, and downward extrapolation would invent a delay
   * smaller than any the library can deliver. */
  CHECK( interpolate_curve( L, v, 1.0f ) == Approx( 100.0 ) );
  CHECK( interpolate_curve( L, v, 10.0f ) == Approx( 100.0 ) );

  /* above the ceiling: extrapolate on the last segment. Clamping here would price an
   * over-loaded driver as if load stopped mattering — optimistic exactly where the cover must
   * be pessimistic. Delay really is linear in load, so the last segment is the honest slope. */
  CHECK( interpolate_curve( L, v, 30.0f ) == Approx( 300.0 ) );
  CHECK( interpolate_curve( L, v, 25.0f ) == Approx( 250.0 ) );

  /* a one-point ladder has no segment to extrapolate on and clamps in both directions */
  load_ladder one{};
  one.size = 1;
  one.points[0] = 10.0f;
  std::array<double, max_load_points> w{ 42.0 };
  CHECK( interpolate_curve( one, w, 1.0f ) == Approx( 42.0 ) );
  CHECK( interpolate_curve( one, w, 1000.0f ) == Approx( 42.0 ) );
}

TEST_CASE( "an empty ladder evaluates to the scalar fallback", "[load_curve]" )
{
  load_ladder const off{};
  std::array<double, max_load_points> v{};
  /* size 0 is the feature being off; the caller's scalar arrival is returned unchanged so a
   * disabled mechanism cannot perturb a single number */
  CHECK( interpolate_curve( off, v, 5.0f, 7.5 ) == Approx( 7.5 ) );
}

TEST_CASE( "one match per ladder point, chosen on a stable key", "[load_curve]" )
{
  /* an empty point takes any candidate */
  CHECK( load_curve_prefer( 5.0, 2.0f, 3, false, 0.0, 0.0f, 0 ) );

  /* at a fixed load the earlier signal wins, whatever it costs */
  CHECK( load_curve_prefer( 4.0, 100.0f, 4, true, 5.0, 1.0f, 1 ) );
  CHECK( !load_curve_prefer( 6.0, 1.0f, 1, true, 5.0, 100.0f, 4 ) );

  /* equal arrival falls to area, then to cut size */
  CHECK( load_curve_prefer( 5.0, 1.0f, 3, true, 5.0, 2.0f, 3 ) );
  CHECK( !load_curve_prefer( 5.0, 2.0f, 3, true, 5.0, 1.0f, 3 ) );
  CHECK( load_curve_prefer( 5.0, 1.0f, 2, true, 5.0, 1.0f, 3 ) );

  /* a full tie keeps the incumbent, so the result depends only on enumeration order */
  CHECK( !load_curve_prefer( 5.0, 1.0f, 3, true, 5.0, 1.0f, 3 ) );

  /* arrivals within the epsilon are a tie, not a win: two drives of one function differ in the
   * last bits of a linearized delay, and that must not outrank a real area difference */
  CHECK( !load_curve_prefer( 5.0 - 1e-12, 2.0f, 3, true, 5.0, 1.0f, 3 ) );
}

TEST_CASE( "store size follows the declared budget formula", "[load_curve]" )
{
  /* ADR-0050 clause 6: the size is computed BEFORE allocating, so a script learns it cannot
   * have the mechanism at this design size instead of discovering it at the OOM killer */
  CHECK( load_curve_store_bytes( 1000, 4, 32 ) == 1000ull * 2ull * 4ull * 32ull );
  CHECK( load_curve_store_bytes( 0, 8, 32 ) == 0ull );
  CHECK( load_curve_store_bytes( 1000, 0, 32 ) == 0ull );

  /* no overflow on a design far larger than anything that fits in memory */
  CHECK( load_curve_store_bytes( 4000000000ull, 16, 64 ) == 4000000000ull * 2ull * 16ull * 64ull );
}
