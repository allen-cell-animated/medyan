#ifdef TESTING

#  define DO_THIS_CUBE_SLICING_TEST
#  ifdef DO_THIS_CUBE_SLICING_TEST

#    include "gtest/gtest.h"

#    include "CuboidSlicing.h"

#    include "MathFunctions.h"
using namespace mathfunc;

namespace {
    bool planeCubeSlicingResultEqual(const PlaneCuboidSlicingResult& r1, const PlaneCuboidSlicingResult& r2, double eps) {
        if(abs(r1.volumeIn - r2.volumeIn) > eps) return false;
        size_t s = r1.areaIn.size();
        for(size_t i = 0; i < r1.areaIn.size(); ++i) {
            if(abs(r1.areaIn[i] - r2.areaIn[i]) > eps) return false;
        }
        return true;
    }

    PlaneCuboidSlicingResult planeUnitCubeSliceByIntersection(double x, double y, double z) {
        auto normal = normalizedVector(std::array<double, 3>{{1.0/x, 1.0/y, 1.0/z}});
        auto point = std::array<double, 3>{{x, 0, 0}};
        return planeUnitCubeSlice(point, normal);
    }
}

TEST(CubeSlicingTest, UnitCubeTransitionContinuity) {
    /**************************************************************************
    Test continuity under different transitions under the simplified
    conditions, i.e. normal has all non-neg components (no normal flip is
    required, however, the reverse might secretly use a whole flip).
    
    Cube slicing result is obtained by calculating intersections on the axes,
    and the formula are different depending on the positions of the
    intersections.
    **************************************************************************/

    double inEps = 1e-6;
    double outEps = 1e-5;

    // From 3 points inside, to 2 points inside
    {
        auto r1 = planeUnitCubeSliceByIntersection(1 - inEps, 0.5, 0.5);
        auto r2 = planeUnitCubeSliceByIntersection(1 + inEps, 0.5, 0.5);
        EXPECT_TRUE(planeCubeSlicingResultEqual(r1, r2, outEps));
    }
    {
        auto r1 = planeUnitCubeSliceByIntersection(0.5, 1 - inEps, 0.5);
        auto r2 = planeUnitCubeSliceByIntersection(0.5, 1 + inEps, 0.5);
        EXPECT_TRUE(planeCubeSlicingResultEqual(r1, r2, outEps));
    }
    {
        auto r1 = planeUnitCubeSliceByIntersection(0.5, 0.5, 1 - inEps);
        auto r2 = planeUnitCubeSliceByIntersection(0.5, 0.5, 1 + inEps);
        EXPECT_TRUE(planeCubeSlicingResultEqual(r1, r2, outEps));
    }

    // From 3 points inside, to 1 point inside
    {
        auto r1 = planeUnitCubeSliceByIntersection(0.5, 1 - inEps, 1 - inEps);
        auto r2 = planeUnitCubeSliceByIntersection(0.5, 1 + inEps, 1 + inEps);
        EXPECT_TRUE(planeCubeSlicingResultEqual(r1, r2, outEps));
    }
    {
        auto r1 = planeUnitCubeSliceByIntersection(1 - inEps, 0.5, 1 - inEps);
        auto r2 = planeUnitCubeSliceByIntersection(1 + inEps, 0.5, 1 + inEps);
        EXPECT_TRUE(planeCubeSlicingResultEqual(r1, r2, outEps));
    }
    {
        auto r1 = planeUnitCubeSliceByIntersection(1 - inEps, 1 - inEps, 0.5);
        auto r2 = planeUnitCubeSliceByIntersection(1 + inEps, 1 + inEps, 0.5);
        EXPECT_TRUE(planeCubeSlicingResultEqual(r1, r2, outEps));
    }

    // Transitions between the scenario with 1 point inside
    {
        auto r1 = planeUnitCubeSliceByIntersection(0.5, 2 - inEps, 2 - inEps);
        auto r2 = planeUnitCubeSliceByIntersection(0.5, 2 + inEps, 2 + inEps);
        EXPECT_TRUE(planeCubeSlicingResultEqual(r1, r2, outEps));
    }
    {
        auto r1 = planeUnitCubeSliceByIntersection(2 - inEps, 0.5, 2 - inEps);
        auto r2 = planeUnitCubeSliceByIntersection(2 + inEps, 0.5, 2 + inEps);
        EXPECT_TRUE(planeCubeSlicingResultEqual(r1, r2, outEps));
    }
    {
        auto r1 = planeUnitCubeSliceByIntersection(2 - inEps, 2 - inEps, 0.5);
        auto r2 = planeUnitCubeSliceByIntersection(2 + inEps, 2 + inEps, 0.5);
        EXPECT_TRUE(planeCubeSlicingResultEqual(r1, r2, outEps));
    }

    // From 3 points inside, to no point inside
    {
        auto r1 = planeUnitCubeSliceByIntersection(1 - inEps, 1 - inEps, 1 - inEps);
        auto r2 = planeUnitCubeSliceByIntersection(1 + inEps, 1 + inEps, 1 + inEps);
        EXPECT_TRUE(planeCubeSlicingResultEqual(r1, r2, outEps));
    }

    // Check reverse condition
    {
        auto r1 = planeUnitCubeSliceByIntersection(1.5 - inEps, 1.5 - inEps, 1.5 - inEps);
        auto r2 = planeUnitCubeSliceByIntersection(1.5 + inEps, 1.5 + inEps, 1.5 + inEps);
        EXPECT_TRUE(planeCubeSlicingResultEqual(r1, r2, outEps));
        EXPECT_FALSE(planeCubeSlicingResultEqual(r1, r2, inEps * 1e-5)); // Check actually changed
    }
    
}

TEST(CubeSlicingTest, UnitCubeFlipping) {
    /**************************************************************************
    Test the flipping of normal vector
    **************************************************************************/

    const double nVal = 1.0 / sqrt(3);
    const double pVal = 1.0 / 6.0;

    auto r1 = planeUnitCubeSlice({{pVal, pVal, pVal}}, {{nVal, nVal, nVal}});

    // Flip x
    {
        std::array<size_t, 6> comp {{1, 0, 2, 3, 4, 5}};
        auto r2 = planeUnitCubeSlice({{1-pVal, pVal, pVal}}, {{-nVal, nVal, nVal}});
        EXPECT_DOUBLE_EQ(r1.volumeIn, r2.volumeIn);
        for(size_t i = 0; i < 6; ++i) {
            EXPECT_DOUBLE_EQ(r1.areaIn[i], r2.areaIn[comp[i]]);
        }
    }

    // Flip y
    {
        std::array<size_t, 6> comp {{0, 1, 3, 2, 4, 5}};
        auto r2 = planeUnitCubeSlice({{pVal, 1-pVal, pVal}}, {{nVal, -nVal, nVal}});
        EXPECT_DOUBLE_EQ(r1.volumeIn, r2.volumeIn);
        for(size_t i = 0; i < 6; ++i) {
            EXPECT_DOUBLE_EQ(r1.areaIn[i], r2.areaIn[comp[i]]);
        }
    }

    // Flip z
    {
        std::array<size_t, 6> comp {{0, 1, 2, 3, 5, 4}};
        auto r2 = planeUnitCubeSlice({{pVal, pVal, 1-pVal}}, {{nVal, nVal, -nVal}});
        EXPECT_DOUBLE_EQ(r1.volumeIn, r2.volumeIn);
        for(size_t i = 0; i < 6; ++i) {
            EXPECT_DOUBLE_EQ(r1.areaIn[i], r2.areaIn[comp[i]]);
        }
    }

    // Flip yz
    {
        std::array<size_t, 6> comp {{0, 1, 3, 2, 5, 4}};
        auto r2 = planeUnitCubeSlice({{pVal, 1-pVal, 1-pVal}}, {{nVal, -nVal, -nVal}});
        EXPECT_DOUBLE_EQ(r1.volumeIn, r2.volumeIn);
        for(size_t i = 0; i < 6; ++i) {
            EXPECT_DOUBLE_EQ(r1.areaIn[i], r2.areaIn[comp[i]]);
        }
    }

    // Flip zx
    {
        std::array<size_t, 6> comp {{1, 0, 2, 3, 5, 4}};
        auto r2 = planeUnitCubeSlice({{1-pVal, pVal, 1-pVal}}, {{-nVal, nVal, -nVal}});
        EXPECT_DOUBLE_EQ(r1.volumeIn, r2.volumeIn);
        for(size_t i = 0; i < 6; ++i) {
            EXPECT_DOUBLE_EQ(r1.areaIn[i], r2.areaIn[comp[i]]);
        }
    }

    // Flip xy
    {
        std::array<size_t, 6> comp {{1, 0, 3, 2, 4, 5}};
        auto r2 = planeUnitCubeSlice({{1-pVal, 1-pVal, pVal}}, {{-nVal, -nVal, nVal}});
        EXPECT_DOUBLE_EQ(r1.volumeIn, r2.volumeIn);
        for(size_t i = 0; i < 6; ++i) {
            EXPECT_DOUBLE_EQ(r1.areaIn[i], r2.areaIn[comp[i]]);
        }
    }

    // Flip xyz
    {
        std::array<size_t, 6> comp {{1, 0, 3, 2, 5, 4}};
        auto r2 = planeUnitCubeSlice({{1-pVal, 1-pVal, 1-pVal}}, {{-nVal, -nVal, -nVal}});
        EXPECT_DOUBLE_EQ(r1.volumeIn, r2.volumeIn);
        for(size_t i = 0; i < 6; ++i) {
            EXPECT_DOUBLE_EQ(r1.areaIn[i], r2.areaIn[comp[i]]);
        }
    }

}

TEST(CubeSlicingTest, TranslationAndScaling) {
    /**************************************************************************
    Test the slicing some cube with certain position and size
    **************************************************************************/

    const double nVal = 1.0 / sqrt(3);
    const double pVal = 1.0 / 6.0;
    const double r0Val = 10.0;
    const double s = 2.0;

    auto r = PlaneCubeSlicer() (
        {{r0Val + s * pVal, r0Val + s * pVal, r0Val + s * pVal}},
        {{-nVal, -nVal, -nVal}},
        {{r0Val, r0Val, r0Val}},
        s
    );

    const double exVolumeIn = 47.0 / 6.0;
    const std::array<double, 6> exAreaIn {{3.5, 4.0, 3.5, 4.0, 3.5, 4.0}};

    EXPECT_NEAR(r.volumeIn, exVolumeIn, abs(exVolumeIn) * 1e-5);
    for(size_t i = 0; i < 6; ++i) {
        EXPECT_NEAR(r.areaIn[i], exAreaIn[i], abs(exAreaIn[i]) * 1e-5);
    }

}

#  endif
#endif //TESTING
