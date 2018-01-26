#ifdef TESTING

#  define DO_THIS_CUBE_SLICING_TEST
#  ifdef DO_THIS_CUBE_SLICING_TEST

#    include "gtest/gtest.h"

#    include "CubeSlicing.h"

#    include "MathFunctions.h"
using namespace mathfunc;

namespace {
    bool planeCubeSlicingResultEqual(const PlaneCubeSlicingResult& r1, const PlaneCubeSlicingResult& r2, double eps) {
        if(abs(r1.volumeIn - r2.volumeIn) > eps) return false;
        size_t s = r1.areaIn.size();
        for(size_t i = 0; i < r1.areaIn.size(); ++i) {
            if(abs(r1.areaIn[i] - r2.areaIn[i]) > eps) return false;
        }
        return true;
    }

    PlaneCubeSlicingResult planeUnitCubeSliceByIntersection(double x, double y, double z) {
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
}

#  endif
#endif //TESTING
