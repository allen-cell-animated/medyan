#include <cstddef> // size_t

#include "catch2/catch.hpp"

#include "util/math/vec.hpp"

TEST_CASE("VecArray tests", "[VecArray]") {

    using std::size_t;
    using namespace mathfunc;

    VecArray< 3, float > v3f;
    VecArray< 4, double > v4d;

    // Appending elements
    const Vec< 3, float > push_v3f[] {
        0.0f, 0.0f, 0.0f,
        1.0f, 1.0f, 1.0f,
        2.0f, 2.0f, 2.0f,
        3.0f, 3.0f, 3.0f,
        4.0f, 4.0f, 4.0f
    };
    const Vec< 4, double > push_v4d[] {
        -4.0, -4.0, -4.0, -4.0,
        -3.0, -3.0, -3.0, -3.0,
        -2.0, -2.0, -2.0, -2.0,
        -1.0, -1.0, -1.0, -1.0
    };

    for(auto& x : push_v3f) v3f.push_back(x);
    for(auto& x : push_v4d) v4d.push_back(x);

    REQUIRE(v3f.size() == 5);
    REQUIRE(v3f.size_raw() == 15);
    REQUIRE(v4d.size() == 4);
    REQUIRE(v4d.size_raw() == 16);

    SECTION("RefVec accessing and iterating") {
        // Accessing
        auto v3f_2 = v3f[2];
        CHECK(v3f_2[0] == Approx(2.0f));
        CHECK(v3f_2[1] == Approx(2.0f));
        CHECK(v3f_2[2] == Approx(2.0f));

        Vec< 3, float > v3f_2_new { 2.1f, 2.1f, 2.1f };
        v3f[2] = v3f_2_new;
        CHECK(v3f_2[0] == Approx(2.1f));
        CHECK(v3f_2[1] == Approx(2.1f));
        CHECK(v3f_2[2] == Approx(2.1f));

        const auto& v4d_cref = v4d;
        auto v4d_1c = v4d_cref[1];
        CHECK(v4d_1c[0] == Approx(-3.0));
        CHECK(v4d_1c[1] == Approx(-3.0));
        CHECK(v4d_1c[2] == Approx(-3.0));
        CHECK(v4d_1c[3] == Approx(-3.0));

        // Iterating
        auto v4d_2 = v4d[2];
        size_t v4d_2_cnt = 0;
        double v4d_2_sum = 0.0;
        for(auto x : v4d_2) {
            ++v4d_2_cnt;
            v4d_2_sum += x;
        }
        REQUIRE(v4d_2_cnt == 4);
        CHECK(v4d_2_sum == Approx(-8.0));

        auto v4d_3 = v4d[3];
        double v4d_3_accu = 0.0;
        for(auto& x : v4d_3) {
            v4d_3_accu += 0.1;
            x += v4d_3_accu;
        }
        CHECK(v4d_3[0] == Approx(-0.9));
        CHECK(v4d_3[1] == Approx(-0.8));
        CHECK(v4d_3[2] == Approx(-0.7));
        CHECK(v4d_3[3] == Approx(-0.6));
    }

    SECTION("Iterator usage") {
        // Iterator pointer, dereference, op[]
        auto v3f_i = v3f.begin() + 1;
        CHECK(v3f_i->size() == 3);
        CHECK((*v3f_i)[0] == Approx(1.0f));
        CHECK(v3f_i[2][0] == Approx(3.0f));

        // Ctor and assignment
        const auto& v3f_cref = v3f;
        auto v3f_i_ccopy = v3f_cref.begin();
        CHECK((*v3f_i_ccopy)[0] == Approx(0.0f));
        v3f_i_ccopy = v3f_i; // Copy assign
        CHECK((*v3f_i_ccopy)[0] == Approx(1.0f));
        decltype(v3f)::const_iterator v3f_ic(v3f_i); // Copy ctor
        CHECK((*v3f_ic)[0] == Approx(1.0f));

        // Check comparison
        auto v3f_ie = v3f.end();
        auto v3f_ice = v3f_cref.end();
        CHECK(v3f_ie == v3f_ice);
        CHECK(v3f_i != v3f_ie);
        CHECK(v3f_i != v3f_ice);
        CHECK(!(v3f_i > v3f_ie));
        CHECK(  v3f_i < v3f_ie );
        CHECK(!(v3f_i >= v3f_ie));
        CHECK(  v3f_i <= v3f_ie );
        CHECK(v3f_ie >= v3f_ice);
        CHECK(v3f_ie <= v3f_ice);

        // Check arithmetics
        CHECK(v3f_i - v3f_ice == -4);
        CHECK(v3f_i + 4 == v3f_ice);
        CHECK(4 + v3f_i == v3f_ice);
        CHECK(v3f_ice - 4 == v3f_i);

        // Check inc/dec
        ++v3f_i;
        CHECK((*v3f_i)[0] == Approx(2.0f));
        --v3f_i;
        CHECK((*v3f_i)[0] == Approx(1.0f));
        v3f_i += 2;
        CHECK((*v3f_i)[0] == Approx(3.0f));
        v3f_i -= 2;
        CHECK((*v3f_i)[0] == Approx(1.0f));
        auto v3f_i_copy = v3f_i++;
        CHECK((*v3f_i_copy)[0] == Approx(1.0f));
        CHECK((*v3f_i     )[0] == Approx(2.0f));
        v3f_i_copy = v3f_i--;
        CHECK((*v3f_i_copy)[0] == Approx(2.0f));
        CHECK((*v3f_i     )[0] == Approx(1.0f));
    }

    SECTION("Element popping") {
        v3f.pop_back();
        CHECK(v3f.size() == 4);
        CHECK(v3f.size_raw() == 12);
    }
}
