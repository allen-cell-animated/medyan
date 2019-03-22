#include <memory> // unique_ptr

#include "catch2/catch.hpp"

#include "Structure/Database.h"

TEST_CASE("Database tests", "[Database]") {
    struct DummyData {
        std::vector< int > d;
        void push_back(int x) { d.push_back(x); }
        void pop_back() { d.pop_back(); }
        void copy_from_back(std::size_t dbIndex) { d[dbIndex] = d.back(); }
    };
    struct Dummy : Database< Dummy, false, DummyData > {
        Dummy(int x) : Database< Dummy, false, DummyData >(x) {}
    };

    // Test adding and removing elements from the vectorized data

    auto dummy0 = std::make_unique<Dummy>(0);
    // d: [0]
    REQUIRE(Dummy::getElements().size() == 1);
    REQUIRE(Dummy::getDbData().d.size() == 1);
    REQUIRE(dummy0->getIndex() == 0);

    auto dummy1 = std::make_unique<Dummy>(1);
    auto dummy2 = std::make_unique<Dummy>(2);
    // d: [0, 1, 2]
    REQUIRE(Dummy::getElements().size() == 3);
    REQUIRE(Dummy::getDbData().d.size() == 3);
    REQUIRE(dummy2->getIndex() == 2);
    REQUIRE(Dummy::getDbData().d[2] == 2);

    dummy1.reset(nullptr);
    // d: [0, 2]
    REQUIRE(Dummy::getElements().size() == 2);
    REQUIRE(Dummy::getDbData().d.size() == 2);
    REQUIRE(dummy2->getIndex() == 1);
    REQUIRE(Dummy::getDbData().d[1] == 2);

    auto dummy3 = std::make_unique<Dummy>(3);
    // d: [0, 2, 3]
    REQUIRE(Dummy::getElements().size() == 3);
    REQUIRE(Dummy::getDbData().d.size() == 3);
    REQUIRE(dummy3->getIndex() == 2);
    REQUIRE(Dummy::getDbData().d[2] == 3);

    dummy3.reset(nullptr);
    // d: [0, 2]
    REQUIRE(Dummy::getElements().size() == 2);
    REQUIRE(Dummy::getDbData().d.size() == 2);
}
