#include <memory> // unique_ptr

#include "catch2/catch.hpp"

#include "Structure/Database.h"

TEST_CASE("Database tests", "[Database]") {
    SECTION("Database with dynamic indexing") {
        struct DummyData {
            std::vector< int > d;
            void push_back(int x) { d.push_back(x); }
            void pop_back() { d.pop_back(); }
            void move_from_back(std::size_t dbIndex) { d[dbIndex] = d.back(); }
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

    SECTION("Database with stable indexing") {
        struct DummyData {
            std::vector< int > d;
            void push_back(int x) { d.push_back(x); }
            void move_content(size_t from, size_t to) { d[to] = d[from]; }
            void resize(size_t size) { d.resize(size); }
        };
        struct Dummy : Database< Dummy, true, DummyData > {
            Dummy(int x) : Database< Dummy, true, DummyData >(x) {}
        };

        // Test adding and removing elements from the vectorized data

        auto dummy0 = std::make_unique<Dummy>(0);
        // d: [0]
        // deleted: []
        REQUIRE(Dummy::getElements().size() == 1);
        REQUIRE(Dummy::getDbData().d.size() == 1);
        REQUIRE(dummy0->getIndex() == 0);
        REQUIRE(Dummy::getStableElement(0) == dummy0.get());
        REQUIRE(dummy0->getStableIndex() == 0);

        auto dummy1 = std::make_unique<Dummy>(1);
        auto dummy2 = std::make_unique<Dummy>(2);
        auto dummy3 = std::make_unique<Dummy>(3);
        auto dummy4 = std::make_unique<Dummy>(4);
        dummy3.reset(nullptr);
        auto dummy5 = std::make_unique<Dummy>(5);
        auto dummy6 = std::make_unique<Dummy>(6);
        dummy5.reset(nullptr);
        dummy1.reset(nullptr);
        dummy2.reset(nullptr);
        // dummies: 0, 6, 4
        // d: [0, 1, 2, 3, 4, 5, 6]
        // deleted: [3, 5, 1, 2]
        REQUIRE(Dummy::getElements().size() == 3);
        REQUIRE(Dummy::getDbData().d.size() == 7);
        REQUIRE(dummy4->getIndex() == 2);
        REQUIRE(Dummy::getStableElement(4) == dummy4.get());
        REQUIRE(dummy4->getStableIndex() == 4);

        Dummy::rearrange();
        // dummies: 0, 6, 4
        // d: [0, 4, 6]
        // deleted: []
        REQUIRE(Dummy::getElements().size() == 3);
        REQUIRE(Dummy::getDbData().d.size() == 3);
        REQUIRE(dummy4->getIndex() == 2);
        REQUIRE(Dummy::getStableElement(1) == dummy4.get());
        REQUIRE(dummy4->getStableIndex() == 1);
        REQUIRE(dummy6->getIndex() == 1);
        REQUIRE(Dummy::getStableElement(2) == dummy6.get());
        REQUIRE(dummy6->getStableIndex() == 2);

    }
}
