#include <memory> // unique_ptr
#include <vector>

#include "catch2/catch.hpp"

#include "Structure/CellList.hpp"

TEST_CASE("Cell list tests", "[CellList]") {
    using namespace cell_list;

    // Types
    //-------------------------------------------------------------------------
    struct DummyElement;
    struct DummyCell;
    using Manager = CellListManager< DummyElement, DummyCell >;

    struct DummyElement {
        CellListElementUser< DummyElement, DummyCell > ele;
    };
    struct DummyCell {
        CellListHeadUser< DummyElement, DummyCell > cell;
    };

    // Test objects
    //-------------------------------------------------------------------------
    Manager m;
    std::vector< std::unique_ptr< DummyElement > > es;
    std::vector< std::unique_ptr< DummyCell    > > cs;

    const auto registerDummyElement = [&](DummyElement& e) { e.ele .manager = &m; };
    const auto registerDummyCell    = [&](DummyCell   & c) { c.cell.manager = &m; };

    const auto buildCells = [&](std::size_t num) {
        for(std::size_t i = 0; i < num; ++i) {
            cs.push_back(std::make_unique< DummyCell >());
            m.addHead(cs.back().get(), cs.back()->cell);
            registerDummyCell(*cs.back());
        }
    };

    SECTION("Cell list operation and view") {
        DummyCell dc[3];
        for(auto& c : dc) { m.addHead(&c, c.cell); registerDummyCell(c); }

        DummyElement de[6];
        for(auto& e : de) registerDummyElement(e);

        // Cell | 0   | 1   | 2   |
        //------|-----|-----|-----|
        // Ele  | 0 4 |     | 1   |
        m.addElement(de + 0, de[0].ele, dc[0].cell);
        m.addElement(de + 1, de[1].ele, dc[2].cell);
        m.addElement(de + 4, de[4].ele, dc[0].cell);

        // Check cell 0 content
        {
            const auto c0View = m.getElementPtrs(dc[0].cell);
            std::vector< DummyElement* > deCell0(c0View.begin(), c0View.end());
            REQUIRE(deCell0.size() == 2);
            REQUIRE(deCell0[0] == de + 0);
            REQUIRE(deCell0[1] == de + 4);
        }
    }
}
