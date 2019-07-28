#ifndef MEDYAN_Structure_CellList_Hpp
#define MEDYAN_Structure_CellList_Hpp

#include <cstddef> // size_t
#include <vector>

namespace internal {

template< typename T >
struct CellListElement {
    T* user;
    std::size_t head;
    std::size_t index;
};

template< typename T >
struct CellListNode {
    using ElementType = CellListElement<T>;

    ElementType* element;
    std::size_t next;
    std::size_t prev;
};

struct CellListHead {
    std::size_t first;
    std::size_t last;
};

// Note:
//   - 0th position for node list is reserved for end() iterator
//   - next/prev index 0 means no next/prev index
template< typename T >
using CellListNodeList = std::vector< CellListNode<T> >;
using CellListHeadList = std::vector< CellListHead >;

template< typename T >
inline void cellListClear(CellListNodeList<T>& nodes, CellListHeadList& heads) {
    nodes.resize(1);
    for(auto& head : heads) {
        head.first = 0;
        head.last  = 0;
    }
} // void cellListClear(...)

template< typename T >
inline void cellListAddElement(CellListElement<T>& element, const CellListHead& head, CellListNodeList<T>& nodes, CellListHeadList& heads) {
    element.head = head;

    // Add new element
    const auto newIndex = nodes.size();
    nodes.push_back({&element});
    auto& lastAddedNode = nodes.back();

    // Reconnect linked list
    const auto last = heads[head].last;
    if(last) {
        // Insert as last
        lastAddedNode.prev = last;
        lastAddedNode.next = 0; // no connection
        nodes[last].next = newIndex;
        heads[head].last = newIndex;
    }
    else {
        // Was empty
        lastAddedNode.prev = 0;
        lastAddedNode.next = 0;
        heads[head].first = newIndex;
        heads[head].last  = newIndex;
    }
} // void cellListAddElement

template< typename T >
inline void cellListUpdateCell(CellListElement<T>& element, const CellListHead& newHead, CellListNodeList<T>& nodes, CellListHeadList& heads) {
    const auto index = element.index;
    const auto oldHead = element.head;

    // Remove from current head
    {
        const auto prev = nodes[index].prev;
        const auto next = nodes[index].next;
        if(prev) nodes[prev].next = next; else heads[oldHead].first = next;
        if(next) nodes[next].prev = prev; else heads[oldHead].last  = prev;
    }

    // Add to new head last
    {
        element.head = newHead;

        const auto last = heads[newHead].last;
        if(last) {
            // Insert as last
            nodes[index].prev = last;
            nodes[index].next = 0; // no connection
            nodes[last].next = index;
            heads[newHead].last = index;
        }
        else {
            // Was empty
            nodes[index].prev = 0;
            nodes[index].next = 0;
            heads[newHead].first = index;
            heads[newHead].last  = index;
        }
    }
} // void cellListUpdateCell(...)

} // namespace internal

#endif
