#ifndef MEDYAN_Structure_CellList_Hpp
#define MEDYAN_Structure_CellList_Hpp

#include <cstddef> // size_t
#include <vector>

namespace internal {

// User of elements, like each molecule in compartments
template< typename T >
struct CellListElementUser {
    T* user;
    std::size_t head;
    std::size_t index;
};

// The necessary structure for element list
template< typename T >
struct CellListElement {
    using UserType = CellListElementUser<T>;

    UserType* elementUser;
    std::size_t next;
    std::size_t prev;
};

// User of heads, like each compartment holding molecules
template< typename T >
struct CellListHeadUser {
    T* user;
    std::size_t head;
};

// The necessary structure for head list
template< typename T >
struct CellListHead {
    using UserType = CellListHeadUser<T>;

    UserType* headUser;
    std::size_t first;
    std::size_t last;
};

// Note:
//   - 0th position for node list is reserved for end() iterator
//   - next/prev index 0 means no next/prev index
template< typename T >
using CellListElementList = std::vector< CellListElement<T> >;
template< typename T >
using CellListHeadList    = std::vector< CellListHead<T>    >;

template< typename TElement, typename THead >
inline void cellListClear(CellListElementList<TElement>& elements, CellListHeadList<THead>& heads) {
    elements.resize(1);
    for(auto& head : heads) {
        head.first = 0;
        head.last  = 0;
    }
} // void cellListClear(...)

template< typename TElement, typename THead >
inline void cellListAddElement(
    CellListElementUser<TElement>& eu, const CellListHeadUser<THead>& hu,
    CellListElementList<TElement>& elements, CellListHeadList<THead>& heads
) {
    const auto head = hu.head;
    eu.head = head;

    // Add new element
    const auto newIndex = elements.size();
    elements.push_back({&eu});
    auto& lastAddedElement = elements.back();

    // Reconnect linked list
    const auto last = heads[head].last;
    if(last) {
        // Insert as last
        lastAddedElement.prev = last;
        lastAddedElement.next = 0; // no connection
        elements[last].next = newIndex;
        heads[head].last = newIndex;
    }
    else {
        // Was empty
        lastAddedElement.prev = 0;
        lastAddedElement.next = 0;
        heads[head].first = newIndex;
        heads[head].last  = newIndex;
    }
} // void cellListAddElement

template< typename TElement, typename THead >
inline void cellListUpdateCell(
    CellListElementUser<TElement>& eu, const CellListHeadUser<THead>& newHu,
    CellListElementList<TElement>& elements, CellListHeadList<THead>& heads
) {
    const auto newHead = newHu.head;
    const auto index = eu.index;
    const auto oldHead = eu.head;

    // Remove from current head
    {
        const auto prev = elements[index].prev;
        const auto next = elements[index].next;
        if(prev) elements[prev].next = next; else heads[oldHead].first = next;
        if(next) elements[next].prev = prev; else heads[oldHead].last  = prev;
    }

    // Add to new head last
    {
        eu.head = newHead;

        const auto last = heads[newHead].last;
        if(last) {
            // Insert as last
            elements[index].prev = last;
            elements[index].next = 0; // no connection
            elements[last].next = index;
            heads[newHead].last = index;
        }
        else {
            // Was empty
            elements[index].prev = 0;
            elements[index].next = 0;
            heads[newHead].first = index;
            heads[newHead].last  = index;
        }
    }
} // void cellListUpdateCell(...)

} // namespace internal

#endif
