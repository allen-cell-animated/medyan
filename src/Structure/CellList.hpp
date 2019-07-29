#ifndef MEDYAN_Structure_CellList_Hpp
#define MEDYAN_Structure_CellList_Hpp

#include <cstddef> // size_t
#include <vector>

namespace cell_list {

// Forward declarations
template< typename TElement, typename THead > class CellListManager;

// Storing redundant information for user of elements,
// like each molecule
template< typename TElement, typename THead >
struct CellListElementUser {
    using ManagerType = CellListManager< TElement, THead >;

    ManagerType* manager = nullptr;
    std::size_t head;
    std::size_t index;
};

// The necessary structure for element list
template< typename T >
struct CellListElement {
    T* ptr;
    std::size_t next;
    std::size_t prev;
    bool hasNext;
    bool hasPrev;
};

// Storing redundant information for user of heads,
// like each compartment storing molecules
template< typename TElement, typename THead >
struct CellListHeadUser {
    using ManagerType = CellListManager< TElement, THead >;

    ManagerType* manager = nullptr;
    std::size_t index;
};

// The necessary structure for head list
template< typename T >
struct CellListHead {
    T* ptr;
    std::size_t size = 0;
    std::size_t first;
    std::size_t last;
};

// The class managing a cell list
template< typename TElement, typename THead >
class CellListManager {
public:
    using ElementList = std::vector< CellListElement<TElement> >;
    using HeadList    = std::vector< CellListHead   <THead   > >;
    using ElementUser = CellListElementUser< TElement, THead >;
    using HeadUser    = CellListHeadUser   < TElement, THead >;

    // Accessors
    //-------------------------------------------------------------------------

    auto getHeadPtr(std::size_t headIndex) const { return headList_[headIndex].ptr; }
    auto getHeadPtr(const ElementUser& eu) const { return getHeadPtr(eu.head); }

    // Element operations with fixed head users
    //-------------------------------------------------------------------------

    void clearElements() {
        elementList_.clear();
        elementDeletedIndices_.clear();
        for(auto& head : headList_) {
            head.size  = 0;
        }
    }

    void addElement(
        TElement* element,
        ElementUser& eu, const HeadUser& headUser
    ) {
        const auto head = headUser.head;

        // Add the new element
        std::size_t newIndex;
        if(elementDeletedIndices_.empty()) {
            elementList_.push_back({element});
            newIndex = elementList_.size() - 1;
        } else {
            newIndex = elementDeletedIndices_.back();
            elementList_[newIndex].ptr = element;
            elementDeletedIndices_.pop_back();
        }

        eu.head = head;
        eu.index = newIndex;

        // Connect the new element
        registerElement_(newIndex, head);

    } // void addElement(...)

    void updateElement(ElementUser& eu, const HeadUser& newHeadUser) {
        const auto newHead = newHeadUser.head;
        const auto index = eu.index;
        const auto oldHead = eu.head;

        // Remove from current head
        unregsiterElement_(index, oldHead);

        // Add to the new head
        eu.head = newHead;
        registerElement_(index, newHead);

    } // void updateElement(...)

    void removeElement(const ElementUser& eu) {
        const auto index = eu.index;
        const auto head = eu.head;

        unregsiterElement_(index, head);

        elementDeletedIndices_.push_back(index);
    } // void removeElement(...)

    // Head operations with fixed element users
    //-------------------------------------------------------------------------

    void addHead(THead* head, HeadUser& hu) {
        headList_.push_back({head});
        hu.index = headList_.size() - 1;
    }

private:
    HeadList                   headList_;
    ElementList                elementList_;
    std::vector< std::size_t > elementDeletedIndices_;

    // Helper function to register an element at an index to a head
    // Note:
    //   - This function causes the head size count to increase by 1.
    //   - This function does not manage allocating the element.
    void registerElement_(const std::size_t index, const std::size_t head) {
        if(headList_[head].size) {
            // Insert as last
            const auto last = headList_[head].last;
            elementList_[index].hasPrev = true;
            elementList_[index].prev = last;
            elementList_[index].hasNext = false;

            elementList_[last].hasNext = true;
            elementList_[last].next = index;
            headList_[head].last = index;
        }
        else {
            // Was empty
            elementList_[index].hasPrev = false;
            elementList_[index].hasNext = false;
            headList_[head].first = index;
            headList_[head].last  = index;
        }

        ++headList_[head].size;
    } // void registerElement_(...)

    // Helper function to unregister an element at an index from a head
    // Note:
    //   - This function causes the head size count to decrease by 1.
    //   - This function does not manage removing the element.
    void unregsiterElement_(const std::size_t index, const std::size_t head) {
        --headList_[head].size;

        // Reconnect linked list
        const auto& e = elementList_[index];

        if(e.hasPrev) {
            elementList_[e.prev].hasNext = e.hasNext;
            if(e.hasNext) elementList_[e.prev].next = e.next;
        } else {
            if(e.hasNext) headList_[head].first = e.next;
        }

        if(e.hasNext) {
            elementList_[e.next].hasPrev = e.hasPrev;
            if(e.hasPrev) elementList_[e.next].prev = e.prev;
        } else {
            if(e.hasPrev) headList_[head].last = e.prev;
        }
    } // void unregisterElement_(...)

}; // template<...> class CellListManager



} // namespace cell_list

#endif
