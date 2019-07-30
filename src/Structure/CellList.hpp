#ifndef MEDYAN_Structure_CellList_Hpp
#define MEDYAN_Structure_CellList_Hpp

#include <cstddef> // ptrdiff_t, size_t
#include <iterator> // tags
#include <vector>

namespace cell_list {

// Forward declarations
template< typename TElement, typename THead > class CellListManager;

// Storing information for user of elements,
// like each molecule.
// Note:
//   - The users must make sure that the manager points to a valid location
template< typename TElement, typename THead >
struct CellListElementUser {
    using ManagerType = CellListManager< TElement, THead >;

    ManagerType* manager = nullptr;
    std::size_t head;
    std::size_t index;
};

// Storing information for user of heads,
// like each compartment storing molecules.
// Note:
//   - The users must make sure that the manager points to a valid location
template< typename TElement, typename THead >
struct CellListHeadUser {
    using ManagerType = CellListManager< TElement, THead >;

    ManagerType* manager = nullptr;
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

    // The class for viewing and iterating the element pointers in a specific cell.
    // The class acts like a double linked list.
    // Only const iterator is offered.
    class CellView {
    private:

        // The const iterator for traversing cell list
        class ConstIterator_ {
        public:
            using iterator_category = std::bidirectional_iterator_tag;
            using value_type        = TElement*;
            using difference_type   = std::ptrdiff_t;
            using pointer           = const value_type*;
            using reference         = const value_type&;

        private:
            const ElementList* el_;
            const CellListHead<THead>* head_;
            std::size_t ei_;
            bool atEnd_;

        public:
            ConstIterator_() = default;
            ConstIterator_(const ElementList* el, const CellListHead<THead>* head, std::size_t elementIndex, bool atEnd) :
                el_(el), head_(head), ei_(elementIndex), atEnd_(atEnd || (head_->size == 0))
            { }
            ConstIterator_(const ConstIterator_&) = default;
            ConstIterator_& operator=(const ConstIterator_&) = default;

            // Dereferencing
            //---------------------------------------------
            reference operator*() const { return (*el_)[ei_].ptr; }
            pointer operator->() const { return &(*el_)[ei_].ptr; }

            // Modification
            //---------------------------------------------
            // Precondition: *this is dereferenceable
            ConstIterator_& operator++() {
                if((*el_)[ei_].hasNext) ei_ = (*el_)[ei_].next;
                else atEnd_ = true;
                return *this;
            }
            // Precondition: *this is decrementable
            ConstIterator_& operator--() {
                if(atEnd_) {
                    ei_ = head_->last; // head_->size == 0 is UB
                    atEnd_ = false;
                }
                else ei_ = (*el_)[ei_].prev; // (*el_)[ei_].hasPrev == false is UB
                return *this;
            }
            ConstIterator_ operator++(int) { ConstIterator_ tmp(*this); ++(*this); return tmp; }
            ConstIterator_ operator--(int) { ConstIterator_ tmp(*this); --(*this); return tmp; }

            // Comparison
            //---------------------------------------------
            bool operator==(const ConstIterator_& rhs) const {
                return (atEnd_ && rhs.atEnd_) || (!atEnd_ && !rhs.atEnd_ && ei_ == rhs.ei_);
            }
            bool operator!=(const ConstIterator_& rhs) const {
                return !(*this == rhs);
            }
        }; // class ConstIterator_

    public:
        using const_iterator = ConstIterator_;
        using size_type = std::size_t;

        CellView(const ElementList* el, const CellListHead<THead>* head) :
            el_(el), head_(head)
        { }

        size_type size() const noexcept { return head_->size; }
        bool empty() const noexcept { return !size(); }

        const_iterator begin() const noexcept { return const_iterator(el_, head_, head_->first, false); }
        const_iterator end() const noexcept { return const_iterator(el_, head_, 0, true); }

    private:
        const ElementList* el_;
        const CellListHead<THead>* head_;
    };

    // Accessors
    //-------------------------------------------------------------------------

    THead* getHeadPtr(std::size_t headIndex) const { return headList_[headIndex].ptr; }
    THead* getHeadPtr(const ElementUser& eu) const { return getHeadPtr(eu.head); }

    CellView getElementPtrs(std::size_t headIndex) const { return CellView(&elementList_, &headList_[headIndex]); }
    CellView getElementPtrs(const HeadUser& hu) const { return getElementPtrs(hu.index); }

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
        ElementUser& eu, const HeadUser& hu
    ) {
        const auto head = hu.index;

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

    void updateElement(ElementUser& eu, const HeadUser& newHu) {
        const auto newHead = newHu.index;
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
