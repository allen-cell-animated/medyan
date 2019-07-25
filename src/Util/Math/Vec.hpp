#ifndef MEDYAN_Util_Math_Vec_Hpp
#define MEDYAN_Util_Math_Vec_Hpp

#include <algorithm> // copy
#include <array>
#include <cmath>
#include <cstddef> // ptrdiff_t
#include <iterator> // tag
#include <ostream>
#include <type_traits> // common_type, conditional, enable_if, is_same, remove_reference

namespace mathfunc {

//-----------------------------------------------------------------------------
// Vec is a simple coordinate type that makes operations easier
// Can be replaced by linalg/tensor libraries in the future
//-----------------------------------------------------------------------------
template< size_t dim, typename Float = double > struct Vec {

    static constexpr size_t vec_size = dim;
    using float_type = Float;

    using storage_type = std::array< Float, dim >;
    using size_type = typename storage_type::size_type;
    using iterator = typename storage_type::iterator;
    using const_iterator = typename storage_type::const_iterator;

    storage_type value;

    Vec& operator=(const Vec& v) = default;
    template< typename VecType, std::enable_if_t< dim == VecType::vec_size > * = nullptr >
    Vec& operator=(const VecType& v) {
        for(size_t i = 0; i < dim; ++i) (*this)[i] = v[i];
        return *this;
    }

    // Conversion operator to another Vec
    template< typename FloatOut >
    explicit operator Vec< vec_size, FloatOut >() const {
        Vec< vec_size, FloatOut > res;
        for(size_t i = 0; i < vec_size; ++i) res[i] = (*this)[i];
        return res;
    }

    constexpr size_type size() const noexcept { return dim; }

    constexpr iterator       begin()       noexcept { return value.begin(); }
    constexpr const_iterator begin() const noexcept { return value.begin(); }
    constexpr iterator       end()       noexcept { return value.end(); }
    constexpr const_iterator end() const noexcept { return value.end(); }

    constexpr       Float& operator[](size_type pos)       { return value[pos]; }
    constexpr const Float& operator[](size_type pos) const { return value[pos]; }
};

// Frequently used type alias
using Vec3  = Vec< 3 >;
using Vec3d = Vec< 3, double >;
using Vec3f = Vec< 3, float >;

namespace internal {

    // Helper functions and types to choose implementation according to whether RefVecBase of raw_ptr_container
    template< typename RefVecType, std::enable_if_t< std::remove_reference_t<RefVecType>::raw_ptr_container >* = nullptr >
    constexpr decltype(auto) refVecBaseGetBegin(RefVecType&& v) { return v.ptr; }
    template< typename RefVecType, std::enable_if_t< !std::remove_reference_t<RefVecType>::raw_ptr_container >* = nullptr >
    constexpr decltype(auto) refVecBaseGetBegin(RefVecType&& v) { return v.ptr->begin(); }
    template< typename RefVecType, std::enable_if_t< std::remove_reference_t<RefVecType>::raw_ptr_container >* = nullptr >
    constexpr decltype(auto) refVecBaseGetContainer(RefVecType&& v) { return v.ptr; }
    template< typename RefVecType, std::enable_if_t< !std::remove_reference_t<RefVecType>::raw_ptr_container >* = nullptr >
    constexpr decltype(auto) refVecBaseGetContainer(RefVecType&& v) { return *v.ptr; }

    template< bool is_const, bool raw_ptr_container, typename Float, typename Container >
    struct RefVecBaseIterator { using type = std::conditional_t< is_const, const Float *, Float * >; };
    template< bool is_const, typename Float, typename Container >
    struct RefVecBaseIterator< is_const, false, Float, Container > {
        using type = std::conditional_t< is_const, typename Container::const_iterator, typename Container::iterator >;
    };
    template< bool is_const, bool raw_ptr_container, typename Float, typename Container >
    using RefVecBaseIterator_t = typename RefVecBaseIterator< is_const, raw_ptr_container, Float, Container >::type;

    template< bool is_const, bool raw_ptr_container, typename Float, typename Container >
    struct RefVecBaseReference { using type = std::conditional_t< is_const, const Float &, Float & >; };
    template< bool is_const, typename Float, typename Container >
    struct RefVecBaseReference< is_const, false, Float, Container > {
        using type = std::conditional_t< is_const, typename Container::const_reference, typename Container::reference >;
    };
    template< bool is_const, bool raw_ptr_container, typename Float, typename Container >
    using RefVecBaseReference_t = typename RefVecBaseReference< is_const, raw_ptr_container, Float, Container >::type;

    // RefVecBase is the base impl for RefVec and ConstRefVec
    // Note:
    //   - There must be a Container::iterator that is a RandomAccessIterator
    //     unless the Container is a Float* raw pointer
    template<
        size_t dim,
        typename Float,
        typename Container,
        bool is_const,
        typename Concrete
    > struct RefVecBase {
        static constexpr size_t vec_size = dim;
        using float_type = Float;

        static constexpr bool raw_ptr_container = std::is_same< Float*, Container >::value;

        using container_type = std::conditional_t< is_const, const Container, Container >;
        using iterator = RefVecBaseIterator_t< is_const, raw_ptr_container, Float, Container >;
        using reference = RefVecBaseReference_t< is_const, raw_ptr_container, Float, Container >;
        using container_ptr_type = std::conditional_t< raw_ptr_container, container_type, container_type * >;

        container_ptr_type ptr;
        size_t pos; // index of first Float

        // Conversion operator to normal Vec
        template< typename FloatOut >
        explicit operator Vec< vec_size, FloatOut >() const {
            Vec< vec_size, FloatOut > res;
            for(size_t i = 0; i < vec_size; ++i) res[i] = (*this)[i];
            return res;
        }

        // Copy/move assignment operators are deleted
        RefVecBase& operator=(const RefVecBase & ) = delete;
        RefVecBase& operator=(      RefVecBase &&) = delete;

        constexpr auto size() const noexcept { return dim; }

        // Iterators
        constexpr iterator begin() const noexcept { return refVecBaseGetBegin(*this) + pos; }
        constexpr iterator end()   const noexcept { return refVecBaseGetBegin(*this) + pos + dim; }

        // sub_pos must be within [0, dim)
        reference operator[](size_t sub_pos) const { return refVecBaseGetContainer(*this)[pos + sub_pos]; }

        // also works like pointer
        Concrete*       operator->()       noexcept { return static_cast<      Concrete*>(this); }
        const Concrete* operator->() const noexcept { return static_cast<const Concrete*>(this); }
    };

} // namespace internal

//-----------------------------------------------------------------------------
// RefVec and ConstRefVec refer to the Vec-like data in an array. They contain
// pointers to the container and relative positions in the array, but do not
// contain data themselves.
//
// They can perform most of the arithmetics just like a normal Vec does. A
// RefVec can also change the data it refers to. If an arithmetic generates
// a new object, then the object should NOT be a RefVec.
//
// Note:
//   - The copy/move constructors here have different semantics from copy/move
//     assignment operators. Both RefVec and ConstRefVec have copy/move
//     constructors to copy the pointer information from another RefVec. On the
//     other hand, assignment operators deal directly with the data they refer
//     to, just like a normal Vec.
//-----------------------------------------------------------------------------
template< size_t dim, typename Float, typename Container >
struct RefVec : internal::RefVecBase< dim, Float, Container, false, RefVec< dim, Float, Container > > {
    using base_type = internal::RefVecBase< dim, Float, Container, false, RefVec >;

    constexpr RefVec(typename base_type::container_ptr_type ptr, size_t pos) : base_type{ptr, pos} {}

    // Copy/move constructors
    constexpr RefVec(const RefVec&  rv) : RefVec(rv.ptr, rv.pos) {}
    constexpr RefVec(      RefVec&& rv) : RefVec(rv.ptr, rv.pos) {}

    RefVec& operator=(const RefVec& v) {
        for(size_t i = 0; i < dim; ++i) (*this)[i] = v[i];
        return *this;
    }

    template< typename VecType, std::enable_if_t< dim == VecType::vec_size > * = nullptr >
    RefVec& operator=(const VecType& v) {
        for(size_t i = 0; i < dim; ++i) (*this)[i] = v[i];
        return *this;
    }
    template< typename VecType, std::enable_if_t< dim == VecType::vec_size > * = nullptr >
    RefVec& operator=(VecType&& v) {
        for(size_t i = 0; i < dim; ++i) (*this)[i] = v[i];
        return *this;
    }
};
template< size_t dim, typename Float, typename Container >
struct ConstRefVec : internal::RefVecBase< dim, Float, Container, true, ConstRefVec< dim, Float, Container > > {
    using base_type = internal::RefVecBase< dim, Float, Container, true, ConstRefVec >;

    constexpr ConstRefVec(typename base_type::container_ptr_type ptr, size_t pos) : base_type{ptr, pos} {}

    constexpr ConstRefVec(const RefVec< dim, Float, Container >&  rv) : ConstRefVec(rv.ptr, rv.pos) {}
    constexpr ConstRefVec(      RefVec< dim, Float, Container >&& rv) : ConstRefVec(rv.ptr, rv.pos) {}

    // Copy/move constructors
    constexpr ConstRefVec(const ConstRefVec&  crv) : ConstRefVec(crv.ptr, crv.pos) {}
    constexpr ConstRefVec(      ConstRefVec&& crv) : ConstRefVec(crv.ptr, crv.pos) {}
};

//-----------------------------------------------------------------------------
// The VecArray type is like an array of Vec's but has flattened storage.
//-----------------------------------------------------------------------------
template<
    size_t dim,
    typename Float = double,
    typename Container = std::vector< Float >
> struct VecArray {

    static constexpr size_t element_vec_size = dim;
    using float_type = Float;
    using container_type = Container;

    using size_type = typename container_type::size_type;
    static_assert(std::is_same<typename container_type::iterator::iterator_category, std::random_access_iterator_tag>::value,
        "The iterator of the VecArray container must be random access iterator.");

    using reference       = RefVec     < dim, Float, Container >;
    using const_reference = ConstRefVec< dim, Float, Container >;

    template< bool is_const > class VecIterator {
        template< bool friend_const > friend class VecIterator;

    public:
        using vec_array_type = VecArray;
        using size_type = vec_array_type::size_type;
        using SolidVec = Vec< dim, Float >;
        using container_type = std::conditional_t< is_const, const vec_array_type::container_type, vec_array_type::container_type >;

        using iterator_category = std::random_access_iterator_tag;
        using value_type = SolidVec; // Not used in class
        using difference_type = std::ptrdiff_t;
        using pointer = std::conditional_t< is_const, VecArray::const_reference, VecArray::reference >; // The pointer type is the reference itself
        using reference = std::conditional_t< is_const, VecArray::const_reference, VecArray::reference >;

    private:
        container_type* _ptr;
        size_type _index;

    public:
        VecIterator() = default;
        VecIterator(container_type* ptr, size_type index) : _ptr(ptr), _index(index) {}
        VecIterator(const VecIterator& rhs) = default;
        template< bool rhs_const >
        VecIterator(const VecIterator<rhs_const>& rhs) : _ptr(rhs._ptr), _index(rhs._index) {}
        VecIterator& operator=(const VecIterator& rhs) = default;
        template< bool rhs_const >
        VecIterator& operator=(const VecIterator<rhs_const>& rhs) { _ptr = rhs._ptr; _index = rhs._index; return *this; }

        reference operator*() const { return reference(_ptr, _index * dim); }
        pointer operator->() const { return pointer(_ptr, _index * dim); }
        reference operator[](difference_type rhs) const { return reference(_ptr, (_index + rhs) * dim); }

        VecIterator& operator+=(difference_type rhs) { _index += rhs; return *this; }
        VecIterator& operator-=(difference_type rhs) { _index -= rhs; return *this; }
        VecIterator& operator++() { ++_index; return *this; }
        VecIterator& operator--() { --_index; return *this; }
        VecIterator operator++(int) { VecIterator tmp(*this); ++_index; return tmp; }
        VecIterator operator--(int) { VecIterator tmp(*this); --_index; return tmp; }

        template< bool rhs_const >
        difference_type operator-(const VecIterator<rhs_const>& rhs) const {
            return difference_type(_index) - difference_type(rhs._index);
        }
        VecIterator operator+(difference_type rhs) const { return VecIterator(_ptr, _index + rhs); }
        VecIterator operator-(difference_type rhs) const { return VecIterator(_ptr, _index - rhs); }
        friend auto operator+(difference_type lhs, const VecIterator& rhs) {
            return VecIterator(rhs._ptr, lhs + rhs._index);
        }

        template< bool rhs_const >
        bool operator==(const VecIterator<rhs_const>& rhs) const { return _ptr == rhs._ptr && _index == rhs._index; }
        template< bool rhs_const >
        bool operator!=(const VecIterator<rhs_const>& rhs) const { return _ptr != rhs._ptr || _index != rhs._index; }
        template< bool rhs_const >
        bool operator>(const VecIterator<rhs_const>& rhs) const { return _ptr > rhs._ptr || (_ptr == rhs._ptr && _index > rhs._index); }
        template< bool rhs_const >
        bool operator<(const VecIterator<rhs_const>& rhs) const { return _ptr < rhs._ptr || (_ptr == rhs._ptr && _index < rhs._index); }
        template< bool rhs_const >
        bool operator>=(const VecIterator<rhs_const>& rhs) const { return _ptr > rhs._ptr || (_ptr == rhs._ptr && _index >= rhs._index); }
        template< bool rhs_const >
        bool operator<=(const VecIterator<rhs_const>& rhs) const { return _ptr < rhs._ptr || (_ptr == rhs._ptr && _index <= rhs._index); }
    };

    using iterator = VecIterator<false>;
    using const_iterator = VecIterator<true>;

    container_type value;

    // returns raw data (Float*, with size equal to value.size())
    Float*       data()       noexcept(noexcept(std::declval<      container_type>().data())) { return value.data(); }
    const Float* data() const noexcept(noexcept(std::declval<const container_type>().data())) { return value.data(); }

    size_type size_raw() const noexcept(noexcept(std::declval<container_type>().data())) { return value.size(); }
    bool empty() const noexcept(noexcept(std::declval<container_type>().empty())) { return value.empty(); }

    // The following considers size as *number of vectors*, which is value.size() / dim
    size_type size() const noexcept(noexcept(std::declval<container_type>().data())) { return value.size() / dim; }

    iterator       begin()       noexcept { return       iterator(&value, 0); }
    const_iterator begin() const noexcept { return const_iterator(&value, 0); }
    iterator       end()       noexcept { return       iterator(&value, size()); }
    const_iterator end() const noexcept { return const_iterator(&value, size()); }

    reference       operator[](size_type index)       { return       reference(&value, index * dim); }
    const_reference operator[](size_type index) const { return const_reference(&value, index * dim); }

    reference       back()       { return       reference(&value, size_raw() - dim); }
    const_reference back() const { return const_reference(&value, size_raw() - dim); }

    template< typename VecType, std::enable_if_t<dim == VecType::vec_size>* = nullptr >
    void push_back(const VecType& v) {
        value.insert(value.end(), v.begin(), v.end());
    }
    void pop_back() {
        value.resize(value.size() - dim);
    }
};

//-----------------------------------------------------------------------------
// Factory functions
//-----------------------------------------------------------------------------
template< size_t dim, typename Float = double >
constexpr auto makeRefVec     (const Float* source) { return ConstRefVec< dim, Float, Float* >(source, 0); }
template< size_t dim, typename Float = double >
constexpr auto makeRefVec     (      Float* source) { return      RefVec< dim, Float, Float* >(source, 0); }
template< size_t dim, typename Float = double >
constexpr auto makeConstRefVec(const Float* source) { return ConstRefVec< dim, Float, Float* >(source, 0); }

template< size_t dim, typename Float = double >
inline auto makeVec(const Float* source) {
    Vec< dim, Float > res;
    std::copy(source, source + dim, res.begin());
    return res;
}

// Formatting
// undefined behavior if size is 0
template< typename VecType, size_t = VecType::vec_size > inline
std::ostream& operator<<(std::ostream& os, const VecType& v) {
    os << '(';
    for(size_t i = 0; i < VecType::vec_size - 1; ++i) os << v[i] << ", ";
    os << v[VecType::vec_size - 1] << ')';
    return os;
}

// Fixed size vector arithmetics

// Magnitude, distance and normalization
template< typename VecType, size_t = VecType::vec_size >
inline auto magnitude2(const VecType& v) {
    typename VecType::float_type mag2 {};
    for(size_t i = 0; i < VecType::vec_size; ++i) mag2 += v[i] * v[i];
    return mag2;
}
template< typename VecType, size_t = VecType::vec_size >
inline auto magnitude(const VecType& v) {
    return std::sqrt(magnitude2(v));
}
template< typename VecType, size_t = VecType::vec_size >
inline void normalize(VecType& v) {
    typename VecType::float_type norm = magnitude(v);
    for(size_t i = 0; i < VecType::vec_size; ++i) v[i] /= norm;
}
template< typename VecType, size_t = VecType::vec_size >
inline auto normalizedVector(const VecType& v) {
    Vec< VecType::vec_size, typename VecType::float_type > res = v;
    normalize(res);
    return res;
}
template< typename VT1, typename VT2, std::enable_if_t<VT1::vec_size == VT2::vec_size>* = nullptr >
inline auto distance2(const VT1& v1, const VT2& v2) {
    std::common_type_t<typename VT1::float_type, typename VT2::float_type> res {};
    for(size_t idx = 0; idx < VT1::vec_size; ++idx) {
        res += (v2[idx] - v1[idx]) * (v2[idx] - v1[idx]);
    }
    return res;
}
template< typename VT1, typename VT2, std::enable_if_t<VT1::vec_size == VT2::vec_size>* = nullptr >
inline auto distance(const VT1& v1, const VT2& v2) {
    return std::sqrt(distance2(v1, v2));
}

// plus, minus, multiply, divide
template< typename VecType, size_t = VecType::vec_size >
inline auto operator-(const VecType& v){
    Vec< VecType::vec_size, typename VecType::float_type > res;
    for(size_t idx = 0; idx < VecType::vec_size; ++idx){
        res[idx] = -v[idx];
    }
    return res;
}
template< typename VT1, typename VT2, std::enable_if_t<VT1::vec_size == VT2::vec_size>* = nullptr >
inline auto operator+(const VT1& v1, const VT2& v2) {
    Vec< VT1::vec_size, std::common_type_t<typename VT1::float_type, typename VT2::float_type> > res;
    for(size_t idx1 = 0; idx1 < VT1::vec_size; ++idx1) {
        res[idx1] = v1[idx1] + v2[idx1];
    }
    return res;
}
template< typename VT1, typename VT2, std::enable_if_t<VT1::vec_size == VT2::vec_size>* = nullptr >
inline auto& operator+=(VT1& v1, const VT2& v2) {
    for(size_t idx = 0; idx < VT1::vec_size; ++idx) {
        v1[idx] += v2[idx];
    }
    return v1;
}
template< typename VT1, typename VT2, std::enable_if_t<VT1::vec_size == VT2::vec_size>* = nullptr >
inline auto operator-(const VT1& v1, const VT2& v2) {
    Vec< VT1::vec_size, std::common_type_t<typename VT1::float_type, typename VT2::float_type> > res;
    for(size_t idx1 = 0; idx1 < VT1::vec_size; ++idx1) {
        res[idx1] = v1[idx1] - v2[idx1];
    }
    return res;
}
template< typename VT1, typename VT2, std::enable_if_t<VT1::vec_size == VT2::vec_size>* = nullptr >
inline auto& operator-=(VT1& v1, const VT2& v2) {
    for(size_t idx = 0; idx < VT1::vec_size; ++idx) {
        v1[idx] -= v2[idx];
    }
    return v1;
}
template< typename VecType, typename Float, size_t = VecType::vec_size >
inline auto operator*(const VecType& v, Float k) {
    Vec< VecType::vec_size, std::common_type_t<typename VecType::float_type, Float> > res;
    for(size_t idx = 0; idx < VecType::vec_size; ++idx){
        res[idx] = v[idx] * k;
    }
    return res;
}
template< typename VecType, typename Float, size_t = VecType::vec_size >
inline auto operator*(Float k, const VecType& v) {
    return v * k;
}
template< typename VecType, typename Float, size_t = VecType::vec_size >
inline auto& operator*=(VecType& v, Float k) {
    for(size_t i = 0; i < VecType::vec_size; ++i) v[i] *= k;
    return v;
}
template< typename VecType, typename Float, size_t = VecType::vec_size >
inline auto operator/(const VecType& v, Float k) {
    return v * (static_cast< std::common_type_t<typename VecType::float_type, Float> >(1.0) / k);
}
template< typename VecType, typename Float, size_t = VecType::vec_size >
inline auto& operator/=(VecType& v, Float k) {
    return v *= (static_cast< std::common_type_t<typename VecType::float_type, Float> >(1.0) / k);
}

// dot product, cross product
template< typename VT1, typename VT2, std::enable_if_t<VT1::vec_size == VT2::vec_size>* = nullptr >
inline auto dot(const VT1& v1, const VT2& v2) {
    std::common_type_t<typename VT1::float_type, typename VT2::float_type> res {};
    for(size_t idx = 0; idx < VT1::vec_size; ++idx)
        res += v1[idx] * v2[idx];
    return res;
}
template<
    typename VT1, typename VT2,
    std::enable_if_t<VT1::vec_size == 3 && VT2::vec_size == 3>* = nullptr
> inline
auto cross(const VT1& v1, const VT2& v2) {
    return Vec< 3, std::common_type_t<typename VT1::float_type, typename VT2::float_type> > {
        v1[1]*v2[2] - v1[2]*v2[1],
        v1[2]*v2[0] - v1[0]*v2[2],
        v1[0]*v2[1] - v1[1]*v2[0]
    };
}

// VecArray arithmetics

// Dot product
// Always using the size of the 1st operand.
// Doing dot product on arrays with different sizes leads to undefined behavior.
template< typename VA1, typename VA2, std::enable_if_t<VA1::element_vec_size == VA2::element_vec_size>* = nullptr >
inline auto dot(const VA1& v1, const VA2& v2) {
    using res_type = std::common_type_t< typename VA1::float_type, typename VA2::float_type >;
    res_type res {};

    const size_t num = v1.value.size();
    for(size_t i = 0; i < num; ++i) {
        res += v1.value[i] * v2.value[i];
    }
    return res;
}

} // namespace mathfunc

#endif