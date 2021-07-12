#ifndef MEDYAN_Util_SExpr_hpp
#define MEDYAN_Util_SExpr_hpp

#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

#include "Util/Io/Log.hpp"
#include "utility.h"

namespace medyan {

struct SExpr {
    using StringType = std::string;
    using ListType   = std::vector< SExpr >;

    std::variant<
        StringType,
        ListType
    > data;
};


//-----------------------------------------------------------------------------
// Common functions of s-expressions
//-----------------------------------------------------------------------------

inline SExpr car(const SExpr& se) {
    return std::visit(
        Overload {
            [](const SExpr::StringType&) -> SExpr {
                LOG(ERROR) << "Expected cons in car, but got a string.";
                throw std::runtime_error("Invalid argument in car");
            },
            [](const SExpr::ListType& l) -> SExpr {
                return SExpr { l[0] };
            }
        },
        se.data
    );
}


} // namespace medyan

#endif
