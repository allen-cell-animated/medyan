#include "util/restart_io.h"

#include <array>
#include <cstdint>
#include <istream>

#include "util/types.h"

namespace restart {

namespace {

template< typename T >
struct HelperIstreamReadAnyType {
    std::istream& operator()(std::istream& is, T& field) const {
        return is.read(reinterpret_cast<char*>(&field), sizeof(T));
    }
};

struct HelperIstreamFillBuffer {
    std::istream& operator()(std::istream& is, vec_byte& buf) const {
        return is.read(reinterpret_cast<char*>(buf.data()), buf.size());
    }
}

} // namespace

void LayerFilament::read() {
    vec_byte buf;

    // Read header
    std::uint32_t numBeads;
    std::uint16_t filType;
    std::uint16_t checksum;
    HelperIstreamReadAnyType()(fs, numBeads);
    HelperIstreamReadAnyType()(fs, filType);
    HelperIstreamReadAnyType()(fs, checksum);

    // Read content
    buf.resize(3 * numBeads * sizeof(mreal));
    HelperIstreamFillBuffer()(fs, buf);
    auto * bufStart = reinterpret_cast<mreal const*>(buf.data());
    std::vector<mreal> beadData(bufStart, 3 * numBeads + bufStart);

    // Collect the data
    // TODO
}

} // namespace restart
