#include "util/restart_io.h"

#include <algorithm>
#include <cstdint>
#include <istream>

#include "util/types.h"

namespace restart {

namespace {

struct HelperIstreamReadAnyType {
    template< typename T >
    std::istream& operator()(std::istream& is, T& field) const {
        return is.read(reinterpret_cast<char*>(&field), sizeof(T));
    }
};
struct HelperOstreamWriteAnyType {
    template< typename T >
    std::ostream& operator()(std::ostream& os, const T& field) const {
        return os.write(reinterpret_cast<char const*>(&field), sizeof(T));
    }
};

struct HelperIstreamFillBuffer {
    std::istream& operator()(std::istream& is, vec_byte& buf) const {
        return is.read(reinterpret_cast<char*>(buf.data()), buf.size());
    }
};

struct HelperOstreamPrintBuffer {
    std::ostream& operator()(std::ostream& os, const vec_byte& buf) const {
        return os.write(reinterpret_cast<char const*>(buf.data()), buf.size());
    }
};

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

    // TODO: validate checksum

    // Collect the data
    auto * bufStart = reinterpret_cast<mreal const*>(buf.data());
    thisFilamentData.beadPosition.resize(numBeads);
    for(auto& eachBead : thisFilamentData.beadPosition) {
        for(auto& eachPos : eachBead) {
            eachPos = *(bufStart++);
        }
    }
}
void LayerFilament::write() {
    vec_byte buf;

    std::uint32_t numBeads = thisFilamentData.beadPosition.size();
    std::uint16_t filType = thisFilamentData.type;
    std::uint16_t checksum;

    // Create content
    buf.resize(3 * numBeads * sizeof(mreal));
    auto * bufStart = reinterpret_cast<mreal*>(buf.data());
    for(const auto& eachBead : thisFilamentData.beadPosition) {
        for(const auto& eachPos : eachBead) {
            *(bufStart++) = eachPos;
        }
    }

    // Create header
    // TODO: create checksum
    HelperOstreamWriteAnyType()(fs, numBeads);
    HelperOstreamWriteAnyType()(fs, filType);
    HelperOstreamWriteAnyType()(fs, checksum);

    // Print content
    HelperOstreamPrintBuffer()(fs, buf);
}

} // namespace restart
