#ifndef MEDYAN_ANALYSIS_IO_READ_SNAPSHOT_H
#define MEDYAN_ANALYSIS_IO_READ_SNAPSHOT_H

#include <string>

namespace medyan {
namespace analysis {

class SnapshotReader {
private:
    std::string _snapshotFilepath;
    std::string _pdbFilepath;
    std::string _psfFilepath;

public:
    SnapshotReader(const std::string& snapshotFilepath, const std::string& pdbFilepath, const std::string& psfFilepath):
        _snapshotFilepath(snapshotFilepath), _pdbFilepath(pdbFilepath), _psfFilepath(psfFilepath) {}

    void readAndConvertToVmd(const size_t maxFrames=0); // 0 means no limit on max frames
};

} // namespace analysis
} // namespace medyan

#endif