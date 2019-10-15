#ifndef MEDYAN_ANALYSIS_IO_ReadSnapshot_hpp
#define MEDYAN_ANALYSIS_IO_ReadSnapshot_hpp

#include <string>

namespace medyan {
namespace analysis {

class SnapshotReader {
private:
    std::string _snapshotFilepath;
    std::string _pdbFilepath;
    std::string psfFileDir_;
    std::string psfFilenameMain_;

public:
    SnapshotReader(const std::string& snapshotFilepath, const std::string& pdbFilepath, const std::string& psfFileDir, const std::string& psfFilenameMain):
        _snapshotFilepath(snapshotFilepath), _pdbFilepath(pdbFilepath), psfFileDir_(psfFileDir), psfFilenameMain_(psfFilenameMain) {}

    void readAndConvertToVmd(const size_t maxFrames=0); // 0 means no limit on max frames
};

} // namespace analysis
} // namespace medyan

#endif