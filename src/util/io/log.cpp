#include "util/io/log.h"

#include <chrono>
#include <ctime>
#include <iomanip>
#include <unordered_map>

namespace medyan {
namespace logger {

namespace {

    std::string timeLiteralGeneration() {
        using namespace std::chrono;

        system_clock::time_point p = system_clock::now();
        milliseconds ms = duration_cast<milliseconds>(p.time_since_epoch());
        seconds s = duration_cast<seconds>(ms);

        std::time_t timeToSec = s.count();
        tm timeinfoToSec;
#ifdef _MSC_VER
        localtime_s(&timeinfoToSec, &timeToSec);
#elif defined(__unix__) || (defined(__APPLE__) && defined(__MACH__))
        localtime_r(&timeToSec, &timeinfoToSec);
#else
        // Not thread safe
        timeinfoToSec = *localtime(&timeToSec);
#endif

        std::size_t msRemain = ms.count() % 1000;

        std::stringstream ss;
        ss << timeinfoToSec.tm_year + 1900 << '-'
            << std::setfill('0') << std::setw(2) << timeinfoToSec.tm_mon + 1 << '-'
            << std::setw(2) << timeinfoToSec.tm_mday << ' '
            << std::setw(2) << timeinfoToSec.tm_hour << ':'
            << std::setw(2) << timeinfoToSec.tm_min << ':'
            << std::setw(2) << timeinfoToSec.tm_sec << '.'
            << std::setw(3) << msRemain;

        return ss.str();
    }

} // namespace

const std::unordered_map<LogLevel, const char*> logLevelLiteral {
    {LogLevel::Debug,   "Debug"},
    {LogLevel::Info,    "Info"},
    {LogLevel::Step,    "Step"},
    {LogLevel::Note,    "Note"},
    {LogLevel::Warning, "Warning"},
    {LogLevel::Error,   "Error"},
    {LogLevel::Fatal,   "Fatal"}
};


bool Logger::defaultLoggerInitialization(const std::string& filepath) {
    Logger& l = getDefaultLogger();

    LoggerOstreamContainer& scrn = l.attachOstream(&std::cout, false);
#ifdef NDEBUG
    scrn.disp.turnOnAtLeast(LogLevel::Info);
#else
    scrn.disp.turnOnAtLeast(LogLevel::Debug);
#endif
    //scrn.dispTime.turnOnAtLeast(LogLevel::Debug);
    scrn.dispFile.turnOnAtLeast(LogLevel::Warning);
    scrn.dispLine.turnOnAtLeast(LogLevel::Error);
    scrn.dispFunc.turnOnAtLeast(LogLevel::Error);
    scrn.dispLevel.turnOnAtLeast(LogLevel::Note);
    scrn.dispLevel.turnOn(LogLevel::Debug);
    scrn.flushLevel.turnOnAtLeast(LogLevel::Debug);

    LoggerOstreamContainer& file = l.addOfstream(filepath);
#ifdef NDEBUG
    scrn.disp.turnOnAtLeast(LogLevel::Info);
#else
    scrn.disp.turnOnAtLeast(LogLevel::Debug);
#endif
    file.dispTime.turnOnAtLeast(LogLevel::Debug);
    file.dispFile.turnOnAtLeast(LogLevel::Note);
    file.dispFile.turnOn(LogLevel::Debug);
    file.dispLine.turnOnAtLeast(LogLevel::Note);
    file.dispLine.turnOn(LogLevel::Debug);
    file.dispFunc.turnOnAtLeast(LogLevel::Note);
    file.dispLevel.turnOnAtLeast(LogLevel::Debug);
    file.flushLevel.turnOnAtLeast(LogLevel::Warning);

    bool fileOpenGood = static_cast<std::ofstream*>(file.os)->is_open();

    if(!fileOpenGood) {
        l.removeOstream(file.os);
        std::cout << "Logger cannot open file " << filepath << std::endl;
    }

    return fileOpenGood;
    
}

namespace internal {

void LogWriter::logDispatch() {
    bool genTime = true; // Generate time literal only once
    std::string strTime;

    const auto& settings = _logger.settings;

    for(auto& eachOs: _logger.getOsContainers()) {
        if(eachOs.disp.isOnWith(_lv)) {
            // The oss for final output
            std::ostringstream finalOss;

            // Prefix generation
            if(eachOs.dispTime.isOnWith(_lv)) {
                if(genTime) {
                    strTime = timeLiteralGeneration();
                    genTime = false;
                }
                finalOss << settings.delimiterBefore << strTime << settings.delimiterAfter
                    << (settings.spaceAfterDelimiter? " ": "");
            }
            if(eachOs.dispLevel.isOnWith(_lv)) {
                finalOss << settings.delimiterBefore << logLevelLiteral.find(_lv)->second << settings.delimiterAfter
                    << (settings.spaceAfterDelimiter? " ": "");
            }
            if(eachOs.dispFile.isOnWith(_lv)) {
                finalOss << settings.delimiterBefore << "File " << _curFile << settings.delimiterAfter
                    << (settings.spaceAfterDelimiter? " ": "");
            }
            if(eachOs.dispLine.isOnWith(_lv)) {
                finalOss << settings.delimiterBefore << "Line " << _curLine << settings.delimiterAfter
                    << (settings.spaceAfterDelimiter? " ": "");
            }
            if(eachOs.dispFunc.isOnWith(_lv)) {
                finalOss << settings.delimiterBefore << "Function " << _curFunc << settings.delimiterAfter
                    << (settings.spaceAfterDelimiter? " ": "");
            }

            // Attach log content
            finalOss << _oss.str();

            // Suffix generation
            if(settings.newLineAfterLog) finalOss << '\n';

            // This output should not cause data races
            (*eachOs.os) << finalOss.str();

            // Flush
            if(eachOs.flushLevel.isOnWith(_lv)) (*eachOs.os) << std::flush;
        }
    }

}
} // namespace internal


} // namespace logger
} // namespace medyan
