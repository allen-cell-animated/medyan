#include "util/io/log.h"

#include <chrono>
#include <ctime>
#include <iomanip>
#include <unordered_map>

#include "util/platform.h"

#ifdef PLATFORM_UNIX_LIKE
    #include <unistd.h>
#endif

namespace medyan {
namespace logger {

namespace {

// Internal configuration
bool stdoutRedirected = false;
bool stderrRedirected = false;

// Helper functions
bool isStdoutRedirected() {
#ifdef PLATFORM_UNIX_LIKE
    return !isatty(STDOUT_FILENO);
#else
    return false;
#endif
}
bool isStderrRedirected() {
#ifdef PLATFORM_UNIX_LIKE
    return !isatty(STDERR_FILENO);
#else
    return false;
#endif
}

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

// Remove by c++14
struct LogLevelHash {
    std::size_t operator()(LogLevel lv)const { return static_cast<std::size_t>(lv); }
};

} // namespace

// Mapping log level to string. Notice that we no longer need to supply hash since c++14.
const std::unordered_map<LogLevel, const char*, LogLevelHash> logLevelLiteral {
    {LogLevel::Debug,   "Debug"},
    {LogLevel::Info,    "Info"},
    {LogLevel::Step,    "Step"},
    {LogLevel::Note,    "Note"},
    {LogLevel::Warning, "Warning"},
    {LogLevel::Error,   "Error"},
    {LogLevel::Fatal,   "Fatal"}
};

// Level color codes. Notice that we no longer need to supply hash since c++14.
const std::unordered_map<LogLevel, const char*, LogLevelHash> logLevelColorAnsi {
    {LogLevel::Debug,   "\033[37m"}, // White
    {LogLevel::Info,    "\033[97m"}, // Bright white
    {LogLevel::Step,    "\033[96m"}, // Bright cyan
    {LogLevel::Note,    "\033[92m"}, // Bright green
    {LogLevel::Warning, "\033[93m"}, // Bright yellow
    {LogLevel::Error,   "\033[91m"}, // Bright red
    {LogLevel::Fatal,   "\033[91m"}  // Bright red
};
constexpr const char * resetAnsi = "\033[0m";


void Logger::defaultLoggerInitialization() {
    Logger& l = getDefaultLogger();

    // Detect output redirection
    stdoutRedirected = isStdoutRedirected();
    stderrRedirected = isStderrRedirected();

    LoggerOstreamContainer& scrn = l.attachOstream(&std::cout, false);
#ifdef NDEBUG
    scrn.disp.turnOnAtLeast(LogLevel::Info);
#else
    scrn.disp.turnOnAtLeast(LogLevel::Debug);
#endif
    if(!(stdoutRedirected && l.settings.supressColorIfRedirected)) scrn.dispColor.turnOnAtLeast(LogLevel::Debug);
    scrn.dispTime.turnOnAtLeast(LogLevel::Warning);
    scrn.dispFile.turnOnAtLeast(LogLevel::Warning);
    scrn.dispLine.turnOnAtLeast(LogLevel::Error);
    scrn.dispFunc.turnOnAtLeast(LogLevel::Error);
    scrn.dispLevel.turnOnAtLeast(LogLevel::Note);
    scrn.dispLevel.turnOn(LogLevel::Debug);
    scrn.flushLevel.turnOnAtLeast(LogLevel::Debug);

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
            if(eachOs.dispColor.isOnWith(_lv)) {
                finalOss << logLevelColorAnsi.find(_lv)->second;
            }
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
            if(eachOs.dispColor.isOnWith(_lv)) {
                finalOss << resetAnsi;
            }
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
