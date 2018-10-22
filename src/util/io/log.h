#ifndef MEDYAN_CORE_IO_LOG_H
#define MEDYAN_CORE_IO_LOG_H

#include <algorithm> // find_if
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>
#include <utility>

namespace medyan {
namespace logger {

/**
 * This file implements a versatile yet easy to use logger.
 * 
 * To use the logger, first configure the output file path via macro
 *     MEDYAN_LOG_DEFAULT_CONFIGURATION(filepath)
 * Then use
 *     LOG(level) << ... << ... ;
 * just as `cout` to invoke logging. Note that `endl` is not needed for default
 * configuration.
 * 
 * The LOG macro actually creates a temporary LogWriter object which will be
 * destroyed at the end of the expression when the log will be printed. The
 * logging is supposed to be thread-safe as long as the underlying `ostream`s
 * are thread-safe.
 * 
 * Note: the logger is designed for easy logging of run-time information with
 * different severity levels, and it should NOT handle the heavy duty output,
 * since it is less efficient than the standard output methods.
 */

enum class LogLevel: int {
    Debug   = 0,
    Info    = 1,
    Step    = 2,
    Note    = 3,
    Warning = 4,
    Error   = 5,
    Fatal   = 6
};

class LoggerLevelFlag {
    int _flag = 0;

    /// Helper function to convert loglevel to flag bit.
    static constexpr int _convert(LogLevel lv) { return 1 << static_cast<int>(lv); }
public:
    int value()const { return _flag; }

    /// Flag manipulations
    bool isOnWith(LogLevel lv)const { return _flag & _convert(lv); }

    void turnOn(LogLevel lv) { _flag |= _convert(lv); }
    void turnOff(LogLevel lv) { _flag &= ~_convert(lv); }
    void turnOnAtLeast(LogLevel lv) {
        switch(lv) {
        case LogLevel::Debug:   turnOn(LogLevel::Debug);
        case LogLevel::Info:    turnOn(LogLevel::Info);
        case LogLevel::Step:    turnOn(LogLevel::Step);
        case LogLevel::Note:    turnOn(LogLevel::Note);
        case LogLevel::Warning: turnOn(LogLevel::Warning);
        case LogLevel::Error:   turnOn(LogLevel::Error);
        case LogLevel::Fatal:   turnOn(LogLevel::Fatal);
        }
    }
};
struct LoggerOstreamContainer {
    std::ostream* os;

    /// For ofstream only
    bool isOfstream;
    std::string filepath;

    /// Display flags
    LoggerLevelFlag disp; ///< Main switch
    LoggerLevelFlag dispTime;
    LoggerLevelFlag dispFile;
    LoggerLevelFlag dispLine;
    LoggerLevelFlag dispFunc;
    LoggerLevelFlag dispLevel;

    /// Other settings
    LoggerLevelFlag flushLevel; ///< Whether to flush after each log.

    LoggerOstreamContainer(std::ostream* os, bool isOfstream): os(os), isOfstream(isOfstream) {}
};

struct LoggerSettings {
    std::string delimiterBefore = "[";
    std::string delimiterAfter = "]";
    bool spaceAfterDelimiter = true;

    bool newLineAfterLog = true;
};

/// Stores configurations of the logger
class Logger {
public:
    /// Logger settings
    LoggerSettings settings;

    /// Getters and setters
    const std::vector<LoggerOstreamContainer>& getOsContainers()const { return _osContainers; }

    /// New ostreams
    LoggerOstreamContainer& attachOstream(std::ostream* os, bool isOfstream) {
        _osContainers.emplace_back(os, isOfstream);
        return _osContainers.back();
    }
    LoggerOstreamContainer& addOfstream(const std::string& filepath) {
        _osManaged.emplace_back(new std::ofstream(filepath));
        _osContainers.emplace_back(_osManaged.back().get(), true);
        _osContainers.back().filepath = filepath;
        return _osContainers.back();
    }
    void removeOstream(std::ostream* os) {
        auto itContainer = std::find_if(
            _osContainers.begin(), _osContainers.end(),
            [&](LoggerOstreamContainer& loc) { return loc.os == os; }
        );
        if(itContainer != _osContainers.end()) _osContainers.erase(itContainer);

        auto itManaged = std::find_if(
            _osManaged.begin(), _osManaged.end(),
            [&](std::unique_ptr<std::ostream>& p) { return p.get() == os; }
        );
        if(itManaged != _osManaged.end()) _osManaged.erase(itManaged);
    }

    /// Default logger
    static Logger& getDefaultLogger() {
        static Logger l;
        return l;
    }
    /// Default initialization. Returns whether the file is successfully opened.
    static void defaultLoggerInitialization();
private:
    /// The actual stringstream
    std::ostringstream _oss;

    /// Ostream containers
    std::vector<LoggerOstreamContainer> _osContainers;

    /// Managed ostreams. will be destroyed when instance of this class is going out of scope
    std::vector<std::unique_ptr<std::ostream>> _osManaged;
};

namespace internal {
/// This is the class that prints the log when destructed.
class LogWriter {
public:
    /// Constructor accepts environmental information.
    LogWriter(const char* curFile, int curLine, const char* curFunc, LogLevel lv, const Logger& logger)
        : _logger(logger), _lv(lv), _curFile(curFile), _curLine(curLine), _curFunc(curFunc) {}

    /// Destructor dispatches the log.
    ~LogWriter() { logDispatch(); }

    /// Copy is not allowed
    LogWriter(const LogWriter&) = delete;
    LogWriter& operator=(const LogWriter&) = delete;

    /// Log generation and dispatch
    void logDispatch();

    /// Overloads operator<< for normal type and user defined class type
    template<typename MsgType>
    LogWriter& operator<<(MsgType&& msg) {
        _oss << std::forward<MsgType>(msg);
        return *this;
    }
    LogWriter& operator<<(std::ostream& (*manip)(std::ostream&)) {
        _oss << manip;
        return *this;
    }

private:
    /// The ref to the logger information
    const Logger& _logger;

    /// The content of the log
    std::ostringstream _oss;
    /// Level of current log
    const LogLevel _lv;
    /// Other information of the log
    const char* _curFile;
    const int _curLine;
    const char* _curFunc;
};
} // namespace internal

} // namespace logger
} // namespace medyan

/// Exposed macro
#define MEDYAN_WRITE_LOG(whichLogger, logLevel) ::medyan::logger::internal::LogWriter(__FILE__, __LINE__, __func__, logLevel, whichLogger)

#define MEDYAN_LOG_GEN_DEBUG(whichLogger)   MEDYAN_WRITE_LOG(whichLogger, ::medyan::logger::LogLevel::Debug)
#define MEDYAN_LOG_GEN_INFO(whichLogger)    MEDYAN_WRITE_LOG(whichLogger, ::medyan::logger::LogLevel::Info)
#define MEDYAN_LOG_GEN_STEP(whichLogger)    MEDYAN_WRITE_LOG(whichLogger, ::medyan::logger::LogLevel::Step)
#define MEDYAN_LOG_GEN_NOTE(whichLogger)    MEDYAN_WRITE_LOG(whichLogger, ::medyan::logger::LogLevel::Note)
#define MEDYAN_LOG_GEN_WARNING(whichLogger) MEDYAN_WRITE_LOG(whichLogger, ::medyan::logger::LogLevel::Warning)
#define MEDYAN_LOG_GEN_ERROR(whichLogger)   MEDYAN_WRITE_LOG(whichLogger, ::medyan::logger::LogLevel::Error)
#define MEDYAN_LOG_GEN_FATAL(whichLogger)   MEDYAN_WRITE_LOG(whichLogger, ::medyan::logger::LogLevel::Fatal)

#define MEDYAN_LOG_GEN(logLevel) MEDYAN_LOG_GEN_##logLevel(::medyan::logger::Logger::getDefaultLogger())

/// User interface
#define LOG(logLevel) MEDYAN_LOG_GEN(logLevel)

#endif
