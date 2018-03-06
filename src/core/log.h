#ifndef MEDYAN_CORE_LOG_H
#define MEDYAN_CORE_LOG_H

#include <algorithm>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <streambuf>
#include <string>
#include <vector>
#include <utility>

namespace medyan {
namespace logger {

enum class LogLevel: int {
    Debug   = 0,
    Info    = 1,
    Step    = 2,
    Note    = 3,
    Warning = 4,
    Error   = 5,
    Fatal   = 6
};

namespace internal {

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
        std::ostringstream& getStream() { return _oss; }
        void setCurLevel(LogLevel lv) { _lv = lv; }

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

        /// Dispatch log
        void generatePrefix(const char* curFile, const char* curLine, const char* curFunc);
        void dispatchStream(); ///< Also clears the stringstream

        /// Default logger
        static Logger* getDefaultLogger() {
            static Logger l;
            return &l;
        }
        /// Default initialization. Returns whether the file is successfully opened.
        static bool defaultLoggerInitialization(const std::string& filepath);
    private:
        /// The actual stringstream
        std::ostringstream _oss;
        /// Level of current log
        LogLevel _lv;

        /// Ostream containers
        std::vector<LoggerOstreamContainer> _osContainers;

        /// Managed ostreams. will be destroyed when instance of this class is going out of scope
        std::vector<std::unique_ptr<std::ostream>> _osManaged;
    };

    /// This is the class that prints the log when destructed.
    class LogWriter {
    public:
        /// Constructor accepts environmental information.
        LogWriter(const char* curFile, const char* curLine, const char* curFunc, LogLevel lv, Logger* logger);
        LogWriter(const LogWriter& lw) = delete;
        ~LogWriter();

        /// Overloads operator<< for normal type and user defined class type
        template<typename MsgType>
        LogWriter& operator<<(MsgType&& msg) {
            _logger->getStream() << std::forward<MsgType>(msg);
            return *this;
        }
        LogWriter& operator<<(std::ostream& (*manip)(std::ostream&)) {
            _logger->getStream() << manip;
            return *this;
        }
    
    private:
        Logger* _logger; ///< The pointer to the actual logger information
    };


} // namespace internal

} // namespace logger
} // namespace medyan

/// Exposed macro
#define MEDYAN_WRITE_LOG(writer, log_level) ::medyan::logger::LogWriter(__FILE__, __LINE__, __func__, log_level, writer)

#define MEDYAN_LOG_GEN_DEBUG(writer)   MEDYAN_WRITE_LOG(writer, ::medyan::logger::LogLevel::Debug)
#define MEDYAN_LOG_GEN_INFO(writer)    MEDYAN_WRITE_LOG(writer, ::medyan::logger::LogLevel::Info)
#define MEDYAN_LOG_GEN_STEP(writer)    MEDYAN_WRITE_LOG(writer, ::medyan::logger::LogLevel::Step)
#define MEDYAN_LOG_GEN_NOTE(writer)    MEDYAN_WRITE_LOG(writer, ::medyan::logger::LogLevel::Note)
#define MEDYAN_LOG_GEN_WARNING(writer) MEDYAN_WRITE_LOG(writer, ::medyan::logger::LogLevel::Warning)
#define MEDYAN_LOG_GEN_ERROR(writer)   MEDYAN_WRITE_LOG(writer, ::medyan::logger::LogLevel::Error)
#define MEDYAN_LOG_GEN_FATAL(writer)   MEDYAN_WRITE_LOG(writer, ::medyan::logger::LogLevel::Fatal)

#define MEDYAN_LOG_GEN(log_level) MEDYAN_LOG_GEN_##log_level(::medyan::logger::internal::Logger::getDefaultLogger())

/// User interface
#define LOG(log_level) MEDYAN_LOG_GEN(log_level)
#define MEDYAN_LOG_DEFAULT_CONFIGURATION(filepath) ::medyan::logger::internal::Logger::defaultLoggerInitialization(filepath)

#endif