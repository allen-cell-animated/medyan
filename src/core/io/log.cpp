#include "core/io/log.h"

#include <chrono>
#include <ctime>
#include <iomanip>
#include <map>

namespace medyan {
namespace logger {

namespace internal {

    static std::string timeLiteralGeneration() {
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

    const std::map<LogLevel, const char*> logLevelLiteral {
        {LogLevel::Debug,   "Debug"},
        {LogLevel::Info,    "Info"},
        {LogLevel::Step,    "Step"},
        {LogLevel::Note,    "Note"},
        {LogLevel::Warning, "Warning"},
        {LogLevel::Error,   "Error"},
        {LogLevel::Fatal,   "Fatal"}
    };
    

    void Logger::generatePrefix(const char* curFile, int curLine, const char* curFunc) {
        bool genTime = false; // Generate time literal only once
        std::string strTime;

        for(auto& eachOs: _osContainers) {
            if(eachOs.disp.isOnWith(_lv)) {
                if(eachOs.dispTime.isOnWith(_lv)) {
                    if(!genTime) {
                        strTime = timeLiteralGeneration();
                        genTime = true;
                    }
                    (*eachOs.os) << settings.delimiterBefore << strTime << settings.delimiterAfter
                        << (settings.spaceAfterDelimiter? " ": "");
                }
                if(eachOs.dispLevel.isOnWith(_lv)) {
                    (*eachOs.os) << settings.delimiterBefore << logLevelLiteral.find(_lv)->second << settings.delimiterAfter
                        << (settings.spaceAfterDelimiter? " ": "");
                }
                if(eachOs.dispFile.isOnWith(_lv)) {
                    (*eachOs.os) << settings.delimiterBefore << "File " << curFile << settings.delimiterAfter
                        << (settings.spaceAfterDelimiter? " ": "");
                }
                if(eachOs.dispLine.isOnWith(_lv)) {
                    (*eachOs.os) << settings.delimiterBefore << "Line " << curLine << settings.delimiterAfter
                        << (settings.spaceAfterDelimiter? " ": "");
                }
                if(eachOs.dispFunc.isOnWith(_lv)) {
                    (*eachOs.os) << settings.delimiterBefore << "Function " << curFunc << settings.delimiterAfter
                        << (settings.spaceAfterDelimiter? " ": "");
                }
            }
        }
    }
    void Logger::dispatchStream() {
        for(auto& eachOs: _osContainers) {
            if(eachOs.disp.isOnWith(_lv)) {
                (*eachOs.os) << _oss.str();
                if(settings.newLineAfterLog) (*eachOs.os) << '\n';
                if(eachOs.flushLevel.isOnWith(_lv)) (*eachOs.os) << std::flush;
            }
        }
        // Clear string stream
        _oss.clear();
        _oss.str(std::string());
    }
    bool Logger::defaultLoggerInitialization(const std::string& filepath) {
        Logger* l = getDefaultLogger();

        LoggerOstreamContainer& scrn = l->attachOstream(&std::cout, false);
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

        LoggerOstreamContainer& file = l->addOfstream(filepath);
        file.disp.turnOnAtLeast(LogLevel::Debug);
        file.dispTime.turnOnAtLeast(LogLevel::Debug);
        file.dispFile.turnOnAtLeast(LogLevel::Note);
        file.dispLine.turnOnAtLeast(LogLevel::Note);
        file.dispFunc.turnOnAtLeast(LogLevel::Note);
        file.dispLevel.turnOnAtLeast(LogLevel::Debug);
        file.flushLevel.turnOnAtLeast(LogLevel::Warning);

        bool fileOpenGood = static_cast<std::ofstream*>(file.os)->is_open();

        if(!fileOpenGood) {
            l->removeOstream(file.os);
            std::cout << "Logger cannot open file " << filepath << std::endl;
        }

        return fileOpenGood;
        
    }

    LogWriter::LogWriter(const char* curFile, int curLine, const char* curFunc, LogLevel lv, Logger* logger):
        _logger(logger) {
        _logger->setCurLevel(lv);
        _logger->generatePrefix(curFile, curLine, curFunc);
    }
    LogWriter::~LogWriter() {
        _logger->dispatchStream();
    }
} // namespace internal


} // namespace logger
} // namespace medyan
