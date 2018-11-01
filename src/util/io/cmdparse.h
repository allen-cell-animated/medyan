#ifndef MEDYAN_UTIL_IO_CMDPARSE_H
#define MEDYAN_UTIL_IO_CMDPARSE_H

#include <functional>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace cmdparse {

// Exceptions
class CommandLogicError : public std::logic_error {
public:
    explicit CommandLogicError( const std::string& what_arg ) : std::logic_error(what_arg) {}
    explicit CommandLogicError( const char* what_arg ) : std::logic_error(what_arg) {}
};
class ParsingError : public std::runtime_error {
public:
    explicit ParsingError( const std::string& what_arg ) : std::runtime_error(what_arg) {}
    explicit ParsingError( const char* what_arg ) : std::runtime_error(what_arg) {}
};
class ValidationError : public std::runtime_error {
public:
    explicit ValidationError( const std::string& what_arg ) : std::runtime_error(what_arg) {}
    explicit ValidationError( const char* what_arg ) : std::runtime_error(what_arg) {}
};

// Positional argument
class PosArg {
private:
    const char * const _name;
    const bool _list;
    const bool _required;

    struct State {
        size_t _occurenceCount = 0;
    } _state;

public:
    bool isList() const { return _list; }
    bool isRequired()const { return _required; }

    const char* getName()const { return _name; }

    const std::function<void(const std::string&)> activate;

    void occur() { ++_state._occurenceCount; }
    size_t getOccurenceCount()const { return _state._occurenceCount; }

};

class Option {
private:

    const char _short = 0; // without "-"
    const std::string _long; // without "--"

    const bool _hasVariable;
    const bool _required;

    struct State {
        size_t _occurenceCount = 0;
    } _state;
public:
    char getShortName()const { return _short; }
    const std::string& getLongName()const { return _long; }
    std::string getReadableName()const {
        std::ostringstream oss;
        if(_short) {
            oss << '-' << _short;
            if(_long.length())
                oss << " (--" << _long << ')';
        } else if(_long.length())
            oss << "--" << _long;
        return oss.str();
    }

    bool hasVariable()const { return _hasVariable; }
    bool isRequired()const { return _required; }

    const std::function<void(const std::string&)> activate;

    void occur() { ++_state._occurenceCount; }
    size_t getOccurenceCount()const { return _state._occurenceCount; }
};


class Command {
private:

    const std::string _name;
    const std::string _inheritedName;

    std::vector<std::unique_ptr<PosArg>> _posArgs;
    std::vector<std::unique_ptr<Option>> _options;
    std::vector<std::unique_ptr<Command>> _subcommands;
    Command* const _parent = nullptr;

    std::function<void()> _userValidation;

    // Check validity of specified rules. Recursive check for subcommands.
    void _ruleCheck()const;
    // Real parsing function
    void _parse(std::vector<std::string>& feed, size_t argp);
    void _parsePosArg(const std::string& arg); // helper
    // Internal validation
    void _validate()const;

    // State variable
    struct State {
        bool _parseAsPosArg = false; // After the delimiter, this will be set to true.
        size_t _posArgCount = 0; // Number of positional argument encountered.
        size_t _posArgIndex = 0; // The index for the next positional argument.
                                 // Normally same with _posArgCount except on arg list.
        size_t _occurenceCount = 0;
    } _state;

public:

    const std::string& getName()const { return _name; }
    std::string getFullName()const { return _inheritedName + ' ' + _name; }

    const std::function<void()> activate;

    void occur() { ++_state._occurenceCount; }
    // The parsing function used only by the root command.
    // States will be changed in this function
    void parse(int argc, char** argv) {
        occur();
        if(activate) activate();

        _ruleCheck();

        std::vector<std::string> inputFeed(argc);
        for(size_t i = 0; i < argc; ++i) inputFeed[i] = argv[i];
        _parse(inputFeed, 1);

        _validate();
    }

    // Auxillary function that prints usage help message
    void printUsage()const;
};

// Helper functions

// Helper function to write to a variable. Throws on error
template< typename T >
struct VariableWrite {
    void operator()(T& var, const std::string& s)const {
        std::istringstream iss(s);
        iss >> var;
        if(iss.fail())
            throw ParsingError("Cannot understand argument value " + s);
    }
};

// Helper function to append to a vector. Throws on error
template< typename T >
struct VectorAppend {
    void operator()(std::vector<T>& var, const std::string& s)const {
        T tmp;
        std::istringstream iss(s);
        iss >> tmp;
        if(iss.fail())
            throw ParsingError("Cannot append the argument value " + s);
        var.emplace_back(std::move(tmp));
    }
};

} // namespace cmdparse

#endif // include guard
