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

class PosArg {
private:
    const char * const _name;
    const bool _list;
    const bool _required;

    size_t _occurenceCount = 0;
public:
    bool isList() const { return _list; }
    bool isRequired()const { return _required; }

    size_t getOccurenceCount()const { return _occurenceCount; }
};

class Option {
protected:
    size_t _occurenceCount = 0;
public:
    size_t getOccurenceCount()const { return _occurenceCount; }
};
// Option that takes no variable
class Option0 : public Option {
};
// Option that takes 1 variable
template< typename T >
class Option1 : public Option {
};

class Command {
private:

    std::vector<std::unique_ptr<PosArg>> _posArgs;
    std::vector<std::unique_ptr<Option>> _options;
    std::vector<std::unique_ptr<Command>> _subcommands;

    std::function<void()> _userValidation;

    // Check validity of specified rules. Recursive check for subcommands.
    void _ruleCheck()const;
    // Real parsing function
    void _parse(int argc, char** argv);
    // Internal validation
    void _validation()const;

public:
    // The parsing function used only by the root command.
    // States will be changed in this function
    void parse(int argc, char** argv) {
        _ruleCheck();
        _parse(argc, argv);
        _validation();
    }
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
