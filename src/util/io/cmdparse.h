#ifndef MEDYAN_UTIL_IO_CMDPARSE_H
#define MEDYAN_UTIL_IO_CMDPARSE_H

#include <cstdlib> // exit
#include <functional>
#include <iostream> // ostream, cout
#include <memory> // unique_ptr
#include <sstream> // string format
#include <stdexcept> // custom exception
#include <string>
#include <utility> // forward
#include <vector>

namespace cmdparse {

// Forward decl
template< typename T > struct VariableWrite;
template< typename T > struct VectorAppend;

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
    const std::string _description;
    const bool _required;
    const bool _list;

    struct State {
        size_t _occurenceCount = 0;
    } _state;

public:

    const std::function<void(const std::string&)> activate;

    PosArg(const char* name, const std::string& description, bool required, bool list, const std::function<void(const std::string&)>& activate) :
        _name(name), _description(description), _required(required), _list(list), activate(activate) {}

    bool isList() const { return _list; }
    bool isRequired()const { return _required; }

    const char* getName()const { return _name; }
    const std::string& getDescription()const { return _description; }

    void occur() { ++_state._occurenceCount; }
    size_t getOccurenceCount()const { return _state._occurenceCount; }

};

class Option {
private:

    const char _short = 0; // without "-". 0 for no short name
    const std::string _long; // without "--". "" for no long name
    const std::string _description;

    const bool _hasVariable;
    const std::string _variableName; // Useless if _hasVariable is false

    const bool _required;

    struct State {
        size_t _occurenceCount = 0;
    } _state;
public:

    const std::function<void(const std::string&)> activate;

    Option(char shortName, const std::string& longName, const std::string& description, bool required, const std::function<void()>& activateWithoutVar) :
        _short(shortName), _long(longName), _description(description), _hasVariable(false), _variableName(), _required(required),
        activate([activateWithoutVar](const std::string&) { activateWithoutVar(); }) {}
    Option(char shortName, const std::string& longName, const std::string& variableName, const std::string& description, bool required, const std::function<void(const std::string&)>& activate):
        _short(shortName), _long(longName), _description(description), _hasVariable(true), _variableName(variableName), _required(required), activate(activate) {}

    char getShortName()const { return _short; }
    const std::string& getLongName()const { return _long; }
    std::string getReadableName()const;
    std::string getUsageName()const;
    const std::string& getDescription()const { return _description; }
    const std::string& getVariableName()const { return _variableName; }

    bool hasVariable()const { return _hasVariable; }
    bool isRequired()const { return _required; }

    void occur() { ++_state._occurenceCount; }
    size_t getOccurenceCount()const { return _state._occurenceCount; }
};


class Command {
private:

    const std::string _name;
    const std::string _inheritedName;
    const std::string _description;

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
    // Internal validation. Recursive check for invoked subcommands.
    void _validate()const;

    // State variable
    struct State {
        bool _parseAsPosArg = false; // After the delimiter, this will be set to true.
        size_t _posArgCount = 0; // Number of positional argument encountered.
        size_t _posArgIndex = 0; // The index for the next positional argument.
                                 // Normally same with _posArgCount except on arg list.
        size_t _occurenceCount = 0;
    } _state;

    // Private constructor, allowing parent assignment and name inheriting
    Command(Command* parent, const std::string& inheritedName, const std::string& name, const std::string& description, const std::function<void()>& activate) :
        _name(name), _inheritedName(inheritedName), _description(description), _parent(parent), activate(activate) {}

public:

    const std::function<void()> activate;

    // Public constructor, where parent and inheritedName are default
    Command(const std::string& name, const std::string& description, const std::function<void()>& activate) :
        _name(name), _description(description), activate(activate) {}
    Command(const std::string& name, const std::string& description) :
        _name(name), _description(description) {}

    const std::string& getName()const { return _name; }
    std::string getFullName()const { return _inheritedName + ' ' + _name; }
    const std::string& getDescription()const { return _description; }

    void occur() { ++_state._occurenceCount; }

    // Specification
    template< typename... Args >
    PosArg* addPosArg(Args&&... args) {
        _posArgs.emplace_back(new PosArg(std::forward<Args>(args)...));
        return _posArgs.back().get();
    }
    template< typename T >
    PosArg* addPosArgForVar(const char* name, const std::string& description, bool required, T& var) {
        // Implies non-list
        return addPosArg(name, description, required, false, [name, &var](const std::string& arg) {
            VariableWrite<T>(std::string(name))(var, arg);
        });
    }
    template< typename T >
    PosArg* addPosArgForVector(const char* name, const std::string& description, bool required, std::vector<T>& var) {
        // Implies list
        return addPosArg(name, description, required, true, [name, &var](const std::string& arg) {
            VectorAppend<T>(std::string(name))(var, arg);
        });
    }
    template< typename... Args >
    Option* addOption(Args&&... args) {
        _options.emplace_back(new Option(std::forward<Args>(args)...));
        return _options.back().get();
    }
    Option* addOptionAsFlag(char shortName, const std::string& longName, const std::string& description, bool required) {
        // Add an option without argument
        return addOption(shortName, longName, description, required, []{});
    }
    template< typename T >
    Option* addOptionWithVar(char shortName, const std::string& longName, const std::string& variableName, const std::string& description, bool required, T& var) {
        return addOption(shortName, longName, variableName, description, required, [variableName, &var](const std::string& arg) {
            VariableWrite<T>{variableName}(var, arg);
            // Caveat: if () is used instead of {}, then T(a)(b, c) will be considered as the definition of variable a with type T initialized with (b, c).
        });
    }
    Option* addHelp() {
        // Auxilliary for
        return addOption('h', "help", "Print usage and exit", false, [this] {
            printUsage();
            std::exit(EXIT_SUCCESS);
        });
    }
    template< typename... Args >
    Command* addCommand(Args&&... args) {
        _subcommands.emplace_back(new Command(
            this,
            (_inheritedName.length() ? _inheritedName + ' ' + _name : _name),
            std::forward<Args>(args)...
        ));
        return _subcommands.back().get();
    }
    Command* addCommandSimple(const std::string& name, const std::string& description) {
        return addCommand(name, description, []{});
    }

    void setValidation(const std::function<void()>& validation) { _userValidation = validation; }

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

    // Auxilliary function that prints usage help message
    void printUsage(std::ostream& os = std::cout)const;
};

// Helper functions

// Helper function to write to a variable. Throws on error
template< typename T >
struct VariableWrite {
    std::string argName;
    VariableWrite() = default;
    VariableWrite(const std::string& argName): argName(argName) {}

    void operator()(T& var, const std::string& s)const {
        std::istringstream iss(s);
        iss >> var;
        if(iss.fail())
            throw ParsingError("Cannot understand argument value " + s + (argName.length() ? " for " + argName : std::string{}));
    }
};

// Helper function to append to a vector. Throws on error
template< typename T >
struct VectorAppend {
    std::string argName;
    VectorAppend() = default;
    VectorAppend(const std::string& argName): argName(argName) {}

    void operator()(std::vector<T>& var, const std::string& s)const {
        T tmp;
        std::istringstream iss(s);
        iss >> tmp;
        if(iss.fail())
            throw ParsingError("Cannot append the argument value " + s + (argName.length() ? " for " + argName : std::string{}));
        var.emplace_back(std::move(tmp));
    }
};

} // namespace cmdparse

#endif // include guard
