#ifndef MEDYAN_CORE_COMMAND_LINE_PARSER_H
#define MEDYAN_CORE_COMMAND_LINE_PARSER_H

#include <functional>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace medyan {
namespace commandline {

enum class ArgType {
    Short, // Short option (eg "-a") or combined short option (eg "-af"). No more '-' is allowed
    Long, // Long option (eg "--analysis").
    ArgOrCmd, // Argument or sub-command, not starting with "-" or "--" (eg "23" or "run")
};

/// Helper function to get arg type
inline ArgType getArgType(const char* argv) {
    int i = 0;
    bool startHyphen = false;
    bool startHyphen2 = false;
    while (argv[i]) {
        if (i == 0 && argv[i] == '-') startHyphen = true;
        if (i == 1 && startHyphen && argv[i] == '-') startHyphen2 = true;
        ++i;
    }
    if (startHyphen2) return ArgType::Long;
    if (startHyphen) return ArgType::Short;
    return ArgType::ArgOrCmd;
}

/**
 * The actual handling of the command line input should be separated into two
 * parts: Parsing and Evaluating.
 * 
 * Parsing is to determine whether the user input is syntactically correct and
 * to select from the multiple groupings. Actual value interpretation and
 * validation do not occur here. Internal errors and command line syntax errors
 * should be captured in this step.
 * 
 * Evaluating is where the value interpretation and validation happens. If any
 * value is invalid, the whole command line reading should fail, without trying
 * to parse in any other ways.
 */

class CommandLineElement {
protected:
    std::string _description;
    CommandLineElement(const std::string& description): _description(description) {}
public:
    const std::string& getDescription()const { return _description; }
};
class OptionBase: public CommandLineElement {
protected:
    /// Input of the option
    const char* _match;
    const bool _takesArg = 0;
    std::string _argName;

    /// Preprocessing
    char _flagShort = 0; // e.g. "-a" w/o "-"
    std::string _flagLong; // e.g. "--analysis" w/o "--"
    std::string _flagCommand; // e.g. "add"
    void _preprocess();

    /// Fail flags and associated info
    bool _inputFailBit = false; // Internal: input parsing failure
    bool _endOfArgListBit = false;
    bool _activateErrorBit = false; // activation fail by the callback
    bool _invalidArgBit = false; // Invalid argument
    std::string _invalidArgContent;

    /// Configuration of the option
    bool _required = false;
    std::vector<OptionBase*> _excluding;

    /// State
    bool _evaluated = false;

    /// Activate callback
    std::function<bool()> _activate;

    /// Constructors
    OptionBase(const std::string& description, const char* match, bool takesArg, const std::string& argName, const std::function<bool()>& activate):
        CommandLineElement(description), _match(match), _takesArg(takesArg), _argName(argName), _activate(activate) {

        _preprocess();
    }

public:
    /// Check state
    operator bool()const {
        return !(_inputFailBit || _endOfArgListBit || _activateErrorBit || _invalidArgBit);
    }

    /// Getters
    const char* getMatch()const { return _match; }

    bool takesArg()const { return _takesArg; }
    const std::string& getArgName()const { return _argName; }

    bool inputFail()const { return _inputFailBit; }
    bool endOfArgList()const { return _endOfArgListBit; }
    bool invalidArg()const { return _invalidArgBit; }
    const std::string& getInvalidArgContent()const { return _invalidArgContent; }

    char getFlagShort()const { return _flagShort; }
    const std::string& getFlagLong()const { return _flagLong; }
    const std::string& getFlagCommand()const { return _flagCommand; }

    bool isRequired()const { return _required; }
    const std::vector<OptionBase*>& getExcluding()const { return _excluding; }

    bool isEvaluated()const { return _evaluated; }

    /// Modify configuration
    virtual OptionBase& require(bool required=true) { _required = required; return *this; }
    virtual OptionBase& exclude(OptionBase* op) { _excluding.push_back(op); return *this; }

    /// Find hit.
    virtual bool findHit(const std::string& arg, ArgType argType);

    /// Evaluate and validate. return how many arguments consumed.
    virtual int evaluate(int argc, char** argv, int argp) = 0;

    /// Print error message
    virtual void printError(std::ostream& out=std::cout)const {
        if(_inputFailBit)
            out << "Internal error: invalid configuration for " << _match << std::endl;
        if(_endOfArgListBit)
            out << "Must specify argument for " << _match << std::endl;
        if(_invalidArgBit)
            out << "Invalid argument for " << _match << ": " << _invalidArgContent << std::endl;
        if(_activateErrorBit)
            out << "The argument value for " << _match << " is not acceptable." << std::endl;
    }
};

/// Command line option which takes exactly one argument
template<typename T>
class Option1: public OptionBase {
private:
    /// The argument value
    T _value;

public:
    Option1(const std::string& description, const char* match, const std::string& argName, T* destination):
        OptionBase(description, match, true, argName, [destination, this]()->bool{ *destination = _value; return true; }) {}

    /// Evaluate and activate
    virtual int evaluate(int argc, char** argv, int argp)override {
        _evaluated = true;

        if(argp + 1 >= argc) {
            _endOfArgListBit = true;
            return 0;
        }

        ++argp;
        std::istringstream iss(argv[argp]);
        iss >> _value;
        if(iss.fail()) {
            _invalidArgBit = true;
            _invalidArgContent = std::string(argv[argp]);
        }

        if(!_activate()) _activateErrorBit = true;

        return 1;
    }

};

/// Command line option which takes no arguments
class Option0: public OptionBase {
public:
    Option0(const std::string& description, const char* match, bool* destination):
        OptionBase(description, match, false, "", [destination]()->bool{ *destination = true; return true; }) {}
    Option0(const std::string& description, const char* match, const std::function<bool()>& activate):
        OptionBase(description, match, false, "", activate) {}

    /// Evaluate
    virtual int evaluate(int argc, char** argv, int argp)override {
        if(!_activate()) _activateErrorBit = true;
        return 0;
    }
};

class PosElement: public CommandLineElement {};
// Note that although argument taken by options are also positional, they are
// handled by the option class and thus not belong here.
class PosArg: public PosElement {
    bool _required;
public:
    bool isRequired()const { return _required; }
};
class Command: public CommandLineElement {
private:
    std::vector<OptionBase*> _op;
    std::vector<Command*> _subcmd;
    std::function<bool()> _activate; ///< Activate callback

    /// Name for the subcommand
    std::string _name;

    /// Fail flags and associated variables
    bool _parseErrorBit = false; // Fail by option reader. Should abort parsing when this is set to true.
    bool _subcmdErrorBit = false; // Fail by subcommand parser. Should abort parsing when this is set to true.
    bool _unusedArgBit = false; // Unused arguments after parsing
    std::vector<std::string> _unusedArgs;
    bool _logicErrorBit = false; // Option requirements not fulfilled
    std::vector<std::string> _logicErrorContent;
    bool _activateErrorBit = false; // Activate error

    /// States
    bool _evaluated = false;

public:

    /// Constructor
    Command(const std::string& description, const std::string& name, const std::vector<OptionBase*>& ops, const std::vector<Command*>& subcmds, const std::function<bool()>& activate) :
        CommandLineElement(description), _name(name), _op(ops), _subcmd(subcmds), _activate(activate) {}

    /// Check state
    operator bool()const {
        return !(_parseErrorBit || _subcmdErrorBit || _unusedArgBit || _logicErrorBit || _activateErrorBit);
    }

    /// Getters
    bool isEvaluated()const { return _evaluated; }

    /// Main parsing function
    bool parse(int argc, char** argv, int argp = 0);

    /// Evaluate
    int evaluate() {
        if(!_activate || !_activate()) _activateErrorBit = true;
        return 1;
    }

    /// Print message
    void printUsage(std::ostream& out = std::cout)const;
    void printError(std::ostream& out = std::cout)const {
        if (_parseErrorBit) {
            for (auto& opPtr : _op) opPtr->printError(out);
        }
        if (_subcmdErrorBit)
            for (auto& cmdPtr : _subcmd) cmdPtr->printError(out);
        if (_logicErrorBit)
            for (auto& content : _logicErrorContent) out << content << std::endl;
        if (_unusedArgBit) {
            out << "Unknown option:" << std::endl;
            for (auto& eachUnusedArg : _unusedArgs) out << eachUnusedArg << " ";
            out << std::endl;
        }
    }
};


} // namespace commandline
} // namespace medyan

#endif
