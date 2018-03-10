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
 * should be captured in this step. Parsing with multiple possibilities is a
 * greedy process.
 * 
 * Evaluating is where the value interpretation and validation happens. If any
 * value is invalid, the whole command line reading should fail, without trying
 * to parse in any other ways.
 */

class CommandLineElement {
protected:
    std::string _description;
    CommandLineElement(const std::string& description): _description(description) {}

    /// Fail bits
    bool _inputFailBit = false; // Interal error: cannot parse input
    std::vector<std::string> _inputFailInfo;

    /// Preprocessing
    virtual void _preprocess() {};

    /// Fields set by Parsing
    bool _offered; ///< Whether the element is offered in input, determined by Parsing

    /// Configuration
    bool _required;

    /// Fields set by Evaluating
    bool _evaluated;

public:
    /// Check state
    operator bool()const {
        return !_inputFailBit;
    }

    /// Getters
    const std::string& getDescription()const { return _description; }

    bool inputFail()const { return _inputFailBit; }

    bool isOffered()const { return _offered; }

    bool isRequired()const { return _required; }

    bool isEvaluated()const { return _evaluated; }

    /// Modifier
    virtual CommandLineElement& require(bool required=true) { _required = required; return *this; }
    virtual CommandLineElement& offer(bool offered=true) { _offered = offered; return *this; }

    /// Evaluate. Returns whether successful
    virtual bool evaluate() = 0;

    /// Error printing
    virtual void printError(std::ostream& os=std::cout)const {
        if(_inputFailBit) {
            for(auto& info: _inputFailInfo)
                os << info << std::endl;
        }
    }
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
    virtual void _preprocess()override;

    /// Fail flags and associated info
    bool _endOfArgListBit = false;

    bool _activateErrorBit = false; // activation fail by the callback
    bool _invalidArgBit = false; // Invalid argument
    std::string _invalidArgInfo;

    /// Configuration of the option
    std::vector<OptionBase*> _excluding;

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

    bool endOfArgList()const { return _endOfArgListBit; }
    bool invalidArg()const { return _invalidArgBit; }
    const std::string& getInvalidArgInfo()const { return _invalidArgInfo; }

    char getFlagShort()const { return _flagShort; }
    const std::string& getFlagLong()const { return _flagLong; }
    const std::string& getFlagCommand()const { return _flagCommand; }

    const std::vector<OptionBase*>& getExcluding()const { return _excluding; }

    /// Modify configuration
    virtual OptionBase& exclude(OptionBase* op) { _excluding.push_back(op); return *this; }
    virtual OptionBase& fillField(const std::string& field) { return *this; }

    /// Find hit.
    virtual bool findHit2(const std::string& arg, ArgType argType);
    virtual bool findHit(char argShort);
    virtual bool findHit(const std::string& argLong);

    /// Evaluate and validate. return how many arguments consumed.
    virtual int evaluate2(int argc, char** argv, int argp) = 0;

    /// Print error message
    virtual void printError(std::ostream& os=std::cout)const {
        CommandLineElement::printError(os);

        if(_endOfArgListBit)
            os << "Must specify argument for " << _match << std::endl;
        if(_invalidArgBit)
            os << "Invalid argument for " << _match << ": " << _invalidArgInfo << std::endl;
        if(_activateErrorBit)
            os << "The argument value for " << _match << " is not acceptable." << std::endl;
    }
};

/// Command line option which takes exactly one argument
template<typename T>
class Option1: public OptionBase {
private:
    /// The argument value
    T _value;

    /// Fields supplied by Parsing
    std::string _field;

public:
    Option1(const std::string& description, const char* match, const std::string& argName, T* destination):
        OptionBase(description, match, true, argName, [destination, this]()->bool{ *destination = _value; return true; }) {}

    /// Modifiers
    virtual Option1& fillField(const std::string& field) { _field = field; return *this; }

    /// Evaluate
    virtual bool evaluate()override {
        _evaluated = true;

        std::istringstream iss(_field);
        iss >> _value;
        if(iss.fail()) {
            _invalidArgBit = true;
            _invalidArgInfo = _field;
        }

        if(!_activate || !_activate()) _activateErrorBit = true;

        return operator bool();
    }
    /// Evaluate and activate
    virtual int evaluate2(int argc, char** argv, int argp)override {
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
            _invalidArgInfo = std::string(argv[argp]);
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
    virtual bool evaluate()override {
        if(!_activate()) _activateErrorBit = true;
        return operator bool();
    }
    virtual int evaluate2(int argc, char** argv, int argp)override {
        if(!_activate()) _activateErrorBit = true;
        return 0;
    }
};


/// Base of positional element
class PosElement: public CommandLineElement {
protected:
    const bool _command; ///< Whether this is a command
    const bool _argument; ///< Whether this is an argument

    PosElement(const std::string& description, bool isCommand, bool isArgument):
        CommandLineElement(description), _command(isCommand), _argument(isArgument) {}
    
public:

    /// Getters
    bool isCommand()const { return _command; }
    bool isArgument()const { return _argument; }
    virtual const char* getCommandName()const { return ""; }

    /// Modifier
    virtual PosElement& require(bool required=true) { _required = required; return *this; }
    virtual PosElement& fillField(const std::string& field) { return *this; }

    /// Main parsing function.
    /// Returns the index just after the last argument it is using, or -1 if anything is wrong.
    virtual int parse(int argc, char** argv, int argp=0) = 0;
};

// Note that although argument taken by options are also positional, they are
// handled by the option class and thus not belong here.
template<typename T>
class PosArg: public PosElement {
    /// Fields set by Parsing
    std::string _field;

    /// Input
    std::string _argName; ///< Name of the argument

    /// Argument value
    T _value;

    /// Activate callback
    std::function<bool()> _activate;

    /// Fail bit
    bool _invalidArgBit = false; // Invalid argument
    std::string _invalidArgInfo;
    bool _activateErrorBit = false;

public:
    PosArg(const std::string& description, const std::string& argName, T* destination):
        PosElement(description, false, true), _argName(argName),
        _activate([destination, this]()->bool{ *destination = _value; return true; }) {}

    /// Check state
    operator bool()const {
        return !(_invalidArgBit || _activateErrorBit);
    }

    /// Modifier
    virtual PosArg& fillField(const std::string& field) { _field = field; return *this; }

    /// Dummy parsing. This should never be called.
    virtual int parse(int argc, char** argv, int argp=0)override { return -1; }

    /// Evaluate
    virtual void evaluate()override {
        _evaluated = true;

        std::istringstream iss(_field);
        iss >> _value;
        if(iss.fail()) {
            _invalidArgBit = true;
            _invalidArgInfo = _field;
        }

        if(!_activate || !_activate()) _activateErrorBit = true;

        return operator bool();
    }

    /// Print error
    virtual void printError(std::ostream& os=std::cout)override {
        PosElement::printError();

        if(_invalidArgBit) os << "Invalid argument for " << _argName << std::endl;
        if(_activateErrorBit) os << "Argument " << _argName << " has incorrect value: " << _value << std::endl;
    }
};

class Command: public PosElement {
private:
    std::vector<PosElement*> _pos;
    std::vector<OptionBase*> _op;
    std::vector<Command*> _subcmd;
    std::function<bool()> _activate; ///< Activate callback

    /// Name for the subcommand
    const char* _name;

    /// Preprocessing
    virtual void _preprocess()override;

    /// Fail flags and associated variables
    bool _parseErrorBit = false; // Syntax error in parsing. Should abort parsing when this is set to true.
    bool _logicErrorBit = false; // Option requirements not fulfilled
    std::vector<std::string> _logicErrorInfo;

    bool _activateErrorBit = false; // Activate error
    bool _subFailBit = false; // Children command/option fail
    std::vector<std::string> _subFailInfo;

    /// States
    bool _evaluated = false;

public:

    /// Constructor
    Command(const std::string& description, const char* name, const std::vector<OptionBase*>& ops, const std::vector<Command*>& subcmds, const std::function<bool()>& activate) :
        PosElement(description, true, false), _name(name), _op(ops), _subcmd(subcmds), _activate(activate) {}

    /// Check state
    operator bool()const {
        return !(_parseErrorBit || _logicErrorBit || _activateErrorBit || _subFailBit);
    }

    /// Getters
    virtual const char* getCommandName()const override { return _name; }

    /// Main parsing function
    virtual int parse(int argc, char** argv, int argp = 0)override;
    int parse2(int argc, char** argv, int argp=0);

    /// Evaluate
    virtual bool evaluate()override {
        _evaluated = true;

        if(!_activate || !_activate()) _activateErrorBit = true;
        for(auto pe: _pos) {
            if(pe->isRequired() || pe->isOffered())
                if(!pe->evaluate()) {
                    _subFailBit = true;
                    std::ostringstream oss;
                    pe->printError(oss);
                    _subFailInfo.emplace_back(oss.str());
                    return false;
                }
        }
        for(auto ob: _op) {
            if(ob->isRequired() || ob->isOffered())
                if(!ob->evaluate()) {
                    _subFailBit = true;
                    std::ostringstream oss;
                    pe->printError(oss);
                    _subFailInfo.emplace_back(oss.str());
                    return false;
                }
        }

        return operator bool();
    }

    /// Print message
    void printUsage(std::ostream& out = std::cout)const;
    virtual void printError(std::ostream& os = std::cout)const override {
        PosElement::printError(os);

        if(_parseErrorBit)
            for(auto& opPtr : _op) opPtr->printError(os);
        if(_logicErrorBit)
            for(auto& info : _logicErrorInfo) os << info << std::endl;
        if(_subFailBit)
            for(auto& info: _subFailInfo) os << info << std::endl;
    }
};


/// Main parsing function
inline bool commandLineParse(int argc, char** argv, Command& cmd) {
    if(!cmd.parse(argc, argv)) {
        cmd.printError();
        return false;
    }
    if(!cmd.evaluate()) {
        cmd.printError();
        return false;
    }
    return true;
}


} // namespace commandline
} // namespace medyan

#endif
