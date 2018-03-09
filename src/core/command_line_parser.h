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
    Argument, // Argument or sub-command, not starting with "-" or "--" (eg "23" or "run")
    Fail // Unsupported syntax (eg "-a-f")
};

/// Helper function to get arg type
inline ArgType getArgType(const char* argv) {
    int i = 0;
    bool startHyphen = false;
    bool startHyphen2 = false;
    while (argv[i]) {
        if (i == 0 && argv[i] == '-') startHyphen = true;
        if (i == 1 && startHyphen && argv[i] == '-') startHyphen2 = true;
        if (i > 1 && startHyphen && !startHyphen2 && argv[i] == '-') return ArgType::Fail;
        ++i;
    }
    if (startHyphen2) return ArgType::Long;
    if (startHyphen) return ArgType::Short;
    return ArgType::Argument;
}


class OptionBase {
protected:
    /// Input of the option
    const char* _match;
    std::string _description;
    const bool _takesArg = 0;

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
    OptionBase(const char* match, const std::string& description, bool takesArg, const std::function<bool()>& activate):
        _match(match), _description(description), _takesArg(takesArg), _activate(activate) {

        _preprocess();
    }

public:
    /// Check state
    operator bool()const {
        return !(_inputFailBit || _endOfArgListBit || _activateErrorBit || _invalidArgBit);
    }

    /// Getters
    const char* getMatch()const { return _match; }
    const std::string& getDescription()const { return _description; }

    bool takesArg()const { return _takesArg; }

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

    /// Activate after successful evaluation.
    virtual bool activate()const = 0; // TODO: const?
    /// Evaluate and validate. return how many arguments consumed.
    virtual int evaluate(int argc, char** argv, int argp) = 0;

    /// Print error message
    virtual void printError(std::ostream& out=std::cout)const {
        if(_inputFailBit)
            out << "Internal error: unrecognized option " << _match << std::endl;
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

    /// Default value
    const bool _hasDefault;
    T _defaultValue;

    /// name
    std::string _argName;

public:
    Option1(const char* match, const std::string& description, const std::string& argName, T* destination):
        OptionBase(match, description, true, [destination, this]()->bool{ *destination = _value; return true; }),
        _argName(argName), _hasDefault(false) {}
    Option1(const char* match, const std::string& description, const std::string& argName, T* destination, T defaultValue):
        OptionBase(match, description, true, [destination, this]()->bool{ *destination = _value; return true; }),
        _argName(argName), _hasDefault(true), _defaultValue(defaultValue) {}

    /// Getters
    const std::string& getArgName()const { return _argName; }

    /// Evaluate and activate
    virtual void evaluate(int argc, char** argv, int argp)override {
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
    Option0(const char* match, const std::string& description, bool* destination):
        OptionBase(match, description, false, [destination, this]()->bool{ *destination = true; return true; }) {}
    Option0(const char* match, const std::string& description, const std::function<bool()>& activate):
        OptionBase(match, description, false, activate) {}

    /// Evaluate
    virtual int evaluate(int argc, char** argv, int argp)override {
        if(!_activate()) _activateErrorBit = true;
        return 0;
    }
};

class Command {
private:
    std::vector<OptionBase*> _op;
    std::vector<Command*> _subcmd;
    OptionBase* _defaultOp = nullptr; ///< The command itself acting as an option

                                      /// Name for the subcommand
    std::string _name;

    /// Fail flags and associated variables
    bool _invalidArgBit = false; // Fail on understanding argument. Should abort parsing when this is set to true.
    std::string _invalidArgContent;
    bool _parseErrorBit = false; // Fail by option reader. Should abort parsing when this is set to true.
    bool _subcmdErrorBit = false; // Fail by subcommand parser. Should abort parsing when this is set to true.
    bool _unusedArgBit = false; // Unused arguments after parsing
    std::vector<std::string> _unusedArgs;
    bool _logicErrorBit = false; // Option requirements not fulfilled
    std::vector<std::string> _logicErrorContent;

    /// States
    bool _evaluated = false;

public:

    /// Constructor
    Command(OptionBase* op) : _defaultOp(op) { if (op) _name = op->getFlagCommand(); }
    Command(const std::string& name, const std::vector<OptionBase*>& ops, const std::vector<Command*>& subcmds) :
        _name(name), _op(ops), _subcmd(subcmds) {}

    /// Check state
    operator bool()const {
        return !(_invalidArgBit || _parseErrorBit || _subcmdErrorBit || _unusedArgBit || _logicErrorBit);
    }

    /// Getters
    OptionBase* getDefaultOp()const { return _defaultOp; }

    /// Main parsing function
    bool parse(int argc, char** argv, int argp = 0);

    /// Print message
    void printUsage(std::ostream& out = std::cout)const;
    void printError(std::ostream& out = std::cout)const {
        if (_invalidArgBit)
            out << "Invalid option: " << _invalidArgContent << std::endl;
        if (_parseErrorBit) {
            _defaultOp->printError(out);
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
