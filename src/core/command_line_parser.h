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

// Forward declarations
class OptionBase;

class Command {
    friend OptionBase;

private:
    std::vector<std::unique_ptr<OptionBase>> _op;
    std::vector<std::unique_ptr<Command>> _subcmd;
    std::unique_ptr<OptionBase> _defaultOp; ///< The command itself acting as an option

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
    enum class ArgType {
        Short, // Short option (eg "-a") or combined short option (eg "-af"). No more '-' is allowed
        Long, // Long option (eg "--analysis").
        Argument, // Argument or sub-command, not starting with "-" or "--" (eg "23" or "run")
        Fail // Unsupported syntax (eg "-a-f")
    };

    operator bool()const {
        return !(_invalidArgBit || _parseErrorBit || _subcmdErrorBit || _unusedArgBit || _logicErrorBit);
    }

    /// Helper function to get arg type
    static ArgType getArgType(const char* argv) {
        int i = 0;
        ArgType t;
        bool startHyphen = false;
        bool startHyphen2 = false;
        while(argv[i]) {
            if(i == 0 && argv[i] == '-') startHyphen = true;
            if(i == 1 && startHyphen && argv[i] == '-') startHyphen2 = true;
            if(i > 1 && startHyphen && !startHyphen2 && argv[i] == '-') return ArgType::Fail;
            ++i;
        }
        if(startHyphen2) return ArgType::Long;
        if(startHyphen) return ArgType::Short;
        return ArgType::Argument;
    }

    bool parse(int argc, char** argv, int argp=0);

    /// Print message
    void printUsage(std::ostream& out=std::cout)const;
    void printError(std::ostream& out=std::cout)const {
        if(_invalidArgBit)
            out << "Invalid option: " << _invalidArgContent << std::endl;
        if(_parseErrorBit) {
            _defaultOp->printError(out);
            for(auto& opPtr: _op) opPtr->printError(out);
        }
        if(_subcmdErrorBit)
            for(auto& cmdPtr: _subcmd) cmdPtr->printError(out);
        if(_logicErrorBit)
            for(auto& content: _logicErrorContent) out << content << std::endl;
        if(_unusedArgBit) {
            out << "Unknown option:" << std::endl;
            for(auto& eachUnusedArg: _unusedArgs) out << eachUnusedArg << " ";
            out << std::endl;
        }
    }
};

class OptionBase {
    friend Command;

protected:
    /// Input of the option
    const char* _match;
    std::string _description;

    /// Preprocessing
    char _flagShort = 0; // e.g. "-a" w/o "-"
    std::string _flagLong; // e.g. "--analysis" w/o "--"
    std::string _flagCommand; // e.g. "add"
    virtual void _preprocess();

    /// Fail flags and associated info
    bool _inputFailBit = false; // input parsing failure
    bool _endOfArgListBit = false; // end of argument list when reading arguments
    bool _invalidArgBit = false; // Invalid argument
    std::string _invalidArgContent;
    bool _validationFailBit = false; // Not passing the validation
    std::string _validationFailContent;

    /// Configuration of the option
    int _numArgs = 0;
    bool _required = false;
    std::vector<OptionBase*> _excluding;

    /// State
    bool _evaluated = false;

    /// Constructors
    OptionBase(const char* match, int numArgs, const std::string& description):
        _match(match), _numArgs(numArgs), _description(description) {}

public:
    /// Check state
    operator bool()const {
        return !(_inputFailBit || _endOfArgListBit || _invalidArgBit || _validationFailBit);
    }

    /// Getters
    const char* getMatch()const { return _match; }
    const std::string& getDescription()const { return _description; }

    int getRequiredArgs()const { return _numArgs; }

    bool inputFail()const { return _inputFailBit; }
    bool endOfArgList()const { return _endOfArgListBit; }
    bool invalidArg()const { return _invalidArgBit; }
    const std::string& getInvalidArgContent()const { return _invalidArgContent; }

    char getFlagShort()const { return _flagShort; }

    bool isRequired()const { return _required; }
    const std::vector<OptionBase*>& getExcluding()const { return _excluding; }

    bool isEvaluated()const { return _evaluated; }

    /// Modify configuration
    virtual OptionBase& require(bool required=true) { _required = required; return *this; }
    virtual OptionBase& exclude(OptionBase* op) { _excluding.push_back(op); return *this; }

    /// Find hit.
    virtual bool findHit(const std::string& arg, Command::ArgType argType);

    /// Evaluate and validate. return how many arguments consumed.
    virtual int evaluate(int argc, char** argv, int argp) = 0;

    /// Print error message
    virtual void printError(std::ostream& out=std::cout)const {
        if(_inputFailBit)
            out << "Internal error: unrecognized flag " << _match << std::endl;
        if(_endOfArgListBit)
            out << "Must specify argument for " << _match << std::endl;
        if(_invalidArgBit)
            out << "Incorrect argument for " << _match << ": " << _invalidArgContent << std::endl;
        if(_validationFailBit)
            out << "The argument value for " << _match << ": " << _validationFailContent << " is not valid." << std::endl;
    }
};

/// Command line option which takes exactly one argument
template<typename T>
class Option1: public OptionBase {
private:
    /// The argument value
    T* _value;

    /// Validator
    std::function<bool(T)> _validator;

public:
    Option1(const char* match, T* destination, const std::string& description):
        OptionBase(match, 1, description), _value(destination)
    {
        _preprocess();
    }
    Option1(const char* match, T* destination, std::function<bool(T)> validator, const std::string& description):
        OptionBase(match, 1, description), _value(destination), _validator(validator)
    {
        _preprocess();
    }

    /// Evaluate
    virtual int evaluate(int argc, char** argv, int argp)override {
        _evaluated = true;

        // _numArgs == 1
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

        if(_validator && !_validator(*_value)) {
            _validationFailBit = true;
            std::ostringstream oss;
            oss << _value;
            _validationFailContent = oss.str();
        }

        return 1;
    }

};

/// Command line option which takes no arguments
class Flag: public OptionBase {
    /// The flag value
    bool* _value;
    bool _turnOn = true; ///< _value is true or false when this flag exists.

public:
    Flag(const char* match, bool* destination, bool turnOn, const std::string& description):
        OptionBase(match, 0, description), _value(destination), _turnOn(turnOn)
    {
        _preprocess();
    }

    /// Evaluate
    virtual int evaluate(int argc, char** argv, int argp)override {
        _evaluated = true;
        *_value = _turnOn;
        return 0;
    }
};

/// Commmand line group logic


} // namespace commandline
} // namespace medyan

#endif
