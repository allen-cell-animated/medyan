#ifndef MEDYAN_CORE_COMMAND_LINE_PARSER_H
#define MEDYAN_CORE_COMMAND_LINE_PARSER_H

#include <functional>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace medyan {

// Forward declarations
class CommandLineOptionBase;

class CommandLineParser {
private:
    std::vector<std::unique_ptr<CommandLineOptionBase>> _op;

    /// Fail flags and associated variables
    bool _invalidOptionBit = false;
    std::string _invalidOptionContent;
    bool _parseErrorBit = false;
    bool _unusedArgBit = false;
    std::vector<std::string> _unusedArgs;

public:
    enum class ArgType {
        Short, // Short option (eg "-a") or combined short option (eg "-af"). No more '-' is allowed
        Long, // Long option (eg "--analysis").
        Argument, // Argument or sub-command, not starting with "-" or "--" (eg "23" or "run")
        Fail // Unsupported syntax (eg "-a-f")
    };

    operator bool()const {
        return !(_invalidOptionBit || _parseErrorBit || _unusedArgBit);
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

    bool parse(int argc, char** argv);

    /// Print error message
    void printError()const {
        if(_invalidOptionBit)
            std::cout << "Invalid option: " << _invalidOptionContent << std::endl;
        if(_parseErrorBit)
            ;// do nothing, because the message is printed by the option class
        if(_unusedArgBit) {
            std::cout << "Unknown option:" << std::endl;
            for(auto& eachUnusedArg: _unusedArgs) std::cout << eachUnusedArg << " ";
            std::cout << std::endl;
        }
    }
};

class CommandLineOptionBase {
protected:
    /// Input of the option
    const char* _match;
    std::string _description;

    CommandLineParser* _p;

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
    bool _required = false;
    std::string _requiredNotFulfilled {""};
    int _numArgs = 0;

    /// Evaluate
    bool _evaluated = false;

    /// Constructors
    CommandLineOptionBase(const std::string& description):
        _description(description) {}
    CommandLineOptionBase(const char* match, const std::string& description):
        _match(match), _description(description) {}

public:
    /// Check state
    virtual operator bool()const {
        return !(_inputFailBit || _endOfArgListBit || _invalidArgBit || _validationFailBit);
    }

    /// Getters
    int getRequiredArgs()const { return _numArgs; }

    bool getInputFailBit()const { return _inputFailBit; }
    bool getEndOfArgListBit()const { return _endOfArgListBit; }
    bool getInvalidArgBit()const { return _invalidArgBit; }
    const std::string& getInvalidArgContent()const { return _invalidArgContent; }

    char getFlagShort()const { return _flagShort; }

    /// Change the option behavior. These manipulators always return the class itself.
    virtual CommandLineOptionBase& required(const std::string& errorMessage) {
        _required = true;
        _requiredNotFulfilled = errorMessage;
        return *this;
    }

    /// Find hit.
    virtual bool findHit(const std::string& arg, CommandLineParser::ArgType argType);

    /// Evaluate. return how many arguments consumed.
    virtual int evaluate(int argc, char** argv, int argp) = 0;

    /// Validate. Options with arguments should override this function
    virtual bool validate() { return true; }

    /// Print error message
    virtual void printError()const {
        if(_inputFailBit)
            std::cout << "Internal error: unrecognized flag " << _match << std::endl;
        if(_endOfArgListBit)
            std::cout << "Must specify argument for " << _match << std::endl;
        if(_invalidArgBit)
            std::cout << "Incorrect argument for " << _match << ": " << _invalidArgContent << std::endl;
        if(_validationFailBit)
            std::cout << "The argument value for " << _match << ": " << _validationFailContent << " is not valid." << std::endl;
    }
};

/// Command line option which takes exactly one argument
template<typename T>
class CommandLineOption: public CommandLineOptionBase {
private:
    /// The argument value
    T* _value;

    /// Validator
    std::function<bool(T)> _validator;

public:
    CommandLineOption(const char* match, T* destination, const std::string& description):
        CommandLineOptionBase(match, description), _value(destination),
        _numArgs(1)
    {
        _preprocess();
    }
    CommandLineOption(const char* match, T* destination, std::function<bool(T)> validator, const std::string& description):
        _match(match), _value(destination), _description(description), _validator(validator),
        _numArgs(1)
    {
        _preprocess();
    }

    /// Evaluate
    virtual int evaluate(int argc, char** argv, int argp)override {
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

        return 1;
    }

    /// Validate.
    virtual bool validate()override {
        bool pass = !_validator || _validator(*_value);
        if(!pass) {
            std::ostringstream oss(_validationFailContent);
            oss << _value;
        }
        return pass;
    }

};

/// Command line option which takes no arguments
class CommandLineFlag: public CommandLineOptionBase {
    /// The flag value
    bool* _value;
    bool _turnOn = true; ///< _value is true or false when this flag exists.

public:
    CommandLineFlag(const char* match, bool* destination, bool turnOn, const std::string& description):
        CommandLineOptionBase(match, description), _value(destination), _turnOn(turnOn)
    {
        _preprocess();
    }

    /// Evaluate
    virtual int evaluate(int argc, char** argv, int argp)override { *_value = _turnOn; return 0; }
};


} // namespace medyan

#endif
