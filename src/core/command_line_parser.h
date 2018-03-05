#ifndef MEDYAN_CORE_COMMAND_LINE_PARSER_H
#define MEDYAN_CORE_COMMAND_LINE_PARSER_H

#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace medyan {

class CommandLineOptionBase;
class CommandLineParser;

class CommandLineOptionBase {
protected:
    /// Input of the option
    const char* _match;
    std::string _description;

    CommandLineParser* _p;

    /// Fail bits
    bool _badinput = false;

    /// Configuration of the option
    bool _required = false;
    std::string _requiredNotFulfilled {""};
    int _requiredArgs = 0;

    /// Evaluate
    bool _evaluated = false;

public:
    /// Check state
    operator bool()const {
        return _badinput;
    }

    /// Getters
    int getRequiredArgs()const { return _requiredArgs; }

    /// Change the option behavior. These manipulators always return the class itself.
    virtual CommandLineOptionBase& required(const std::string& errorMessage) {
        _required = true;
        _requiredNotFulfilled = errorMessage;
        return *this;
    }

    /// Validate. Options with arguments should override this function
    virtual bool validate() { return true; }

    /// Find hit.
    virtual bool findHit(const std::string& arg, CommandLineParser::ArgType argType);
};

template<typename T>
class CommandLineOption: public CommandLineOptionBase {
private:
    /// The argument value
    T* _value;

    /// Validator
    std::function<bool(T)> _validator;

public:
    CommandLineOption(const char* match, T* destination, const std::string& description):
        _match(match), _value(destination), _description(description),
        _requiredArguments(1) {}
    CommandLineOption(const char* match, T* destination, std::function<bool(T)> validator, const std::string& description):
        _match(match), _value(destination), _description(description), _validator(validator),
        _requiredArguments(1) {}

    /// Validate.
    virtual bool validate()override { return !_validator || _validator(*_value); }

};

class CommandLineFlag: public CommandLineOptionBase {
    /// The flag value
    bool* _value;
public:
    CommandLineFlag(const char* match, bool* destination, const std::string& description);
};

/* Broken type... Do not use it for now. */
class CommandLineGroup: public CommandLineOptionBase {
    /// Pointers to other options
    std::vector<std::unique_ptr<CommandLineOptionBase>> _op;
    /// Aggregate type
public:
    CommandLineGroup(CommandLineParser* p): _p(p) {}

    /// Validate.
    virtual bool validate()override {
        for(auto& opPtr: _op) {
            if(!opPtr->validate()) return false;
        }
        return true;
    }
};

class CommandLineParser {
public:
    enum class ArgType {
        Short, // Short option (eg "-a") or combined short option (eg "-af"). No more '-' is allowed
        Long, // Long option (eg "--analysis").
        Argument, // Argument or sub-command, not starting with "-" or "--" (eg "23" or "run")
        Fail // Unsupported syntax (eg "-a-f")
    };

private:
    std::vector<std::unique_ptr<CommandLineOptionBase>> _op;

    /// Helper function to get arg type
    static ArgType getArgType(char* argv) {
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

public:
    void parse(int argc, char** argv);
};

} // namespace medyan

#endif