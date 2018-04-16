#ifndef MEDYAN_CORE_IO_COMMAND_LINE_OPT_H
#define MEDYAN_CORE_IO_COMMAND_LINE_OPT_H

#include <functional>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_set>
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
    bool _offered = false; ///< Whether the element is offered in input, determined by Parsing

    /// Configuration
    bool _required = false;

    /// Fields set by Evaluating
    bool _evaluated = false;

public:
    /// Check state
    virtual operator bool()const {
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
                os << "Internal error: " << info << std::endl;
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
    virtual void _preprocess()override;

    /// Fail flags and associated info
    bool _activateErrorBit = false; // activation fail by the callback
    bool _invalidArgBit = false; // Invalid argument
    std::string _invalidArgInfo;

    /// Activate callback
    std::function<bool()> _activate;

    /// Constructors
    OptionBase(const std::string& description, const char* match, bool takesArg, const std::string& argName, const std::function<bool()>& activate):
        CommandLineElement(description), _match(match), _takesArg(takesArg), _argName(argName), _activate(activate) {

        _preprocess();
    }

public:
    /// Check state
    virtual operator bool()const override {
        return CommandLineElement::operator bool() &&
            !(_activateErrorBit || _invalidArgBit);
    }

    /// Getters
    const char* getMatch()const { return _match; }

    bool takesArg()const { return _takesArg; }
    const std::string& getArgName()const { return _argName; }

    bool invalidArg()const { return _invalidArgBit; }
    const std::string& getInvalidArgInfo()const { return _invalidArgInfo; }

    char getFlagShort()const { return _flagShort; }
    const std::string& getFlagLong()const { return _flagLong; }

    /// Modify configuration
    virtual OptionBase& fillField(const std::string& field) { return *this; }

    /// Find hit.
    virtual bool findHit(char argShort);
    virtual bool findHit(const std::string& argLong);

    /// Print error message
    virtual void printError(std::ostream& os=std::cout)const override {
        CommandLineElement::printError(os);

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
    Option1(const std::string& description, const char* match, const std::string& argName, T* destination, const std::function<bool()>& activate):
        OptionBase(description, match, true, argName, [destination, activate, this]()->bool{ *destination = _value; return activate(); }) {}

    /// Modifiers
    virtual Option1& fillField(const std::string& field)override { _field = field; return *this; }

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
};


/// Base of positional element
class PosElement: public CommandLineElement {
public:
    enum class PosType {
        Command, Arg, Holder, MutuallyExclusive
    };

protected:
    const PosType _posType; ///< The element type

    /// Fail bits
    bool _parseErrorBit = false; ///< Syntax error in parsing. Should abort parsing when this is set to true.

    PosElement(const std::string& description, PosType posType):
        CommandLineElement(description), _posType(posType) {}

public:

    /// Check state
    virtual operator bool()const override {
        return CommandLineElement::operator bool() &&
            !_parseErrorBit;
    }

    /// Getters
    PosType getPosType()const { return _posType; }
    bool isCommand()const { return _posType == PosType::Command; }
    bool isArgument()const { return _posType == PosType::Arg; }
    bool isHolder()const { return _posType == PosType::Holder; }
    bool isMutuallyExclusive()const { return _posType == PosType::MutuallyExclusive; }

    virtual const char* getCommandName()const { return ""; }

    /// Modifier
    virtual PosElement& require(bool required=true)override { _required = required; return *this; }
    virtual PosElement& fillField(const std::string& field) { return *this; }

    /// Helper function used by parent. This function deals with requiredness and offering
    virtual int parseThis(int argc, char** argv, int argp) {
        int newArgp = -1;

        // Un-offer this
        offer(false);

        if(_required) {
            // We don't care about the offer for required positional elements.
            newArgp = parse(argc, argv, argp);
        } else { // If it is optional
            bool offerThisTime = true;
            do {
                offer(offerThisTime);

                newArgp = parse(argc, argv, argp);
                if(newArgp >= 0) break; // else we wait

                offerThisTime = !offerThisTime; // false for the 2nd time
            } while(!offerThisTime && newArgp < 0); // This will exit for at most 2 iterates.

            // if newArgp < 0 here, it must be in the un-offered state.
        }

        return newArgp;
    }
    /// Main parsing function.
    /// Returns the index just after the last argument it is using, or -1 if anything is wrong.
    virtual int parse(int argc, char** argv, int argp=0) = 0;

    /// Helper function for printUsage
    virtual void printContent(std::ostream& os=std::cout)const = 0;
    virtual void addToCmdOpSet(std::unordered_set<const PosElement*>& pos, std::unordered_set<const OptionBase*>& op)const = 0;
    
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
        PosElement(description, PosType::Arg), _argName(argName),
        _activate([destination, this]()->bool{ *destination = _value; return true; }) {}

    /// Check state
    virtual operator bool()const override {
        return PosElement::operator bool() &&
            !(_invalidArgBit || _activateErrorBit);
    }

    /// Modifier
    virtual PosArg& fillField(const std::string& field)override { _field = field; return *this; }

    virtual int parse(int argc, char** argv, int argp=0)override {
        if(argp >= argc) {
            _parseErrorBit = true;
            return -1;
        }

        _parseErrorBit = false;

        // Fill itself
        fillField(std::string(argv[argp]));

        return argp + 1;
    }

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

    /// Helper function for printUsage
    virtual void printContent(std::ostream& os=std::cout)const override {
        os << '<' << _argName << '>';
    }
    virtual void addToCmdOpSet(std::unordered_set<const PosElement*>& pos, std::unordered_set<const OptionBase*>& op)const {}

    /// Print error
    virtual void printError(std::ostream& os=std::cout)const override {
        PosElement::printError();

        if(_invalidArgBit) os << "Invalid argument for " << _argName << std::endl;
        if(_activateErrorBit) os << "Argument " << _argName << " has incorrect value: " << _value << std::endl;
    }
};

/// Class PosHolder is a container of PosElements and Options
class PosHolder: public PosElement {
private:
    std::vector<PosElement*> _pos;
    std::vector<OptionBase*> _op;

    /// Fail flags and associated variables
    bool _logicErrorBit = false; // Option requirements not fulfilled
    std::vector<std::string> _logicErrorInfo;

    bool _subFailBit = false; // Children command/option fail
    std::vector<std::string> _subFailInfo;

public:

    /// Constructor
    PosHolder(const std::vector<OptionBase*>& op, const std::vector<PosElement*>& pos) :
        PosElement("", PosType::Holder), _op(op), _pos(pos) {}

    /// Check state
    virtual operator bool()const override {
        return PosElement::operator bool() &&
            !(_logicErrorBit || _subFailBit);
    }

    /// Getters
    const std::vector<PosElement*>& getPos()const { return _pos; }
    const std::vector<OptionBase*>& getOp()const { return _op; }

    /// Main parsing function
    virtual int parse(int argc, char** argv, int argp=0)override;

    /// Evaluate
    virtual bool evaluate()override {
        _evaluated = true;

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
                    ob->printError(oss);
                    _subFailInfo.emplace_back(oss.str());
                    return false;
                }
        }

        return operator bool();
    }

    /// Helper function for printUsage.
    virtual void printContent(std::ostream& os=std::cout)const override;
    virtual void addToCmdOpSet(std::unordered_set<const PosElement*>& pos, std::unordered_set<const OptionBase*>& op)const {
        for(auto& pe: _pos) pe->addToCmdOpSet(pos, op);
        for(auto& ob: _op) op.insert(ob);
    }

    /// Print errors
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

class PosMutuallyExclusive: public PosElement {
    std::vector<PosElement*> _pos; ///< Store a list of mutually exclusive poses.

    /// Fields set by parsing
    int _index = -1; ///< which one to choose

    /// Fail bits
    bool _logicErrorBit = false; // Parsing logic error
    std::vector<std::string> _logicErrorInfo;

    bool _unknownBit = false; // Unknown internal error
    bool _subFailBit = false; // Children command/option fail
    std::vector<std::string> _subFailInfo;

    /// Preprocessing. This class only accepts required elements
    virtual void _preprocess()override;

public:
    PosMutuallyExclusive(const std::vector<PosElement*>& pos):
        PosElement("", PosType::MutuallyExclusive), _pos(pos) {}

    /// Check state
    virtual operator bool()const override {
        return PosElement::operator bool() &&
            !(_logicErrorBit || _unknownBit || _subFailBit);
    }

    /// Getters
    const std::vector<PosElement*>& getPos()const { return _pos; }

    /// Main parsing function
    virtual int parse(int argc, char** argv, int argp=0)override;

    /// Evaluate
    virtual bool evaluate()override {
        _evaluated = true;
        if(_index < 0 || _index >= _pos.size()) {
            _unknownBit = true;
            return false;
        }
        if(!_pos[_index]->evaluate()) {
            _subFailBit = true;
            std::ostringstream oss;
            _pos[_index]->printError(oss);
            _subFailInfo.emplace_back(oss.str());
            return false;
        }

        return operator bool();
    }

    /// Helper function for printUsage
    virtual void printContent(std::ostream& os=std::cout)const override {
        for(auto it = _pos.begin(); it != _pos.end(); ++it) {
            if(it != _pos.begin()) os << '|';
            (*it)->printContent(os);
        }
    }
    virtual void addToCmdOpSet(std::unordered_set<const PosElement*>& pos, std::unordered_set<const OptionBase*>& op)const {
        for(auto& pe: _pos) pe->addToCmdOpSet(pos, op);
    }

    /// Print error
    virtual void printError(std::ostream& os=std::cout)const override {
        PosElement::printError(os);

        if(_logicErrorBit) for(auto& info: _logicErrorInfo) os << info << std::endl;
        if(_unknownBit) os << "Something is wrong with command line mutual exclusive groups." << std::endl;
        if(_subFailBit) for(auto& info: _subFailInfo) os << info << std::endl;
    }
};

class Command: public PosElement {
private:
    PosElement* _content;
    std::function<bool()> _activate; ///< Activate callback

    /// Name for the subcommand
    const char* _name;

    /// Preprocessing
    virtual void _preprocess()override;

    /// Fail flags and associated variables
    bool _activateErrorBit = false; // Activate error
    bool _subFailBit = false; // Children command/option fail
    std::vector<std::string> _subFailInfo;

    /// Configurations
    bool _main = false;

public:

    /// Constructor
    Command(const std::string& description, const char* name, PosElement* content, const std::function<bool()>& activate) :
        PosElement(description, PosType::Command), _name(name), _content(content), _activate(activate) {}

    /// Check state
    virtual operator bool()const override {
        return PosElement::operator bool() &&
            !(_activateErrorBit || _subFailBit);
    }

    /// Getters
    virtual const char* getCommandName()const override { return _name; }

    /// Modifiers
    virtual Command& setMain(bool isMain=true) { _main = isMain; return *this; }

    /// Main parsing function
    virtual int parse(int argc, char** argv, int argp=0)override;

    /// Evaluate
    virtual bool evaluate()override {
        _evaluated = true;

        if(!_activate || !_activate()) _activateErrorBit = true;
        if(_content && !_content->evaluate()) {
            _subFailBit = true;
            std::ostringstream oss;
            _content->printError(oss);
            _subFailInfo.emplace_back(oss.str());
            return false;
        }

        return operator bool();
    }

    /// Helper function for printUsage
    virtual void printContent(std::ostream& os=std::cout)const override {
        os << _name;
    }
    virtual void addToCmdOpSet(std::unordered_set<const PosElement*>& pos, std::unordered_set<const OptionBase*>& op)const {
        pos.insert(this);
    }
    /// Print help message
    void printUsage(std::ostream& os=std::cout)const;

    /// Print errors
    virtual void printError(std::ostream& os = std::cout)const override {
        PosElement::printError(os);

        if(_parseErrorBit)
            if(_content) _content->printError(os);
        if(_subFailBit)
            for(auto& info: _subFailInfo) os << info << std::endl;
    }
};


/// Main parsing function
inline bool commandLineParse(int argc, char** argv, Command& cmd) {
    if(cmd.parse(argc, argv) < 0) {
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
