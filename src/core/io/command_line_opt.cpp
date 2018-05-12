#include "core/io/command_line_opt.h"

#include <algorithm>
#include <cstring>

namespace medyan {
namespace commandline {

namespace {

    void usagePairFormatter(const std::string& key, const std::string& description, std::ostream& os) {
        static const size_t leftMargin = 2;
        static const size_t maxKeySize = 21;
        static const size_t midMargin = 1;

        os << std::string(leftMargin, ' ');

        size_t len = key.length();
        os << key;
        if(len > maxKeySize) os << '\n' << std::string(leftMargin + maxKeySize + midMargin, ' ');
        else os << std::string(maxKeySize + midMargin - len, ' ');

        os << description << std::endl;
    }
}

void OptionBase::_preprocess() {
    std::string eachStr;
    std::istringstream iss{std::string(_match)};
    while(std::getline(iss, eachStr, ',')) {
        size_t len = eachStr.length();
        if(len == 2 && eachStr[0] == '-' && eachStr[1] != '-') {
            // Short flag
            _flagShort = eachStr[1];
        } else if(len > 2 && std::string(eachStr, 0, 2) == "--") {
            // Long flag
            _flagLong = std::string(eachStr, 2);
        } else {
            // Failed
            _inputFailBit = true;
        }
    }
    
    // Callback error
    if(!_activate) {
        _inputFailBit = true;
        _inputFailInfo.emplace_back("Option activate callback not set.");
    }
}

bool OptionBase::findHit(char argShort)              { return _flagShort == argShort; }
bool OptionBase::findHit(const std::string& argLong) { return _flagLong  == argLong;  }

int PosHolder::parse(int argc, char** argv, int argp) {
    // Uncheck error bit
    _parseErrorBit = false;
    _logicErrorBit = false;
    _logicErrorInfo.clear();

    std::vector<OptionBase*> op = _op; // Copy
    for(auto ob: op) {
        if(!ob->isRequired()) ob->offer(false); // Un-offer optional options
    }
    auto itPos = _pos.begin();

    if(!*this) return -1;

    // Try to read to the end of the arguments
    bool shouldEnd = false;
    while(argp < argc && !shouldEnd) {
        ArgType t = getArgType(argv[argp]);
        std::string arg {argv[argp]};

        // parse arg
        switch(t) {
        case ArgType::Long:
            {
                size_t eq = arg.find('=');
                std::string flag(arg, 2, eq-2);
                auto it = std::find_if(op.begin(), op.end(),
                    [&flag](OptionBase* ob) { return ob->findHit(flag); });
                if(it == op.end()) { // Not found
                    shouldEnd = true;
                } else { // Found
                    if((*it)->takesArg()) {
                        if(eq == std::string::npos) { // no equal sign
                            ++argp;
                            if(argp >= argc) {
                                _parseErrorBit = true;
                                return -1;
                            } else {
                                // take the whole argument as a parameter
                                (*it)->fillField(std::string(argv[argp]));
                            }
                        } else { // with equal sign
                            (*it)->fillField(std::string(arg, eq+1));
                        }
                    } // else good. Not taking arguments
                    // Offer it and remove from list
                    (*it)->offer();
                    op.erase(it);
                }
                ++argp;
            }
            break;
        case ArgType::Short:
            {
                // Offering is set at the end.

                while(arg.length() >= 2 && !shouldEnd) {
                    char flag = arg[1];
                    auto it = std::find_if(op.begin(), op.end(),
                        [&flag](OptionBase* ob) { return ob->findHit(flag); });
                    if(it == op.end()) { // Not found
                        shouldEnd = true;
                    } else { // Found
                        if((*it)->takesArg()) {
                            if(arg.length() == 2) { // Single short option
                                ++argp;
                                if(argp >= argc) {
                                    _parseErrorBit = true;
                                    return -1;
                                } else {
                                    // take the whole argument as a parameter
                                    (*it)->fillField(std::string(argv[argp]));
                                }
                            } else {
                                // Use the rest as parameter
                                (*it)->fillField(std::string(arg, 2));
                            }
                            // Offer it and remove from the list
                            (*it)->offer();
                            op.erase(it);
                            break;
                        } else { // Not taking arguments
                            arg.erase(1); // Remove short option
                            // Offer it and remove from the list
                            (*it)->offer();
                            op.erase(it);
                        }
                    }
                } // End of cropping short argument
                ++argp;
            }
            break;
        case ArgType::ArgOrCmd:
            {
                if(itPos == _pos.end()) {
                    shouldEnd = true;
                    break;
                }

                int newArgp = (*itPos)->parseThis(argc, argv, argp);
                if(newArgp < 0) {
                    _parseErrorBit = true;
                    return -1;
                } else {
                    argp = newArgp;
                }

                // Now update itPos
                ++itPos;
            }
            break;
        }
    }

    // Check required options
    for(auto ob: op) {
        if(ob->isRequired()) {
            _logicErrorBit = true;
            std::ostringstream oss;
            oss << ob->getMatch() << " (" << ob->getDescription() << ") is required.";
            _logicErrorInfo.emplace_back(oss.str());
            return -1;
        }
    }

    // Check required pos eles
    for(; itPos != _pos.end(); ++itPos) {
        if((*itPos)->isRequired()) {
            _logicErrorBit = true;
            _logicErrorInfo.emplace_back("Insufficient arguments.");
            return -1;
        }
    }
        
    return argp;

}

void PosHolder::printContent(std::ostream& os)const {
    for(auto opPtr: _op) {
        bool required = opPtr->isRequired();
        os << (required? " ": " [");

        if(opPtr->getFlagLong().length()) {
            os << "--" << opPtr->getFlagLong();
            if(opPtr->takesArg()) os << "=" << opPtr->getArgName();
        }
        else if(opPtr->getFlagShort()) os << '-' << opPtr->getFlagShort();

        os << (required? "": "]");
    }
    for(auto pe: _pos) {
        bool required = pe->isRequired();
        bool group = pe->isMutuallyExclusive();
        os << (required? (group? " (": " "): " [");
        pe->printContent(os);
        os << (required? (group? ")": ""): "]");
    }
}

void PosMutuallyExclusive::_preprocess() {
    if(_pos.empty()) {
        _inputFailBit = true;
        _inputFailInfo.emplace_back("Mutually exclusive group cannot be empty.");
    }

    for(auto pe: _pos) {
        if(!pe->isRequired()) {
            _inputFailBit = true;
            _inputFailInfo.emplace_back("Elements inside mutually exclusive group must be required.");
            break;
        }
    }
}

int PosMutuallyExclusive::parse(int argc, char** argv, int argp) {
    _parseErrorBit = false;
    _logicErrorBit = false;
    _logicErrorInfo.clear();

    if(!*this) return -1;

    size_t successCnt = 0;
    int lastGoodIndex = 0;
    int lastGoodArgp = -1;

    int newArgp;
    size_t n = _pos.size();
    for(size_t i = 0; i < n; ++i) {
        newArgp = _pos[i]->parseThis(argc, argv, argp);
        if(newArgp >= 0) {
            ++successCnt;
            lastGoodIndex = i;
            lastGoodArgp = newArgp;
        }
    }

    if(successCnt == 0) {
        _parseErrorBit = true;
        return -1;
    } else if(successCnt > 1) {
        _logicErrorBit = true;
        _logicErrorInfo.emplace_back("More than 1 options that are mutually exclusive are supplied.");
        return -1;
    }

    _index = lastGoodIndex;
    return lastGoodArgp;
}

void Command::_preprocess() {
    if(getArgType(_name) != ArgType::ArgOrCmd) {
        _inputFailBit = true;
        std::ostringstream oss;
        oss << "Command name " << _name << " is invalid.";
        _inputFailInfo.emplace_back(oss.str());
    }
    
    // Callback error
    if(!_activate) {
        _inputFailBit = true;
        _inputFailInfo.emplace_back("Command activate callback not set.");
    }
}

int Command::parse(int argc, char** argv, int argp) {
    if(argp >= argc || (!_main && std::strcmp(argv[argp], _name) != 0)) {
        _parseErrorBit = true;
        return -1;
    }

    // Uncheck error bit
    _parseErrorBit = false;

    if(!*this) return -1;

    ++argp; // Skip the "name" of this command

    if(!_content) return argp;

    int newArgp = _content->parseThis(argc, argv, argp);
    if(newArgp < 0) {
        _parseErrorBit = true;
        return -1;
    }

    // For main command, report error if there are unused arguments
    if(_main && newArgp < argc) {
        _unusedArgumentBit = true;
        for(int idx = newArgp; idx < argc; ++idx) _unusedArguments.emplace_back(argv[idx]);
        return -1;
    }
        
    return newArgp;
}

void Command::printUsage(std::ostream& os)const {
    os << "Usage:" << '\n';
    if(!_content) return;

    std::unordered_set<const PosElement*> cmd;
    std::unordered_set<const OptionBase*> op;
    _content->addToCmdOpSet(cmd, op);

    switch(_content->getPosType()) {
    case PosType::Command:
    case PosType::Arg:
    case PosType::Holder:
        os << "  " << _name;
        _content->printContent(os);
        os << '\n';
        break;
    case PosType::MutuallyExclusive:
        for(auto pe: static_cast<PosMutuallyExclusive*>(_content)->getPos()) {
            os << "  " << _name;
            pe->printContent(os);
            os << '\n';
        }
        break;
    }
    os << std::endl;

    if(!cmd.empty()) {
        os << "Commands:\n";
        for(const PosElement* pe: cmd) {
            usagePairFormatter(std::string(pe->getCommandName()), pe->getDescription(), os);
        }
        os << std::endl;
    }
    if(!op.empty()) {
        os << "Options:\n";
        for(const OptionBase* ob: op) {
            std::ostringstream oss;

            if(ob->getFlagShort()) oss << '-' << ob->getFlagShort();
            else oss << "    ";

            if(!ob->getFlagLong().empty()) {
                if(ob->getFlagShort()) oss << ", ";
                oss << "--" << ob->getFlagLong();
                if(ob->takesArg()) oss << '=' << ob->getArgName();
            } else if(ob->takesArg()) oss << " " << ob->getArgName();

            usagePairFormatter(oss.str(), ob->getDescription(), os);
        }
        os << std::endl;
    }

    os << std::endl;
}



} // namespace commandline
} // namespace medyan
