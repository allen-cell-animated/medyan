#include "core/command_line_parser.h"

#include <algorithm>
#include <cstring>

namespace medyan {
namespace commandline {

namespace {
    void usagePairFormatter(const std::string& key, const std::string& description, std::ostream& os) {
        static const size_t leftMargin = 2;
        static const size_t maxKeySize = 13;
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
    if(!_activate) _inputFailBit = true;
}

bool OptionBase::findHit(char argShort)              { return _flagShort == argShort; }
bool OptionBase::findHit(const std::string& argLong) { return _flagLong  == argLong;  }

void Command::_preprocess() {
    if(getArgType(_name) != ArgType::ArgOrCmd) _inputFailBit = true;
    
    // Callback error
    if(!_activate) _inputFailBit = true;
}

int Command::parse(int argc, char** argv, int argp) {
    // This command should have already passed name check before entering
    // So no self-check is performed

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
    ++argp; // Skip the "name" of this command
    while(argp < argc && !shouldEnd) {
        ArgType t = getArgType(argv[argp]);
        std::string arg {argv[argp]};

        // parse arg
        switch(t) {
        case ArgType::Long:
            {
                size_t eq = arg.find('=');
                std::string flag(arg, 2, eq);
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

                while(arg.length() >= 2) {
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

                // Un-offering current argument
                (*itPos)->offer(false);

                if((*itPos)->isRequired()) {
                    // We don't care about the offer for required positional elements.
                    if((*itPos)->isCommand()) {
                        if(arg == (*itPos)->getCommandName()) { // Name matches
                            int newArgp = (*itPos)->parse(argc, argv, argp);
                            if(newArgp < 0) {
                                _parseErrorBit = true;
                                return -1;
                            } else {
                                argp = newArgp; // Update array index
                            }
                        } else { // Name does not match
                            _parseErrorBit = true;
                            return -1;
                        }
                    } else if((*itPos)->isArgument()) {
                        (*itPos)->fillField(arg);
                        ++argp;
                    } else { // Other cases, directly parse them.
                        int newArgp = (*itPos)->parse(argc, argv, argp);
                        if(newArgp < 0) {
                            _parseErrorBit = true;
                            return -1;
                        } else {
                            argp = newArgp; // Update array index
                        }
                    }
                } else { // If it is optional
                    bool offerThisTime = true;
                    int newArgp = -1;
                    do {
                        (*itPos)->offer(offerThisTime);
                        if((*itPos)->isCommand()) {
                            if(arg == (*itPos)->getCommandName()) { // Name matches
                                newArgp = (*itPos)->parse(argc, argv, argp);
                                if(newArgp < 0) {
                                    // Things aren't correct, but we'll wait
                                } else {
                                    argp = newArgp; // Update array index
                                    break;
                                }
                            } else { // Name does not match
                                // Things aren't correct, but we'll wait
                            }
                        } else if((*itPos)->isArgument()) {
                            (*itPos)->fillField(arg); // Fill it anyway
                            ++argp;
                            newArgp = argp; // Set positive value, because we need to check it later
                            break;
                        } else { //  Other cases
                            newArgp = (*itPos)->parse(argc, argv, argp);
                            if(newArgp < 0) {
                                // Things aren't correct, but we'll wait
                            } else {
                                argp = newArgp; // Update array index
                                break;
                            }
                        }
                        offerThisTime = !offerThisTime; // false for the 2nd time
                    } while(!offerThisTime && newArgp < 0); // This will exit for at most 2 iterates.
                    if(newArgp < 0) {
                        _parseErrorBit = true;
                        return -1;
                    } // Else the argp has already been updated inside the loop, so good!
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
            std::ostringstream oss;
            oss << "Insufficient arguments for " << _name;
            _logicErrorInfo.emplace_back(oss.str());
            return -1;
        }
    }
        
    return argp;
}

void Command::printUsage(std::ostream& os)const {
    os << "Usage: " << _name;
    for(auto& opPtr: _op) {
        bool required = opPtr->isRequired();
        os << (required? " ": " [");

        if(opPtr->getFlagLong().length()) os << "--" << opPtr->getFlagLong();
        else if(opPtr->getFlagShort()) os << '-' << opPtr->getFlagShort();

        if (opPtr->takesArg()) os << ' ' << opPtr->getArgName();

        os << (required? "": "]");
    }
    for(auto pe: _pos) {
        bool required = pe->isRequired();
        os << (required? " ": " [");
        pe->printContent(os);
        if(!required) os << "]";
    }
    os << '\n' << std::endl;

    std::vector<PosElement*> commands;
    std::copy_if(_pos.begin(), _pos.end(), std::back_inserter(commands),
        [](PosElement* pe) { return pe->isCommand(); });
    if(!commands.empty()) {
        os << "Commands:\n";
        for(PosElement* pe: commands) {
            usagePairFormatter(std::string(static_cast<Command*>(pe)->_name), pe->getDescription(), os);
        }
        os << std::endl;
    }
    if(!_op.empty()) {
        os << "Options:\n";
        for(auto& opPtr: _op) {
            usagePairFormatter(std::string(opPtr->getMatch()), opPtr->getDescription(), os);
        }
        os << std::endl;
    }
}



} // namespace commandline
} // namespace medyan
