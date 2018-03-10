#include "core/command_line_parser.h"

#include <algorithm>
#include <cstring>

namespace medyan {
namespace commandline {

int Command::parse2(int argc, char** argv, int argp) {
    _unusedArgs.clear();
    _unusedArgBit = false;

    // Evaluate itself first
    if(!_evaluated) {
        evaluate();
        if(!*this) {
            _parseErrorBit = true;
            return -1;
        }
        argp += 1;
    }
    _evaluated = true;

    // Try to read to the end of the arguments
    int idx = argp;
    for(; idx < argc; ++idx) {
        ArgType t = getArgType(argv[idx]);
        std::string arg {argv[idx]};

        int maxMove = 0; // How many arguments should idx move forward

        // Evaluate and parse arg
        switch(t) {
        case ArgType::Long:
        case ArgType::Short:
            for(auto& opPtr: _op) {
                if(opPtr->findHit2(arg, t)) {
                    int iMove = opPtr->evaluate2(argc, argv, idx);
                    if(!*opPtr) {
                        _parseErrorBit = true;
                        return -1;
                    }
                    if(iMove > maxMove) maxMove = iMove;

                    // Short args may be evaluated by other options, e.g. "-af" would be evaluated by "-a" and "-f"
                    if(t != ArgType::Short || arg.length() == 2) {
                        arg.clear();
                        break;
                    } else {
                        char f = opPtr->getFlagShort();
                        arg.erase(std::remove(arg.begin(), arg.end(), f), arg.end());
                    }
                }

            }

            if(!arg.empty()) _unusedArgs.push_back(arg);

            // Increase idx by required
            idx += maxMove;

            break;
        case ArgType::ArgOrCmd:
            // Currently only one subcommand in each command is allowed,
            // because any command will try to read to the end of the arg list
            for(auto& cmdPtr: _subcmd) {
                if(cmdPtr->_name == arg) {
                    if(cmdPtr->parse(argc, argv, idx) == -1) {
                        _subcmdErrorBit = true;
                        return -1;
                    }
                    break;
                }
            }
            break;
        }
    }

    if(_unusedArgs.size()) {
        _unusedArgBit = true;
        return -1;
    }

    // Check for option logic
    for(auto& opPtr: _op) {
        // Required option
        if(opPtr->isRequired() && !opPtr->isEvaluated()) {
            _logicErrorBit = true;
            std::ostringstream oss;
            oss << opPtr->getMatch() << " (" << opPtr->getDescription() << ") is required.";
            _logicErrorContent.emplace_back(oss.str());
            return -1;
        }
        // Excluding
        if(opPtr->isEvaluated())
            for(auto ex: opPtr->getExcluding()) {
                if(ex->isEvaluated()) {
                    _logicErrorBit = true;
                    std::ostringstream oss;
                    oss << ex->getMatch() << " cannot be specified when " << opPtr->getMatch() << " is present.";
                    _logicErrorContent.emplace_back(oss.str());
                }
            }
    }
        
    return idx;
}

void Command::printUsage(std::ostream& out)const {
    out << "Usage: " << _name;
    for(auto& opPtr: _op) {
        bool required = opPtr->isRequired();
        out << (required? " ": " [");

        if(opPtr->getFlagLong().length()) out << "--" << opPtr->getFlagLong();
        else if(opPtr->getFlagShort()) out << '-' << opPtr->getFlagShort();

        if (opPtr->takesArg()) out << '=' << opPtr->getArgName();

        out << (required? "": "]");
    }
    if(!_subcmd.empty()) {
        //bool required = defaultOp && defaultOp->isRequired();
        //if(!required) out << "[";
        out << " command";
        //if(!required) out << "]";
    }
    out << '\n' << std::endl;

    if(!_subcmd.empty()) {
        out << "Commands:\n";
        for(auto& cmdPtr: _subcmd) {
            out << "  ";
            if(std::strlen(cmdPtr->_name) > 13) {
                out << cmdPtr->_name << '\n' << std::string(' ', 14);
            } else {
                std::string tmp = cmdPtr->_name;
                tmp.resize(16, ' ');
                out << tmp;
            }
            out << cmdPtr->getDescription() << '\n';
        }
        out << std::endl;
    }
    if(!_op.empty()) {
        out << "Options:\n";
        for(auto& opPtr: _op) {
            out << "  ";
            if(std::strlen(opPtr->getMatch()) > 13) {
                out << opPtr->getMatch() << '\n' << std::string(' ', 14);
            } else {
                std::string tmp {opPtr->getMatch()};
                tmp.resize(16, ' ');
                out << tmp;
            }
            out << opPtr->getDescription() << '\n';
        }
        out << std::endl;
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
        } else if(len && eachStr[0] != '-') {
            // Command flag
            _flagCommand = eachStr;
        } else {
            // Failed
            _inputFailBit = true;
        }
    }
    
    // Callback error
    if(!_activate) _inputFailBit = true;
}

bool OptionBase::findHit2(const std::string& arg, ArgType argType) {
    switch(argType) {
    case ArgType::Short:
        if(_flagShort) {
            size_t pos = arg.find(_flagShort);
            if(pos != std::string::npos) {
                if(_takesArg) {
                    if(pos == arg.length() - 1) return true;
                }
                else return true;
            }
        }
        break;
    case ArgType::Long:
        if(_flagLong == std::string(arg, 2)) return true;
        break;
    case ArgType::ArgOrCmd:
        if(_flagCommand == arg) return true;
        break;
    default:
        break;
    }

    return false;
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
                    if(it->takesArg()) {
                        if(eq == std::string::npos) { // no equal sign
                            ++argp;
                            if(argp >= argc) {
                                _parseErrorBit = true;
                                return -1;
                            } else {
                                // take the whole argument as a parameter
                                it->fillField(std::string(argv[argp]));
                            }
                        } else { // with equal sign
                            it->fillField(std::string(arg, eq+1));
                        }
                    } // else good. Not taking arguments
                    // Offer it and remove from list
                    it->offer();
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
                        if(it->takesArg()) {
                            if(arg.length() == 2) { // Single short option
                                ++argp;
                                if(argp >= argc) {
                                    _parseErrorBit = true;
                                    return -1;
                                } else {
                                    // take the whole argument as a parameter
                                    it->fillField(std::string(argv[argp]));
                                }
                            } else {
                                // Use the rest as parameter
                                it->fillField(std::string(arg, 2));
                            }
                        } else { // Not taking arguments
                            arg.erase(1); // Remove short option
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
                itPos->offer(false);

                if(itPos->isRequired()) {
                    // We don't care about the offer for required positional elements.
                    if(itPos->isCommand()) {
                        if(arg == itPos->getCommandName()) { // Name matches
                            int newArgp = itPos->parse(argc, argv, argp);
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
                    } else if(itPos->isArgument()) {
                        it->fillField(arg);
                        ++argp;
                    } else { // Other cases, directly parse them.
                        int newArgp = itPos->parse(argc, argv, argp);
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
                        itPos->offer(offerThisTime);
                        if(itPos->isCommand()) {
                            if(arg == itPos->getCommandName()) { // Name matches
                                newArgp = itPos->parse(argc, argv, argp);
                                if(newArgp < 0) {
                                    // Things aren't correct, but we'll wait
                                } else {
                                    argp = newArgp; // Update array index
                                    break;
                                }
                            } else { // Name does not match
                                // Things aren't correct, but we'll wait
                            }
                        } else if(itPos->isArgument()) {
                            it->fillField(arg); // Fill it anyway
                            ++argp;
                            newArgp = argp; // Set positive value, because we need to check it later
                            break;
                        } else { //  Other cases
                            newArgp = itPos->parse(argc, argv, argp);
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
            oss << opPtr->getMatch() << " (" << opPtr->getDescription() << ") is required.";
            _logicErrorInfo.emplace_back(oss.str());
            return -1;
        }
    }

    // Check required pos eles
    for(; itPos != _pos.end(); ++itPos) {
        if(itPos->isRequired()) {
            _logicErrorBit = true;
            std::ostringstream oss;
            oss << "Insufficient arguments for " << _name;
            _logicErrorInfo.emplace_back(oss.str());
            return -1;
        }
    }
        
    return argp;
}



} // namespace commandline
} // namespace medyan
