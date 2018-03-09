#include "core/command_line_parser.h"

#include <algorithm>
#include <cstring>

namespace medyan {
namespace commandline {

bool Command::parse(int argc, char** argv, int argp) {
    _unusedArgs.clear();
    _unusedArgBit = false;

    _evaluated = true;

    // Evaluate default option first
    if(_defaultOp && !_defaultOp->isEvaluated()) {
        ++argp;
        int iMove = _defaultOp->evaluate(argc, argv, argp);
        if(!*_defaultOp) {
            _parseErrorBit = true;
            return false;
        }
        argp += iMove;
    }
    // Try to read to the end of the arguments
    for(int idx = argp + 1; idx < argc; ++idx) {
        ArgType t = getArgType(argv[idx]);
        std::string arg {argv[idx]};
        if(t == ArgType::Fail) {
            _invalidArgContent = arg;
            _invalidArgBit = true;
            return false;
        }

        int maxMove = 0; // How many arguments should idx move forward

        // Evaluate and parse arg
        switch(t) {
        case ArgType::Long:
        case ArgType::Short:
            for(auto& opPtr: _op) {
                if(opPtr->findHit(arg, t)) {
                    int iMove = opPtr->evaluate(argc, argv, idx);
                    if(!*opPtr) {
                        _parseErrorBit = true;
                        return false;
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
        case ArgType::Argument:
            // Currently only one subcommand in each command is allowed,
            // because any command will try to read to the end of the arg list
            for(auto& cmdPtr: _subcmd) {
                if(cmdPtr->_name == arg) {
                    if(!cmdPtr->parse(argc, argv, idx)) {
                        _subcmdErrorBit = true;
                        return false;
                    }
                    break;
                }
            }
            break;
        }
    }

    if(_unusedArgs.size()) {
        _unusedArgBit = true;
        return false;
    }

    // Check for option logic
    for(auto& opPtr: _op) {
        // Required option
        if(opPtr->isRequired() && !opPtr->isEvaluated()) {
            _logicErrorBit = true;
            std::ostringstream oss;
            oss << opPtr->getMatch() << " (" << opPtr->getDescription() << ") is required.";
            _logicErrorContent.emplace_back(oss.str());
            return false;
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
        
    return true;
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
        OptionBase* defaultOp = _subcmd[0]->getDefaultOp();
        bool required = defaultOp && defaultOp->isRequired();
        if(!required) out << "[";
        out << " command";
        if(!required) out << "]";
    }
    out << '\n' << std::endl;

    if(!_subcmd.empty()) {
        out << "Commands:\n";
        for(auto& cmdPtr: _subcmd) {
            out << "  ";
            if(cmdPtr->_name.length() > 13) {
                out << cmdPtr->_name << '\n' << std::string(' ', 14);
            } else {
                std::string tmp = cmdPtr->_name;
                tmp.resize(16, ' ');
                out << tmp;
            }
            out << cmdPtr->_defaultOp->getDescription() << '\n';
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
}

bool OptionBase::findHit(const std::string& arg, ArgType argType) {
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
    case ArgType::Argument:
        if(_flagCommand == arg) return true;
        break;
    default:
        break;
    }

    return false;
}


} // namespace commandline
} // namespace medyan
