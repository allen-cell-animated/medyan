#include "core/command_line_parser.h"

#include <algorithm>

#include "core/globals.h"

namespace medyan {

bool CommandLineParser::parse(int argc, char** argv) {
    _unusedArgs.clear();
    _unusedArgBit = false;

    for(int idx = 1; idx < argc; ++idx) {
        ArgType t = getArgType(argv[idx]);
        std::string arg {argv[idx]};
        if(t == ArgType::Fail) {
            _invalidOptionContent = arg;
            _invalidOptionBit = true;
            return false;
        }

        int maxMove = 0; // How many arguments should idx move forward
        
        for(auto& opPtr: _op) {
            if(opPtr->findHit(arg, t)) {
                if(!*opPtr) {
                    opPtr->printError();
                    _parseErrorBit = true;
                    return false;
                }

                int iMove = opPtr->evaluate(argc, argv, idx); // idx may be changed in this function
                if(!*opPtr) {
                    opPtr->printError();
                    _parseErrorBit = true;
                    return false;
                }
                if(iMove > maxMove) maxMove = iMove;

                opPtr->validate();
                if(!*opPtr) {
                    opPtr->printError();
                    _parseErrorBit = true;
                    return false;
                }

                // Short args may be evaluated by other options, e.g. "-af" would be evaluated by "-a" and "-f"
                if(t != ArgType::Short || arg.length() == 2) {
                    arg = "";
                }
                else {
                    char f = opPtr->getFlagShort();
                    arg.erase(std::remove(arg.begin(), arg.end(), f), arg.end());
                }
            }

        }

        if(!arg.empty()) _unusedArgs.push_back(arg);

        // Increase idx by required
        idx += maxMove;
    }

    if(_unusedArgs.size()) {
        _unusedArgBit = true;
        return false;
    }
        
    return true;
}

void CommandLineParser::printUsage()const {
    // TODO:
}


void CommandLineOptionBase::_preprocess() {
    int prevLoc = 0;
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

bool CommandLineOptionBase::findHit(const std::string& arg, CommandLineParser::ArgType argType) {
    switch(argType) {
    case CommandLineParser::ArgType::Short:
        if(_flagShort) {
            size_t pos = arg.find(_flagShort);
            if(pos != std::string::npos) {
                if(_numArgs) {
                    if(pos == arg.length() - 1) return true;
                }
                else return true;
            }
        }
        break;
    case CommandLineParser::ArgType::Long:
        if(_flagLong == std::string(arg, 2)) return true;
        break;
    case CommandLineParser::ArgType::Argument:
        if(_flagCommand == arg) return true;
        break;
    default:
        break;
    }

    return false;
}


} // namespace medyan
