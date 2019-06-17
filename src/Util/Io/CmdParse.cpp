#include "Util/Io/CmdParse.hpp"

#include <algorithm>

namespace cmdparse {

namespace {

void usagePairFormatter(const std::string& key, const std::string& description, std::ostream& os) {
    static const size_t leftMargin = 4;
    static const size_t maxKeySize = 21;
    static const size_t midMargin = 1;

    os << std::string(leftMargin, ' ');

    size_t len = key.length();
    os << key;
    if(len > maxKeySize) os << '\n' << std::string(leftMargin + maxKeySize + midMargin, ' ');
    else os << std::string(maxKeySize + midMargin - len, ' ');

    os << description << '\n';
}

} // namespace


std::string Option::getReadableName()const {
    std::ostringstream oss;
    if(_short) {
        oss << '-' << _short;
        if(_long.length())
            oss << " (--" << _long << ')';
    } else if(_long.length())
        oss << "--" << _long;
    return oss.str();
}

std::string Option::getUsageName()const {
    std::ostringstream oss;
    if(_short) {
        oss << '-' << _short;
        if(_long.length()) {
            oss << ", --" << _long;
            if(_hasVariable) {
                oss << "=<" << _variableName << '>';
            }
        } else {
            if(_hasVariable) {
                oss << " <" << _variableName << '>';
            }
        }
    } else if(_long.length()) {
        oss << "--" << _long;
        if(_hasVariable) {
            oss << "=<" << _variableName << '>';
        }
    }
    return oss.str();
}


void Command::_ruleCheck()const {

    // Check positional argument positions
    {
        size_t num = _posArgs.size();
        size_t numOpt = 0;
        size_t numList = 0;
        for(size_t i = 0; i < num; ++i) {
            if(_posArgs[i]->isRequired()) {
                if(numOpt)
                    throw CommandLogicError("Required positional argument should not appear after optional positional argument.");
            } else {
                ++numOpt;
            }
            if(_posArgs[i]->isList()) {
                ++numList;
                if(numList > 1)
                    throw CommandLogicError("Command should not contain more than one positional argument list.");
            }
        }
    }

    // Check option names
    for(auto& op : _options) {
        if(op->getShortName() == 0 && op->getLongName() == "")
            throw CommandLogicError("Command should not contain options without name.");
    }

    // Recursively check all subcommands
    for(auto& sc : _subcommands) {
        sc->_ruleCheck();
    }
}

void Command::_parsePosArg(const std::string& arg) {
    if(_state._posArgIndex >= _posArgs.size())
        throw ParsingError("More positional argument specified than required: " + arg);

    ++_state._posArgCount;

    _posArgs[_state._posArgIndex]->occur();
    _posArgs[_state._posArgIndex]->activate(arg);
    if(!_posArgs[_state._posArgIndex]->isList())
        ++_state._posArgIndex;

}

void Command::_parse(std::vector<std::string>& feed, size_t argp) {
    for(; argp < feed.size(); ++argp) {
        std::string& thisArg = feed[argp];

        // Deduce the argument type
        if(_state._parseAsPosArg) {
            _parsePosArg(thisArg);

            continue;
        }
        
        if (thisArg == "--") {
            _state._parseAsPosArg = true;
            continue;
        }
        
        {
            // if no positional argument shows up, check if the argument is a subcommand
            auto cmdMatch = std::find_if(
                _subcommands.begin(), _subcommands.end(),
                [&thisArg](const std::unique_ptr<Command>& sc) { return sc->getName() == thisArg; }
            );
            if(_state._posArgCount == 0 && cmdMatch != _subcommands.end()) {
                // hand over the parsing to the subcommand
                (*cmdMatch)->occur();
                if((*cmdMatch)->activate) (*cmdMatch)->activate();
                return (*cmdMatch)->_parse(feed, argp + 1);
            }
        }

        // Check if the argument is a short option.
        if(thisArg.length() >= 2 && thisArg[0] == '-' && thisArg[1] != '-') {
            char shortName = thisArg[1];
            Command *p = this;
            auto shortNameMatch = [shortName](const std::unique_ptr<Option>& op) { return op->getShortName() == shortName; };
            auto shortOpMatch = std::find_if(
                p->_options.begin(), p->_options.end(),
                shortNameMatch
            );
            while(shortOpMatch == p->_options.end() && p->_parent) {
                p = p->_parent;
                shortOpMatch = std::find_if(
                    p->_options.begin(), p->_options.end(),
                    shortNameMatch
                );
            }

            if(shortOpMatch == p->_options.end())
                throw ParsingError(std::string("Unrecognized short option -") + shortName);

            // Short Option found
            (*shortOpMatch)->occur();

            if((*shortOpMatch)->hasVariable()) {
                if(thisArg.length() > 2) {
                    // Use the rest as argument
                    std::string content = thisArg.substr(2);
                    (*shortOpMatch)->activate(content);
                } else {
                    // Use the next as argument
                    if(argp + 1 == feed.size()) {
                        // No next argument
                        throw ParsingError(std::string("Argument required for option -") + shortName);
                    } else {
                        ++argp;
                        (*shortOpMatch)->activate(feed[argp]);
                    }
                }
            } else {
                (*shortOpMatch)->activate("");
                if(thisArg.length() > 2) {
                    // Prepare the rest for next round
                    thisArg = '-' + thisArg.substr(2);
                    --argp; // To keep argp unchanged in the next iteration
                } // else do nothing
            }

            continue;
        }
        
        // Check if the argument is a long option.
        if(thisArg.length() > 2 && thisArg[0] == '-' && thisArg[1] == '-') {
            size_t eqIdx = thisArg.find('=');
            std::string longName = (
                eqIdx == std::string::npos ?
                thisArg.substr(2) :
                thisArg.substr(2, eqIdx - 2)
            );

            if(longName.length() == 0)
                throw ParsingError("Invalid option " + thisArg);
            
            Command *p = this;
            auto longNameMatch = [&longName](const std::unique_ptr<Option>& op) { return op->getLongName() == longName; };
            auto longOpMatch = std::find_if(
                p->_options.begin(), p->_options.end(),
                longNameMatch
            );
            while(longOpMatch == p->_options.end() && p->_parent) {
                p = p->_parent;
                longOpMatch = std::find_if(
                    p->_options.begin(), p->_options.end(),
                    longNameMatch
                );
            }

            if(longOpMatch == p->_options.end())
                throw ParsingError(std::string("Unrecognized long option --") + longName);

            // Long Option found
            (*longOpMatch)->occur();

            if((*longOpMatch)->hasVariable()) {
                if(eqIdx != std::string::npos) {
                    // Use everything after '=' as the argument
                    std::string content = thisArg.substr(eqIdx + 1);
                    (*longOpMatch)->activate(content);
                } else {
                    // Use the next as argument
                    if(argp + 1 == feed.size()) {
                        // No next argument
                        throw ParsingError(std::string("Argument required for option --") + longName);
                    } else {
                        ++argp;
                        (*longOpMatch)->activate(feed[argp]);
                    }
                }
            } else {
                (*longOpMatch)->activate("");
                if(eqIdx != std::string::npos)
                    throw ParsingError(std::string("Option --") + longName + std::string(" should not take variable"));
            }

            continue;
        }

        // Default as positional argument
        _parsePosArg(thisArg);
    }
}


void Command::_validate()const {
    // Check for required positional argument
    for(auto& pa : _posArgs) {
        if(pa->isRequired() && pa->getOccurenceCount() == 0)
            throw ValidationError(std::string("Must specify positional argument ") + pa->getName());
    }

    // Check for required option
    for(auto& op : _options) {
        if(op->isRequired() && op->getOccurenceCount() == 0)
            throw ValidationError(std::string("Must specify option ") + op->getReadableName());
    }

    // Run user validation
    if(_userValidation) _userValidation();

    // Recursively validate subcommands
    for(auto& sc : _subcommands) {
        if(sc->_state._occurenceCount)
            sc->_validate();
    }
}

void Command::printUsage(std::ostream& os)const {
    os << "Usage:\n";

    std::ostringstream ossFullName;
    if(_inheritedName.length()) ossFullName << _inheritedName << ' ';
    ossFullName << _name;

    size_t numReqOp = std::count_if(
        _options.begin(), _options.end(),
        [](const std::unique_ptr<Option>& op) { return op->isRequired(); }
    );
    size_t numOptOp = _options.size() - numReqOp;

    // Usage of command
    os << "    " << ossFullName.str();
    if(numReqOp)
        os << " <options>";
    if(numOptOp)
        os << " [<options>]";
    if(_posArgs.size()) {
        os << " [--]";
        for(auto& pa : _posArgs) {
            os << ' ';
            if(!pa->isRequired()) os << '[';
            os << '<' << pa->getName() << '>';
            if(pa->isList()) os << "...";
            if(!pa->isRequired()) os << ']';
        }
    }
    os << '\n';

    // Usage with subcommand
    if(_subcommands.size())
        os << "    " << ossFullName.str() << " <command>\n";

    // Section of subcommand
    if(_subcommands.size()) {
        os << "\nCommands:\n";
        for(auto& sc : _subcommands)
            usagePairFormatter(sc->getName(), sc->getDescription(), os);
    }

    // Section of option
    if(_options.size()) {
        os << "\nOptions:\n";
        for(auto& op : _options)
            if(op->isRequired())
                usagePairFormatter(op->getUsageName(), "[Required] " + op->getDescription(), os);
        for(auto& op : _options)
            if(!op->isRequired())
                usagePairFormatter(op->getUsageName(), op->getDescription(), os);
    }

    os << std::endl;
}

} // namespace cmdparse
