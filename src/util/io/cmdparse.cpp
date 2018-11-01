#include "util/io/cmdparse.h"

#include <algorithm>

namespace cmdparse {

namespace {

void commandParsePosArg(const std::string& arg, PosArg& pa) {
}
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

    // Recursively check all subcommands
    for(auto& sc : _subcommands) {
        sc->_ruleCheck();
    }
}


void Command::_parse(std::vector<std::string>& feed, size_t argp) {
    for(; argp < feed.size(); ++argp) {
        std::string& thisArg = feed[argp];

        // Deduce the argument type
        if(_state._parseAsPosArg) {
            if(_state._posArgIndex >= _posArgs.size())
                throw ParsingError("More positional argument specified than required.");

            ++_state._posArgCount;

            _posArgs[_state._posArgIndex]->occur();
            _posArgs[_state._posArgIndex]->activate(thisArg);
            if(!_posArgs[_state._posArgIndex]->isList())
                ++_state._posArgIndex;

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
                return (*cmdMatch)->_parse(feed, argp + 1);
            }
        }

        {
            // Check if the argument is a short option.
            if(thisArg.length() >= 2 && thisArg[0] == '-' && thisArg[1] != '-') {
                char shortName = thisArg[1];
                bool foundShortName = false;
                Command *p = this;
                auto shortOpMatch = std::find_if(
                    p->_options.begin(), p->_options.end(),
                    [shortName](const std::unique_ptr<Option>& op) { return op->getShortName() == shortName; }
                );
                while(shortOpMatch == p->_options.end() && p->_parent) {
                    p = p->_parent;
                    shortOpMatch = std::find_if(
                        p->_options.begin(), p->_options.end(),
                        [shortName](const std::unique_ptr<Option>& op) { return op->getShortName() == shortName; }
                    );
                }

                if(shortOpMatch == p->_options.end())
                    throw ParsingError(std::string("Unrecognized short option -") + shortName);

                // Short Option found
                (*shortOpMatch)->occur();

                if((*shortOpMatch)->hasVariable()) {
                    // TODO having variable
                } else {
                    // TODO without variable
                }
            } else {
                // TODO long option
            }
        }

        // TODO default pos args
    }
}


void Command::_validate()const {

}

} // namespace cmdparse
