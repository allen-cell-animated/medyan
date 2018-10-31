#include "util/io/cmdparse.h"

namespace cmdparse {

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

/*
void Command::_parse() {
}
*/

void Command::_validation()const {

}

} // namespace cmdparse
