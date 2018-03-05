#include "core/command_line_parser.h"

namespace medyan {

namespace
bool CommandLineOptionBase::findHit(const std::string& arg, CommandLineParser::ArgType argType) {

}
void CommandLineParser::parse(int argc, char** argv) {
    std::vector<char*> unusedArgs;

    for(int idx = 1; idx < argc; ++idx) {
        std::string arg {argv[idx]};
        
        for(auto& opPtr: _op) {
            int requireArgs = 0;
            opPtr->findHit(arg, requireArgs);
        }
    }
};

} // namespace medyan
