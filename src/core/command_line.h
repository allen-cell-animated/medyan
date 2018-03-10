#ifndef MEDYAN_CORE_COMMAND_LINE_H
#define MEDYAN_CORE_COMMAND_LINE_H

namespace medyan {
namespace commandline {

/**
 * MEDYAN - read from command line input.
 * 
 * Some other initializations may also happen inside this function, such as
 *   - Logger initialization
 * 
 * Note that this function may abort the program on error.
 */
void initializeFromCommandLine(int argc, char** argv);


} // namespace commandline
} // namespace medyan

#endif
