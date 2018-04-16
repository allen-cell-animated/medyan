#ifndef MEDYAN_CORE_IO_COMMAND_LINE_H
#define MEDYAN_CORE_IO_COMMAND_LINE_H

namespace medyan {
namespace commandline {

/**
 * MEDYAN - read from command line input.
 * 
 * Note that this function may abort the program on error.
 */
void readFromCommandLine(int argc, char** argv);


} // namespace commandline
} // namespace medyan

#endif
