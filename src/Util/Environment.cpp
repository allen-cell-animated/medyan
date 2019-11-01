#include "Environment.hpp"

#ifdef PLATFORM_UNIX_LIKE
    #include <unistd.h>
#elif defined(PLATFORM_WINDOWS)
    #include <Windows.h>
#endif

IoEnv::IoEnv() {

#ifdef PLATFORM_UNIX_LIKE
    // Detect redirection
    stdoutRedirected = !isatty(STDOUT_FILENO);
    stderrRedirected = !isatty(STDERR_FILENO);

#elif defined(PLATFORM_WINDOWS)
    // Detect redirection and set VT mode
    HANDLE hOut = GetStdHandle(STD_OUTPUT_HANDLE);
    if(hOut == INVALID_HANDLE_VALUE) {
        stdoutRedirected = true;
    } else {
        DWORD mode = 0;
        stdoutRedirected = !GetConsoleMode(hOut, &mode);
        if(!stdoutRedirected) {
            mode |= ENABLE_VIRTUAL_TERMINAL_PROCESSING;
            SetConsoleMode(hOut, mode);
        }
    }
    HANDLE hErr = GetStdHandle(STD_ERROR_HANDLE);
    if(hErr == INVALID_HANDLE_VALUE) {
        stderrRedirected = true;
    } else {
        DWORD mode = 0;
        stderrRedirected = !GetConsoleMode(hErr, &mode);
        if(!stderrRedirected) {
            mode |= ENABLE_VIRTUAL_TERMINAL_PROCESSING;
            SetConsoleMode(hErr, mode);
        }
    }

#endif

} // IoEnv::IoEnv()
