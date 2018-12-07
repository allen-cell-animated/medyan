#ifndef MEDYAN_UTIL_PLATFORM_H
#define MEDYAN_UTIL_PLATFORM_H

#if defined(__unix__) || defined(__APPLE__) && defined(__MACH__)
    #define PLATFORM_UNIX_LIKE
#elif defined(_WIN32)
    #define PLATFORM_WINDOWS
#endif

#endif
