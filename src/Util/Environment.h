#ifndef MEDYAN_UTIL_ENVIRONMENT_H
#define MEDYAN_UTIL_ENVIRONMENT_H

#if defined(__unix__) || defined(__APPLE__) && defined(__MACH__)
    #define PLATFORM_UNIX_LIKE
#elif defined(_WIN32)
    #define PLATFORM_WINDOWS
#endif

#if defined(__clang__)
    #define COMPILER_CLANG
#elif defined(__GNUC__)
    #define COMPILER_GCC
#elif defined(_MSC_VER)
    #define COMPILER_MSVC
#endif

#endif
