#ifndef MEDYAN_Util_Environment_Hpp
#define MEDYAN_Util_Environment_Hpp

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
