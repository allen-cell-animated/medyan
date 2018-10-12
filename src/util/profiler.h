#ifndef MEDYAN_UTIL_PROFILER_H
#define MEDYAN_UTIL_PROFILER_H

#include <chrono>
#include <ostream>
#include <string>
#include <type_traits>

#include "util/io/log.h"

namespace profiler {

#ifdef USE_PROFILER
constexpr bool useProfiler = true;
#else
constexpr bool useProfiler = false;
#endif

namespace internal {

template<
    bool enable,
    bool raii
> class SimpleTimerImpl;

// Specialization
template<bool raii> class SimpleTimerImpl< true, raii > {
private:
    using _time_point_t = std::chrono::time_point< std::chrono::steady_clock >;
    using _time_diff_t = decltype(_time_point_t() - _time_point_t());

    _time_point_t _start;
    _time_diff_t _elapsed;

    std::string _name;
public:
    // Constructor for block timer
    template< typename = typename std::enable_if< raii >::type >
    SimpleTimerImpl(const std::string& name) : _name(name) {
        start();
    }
    // Constructor for manual timer
    template< typename = typename std::enable_if< !raii >::type >
    SimpleTimerImpl(const std::string& name) : _name(name) {}
    // Destructor. Unfortunately SFINAE is not an option here.
    ~SimpleTimerImpl() {
        if /* constexpr since c++17 */ (raii) {
            elapse();
            report();
        }
    }

    // Record the current time as start time.
    void start() {
        _start = std::chrono::steady_clock::now();
    }
    // Calculate the time elapsed from start to current time.
    void elapse() {
        _elapsed = std::chrono::steady_clock::now() - _start;
    }
    // Output the timer result.
    void report() {
        LOG(INFO) << "Time elapsed for " << _name << ": "
            << std::chrono::duration_cast< std::chrono::duration< double > >(_elapsed).count()
            << "s.";
    }
};

template<bool raii> class SimpleTimerImpl< false, raii > {
public:
    // General constructor which does nothing
    template< typename... Args >
    SimpleTimerImpl(Args&&... args) {}

    void start() {}
    void elapse() {}
    void report() {}
};

} // namespace internal

using Timer = internal::SimpleTimerImpl< useProfiler, false >;
using ScopeTimer = internal::SimpleTimerImpl< useProfiler, true >;

} // namespace profiler

#endif
