#ifndef MEDYAN_UTIL_PROFILER_H
#define MEDYAN_UTIL_PROFILER_H

#include <atomic>
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

using time_count_float_t = double;

namespace internal {

/**
 * SimpleTimerImpl is used individually to measure a certain time elapsed.
 * 
 * Template parameters
 *   - enable: If set false, no member variables should exist. All function
 *     calls should be valid but empty.
 *   - raii: If set true, timer will automatically measure time between
 *     initialization and destruction (normally at the end of the scope).
 *   - print: If set true, calling report() will generate a log message
 *     containing the name and elapsed time.
 *   - worker: If set true, the timer must specify a manager obejct. The
 *     elapsed time will be automatically submitted to the manager.
 * 
 * Note:
 *   - At least one in print and worker must be true.
 */
template< bool enable, bool raii, bool print, bool worker> class SimpleTimerImpl;

/**
 * TimerManagerImpl is used to collect timer results from worker timers.
 * 
 * Template parameters
 *   - enable: If set false, no member variables should exist. All function
 *     calls should be valid but empty.
 */
template< bool enable > class TimerManagerImpl;


// Helper classes for SimpleTimerImpl
template< bool enable > struct SimpleTimerMember {};
template<> struct SimpleTimerMember<true> {
    using time_point_t = std::chrono::time_point< std::chrono::steady_clock >;
    using time_diff_t = decltype(time_point_t() - time_point_t());
    time_point_t start;
    time_diff_t elapsed;
};

template< bool enable, bool print > struct SimpleTimerPrintMember {
    SimpleTimerPrintMember(const std::string& name) {} // Discard input
};
template<> struct SimpleTimerPrintMember< true, true > {
    std::string name;
    SimpleTimerPrintMember(const std::string& name) : name(name) {}
};

template< bool enable, bool worker > struct SimpleTimerManagerMember {
    using timer_manager_t = TimerManager< enable >;
    SimpleTimerManagerMember(timer_manager_t& manager) {} // Discard input
};
template<> struct SimpleTimerManagerMember< true, true > {
    using timer_manager_t = TimerManager< true >;
    timer_manager_t& manager;
    SimpleTimerManagerMember(timer_manager_t& manager) : manager(manager) {}
};

// Definition
template< bool enable, bool raii, bool print, bool worker > class SimpleTimerImpl
    : public SimpleTimerMember< enable >,
      public SimpleTimerPrintMember< enable, print >,
      public SimpleTimerManagerMember< enable, worker > {
public:
    // Constructors
    template< typename T = void, typename std::enable_if< print && worker, T >::type* = nullptr >
    SimpleTimerImpl(const std::string& name, typename SimpleTimerManagerMember<enable, worker>::timer_manager_t& manager)
        : SimpleTimerPrintMember<enable, print>(name), SimpleTimerManagerMember<enable, worker>(manager) {
        if /* constexpr since c++17 */ (raii) start();
    }
    template< typename T = void, typename std::enable_if< print && !worker, T >::type* = nullptr >
    SimpleTimerImpl(const std::string& name)
        : SimpleTimerPrintMember<enable, print>(name) {
        if /* constexpr since c++17 */ (raii) start();
    }
    template< typename T = void, typename std::enable_if< !print && worker, T >::type* = nullptr >
    SimpleTimerImpl(typename SimpleTimerManagerMember<enable, worker>::timer_manager_t& manager)
        : SimpleTimerManagerMember<enable, worker>(manager) {
        if /* constexpr since c++17 */ (raii) start();
    }
    // Non-printing and non-worker timer is not allowed

    // Destructor. SFINAE cannot apply here.
    ~SimpleTimerImpl() {
        if /* constexpr since c++17 */ (raii) {
            elapse();
            report();
        }
    }

    // Record the current time as start time.
    template< typename T = void, typename std::enable_if< enable, T >::type* = nullptr >
    void start() { this->start = std::chrono::steady_clock::now(); }
    template< typename T = void, typename std::enable_if< !enable, T >::type* = nullptr >
    void start() {}

    // Calculate the time elapsed from start to current time.
    template< typename T = void, typename std::enable_if< enable, T >::type* = nullptr >
    void elapse() {
        this->elapsed = std::chrono::steady_clock::now() - this->start;
        _submit();
    }
    template< typename T = void, typename std::enable_if< !enable, T >::type* = nullptr >
    void elapse() {}

private:
    // _submit will submit the current result to the manager.
    // Will be called by elapse()
    template< typename T = void, typename std::enable_if< enable && worker, T >::type* = nullptr >
    void _submit()const {
        this->manager.submit(std::chrono::duration_cast< std::chrono::duration< time_count_float_t > >(this->elapsed).count());
    }
    template< typename T = void, typename std::enable_if< !(enable && worker), T >::type* = nullptr >
    void _submit()const {}

public:
    // Output the timer result.
    template< typename T = void, typename std::enable_if< enable && print, T >::type* = nullptr >
    void report()const {
        LOG(INFO) << "Time elapsed for " << this->name << ": "
            << std::chrono::duration_cast< std::chrono::duration< time_count_float_t > >(this->elapsed).count()
            << "s.";
    }
    template< typename T = void, typename std::enable_if< !(enable && print), T >::type* = nullptr >
    void report()const {}
};


// Helper classes for TimerManager
template< bool enable > struct TimerManagerMember {
    TimerManagerMember() = default;
};
template<> struct TimerManagerMember<true> {
    std::atomic<unsigned long long> count;
    std::atomic<time_count_float_t> elapsed;
    TimerManagerMember() : count(0ull), elapsed(0.0) {}
};

template< bool enable > struct TimerManagerPrintMember {
    TimerManagerPrintMember(const std::string& name) {} // Discard input
};
template<> struct TimerManagerPrintMember<true> {
    std::string name;
    TimerManagerPrintMember(const std::string& name) : name(name) {}
};

// Definition
template< bool enable > class TimerManagerImpl
    : public TimerManagerMember<enable> {
public:
    // Constructors
    TimerManagerImpl(const std::string& name)
        : TimerManagerPrintMember<enable>(name) {}

    // Reset the counter and timer.
    template< typename T = void, typename std::enable_if< enable, T >::type* = nullptr >
    void reset() { this->elapsed = 0.0; this->count = 0ull; }
    template< typename T = void, typename std::enable_if< !enable, T >::type* = nullptr >
    void reset() {}

    // submit(elapsed) will add the time to the total elapsed time, and will
    // increase the timer fire count by 1.
    // Normally it should only be called by worker timer.
    template< typename T = void, typename std::enable_if< enable, T >::type* = nullptr >
    void submit(time_count_float_t elapsedSingle) {
        ++(this->count);
        this->elapsed += elapsedSingle;
    }
    template< typename T = void, typename std::enable_if< !enable, T >::type* = nullptr >
    void submit(time_count_float_t elapsedSingle) {} // Do nothing

    // print the total time information
    template< typename T = void, typename std::enable_if< enable, T >::type* = nullptr >
    void report()const {
        LOG(INFO)
            << "Time elapsed for " << this->count.load() << " occurrences for " << this->name
            << this->elapsed.load() << "s.";
    }
    template< typename T = void, typename std::enable_if< enable, T >::type* = nullptr >
    void report()const {}
};

} // namespace internal

using Timer            = internal::SimpleTimerImpl< useProfiler, false, true,  false >;
using ScopeTimer       = internal::SimpleTimerImpl< useProfiler, true,  true,  false >;
using TimerWorker      = internal::SimpleTimerImpl< useProfiler, false, false, true  >;
using ScopeTimerWorker = internal::SimpleTimerImpl< useProfiler, true,  false, true  >;

using TimerManager = internal::TimerManagerImpl< useProfiler >;

} // namespace profiler

#endif
