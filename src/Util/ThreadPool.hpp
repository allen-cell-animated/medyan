#ifndef MEDYAN_Util_ThreadPool_hpp
#define MEDYAN_Util_ThreadPool_hpp

#include <atomic>
#include <condition_variable>
#include <cstddef> // size_t
#include <functional>
#include <future>
#include <memory> // unique_ptr
#include <mutex>
#include <queue>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

// The implementation for thread pooling in MEDYAN.
//
// Purpose: Increase performance by running computational-heavy work on
//   multiple threads (possibly multiple processor cores). Thread pooling helps
//   balance work loads according to the avaiable processors, and reduces
//   the overhead of creating/destroying threads.
//
// Interface:
//   - The thread pool should be globally accessible in MEDYAN, so all the
//     interface functions must be thread-safe.
//   - Can push tasks to the queue.
//   - Can obtain the return value via a future.
//   - Can explicitly wait for specific tasks to complete.
//   - Tasks waiting for other tasks can work on pending tasks.
//   - Can terminate running tasks (via setting and polling shared states).
//
// Possible Additional Features
//   [ ] priority queues
//   [ ] thread idle/contention stats and recommendations
//   [ ] automatic load balancing

class ThreadPool {
private:

    // An implementation of function wrapper to store movable objects
    // It can be also used to store additional task information
    class FuncWrapper_ {
    private:
        struct Base_ {
            virtual ~Base_() = default;
            virtual void operator()() = 0;
        };

        template< typename F >
        struct Concrete_ : Base_ {
            F f;
            Concrete_(F&& f) : f(std::move(f)) {}
            void operator()() { f(); }
        };

    public:
        // Default and move constructor
        FuncWrapper_() = default;
        FuncWrapper_(FuncWrapper_&&) = default;

        // Constructor from callable
        template< typename F >
        FuncWrapper_(F&& f) : f_(new Concrete_<F>(std::forward<F>(f))) {}

        // Move assignment operator
        FuncWrapper_& operator=(FuncWrapper_&&) = default;

        void operator()() const { (*f_)(); }

    private:
        std::unique_ptr< Base_ > f_;
    };

    template< typename IntegralType >
    struct ScopeCounter_ {
        std::atomic< IntegralType >& c;
        ScopeCounter_(std::atomic< IntegralType >& c) : c(c) { ++c; }
        ~ScopeCounter_() { --c; }
    };

public:

    // Constructor
    ThreadPool(std::size_t numThreads) {
        // Create working threads
        threads_.reserve(numThreads);
        for(int i = 0; i < numThreads; ++i) {
            threads_.emplace_back(&ThreadPool::work_, this);
        }
    }

    // Destructor
    ~ThreadPool() {
        done_ = true;
        cvWork_.notify_all();
        for(auto& t : threads_) t.join();
    }

    // Submit a new task
    //
    // Note:
    //   - When arguments include references or pointers, it is the caller's
    //     job to ensure that the ref or ptr is valid when the job is running.
    template< typename F, typename... Args >
    auto submit(F&& f, Args&&... args) {
        using ReturnType = std::result_of_t< F(Args...) >;

        // Create a task containing the callable and the arguments
        std::packaged_task< ReturnType() > task(
            std::bind(std::forward<F>(f), std::forward<Args>(args)...) // C++20: change to perfect forwarding lambda
        );
        auto res = task.get_future();

        {
            std::lock_guard< std::mutex > guard(meQueue_);

            queue_.push(
                [task{std::move(task)}]() mutable { task(); }
            );
        }
        cvWork_.notify_one();

        return res;
    }

    // Accessors
    auto numThreads()        const noexcept { return threads_.size(); }
    auto numWorkingThreads() const noexcept { return numWorkingThreads_.load(); }
    auto numIdleThreads()    const noexcept { return numThreads() - numWorkingThreads(); }

    auto numTasks() {
        std::lock_guard< std::mutex > guard(meQueue_);
        return queue_.size();
    }

private:

    // Working thread
    void work_() {
        while(true) {
            FuncWrapper_ f;
            bool queuePopped = false;

            {
                std::unique_lock< std::mutex > lk(meQueue_);
                cvWork_.wait(
                    lk,
                    [&, this] {
                        queuePopped = tryPop_(f);
                        return queuePopped || done_;
                    }
                );
            }

            if(done_) return;
            if(queuePopped) {
                ScopeCounter_<int> c(numWorkingThreads_);
                f();
            }

        }
    }

    // Utility function for popping the queue:
    // If the queue is empty, then return false
    // Else return true and the popped value is written in the parameter x
    //
    // Note:
    //   - This function is NOT thread-safe
    bool tryPop_(FuncWrapper_& x) {
        if(queue_.empty()) return false;
        x = std::move(queue_.front());
        queue_.pop();
        return true;
    }

    std::atomic_bool done_ {false};

    std::condition_variable cvWork_;
    std::mutex              meQueue_;

    std::queue< FuncWrapper_ > queue_; // not thread-safe
    std::vector< std::thread > threads_; // not thread-safe

    // Stats
    std::atomic_int numWorkingThreads_ {0};

};

#endif
