#ifndef MEDYAN_Util_ThreadPool_hpp
#define MEDYAN_Util_ThreadPool_hpp

#include <atomic>
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
public:
    // An implementation of a thread-safe queue
    template< typename T >
    class SafeQueue {
    public:
        SafeQueue() = default;

        void push(const T& x) {
            std::lock_guard< std::mutex > guard(me_);
            q_.push(x);
        }
        void push(T&& x) {
            std::lock_guard< std::mutex > guard(me_);
            q_.push(std::move(x));
        }

        // This function tries to pop the queue and get the value
        // If the queue is empty, then return false
        // Else return true and the popped value is written in the parameter x
        bool tryPop(T& x) {
            std::lock_guard< std::mutex > guard(me_);
            if(q_.empty()) return false;
            x = std::move(q_.front());
            q_.pop();
            return true;
        }

    private:
        std::queue<T> q_;
        std::mutex me_;
    };

private:

    // An implementation of function wrapper to store movable objects, and to store task information
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

        queue_.push(
            [task{std::move(task)}]() mutable { task(); }
        );

        return res;
    }

private:

    // Working thread
    void work_() {
        while(!done_) {
            FuncWrapper_ f;
            if(queue_.tryPop(f)) f();

            // TODO
        }
    }

    std::atomic_bool done_ {false};

    SafeQueue< FuncWrapper_ >  queue_;
    std::vector< std::thread > threads_;

};

#endif
