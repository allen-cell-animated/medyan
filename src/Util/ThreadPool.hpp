#ifndef MEDYAN_Util_ThreadPool_hpp
#define MEDYAN_Util_ThreadPool_hpp

#include <cstddef> // size_t
#include <functional>
#include <future>
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
        // TODO
    }

    // Submit a new task
    template< typename F, typename... Args >
    auto submit(F&& f, Args&&... args) {
        using ReturnType = std::result_of_t< F(Args...) >;
        std::packaged_task< ReturnType() > task(std::forward<F>(f));
        queue_.push(
            /*[=, task{std::move(task)}] {
                return task()
            }*/
        )
        // TODO
    }

private:

    // Working thread
    void work_() {
        while(true) {
            // TODO
        }
    }

    std::queue< std::function< void() > > queue_; // TODO: use thread-safe queue
    std::vector< std::thread >            threads_;
};

#endif
