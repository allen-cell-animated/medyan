#include "util/profiler.h"

#include "util/io/log.h"

namespace profiler {

namespace internal {
void SimpleTimerImpl< true >::report() {
    LOG(DEBUG) << "Time elapsed for " << _name << ": "
        << std::chrono::duration_cast< std::chrono::duration< double > >(_elapsed).count()
        << "s.";
}

} // namespace internal

} // namespace profiler
