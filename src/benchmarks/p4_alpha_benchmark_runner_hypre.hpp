#pragma once

#include "benchmarks/p4_alpha_benchmark_runner.hpp"

#include <mpi.h>

namespace dec3d::benchmarks {

[[nodiscard]] P4AlphaBenchmarkRunResult RunP4AlphaBenchmarkCaseDistributed(
    MPI_Comm communicator,
    const P4AlphaBenchmarkDescriptor& descriptor,
    const P4AlphaBenchmarkRunOptions& options) noexcept;

[[nodiscard]] P4AlphaBenchmarkRunResult
RunP4AlphaBenchmarkSerialDistributedComparison(
    MPI_Comm communicator,
    const P4AlphaBenchmarkDescriptor& descriptor,
    const P4AlphaBenchmarkRunOptions& options) noexcept;

}  // namespace dec3d::benchmarks
