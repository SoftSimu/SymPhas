
#pragma once

#include "solver.h"
#include "spsmpi.h"

#ifdef USING_MPI

namespace symphas::parallel {

template <size_t D>
grid::region_interval<D> local_iterable_domain(const domain_info<D>& dinfo) {
  len_type dims[D];
  len_type intervals[D][2];

  for (iter_type i = 0; i < D; ++i) {
    dims[i] = dinfo.global_dims[i];
  }

  intervals[0][0] = dinfo.boundary_depth + dinfo.local_x_start;
  intervals[0][1] = dinfo.boundary_depth + dinfo.local_x_end;

  if constexpr (D >= 2) {
    intervals[1][0] = dinfo.boundary_depth + dinfo.local_y_start;
    intervals[1][1] = dinfo.boundary_depth + dinfo.local_y_end;
  }

  if constexpr (D >= 3) {
    intervals[2][0] = dinfo.boundary_depth;
    intervals[2][1] = dinfo.global_dims[2] - dinfo.boundary_depth;
  }

  return grid::region_interval<D>(dims, intervals);
}

}  // namespace symphas::parallel

#endif
