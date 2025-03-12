
/* ***************************************************************************
 * This file is part of the SymPhas library, a framework for implementing
 * solvers for phase-field problems with compile-time symbolic algebra.
 *
 * Copyright (c) 2018-2021 by Steven A. Silber and Mikko Karttunen
 *
 * SymPhas is free software, which can be redistributed or modified under
 * the terms of the GNU Lesser General Public License (LGPL) as published
 * by the Free Software Foundation; LGPL version 3, or later versions at
 * your choice.
 *
 * SymPhas is distributed with the faith that it will be helpful and
 * practical but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser
 * General Public License for more details.
 *
 * ***************************************************************************
 *
 * MODULE:  datatypes
 * PURPOSE: Defines functions for using grids, such as scaling grid values.
 *
 * ***************************************************************************
 */

#pragma once

#include "dataiterator.h"
#include "dft.h"

namespace grid {
//! Specialization based on symphas::dft::scale() using ::Grid type.
/*!
 * Scales the values of the grid type based on the number of elements
 * in the grid.
 *
 * \param grid The grid containing the values to be scaled.
 */
template <typename T>
void scale(Grid<T, 1>& grid) {
  scale(grid.values, grid.dims[0]);
}

//! Specialization based on symphas::dft::scale()
/*!
 * Scales the values of the grid type based on the number of elements
 * in the grid.
 *
 * \param grid The grid containing the values to be scaled.
 */
template <typename T>
void scale(Grid<T, 2>& grid) {
  scale(grid.values, grid.dims[0], grid.dims[1]);
}

//! Specialization based on symphas::dft::scale()
/*!
 * Scales the values of the grid type based on the number of elements
 * in the grid.
 *
 * \param grid The grid containing the values to be scaled.
 */
template <typename T>
void scale(Grid<T, 3>& grid) {
  scale(grid.values, grid.dims[0], grid.dims[1], grid.dims[2]);
}

//! Specialization based on symphas::dft::scale()
/*!
 * Scales the values of the grid type based on the number of elements
 * in the grid.
 *
 * \param grid The grid containing the values to be scaled.
 */
template <typename T>
void scale(Grid<any_vector_t<T, 1>, 1>& grid) {
  scale(grid.axis(Axis::X), grid.dims[0]);
}

//! Specialization based on symphas::dft::scale()
/*!
 * Scales the values of the grid type based on the number of elements
 * in the grid.
 *
 * \param grid The grid containing the values to be scaled.
 */
template <typename T>
void scale(Grid<any_vector_t<T, 2>, 2>& grid) {
  scale(grid.axis(Axis::X), grid.dims[0], grid.dims[1]);
  scale(grid.axis(Axis::Y), grid.dims[0], grid.dims[1]);
}

//! Specialization based on symphas::dft::scale()
/*!
 * Scales the values of the grid type based on the number of elements
 * in the grid.
 *
 * \param grid The grid containing the values to be scaled.
 */
template <typename T>
void scale(Grid<any_vector_t<T, 3>, 3>& grid) {
  scale(grid.axis(Axis::X), grid.dims[0], grid.dims[1], grid.dims[2]);
  scale(grid.axis(Axis::Y), grid.dims[0], grid.dims[1], grid.dims[2]);
  scale(grid.axis(Axis::Z), grid.dims[0], grid.dims[1], grid.dims[2]);
}

#ifdef USING_CUDA

template <typename T, size_t D>
void scale(GridCUDA<T, D>& grid) {
  scale_cuda(grid.values, grid::length<D>(grid.dims));
}

template <typename T, size_t D>
void scale(GridCUDA<any_vector_t<T, D>, D>& grid) {
  for (iter_type i = 0; i < D; ++i) {
    scale_cuda(grid.axis(symphas::index_to_axis(i)),
               grid::length<D>(grid.dims));
  }
}

#endif

}  // namespace grid

/*!
 * \addtogroup grid
 * @{
 */

namespace grid {

template <typename T, size_t D, size_t... Is>
decltype(auto) value_at(Grid<T, D> const& grid, iter_type (&pos)[D],
                        std::index_sequence<Is...>) {
  return grid(pos[Is]...);
}

template <typename T, size_t D>
decltype(auto) value_at(Grid<T, D> const& grid, iter_type (&pos)[D]) {
  return value_at(grid, pos, std::make_index_sequence<D>{});
}

template <typename T, size_t D, size_t... Is>
decltype(auto) value_at(RegionalGrid<T, D> const& grid, iter_type (&pos)[D],
                        std::index_sequence<Is...>) {
  return grid(pos[Is]...);
}

template <typename T, size_t D>
decltype(auto) value_at(RegionalGrid<T, D> const& grid, iter_type (&pos)[D]) {
  return value_at(grid, pos, std::make_index_sequence<D>{});
}

template <typename T, size_t D, size_t... Is>
decltype(auto) value_at(BoundaryGrid<T, D> const& grid, iter_type (&pos)[D],
                        std::index_sequence<Is...>) {
  return grid(pos[Is]...);
}

template <typename T, size_t D>
decltype(auto) value_at(BoundaryGrid<T, D> const& grid, iter_type (&pos)[D]) {
  return value_at(grid, pos, std::make_index_sequence<D>{});
}

#ifdef USING_CUDA

template <typename T, size_t D, size_t... Is>
decltype(auto) value_at(BoundaryGridCUDA<T, D> const& grid, iter_type (&pos)[D],
                        std::index_sequence<Is...>) {
  return grid(pos[Is]...);
}

template <typename T, size_t D>
decltype(auto) value_at(BoundaryGridCUDA<T, D> const& grid,
                        iter_type (&pos)[D]) {
  return value_at(grid, pos, std::make_index_sequence<D>{});
}

template <typename T, size_t D, size_t... Is>
decltype(auto) value_at(GridCUDA<T, D> const& grid, iter_type (&pos)[D],
                        std::index_sequence<Is...>) {
  return grid(pos[Is]...);
}

template <typename T, size_t D>
decltype(auto) value_at(GridCUDA<T, D> const& grid, iter_type (&pos)[D]) {
  return value_at(grid, pos, std::make_index_sequence<D>{});
}

template <typename T, size_t D, size_t... Is>
decltype(auto) value_at(RegionalGridCUDA<T, D> const& grid, iter_type (&pos)[D],
                        std::index_sequence<Is...>) {
  return grid(pos[Is]...);
}

template <typename T, size_t D>
decltype(auto) value_at(RegionalGridCUDA<T, D> const& grid,
                        iter_type (&pos)[D]) {
  return value_at(grid, pos, std::make_index_sequence<D>{});
}

#endif
template <typename interval_type, size_t D>
auto break_region(const iter_type (&intervals)[D][2],
                  region_interval<D> const& region) {
  if (grid::length<D>(region) == 0) {
    return std::vector<interval_type>{};
  } else if (grid::is_in_region(region.intervals, intervals)) {
    std::vector<interval_type> new_regions;
    new_regions.push_back(region.intervals);
    return new_regions;
  } else {
    std::vector<interval_type> new_regions(1);
    new_regions.reserve(fixed_pow<2, D>);

    for (iter_type i = 0; i < D; ++i) {
      new_regions.front()[i][0] = region[i][0];
      new_regions.front()[i][1] = region[i][1];
    }

    for (iter_type i = 0; i < D; ++i) {
      for (iter_type n = 0; n < new_regions.size(); ++n) {
        auto it = new_regions.begin() + n;
        if ((*it)[i][1] > intervals[i][1]) {
          if ((*it)[i][0] >= intervals[i][1]) {
            (*it)[i][0] = intervals[i][0] + ((*it)[i][0] - intervals[i][1]);
            (*it)[i][1] = intervals[i][0] + ((*it)[i][1] - intervals[i][1]);
          } else {
            interval_type add_region(*it);
            add_region[i][0] = intervals[i][0];
            add_region[i][1] = intervals[i][0] + (*it)[i][1] - intervals[i][1];

            (*it)[i][1] = intervals[i][1];
            if (grid::is_contact_overlapping(*it, add_region)) {
              grid::get_region_union(*it, *it, add_region);
            } else {
              new_regions.push_back(add_region);
            }
          }
        } else if ((*it)[i][0] < intervals[i][0]) {
          if ((*it)[i][1] <= intervals[i][0]) {
            (*it)[i][0] = intervals[i][1] - (intervals[i][0] - (*it)[i][0]);
            (*it)[i][1] = intervals[i][1] - (intervals[i][0] - (*it)[i][1]);
          } else {
            interval_type add_region(*it);
            add_region[i][1] = intervals[i][1];
            add_region[i][0] = intervals[i][1] + (*it)[i][0] - intervals[i][0];

            (*it)[i][0] = intervals[i][0];
            if (grid::is_contact_overlapping(*it, add_region)) {
              grid::get_region_union(*it, *it, add_region);
            } else {
              new_regions.push_back(add_region);
            }
          }
        }
      }
    }
    return new_regions;
  }
}

//! Creates a sub-grid centered around the given index.
/*!
 * Given a grid, the get_subgrid function will copy a small part of the Grid
 * "src" into another grid and then return that grid.
 *
 * The part that is copied is centered on the raw index n, and the given
 * len_type array "extent" indicates how far in each axis to copy around the
 * given point.
 *
 * The copy algorithm uses periodic boundaries, so sub-grids where the source
 * point is near the edges will be continued onto the opposite side.
 *
 * One restriction is that the individual extent dimensions can't be larger than
 * the grid dimensions.
 *
 * \param src The grid from which is taken the subset (subgrid).
 * \param extent The lengths beyond the center point from which to take
 * the subgrid.
 * \param n The index at which the subgrid is centered.
 */
template <template <typename, size_t> typename grid_type, typename T>
Grid<T, 1> get_subgrid(grid_type<T, 1> const& src, len_type* extent,
                       iter_type n) {
  len_type ds = extent[0] + extent[0] + 1;
  len_type s2 = extent[0];
  len_type cdims[] = {ds};

  Grid<T, 1> grid(cdims);
  T* dest = grid.values;

  if (n < s2) {
    iter_type offset = src.dims[0] - (s2 - n);
    for (iter_type i = 0; i < s2 - n; ++i) {
      dest[i] = src[offset + i];
    }

    for (iter_type i = 0; i < ds - (s2 - n); ++i) {
      dest[i + s2 - n] = src[i];
    }
  } else if (n > src.dims[0] - s2) {
    iter_type offset = n - s2;
    for (iter_type i = 0; i < src.dims[0] - (n - s2); ++i) {
      dest[i] = src[offset + i];
    }
    offset = src.dims[0] - n + s2;
    for (iter_type i = 0; i < ds - offset; ++i) {
      dest[offset + i] = src[i];
    }
  } else {
    iter_type index = n - s2;
    for (iter_type i = 0; i < ds; ++i) {
      dest[i] = src[index + i];
    }
  }

  return grid;
}

//! Specialization based on grid::get_subgrid().
template <template <typename, size_t> typename grid_type, typename T>
Grid<T, 2> get_subgrid(grid_type<T, 2> const& src, len_type* extent,
                       iter_type n) {
  len_type dsx = extent[0] + extent[0] + 1;
  len_type dsy = extent[1] + extent[1] + 1;
  len_type s2x = extent[0];
  len_type s2y = extent[1];
  len_type cdims[] = {dsx, dsy};

  Grid<T, 2> grid(cdims);
  T* dest = grid.values;

  iter_type x = n % src.dims[0] - s2x;
  iter_type y = (n / src.dims[0]) - s2y;

  iter_type index = 0;
  for (iter_type j = 0; j < dsy; ++j) {
    y %= src.dims[1];
    for (iter_type i = 0; i < dsx; ++i) {
      x %= src.dims[0];
      dest[index++] = src[x + y * src.dims[0]];
      ++x;
    }
    ++y;
  }

  return grid;
}

//! Specialization based on grid::get_subgrid().
template <template <typename, size_t> typename grid_type, typename T>
Grid<T, 3> get_subgrid(grid_type<T, 3> const& src, len_type* extent,
                       iter_type n) {
  len_type dsx = extent[0] + extent[0] + 1;
  len_type dsy = extent[1] + extent[1] + 1;
  len_type dsz = extent[2] + extent[2] + 1;
  len_type s2x = extent[0];
  len_type s2y = extent[1];
  len_type s2z = extent[2];
  len_type cdims[] = {dsx, dsy, dsz};

  Grid<T, 3> grid(cdims);
  T* dest = grid.values;

  iter_type x = n % src.dims[0] - s2x;
  iter_type y = (n / src.dims[0]) % src.dims[1] - s2y;
  iter_type z = n / (src.dims[0] * src.dims[1]) - s2z;
  iter_type index = 0;
  for (iter_type k = 0; k < dsz; ++k) {
    z %= src.dims[2];
    for (iter_type j = 0; j < dsy; ++j) {
      y %= src.dims[1];
      for (iter_type i = 0; i < dsx; ++i) {
        x %= src.dims[0];
        dest[index++] =
            src[x + y * src.dims[0] + z * src.dims[0] * src.dims[1]];
        ++x;
      }
      ++y;
    }
    ++z;
  }

  return grid;
}

template <typename T>
auto length(Block<T> const& grid) {
  return grid.len;
}

template <typename T, size_t D>
auto length(RegionalGrid<T, D> const& grid) {
  return length<D>(grid.region.dims);
}

template <typename T, size_t D>
auto length(BoundaryGrid<T, D> const& grid) {
  return length<D>(grid.dims);
}

template <typename T, size_t D>
auto length_interior(RegionalGrid<T, D> const& grid) {
  return length_interior<D>(grid.region.dims, grid.region.boundary_size);
}

template <typename T, size_t D>
auto length_interior(BoundaryGrid<T, D> const& grid) {
  return length_interior<D>(grid.dims, BOUNDARY_DEPTH);
}

template <typename T, size_t D>
inline bool is_interior_point(const iter_type (&pos)[D],
                              BoundaryGrid<T, D> const& grid) {
  return is_in_region(pos, grid.dims, BOUNDARY_DEPTH);
}

template <typename T, size_t D>
inline bool is_interior_point(iter_type n, BoundaryGrid<T, D> const& grid) {
  iter_type pos[D];
  grid::get_grid_position(pos, grid.dims, n);
  return is_interior_point(pos, grid);
}

template <typename T, size_t D, size_t... Is>
inline bool is_interior_point(const iter_type (&pos)[D],
                              RegionalGrid<T, D> const& grid,
                              std::index_sequence<Is...>) {
  iter_type pos0[]{
      ((pos[Is] >= grid.region.origin[Is] + grid.region.boundary_size)
           ? pos[Is] - grid.region.origin[Is]
           : pos[Is] - grid.region.origin[Is] +
                 (grid.dims[Is] - 2 * grid.region.boundary_size))...};
  return is_in_region(pos0, grid.region.dims, grid.region.boundary_size);
}

template <typename T, size_t D>
inline bool is_interior_point(const iter_type (&pos)[D],
                              RegionalGrid<T, D> const& grid) {
  return is_interior_point(pos, grid, std::make_index_sequence<D>{});
}

template <typename T, size_t D>
inline bool is_interior_point(iter_type n, RegionalGrid<T, D> const& grid) {
  iter_type pos[D];
  grid::get_grid_position(pos, grid.dims, n);
  return is_interior_point(pos, grid);
}

#ifdef USING_CUDA

template <typename T>
auto length(BlockCUDA<T> const& grid) {
  return grid.len;
}

template <typename T, size_t D>
auto length(RegionalGridCUDA<T, D> const& grid) {
  return length<D>(grid.region.dims);
}

template <typename T, size_t D>
auto length(BoundaryGridCUDA<T, D> const& grid) {
  return length<D>(grid.dims);
}

template <typename T, size_t D>
auto length_interior(RegionalGridCUDA<T, D> const& grid) {
  return length_interior<D>(grid.region.dims, grid.region.boundary_size);
}

template <typename T, size_t D>
auto length_interior(BoundaryGridCUDA<T, D> const& grid) {
  return length_interior<D>(grid.dims, BOUNDARY_DEPTH);
}

template <typename T, size_t D>
inline bool is_interior_point(const iter_type (&pos)[D],
                              BoundaryGridCUDA<T, D> const& grid) {
  return is_in_region(pos, grid.dims, BOUNDARY_DEPTH);
}

template <typename T, size_t D>
inline bool is_interior_point(iter_type n, BoundaryGridCUDA<T, D> const& grid) {
  iter_type pos[D];
  grid::get_grid_position(pos, grid.dims, n);
  return is_interior_point(pos, grid);
}

template <typename T, size_t D, size_t... Is>
inline bool is_interior_point(const iter_type (&pos)[D],
                              RegionalGridCUDA<T, D> const& grid,
                              std::index_sequence<Is...>) {
  iter_type pos0[]{
      ((pos[Is] >= grid.region.origin[Is] + grid.region.boundary_size)
           ? pos[Is] - grid.region.origin[Is]
           : pos[Is] - grid.region.origin[Is] +
                 (grid.dims[Is] - 2 * grid.region.boundary_size))...};
  return is_in_region(pos0, grid.region.dims, grid.region.boundary_size);
}

template <typename T, size_t D>
inline bool is_interior_point(const iter_type (&pos)[D],
                              RegionalGridCUDA<T, D> const& grid) {
  return is_interior_point(pos, grid, std::make_index_sequence<D>{});
}

template <typename T, size_t D>
inline bool is_interior_point(iter_type n, RegionalGridCUDA<T, D> const& grid) {
  iter_type pos[D];
  grid::get_grid_position(pos, grid.dims, n);
  return is_interior_point(pos, grid);
}

#endif

inline bool has_subdomain(symphas::interval_element_type const& interval) {
  return (interval.left() > interval.domain_left() ||
          interval.right() < interval.domain_right());
}

inline bool has_subdomain(symphas::grid_info const& info) {
  bool flag = false;
  for (auto const& [axis, interval] : info) {
    flag = (has_subdomain(interval)) ? true : flag;
  }
  return flag;
}

inline auto interior_dimensions(symphas::grid_info const& info) {
  grid::dim_list dims(nullptr, info.dimension());
  for (iter_type i = 0; i < info.dimension(); ++i) {
    dims[i] = info.at(symphas::index_to_axis(i)).get_interval_count();
  }
  return dims;
}

inline auto dimensions(symphas::grid_info const& info) {
  return info.get_dims();
}
}  // namespace grid

namespace grid {

#ifdef USING_CUDA

void fill_interior_cuda_1d(const scalar_t* srcHost, scalar_t* destDevice,
                           int srcLength, int destLength);

void fill_interior_cuda_2d(const scalar_t* srcHost, scalar_t* destDevice,
                           int srcWidth, int srcHeight, int destWidth,
                           int destHeight);

void fill_interior_cuda_3d(const scalar_t* srcHost, scalar_t* destDevice,
                           int srcWidth, int srcHeight, int srcDepth,
                           int destWidth, int destHeight, int destDepth);

void copy_interior_cuda_1d(const scalar_t* srcDevice, scalar_t* destHost,
                           int srcLength, int destLength);

void copy_interior_cuda_2d(const scalar_t* srcDevice, scalar_t* destHost,
                           int srcWidth, int srcHeight, int destWidth,
                           int destHeight);

void copy_interior_cuda_3d(const scalar_t* srcDevice, scalar_t* destHost,
                           int srcWidth, int srcHeight, int srcDepth,
                           int destWidth, int destHeight, int destDepth);

void fill_interior_cuda_1d(const complex_t* srcHost, complex_t* destDevice,
                           int srcLength, int destLength);

void fill_interior_cuda_2d(const complex_t* srcHost, complex_t* destDevice,
                           int srcWidth, int srcHeight, int destWidth,
                           int destHeight);

void fill_interior_cuda_3d(const complex_t* srcHost, complex_t* destDevice,
                           int srcWidth, int srcHeight, int srcDepth,
                           int destWidth, int destHeight, int destDepth);

void copy_interior_cuda_1d(const complex_t* srcDevice, complex_t* destHost,
                           int srcLength, int destLength);

void copy_interior_cuda_2d(const complex_t* srcDevice, complex_t* destHost,
                           int srcWidth, int srcHeight, int destWidth,
                           int destHeight);

void copy_interior_cuda_3d(const complex_t* srcDevice, complex_t* destHost,
                           int srcWidth, int srcHeight, int srcDepth,
                           int destWidth, int destHeight, int destDepth);

void fill_interior_cuda_1d(const any_vector_t<scalar_t, 1>* srcHost,
                           scalar_t* (&destDevice)[1], int srcLength,
                           int destLength);

void fill_interior_cuda_2d(const any_vector_t<scalar_t, 2>* srcHost,
                           scalar_t* (&destDevice)[2], int srcWidth,
                           int srcHeight, int destWidth, int destHeight);

void fill_interior_cuda_3d(const any_vector_t<scalar_t, 3>* srcHost,
                           scalar_t* (&destDevice)[3], int srcWidth,
                           int srcHeight, int srcDepth, int destWidth,
                           int destHeight, int destDepth);

void copy_interior_cuda_1d(scalar_t* const (&srcDevice)[1],
                           any_vector_t<scalar_t, 1>* destHost, int srcLength,
                           int destLength);

void copy_interior_cuda_2d(scalar_t* const (&srcDevice)[2],
                           any_vector_t<scalar_t, 2>* destHost, int srcWidth,
                           int srcHeight, int destWidth, int destHeight);

void copy_interior_cuda_3d(scalar_t* const (&srcDevice)[3],
                           any_vector_t<scalar_t, 3>* destHost, int srcWidth,
                           int srcHeight, int srcDepth, int destWidth,
                           int destHeight, int destDepth);

void fill_interior_cuda_1d(const any_vector_t<complex_t, 1>* srcHost,
                           complex_t* (&destDevice)[1], int srcLength,
                           int destLength);

void fill_interior_cuda_2d(const any_vector_t<complex_t, 2>* srcHost,
                           complex_t* (&destDevice)[2], int srcWidth,
                           int srcHeight, int destWidth, int destHeight);

void fill_interior_cuda_3d(const any_vector_t<complex_t, 3>* srcHost,
                           complex_t* (&destDevice)[3], int srcWidth,
                           int srcHeight, int srcDepth, int destWidth,
                           int destHeight, int destDepth);

void copy_interior_cuda_1d(complex_t* const (&srcDevice)[1],
                           any_vector_t<complex_t, 1>* destHost, int srcLength,
                           int destLength);

void copy_interior_cuda_2d(complex_t* const (&srcDevice)[2],
                           any_vector_t<complex_t, 2>* destHost, int srcWidth,
                           int srcHeight, int destWidth, int destHeight);

void copy_interior_cuda_3d(complex_t* const (&srcDevice)[3],
                           any_vector_t<complex_t, 3>* destHost, int srcWidth,
                           int srcHeight, int srcDepth, int destWidth,
                           int destHeight, int destDepth);

#endif
namespace {

//! Copy the interior values of the grid into an array.
/*!
 * Fills the interior values of an array using ghost cells for the boundary.
 *
 * \tparam T Type of the array.
 * \tparam D The dimensions of the array.
 */
template <typename T, size_t D>
struct fill_interior_apply;

template <typename T>
struct fill_interior_apply<T, 1> {
  //! Copy an array into the interior values of another.
  /*!
   * \param input The grid from which the sequential values are copied.
   * \param output The array into which the values are transcribed in
   * interior-based sequencing. \param dims The dimension of the interior
   * region.
   */
  template <template <typename, size_t> typename G>
  void operator()(const T* input, G<T, 1>& output, const len_type* dims) {
    ITER_GRID1(output[INDEX] = input[ENTRY], dims[0])
  }

#ifdef USING_CUDA
  //! Copy an array into the interior values of a CUDA array.
  /*!
   * \param input The grid from which the sequential values are copied.
   * \param output The array into which the values are transcribed in
   * interior-based sequencing. \param dims The dimension of the interior
   * region.
   */

  void operator()(const T* input, GridCUDA<T, 1>& output,
                  const len_type* dims) {
    fill_interior_cuda_1d(input, output.values, dims[0],
                          dims[0] + BOUNDARY_DEPTH * 2);
  }
#endif
};

template <typename T>
struct fill_interior_apply<T, 2> {
  //! Copy an array into the interior values of another.
  /*!
   * \param input The grid from which the sequential values are copied.
   * \param output The array into which the values are transcribed in
   * interior-based sequencing. \param dims The dimension of the interior
   * region.
   */
  template <template <typename, size_t> typename G>
  void operator()(const T* input, G<T, 2>& output, const len_type* dims) {
    ITER_GRID2(output[INDEX] = input[ENTRY], dims[0], dims[1])
  }

#ifdef USING_CUDA
  //! Copy an array into the interior values of a CUDA array.
  /*!
   * \param input The grid from which the sequential values are copied.
   * \param output The array into which the values are transcribed in
   * interior-based sequencing. \param dims The dimension of the interior
   * region.
   */
  void operator()(const T* input, GridCUDA<T, 1>& output,
                  const len_type* dims) {
    fill_interior_cuda_2d(input, output.values, dims[0], dims[1],
                          dims[0] + BOUNDARY_DEPTH * 2,
                          dims[1] + BOUNDARY_DEPTH * 2);
  }
#endif
};

template <typename T>
struct fill_interior_apply<T, 3> {
  //! Copy an array into the interior values of another.
  /*!
   * \param input The grid from which the sequential values are copied.
   * \param output The array into which the values are transcribed in
   * interior-based sequencing. \param dims The dimension of the interior
   * region.
   */
  template <template <typename, size_t> typename G>
  void operator()(const T* input, G<T, 3>& output, const len_type* dims) {
    ITER_GRID3(output[INDEX] = input[ENTRY], dims[0], dims[1], dims[2])
  }

#ifdef USING_CUDA
  //! Copy an array into the interior values of a CUDA array.
  /*!
   * \param input The grid from which the sequential values are copied.
   * \param output The array into which the values are transcribed in
   * interior-based sequencing. \param dims The dimension of the interior
   * region.
   */
  void operator()(const T* input, GridCUDA<T, 1>& output,
                  const len_type* dims) {
    fill_interior_cuda_3d(input, output.values, dims[0], dims[1], dims[2],
                          dims[0] + BOUNDARY_DEPTH * 2,
                          dims[1] + BOUNDARY_DEPTH * 2,
                          dims[2] + BOUNDARY_DEPTH * 2);
  }
#endif
};

template <typename T>
struct fill_interior_apply<any_vector_t<T, 1>, 1> {
  //! Copy an array into the interior values of another.
  /*!
   * \param input The grid from which the sequential values are copied.
   * \param output The array into which the values are transcribed in
   * interior-based sequencing. \param dims The dimension of the interior
   * region.
   */
  template <template <typename, size_t> typename G>
  void operator()(const any_vector_t<T, 1>* input,
                  G<any_vector_t<T, 1>, 1>& output, const len_type* dims) {
    ITER_GRID1(output.axis(Axis::X).values[INDEX] = input[ENTRY][0], dims[0])
  }
};

template <typename T>
struct fill_interior_apply<any_vector_t<T, 2>, 2> {
  //! Copy an array into the interior values of another.
  /*!
   * \param input The grid from which the sequential values are copied.
   * \param output The array into which the values are transcribed in
   * interior-based sequencing. \param dims The dimension of the interior
   * region.
   */
  template <template <typename, size_t> typename G>
  void operator()(const any_vector_t<T, 2>* input,
                  G<any_vector_t<T, 2>, 2>& output, const len_type* dims) {
    ITER_GRID2(output.axis(Axis::X).values[INDEX] = input[ENTRY][0], dims[0],
               dims[1])
    ITER_GRID2(output.axis(Axis::Y).values[INDEX] = input[ENTRY][1], dims[0],
               dims[1])
  }
};

template <typename T>
struct fill_interior_apply<any_vector_t<T, 3>, 3> {
  //! Copy an array into the interior values of another.
  /*!
   * \param input The grid from which the sequential values are copied.
   * \param output The array into which the values are transcribed in
   * interior-based sequencing. \param dims The dimension of the interior
   * region.
   */
  template <template <typename, size_t> typename G>
  void operator()(const any_vector_t<T, 3>* input,
                  G<any_vector_t<T, 3>, 3>& output, const len_type* dims) {
    ITER_GRID3(output.axis(Axis::X).values[INDEX] = input[ENTRY][0], dims[0],
               dims[1], dims[2])
    ITER_GRID3(output.axis(Axis::Y).values[INDEX] = input[ENTRY][1], dims[0],
               dims[1], dims[2])
    ITER_GRID3(output.axis(Axis::Z).values[INDEX] = input[ENTRY][2], dims[0],
               dims[1], dims[2])
  }
};

}  // namespace

//! Copies the system data into the given array.
/*!
 * The values of the system data block are copied into a new one. The copy
 * is performed point-wise for all data points.
 */
template <typename T, typename S>
void copy(Block<T> const& from, S* to) {
  std::copy(
#ifdef EXECUTION_HEADER_AVAILABLE
      std::execution::par_unseq,
#endif
      from.values, from.values + from.len, to);
}

//! Copies the system data into the given array.
/*!
 * The values of the system data block are copied into a new one. The copy
 * is performed point-wise for all data points.
 */
template <typename T, typename S>
void copy(Block<T> const& from, Block<S>& to) {
  copy(from, to.values);
}

//! Copies the system data into the given array.
/*!
 * The values of the system data block are copied into a new one. The copy
 * is performed point-wise for all data points.
 */
template <typename T, typename S>
void copy(MultiBlock<1, T> const& from, any_vector_t<S, 1>* to) {
  std::for_each(
#ifdef EXECUTION_HEADER_AVAILABLE
      std::execution::par_unseq,
#endif
      to, to + from.len, [&](auto& e) {
        size_t i = &e - to;
        e = any_vector_t<T, 1>{from(Axis::X)[i]};
      });
}

//! Copies the system data into the given array.
/*!
 * The values of the system data block are copied into a new one. The copy
 * is performed point-wise for all data points.
 */
template <typename T, typename S>
void copy(MultiBlock<2, T> const& from, any_vector_t<S, 2>* to) {
  std::for_each(
#ifdef EXECUTION_HEADER_AVAILABLE
      std::execution::par_unseq,
#endif
      to, to + from.len, [&](auto& e) {
        size_t i = &e - to;
        e = any_vector_t<T, 2>{from(Axis::X)[i], from(Axis::Y)[i]};
      });
}

//! Copies the system data into the given array.
/*!
 * The values of the system data block are copied into a new one. The copy
 * is performed point-wise for all data points.
 */
template <typename T, typename S>
void copy(MultiBlock<3, T> const& from, any_vector_t<S, 3>* to) {
  std::for_each(
#ifdef EXECUTION_HEADER_AVAILABLE
      std::execution::par_unseq,
#endif
      to, to + from.len, [&](auto& e) {
        size_t i = &e - to;
        e = any_vector_t<S, 3>{from(Axis::X)[i], from(Axis::Y)[i],
                               from(Axis::Z)[i]};
      });
}

//! Fills the block-type data with values from the array.
/*!
 * The values of the system data block are initialized from the
 * given values, correctly transcribing all values.
 */
template <typename S, typename T>
void copy(const S* from, Block<T>& to) {
  std::copy(
#ifdef EXECUTION_HEADER_AVAILABLE
      std::execution::par_unseq,
#endif
      from, from + to.len, to.values);
}

//! Fills the block-type data with values from the array.
/*!
 * The values of the system data block are copied into a new one. The copy
 * is performed point-wise for all data points.
 */
template <typename S, typename T>
void copy(const any_vector_t<S, 1>* from, MultiBlock<1, T>& to) {
  std::for_each(
#ifdef EXECUTION_HEADER_AVAILABLE
      std::execution::par_unseq,
#endif
      from, from + to.len, [&](auto& e) {
        size_t i = &e - from;
        to(Axis::X)[i] = e[0];
      });
}

//! Fills the block-type data with values from the array.
/*!
 * The values of the system data block are copied into a new one. The copy
 * is performed point-wise for all data points.
 */
template <typename S, typename T>
void copy(const any_vector_t<S, 2>* from, MultiBlock<2, T>& to) {
  std::for_each(
#ifdef EXECUTION_HEADER_AVAILABLE
      std::execution::par_unseq,
#endif
      from, from + to.len, [&](auto& e) {
        size_t i = &e - from;
        to(Axis::X)[i] = e[0];
        to(Axis::Y)[i] = e[1];
      });
}

//! Copies the system data into the given array.
/*!
 * The values of the system data block are copied into a new one. The copy
 * is performed point-wise for all data points.
 */
template <typename S, typename T>
void copy(const any_vector_t<S, 3>* from, MultiBlock<3, T>& to) {
  std::for_each(
#ifdef EXECUTION_HEADER_AVAILABLE
      std::execution::par_unseq,
#endif
      from, from + to.len, [&](auto& e) {
        size_t i = &e - from;
        to(Axis::X)[i] = e[0];
        to(Axis::Y)[i] = e[1];
        to(Axis::Z)[i] = e[2];
      });
}

//! Fills the block-type data with values from the array.
/*!
 * The values of the system data block are initialized from the
 * given values, correctly transcribing all values.
 */
template <typename S, typename T, size_t D>
void copy(const S* from, RegionalGrid<T, D>& to) {
  auto intervals = grid::get_iterable_domain(to);
  symphas::data_iterator_region from_it(from, intervals);
  symphas::data_iterator_region to_it(to, intervals);
  auto len = grid::length(intervals);

  std::copy(
#ifdef EXECUTION_HEADER_AVAILABLE
      std::execution::par_unseq,
#endif
      from_it, from_it + len, to_it);
}

//! Copies the system data into the given array.
/*!
 * The values of the system data block are copied into a new one. The copy
 * is performed point-wise for all data points.
 */
template <typename T, typename S, size_t D>
void copy(RegionalGrid<T, D> const& from, RegionalGrid<S, D>& to) {
  auto intervals = grid::get_iterable_domain(from);
  symphas::data_iterator_region from_it(from, intervals);
  symphas::data_iterator_region to_it(to, intervals);
  auto len = grid::length(intervals);

  std::copy(
#ifdef EXECUTION_HEADER_AVAILABLE
      std::execution::par_unseq,
#endif
      from_it, from_it + len, to_it);
}

//! Copies the system data into the given array.
/*!
 * The values of the grid are filled with their index.
 */
template <typename T>
void fill(T* data, len_type len, T value) {
  std::fill(
#ifdef EXECUTION_HEADER_AVAILABLE
      std::execution::par_unseq,
#endif
      data, data + len, value);
}

//! Copies the system data into the given array.
/*!
 * The values of the grid are filled with their index.
 */
template <typename T, typename S>
void fill(Block<T>& grid, S value) {
  std::fill(
#ifdef EXECUTION_HEADER_AVAILABLE
      std::execution::par_unseq,
#endif
      grid.values, grid.values + grid.len, value);
}

//! Copies the system data into the given array.
/*!
 * The values of the grid are filled with their index.
 */
template <typename T, size_t N, typename S>
void fill(MultiBlock<N, T>& grid, any_vector_t<S, N> value) {
  std::fill(
#ifdef EXECUTION_HEADER_AVAILABLE
      std::execution::par_unseq,
#endif
      grid.values, grid.values + grid.len, value);
}

//! Copies the system data into the given array.
/*!
 * The values of the system data block are copied into a new one. The copy
 * is performed point-wise for all data points.
 */
template <typename T>
void fill_random(T* data, len_type len, scalar_t min = 0, scalar_t max = 1) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> dist(min, max);

  std::for_each(
#ifdef EXECUTION_HEADER_AVAILABLE
      std::execution::par_unseq,
#endif
      data, data + len, [&](auto& e) { e = dist(gen); });
}

//! Copies the system data into the given array.
/*!
 * The values of the system data block are copied into a new one. The copy
 * is performed point-wise for all data points.
 */
template <typename T>
void fill_random(Block<T>& grid, scalar_t min = 0, scalar_t max = 1) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> dist(min, max);

  std::for_each(
#ifdef EXECUTION_HEADER_AVAILABLE
      std::execution::par_unseq,
#endif
      grid.values, grid.values + grid.len, [&](auto& e) { e = dist(gen); });
}

//! Copies the system data into the given array.
/*!
 * The values of the system data block are copied into a new one. The copy
 * is performed point-wise for all data points.
 */
inline void fill_random(Block<complex_t>& grid, scalar_t min = 0,
                        scalar_t max = 1) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> dist(min, max);

  std::for_each(
#ifdef EXECUTION_HEADER_AVAILABLE
      std::execution::par_unseq,
#endif
      grid.values, grid.values + grid.len,
      [&](auto& e) { e = {dist(gen), dist(gen)}; });
}

//! Copies the system data into the given array.
/*!
 * The values of the system data block are copied into a new one. The copy
 * is performed point-wise for all data points.
 */
template <typename T, size_t N>
void fill_random(MultiBlock<N, T>& grid, scalar_t min = 0, scalar_t max = 1) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> dist(min, max);

  for (iter_type i = 0; i < N; ++i) {
    std::for_each(
#ifdef EXECUTION_HEADER_AVAILABLE
        std::execution::par_unseq,
#endif
        grid.values[i], grid.values[i] + grid.len,
        [&](auto& e) { e = dist(gen); });
  }
}

//! Copies the system data into the given array.
/*!
 * The values of the region within the RegionalGrid are initialized to random
 * values in the given range.
 */
template <typename T, size_t D>
void fill_random(RegionalGrid<T, D>& grid, scalar_t min = 0, scalar_t max = 1) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> dist(min, max);

  auto intervals = grid::get_iterable_domain(grid);
  symphas::data_iterator_region it(grid, intervals);
  auto len = grid::length(intervals);

  std::for_each(
#ifdef EXECUTION_HEADER_AVAILABLE
      std::execution::par_unseq,
#endif
      it, it + len, [&](auto e) { e = dist(gen); });
}

//! Copies the system data into the given array.
/*!
 * The values of the system data block are copied into a new one. The copy
 * is performed point-wise for all data points.
 */
template <size_t D>
inline void fill_random(RegionalGrid<complex_t, D>& grid, scalar_t min = 0,
                        scalar_t max = 1) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> dist(min, max);

  auto intervals = grid::get_iterable_domain(grid);
  symphas::data_iterator_region it(grid, intervals);
  auto len = grid::length(intervals);

  std::for_each(
#ifdef EXECUTION_HEADER_AVAILABLE
      std::execution::par_unseq,
#endif
      it, it + len, [&](auto e) { e = {dist(gen), dist(gen)}; });
}

//! Copies the system data into the given array.
/*!
 * The values of the system data block are copied into a new one. The copy
 * is performed point-wise for all data points.
 */
template <typename T, size_t D>
void fill_random(RegionalGrid<any_vector_t<T, D>, D>& grid, scalar_t min = 0,
                 scalar_t max = 1) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> dist(min, max);

  auto intervals = grid::get_iterable_domain(grid);
  symphas::data_iterator_region it(grid, intervals);
  auto len = grid::length(intervals);

  for (iter_type i = 0; i < D; ++i) {
    std::for_each(
#ifdef EXECUTION_HEADER_AVAILABLE
        std::execution::par_unseq,
#endif
        it, it + len, [&](auto e) {
          for (iter_type i = 0; i < D; ++i) {
            e[i] = dist(gen);
          }
        });
  }
}

//! Copies the system data into the given array.
/*!
 * The values of the grid are filled with their index.
 */
template <typename T>
void fill_n(T* data, len_type len) {
  std::for_each(
#ifdef EXECUTION_HEADER_AVAILABLE
      std::execution::par_unseq,
#endif
      data, data + len, [&](auto& e) { e = T(&e - data); });
}

//! Copies the system data into the given array.
/*!
 * The values of the grid are filled with their index.
 */
template <typename T>
void fill_n(Block<T>& grid) {
  std::for_each(
#ifdef EXECUTION_HEADER_AVAILABLE
      std::execution::par_unseq,
#endif
      grid.values, grid.values + grid.len,
      [&](auto& e) { e = T(&e - grid.values); });
}

//! Copies the system data into the given array.
/*!
 * The values of the grid are filled with their index.
 */
template <typename T, size_t N>
void fill_n(MultiBlock<N, T>& grid) {
  for (iter_type i = 0; i < N; ++i) {
    std::for_each(
#ifdef EXECUTION_HEADER_AVAILABLE
        std::execution::par_unseq,
#endif
        grid.values[i], grid.values[i] + grid.len,
        [&](auto& e) { e = T(&e - grid.values[i]); });
  }
}

//! Copy the interior values of the grid into an array.
/*!
 * The interior values of the given grid are copied into an array.
 * It is assumed that the array has enough space to store all the interior
 * values. For computing the number of interior points, see
 * grid::length_interior(len_type const*). The grid is 1-dimensional.
 *
 * \param input The grid from which the interior values are copied.
 * \param output The array into which the values are copied.
 */
template <typename T>
void copy_interior(Grid<T, 1> const& input, T* output) {
  ITER_GRID1(output[ENTRY] = input.values[INDEX], input.dims[0])
}

//! Copy the interior values of the grid into an array.
/*!
 * Implementation of copying interior values for a 2-dimensional Grid, see
 * grid::copy_interior(Grid<T, 1> const&, T*).
 *
 * \param input The grid from which the interior values are copied.
 * \param output The array into which the values are copied.
 */
template <typename T>
void copy_interior(Grid<T, 2> const& input, T* output) {
  ITER_GRID2(output[ENTRY] = input.values[INDEX], input.dims[0], input.dims[1]);
}

//! Copy the interior values of the grid into an array.
/*!
 * Implementation of copying interior values for a 3-dimensional Grid, see
 * grid::copy_interior(Grid<T, 1> const&, T*).
 *
 * \param input The grid from which the interior values are copied.
 * \param output The array into which the values are copied.
 */
template <typename T>
void copy_interior(Grid<T, 3> const& input, T* output) {
  ITER_GRID3(output[ENTRY] = input.values[INDEX], input.dims[0], input.dims[1],
             input.dims[2]);
}

//! Copy the interior values of the grid into an array.
/*!
 * The interior values of the given grid are copied into an array.
 * It is assumed that the array has enough space to store all the interior
 * values. For computing the number of interior points, see
 * grid::length_interior(len_type const*). The grid is 1-dimensional.
 *
 * \param input The grid from which the interior values are copied.
 * \param output The array into which the values are copied.
 */
template <typename T>
void copy_interior(Grid<any_vector_t<T, 1>, 1> const& input,
                   any_vector_t<T, 1>* output) {
  ITER_GRID1(output[ENTRY][0] = input.axis(Axis::X)[INDEX], input.dims[0])
}

//! Copy the interior values of the grid into an array.
/*!
 * Implementation of copying interior values for a 2-dimensional Grid, see
 * grid::copy_interior(Grid<T, 1> const&, T*).
 *
 * \param input The grid from which the interior values are copied.
 * \param output The array into which the values are copied.
 */
template <typename T>
void copy_interior(Grid<any_vector_t<T, 2>, 2> const& input,
                   any_vector_t<T, 2>* output) {
  ITER_GRID2(output[ENTRY][0] = input.axis(Axis::X)[INDEX], input.dims[0],
             input.dims[1]);
  ITER_GRID2(output[ENTRY][1] = input.axis(Axis::Y)[INDEX], input.dims[0],
             input.dims[1]);
}

//! Copy the interior values of the grid into an array.
/*!
 * Implementation of copying interior values for a 3-dimensional Grid, see
 * grid::copy_interior(Grid<T, 1> const&, T*).
 *
 * \param input The grid from which the interior values are copied.
 * \param output The array into which the values are copied.
 */
template <typename T>
void copy_interior(Grid<any_vector_t<T, 3>, 3> const& input,
                   any_vector_t<T, 3>* output) {
  ITER_GRID3(output[ENTRY][0] = input.axis(Axis::X)[INDEX], input.dims[0],
             input.dims[1], input.dims[2]);
  ITER_GRID3(output[ENTRY][1] = input.axis(Axis::Y)[INDEX], input.dims[0],
             input.dims[1], input.dims[2]);
  ITER_GRID3(output[ENTRY][2] = input.axis(Axis::Z)[INDEX], input.dims[0],
             input.dims[1], input.dims[2]);
}

//! Copy the interior values of the regional grid into an array.
/*!
 * The interior values of the given grid are copied into an array.
 * It is assumed that the array has enough space to store all the interior
 * values. For computing the number of interior points, see
 * grid::length_interior(len_type const*). The grid is 1-dimensional.
 *
 * \param input The grid from which the interior values are copied.
 * \param output The array into which the values are copied.
 */
template <typename T>
void copy_interior(RegionalGrid<T, 1> const& input, T* output) {
  ITER_GRID1(output[ENTRY] = input[INDEX], input.dims[0])
}

//! Copy the interior values of the grid into an array.
/*!
 * Implementation of copying interior values for a 2-dimensional Grid, see
 * grid::copy_interior(RegionalGrid<T, 1> const&, T*).
 *
 * \param input The grid from which the interior values are copied.
 * \param output The array into which the values are copied.
 */
template <typename T>
void copy_interior(RegionalGrid<T, 2> const& input, T* output) {
  ITER_GRID2(output[ENTRY] = input[INDEX], input.dims[0], input.dims[1]);
}

//! Copy the interior values of the grid into an array.
/*!
 * Implementation of copying interior values for a 3-dimensional Grid, see
 * grid::copy_interior(RegionalGrid<T, 1> const&, T*).
 *
 * \param input The grid from which the interior values are copied.
 * \param output The array into which the values are copied.
 */
template <typename T>
void copy_interior(RegionalGrid<T, 3> const& input, T* output) {
  ITER_GRID3(output[ENTRY] = input[INDEX], input.dims[0], input.dims[1],
             input.dims[2]);
}

//! Copy the interior values of the grid into an array.
/*!
 * The interior values of the given grid are copied into an array.
 * It is assumed that the array has enough space to store all the interior
 * values. For computing the number of interior points, see
 * grid::length_interior(len_type const*). The grid is 1-dimensional.
 *
 * \param input The grid from which the interior values are copied.
 * \param output The array into which the values are copied.
 */
template <typename T>
void copy_interior(RegionalGrid<any_vector_t<T, 1>, 1> const& input,
                   any_vector_t<T, 1>* output) {
  ITER_GRID1(output[ENTRY][0] = input.axis(Axis::X)[INDEX], input.dims[0])
}

//! Copy the interior values of the grid into an array.
/*!
 * Implementation of copying interior values for a 2-dimensional Grid, see
 * grid::copy_interior(RegionalGrid<T, 1> const&, T*).
 *
 * \param input The grid from which the interior values are copied.
 * \param output The array into which the values are copied.
 */
template <typename T>
void copy_interior(RegionalGrid<any_vector_t<T, 2>, 2> const& input,
                   any_vector_t<T, 2>* output) {
  ITER_GRID2(output[ENTRY][0] = input.axis(Axis::X)[INDEX], input.dims[0],
             input.dims[1]);
  ITER_GRID2(output[ENTRY][1] = input.axis(Axis::Y)[INDEX], input.dims[0],
             input.dims[1]);
}

//! Copy the interior values of the grid into an array.
/*!
 * Implementation of copying interior values for a 3-dimensional Grid, see
 * grid::copy_interior(RegionalGrid<T, 1> const&, T*).
 *
 * \param input The grid from which the interior values are copied.
 * \param output The array into which the values are copied.
 */
template <typename T>
void copy_interior(RegionalGrid<any_vector_t<T, 3>, 3> const& input,
                   any_vector_t<T, 3>* output) {
  ITER_GRID3(output[ENTRY][0] = input.axis(Axis::X)[INDEX], input.dims[0],
             input.dims[1], input.dims[2]);
  ITER_GRID3(output[ENTRY][1] = input.axis(Axis::Y)[INDEX], input.dims[0],
             input.dims[1], input.dims[2]);
  ITER_GRID3(output[ENTRY][2] = input.axis(Axis::Z)[INDEX], input.dims[0],
             input.dims[1], input.dims[2]);
}

//! Copy the interior values of the grid into an array.
/*!
 * Fills the interior values of an array using ghost cells for the boundary.
 *
 * \param input The grid from which the sequential values are copied.
 * \param output The array into which the values are transcribed in
 * interior-based sequencing. \param dims The dimension of the interior region.
 */
template <typename T, size_t D>
void fill_interior(const T* input, Grid<T, D>& output,
                   const len_type (&dims)[D]) {
  fill_interior_apply<T, D>{}(input, output, dims);
}

//! Copy the interior values of the grid into an array.
/*!
 * Fills the interior values of an array using ghost cells for the boundary.
 *
 * \param input The grid from which the sequential values are copied.
 * \param output The array into which the values are transcribed in
 * interior-based sequencing. \param dims The dimension of the interior region.
 */
template <typename T, size_t D>
void fill_interior(const T* input, RegionalGrid<T, D>& output,
                   const len_type (&dims)[D]) {
  fill_interior_apply<T, D>{}(input, output, dims);
}

#ifdef USING_CUDA
//! Copy the interior values of the grid into an array.
/*!
 * Fills the interior values of an array using ghost cells for the boundary.
 *
 * \param input The grid from which the sequential values are copied.
 * \param output The array into which the values are transcribed in
 * interior-based sequencing. \param dims The dimension of the interior region.
 */
template <typename T, size_t D>
void fill_interior(const T* input, GridCUDA<T, D>& output,
                   const len_type (&dims)[D]) {
  fill_interior_apply<T, D>{}(input, output, dims);
}

//! Copy the interior values of the grid into an array.
/*!
 * Fills the interior values of an array using ghost cells for the boundary.
 *
 * \param input The grid from which the sequential values are copied.
 * \param output The array into which the values are transcribed in
 * interior-based sequencing. \param dims The dimension of the interior region.
 */
template <typename T, size_t D>
void fill_interior(const T* input, RegionalGridCUDA<T, D>& output,
                   const len_type (&dims)[D]) {
  fill_interior_apply<T, D>{}(input, output, dims);
}
#endif

template <typename T>
__host__ __device__ bool compare_cutoff(T const& left, T const& right) {
  return left >= right;
}

template <typename T>
__host__ __device__ bool compare_cutoff(const T* left, iter_type index,
                                        T const& right) {
  return compare_cutoff(left[index], right);
}

template <typename T, size_t D>
__host__ __device__ bool compare_cutoff(T* const (&left)[D], iter_type index,
                                        const T (&right)[D]) {
  T result{};
  for (iter_type i = 0; i < D; ++i) {
    result += left[i][index] * left[i][index];
  }
  using std::abs;
  using std::sqrt;
  using symphas::math::abs;
  using symphas::math::sqrt;
  return compare_cutoff(sqrt(result), abs(right));
}

template <typename T, size_t D>
__host__ __device__ bool compare_cutoff(const T (&left)[D],
                                        any_vector_t<T, D> const& right) {
  T result{};
  for (iter_type i = 0; i < D; ++i) {
    result += left[i] * left[i];
  }
  using std::sqrt;
  using symphas::math::abs;
  using symphas::math::sqrt;
  return compare_cutoff(sqrt(result), abs(right));
}

template <typename T, size_t D>
__host__ __device__ bool compare_cutoff(T* const (&left)[D], iter_type index,
                                        any_vector_t<T, D> const& right) {
  T result{};
  for (iter_type i = 0; i < D; ++i) {
    result += left[i][index] * left[i][index];
  }
  using std::sqrt;
  using symphas::math::abs;
  using symphas::math::sqrt;
  return compare_cutoff(sqrt(result), abs(right));
}

template <size_t D>
struct RegionAdjustParams {
  len_type new_interior_dims[D];
  len_type old_interior_dims[D];
  iter_type new_stride[D];
  iter_type old_stride[D];
  len_type wrap_dims[D];
  len_type origin_delta[D];

  RegionAdjustParams(const iter_type (&new_origin)[D],
                     const len_type (&new_dims)[D],
                     const iter_type (&old_origin)[D],
                     const len_type (&old_dims)[D],
                     const len_type (&global_dims)[D], len_type boundary_size)
      : new_interior_dims{},
        old_interior_dims{},
        new_stride{},
        old_stride{},
        wrap_dims{},
        origin_delta{} {
    init_region_stride(new_stride, new_dims);
    init_region_stride(old_stride, old_dims);

    for (iter_type i = 0; i < D; ++i) {
      wrap_dims[i] = global_dims[i] - 2 * boundary_size;

      len_type old_origin_normalized_i =
          old_origin[i] -
          wrap_dims[i] * std::lround(double(old_origin[i]) / wrap_dims[i]);
      len_type new_origin_normalized_i =
          new_origin[i] -
          wrap_dims[i] * std::lround(double(new_origin[i]) / wrap_dims[i]);
      origin_delta[i] = new_origin_normalized_i - old_origin_normalized_i;

      new_interior_dims[i] = new_dims[i] - 2 * boundary_size;
      old_interior_dims[i] = old_dims[i] - 2 * boundary_size;
    }
  }

  __device__ __host__ void compute_old_position(iter_type (&pos)[D],
                                                iter_type n) const {
    for (iter_type i = 0; i < D; ++i) {
      pos[i] += origin_delta[i];
      pos[i] = (pos[i] < 0)               ? pos[i] + wrap_dims[i]
               : (pos[i] >= wrap_dims[i]) ? pos[i] - wrap_dims[i]
                                          : pos[i];
    }
  }
};

template <size_t D>
struct OriginAdjustParams {
  len_type overlap_dims[D];
  len_type interior_dims[D];
  len_type shift_pos[D];
  iter_type offset;
  iter_type boundary_offset;

  OriginAdjustParams(const iter_type (&new_origin)[D],
                     const iter_type (&old_origin)[D],
                     const len_type (&dims)[D], const len_type (&direction)[D],
                     len_type boundary_size)
      : overlap_dims{},
        interior_dims{},
        shift_pos{},
        offset{},
        boundary_offset{} {
    iter_type stride[D]{};
    init_region_stride(stride, dims);

    for (iter_type i = 0; i < D; ++i) {
      interior_dims[i] = dims[i] - boundary_size * 2;
      shift_pos[i] = new_origin[i] - old_origin[i];
      overlap_dims[i] = interior_dims[i] - std::abs(shift_pos[i]);
    }

    // iter_type pos[D]{};
    offset = grid::index_from_position(shift_pos, stride);

    for (iter_type i = 0; i < D; ++i) {
      // if (direction[i] < 0) {
      //   pos[i] = interior_dims[i] - 1;
      // }

      boundary_offset += stride[i] * boundary_size;
    }
  }
};

template <size_t D>
struct MinimalRegionParams {
  iter_type stride[D];
  iter_type dims[D];
  MinimalRegionParams(const len_type (&dims)[D]) : stride{} {
    for (iter_type d = 0; d < D; ++d) {
      this->dims[d] = dims[d];
    }
    init_region_stride(stride, dims);
  }
};

//! Adjust the floating region inside the grid to a new position and
//! dimensions
/*!
 * Adjust the floating region inside the grid to a new position and
 * dimensions.
 *
 * \param new_values The new region to fit the old region inside. This needs
 * to have as much space as the given dimensions with an added boundary
 * layer.
 *
 * \param new_origin The new starting point of the new region inside the
 * global domain. This origin point includes the boundary layer as well,
 * meaning that the interior data starts after this origin point.
 *
 * \param new_dimensions The dimensions of the new region, also counting the
 * boundary layer.
 *
 * \param old_values The old region which contains data, from which the data
 * that is overlapping with the new region will be copied.
 *
 * \param old_origin The starting point of the old region inside the global
 * domain. The origin point includes the boundary layer as well, meaning
 * that the interior data starts after this origin point.
 *
 * \param old_dims The dimensions of the old region, also counting the
 * boundary.
 *
 * \param global_dims The dimensions of the global domain.
 *
 * \param boundary_size The width of the boundary region inside both the old
 * and new regions.
 */
template <typename T, size_t D>
void adjust_region_to_from(T* new_values, const iter_type (&new_origin)[D],
                           const len_type (&new_dims)[D], const T* old_values,
                           const iter_type (&old_origin)[D],
                           const len_type (&old_dims)[D],
                           const len_type (&global_dims)[D], const T empty,
                           len_type boundary_size) {
  RegionAdjustParams<D> rp{new_origin, new_dims,    old_origin,
                           old_dims,   global_dims, boundary_size};

  for (iter_type n = 0; n < grid::length<D>(rp.new_interior_dims); ++n) {
    iter_type pos[D];

    grid::get_grid_position(pos, rp.new_interior_dims, n);
    iter_type new_index =
        index_from_position(pos, rp.new_stride, boundary_size);
    rp.compute_old_position(pos, n);

    bool is_in_old = grid::is_in_region(pos, rp.old_interior_dims);
    if (is_in_old) {
      iter_type old_index =
          index_from_position(pos, rp.old_stride, boundary_size);
      new_values[new_index] = old_values[old_index];
    }
  }
}

//! Adjust the floating region inside the grid to a new position and
//! dimensions
/*!
 * Adjust the floating region inside the grid to a new position and
 * dimensions.
 *
 * \param values The existing region which will have its values moved in
 * place to accommodate a new origin.
 *
 * \param new_origin The new starting point of the new region inside the
 * global domain. This origin point includes the boundary layer as well,
 * meaning that the interior data starts after this origin point.
 *
 * \param old_origin The starting point of the old region inside the global
 * domain. The origin point includes the boundary layer as well, meaning
 * that the interior data starts after this origin point.
 *
 * \param dims The dimensions of the existing region, counting the boundary
 * layer.
 *
 * \param global_dims The dimensions of the global domain.
 *
 * \param boundary_size The width of the boundary region inside both the old
 * and new regions.
 */
template <typename T, size_t D>
void adjust_origin_to_from_replace(T*(&values),
                                   const iter_type (&new_origin)[D],
                                   const iter_type (&old_origin)[D],
                                   const len_type (&dims)[D],
                                   const len_type (&global_dims)[D], T empty,
                                   len_type boundary_size) {
  T* new_values = new T[grid::length<D>(dims)]{};
  std::fill(new_values, new_values + grid::length<D>(dims), empty);

  RegionAdjustParams<D> rp{new_origin, dims,        old_origin,
                           dims,       global_dims, boundary_size};

  for (iter_type n = 0; n < grid::length<D>(rp.new_interior_dims); ++n) {
    iter_type pos[D];

    grid::get_grid_position(pos, rp.new_interior_dims, n);
    iter_type new_index =
        index_from_position(pos, rp.new_stride, boundary_size);
    rp.compute_old_position(pos, n);

    bool is_in_old = grid::is_in_region(pos, rp.new_interior_dims);
    if (is_in_old) {
      iter_type old_index =
          index_from_position(pos, rp.new_stride, boundary_size);
      new_values[new_index] = values[old_index];
    }
  }
  std::swap(new_values, values);
  delete[] new_values;
}

//! Adjust the floating region inside the grid to a new position and
//! dimensions
/*!
 * Adjust the floating region inside the grid to a new position and
 * dimensions.
 *
 * \param values The existing region which will have its values moved in
 * place to accommodate a new origin.
 *
 * \param new_origin The new starting point of the new region inside the
 * global domain. This origin point includes the boundary layer as well,
 * meaning that the interior data starts after this origin point.
 *
 * \param old_origin The starting point of the old region inside the global
 * domain. The origin point includes the boundary layer as well, meaning
 * that the interior data starts after this origin point.
 *
 * \param dims The dimensions of the existing region, counting the boundary
 * layer.
 *
 * \param global_dims The dimensions of the global domain.
 *
 * \param boundary_size The width of the boundary region inside both the old
 * and new regions.
 */
template <typename T, size_t D>
void adjust_origin_to_from(T*(&values), const iter_type (&new_origin)[D],
                           const iter_type (&old_origin)[D],
                           const len_type (&dims)[D],
                           const len_type (&global_dims)[D], T empty,
                           len_type boundary_size) {
  iter_type direction[D];

  bool skip = true;
  for (iter_type i = 0; i < D; ++i) {
    iter_type delta = new_origin[i] - old_origin[i];
    direction[i] = (delta >= 0) ? 1 : -1;
    skip = skip && (delta == 0);
  }

  if (skip) {
    return;
  }

  for (iter_type i = 0; i < D; ++i) {
    iter_type delta = std::abs(new_origin[i] - old_origin[i]);
    skip = skip || (delta + dims[i] > global_dims[i]);
  }

  if (skip) {
    adjust_origin_to_from_replace(values, new_origin, old_origin, dims,
                                  global_dims, empty, boundary_size);
  } else {
    OriginAdjustParams<D> op{new_origin, old_origin, dims, direction,
                             boundary_size};

    iter_type stride[D]{};
    init_region_stride(stride, dims);

    iter_type overlap_dims[D]{};
    iter_type interior_dims[D]{};
    iter_type shift_pos[D]{};
    for (iter_type i = 0; i < D; ++i) {
      interior_dims[i] = dims[i] - boundary_size * 2;
      shift_pos[i] = new_origin[i] - old_origin[i];
      overlap_dims[i] = interior_dims[i] - std::abs(shift_pos[i]);
    }

    iter_type pos[D]{};
    iter_type offset = grid::index_from_position(shift_pos, stride);
    iter_type boundary_offset = 0;

    for (iter_type i = 0; i < D; ++i) {
      if (direction[i] < 0) {
        pos[i] = interior_dims[i] - 1;
      }

      boundary_offset += stride[i] * boundary_size;
    }

    for (iter_type n = 0; n < grid::length<D>(interior_dims); ++n) {
      bool set_empty = false;
      for (iter_type i = 0; i < D; ++i) {
        if (pos[i] >= interior_dims[i]) {
          pos[i + 1] += direction[i + 1];
          pos[i] = 0;
        } else if (pos[i] < 0) {
          pos[i + 1] += direction[i + 1];
          pos[i] = interior_dims[i] - 1;
        }
        set_empty = (direction[i] > 0)
                        ? (pos[i] >= overlap_dims[i])
                        : (pos[i] < interior_dims[i] - overlap_dims[i]);
      }

      iter_type index = grid::index_from_position(pos, stride);
      if (set_empty) {
        values[index + boundary_offset] = empty;
      } else {
        values[index + boundary_offset] =
            values[index + offset + boundary_offset];
      }

      pos[0] += direction[0];
    }
  }
}

//! Construct a new view around a region with values greater than a cutoff.
/*!
 * Given the regional grid, populate the origin and dimensions for the
 * minimal new region where all values outside are less than cutoff. The
 * regional grid dimensions give the total region size including boundaries,
 * and this algorithm will compute the dimensions without extending them to
 * boundaries.
 *
 * \param grid The regional grid from which to find a new view.
 * \param origin The values of the new region origin are populated here.
 * \param dims The values of the dimensions of the new region are populated
 * here. \param boundary_size The width of the boundary used in the
 * regional_grid.
 */
template <typename T, size_t D>
void get_view_resized(RegionalGrid<T, D>& grid, T cutoff_value,
                      iter_type (&origin)[D], len_type (&dims)[D]) {
  iter_type intervals[D][2];
  iter_type offset[D]{};
  iter_type stride[D]{};

  for (iter_type i = 0; i < D; ++i) {
    intervals[i][0] = grid.region.dims[i];
    intervals[i][1] = 0;

    dims[i] = grid.region.dims[i] - grid.region.boundary_size * 2;
    offset[i] = grid.region.boundary_size;
  }
  init_region_stride(stride, grid.region.dims);

  omp_lock_t interval_lock;
  omp_init_lock(&interval_lock);

#pragma omp parallel for
  for (iter_type n = 0; n < grid::length<D>(dims); ++n) {
    iter_type pos[D];
    get_grid_position_offset(pos, dims, offset, n);
    iter_type index = index_from_position(pos, stride);

    using std::abs;
    if (compare_cutoff(grid.values, index, cutoff_value) &&
        !grid::is_in_region(pos, intervals)) {
      omp_set_lock(&interval_lock);
      for (iter_type i = 0; i < D; ++i) {
        if (pos[i] >= intervals[i][0]) {
          intervals[i][1] = std::max(pos[i], intervals[i][1]);
        } else {
          intervals[i][0] = std::min(pos[i], intervals[i][0]);
        }
      }
      omp_unset_lock(&interval_lock);
    }
  }
  omp_destroy_lock(&interval_lock);

  for (iter_type i = 0; i < D; ++i) {
    origin[i] =
        intervals[i][0] - grid.region.boundary_size + grid.region.origin[i];
    dims[i] = std::max(0, intervals[i][1] - intervals[i][0] + 1 +
                              grid.region.boundary_size * 2);
  }
}

template <typename T, size_t D>
void get_view_intervals_old(iter_type (&intervals)[D][2],
                            iter_type (&stride)[D], iter_type (&dims)[D],
                            RegionalGrid<T, D>& grid, T cutoff_value,
                            len_type boundary_size) {
  for (iter_type ii = 0; ii < D; ++ii) {
    iter_type pos[D]{};

    for (iter_type n = 0;
         n < grid::length<D>(dims) / dims[ii] &&
         (intervals[ii][0] > 0 || intervals[ii][1] < dims[ii] - 1);
         ++n) {
      for (iter_type m = 0; m < intervals[ii][0]; ++m) {
        iter_type index = index_from_position(pos, stride, boundary_size);
        if (compare_cutoff(grid.values, index, cutoff_value)) {
          intervals[ii][0] = std::min(pos[ii], intervals[ii][0]);
        }
        ++pos[ii];
      }
      pos[ii] = dims[ii];
      for (iter_type m = dims[ii] - 1; m >= intervals[ii][1]; --m) {
        --pos[ii];
        iter_type index = index_from_position(pos, stride, boundary_size);
        if (compare_cutoff(grid.values, index, cutoff_value)) {
          intervals[ii][1] = std::max(pos[ii], intervals[ii][1]);
        }
      }

      pos[ii] = 0;
      for (iter_type ij = 0; ij < D; ++ij) {
        if (ij != ii) {
          ++pos[ij];
          if (pos[ij] == dims[ij]) {
            ++pos[ij + 1 + ((ij + 1 == ii) ? 1 : 0)];
            pos[ij] = 0;
          }
        }
      }
    }
  }
}

template <typename T, size_t D>
void get_view_intervals(iter_type (&intervals)[D][2], iter_type (&stride)[D],
                        iter_type (&dims)[D], RegionalGrid<T, D>& grid,
                        T cutoff_value, len_type boundary_size) {
  iter_type pos[D]{};
  for (iter_type n = 0; n < grid::length<D>(dims); ++n) {
    for (iter_type d = 0; d < D; ++d) {
      pos[d] = (n / stride[d]) % dims[d];
    }
    if (compare_cutoff(grid.values, n, cutoff_value)) {
      for (iter_type d = 0; d < D; ++d) {
        intervals[d][0] = std::min(pos[d], intervals[d][0]);
        intervals[d][1] = std::max(pos[d], intervals[d][1]);
      }
    }
  }
}

//! Construct a new view around a region with values greater than a cutoff.
/*!
 * Given the regional grid, populate the origin for the new region
 * where all values outside are less than cutoff. The dimensions of the
 * regional grid will be the same, since only the position of the view
 * itself is adjusted.
 *
 * \param grid The regional grid from which to find a new view.
 * \param origin The values of the new region origin are populated here.
 * \param boundary_size The width of the boundary used in the regional_grid.
 */
template <typename T, size_t D>
void get_view(RegionalGrid<T, D>& grid, T cutoff_value,
              iter_type (&origin)[D]) {
  iter_type intervals[D][2];
  iter_type stride[D]{};
  iter_type dims[D]{};
  for (iter_type i = 0; i < D; ++i) {
    dims[i] = grid.region.dims[i] - grid.region.boundary_size * 2;
    intervals[i][0] = dims[i];
    intervals[i][1] = 0;
  }
  init_region_stride(stride, grid.region.dims);
  get_view_intervals(intervals, stride, dims, grid, cutoff_value,
                     grid.region.boundary_size);

  for (iter_type i = 0; i < D; ++i) {
    dims[i] =
        intervals[i][1] - intervals[i][0] + 1 + grid.region.boundary_size * 2;
    iter_type delta = grid.region.dims[i] - dims[i];
    origin[i] = grid.region.origin[i] + (intervals[i][0] - (delta / 2));
  }
}

template <typename T, size_t D>
void get_view_resized_periodic(RegionalGrid<T, D>& grid, T cutoff_value,
                               iter_type (&origin)[D], len_type (&dims)[D]) {
  iter_type offset[D]{};
  iter_type stride[D]{};

  for (iter_type i = 0; i < D; ++i) {
    dims[i] = grid.region.dims[i] - grid.region.boundary_size * 2;
    offset[i] = grid.region.boundary_size;
  }
  init_region_stride(stride, grid.region.dims);

  iter_type start_n = -1;

  for (iter_type n = 0; n < grid::length<D>(dims) && start_n < 0; ++n) {
    iter_type pos[D];
    get_grid_position_offset(pos, dims, offset, n);
    iter_type index = index_from_position(pos, stride);

    using std::abs;
    if (compare_cutoff(grid.values, index, cutoff_value)) {
      start_n = n;
    }
  }

  if (start_n < 0) {
    for (iter_type i = 0; i < D; ++i) {
      dims[i] = 0;
    }
  } else {
    iter_type intervals[D][2];
    iter_type intervals_updated[D][2];
    iter_type start_pos[D];
    get_grid_position(start_pos, dims, start_n);
    for (iter_type i = 0; i < D; ++i) {
      intervals_updated[i][0] = start_pos[i];
      intervals_updated[i][1] = start_pos[i];
    }

    do {
      for (iter_type i = 0; i < D; ++i) {
        intervals[i][0] = intervals_updated[i][0];
        intervals[i][1] = intervals_updated[i][1];
      }

      for (iter_type i = 0; i < D; ++i) {
        bool offset = true;

        // set the starting position
        iter_type pos[D]{};
        for (iter_type j = 0; j < D; ++j) {
          pos[j] = (j == i) ? (intervals_updated[j][1] + 1)
                            : intervals_updated[j][0];
          pos[j] = (pos[j] >= dims[j]) ? pos[j] - dims[j]
                   : (pos[j] < 0)      ? dims[i] + pos[j]
                                       : pos[j];
        }

        // search within the intervals of the existing region
        len_type search_dims[D - 1]{};
        len_type search_stride[D - 1]{};
        for (iter_type j = 0; j < i; ++j) {
          search_dims[j] =
              intervals_updated[j][1] - intervals_updated[j][0] + 1;
        }
        for (iter_type j = i + 1; j < D; ++j) {
          search_dims[j - 1] =
              intervals_updated[j][1] - intervals_updated[j][0] + 1;
        }

        init_region_stride(search_stride, search_dims);

        len_type len =
            dims[i] - (intervals_updated[i][1] - intervals_updated[i][0] + 1);
        bool found_value[]{true, false};
        iter_type last_pos[]{intervals_updated[i][0], intervals_updated[i][1]};
        for (iter_type n = 0; n < len; ++n) {
          for (iter_type m = 0; m < grid::length<D - 1>(search_dims) &&
                                !found_value[(offset) ? 1 : 0];
               ++m) {
            len_type search_pos[D - 1]{};
            get_grid_position(search_pos, search_dims, m);
            for (iter_type j = 0; j < i; ++j) {
              pos[j] = search_pos[j] + intervals_updated[j][0];
              pos[j] = (pos[j] < 0)          ? pos[j] + dims[j]
                       : (pos[j] >= dims[j]) ? pos[j] - dims[j]
                                             : pos[j];
            }
            for (iter_type j = i + 1; j < D; ++j) {
              pos[j] = search_pos[j - 1] + intervals_updated[j][0];
              pos[j] = (pos[j] < 0)          ? pos[j] + dims[j]
                       : (pos[j] >= dims[j]) ? pos[j] - dims[j]
                                             : pos[j];
            }

            iter_type index = grid::index_from_position(
                pos, stride, grid.region.boundary_size);
            if (compare_cutoff(grid.values, index, cutoff_value)) {
              found_value[(offset) ? 1 : 0] = true;
              // if (offset)
              //{
              //	intervals_updated[i][1]++;
              //	found_value[1] = true;
              // }
              // else
              //{
              //	intervals_updated[i][0]--;
              //	found_value[0] = true;
              // }
            }
          }

          if (offset) {
            if (found_value[1]) {
              intervals_updated[i][1] = (pos[i] < intervals_updated[i][1])
                                            ? pos[i] + dims[i]
                                            : pos[i];
              ++last_pos[1];
            }

            pos[i] = (last_pos[0] <= 0) ? dims[i] + last_pos[0] - 1
                                        : last_pos[0] - 1;
            found_value[0] = false;
            offset = !offset;
          } else {
            if (found_value[0]) {
              intervals_updated[i][0] = (pos[i] > intervals_updated[i][0])
                                            ? pos[i] - dims[i]
                                            : pos[i];
              --last_pos[0];
            }

            pos[i] = (last_pos[1] >= dims[i] - 1) ? last_pos[1] + 1 - dims[i]
                                                  : last_pos[1] + 1;
            found_value[1] = false;
            offset = !offset;
          }
        }
      }

    } while (!grid::is_fully_overlapping(intervals, intervals_updated));

    for (iter_type i = 0; i < D; ++i) {
      intervals[i][0] = intervals_updated[i][0];
      intervals[i][1] = intervals_updated[i][1];
    }

    for (iter_type i = 0; i < D; ++i) {
      origin[i] = intervals[i][0] + grid.region.origin[i];
      dims[i] = std::max(0, intervals[i][1] - intervals[i][0] + 1 +
                                grid.region.boundary_size * 2);
    }
  }
}

template <typename T, size_t D>
void adjust_region(RegionalGrid<T, D>& grid, T cutoff) {
  iter_type origin[D];
  get_view(grid, cutoff, origin);

  for (iter_type i = 0; i < D; ++i) {
    len_type delta = (grid.dims[i] - grid.region.boundary_size * 2);
    origin[i] += (origin[i] < 0) ? delta : (origin[i] >= delta) ? -delta : 0;
  }

  grid.adjust(origin);
}

//! Adjust the region by finding a new origin and dimensions.
/*!
 * Adjust the region by finding a new origin and dimensions.
 *
 * \param grid The grid to resize.
 * \param cutoff The cutoff value to use when determining the region.
 * \param padding_factor Increases the smallest determined dimensions by the
 * given ratio. \param dims_relative_eps Computes the relative distance
 * between new and current dimensions, and uses the new dimensions if the
 * relative distance exceeds this value. The relative distance is defined
 * as:
 * $\f(\\text{dims}_0 - \\text{dims}_1) / \\text{dims}_0\f$, where subscript
 * 0 represents current dimensions and subscript 1 represents the newly
 * determined dimensions.
 */
template <typename T, size_t D>
void resize_adjust_region(RegionalGrid<T, D>& grid, T cutoff,
                          double padding_factor = 1.0,
                          double dims_relative_eps = 0.0) {
  iter_type origin[D]{};
  len_type dims[D]{};

  len_type dims_set[D]{};
  len_type origin_set[D]{};

  if (grid::is_same_point(grid.dims, grid.region.dims)) {
    get_view_resized_periodic(grid, cutoff, origin, dims);
  } else {
    get_view_resized(grid, cutoff, origin, dims);
  }

  bool same_dims_flag = true;
  for (iter_type i = 0; i < D; ++i) {
    auto dim0 = std::min(grid.dims[i], iter_type(dims[i] * padding_factor));
    origin_set[i] = origin[i] - (dim0 - dims[i]) / 2;
    origin_set[i] += (origin_set[i] < 0)
                         ? (grid.dims[i] - grid.region.boundary_size * 2)
                         : 0;
    dims_set[i] = dim0;

    // if the new dimensions are close to the current dimensions, don't
    // change it
    if (std::abs(dims_set[i] - grid.region.dims[i]) /
            double(grid.region.dims[i]) >
        dims_relative_eps) {
      same_dims_flag = false;
    }
  }

  if (same_dims_flag) {
    for (iter_type i = 0; i < D; ++i) {
      origin_set[i] = origin[i] - (grid.region.dims[i] - dims[i]) / 2;
      origin_set[i] += (origin_set[i] < 0)
                           ? (grid.dims[i] - grid.region.boundary_size * 2)
                           : 0;
    }
    grid.adjust(grid::select_region<D>(origin_set, grid.region.dims));
  } else {
    grid.adjust(grid::select_region<D>(origin_set, dims_set));
  }
}

template <typename T, size_t D>
void resize_adjust_region(RegionalGrid<T, D>& grid,
                          symphas::grid_info const& info) {
  len_type origin[D]{};
  len_type dims[D]{};
  for (auto const& [axis, interval] : info) {
    dims[symphas::axis_to_index(axis)] = interval.get_interval_count();
    origin[symphas::axis_to_index(axis)] =
        len_type((interval.left() - interval.domain_left()) / interval.width());
  }
  grid.adjust(grid::select_region<D>(origin, dims));
}

template <typename T, size_t D>
void resize_adjust_region(RegionalGrid<T, D>& grid,
                          const len_type (&intervals)[D][2]) {
  len_type origin[D]{};
  len_type dims[D]{};
  for (iter_type i = 0; i < D; ++i) {
    dims[i] = intervals[i][1] - intervals[i][0];
    origin[i] = intervals[i][0];
  }
  grid.adjust(grid::select_region<D>(origin, dims));
}

template <typename T, size_t D>
void resize_adjust_region(RegionalGrid<T, D>& grid,
                          grid::region_interval<D> const& interval) {
  resize_adjust_region(grid, interval.intervals);
}

template <typename T, size_t D>
void resize_adjust_region(RegionalGrid<T, D>& grid,
                          grid::region_interval_multiple<D> const& regions) {
  len_type intervals[D][2]{};

  for (iter_type i = 0; i < D; ++i) {
    intervals[i][0] = grid.dims[i];
    intervals[i][1] = 0;
  }

  for (grid::region_interval<D> region : regions) {
    for (iter_type i = 0; i < D; ++i) {
      intervals[i][0] = std::min(intervals[i][0], region[i][0]);
      intervals[i][1] = std::max(intervals[i][1], region[i][1]);
    }
  }
  resize_adjust_region(grid, intervals);
}

template <typename T, size_t D>
void resize_adjust_region(RegionalGrid<T, D>& grid, T cutoff,
                          const len_type (&minimum_dims)[D]) {
  iter_type origin[D];
  len_type dims[D];
  if (grid::is_same_point(grid.dims, grid.region.dims)) {
    get_view_resized_periodic(grid, cutoff, origin, dims);
  } else {
    get_view_resized(grid, cutoff, origin, dims);
  }

  for (iter_type i = 0; i < D; ++i) {
    iter_type delta = minimum_dims[i] - dims[i];
    if (delta > 0) {
      origin[i] -= delta / 2;
      dims[i] += delta - (delta / 2);
    }
  }
  grid.adjust(grid::select_region<D>(origin, dims));
}

}  // namespace grid

//! @}
