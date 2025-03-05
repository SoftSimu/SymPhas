
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
 * MODULE:  expr
 * PURPOSE: Defines an iterator that is used for iterating over data so that
 * expressions can be evaluated into the data.
 *
 * ***************************************************************************
 */

#pragma once

#include "grid.h"
#include "gridcudaincludes.h"

namespace grid {
template <size_t D>
struct region_interval_multiple;
template <size_t D>
struct region_interval;
template <size_t D>
struct region_extent;
struct region_size;

inline len_type length(grid::region_interval<0> const& interval) { return 1; }

template <size_t D>
len_type length(grid::region_interval<D> const& interval) {
  len_type dims[D];
  for (iter_type i = 0; i < D; ++i) {
    dims[i] = interval[i][1] - interval[i][0];
  }
  return length<D>(dims);
}

template <size_t D>
len_type length(grid::region_interval_multiple<D> const& interval) {
  len_type len = 0;
  for (grid::region_interval<D> region : interval) {
    len += length<D>(region);
  }
  return len;
}

inline len_type length(grid::region_size const& interval);

template <size_t D>
len_type length(grid::region_extent<D> const& region) {
  return length<D>(region.dims);
}

template <size_t D>
len_type length(const iter_type (&intervals)[D][2]) {
  return length<D>(region_extent(intervals));
}

template <size_t D>
bool is_empty(grid::region_interval<D> const& interval) {
  return length<D>(interval) == 0;
}

template <size_t D>
bool is_empty(grid::region_interval_multiple<D> const& interval) {
  return interval.regions.size() == 0;
}
}  // namespace grid

namespace symphas::lib {
template <size_t D, typename region_iterator_data_type>
struct region_iterator_container_impl {
  using arr_entry_t = iter_type[2];
  using interval_entry_t = iter_type[D][2];

  template <size_t... Is>
  region_iterator_container_impl(const iter_type (&dims)[D],
                                 region_iterator_data_type data, iter_type pos,
                                 std::index_sequence<Is...>)
      : pos{pos}, dims{dims[Is]...}, data{data} {}
  region_iterator_container_impl(const iter_type (&dims)[D],
                                 region_iterator_data_type data, iter_type pos)
      : region_iterator_container_impl(dims, data, pos,
                                       std::make_index_sequence<D>{}) {}
  template <typename T>
  region_iterator_container_impl(T&& regions, iter_type pos = 0)
      : region_iterator_container_impl(std::forward<T>(regions).dims,
                                       std::forward<T>(regions).regions.begin(),
                                       pos) {}

  operator grid::region_interval<D>() const {
    return grid::region_interval(
        dims, static_cast<interval_entry_t const&>(data[pos]));
  }

  arr_entry_t& operator[](iter_type i) { return (data[pos])[i]; }

  const arr_entry_t& operator[](iter_type i) const { return (data[pos])[i]; }

  operator const interval_entry_t&() const { return (data[pos]); }

  operator interval_entry_t&() { return (data[pos]); }

  iter_type pos;
  len_type dims[D];

 protected:
  region_iterator_data_type data;
};

template <size_t D>
struct basic_forward_iterator_container<grid::region_interval_multiple<D>>
    : region_iterator_container_impl<
          D, typename std::vector<typename grid::region_interval_multiple<
                 D>::subdomain_t>::iterator> {
  using parent_type = region_iterator_container_impl<
      D, typename std::vector<typename grid::region_interval_multiple<
             D>::subdomain_t>::iterator>;
  using parent_type::parent_type;
};

template <size_t D>
struct basic_forward_iterator_container<const grid::region_interval_multiple<D>>
    : region_iterator_container_impl<
          D, typename std::vector<typename grid::region_interval_multiple<
                 D>::subdomain_t>::const_iterator> {
  using parent_type = region_iterator_container_impl<
      D, typename std::vector<typename grid::region_interval_multiple<
             D>::subdomain_t>::const_iterator>;
  using parent_type::parent_type;
};

}  // namespace symphas::lib

namespace grid {
template <size_t D, size_t I0, size_t... Is>
iter_type stride_element(const len_type (&dims)[D],
                         std::index_sequence<I0, Is...>) {
  return (dims[I0] * ... * dims[Is]);
}

template <size_t D>
iter_type stride_element(const len_type (&dims)[D], std::index_sequence<>) {
  return 1;
}

template <size_t I, size_t D>
iter_type stride_element(const len_type (&dims)[D]) {
  return stride_element(dims, std::make_index_sequence<I>{});
}

template <size_t D>
struct region_interval {
  iter_type intervals[D][2];  //!< Intervals of the region within the global
                              //!< dimensions.
  iter_type dims[D];          //!< Global dimensions.

  template <size_t... Is>
  region_interval(const len_type (&dims)[D], len_type boundary_size,
                  std::index_sequence<Is...>)
      : intervals{{boundary_size, dims[Is] - boundary_size}...},
        dims{dims[Is]...} {}
  region_interval(const len_type (&dims)[D], len_type boundary_size = 0)
      : region_interval(dims, boundary_size, std::make_index_sequence<D>{}) {}

  template <size_t... Is>
  region_interval(const len_type (&origin)[D], const len_type (&dims)[D],
                  std::index_sequence<Is...>)
      : intervals{{origin[Is], origin[Is] + dims[Is]}...}, dims{dims[Is]...} {}
  region_interval(const len_type (&origin)[D], const len_type (&dims)[D])
      : region_interval(origin, dims, std::make_index_sequence<D>{}) {}

  template <size_t... Is>
  region_interval(const len_type (&dims)[D], const len_type (&intervals)[D][2],
                  std::index_sequence<Is...>)
      : intervals{{intervals[Is][0], intervals[Is][1]}...}, dims{dims[Is]...} {}
  region_interval(const len_type (&dims)[D], const len_type (&intervals)[D][2])
      : region_interval(dims, intervals, std::make_index_sequence<D>{}) {}

  template <size_t... Is>
  region_interval(const len_type (&dims)[D],
                  symphas::interval_data_type const& intervals,
                  std::index_sequence<Is...>)
      : region_interval(
            symphas::arr_l_t<D>{dims[Is]...},
            symphas::arr2_l_t<D>{symphas::arr_l_t<2>{
                len_type(intervals.at(symphas::index_to_axis(Is)).left() /
                         intervals.at(symphas::index_to_axis(Is)).width()),
                len_type(intervals.at(symphas::index_to_axis(Is)).right() /
                         intervals.at(symphas::index_to_axis(Is)).width()) +
                    1}...}) {}

  region_interval(const len_type (&dims)[D],
                  symphas::interval_data_type const& intervals)
      : region_interval(dims, intervals, std::make_index_sequence<D>{}) {}

  template <size_t... Is>
  region_interval(grid::dim_list const& dims,
                  symphas::interval_data_type const& intervals,
                  std::index_sequence<Is...>)
      : region_interval(symphas::arr_l_t<D>{dims[Is]...}, intervals) {}
  region_interval(grid::dim_list const& dims,
                  symphas::interval_data_type const& intervals)
      : region_interval(dims, intervals, std::make_index_sequence<D>{}) {}
  region_interval(symphas::grid_info const& info)
      : region_interval(info.get_dims(), info.intervals) {}

  region_interval() : intervals{}, dims{} {}

  using arr_entry_t = iter_type[2];
  arr_entry_t& operator[](iter_type i) { return intervals[i]; }
  const arr_entry_t& operator[](iter_type i) const { return intervals[i]; }

  region_interval& operator+=(region_interval<D> const& other) {
    for (iter_type i = 0; i < D; ++i) {
      intervals[i][0] = std::min(intervals[i][0], other[i][0]);
      intervals[i][1] = std::max(intervals[i][1], other[i][1]);
    }
    return *this;
  }

  region_interval& operator+=(region_interval<0> const& other) { return *this; }

  region_interval& operator/=(region_interval<D> const& other) {
    for (iter_type i = 0; i < D; ++i) {
      intervals[i][0] = std::max(intervals[i][0], other[i][0]);
      intervals[i][1] = std::min(intervals[i][1], other[i][1]);

      if (intervals[i][1] <= intervals[i][0]) {
        intervals[i][0] = 0;
        intervals[i][1] = 0;
      }
    }
    return *this;
  }

  // returns the local index given the index in the global region
  template <size_t... Is>
  iter_type to_local_index(iter_type n, std::index_sequence<Is...>) {
    return grid::index_from_position(
        n, dims, symphas::arr_n_t<D>{-intervals[Is][0]...},
        symphas::arr_n_t<D>{stride_element<Is>(
            symphas::arr_n_t<D>{(intervals[Is][1] - intervals[Is][0])...})...});
  }

  iter_type to_local_index(iter_type n) {
    return to_local_index(n, std::make_index_sequence<D>{});
  }

  // returns the global index given the index in the local region
  template <size_t... Is>
  iter_type to_global_index(iter_type n, std::index_sequence<Is...>) {
    return grid::index_from_position(
        n, symphas::arr_n_t<D>{(intervals[Is][1] - intervals[Is][0])...},
        symphas::arr_n_t<D>{intervals[Is][0]...},
        symphas::arr_n_t<D>{stride_element<Is>(dims)...});
  }

  iter_type to_global_index(iter_type n) {
    return to_global_index(n, std::make_index_sequence<D>{});
  }

  operator symphas::grid_info() const { return intervals; }
};

template <>
struct region_interval<0> {
  region_interval() {}

  region_interval& operator+=(region_interval<0> const& other) { return *this; }
};

template <size_t D>
region_interval(len_type (&)[D]) -> region_interval<D>;

template <typename interval_type, size_t D>
auto break_region(const iter_type (&intervals)[D][2],
                  region_interval<D> const& region);

template <size_t D>
struct region_interval_multiple : region_interval<D> {
  using parent_type = region_interval<D>;
  using parent_type::dims;
  using parent_type::intervals;
  using parent_type::operator[];

  struct subdomain_t {
    iter_type intervals[D][2];  //!< Intervals of the region within the global
                                //!< dimensions.
    using arr_entry_t = iter_type[2];
    using interval_entry_t = iter_type[D][2];

    template <size_t... Is>
    subdomain_t(const iter_type (&intervals)[D][2], std::index_sequence<Is...>)
        : intervals{{intervals[Is][0], intervals[Is][1]}...} {}
    subdomain_t(const iter_type (&intervals)[D][2])
        : subdomain_t(intervals, std::make_index_sequence<D>{}) {}
    subdomain_t() : intervals{} {}

    arr_entry_t& operator[](iter_type i) { return intervals[i]; }
    const arr_entry_t& operator[](iter_type i) const { return intervals[i]; }

    operator const interval_entry_t&() const { return intervals; }

    operator interval_entry_t&() { return intervals; }
  };

  using interval_type = subdomain_t;

  region_interval_multiple(const iter_type (&intervals)[D][2])
      : parent_type(intervals) {}
  region_interval_multiple(region_interval<D> const& region)
      : parent_type(region.dims),
        regions{break_region<interval_type>(intervals, region)} {}

  region_interval_multiple(const len_type (&dims)[D], len_type boundary_size)
      : parent_type(dims), regions{} {
    for (iter_type i = 0; i < D; ++i) {
      intervals[i][0] = boundary_size;
      intervals[i][1] = dims[i] - boundary_size;
    }
  }

  region_interval_multiple() : parent_type(), regions{} {}

  region_interval_multiple& operator+=(
      std::vector<interval_type> const& add_regions) {
    TIME_THIS_CONTEXT_LIFETIME(region_interval_multiple_union)
    for (const auto& add_region : add_regions) {
      if (grid::length<D>(add_region) > 0) {
        auto last_union = regions.end();
        for (auto it = regions.begin(); it < regions.end();) {
          if (last_union == regions.end()) {
            if (grid::is_contact_overlapping(add_region, *it)) {
              grid::get_region_union(*it, *it, add_region);
              last_union = it;
            }
            ++it;
          } else {
            if (grid::is_contact_overlapping(*last_union, *it)) {
              grid::get_region_union(*last_union, *last_union, *it);
              auto diff = it - regions.begin();
              regions.erase(it, it + 1);
              it = regions.begin() + diff;
            } else {
              ++it;
            }
          }
        }
        if (last_union == regions.end()) {
          regions.push_back(add_region);
        }
      }
    }

    return *this;
  }

  region_interval_multiple& operator-=(
      std::vector<interval_type> const& diff_regions) {
    TIME_THIS_CONTEXT_LIFETIME(region_interval_multiple_difference)
    std::vector<interval_type> result_regions;
    for (const auto& region : regions) {
      std::vector<interval_type> work_list;
      work_list.push_back(region);

      for (iter_type i = 0; i < D; ++i) {
        for (auto it = work_list.begin(); it < work_list.end();) {
          bool split = false;
          for (const auto& diff_region : diff_regions) {
            if (grid::is_fully_overlapping(diff_region, *it)) {
              auto diff = it - work_list.begin();
              work_list.erase(it, it + 1);
              it = work_list.begin() + diff;

              split = true;
            }
            if (grid::is_overlapping(diff_region, *it)) {
              auto intersection(*it);
              grid::get_region_intersection(intersection, *it, diff_region);

              auto first(*it);
              first[i][1] = intersection[i][0];

              auto second(*it);
              second[i][0] = intersection[i][0];
              second[i][1] = intersection[i][1];

              auto third(*it);
              third[i][0] = intersection[i][1];

              auto diff = it - work_list.begin();
              work_list.erase(work_list.begin() + diff,
                              work_list.begin() + diff + 1);
              for (auto check_region : {first, second, third}) {
                if (grid::length<D>(check_region) > 0) {
                  if (i < D - 1) {
                    work_list.insert(work_list.begin() + diff++, check_region);
                  } else if (!grid::is_fully_overlapping(diff_region,
                                                         check_region)) {
                    result_regions.insert(result_regions.begin(), check_region);
                  }
                }
              }
              it = work_list.begin() + diff;

              split = true;
            }

            if (split) {
              break;
            }
          }

          if (!split) {
            result_regions.insert(result_regions.end(), *it);
            auto diff = it - work_list.begin();
            work_list.erase(it, it + 1);
            it = work_list.begin() + diff;
          }
        }
      }
    }
    std::swap(regions, result_regions);
    return *this;
  }

  region_interval_multiple& operator+=(
      region_interval_multiple<D> const& add_regions) {
    return operator+=(add_regions.regions);
  }

  region_interval_multiple& operator+=(region_interval<D> const& region) {
    return operator+=(break_region<interval_type>(intervals, region));
  }

  region_interval_multiple& operator-=(
      region_interval_multiple<D> const& diff_regions) {
    return operator-=(diff_regions.regions);
  }

  region_interval_multiple& operator-=(region_interval<D> const& region) {
    return operator-=(break_region<interval_type>(intervals, region));
  }

  region_interval_multiple& operator/=(
      std::vector<interval_type> const& add_regions) {
    TIME_THIS_CONTEXT_LIFETIME(region_interval_multiple_intersection)
    std::vector<interval_type> work_list;
    work_list.reserve(regions.size());

    for (auto it = regions.begin(); it < regions.end(); ++it) {
      for (const auto& add_region : add_regions) {
        if (grid::is_overlapping(add_region, *it)) {
          work_list.push_back({});
          grid::get_region_intersection(work_list.back(), *it, add_region);
        }
      }
    }
    std::swap(work_list, regions);

    return *this;
  }

  region_interval_multiple& operator/=(
      region_interval_multiple<D> const& add_regions) {
    return operator/=(add_regions.regions);
  }

  region_interval_multiple& operator/=(region_interval<D> const& region) {
    return operator/=(break_region<interval_type>(intervals, region));
  }

  auto begin() const;
  auto begin();
  auto end() const;
  auto end();

  //! Convert the multiple region into a single unified region.
  /*!
   * The multiple region will be converted to the smallest possible
   * region_interval, assuming periodic boundary conditions. Therefore, the
   * interval of the result might be outside the interval of the multiple
   * region.
   *
   * The point is to have a minimal region_interval which can contain this
   * multiple region assuming periodic boundary conditions (this is the default
   * assumption made when breaking up an input region when performing the union.
   */
  grid::region_interval<D> operator+() const {
    if (regions.size() == 0) {
      return *this;
    }

    grid::region_interval<D> combined;
    for (iter_type i = 0; i < D; ++i) {
      combined.dims[i] = dims[i];
      combined[i][0] = dims[i];
    }
    for (grid::region_interval<D> region : *this) {
      for (iter_type i = 0; i < D; ++i) {
        if (region.intervals[i][0] < 0) {
          region.intervals[i][0] += combined.dims[i];
          region.intervals[i][1] += combined.dims[i];
        }

        combined[i][0] = std::min(combined[i][0], region.intervals[i][0]);
        combined[i][1] = std::max(combined[i][1], region.intervals[i][1]);
      }
    }
    return combined;

    // grid::region_interval<D> combined;
    // for (iter_type i = 0; i < D; ++i) {
    //   combined.dims[i] = dims[i];
    //   combined[i][0] = dims[i];
    // }

    // for (grid::region_interval<D> region : *this) {
    //   for (iter_type i = 0; i < D; ++i) {
    //     len_type width = intervals[i][1] - intervals[i][0];
    //     iter_type left = region[i][0] - intervals[i][0];
    //     iter_type right = region[i][1] - intervals[i][0];
    //     left -= width * std::round(double(left) / width);

    //    if (left < 0) {
    //      left += intervals[i][0];
    //      right -= width - intervals[i][0];
    //    } else {
    //      left += intervals[i][0];
    //      right += intervals[i][0];
    //    }

    //    combined[i][0] = std::min(combined[i][0], left);
    //    combined[i][1] = std::max(combined[i][1], right);
    //  }
    //}
    // return combined;
  }

  std::vector<interval_type> regions;
};

struct region_size {
  region_size() : len{0} {}
  region_size(len_type len) : len{len} {}
  len_type len;

  operator len_type() const { return len; }

  region_size& operator+=(region_size const& other) {
    len = std::max(len, other.len);
    return *this;
  }

  region_size& operator/=(region_size const& other) {
    len = std::min(len, other.len);
    return *this;
  }

  region_size& operator-=(region_size const& other) {
    len = std::max(len, other.len) - std::min(len, other.len);
    return *this;
  }
};

len_type length(grid::region_size const& interval) { return interval; }

template <size_t D>
struct region_extent {
  iter_type dims[D];  //!< Dimensions of the local region within global region.
  iter_type stride[D];  //!< Stride of the global region.
  iter_type shift[D];
  iter_type delta;

  template <size_t... Is>
  region_extent(region_interval<D> const& interval, std::index_sequence<Is...>)
      : dims{(interval[Is][1] - interval[Is][0])...},
        stride{stride_element<Is>(interval.dims)...},
        shift{(interval[Is][0] * stride[Is])...},
        delta{compute_delta(dims, stride)} {}

  region_extent(region_interval<D> const& interval)
      : region_extent(interval, std::make_index_sequence<D>{}) {}

  template <size_t... Is>
  region_extent(const iter_type (&intervals)[D][2], std::index_sequence<Is...>)
      : dims{(intervals[Is][1] - intervals[Is][0])...},
        stride{stride_element<Is>(dims)...},
        shift{(intervals[Is][0] * stride[Is])...},
        delta{compute_delta(dims, stride)} {}

  region_extent(const iter_type (&intervals)[D][2])
      : region_extent(intervals, std::make_index_sequence<D>{}) {}

  template <size_t... Is>
  region_extent(const len_type (&origin)[D], const len_type (&dims)[D],
                len_type boundary_size, std::index_sequence<Is...>)
      : dims{dims[Is]...},
        stride{stride_element<Is>(dims)...},
        shift{((origin[Is] + boundary_size) * stride[Is])...},
        delta{compute_delta(dims, stride)} {}

  region_extent(const len_type (&origin)[D], const len_type (&dims)[D],
                len_type boundary_size = 0)
      : region_extent(origin, dims, boundary_size,
                      std::make_index_sequence<D>{}) {}

  template <size_t... Is>
  region_extent(const len_type (&dims)[D], len_type boundary_size,
                std::index_sequence<Is...>)
      : dims{dims[Is]...},
        stride{stride_element<Is>(dims)...},
        shift{(boundary_size * stride[Is])...},
        delta{compute_delta(dims, stride)} {}

  region_extent(const len_type (&dims)[D], len_type boundary_size = 0)
      : region_extent(dims, boundary_size, std::make_index_sequence<D>{}) {}

  region_extent(grid::select_region<D> const& region)
      : region_extent(region.origin, region.dims, region.boundary_size) {}

  region_extent(symphas::grid_info const& info)
      : region_extent(region_interval<D>(info)) {}

  iter_type operator[](iter_type n) const {
    return grid::index_from_position(n, dims, stride, shift, delta);
  }

  region_extent() : dims{}, stride{}, shift{}, delta{} {}
};

template <size_t D>
region_extent(const iter_type (&intervals)[D][2]) -> region_extent<D>;
template <size_t D>
region_extent(const iter_type (&dims)[D]) -> region_extent<D>;
template <size_t D>
region_extent(region_interval<D>) -> region_extent<D>;

template <size_t D>
using region_group = region_extent<D>;

// template<size_t D>
// struct region_group : region_extent<D>
//{
//	using parent_type = region_extent<D>;
//	using parent_type::dims;
//	using parent_type::stride;
//	using parent_type::shift;
//	using parent_type::delta;
//	using parent_type::operator[];

//	iter_type global_dims[D];		//!< Dimensions of the global
// domain. 	iter_type origin[D];			//!< Origin of the
// region, to check for bounds wrapping.

//	region_group(const len_type(&dims)[D]) :
//		parent_type(dims)
//	{
//		for (iter_type i = 0; i < D; ++i)
//		{
//			global_dims[i] = dims[i];
//		}
//	}

//	template<size_t... Is>
//	region_group(region_interval<D> const& interval,
// std::index_sequence<Is...>) : 		parent_type(interval,
// std::index_sequence<Is...>{}), 		global_dims{
// interval.dims[Is]... }, 		origin{ interval[Is][0]... } {}

//	region_group(region_interval<D> const& interval) :
// region_group(interval, std::make_index_sequence<D>{}) {}

//	region_group() : parent_type() {}

//};

template <size_t D>
struct region_index_list {
  iter_type* iters;
  len_type len;

  region_index_list(region_interval<D> const& intervals)
      : iters{(grid::length<D>(intervals) > 0)
                  ? new iter_type[grid::length<D>(intervals)]
                  : nullptr},
        len{grid::length<D>(intervals)} {
    region_extent region(intervals);
    iter_type origin[D];
    for (iter_type i = 0; i < D; ++i) {
      origin[i] = intervals[i][0];
    }
    for (iter_type n = 0; n < len; ++n) {
      iter_type pos[D];
      iters[n] =
          grid::index_from_position(n, region.dims, origin, region.stride);
    }
  }

  region_index_list(const len_type* dims, len_type boundary_size)
      : iters{(grid::length_interior<D>(dims, boundary_size) > 0)
                  ? new iter_type[grid::length_interior<D>(dims, boundary_size)]
                  : nullptr},
        len{grid::length_interior<D>(dims, boundary_size)} {
    iter_type origin[D];
    iter_type interior_dims[D];
    for (iter_type i = 0; i < D; ++i) {
      origin[i] = boundary_size;
      interior_dims[i] = dims[i] - boundary_size * 2;
    }

    iter_type stride[D];
    grid::get_stride(stride, dims);

    for (iter_type n = 0; n < len; ++n) {
      iter_type pos[D];
      iters[n] = grid::index_from_position(n, interior_dims, origin, stride);
    }
  }

  region_index_list(region_index_list<D> const& other)
      : iters{new iter_type[other.len]}, len{other.len} {
    std::copy(other.iters, other.iters + other.len, iters);
  }

  region_index_list(region_index_list<D>&& other) : region_index_list() {
    swap(*this, other);
  }

  region_index_list& operator=(region_index_list<D> other) {
    swap(*this, other);
    return *this;
  }

  friend void swap(region_index_list& first, region_index_list& second) {
    using std::swap;
    swap(first.iters, second.iters);
    swap(first.len, second.len);
  }

  iter_type operator[](iter_type n) const { return iters[n]; }

  region_index_list() : iters{nullptr}, len{0} {}

  ~region_index_list() { delete[] iters; }
};

template <size_t D>
using region_const_iterator = symphas::lib::basic_forward_iterator<
    const grid::region_interval_multiple<D>>;
template <size_t D>
using region_iterator =
    symphas::lib::basic_forward_iterator<grid::region_interval_multiple<D>>;

template <size_t D>
auto region_interval_multiple<D>::begin() const {
  return region_const_iterator<D>(*this);
}

template <size_t D>
auto region_interval_multiple<D>::begin() {
  return region_iterator<D>(*this);
}

template <size_t D>
auto region_interval_multiple<D>::end() const {
  iter_type n = (iter_type)regions.size();
  return region_const_iterator<D>(*this, n);
}

template <size_t D>
auto region_interval_multiple<D>::end() {
  return region_iterator<D>(*this, regions.size());
}

// an object representing no iterable region.
struct region_empty {};

template <size_t D>
region_interval_multiple<D> operator-(region_interval<D> const& first,
                                      region_interval<D> const& second) {
  region_interval_multiple<D> region(first);
  return region -= second;
}

template <size_t D>
region_interval_multiple<D> operator-(
    region_interval<D> const& first,
    region_interval_multiple<D> const& second) {
  region_interval_multiple<D> region(first);
  return region -= second;
}

template <size_t D>
region_interval_multiple<D> operator-(region_interval_multiple<D> const& first,
                                      region_interval<D> const& second) {
  region_interval_multiple<D> region(first);
  return region -= second;
}

template <size_t D>
region_interval_multiple<D> operator-(
    region_interval_multiple<D> const& first,
    region_interval_multiple<D> const& second) {
  region_interval_multiple<D> region(first);
  return region -= second;
}

template <size_t D>
region_interval<D> operator-(region_interval<D> const& first,
                             region_empty const& second) {
  return first;
}

template <size_t D>
region_interval_multiple<D> operator-(region_interval_multiple<D> const& first,
                                      region_empty const& second) {
  return first;
}

template <size_t D>
region_empty operator-(region_empty const& first,
                       region_interval<D> const& second) {
  return first;
}

template <size_t D>
region_empty operator-(region_empty const& first,
                       region_interval_multiple<D> const& second) {
  return first;
}

inline auto operator+(region_interval<0> const& first,
                      region_interval<0> const& second) {
  return first;
}

inline auto operator/(region_interval<0> const& first,
                      region_interval<0> const& second) {
  return first;
}

template <size_t D>
auto operator+(region_interval<D> const& first,
               region_interval<D> const& second) {
  auto result(first);
  result += second;
  return result;
}

template <size_t D>
auto operator/(region_interval<D> const& first,
               region_interval<D> const& second) {
  auto result(first);
  result /= second;
  return result;
}

inline auto operator+(region_size const& first, region_size const& second) {
  auto result(first);
  result += second;
  return result;
}

inline auto operator/(region_size const& first, region_size const& second) {
  auto result(first);
  result /= second;
  return result;
}

template <size_t D>
auto operator+(region_interval_multiple<D> const& first,
               region_interval_multiple<D> const& second) {
  auto result(first);
  result += second;
  return result;
}

template <size_t D>
auto operator/(region_interval_multiple<D> const& first,
               region_interval_multiple<D> const& second) {
  auto result(first);
  result /= second;
  return result;
}

inline auto operator+(region_size const& first, region_empty const& second) {
  return first;
}

inline auto operator+(region_empty const& first, region_size const& second) {
  return second + first;
}

inline auto operator/(region_size const& first, region_empty const& second) {
  auto result(first);
  result.len = 0;
  return result;
}

inline auto operator/(region_empty const& first, region_size const& second) {
  return second / first;
}

inline auto operator+(region_interval<0> const& first,
                      region_empty const& second) {
  return first;
}

inline auto operator+(region_empty const& first,
                      region_interval<0> const& second) {
  return second + first;
}

inline auto operator/(region_interval<0> const& first,
                      region_empty const& second) {
  return first;
}

inline auto operator/(region_empty const& first,
                      region_interval<0> const& second) {
  return second / first;
}

template <size_t D>
auto operator+(region_interval<D> const& first, region_size const& second) {
  return first;
}

template <size_t D>
auto operator+(region_size const& first, region_interval<D> const& second) {
  return second + first;
}

template <size_t D>
auto operator+(region_interval<D> const& first, region_empty const& second) {
  return first;
}

template <size_t D>
auto operator+(region_empty const& first, region_interval<D> const& second) {
  return second + first;
}

template <size_t D>
auto operator/(region_interval<D> const& first, region_size const& second) {
  return first;
}

template <size_t D>
auto operator/(region_size const& first, region_interval<D> const& second) {
  return second / first;
}

template <size_t D>
auto operator/(region_interval<D> const& first, region_empty const& second) {
  auto result(first);
  for (iter_type i = 0; i < D; ++i) {
    result[i][0] = 0;
    result[i][1] = 0;
  }
  return result;
}

template <size_t D>
auto operator/(region_empty const& first, region_interval<D> const& second) {
  return second / first;
}

template <size_t D>
auto operator+(region_interval<D> const& first,
               region_interval<0> const& second) {
  return first;
}

template <size_t D>
auto operator+(region_interval<0> const& first,
               region_interval<D> const& second) {
  return second + first;
}

template <size_t D>
auto operator/(region_interval<D> const& first,
               region_interval<0> const& second) {
  return first;
}

template <size_t D>
auto operator/(region_interval<0> const& first,
               region_interval<D> const& second) {
  return second / first;
}

inline auto operator+(region_interval<0> const& first,
                      region_size const& second) {
  return second;
}

inline auto operator+(region_size const& first,
                      region_interval<0> const& second) {
  return second + first;
}

inline auto operator/(region_interval<0> const& first,
                      region_size const& second) {
  return second;
}

inline auto operator/(region_size const& first,
                      region_interval<0> const& second) {
  return second / first;
}

template <size_t D>
auto operator+(region_interval_multiple<D> const& first,
               region_size const& second) {
  return first;
}

template <size_t D>
auto operator+(region_size const& first,
               region_interval_multiple<D> const& second) {
  return second + first;
}

template <size_t D>
auto operator+(region_interval_multiple<D> const& first,
               region_empty const& second) {
  return first;
}

template <size_t D>
auto operator+(region_empty const& first,
               region_interval_multiple<D> const& second) {
  return second + first;
}

template <size_t D>
auto operator/(region_interval_multiple<D> const& first,
               region_size const& second) {
  return first;
}

template <size_t D>
auto operator/(region_size const& first,
               region_interval_multiple<D> const& second) {
  return second / first;
}

template <size_t D>
auto operator/(region_interval_multiple<D> const& first,
               region_empty const& second) {
  auto result(first);
  result.regions.clear();
  return result;
}

template <size_t D>
auto operator/(region_empty const& first,
               region_interval_multiple<D> const& second) {
  return second / first;
}

template <size_t D>
auto operator+(region_interval_multiple<D> const& first,
               region_interval<D> const& second) {
  auto result(first);
  result += second;
  return result;
}

template <size_t D>
auto operator+(region_interval<D> const& first,
               region_interval_multiple<D> const& second) {
  return second + first;
}

template <size_t D>
auto operator/(region_interval_multiple<D> const& first,
               region_interval<D> const& second) {
  auto result(first);
  result /= second;
  return result;
}

template <size_t D>
auto operator/(region_interval<D> const& first,
               region_interval_multiple<D> const& second) {
  return second / first;
}

template <size_t D>
auto operator+(region_interval<0> const& first,
               region_interval_multiple<D> const& second) {
  return second;
}

template <size_t D>
auto operator+(region_interval_multiple<D> const& first,
               region_interval<0> const& second) {
  return second + first;
}

template <size_t D>
auto operator/(region_interval<0> const& first,
               region_interval_multiple<D> const& second) {
  return second;
}

template <size_t D>
auto operator/(region_interval_multiple<D> const& first,
               region_interval<0> const& second) {
  return second / first;
}

// template<size_t D>
// iter_type add_grid_position(iter_type offset, iter_type(&pos)[D], const
// iter_type (&dims)[D], const iter_type (&stride)[D])
//{
//	for (iter_type i = 0; i < D; ++i)
//	{
//		pos[i] += grid::position(offset, dims[i], stride[i]);
//	}
// }

// template<size_t D>
// struct region_position
//{
//	iter_type dims[D];			//!< Dimensions of the local
// region within global region. 	iter_type stride[D];		//!<
// Stride of the global region. 	iter_type pos[D];
// //!< Position inside the region represented by n. 	iter_type n;

//	template<size_t... Is>
//	region_position(region_interval<D> const& interval, std::ptrdiff_t pos,
// std::index_sequence<Is...>) : 		dims{ (interval[Is][1] -
// interval[Is][0])... }, 		stride{
// stride_element<Is>(interval.dims)... }, 		pos{
// position((iter_type)pos, (interval[Is][1] - interval[Is][0]),
// stride[Is])...},
// n{ grid::index_from_position((iter_type)pos, dims, symphas::arr_n_t<D> {
// interval[Is][0]... }, stride) }
//	{}

//	region_position(region_interval<D> const& interval, std::ptrdiff_t pos)
//: region_position(interval, pos, std::make_index_sequence<D>{}) {}

//	//! Dereference the iterator.
//	region_position& operator+=(region_position const& other)
//	{
//		n += other.n;
//		//region += other.region;
//		return *this;
//	};

//	//! Dereference the iterator.
//	region_position& operator-=(region_position const& other)
//	{
//		//region -= other.region;
//		n -= other.n;
//		return *this;
//	};

//	//! Dereference the iterator.
//	region_position& operator+=(iter_type offset)
//	{
//		add_grid_position(offset, pos, dims);
//		n = grid::index_from_position(pos, stride);
//		return *this;
//	};

//	//! Dereference the iterator.
//	region_position& operator-=(iter_type offset)
//	{
//		add_grid_position(-offset, pos, dims);
//		n = grid::index_from_position(pos, stride);
//		return *this;
//	};

//	//! Prefix decrement, returns itself.
//	region_position& operator--()
//	{
//		for (iter_type i = 0; i < D; ++i)
//		{
//			if (pos[i] == 0)
//			{
//				n -= stride[i];
//				pos[i] = dims[i] - 1;
//				if (i != 0) --pos[i + 1];
//			}
//		}
//		--n;
//		return *this;
//	}

//	//! Prefix decrement, returns itself.
//	region_position& operator++()
//	{
//		for (iter_type i = 0; i < D; ++i)
//		{
//			if (pos[i] == dims[i] - 1)
//			{
//				n += stride[i];
//				pos[i] = 0;
//				if (i != D - 1) ++pos[i + 1];
//			}
//		}
//		++n;
//		return *this;
//	}

//	iter_type operator[](iter_type offset) const
//	{
//		return n;
//	}

//	operator iter_type() const
//	{
//		return n;
//	}

//	region_position() : dims{}, stride{}, n{} {}

//};

template <size_t D, typename T>
bool is_fully_overlapping(select_region<D> const& region, T&& other) {
  return is_fully_overlapping<D>(
      region_interval<D>(region.origin, region.dims).intervals, other);
}

template <size_t D, typename T>
bool is_fully_overlapping(T const& other,
                          grid::region_interval_multiple<D> const& region) {
  return is_fully_overlapping<D>(other, (+region).intervals);
}

template <size_t D, typename T>
bool is_fully_overlapping(grid::region_interval_multiple<D> const& region,
                          T const& other) {
  return is_fully_overlapping<D>((+region).intervals, other);
}

template <size_t D>
bool is_fully_overlapping(select_region<D> const& region0,
                          grid::region_interval_multiple<D> const& region1) {
  return is_fully_overlapping<D>(
      region_interval<D>(region0.origin, region0.dims).intervals,
      (+region1).intervals);
}

template <size_t D>
bool is_fully_overlapping(grid::region_interval_multiple<D> const& region0,
                          grid::region_interval_multiple<D> const& region1) {
  return is_fully_overlapping<D>((+region0).intervals, (+region1).intervals);
}

template <size_t D>
bool is_fully_overlapping(grid::region_interval<D> const& region0,
                          grid::region_interval_multiple<D> const& region1) {
  return is_fully_overlapping<D>(region0.intervals, (+region1).intervals);
}

template <size_t D>
bool is_fully_overlapping(grid::region_interval_multiple<D> const& region0,
                          grid::region_interval<D> const& region1) {
  return is_fully_overlapping<D>((+region0).intervals, region1.intervals);
}

template <size_t D>
bool is_fully_overlapping(grid::region_interval<D> const& region0,
                          grid::region_interval<D> const& region1) {
  return is_fully_overlapping<D>(region0.intervals, region1.intervals);
}

template <size_t D>
bool is_overlapping(grid::region_interval<D> const& region0,
                    grid::region_interval<D> const& region1) {
  return is_overlapping<D>(region0.intervals, region1.intervals);
}
}  // namespace grid

namespace symphas::internal {
template <typename T>
struct data_value_type {
  using type = T;
  using ref = type&;

  ref operator()(T* data, iter_type n) { return data[n]; }
};

template <template <typename, size_t> typename G, typename T, size_t D>
struct data_value_type<G<T, D>> {
  using type = T;
  using ref = type&;

  ref operator()(G<T, D>* data, iter_type n) { return (*data)[n]; }
};

template <typename T, size_t D>
struct data_value_type<RegionalGrid<any_vector_t<T, D>, D>> {
  using type = any_vector_t<T, D>;
  using ref = multi_value<D, T>;

  ref operator()(RegionalGrid<any_vector_t<T, D>, D>* data, iter_type n) {
    return (*data)[n];
  }
};

template <typename T, size_t D>
struct data_value_type<RegionalGrid<T, D>> {
  using type = T;
  using ref = carry_value<T>;

  ref operator()(RegionalGrid<T, D>* data, iter_type n) { return (*data)[n]; }
};

#ifdef USING_MPI

template <typename T, size_t D>
struct data_value_type<RegionalGridMPI<any_vector_t<T, D>, D>> {
  using type = any_vector_t<T, D>;
  using ref = multi_value<D, T>;

  ref operator()(RegionalGridMPI<any_vector_t<T, D>, D>* data, iter_type n) {
    return (*data)[n];
  }
};

template <typename T, size_t D>
struct data_value_type<RegionalGridMPI<T, D>> {
  using type = T;
  using ref = carry_value<T>;

  ref operator()(RegionalGridMPI<T, D>* data, iter_type n) {
    return (*data)[n];
  }
};

#endif

template <typename T>
struct data_value_type<Block<T>> {
  using type = T;
  using ref = type&;

  ref operator()(Block<T>* data, iter_type n) { return (*data)[n]; }
};

template <size_t N, typename T>
struct data_value_type<MultiBlock<N, T>> {
  using type = any_vector_t<T, N>;
  using ref = multi_value<N, T>;

  ref operator()(MultiBlock<N, T>* data, iter_type n) { return (*data)[n]; }
};

template <typename T, size_t D>
struct data_value_type<any_vector_t<T, D>> {
  using type = any_vector_t<T, D>;
  using ref = type&;

  ref operator()(any_vector_t<T, D>* data, iter_type n) { return data[n]; }
};

template <template <typename, size_t> typename G, typename T, size_t D>
struct data_value_type<G<any_vector_t<T, D>, D>>
    : data_value_type<MultiBlock<D, T>> {
  using parent_type = data_value_type<MultiBlock<D, T>>;
  using typename parent_type::ref;
  using typename parent_type::type;
};

template <typename T>
struct data_value_type<const T> {
  using type = const typename data_value_type<T>::type;
  using ref = const typename data_value_type<T>::ref;

  ref operator()(const T* data, iter_type n) {
    return data_value_type<T>{}(const_cast<T*>(data), n);
  }
};

template <typename specialized_iterator>
struct iterator_difference_type_impl {
  using difference_type = std::ptrdiff_t;

  iterator_difference_type_impl(difference_type pos) : pos{pos} {}

  bool operator<(specialized_iterator const& other) const {
    return pos < other.pos;
  }

  bool operator<(size_t value) const { return pos < value; }

  //! Dereference the iterator.
  inline decltype(auto) operator+(specialized_iterator const& other) const {
    specialized_iterator diff(cast());
    diff.pos += other.pos;
    return diff;
  };

  //! Dereference the iterator.
  inline decltype(auto) operator-(specialized_iterator const& other) const {
    specialized_iterator diff(cast());
    diff.pos -= other.pos;
    return diff;
  };

  //! Dereference the iterator.
  inline decltype(auto) operator+(difference_type offset) const {
    specialized_iterator diff(cast());
    diff.pos += offset;
    return diff;
  };

  //! Dereference the iterator.
  inline decltype(auto) operator-(difference_type offset) const {
    specialized_iterator diff(cast());
    diff.pos += offset;
    return diff;
  };

  //! Dereference the iterator.
  inline decltype(auto) operator/(size_t div) const {
    specialized_iterator diff(cast());
    diff.pos /= div;
    return diff;
  };

  //! Dereference the iterator.
  inline decltype(auto) operator/(specialized_iterator const& div) const {
    specialized_iterator diff(cast());
    diff.pos /= div.pos;
    return diff;
  };

  //! Dereference the iterator.
  inline decltype(auto) operator%(size_t div) const {
    specialized_iterator diff(cast());
    diff.pos %= div;
    return diff;
  };

  //! Dereference the iterator.
  inline decltype(auto) operator%(specialized_iterator const& div) const {
    specialized_iterator diff(cast());
    diff.pos %= div.pos;
    return diff;
  };

  //! Dereference the iterator.
  inline decltype(auto) operator*(size_t div) const {
    specialized_iterator diff(cast());
    diff.pos *= div;
    return diff;
  };

  //! Dereference the iterator.
  inline decltype(auto) operator*(specialized_iterator const& div) const {
    specialized_iterator diff(cast());
    diff.pos *= div.pos;
    return diff;
  };

  //! Dereference the iterator.
  specialized_iterator& operator+=(specialized_iterator const& other) {
    pos += other.pos;
    return cast();
  };

  //! Dereference the iterator.
  specialized_iterator& operator-=(specialized_iterator const& other) {
    pos -= other.pos;
    return cast();
  };

  //! Dereference the iterator.
  specialized_iterator& operator+=(difference_type offset) {
    pos += offset;
    return cast();
  };

  //! Dereference the iterator.
  specialized_iterator& operator-=(difference_type offset) {
    pos -= offset;
    return cast();
  };

  //! Prefix decrement, returns itself.
  specialized_iterator& operator--() {
    --pos;
    return cast();
  }

  //! Prefix decrement, returns itself.
  specialized_iterator& operator++() {
    ++pos;
    return cast();
  }

  specialized_iterator operator-() const {
    specialized_iterator result(cast());
    result.pos = -result.pos;
    return result;
  }

  explicit operator size_t() const { return size_t(pos); }

  explicit operator int() const { return pos; }

  explicit operator difference_type() const { return difference_type(pos); }

  specialized_iterator& cast() {
    return *static_cast<specialized_iterator*>(this);
  }

  const specialized_iterator& cast() const {
    return *static_cast<specialized_iterator const*>(this);
  }

  difference_type pos;
};

template <typename G>
struct iterator_difference_type
    : iterator_difference_type_impl<iterator_difference_type<G>> {
  using parent_type =
      iterator_difference_type_impl<iterator_difference_type<G>>;
  using difference_type = typename parent_type::difference_type;
  using parent_type::pos;

  iterator_difference_type(difference_type pos)
      : iterator_difference_type(nullptr, pos) {}
  iterator_difference_type(G* ptr, difference_type pos)
      : parent_type(pos), ptr{ptr} {}
  iterator_difference_type() : iterator_difference_type(nullptr, 0) {}

  bool operator==(iterator_difference_type<G> const& other) const {
    return ptr == other.ptr && pos == other.pos;
  }

  bool operator==(size_t value) const { return pos == value; }

  //! Comparison with another iterator.
  /*!
   * Greater than comparison with another iterator.
   * Compares the current position.
   */
  bool operator>(iterator_difference_type<G> const& other) const {
    return pos > other.pos && ptr == other.ptr;
  }

  bool operator>(size_t value) const { return pos > value; }

  G* ptr;
};

template <typename G>
struct iterator_selection_difference_type
    : iterator_difference_type_impl<iterator_selection_difference_type<G>> {
  using parent_type =
      iterator_difference_type_impl<iterator_selection_difference_type<G>>;
  using difference_type = typename parent_type::difference_type;
  using parent_type::pos;

  iterator_selection_difference_type(difference_type pos)
      : iterator_selection_difference_type(nullptr, nullptr, pos) {}
  iterator_selection_difference_type(G* ptr, const iter_type* iters,
                                     difference_type pos)
      : parent_type(pos), ptr{ptr}, iters{iters} {}
  iterator_selection_difference_type()
      : iterator_selection_difference_type(nullptr, nullptr, 0) {}

  bool operator==(iterator_selection_difference_type<G> const& other) const {
    return ptr == other.ptr && pos == other.pos;
  }

  bool operator==(size_t value) const { return pos == value; }

  //! Comparison with another iterator.
  /*!
   * Greater than comparison with another iterator.
   * Compares the current position.
   */
  bool operator>(iterator_selection_difference_type<G> const& other) const {
    return pos > other.pos && ptr == other.ptr;
  }

  bool operator>(size_t value) const { return pos > value; }

  G* ptr;
  const iter_type* iters;
};

template <typename G, size_t D>
struct iterator_region_difference_type
    : iterator_difference_type_impl<iterator_region_difference_type<G, D>> {
  using parent_type =
      iterator_difference_type_impl<iterator_region_difference_type<G, D>>;
  using difference_type = typename parent_type::difference_type;
  using parent_type::pos;

  iterator_region_difference_type(difference_type pos)
      : iterator_region_difference_type(nullptr, grid::region_extent<D>(),
                                        pos) {}
  iterator_region_difference_type(G* ptr, grid::region_interval<D> intervals,
                                  difference_type pos)
      : parent_type(pos), ptr{ptr}, region{intervals} {}
  iterator_region_difference_type(G* ptr, grid::region_extent<D> region,
                                  difference_type pos)
      : parent_type(pos), ptr{ptr}, region{region} {}
  iterator_region_difference_type()
      : iterator_region_difference_type(nullptr, {}, 0) {}

  bool operator==(iterator_region_difference_type<G, D> const& other) const {
    return ptr == other.ptr && pos == other.pos;
  }

  bool operator==(size_t value) const { return pos == value; }

  //! Comparison with another iterator.
  /*!
   * Greater than comparison with another iterator.
   * Compares the current position.
   */
  bool operator>(iterator_region_difference_type<G, D> const& other) const {
    return pos > other.pos && ptr == other.ptr;
  }

  bool operator>(size_t value) const { return pos > value; }

  iter_type operator[](iter_type offset) const { return region[pos + offset]; }

  iter_type operator*() const {
    iter_type n = iter_type(pos);
    return region[n];
  }

  G* ptr;
  grid::region_extent<D> region;
};

template <typename G, size_t D>
struct iterator_group_difference_type
    : iterator_difference_type_impl<iterator_group_difference_type<G, D>> {
  using parent_type =
      iterator_difference_type_impl<iterator_group_difference_type<G, D>>;
  using difference_type = typename parent_type::difference_type;
  using parent_type::pos;

  iterator_group_difference_type(difference_type pos)
      : iterator_group_difference_type(nullptr, {}, pos) {}
  iterator_group_difference_type(G* ptr, grid::region_interval<D> intervals,
                                 difference_type pos)
      : parent_type(pos), ptr{ptr}, region{intervals} {}
  iterator_group_difference_type()
      : iterator_group_difference_type(nullptr, {}, 0) {}

  bool operator==(iterator_group_difference_type<G, D> const& other) const {
    return ptr == other.ptr && pos == other.pos;
  }

  bool operator==(size_t value) const { return pos == value; }

  //! Comparison with another iterator.
  /*!
   * Greater than comparison with another iterator.
   * Compares the current position.
   */
  bool operator>(iterator_group_difference_type<G, D> const& other) const {
    return pos > other.pos && ptr == other.ptr;
  }

  bool operator>(size_t value) const { return pos > value; }

  iter_type operator[](iter_type offset) const { return region[pos + offset]; }

  iter_type operator*() const { return region[pos]; }

  G* ptr;
  grid::region_group<D> region;
};

template <typename... iterator_ts>
struct reduce_iterator_difference_type
    : iterator_difference_type_impl<
          reduce_iterator_difference_type<iterator_ts...>> {
  using parent_type = iterator_difference_type_impl<
      reduce_iterator_difference_type<iterator_ts...>>;
  using difference_type = typename parent_type::difference_type;
  using parent_type::pos;

  reduce_iterator_difference_type(iterator_ts const&... iterators)
      : parent_type(get_position(iterators...)), iterators{iterators...} {}
  reduce_iterator_difference_type(difference_type pos = {})
      : parent_type(pos), iterators{} {}

  auto get_position(iterator_ts const&... iterators) {
    return std::min({iterators.ptr.pos...});
  }

  bool operator==(
      reduce_iterator_difference_type<iterator_ts...> const& other) const {
    return pos == other.pos;
  }

  bool operator==(size_t value) const { return pos == value; }

  //! Comparison with another iterator.
  /*!
   * Greater than comparison with another iterator.
   * Compares the current position.
   */
  bool operator>(
      reduce_iterator_difference_type<iterator_ts...> const& other) const {
    return pos > other.pos;
  }

  bool operator>(size_t value) const { return pos > value; }

  iter_type operator[](iter_type offset) const { return pos + offset; }

  iter_type operator*() const { return pos; }

  std::tuple<iterator_ts...> iterators;
};
}  // namespace symphas::internal

namespace symphas::internal {
template <typename E, size_t D>
struct iterator_group_expression;

template <typename G, size_t D>
struct iterator_group;
}  // namespace symphas::internal

namespace symphas {
enum class iterator_key { SEQ, SEL, GRP, REG };

template <iterator_key type>
struct iterator_policy {};

using sequential_iterator_policy = iterator_policy<iterator_key::SEQ>;
using selection_iterator_policy = iterator_policy<iterator_key::SEL>;
using group_iterator_policy = iterator_policy<iterator_key::GRP>;
using region_iterator_policy = iterator_policy<iterator_key::REG>;

inline constexpr auto it_seq = sequential_iterator_policy{};
inline constexpr auto it_sel = selection_iterator_policy{};
inline constexpr auto it_grp = group_iterator_policy{};
inline constexpr auto it_reg = region_iterator_policy{};

template <typename specialized_iterator, typename iterator_category,
          typename value_type, typename reference_type,
          typename difference_type, typename pointer, typename reference>
struct iterator_type_impl;

template <typename T>
struct stride_type
    : iterator_type_impl<stride_type<T>, std::forward_iterator_tag, T, T&,
                         std::ptrdiff_t, T*, T&> {
 protected:
  stride_type(T* data, len_type stride, len_type len, std::ptrdiff_t ptr)
      : data{data}, stride{stride}, len{len}, ptr{ptr} {}

 public:
  stride_type(T* data, len_type stride, len_type len)
      : data{data}, stride{stride}, len{len}, ptr{0} {}
  stride_type(T* data, len_type len)
      : data{data}, stride{1}, len{len}, ptr{0} {}

  stride_type(T* data, len_type stride, size_t len)
      : stride_type(data, stride, len_type(len)) {}
  stride_type(T* data, size_t len) : stride_type(data, len_type(len)) {}

  auto begin() const { return stride_type<T>(data, stride, len, 0); }

  auto end() const { return stride_type<T>(data, stride, len, len); }

  const auto& operator[](iter_type i) const { return data[(ptr + i) * stride]; }

  auto& operator[](iter_type i) { return data[(ptr + i) * stride]; }

  const auto& operator*() const { return data[ptr * stride]; }

  auto& operator*() { return data[ptr * stride]; }

  T* data;
  len_type stride;
  len_type len;
  std::ptrdiff_t ptr;
};

template <typename T>
stride_type(T*, size_t) -> stride_type<T>;
template <typename T>
stride_type(T*, len_type) -> stride_type<T>;

//! An iterator used to evaluate an expression on its underlying data.
/*!
 * Implements the forward iterator functionality, in order
 * to support parallelization using the standard library.
 *
 * \tparam G The expression type which is evaluated.
 */
template <typename G>
struct data_iterator
    : iterator_type_impl<data_iterator<G>, std::random_access_iterator_tag,
                         typename symphas::internal::data_value_type<G>::type,
                         typename symphas::internal::data_value_type<G>::ref,
                         std::ptrdiff_t,
                         typename symphas::internal::data_value_type<G>::type*,
                         typename symphas::internal::data_value_type<G>::ref> {
  using difference_type = std::ptrdiff_t;

  //! Create an iterator starting at the given position.
  /*!
   * Create an iterator over an expression starting at the given
   * position.
   *
   * \param pos The index of the underlying data in the expression
   * which is the first index in the iterator.
   */
  template <typename specialized_difference =
                symphas::internal::iterator_difference_type<G>>
  data_iterator(symphas::internal::iterator_difference_type_impl<
                specialized_difference> const& ptr)
      : ptr{*static_cast<specialized_difference const*>(&ptr)} {}

  data_iterator(difference_type pos = {}) : ptr{pos} {}

  //! Create an iterator starting at the given position.
  /*!
   * Create an iterator over an expression starting at the given
   * position. The expression is explicitly given.
   *
   * \param data The expression for this iterator.
   * \param pos The index of the underlying data in the expression
   * which is the first index in the iterator.
   */
  explicit data_iterator(G& data, difference_type pos = {}) : ptr{&data, pos} {}

  //! Create an iterator starting at the given position.
  /*!
   * Create an iterator over an expression starting at the given
   * position. The expression is explicitly given.
   *
   * \param data The expression for this iterator.
   * \param pos The index of the underlying data in the expression
   * which is the first index in the iterator.
   */
  explicit data_iterator(G* data, difference_type pos = {}) : ptr{data, pos} {}

  data_iterator(data_iterator<G> const& other) : data_iterator(other.ptr) {}
  data_iterator(data_iterator<G>&& other) : data_iterator(other.ptr) {}
  data_iterator<G>& operator=(data_iterator<G> other) {
    using std::swap;
    swap(ptr, other.ptr);
    return *this;
  }

  //! Dereference the iterator.
  inline decltype(auto) operator*() {
    iter_type n = (iter_type)ptr.pos;
    return symphas::internal::data_value_type<G>{}(ptr.ptr, n);
  };

  //! Dereference past the iterator.
  inline decltype(auto) operator[](difference_type given_pos) {
    return symphas::internal::data_value_type<G>{}(
        ptr.ptr, iter_type(ptr.pos + given_pos));
  }

  //! Member access of the iterated expression.
  inline G* operator->() { return ptr.ptr; };

  symphas::internal::iterator_difference_type<G> ptr;
};

template <typename G>
data_iterator(G*) -> data_iterator<G>;

template <typename G>
data_iterator(G*, int) -> data_iterator<G>;

//! An iterator used to evaluate an expression on its underlying data.
/*!
 * Implements the forward iterator functionality, in order
 * to support parallelization using the standard library.
 *
 * \tparam G The expression type which is evaluated.
 */
template <typename G>
struct data_iterator_selection
    : iterator_type_impl<
          data_iterator_selection<G>, std::random_access_iterator_tag,
          typename symphas::internal::data_value_type<G>::type,
          typename symphas::internal::data_value_type<G>::ref, std::ptrdiff_t,
          typename symphas::internal::data_value_type<G>::type*,
          typename symphas::internal::data_value_type<G>::ref> {
  using difference_type =
      std::ptrdiff_t;  // iterator_selection_difference_type<G>;

  //! Create an iterator starting at the given position.
  /*!
   * Create an iterator over an expression starting at the given
   * position.
   *
   * \param pos The index of the underlying data in the expression
   * which is the first index in the iterator.
   */
  template <typename specialized_difference =
                symphas::internal::iterator_selection_difference_type<G>>
  data_iterator_selection(symphas::internal::iterator_difference_type_impl<
                          specialized_difference> const& ptr)
      : ptr{*static_cast<specialized_difference const*>(&ptr)} {}

  data_iterator_selection(difference_type pos = {}) : ptr{pos} {}

  //! Create an iterator starting at the given position.
  /*!
   * Create an iterator over an expression starting at the given
   * position. The expression is explicitly given.
   *
   * \param data The expression for this iterator.
   * \param pos The index of the underlying data in the expression
   * which is the first index in the iterator.
   */
  explicit data_iterator_selection(G& data, const iter_type* iters,
                                   difference_type pos = 0)
      : ptr{&data, iters, pos} {}

  //! Create an iterator starting at the given position.
  /*!
   * Create an iterator over an expression starting at the given
   * position. The expression is explicitly given.
   *
   * \param data The expression for this iterator.
   * \param pos The index of the underlying data in the expression
   * which is the first index in the iterator.
   */
  explicit data_iterator_selection(G* data, const iter_type* iters,
                                   difference_type pos = 0)
      : ptr{data, iters, pos} {}

  //! Create an iterator starting at the given position.
  /*!
   * Create an iterator over an expression starting at the given
   * position. The expression is explicitly given.
   *
   * \param data The expression for this iterator.
   * \param pos The index of the underlying data in the expression
   * which is the first index in the iterator.
   */
  template <size_t D>
  explicit data_iterator_selection(G& data,
                                   grid::region_index_list<D> const& list,
                                   difference_type pos = 0)
      : data_iterator_selection(data, list.iters, pos) {}

  //! Create an iterator starting at the given position.
  /*!
   * Create an iterator over an expression starting at the given
   * position. The expression is explicitly given.
   *
   * \param data The expression for this iterator.
   * \param pos The index of the underlying data in the expression
   * which is the first index in the iterator.
   */
  template <size_t D>
  explicit data_iterator_selection(G* data,
                                   grid::region_index_list<D> const& list,
                                   difference_type pos = 0)
      : data_iterator_selection(data, list.iters, pos) {}

  data_iterator_selection(data_iterator_selection<G> const& other)
      : data_iterator_selection(other.ptr) {}
  data_iterator_selection(data_iterator_selection<G>&& other)
      : data_iterator_selection(other.ptr) {}
  data_iterator_selection<G>& operator=(data_iterator_selection<G> other) {
    using std::swap;
    swap(ptr, other.ptr);
    return *this;
  }

  //! Dereference the iterator.
  inline decltype(auto) operator*() {
    return symphas::internal::data_value_type<G>{}(ptr.ptr, ptr.iters[ptr.pos]);
  };

  //! Dereference past the iterator.
  inline decltype(auto) operator[](difference_type given_pos) {
    return symphas::internal::data_value_type<G>{}(
        ptr.ptr, ptr.iters[ptr.pos + given_pos]);
  }

  //! Member access of the iterated expression.
  inline G* operator->() { return ptr.ptr; };

  symphas::internal::iterator_selection_difference_type<G> ptr;
};

//! An iterator used to evaluate an expression on its underlying data.
/*!
 * Implements the forward iterator functionality, in order
 * to support parallelization using the standard library.
 *
 * \tparam G The expression type which is evaluated.
 */
template <typename G, size_t D>
struct data_iterator_region
    : iterator_type_impl<
          data_iterator_region<G, D>, std::random_access_iterator_tag,
          typename symphas::internal::data_value_type<G>::type,
          typename symphas::internal::data_value_type<G>::ref, std::ptrdiff_t,
          typename symphas::internal::data_value_type<G>::type*,
          typename symphas::internal::data_value_type<G>::ref> {
  // using iterator_category = typename iterator_policy_impl::iterator_category;
  // using value_type = typename iterator_policy_impl::value_type;
  // using reference_type = typename iterator_policy_impl::reference_type;
  // using difference_type = typename iterator_policy_impl::difference_type;
  // using pointer = typename iterator_policy_impl::pointer;
  // using reference = typename iterator_policy_impl::reference;

  using difference_type = std::ptrdiff_t;

  //! Create an iterator starting at the given position.
  /*!
   * Create an iterator over an expression starting at the given
   * position.
   *
   * \param pos The index of the underlying data in the expression
   * which is the first index in the iterator.
   */
  template <typename specialized_difference =
                symphas::internal::iterator_region_difference_type<G, D>>
  data_iterator_region(symphas::internal::iterator_difference_type_impl<
                       specialized_difference> const& ptr)
      : ptr{*static_cast<specialized_difference const*>(&ptr)} {}

  data_iterator_region(difference_type pos = {}) : ptr{pos} {}

  //! Create an iterator starting at the given position.
  /*!
   * Create an iterator over an expression starting at the given
   * position. The expression is explicitly given.
   *
   * \param data The expression for this iterator.
   * \param pos The index of the underlying data in the expression
   * which is the first index in the iterator.
   */
  explicit data_iterator_region(G& data,
                                grid::region_interval<D> const& interval,
                                difference_type pos = 0)
      : ptr{&data, interval, pos} {}

  //! Create an iterator starting at the given position.
  /*!
   * Create an iterator over an expression starting at the given
   * position. The expression is explicitly given.
   *
   * \param data The expression for this iterator.
   * \param pos The index of the underlying data in the expression
   * which is the first index in the iterator.
   */
  explicit data_iterator_region(G* data,
                                grid::region_interval<D> const& interval,
                                difference_type pos = 0)
      : ptr{data, interval, pos} {}

  explicit data_iterator_region(G* data, grid::region_extent<D> const& region,
                                difference_type pos = 0)
      : ptr{data, region, pos} {}
  explicit data_iterator_region(G& data, grid::region_extent<D> const& region,
                                difference_type pos = 0)
      : ptr{data, region, pos} {}

  data_iterator_region(data_iterator_region<G, D> const& other)
      : data_iterator_region(other.ptr) {}
  data_iterator_region(data_iterator_region<G, D>&& other)
      : data_iterator_region(other.ptr) {}
  data_iterator_region<G, D>& operator=(data_iterator_region<G, D> other) {
    using std::swap;
    swap(ptr, other.ptr);
    return *this;
  }

  //! Dereference the iterator.
  inline decltype(auto) operator*() {
    return symphas::internal::data_value_type<G>{}(ptr.ptr, *ptr);
  };

  //! Dereference past the iterator.
  inline decltype(auto) operator[](difference_type given_pos) {
    return symphas::internal::data_value_type<G>{}(ptr.ptr, ptr[given_pos]);
  }

  //! Member access of the iterated expression.
  inline G* operator->() { return ptr.ptr; };

  symphas::internal::iterator_region_difference_type<G, D> ptr;
};

//! An iterator used to evaluate an expression on its underlying data.
/*!
 * Implements the forward iterator functionality, in order
 * to support parallelization using the standard library.
 *
 * \tparam G The expression type which is evaluated.
 */
template <typename G, size_t D>
struct data_iterator_group
    : iterator_type_impl<
          data_iterator_group<G, D>, std::random_access_iterator_tag,
          typename symphas::internal::data_value_type<G>::type,
          typename symphas::internal::data_value_type<G>::ref, std::ptrdiff_t,
          typename symphas::internal::data_value_type<G>::type*,
          symphas::internal::iterator_group<G, D>> {
  using difference_type = std::ptrdiff_t;

  //! Create an iterator starting at the given position.
  /*!
   * Create an iterator over an expression starting at the given
   * position.
   *
   * \param pos The index of the underlying data in the expression
   * which is the first index in the iterator.
   */
  template <typename specialized_difference =
                symphas::internal::iterator_region_difference_type<G, D>>
  data_iterator_group(symphas::internal::iterator_difference_type_impl<
                      specialized_difference> const& ptr)
      : ptr{*static_cast<specialized_difference const*>(&ptr)} {}

  data_iterator_group(difference_type pos = {}) : ptr{pos} {}

  //! Create an iterator starting at the given position.
  /*!
   * Create an iterator over an expression starting at the given
   * position. The expression is explicitly given.
   *
   * \param data The expression for this iterator.
   * \param pos The index of the underlying data in the expression
   * which is the first index in the iterator.
   */
  explicit data_iterator_group(G& data,
                               grid::region_interval<D> const& interval,
                               difference_type pos = 0)
      : ptr{&data, interval, pos} {}

  //! Create an iterator starting at the given position.
  /*!
   * Create an iterator over an expression starting at the given
   * position. The expression is explicitly given.
   *
   * \param data The expression for this iterator.
   * \param pos The index of the underlying data in the expression
   * which is the first index in the iterator.
   */
  explicit data_iterator_group(G* data,
                               grid::region_interval<D> const& interval,
                               difference_type pos = 0)
      : ptr{data, interval, pos} {}

  data_iterator_group(data_iterator_group<G, D> const& other)
      : data_iterator_group(other.ptr) {}
  data_iterator_group(data_iterator_group<G, D>&& other)
      : data_iterator_group(other.ptr) {}
  data_iterator_group<G, D>& operator=(data_iterator_group<G, D> other) {
    using std::swap;
    swap(ptr, other.ptr);
    return *this;
  }

  //! Dereference the iterator.
  inline decltype(auto) operator*() {
    return symphas::internal::iterator_group(ptr);
  };

  //! Dereference past the iterator.
  inline decltype(auto) operator[](difference_type given_pos) {
    return symphas::internal::iterator_group(ptr, given_pos);
  }

  //! Member access of the iterated expression.
  inline G* operator->() { return ptr.ptr; };

  symphas::internal::iterator_group_difference_type<G, D> ptr;
};

}  // namespace symphas

namespace symphas {

template <typename specialized_iterator, typename iterator_category_t,
          typename value_type_t, typename reference_type_t,
          typename difference_type_t, typename pointer_t, typename reference_t>
struct iterator_type_impl {
  using iterator_category = iterator_category_t;
  using value_type = value_type_t;
  using reference_type = reference_type_t;
  using difference_type = difference_type_t;
  using pointer = pointer_t;
  using reference = reference_t;

  //! Prefix increment, returns itself.
  specialized_iterator& operator++() {
    ++(cast().ptr);
    return cast();
  }

  //! Postfix increment, return a copy before the increment.
  specialized_iterator operator++(int) {
    specialized_iterator it = cast();
    ++(cast().ptr);
    return it;
  }

  //! Prefix decrement, returns itself.
  specialized_iterator& operator--() {
    --(cast().ptr);
    return cast();
  }

  //! Postfix decrement, return a copy before the increment.
  specialized_iterator operator--(int) {
    specialized_iterator it = cast();
    --(cast().ptr);
    return it;
  }

  specialized_iterator& operator+=(difference_type offset) {
    cast().ptr += offset;
    return cast();
  }

  specialized_iterator& operator-=(difference_type offset) {
    cast().ptr -= offset;
    return cast();
  }

  specialized_iterator& operator+=(int offset) {
    return this->operator+=(difference_type(offset));
  }

  specialized_iterator& operator-=(int offset) {
    return this->operator-=(difference_type(offset));
  }

  //! Equality comparison with another iterator.
  /*!
   * Equality comparison with another iterator.
   * Compares the current position.
   */
  bool operator==(specialized_iterator const& other) const {
    return cast().ptr == other.ptr;
  }

  //! Inequality comparison with another iterator.
  /*!
   * Inequality comparison with another iterator.
   * Compares the current position.
   */
  bool operator!=(specialized_iterator const& other) const {
    return !(cast() == other);
  }

  //! Comparison with another iterator.
  /*!
   * Greater than comparison with another iterator.
   * Compares the current position.
   */
  bool operator>(specialized_iterator const& other) const {
    return cast().ptr > other.ptr;
  }

  //! Comparison with another iterator.
  /*!
   * Less than comparison with another iterator.
   * Compares the current position.
   */
  bool operator<(specialized_iterator const& other) const {
    return other > cast();
  }

  //! Comparison with another iterator.
  /*!
   * Greater than or equal to comparison with another iterator.
   * Compares the current position.
   */
  bool operator>=(specialized_iterator const& other) const {
    return !(cast() < other);
  }

  //! Comparison with another iterator.
  /*!
   * Less than or equal to comparison with another iterator.
   * Compares the current position.
   */
  bool operator<=(specialized_iterator const& other) const {
    return !(cast() > other);
  }

  //! Convertible to the difference type of two iterators.
  operator difference_type() const { return cast().ptr; }

  //! Convertible to the difference type of two iterators.
  operator difference_type&() { return cast().ptr; }

  //! Convertible to the difference type of two iterators.
  operator const difference_type&() const { return cast().ptr; }

  //! Add two iterators.
  difference_type operator+(specialized_iterator const& rhs) const {
    return difference_type(cast().ptr + rhs.ptr);
  }

  //! Subtract two iterators.
  difference_type operator-(specialized_iterator const& rhs) const {
    return difference_type(cast().ptr - rhs.ptr);
  }

  //! Add an offset from the iterator.
  // template<typename G0>
  specialized_iterator operator+(difference_type offset) const {
    specialized_iterator it = cast();
    return it += offset;
  }

  //! Subtract an offset from the iterator.
  // template<typename G0>
  specialized_iterator operator-(difference_type offset) const {
    specialized_iterator it = cast();
    return it -= offset;
  }

  //! Add an offset from the iterator.
  specialized_iterator operator+(iter_type offset) const {
    return (cast()) + difference_type(offset);
  }

  //! Subtract an offset from the iterator.
  specialized_iterator operator-(iter_type offset) const {
    return (cast())-difference_type(offset);
  }

  //! Add an offset from the left hand side to an iterator.
  friend specialized_iterator operator+(difference_type const& offset,
                                        specialized_iterator const& rhs) {
    auto result(rhs);
    rhs.ptr = rhs.ptr + offset;
    return result;
  }

  //! Subtract an offset from the left hand side to an iterator.
  friend specialized_iterator operator-(difference_type const& offset,
                                        specialized_iterator const& rhs) {
    auto result(rhs);
    rhs.ptr = -rhs.ptr + offset;
    return result;
  }

  specialized_iterator& cast() {
    return *static_cast<specialized_iterator*>(this);
  }

  const specialized_iterator& cast() const {
    return *static_cast<specialized_iterator const*>(this);
  }
};

}  // namespace symphas

namespace grid {

//! Obtains the iterable_domain from the Block compatible instance.
template <typename T>
auto get_iterable_domain(Block<T> const& data) {
  return region_size(data.len);
}

//! Obtains the iterable_domain from the Block compatible instance.
template <typename T, size_t D>
auto get_iterable_domain(Grid<T, D> const& data) {
  return region_size(data.len);
}

//! Obtains the iterable_domain from the Block compatible instance.
template <typename T, size_t D>
auto get_iterable_domain(BoundaryGrid<T, D> const& data) {
  region_interval<D> region(data.dims);
  for (iter_type i = 0; i < D; ++i) {
    region[i][0] = BOUNDARY_DEPTH;
    region[i][1] = data.dims[i] - BOUNDARY_DEPTH;
  }
  return region;
}

//! Obtains the iterable_domain from the Block compatible instance.
template <typename T, size_t D>
auto get_iterable_domain(RegionalGrid<T, D> const& data) {
  region_interval<D> region(data.dims);
  if (data.region.len > 0) {
    for (iter_type i = 0; i < D; ++i) {
      region[i][0] = data.region.origin[i] + data.region.boundary_size;
      region[i][1] = data.region.origin[i] + data.region.dims[i] -
                     data.region.boundary_size;
    }
  } else {
    for (iter_type i = 0; i < D; ++i) {
      region[i][0] = data.region.origin[i];
      region[i][1] = data.region.origin[i];
    }
  }
  return region;
}

//! Obtains the iterable_domain from the Block compatible instance.
template <typename T, size_t D>
auto get_data_domain(Grid<T, D> const& data) {
  return region_size(data.len);
}

//! Obtains the iterable_domain from the Block compatible instance.
template <typename T, size_t D>
auto get_data_domain(BoundaryGrid<T, D> const& data) {
  region_interval<D> region(data.dims);
  for (iter_type i = 0; i < D; ++i) {
    region[i][0] = 0;
    region[i][1] = data.dims[i];
  }
  return region;
}

//! Obtains the iterable_domain from the Block compatible instance.
template <typename T, size_t D>
auto get_data_domain(RegionalGrid<T, D> const& data) {
  region_interval<D> region(data.dims);
  for (iter_type i = 0; i < D; ++i) {
    region[i][0] = data.region.origin[i];
    region[i][1] = data.region.origin[i] + data.region.dims[i];
  }
  return region;
}

}  // namespace grid

#ifdef USING_CUDA
namespace grid {

template <typename T, size_t D>
auto get_iterable_domain(GridCUDA<T, D> const& data);

template <typename T, size_t D>
auto get_iterable_domain(BoundaryGridCUDA<T, D> const& data);

template <typename T, size_t D>
auto get_iterable_domain(RegionalGridCUDA<T, D> const& data);

template <typename T, size_t D>
auto get_data_domain(GridCUDA<T, D> const& data);

template <typename T, size_t D>
auto get_data_domain(BoundaryGridCUDA<T, D> const& data);

template <typename T, size_t D>
auto get_data_domain(RegionalGridCUDA<T, D> const& data);

}  // namespace grid
#endif
