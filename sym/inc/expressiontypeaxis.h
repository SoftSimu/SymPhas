
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
 * MODULE:  lib
 * PURPOSE: Defines the wavevector object, used for example, to construct
 * a Gaussian smoothing kernel.
 *
 * ***************************************************************************
 */

#pragma once

#include "expressions.h"

struct TimeValue {
  TimeValue() : time{nullptr} {}
  TimeValue(const double* time) : time{time} {}

  double get_time() const { return *time; }

 protected:
  const double* time;
};

template <size_t D, Axis ax>
struct GridAxis {
  len_type dims0[D];
  len_type dims[D];
  symphas::interval_element_type v;

  //! Construct data for computing the axis position in a grid.
  /*!
   * Construct data for computing the axis position in a grid. Based on using
   * the spatial discretization and the system dimensions to determine
   * the position in the grid from a flattened index.
   *
   * \param dims0 The interior dimensions which are computed by the
   * expression evaluation.
   * \param dims The whole system dimensions.
   * \param left The first endpoint of the axis.
   * \param right The last endpoint of the axis.
   */
  GridAxis(len_type const (&dims0)[D], len_type const (&dims)[D], double left,
           double right)
      : GridAxis() {
    std::copy(dims0, dims0 + D, this->dims0);
    std::copy(dims, dims + D, this->dims);
    v.set_count(left, right, dims0[symphas::axis_to_index(ax)]);

    for (iter_type i = 0; i < D; ++i) {
      this->dims0[i] = (dims[i] - dims0[i]) / 2;
    }
  }

  GridAxis(len_type const (&dims)[D]) : GridAxis(dims, dims, 0, 1) {}

  GridAxis(len_type const (&dims)[D], symphas::interval_data_type const& v)
      : GridAxis(dims, dims, v.at(ax).left(), v.at(ax).right()) {}

 protected:
  GridAxis() : dims{0}, v{} {}
};

namespace expr {

template <Axis ax>
using axis_var_t = OpTerms<OpIdentity, Term<GridAxis<1, ax>, 1>>;

//! Specialization based on expr::BaseData.
template <>
struct BaseData<TimeValue> {
  static auto get(TimeValue const& data, iter_type) { return data.get_time(); }
  static auto get(TimeValue const& data) { return data; }
  static auto get(TimeValue& data, iter_type) { return data.get_time(); }
  static auto get(TimeValue& data) { return data; }
};

//! Specialization based on expr::BaseData.
template <>
struct BaseData<GridAxis<1, Axis::X>> {
  static axis_coord_t get(GridAxis<1, Axis::X> const& data, iter_type n) {
    return data.v.left() + (n - data.dims0[0]) * data.v.width();
  }
  static auto get(GridAxis<1, Axis::X> const& data) { return data; }
  static axis_coord_t get(GridAxis<1, Axis::X>& data, iter_type n) {
    return data.v.left() + (n - data.dims0[0]) * data.v.width();
  }
  static auto get(GridAxis<1, Axis::X>& data) { return data; }
};

//! Specialization based on expr::BaseData.
template <>
struct BaseData<GridAxis<2, Axis::X>> {
  static axis_coord_t get(GridAxis<2, Axis::X> const& data, iter_type n) {
    iter_type x = (n % data.dims[0]) - data.dims0[0];
    return data.v.left() + x * data.v.width();
  }
  static auto get(GridAxis<2, Axis::X> const& data) { return data; }
  static axis_coord_t get(GridAxis<2, Axis::X>& data, iter_type n) {
    iter_type x = (n % data.dims[0]) - data.dims0[0];
    return data.v.left() + x * data.v.width();
  }
  static auto get(GridAxis<2, Axis::X>& data) { return data; }
};

//! Specialization based on expr::BaseData.
template <>
struct BaseData<GridAxis<3, Axis::X>> {
  static axis_coord_t get(GridAxis<3, Axis::X> const& data, iter_type n) {
    iter_type x = (n % data.dims[0]) - data.dims0[0];
    return data.v.left() + x * data.v.width();
  }
  static auto get(GridAxis<3, Axis::X> const& data) { return data; }
  static axis_coord_t get(GridAxis<3, Axis::X>& data, iter_type n) {
    iter_type x = (n % data.dims[0]) - data.dims0[0];
    return data.v.left() + x * data.v.width();
  }
  static auto get(GridAxis<3, Axis::X>& data) { return data; }
};

//! Specialization based on expr::BaseData.
template <>
struct BaseData<GridAxis<2, Axis::Y>> {
  static axis_coord_t get(GridAxis<2, Axis::Y> const& data, iter_type n) {
    iter_type y = (n / data.dims[0]) - data.dims0[1];
    return data.v.left() + y * data.v.width();
  }
  static auto get(GridAxis<2, Axis::Y> const& data) { return data; }
  static axis_coord_t get(GridAxis<2, Axis::Y>& data, iter_type n) {
    iter_type y = (n / data.dims[0]) - data.dims0[1];
    return data.v.left() + y * data.v.width();
  }
  static auto get(GridAxis<2, Axis::Y>& data) { return data; }
};

//! Specialization based on expr::BaseData.
template <>
struct BaseData<GridAxis<3, Axis::Y>> {
  static axis_coord_t get(GridAxis<3, Axis::Y> const& data, iter_type n) {
    iter_type y = ((n / data.dims[0]) % data.dims[1]) - data.dims0[1];
    return data.v.left() + y * data.v.width();
  }
  static auto get(GridAxis<3, Axis::Y> const& data) { return data; }
  static axis_coord_t get(GridAxis<3, Axis::Y>& data, iter_type n) {
    iter_type y = ((n / data.dims[0]) % data.dims[1]) - data.dims0[1];
    return data.v.left() + y * data.v.width();
  }
  static auto get(GridAxis<3, Axis::Y>& data) { return data; }
};

//! Specialization based on expr::BaseData.
template <>
struct BaseData<GridAxis<3, Axis::Z>> {
  static axis_coord_t get(GridAxis<3, Axis::Z> const& data, iter_type n) {
    iter_type z = (n / (data.dims[0] * data.dims[1])) - data.dims0[2];
    return data.v.left() + z * data.v.width();
  }
  static auto get(GridAxis<3, Axis::Z> const& data) { return data; }
  static axis_coord_t get(GridAxis<3, Axis::Z>& data, iter_type n) {
    iter_type z = (n / (data.dims[0] * data.dims[1])) - data.dims0[2];
    return data.v.left() + z * data.v.width();
  }
  static auto get(GridAxis<3, Axis::Z>& data) { return data; }
};

//! Helper in order to construct the x-axis variable.
/*!
 * The dimensions and interval information from the given data is
 * copied to the axis data object.
 *
 * The endpoints of the interval are defaultly initialized by the
 * GridAxis object.
 *
 * \param data The grid object from which to derive the side length.
 */
template <typename T, size_t D>
auto make_varx(Grid<T, D> const& data) {
  return expr::make_term(GridAxis<D, Axis::X>{data});
}

//! Helper in order to construct the y-axis variable.
/*!
 * The dimensions and interval information from the given data is
 * copied to the axis data object.
 *
 * The endpoints of the interval are defaultly initialized by the
 * GridAxis object.
 *
 * \param data The grid object from which to derive the side length.
 */
template <typename T, size_t D>
auto make_vary(Grid<T, D> const& data) {
  return expr::make_term(GridAxis<D, Axis::Y>{data});
}

//! Helper in order to construct the z-axis variable.
/*!
 * The dimensions and interval information from the given data is
 * copied to the axis data object.
 *
 * The endpoints of the interval are defaultly initialized by the
 * GridAxis object.
 *
 * \param data The grid object from which to derive the side length.
 */
template <typename T, size_t D>
auto make_varz(Grid<T, D> const& data) {
  return expr::make_term(GridAxis<D, Axis::Z>{data});
}

//! Helper in order to construct the x-axis variable.
/*!
 * The dimensions and interval information from the given data is
 * copied to the axis data object. Used to initialize the x axis variable
 * for grids with boundaries, in order to correctly generate the interior
 * axis values.
 *
 * \param dims0 The interior dimensions of the system.
 * \param dims The number of points on each axis.
 * \param left The left endpoint of the interval.
 */
inline auto make_varx(len_type dims0, len_type dims, double left,
                      double right) {
  len_type lenarr0[]{dims0};
  len_type lenarr[]{dims};
  return expr::make_term(GridAxis<1, Axis::X>{lenarr0, lenarr, left, right});
}

//! Helper in order to construct the x-axis variable.
/*!
 * The dimensions and interval information from the given data is
 * copied to the axis data object. Used to initialize the x axis variable
 * for grids with boundaries, in order to correctly generate the interior
 * axis values.
 *
 * \param dims0 The interior dimensions of the system.
 * \param dims The number of points on each axis.
 * \param left The left endpoint of the interval.
 */
inline auto make_varx(len_type const (&dims0)[2], len_type const (&dims)[2],
                      double left, double right) {
  return expr::make_term(GridAxis<2, Axis::X>{dims0, dims, left, right});
}

//! Helper in order to construct the x-axis variable.
/*!
 * The dimensions and interval information from the given data is
 * copied to the axis data object. Used to initialize the x axis variable
 * for grids with boundaries, in order to correctly generate the interior
 * axis values.
 *
 * \param dims0 The interior dimensions of the system.
 * \param dims The number of points on each axis.
 * \param left The left endpoint of the interval.
 */
inline auto make_varx(len_type const (&dims0)[3], len_type const (&dims)[3],
                      double left, double right) {
  return expr::make_term(GridAxis<3, Axis::X>{dims0, dims, left, right});
}

//! Helper in order to construct the y-axis variable.
/*!
 * The dimensions and interval information from the given data is
 * copied to the axis data object. Used to initialize the y axis variable
 * for grids with boundaries, in order to correctly generate the interior
 * axis values.
 *
 * \param dims0 The interior dimensions of the system.
 * \param dims The number of points on each axis.
 * \param left The left endpoint of the interval.
 */
inline auto make_vary(len_type const (&dims0)[2], len_type const (&dims)[2],
                      double left, double right) {
  return expr::make_term(GridAxis<2, Axis::Y>{dims0, dims, left, right});
}

//! Helper in order to construct the y-axis variable.
/*!
 * The dimensions and interval information from the given data is
 * copied to the axis data object. Used to initialize the y-axis variable
 * for grids with boundaries, in order to correctly generate the interior
 * axis values.
 *
 * \param dims0 The interior dimensions of the system.
 * \param dims The number of points on each axis.
 * \param left The left endpoint of the interval.
 */
inline auto make_vary(len_type const (&dims0)[3], len_type const (&dims)[3],
                      double left, double right) {
  return expr::make_term(GridAxis<3, Axis::Y>{dims0, dims, left, right});
}

//! Helper in order to construct the z-axis variable.
/*!
 * The dimensions and interval information from the given data is
 * copied to the axis data object. Used to initialize the z-axis variable
 * for grids with boundaries, in order to correctly generate the interior
 * axis values.
 *
 * \param dims0 The interior dimensions of the system.
 * \param dims The number of points on each axis.
 * \param left The left endpoint of the interval.
 */
inline auto make_varz(len_type const (&dims0)[3], len_type const (&dims)[3],
                      double left, double right) {
  return expr::make_term(GridAxis<3, Axis::Z>{dims0, dims, left, right});
}

//! Helper in order to construct the x-axis variable.
/*!
 * The dimensions and interval information from the given data is
 * copied to the axis data object.
 *
 * \param dims The number of points on each axis.
 * \param left The left endpoint of the interval.
 */
inline auto make_varx(len_type dims, double left, double right) {
  return make_varx(dims, dims, left, right);
}

//! Helper in order to construct the x-axis variable.
/*!
 * The dimensions and interval information from the given data is
 * copied to the axis data object.
 *
 * \param dims The number of points on each axis.
 * \param left The left endpoint of the interval.
 */
inline auto make_varx(len_type const (&dims)[2], double left, double right) {
  return make_varx(dims, dims, left, right);
}

//! Helper in order to construct the x-axis variable.
/*!
 * The dimensions and interval information from the given data is
 * copied to the axis data object.
 *
 * \param dims The number of points on each axis.
 * \param left The left endpoint of the interval.
 */
inline auto make_varx(len_type const (&dims)[3], double left, double right) {
  return make_varx(dims, dims, left, right);
}

//! Helper in order to construct the y-axis variable.
/*!
 * The dimensions and interval information from the given data is
 * copied to the axis data object.
 *
 * \param dims The number of points on each axis.
 * \param left The left endpoint of the interval.
 */
inline auto make_vary(len_type const (&dims)[2], double left, double right) {
  return make_vary(dims, dims, left, right);
}

//! Helper in order to construct the y-axis variable.
/*!
 * The dimensions and interval information from the given data is
 * copied to the axis data object.
 *
 * \param dims The number of points on each axis.
 * \param left The left endpoint of the interval.
 */
inline auto make_vary(len_type const (&dims)[3], double left, double right) {
  return make_vary(dims, dims, left, right);
}

//! Helper in order to construct the z-axis variable.
/*!
 * The dimensions and interval information from the given data is
 * copied to the axis data object.
 *
 * \param dims The number of points on each axis.
 * \param left The left endpoint of the interval.
 */
inline auto make_varz(len_type const (&dims)[3], double left, double right) {
  return make_varz(dims, dims, left, right);
}

//! Helper in order to construct the generalized axis variable.
/*!
 * The dimensions and interval information from the given data is
 * copied to the axis data object. Specializations of this function
 * refer to the make_var* functions. Used to initialize the axis variable
 * for grids with boundaries, in order to correctly generate the interior
 * axis values.
 *
 * \param dims0 The interior dimensions of the system.
 * \param dims The number of points on each axis.
 * \param left The left endpoint of the interval.
 *
 * \tparam N The index of the axis.
 * \tparam D The dimension of the system.
 */
template <Axis ax, size_t D>
auto make_var(len_type const* dims0, len_type const* dims, double left,
              double right) {
  return OpVoid{};
}

//! Helper in order to construct the generalized axis variable.
/*!
 * The dimensions and interval information from the given data is
 * copied to the axis data object. Specializations of this function
 * refer to the make_var* functions.
 *
 * \param dims The number of points on each axis.
 * \param left The left endpoint of the interval.
 *
 * \tparam N The index of the axis.
 * \tparam D The dimension of the system.
 */
template <Axis ax, size_t D>
auto make_var(len_type const* dims, double left, double right) {
  return make_var<ax, D>(dims, dims, left, right);
}

template <>
inline auto make_var<Axis::X, 1>(len_type const* dims0, len_type const* dims,
                                 double left, double right) {
  return make_varx(dims0[0], dims[0], left, right);
}

template <>
inline auto make_var<Axis::X, 2>(len_type const* dims0, len_type const* dims,
                                 double left, double right) {
  len_type dims0_arr[]{dims0[0], dims0[1]};
  len_type dims_arr[]{dims[0], dims[1]};
  return make_varx(dims0_arr, dims_arr, left, right);
}

template <>
inline auto make_var<Axis::X, 3>(len_type const* dims0, len_type const* dims,
                                 double left, double right) {
  len_type dims0_arr[]{dims0[0], dims0[1], dims0[2]};
  len_type dims_arr[]{dims[0], dims[1], dims[2]};
  return make_varx(dims0_arr, dims_arr, left, right);
}

template <>
inline auto make_var<Axis::Y, 2>(len_type const* dims0, len_type const* dims,
                                 double left, double right) {
  len_type dims0_arr[]{dims0[0], dims0[1]};
  len_type dims_arr[]{dims[0], dims[1]};
  return make_vary(dims0_arr, dims_arr, left, right);
}

template <>
inline auto make_var<Axis::Y, 3>(len_type const* dims0, len_type const* dims,
                                 double left, double right) {
  len_type dims0_arr[]{dims0[0], dims0[1], dims0[2]};
  len_type dims_arr[]{dims[0], dims[1], dims[2]};
  return make_vary(dims0_arr, dims_arr, left, right);
}

template <>
inline auto make_var<Axis::Z, 3>(len_type const* dims0, len_type const* dims,
                                 double left, double right) {
  len_type dims0_arr[]{dims0[0], dims0[1], dims0[2]};
  len_type dims_arr[]{dims[0], dims[1], dims[2]};
  return make_varz(dims0_arr, dims_arr, left, right);
}

//! Helper in order to construct the x-axis variable.
/*!
 * The dimensions and interval information from the given data is
 * copied to the axis data object.
 *
 * \param dims The number of points on each axis.
 * \param left The left endpoint of the interval.
 */
inline auto make_varx(len_type dims, symphas::interval_data_type const& vdata) {
  len_type dims_arr[]{dims};
  return make_var<Axis::X, 1>(symphas::grid_info{vdata}.get_dims(), dims_arr,
                              vdata.at(Axis::X).left(),
                              vdata.at(Axis::X).right());
}

//! Helper in order to construct the x-axis variable.
/*!
 * The dimensions and interval information from the given data is
 * copied to the axis data object.
 *
 * \param dims The number of points on each axis.
 * \param left The left endpoint of the interval.
 */
inline auto make_varx(len_type const (&dims)[2],
                      symphas::interval_data_type const& vdata) {
  return make_var<Axis::X, 2>(symphas::grid_info{vdata}.get_dims(), dims,
                              vdata.at(Axis::X).left(),
                              vdata.at(Axis::X).right());
}

//! Helper in order to construct the x-axis variable.
/*!
 * The dimensions and interval information from the given data is
 * copied to the axis data object.
 *
 * \param dims The number of points on each axis.
 * \param left The left endpoint of the interval.
 */
inline auto make_varx(len_type const (&dims)[3],
                      symphas::interval_data_type const& vdata) {
  return make_var<Axis::X, 3>(symphas::grid_info{vdata}.get_dims(), dims,
                              vdata.at(Axis::X).left(),
                              vdata.at(Axis::X).right());
}

//! Helper in order to construct the x-axis variable.
/*!
 * The dimensions and interval information from the given data is
 * copied to the axis data object.
 *
 * \param dims The number of points on each axis.
 * \param left The left endpoint of the interval.
 */
inline auto make_vary(len_type const (&dims)[2],
                      symphas::interval_data_type const& vdata) {
  return make_var<Axis::Y, 2>(symphas::grid_info{vdata}.get_dims(), dims,
                              vdata.at(Axis::Y).left(),
                              vdata.at(Axis::Y).right());
}

//! Helper in order to construct the x-axis variable.
/*!
 * The dimensions and interval information from the given data is
 * copied to the axis data object.
 *
 * \param dims The number of points on each axis.
 * \param left The left endpoint of the interval.
 */
inline auto make_vary(len_type const (&dims)[3],
                      symphas::interval_data_type const& vdata) {
  return make_var<Axis::Y, 3>(symphas::grid_info{vdata}.get_dims(), dims,
                              vdata.at(Axis::Y).left(),
                              vdata.at(Axis::Y).right());
}

//! Helper in order to construct the x-axis variable.
/*!
 * The dimensions and interval information from the given data is
 * copied to the axis data object.
 *
 * \param dims The number of points on each axis.
 * \param left The left endpoint of the interval.
 */
inline auto make_varz(len_type const (&dims)[3],
                      symphas::interval_data_type const& vdata) {
  return make_var<Axis::Z, 3>(symphas::grid_info{vdata}.get_dims(), dims,
                              vdata.at(Axis::Z).left(),
                              vdata.at(Axis::Z).right());
}

//! Helper in order to construct the generalized axis variable.
/*!
 * The dimensions and interval information from the given data is
 * copied to the axis data object. Specializations of this function
 * refer to the make_var* functions.
 *
 * \param dims The number of points on each axis.
 * \param left The left endpoint of the interval.
 *
 * \tparam N The index of the axis.
 * \tparam D The dimension of the system.
 */
template <Axis ax, size_t D>
auto make_var(len_type const* dims, symphas::interval_data_type const& vdata) {
  return make_var<ax, D>(symphas::grid_info{vdata}.get_dims(), dims,
                         vdata.at(ax).left(), vdata.at(ax).right());
}

template <size_t D>
auto make_coords(len_type const* dims,
                 symphas::interval_data_type const& vdata) {
  using namespace expr::symbols;
  if constexpr (D == 1) {
    return std::make_tuple(make_var<Axis::X, 1>(dims, vdata), zero, zero);
  } else if constexpr (D == 2) {
    return std::make_tuple(make_var<Axis::X, 2>(dims, vdata),
                           make_var<Axis::Y, 2>(dims, vdata), zero);
  } else if constexpr (D == 3) {
    return std::make_tuple(make_var<Axis::X, 3>(dims, vdata),
                           make_var<Axis::Y, 3>(dims, vdata),
                           make_var<Axis::Z, 3>(dims, vdata));
  } else {
    return std::make_tuple(zero, zero, zero);
  }
}
}  // namespace expr

//! \cond

ALLOW_COMBINATION((size_t D, Axis ax), (GridAxis<D, ax>))

#ifdef LATEX_PLOT

#define STR_AXIS(C) #C

#else

#define STR_AXIS(C) #C

#endif

DEFINE_SYMBOL_ID((), (GridAxis<1, Axis::X>), return STR_AXIS(x))
DEFINE_SYMBOL_ID((), (GridAxis<2, Axis::X>), return STR_AXIS(x))
DEFINE_SYMBOL_ID((), (GridAxis<2, Axis::Y>), return STR_AXIS(y))
DEFINE_SYMBOL_ID((), (GridAxis<3, Axis::X>), return STR_AXIS(x))
DEFINE_SYMBOL_ID((), (GridAxis<3, Axis::Y>), return STR_AXIS(y))
DEFINE_SYMBOL_ID((), (GridAxis<3, Axis::Z>), return STR_AXIS(z))

DEFINE_SYMBOL_ID((), (TimeValue), return STR_AXIS(t))

//! \endcond
