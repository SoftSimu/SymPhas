
/* ***************************************************************************
 * This file is part of the SymPhas library, a framework for implementing
 * solvers for phase-field problems with compile-time symbolic algebra.
 *
 * Copyright (c) 2018-2020 by Steven A. Silber and Mikko Karttunen
 *
 * SymPhas is free software, which can be redistributed or modified under
 * the terms of the GNU Lesser General Public License (LGPL) as published
 * by the Free Software Foundation; LGPL version 3, or later versions at
 * your choice
 *
 * SymPhas is distributed with the faith that it will be helpful and
 * practical but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser
 * General Public License for more details
 *
 * ***************************************************************************
 *
 * MODULE:  sol
 * PURPOSE:
 *
 * The stencil object defines the three finite difference functions, and
 * uses templates to specialize each one by dimension and order, for the
 * gradient, laplacian, bilaplacian and gradlaplacian.
 *
 * A specialized template for the stencil object in turn specializes each
 * function to take advantage of generic programming.
 *
 * A generic stencil can also be specialized for any stencils which are not
 * included in the specialization of first to fourth order derivatives, and
 * allows higher order derivatives to be programmatically or manually computed
 * by specializing the dimension, accuracy and derivative order.
 *
 * ***************************************************************************
 */

#pragma once

#include "expressionstencils.h"
#include "gridcudaincludes.h"
#include "gridinfo.h"
#include "stencildefs.h"

#define DEFAULT_STENCIL_ACCURACY 2

namespace symphas::internal {
template <template <typename, size_t> typename grid_type, typename data_type,
          size_t dimension>
constexpr bool is_convertible_to_block =
    std::is_convertible<grid_type<data_type, dimension>,
                        Block<data_type>>::value ||
    std::is_same_v<grid_type<data_type, dimension>,
                   GridData<data_type, dimension>>
#ifdef USING_CUDA
    || std::is_convertible<grid_type<data_type, dimension>,
                           BlockCUDA<data_type>>::value ||
    std::is_same_v<grid_type<data_type, dimension>,
                   GridDataCUDA<data_type, dimension>>
#endif
    ;
}  // namespace symphas::internal

/*!
 * \defgroup stencil Finite Difference Stencils
 * @{
 */

// ******************************************************************************************

//!
/*!
 * Basic stencil object on which concrete stencils are inherited. This stencil
 * is used for uniform grids where the spatial discretization is identical
 * across dimensions.
 */
template <typename Sp>
struct Stencil {
  template <typename T, size_t D>
  __host__ __device__ auto laplacian(T *const v,
                                     const len_type (&stride)[D]) const {
    return cast().laplacian(v, stride);
  }
  template <typename T, size_t D>
  __host__ __device__ auto bilaplacian(T *const v,
                                       const len_type (&stride)[D]) const {
    return cast().bilaplacian(v, stride);
  }
  template <typename T, size_t D>
  __host__ __device__ auto gradlaplacian(T *const v,
                                         const len_type (&stride)[D]) const {
    return cast().gradlaplacian(v, stride);
  }

  template <typename T, size_t D>
  __host__ __device__ auto gradient(T *const v,
                                    const len_type (&stride)[D]) const {
    return cast().gradient(v, stride);
  }

  template <Axis ax, size_t O, template <typename, size_t> typename G,
            typename T, size_t D,
            typename = std::enable_if_t<
                symphas::internal::is_convertible_to_block<G, T, D>, int>>
  __host__ __device__ auto applied_generalized_directional_derivative(
      G<T, D> const &grid, iter_type n) const {
    return cast().template apply_directional<ax, O>(&grid[0] + n);
  }

  template <size_t... Os, template <typename, size_t> typename G, typename T,
            size_t D,
            typename = std::enable_if_t<
                symphas::internal::is_convertible_to_block<G, T, D>, int>>
  __host__ __device__ auto applied_generalized_mixed_derivative(
      G<T, D> const &grid, iter_type n) const {
    return cast().template apply_mixed<Os...>(&grid[0] + n);
  }

  template <Axis ax, size_t O, template <typename, size_t> typename G,
            typename T, size_t D,
            typename = std::enable_if_t<
                symphas::internal::is_convertible_to_block<G, T, D>, int>>
  __host__ __device__ auto applied_generalized_derivative(G<T, D> const &grid,
                                                          iter_type n) const {
    len_type stride[D];
    grid::get_stride<ax>(stride, cast().dims);
    return cast().template apply<O>(&grid[0] + n, stride);
  }

  template <template <typename, size_t> typename G, typename T, size_t D,
            typename = std::enable_if_t<
                symphas::internal::is_convertible_to_block<G, T, D>, int>>
  __host__ __device__ auto applied_laplacian(G<T, D> const &grid,
                                             iter_type n) const {
    len_type stride[D];
    grid::get_stride<Axis::X>(stride, cast().dims);
    return laplacian(&grid[0] + n, stride);
  }

  template <template <typename, size_t> typename G, typename T, size_t D,
            typename = std::enable_if_t<
                symphas::internal::is_convertible_to_block<G, T, D>, int>>
  __host__ __device__ auto applied_bilaplacian(G<T, D> const &grid,
                                               iter_type n) const {
    len_type stride[D];
    grid::get_stride<Axis::X>(stride, cast().dims);
    return bilaplacian(&grid[0] + n, stride);
  }

  template <Axis ax, template <typename, size_t> typename G, typename T,
            size_t D,
            typename = std::enable_if_t<
                symphas::internal::is_convertible_to_block<G, T, D>, int>>
  __host__ __device__ auto applied_gradlaplacian(G<T, D> const &grid,
                                                 iter_type n) const {
    len_type stride[D];
    grid::get_stride<ax>(stride, cast().dims);
    return gradlaplacian(&grid[0] + n, stride);
  }

  template <Axis ax, template <typename, size_t> typename G, typename T,
            size_t D,
            typename = std::enable_if_t<
                symphas::internal::is_convertible_to_block<G, T, D>, int>>
  __host__ __device__ auto applied_gradient(G<T, D> const &grid,
                                            iter_type n) const {
    len_type stride[D];
    grid::get_stride<ax>(stride, cast().dims);
    return gradient(&grid[0] + n, stride);
  }

  template <Axis ax, size_t O, typename T, size_t D>
  __device__ __host__ auto applied_generalized_directional_derivative(
      RegionalGrid<T, D> const &grid, iter_type n) const {
    auto v = grid[n];
    if (v.is_valid()) {
      return cast().template apply_directional<ax, O>(v.value,
                                                      grid.region.dims);
    } else {
      return *v.value;
    }
  }

  template <size_t... Os, typename T, size_t D>
  __device__ __host__ auto applied_generalized_mixed_derivative(
      RegionalGrid<T, D> const &grid, iter_type n) const {
    auto v = grid[n];
    if (v.is_valid()) {
      return cast().template apply_mixed<Os...>(v.value, grid.region.dims);
    } else {
      return *v.value;
    }
  }

  template <Axis ax, size_t O, typename T, size_t D>
  __device__ __host__ auto applied_generalized_derivative(
      RegionalGrid<T, D> const &grid, iter_type n) const {
    auto v = grid[n];
    if (v.is_valid()) {
      len_type stride[D];
      grid::get_stride<ax>(stride, grid.region.dims);
      return cast().template apply<O>(v.value, stride);
    } else {
      return *v.value;
    }
  }

  template <typename T, size_t D>
  __device__ __host__ auto applied_laplacian(RegionalGrid<T, D> const &grid,
                                             iter_type n) const {
    auto v = grid[n];
    if (v.is_valid()) {
      return laplacian(v.value, grid.region.stride);
    } else {
      return *v.value;
    }
  }

  template <typename T, size_t D>
  __device__ __host__ auto applied_bilaplacian(RegionalGrid<T, D> const &grid,
                                               iter_type n) const {
    auto v = grid[n];
    if (v.is_valid()) {
      return bilaplacian(v.value, grid.region.stride);
    } else {
      return *v.value;
    }
  }

  template <Axis ax, typename T, size_t D>
  __device__ __host__ auto applied_gradlaplacian(RegionalGrid<T, D> const &grid,
                                                 iter_type n) const {
    auto v = grid[n];
    if (v.is_valid()) {
      len_type stride[D];
      grid::get_stride<ax>(stride, grid.region.dims);
      return gradlaplacian(v.value, stride);
    } else {
      return *v.value;
    }
  }

  template <Axis ax, typename T, size_t D>
  __device__ __host__ auto applied_gradient(RegionalGrid<T, D> const &grid,
                                            iter_type n) const {
    auto v = grid[n];
    if (v.is_valid()) {
      len_type stride[D];
      grid::get_stride<ax>(stride, grid.region.dims);
      return gradient(v.value, stride);
    } else {
      return *v.value;
    }
  }

#ifdef USING_CUDA
  template <Axis ax, size_t O, typename T, size_t D>
  __host__ __device__ auto applied_generalized_directional_derivative(
      RegionalGridCUDA<T, D> const &grid, iter_type n) const;

  template <size_t... Os, typename T, size_t D>
  __host__ __device__ auto applied_generalized_mixed_derivative(
      RegionalGridCUDA<T, D> const &grid, iter_type n) const;

  template <Axis ax, size_t O, typename T, size_t D>
  __host__ __device__ auto applied_generalized_derivative(
      RegionalGridCUDA<T, D> const &grid, iter_type n) const;

  template <typename T, size_t D>
  __host__ __device__ auto applied_laplacian(RegionalGridCUDA<T, D> const &grid,
                                             iter_type n) const;

  template <typename T, size_t D>
  __host__ __device__ auto applied_bilaplacian(
      RegionalGridCUDA<T, D> const &grid, iter_type n) const;

  template <Axis ax, typename T, size_t D>
  __host__ __device__ auto applied_gradlaplacian(
      RegionalGridCUDA<T, D> const &grid, iter_type n) const;

  template <Axis ax, typename T, size_t D>
  __host__ __device__ auto applied_gradient(RegionalGridCUDA<T, D> const &grid,
                                            iter_type n) const;

  template <Axis ax, size_t O, typename T, size_t D>
  __host__ __device__ auto applied_generalized_directional_derivative(
      RegionGridDataCUDA<T, D> const &grid, iter_type n) const;

  template <size_t... Os, typename T, size_t D>
  __host__ __device__ auto applied_generalized_mixed_derivative(
      RegionGridDataCUDA<T, D> const &grid, iter_type n) const;

  template <Axis ax, size_t O, typename T, size_t D>
  __host__ __device__ auto applied_generalized_derivative(
      RegionGridDataCUDA<T, D> const &grid, iter_type n) const;

  template <typename T, size_t D>
  __host__ __device__ auto applied_laplacian(
      RegionGridDataCUDA<T, D> const &grid, iter_type n) const;

  template <typename T, size_t D>
  __host__ __device__ auto applied_bilaplacian(
      RegionGridDataCUDA<T, D> const &grid, iter_type n) const;

  template <Axis ax, typename T, size_t D>
  __host__ __device__ auto applied_gradlaplacian(
      RegionGridDataCUDA<T, D> const &grid, iter_type n) const;

  template <Axis ax, typename T, size_t D>
  __host__ __device__ auto applied_gradient(
      RegionGridDataCUDA<T, D> const &grid, iter_type n) const;
#endif

#ifdef USING_MPI

  template <Axis ax, size_t O, typename T, size_t D>
  auto applied_generalized_directional_derivative(
      RegionalGridMPI<T, D> const &grid, iter_type n) const {
    if (grid.values != nullptr) {
      return applied_generalized_directional_derivative<ax, O>(
          static_cast<RegionalGrid<T, D> const &>(grid), n);
    } else {
      return grid.empty;
    }
  }

  template <size_t... Os, typename T, size_t D>
  auto applied_generalized_mixed_derivative(RegionalGridMPI<T, D> const &grid,
                                            iter_type n) const {
    if (grid.values != nullptr) {
      return applied_generalized_mixed_derivative<Os...>(
          static_cast<RegionalGrid<T, D> const &>(grid), n);
    } else {
      return grid.empty;
    }
  }

  template <Axis ax, size_t O, typename T, size_t D>
  auto applied_generalized_derivative(RegionalGridMPI<T, D> const &grid,
                                      iter_type n) const {
    if (grid.values != nullptr) {
      return applied_generalized_derivative<ax, O>(
          static_cast<RegionalGrid<T, D> const &>(grid), n);
    } else {
      return grid.empty;
    }
  }

  template <typename T, size_t D>
  auto applied_laplacian(RegionalGridMPI<T, D> const &grid, iter_type n) const {
    if (grid.values != nullptr) {
      return applied_laplacian(static_cast<RegionalGrid<T, D> const &>(grid),
                               n);
    } else {
      return grid.empty;
    }
  }

  template <typename T, size_t D>
  auto applied_bilaplacian(RegionalGridMPI<T, D> const &grid,
                           iter_type n) const {
    if (grid.values != nullptr) {
      return applied_bilaplacian(static_cast<RegionalGrid<T, D> const &>(grid),
                                 n);
    } else {
      return grid.empty;
    }
  }

  template <Axis ax, typename T, size_t D>
  auto applied_gradlaplacian(RegionalGridMPI<T, D> const &grid,
                             iter_type n) const {
    if (grid.values != nullptr) {
      return applied_gradlaplacian<ax>(
          static_cast<RegionalGrid<T, D> const &>(grid), n);
    } else {
      return grid.empty;
    }
  }

  template <Axis ax, typename T, size_t D>
  auto applied_gradient(RegionalGridMPI<T, D> const &grid, iter_type n) const {
    if (grid.values != nullptr) {
      return applied_gradient<ax>(static_cast<RegionalGrid<T, D> const &>(grid),
                                  n);
    } else {
      return grid.empty;
    }
  }

#endif

  template <Axis ax, size_t O, typename T, size_t D>
  __host__ __device__ auto applied_generalized_directional_derivative(
      VectorComponentData<ax, T *, D> const &data, iter_type n) const {
    return cast().template apply_directional<ax, O>(data.values + n);
  }

  template <size_t... Os, Axis ax, typename T, size_t D>
  __host__ __device__ auto applied_generalized_mixed_derivative(
      VectorComponentData<ax, T *, D> const &grid, iter_type n) const {
    return cast().template apply_mixed<Os...>(grid.values + n);
  }

  template <Axis ax, size_t O, typename T, size_t D>
  __host__ __device__ auto applied_generalized_derivative(
      VectorComponentData<ax, T *, D> const &data, iter_type n) const {
    len_type stride[D];
    grid::get_stride<ax>(stride, cast().dims);
    return cast().template apply<O>(data.values + n, stride);
  }

  template <Axis ax, typename T, size_t D>
  __host__ __device__ auto applied_laplacian(
      VectorComponentData<ax, T *, D> const &data, iter_type n) const {
    len_type stride[D];
    grid::get_stride<Axis::X>(stride, cast().dims);
    return laplacian(data.values + n, stride);
  }

  template <Axis ax, typename T, size_t D>
  __host__ __device__ auto applied_bilaplacian(
      VectorComponentData<ax, T *, D> const &data, iter_type n) const {
    len_type stride[D];
    grid::get_stride<Axis::X>(stride, cast().dims);
    return bilaplacian(data.values + n, stride);
  }

  template <Axis ax, typename T, size_t D>
  __host__ __device__ auto applied_gradlaplacian(
      VectorComponentData<ax, T *, D> const &data, iter_type n) const {
    len_type stride[D];
    grid::get_stride<ax>(stride, cast().dims);
    return gradlaplacian(data.values + n, stride);
  }

  template <Axis ax, typename T, size_t D>
  __host__ __device__ auto applied_gradient(
      VectorComponentData<ax, T *, D> const &data, iter_type n) const {
    len_type stride[D];
    grid::get_stride<ax>(stride, cast().dims);
    return gradient(data.values + n, stride);
  }

  // -------------------------------------------------------------------------
  // Cross-axis VectorComponentData overloads. The derivative axis `axd`
  // is independent of the component axis `axc`. The original overloads
  // above tied them together (single template param `ax`), which silently
  // refused to match expressions like `gradx(my)` (∂/∂x of the y-component
  // of m); SFINAE failure there fell through to a fallback that returned
  // wrong values. Provide explicit cross-axis overloads.
  // The component axis selects WHICH array (data.values is the right one
  // already, picked by resolve_axis_component<axc>). The derivative axis
  // controls the stride used by the stencil.
  // -------------------------------------------------------------------------
  template <Axis axd, size_t O, Axis axc, typename T, size_t D,
            typename = std::enable_if_t<axd != axc, int>>
  __host__ __device__ auto applied_generalized_directional_derivative(
      VectorComponentData<axc, T *, D> const &data, iter_type n) const {
    return cast().template apply_directional<axd, O>(data.values + n);
  }

  template <Axis axd, size_t O, Axis axc, typename T, size_t D,
            typename = std::enable_if_t<axd != axc, int>>
  __host__ __device__ auto applied_generalized_derivative(
      VectorComponentData<axc, T *, D> const &data, iter_type n) const {
    len_type stride[D];
    grid::get_stride<axd>(stride, cast().dims);
    return cast().template apply<O>(data.values + n, stride);
  }

  template <Axis axd, Axis axc, typename T, size_t D,
            typename = std::enable_if_t<axd != axc, int>>
  __host__ __device__ auto applied_gradlaplacian(
      VectorComponentData<axc, T *, D> const &data, iter_type n) const {
    len_type stride[D];
    grid::get_stride<axd>(stride, cast().dims);
    return gradlaplacian(data.values + n, stride);
  }

  template <Axis axd, Axis axc, typename T, size_t D,
            typename = std::enable_if_t<axd != axc, int>>
  __host__ __device__ auto applied_gradient(
      VectorComponentData<axc, T *, D> const &data, iter_type n) const {
    len_type stride[D];
    grid::get_stride<axd>(stride, cast().dims);
    return gradient(data.values + n, stride);
  }

  template <Axis ax, size_t O, typename T, size_t D>
  __host__ __device__ auto applied_generalized_directional_derivative(
      VectorComponentRegionData<ax, T *, D> const &data, iter_type n) const {
    auto v = data[n];
    if (v.is_valid()) {
      return cast().template apply_directional<ax, O>(v.value,
                                                      data.region.dims);
    } else {
      return *v.value;
    }
  }

  template <size_t... Os, Axis ax, typename T, size_t D>
  __host__ __device__ auto applied_generalized_mixed_derivative(
      VectorComponentRegionData<ax, T *, D> const &grid, iter_type n) const {
    auto v = grid[n];
    if (v.is_valid()) {
      return cast().template apply_mixed<Os...>(v.value, grid.region.dims);
    } else {
      return *v.value;
    }
  }

  template <Axis ax, size_t O, typename T, size_t D>
  __host__ __device__ auto applied_generalized_derivative(
      VectorComponentRegionData<ax, T *, D> const &data, iter_type n) const {
    auto v = data[n];
    if (v.is_valid()) {
      len_type stride[D];
      grid::get_stride<ax>(stride, data.region.dims);
      return cast().template apply<O>(v.value, stride);
    } else {
      return *v.value;
    }
  }

  template <Axis ax, typename T, size_t D>
  __host__ __device__ auto applied_laplacian(
      VectorComponentRegionData<ax, T *, D> const &data, iter_type n) const {
    auto v = data[n];
    if (v.is_valid()) {
      return laplacian(v.value, data.region.stride);
    } else {
      return *v.value;
    }
  }

  template <Axis ax, typename T, size_t D>
  __host__ __device__ auto applied_bilaplacian(
      VectorComponentRegionData<ax, T *, D> const &data, iter_type n) const {
    auto v = data[n];
    if (v.is_valid()) {
      return bilaplacian(v.value, data.region.stride);
    } else {
      return *v.value;
    }
  }

  template <Axis ax, typename T, size_t D>
  __host__ __device__ auto applied_gradlaplacian(
      VectorComponentRegionData<ax, T *, D> const &data, iter_type n) const {
    auto v = data[n];
    if (v.is_valid()) {
      len_type stride[D];
      grid::get_stride<ax>(stride, data.region.dims);
      return gradlaplacian(v.value, stride);
    } else {
      return *v.value;
    }
  }

  template <Axis ax, typename T, size_t D>
  __host__ __device__ auto applied_gradient(
      VectorComponentRegionData<ax, T *, D> const &data, iter_type n) const {
    auto v = data[n];
    if (v.is_valid()) {
      len_type stride[D];
      grid::get_stride<ax>(stride, data.region.dims);
      return gradient(v.value, stride);
    } else {
      return *v.value;
    }
  }

  template <Axis ax, size_t O, typename T, size_t D>
  auto applied_generalized_directional_derivative(GridSymbol<T, D> const &grid,
                                                  iter_type n) const {
    return T{};
  }

  template <size_t... Os, typename T, size_t D>
  auto applied_generalized_mixed_derivative(GridSymbol<T, D> const &grid,
                                            iter_type n) const {
    return T{};
  }

  template <Axis ax, size_t O, typename T, size_t D>
  auto applied_generalized_derivative(GridSymbol<T, D> const &grid,
                                      iter_type n) const {
    return T{};
  }

  template <typename T, size_t D>
  auto applied_laplacian(GridSymbol<T, D> const &grid, iter_type n) const {
    return T{};
  }

  template <typename T, size_t D>
  auto applied_bilaplacian(GridSymbol<T, D> const &grid, iter_type n) const {
    return T{};
  }

  template <Axis ax, typename T, size_t D>
  auto applied_gradlaplacian(GridSymbol<T, D> const &grid, iter_type n) const {
    return T{};
  }

  template <Axis ax, typename T, size_t D>
  auto applied_gradient(GridSymbol<T, D> const &grid, iter_type n) const {
    return T{};
  }

  template <Axis ax, size_t O>
  __host__ __device__ auto applied_generalized_directional_derivative(
      ...) const {
    return symphas::lib::get_identity<scalar_t>();
  }

  template <size_t... Os>
  __host__ __device__ auto applied_generalized_mixed_derivative(...) const {
    return symphas::lib::get_identity<scalar_t>();
  }

  template <Axis ax, size_t O>
  __host__ __device__ auto applied_generalized_derivative(...) const {
    return symphas::lib::get_identity<scalar_t>();
  }

  __host__ __device__ auto applied_laplacian(...) const {
    return symphas::lib::get_identity<scalar_t>();
  }

  __host__ __device__ auto applied_bilaplacian(...) const {
    return symphas::lib::get_identity<scalar_t>();
  }

  template <Axis ax>
  __host__ __device__ auto applied_gradlaplacian(...) const {
    return symphas::lib::get_identity<scalar_t>();
  }

  template <Axis ax>
  __host__ __device__ auto applied_gradient(...) const {
    return symphas::lib::get_identity<scalar_t>();
  }

  __host__ __device__ const Sp &cast() const {
    return (*static_cast<const Sp *>(this));
  }

  __host__ __device__ Sp &cast() { return (*static_cast<Sp *>(this)); }
};

namespace symphas::internal {

#ifdef GENERATE_UNDEFINED_STENCILS_ON

template <size_t OD, size_t OA, size_t D>
using stencil_type =
    std::invoke_result_t<decltype(&expr::get_central_space_stencil<OD, OA, D>)>;

template <size_t OA, typename Seq>
__host__ __device__ auto get_mixed_stencil_type(Seq) {
  return expr::get_central_space_mixed_stencil<OA>(Seq{});
};

template <size_t OA, typename Seq>
using mixed_stencil_type =
    std::invoke_result_t<decltype(&get_mixed_stencil_type<OA, Seq>), Seq>;

template <typename stencil_type>
__host__ __device__ auto get_stencil_apply_type(stencil_type) {
  return symphas::internal::GeneratedStencilApply{stencil_type{}};
};

template <size_t OD, size_t OA, size_t D>
struct stencil_apply_type_impl {
  using type = std::invoke_result_t<
      decltype(&get_stencil_apply_type<stencil_type<OD, OA, D>>),
      stencil_type<OD, OA, D>>;
};

template <size_t OD, size_t OA, size_t D>
using stencil_apply_type = typename stencil_apply_type_impl<OD, OA, D>::type;

template <size_t OA, typename Seq>
struct mixed_stencil_apply_type_impl {
  using type = std::invoke_result_t<
      decltype(&get_stencil_apply_type<mixed_stencil_type<OA, Seq>>),
      mixed_stencil_type<OA, Seq>>;
};

template <size_t OA, typename Seq>
using mixed_stencil_apply_type =
    typename mixed_stencil_apply_type_impl<OA, Seq>::type;

template <size_t OD, size_t OA, typename T>
__host__ __device__ auto apply_generalized_derivative(
    T *const v, double divh, const len_type (&stride)[1]) {
  return stencil_apply_type<OD, OA, 1>{}(v, stride, divh);
}

template <size_t OD, size_t OA, typename T>
__host__ __device__ auto apply_generalized_derivative(
    T *const v, double divh, const len_type (&stride)[2]) {
  return stencil_apply_type<OD, OA, 2>{}(v, stride, divh);
}

template <size_t OD, size_t OA, typename T>
__host__ __device__ auto apply_generalized_derivative(
    T *const v, double divh, const len_type (&stride)[3]) {
  return stencil_apply_type<OD, OA, 3>{}(v, stride, divh);
}

#else

struct no_derivative_message_printed {
  __host__ __device__ no_derivative_message_printed(size_t OD, size_t OA,
                                                    size_t D) {
#ifndef __CUDA_ARCH__
    fprintf(
        SYMPHAS_ERR,
        "no derivative of order %zd accuracy %zd available in dimension %zd\n",
        OD, OA, D);
    exit(565);
#endif
  }
};

template <size_t OD, size_t OA, size_t D>
__host__ __device__ void print_no_derivative_message() {
  no_derivative_message_printed message{OD, OA, D};
}

template <size_t OD, size_t OA, size_t D>
struct stencil_apply_type;

template <size_t OD, size_t OA, typename T>
__host__ __device__ auto apply_generalized_derivative(
    T *const v, double divh, const len_type (&stride)[1]) {
  // Delegate to stencil_apply_type so a static specialization (e.g. for
  // OD=1 axial gradients) can intercept; the primary template falls
  // through to print_no_derivative_message on its own.
  return stencil_apply_type<OD, OA, 1>{}(v, stride, divh);
}

template <size_t OD, size_t OA, typename T>
__host__ __device__ auto apply_generalized_derivative(
    T *const v, double divh, const len_type (&stride)[2]) {
  return stencil_apply_type<OD, OA, 2>{}(v, stride, divh);
}

template <size_t OD, size_t OA, typename T>
__host__ __device__ auto apply_generalized_derivative(
    T *const v, double divh, const len_type (&stride)[3]) {
  return stencil_apply_type<OD, OA, 3>{}(v, stride, divh);
}

template <size_t OD, size_t OA, size_t D>
struct stencil_apply_type {
  template <typename T>
  __host__ __device__ auto operator()(const T *, ...) const {
    print_no_derivative_message<OD, OA, D>();
    return T{};
  }
};

// Statically-registered 1D central-difference stencils (1 axis through
// any-D grid).  Without these, models that decompose grad / div / curl
// per-axis at apply-time fail with "no derivative of order N accuracy 2
// available in dimension 1" unless AVAILABLE_STENCILS_AUTOGENERATION=ON
// (which has its own pre-existing internal bugs).
//
// The 1D stencil's `apply_directional<ax, OD>` calls
// `stencil_apply_type<OD, OA, 1>(v, stride, divh)` where `v` points at
// a sample, `stride` is the linear index step along the chosen axis,
// and `divh = 1/h` (so divh^k = 1/h^k).  These specializations use
// standard 2nd-order accurate central differences.

// Order 1: (f(+h) - f(-h)) / (2h)
template <>
struct stencil_apply_type<1, 2, 1> {
  template <typename T>
  __host__ __device__ auto operator()(T *const v, len_type stride,
                                      double divh) const {
    return 0.5 * divh * (v[stride] - v[-stride]);
  }
  template <typename T>
  __host__ __device__ auto operator()(T *const v,
                                      const len_type (&stride)[1],
                                      double divh) const {
    return 0.5 * divh * (v[stride[0]] - v[-stride[0]]);
  }
};

// Order 2: (f(+h) - 2 f(0) + f(-h)) / h^2
template <>
struct stencil_apply_type<2, 2, 1> {
  template <typename T>
  __host__ __device__ auto operator()(T *const v, len_type stride,
                                      double divh) const {
    const double divh2 = divh * divh;
    return divh2 * (v[stride] - 2.0 * v[0] + v[-stride]);
  }
  template <typename T>
  __host__ __device__ auto operator()(T *const v,
                                      const len_type (&stride)[1],
                                      double divh) const {
    return (*this)(v, stride[0], divh);
  }
};

// Order 3: (-f(-2h) + 2f(-h) - 2f(+h) + f(+2h)) / (2 h^3)
template <>
struct stencil_apply_type<3, 2, 1> {
  template <typename T>
  __host__ __device__ auto operator()(T *const v, len_type stride,
                                      double divh) const {
    const double divh3 = divh * divh * divh;
    return 0.5 * divh3 *
           (-v[-2 * stride] + 2.0 * v[-stride] - 2.0 * v[stride] +
            v[2 * stride]);
  }
  template <typename T>
  __host__ __device__ auto operator()(T *const v,
                                      const len_type (&stride)[1],
                                      double divh) const {
    return (*this)(v, stride[0], divh);
  }
};

// Order 4: (f(-2h) - 4f(-h) + 6f(0) - 4f(+h) + f(+2h)) / h^4
template <>
struct stencil_apply_type<4, 2, 1> {
  template <typename T>
  __host__ __device__ auto operator()(T *const v, len_type stride,
                                      double divh) const {
    const double divh4 = divh * divh * divh * divh;
    return divh4 * (v[-2 * stride] - 4.0 * v[-stride] + 6.0 * v[0] -
                    4.0 * v[stride] + v[2 * stride]);
  }
  template <typename T>
  __host__ __device__ auto operator()(T *const v,
                                      const len_type (&stride)[1],
                                      double divh) const {
    return (*this)(v, stride[0], divh);
  }
};

// 2D-stride single-axis stencils.  These are reached when the generalized
// dispatch routes `applied_generalized_derivative<ax, O>(grid_2D, n)` ->
// `apply<O>(v, stride[2])` -> `stencil_apply_type<O, 2, 2>(v, stride[2], divh)`.
// `grid::get_stride<ax>(stride, dims)` puts the relevant axis step in
// stride[0]; the other entries are inactive for an axial derivative.
template <>
struct stencil_apply_type<1, 2, 2> {
  template <typename T>
  __host__ __device__ auto operator()(T *const v,
                                      const len_type (&stride)[2],
                                      double divh) const {
    return 0.5 * divh * (v[stride[0]] - v[-stride[0]]);
  }
};

template <>
struct stencil_apply_type<2, 2, 2> {
  template <typename T>
  __host__ __device__ auto operator()(T *const v,
                                      const len_type (&stride)[2],
                                      double divh) const {
    const double divh2 = divh * divh;
    const len_type s = stride[0];
    return divh2 * (v[s] - 2.0 * v[0] + v[-s]);
  }
};

template <>
struct stencil_apply_type<3, 2, 2> {
  template <typename T>
  __host__ __device__ auto operator()(T *const v,
                                      const len_type (&stride)[2],
                                      double divh) const {
    const double divh3 = divh * divh * divh;
    const len_type s = stride[0];
    return 0.5 * divh3 *
           (-v[-2 * s] + 2.0 * v[-s] - 2.0 * v[s] + v[2 * s]);
  }
};

template <>
struct stencil_apply_type<4, 2, 2> {
  template <typename T>
  __host__ __device__ auto operator()(T *const v,
                                      const len_type (&stride)[2],
                                      double divh) const {
    const double divh4 = divh * divh * divh * divh;
    const len_type s = stride[0];
    return divh4 * (v[-2 * s] - 4.0 * v[-s] + 6.0 * v[0] -
                    4.0 * v[s] + v[2 * s]);
  }
};

template <size_t OA, typename Seq>
struct mixed_stencil_apply_type;

template <size_t OA, size_t... Os>
struct mixed_stencil_apply_type<OA, std::index_sequence<Os...>> {
  template <typename T, size_t D>
  __host__ __device__ auto operator()(const T *, const len_type (&stride)[D],
                                      ...) const {
    print_no_derivative_message<(0 + ... + Os), OA, D>();
    return T{};
  }
};

// 2D axial-as-mixed: the symbolic layer encodes directional_derivative<X,O>
// as mixed (O, 0) and directional_derivative<Y,O> as mixed (0, O). These
// reduce to plain 1D central differences applied along the appropriate
// stride. OA = 2, axial order O in {1, 2, 3, 4}.
template <>
struct mixed_stencil_apply_type<2, std::index_sequence<1, 0>> {
  template <typename T>
  __host__ __device__ auto operator()(T *const v,
                                      const len_type (&stride)[2],
                                      double divh) const {
    const len_type sx = stride[0];
    return 0.5 * divh * (v[sx] - v[-sx]);
  }
};

template <>
struct mixed_stencil_apply_type<2, std::index_sequence<0, 1>> {
  template <typename T>
  __host__ __device__ auto operator()(T *const v,
                                      const len_type (&stride)[2],
                                      double divh) const {
    const len_type sy = stride[1];
    return 0.5 * divh * (v[sy] - v[-sy]);
  }
};

template <>
struct mixed_stencil_apply_type<2, std::index_sequence<2, 0>> {
  template <typename T>
  __host__ __device__ auto operator()(T *const v,
                                      const len_type (&stride)[2],
                                      double divh) const {
    const double divh2 = divh * divh;
    const len_type sx = stride[0];
    return divh2 * (v[sx] - 2.0 * v[0] + v[-sx]);
  }
};

template <>
struct mixed_stencil_apply_type<2, std::index_sequence<0, 2>> {
  template <typename T>
  __host__ __device__ auto operator()(T *const v,
                                      const len_type (&stride)[2],
                                      double divh) const {
    const double divh2 = divh * divh;
    const len_type sy = stride[1];
    return divh2 * (v[sy] - 2.0 * v[0] + v[-sy]);
  }
};

template <>
struct mixed_stencil_apply_type<2, std::index_sequence<3, 0>> {
  template <typename T>
  __host__ __device__ auto operator()(T *const v,
                                      const len_type (&stride)[2],
                                      double divh) const {
    const double divh3 = divh * divh * divh;
    const len_type sx = stride[0];
    return 0.5 * divh3 *
           (-v[-2 * sx] + 2.0 * v[-sx] - 2.0 * v[sx] + v[2 * sx]);
  }
};

template <>
struct mixed_stencil_apply_type<2, std::index_sequence<0, 3>> {
  template <typename T>
  __host__ __device__ auto operator()(T *const v,
                                      const len_type (&stride)[2],
                                      double divh) const {
    const double divh3 = divh * divh * divh;
    const len_type sy = stride[1];
    return 0.5 * divh3 *
           (-v[-2 * sy] + 2.0 * v[-sy] - 2.0 * v[sy] + v[2 * sy]);
  }
};

template <>
struct mixed_stencil_apply_type<2, std::index_sequence<4, 0>> {
  template <typename T>
  __host__ __device__ auto operator()(T *const v,
                                      const len_type (&stride)[2],
                                      double divh) const {
    const double divh4 = divh * divh * divh * divh;
    const len_type sx = stride[0];
    return divh4 * (v[-2 * sx] - 4.0 * v[-sx] + 6.0 * v[0] -
                    4.0 * v[sx] + v[2 * sx]);
  }
};

template <>
struct mixed_stencil_apply_type<2, std::index_sequence<0, 4>> {
  template <typename T>
  __host__ __device__ auto operator()(T *const v,
                                      const len_type (&stride)[2],
                                      double divh) const {
    const double divh4 = divh * divh * divh * divh;
    const len_type sy = stride[1];
    return divh4 * (v[-2 * sy] - 4.0 * v[-sy] + 6.0 * v[0] -
                    4.0 * v[sy] + v[2 * sy]);
  }
};

// Mixed (1,1) on a 2D grid: cross-derivative d^2 f / (dx dy)
// = (f(+h,+h) - f(+h,-h) - f(-h,+h) + f(-h,-h)) / (4 h^2)
template <>
struct mixed_stencil_apply_type<2, std::index_sequence<1, 1>> {
  template <typename T>
  __host__ __device__ auto operator()(T *const v,
                                      const len_type (&stride)[2],
                                      double divh) const {
    const double divh2 = divh * divh;
    const len_type sxy_pp = stride[0] + stride[1];
    const len_type sxy_pm = stride[0] - stride[1];
    return 0.25 * divh2 *
           (v[sxy_pp] - v[sxy_pm] - v[-sxy_pm] + v[-sxy_pp]);
  }
};

// Mixed (2,2) on a 2D grid: d^4 f / (dx^2 dy^2) via tensor product of
// the 1D 3-point second derivatives.
// = sum_{i,j in {-1,0,+1}} w_i w_j * f(i h_x, j h_y) / (h_x^2 h_y^2)
// with w = (1, -2, 1).
template <>
struct mixed_stencil_apply_type<2, std::index_sequence<2, 2>> {
  template <typename T>
  __host__ __device__ auto operator()(T *const v,
                                      const len_type (&stride)[2],
                                      double divh) const {
    const double divh4 = divh * divh * divh * divh;
    const len_type sx = stride[0];
    const len_type sy = stride[1];
    return divh4 *
           (v[sx + sy] - 2.0 * v[sy] + v[-sx + sy] -
            2.0 * v[sx] + 4.0 * v[0] - 2.0 * v[-sx] +
            v[sx - sy] - 2.0 * v[-sy] + v[-sx - sy]);
  }
};

// Mixed (1,2) on a 2D grid: d^3 f / (dx dy^2)
// = (1/(2 h_x h_y^2)) * [ (f(+x,+y) - 2 f(+x,0) + f(+x,-y))
//                       - (f(-x,+y) - 2 f(-x,0) + f(-x,-y)) ]
template <>
struct mixed_stencil_apply_type<2, std::index_sequence<1, 2>> {
  template <typename T>
  __host__ __device__ auto operator()(T *const v,
                                      const len_type (&stride)[2],
                                      double divh) const {
    const double divh3 = divh * divh * divh;
    const len_type sx = stride[0];
    const len_type sy = stride[1];
    return 0.5 * divh3 *
           ((v[sx + sy] - 2.0 * v[sx] + v[sx - sy]) -
            (v[-sx + sy] - 2.0 * v[-sx] + v[-sx - sy]));
  }
};

// Mixed (2,1) on a 2D grid: d^3 f / (dx^2 dy)
// = (1/(h_x^2 * 2 h_y)) * [ (f(+x,+y) - f(+x,-y))
//                         - 2 (f(0,+y) - f(0,-y))
//                         + (f(-x,+y) - f(-x,-y)) ]
template <>
struct mixed_stencil_apply_type<2, std::index_sequence<2, 1>> {
  template <typename T>
  __host__ __device__ auto operator()(T *const v,
                                      const len_type (&stride)[2],
                                      double divh) const {
    const double divh3 = divh * divh * divh;
    const len_type sx = stride[0];
    const len_type sy = stride[1];
    return 0.5 * divh3 *
           ((v[sx + sy] - v[sx - sy]) -
            2.0 * (v[sy] - v[-sy]) +
            (v[-sx + sy] - v[-sx - sy]));
  }
};

// Mixed (1,3) on a 2D grid: d^4 f / (dx dy^3) via tensor product of
// the 1D 2-point first derivative (weights +/- 1/2) with the 1D 5-point
// third derivative (weights (-1, +2, 0, -2, +1)/2).
// = (1/(4 h_x h_y^3)) * sum_{i,j} w1_i w3_j f(i h_x, j h_y).
template <>
struct mixed_stencil_apply_type<2, std::index_sequence<1, 3>> {
  template <typename T>
  __host__ __device__ auto operator()(T *const v,
                                      const len_type (&stride)[2],
                                      double divh) const {
    const double divh4 = divh * divh * divh * divh;
    const len_type sx = stride[0];
    const len_type sy = stride[1];
    return 0.25 * divh4 *
           ((-v[sx - 2 * sy] + 2.0 * v[sx - sy] -
             2.0 * v[sx + sy] + v[sx + 2 * sy]) -
            (-v[-sx - 2 * sy] + 2.0 * v[-sx - sy] -
             2.0 * v[-sx + sy] + v[-sx + 2 * sy]));
  }
};

// Mixed (3,1) on a 2D grid: d^4 f / (dx^3 dy) via tensor product of
// the 1D 5-point third derivative (weights (-1, +2, 0, -2, +1)/2) with
// the 1D 2-point first derivative (weights +/- 1/2).
// = (1/(4 h_x^3 h_y)) * sum_{i,j} w3_i w1_j f(i h_x, j h_y).
template <>
struct mixed_stencil_apply_type<2, std::index_sequence<3, 1>> {
  template <typename T>
  __host__ __device__ auto operator()(T *const v,
                                      const len_type (&stride)[2],
                                      double divh) const {
    const double divh4 = divh * divh * divh * divh;
    const len_type sx = stride[0];
    const len_type sy = stride[1];
    return 0.25 * divh4 *
           ((-(v[-2 * sx + sy] - v[-2 * sx - sy])) +
            2.0 * (v[-sx + sy] - v[-sx - sy]) -
            2.0 * (v[sx + sy] - v[sx - sy]) +
            (v[2 * sx + sy] - v[2 * sx - sy]));
  }
};

#endif
}  // namespace symphas::internal

namespace symphas::internal {}

//! A generalized derivative stencil.
/*!
 * A generalized derivative object that will produce stencils for any order of
 * derivative with a maximum boundary expansion of 3.
 *
 * \tparam DD Dimension the stencil is applied in.
 * \tparam OA The order of accuracy for the stencil.
 */
template <size_t DD, size_t OA = DEFAULT_STENCIL_ACCURACY>
struct GeneralizedStencil {
  len_type dims[DD];
  double divh;

  GeneralizedStencil() : dims{}, divh{0} {}

  GeneralizedStencil(const len_type *dims, double h) : divh{1. / h} {
    std::copy(dims, dims + DD, this->dims);
  }

  //! Determine the finite difference approximation to the derivative.
  /*!
   * Determine the finite difference approximation to the derivative of
   * the given order.
   *
   * \tparam OD Order of the derivative.
   */
  template <size_t OD, typename T>
  __host__ __device__ auto apply(T *const v,
                                 const len_type (&stride)[DD]) const {
    return symphas::internal::apply_generalized_derivative<OD, OA>(v, divh,
                                                                   stride);
  }

  //! Determine the finite difference approximation to the directional
  //! derivative.
  /*!
   * Determine the finite difference approximation to the directional derivative
   * of the given order.
   *
   * \tparam OD Order of the derivative.
   * \tparam ax Axis of the directional derivative.
   */
  template<Axis ax, size_t OD, typename T/*, typename Stt = std::invoke_result_t<decltype(&expr::get_central_space_stencil<OD, OA, 1>)>,
		typename = std::enable_if_t<((ax == Axis::Z) ? (DD == 3) : (ax == Axis::Y) ? (DD >= 2) : (ax == Axis::X) ? (DD >= 1) : false), int>*/>
	__host__ __device__ auto apply_directional(T* const v) const
	{
    len_type stride = (ax == Axis::Z)   ? (dims[0] * dims[1])
                      : (ax == Axis::Y) ? dims[0]
                                        : 1;
    return symphas::internal::stencil_apply_type<OD, OA, 1>{}(v, stride, divh);
  }

  //! Determine the finite difference approximation to the directional
  //! derivative.
  /*!
   * Determine the finite difference approximation to the directional derivative
   * of the given order.
   *
   * \tparam OD Order of the derivative.
   * \tparam ax Axis of the directional derivative.
   */
  template<Axis ax, size_t OD, typename T/*, typename Stt = std::invoke_result_t<decltype(&expr::get_central_space_stencil<OD, OA, 1>)>,
		typename = std::enable_if_t<((ax == Axis::Z) ? (DD == 3) : (ax == Axis::Y) ? (DD >= 2) : (ax == Axis::X) ? (DD >= 1) : false), int>*/>
	__host__ __device__ auto apply_directional(T* const v, len_type const (&override_dims)[DD]) const
	{
    len_type stride = (ax == Axis::Z)   ? (override_dims[0] * override_dims[1])
                      : (ax == Axis::Y) ? override_dims[0]
                                        : 1;
    return symphas::internal::stencil_apply_type<OD, OA, 1>{}(v, stride, divh);
  }

  //! Determine the finite difference approximation to the mixed derivative.
  /*!
   * Determine the finite difference approximation to the mixed derivative of
   * the given order.
   *
   * \tparam OD Order of the derivative.
   * \tparam ax Axis of the directional derivative.
   */
  template<size_t OD1, size_t OD2, size_t OD3, typename T>
	__host__ __device__ auto apply_mixed(T* const v) const
	{
    len_type stride[DD];
    grid::get_stride<Axis::X>(stride, dims);
    return symphas::internal::mixed_stencil_apply_type<
        OA, std::index_sequence<OD1, OD2, OD3>>{}(v, stride, divh);
  }

  //! Determine the finite difference approximation to the mixed derivative.
  /*!
   * Determine the finite difference approximation to the mixed derivative of
   * the given order.
   *
   * \tparam OD Order of the derivative.
   * \tparam ax Axis of the directional derivative.
   */
  template<size_t OD1, size_t OD2, typename T>
	__host__ __device__ auto apply_mixed(T* const v) const
	{
    len_type stride[DD];
    grid::get_stride<Axis::X>(stride, dims);
    return symphas::internal::mixed_stencil_apply_type<
        OA, std::index_sequence<OD1, OD2>>{}(v, stride, divh);
  }

  //! Determine the finite difference approximation to the mixed derivative.
  /*!
   * Determine the finite difference approximation to the mixed derivative of
   * the given order.
   *
   * \tparam OD Order of the derivative.
   * \tparam ax Axis of the directional derivative.
   */
  template<size_t OD1, typename T>
	__host__ __device__ auto apply_mixed(T* const v) const
	{
    len_type stride[DD];
    grid::get_stride<Axis::X>(stride, dims);
    static symphas::internal::GeneratedStencilApply stencil(
        expr::get_central_space_stencil<OD1, OA, 1>());

    return stencil(v, stride, divh);
    // return symphas::internal::GeneratedStencilApply<Stt>{ stride, divh }(v);
  }

  //! Determine the finite difference approximation to the mixed derivative.
  /*!
   * Determine the finite difference approximation to the mixed derivative of
   * the given order.
   *
   * \tparam OD Order of the derivative.
   * \tparam ax Axis of the directional derivative.
   */
  template <size_t OD1, size_t OD2, size_t OD3, typename T>
  __host__ __device__ auto apply_mixed(
      T *const v, len_type const (&override_dims)[DD]) const {
    len_type stride[DD];
    grid::get_stride<Axis::X>(stride, override_dims);
    return symphas::internal::mixed_stencil_apply_type<
        OA, std::index_sequence<OD1, OD2, OD3>>{}(v, stride, divh);
  }

  //! Determine the finite difference approximation to the mixed derivative.
  /*!
   * Determine the finite difference approximation to the mixed derivative of
   * the given order.
   *
   * \tparam OD Order of the derivative.
   * \tparam ax Axis of the directional derivative.
   */
  template <size_t OD1, size_t OD2, typename T>
  __host__ __device__ auto apply_mixed(
      T *const v, len_type const (&override_dims)[DD]) const {
    len_type stride[DD];
    grid::get_stride<Axis::X>(stride, override_dims);
    return symphas::internal::mixed_stencil_apply_type<
        OA, std::index_sequence<OD1, OD2>>{}(v, stride, divh);
  }

  //! Determine the finite difference approximation to the mixed derivative.
  /*!
   * Determine the finite difference approximation to the mixed derivative of
   * the given order.
   *
   * \tparam OD Order of the derivative.
   * \tparam ax Axis of the directional derivative.
   */
  template <size_t OD1, typename T>
  __host__ __device__ auto apply_mixed(
      T *const v, len_type const (&override_dims)[DD]) const {
    len_type stride[DD];
    grid::get_stride<Axis::X>(stride, override_dims);
    static symphas::internal::GeneratedStencilApply stencil(
        expr::get_central_space_stencil<OD1, OA, 1>());

    return stencil(v, stride, divh);
    // return symphas::internal::GeneratedStencilApply<Stt>{ stride, divh }(v);
  }

  //! Laplacian (2nd order derivative) of the field.
  template <typename T>
  __host__ __device__ inline auto laplacian(
      T *const v, const len_type (&stride)[DD]) const {
    return apply<2>(v, stride);
  }

  //! Bilaplacian (4th order derivative) of the field.
  template <typename T>
  __host__ __device__ inline auto bilaplacian(
      T *const v, const len_type (&stride)[DD]) const {
    return apply<4>(v, stride);
  }

  //! Gradlaplacian (gradient of the laplacian) of the field.
  template <typename T>
  __host__ __device__ inline auto gradlaplacian(
      T *const v, const len_type (&stride)[DD]) const {
    return apply<3>(v, stride);
  }

  template <typename T>
  __host__ __device__ auto gradient(T *const v,
                                    const len_type (&stride)[DD]) const {
    return apply<1>(v, stride);
  }
};

//! A specialization of the generalized derivative which returns the original
//! value.
/*!
 * This specialization is meant to match the return value type for the
 * derivative of the specified order. It does not actually compute the
 * derivative and is primarily used for testing by substituting it where a
 * functional stencil is intended.
 */
template <>
struct GeneralizedStencil<0, DEFAULT_STENCIL_ACCURACY> {
  GeneralizedStencil(...) {}

  template <size_t OD, typename T>
  auto apply(T *const v) const {
    using namespace std;
    using namespace symphas::math;
    if constexpr (OD % 2 == 1) {
      return *v;
    } else {
      return abs(*v);
    }
  }
};

// ******************************************************************************************

//! @}

// \cond

/* translations
 */

#define __DZ stride[2]
#define __DY stride[1]
#define __DX stride[0]

#define __2DZ (2 * __DZ)
#define __2DY (2 * __DY)
#define __2DX (2 * __DX)

#define __3DZ (3 * __DZ)
#define __3DY (3 * __DY)
#define __3DX (3 * __DX)

// primary index
#define v0 v[0]

/* single combinations
 */

#define vx v[__DX]
#define vx_ v[-__DX]
#define vy v[__DY]
#define vy_ v[-__DY]
#define vz v[__DZ]
#define vz_ v[-__DZ]

#define vx2 v[__2DX]
#define vx2_ v[-__2DX]
#define vy2 v[__2DY]
#define vy2_ v[-__2DY]
#define vz2 v[__2DZ]
#define vz2_ v[-__2DZ]

#define vx3 v[__3DX]
#define vx3_ v[-__3DX]
#define vy3 v[__3DY]
#define vy3_ v[-__3DY]
#define vz3 v[__3DZ]
#define vz3_ v[-__3DZ]

/*
 * all xy combinations
 */

#define vxy v[__DX + __DY]
#define vxy_ v[__DX - __DY]
#define vx_y v[-__DX + __DY]
#define vx_y_ v[-__DX - __DY]

#define vx2y v[__2DX + __DY]
#define vx2y_ v[__2DX - __DY]
#define vx2_y v[-__2DX + __DY]
#define vx2_y_ v[-__2DX - __DY]

#define vxy2 v[__DX + __2DY]
#define vxy2_ v[__DX - __2DY]
#define vx_y2 v[-__DX + __2DY]
#define vx_y2_ v[-__DX - __2DY]

#define vx3y v[__3DX + __DY]
#define vx3y_ v[__3DX - __DY]
#define vx3_y v[-__3DX + __DY]
#define vx3_y_ v[-__3DX - __DY]

#define vxy3 v[__DX + __3DY]
#define vxy3_ v[__DX - __3DY]
#define vx_y3 v[-__DX + __3DY]
#define vx_y3_ v[-__DX - __3DY]

#define vx2y2 v[__2DX + __2DY]
#define vx2y2_ v[__2DX - __2DY]
#define vx2_y2 v[-__2DX + __2DY]
#define vx2_y2_ v[-__2DX - __2DY]

#define vx3y3 v[__3DX + __3DY]
#define vx3y3_ v[__3DX - __3DY]
#define vx3_y3 v[-__3DX + __3DY]
#define vx3_y3_ v[-__3DX - __3DY]

#define vx3y2 v[__3DX + __2DY]
#define vx3y2_ v[__3DX - __2DY]
#define vx3_y2 v[-__3DX + __2DY]
#define vx3_y2_ v[-__3DX - __2DY]

#define vx2y3 v[__2DX + __3DY]
#define vx2y3_ v[__2DX - __3DY]
#define vx2_y3 v[-__2DX + __3DY]
#define vx2_y3_ v[-__2DX - __3DY]

/*
 * all xz combinations
 */

#define vxz v[__DX + __DZ]
#define vxz_ v[__DX - __DZ]
#define vx_z v[-__DX + __DZ]
#define vx_z_ v[-__DX - __DZ]

#define vx2z v[__2DX + __DZ]
#define vx2z_ v[__2DX - __DZ]
#define vx2_z v[-__2DX + __DZ]
#define vx2_z_ v[-__2DX - __DZ]

#define vxz2 v[__DX + __2DZ]
#define vxz2_ v[__DX - __2DZ]
#define vx_z2 v[-__DX + __2DZ]
#define vx_z2_ v[-__DX - __2DZ]

#define vx2z2 v[__2DX + __2DZ]
#define vx2z2_ v[__2DX - __2DZ]
#define vx2_z2 v[-__2DX + __2DZ]
#define vx2_z2_ v[-__2DX - __2DZ]

/*
 * all yz combinations
 */

#define vyz v[__DY + __DZ]
#define vyz_ v[__DY - __DZ]
#define vy_z v[-__DY + __DZ]
#define vy_z_ v[-__DY - __DZ]

#define vyz2 v[__DY + __2DZ]
#define vyz2_ v[__DY - __2DZ]
#define vy_z2 v[-__DY + __2DZ]
#define vy_z2_ v[-__DY - __2DZ]

#define vy2z v[__2DY + __DZ]
#define vy2z_ v[__2DY - __DZ]
#define vy2_z v[-__2DY + __DZ]
#define vy2_z_ v[-__2DY - __DZ]

#define vy2z2 v[__2DY + __2DZ]
#define vy2z2_ v[__2DY - __2DZ]
#define vy2_z2 v[-__2DY + __2DZ]
#define vy2_z2_ v[-__2DY - __2DZ]

/*
 * xyz combinations
 */

#define vxyz v[__DX + __DY + __DZ]
#define vxyz_ v[__DX + __DY - __DZ]
#define vxy_z v[__DX - __DY + __DZ]
#define vxy_z_ v[__DX - __DY - __DZ]
#define vx_yz v[-__DX + __DY + __DZ]
#define vx_yz_ v[-__DX + __DY - __DZ]
#define vx_y_z v[-__DX - __DY + __DZ]
#define vx_y_z_ v[-__DX - __DY - __DZ]

#define vx2yz v[__2DX + __DY + __DZ]
#define vx2yz_ v[__2DX + __DY - __DZ]
#define vx2y_z v[__2DX - __DY + __DZ]
#define vx2y_z_ v[__2DX - __DY - __DZ]
#define vx2_yz v[-__2DX + __DY + __DZ]
#define vx2_yz_ v[-__2DX + __DY - __DZ]
#define vx2_y_z v[-__2DX - __DY + __DZ]
#define vx2_y_z_ v[-__2DX - __DY - __DZ]

#define vxy2z v[__DX + __2DY + __DZ]
#define vxy2z_ v[__DX + __2DY - __DZ]
#define vxy2_z v[__DX - __2DY + __DZ]
#define vxy2_z_ v[__DX - __2DY - __DZ]
#define vx_y2z v[-__DX + __2DY + __DZ]
#define vx_y2z_ v[-__DX + __2DY - __DZ]
#define vx_y2_z v[-__DX - __2DY + __DZ]
#define vx_y2_z_ v[-__DX - __2DY - __DZ]

#define vxyz2 v[__DX + __DY + __2DZ]
#define vxyz2_ v[__DX + __DY - __2DZ]
#define vxy_z2 v[__DX - __DY + __2DZ]
#define vxy_z2_ v[__DX - __DY - __2DZ]
#define vx_yz2 v[-__DX + __DY + __2DZ]
#define vx_yz2_ v[-__DX + __DY - __2DZ]
#define vx_y_z2 v[-__DX - __DY + __2DZ]
#define vx_y_z2_ v[-__DX - __DY - __2DZ]

#define vx2y2z2 v[__2DX + __2DY + __2DZ]
#define vx2y2z2_ v[__2DX + __2DY - __2DZ]
#define vx2y2_z2 v[__2DX - __2DY + __2DZ]
#define vx2y2_z2_ v[__2DX - __2DY - __2DZ]
#define vx2_y2z2 v[-__2DX + __2DY + __2DZ]
#define vx2_y2z2_ v[-__2DX + __2DY - __2DZ]
#define vx2_y2_z2 v[-__2DX - __2DY + __2DZ]
#define vx2_y2_z2_ v[-__2DX - __2DY - __2DZ]

// \endcond
