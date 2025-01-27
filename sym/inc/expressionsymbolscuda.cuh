
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
 * MODULE:  sym
 * PURPOSE: Defines all the expression types used in the symbolic
 * algebra for cuda functionality.
 *
 * ***************************************************************************
 */

#include "expressiontypeincludes.h"

#ifdef USING_CUDA

#include <cuda_runtime.h>

namespace symphas::internal {

template <size_t D>
__host__ __device__ iter_type get_index(const len_type (&dims)[D],
                                        const len_type (&pos)[D]) {
  iter_type n = 0;
  len_type stride = 1;
  for (iter_type i = 0; i < D; ++i) {
    n += stride * pos[i];
    stride *= dims[i];
  }
  return n;
}

}  // namespace symphas::internal

template <size_t D>
struct position_type {
  using arr_t = iter_type[D];

  len_type dims[D];
  iter_type pos[D];
  iter_type index;
  __host__ __device__ position_type() : dims{}, pos{}, index{} {}
  template <size_t... Is>
  __host__ __device__ position_type(const len_type (&dims)[D],
                                    const unsigned (&pos)[D],
                                    std::index_sequence<Is...>)
      : position_type(dims, arr_t{static_cast<iter_type>(pos[Is])...}) {}
  __host__ __device__ position_type(const len_type (&dims)[D],
                                    const unsigned (&pos)[D])
      : position_type(dims, pos, std::make_index_sequence<D>{}) {}

  template <size_t... Is>
  __host__ __device__ position_type(const len_type (&dims)[D],
                                    const iter_type (&pos)[D],
                                    std::index_sequence<Is...>)
      : dims{dims[Is]...},
        pos{pos[Is]...},
        index{symphas::internal::get_index<D>(dims, pos)} {}
  __host__ __device__ position_type(const len_type (&dims)[D],
                                    const iter_type (&pos)[D])
      : position_type(dims, pos, std::make_index_sequence<D>{}) {}

  __host__ __device__ position_type(position_type<D> const &other)
      : position_type(other.dims, other.pos, std::make_index_sequence<D>{}) {}
  __host__ __device__ position_type(position_type<D> &&other)
      : position_type(other.dims, other.pos, std::make_index_sequence<D>{}) {}
};

template <typename T, expr::exp_key_t X>
struct CuTerm : CuEvaluable<CuTerm<T, X>> {
  T *values;

  CuTerm(T *data) : values{data} {}

  __host__ __device__ T result(T const &value) const {
    if constexpr (expr::_Xk<X> == 1) {
      return value;
    } else {
      using symphas::math::pow;

      if constexpr (expr::_Xk_t<X>::D == 1) {
        auto result = pow<size_t(expr::_Xk_t<X>::N)>(value);
        if constexpr (expr::_Xk_t<X>::sign) {
          return 1.0 / result;
        } else {
          return result;
        }
      } else {
        return pow(value, expr::_Xk<X>);
      }
    }
  }

  template <size_t D>
  __host__ __device__ T eval(position_type<D> const &n) const {
    return result(values[n.index]);
  }
  __host__ __device__ T eval(iter_type n) const { return result(values[n]); }
};

template <typename T, size_t D, expr::exp_key_t X>
struct CuTerm<any_vector_t<T, D>, X>
    : CuEvaluable<CuTerm<any_vector_t<T, D>, X>> {
  T *values[D];

  CuTerm(T *const (&data)[D]) : values{} {
    for (iter_type i = 0; i < D; ++i) {
      values[i] = data[i];
    }
  }

  __host__ __device__ any_vector_t<T, D> result(
      any_vector_t<T, D> const &value) const {
    if constexpr (expr::_Xk<X> == 1) {
      return value;
    } else {
      using symphas::math::pow;

      if constexpr (expr::_Xk_t<X>::D == 1) {
        auto result = pow<size_t(expr::_Xk_t<X>::N)>(value);
        if constexpr (expr::_Xk_t<X>::sign) {
          return 1.0 / result;
        } else {
          return result;
        }
      } else {
        return pow(value, expr::_Xk<X>);
      }
    }
  }

  __host__ __device__ any_vector_t<T, D> eval(position_type<D> const &n) const {
    any_vector_t<T, D> value;
    for (iter_type i = 0; i < D; ++i) {
      value[i] = values[i][n.index];
    }
    return result(value);
  }
  __host__ __device__ any_vector_t<T, D> eval(iter_type n) const {
    any_vector_t<T, D> value;
    for (iter_type i = 0; i < D; ++i) {
      value[i] = values[i][n];
    }
    return result(value);
  }
};

template <size_t D, typename T, expr::exp_key_t X>
struct CuTermRegion : CuEvaluable<CuTermRegion<D, T, X>> {
  CuTerm<T, X> term;
  grid::select_region_cuda<D> region;
  len_type dims[D];
  T empty;

  template <size_t... Is>
  CuTermRegion(T *data, grid::select_region_cuda<D> const &region,
               const len_type (&dims)[D], T const &empty,
               std::index_sequence<Is...>)
      : term{data}, region{region}, dims{dims[Is]...}, empty{empty} {}
  CuTermRegion(T *data, grid::select_region_cuda<D> const &region,
               const len_type (&dims)[D], T const &empty)
      : CuTermRegion(data, region, dims, empty, std::make_index_sequence<D>{}) {
  }

  __host__ __device__ T eval(position_type<D> const &n) const {
    return term.result(region(n.pos, term.values, dims, empty).dev());
  }

  __host__ __device__ T eval(iter_type n) const {
    iter_type pos[D]{};
    grid::get_grid_position(pos, dims, n);
    return term.result(region(pos, term.values, dims, empty).dev());
  }
};

template <size_t D, typename T, expr::exp_key_t X>
struct CuTermRegion<D, any_vector_t<T, D>, X>
    : CuEvaluable<CuTermRegion<D, any_vector_t<T, D>, X>> {
  CuTerm<any_vector_t<T, D>, X> term;
  grid::select_region_cuda<D> region;
  len_type dims[D];
  T empty[D];

  template <size_t... Is>
  CuTermRegion(T *const (&data)[D], grid::select_region_cuda<D> const &region,
               const len_type (&dims)[D], T const (&empty)[D],
               std::index_sequence<Is...>)
      : term{data}, region{region}, dims{dims[Is]...}, empty{empty[Is]...} {}
  CuTermRegion(T *const (&data)[D], grid::select_region_cuda<D> const &region,
               const len_type (&dims)[D], T const (&empty)[D])
      : CuTermRegion(data, region, dims, empty, std::make_index_sequence<D>{}) {
  }

  __host__ __device__ any_vector_t<T, D> eval(position_type<D> const &n) const {
    return term.result(region(n.pos, term.values, dims, empty).dev());
  }

  __host__ __device__ any_vector_t<T, D> eval(iter_type n) const {
    iter_type pos[D]{};
    grid::get_grid_position(pos, dims, n);
    return term.result(region(pos, term.values, dims, empty).dev());
  }
};

template <size_t D, typename T, expr::exp_key_t X>
struct CuTermRegionHost : CuEvaluable<CuTermRegionHost<D, T, X>> {
  CuTerm<T, X> term;
  grid::select_region_cuda<D> region;
  len_type dims[D];
  T empty;

  template <size_t... Is>
  CuTermRegionHost(T *data, grid::select_region_cuda<D> const &region,
                   const len_type (&dims)[D], T const &empty,
                   std::index_sequence<Is...>)
      : term{data}, region{region}, dims{dims[Is]...}, empty{empty} {}
  CuTermRegionHost(T *data, grid::select_region_cuda<D> const &region,
                   const len_type (&dims)[D], T const &empty)
      : CuTermRegionHost(data, region, dims, empty,
                         std::make_index_sequence<D>{}) {}

  __host__ __device__ T eval(position_type<D> const &n) const {
    return T{};
    // return term.result(region(n.pos, term.values, dims, empty));
  }

  __host__ __device__ T eval(iter_type n) const {
    return T{};
    // iter_type pos[D]{};
    // grid::get_grid_position(pos, dims, n);
    // return term.result(region(pos, term.values, dims, empty).dev());
  }
};

template <size_t D, typename T, expr::exp_key_t X>
struct CuTermRegionHost<D, any_vector_t<T, D>, X>
    : CuEvaluable<CuTermRegionHost<D, any_vector_t<T, D>, X>> {
  CuTerm<any_vector_t<T, D>, X> term;
  grid::select_region_cuda<D> region;
  len_type dims[D];
  T empty[D];

  template <size_t... Is>
  CuTermRegionHost(T *data, grid::select_region_cuda<D> const &region,
                   const len_type (&dims)[D], T const (&empty)[D],
                   std::index_sequence<Is...>)
      : term{data}, region{region}, dims{dims[Is]...}, empty{empty[Is]...} {}
  CuTermRegionHost(T *data, grid::select_region_cuda<D> const &region,
                   const len_type (&dims)[D], T const (&empty)[D])
      : CuTermRegionHost(data, region, dims, empty,
                         std::make_index_sequence<D>{}) {}

  __host__ __device__ any_vector_t<T, D> eval(position_type<D> const &n) const {
    return T{};
    // return term.result(region(n.pos, term.values, dims, empty).dev());
  }

  __host__ __device__ any_vector_t<T, D> eval(iter_type n) const {
    return T{};
    /*iter_type pos[D]{};
    grid::get_grid_position(pos, dims, n);
    return term.result(region(pos, term.values, dims, empty).dev());*/
  }
};

template <typename... Es>
struct CuTermList;

template <typename E1, typename E2, typename... Es>
struct CuTermList<E1, E2, Es...> : CuTermList<E2, Es...> {
  E1 e;
  CuTermList(E1 const &e1, E2 const &e2, Es const &...es)
      : CuTermList<E2, Es...>{e2, es...}, e{e1} {}
  template <size_t I>
  __device__ __host__ const symphas::lib::type_at_index<I, E1, E2, Es...> &get()
      const {
    if constexpr (I == 0) {
      return e;
    } else {
      return CuTermList<E2, Es...>::template get<I - 1>();
    }
  }
};

template <typename E1>
struct CuTermList<E1> {
  E1 e;

  CuTermList(E1 const &e) : e{e} {}

  template <size_t I>
  __device__ __host__ const E1 &get() const {
    if constexpr (I == 0) {
      return e;
    }
  }
};
//
// template <typename T, size_t... Ns>
// struct CuTensor {
//  explicit CuTensor(T const& entry) : value{entry} {}
//
//  __host__ __device__ auto eval(iter_type n = 0) const {
//    return symphas::internal::tensor_as_coeff<Ns...>(value);
//  }
//
//  T value;
//};

struct CuIdentity : CuEvaluable<CuIdentity> {
  template <typename index_type>
  __host__ __device__ scalar_t eval(index_type &&) const {
    return 1;
  }
};

struct CuNegIdentity : CuIdentity {
  template <typename index_type>
  __host__ __device__ scalar_t eval(index_type &&) const {
    return -1;
  }
};

struct CuVoid : CuEvaluable<CuVoid> {
  template <typename index_type>
  __host__ __device__ CuVoid eval(index_type &&) const {
    return {};
  }
  __host__ __device__ operator scalar_t() const { return 0; }
  template <size_t D>
  __host__ __device__ operator any_vector_t<scalar_t, D>() const {
    return {};
  }
  __device__ __host__ CuVoid operator+(CuVoid) { return {}; }
  __device__ __host__ CuVoid operator-(CuVoid) { return {}; }
  __device__ __host__ CuVoid operator*(CuVoid) { return {}; }
};

template <typename T>
struct CuLiteral : CuEvaluable<CuLiteral<T>> {
  T value;
  CuLiteral(T const &value) : value{value} {}
  template <typename index_type>
  __host__ __device__ T eval(index_type &&) const {
    return value;
  }
};

// template <typename T, typename I = void>
// struct CuCoeff;
//
// template <typename T>
// struct CuCoeff<T, void>;
//
// template <typename T, typename I>
// struct CuCoeff : CuEvaluable<CuCoeff<T, I>> {
//   CuCoeff(symphas::internal::coeff_data<T> const &data) : data{data} {}
//
//   template <typename index_type>
//   __host__ __device__ auto eval(index_type &&) const {
//     if (data.len > 0) {
//       return data[0];
//     } else {
//       return T{};
//     }
//   }
//
//   symphas::internal::coeff_data_cuda<T> data;
// };
//
// template <typename T>
// struct CuCoeff<T, DynamicIndex> : CuEvaluable<CuCoeff<T, DynamicIndex>> {
//   CuCoeff(DynamicIndex const &index,
//           symphas::internal::coeff_data<T> const &data)
//       : data{data}, index{index} {}
//
//   template <typename index_type>
//   __host__ __device__ auto eval(index_type &&) const {
//     if ((int)index < data.len) {
//       return data[(int)index];
//     } else {
//       return T{};
//     }
//   }
//
//   symphas::internal::coeff_data<T> data;
//   DynamicIndex index;
// };
//
// template <typename I>
// struct CuCoeff<void, I>;
//
// template <>
// struct CuCoeff<void, DynamicIndex>;

template <size_t N, size_t D>
struct CuFractionLiteral : CuEvaluable<CuFractionLiteral<N, D>> {
  template <typename index_type>
  __host__ __device__ scalar_t eval(index_type &&) const {
    return scalar_t(N) / scalar_t(D);
  }
};

template <size_t N, size_t D>
struct CuNegFractionLiteral : CuEvaluable<CuNegFractionLiteral<N, D>> {
  template <typename index_type>
  __host__ __device__ scalar_t eval(index_type &&n) const {
    return -CuFractionLiteral<N, D>().eval(std::forward<index_type>(n));
  }
};

__device__ void print_value(const char *str, scalar_t value, iter_type i,
                            iter_type n) {
  if (value > 100 || value < -1) {
    printf("%s term %d of %d > wrong value\n", str, i, n);
  } else {
    printf("%s term %d of %d > %lf\n", str, i, n, value);
  }
}

__device__ void print_value(const char *str, any_vector_t<scalar_t, 2> value,
                            iter_type i, iter_type n) {
  if (value[0] > 100 || value[0] < -1 || value[1] > 100 || value[1] < -1) {
    printf("%s term %d of %d > wrong value\n", str, i, n);
  } else {
    printf("%s term %d of %d > {%lf,%lf}\n", str, i, n, value[0], value[1]);
  }
}

template <typename... Es, typename index_type>
__device__ void print_value_next(const char *str,
                                 CuTermList<Es...> const &terms, index_type &&i,
                                 iter_type n, std::index_sequence<>) {}

template <typename... Es, typename index_type, size_t I0, size_t... Is>
__device__ void print_value_next(const char *str,
                                 CuTermList<Es...> const &terms, index_type &&i,
                                 iter_type n, std::index_sequence<I0, Is...>) {
  print_value(str, terms.template get<I0>().eval(std::forward<index_type>(i)),
              int(I0 + 1), int(sizeof...(Es)));
  print_value_next(str, terms, std::forward<index_type>(i), n,
                   std::index_sequence<Is...>{});
}

template <auto f, typename E>
struct CuFunctionApply : CuEvaluable<CuFunctionApply<f, E>> {
  E e;
  CuFunctionApply(E const &e) : e{e} {}
  template <typename index_type>
  __host__ __device__ scalar_t eval(index_type &&n) const {
    return f(e.eval(std::forward<index_type>(n)));
  }
};

template <typename... Es>
struct CuAdd : CuEvaluable<CuAdd<Es...>> {
  CuTermList<Es...> terms;
  CuAdd(Es const &...es) : terms{es...} {}

  template <typename index_type, size_t... Is>
  __host__ __device__ auto eval(index_type &&i,
                                std::index_sequence<Is...>) const {
    // print_value_next("add", terms, std::forward<index_type>(i),
    //                  int(sizeof...(Is)), std::index_sequence<Is...>{});
    return (terms.template get<Is>().eval(std::forward<index_type>(i)) + ...);
  }

  template <typename index_type>
  __host__ __device__ auto eval(index_type &&i) const {
    return this->eval(std::forward<index_type>(i),
                      std::make_index_sequence<sizeof...(Es)>{});
  }
};

template <typename... Es>
struct CuMul : CuEvaluable<CuMul<Es...>> {
  CuTermList<Es...> terms;
  CuMul(Es const &...es) : terms{es...} {}
  template <typename index_type, size_t... Is>
  __host__ __device__ auto eval(index_type &&i,
                                std::index_sequence<Is...>) const {
    // print_value_next("mul", terms, std::forward<index_type>(i),
    //                  int(sizeof...(Is)), std::index_sequence<Is...>{});
    return (terms.template get<Is>().eval(std::forward<index_type>(i)) * ...);
  }

  template <typename index_type>
  __host__ __device__ auto eval(index_type &&i) const {
    return this->eval(std::forward<index_type>(i),
                      std::make_index_sequence<sizeof...(Es)>{});
  }
};

template <typename Dd, typename T, size_t D, typename Sp>
struct CuDerivative : CuEvaluable<CuDerivative<Dd, T, D, Sp>> {
  T *values;
  Sp *solver;

  Sp *copySolverToDevice(Sp const &hostSolver) const {
    Sp *deviceSolver;
    CHECK_CUDA_ERROR(cudaMalloc(&deviceSolver, sizeof(Sp)));
    CHECK_CUDA_ERROR(cudaMemcpy(deviceSolver, &hostSolver, sizeof(Sp),
                                cudaMemcpyHostToDevice));
    return deviceSolver;
  }

  Sp *copySolverFromDevice(Sp const *otherDeviceSolver) const {
    Sp *deviceSolver;
    CHECK_CUDA_ERROR(cudaMalloc(&deviceSolver, sizeof(Sp)));
    CHECK_CUDA_ERROR(cudaMemcpy(deviceSolver, otherDeviceSolver, sizeof(Sp),
                                cudaMemcpyDeviceToDevice));
    return deviceSolver;
  }

  CuDerivative() : values{nullptr}, solver{nullptr} {}

  CuDerivative(T *values, Sp const &hostSolver)
      : values{values}, solver{copySolverToDevice(hostSolver)} {}

  CuDerivative(T *values, Sp *solver) : values{values}, solver{solver} {}

  CuDerivative(CuDerivative<Dd, T, D, Sp> const &other)
      : values{other.values}, solver{copySolverFromDevice(other.solver)} {}

  CuDerivative(CuDerivative<Dd, T, D, Sp> &&other) : CuDerivative() {
    swap(*this, other);
  }

  friend void swap(CuDerivative &first, CuDerivative &second) {
    using std::swap;
    swap(first.values, second.values);
    swap(first.solver, second.solver);
  }

  template <typename value_type>
  __host__ __device__ auto eval(value_type &&value, iter_type i) const {
    return Dd{}(*solver, std::forward<value_type>(value), i);
  }

  __host__ __device__ auto eval(iter_type i) const {
    return eval(GridDataCUDA<T, D>(values), i);
  }

  template <size_t D>
  __host__ __device__ auto eval(position_type<D> const &n) const {
    return eval(n.index);
  }

  ~CuDerivative() { CHECK_CUDA_ERROR(cudaFree(solver)); }
};

template <typename Dd, typename T, size_t D, typename Sp>
struct CuDerivativeRegion : CuEvaluable<CuDerivativeRegion<Dd, T, D, Sp>> {
  CuDerivative<Dd, T, D, Sp> deriv;
  grid::select_region_cuda<D> region;
  len_type dims[D];
  T empty;

  CuDerivativeRegion() : deriv{}, region{}, dims{}, empty{} {}

  template <size_t... Is>
  CuDerivativeRegion(T *values, grid::select_region_cuda<D> const &region,
                     const len_type (&dims)[D], T const &empty,
                     Sp const &hostSolver, std::index_sequence<Is...>)
      : deriv{values, hostSolver},
        region{region},
        dims{dims[Is]...},
        empty{empty} {}
  template <size_t... Is>
  CuDerivativeRegion(T *values, grid::select_region_cuda<D> const &region,
                     const len_type (&dims)[D], T const &empty,
                     Sp const &hostSolver)
      : CuDerivativeRegion(values, region, dims, empty, hostSolver,
                           std::make_index_sequence<D>{}) {}

  template <size_t... Is>
  CuDerivativeRegion(T *values, grid::select_region_cuda<D> const &region,
                     len_type (&dims)[D], T const &empty, Sp *solver,
                     std::index_sequence<Is...>)
      : deriv{values, hostSolver},
        region{region},
        dims{dims[Is]...},
        empty{empty} {}

  CuDerivativeRegion(T *values, grid::select_region_cuda<D> const &region,
                     len_type (&dims)[D], T const &empty, Sp *solver)
      : CuDerivativeRegion(values, region, dims, empty, solver,
                           std::make_index_sequence<D>{}) {}

  template <size_t... Is>
  CuDerivativeRegion(CuDerivativeRegion<Dd, T, D, Sp> const &other,
                     std::index_sequence<Is...>)
      : deriv{other.deriv.values,
              other.deriv.copySolverFromDevice(other.deriv.solver)},
        region{other.region},
        dims{other.dims[Is]...},
        empty{other.empty} {}
  CuDerivativeRegion(CuDerivativeRegion<Dd, T, D, Sp> const &other)
      : CuDerivativeRegion(other, std::make_index_sequence<D>{}) {}

  CuDerivativeRegion(CuDerivativeRegion<Dd, T, D, Sp> &&other)
      : CuDerivativeRegion() {
    swap(*this, other);
  }

  friend void swap(CuDerivativeRegion &first, CuDerivativeRegion &second) {
    using std::swap;
    swap(first.deriv, second.deriv);
    swap(first.region, second.region);
    swap(first.dims, second.dims);
    swap(first.empty, second.empty);
  }

  __host__ __device__ auto eval(iter_type n) const {
    iter_type pos[D]{};
    grid::get_grid_position(pos, dims, n);
    return eval(position_type<D>(dims, pos));
  }

  template <size_t D>
  __host__ __device__ auto eval(position_type<D> const &n) const {
    auto value = region(n.pos, deriv.values, dims, empty);
    if (value.is_valid()) {
      auto r = deriv.eval(
          RegionGridDataCUDA<T, D>(value.value, region.dims, region.stride), 0);
      // print_value("deriv", r, n.index, grid::length<D>(dims));
      return r;
    }
    // print_value("deriv", value.dev(), n.index, grid::length<D>(dims));
    return value.dev();
  }
};

template <typename Dd, typename T, size_t D, typename Sp>
struct CuDerivativeHost : CuEvaluable<CuDerivativeHost<Dd, T, D, Sp>> {
  CuDerivative<Dd, T, D, Sp> derivative;
  len_type len;

  T *copyValuesToDevice(const T *hostValues, len_type len) const {
    T *deviceValues;
    CHECK_CUDA_ERROR(cudaMalloc(&deviceValues, len * sizeof(T)));
    CHECK_CUDA_ERROR(cudaMemcpy(deviceValues, &hostValues, len * sizeof(T),
                                cudaMemcpyHostToDevice));
    return deviceValues;
  }

  T *copyValuesFromDevice(const T *otherDeviceValues, len_type len) const {
    T *deviceValues;
    CHECK_CUDA_ERROR(cudaMalloc(&deviceValues, len * sizeof(T)));
    CHECK_CUDA_ERROR(cudaMemcpy(deviceValues, otherDeviceValues,
                                len * sizeof(T), cudaMemcpyDeviceToDevice));
    return deviceValues;
  }

  CuDerivativeHost() : derivative{}, len{} {}

  CuDerivativeHost(T *hostValues, len_type len, Sp const &hostSolver)
      : derivative{copyValuesToDevice(hostValues, len), hostSolver}, len{len} {}

  CuDerivativeHost(CuDerivativeHost<Dd, T, D, Sp> const &other)
      : derivative{copyValuesFromDevice(other.derivative.values, other.len),
                   other.derivative.copySolverFromDevice(
                       other.derivative.solver)},
        len{other.len} {}

  CuDerivativeHost(CuDerivativeHost<Dd, T, D, Sp> &&other)
      : CuDerivativeHost() {
    swap(*this, other);
  }

  friend void swap(CuDerivativeHost &first, CuDerivativeHost &second) {
    using std::swap;
    swap(first.derivative, second.derivative);
    swap(first.len, second.len);
  }

  __host__ __device__ auto eval(iter_type i) const {
    return derivative.eval(i);
  }

  template <size_t D>
  __host__ __device__ auto eval(position_type<D> const &i) const {
    return eval(i.index);
  }

  ~CuDerivativeHost() { CHECK_CUDA_ERROR(cudaFree(derivative.values)); }
};

template <typename T, typename E>
struct CuSeries : CuEvaluable<CuSeries<T, E>> {
  CuSeries(E *persistent, len_type len) : persistent{}, len{len} {
    CHECK_CUDA_ERROR(cudaMalloc(&this->persistent, sizeof(E) * len));
    CHECK_CUDA_ERROR(cudaMemcpy(this->persistent, persistent, sizeof(E) * len,
                                cudaMemcpyHostToDevice));
  }
  CuSeries(E *const *persistent, len_type len) : persistent{}, len{len} {
    CHECK_CUDA_ERROR(cudaMalloc(&this->persistent, sizeof(E) * len));
    CHECK_CUDA_ERROR(cudaMemcpy(this->persistent, *persistent, sizeof(E) * len,
                                cudaMemcpyDeviceToDevice));
  }

  CuSeries() : persistent{}, len{} {}

  CuSeries(CuSeries<T, E> const &other)
      : CuSeries{&other.persistent, other.len} {}
  CuSeries(CuSeries<T, E> &&other) : CuSeries() { swap(*this, other); }

  friend void swap(CuSeries &first, CuSeries &second) {
    using std::swap;
    swap(first.persistent, second.persistent);
    swap(first.len, second.len);
  }

  template <typename index_type>
  __host__ __device__ auto eval(index_type &&n) const {
    using namespace symphas::internal;

    if (len > 0) {
      auto result = persistent[0].eval(std::forward<index_type>(n));

      for (iter_type i = 1; i < len; ++i) {
        result = result + persistent[i].eval(std::forward<index_type>(n));
      }
      return result;
    } else {
      return T{};
    }
    // return T{};
  }

  ~CuSeries() { CHECK_CUDA_ERROR(cudaFree(persistent)); }

  E *persistent;
  len_type len;
};

//
// template <typename Dd, typename T, size_t D, typename Sp>
// struct CuDerivative : CuEvaluable<CuDerivative<Dd, T, D, Sp>> {
//  GridData<T, D> data;
//  T* values;
//  Sp solver;
//  bool clear_data;
//
//  T* copy_data_to_cuda(T* srcHost, len_type const* dims) {
//    T* destDevice;
//    len_type size = grid::length<D>(dims);
//    CHECK_CUDA_ERROR(cudaMalloc(&destDevice, sizeof(T) * size));
//    CHECK_CUDA_ERROR(cudaMemcpy(destDevice, srcHost, sizeof(T) * size,
//                                cudaMemcpyHostToDevice));
//    return destDevice;
//  }
//
//  T* cast_to_non_const(T const* src) { return const_cast<T*>(src); }
//
//  CuDerivative()
//      : data{nullptr, nullptr}, values{nullptr}, solver{}, clear_data{false}
//      {}
//  CuDerivative(CuDerivative<Dd, T, D, Sp> const& other)
//      : data{},
//        solver{other.solver},
//        values{other.values},
//        clear_data{other.clear_data} {
//    if (other.clear_data) {
//      T* values0;
//      len_type len = grid::length<D>(other.data.dims);
//      symphas::cuda::allocate(&values0, sizeof(T) * len);
//      CHECK_CUDA_ERROR(cudaMemcpy(values0, &other.data[0], sizeof(T) * len,
//                            cudaMemcpyDeviceToDevice);
//      data = GridData<T, D>(values0, other.data.dims);
//    } else {
//      data = other.data;
//    }
//  }
//  CuDerivative(CuDerivative<Dd, T, D, Sp>&& other) : CuDerivative() {
//    swap(*this, other);
//  }
//  friend void swap(CuDerivative& first, CuDerivative& second) {
//    using std::swap;
//    swap(first.data, second.data);
//    swap(first.values, second.values);
//    swap(first.solver, second.solver);
//    swap(first.clear_data, second.clear_data);
//  }
//
//  CuDerivative(BlockCUDA<T> const& data, const len_type* dims, Sp const&
//  solver)
//      : data{cast_to_non_const(data.values), dims},
//        values{cast_to_non_const(data.values)},
//        solver{solver},
//        clear_data{false} {}
//  template <Axis ax>
//  CuDerivative(VectorComponentData<ax, T*, D> const& data, const len_type*
//  dims,
//               Sp const& solver)
//      : data{cast_to_non_const(data.values), dims},
//        values{cast_to_non_const(data.values)},
//        solver{solver},
//        clear_data{false} {}
//  CuDerivative(Grid<T, D> const& data, const len_type* dims, Sp const&
//  solver)
//      : data{copy_data_to_cuda(data.values, dims), dims},
//        values{&this->data[0]},
//        solver{solver},
//        clear_data{true} {}
//
//  __device__ auto eval(iter_type i) const {
//    printf("%d %p %p\n", i, &values[0], &values[i]);
//    return Dd{}(solver, data, i);
//  }
//
//  __device__ __host__ ~CuDerivative() {
//    if (clear_data) {
//      cudaFree(&data[0]);
//    }
//  }
//};

// CuDerivative<Dd, CuTerm<T, 1>, Sp>

namespace expr {

template <typename E1, typename E2>
auto operator*(CuEvaluable<E1> const &a, CuEvaluable<E2> const &b) {
  return CuMul{*static_cast<E1 const *>(&a), *static_cast<E2 const *>(&b)};
}

template <typename... As, typename E2, size_t... Is>
auto to_cu_mul(CuMul<As...> const &a, CuEvaluable<E2> const &b,
               std::index_sequence<Is...>) {
  return CuMul{a.terms.get<Is>()..., *static_cast<E2 const *>(&b)};
}

template <typename E1, typename... Bs, size_t... Is>
auto to_cu_mul(CuEvaluable<E1> const &a, CuMul<Bs...> const &b,
               std::index_sequence<Is...>) {
  return CuMul{*static_cast<E1 const *>(&a), b.terms.get<Is>()...};
}

template <typename... As, typename... Bs, size_t... Is, size_t... Js>
auto to_cu_mul(CuMul<As...> const &a, CuMul<Bs...> const &b,
               std::index_sequence<Is...>, std::index_sequence<Js...>) {
  return CuMul{a.terms.get<Is>()..., b.terms.get<Js>()...};
}

template <typename... As, typename E2>
auto operator*(CuMul<As...> const &a, CuEvaluable<E2> const &b) {
  return to_cu_mul(a, *static_cast<E2 const *>(&b),
                   std::make_index_sequence<sizeof...(As)>{});
}

template <typename E1, typename... Bs>
auto operator*(CuEvaluable<E1> const &a, CuMul<Bs...> const &b) {
  return to_cu_mul(*static_cast<E1 const *>(&a), b,
                   std::make_index_sequence<sizeof...(Bs)>{});
}
template <typename... As, typename... Bs>
auto operator*(CuMul<As...> const &a, CuMul<Bs...> const &b) {
  return to_cu_mul(a, b, std::make_index_sequence<sizeof...(As)>{},
                   std::make_index_sequence<sizeof...(Bs)>{});
}

template <expr::exp_key_t X, typename T>
auto to_cuda_expr(T *data) {
  return CuTerm<T, X>(data);
}

template <expr::exp_key_t X, typename T, size_t D>
auto to_cuda_expr(GridCUDA<T, D> const &data) {
  return CuTerm<T, X>(data.values);
}

template <expr::exp_key_t X, typename T, size_t D>
auto to_cuda_expr(BoundaryGridCUDA<T, D> const &data) {
  return CuTerm<T, X>(data.values);
}

template <expr::exp_key_t X, typename T, size_t D>
auto to_cuda_expr(RegionalGridCUDA<T, D> const &data) {
  return CuTermRegion<D, T, X>(data.values, data.region, data.dims, data.empty);
}

template <expr::exp_key_t X, typename T, size_t D>
auto to_cuda_expr(RegionalGrid<T, D> const &data) {
  return CuTermRegionHost<D, T, X>(data.values, data.region, data.dims,
                                   data.empty);
}

template <expr::exp_key_t X, Axis ax, typename T, size_t D>
auto to_cuda_expr(VectorComponentData<ax, T *, D> const &data) {
  return CuTerm<T, X>(data.values);
}

template <typename G, expr::exp_key_t X>
auto to_cuda_expr(Term<G, X> const &e) {
  return to_cuda_expr<X>(BaseData<G>::get(e));
}

template <typename A, typename B>
auto to_cuda_expr(OpBinaryMul<A, B> const &e) {
  return to_cuda_expr(e.a) * to_cuda_expr(e.b);
}

template <typename V, typename T0, size_t I0>
auto to_cuda_expr(OpTerms<V, T0> const &e, std::index_sequence<I0>) {
  return CuMul{to_cuda_expr(expr::get<0>(e)),
               to_cuda_expr(expr::get<I0 + 1>(e))};
}

template <typename T0, size_t I0>
auto to_cuda_expr(OpTerms<OpIdentity, T0> const &e, std::index_sequence<I0>) {
  return to_cuda_expr(expr::get<I0 + 1>(e));
}

template <typename V, typename... Ts, size_t I0, size_t I1, size_t... Is>
auto to_cuda_expr(OpTerms<V, Ts...> const &e,
                  std::index_sequence<I0, I1, Is...>) {
  return CuMul{to_cuda_expr(expr::get<0>(e)),
               to_cuda_expr(expr::get<I0 + 1>(e)),
               to_cuda_expr(expr::get<I1 + 1>(e)),
               to_cuda_expr(expr::get<Is + 1>(e))...};
}

template <typename... Ts, size_t I0, size_t I1, size_t... Is>
auto to_cuda_expr(OpTerms<OpIdentity, Ts...> const &e,
                  std::index_sequence<I0, I1, Is...>) {
  return CuMul{to_cuda_expr(expr::get<I0 + 1>(e)),
               to_cuda_expr(expr::get<I1 + 1>(e)),
               to_cuda_expr(expr::get<Is + 1>(e))...};
}

template <typename V, typename... Ts>
auto to_cuda_expr(OpTerms<V, Ts...> const &e) {
  return to_cuda_expr(e, std::make_index_sequence<sizeof...(Ts)>{});
}

template <typename Dd, typename T, size_t D, typename Sp>
auto to_cuda_expr_deriv(GridCUDA<T, D> const &data, Sp const &solver) {
  return CuDerivative<Dd, T, D, Sp>{data.values, solver};
}

template <typename Dd, typename T, size_t D, typename Sp>
auto to_cuda_expr_deriv(RegionalGridCUDA<T, D> const &data, Sp const &solver) {
  return CuDerivativeRegion<Dd, T, D, Sp>{data.values, data.region, data.dims,
                                          data.empty, solver};
}

template <typename Dd, Axis ax, typename T, size_t D, typename Sp>
auto to_cuda_expr_deriv(VectorComponentData<ax, T *, D> const &data,
                        Sp const &solver) {
  return CuDerivative<Dd, T, D, Sp>{data.values, solver};
}

template <typename Dd, typename T, size_t D, typename Sp>
auto to_cuda_expr_deriv(Grid<T, D> const &data, Sp const &solver) {
  return CuDerivativeHost<Dd, T, D, Sp>{data.values, data.len, solver};
}

template <typename Dd, typename G, typename Sp>
auto to_cuda_expr(
    OpDerivative<Dd, OpIdentity, OpTerm<OpIdentity, G>, Sp> const &e) {
  auto term = expr::get_enclosed_expression(e);
  auto &data = expr::BaseData<G>::get(expr::data(term));
  auto dimensions = expr::data_dimensions(term);

  using T = expr::eval_type_t<OpTerm<OpIdentity, G>>;

  return to_cuda_expr_deriv<Dd>(data, e.solver);
}

template <typename Dd, typename E, typename Sp>
auto to_cuda_expr(OpDerivative<Dd, OpIdentity, E, Sp> const &e) {
  using T = expr::eval_type_t<E>;
  static const size_t D = expr::grid_dim<E>::value;

  auto &result_data = expr::get_result_data(e);
  auto dimensions = expr::data_dimensions(result_data);

  return to_cuda_expr_deriv<Dd>(result_data, e.solver);
}

template <typename Dd, typename V, typename G, typename Sp>
auto to_cuda_expr(OpDerivative<Dd, V, OpTerm<OpIdentity, G>, Sp> const &e) {
  auto term = expr::get_enclosed_expression(e);
  auto &data = expr::BaseData<G>::get(expr::data(term));
  auto dimensions = expr::data_dimensions(term);

  using T = expr::eval_type_t<OpTerm<OpIdentity, G>>;

  return CuMul{to_cuda_expr(e.value), to_cuda_expr_deriv<Dd>(data, e.solver)};
}

template <typename Dd, typename V, typename E, typename Sp>
auto to_cuda_expr(OpDerivative<Dd, V, E, Sp> const &e) {
  using T = expr::eval_type_t<E>;
  static const size_t D = expr::grid_dim<E>::value;

  auto &result_data = expr::get_result_data(e);
  auto dimensions = expr::data_dimensions(result_data);

  return CuMul{to_cuda_expr(e.value),
               to_cuda_expr_deriv<Dd>(result_data, e.solver)};
}

template <typename V, typename E, typename Sp>
auto to_cuda_expr(OpIntegral<V, E, Sp> const &e) {
  return to_cuda_expr(e.value) * CuLiteral{expr::get_result_data(e)};
}

template <typename E>
auto to_cuda_expr(OpOperatorChain<OpIdentity, E> const &e) {
  return to_cuda_expr(e.g);
}

inline auto to_cuda_expr(scalar_t const &value) {
  return CuLiteral<scalar_t>{value};
}

inline auto to_cuda_expr(complex_t const &value) {
  return CuLiteral<complex_t>{value};
}

template <typename T, size_t D>
auto to_cuda_expr(any_vector_t<T, D> const &value) {
  return CuLiteral<any_vector_t<T, D>>{value};
}

auto to_cuda_expr(OpIdentity const &e) { return CuIdentity(); }

auto to_cuda_expr(OpNegIdentity const &e) { return CuNegIdentity(); }

auto to_cuda_expr(OpVoid const &e) { return CuVoid(); }

template <typename T>
auto to_cuda_expr(OpLiteral<T> const &e) {
  return to_cuda_expr(e.eval(0));
}

template <typename T, typename I>
auto to_cuda_expr(OpCoeff<T, I> const &e) {
  return to_cuda_expr(e.eval(0));
}

template <typename T>
auto to_cuda_expr(OpCoeff<T, DynamicIndex> const &e) {
  return to_cuda_expr(e.eval(0));
}

template <typename T, size_t... Ns>
auto to_cuda_expr(OpTensor<T, Ns...> const &e) {
  return to_cuda_expr(e.eval(0));
}

template <size_t N, size_t D>
auto to_cuda_expr(OpNegFractionLiteral<N, D> const &e) {
  return to_cuda_expr(e.eval(0));
}

template <size_t N, size_t D>
auto to_cuda_expr(OpFractionLiteral<N, D> const &e) {
  return to_cuda_expr(e.eval(0));
}

template <typename V, typename E, typename I, typename F, size_t... Ns,
          typename... Ts>
auto to_cuda_expr(SymbolicListIndex<E, I> const &list_index,
                  SymbolicFunction<F, Variable<Ns, Ts>...> const &func) {
  return CuVoid{};
}

// SymbolicListIndex<OpAdd<DynamicIndex, DynamicIndex, OpIdentity>,
// expr::symbols::i_<0, 0>>
template <typename V, typename E, typename I, typename F, size_t... Ns,
          typename... Ts>
auto to_cuda_expr(
    OpSymbolicEval<V, SymbolicListIndex<E, I>,
                   SymbolicFunction<F, Variable<Ns, Ts>...>> const &e) {
  return to_cuda_expr(expr::coeff(e)) * to_cuda_expr(e.data, e.f);
}

template <expr::NoiseType nt, typename T, size_t D, typename E, typename... Ts>
auto to_cuda_expr(NoiseData<nt, T, D> const &noise_data,
                  SymbolicFunction<E, Ts...> const &func) {
  return CuVoid{};
}

template <typename V, expr::NoiseType nt, typename T, size_t D, typename E,
          typename... Ts>
auto to_cuda_expr(OpSymbolicEval<V, NoiseData<nt, T, D>,
                                 SymbolicFunction<E, Ts...>> const &e) {
  return to_cuda_expr(expr::coeff(e)) * to_cuda_expr(e.data, e.f);
}

template <auto f, typename E>
auto to_cuda_fn_expr(CuEvaluable<E> const &e) {
  return CuFunctionApply<f, E>{*static_cast<E const *>(&e)};
}

template <auto f, typename E>
auto to_cuda_expr(OpFunctionApply<f, OpIdentity, E> const &e) {
  return to_cuda_fn_expr<f>(to_cuda_expr(e.e));
}

template <auto f, typename V, typename E>
auto to_cuda_expr(OpFunctionApply<f, V, E> const &e) {
  return to_cuda_expr(expr::coeff(e)) * to_cuda_fn_expr<f>(to_cuda_expr(e.e));
}

template <typename E>
auto to_cuda_expr(OpOptimized<E> const &e) {
  return CuVoid{};
  // return to_cuda_expr(e.e) * to_cuda_expr<1>(e.working);
}

template <typename... Es, size_t... Is>
auto to_cuda_expr(OpAdd<Es...> const &e, std::index_sequence<Is...>) {
  return CuAdd{to_cuda_expr(get<Is>(e))...};
}

template <typename... Es>
auto to_cuda_expr(OpAdd<Es...> const &e) {
  return to_cuda_expr(e, std::make_index_sequence<sizeof...(Es)>{});
}

template <typename E0, typename... Ts>
auto to_cuda_expr(SymbolicFunctionArray<E0, Ts...> const &persistent) {
  using cuda_expr_t = decltype(to_cuda_expr(persistent.data[0].e));
  using eval_type_t = expr::eval_type_t<E0>;
  cuda_expr_t *e;
  if (persistent.len == 0) {
    return CuSeries<eval_type_t, cuda_expr_t>{};
  }
  CHECK_CUDA_ERROR(cudaMalloc(&e, sizeof(cuda_expr_t) * persistent.len));
  for (iter_type i = 0; i < persistent.len; ++i) {
    auto expr = to_cuda_expr(persistent.data[i].e);
    CHECK_CUDA_ERROR(
        cudaMemcpy(&e[i], &expr, sizeof(cuda_expr_t), cudaMemcpyHostToDevice));
  }
  auto series = CuSeries<eval_type_t, cuda_expr_t>{&e, persistent.len};
  CHECK_CUDA_ERROR(cudaFree(&e[0]));
  return series;
}

template <typename... Ts, typename E, size_t I0, size_t P0, size_t I1,
          size_t P1, size_t... Is, size_t... Ps, typename... Rest>
auto to_cuda_expr(
    SymbolicSeries<expr::sum_op, Substitution<SymbolicDataArray<Ts>...>,
                   symphas::lib::types_list<
                       E,
                       symphas::lib::types_list<expr::symbols::i_<I0, P0>,
                                                expr::symbols::i_<I1, P1>,
                                                expr::symbols::i_<Is, Ps>...>,
                       Rest...>> const &series) {
  return to_cuda_expr(series.persistent);
}

template <typename T, typename E, size_t I0, size_t P0, typename... Rest>
auto to_cuda_expr(
    SymbolicSeries<expr::sum_op, Substitution<SymbolicDataArray<T>>,
                   symphas::lib::types_list<
                       E, symphas::lib::types_list<expr::symbols::i_<I0, P0>>,
                       Rest...>> const &series) {
  return to_cuda_expr(series.persistent);
}

template <typename V, typename... Ss, typename F>
auto to_cuda_expr(OpSymbolicEval<V, SymbolicSeries<Ss...>, F> const &e) {
  return to_cuda_expr(e.value) * to_cuda_expr(e.data);
}

}  // namespace expr

#endif
