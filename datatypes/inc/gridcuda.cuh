
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
 * PURPOSE: Defines the grid types used in the finite difference and other
 * semi-discrete numerical solvers.
 *
 * ***************************************************************************
 */

#pragma once

#include "grid.h"

#ifdef USING_CUDA

#include <cuda_runtime.h>

namespace symphas {
namespace cuda {}
}  // namespace symphas

namespace symphas::cuda {

template <typename T>
__global__ void initializeArray(T* array, T value, int size) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < size) {
    array[idx] = value;
  }
}

template <typename T, size_t D>
__global__ void initializeArray(any_vector_t<T, D>* array, T value, int size) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < size) {
    for (iter_type i = 0; i < D; ++i) {
      array[idx][i] = value;
    }
  }
}
template <typename T, size_t D>
__global__ void initializeArray(any_vector_t<T, D>* array,
                                any_vector_t<T, D> value, int size) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < size) {
    for (iter_type i = 0; i < D; ++i) {
      array[idx][i] = value[i];
    }
  }
}
void allocate(void** devPtr, size_t size);
template <typename T>
void allocate(T** devPtr, size_t size) {
  allocate(reinterpret_cast<void**>(devPtr), size);
}

template <typename T>
void allocate(T** devPtr, size_t size, T value) {
  allocate(devPtr, size);
  int numBlocks = (size + BLOCK_SIZE - 1) / BLOCK_SIZE;
  initializeArray CUDA_KERNEL(numBlocks, BLOCK_SIZE)(*devPtr, value, size);
  CHECK_CUDA_ERROR(cudaPeekAtLastError());
  CHECK_CUDA_ERROR(cudaDeviceSynchronize());
}
void memcpy(void* dest, const void* src, size_t count, cudaMemcpyKind kind);
void free(void* devPtr);
}  // namespace symphas::cuda

template <typename T>
struct carry_value_cuda {
  __device__ __host__ carry_value_cuda() : value{&fallback}, fallback{0} {}
  __device__ __host__ carry_value_cuda(T fallback)
      : value{&this->fallback}, fallback{fallback} {}
  __device__ __host__ carry_value_cuda(T* value, T fallback)
      : value{value}, fallback{fallback} {}
  __device__ __host__ carry_value_cuda(T* value, T fallback, bool valid)
      : value{this->value = (valid) ? value : &this->fallback},
        fallback{fallback} {}
  __device__ __host__ carry_value_cuda(carry_value_cuda<T> const& other)
      : carry_value_cuda{other.value, other.fallback, other.is_valid()} {}
  __device__ __host__ carry_value_cuda(carry_value_cuda<T>&& other)
      : carry_value_cuda() {
    swap(*this, other);
  }

  __host__ friend void swap(carry_value_cuda<T>& first,
                            carry_value_cuda<T>& second) {
    using std::swap;
    T* tmp = second.value;
    if (first.is_valid()) {
      second.value = first.value;
    } else {
      second.value = &second.fallback;
    }

    if (second.is_valid()) {
      first.value = tmp;
    } else {
      first.value = &first.fallback;
    }
    swap(first.fallback, second.fallback);
  }

  __host__ carry_value_cuda& operator=(T const& other) {
    if (is_valid()) {
      CHECK_CUDA_ERROR(
          cudaMemcpy(value, &other, sizeof(T), cudaMemcpyHostToDevice));
    }
    return *this;
  }

  __host__ carry_value_cuda& operator=(carry_value_cuda<T> const& other) {
    if (is_valid()) {
      if (other.is_valid()) {
        CHECK_CUDA_ERROR(cudaMemcpy(value, other.value, sizeof(T),
                                    cudaMemcpyDeviceToDevice));
      }
      CHECK_CUDA_ERROR(
          cudaMemcpy(value, other.value, sizeof(T), cudaMemcpyHostToDevice));
    }
    return *this;
  }

  __host__ __device__ T* operator&() const { return nullptr; }
  __host__ operator T() const {
    T v;
    if (is_valid()) {
      CHECK_CUDA_ERROR(
          cudaMemcpy(&v, value, sizeof(T), cudaMemcpyDeviceToHost));
    } else {
      v = fallback;
    }
    return v;
  }
  __host__ __device__ inline bool is_valid() const {
    return value != &fallback;
  }

  __device__ T dev() const { return (is_valid()) ? *value : fallback; }

  T* value;
  T fallback;
};

template <size_t N, typename T>
struct multi_value_cuda;

namespace grid {

template <typename T>
struct selected_entry_cuda {
  T* value;
  operator T() const {
    T host_value;
    CHECK_CUDA_ERROR(
        cudaMemcpy(&host_value, value, sizeof(T), cudaMemcpyDeviceToHost));
    return host_value;
  }

  selected_entry_cuda<T>& operator=(T const& new_value) {
    CHECK_CUDA_ERROR(
        cudaMemcpy(value, &new_value, sizeof(T), cudaMemcpyHostToDevice));
    return *this;
  }

  T operator!() const { return *this; }
};

/**
 * @brief A structure to manage CUDA device memory for arrays with a specific
 * stride.
 *
 * This structure provides functionality to copy values between host and device
 * memory with a specific stride, useful for multi-dimensional arrays.
 *
 * @tparam T The type of the elements.
 * @tparam D The dimension of the array.
 */
template <typename T, size_t D>
struct selected_entry_cuda<T[D]> {
  T* value[D];  ///< Pointer to the array in device memory.

  selected_entry_cuda(T* (&value)[D]) : value{} {
    for (iter_type i = 0; i < D; ++i) {
      this->value[i] = value[i];
    }
  }

  /**
   * @brief Converts the CUDA device memory to a host vector.
   *
   * This operator retrieves the values from the CUDA device memory and returns
   * them as a host vector. It performs a device-to-host memory copy to achieve
   * this.
   *
   * @return A host vector containing the values stored in the CUDA device
   * memory.
   */
  operator any_vector_t<T, D>() const {
    any_vector_t<T, D> host_value;
    for (iter_type i = 0; i < D; ++i) {
      CHECK_CUDA_ERROR(cudaMemcpy(&host_value[i], value[i], sizeof(T),
                                  cudaMemcpyDeviceToHost));
    }
    return host_value;
  }

  /**
   * @brief Assigns a new value to the CUDA device memory.
   *
   * This operator assigns a new value to the CUDA device memory by performing a
   * host-to-device memory copy.
   *
   * @param new_value The new value to be assigned to the CUDA device memory.
   * @return A reference to the updated selected_entry_cuda object.
   */
  selected_entry_cuda<T>& operator=(T const& new_value) {
    CHECK_CUDA_ERROR(
        cudaMemcpy(value, &new_value, sizeof(T), cudaMemcpyHostToDevice));
    return *this;
  }

  /**
   * @brief Returns the value stored in the CUDA device memory.
   *
   * This operator retrieves the value from the CUDA device memory and returns
   * it as a host value. It performs a device-to-host memory copy to achieve
   * this.
   *
   * @return The value stored in the CUDA device memory.
   */
  any_vector_t<T, D> operator!() const { return *this; }
};

template <size_t D>
struct select_grid_index_cuda {
  len_type stride[D];
  select_grid_index_cuda(const len_type (&dims)[D]) {
    stride[D - 1] = 1;
    for (iter_type i = D - 1; i > 0; --i) {
      stride[i - 1] = stride[i] * dims[i];
    }
  }

 private:
  template <typename T>
  selected_entry_cuda<T> get_value(T* values, iter_type n) {
    selected_entry_cuda<T> value{values + n};
    return value;
  }

  template <typename T>
  selected_entry_cuda<T[D]> get_value_arr(T* (&values)[D], iter_type n) {
    selected_entry_cuda<T[D]> value{values + n};
    return value;
  }

  template <typename T>
  selected_entry_cuda<T> get_value(const T* values, iter_type n) {
    selected_entry_cuda<T> value{values + n};
    return value;
  }

  template <typename T>
  selected_entry_cuda<T[D]> get_value_arr(T* const (&values)[D], iter_type n) {
    T* value_entry[D];
    for (iter_type i = 0; i < D; ++i) {
      value_entry[i] = values[i] + n;
    }
    selected_entry_cuda<T[D]> value{value_entry};
    return value;
  }

 public:
  template <typename T, typename... iter_types>
  decltype(auto) entry(const T* values, iter_types... indices) {
    iter_type n = 0;
    iter_type indices_array[] = {indices...};
    for (iter_type i = 0; i < D; ++i) {
      n += indices_array[i] * stride[i];
    }
    return get_value(values, n);
  }

  template <typename T, typename... iter_types>
  decltype(auto) entry_arr(T* const (&values)[D], iter_types... indices) {
    iter_type n = 0;
    iter_type indices_array[] = {indices...};
    for (iter_type i = 0; i < D; ++i) {
      n += indices_array[i] * stride[i];
    }
    return get_value_arr(values, n);
  }

  template <typename T, typename... iter_types>
  decltype(auto) entry(T* values, iter_types... indices) {
    iter_type n = 0;
    iter_type indices_array[] = {indices...};
    for (iter_type i = 0; i < D; ++i) {
      n += indices_array[i] * stride[i];
    }
    return get_value(values, n);
  }

  template <typename T, typename... iter_types>
  decltype(auto) entry_arr(T* (&values)[D], iter_types... indices) {
    iter_type n = 0;
    iter_type indices_array[] = {indices...};
    for (iter_type i = 0; i < D; ++i) {
      n += indices_array[i] * stride[i];
    }
    return get_value_arr(values, n);
  }
};

template <size_t D>
select_grid_index_cuda(const len_type (&)[D]) -> select_grid_index_cuda<D>;

template <size_t D>
struct select_region_cuda : select_region<D> {
  using parent_type = select_region<D>;
  using parent_type::parent_type;

  using parent_type::boundary_size;
  using parent_type::dims;
  using parent_type::len;
  using parent_type::offset;
  using parent_type::origin;
  using parent_type::stride;

  template <typename T>
  __device__ __host__ inline carry_value_cuda<T> operator()(
      const iter_type (&pos0)[D], T* values, const iter_type (&domain_dims)[D],
      const T& empty) const {
    iter_type pos[D];
    for (iter_type i = 0; i < D; ++i) {
      pos[i] = (pos0[i] >= origin[i] + boundary_size)
                   ? pos0[i] - origin[i]
                   : pos0[i] - origin[i] + domain_dims[i] - 2 * boundary_size;
    }
    if (is_in_region(pos, dims, boundary_size)) {
      return {&values[grid::index_from_position(pos, stride)], empty};
    } else {
      return {empty};
    }
  }

  template <typename T>
  __device__ __host__ inline carry_value_cuda<const T> operator()(
      const iter_type (&pos0)[D], const T* values,
      const iter_type (&domain_dims)[D], const T& empty) const {
    iter_type pos[D];
    for (iter_type i = 0; i < D; ++i) {
      pos[i] = (pos0[i] >= origin[i] + boundary_size)
                   ? pos0[i] - origin[i]
                   : pos0[i] - origin[i] + domain_dims[i] - 2 * boundary_size;
    }
    if (is_in_region(pos, dims, boundary_size)) {
      return {&values[grid::index_from_position(pos, stride)], empty};
    } else {
      return {empty};
    }
  }

  template <typename T>
  __device__ __host__ inline T& operator()(iter_type n, T* values,
                                           const iter_type (&domain_dims)[D],
                                           const T& empty) const {
    return values[(n - offset) % grid::length<D>(domain_dims)];
  }

  template <typename T>
  __device__ __host__ inline const T& operator()(
      iter_type n, const T* values, const iter_type (&domain_dims)[D],
      const T& empty) const {
    return values[(n - offset) % grid::length<D>(domain_dims)];
  }

  template <typename T>
  __device__ __host__ inline decltype(auto) operator()(
      T* values, const iter_type (&domain_dims)[D], const T& empty,
      const iter_type (&pos)[D]) const {
    iter_type pos0[D]{};
    for (iter_type i = 0; i < D; ++i) pos0[i] = pos[i];
    return operator()(pos0, values, domain_dims, empty);
  }

  template <typename T>
  __device__ __host__ inline decltype(auto) operator()(
      const T* values, const iter_type (&domain_dims)[D], const T& empty,
      const iter_type (&pos)[D]) const {
    iter_type pos0[D]{};
    for (iter_type i = 0; i < D; ++i) pos0[i] = pos[i];
    return operator()(pos0, values, domain_dims, empty);
  }

  template <typename T>
  __device__ __host__ inline multi_value_cuda<D, T> operator()(
      iter_type n, T* (&values)[D], const iter_type (&domain_dims)[D],
      T (&empty)[D]) const {
    multi_value_cuda<D, T> value;
    for (iter_type i = 0; i < D; ++i) {
      value[i] = operator()(n, values[i], domain_dims, empty[i]);
    }
    return value;
  }

  template <typename T>
  __device__ __host__ inline multi_value_cuda<D, T> operator()(
      iter_type n, T* const (&values)[D], const iter_type (&domain_dims)[D],
      const T (&empty)[D]) const {
    multi_value_cuda<D, T> value;
    for (iter_type i = 0; i < D; ++i) {
      value[i] = operator()(n, values[i], domain_dims, empty[i]);
    }
    return value;
  }

  template <typename T>
  __device__ __host__ inline multi_value_cuda<D, T> operator()(
      const iter_type (&pos0)[D], T* (&values)[D],
      const iter_type (&domain_dims)[D], const T (&empty)[D]) const {
    multi_value_cuda<D, T> value;
    iter_type pos[D];
    for (iter_type i = 0; i < D; ++i) {
      pos[i] = (pos0[i] >= origin[i] + boundary_size)
                   ? pos0[i] - origin[i]
                   : pos0[i] - origin[i] + domain_dims[i] - 2 * boundary_size;
    }
    if (is_in_region(pos, dims, boundary_size)) {
      return multi_value_cuda<D, T>(values,
                                    grid::index_from_position(pos, stride));
    } else {
      return {empty};
    }
  }

  template <typename T>
  __device__ __host__ inline multi_value_cuda<D, T> operator()(
      const iter_type (&pos0)[D], T* const (&values)[D],
      const iter_type (&domain_dims)[D], const T (&empty)[D]) const {
    iter_type pos[D];
    for (iter_type i = 0; i < D; ++i) {
      pos[i] = (pos0[i] >= origin[i] + boundary_size)
                   ? pos0[i] - origin[i]
                   : pos0[i] - origin[i] + domain_dims[i] - 2 * boundary_size;
    }
    if (is_in_region(pos, dims, boundary_size)) {
      return multi_value_cuda<D, T>(values,
                                    grid::index_from_position(pos, stride));
    } else {
      return {empty};
    }
  }

  template <typename T>
  __device__ __host__ inline multi_value_cuda<D, T> operator()(
      T* (&values)[D], const iter_type (&domain_dims)[D], const T (&empty)[D],
      const iter_type (&pos)[D]) const {
    iter_type pos0[D]{};
    for (iter_type i = 0; i < D; ++i) pos0[i] = pos[i];
    return operator()(pos0, values, domain_dims, empty);
  }

  template <typename T>
  __device__ __host__ inline multi_value_cuda<D, T> operator()(
      const T* (&values)[D], const iter_type (&domain_dims)[D],
      const T (&empty)[D], const iter_type (&pos)[D]) const {
    iter_type pos0[D]{};
    for (iter_type i = 0; i < D; ++i) pos0[i] = pos[i];
    return operator()(pos0, values, domain_dims, empty);
  }

  template <
      typename T, typename T0, typename... Ts,
      std::enable_if_t<symphas::are_all_same_v<iter_type, Ts...>, int> = 0>
  __device__ __host__ decltype(auto) operator()(
      T&& values, const iter_type (&domain_dims)[D], T0 empty, iter_type coord0,
      Ts&&... coords) const {
    iter_type pos[D]{coord0, std::forward<Ts>(coords)...};
    return operator()(pos, std::forward<T>(values), domain_dims, empty);
  }
};

template <typename T, size_t D>
void adjust_region_to_from_cuda(T* new_values, const iter_type (&new_origin)[D],
                                const len_type (&new_dims)[D],
                                const T* old_values,
                                const iter_type (&old_origin)[D],
                                const len_type (&old_dims)[D],
                                const len_type (&global_dims)[D], T empty,
                                len_type boundary_size = 0);

template <typename T, size_t D>
void adjust_origin_to_from_cuda(T*(&values), const iter_type (&new_origin)[D],
                                const iter_type (&old_origin)[D],
                                const len_type (&dims)[D],
                                const len_type (&global_dims)[D], T empty,
                                len_type boundary_size = 0);

}  // namespace grid

template <size_t N, typename T>
struct multi_value_cuda {
  T* value[N];

  __host__ __device__ multi_value_cuda(T* (&value)[N], iter_type offset = 0)
      : multi_value_cuda() {
    for (iter_type i = 0; i < N; ++i) {
      this->value[i] = value[i] + offset;
    }
  }

  __host__ __device__ multi_value_cuda(T* const (&value)[N],
                                       iter_type offset = 0)
      : multi_value_cuda() {
    for (iter_type i = 0; i < N; ++i) {
      this->value[i] = value[i] + offset;
    }
  }

  __host__ __device__ multi_value_cuda(T* const value) : multi_value_cuda() {
    for (iter_type i = 0; i < N; ++i) {
      this->value[i] = &value[i];
    }
  }

  __host__ __device__ multi_value_cuda(T (&value)[N]) : multi_value_cuda() {
    for (iter_type i = 0; i < N; ++i) {
      this->value[i] = &value[i];
    }
  }

  __host__ __device__ multi_value_cuda(const T (&value)[N])
      : multi_value_cuda() {
    for (iter_type i = 0; i < N; ++i) {
      this->value[i] = &const_cast<T&>(value[i]);
    }
  }

  __device__ const T& operator[](iter_type i) const { return *(value[i]); }

  __device__ T& operator[](iter_type i) { return *(value[i]); }

  //__host__ __device__ T operator[](iter_type i) const {
  //  T v;
  //  symphas::cuda::memcpy(&v, value[i], sizeof(T), cudaMemcpyDeviceToHost);
  //  return v;
  //}

  __host__ __device__ bool is_valid() { return true; }

  //! Return a multi_value_cuda as a vector for compatibility.
  __host__ operator any_vector_t<T, N>() const;

  __device__ any_vector_t<T, N> dev() const {
    any_vector_t<T, N> v;
    for (iter_type i = 0; i < N; ++i) {
      v[i] = *(value[i]);
    }
    return v;
  }

  //! Set the values of the multi_value_cuda from a vector.
  /*!
   * Assigning the value from a vector can update the values
   * in the MultiBlock list, using the pointers stored by
   * the multi_value_cuda object.
   */
  __host__ multi_value_cuda& operator=(any_vector_t<T, N> const& vector);
  __host__ multi_value_cuda& operator=(any_vector_t<T, N>&& vector);
  __host__ multi_value_cuda& operator=(multi_value_cuda<N, T> const& other);
  __host__ multi_value_cuda& operator=(multi_value_cuda<N, T>&& other);
  __host__ multi_value_cuda& operator=(T const& other);

  __host__ multi_value_cuda(multi_value_cuda<N, T> const& other)
      : multi_value_cuda(other.value) {}
  __host__ multi_value_cuda(multi_value_cuda<N, T>&& other)
      : multi_value_cuda(other.value) {}

  //! Set the values of the multi_value_cuda from a vector.
  /*!
   * Assigning the value from a vector can update the values
   * in the MultiBlock list, using the pointers stored by
   * the multi_value_cuda object.
   */
  // multi_value_cuda& operator=(multi_value_cuda other);

  __host__ __device__ multi_value_cuda() : value{} {}
};

template <size_t N, typename T>
__host__ multi_value_cuda<N, T>::operator any_vector_t<T, N>() const {
  any_vector_t<T, N> vector;

  for (iter_type i = 0; i < N; ++i) {
    CHECK_CUDA_ERROR(
        cudaMemcpy(&vector[i], value[i], sizeof(T), cudaMemcpyDeviceToHost));
  }

  return vector;
}

template <size_t N, typename T>
__host__ multi_value_cuda<N, T>& multi_value_cuda<N, T>::operator=(
    any_vector_t<T, N> const& vector) {
  for (iter_type i = 0; i < N; ++i) {
    CHECK_CUDA_ERROR(
        cudaMemcpy(value[i], &vector[i], sizeof(T), cudaMemcpyHostToDevice));
  }

  return *this;
}

template <size_t N, typename T>
__host__ multi_value_cuda<N, T>& multi_value_cuda<N, T>::operator=(
    any_vector_t<T, N>&& vector) {
  for (iter_type i = 0; i < N; ++i) {
    CHECK_CUDA_ERROR(
        cudaMemcpy(value[i], &vector[i], sizeof(T), cudaMemcpyHostToDevice));
  }

  return *this;
}

template <size_t N, typename T>
__host__ multi_value_cuda<N, T>& multi_value_cuda<N, T>::operator=(
    multi_value_cuda<N, T> const& other) {
  any_vector_t<T, N> vector = other;
  *this = vector;
  return *this;
}

template <size_t N, typename T>
__host__ multi_value_cuda<N, T>& multi_value_cuda<N, T>::operator=(
    multi_value_cuda<N, T>&& other) {
  any_vector_t<T, N> vector = other;
  *this = vector;
  return *this;
}

template <size_t N, typename T>
__host__ multi_value_cuda<N, T>& multi_value_cuda<N, T>::operator=(
    T const& other) {
  any_vector_t<T, N> vector;
  ;
  for (iter_type i = 0; i < N; ++i) {
    vector[i] = other;
  }
  *this = vector;
  return *this;
}

template <size_t N, typename T, typename V>
__host__ auto operator*(multi_value_cuda<N, T> const& a, V&& b) {
  return any_vector_t<T, N>(a) * std::forward<V>(b);
}

template <size_t N, typename T, typename V>
__host__ auto operator*(V&& a, multi_value_cuda<N, T> const& b) {
  return std::forward<V>(a) * any_vector_t<T, N>(b);
}

template <size_t N, typename T, typename V>
__host__ auto operator+(multi_value_cuda<N, T> const& a,
                        multi_value_cuda<N, V> const& b) {
  return any_vector_t<T, N>(a) + any_vector_t<V, N>(a);
}

template <size_t N, typename T, typename V>
__host__ auto operator-(multi_value_cuda<N, T> const& a,
                        multi_value_cuda<N, V> const& b) {
  return any_vector_t<T, N>(a) - any_vector_t<V, N>(a);
}

//! Manages an array of values of arbitrary type.
/*!
 * Basic array type object used in constructing the finite difference grid.
 * Values are always initialized to the empty value.
 *
 * \tparam T The value type of the underlying array.
 */
template <typename T>
struct BlockCUDA {
  T* values;     //!< The list of values managed by this object.
  len_type len;  //!< The number of values in the list.

  //! Create this object with \p len values.
  /*!
   * Create this object with \p len values. The length can be 0, in which case
   * no memory will be allocated, no values can be accessed, but the
   * object can still be constructed.
   *
   * \param len The number of values to create.
   */
  BlockCUDA(len_type len) : values{nullptr}, len{len} {
    if (len > 0) {
      CHECK_CUDA_ERROR(cudaMalloc(&values, len * sizeof(T)));
    }
  }

  explicit BlockCUDA(const len_type* len, size_t dimensions = 1)
      : BlockCUDA((len != nullptr) ? grid::length(len, dimensions) : 0) {}
  explicit BlockCUDA(grid::dim_list const& dims)
      : BlockCUDA(dims.data, dims.n) {}

  BlockCUDA(BlockCUDA<T> const& other) : BlockCUDA(other.len) {
    // std::copy(other.values, other.values + other.len, values);
    CHECK_CUDA_ERROR(cudaMemcpy(values, other.values, other.len * sizeof(T),
                                cudaMemcpyDeviceToDevice));
  }

  BlockCUDA(BlockCUDA<T>&& other) noexcept : BlockCUDA() { swap(*this, other); }

  BlockCUDA& operator=(BlockCUDA<T> other) {
    swap(*this, other);
    return *this;
  }

  friend void swap(BlockCUDA<T>& first, BlockCUDA<T>& second) {
    using std::swap;
    swap(first.len, second.len);
    swap(first.values, second.values);
  }

  //! Return the value at the \p i index in the list.
  /*!
   * Return the value from the list from index \p i.
   *
   * \param i The index of the value to return.
   */
  __host__ __device__ const T& operator[](iter_type i) const {
    return values[i];
  }

  //! Return the value at the \p i index in the list.
  /*!
   * Return the value from the list from index \p i.
   *
   * \param i The index of the value to return.
   */
  __host__ __device__ T& operator[](iter_type i) { return values[i]; }

  __device__ operator const T*() const { return values; }

  __device__ operator T*() { return values; }

  operator Block<T>() {
    Block<T> block(len);
    CHECK_CUDA_ERROR(cudaMemcpy(block.values, values, len * sizeof(T),
                                cudaMemcpyDeviceToHost));
    return block;
  }

  decltype(auto) operator()(iter_type i) {
    return grid::selected_entry_cuda<T>{values + i};
  }

  ~BlockCUDA() { CHECK_CUDA_ERROR(cudaFree(values)); }

  BlockCUDA() : values{nullptr}, len{0} {}
};

// ***********************************************************************************************

//! Manages an array of values of arbitrary type.
/*!
 * Basic array type object used in constructing the finite difference grid.
 * Values are always initialized to the empty value.
 *
 * \tparam T The value type of the underlying array.
 */
template <size_t N, typename T>
struct MultiBlockCUDA {
  T* values[N];  //!< The list of values managed by this object.
  len_type len;  //!< The number of values in the list.

  //! Create this object with \p len values.
  /*!
   * Create this object with \p len values. The length can be 0, in which case
   * no memory will be allocated, no values can be accessed, but the
   * object can still be constructed.
   *
   * \param len The number of values to create.
   */
  MultiBlockCUDA(len_type len) : values{}, len{len} {
    for (iter_type n = 0; n < N; ++n) {
      if (len > 0) {
        CHECK_CUDA_ERROR(cudaMalloc(&values[n], len * sizeof(T)));
      } else {
        values[n] = nullptr;
      }
    }
  }

  explicit MultiBlockCUDA(const len_type* len, size_t dimensions = 1)
      : MultiBlockCUDA((len != nullptr) ? grid::length(len, dimensions) : 0) {}

  MultiBlockCUDA(MultiBlockCUDA<N, T> const& other)
      : MultiBlockCUDA(other.len) {
    for (iter_type n = 0; n < N; ++n) {
      // std::copy(other.values[n], other.values[n] + other.len, values[n]);
      CHECK_CUDA_ERROR(cudaMemcpy(values[n], other.values[n],
                                  other.len * sizeof(T),
                                  cudaMemcpyDeviceToDevice));
    }
  }

  MultiBlockCUDA(MultiBlockCUDA<N, T>&& other) noexcept : MultiBlockCUDA() {
    swap(*this, other);
  }

  MultiBlockCUDA& operator=(MultiBlockCUDA<N, T> other) {
    swap(*this, other);
    return *this;
  }

  friend void swap(MultiBlockCUDA<N, T>& first, MultiBlockCUDA<N, T>& second) {
    using std::swap;
    swap(first.len, second.len);
    swap(first.values, second.values);
  }

  //! Return the value at the \p i index in the list.
  /*!
   * Return the value from the list from index \p i.
   *
   * \param i The index of the value to return.
   */
  __host__ __device__ multi_value_cuda<N, T> operator[](iter_type i) const {
    multi_value_cuda<N, T> value;
    for (iter_type n = 0; n < N; ++n) {
      value.value[n] = values[n] + i;
    }
    return value;
  }

  __device__ T* operator()(Axis ax) const {
    return values[symphas::axis_to_index(ax)];
  }

  operator MultiBlock<N, T>() {
    MultiBlock<N, T> block(len);
    for (iter_type n = 0; n < N; ++n) {
      CHECK_CUDA_ERROR(cudaMemcpy(block.values[n], values[n], len * sizeof(T),
                                  cudaMemcpyDeviceToHost));
    }
    return block;
  }

  ~MultiBlockCUDA() {
    for (iter_type n = 0; n < N; ++n) {
      CHECK_CUDA_ERROR(cudaFree(values[n]));
    }
  }

  MultiBlockCUDA() : values{0}, len{0} {}
};

// ***********************************************************************************************

//! A grid object of arbitrary dimension and arbitrary value type.
/*!
 * A grid object of arbitrary dimension and arbitrary value type. The grid
 * forms the basis of a finite difference grid. Only manages its own
 * dimensions and list of values. The values are inherited from ::Block, meaning
 * that the data is flattened and is not `D`-dimensional in memory.
 *
 * \tparam T The value type of the underlying array.
 * \tparam D The dimension of the grid.
 */
template <typename T, size_t D>
struct GridCUDA : BlockCUDA<T> {
 private:
  void set_dimensions(const len_type* dimensions) {
    if (dimensions == nullptr) {
      std::fill(dims, dims + D, 0);
    } else {
      for (iter_type i = 0; i < D; ++i) {
        dims[i] = dimensions[i];
      }
    }
  }

  void set_dimensions(grid::dim_list dimensions) {
    std::fill(dims, dims + D, 0);
    for (iter_type i = 0; i < dimensions.n; ++i) {
      dims[i] = dimensions[i];
    }
  }

 public:
  len_type dims[D];  //!< Dimensions of the grid, arranged beginning from the
                     //!< horizontal coordinate.

  //! Create a grid of the prescribed dimensions.
  /*!
   * Creates a new grid using the given dimensions. The number of values in
   * the grid directly correspond to the dimensions, equal to the product
   * of all dimensions.
   *
   * \param dimensions The dimensions of the grid.
   */
  GridCUDA(grid::dim_list dimensions)
      : BlockCUDA<T>{grid::length<D>(dimensions)}, dims{0} {
    set_dimensions(dimensions);
  }

  //! Create a grid of the prescribed dimensions.
  /*!
   * Creates a new grid using the given dimensions. The number of values in
   * the grid directly correspond to the dimensions, equal to the product
   * of all dimensions.
   *
   * \param dimensions The dimensions of the grid.
   */
  GridCUDA(const len_type* dimensions)
      : BlockCUDA<T>{grid::length<D>(dimensions)}, dims{0} {
    set_dimensions(dimensions);
  }

  const GridCUDA<T, D>& as_grid() const { return *this; }

  GridCUDA<T, D>& as_grid() { return *this; }

  operator Grid<T, D>() {
    Grid<T, D> grid(dims);
    CHECK_CUDA_ERROR(cudaMemcpy(grid.values, BlockCUDA<T>::values,
                                grid::length<D>(dims) * sizeof(T),
                                cudaMemcpyDeviceToHost));
    return grid;
  }

  template <typename... Ts,
            std::enable_if_t<
                (sizeof...(Ts) == D &&
                 std::conjunction_v<std::is_convertible<Ts, iter_type>...>),
                int> = 0>
  decltype(auto) operator()(Ts&&... indices) const {
    return grid::select_grid_index_cuda(dims).entry(
        BlockCUDA<T>::values, std::forward<Ts>(indices)...);
  }

 protected:
  constexpr GridCUDA() : GridCUDA(nullptr) {}
};

template <typename T, size_t D>
struct GridCUDA<any_vector_t<T, D>, D> : MultiBlockCUDA<D, T> {
  using parent_type = MultiBlockCUDA<D, T>;
  using parent_type::parent_type;

 private:
  void set_dimensions(const len_type* dimensions) {
    if (dimensions == nullptr) {
      std::fill(dims, dims + D, 0);
    } else {
      for (iter_type i = 0; i < D; ++i) {
        dims[i] = dimensions[i];
      }
    }
  }

 public:
  len_type dims[D];  //!< Dimensions of the grid, arranged beginning from the
                     //!< horizontal coordinate.

  //! Create a grid of the prescribed dimensions.
  /*!
   * Creates a new grid using the given dimensions. The number of values in
   * the grid directly correspond to the dimensions, equal to the product
   * of all dimensions.
   *
   * \param dimensions The dimensions of the grid.
   */
  GridCUDA(grid::dim_list dimensions)
      : parent_type{grid::length<D>(dimensions)}, dims{0} {
    set_dimensions(dimensions);
  }

  //! Create a grid of the prescribed dimensions.
  /*!
   * Creates a new grid using the given dimensions. The number of values in
   * the grid directly correspond to the dimensions, equal to the product
   * of all dimensions.
   *
   * \param dimensions The dimensions of the grid.
   */
  GridCUDA(const len_type* dimensions)
      : parent_type{grid::length<D>(dimensions)}, dims{0} {
    set_dimensions(dimensions);
  }

  const GridCUDA<any_vector_t<T, D>, D>& as_grid() const { return *this; }

  GridCUDA<any_vector_t<T, D>, D>& as_grid() { return *this; }

  operator Grid<any_vector_t<T, D>, D>() {
    Grid<any_vector_t<T, D>, D> grid(dims);
    for (iter_type n = 0; n < D; ++n) {
      CHECK_CUDA_ERROR(cudaMemcpy(grid.values[n], parent_type::values[n],
                                  grid::length<D>(dims) * sizeof(T),
                                  cudaMemcpyDeviceToHost));
    }
    return grid;
  }

  const T* axis(Axis ax) const {
    return parent_type::values[symphas::axis_to_index(ax)];
  }

  T* axis(Axis ax) { return parent_type::values[symphas::axis_to_index(ax)]; }

  template <typename... Ts,
            std::enable_if_t<
                (sizeof...(Ts) == D &&
                 std::conjunction_v<std::is_convertible<Ts, iter_type>...>),
                int> = 0>
  decltype(auto) operator()(Ts&&... rest) const {
    return grid::select_grid_index_cuda(dims).entry_arr(
        MultiBlockCUDA<D, T>::values, std::forward<Ts>(rest)...);
  }

 protected:
  constexpr GridCUDA() : GridCUDA(nullptr) {}
};

//! A grid object of arbitrary dimension and arbitrary value type.
/*!
 * A grid object of arbitrary dimension and arbitrary value type.
 * This grid implementation is an extension of the base ::Grid, but it
 * contains virtual boundary conditions.
 * The grid
 * forms the basis of a finite difference grid. Only manages its own
 * dimensions and list of values. The values are inherited from ::BlockCUDA,
 * meaning that the data is flattened and is not `D`-dimensional in memory.
 *
 * \tparam T The value type of the underlying array.
 * \tparam D The dimension of the grid.
 */
template <typename T, size_t D>
struct BoundaryGridCUDA : GridCUDA<T, D> {
  using parent_type = GridCUDA<T, D>;
  using parent_type::dims;

  BoundaryGridCUDA(grid::dim_list dimensions)
      : BoundaryGridCUDA(static_cast<const len_type*>(dimensions)) {}
  BoundaryGridCUDA(const len_type* dimensions) : GridCUDA<T, D>{dimensions} {}
  BoundaryGridCUDA(GridCUDA<T, D> const& other) : GridCUDA<T, D>{other} {}
  BoundaryGridCUDA(GridCUDA<T, D>&& other) noexcept : GridCUDA<T, D>{other} {}

  const BoundaryGridCUDA<T, D>& as_grid() const { return *this; }

  BoundaryGridCUDA<T, D>& as_grid() { return *this; }

  template <typename... Ts,
            std::enable_if_t<
                (sizeof...(Ts) == D &&
                 std::conjunction_v<std::is_convertible<Ts, iter_type>...>),
                int> = 0>
  decltype(auto) operator()(Ts&&... rest) const {
    return GridCUDA<T, D>::operator()(
        (std::forward<Ts>(rest) + BOUNDARY_DEPTH)...);
  }

 protected:
  constexpr BoundaryGridCUDA() : GridCUDA<T, D>() {}
};

template <typename T>
struct BoundaryGridCUDA<T, 3> : GridCUDA<T, 3> {
  using parent_type = GridCUDA<T, 3>;
  using parent_type::dims;

  BoundaryGridCUDA(grid::dim_list dimensions)
      : BoundaryGridCUDA(static_cast<const len_type*>(dimensions)) {}
  BoundaryGridCUDA(const len_type* dimensions) : GridCUDA<T, 3>(dimensions) {}
  BoundaryGridCUDA(GridCUDA<T, 3> const& other) : GridCUDA<T, 3>(other) {}
  BoundaryGridCUDA(GridCUDA<T, 3>&& other) noexcept : GridCUDA<T, 3>(other) {}

  const BoundaryGridCUDA<T, 3>& as_grid() const { return *this; }

  BoundaryGridCUDA<T, 3>& as_grid() { return *this; }

  decltype(auto) operator()(iter_type x, iter_type y, iter_type z) const {
    return GridCUDA<T, 3>::operator()(x + BOUNDARY_DEPTH, y + BOUNDARY_DEPTH,
                                      z + BOUNDARY_DEPTH);
  }

 protected:
  constexpr BoundaryGridCUDA() : GridCUDA<T, 3>() {}
};

template <typename T>
struct BoundaryGridCUDA<T, 2> : GridCUDA<T, 2> {
  using parent_type = GridCUDA<T, 2>;
  using parent_type::dims;

  BoundaryGridCUDA(grid::dim_list dimensions) : GridCUDA<T, 2>(dimensions) {}
  BoundaryGridCUDA(const len_type* dimensions) : GridCUDA<T, 2>(dimensions) {}
  BoundaryGridCUDA(GridCUDA<T, 2> const& other) : GridCUDA<T, 2>(other) {}
  BoundaryGridCUDA(GridCUDA<T, 2>&& other) noexcept : GridCUDA<T, 2>(other) {}

  const BoundaryGridCUDA<T, 2>& as_grid() const { return *this; }

  BoundaryGridCUDA<T, 2>& as_grid() { return *this; }

  decltype(auto) operator()(iter_type x, iter_type y) const {
    return GridCUDA<T, 2>::operator()(x + BOUNDARY_DEPTH, y + BOUNDARY_DEPTH);
  }

 protected:
  constexpr BoundaryGridCUDA() : GridCUDA<T, 2>() {}
};

template <typename T>
struct BoundaryGridCUDA<T, 1> : GridCUDA<T, 1> {
  using parent_type = GridCUDA<T, 1>;
  using parent_type::dims;

  BoundaryGridCUDA(grid::dim_list dimensions)
      : BoundaryGridCUDA(static_cast<const len_type*>(dimensions)) {}
  BoundaryGridCUDA(const len_type* dimensions) : GridCUDA<T, 1>(dimensions) {}
  BoundaryGridCUDA(GridCUDA<T, 1> const& other) : GridCUDA<T, 1>(other) {}
  BoundaryGridCUDA(GridCUDA<T, 1>&& other) noexcept : GridCUDA<T, 1>(other) {}

  const BoundaryGridCUDA<T, 1>& as_grid() const { return *this; }

  BoundaryGridCUDA<T, 1>& as_grid() { return *this; }

  decltype(auto) operator()(iter_type x) const {
    return GridCUDA<T, 1>::operator()(x + BOUNDARY_DEPTH);
  }

 protected:
  constexpr BoundaryGridCUDA() : GridCUDA<T, 1>() {}
};

//! A grid object of arbitrary dimension and arbitrary value type.
/*!
 * A grid object of arbitrary dimension and arbitrary value type.
 * This grid implementation is an extension of the base ::Grid, but it
 * is meant to sub-domain a larger domain when there are many fields.
 * The grid
 * forms the basis of a finite difference grid. Only manages its own
 * dimensions and list of values. The values are inherited from ::Block, meaning
 * that the data is flattened and is not `D`-dimensional in memory.
 *
 * \tparam T The value type of the underlying array.
 * \tparam D The dimension of the grid.
 */
template <typename T, size_t D>
struct RegionalGridCUDA : GridCUDA<T, D> {
  using parent_type = GridCUDA<T, D>;
  using parent_type::dims;

  grid::select_region_cuda<D> region;
  T empty;

 public:
  RegionalGridCUDA(grid::dim_list dimensions,
                   T empty = REGIONAL_GRID_OUTER_VALUE,
                   len_type boundary_size = BOUNDARY_DEPTH)
      : RegionalGridCUDA(static_cast<const len_type*>(dimensions), empty,
                         boundary_size) {}
  RegionalGridCUDA(const len_type* dimensions,
                   T empty = REGIONAL_GRID_OUTER_VALUE,
                   len_type boundary_size = BOUNDARY_DEPTH)
      : parent_type{dimensions},
        region{parent_type::dims, boundary_size},
        empty{empty} {}

  RegionalGridCUDA(GridCUDA<T, D> const& other)
      : parent_type{other},
        region{parent_type::dims, BOUNDARY_DEPTH},
        empty{REGIONAL_GRID_OUTER_VALUE} {}
  RegionalGridCUDA(GridCUDA<T, D>&& other) noexcept
      : parent_type{other},
        region{parent_type::dims, BOUNDARY_DEPTH},
        empty{REGIONAL_GRID_OUTER_VALUE} {}

  RegionalGridCUDA(RegionalGridCUDA<T, D> const& other)
      : RegionalGridCUDA(nullptr, other.empty, other.region.boundary_size) {
    using std::swap;

    parent_type::len = other.len;
    std::copy(other.dims, other.dims + D, parent_type::dims);
    region = other.region;

    BlockCUDA<T> tmp(other.region.len);
    CHECK_CUDA_ERROR(cudaMemcpy(tmp.values, other.values,
                                other.region.len * sizeof(T),
                                cudaMemcpyDeviceToDevice));
    swap(tmp.values, parent_type::values);
  }

  RegionalGridCUDA(RegionalGridCUDA<T, D>&& other) noexcept
      : RegionalGridCUDA() {
    swap(*this, other);
  }
  RegionalGridCUDA& operator=(RegionalGridCUDA<T, D> other) {
    swap(*this, other);
    return *this;
  }

  friend void swap(RegionalGridCUDA<T, D>& first,
                   RegionalGridCUDA<T, D>& second) {
    using std::swap;
    swap(static_cast<parent_type&>(first), static_cast<parent_type&>(second));
    swap(first.empty, second.empty);
    swap(first.region, second.region);
  }

  const RegionalGridCUDA<T, D>& as_grid() const { return *this; }

  RegionalGridCUDA<T, D>& as_grid() { return *this; }

  template <typename... Ts,
            std::enable_if_t<
                (sizeof...(Ts) == D &&
                 std::conjunction_v<std::is_convertible<Ts, iter_type>...>),
                int> = 0>
  __device__ __host__ decltype(auto) operator()(Ts&&... rest) const {
    return region(BlockCUDA<T>::values, parent_type::dims, empty,
                  (std::forward<Ts>(rest) + region.boundary_size)...);
  }

  __device__ __host__ decltype(auto) operator[](iter_type n) const {
    iter_type pos[D]{};
    grid::get_grid_position(pos, parent_type::dims, n);
    return region(pos, BlockCUDA<T>::values, parent_type::dims, empty);
  }

  __device__ __host__ T const& get_unsafe(iter_type n) const {
    return region(n, BlockCUDA<T>::values, parent_type::dims, empty);
  }

  __device__ __host__ T& get_unsafe(iter_type n) {
    return region(n, BlockCUDA<T>::values, parent_type::dims, empty);
  }

  void adjust(const iter_type (&new_origin)[D]) {
    grid::adjust_origin_to_from_cuda(
        parent_type::values, new_origin, region.origin, region.dims,
        parent_type::dims, empty, region.boundary_size);
    region.update(new_origin);
  }

  void adjust(const iter_type (&new_origin)[D], const len_type (&new_dims)[D]) {
    T* new_values;

    CHECK_CUDA_ERROR(
        cudaMalloc(&new_values, grid::length<D>(new_dims) * sizeof(T)));
    int numBlocks = (grid::length<D>(new_dims) + BLOCK_SIZE - 1) / BLOCK_SIZE;
    symphas::cuda::initializeArray CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
        new_values, empty, grid::length<D>(new_dims));
    CHECK_CUDA_ERROR(cudaPeekAtLastError());
    CHECK_CUDA_ERROR(cudaDeviceSynchronize());

    grid::adjust_region_to_from_cuda(
        new_values, new_origin, new_dims, parent_type::values, region.origin,
        region.dims, parent_type::dims, empty, region.boundary_size);

    using std::swap;
    swap(parent_type::values, new_values);
    CHECK_CUDA_ERROR(cudaFree(new_values));

    region.update(new_origin, new_dims);
  }

  void adjust(grid::select_region_cuda<D> const& other) {
    if (!grid::is_same(region.dims, other.dims) ||
        !grid::is_same(region.origin, other.origin)) {
      if (grid::is_same(region.dims, other.dims)) {
        adjust(other.origin);
      } else {
        adjust(other.origin, other.dims);
      }
    }
  }

 protected:
  constexpr RegionalGridCUDA() : parent_type(), region{}, empty{} {}
};

template <typename T, size_t D>
struct RegionalGridCUDA<any_vector_t<T, D>, D>
    : GridCUDA<any_vector_t<T, D>, D> {
  using parent_type = GridCUDA<any_vector_t<T, D>, D>;

  grid::select_region_cuda<D> region;
  T empty[D];

 public:
  template <size_t... Is>
  RegionalGridCUDA(const len_type* dimensions, const T (&empty)[D],
                   len_type boundary_size, std::index_sequence<Is...>)
      : parent_type{dimensions},
        region{parent_type::dims, boundary_size},
        empty{empty[Is]...} {}
  RegionalGridCUDA(const len_type* dimensions, const T (&empty)[D],
                   len_type boundary_size = BOUNDARY_DEPTH)
      : RegionalGridCUDA(static_cast<const len_type*>(dimensions), empty,
                         boundary_size, std::make_index_sequence<D>{}) {}
  RegionalGridCUDA(grid::dim_list dimensions, const T (&empty)[D],
                   len_type boundary_size = BOUNDARY_DEPTH)
      : RegionalGridCUDA(static_cast<const len_type*>(dimensions), empty,
                         boundary_size) {}

  template <size_t... Is>
  RegionalGridCUDA(const len_type* dimensions, T empty, len_type boundary_size,
                   std::index_sequence<Is...>)
      : RegionalGridCUDA(dimensions,
                         symphas::array_values_type<T, D>{
                             symphas::lib::repeat_value<Is>(empty)...},
                         boundary_size, std::index_sequence<Is...>{}) {}
  RegionalGridCUDA(const len_type* dimensions,
                   T empty = REGIONAL_GRID_OUTER_VALUE,
                   len_type boundary_size = BOUNDARY_DEPTH)
      : RegionalGridCUDA(dimensions, empty, boundary_size,
                         std::make_index_sequence<D>{}) {}
  RegionalGridCUDA(grid::dim_list dimensions,
                   T empty = REGIONAL_GRID_OUTER_VALUE,
                   len_type boundary_size = BOUNDARY_DEPTH)
      : RegionalGridCUDA(static_cast<const len_type*>(dimensions), empty,
                         boundary_size) {}

  template <size_t... Is>
  RegionalGridCUDA(GridCUDA<any_vector_t<T, D>, D> const& other,
                   std::index_sequence<Is...>)
      : parent_type{other},
        region{parent_type::dims, BOUNDARY_DEPTH},
        empty{symphas::lib::repeat_value<Is>(REGIONAL_GRID_OUTER_VALUE)...} {}
  template <size_t... Is>
  RegionalGridCUDA(GridCUDA<any_vector_t<T, D>, D>&& other,
                   std::index_sequence<Is...>) noexcept
      : parent_type{other},
        region{parent_type::dims, BOUNDARY_DEPTH},
        empty{symphas::lib::repeat_value<Is>(REGIONAL_GRID_OUTER_VALUE)...} {}

  RegionalGridCUDA(GridCUDA<any_vector_t<T, D>, D> const& other)
      : RegionalGridCUDA(other, std::make_index_sequence<D>{}) {}
  RegionalGridCUDA(GridCUDA<any_vector_t<T, D>, D>&& other) noexcept
      : RegionalGridCUDA(other, std::make_index_sequence<D>{}) {}

  RegionalGridCUDA(RegionalGridCUDA<any_vector_t<T, D>, D> const& other)
      : RegionalGridCUDA(nullptr, other.empty, other.region.boundary_size) {
    using std::swap;

    parent_type::len = other.len;
    std::copy(other.dims, other.dims + D, parent_type::dims);
    region = other.region;

    MultiBlockCUDA<D, T> tmp(other.region.len);
    for (iter_type i = 0; i < D; ++i) {
      /*std::copy(other.values[i], other.values[i] + other.region.len,
                tmp.values[i]);*/
      CHECK_CUDA_ERROR(cudaMemcpy(tmp.values[i], other.values[i],
                                  other.region.len * sizeof(T),
                                  cudaMemcpyDeviceToDevice));
    }
    swap(tmp.values, parent_type::values);
  }

  RegionalGridCUDA(RegionalGridCUDA<any_vector_t<T, D>, D>&& other) noexcept
      : RegionalGridCUDA() {
    swap(*this, other);
  }
  RegionalGridCUDA& operator=(RegionalGridCUDA<any_vector_t<T, D>, D> other) {
    swap(*this, other);
    return *this;
  }

  const RegionalGridCUDA<any_vector_t<T, D>, D>& as_grid() const {
    return *this;
  }

  RegionalGridCUDA<any_vector_t<T, D>, D>& as_grid() { return *this; }

  __device__ __host__ const T* axis(Axis ax) const {
    return parent_type::values[symphas::axis_to_index(ax)];
  }

  __device__ __host__ T* axis(Axis ax) {
    return parent_type::values[symphas::axis_to_index(ax)];
  }

  template <typename... Ts,
            std::enable_if_t<
                (sizeof...(Ts) == D &&
                 std::conjunction_v<std::is_convertible<Ts, iter_type>...>),
                int> = 0>
  __device__ __host__ decltype(auto) operator()(Ts&&... rest) const {
    return region(MultiBlockCUDA<D, T>::values, parent_type::dims, empty,
                  (std::forward<Ts>(rest) + region.boundary_size)...);
  }

  __device__ __host__ decltype(auto) operator[](iter_type n) const {
    iter_type pos[D]{};
    grid::get_grid_position(pos, parent_type::dims, n);
    return region(pos, MultiBlockCUDA<D, T>::values, parent_type::dims, empty);
  }

  __device__ __host__ T const& get_unsafe(iter_type n) const {
    return region(n, MultiBlockCUDA<D, T>::values, parent_type::dims, empty);
  }

  __device__ __host__ T& get_unsafe(iter_type n) {
    return region(n, MultiBlockCUDA<D, T>::values, parent_type::dims, empty);
  }

  void adjust(const iter_type (&new_origin)[D]) {
    for (iter_type i = 0; i < D; ++i) {
      grid::adjust_origin_to_from_cuda(
          MultiBlockCUDA<D, T>::values[i], new_origin, region.origin,
          region.dims, parent_type::dims, empty[i], region.boundary_size);
    }
    region.update(new_origin);
  }

  void adjust(const iter_type (&new_origin)[D], const len_type (&new_dims)[D]) {
    for (iter_type i = 0; i < D; ++i) {
      T* new_values;
      CHECK_CUDA_ERROR(
          cudaMalloc(&new_values, grid::length<D>(new_dims) * sizeof(T)));
      int numBlocks = (grid::length<D>(new_dims) + BLOCK_SIZE - 1) / BLOCK_SIZE;
      symphas::cuda::initializeArray CUDA_KERNEL(numBlocks, BLOCK_SIZE)(
          new_values, empty[i], grid::length<D>(new_dims));
      CHECK_CUDA_ERROR(cudaPeekAtLastError());
      CHECK_CUDA_ERROR(cudaDeviceSynchronize());

      grid::adjust_region_to_from_cuda(
          new_values, new_origin, new_dims, MultiBlockCUDA<D, T>::values[i],
          region.origin, region.dims, parent_type::dims, empty[i],
          region.boundary_size);

      using std::swap;
      swap(MultiBlockCUDA<D, T>::values[i], new_values);
      CHECK_CUDA_ERROR(cudaFree(new_values));
    }

    region.update(new_origin, new_dims);
  }

  void adjust(grid::select_region_cuda<D> const& other) {
    if (!grid::is_same(region.dims, other.dims) ||
        !grid::is_same(region.origin, other.origin)) {
      if (grid::is_same(region.dims, other.dims)) {
        adjust(other.origin);
      } else {
        adjust(other.origin, other.dims);
      }
    }
  }

 protected:
  constexpr RegionalGridCUDA() : parent_type(), region{}, empty{} {}
};

namespace grid {
template <typename T, size_t D>
Grid<T, D>& copy(Grid<T, D>& to, GridCUDA<T, D> const& from) {
  CHECK_CUDA_ERROR(cudaMemcpy(to.values, from.values,
                              grid::length<D>(from.dims) * sizeof(T),
                              cudaMemcpyDeviceToHost));
  return to;
}

template <typename T, size_t D>
Grid<any_vector_t<T, D>, D>& copy(Grid<any_vector_t<T, D>, D>& to,
                                  GridCUDA<any_vector_t<T, D>, D> const& from) {
  for (iter_type i = 0; i < D; ++i) {
    CHECK_CUDA_ERROR(cudaMemcpy(to.values[i], from.values[i],
                                grid::length<D>(from.dims) * sizeof(T),
                                cudaMemcpyDeviceToHost));
  }
  return to;
}

template <typename T, size_t D>
GridCUDA<T, D>& copy(GridCUDA<T, D>& to, Grid<T, D> const& from) {
  CHECK_CUDA_ERROR(cudaMemcpy(to.values, from.values,
                              grid::length<D>(from.dims) * sizeof(T),
                              cudaMemcpyHostToDevice));
  return to;
}

template <typename T, size_t D>
GridCUDA<any_vector_t<T, D>, D>& copy(GridCUDA<any_vector_t<T, D>, D>& to,
                                      Grid<any_vector_t<T, D>, D> const& from) {
  for (iter_type i = 0; i < D; ++i) {
    CHECK_CUDA_ERROR(cudaMemcpy(to.values[i], from.values[i],
                                grid::length<D>(from.dims) * sizeof(T),
                                cudaMemcpyHostToDevice));
  }
  return to;
}

template <typename T, size_t D>
BoundaryGrid<T, D>& copy(BoundaryGrid<T, D>& to,
                         BoundaryGridCUDA<T, D> const& from) {
  CHECK_CUDA_ERROR(cudaMemcpy(to.values, from.values,
                              grid::length<D>(from.dims) * sizeof(T),
                              cudaMemcpyDeviceToHost));
  return to;
}

template <typename T, size_t D>
BoundaryGrid<any_vector_t<T, D>, D>& copy(
    BoundaryGrid<any_vector_t<T, D>, D>& to,
    BoundaryGridCUDA<any_vector_t<T, D>, D> const& from) {
  for (iter_type i = 0; i < D; ++i) {
    CHECK_CUDA_ERROR(cudaMemcpy(to.values[i], from.values[i],
                                grid::length<D>(from.dims) * sizeof(T),
                                cudaMemcpyDeviceToHost));
  }
  return to;
}

template <typename T, size_t D>
BoundaryGridCUDA<T, D>& copy(BoundaryGridCUDA<T, D>& to,
                             BoundaryGrid<T, D> const& from) {
  CHECK_CUDA_ERROR(cudaMemcpy(to.values, from.values,
                              grid::length<D>(from.dims) * sizeof(T),
                              cudaMemcpyHostToDevice));
  return to;
}

template <typename T, size_t D>
BoundaryGridCUDA<any_vector_t<T, D>, D>& copy(
    BoundaryGridCUDA<any_vector_t<T, D>, D>& to,
    BoundaryGrid<any_vector_t<T, D>, D> const& from) {
  for (iter_type i = 0; i < D; ++i) {
    CHECK_CUDA_ERROR(cudaMemcpy(to.values[i], from.values[i],
                                grid::length<D>(from.dims) * sizeof(T),
                                cudaMemcpyHostToDevice));
  }
  return to;
}

template <typename T, size_t D>
RegionalGrid<T, D>& copy(RegionalGrid<T, D>& to,
                         RegionalGridCUDA<T, D> const& from) {
  CHECK_CUDA_ERROR(cudaMemcpy(to.values, from.values,
                              grid::length<D>(to.region.dims) * sizeof(T),
                              cudaMemcpyDeviceToHost));
  return to;
}

template <typename T, size_t D>
RegionalGrid<any_vector_t<T, D>, D>& copy(
    RegionalGrid<any_vector_t<T, D>, D>& to,
    RegionalGridCUDA<any_vector_t<T, D>, D> const& from) {
  for (iter_type i = 0; i < D; ++i) {
    CHECK_CUDA_ERROR(cudaMemcpy(to.values[i], from.values[i],
                                grid::length<D>(to.region.dims) * sizeof(T),
                                cudaMemcpyDeviceToHost));
  }
  return to;
}

template <typename T, size_t D>
RegionalGridCUDA<T, D>& copy(RegionalGridCUDA<T, D>& to,
                             RegionalGrid<T, D> const& from) {
  to.adjust(from.region.origin, from.region.dims);
  CHECK_CUDA_ERROR(cudaMemcpy(to.values, from.values,
                              grid::length<D>(from.region.dims) * sizeof(T),
                              cudaMemcpyHostToDevice));
  return to;
}

template <typename T, size_t D>
RegionalGridCUDA<any_vector_t<T, D>, D>& copy(
    RegionalGridCUDA<any_vector_t<T, D>, D>& to,
    RegionalGrid<any_vector_t<T, D>, D> const& from) {
  to.adjust(from.region.origin, from.region.dims);
  for (iter_type i = 0; i < D; ++i) {
    CHECK_CUDA_ERROR(cudaMemcpy(to.values[i], from.values[i],
                                grid::length<D>(to.region.dims) * sizeof(T),
                                cudaMemcpyHostToDevice));
  }
  return to;
}

}  // namespace grid

#endif
