
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
 * PURPOSE: Defines functions for basic functionality throughout SymPhas.
 * This includes functions which perform common tasks.
 *
 * ***************************************************************************
 */

#pragma once

#include <algorithm>
#include <cstring>
#include <numeric>

#include "spslibseq.h"

namespace symphas {
namespace lib {

//! Collect and count like types from a tuple list.
/*!
 * Give the list of types, the list is made unique and number of times
 * each type appears in the list is repeated.
 */
template <typename... Gs>
struct filter_types_impl;

template <typename... G0s>
struct filter_types_impl<types_list<G0s...>> {
  using type = types_list<G0s...>;
};

template <typename... G0s>
struct filter_types_impl<types_list<G0s...>, types_list<>> {
  using type = types_list<G0s...>;
};

template <bool... flags, typename... G0s>
struct filter_types_impl<types_list<G0s...>,
                         std::integer_sequence<bool, flags...>> {
  using type = expand_types_list<
      std::conditional_t<flags, types_list<G0s>, types_list<>>...>;
};

template <typename... G0s, typename G1, typename... G1s>
struct filter_types_impl<types_list<G0s...>, types_list<G1, G1s...>> {
 protected:
  using filter_mask =
      std::integer_sequence<bool, !std::is_same<G1, G0s>::value...>;
  using filtered_1 =
      typename filter_types_impl<types_list<G0s...>, filter_mask>::type;

 public:
  using type = typename filter_types_impl<filtered_1, types_list<G1s...>>::type;
};

template <typename... G0s, typename... G1s, typename... Rest>
struct filter_types_impl<types_list<G0s...>, types_list<G1s...>, Rest...> {
  using type = typename filter_types_impl<
      typename filter_types_impl<types_list<G0s...>, types_list<G1s...>>::type,
      Rest...>::type;
};

//! Filter types from the first type list.
/*!
 * Filter from the first list everything that's in the subsequent lists. That
 * is, any element in the lists after the first will be removed from the first.
 *
 * If an std::integer_sequence<bool> type is passed as the second argument, it
 * will act as a mask.
 */
template <typename... Gs>
using filter_types = typename filter_types_impl<Gs...>::type;

template <typename Seq, typename... Ts>
struct filter_types_on_index_impl;

template <size_t... Is, typename... Ts>
struct filter_types_on_index_impl<std::index_sequence<Is...>, Ts...> {
  using type = select_types<
      symphas::lib::filter_seq_t<std::make_index_sequence<sizeof...(Ts)>,
                                 std::index_sequence<Is...>>,
      Ts...>;
};

template <size_t... Is, typename... Ts>
struct filter_types_on_index_impl<std::index_sequence<Is...>,
                                  unroll_types_list<types_list<Ts...>>> {
  using type = select_types<
      symphas::lib::filter_seq_t<std::make_index_sequence<sizeof...(Ts)>,
                                 std::index_sequence<Is...>>,
      Ts...>;
};

template <typename Seq, typename... Ts>
using filter_types_on_index =
    typename filter_types_on_index_impl<Seq, Ts...>::type;

//! Return the index of the given type in the type list.
/*!
 * Given the type list, identify the index of the `I`-th chosen type. If
 * the type is not identified, then `-1` will be returned.
 *
 * \tparam Type The type to be found in the list.
 * \tparam I The index of the type to find out of those that match.
 * \tparam Ss... The types in the list.
 */
template <typename Type, size_t I, typename Seq, typename... Ss>
struct index_of_type_impl;

template <typename Type, size_t I, size_t... Is, typename... Ss>
struct index_of_type_impl<Type, I, std::index_sequence<Is...>, Ss...> {
  using mask_t =
      seq_join_t<std::integer_sequence<int>,
                 std::conditional_t<(std::is_same_v<Ss, Type>),
                                    std::integer_sequence<int, int(Is)>,
                                    std::integer_sequence<int>>...>;

  static constexpr int value =
      (I < mask_t::size())
          ? seq_index_value_any_default<I, mask_t,
                                        std::integer_sequence<int, -1>>::value
          : -1l;
};

//! Return the index of the given type in the type list.
/*!
 * Given the type list, identify the index of the `I`-th chosen type. If
 * the type is not identified, then `-1` will be returned.
 *
 * \tparam Type The type to be found in the list.
 * \tparam I The index of the type to find out of those that match.
 * \tparam S0 The next type in the type list to compare.
 * \tparam Ss... The remaining types in the list.
 */
template <typename Type, size_t I, typename... Ss>
constexpr int nth_index_of_type =
    index_of_type_impl<Type, I, std::make_index_sequence<sizeof...(Ss)>,
                       Ss...>::value;

//! Return the index of the given type in the type list.
/*!
 * Given the type list, identify the index of the `I`-th chosen type. If
 * the type is not identified, then `-1` will be returned.
 *
 * \tparam Type The type to be found in the list.
 * \tparam S0 The next type in the type list to compare.
 * \tparam Ss... The remaining types in the list.
 */
template <typename Type, typename... Ss>
constexpr int index_of_type = nth_index_of_type<Type, 0, Ss...>;

template <typename Type>
constexpr int index_of_type<Type> = -1;

template <typename Type, Type N, Type... Is>
constexpr int index_of_value =
    index_of_type<std::integer_sequence<Type, N>,
                  std::integer_sequence<Type, Is>...>;

template <typename Type, size_t N, typename S0, typename... Ss>
constexpr int nth_index_of_type<Type, N, types_list<S0, Ss...>> =
    nth_index_of_type<Type, N, S0, Ss...>;

template <typename Type, typename S0, typename... Ss>
constexpr int index_of_type<Type, types_list<S0, Ss...>> =
    index_of_type<Type, S0, Ss...>;

template <typename Type, typename S0, typename... Ss>
constexpr int index_of_type<Type, unroll_types_list<types_list<S0, Ss...>>> =
    index_of_type<Type, S0, Ss...>;

template <typename... Gs>
struct combine_collect_like_types;

template <>
struct combine_collect_like_types<> {
  using type = types_list<>;
  using count_seq = std::index_sequence<>;
};

template <typename G0, typename... Gs>
struct combine_collect_like_types<G0, Gs...> {
 protected:
  using filtered =
      symphas::lib::filter_types<types_list<G0, Gs...>, types_list<G0>>;
  using like_types_rest = combine_collect_like_types<filtered>;
  static const size_t count =
      (((std::is_same<G0, Gs>::value) ? 1 : 0) + ... + 1);

 public:
  using type = expand_types_list<G0, typename like_types_rest::type>;
  using count_seq = seq_join_t<std::index_sequence<count>,
                               typename like_types_rest::count_seq>;
};

template <typename... Gs>
struct combine_collect_like_types<types_list<Gs...>> {
  using type = typename combine_collect_like_types<Gs...>::type;
  using count_seq = typename combine_collect_like_types<Gs...>::count_seq;
};

//! Returns the unique list of types from the list.
/*!
 * Returns the unique list of types is returned from joining together
 * the list of types.
 */
template <typename... Gs>
struct combine_types_unique {
  using type =
      typename combine_collect_like_types<expand_types_list<Gs...>>::type;
};

/*!
 * @brief A container for a dynamically allocated array.
 *
 * This struct represents a container for a dynamically allocated array of type
 * `T`. It manages the memory for the array, including allocation and
 * deallocation. It also provides several utility functions for working with the
 * array, such as accessing elements, getting the size, and expanding the array.
 *
 * @tparam T The type of the elements in the array.
 */
template <typename T>
struct array_container {
  T* data;
  size_t n;

  __host__ __device__ array_container() : array_container(0) {}
  __host__ __device__ array_container(size_t n)
      : data{(n > 0) ? new T[n]{} : nullptr}, n{n} {}
  __host__ __device__ array_container(const T* data, size_t n)
      : array_container(n) {
    if (data) {
      for (iter_type i = 0; i < n; ++i) {
        this->data[i] = data[i];
      }
    }
  }

  template <
      typename... Ts,
      std::enable_if_t<std::conjunction_v<std::is_same<Ts, T>...>, int> = 0>
  __host__ __device__ array_container(T data0, Ts... datas)
      : data{new T[sizeof...(Ts) + 1]{data0, static_cast<T>(datas)...}},
        n{sizeof...(Ts) + 1} {}

  __host__ __device__ array_container(array_container<T> const& other) noexcept
      : array_container(other.data, other.n) {}

  __host__ __device__ array_container(array_container<T>&& other) noexcept
      : array_container() {
    swap(*this, other);
  }

  __host__ __device__ array_container& operator=(array_container other) {
    swap(*this, other);
    return *this;
  }

  __host__ __device__ ~array_container() { delete[] data; }

  template <typename... Ts,
            std::enable_if_t<std::conjunction_v<std::is_convertible<Ts, T>...>,
                             int> = 0>
  static array_container<T> create(Ts&&... ts) {
    return array_container<T>(static_cast<T>(std::forward<Ts>(ts))...);
  }

  __host__ __device__ operator const T*() const { return data; }

  __host__ __device__ operator T*() { return data; }

  __host__ __device__ T& operator[](iter_type i) { return data[i]; }

  __host__ __device__ const T& operator[](iter_type i) const { return data[i]; }

  __host__ __device__ friend void swap(array_container<T>& first,
                                       array_container<T>& second) {
    auto data = first.data;
    first.data = second.data;
    second.data = data;

    auto n = first.n;
    first.n = second.n;
    second.n = n;
  }

  __host__ __device__ const T* begin() const { return data; }

  __host__ __device__ const T* end() const { return data + n; }

  __host__ __device__ T* begin() { return data; }

  __host__ __device__ T* end() { return data + n; }

  __host__ __device__ array_container<T> operator*(size_t N) {
    array_container<T> expanded(n * N);
    for (iter_type i = 0; i < N * n; ++i) {
      expanded[i] = data[i % n];
    }
    return expanded;
  }

  __host__ __device__ array_container<T>& operator*=(size_t N) {
    array_container<T> expanded(n * N);
    for (iter_type i = 0; i < N * n; ++i) {
      expanded[i] = data[i % n];
    }
    swap(*this, expanded);
    return *this;
  }
};

template <typename T, size_t N>
struct fixed_array_container {
  T data[N];

  __host__ __device__ fixed_array_container() : data{} {}

  __host__ __device__ fixed_array_container(const T* data, size_t n) : data{} {
    if (data) {
      for (iter_type i = 0; i < n; ++i) {
        this->data[i] = data[i];
      }
    }
  }

  template <
      typename... Ts,
      std::enable_if_t<std::conjunction_v<std::is_same<Ts, T>...>, int> = 0>
  __host__ __device__ fixed_array_container(T data0, Ts... datas)
      : data{data0, static_cast<T>(datas)...} {}

  __host__ __device__ operator const T*() const { return data; }

  __host__ __device__ operator T*() { return data; }

  __host__ __device__ T& operator[](iter_type i) { return data[i]; }

  __host__ __device__ const T& operator[](iter_type i) const { return data[i]; }

  __host__ __device__ const T* begin() const { return data; }

  __host__ __device__ const T* end() const { return data + N; }

  __host__ __device__ T* begin() { return data; }

  __host__ __device__ T* end() { return data + N; }
};

//! Generalized function to expand the given list.
/*!
 * All the elements from the input are added to the out array, which is
 * appropriately resized in order to accommodate the appended elements.
 * Since the array is resized, the array needs to be taken by reference to
 * be properly reassigned.
 *
 * \param[in] in The buffer list which will be copied to out.
 * \param[out] out The elements from the buffer parameter `in` are added to
 * the end of this list and its memory is reallocated to fit the new size.
 * \param in_len The total length of the input, which is how many elements
 * will be appended to the output.
 * \param out_len The current length of out, so that the memory can be
 * appropriately reallocated with a larger length.
 */
template <typename T>
void expand_append_array(const T* in, T*& out, size_t in_len, size_t out_len) {
  size_t total_len = out_len + in_len;
  T* out_buffer = new T[total_len];
  std::copy(out, out + out_len, out_buffer);
  std::copy(in, in + in_len, out_buffer + out_len);

  delete[] out;
  out = out_buffer;
}

//! Generalized function to expand the given list.
/*!
 * A given number of copies of the provided input value are added to the out
 * array, which is appropriately resized in order to accommodate the
 * appended elements. Since the array is resized, the array needs to be
 * taken by reference to be properly reassigned.
 *
 * \param[in] in The element of which copies will be appended to the output.
 * \param[out] out The elements from the buffer parameter `in` are added to
 * the end of this list and its memory is reallocated to fit the new size.
 * \param in_len The total length of the input, which is how many elements
 * will be appended to the output.
 * \param out_len The current length of out, so that the memory can be
 * appropriately reallocated with a larger length.
 */
template <typename T>
void expand_append_array(T in, T*& out, size_t in_len, size_t out_len = 0) {
  size_t total_len = out_len + in_len;
  T* out_buffer = new T[total_len];
  if (out_len > 0) {
    std::copy(out, out + out_len, out_buffer);
  }
  std::fill(out_buffer + out_len, out_buffer + total_len, in);

  delete[] out;
  out = out_buffer;
}

template <typename array_type>
void clear_array(array_type*& arr, size_t const& len) {
  delete[] arr;
  //for (iter_type i = 0; i < len; ++i) {
  //  arr[i].~array_type();
  //}
  //operator delete[](arr);
}

template <typename array_type>
void clear_array(array_type**& arr, size_t const& len, size_t element_len) {
  for (iter_type i = 0; i < len; ++i) {
    for (iter_type n = 0; n < element_len; ++n) {
      arr[i][n].~array_type();
    }
    operator delete[](arr[i]);
  }
  delete[] arr;
}

inline void clear_array(char**& arr, size_t const& len) {
  for (iter_type i = 0; i < len; ++i) {
    delete[] arr[i];
  }
  delete[] arr;
}

/*!
 * @brief Resizes an array to a new length.
 *
 * This function resizes an array to a new length. If the new length is greater
 * than the current length, new elements are added to the end of the array and
 * initialized with a provided default value. If the new length is less than the
 * current length, elements at the end of the array are discarded. The array is
 * resized in-place, and the original array is deleted.
 *
 * @tparam array_type The type of the elements in the array.
 * @param new_len The new length for the array.
 * @param arr A pointer to the array to be resized. This should be a pointer to
 * a dynamically-allocated array.
 * @param len The current length of the array.
 * @param copy The value to initialize new elements with, if the array is being
 * expanded. Defaults to a default-constructed instance of array_type.
 */
template <typename array_type>
void resize_array_impl(size_t new_len, array_type*& arr, size_t const& len,
                       array_type const& copy = array_type{}) {
  if (new_len != len) {
    array_type* ptr = new array_type[new_len]{};
    //void* extend = operator new[](new_len * sizeof(array_type));
    //array_type* ptr = static_cast<array_type*>(extend);

    for (iter_type i = 0; i < std::min(len, new_len); ++i) {
      using std::swap;
      swap(ptr[i], arr[i]);
    }

    for (iter_type i = (iter_type)len; i < new_len; ++i) {
      array_type copy2(copy);
      swap(ptr[i], copy2);
    }

    std::swap(ptr, arr);
    clear_array(ptr, len);
  }
}

template <typename array_type>
void resize_array(size_t new_len, array_type*& arr, size_t const& len,
                  array_type const& copy = array_type{}) {
  resize_array_impl(new_len, arr, len, copy);
}

/*!
 * @brief Resizes an array to a new length.
 *
 * This function resizes an array to a new length. If the new length is greater
 * than the current length, new elements are added to the end of the array and
 * initialized with a provided default value. If the new length is less than the
 * current length, elements at the end of the array are discarded. The array is
 * resized in-place, and the original array is deleted.
 *
 * @tparam array_type The type of the elements in the array.
 * @param new_len The new length for the array.
 * @param arr A pointer to the array to be resized. This should be a pointer to
 * a dynamically-allocated array.
 * @param len The current length of the array. This will be updated to the new
 * length.
 * @param copy The value to initialize new elements with, if the array is being
 * expanded. Defaults to a default-constructed instance of array_type.
 */
template <typename array_type>
void resize_array(size_t new_len, array_type*& arr, size_t& len,
                  array_type const& copy = array_type{}) {
  resize_array_impl(new_len, arr, len, copy);
  len = new_len;
}

/*!
 * \brief Resizes a 2D array to a new length.
 *
 * This function resizes a 2D array to a new length. If the new length is
 * greater than the current length, new elements are added to the end of the
 * array. If the new length is less than the current length, elements at the end
 * of the array are discarded. The array is resized in-place, and the original
 * array is deleted.
 *
 * \tparam array_type The type of the elements in the array.
 * \param new_len The new length for the array.
 * \param arr A pointer to the array to be resized. This should be a pointer to
 * a dynamically-allocated array.
 * \param len The current length of the array.
 * \param element_len The length of each element in the 2D array.
 */
template <typename array_type>
void resize_array_impl(size_t new_len, array_type**& arr, size_t const& len,
                       size_t element_len) {
  if (new_len != len) {
    len_type** extend = new array_type* [new_len] {};

    for (iter_type i = 0; i < std::min(len, new_len); ++i) {
      extend[i] = nullptr;
      std::swap(arr[i], extend[i]);
      void* element = operator new[](element_len * sizeof(array_type));
    }
    for (iter_type i = iter_type(len); i < new_len; ++i) {
      void* element = operator new[](element_len * sizeof(array_type));
      array_type* ptr = static_cast<array_type*>(element);
      for (iter_type n = 0; n < element_len; ++n) {
        new (&ptr[n]) array_type{};
      }
      extend[i] = ptr;
    }
    std::swap(extend, arr);
    clear_array(extend, len, element_len);
  }
}

/*!
 * \brief Resizes a 2D array to a new length.
 *
 * This function resizes a 2D array to a new length. If the new length is
 * greater than the current length, new elements are added to the end of the
 * array. If the new length is less than the current length, elements at the end
 * of the array are discarded. The array is resized in-place, and the original
 * array is deleted.
 *
 * \tparam array_type The type of the elements in the array.
 * \param new_len The new length for the array.
 * \param arr A pointer to the array to be resized. This should be a pointer to
 * a dynamically-allocated array.
 * \param len The current length of the array.
 * \param element_len The length of each element in the 2D array.
 */
template <typename array_type>
void resize_array(size_t new_len, array_type**& arr, size_t const& len,
                  size_t element_len) {
  resize_array_impl(new_len, arr, len, element_len);
}

/*!
 * @brief Resizes a 2D array to a new length.
 *
 * This function resizes a 2D array to a new length. If the new length is
 * greater than the current length, new elements are added to the end of the
 * array. If the new length is less than the current length, elements at the
 * end of the array are discarded. The array is resized in-place, and the
 * original array is deleted.
 *
 * @tparam array_type The type of the elements in the array.
 * @param new_len The new length for the array.
 * @param arr A pointer to the array to be resized. This should be a pointer
 * to a dynamically-allocated array.
 * @param len The current length of the array. This will be updated to the new
 * length.
 * @param element_len The length of each element in the 2D array.
 */
template <typename array_type>
void resize_array(size_t new_len, array_type**& arr, size_t& len,
                  size_t element_len) {
  resize_array_impl(new_len, arr, len, element_len);
  len = new_len;
}

/*!
 * @brief Resizes a 2D array of char pointers to a new length.
 *
 * This function resizes a 2D array of char pointers to a new length. If the new
 * length is greater than the current length, new elements are added to the end
 * of the array and initialized to nullptr. If the new length is less than the
 * current length, elements at the end of the array are discarded. The array is
 * resized in-place, and the original array is deleted.
 *
 * @param new_len The new length for the array.
 * @param arr A pointer to the array to be resized. This should be a pointer
 * to a dynamically-allocated array.
 * @param len The current length of the array.
 */

inline void resize_array_impl(size_t new_len, char**& arr, size_t const& len) {
  if (new_len != len) {
    char** extend = new char* [new_len] {};

    for (iter_type i = 0; i < std::min(len, new_len); ++i) {
      extend[i] = nullptr;
      std::swap(arr[i], extend[i]);
    }

    for (iter_type i = iter_type(len); i < new_len; ++i) {
      extend[i] = nullptr;
    }
    std::swap(extend, arr);
    clear_array(extend, len);
  }
}

inline void resize_array(size_t new_len, char**& arr, size_t const& len) {
  resize_array_impl(new_len, arr, len);
}

inline void resize_array(size_t new_len, char**& arr, size_t& len) {
  resize_array_impl(new_len, arr, len);
  len = new_len;
}

template <typename value_type, typename array_type>
void set_array_value(value_type const& value, size_t n, array_type*& arr,
                     size_t& len) {
  if (n >= len) {
    resize_array(n + 1, arr, len);
  }

  arr[n] = value;
}

template <>
inline void set_array_value<>(const char* const& value, size_t n, char**& arr,
                              size_t& len) {
  if (n >= len) {
    resize_array(n + 1, arr, len);
  }

  delete[] arr[n];
  arr[n] = new char[std::strlen(value) + 1]{};
  std::strcpy(arr[n], value);
}

}  // namespace lib
}  // namespace symphas