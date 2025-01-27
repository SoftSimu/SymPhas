
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
#include <cmath>
#include <cstring>
#include <iostream>
#include <numeric>
#include <tuple>
#include <vector>

namespace symphas {
namespace lib {

template <typename data_type>
struct basic_forward_iterator_container;

template <typename data_type>
struct basic_forward_iterator {
  using iterator_category = std::forward_iterator_tag;
  using difference_type = int;
  using value_type = basic_forward_iterator_container<data_type>;
  using pointer = value_type*;
  using reference = int;

  template <typename T>
  basic_forward_iterator(T&& data, iter_type pos = 0)
      : ptr{std::forward<T>(data), pos} {}

  basic_forward_iterator_container<data_type> operator*() const { return ptr; }

  // Prefix increment
  basic_forward_iterator& operator++() {
    ++ptr.pos;
    return *this;
  }

  // Postfix increment
  basic_forward_iterator operator++(int) {
    basic_forward_iterator tmp = *this;
    ++(*this);
    return tmp;
  }

  friend bool operator==(basic_forward_iterator const& a,
                         basic_forward_iterator const& b) {
    return a.ptr.pos == b.ptr.pos;
  }

  friend bool operator!=(basic_forward_iterator const& a,
                         basic_forward_iterator const& b) {
    return !(a == b);
  }

  basic_forward_iterator_container<data_type> ptr;
};

template <typename... Ts>
struct zip_iterator {
  using iterator_category = std::forward_iterator_tag;
  using difference_type = int;
  using value_type = std::tuple<Ts...>;
  using pointer = value_type*;
  using reference = int;

  using seq = std::make_index_sequence<sizeof...(Ts)>;

  zip_iterator(Ts const&... ts)
      : data{std::begin(ts)...},
        n{0}  //, end{ std::end(ts)... }
  {}

  decltype(auto) operator*() const { return get_data(seq{}); }

  decltype(auto) operator[](iter_type offset) const {
    return get_data(seq{}, offset);
  }

  // Prefix increment
  zip_iterator<Ts...>& operator++() {
    n += 1;
    return *this;
  }

  // Postfix increment
  zip_iterator<Ts...> operator++(int) {
    zip_iterator<Ts...> tmp = *this;
    ++(*this);
    return tmp;
  }

  friend bool operator==(zip_iterator<Ts...> const& a,
                         zip_iterator<Ts...> const& b) {
    return a.check_all_same(b, seq{});
  }

  friend bool operator!=(zip_iterator<Ts...> const& a,
                         zip_iterator<Ts...> const& b) {
    return a.check_different(b, seq{});
  }

  template <size_t... Is>
  bool check_all_same(zip_iterator<Ts...> const& b,
                      std::index_sequence<Is...>) const {
    return (check_one_same<Is>(b) && ...);
  }

  template <size_t... Is>
  bool check_different(zip_iterator<Ts...> const& b,
                       std::index_sequence<Is...>) const {
    return (!check_one_same<Is>(b) || ...);
  }

  template <size_t I>
  bool check_one_same(zip_iterator<Ts...> const& b) const {
    if (n == b.n) {
      return std::get<I>(data) == std::get<I>(b.data);
    } else {
      return std::get<I>(b.data) - std::get<I>(data) == b.n - n;
    }
  }

  template <size_t... Is>
  decltype(auto) get_data(std::index_sequence<Is...>,
                          iter_type offset = 0) const {
    return std::tie(*(std::get<Is>(data) + std::ptrdiff_t(n + offset))...);
  }

  template <size_t... Is>
  decltype(auto) get_data(std::index_sequence<Is...>, iter_type offset = 0) {
    return std::tie(std::get<Is>(data)[n + offset]...);
  }

 protected:
  template <typename T>
  static decltype(auto) begin(T const& t) {
    return std::begin(t);
  }

 public:
  std::tuple<
      std::invoke_result_t<decltype(&zip_iterator<Ts...>::begin<Ts>), Ts>...>
      data;
  iter_type n;
};

template <typename T0, typename... Ts>
struct zip_container {
  zip_container(T0 const& t0, Ts const&... ts)
      : iter{t0, ts...},
        len{static_cast<len_type>(std::end(t0) - std::begin(t0))} {}

  auto& begin() { return iter; }

  auto end() {
    zip_iterator end(iter);
    end.n = len;
    return end;
  }

  auto& begin() const { return iter; }

  auto end() const {
    zip_iterator end(iter);
    end.n = len;
    return end;
  }

  zip_iterator<T0, Ts...> iter;
  len_type len;
};

template <typename T0, typename... Ts>
struct zip_container<T0&&, Ts&&...> {
  zip_container(T0&& t0, Ts&&... ts)
      : data{t0, ts...},
        iter{get_iter()},
        len{static_cast<len_type>(std::end(std::get<0>(data)) -
                                  std::begin(std::get<0>(data)))} {}

  auto begin() { return iter; }

  auto end() {
    zip_iterator<T0, Ts...> end(iter);
    end.n = len;
    return end;
  }

  auto begin() const { return iter; }

  auto end() const {
    zip_iterator end(iter);
    end.n = len;
    return end;
  }

  std::tuple<T0, Ts...> data;
  zip_iterator<T0, Ts...> iter;
  len_type len;

 protected:
  auto get_iter() const {
    return get_iter(std::make_index_sequence<1 + sizeof...(Ts)>{});
  }

  template <size_t... Is>
  auto get_iter(std::index_sequence<Is...>) const {
    return zip_iterator(std::get<Is>(data)...);
  }
};

template <typename... Ts>
zip_container(Ts const&...) -> zip_container<Ts...>;

template <typename... Ts>
zip_container(Ts&&...) -> zip_container<Ts&&...>;

template <typename... Ts>
auto zip(Ts&&... ts) {
  return zip_container(std::forward<Ts>(ts)...);
}
}  // namespace lib
}  // namespace symphas
