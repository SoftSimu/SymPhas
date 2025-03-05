
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
 * PURPOSE: Defines elements that can allow definition of additional symbols
 * that can be used in the symbolic algebra framework, as well as type traits
 * that define rules and factoring of expressions.
 *
 * ***************************************************************************
 */

#pragma once

#include "expressionlogic.h"
#include "expressionproperties.cuh"

// **************************************************************************************

//! Helper for accessing data encapsulated by different types.
/*!
 * For the different types of objects that be data to an OpTerm, a
 * method is implemented to correctly obtain the value (the data point) at a
 * given index. Additionally, the entire data itself can be obtained.
 */
template <typename A>
template <typename T>
const auto& expr::BaseData<A>::get(BlockCUDA<T> const& data, iter_type n) {
  return data[n];
}

template <typename A>
template <typename T>
const auto& expr::BaseData<A>::get(BlockCUDA<T> const& data) {
  return data;
}

template <typename A>
template <typename T>
auto& expr::BaseData<A>::get(BlockCUDA<T>& data, iter_type n) {
  return data[n];
}

template <typename A>
template <typename T>
auto& expr::BaseData<A>::get(BlockCUDA<T>& data) {
  return data;
}

//! Specialization based on expr::BaseData.
template <typename T>
struct expr::BaseData<BlockCUDA<T>> {
  static auto const& get(BlockCUDA<T> const& data, iter_type n) {
    return data[n];
  }
  static auto const& get(BlockCUDA<T> const& data) { return data; }
  static auto& get(BlockCUDA<T>& data, iter_type n) { return data[n]; }
  static auto& get(BlockCUDA<T>& data) { return data; }
};

//! Specialization based on expr::BaseData.
template <size_t N, typename T>
struct expr::BaseData<MultiBlockCUDA<N, T>> {
  static decltype(auto) get(MultiBlockCUDA<N, T> const& data, iter_type n) {
    return data[n];
  }
  static auto const& get(MultiBlockCUDA<N, T> const& data) { return data; }
  static decltype(auto) get(MultiBlockCUDA<N, T>& data, iter_type n) {
    return data[n];
  }
  static auto& get(MultiBlockCUDA<N, T>& data) { return data; }
};

//! Specialization based on expr::BaseData.
template <typename T, size_t N>
struct expr::BaseData<GridData<MultiBlockCUDA<N, T>, N>> {
  static decltype(auto) get(GridData<MultiBlockCUDA<N, T>, N> const& data,
                            iter_type n) {
    return BaseData<MultiBlockCUDA<N, T>>::get(data, n);
  }
  static decltype(auto) get(GridData<MultiBlockCUDA<N, T>, N> const& data) {
    return BaseData<MultiBlockCUDA<N, T>>::get(data);
  }
  static decltype(auto) get(GridData<MultiBlockCUDA<N, T>, N>& data,
                            iter_type n) {
    return BaseData<MultiBlockCUDA<N, T>>::get(data, n);
  }
  static decltype(auto) get(GridData<MultiBlockCUDA<N, T>, N>& data) {
    return BaseData<MultiBlockCUDA<N, T>>::get(data);
  }
};

template <typename G>
template <size_t N, typename T>
decltype(auto) expr::construct_result_data<G>::get_constructor(
    MultiBlockCUDA<N, T>) const {
  return construct_result_data<MultiBlockCUDA<N, T>>{};
}
template <typename G>
template <typename T, size_t D>
decltype(auto) expr::construct_result_data<G>::get_constructor(
    GridCUDA<T, D>) const {
  return construct_result_data<GridCUDA<T, D>>{};
}
template <typename G>
template <typename T, size_t D>
decltype(auto) expr::construct_result_data<G>::get_constructor(
    GridCUDA<any_vector_t<T, D>, D>) const {
  return construct_result_data<MultiBlockCUDA<D, T>>{};
}
template <typename G>
template <typename T>
decltype(auto) expr::construct_result_data<G>::get_constructor(
    BlockCUDA<T>) const {
  return construct_result_data<BlockCUDA<T>>{};
}

DEFINE_BASE_DATA((typename T, size_t D), (GridCUDA<T, D>), (T)data[n], data)
DEFINE_BASE_DATA((typename T, size_t D), (BoundaryGridCUDA<T, D>), (T)data[n],
                 data)
DEFINE_BASE_DATA((typename T, size_t D), (RegionalGridCUDA<T, D>), (T)data[n],
                 data)