
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
 * MODULE:  sol
 * PURPOSE: Defines a group of boundaries, which is used for a whole
 * system in applying the boundary conditions to a grid.
 *
 * ***************************************************************************
 */

#pragma once

#include "boundaryupdatecuda1d.cuh"
#include "boundaryupdatecuda2d.cuh"
#include "boundaryupdatecuda3dboundary.cuh"
#include "boundaryupdatecuda3ddefault.cuh"
#include "boundaryupdatecuda3dregional.cuh"

//! Manages the boundaries across the whole `D`-dimensional system.
/*!
 * The boundary group manages and controls the update process of the boundaries
 * in a `D`-dimensional system. The boundaries have corresponding dimension
 * `D-1`.
 *
 * \tparam T The system value type.
 * \tparam D The dimension of the system.
 */
template <typename T, size_t D>
struct BoundaryGroup {
  BoundaryType types[D * 2];  //!< The system boundary types.
  grid::Boundary<T, D - 1>*
      boundaries[D * 2];  //!< Boundaries, indexing based on
                          //!< symphas::index_to_side().

  //! Creates boundaries from the boundary and interval data.
  /*!
   * Creates boundaries from the boundary and interval data.
   *
   * \param vdata The interval data used to initialize the boundaries. The
   * range of the boundaries will be taken from this information.
   * \param bdata The information about the system boundaries which is used to
   * initialize each boundary in the group.
   */
  BoundaryGroup(symphas::interval_data_type const& vdata,
                symphas::b_data_type const& bdata);
  BoundaryGroup(BoundaryGroup<T, D> const&);
  BoundaryGroup(BoundaryGroup<T, D>&&) noexcept;
  BoundaryGroup<T, D>& operator=(BoundaryGroup<T, D>);

  template <template <typename, size_t> typename G, size_t... Is>
  void update_boundaries(G<T, D>& grid, iter_type index, double time,
                         std::index_sequence<Is...>);

  //! Update the boundaries of a grid.
  /*!
   * Iterate over the sides and update all the boundaries
   * in this object based on the type of the boundary.
   *
   * \param grid The grid to which the boundaries apply.
   * \param index The current solution index.
   * \param time The current time of the solution data.
   */
  void update_boundaries(Grid<T, D>& grid, iter_type index, double time) {
    update_boundaries(grid, index, time, std::make_index_sequence<D * 2>{});
  }
  //! Update the boundaries of a regional grid.
  /*!
   * Iterate over the sides and update all the boundaries
   * in this object based on the type of the boundary.
   *
   * \param grid The grid to which the boundaries apply.
   * \param index The current solution index.
   * \param time The current time of the solution data.
   */
  void update_boundaries(RegionalGrid<T, D>& grid, iter_type index,
                         double time) {
    update_boundaries(grid, index, time, std::make_index_sequence<D * 2>{});
  }

#ifdef USING_CUDA
  //! Update the boundaries of a grid.
  /*!
   * Iterate over the sides and update all the boundaries
   * in this object based on the type of the boundary.
   *
   * \param grid The grid to which the boundaries apply.
   * \param index The current solution index.
   * \param time The current time of the solution data.
   */
  void update_boundaries(GridCUDA<T, D>& grid, iter_type index, double time) {
    update_boundaries(grid, index, time, std::make_index_sequence<D * 2>{});
  }

  //! Update the boundaries of a regional grid.
  /*!
   * Iterate over the sides and update all the boundaries
   * in this object based on the type of the boundary.
   *
   * \param grid The grid to which the boundaries apply.
   * \param index The current solution index.
   * \param time The current time of the solution data.
   */
  void update_boundaries(RegionalGridCUDA<T, D>& grid, iter_type index,
                         double time) {
    update_boundaries(grid, index, time, std::make_index_sequence<D * 2>{});
  }
#endif

  friend void swap(BoundaryGroup<T, D>& first, BoundaryGroup<T, D>& second) {
    using std::swap;
    swap(first.types, second.types);
    swap(first.boundaries, second.boundaries);
  }

  ~BoundaryGroup();

 protected:
  BoundaryGroup();
};

namespace symphas::internal {

template <size_t D>
void validate_periodic_boundaries(BoundaryType (&types)[D]);

template <size_t N>
bool all_periodic(BoundaryType (&types)[N]) {
  bool flag = true;
  for (iter_type i = 0; i < N; ++i) {
    flag = flag && (types[i] == BoundaryType::PERIODIC);
  }
  return flag;
}

template <>
inline void validate_periodic_boundaries<2>(BoundaryType (&types)[2]) {
  if (types[0] == BoundaryType::PERIODIC &&
      types[1] != BoundaryType::PERIODIC) {
    types[0] = types[1];
  } else if (types[1] == BoundaryType::PERIODIC &&
             types[0] != BoundaryType::PERIODIC) {
    types[1] = types[0];
  }
}

template <>
inline void validate_periodic_boundaries<4>(BoundaryType (&types)[4]) {
  auto [left, right, top, bottom] = std::make_tuple(
      symphas::side_to_index(Side::LEFT), symphas::side_to_index(Side::RIGHT),
      symphas::side_to_index(Side::TOP), symphas::side_to_index(Side::BOTTOM));

  for (auto [side0, side1] :
       {std::make_tuple(left, right), std::make_tuple(top, bottom)}) {
    for (auto [check0, check1] :
         {std::make_tuple(side0, side1), std::make_tuple(side1, side0)}) {
      if (types[check0] == BoundaryType::PERIODIC &&
          types[check1] != BoundaryType::PERIODIC) {
        types[check0] = types[check1];
      }
    }
  }

  if (!all_periodic(types)) {
    if (types[left] == BoundaryType::PERIODIC) {
      types[left] = BoundaryType::PERIODIC0;
      types[right] = BoundaryType::PERIODIC0;
    } else if (types[top] == BoundaryType::PERIODIC) {
      types[top] = BoundaryType::PERIODIC0;
      types[bottom] = BoundaryType::PERIODIC0;
    }
  }
}

template <>
inline void validate_periodic_boundaries<6>(BoundaryType (&types)[6]) {
  auto [left, right, top, bottom, front, back] = std::make_tuple(
      symphas::side_to_index(Side::LEFT), symphas::side_to_index(Side::RIGHT),
      symphas::side_to_index(Side::TOP), symphas::side_to_index(Side::BOTTOM),
      symphas::side_to_index(Side::FRONT), symphas::side_to_index(Side::BACK));

  for (auto [side0, side1] :
       {std::make_tuple(left, right), std::make_tuple(top, bottom),
        std::make_tuple(front, back)}) {
    for (auto [check0, check1] :
         {std::make_tuple(side0, side1), std::make_tuple(side1, side0)}) {
      if (types[check0] == BoundaryType::PERIODIC &&
          types[check1] != BoundaryType::PERIODIC) {
        types[check0] = types[check1];
      }
    }
  }

  if (!all_periodic(types)) {
    if (types[left] == BoundaryType::PERIODIC) {
      BoundaryType select_type =
          (types[top] == BoundaryType::PERIODIC)     ? BoundaryType::PERIODIC3XY
          : (types[front] == BoundaryType::PERIODIC) ? BoundaryType::PERIODIC3XZ
                                                     : BoundaryType::PERIODIC0;
      types[left] = select_type;
      types[right] = select_type;
    } else if (types[top] == BoundaryType::PERIODIC) {
      BoundaryType select_type =
          (types[front] == BoundaryType::PERIODIC)  ? BoundaryType::PERIODIC3XY
          : (types[left] == BoundaryType::PERIODIC) ? BoundaryType::PERIODIC3XZ
                                                    : BoundaryType::PERIODIC0;
      types[top] = select_type;
      types[bottom] = select_type;
    } else if (types[front] == BoundaryType::PERIODIC) {
      BoundaryType select_type =
          (types[top] == BoundaryType::PERIODIC)    ? BoundaryType::PERIODIC3XY
          : (types[left] == BoundaryType::PERIODIC) ? BoundaryType::PERIODIC3XZ
                                                    : BoundaryType::PERIODIC0;
      types[front] = select_type;
      types[back] = select_type;
    }
  }
}

}  // namespace symphas::internal

template <typename T, size_t D>
BoundaryGroup<T, D>::BoundaryGroup(symphas::interval_data_type const& vdata,
                                   symphas::b_data_type const& bdata) {
  for (iter_type i = 0; i < D * 2; ++i) {
    boundaries[i] = nullptr;
  }

  for (auto&& [side, data] : bdata) {
    types[symphas::side_to_index(side)] = data.type;
  }

  symphas::internal::validate_periodic_boundaries(types);

  for (auto&& [side, data] : bdata) {
    boundaries[symphas::side_to_index(side)] =
        symphas::internal::new_boundary<T, D - 1>(data);
  }

  symphas::boundary_setup<D - 1>{}(boundaries, bdata, vdata);
}

template <typename T, size_t D>
BoundaryGroup<T, D>::BoundaryGroup() : types{}, boundaries{} {
  for (iter_type i = 0; i < D * 2; ++i) {
    types[i] = BoundaryType::DEFAULT;
    boundaries[i] = nullptr;
  }
}

template <typename T, size_t D>
BoundaryGroup<T, D>::BoundaryGroup(BoundaryGroup<T, D> const& other) {
  for (iter_type i = 0; i < D * 2; ++i) {
    types[i] = other.types[i];
    if (other.boundaries[i] != nullptr) {
      boundaries[i] = other.boundaries[i]->clone();
    } else {
      boundaries[i] = nullptr;
    }
  }
}

template <typename T, size_t D>
BoundaryGroup<T, D>::BoundaryGroup(BoundaryGroup<T, D>&& other) noexcept
    : BoundaryGroup() {
  swap(*this, other);
}

template <typename T, size_t D>
BoundaryGroup<T, D>& BoundaryGroup<T, D>::operator=(BoundaryGroup<T, D> other) {
  swap(*this, other);
  return *this;
}

template <typename T, size_t D>
BoundaryGroup<T, D>::~BoundaryGroup() {
  for (iter_type i = 0; i < D * 2; ++i) {
    delete boundaries[i];
  }
}

namespace symphas::internal {
template <size_t I>
struct update_boundary_call {
  template <typename B, template <typename, size_t> typename G, typename T,
            size_t D>
  void operator()(B* boundaries, BoundaryType types[D * 2], G<T, D>& grid,
                  iter_type, double time) {
    if (types[I] == BoundaryType::OPEN) {
      throw;
    } else if (types[I] == BoundaryType::DEFAULT) {
      symphas::internal::update_boundary<BoundaryType::DEFAULT,
                                         symphas::index_to_side(I), D - 1>{}(
          boundaries[I], grid, time);
    } else if (types[I] != BoundaryType::NONE) {
      if constexpr (D == 3) {
        switch (types[I]) {
          case BoundaryType::PERIODIC3XZ:
            symphas::internal::update_boundary<
                BoundaryType::PERIODIC3XZ, symphas::index_to_side(I), D - 1>{}(
                boundaries[I], grid);
            break;
          case BoundaryType::PERIODIC3XY:
            symphas::internal::update_boundary<
                BoundaryType::PERIODIC3XY, symphas::index_to_side(I), D - 1>{}(
                boundaries[I], grid);
            break;
          case BoundaryType::PERIODIC0:
            symphas::internal::update_boundary<
                BoundaryType::PERIODIC0, symphas::index_to_side(I), D - 1>{}(
                boundaries[I], grid);
            break;
          case BoundaryType::PERIODIC:
            symphas::internal::update_boundary<
                BoundaryType::PERIODIC, symphas::index_to_side(I), D - 1>{}(
                boundaries[I], grid);
            break;
          default:
            throw;
        }
      } else if constexpr (D == 2) {
        switch (types[I]) {
          case BoundaryType::PERIODIC:
            symphas::internal::update_boundary<
                BoundaryType::PERIODIC, symphas::index_to_side(I), D - 1>{}(
                boundaries[I], grid);
            break;
          case BoundaryType::PERIODIC0:
            symphas::internal::update_boundary<
                BoundaryType::PERIODIC0, symphas::index_to_side(I), D - 1>{}(
                boundaries[I], grid);
            break;
          default:
            throw;
        }
      } else {
        symphas::internal::update_boundary<BoundaryType::PERIODIC,
                                           symphas::index_to_side(I), D - 1>{}(
            boundaries[I], grid);
      }
    }
  }
};

}  // namespace symphas::internal

template <typename T, size_t D>
template <template <typename, size_t> typename G, size_t... Is>
void BoundaryGroup<T, D>::update_boundaries(G<T, D>& grid, iter_type index,
                                            double time,
                                            std::index_sequence<Is...>) {
  (symphas::internal::update_boundary_call<Is>{}(boundaries, types, grid, index,
                                                 time),
   ...);
}
