
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
 * PURPOSE: Elements related to defining a phase field system, including
 * data persistence and representing information. Also contains methods of
 * getting the correct grid data based on different grid implementations.
 *
 * ***************************************************************************
 */

#pragma once

#include "boundary.cuh"
#include "systemlib.h"

#ifdef USING_CUDA

//! Maintains grid data and associated information.
/*!
 * Maintains grid data and parameters, and defines the ability to use data
 * persistence by saving a snapshot. Specialization implementing
 * a system with boundaries.
 *
 * See SystemData<Grid<T, D>>.
 */
template <typename T, size_t D>
struct SystemData<GridCUDA<T, D>> : GridCUDA<T, D>, SystemInfo {
  using SystemInfo::id;
  using SystemInfo::info;

  using GridCUDA<T, D>::dims;
  using GridCUDA<T, D>::len;
  using GridCUDA<T, D>::as_grid;

  //! Create a system data instance.
  /*!
   * Data for a phase field system is generated, resulting in a new
   * grid instance and information defining that grid.
   */
  SystemData(symphas::interval_data_type const& vdata, size_t id)
      : GridCUDA<T, D>{grid::construct<::GridCUDA, T, D>(vdata)},
        SystemInfo{{vdata}, id} {
    if (info.dimension() > 0) {
      for (iter_type i = info.dimension() - 1; i >= D; --i) {
        info.intervals.erase(symphas::index_to_axis(i));
      }
    }
  }

  //! Get the grid representation of this system.
  /*!
   * Get the grid representation of this system, which will return the
   * underlying grid with boundaries.
   */
  operator const GridCUDA<T, D>&() const { return as_grid(); }

  //! Get the grid representation of this system.
  /*!
   * Get the grid representation of this system, which will return the
   * underlying grid with boundaries.
   */
  operator GridCUDA<T, D>&() { return as_grid(); }

  //! Get the number of elements of true data.
  /*!
   * The elements which are considered true data are the interior domain,
   * that is, all non-boundary elements. Return the count of these elements.
   */
  auto length() const { return len; }

  //! Copies the system data into the given array.
  /*!
   * The values of the system data block are copied into a new one. The copy
   * is performed point-wise for all interior points, meaning points not on
   * the boundary.
   */
  void persist(T* out) const {
    CHECK_CUDA_ERROR(cudaMemcpy(out, GridCUDA<T, D>::values, len * sizeof(T),
                                cudaMemcpyDeviceToHost));
  }

  //! Copies the input data into the system values.
  /*!
   * The values of the system data block are initialized from the
   * given values, correctly transcribing all values.
   */
  void fill(const T* in) const {
    CHECK_CUDA_ERROR(cudaMemcpy(GridCUDA<T, D>::values, in, len * sizeof(T),
                                cudaMemcpyHostToDevice));
  }

 protected:
  SystemData() : GridCUDA<T, D>{}, SystemInfo{} {}
};

//! Maintains grid data and associated information.
/*!
 * Maintains grid data and parameters, and defines the ability to use data
 * persistence by saving a snapshot. Specialization implementing
 * a system with boundaries.
 *
 * See SystemData<Grid<T, D>>.
 */
template <typename T, size_t D>
struct SystemData<BoundaryGridCUDA<T, D>> : BoundaryGridCUDA<T, D>, SystemInfo {
  using SystemInfo::id;
  using SystemInfo::info;

  using GridCUDA<T, D>::dims;
  using GridCUDA<T, D>::len;
  using BoundaryGridCUDA<T, D>::as_grid;

  //! Create a system data instance.
  /*!
   * Data for a phase field system is generated, resulting in a new
   * grid instance and information defining that grid.
   */
  SystemData(symphas::interval_data_type const& vdata, size_t id)
      : BoundaryGridCUDA<T, D>{grid::construct<::BoundaryGridCUDA, T, D>(
            vdata)},
        SystemInfo{{vdata}, id} {
    if (!info.intervals.empty()) {
      for (iter_type i = 0; i < D; ++i) {
        Axis ax = symphas::index_to_axis(i);
        auto& interval = info.intervals.at(ax);

        interval.set_count(dims[i] - 2 * BOUNDARY_DEPTH);
        interval.interval_to_domain();
      }
    }
  }

  //! Get the number of elements of true data.
  /*!
   * The elements which are considered true data are the interior domain,
   * that is, all non-boundary elements. Return the count of these elements.
   */
  auto length() const { return grid::length_interior<D>(dims); }

  //! Copies the system data into the given array.
  /*!
   * The values of the system data block are copied into a new one. The copy
   * is performed point-wise for all interior points, meaning points not on
   * the boundary.
   */
  void persist(T* out) const { grid::copy_interior(*this, out); }

  //! Copies the input data into the system values.
  /*!
   * The values of the system data block are initialized from the
   * given values, correctly transcribing all values.
   */
  void fill(const T* in) const { grid::fill_interior(in, *this, dims); }

 protected:
  SystemData() : BoundaryGridCUDA<T, D>{}, SystemInfo{} {}
};

//! Maintains grid data and associated information.
/*!
 * Maintains grid data and parameters, and defines the ability to use data
 * persistence by saving a snapshot. Specialization implementing
 * a system with boundaries.
 *
 * See SystemData<Grid<T, D>>.
 */
template <typename T, size_t D>
struct SystemData<RegionalGridCUDA<T, D>> : RegionalGridCUDA<T, D>, SystemInfo {
  using SystemInfo::id;
  using SystemInfo::info;

  using GridCUDA<T, D>::dims;
  using GridCUDA<T, D>::len;
  using RegionalGridCUDA<T, D>::as_grid;
  using RegionalGridCUDA<T, D>::region;

  //! Create a system data instance.
  /*!
   * Data for a phase field system is generated, resulting in a new
   * grid instance and information defining that grid.
   */
  SystemData(symphas::interval_data_type const& vdata, size_t id)
      : RegionalGridCUDA<T, D>{grid::construct<::RegionalGridCUDA, T, D>(
            vdata)},
        SystemInfo{{vdata}, id} {
    if (!info.intervals.empty()) {
      for (iter_type i = 0; i < D; ++i) {
        Axis ax = symphas::index_to_axis(i);
        auto& interval = info.intervals.at(ax);
        double boundary = region.boundary_size * interval.width();

        interval.set_domain_count(dims[i] - 2 * BOUNDARY_DEPTH);
        interval.set_interval(interval.left() + boundary,
                              interval.right() - boundary);
      }
    }
  }

  //! Get the number of elements of true data.
  /*!
   * The elements which are considered true data are the interior domain,
   * that is, all non-boundary elements. Return the count of these elements.
   */
  auto length() const { return grid::length_interior<D>(dims); }

  //! Copies the system data into the given array.
  /*!
   * The values of the system data block are copied into a new one. The copy
   * is performed point-wise for all interior points, meaning points not on
   * the boundary.
   */
  void persist(T* out) const { grid::copy_region(*this, out); }

  //! Copies the input data into the system values.
  /*!
   * The values of the system data block are initialized from the
   * given values, correctly transcribing all values.
   */
  void fill(const T* in) const { grid::fill_interior(in, *this, dims); }

 protected:
  SystemData() : RegionalGridCUDA<T, D>{}, SystemInfo{} {}
};

//! Maintains grid data parameters as well as a snapshot.
/*!
 * Maintains grid data parameters, and contains implementations used in data
 * persistence.
 *
 * The snapshot may be written to disk only if IO is enabled,
 * otherwise no write utility would be available. The snapshot and write utility
 * are available through the object WriteParallel, which is conditionally
 * compiled if IO is enabled.
 */
template <typename T, size_t D>
struct PersistentSystemData<GridCUDA<T, D>> : SystemData<GridCUDA<T, D>> {
  using parent_type = SystemData<GridCUDA<T, D>>;

  using parent_type::dims;
  using parent_type::id;
  using parent_type::info;
  using parent_type::len;
  using parent_type::length;
  using parent_type::persist;

  //! Create a system that can persist data to disk.
  /*!
   * Data for a phase field system is generated, resulting in a new
   * grid instance and information defining that grid.
   */
  PersistentSystemData(symphas::interval_data_type const& vdata, size_t id)
      : parent_type(vdata, id), writer{length()} {}

  void persist() const { persist(writer.get_snapshot().values); }

  //! Writes the current snapshot to the disk if IO is enabled.
  /*!
   * The current snapshot is written to disk if the IO macro keyword is
   * defined. The snapshot should be updated using persist before writing.
   */
  void write(const char* dir, iter_type index) const {
    symphas::io::write_info w{dir, index, id, DataFileType::SOLUTION_DATA};
    symphas::grid_info g{info.intervals};
    writer.write(w, g);
  }

  //! Writes the current snapshot to the disk if IO is enabled.
  /*!
   * The current snapshot is written to disk if the IO macro keyword is
   * defined. The snapshot should be updated using persist before writing.
   *
   * \param dir The directory to which the data file is saved.
   * \param name The name of the file to save to disk.
   */
  void write(const char* dir, const char* name, iter_type index) const {
    char* write_loc = new char[std::strlen(dir) + std::strlen(name) + 2];
    snprintf(write_loc, BUFFER_LENGTH, "%s/%s", dir, name);

    symphas::io::write_info w{write_loc, index, id, DataFileType::NAMED_DATA};
    symphas::grid_info g{info.intervals};
    writer.write(w, g);

    delete[] write_loc;
  }

  //! Save the data information to disk.
  /*!
   * The system data is saved to a disk.
   *
   * \param dir The directory to which the data file is saved.
   */
  inline void save(const char* dir, iter_type index) const {
    persist();
    write(dir, index);
  }

  //! Save the data information to disk with a name.
  /*!
   * The system data is saved to a disk.
   *
   * \param dir The directory to which the data file is saved.
   * \param name The name of the file to save to disk.
   */
  inline void save(const char* dir, const char* name, iter_type index) const {
    persist();
    write(dir, name, index);
  }

  //! Get a copy of the snapshot.
  /*!
   * This updates the member snapshot data with the current system values,
   * and copies the snapshot n the snapshot and returns it.
   */
  Block<T> get_snapshot() const {
    persist();
    return writer.get_snapshot();
  }

  using field_type = symphas::FieldAxis<D, T*>;
  operator field_type() const {
    auto values = std::shared_ptr<T[]>(new T[length()]);
    persist(values.get());

    return {std::shared_ptr<axis_nd_t<D>[]>(std::move(
                symphas::lib::new_system_axis_list<D>(info.intervals))),
            values, length()};
  }

  field_type as_field() const { return *this; }

 protected:
  PersistentSystemData() : PersistentSystemData{{}, 0} {}

  WriteParallel<T> writer;
};

//! Maintains grid data parameters as well as a snapshot.
/*!
 * Maintains grid data parameters, and contains implementations used in data
 * persistence.
 *
 * The snapshot may be written to disk only if IO is enabled,
 * otherwise no write utility would be available. The snapshot and write utility
 * are available through the object WriteParallel, which is conditionally
 * compiled if IO is enabled.
 */
template <typename T, size_t D>
struct PersistentSystemData<RegionalGridCUDA<T, D>>
    : SystemData<RegionalGridCUDA<T, D>> {
  using parent_type = SystemData<RegionalGridCUDA<T, D>>;

  using parent_type::dims;
  using parent_type::id;
  using parent_type::info;
  using parent_type::len;
  using parent_type::length;
  using parent_type::persist;
  using parent_type::region;
  using parent_type::values;

  //! Create a system that can persist data to disk.
  /*!
   * Data for a phase field system is generated, resulting in a new
   * grid instance and information defining that grid.
   */
  PersistentSystemData(symphas::interval_data_type const& vdata, size_t id)
      : parent_type(vdata, id) {}

  void write(symphas::io::write_info w, symphas::grid_info g,
             WriteParallel<T> const& writer) const {
    for (auto& [axis, interval] : g) {
      w.intervals[axis][0] = interval.domain_left();
      w.intervals[axis][1] = interval.domain_right();
    }

    for (auto& [axis, interval] : g) {
      interval.interval_to_domain();

      // double offset = region.boundary_size * interval.width();
      // interval.set_domain(interval.domain_left() - offset,
      //                     interval.domain_right() + offset);
    }

    g.update_strides();

    writer.write(w, g);
  }

  WriteParallel<T> persist() const {
    WriteParallel<T> writer(length());
    persist(writer.get_snapshot().values);
    return writer;
  }

  //! Writes the current snapshot to the disk if IO is enabled.
  /*!
   * The current snapshot is written to disk if the IO macro keyword is
   * defined. The snapshot should be updated using persist before writing.
   */
  void write(const char* dir, iter_type index,
             WriteParallel<T> const& writer) const {
    symphas::io::write_info w{dir, index, id, DataFileType::SOLUTION_DATA};
    symphas::grid_info g{info};
    write(w, g, writer);
  }

  //! Writes the current snapshot to the disk if IO is enabled.
  /*!
   * The current snapshot is written to disk if the IO macro keyword is
   * defined. The snapshot should be updated using persist before writing.
   *
   * \param dir The directory to which the data file is saved.
   * \param name The name of the file to save to disk.
   */
  void write(const char* dir, const char* name, iter_type index,
             WriteParallel<T> const& writer) const {
    char* write_loc = new char[std::strlen(dir) + std::strlen(name) + 2];
    snprintf(write_loc, BUFFER_LENGTH, "%s/%s", dir, name);

    symphas::io::write_info w{write_loc, index, id, DataFileType::NAMED_DATA};
    symphas::grid_info g{info};
    write(w, g, writer);

    delete[] write_loc;
  }

  //! Save the data information to disk.
  /*!
   * The system data is saved to a disk.
   *
   * \param dir The directory to which the data file is saved.
   */
  inline void save(const char* dir, iter_type index) const {
    auto writer = persist();
    write(dir, index, writer);
  }

  //! Save the data information to disk with a name.
  /*!
   * The system data is saved to a disk.
   *
   * \param dir The directory to which the data file is saved.
   * \param name The name of the file to save to disk.
   */
  inline void save(const char* dir, const char* name, iter_type index) const {
    auto writer = persist();
    write(dir, name, index, writer);
  }

  //! Get a copy of the snapshot.
  /*!
   * This updates the member snapshot data with the current system values,
   * and copies the snapshot n the snapshot and returns it.
   */
  Block<T> get_snapshot() const {
    Block<T> snapshot(length());
    persist(snapshot.values);
    return snapshot;
  }

  using field_type = symphas::FieldAxis<D, T*>;
  operator field_type() const {
    auto values = std::shared_ptr<T[]>(new T[length()]);
    persist(values.get());

    return {std::shared_ptr<axis_nd_t<D>[]>(std::move(
                symphas::lib::new_system_axis_list<D>(info.intervals))),
            values, length()};
  }

  field_type as_field() const { return *this; }

 protected:
  PersistentSystemData() : PersistentSystemData{{}, 0} {}
};

#endif
