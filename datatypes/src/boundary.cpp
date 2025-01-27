
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
 */

#include "boundary.h"

symphas::b_element_type::b_element_type(BoundaryType type,
                                        std::initializer_list<BoundaryTag> tags,
                                        const double* params, int argc)
    : type{type},
      tag{tags.begin()[0], tags.begin()[1]},
      params{(argc > 0) ? new double[argc] : nullptr},
      argc{argc} {
  if (type == BoundaryType::PERIODIC0 || type == BoundaryType::PERIODIC3XZ ||
      type == BoundaryType::PERIODIC3XY) {
    this->type = BoundaryType::PERIODIC;
  }

  std::copy(params, params + argc, this->params);
}

symphas::b_element_type::b_element_type(b_element_type const& other)
    : b_element_type(other.type, other.tag, other.params, other.argc) {}

symphas::b_element_type::b_element_type(b_element_type&& other) noexcept
    : b_element_type() {
  swap(*this, other);
}

symphas::b_element_type& symphas::b_element_type::operator=(
    b_element_type other) {
  symphas::b_element_type b{other};
  swap(*this, b);
  return *this;
}

symphas::b_element_type::~b_element_type() { delete[] params; }

void symphas::b_element_type::set_parameter(double value, iter_type n) {
  if (argc <= n) {
    double* params_append = new double[n + 1]{};
    for (iter_type i = 0; i < argc; ++i) {
      params_append[i] = params[i];
    }
    for (iter_type i = argc; i < n; ++i) {
      params_append[i] = 0.0;
    }

    delete[] params;
    params = params_append;
    argc = n + 1;
  }
  params[n] = value;
}

std::map<const char*, BoundaryType, symphas::lib::any_case_comparator>
    symphas::internal::boundary_key_map = {{"DEFAULT", BoundaryType::DEFAULT},
                                           {"PERIODIC", BoundaryType::PERIODIC},
                                           {"OPEN", BoundaryType::OPEN}};

std::map<const char*, BoundaryTag, symphas::lib::any_case_comparator>
    symphas::internal::boundary_tag_key_map = {
        {"GAUSSIAN", BoundaryTag::GAUSSIAN},  {"GA", BoundaryTag::GAUSSIAN},
        {"RANDOM", BoundaryTag::RANDOM},      {"RA", BoundaryTag::RANDOM},
        {"TRIGONOMETRIC", BoundaryTag::TRIG}, {"TR", BoundaryTag::TRIG},
        {"CONSTANT", BoundaryTag::CONSTANT},  {"CO", BoundaryTag::CONSTANT},
        {"LINEAR", BoundaryTag::LINEAR},      {"LI", BoundaryTag::LINEAR},
        {"TIME", BoundaryTag::TIME},          {"T", BoundaryTag::TIME}};

BoundaryType symphas::boundary_from_str(const char* type) {
  auto condition = internal::boundary_key_map.find(type);
  if (condition != internal::boundary_key_map.end()) {
    return condition->second;
  } else {
    if (boundary_tag_from_str(type) != BoundaryTag::NONE) {
      return BoundaryType::DEFAULT;
    }
  }
  return BoundaryType::NONE;
}

const char* symphas::str_from_boundary(BoundaryType in) {
  for (auto& [k, v] : internal::boundary_key_map) {
    if (v == in) {
      return k;
    }
  }
  return nullptr;
}

BoundaryTag symphas::boundary_tag_from_str(const char* type) {
  auto tag = internal::boundary_tag_key_map.find(type);
  if (std::strlen(type) == 0) {
    return BoundaryTag::NONE;
  }

  if (tag != internal::boundary_tag_key_map.end()) {
    return tag->second;
  }
  return BoundaryTag::NONE;
}

const char* symphas::str_from_boundary_tag(BoundaryTag tag) {
  for (auto& [k, v] : internal::boundary_tag_key_map) {
    if (v == tag) {
      return k;
    }
  }
  return nullptr;
}
