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
 * PURPOSE: Implements string conversion functions for Side and Axis enums.
 *
 * ***************************************************************************
 */

#include "gridinfo.h"
#include "spslib.h"
#include <map>
#include <cstring>

namespace symphas::internal {

// Map from string to Side enum values
std::map<const char*, Side, symphas::lib::any_case_comparator> side_key_map = {
    {"LEFT", Side::LEFT},
    {"L", Side::LEFT},
    {"RIGHT", Side::RIGHT}, 
    {"R", Side::RIGHT},
    {"TOP", Side::TOP},
    {"T", Side::TOP},
    {"BOTTOM", Side::BOTTOM},
    {"B", Side::BOTTOM},
    {"FRONT", Side::FRONT},
    {"F", Side::FRONT},
    {"BACK", Side::BACK},
    {"K", Side::BACK}
};

// Map from string to Axis enum values
std::map<const char*, Axis, symphas::lib::any_case_comparator> axis_key_map = {
    {"X", Axis::X},
    {"Y", Axis::Y},
    {"Z", Axis::Z},
    {"T", Axis::T},
    {"S", Axis::S}, 
    {"R", Axis::R}
};

} // namespace symphas::internal

namespace symphas {

Side side_from_str(const char* str) {
    auto it = internal::side_key_map.find(str);
    return (it != internal::side_key_map.end()) ? it->second : static_cast<Side>(-1);
}

const char* str_from_side(Side side) {
    for (const auto& [key, value] : internal::side_key_map) {
        if (value == side) {
            return key;
        }
    }
    return nullptr;
}

Axis axis_from_str(const char* str) {
    if (std::strlen(str) == 0) {
        return Axis::NONE;
    }
    
    auto it = internal::axis_key_map.find(str);
    return (it != internal::axis_key_map.end()) ? it->second : Axis::NONE;
}

const char* str_from_axis(Axis axis) {
    for (const auto& [key, value] : internal::axis_key_map) {
        if (value == axis) {
            return key;
        }
    }
    return nullptr;
}

} // namespace symphas
