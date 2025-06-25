
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

#ifdef _MSC_VER
#define strncasecmp _strnicmp
#define strcasecmp _stricmp
#endif

#include <cstring>

#include "definitions.h"
#include "spsliblist.h"

namespace symphas {
namespace lib {

//! Transform a string into a file name format.
/*
 * Transform a string to lowercase and remove spaces. Additional
 * transformations are also made to convert any string to an appropriate
 * file name. The string is not changed in place, but is copied
 * into the given output, which is expected to be at least as long as
 * input.
 *
 * \param[in] in The input string which is transformed.
 * \param[out] out The result of converting the input to a file name format.
 * The function expects that this string is already allocated to a length
 * at least as long as the input.
 */
void to_file_name(char const* in, char* out, size_t count = 0);

//! Transform the string in place to lowercase.
/*!
 * The given string is manipulated so that all characters are converted
 * into lowercase.
 *
 * \param[inout] str The input string which is manipulated in place.
 */
void to_lower(char* str);

//! Transform the string in place to lowercase.
/*!
 * The given string is transformed into the output string so that all
 * characters are converted into lowercase.
 *
 * \param[in] str The input string which is converted to lowercase.
 * \param[out] out The result of converting the input to lowercase.
 */
void to_lower(const char* str, char* out);

//! Transform the string in place to uppercase.
/*!
 * The given string is manipulated so that all characters are converted
 * into uppercase.
 *
 * \param[inout] str The input string which is manipulated in place.
 */
void to_upper(char* str);

//! Transform the string in place to uppercase.
/*!
 * The given string is transformed into the output string so that all
 * characters are converted into uppercase.
 *
 * \param[in] str The input string which is converted to uppercase.
 * \param[out] out The result of converting the input to uppercase.
 */
void to_upper(const char* str, char* out);

//! Trim whitespace from the beginning and end.
/*!
 * The given string is manipulated so that all whitespace is removed from
 * the beginning and end. The whitespace characters includes `\n` `\t` and
 * `\r` in addition to space, ` `.
 *
 * \param[inout] str The input string which is manipulated in place.
 */
void str_trim(char* str);

//! Trim whitespace from the beginning and end.
/*!
 * The given string is manipulated so that all whitespace is removed from
 * the beginning and end. The whitespace characters includes `\n` `\t` and
 * `\r` in addition to space, ` `.
 *
 * \param[in] str The input string which is trimmed.
 * \param[out] out The result of trimming the input.
 */
void str_trim(const char* str, char* out);

//! Returns offset to next non-digit character, which may the first one.
/*!
 * Returns offset to a position in the given name string to the next
 * character which is not a digit, which can also be the end of the string.
 * In the case that the offset is the end of the string, then the offset
 * will be equal to the result of `std::strlen`. So this case, the result
 * of dereferencing the input by the result would be the terminating
 * character.
 *
 * \param content The pointer to the start of the string that will be iterated
 * over until the current index no longer refers to a digit.
 */
iter_type pos_after_digit(const char* content);

//! Returns offset to next character after matching the start of the string.
/*!
 * Returns offset to a position in the given string to the next
 * character that follows the provided token, if the token exists. This
 * can also be the end of the string, see pos_after_digit. In the case the
 * input does not start with the given token, the function will return
 * 0. Likewise if there is no complete tokens, an offset to a partial match
 * (where the token no longer matches) is still considered 0. The match
 * is performed index by index.
 *
 * \param content the pointer to the start of the string that will be compared
 * with the token index by index to either return the position or return 0 if
 * the token is non-matching.
 * \param token The token which is matched against the start of the given
 * string.
 */
iter_type pos_after_token(const char* content, const char* token);

//! Returns offset to next character after a token in the set if matching.
/*!
 * Returns offset to a position in the given name string to the next
 * character which follows the end of the token string for one of the tokens
 * in the provided set. Regarding return value and end of string behaviour,
 * see pos_after_digit.
 *
 * \param content the pointer to the start of the string that will be compared
 * with the token index by index to either return the position or return 0 if
 * the token is non-matching.
 * \param tokens The tokens which are matched against the start of the given
 * string.
 */
iter_type pos_after_token(const char* content,
                          std::initializer_list<const char*> tokens);

namespace {

template <char...>
struct char_list {};

template <typename A, typename B>
struct filter_decimal_left;

template <char... cs>
struct filter_decimal_left<char_list<cs...>, char_list<>> {
  using type = char_list<cs...>;
};

template <char... cs, char c0, char... c1s>
struct filter_decimal_left<char_list<cs...>, char_list<c0, c1s...>> {
  using type = typename filter_decimal_left<char_list<cs..., c0>,
                                            char_list<c1s...>>::type;
};

template <char... cs, char... c1s>
struct filter_decimal_left<char_list<cs...>, char_list<'.', c1s...>> {
  using type = char_list<cs...>;
};

template <typename A, typename B>
struct filter_decimal_right;

template <>
struct filter_decimal_right<char_list<>, char_list<>> {
  using type = char_list<>;
};

template <char... cs>
struct filter_decimal_right<char_list<'.', cs...>, char_list<>> {
  using type = char_list<'.', cs...>;
};

template <char... cs, char c0, char... c1s>
struct filter_decimal_right<char_list<'.', cs...>, char_list<c0, c1s...>> {
  using type = typename filter_decimal_right<char_list<'.', cs..., c0>,
                                             char_list<c1s...>>::type;
};

template <char... c1s>
struct filter_decimal_right<char_list<>, char_list<'.', c1s...>> {
  using type =
      typename filter_decimal_right<char_list<'.'>, char_list<c1s...>>::type;
};

template <char c0, char... c1s>
struct filter_decimal_right<char_list<>, char_list<c0, c1s...>> {
  using type =
      typename filter_decimal_right<char_list<>, char_list<c1s...>>::type;
};

template <char n>
constexpr size_t get_value() {
  return n - '0';
}

}  // namespace

template <char... n>
struct char_to_number {
  template <char... ln, size_t... Is>
  constexpr int left(char_list<ln...>, std::index_sequence<Is...>) const {
    return (int)((get_value<ln>() *
                  fixed_pow<10, sizeof...(Is) - 1 - int(Is)>)+... +
                 0);
  }

  template <char... ln>
  constexpr int left(char_list<ln...>) const {
    return left(char_list<ln...>{}, std::make_index_sequence<sizeof...(ln)>{});
  }

  template <char... rn, size_t... Is>
  constexpr double right(char_list<'.', rn...>,
                         std::index_sequence<Is...>) const {
    return ((get_value<rn>() / double(fixed_pow<10, Is + 1>)) + ... + 0);
  }

  template <char... rn>
  constexpr double right(char_list<'.', rn...>) const {
    return right(char_list<'.', rn...>{},
                 std::make_index_sequence<sizeof...(rn)>{});
  }

  constexpr int right(char_list<>) const { return 0; }

  constexpr int left() const {
    return left(
        typename filter_decimal_left<char_list<>, char_list<n...>>::type{});
  }

  constexpr auto right() const {
    return right(
        typename filter_decimal_right<char_list<>, char_list<n...>>::type{});
  }

  constexpr auto operator()() const { return left() + right(); }
};

struct string : array_container<char> {
  using array_container<char>::array_container;
  string() : string("") {}  // Default constructor initializes to empty string
  string(const char* str) : array_container(std::strlen(str) + 1) {
    std::strcpy(array_container<char>::data, str);
  }
  bool operator==(string const& other) const {
    return std::strcmp(array_container<char>::data, other.data) == 0;
  }
  bool empty() const { return array_container<char>::data[0] == '\0'; }
};

inline bool operator<(string const& a, string const& b) {
  return std::strcmp(a, b) < 0;
}

struct any_case_comparator {
  bool operator()(std::string const& s1, std::string const& s2) const {
    return strcasecmp(s1.c_str(), s2.c_str()) < 0;
  }
  bool operator()(std::string const& s1, const char* s2) const {
    return strcasecmp(s1.c_str(), s2) < 0;
  }
  bool operator()(const char* s1, std::string const& s2) const {
    return strcasecmp(s1, s2.c_str()) < 0;
  }
  bool operator()(const char* s1, const char* s2) const {
    return strcasecmp(s1, s2) < 0;
  }
};

}  // namespace lib
}  // namespace symphas