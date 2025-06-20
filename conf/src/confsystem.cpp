
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

#include "confsystem.h"

#include <regex>

namespace symphas::internal {

/* the names of all the configuration parameters
 */

char C_DIM[] = "DIM";
char C_SOLVERVAR[] = "SOLVERVAR";
char C_ORDER[] = "ORDER";
char C_PTL[] = "PTL";
char C_PTB[] = "PTB";
char C_PTG[] = "PTG";
char C_WIDTH[] = "WIDTH";
char C_FORM[] = "FORM";
char C_RNGX[] = "RNGX";
char C_RNGY[] = "RNGY";
char C_RNGZ[] = "RNGZ";
char C_RNGR[] = "RNGR";
char C_RNGT[] = "RNGT";
char C_RNGF[] = "RNGF";
char C_BNDLT[] = "BNDLT";
char C_BNDRT[] = "BNDRT";
char C_BNDTP[] = "BNDTP";
char C_BNDBT[] = "BNDBT";
char C_BNDFT[] = "BNDFT";
char C_BNDBK[] = "BNDBK";
char C_INSIDE[] = "INSIDE";
char C_DELTA[] = "DELTA";
char C_MODEL[] = "MODEL";
char C_NAMES[] = "NAMES";
char C_FRAMES[] = "FRAMES";
char C_SAVE[] = "SAVE";
char C_SAVEINIT[] = "SAVEINIT";
char C_RUNS[] = "RUNS";
char C_DIR[] = "DIR";

//! The opening and closing brackets of a coefficient index.
/*!
 * Defines the allowable brackets which are used in the grammar of
 * directly specifying a coefficient index in the coefficient file.
 */
constexpr std::initializer_list<std::pair<char, char>> COEFF_ID_BRACKETS = {
    {'[', ']'}, {'(', ')'}};

//! Check whether the coefficient name is a coefficient id.
/*!
 * Check whether the start and end of the coefficient name is bounded by the
 * coefficient bounding characters defined by #COEFF_ID_BRACKETS. If the name is
 * bounded, then the string between the bounding characters refers to a
 * coefficient id and the function returns true. The content between the
 * brackets is not checked and might be invalid, but as long as the string is
 * bounded the function will return true. If the brackets types are not
 * correctly matched, i.e. using `[xxx)`, the function will return false.
 *
 * The given input is expected to be trimmed (no whitespace at the start or end
 * of the string).
 *
 * \param name The coefficient name that is checked to be in the format
 * corresponding to direct index.
 */
bool is_name_coeff_id(const char* name) {
  size_t len = std::strlen(name);

  for (auto [openc, closec] : COEFF_ID_BRACKETS) {
    if (name[0] == openc && name[len - 1] == closec) {
      return true;
    }
  }

  return false;
}

constexpr std::initializer_list<const char*> COEFF_ID_RANGE_STRS = {"..", "-"};
constexpr std::initializer_list<const char*> COEFF_ID_SEP_STRS = {","};

void assign_coeff_by_index(const char* value_str, double*& coeff,
                           size_t& coeff_len) {
  if (std::all_of(value_str, value_str + std::strlen(value_str), ::isspace)) {
    return;
  }

  std::string value_pat = "([-\\.\\deE\\+]+)";
  std::regex r_range("\\[(\\d+)..(\\d+)\\]=" + value_pat);
  std::regex r_single("\\[(\\d+)\\]=" + value_pat);
  std::regex c_single("c(\\d+)=" + value_pat);
  std::smatch match;

  std::string s(value_str);
  if (std::regex_search(s, match, r_range) && match.size() > 3) {
    int start_index = std::stoi(match.str(1));
    int end_index = std::stoi(match.str(2));
    double value = std::stod(match.str(3));

    if (end_index >= coeff_len) {
      symphas::lib::expand_append_array(DEFAULT_COEFF_VALUE, coeff,
                                        end_index + 1 - coeff_len, coeff_len);
      coeff_len = end_index + 1;
    }

    if (start_index >= 0) {
      for (int i = start_index; i <= end_index; ++i) {
        coeff[i] = value;
      }
    }
  } else if (std::regex_search(s, match, r_single) && match.size() > 2) {
    int index = std::stoi(match.str(1));
    double value = std::stod(match.str(2));

    if (index >= coeff_len) {
      symphas::lib::expand_append_array(DEFAULT_COEFF_VALUE, coeff,
                                        index + 1 - coeff_len, coeff_len);
      coeff_len = index + 1;
    }

    coeff[index] = value;
  } else if (std::regex_search(s, match, c_single) && match.size() > 2) {
    int index = std::stoi(match.str(1)) - 1;
    double value = std::stod(match.str(2));

    if (index >= coeff_len) {
      symphas::lib::expand_append_array(DEFAULT_COEFF_VALUE, coeff,
                                        index + 1 - coeff_len, coeff_len);
      coeff_len = index + 1;
    }

    coeff[index] = value;
  } else {
    fprintf(SYMPHAS_WARN,
            "coefficient format '%s' to assign isn't recognized,"
            " line is skipped\n",
            value_str);
  }
}

// parse one line from a coefficient initialization file
void coeff_read_entry(const char* line, int i, double*& coeff,
                      size_t& coeff_len) {
  char* line_buffer = new char[std::strlen(line) + 1];
  std::strcpy(line_buffer, line);
  symphas::lib::str_trim(line_buffer);

  // check if the line is not a comment, or empty
  if (*line_buffer && *line_buffer != CONFIG_COMMENT_PREFIX_C) {
    char* endptr = nullptr;
    double value = strtod(line_buffer, &endptr);
    if (endptr > line_buffer) {
      // The line is a simple number
      if (i >= coeff_len) {
        symphas::lib::expand_append_array(DEFAULT_COEFF_VALUE, coeff,
                                          i + 1 - coeff_len, coeff_len);
        coeff_len = i + 1;
      }
      coeff[i] = value;
    } else {
      // The line is in the format [a]=value or [a..b]=value
      assign_coeff_by_index(line_buffer, coeff, coeff_len);
    }
  }
  delete[] line_buffer;
}

// set the coefficients using the file
void coeff_from_file(FILE* f, double*& coeff, size_t& coeff_len) {
  char line[BUFFER_LENGTH];
  while (fgets(line, BUFFER_LENGTH, f) != NULL) {
    symphas::lib::str_trim(line);
    if (*line && *line != CONFIG_COMMENT_PREFIX_C) {
      assign_coeff_by_index(line, coeff, coeff_len);
    }
  }
}

// Parse the coefficients list, optionally looking at the given directory
// if there is a coefficients file provided.
static void parse_coeff_list(const char* str, const char* dir, double*& coeff,
                             size_t& coeff_len) {
  char constants[LINE_READ_BUFFER];
  symphas::lib::str_trim(str, constants);

  /*!
   * If the first character of the constants is #CONFIG_OPTION_PREFIX_C, then
   * we will read the constants from a file instead. All the constants which are
   * not mentioned in the file, but which should be set because a coefficient of
   * a higher index is initialized, will be initialized to #DEFAULT_COEFF_VALUE.
   */
  if (constants[0] == CONFIG_OPTION_PREFIX_C) {
    char fname[BUFFER_LENGTH_R2];

    if (sscanf(constants, STR(CONFIG_OPTION_PREFIX) "%s", fname) != 1) {
      fprintf(SYMPHAS_ERR,
              "a '%c' was given in the coefficients but an "
              "incorrect format provided: '%s'\n",
              CONFIG_OPTION_PREFIX_C, fname);
      exit(824);
    }

    // add the extension to the constants file name if it is not already
    // appended, even if there is another string after an occurrence of a dot
    char* dot = std::strrchr(constants, '.');
    if (dot == NULL || std::strcmp(dot + 1, COEFF_SPEC_EXTENSION) != 0) {
      std::strcat(fname, "." COEFF_SPEC_EXTENSION);
    }

    FILE* f;
    fname[BUFFER_LENGTH_R2 - 1] = '\0';
    if ((f = fopen(fname, "r")) == 0) {
      char* fname_with_dir =
          new char[std::strlen(dir) + std::strlen(fname) + 1];
      sprintf(fname_with_dir, "%s%s", dir, fname);
      if ((f = fopen(fname_with_dir, "r")) == 0) {
        fprintf(SYMPHAS_ERR, "constants file named '%s' could not be read\n",
                fname);
        exit(825);
      }
      delete[] fname_with_dir;
    }

    coeff_from_file(f, coeff, coeff_len);
  } else if (std::strlen(str) > 0) {
    std::istringstream iss(constants);
    std::string constant;
    iter_type i = 0;
    while (std::getline(iss, constant, ' ')) {
      coeff_read_entry(constant.c_str(), i++, coeff, coeff_len);
    }
  }
}

// Parse the simple coefficients list
static void parse_simple_coeff_list(const char* str, double*& coeff,
                                    size_t& coeff_len) {
  coeff_len = 0;
  delete[] coeff;
  coeff = nullptr;

  char* strcpy = new char[std::strlen(str) + 1];
  std::strcpy(strcpy, str);
  char* tok = std::strtok(strcpy, " ");

  if (tok) {
    do {
      char* endptr;
      double parameter = strtod(tok, &endptr);

      if (!*endptr) {
        double* new_coeff = new double[++coeff_len];
        for (iter_type i = 0; i < coeff_len - 1; ++i) {
          new_coeff[i] = coeff[i];
        }
        new_coeff[coeff_len - 1] = parameter;

        delete[] coeff;
        coeff = new_coeff;
      } else {
        fprintf(SYMPHAS_WARN,
                "the value '%s' given to the initial conditions "
                "could not be interpreted, skipping...\n",
                tok);
      }
    } while ((tok = std::strtok(NULL, " ")) != 0);
  }
  delete[] strcpy;
}

//! Delimiter list used for configuration options.
/*!
 * Some parameters can be given as a list of options with delimiters combining
 * individual elements. There are several types of delimiters, and currently
 * all of them are valid for packaging a full option inside a config parameter.
 */
inline std::initializer_list<char> OPTION_DELIMITERS = {'\'', '"'};

//! Brackets used to group configuration options.
/*!
 * Brackets are used in the configuration to group together an individual
 * configuration item, rather than individual parameters for one configuration.
 * For instance, if a name requires spaces, then the characters specified in
 * #OPTION_DELIMITERS are used. In order to group together an entire
 * configuration input for parameters which support multiple configuration
 * (as an example, see ConfSystem::intervals), then brackets are used for each
 * configuration.
 *
 * Everything inside the brackets will be matched until appearance of the
 * closing one.
 */
inline std::initializer_list<std::pair<char, char>> OPTION_BRACKETS = {
    {'{', '}'}, {'(', ')'}};

//! Returns open and closed brackets as their own lists.
/*!
 * Using the bracket data in #OPTION_BRACKETS, all open brackets are aggregated
 * into one list and all close brackets are aggregated in one list. Both of
 * these lists are returned in a pair object.
 */
std::pair<std::vector<char>, std::vector<char>> bracket_lists() {
  std::vector<char> out_open, out_close;
  for (auto [open, close] : OPTION_BRACKETS) {
    out_open.push_back(open);
    out_close.push_back(close);
  }
  return {out_open, out_close};
}

char option_open_bracket() {
  return symphas::internal::OPTION_BRACKETS.begin()->first;
}

char option_close_bracket() {
  return symphas::internal::OPTION_BRACKETS.begin()->second;
}
}  // namespace symphas::internal

//! Parse the options provided and return a list of configurations.
/*!
 * A given string is split by option
 *
 * If none of the delimiters are found, then the algorithm will assume there is
 * only a single option in particular, if `spaces_are_delimiters' is false, then
 * a configuration parameter with spaces will be one option.
 *
 * this will put each individual option in an std::string
 */
std::vector<std::string> symphas::conf::parse_options(
    const char* options, bool spaces_are_delimiters,
    const char* extra_delimiters) {
  std::vector<std::string> out;
  char* option_buffer = new char[std::strlen(options) + 1];
  std::strcpy(option_buffer, options);
  symphas::lib::str_trim(option_buffer);

  std::vector<char> delimiters{symphas::internal::OPTION_DELIMITERS};
  char add_delimiter;
  while ((add_delimiter = *extra_delimiters++)) {
    delimiters.push_back(add_delimiter);
  }

  auto open_list = symphas::internal::bracket_lists().first;
  auto close_list = symphas::internal::bracket_lists().second;
  for (const char* c = option_buffer; *c;) {
    /*
     * copy the names into the preallocated buffers, from
     * input between double or single quotes, or spaces
     */

    char buffer[BUFFER_LENGTH];
    bool add_option = false;

    for (char chk : delimiters) {
      if (*c == chk) {
        char* b = buffer;
        while (*(++c) && *c != chk) {
          *(b++) = *c;
        }
        *b = '\0';
        add_option = true;
      }
    }

    if (!add_option) {
      /* to be general, this determines the last closing bracket, even if
       * there is more than one nested bracket type
       */
      for (auto [open, close] : symphas::internal::OPTION_BRACKETS) {
        if (*c == open) {
          std::stack<std::ptrdiff_t> bracket_nests;

          // continue gobbling until the final bracket is found
          char* b = buffer;
          while (!(*(++c) == close && bracket_nests.empty())) {
            /* add a bracket index if an opening bracket is found
             */
            auto it_open = std::find(open_list.begin(), open_list.end(), *c);
            if (it_open != open_list.end()) {
              bracket_nests.push(std::distance(open_list.begin(), it_open));
            }

            /* remove a bracket index if a closing bracket is found
             * if the found closing bracket doesn't match the previous
             * opening bracket, give an error
             */
            auto it_close = std::find(close_list.begin(), close_list.end(), *c);
            if (it_close != close_list.end()) {
              if (std::distance(open_list.begin(), it_close) ==
                  bracket_nests.top()) {
                bracket_nests.pop();
              } else {
                size_t pos = static_cast<size_t>(c - option_buffer);
                fprintf(
                    SYMPHAS_WARN,
                    "there is an unmatched bracket in the options starting at "
                    "character %zd in options '%s'\n",
                    pos, option_buffer);
              }
            }

            *(b++) = *c;
          }
          *b = '\0';
          ++c;
          add_option = true;
        }
      }
    }

    if (!add_option) {
      /* if each delimiter has been checked but nothing has been copied,
       * then copy in everything until a delimiter has been found
       */
      if (*c != ' ') {
        char* b = buffer;

        if (spaces_are_delimiters) {
          while (*c && *c != ' ') {
            *(b++) = *(c++);
          }

        } else {
          auto a0 = open_list.begin();
          auto a1 = open_list.end();
          auto b0 = delimiters.begin();
          auto b1 = delimiters.end();

          while (*c && (std::find(a0, a1, *c) == a1) &&
                 (std::find(b0, b1, *c) == b1)) {
            *(b++) = *(c++);
          }
        }

        *b = '\0';
        symphas::lib::str_trim(buffer);
        add_option = true;
      } else {
        for (; *c == ' '; ++c)
          ;  // skip all the next spaces
      }
    }

    /* copy the buffer into an option string if an option
     * was copied
     */
    if (add_option) {
      out.emplace_back(buffer);
    } else {
      //++c;
    }
  }

  delete[] option_buffer;
  return out;
}

symphas::lib::string handle_substitutions(const char* value) {
  const char* pos_sub = std::strchr(value, CONFIG_SUBSTITUTION_PREFIX_C);
  if (pos_sub == NULL) {
    symphas::lib::string result(value, std::strlen(value) + 1);
    return result;
  } else {
    const char* start = value;
    const char* end = NULL;
    char* buffer = new char[1]{};
    while (pos_sub != NULL) {
      int index = -1;
      const char* it;
      for (auto [open, close] : symphas::internal::OPTION_BRACKETS) {
        it = pos_sub + 1;
        if (*it == open) {
          ++it;
          while (*it >= '0' && *it <= '9') {
            ++it;
          }
          if (*it == close) {
            index = atoi(pos_sub + 2);
            end = pos_sub;
            break;
          }
        }
      }

      if (index > 0) {
        const char* default_substitution = "";
        const char* substitution;
        if (index < params::config_key_values.size()) {
          auto& config_entry = params::config_key_values[index - 1];
          substitution = config_entry.data;
        } else {
          substitution = default_substitution;
        }

        len_type size = end - start;
        len_type substitution_len = std::strlen(substitution);
        len_type buffer_len = std::strlen(buffer);

        char* append_buffer =
            new char[buffer_len + size + substitution_len + 1]{};
        sprintf(append_buffer, "%s", buffer);
        std::copy(start, end, append_buffer + buffer_len);
        sprintf(append_buffer + buffer_len + size, "%s", substitution);

        std::swap(buffer, append_buffer);
        delete[] append_buffer;

        end = it + 1;
        start = end;
        pos_sub = std::strchr(end + 1, CONFIG_SUBSTITUTION_PREFIX_C);
      }
    }

    end = (value + std::strlen(value));
    len_type size = end - start;
    len_type buffer_len = std::strlen(buffer);

    symphas::lib::string result(buffer_len + size + 1);
    sprintf(result.data, "%s", buffer);
    std::copy(start, end, result.data + buffer_len);
    result.data[buffer_len + size] = '\0';

    delete[] buffer;
    return result;
  }
}

void SymPhasSettings::init_work_dirs() const {
  // #ifdef FILESYSTEM_HEADER_AVAILABLE
  std::filesystem::path p_dir(".");
  p_dir.append(directory_settings.result_dir);

  symphas::lib::make_directory(p_dir / CHECKPOINT_DIR, 882);
  symphas::lib::make_directory(p_dir / DATA_DIR, 883);
  symphas::lib::make_directory(p_dir / PLOT_DIR, 884);
  //
  // #else
  //
  //  size_t p_dir_len = std::strlen(directory_settings.result_dir) + 3;
  //  char* p_dir = new char[p_dir_len];
  //  snprintf(p_dir, p_dir_len, "./%s", directory_settings.result_dir);
  //
  //  size_t sub_dir_len = p_dir_len + BUFFER_LENGTH;
  //  char* sub_dir = new char[sub_dir_len];
  //
  //  snprintf(sub_dir, p_dir_len + sizeof(CHECKPOINT_DIR) / sizeof(char) + 1,
  //           "%s/%s", p_dir, CHECKPOINT_DIR);
  //  symphas::lib::make_directory(sub_dir, 882);
  //  snprintf(sub_dir, p_dir_len + sizeof(DATA_DIR) / sizeof(char) + 1,
  //  "%s/%s",
  //           p_dir, DATA_DIR);
  //  symphas::lib::make_directory(sub_dir, 882);
  //  snprintf(sub_dir, p_dir_len + sizeof(PLOT_DIR) / sizeof(char) + 1,
  //  "%s/%s",
  //           p_dir, PLOT_DIR);
  //  symphas::lib::make_directory(sub_dir, 882);
  //
  //  delete[] p_dir;
  //  delete[] sub_dir;
  //
  // #endif
}

void SimulationSettings::parse_dimension(const char* value) {
  len_type available_dimensions[]{AVAILABLE_DIMENSIONS};
  if (value != NULL) {
    if (value[0] == CONFIG_OPTION_PREFIX_C) {
      if (sizeof(available_dimensions) / sizeof(len_type) > 0) {
        dimension = available_dimensions[0];
      } else {
        dimension = 0;
      }
    } else {
      dimension = atoi(value);
    }
  } else {
    if (sizeof(available_dimensions) / sizeof(len_type) > 0) {
      dimension = available_dimensions[0];
    } else {
      dimension = 0;
    }
  }
}

void DomainSettings::set_dimensions(size_t n) {
  const char msg[] = "the number of grid points cannot be less than 1\n";
  delete[] dims[n];

  size_t dimension = intervals[n].size();
  switch (dimension) {
    case 1: {
      iter_type L = intervals[n].at(Axis::X).get_domain_count();

      if (L < 1) {
        fprintf(SYMPHAS_LOG, msg);
        exit(802);
      }

      dims[n] = new len_type[1]{L};

      break;
    }
    case 2: {
      iter_type L = intervals[n].at(Axis::X).get_domain_count();
      iter_type M = intervals[n].at(Axis::Y).get_domain_count();

      if (L < 1 || M < 1) {
        fprintf(SYMPHAS_LOG, msg);
        exit(803);
      }

      dims[n] = new len_type[2]{L, M};

      break;
    }
    case 3: {
      iter_type L = intervals[n].at(Axis::X).get_domain_count();
      iter_type M = intervals[n].at(Axis::Y).get_domain_count();
      iter_type N = intervals[n].at(Axis::Z).get_domain_count();

      if (L < 1 || M < 1 || N < 1) {
        fprintf(SYMPHAS_LOG, msg);
        exit(804);
      }

      dims[n] = new len_type[3]{L, M, N};

      break;
    }
    default:
      fprintf(SYMPHAS_ERR, "incompatible dimensions\n");
      exit(805);
  }
}

void DomainSettings::parse_boundaries_array(const char* const* b_specs,
                                            size_t dimension) {
  std::vector<std::string>* boundary_lists =
      new std::vector<std::string>[dimension * 2];
  size_t* boundary_lens = new size_t[dimension * 2];

  // compile all the boundaries given to the configuration
  for (iter_type i = 0; i < dimension * 2; ++i) {
    boundary_lists[i] = symphas::conf::parse_options(b_specs[i]);
    boundary_lens[i] = boundary_lists[i].size();
  }

  // set the number of boundary options given and top off any boundaries
  // which haven't been fully specified
  size_t len = *std::max_element(boundary_lens, boundary_lens + dimension * 2);
  if (len > bdata_len) {
    delete[] bdata;

    bdata_len = len;
    bdata = new symphas::b_data_type[bdata_len];
  }

  for (iter_type i = 0; i < dimension * 2; ++i) {
    // if the boundary specification for a given side is not "full", i.e. not
    // of all systems have boundaries specified, then append copies of the first
    // boundary to the rest.
    if (boundary_lens[i] < len) {
      std::vector<std::string> append{len - boundary_lens[i],
                                      boundary_lists[i][0]};
      boundary_lists[i].insert(boundary_lists[i].end(),
                               std::make_move_iterator(append.begin()),
                               std::make_move_iterator(append.end()));
    }
  }

  // set the boundaries for each of the options by copying the spec string
  // into the temporary array
  const char** b_spec = new const char*[dimension * 2];
  for (iter_type n = 0; n < len; ++n) {
    for (iter_type i = 0; i < dimension * 2; ++i) {
      b_spec[i] = boundary_lists[i][n].c_str();
    }

    parse_boundaries(b_spec, n, dimension);
  }

  delete[] boundary_lists;
  delete[] boundary_lens;
  delete[] b_spec;
}

void DomainSettings::parse_boundaries_array(Side side, const char* b_spec,
                                            size_t dimension) {
  auto boundary_list = symphas::conf::parse_options(b_spec);
  size_t boundary_len = boundary_list.size();
  if (boundary_len > bdata_len) {
    symphas::lib::resize_array(boundary_len, bdata, bdata_len);
  }

  for (iter_type i = 0; i < boundary_len; ++i) {
    parse_boundaries(side, boundary_list[i].c_str(), i);
  }
}

void DomainSettings::parse_boundaries(const char* const* b_spec, size_t n,
                                      size_t dimension) {
  /* set the boundary information
   *
   * in order to set the correct intervals, the index n is maintained
   * moreover, the boundaries are filled in reverse order, since the intervals
   * correspond with the boundaries this way
   */

  switch (dimension) {
    case 1:
      for (auto side : {Side::LEFT, Side::RIGHT}) {
        parse_boundaries(side, b_spec[symphas::side_to_index(side)], n);
      }
      break;

    case 2:
      for (auto side : {Side::LEFT, Side::RIGHT, Side::TOP, Side::BOTTOM}) {
        parse_boundaries(side, b_spec[symphas::side_to_index(side)], n);
      }
      break;

    case 3:
      for (auto side : {Side::LEFT, Side::RIGHT, Side::TOP, Side::BOTTOM,
                        Side::FRONT, Side::BACK}) {
        parse_boundaries(side, b_spec[symphas::side_to_index(side)], n);
      }
      break;

    default:
      break;
  }
}

void DomainSettings::parse_boundaries(Side side, const char* b_spec, size_t n) {
  /* set the boundary information
   *
   * in order to set the correct intervals, the index n is maintained
   * moreover, the boundaries are filled in reverse order, since the intervals
   * correspond with the boundaries this way
   */
  auto [it, flag] = bdata[n].emplace(side, symphas::b_element_type{});
  if (!flag) {
    return;
  }

  symphas::b_element_type& b = it->second;

  char typestr[BUFFER_LENGTH_R2], str[BUFFER_LENGTH];
  int end;

  int args = sscanf(b_spec, "%s%n", typestr, &end);

  b.type = symphas::boundary_from_str(typestr);
  BoundaryTag tag = symphas::boundary_tag_from_str(typestr);

  bool type_as_tag =
      (b.type == BoundaryType::DEFAULT && tag != BoundaryTag::NONE);
  if (type_as_tag) {
    b.tag[0] = tag;
  }

  symphas::lib::str_trim(b_spec + end, str);
  if (end < std::strlen(b_spec)) {
    char* ptr;
    char* tok;

    if ((tok = std::strtok(str, " ")) != NULL) {
      double value = std::strtod(tok, &ptr);
      if (ptr == tok) {
        iter_type i = (type_as_tag) ? 1 : 0;
        do {
          if (tok != NULL) {
            if (i < sizeof(b.tag) / sizeof(BoundaryTag)) {
              b.tag[i++] = symphas::boundary_tag_from_str(tok);
            }
            tok = std::strtok(NULL, " ");
            value = std::strtod(tok, &ptr);
          }
        } while (ptr == tok);
      }

      iter_type i = 0;
      while (ptr > tok) {
        b.set_parameter(value, i++);
        if ((tok = std::strtok(NULL, " ")) != NULL) {
          value = std::strtod(tok, &ptr);
        } else {
          ptr = tok;
        }
      }
    }
  }
}

void DomainSettings::parse_initial_condition_array(
    const char* value, ModelSettings* model_settings) {
  auto inits = symphas::conf::parse_options(value);

  symphas::lib::resize_array(inits.size(), tdata, tdata_len);
  for (iter_type n = 0; n < inits.size(); ++n) {
    parse_initial_condition(inits[n].c_str(), n, model_settings);
  }
}

symphas::init_entry_type get_initial_condition_entry(char* input,
                                                     size_t dimension,
                                                     double* coeff,
                                                     size_t coeff_len,
                                                     bool& init_coeff_copied) {
  symphas::init_entry_type init;

  // char* pos0 = input;
  char* tok = input;
  while (*++tok && *tok != ' ')
    ;

  char buffer[BUFFER_LENGTH]{0};
  char name[BUFFER_LENGTH]{0};
  std::copy(input, tok, name);
  std::copy(tok, input + std::strlen(input) + 1, buffer);

  init.in = symphas::in_from_str(name);

  if (init.in == Inside::EXPRESSION) {
    // char* pos0 = input;
    tok = std::strtok(buffer, " ");
    char* expression_name = new char[std::strlen(tok) + 1];

    if (sscanf(tok, "%s", expression_name) < 1) {
      fprintf(SYMPHAS_ERR,
              "unable to read the expression name given to the initial "
              "conditions, '%s'\n",
              input);
      exit(841);
    } else {
      init.expr_data = expression_name;
    }

    delete[] expression_name;
    tok = std::strtok(NULL, " ");
  } else {
    switch (init.in) {
      case Inside::GAUSSIAN: {
        init.data.gp[0] = 0;     // mean
        init.data.gp[1] = 0.25;  // std deviation
        break;
      }
      case Inside::UNIFORM: {
        init.data.gp[0] = -1;  // minimum value
        init.data.gp[1] = 1;   // maximum value
        break;
      }
      case Inside::CAPPED: {
        init.data.gp[0] = -1;  // minimum value
        init.data.gp[1] = 1;   // maximum value
        break;
      }
      case Inside::CONSTANT: {
        init.data.gp[0] = 0;  // the constant value
        break;
      }
      case Inside::CIRCLE: {
        init.data.gp[0] = 2;  // ratio of the x width
        init.data.gp[1] = 2;  // ratio of the y width
        if (dimension == 3) {
          init.data.gp[2] = 2;  // ratio of z width
        }
        break;
      }
      case Inside::HEXAGONAL: {
        init.data.gp[0] = 1;  // ratio of density along x
        init.data.gp[1] = 1;  // ratio of density along y
        if (dimension == 2) {
          init.data.gp[2] = 2;  // ratio of size of nodes
        } else if (dimension == 3) {
          init.data.gp[2] = 1;  // ratio of density along z
          init.data.gp[3] = 2;  // ratio of density of nodes
        }
        break;
      }
      case Inside::CUBIC: {
        init.data.gp[0] = 1;  // same as HEX
        init.data.gp[1] = 1;
        if (dimension == 2) {
          init.data.gp[2] = 2;
        } else if (dimension == 3) {
          init.data.gp[2] = 1;
          init.data.gp[3] = 2;
        }
        break;
      }
      case Inside::SEEDSSQUARE: {
        init.data.gp[0] = 10;  // the number of seeds in the system
        init.data.gp[1] = 1;   // the scale factor of the size
        // the value inside the seed
        init.data.gp[2] = params::init_inside_val;
        // the value outside the seed
        init.data.gp[3] = params::init_outside_val;
        break;
      }
      case Inside::SEEDSCIRCLE: {
        init.data.gp[0] = 10;
        init.data.gp[1] = 1;
        init.data.gp[2] = params::init_inside_val;
        init.data.gp[3] = params::init_outside_val;
        break;
      }
      case Inside::VORONOI: {
        init.data.gp[0] = 10;  // number of crystals
        init.data.gp[1] = -1;  // lower range
        init.data.gp[2] = 1;   // upper range
        break;
      }
      case Inside::BUBBLE: {
        init.data.gp[0] = 1;    // The number of bubbles to fill
        init.data.gp[1] = -1;   // lower range
        init.data.gp[2] = 1;    // upper range
        init.data.gp[3] = .75;  // The filling ratio
        init.data.gp[4] = 1;    // THe overlap ratio
        break;
      }
      case Inside::SPIRALHEX: {
        init.data.gp[0] = 1;   // The number of hexes to fill
        init.data.gp[1] = -1;  // lower range
        init.data.gp[2] = 1;   // upper range
        init.data.gp[3] = 1;   // size of the circle
        break;
      }
      default:
        init.data = symphas::init_data_parameters::one();
        break;
    }

    input = tok;
    while (true) {
      if (*tok == ' ')
        while (*++tok == ' ')
          ;
      while (*tok != ' ' && *++tok)
        ;

      std::copy(input, tok, buffer);
      buffer[tok - input] = '\0';

      symphas::lib::str_trim(buffer);
      InsideTag tag = symphas::in_tag_from_str(buffer);

      if (tag != InsideTag::NONE) {
        init.intag = symphas::build_intag(init.intag, tag);
        input = tok;
      } else {
        break;
      }
    }
  }

  size_t gp_count = 0;
  if (tok) {
    if (*tok == CONFIG_OPTION_PREFIX_C) {
      if (init.in == Inside::EXPRESSION) {
        init.expr_data.set_coeff(coeff, coeff_len);
      } else {
        std::copy(coeff,
                  coeff + std::min(coeff_len, size_t(NUM_INIT_CONSTANTS)),
                  init.data.gp);

        if (coeff_len > NUM_INIT_CONSTANTS) {
          fprintf(SYMPHAS_WARN,
                  "using only %d numeric arguments from model coefficients\n",
                  NUM_INIT_CONSTANTS);
        }
      }
      init_coeff_copied = true;
    } else {
      double* init_coeff = nullptr;
      symphas::internal::parse_simple_coeff_list(input, init_coeff, gp_count);

      if (init.in == Inside::EXPRESSION) {
        init.expr_data.set_coeff(init_coeff, gp_count);
      } else {
        std::copy(init_coeff,
                  init_coeff + std::min(gp_count, size_t(NUM_INIT_CONSTANTS)),
                  init.data.gp);
        if (gp_count > NUM_INIT_CONSTANTS) {
          fprintf(SYMPHAS_WARN,
                  "only %d numeric arguments may be provided "
                  "to the initial conditions\n",
                  NUM_INIT_CONSTANTS);
        }
      }

      delete[] init_coeff;
    }
  }

  auto message_unset = [&](size_t m) {
    if (gp_count < m) {
      fprintf(SYMPHAS_LOG,
              "using fewer than possible number of parameters for the initial "
              "condition '%s', the remaining %zd parameters are given the "
              "default values: \n",
              name, m - gp_count);
      for (size_t i = gp_count; i < m;
           fprintf(SYMPHAS_LOG, "\t%.2f \n", init.data.gp[i++]))
        ;
    } else if (gp_count > m) {
      fprintf(SYMPHAS_LOG,
              "more initial condition parameters than required, the initial "
              "condition '%s' uses up to %zd parameters\n",
              name, m);
    } else {
    }
  };

  /* if there are fewer parameters than what is given,
   */

  switch (init.in) {
    case Inside::GAUSSIAN: {
      message_unset(2);
      break;
    }
    case Inside::UNIFORM: {
      message_unset(2);
      break;
    }
    case Inside::CAPPED: {
      message_unset(2);
      break;
    }
    case Inside::CONSTANT: {
      message_unset(1);
      break;
    }
    case Inside::CIRCLE: {
      if (dimension == 2) {
        message_unset(2);
      } else if (dimension == 3) {
        message_unset(3);
      }
      break;
    }
    case Inside::HEXAGONAL: {
      if (dimension == 2) {
        message_unset(3);
      } else if (dimension == 3) {
        message_unset(4);
      }
      break;
    }
    case Inside::CUBIC: {
      if (dimension == 2) {
        message_unset(3);
      } else if (dimension == 3) {
        message_unset(4);
      }
      break;
    }
    case Inside::SEEDSSQUARE: {
      message_unset(4);
      break;
    }
    case Inside::SEEDSCIRCLE: {
      message_unset(4);
      break;
    }
    case Inside::VORONOI: {
      message_unset(3);
      break;
    }
    case Inside::BUBBLE: {
      message_unset(5);
      break;
    }
    case Inside::SPIRALHEX: {
      message_unset(4);
      break;
    }
    default:
      break;
  }

  return init;
}

void DomainSettings::parse_initial_condition_expression(
    char* input, size_t n, ModelSettings* model_settings) {
  Axis ax = Axis::NONE;

  // go to the next exclamation mark
  char* input0 = input;
  char* last = input0;

  do {
    while (*++last && *last != CONFIG_TITLE_PREFIX_C)
      ;

    if (*input0 == CONFIG_TITLE_PREFIX_C) {
      char ax_str[6]{0};
      size_t nn = sscanf(input0 + 1, "%5s", ax_str);
      input0 += nn + 1;

      if (std::strcmp(ax_str, "X") == 0) {
        ax = Axis::X;
      } else if (std::strcmp(ax_str, "Y") == 0) {
        ax = Axis::Y;
      } else if (std::strcmp(ax_str, "Z") == 0) {
        ax = Axis::Z;
      } else if (std::strcmp(ax_str, "R") == 0) {
        ax = Axis::R;
      } else if (std::strcmp(ax_str, "T") == 0) {
        ax = Axis::T;
      } else if (std::strcmp(ax_str, "S") == 0) {
        ax = Axis::S;
      } else {
        fprintf(SYMPHAS_ERR,
                "incorrect axis specification for initial condition, '%s'",
                ax_str);
      }
    }
    char c = *last;
    *last = '\0';

    symphas::lib::str_trim(input0);
    tdata[n][ax] = get_initial_condition_entry(
        input0, intervals[n].size(), model_settings->coeff,
        model_settings->coeff_len, model_settings->init_coeff_copied);

    *last = c;
    input0 = last;

  } while (*input0);

  if (tdata[n].find(Axis::NONE) == tdata[n].end()) {
    tdata[n][Axis::NONE] = tdata[n].begin()->second;
  }
}

void DomainSettings::parse_initial_condition_file(char* input, size_t n) {
  char* file_name = new char[std::strlen(input) + 1];
  int index = 0;

  if (sscanf(input, "%s %d", file_name, &index) < 2) {
    fprintf(SYMPHAS_WARN, "using index %d from intial conditions input file\n",
            index);
    if (sscanf(input, "%s", file_name) < 1) {
      fprintf(SYMPHAS_WARN,
              "unable to read the filename given to the initial conditions, "
              "'%s'\n",
              input);
    }
  }

  if (file_name) {
    tdata[n][Axis::NONE].file = {file_name, index};
    tdata[n][Axis::NONE].in = Inside::FILE;
  }

  delete[] file_name;
}

void DomainSettings::parse_initial_condition_params(char* input, size_t n) {
  if (params::input_data_file) {
    auto inits = symphas::conf::parse_options(params::input_data_file, ",");
    tdata[n][Axis::NONE].file = {inits[n % inits.size()].c_str(), 0};
    tdata[n][Axis::NONE].in = Inside::FILE;
  } else {
    tdata[n][Axis::NONE] = Inside::CONSTANT;
  }
}

void DomainSettings::parse_initial_condition(const char* value, size_t n,
                                             ModelSettings* model_settings) {
  char input[BUFFER_LENGTH];
  std::strncpy(input, value, sizeof(input) / sizeof(char) - 1);
  symphas::lib::str_trim(input);

  if (*input == CONFIG_TITLE_PREFIX_C && std::strlen(input) == 1) {
    parse_initial_condition_params(input, n);
  } else if (*input == CONFIG_OPTION_PREFIX_C) {
    char* tok = std::strtok(input, STR(CONFIG_OPTION_PREFIX));
    if (tok) {
      parse_initial_condition_file(tok, n);
    }
  } else {
    parse_initial_condition_expression(input, n, model_settings);
  }
}

void DomainSettings::parse_interval_array(const char* const* r_specs,
                                          size_t dimension) {
  for (iter_type d = 0; d < dimension; ++d) {
    Axis ax = symphas::index_to_axis(d);
    const char* r_spec = r_specs[d];

    size_t dims_len = intervals_len;
    auto ranges = symphas::conf::parse_options(r_spec);
    if (ranges.size() > intervals_len) {
      symphas::lib::resize_array(ranges.size(), intervals, intervals_len);
      symphas::lib::resize_array(ranges.size(), dims, dims_len, dimension);
    }

    for (iter_type i = 0; i < intervals_len; ++i) {
      parse_interval(ranges[i].c_str(), ax, i);
    }
  }

  for (iter_type i = 0; i < intervals_len; ++i) {
    set_dimensions(i);
  }
}

void DomainSettings::parse_interval_array(Axis ax, const char* r_spec,
                                          size_t dimension) {
  auto ranges = symphas::conf::parse_options(r_spec);
  size_t dims_len = intervals_len;
  if (ranges.size() > intervals_len) {
    symphas::lib::resize_array(ranges.size(), intervals, intervals_len);
    symphas::lib::resize_array(ranges.size(), dims, dims_len, dimension);
  }

  for (iter_type i = 0; i < intervals_len; ++i) {
    parse_interval(ranges[i].c_str(), ax, i);
    set_dimensions(i);
  }
}

void DomainSettings::parse_interval(const char* value, Axis ax, size_t n) {
  if (value[0] == CONFIG_OPTION_PREFIX_C) {
    iter_type count = 1;
    iter_type stop = 0;

    if (sscanf(value, STR(CONFIG_OPTION_PREFIX) " %d%n", &count, &stop) != 1) {
      fprintf(SYMPHAS_WARN,
              "the grid points given by '%s' is not in the correct "
              "format, expecting the number of grid points as `@ L`\n",
              value);
    }

    intervals[n][ax].set_count_from_r(count, intervals[n][ax].width(), 0.0);

    const char* split = strchr(value, CONFIG_TITLE_PREFIX_C);
    if (split != NULL) {
      iter_type left, right;
      if (sscanf(split + 1, "%d %d", &left, &right) != 2) {
        fprintf(
            SYMPHAS_WARN,
            "the subinterval specification given by '%s' is not in the correct "
            "format, expecting the format to be provided as `@ L ! A B`\n",
            split + 1);
      }
      intervals[n][ax].set_interval_fraction((double)left / count,
                                             (double)right / count);
    } else {
      intervals[n][ax].interval_to_domain();
    }
  } else {
    double left = 0, right = 1;
    if (sscanf(value, "%lf %lf", &left, &right) != 2) {
      fprintf(SYMPHAS_WARN,
              "the interval given by '%s' is not in the correct "
              "format, expecting `A B`\n",
              value);
    }
    intervals[n][ax].set_domain(left, right, intervals[n][ax].width());
    intervals[n][ax].interval_to_domain();
  }
}

void NameSettings::parse_names(const char* value) {
  auto name_list = symphas::conf::parse_options(value, true);

  if (name_list.size() > 0) {
    symphas::lib::resize_array(name_list.size(), names, names_len);

    char** name_it = names;
    for (auto&& name : name_list) {
      *name_it = new char[std::strlen(name.c_str()) + 1];
      std::strcpy(*name_it++, name.c_str());
    }
  }
}

symphas::lib::array_container<symphas::lib::string> NameSettings::get_names(
    len_type desired_len) const {
  desired_len = (desired_len < 1) ? len_type(names_len) : desired_len;
  symphas::lib::array_container<symphas::lib::string> name_list(desired_len);

  iter_type i = 0;
  for (auto& name : name_list) {
    if (i >= names_len) {
      name = symphas::lib::string(STR_ARR_LEN(DEFAULT_FIELD_NAME) +
                                  symphas::lib::num_digits(i));
      sprintf(name.begin(), DEFAULT_FIELD_NAME "%d", i);
    } else {
      name = symphas::lib::string(get_name(i), std::strlen(get_name(i)) + 1);
    }
    ++i;
  }
  return name_list;
}

symphas::lib::string NameSettings::get_name(int n) const {
  if (n < names_len) {
    return names[n];
  } else {
    return "";
  }
}

void DirectorySettings::set_directory(const char* directory,
                                      const char* title) {
  delete[] result_dir;
  delete[] root_dir;

  root_dir = new char[std::strlen(directory) + 1];
  std::strcpy(root_dir, directory);

  char append_dir[BUFFER_LENGTH_R1];
  if (title != nullptr) {
    char title_dir[BUFFER_LENGTH_R2];
    symphas::lib::to_file_name(title, title_dir, BUFFER_LENGTH_R2);
    snprintf(append_dir, BUFFER_LENGTH_R1, "%s/%s", root_dir, title_dir);
  } else {
    std::strcpy(append_dir, root_dir);
  }

  // if we're using a timestamp, then the created directory using that timestamp
  // must be unique, add an index to dir if one already exists
  if (params::use_timestamp) {
    /* add timestamp to the directory tree
     */

    char ts_buffer[BUFFER_LENGTH_R2];

#ifdef USING_MPI
    MPI_Request* rqst = new MPI_Request[symphas::parallel::get_num_nodes()]{};

    if (symphas::parallel::is_host_node()) {
      symphas::lib::write_ts_str(ts_buffer);
      for (iter_type i = 0; i < symphas::parallel::get_num_nodes(); ++i) {
        if (!symphas::parallel::is_host_node(i)) {
          rqst[i] = MPI_REQUEST_NULL;
          MPI_Send(&ts_buffer[0], BUFFER_LENGTH_R2, MPI_CHAR, i, 0,
                   MPI_COMM_WORLD);
        }
      }
    } else {
      MPI_Status status;
      MPI_Recv(&ts_buffer[0], BUFFER_LENGTH_R2, MPI_CHAR, SYMPHAS_MPI_HOST_RANK,
               0, MPI_COMM_WORLD, &status);
    }
#else
    symphas::lib::write_ts_str(ts_buffer);
#endif

    size_t dlen = std::strlen(append_dir) + std::strlen(ts_buffer);

#ifdef USING_MPI
    if (symphas::parallel::is_host_node()) {
      result_dir = new char[dlen + 2];
      sprintf(result_dir, "%s/%s", append_dir, ts_buffer);
    } else {
      int rank = symphas::parallel::get_node_rank();
      result_dir = new char[dlen + symphas::lib::num_digits(rank) + 3];
      sprintf(result_dir, "%s/%s/%d", append_dir, ts_buffer, rank);
    }
#else
    result_dir = new char[dlen + 2];
    sprintf(result_dir, "%s/%s", append_dir, ts_buffer);
#endif

    /*
     * decide if index should be appended
     *
     * if multiple runs are initiated at the same time, then the program will
     * detect that the directory with that timestamp name already exists and it
     * will append an index to the end
     */

    iter_type i = 0;
    struct stat info;
    while (stat(result_dir, &info) == 0 && (info.st_mode & S_IFMT) == S_IFDIR &&
           i < 10) {
#ifdef USING_MPI
      char* new_dir = nullptr;
      if (symphas::parallel::is_host_node()) {
        char* new_dir = new char[dlen + STR_ARR_LEN(TIMESTAMP_ID_APPEND) + 2]{};
        sprintf(new_dir, "%s/%s" TIMESTAMP_ID_APPEND "", append_dir, ts_buffer,
                ++i);
      } else {
        int rank = symphas::parallel::get_node_rank();
        char* new_dir = new char[dlen + symphas::lib::num_digits(rank) +
                                 STR_ARR_LEN(TIMESTAMP_ID_APPEND) + 2];
        sprintf(new_dir, "%s/%s" TIMESTAMP_ID_APPEND "/%d", append_dir,
                ts_buffer, ++i, rank);
      }
#else
      char* new_dir = new char[dlen + STR_ARR_LEN(TIMESTAMP_ID_APPEND) + 2]{};
      sprintf(new_dir, "%s/%s" TIMESTAMP_ID_APPEND, append_dir, ts_buffer, ++i);
#endif

      std::swap(result_dir, new_dir);
      delete[] new_dir;
    }
  } else {
    size_t dlen = std::strlen(append_dir);
#ifdef USING_MPI
    if (symphas::parallel::is_host_node()) {
      result_dir = new char[dlen + 1];
      sprintf(result_dir, "%s", append_dir);
    } else {
      int rank = symphas::parallel::get_node_rank();
      result_dir = new char[dlen + symphas::lib::num_digits(rank) + 2];
      sprintf(result_dir, "%s/%d", append_dir, rank);
    }
#else
    result_dir = new char[dlen + 1];
    std::strcpy(result_dir, append_dir);
#endif
  }
}

void SimulationSettings::parse_stencil(size_t order, const char* str) {
  const char msg[] =
      "derivative order '%zd' is invalid in selection of stencil point "
      "values\n";
  if (*str == CONFIG_OPTION_PREFIX_C) {
    StencilParams default_stp = DefaultStencil{dimension, stp.ord}();
    if (default_stp != StencilParams{}) {
      switch (order) {
        case 2:
          stp.ptl = default_stp.ptl;
          break;
        case 3:
          stp.ptg = default_stp.ptg;
          break;
        case 4:
          stp.ptb = default_stp.ptb;
          break;
        default:
          fprintf(SYMPHAS_WARN, msg, order);
      }
    } else {
      fprintf(SYMPHAS_WARN,
              "not able to select a default stencil with the given "
              "dimension (%zd) and order of accuracy (%hu)\n",
              dimension, stp.ord);
    }
  } else {
    unsigned short value =
        static_cast<unsigned short>(std::strtoul(str, NULL, 10));
    if (value > 0 && errno != ERANGE) {
      switch (order) {
        case 2:
          stp.ptl = value;
          break;
        case 3:
          stp.ptg = value;
          break;
        case 4:
          stp.ptb = value;
          break;
        default:
          fprintf(SYMPHAS_WARN, msg, order);
      }
    }
  }
}

void SimulationSettings::parse_stencil_accuracy(const char* str) {
  if (*str == CONFIG_OPTION_PREFIX_C) {
    StencilParams default_stp = DefaultStencil{dimension, 0}();
    if (default_stp != StencilParams{}) {
      stp.ord = default_stp.ord;
    } else {
      fprintf(SYMPHAS_WARN,
              "not able to select a default stencil "
              "accuracy with the given dimension (%zd) \n",
              dimension);
    }
  } else {
    unsigned short value =
        static_cast<unsigned short>(std::strtoul(str, NULL, 10));
    stp.ord = value;
  }
}

void SimulationSettings::parse_solver_type(const char* str) {
  try {
    stp.type = std::stoi(str);
  } catch (const std::invalid_argument&) {
    fprintf(SYMPHAS_ERR,
            "could not convert solver variation configuration '%s'\n", str);
    stp.type = 0;
  }
}

void ModelSettings::parse_model_spec(const char* value, const char* dir) {
  char type[BUFFER_LENGTH + 1]{};
  int model_end_pos;

  if (std::strlen(value) > 1 && value[0] == CONFIG_TITLE_PREFIX_C) {
    iter_type start, end;
    size_t n =
        sscanf(value, CONFIG_TITLE_PREFIX " %" STR(BUFFER_LENGTH) "s %d %d",
               type, &start, &end);

    if (n == 0) {
      fprintf(SYMPHAS_ERR,
              "expected file name after '%c' in model specification\n",
              CONFIG_TITLE_PREFIX_C);
    } else {
      size_t type_len = std::strlen(type);
      size_t keyword_len = sizeof(STR(VIRTUAL_MODEL_KEYWORD)) / sizeof(char);
      size_t index_len = size_t(n > 1) * symphas::lib::num_digits(start) +
                         size_t(n > 2) * symphas::lib::num_digits(end) +
                         (n - 1);
      char* model_spec = new char[keyword_len + type_len + 1 + index_len]{};

      sprintf(model_spec, "%s%c%s", STR(VIRTUAL_MODEL_KEYWORD),
              VIRTUAL_MODEL_SEP_KEY_C, type);
      if (n > 1) {
        sprintf(model_spec + keyword_len + type_len, "%c%d",
                VIRTUAL_MODEL_SEP_KEY_C, start);
      }
      if (n > 2) {
        sprintf(model_spec + keyword_len + type_len +
                    symphas::lib::num_digits(start) + 1,
                "%c%d", VIRTUAL_MODEL_SEP_KEY_C, end);
      }

      std::swap(model, model_spec);
      delete[] model_spec;
    }
  } else {
    // set the model and list of coefficients to read them, the
    // default is given by #DEFAULT_COEFF_VALUE.
    size_t n =
        sscanf(value, " %" STR(BUFFER_LENGTH) "s %n", type, &model_end_pos);

    if (n == 1) {
      set_model_name(type);

      if (model_end_pos > std::strlen(model)) {
        const char* iter = value + model_end_pos - 1;
        while (*++iter && *iter != CONFIG_TITLE_PREFIX_C)
          ;

        if (*iter == CONFIG_TITLE_PREFIX_C) {
          auto options = symphas::conf::parse_options(iter + 1, ",");
          delete[] num_fields;
          for (iter_type i = 0; i < num_fields_len; ++i) {
            delete[] modifiers[i];
          }
          delete[] modifiers;

          num_fields_len = options.size();
          num_fields = new len_type[num_fields_len]{};
          modifiers = new char* [num_fields_len] {};

          iter_type i = 0;
          for (auto option : options) {
            *type = '\0';
            sscanf(option.c_str(), "%d %s", num_fields + i, type);
            if (*type) {
              modifiers[i] = new char[std::strlen(type) + 1];
              std::strcpy(modifiers[i], type);
            } else {
              modifiers[i] = nullptr;
            }
            ++i;
          }
        }

        std::copy(value + model_end_pos, iter, type);
        type[iter - (value + model_end_pos)] = '\0';
        symphas::internal::parse_coeff_list(type, dir, coeff, coeff_len);
      } else {
        if (coeff) {
          coeff_len = 0;
        }
      }
    } else {
      fprintf(SYMPHAS_ERR,
              "model name was not specified at the configuration "
              "specification for %s\n",
              symphas::internal::C_MODEL);
      model = nullptr;
    }
  }
}

void DomainSettings::parse_width(const char* str, size_t dimension) {
  auto widths = symphas::conf::parse_options(str);
  if (widths.size() > intervals_len) {
    size_t dims_len = intervals_len;
    symphas::lib::resize_array(widths.size(), intervals, intervals_len);
    symphas::lib::resize_array(widths.size(), dims, dims_len, dimension);
  }

  for (iter_type i = 0; i < widths.size(); ++i) {
    auto values = symphas::conf::parse_options(widths[i].c_str(), true);

    iter_type len = static_cast<iter_type>(values.size());
    if (len > 0) {
      for (iter_type d = 0; d < len; ++d) {
        char* end;
        double width = strtod(values[d].c_str(), &end);

        if (*end) {
          fprintf(SYMPHAS_WARN,
                  "there was an error interpreting the width "
                  "in '%s'",
                  widths[i].c_str());
          width = 1.0;
        }

        Axis ax = symphas::index_to_axis(d);
        if (intervals[i].count(ax) > 0) {
          auto& interval = intervals[i].at(ax);
          interval.set_domain(interval.left(), interval.right(), width);
        } else {
          intervals[i][ax] = symphas::interval_element_type(0, width, width);
        }
      }
      for (iter_type d = len; d < dimension; ++d) {
        auto& interval0 = intervals[i].at(symphas::index_to_axis(0));
        Axis ax = symphas::index_to_axis(d);
        if (intervals[i].count(ax) > 0) {
          auto& interval = intervals[i].at(symphas::index_to_axis(d));
          interval.set_domain(interval.left(), interval.right(),
                              interval0.width());
        } else {
          intervals[i][ax] = symphas::interval_element_type(
              0, interval0.width(), interval0.width());
        }
      }
    }
  }
}

void SimulationSettings::parse_dt(const char* str) {
  auto dts = symphas::conf::parse_options(str, ",");
  if (dts.size() > 0) {
    if (dts.size() == 1) {
      double dt;
      if (sscanf(dts.front().c_str(), "%lf", &dt) != 1) {
        fprintf(SYMPHAS_WARN,
                "ignoring the time step specification '%s' that is not in the "
                "correct "
                "format, expected `dt`\n",
                dts.front().c_str());
      } else {
        dt_list.set_time_step(dt, 0, false);
      }
    } else {
      dt_list.clear_time_steps(1.0);
      double dt, time;
      for (const auto& spec : dts) {
        if (sscanf(spec.c_str(), "%lf " STR(CONFIG_OPTION_PREFIX) " %lf", &dt,
                   &time) != 2) {
          fprintf(SYMPHAS_WARN,
                  "ignoring the time step specification '%s' that is not in "
                  "the correct "
                  "format, expected `dt`\n",
                  spec.c_str());
        } else {
          dt_list.set_time_step(dt, time);
        }
      }
    }
  }
}

void ModelSettings::set_model_name(const char* str) {
  delete[] model;
  model = new char[std::strlen(str) + 1];
  symphas::lib::to_upper(str, model);
}

SymPhasSettings::SymPhasSettings(settings_spec_type const& settings_list,
                                 const char* title, const char* dir)
    : SymPhasSettings() {
  std::map<symphas::lib::string, symphas::lib::string> settings_map;
  for (auto setting : settings_list) {
    auto& [key, value] = setting;
    settings_map[key.c_str()] = handle_substitutions(value.c_str());
  }

  using namespace symphas::internal;
  simulation_settings.parse_dimension(settings_map[C_DIM]);
  simulation_settings.parse_stencil_accuracy(settings_map[C_ORDER]);
  simulation_settings.parse_dt(settings_map[C_DELTA]);
  domain_settings.parse_width(settings_map[C_WIDTH],
                              simulation_settings.dimension);

  if (settings_map.find(C_SOLVERVAR) != settings_map.end()) {
    simulation_settings.parse_solver_type(settings_map[C_SOLVERVAR]);
  }

  for (auto [key, order] : {std::make_pair(C_PTL, 2), std::make_pair(C_PTG, 3),
                            std::make_pair(C_PTB, 4)}) {
    if (settings_map.find(key) != settings_map.end())
      simulation_settings.parse_stencil(order, settings_map[key]);
    else
      simulation_settings.parse_stencil(order, STR(CONFIG_OPTION_PREFIX));
  }

  for (auto& [key, ax] :
       {std::make_pair(C_RNGX, Axis::X), std::make_pair(C_RNGY, Axis::Y),
        std::make_pair(C_RNGZ, Axis::Z)}) {
    if (settings_map.find(key) != settings_map.end()) {
      domain_settings.parse_interval_array(ax, settings_map[key],
                                           simulation_settings.dimension);
    }
  }

  if (domain_settings.get_dimension() != simulation_settings.dimension) {
    fprintf(SYMPHAS_ERR,
            "the dimension of the domain and the simulation "
            "settings do not match\n");
    throw;
  }

  for (auto& [key, side] : {std::make_pair(C_BNDLT, Side::LEFT),
                            std::make_pair(C_BNDRT, Side::RIGHT),
                            std::make_pair(C_BNDTP, Side::TOP),
                            std::make_pair(C_BNDBT, Side::BOTTOM),
                            std::make_pair(C_BNDFT, Side::FRONT),
                            std::make_pair(C_BNDBK, Side::BACK)}) {
    if (settings_map.find(key) != settings_map.end()) {
      domain_settings.parse_boundaries_array(side, settings_map[key],
                                             simulation_settings.dimension);
    }
  }

  domain_settings.parse_initial_condition_array(settings_map[C_INSIDE],
                                                &model_settings);
  model_settings.parse_model_spec(settings_map[C_MODEL], dir);

  if (settings_map.find(C_NAMES) != settings_map.end())
    name_settings.parse_names(settings_map[C_NAMES]);
  name_settings.set_title(title);
  simulation_settings.parse_save(settings_map[C_SAVE]);
  simulation_settings.parse_frame_count(settings_map[C_FRAMES]);

  if (settings_map.find(C_SAVEINIT) != settings_map.end())
    simulation_settings.parse_save_init(settings_map[C_SAVEINIT]);

  if (settings_map.find(C_RUNS) != settings_map.end())
    runs = atoi(settings_map[C_RUNS]);

  if (settings_map.find(C_DIR) != settings_map.end())
    directory_settings.set_directory(settings_map[C_DIR], name_settings.title);

  make_directories();
}

symphas::problem_parameters_type SymPhasSettings::get_problem_parameters()
    const {
  size_t sys_len =
      std::max({domain_settings.intervals_len, domain_settings.bdata_len,
                domain_settings.tdata_len});
  symphas::problem_parameters_type pp(0);
  if (model_settings.num_fields_len > 0) {
    len_type field_len =
        std::reduce(model_settings.num_fields,
                    model_settings.num_fields + model_settings.num_fields_len);

    if (sys_len > field_len) {
      len_type* num_fields_expand =
          new len_type[model_settings.num_fields_len + 1];
      char** modifiers_expand =
          new char* [model_settings.num_fields_len + 1] {};
      for (iter_type i = 0; i < model_settings.num_fields_len; ++i) {
        num_fields_expand[i] = model_settings.num_fields[i];
        modifiers_expand[i] =
            new char[std::strlen(model_settings.modifiers[i]) + 1]{};
        std::strcpy(modifiers_expand[i], model_settings.modifiers[i]);
      }
      num_fields_expand[model_settings.num_fields_len] =
          len_type(sys_len) - field_len;
      modifiers_expand[model_settings.num_fields_len] = nullptr;

      pp = symphas::problem_parameters_type(num_fields_expand,
                                            model_settings.num_fields_len + 1);
      pp.set_modifiers(modifiers_expand, model_settings.num_fields_len + 1);
      delete[] num_fields_expand;
      delete[] modifiers_expand;

    } else {
      pp = symphas::problem_parameters_type(model_settings.num_fields,
                                            model_settings.num_fields_len);
      pp.set_modifiers(model_settings.modifiers, model_settings.num_fields_len);
    }
  } else {
    pp = symphas::problem_parameters_type(sys_len);
  }

  pp.set_boundary_data(domain_settings.bdata, domain_settings.bdata_len);
  pp.set_initial_data(domain_settings.tdata, domain_settings.tdata_len);
  pp.set_interval_data(domain_settings.intervals,
                       domain_settings.intervals_len);
  pp.set_time_steps(simulation_settings.dt_list);

  return pp;
}

inline void print_model_coeff(char* out, const double* coeff, iter_type i) {
  if (coeff[i] * 100 == int(coeff[i] * 100)) {
    sprintf(out, "[%d]=%.2lf ", i, coeff[i]);
  } else if (coeff[i] / 100. >= 1. || 1. / coeff[i] >= 100.) {
    sprintf(out, "[%d]=%" DATA_OUTPUT_ACCURACY_STR "E ", i, coeff[i]);
  } else {
    sprintf(out, "[%d]=%" DATA_OUTPUT_ACCURACY_STR "lf ", i, coeff[i]);
  }
}

inline void print_coeff(FILE* out, const double* coeff, iter_type i) {
  if (coeff[i] == int(coeff[i])) {
    fprintf(out, "%d ", (int)coeff[i]);
  } else if (coeff[i] / 100. >= 1. || 1. / coeff[i] >= 100.) {
    fprintf(out, "%" DATA_OUTPUT_ACCURACY_STR "E ", coeff[i]);
  } else {
    fprintf(out, "%" DATA_OUTPUT_ACCURACY_STR "lf ", coeff[i]);
  }
}

void SymPhasSettings::write(const char* savedir, const char* name) const {
  char* model_name;

  if (model_settings.model != nullptr) {
    auto sep_it = std::strchr(model_settings.model, VIRTUAL_MODEL_SEP_KEY_C);

    if (sep_it == NULL) {
      model_name = new char[std::strlen(model_settings.model) + 1]{};
      std::strcpy(model_name, model_settings.model);
    } else {
      model_name = new char[sep_it - model_settings.model + 1]{};
      std::copy(model_settings.model, sep_it, model_name);
      model_name[sep_it - model_settings.model] = '\0';
    }
  } else {
    model_name = new char[1]{};
  }

  char param_file[BUFFER_LENGTH]{};
  snprintf(param_file, BUFFER_LENGTH, "%s/%s.constants", savedir, model_name);

  FILE* pf;
  if ((pf = fopen(param_file, "a")) == 0) {
    fprintf(SYMPHAS_ERR,
            "error opening write configuration parameters file '%s'\n",
            param_file);
    exit(ERR_CODE_FILE_OPEN);
  }

  char model_spec[BUFFER_LENGTH_L4]{};
  snprintf(model_spec, BUFFER_LENGTH_L4, "%s ", model_name);
  delete[] model_name;

  // print to the parameters file and the model specification line
  for (iter_type i = 0; i < model_settings.coeff_len; ++i) {
    char coeff_spec[BUFFER_LENGTH_R2];
    print_model_coeff(coeff_spec, model_settings.coeff, i);

    fprintf(pf, "%s\n", coeff_spec);
    if (std::strlen(model_spec) + std::strlen(coeff_spec) < BUFFER_LENGTH_L4) {
      std::strcat(model_spec, coeff_spec);
    } else {
      fprintf(
          SYMPHAS_WARN,
          "printing configuration could not include coefficient index '%d'\n",
          i);
    }
  }
  fclose(pf);

  if (model_settings.num_fields_len > 0) {
    std::strcat(model_spec, CONFIG_TITLE_PREFIX);
    for (iter_type i = 0; i < model_settings.num_fields_len; ++i) {
      char modifier_spec[BUFFER_LENGTH_R2];
      if (i == model_settings.num_fields_len - 1) {
        sprintf(modifier_spec, " %d %s", model_settings.num_fields[i],
                model_settings.modifiers[i]);
      } else {
        sprintf(modifier_spec, " %d %s,", model_settings.num_fields[i],
                model_settings.modifiers[i]);
      }
      std::strcat(model_spec, modifier_spec);
    }
  }

  FILE* f;
  char bname[BUFFER_LENGTH]{};
  sprintf(bname,
          "%s/%s"
          "." CONFIG_EXTENSION,
          savedir, name);
  if ((f = fopen(bname, "w")) == 0) {
    fprintf(SYMPHAS_ERR, "error opening write configuration file '%s'\n",
            bname);
    exit(ERR_CODE_FILE_OPEN);
  }

  if (std::strlen(name_settings.title) > 0) {
    fprintf(f, "%c%s\n", CONFIG_TITLE_PREFIX_C, name_settings.title);
  }

  char open = symphas::internal::option_open_bracket();
  char close = symphas::internal::option_close_bracket();

  fprintf(f, CONFIG_NAME_FMT "%zd\n", symphas::internal::C_DIM,
          simulation_settings.dimension);
  fprintf(f, CONFIG_NAME_FMT "%hu\n", symphas::internal::C_ORDER,
          simulation_settings.stp.ord);
  fprintf(f, CONFIG_NAME_FMT "%hu\n", symphas::internal::C_PTL,
          simulation_settings.stp.ptl);
  fprintf(f, CONFIG_NAME_FMT "%hu\n", symphas::internal::C_PTB,
          simulation_settings.stp.ptb);
  fprintf(f, CONFIG_NAME_FMT "%hu\n", symphas::internal::C_PTG,
          simulation_settings.stp.ptg);
  fprintf(f, CONFIG_NAME_FMT "%s\n", symphas::internal::C_FORM, "");
  fprintf(f, CONFIG_NAME_FMT "%hu\n", symphas::internal::C_SOLVERVAR,
          simulation_settings.stp.type);

  fprintf(f, CONFIG_NAME_FMT, symphas::internal::C_WIDTH);
  for (iter_type i = 0; i < domain_settings.intervals_len; ++i) {
    fprintf(f, "%c ", open);
    for (iter_type d = 0; d < simulation_settings.dimension; ++d) {
      auto& interval =
          domain_settings.intervals[i].at(symphas::index_to_axis(d));
      double width = interval.width();
      fprintf(f, "%lf ", width);
    }
    fprintf(f, "%c", close);
  }
  fprintf(f, "\n");

  std::map<Side, const char*> side_key_map = {
      {Side::LEFT, symphas::internal::C_BNDLT},
      {Side::RIGHT, symphas::internal::C_BNDRT},
      {Side::TOP, symphas::internal::C_BNDTP},
      {Side::BOTTOM, symphas::internal::C_BNDBT},
      {Side::FRONT, symphas::internal::C_BNDFT},
      {Side::BACK, symphas::internal::C_BNDBK}};

  for (iter_type side_index = 0; side_index < simulation_settings.dimension * 2;
       ++side_index) {
    Side side = symphas::index_to_side(side_index);

    fprintf(f, CONFIG_NAME_FMT, side_key_map[side]);
    for (iter_type i = 0; i < domain_settings.bdata_len; ++i) {
      if (domain_settings.bdata[i].find(side) ==
          domain_settings.bdata[i].end()) {
        throw std::runtime_error(
            "The dimension of the domain and the simulation settings do not "
            "match.");
      }
      fprintf(
          f, "%c %s ", open,
          symphas::str_from_boundary(domain_settings.bdata[i].at(side).type));

      for (iter_type tag_index = 0; tag_index < 2; ++tag_index) {
        const char* tag_name = symphas::str_from_boundary_tag(
            domain_settings.bdata[i].at(side).tag[tag_index]);
        if (tag_name) {
          fprintf(f, "%s ", tag_name);
        }
      }
      for (iter_type n = 0; n < domain_settings.bdata[i].at(side).argc; ++n) {
        fprintf(f, "%.8E ", domain_settings.bdata[i].at(side).params[n]);
      }
      fprintf(f, "%c", close);
    }
    fprintf(f, "\n");
  }

  const char interval_fmt[] = "%c %lf %lf %c";
  if (simulation_settings.dimension > 0) {
    fprintf(f, CONFIG_NAME_FMT, symphas::internal::C_RNGX);
    for (iter_type i = 0; i < domain_settings.intervals_len; ++i) {
      fprintf(f, interval_fmt, open, domain_settings.DOMAIN_X0_AT(i),
              domain_settings.DOMAIN_Xn_AT(i), close);
    }
    fprintf(f, "\n");
  }
  if (simulation_settings.dimension > 1) {
    fprintf(f, CONFIG_NAME_FMT, symphas::internal::C_RNGY);
    for (iter_type i = 0; i < domain_settings.intervals_len; ++i) {
      fprintf(f, interval_fmt, open, domain_settings.DOMAIN_Y0_AT(i),
              domain_settings.DOMAIN_Yn_AT(i), close);
    }
    fprintf(f, "\n");
  }
  if (simulation_settings.dimension > 2) {
    fprintf(f, CONFIG_NAME_FMT, symphas::internal::C_RNGZ);
    for (iter_type i = 0; i < domain_settings.intervals_len; ++i) {
      fprintf(f, interval_fmt, open, domain_settings.DOMAIN_Z0_AT(i),
              domain_settings.DOMAIN_Zn_AT(i), close);
    }
    fprintf(f, "\n");
  }

  fprintf(f, CONFIG_NAME_FMT, symphas::internal::C_INSIDE);
  for (iter_type i = 0; i < domain_settings.tdata_len; ++i) {
    fprintf(f, "%c ", open);
    for (auto const& [key, entry] : domain_settings.tdata[i]) {
      if (domain_settings.tdata[i].size() > 1 && key != Axis::NONE) {
        fprintf(f, "%c%s ", CONFIG_TITLE_PREFIX_C,
                (key == Axis::X)   ? "X"
                : (key == Axis::Y) ? "Y"
                : (key == Axis::Z) ? "Z"
                : (key == Axis::R) ? "R"
                : (key == Axis::S) ? "S"
                : (key == Axis::T) ? "T"
                                   : "?");
      }

      if (entry.in != Inside::NONE &&
          ((domain_settings.tdata[i].size() > 1 && key != Axis::NONE) ||
           (domain_settings.tdata[i].size() == 1))) {
        if (entry.in == Inside::EXPRESSION) {
          fprintf(f, "%s %s ", symphas::str_from_in(entry.in),
                  entry.expr_data.get_name());
          if (!model_settings.init_coeff_copied) {
            for (iter_type n = 0; n < entry.expr_data.get_num_coeff(); ++n) {
              print_coeff(f, entry.expr_data.get_coeff(), n);
            }
          } else {
            fprintf(f, STR(CONFIG_OPTION_PREFIX));
          }
        } else if (entry.in == Inside::FILE || entry.in == Inside::CHECKPOINT) {
          fprintf(f, STR(CONFIG_OPTION_PREFIX) " %s %d", entry.file.get_name(),
                  entry.file.get_index());
        } else {
          fprintf(f, "%s ", symphas::str_from_in(entry.in));

          size_t tag = entry.intag;
          size_t pos = 0;
          while (tag > 0) {
            bool is_bit_set = (tag & (1ull << pos)) > 0;
            if (is_bit_set) {
              tag = (tag & ~(1ull << pos));
              const char* tag_name =
                  symphas::str_from_in_tag(static_cast<InsideTag>(pos));
              if (tag_name) {
                fprintf(f, "%s ", tag_name);
              }
            }
            pos += 1;
          }
          for (iter_type n = 0; n < NUM_INIT_CONSTANTS; ++n) {
            print_coeff(f, entry.data.gp, n);
          }
        }
      }
    }

    fprintf(f, "%c", close);
  }
  fprintf(f, "\n");

  fprintf(f, CONFIG_NAME_FMT, symphas::internal::C_DELTA);
  if (simulation_settings.dt_list.get_num_time_steps() > 1) {
    for (auto [time, dt] : simulation_settings.dt_list) {
      fprintf(f, "%c %lf " STR(CONFIG_OPTION_PREFIX) " %lf %c ", open, dt, time,
              close);
    }
    fprintf(f, "\n");
  } else {
    fprintf(f, "%lf\n", simulation_settings.dt_list.get_time_step());
  }

  if (model_settings.model != nullptr) {
    char* sep_it;
    if ((sep_it = std::strchr(model_settings.model, VIRTUAL_MODEL_SEP_KEY_C)) ==
        NULL) {
      fprintf(f, CONFIG_NAME_FMT "%s\n", symphas::internal::C_MODEL,
              model_spec);
    } else {
      char* model_cpy = new char[std::strlen(sep_it)]{};
      auto index_it = std::strchr(sep_it + 1, VIRTUAL_MODEL_SEP_KEY_C);
      auto end_it =
          (index_it == NULL) ? sep_it + std::strlen(sep_it) : index_it;

      std::copy(sep_it + 1, end_it, model_cpy);
      model_cpy[end_it - sep_it - 1] = '\0';
      fprintf(f, CONFIG_NAME_FMT "%c %s", symphas::internal::C_MODEL,
              CONFIG_TITLE_PREFIX_C, model_cpy);

      if (index_it != NULL) {
        iter_type index;
        sscanf(index_it + 1, "%d", &index);
        fprintf(f, " %d", index);
        if ((index_it = std::strchr(index_it + 1, VIRTUAL_MODEL_SEP_KEY_C)) !=
            NULL) {
          sscanf(index_it + 1, "%d", &index);
          fprintf(f, " %d", index);
        }
      }
      fprintf(f, "\n");

      delete[] model_cpy;
    }
  }

  fprintf(f, CONFIG_NAME_FMT "%d\n", symphas::internal::C_FRAMES,
          simulation_settings.save.get_stop());

  char save_spec[BUFFER_LENGTH];
  symphas::io::save_as_str(&simulation_settings.save, save_spec, BUFFER_LENGTH);
  fprintf(f, CONFIG_NAME_FMT "%s\n", symphas::internal::C_SAVE, save_spec);

  fprintf(f, CONFIG_NAME_FMT "%s\n", symphas::internal::C_SAVEINIT,
          simulation_settings.save.get_init_flag() ? "YES" : "NO");
  fprintf(f, CONFIG_NAME_FMT "%zd\n", symphas::internal::C_RUNS, runs);
  fprintf(f, CONFIG_NAME_FMT "%s\n", symphas::internal::C_DIR,
          directory_settings.root_dir);

  fprintf(f, "\n");

  /* also copy into the backup configuration the parameters that were given
   * in the last run
   */
  for (auto s : params::rawparams) {
    if (s.second.size() == 0) {
      fprintf(f, CONFIG_PARAM_PREFIX "%s=\n", s.first.c_str());
    } else {
      fprintf(f, CONFIG_PARAM_PREFIX "%s=%s\n", s.first.c_str(),
              s.second.c_str());
    }
  }
  fprintf(f, "\n");
  fclose(f);
}
