
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

namespace symphas::internal
{



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






	/* utility functions for converting the boundary label into the enum boundary type/tag
	 */

	/* initializes the boundary based on the given specification
	 */
	inline void fill_boundary(symphas::b_element_type& b, const char* spec)
	{
		char typestr[BUFFER_LENGTH_R2], str[BUFFER_LENGTH];
		int end;

		int args = sscanf(spec, "%s%n", typestr, &end);
		
		b.type = boundary_from_str(typestr);
		BoundaryTag tag = boundary_tag_from_str(typestr);

		bool type_as_tag = (b.type == BoundaryType::DEFAULT && tag != BoundaryTag::NONE);
		if (type_as_tag)
		{
			b.tag[0] = tag;
		}

		symphas::lib::str_trim(spec + end, str);
		if (end < std::strlen(spec))
		{
			char* ptr;
			char* tok;

			if ((tok = std::strtok(str, " ")) != NULL)
			{
				double value = std::strtod(tok, &ptr);
				if (ptr == tok)
				{
					iter_type i = (type_as_tag) ? 1 : 0;
					do
					{
						if (tok != NULL)
						{
							if (i < sizeof(b.tag) / sizeof(BoundaryTag))
							{
								b.tag[i++] = boundary_tag_from_str(tok);
							}
							tok = std::strtok(NULL, " ");
							value = std::strtod(tok, &ptr);
						}
					} while (ptr == tok);
				}

				iter_type i = 0;
				while (ptr > tok)
				{
					b.set_parameter(value, i++);
					if ((tok = std::strtok(NULL, " ")) != NULL)
					{
						value = std::strtod(tok, &ptr);
					}
					else
					{
						ptr = tok;
					}
				}
			}
		}


	}



	//! The letter range used to refer to coefficients.
	/*!
	 * Defines the first and last letters which define the range within
	 * the alphabet that are used to refer to coefficients. These are
	 * the same as the coefficient macro definitions named in the
	 * model definitions.
	 */
	constexpr int NAMED_RANGE[] = { 'A', 'H' };

	//! The number of letters in #NAMED_RANGE.
	constexpr iter_type RANGE_LEN = NAMED_RANGE[1] - NAMED_RANGE[0] + 1;

	//! The opening and closing brackets of a coefficient index.
	/*!
	 * Defines the allowable brackets which are used in the grammar of
	 * directly specifying a coefficient index in the coefficient file.
	 */
	constexpr std::initializer_list<std::pair<char, char>> COEFF_ID_BRACKETS = { {'[', ']'}, {'(', ')'} };

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
	bool is_name_coeff_id(const char* name)
	{
		size_t len = std::strlen(name);

		for (auto [openc, closec] : COEFF_ID_BRACKETS)
		{
			if (name[0] == openc && name[len - 1] == closec)
			{
				return true;
			}
		}

		return false;
	}

	//! Check whether the coefficient name contains letters in #NAMED_RANGE.
	/*!
	 * Returns true if each of the characters of name are characters in the
	 * range given by #NAMED_RANGE, which istypically between 'A' and 'H'. The
	 * string is case sensitive (the range is only between the capital letters 'A'
	 * and 'H', and not the small letter equivalents).
	 *
	 * The function coeff_index(const char*) processes the name under similar
	 * requirements (not necessarily identical) to extract the index represented
	 * by the name.
	 *
	 * \param name The string that is checked to satisfy the condition that each
	 * letter belongs to the range specified by #NAMED_RANGE.
	 */
	bool is_name_coeff_alph(const char* name)
	{
		int value_check;
		if (sscanf(name, COEFF_STR_PREPEND_NAME "%d", &value_check) == 1)
		{
			return (value_check > 0);
		}
		else
		{
			for (const char* it = name; *it; ++it)
			{
				if (*it < NAMED_RANGE[0] ||
					*it > NAMED_RANGE[1])
				{
					return false;
				}
			}
			return true;
		}
	}

	//! Check whether the name consists of all the same coefficient letter.
	/*!
	 * The given name must satisfy the condition checked by
	 * is_name_coeff_alph(const char*) and then the additional condition that all
	 * letters in the name are the same.
	 */
	bool is_name_coeff_alph_same(const char* name)
	{
		if (is_name_coeff_alph(name))
		{
			iter_type len = static_cast<iter_type>(std::strlen(name));
			for (size_t i = 1; i < len; ++i)
			{
				if (name[0] != name[i])
				{
					return false;
				}
			}

			return true;
		}
		else
		{
			return false;
		}
	}


	//! Return coefficient index from the name.
	/*!
	 * Given the coefficient name, return the index of the coefficient in a
	 * coefficient list. The coefficients are ordered alphabetically by capital
	 * letters from #NAMED_RANGE[0] to #NAMED_RANGE[1], typically `A` to `H`.
	 *
	 * For a name to refer to a coefficient of index larger than the index specified
	 * by #NAMED_RANGE[1], additional letters from the range are appended. Only
	 * the same letter as the starting letter may be appended, otherwise the name
	 * is not a valid index. The length of the string determines the interval of
	 * index value specified by the name.
	 *
	 * If the name is valid, then the index is related to the letter index `N` of
	 * the name, the length of the name `l` and the number of letters in
	 *  #NAMED_RANGE `dl` by the relation `N + l * dl`.
	 */
	iter_type index_coeff_alph_same(const char* name)
	{
		if (is_name_coeff_alph_same(name))
		{
			iter_type letter_index = static_cast<iter_type>(*name - NAMED_RANGE[0]);
			iter_type len_offset = static_cast<iter_type>(std::strlen(name) - 1) * RANGE_LEN;

			return letter_index + len_offset;
		}
		else
		{
			int value_check;
			if (sscanf(name, COEFF_STR_PREPEND_NAME "%d", &value_check) == 1)
			{
				if (value_check > 0)
				{
					return value_check - 1;
				}
				else
				{
					fprintf(SYMPHAS_WARN, "the value given to '%s' "
						"must be a positive integer\n", name);
					return -1;
				}
			}
			else if (is_name_coeff_alph(name))
			{
				fprintf(SYMPHAS_WARN, "the coefficient name '%s' is invalid "
					"because all characters must match '%c'\n", name, name[0]);
				return -1;
			}
		}
		fprintf(SYMPHAS_WARN, "the given coefficient name '%s' is not "
			"in a valid format\n", name);
		return -1;
	}


	constexpr std::initializer_list<const char*>
		COEFF_ID_RANGE_STRS = { "..", "-" };
	constexpr std::initializer_list<const char*>
		COEFF_ID_SEP_STRS = { "," };

	//! Parse the given string for coefficients in list format
	/*!
	 * Read the coefficient indices from the given string. If `brackets` is `true`,
	 * then the coefficients are provided in brackets and the parser takes this
	 * into account. Otherwise, the list format can be directly parsed.
	 */
	std::vector<iter_type> index_list_coeff_id(const char* name, bool brackets = true)
	{
		std::vector<iter_type> indices;
		if (!brackets || is_name_coeff_id(name))
		{
			char buffer[BUFFER_LENGTH_R4];

			size_t content_len;
			char* content;

			if (brackets)
			{
				// read in the content between the brackets
				content_len = std::strlen(name) - 2;
				content = new char[content_len + 1];
				for (iter_type i = 0; i < content_len; ++i)
				{
					content[i] = name[i + 1];
				}
				content[content_len] = '\0';
			}
			else
			{
				content_len = std::strlen(name);
				content = new char[content_len + 1];
				std::strcpy(content, name);
			}

			iter_type delim_pos = symphas::lib::pos_after_digit(content);
			if (delim_pos > 0)
			{
				if (delim_pos == content_len)
				{
					std::strncpy(buffer, content, sizeof(buffer) / sizeof(char) - 1);
					indices.emplace_back(atoi(buffer));
				}
				else
				{
					char* delim_it = content + delim_pos;
					std::copy(content, delim_it, buffer);
					buffer[delim_pos] = '\0';

					indices.emplace_back(atoi(buffer) - 1);
					iter_type n = symphas::lib::pos_after_token(delim_it, COEFF_ID_RANGE_STRS);
					bool add_interval = (n != 0);

					if (!add_interval)
					{
						n = symphas::lib::pos_after_token(delim_it, COEFF_ID_SEP_STRS);
					}

					auto rest = index_list_coeff_id(delim_it + n, false);
					if (rest.size() == 0)
					{
						std::copy(delim_it, delim_it + n, buffer);
						buffer[n] = '\0';
						fprintf(SYMPHAS_WARN, "the key '%s' in '%s' must "
							"be followed by an index \n", buffer, name);
					}
					else
					{
						iter_type next_index = rest[0];

						// if a range key was specified, all the indices between
						// the currently parsed index and the next parsed index
						if (add_interval)
						{
							iter_type current_index = indices[0];
							indices.reserve(indices.size() + (next_index - current_index));
							for (iter_type i = current_index + 1; i < next_index; ++i)
							{
								indices.emplace_back(i);
							}
						}

						indices.insert(indices.end(), rest.begin(), rest.end());
					}
				}
			}
			else
			{
				fprintf(SYMPHAS_WARN, "the content of the indexed coefficient "
					"specification '%s' must begin with an index\n", name);
			}

			delete[] content;

		}
		else
		{
			fprintf(SYMPHAS_WARN, "the coefficient name '%s' is invalid because it "
				"must be surrounded by the correct brackets\n", name);
		}
		return indices;
	}


	//! Assign the value taken from the value string to the given index
	/*!
	 * The `value_str` is parsed to identify the type of value it represents, which
	 * can be numeric or a name, in which case the value corresponding to that name
	 * is taken instead. The value is then assigned to the given index in the
	 * coefficient array. Depending on where the value string
	 */
	void assign_coeff(const char* value_str, iter_type assign_index, double*& coeff, size_t& coeff_len)
	{
		char coeff_value_str[BUFFER_LENGTH_R2];
		symphas::lib::str_trim(value_str, coeff_value_str);


		if (assign_index >= coeff_len)
		{
			symphas::lib::expand_append_list(DEFAULT_COEFF_VALUE, coeff, assign_index + 1 - coeff_len, coeff_len);
			coeff_len = assign_index + 1;
		}

		iter_type access_index;
		if (is_name_coeff_alph(coeff_value_str))
		{
			access_index = index_coeff_alph_same(coeff_value_str);
		}
		else if (is_name_coeff_id(coeff_value_str))
		{
			auto indices = index_list_coeff_id(coeff_value_str);
			if (indices.size() > 1)
			{
				fprintf(SYMPHAS_WARN, "a coefficient list may not appear on the left "
					"of an equals sign, incorrectly given in '%s'.\n", coeff_value_str);
				access_index = indices[0];
			}
			else if (indices.empty())
			{
				access_index = static_cast<iter_type>(coeff_len);
			}
			else
			{
				access_index = indices[0];
			}
		}
		else
		{
			double value;
			if (sscanf(value_str, "%lf", &value) == 1)
			{
				coeff[assign_index] = value;
				return;
			}
			else
			{
				access_index = static_cast<iter_type>(coeff_len);
			}
		}


		if (access_index >= coeff_len)
		{
			fprintf(SYMPHAS_WARN, "the coefficient '%d' is being set equal to the undefined "
				"coefficient '%s', so will instead be set equal to %.2f.\n",
				assign_index, value_str, DEFAULT_COEFF_VALUE);
			coeff[assign_index] = DEFAULT_COEFF_VALUE;
		}
		else if (access_index >= 0)
		{
			coeff[assign_index] = coeff[access_index];
		}
	}

	// parse one line from a coefficient initialization file
	void coeff_read_line(const char* line, double*& coeff, size_t& coeff_len)
	{
		char* line_buffer = new char[std::strlen(line) + 1];
		std::strcpy(line_buffer, line);

		symphas::lib::str_trim(line_buffer);

		// check if the line is not a comment, or empty
		if (*line_buffer && *line_buffer != CONFIG_COMMENT_PREFIX_C)
		{
			char coeff_name[BUFFER_LENGTH];
			int end; // stores the position of the last scanned character
			size_t n = sscanf(line_buffer, "%[^=]=%n", coeff_name, &end);

			if (n == 1)
			{
				symphas::lib::str_trim(coeff_name);
				if (is_name_coeff_alph(coeff_name))
				{
					iter_type index = index_coeff_alph_same(coeff_name);
					if (index >= 0)
					{
						assign_coeff(line_buffer + end, index, coeff, coeff_len);
					}
				}
				else if (is_name_coeff_id(coeff_name))
				{
					auto indices = index_list_coeff_id(coeff_name);
					for (auto index : indices)
					{
						assign_coeff(line_buffer + end, index, coeff, coeff_len);
					}
				}
				else
				{
					fprintf(SYMPHAS_WARN, "coefficient name '%s' to assign isn't recognized,"
						" line is skipped\n", coeff_name);
				}
			}
			else
			{
				fprintf(SYMPHAS_WARN, "coefficient spec '%s' is not in a"
					" valid format\n", line_buffer);
			}
		}
	}

	// set the coefficients using the file
	void coeff_from_file(FILE* f, double*& coeff, size_t& coeff_len)
	{
		char line[BUFFER_LENGTH];
		while (fgets(line, BUFFER_LENGTH, f) != NULL)
		{
			coeff_read_line(line, coeff, coeff_len);
		}
	}

	// Parse the coefficients list, optionally looking at the given directory
	// if there is a coefficients file provided.
	void parse_coeff_list(const char* str, const char* dir, double*& coeff, size_t& coeff_len)
	{
		char constants[LINE_READ_BUFFER];
		symphas::lib::str_trim(str, constants);

		/*!
		 * If the first character of the constants is #CONFIG_OPTION_PREFIX_C, then
		 * we will read the constants from a file instead. All the constants which are
		 * not mentioned in the file, but which should be set because a coefficient of
		 * a higher index is initialized, will be initialized to #DEFAULT_COEFF_VALUE.
		 */
		if (constants[0] == CONFIG_OPTION_PREFIX_C)
		{
			char fname[BUFFER_LENGTH_R2];

			if (sscanf(constants, STR(CONFIG_OPTION_PREFIX) "%s", fname) != 1)
			{
				fprintf(SYMPHAS_ERR, "a '%c' was given in the coefficients but an "
					"incorrect format provided: '%s'\n", CONFIG_OPTION_PREFIX_C, fname);
				exit(824);
			}

			// add the extension to the constants file name if it is not already
			// appended, even if there is another string after an occurrence of a dot
			char* dot = std::strrchr(constants, '.');
			if (dot == NULL || std::strcmp(dot + 1, COEFF_SPEC_EXTENSION) != 0)
			{
				std::strcat(fname, "." COEFF_SPEC_EXTENSION);
			}

			FILE* f;
			fname[BUFFER_LENGTH_R2 - 1] = '\0';
			if ((f = fopen(fname, "r")) == 0)
			{
				char* fname_with_dir = new char[std::strlen(dir) + std::strlen(fname) + 1];
				sprintf(fname_with_dir, "%s%s", dir, fname);
				if ((f = fopen(fname_with_dir, "r")) == 0)
				{
					fprintf(SYMPHAS_ERR, "constants file named '%s' could not be read\n", fname);
					exit(825);
				}
				delete[] fname_with_dir;
			}

			coeff_from_file(f, coeff, coeff_len);
		}
		else if (std::strlen(str) > 0)
		{
			size_t i = 0;
			char* constant = std::strtok(constants, " ");
			do
			{
				char* endptr = nullptr;
				double converted = strtod(constant, &endptr);
				if (endptr > constant)
				{
					if (i == coeff_len)
					{
						symphas::lib::expand_append_list(converted, coeff, 1, coeff_len);
						coeff_len += 1;
					}
					else
					{
						coeff[i] = converted;
					}
					i += 1;
				}
				else
				{
					coeff_read_line(constant, coeff, coeff_len);
				}
			} while ((constant = std::strtok(NULL, " ")) != 0);
		}
	}


	// Parse the simple coefficients list
	void parse_simple_coeff_list(const char* str, double*& coeff, size_t& coeff_len)
	{
		coeff_len = 0;
		delete[] coeff;
		coeff = nullptr;

		char* strcpy = new char[std::strlen(str) + 1];
		std::strcpy(strcpy, str);
		char* tok = std::strtok(strcpy, " ");

		if (tok)
		{
			do
			{
				char* endptr;
				double parameter = strtod(tok, &endptr);

				if (!*endptr)
				{
					double* new_coeff = new double[++coeff_len];
					for (iter_type i = 0; i < coeff_len - 1; ++i)
					{
						new_coeff[i] = coeff[i];
					}
					new_coeff[coeff_len - 1] = parameter;
					
					delete[] coeff;
					coeff = new_coeff;
				}
				else
				{
					fprintf(SYMPHAS_WARN, "the value '%s' given to the initial conditions "
						"could not be interpreted, skipping...\n", tok);
				}
			} while ((tok = std::strtok(NULL, " ")) != 0);
		}
	}



	//! Delimiter list used for configuration options.
	/*!
	 * Some parameters can be given as a list of options with delimiters combining
	 * individual elements. There are several types of delimiters, and currently
	 * all of them are valid for packaging a full option inside a config parameter.
	 */
	inline std::initializer_list<char>
		OPTION_DELIMITERS = { '\'', '"' };

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
	inline std::initializer_list<std::pair<char, char>>
		OPTION_BRACKETS = { {'{', '}'}, {'(', ')'} };

	//! Returns open and closed brackets as their own lists.
	/*!
	 * Using the bracket data in #OPTION_BRACKETS, all open brackets are aggregated
	 * into one list and all close brackets are aggregated in one list. Both of
	 * these lists are returned in a pair object.
	 */
	std::pair<std::vector<char>, std::vector<char>> bracket_lists()
	{
		std::vector<char> out_open, out_close;
		for (auto [open, close] : OPTION_BRACKETS)
		{
			out_open.push_back(open);
			out_close.push_back(close);
		}
		return { out_open, out_close };
	}


	char option_open_bracket()
	{
		return symphas::internal::OPTION_BRACKETS.begin()->first;
	}

	char option_close_bracket()
	{
		return symphas::internal::OPTION_BRACKETS.begin()->second;
	}
}



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
	const char* options, bool spaces_are_delimiters, const char* extra_delimiters)
{
	std::vector<std::string> out;
	char* option_buffer = new char[std::strlen(options) + 1];
	std::strcpy(option_buffer, options);
	symphas::lib::str_trim(option_buffer);


	std::vector<char> delimiters{ symphas::internal::OPTION_DELIMITERS };
	char add_delimiter;
	while ((add_delimiter = *extra_delimiters++))
	{
		delimiters.push_back(add_delimiter);
	}

	auto open_list = symphas::internal::bracket_lists().first;
	auto close_list = symphas::internal::bracket_lists().second;
	for (const char* c = option_buffer; *c;)
	{
		/*
		 * copy the names into the preallocated buffers, from
		 * input between double or single quotes, or spaces
		 */

		char buffer[BUFFER_LENGTH];
		bool add_option = false;

		for (char chk : delimiters)
		{
			if (*c == chk)
			{
				char* b = buffer;
				while (*(++c) && *c != chk)
				{
					*(b++) = *c;
				}
				*b = '\0';
				add_option = true;
			}
		}

		if (!add_option)
		{
			/* to be general, this determines the last closing bracket, even if
			 * there is more than one nested bracket type
			 */
			for (auto [open, close] : symphas::internal::OPTION_BRACKETS)
			{
				if (*c == open)
				{
					std::stack<std::ptrdiff_t> bracket_nests;

					// continue gobbling until the final bracket is found
					char* b = buffer;
					while (!(*(++c) == close && bracket_nests.empty()))
					{
						/* add a bracket index if an opening bracket is found
							*/
						auto it_open = std::find(open_list.begin(), open_list.end(), *c);
						if (it_open != open_list.end())
						{
							bracket_nests.push(std::distance(open_list.begin(), it_open));
						}

						/* remove a bracket index if a closing bracket is found
							* if the found closing bracket doesn't match the previous
							* opening bracket, give an error
							*/
						auto it_close = std::find(close_list.begin(), close_list.end(), *c);
						if (it_close != close_list.end())
						{
							if (std::distance(open_list.begin(), it_close) == bracket_nests.top())
							{
								bracket_nests.pop();
							}
							else
							{
								size_t pos = static_cast<size_t>(c - option_buffer);
								fprintf(SYMPHAS_WARN, "there is an unmatched bracket in the options starting at "
									"character %zd in options '%s'\n", pos, option_buffer);
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

		if (!add_option)
		{
			/* if each delimiter has been checked but nothing has been copied,
			 * then copy in everything until a delimiter has been found
			 */
			if (*c != ' ')
			{
				char* b = buffer;

				if (spaces_are_delimiters)
				{
					while (
						*c &&
						*c != ' ')
					{
						*(b++) = *(c++);
					}

				}
				else
				{
					auto a0 = open_list.begin();
					auto a1 = open_list.end();
					auto b0 = delimiters.begin();
					auto b1 = delimiters.end();

					while (*c &&
						(std::find(a0, a1, *c) == a1) &&
						(std::find(b0, b1, *c) == b1))
					{
						*(b++) = *(c++);
					}
				}

				*b = '\0';
				symphas::lib::str_trim(buffer);
				add_option = true;
			}
			else
			{
				for (; *c == ' '; ++c); // skip all the next spaces
			}
		}

		/* copy the buffer into an option string if an option
		 * was copied
		 */
		if (add_option)
		{
			out.emplace_back(buffer);
		}
		else
		{
			//++c;
		}
	}

	delete[] option_buffer;
	return out;
}





SystemConf::SystemConf() :
	dimension{ 1 },
	dt_list{},
	stp{ StencilParams{} },
	runs{ 1 },
	save{ SaveParams{} },
	g{ Geometry::CARTESIAN },
	root_dir{ nullptr },
	result_dir{ nullptr },
	title{ nullptr },
	model{ nullptr },
	names{ nullptr },
	coeff{ 0 },
	intervals{ nullptr },
	bdata{ nullptr },
	tdata{ nullptr },
	num_fields{ nullptr },
	dims{ nullptr },
	intervals_len{ 0 },
	bdata_len{ 0 },
	tdata_len{ 0 },
	names_len{ 0 },
	coeff_len{ 0 },
	num_fields_len{ 0 },
	modifiers{ nullptr },
	init_coeff_copied{ false }
{}

SystemConf::SystemConf(SystemConf const& other) : SystemConf()
{
	g = other.g;
	dimension = other.dimension;
	stp = other.stp;
	runs = other.runs;
	save = other.save;

	dt_list = other.dt_list;

	if (other.title != nullptr)
	{
		title = new char[std::strlen(other.title) + 1];
		std::strcpy(title, other.title);
	}
	else
	{
		title = nullptr;
	}

	if (other.model != nullptr)
	{
		model = new char[std::strlen(other.model) + 1];
		std::strcpy(model, other.model);
	}
	else
	{
		model = nullptr;
	}

	if (other.root_dir != nullptr)
	{
		root_dir = new char[std::strlen(other.root_dir) + 1];
		std::strcpy(root_dir, other.root_dir);
	}
	else
	{
		root_dir = nullptr;
	}

	if (other.result_dir != nullptr)
	{
		result_dir = new char[std::strlen(other.result_dir) + 1];
		std::strcpy(result_dir, other.result_dir);
	}
	else
	{
		result_dir = nullptr;
	}

	intervals_len = other.intervals_len;
	bdata_len = other.bdata_len;
	tdata_len = other.tdata_len;
	names_len = other.names_len;
	coeff_len = other.coeff_len;
	num_fields_len = other.num_fields_len;

	intervals = new symphas::interval_data_type[other.intervals_len]{};
	bdata = new symphas::b_data_type[other.bdata_len]{};
	tdata = new symphas::init_data_type[other.tdata_len]{};
	coeff = new double[other.coeff_len]{};
	names = new char* [other.names_len] {};
	num_fields = new len_type[other.num_fields_len]{};
	dims = new len_type*[other.intervals_len]{};
	
	if (other.num_fields_len > 0)
	{
		modifiers = new char* [other.num_fields_len] {};
		for (iter_type i = 0; i < other.num_fields_len; ++i)
		{
			if (other.modifiers[i])
			{
				modifiers[i] = new char[std::strlen(other.modifiers[i]) + 1];
				std::strcpy(modifiers[i], other.modifiers[i]);
			}
			num_fields[i] = other.num_fields[i];
		}
	}

	std::copy(other.intervals, other.intervals + other.intervals_len, intervals);
	std::copy(other.bdata, other.bdata + other.bdata_len, bdata);
	std::copy(other.tdata, other.tdata + other.tdata_len, tdata);
	std::copy(other.coeff, other.coeff + other.coeff_len, coeff);
	std::copy(other.dims, other.dims + other.intervals_len, dims);

	for (iter_type i = 0; i < other.intervals_len; ++i)
	{
		dims[i] = new len_type[other.dimension];
		std::copy(other.dims[i], other.dims[i] + other.dimension, dims[i]);
	}

	for (iter_type i = 0; i < other.names_len; ++i)
	{
		names[i] = new char[std::strlen(other.names[i]) + 1];
		std::strcpy(names[i], other.names[i]);
	}

	init_coeff_copied = other.init_coeff_copied;

}

SystemConf::SystemConf(SystemConf&& other) noexcept : SystemConf()
{
	swap(*this, other);
}

SystemConf& SystemConf::operator=(SystemConf other)
{
	swap(*this, other);
	return *this;
}

SystemConf::SystemConf(symphas::problem_parameters_type const& parameters, const char* title, const char* dir) : SystemConf()
{
	this->title = new char[std::strlen(title) + 1];
	std::strcpy(this->title, title);

	set_directory(".");

	g = Geometry::CARTESIAN;
	set_time_steps(parameters.get_time_step_list());
	dimension = parameters.get_dimension();

	tdata_len = parameters.length();
	bdata_len = parameters.length();
	intervals_len = parameters.length();

	tdata = new symphas::init_data_type[tdata_len]{};
	bdata = new symphas::b_data_type[bdata_len]{};
	intervals = new symphas::interval_data_type[intervals_len]{};
	dims = new len_type*[intervals_len]{};

	for (iter_type i = 0; i < parameters.length(); ++i)
	{
		set_initial_condition(parameters.get_initial_data()[i], i);

		for (iter_type d = 0; d < dimension * 2; ++d)
		{
			Side side = symphas::index_to_side(d);
			set_boundary(parameters.get_boundary_data()[i].at(side), side, i);
		}

		for (iter_type d = 0; d < dimension; ++d)
		{
			Axis ax = symphas::index_to_axis(d);
			set_interval(parameters.get_interval_data()[i].at(ax), ax, i);
		}
	}
}

symphas::lib::string handle_substitutions(const char* value)
{
	const char* pos_sub = std::strchr(value, CONFIG_SUBSTITUTION_PREFIX_C);
	if (pos_sub == NULL)
	{
		symphas::lib::string result(value, std::strlen(value) + 1);
		return result;
	}
	else
	{
		const char* start = value;
		const char* end = NULL;
		char* buffer = new char[1] {};
		while (pos_sub != NULL)
		{
			int index = -1;
			const char* it;
			for (auto [open, close] : symphas::internal::OPTION_BRACKETS)
			{
				it = pos_sub + 1;
				if (*it == open)
				{
					++it;
					while (*it >= '0' && *it <= '9')
					{
						++it;
					}
					if (*it == close)
					{
						index = atoi(pos_sub + 2);
						end = pos_sub;
						break;
					}
				}
			}

			if (index > 0)
			{
				const char* default_substitution = "";
				const char* substitution;
				if (index < params::config_key_values.size())
				{
					auto& config_entry = params::config_key_values[index - 1];
					substitution = config_entry.data;
				}
				else
				{
					substitution = default_substitution;
				}

				len_type size = end - start;
				len_type substitution_len = std::strlen(substitution);
				len_type buffer_len = std::strlen(buffer);

				char* append_buffer = new char[buffer_len + size + substitution_len + 1] {};
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

SystemConf::SystemConf(std::vector<std::pair<std::string, std::string>> params, const char* title, const char* dir) : SystemConf()
{
	this->title = new char[std::strlen(title) + 1];
	std::strcpy(this->title, title);


	/* raw boundary information which will be placed in a map later
	 */
	auto b_spec = new char* [6];
	auto r_spec = new char* [3];
	char in_spec[BUFFER_LENGTH]{};
	char save_spec[BUFFER_LENGTH]{};
	char model_spec[BUFFER_LENGTH_L3]{};
	char original_dir[BUFFER_LENGTH_L4]{};
	char width_spec[BUFFER_LENGTH]{};
	char dt_spec[BUFFER_LENGTH]{};
	int stop_index = 0;

	for (iter_type i = 0; i < 6; ++i)
	{
		b_spec[i] = new char[BUFFER_LENGTH];
	}

	for (iter_type i = 0; i < 3; ++i)
	{
		r_spec[i] = new char[BUFFER_LENGTH];
	}



	auto copy_range_f = [&](const char* v, Axis ax) 
	{
		std::strncpy(
			r_spec[symphas::axis_to_index(ax)], v,
			BUFFER_LENGTH - 1);
	};

	auto copy_boundary_f = [&](const char* v, Side side)
	{
		std::strncpy(
			b_spec[symphas::side_to_index(side)], v,
			BUFFER_LENGTH - 1);
	};

	auto set_dimension = [](const char* v)
	{
		len_type available_dimensions[]{ AVAILABLE_DIMENSIONS };
		if (v != NULL)
		{
			if (v[0] == CONFIG_OPTION_PREFIX_C)
			{
				if (sizeof(available_dimensions) / sizeof(len_type) > 0)
				{
					return available_dimensions[0];
				}
				else
				{
					return 0;
				}
			}
			else
			{
				return atoi(v);
			}
		}
		else
		{
			if (sizeof(available_dimensions) / sizeof(len_type) > 0)
			{
				return available_dimensions[0];
			}
			else
			{
				return 0;
			}
		}
	};


	std::map<std::string, std::function<void(const char*)>, symphas::internal::any_case_comparator> system_map =
	{
		{ symphas::internal::C_DIM,			[&](const char* v) { dimension = set_dimension(v); } },
		{ symphas::internal::C_SOLVERVAR,	[&](const char* v) { stp.type = atoi(v); } },
		{ symphas::internal::C_ORDER,		[&](const char* v) { select_stencil_accuracy(dimension, v); }},
		{ symphas::internal::C_PTL,			[&](const char* v) { select_stencil(2, v); } },
		{ symphas::internal::C_PTB,			[&](const char* v) { select_stencil(4, v); } },
		{ symphas::internal::C_PTG,			[&](const char* v) { select_stencil(3, v); } },
		{ symphas::internal::C_WIDTH,		[&](const char* v) { std::strncpy(width_spec, v, sizeof(width_spec) / sizeof(char) - 1); } },
		{ symphas::internal::C_FORM,		[&](const char* v) { parse_form(v); } },
		{ symphas::internal::C_RNGX,		[&](const char* v) { copy_range_f(v, Axis::X); } },
		{ symphas::internal::C_RNGY,		[&](const char* v) { copy_range_f(v, Axis::Y); } },
		{ symphas::internal::C_RNGZ,		[&](const char* v) { copy_range_f(v, Axis::Z); } },
		{ symphas::internal::C_RNGR,		[&](const char* v) { copy_range_f(v, Axis::X); } },
		{ symphas::internal::C_RNGT,		[&](const char* v) { copy_range_f(v, Axis::Y); } },
		{ symphas::internal::C_RNGF,		[&](const char* v) { copy_range_f(v, Axis::Z); } },
		{ symphas::internal::C_BNDLT,		[&](const char* v) { copy_boundary_f(v, Side::LEFT); } },
		{ symphas::internal::C_BNDRT,		[&](const char* v) { copy_boundary_f(v, Side::RIGHT); } },
		{ symphas::internal::C_BNDTP,		[&](const char* v) { copy_boundary_f(v, Side::TOP); } },
		{ symphas::internal::C_BNDBT,		[&](const char* v) { copy_boundary_f(v, Side::BOTTOM); } },
		{ symphas::internal::C_BNDFT,		[&](const char* v) { copy_boundary_f(v, Side::FRONT); } },
		{ symphas::internal::C_BNDBK,		[&](const char* v) { copy_boundary_f(v, Side::BACK); } },
		{ symphas::internal::C_INSIDE,		[&](const char* v) { std::strncpy(in_spec, v, sizeof(in_spec) / sizeof(char) - 1); } },
		{ symphas::internal::C_DELTA,		[&](const char* v) { std::strncpy(dt_spec, v, sizeof(dt_spec) / sizeof(char) - 1); } },
		{ symphas::internal::C_MODEL,		[&](const char* v) { std::strncpy(model_spec, v, sizeof(model_spec) / sizeof(char) - 1); } },
		{ symphas::internal::C_NAMES,		[&](const char* v) { parse_names(v); } },
		{ symphas::internal::C_FRAMES,		[&](const char* v) { stop_index = atoi(v); } },
		{ symphas::internal::C_SAVE,		[&](const char* v) { std::strncpy(save_spec, v, sizeof(save_spec) / sizeof(char) - 1); } },
		{ symphas::internal::C_SAVEINIT,	[&](const char* v) { save.set_init_flag(std::strcmp(v, "NO") != 0); } },
		{ symphas::internal::C_RUNS,		[&](const char* v) { runs = atoi(v); } },
		{ symphas::internal::C_DIR,			[&](const char* v) { std::strncpy(original_dir, v, sizeof(original_dir) / sizeof(char) - 1); } }
	};

	for (auto& [k, v] : params)
	{
		if (system_map.count(k) != 0)
		{
			system_map[k](v.c_str());
		}
	}

	save.set_stop(stop_index);
	symphas::io::parse_save_str(save_spec, &save);

	if (stp.ord == 0)
	{
		fprintf(SYMPHAS_ERR, "an order of accuracy greater than 0 must be selected\n");
		exit(8010);
	}


	if (stp.ptl == StencilParams{}.ptl)
	{
		select_stencil(2, STR(CONFIG_OPTION_PREFIX));
	}
	if (stp.ptg == StencilParams{}.ptg)
	{
		select_stencil(3, STR(CONFIG_OPTION_PREFIX));
	}
	if (stp.ptb == StencilParams{}.ptb)
	{
		select_stencil(4, STR(CONFIG_OPTION_PREFIX));
	}

	auto dts = dt_list.get_time_steps();

	// check integrity of parameters
	if (stp.ptl < 1 || stp.ptb < 1 || stp.ptg < 1)
	{
		fprintf(SYMPHAS_ERR, "a stencil parameter is out of range\n");
		exit(8011);
	}
	else if (runs < 1)
	{
		fprintf(SYMPHAS_ERR, "the number of runs must be greater than 0\n");
		exit(8012);
	}
	else if (*std::min_element(dts.begin(), dts.end()) <= 0)
	{
		fprintf(SYMPHAS_ERR, "the time step must be greater than 0\n");
		exit(8013);
	}
	else if (stop_index < INDEX_INIT || save.get_start() < INDEX_INIT)
	{
		fprintf(SYMPHAS_ERR, "a save parameter is out of range\n");
		exit(8014);
	}
	else if (dimension > 3)
	{
		fprintf(SYMPHAS_ERR, "dimension must be either 1, 2 or 3\n");
		exit(8015);
	}



	// set the final configuration parameters
	parse_model_spec(handle_substitutions(model_spec), dir);
	parse_width(width_spec);
	parse_interval_array(r_spec);
	parse_boundaries_array(b_spec);
	parse_initial_condition_array(in_spec);
	parse_dt(dt_spec);
	set_directory(original_dir);

#ifdef FILESYSTEM_HEADER_AVAILABLE
	std::filesystem::path p_dir(".");
	p_dir.append(result_dir);

	symphas::lib::make_directory(p_dir / CHECKPOINT_DIR, 882);
	symphas::lib::make_directory(p_dir / DATA_DIR, 883);
	symphas::lib::make_directory(p_dir / PLOT_DIR, 884);

#else

	size_t p_dir_len = std::strlen(result_dir) + 3;
	char* p_dir = new char[p_dir_len];
	snprintf(p_dir, p_dir_len, "./%s", result_dir);


	size_t sub_dir_len = p_dir_len + BUFFER_LENGTH;
	char* sub_dir = new char[sub_dir_len];

	snprintf(sub_dir, p_dir_len + sizeof(CHECKPOINT_DIR) / sizeof(char) + 1, "%s/%s", p_dir, CHECKPOINT_DIR);
	symphas::lib::make_directory(sub_dir, 882);
	snprintf(sub_dir, p_dir_len + sizeof(DATA_DIR) / sizeof(char) + 1, "%s/%s", p_dir, DATA_DIR);
	symphas::lib::make_directory(sub_dir, 882);
	snprintf(sub_dir, p_dir_len + sizeof(PLOT_DIR) / sizeof(char) + 1, "%s/%s", p_dir, PLOT_DIR);
	symphas::lib::make_directory(sub_dir, 882);


	delete[] p_dir;
	delete[] sub_dir;

#endif


	/* create a backup configuration under the parent output directory
	 */


	for (iter_type i = 0; i < 6; ++i)
	{
		delete[] b_spec[i];
	}

	for (iter_type i = 0; i < 3; ++i)
	{
		delete[] r_spec[i];
	}

	delete[] b_spec;
	delete[] r_spec;
}

void SystemConf::set_dimensions(size_t n)
{
	const char msg[] = "the number of grid points cannot be less than 1\n";
	delete[] dims[n];

	switch (dimension)
	{
	case 1:
	{
		iter_type L = intervals[n].at(Axis::X).get_domain_count();

		if (L < 1)
		{
			fprintf(SYMPHAS_LOG, msg);
			exit(802);
		}

		dims[n] = new len_type[1]{ L };

		break;
	}
	case 2:
	{
		iter_type L = intervals[n].at(Axis::X).get_domain_count();
		iter_type M = intervals[n].at(Axis::Y).get_domain_count();

		if (L < 1 || M < 1)
		{
			fprintf(SYMPHAS_LOG, msg);
			exit(803);
		}

		dims[n] = new len_type[2]{ L, M };

		break;
	}
	case 3:
	{
		iter_type L = intervals[n].at(Axis::X).get_domain_count();
		iter_type M = intervals[n].at(Axis::Y).get_domain_count();
		iter_type N = intervals[n].at(Axis::Z).get_domain_count();

		if (L < 1 || M < 1 || N < 1)
		{
			fprintf(SYMPHAS_LOG, msg);
			exit(804);
		}

		dims[n] = new len_type[3]{ L, M, N };

		break;
	}
	default:
		fprintf(SYMPHAS_ERR, "incompatible dimensions\n");
		exit(805);
	}

}

void SystemConf::parse_boundaries_array(const char* const* b_specs)
{
	std::vector<std::string> *boundary_lists = new std::vector<std::string>[dimension * 2];
	size_t* boundary_lens = new size_t[dimension * 2];

	// compile all the boundaries given to the configuration
	for (iter_type i = 0; i < dimension * 2; ++i)
	{
		boundary_lists[i] = symphas::conf::parse_options(b_specs[i]);
		boundary_lens[i] = boundary_lists[i].size();
	}

	// set the number of boundary options given and top off any boundaries
	// which haven't been fully specified
	size_t len = *std::max_element(boundary_lens, boundary_lens + dimension * 2);
	if (len > bdata_len)
	{
		delete[] bdata;

		bdata_len = len;
		bdata = new symphas::b_data_type[bdata_len];
	}

	for (iter_type i = 0; i < dimension * 2; ++i)
	{
		// if the boundary specification for a given side is not "full", i.e. not
		// of all systems have boundaries specified, then append copies of the first
		// boundary to the rest.
		if (boundary_lens[i] < len)
		{
			std::vector<std::string> append{ len - boundary_lens[i], boundary_lists[i][0] };
			boundary_lists[i].insert(
				boundary_lists[i].end(), 
				std::make_move_iterator(append.begin()), std::make_move_iterator(append.end()));
		}
	}
	

	// set the boundaries for each of the options by copying the spec string
	// into the temporary array
	const char** b_spec = new const char* [dimension * 2];
	for (iter_type n = 0; n < len; ++n)
	{
		for (iter_type i = 0; i < dimension * 2; ++i)
		{
			b_spec[i] = boundary_lists[i][n].c_str();
		}

		parse_boundaries(b_spec, n);
	}



	delete[] boundary_lists;
	delete[] boundary_lens;
	delete[] b_spec;

}

void SystemConf::parse_boundaries(const char* const* b_spec, size_t n)
{
	/* set the boundary information
	 *
	 * in order to set the correct intervals, the index n is maintained
	 * moreover, the boundaries are filled in reverse order, since the intervals
	 * correspond with the boundaries this way
	 */

	switch (dimension)
	{
	case 1:
		for (auto side : { Side::LEFT, Side::RIGHT })
		{
			auto [it, flag] = bdata[n].emplace(side, symphas::b_element_type{});
			if (flag)
			{
				symphas::internal::fill_boundary(it->second, b_spec[symphas::side_to_index(side)]);
			}
		}
		break;

	case 2:
		for (auto side : { Side::LEFT, Side::RIGHT, Side::TOP, Side::BOTTOM })
		{
			auto [it, flag] = bdata[n].emplace(side, symphas::b_element_type{});
			if (flag)
			{
				symphas::internal::fill_boundary(it->second, b_spec[symphas::side_to_index(side)]);
			}
		}
		break;

	case 3:
		for (auto side : { Side::LEFT, Side::RIGHT, Side::TOP, Side::BOTTOM, Side::FRONT, Side::BACK })
		{
			auto [it, flag] = bdata[n].emplace(side, symphas::b_element_type{});
			if (flag)
			{
				symphas::internal::fill_boundary(it->second, b_spec[symphas::side_to_index(side)]);
			}
		}
		break;

	default:
		break;
	}


}


void SystemConf::parse_form(const char* form)
{
	if (std::strcmp(form, "CARTESIAN") == 0)
	{
		g = Geometry::CARTESIAN;
	}
	else if (std::strcmp(form, "POLAR") == 0)
	{
		g = Geometry::POLAR;
	}
	else
	{
		g = Geometry::CARTESIAN;
	}
}

void SystemConf::parse_initial_condition_array(const char* value)
{
	auto inits = symphas::conf::parse_options(value);

	if (inits.size() > tdata_len)
	{
		delete[] tdata;

		tdata_len = inits.size();
		tdata = new symphas::init_data_type[tdata_len];
	}

	for (iter_type n = 0; n < inits.size(); ++n)
	{
		parse_initial_condition(inits[n].c_str(), n);
	}
}

symphas::init_entry_type get_initial_condition_entry(char* input, size_t dimension, double* coeff, size_t coeff_len, bool& init_coeff_copied)
{
	symphas::init_entry_type init;

	//char* pos0 = input;
	char* tok = input;
	while (*++tok && *tok != ' ');

	char buffer[BUFFER_LENGTH]{ 0 };
	char name[BUFFER_LENGTH]{ 0 };
	std::copy(input, tok, name);
	std::copy(tok, input + std::strlen(input) + 1, buffer);

	init.in = symphas::in_from_str(name);

	if (init.in == Inside::EXPRESSION)
	{
		//char* pos0 = input;
		tok = std::strtok(buffer, " ");
		char* expression_name = new char[std::strlen(tok) + 1];

		if (sscanf(tok, "%s", expression_name) < 1)
		{
			fprintf(SYMPHAS_ERR, "unable to read the expression name given to the initial conditions, '%s'\n", input);
			exit(841);
		}
		else
		{
			init.expr_data = expression_name;
		}

		delete[] expression_name;
		tok = std::strtok(NULL, " ");
	}
	else
	{
		switch (init.in)
		{
		case Inside::GAUSSIAN:
		{
			init.data.gp[0] = 0;		// mean
			init.data.gp[1] = 0.25;		// std deviation
			break;
		}
		case Inside::UNIFORM:
		{
			init.data.gp[0] = -1;		// minimum value
			init.data.gp[1] = 1;		// maximum value
			break;
		}
		case Inside::CAPPED:
		{
			init.data.gp[0] = -1;		// minimum value
			init.data.gp[1] = 1;		// maximum value
			break;
		}
		case Inside::CONSTANT:
		{
			init.data.gp[0] = 0;		// the constant value
			break;
		}
		case Inside::CIRCLE:
		{
			init.data.gp[0] = 2;		// ratio of the x width
			init.data.gp[1] = 2;		// ratio of the y width
			if (dimension == 3)
			{
				init.data.gp[2] = 2;	// ratio of z width
			}
			break;
		}
		case Inside::HEXAGONAL:
		{
			init.data.gp[0] = 1;		// ratio of density along x
			init.data.gp[1] = 1;		// ratio of density along y
			if (dimension == 2)
			{
				init.data.gp[2] = 2;	// ratio of size of nodes
			}
			else if (dimension == 3)
			{
				init.data.gp[2] = 1;	// ratio of density along z
				init.data.gp[3] = 2;	// ratio of density of nodes
			}
			break;
		}
		case Inside::CUBIC:
		{
			init.data.gp[0] = 1;		// same as HEX
			init.data.gp[1] = 1;
			if (dimension == 2)
			{
				init.data.gp[2] = 2;
			}
			else if (dimension == 3)
			{
				init.data.gp[2] = 1;
				init.data.gp[3] = 2;
			}
			break;
		}
		case Inside::SEEDSSQUARE:
		{
			init.data.gp[0] = 10;		// the number of seeds in the system
			init.data.gp[1] = 1;		// the scale factor of the size
			// the value inside the seed
			init.data.gp[2] = params::init_inside_val;
			// the value outside the seed
			init.data.gp[3] = params::init_outside_val;
			break;
		}
		case Inside::SEEDSCIRCLE:
		{
			init.data.gp[0] = 10;
			init.data.gp[1] = 1;
			init.data.gp[2] = params::init_inside_val;
			init.data.gp[3] = params::init_outside_val;
			break;
		}
		case Inside::VORONOI:
		{
			init.data.gp[0] = 10;		// number of crystals
			init.data.gp[1] = -1;		// lower range
			init.data.gp[2] = 1;		// upper range
			break;
		}
		case Inside::BUBBLE:
		{
			init.data.gp[0] = 1;		// The number of bubbles to fill
			init.data.gp[1] = -1;		// lower range
			init.data.gp[2] = 1;		// upper range
			init.data.gp[3] = .75;		// The filling ratio
			init.data.gp[4] = 1;		// THe overlap ratio
			break;
		}
		case Inside::SPIRALHEX:
		{
			init.data.gp[0] = 1;		// The number of hexes to fill
			init.data.gp[1] = -1;		// lower range
			init.data.gp[2] = 1;		// upper range
			init.data.gp[3] = 1;		// size of the circle
			break;
		}
		default:
			init.data = symphas::init_data_parameters::one();
			break;
		}

		input = tok;
		while (true)
		{
			if (*tok == ' ') while (*++tok == ' ');
			while (*tok != ' ' && *++tok);

			std::copy(input, tok, buffer);
			buffer[tok - input] = '\0';

			symphas::lib::str_trim(buffer);
			InsideTag tag = symphas::in_tag_from_str(buffer);

			if (tag != InsideTag::NONE)
			{
				init.intag = symphas::build_intag(init.intag, tag);
				input = tok;
			}
			else
			{
				break;
			}
		}
	}

	size_t gp_count = 0;
	if (tok)
	{
		if (*tok == CONFIG_OPTION_PREFIX_C)
		{
			if (init.in == Inside::EXPRESSION)
			{
				init.expr_data.set_coeff(coeff, coeff_len);
			}
			else
			{
				std::copy(coeff, coeff + std::min(coeff_len, size_t(NUM_INIT_CONSTANTS)), init.data.gp);

				if (coeff_len > NUM_INIT_CONSTANTS)
				{
					fprintf(SYMPHAS_WARN, "using only %d numeric arguments from model coefficients\n",
						NUM_INIT_CONSTANTS);
				}
			}
			init_coeff_copied = true;
		}
		else
		{
			double* init_coeff = nullptr;
			symphas::internal::parse_simple_coeff_list(input, init_coeff, gp_count);

			if (init.in == Inside::EXPRESSION)
			{
				init.expr_data.set_coeff(init_coeff, gp_count);
			}
			else
			{
				std::copy(init_coeff, init_coeff + std::min(gp_count, size_t(NUM_INIT_CONSTANTS)), init.data.gp);
				if (gp_count > NUM_INIT_CONSTANTS)
				{
					fprintf(SYMPHAS_WARN, "only %d numeric arguments may be provided "
						"to the initial conditions\n", NUM_INIT_CONSTANTS);
				}
			}

			delete[] init_coeff;
		}
	}

	auto message_unset = [&] (size_t m)
	{
		if (gp_count < m)
		{
			fprintf(SYMPHAS_LOG, "using fewer than possible number of parameters for the initial "
				"condition '%s', the remaining %zd parameters are given the default values: \n",
				name, m - gp_count);
			for (size_t i = gp_count; i < m; fprintf(SYMPHAS_LOG, "\t%.2f \n", init.data.gp[i++]));
		}
		else if (gp_count > m)
		{
			fprintf(SYMPHAS_LOG, "more initial condition parameters than required, the initial "
				"condition '%s' uses up to %zd parameters\n", name, m);
		}
		else {}
	};


	/* if there are fewer parameters than what is given,
	 */

	switch (init.in)
	{
	case Inside::GAUSSIAN:
	{
		message_unset(2);
		break;
	}
	case Inside::UNIFORM:
	{
		message_unset(2);
		break;
	}
	case Inside::CAPPED:
	{
		message_unset(2);
		break;
	}
	case Inside::CONSTANT:
	{
		message_unset(1);
		break;
	}
	case Inside::CIRCLE:
	{
		if (dimension == 2)
		{
			message_unset(2);
		}
		else if (dimension == 3)
		{
			message_unset(3);
		}
		break;
	}
	case Inside::HEXAGONAL:
	{
		if (dimension == 2)
		{
			message_unset(3);
		}
		else if (dimension == 3)
		{
			message_unset(4);
		}
		break;
	}
	case Inside::CUBIC:
	{
		if (dimension == 2)
		{
			message_unset(3);
		}
		else if (dimension == 3)
		{
			message_unset(4);
		}
		break;
	}
	case Inside::SEEDSSQUARE:
	{
		message_unset(4);
		break;
	}
	case Inside::SEEDSCIRCLE:
	{
		message_unset(4);
		break;
	}
	case Inside::VORONOI:
	{
		message_unset(3);
		break;
	}
	case Inside::BUBBLE:
	{
		message_unset(5);
		break;
	}
	case Inside::SPIRALHEX:
	{
		message_unset(4);
		break;
	}
	default:
		break;
	}


	return init;
}

void SystemConf::parse_initial_condition(const char* value, size_t n)
{
	char input[BUFFER_LENGTH];
	std::strncpy(input, value, sizeof(input) / sizeof(char) - 1);
	symphas::lib::str_trim(input);

	if (*input == CONFIG_TITLE_PREFIX_C && std::strlen(input) == 1)
	{
		if (params::input_data_file)
		{
			auto inits = symphas::conf::parse_options(params::input_data_file, ",");

			tdata[n][Axis::NONE].file = { inits[n % inits.size()].c_str(), 0};
			tdata[n][Axis::NONE].in = Inside::FILE;
		}
		else
		{
			tdata[n][Axis::NONE] = Inside::CONSTANT;
		}
	}
	else if (*input == CONFIG_OPTION_PREFIX_C)
	{
		char* tok = std::strtok(input, STR(CONFIG_OPTION_PREFIX));
		if (tok)
		{
			char* file_name = new char[std::strlen(tok) + 1] { 0 };
			int index = 0;

			if (sscanf(tok, "%s %d", file_name, &index) < 2)
			{
				fprintf(SYMPHAS_WARN, "using index %d from intial conditions input file\n", index);
				if (sscanf(tok, "%s", file_name) < 1)
				{
					fprintf(SYMPHAS_WARN, "unable to read the filename given to the initial conditions, '%s'\n", value);
				}
			}

			if (file_name)
			{
				tdata[n][Axis::NONE].file = { file_name, index };
				tdata[n][Axis::NONE].in = Inside::FILE;
			}

			delete[] file_name;
		}
	}
	else
	{
		Axis ax = Axis::NONE;

		// go to the next exclamation mark
		char* input0 = input;
		char* last = input0;

		do
		{
			while (*++last && *last != CONFIG_TITLE_PREFIX_C);

			if (*input0 == CONFIG_TITLE_PREFIX_C)
			{

				char ax_str[6]{ 0 };
				size_t nn = sscanf(input0 + 1, "%5s", ax_str);
				input0 += nn + 1;

				if (std::strcmp(ax_str, "X") == 0)
				{
					ax = Axis::X;
				}
				else if (std::strcmp(ax_str, "Y") == 0)
				{
					ax = Axis::Y;
				}
				else if (std::strcmp(ax_str, "Z") == 0)
				{
					ax = Axis::Z;
				}
				else if (std::strcmp(ax_str, "R") == 0)
				{
					ax = Axis::R;
				}
				else if (std::strcmp(ax_str, "T") == 0)
				{
					ax = Axis::T;
				}
				else if (std::strcmp(ax_str, "S") == 0)
				{
					ax = Axis::S;
				}
				else
				{
					fprintf(SYMPHAS_ERR, "incorrect axis specification for initial condition, '%s'", ax_str);
				}

			}
			char c = *last;
			*last = '\0';

			symphas::lib::str_trim(input0);
			tdata[n][ax] = get_initial_condition_entry(input0, dimension, coeff, coeff_len, init_coeff_copied);

			*last = c;
			input0 = last;
			
		} while (*input0);

		if (tdata[n].find(Axis::NONE) == tdata[n].end())
		{
			tdata[n][Axis::NONE] = tdata[n].begin()->second;
		}
	}

}

void SystemConf::parse_interval_array(const char* const* r_specs)
{
	for (iter_type d = 0; d < dimension; ++d)
	{
		Axis ax = symphas::index_to_axis(d);
		const char* r_spec = r_specs[d];

		auto ranges = symphas::conf::parse_options(r_spec);
		if (ranges.size() > intervals_len)
		{
			extend_intervals(ranges.size());
		}

		for (iter_type i = 0; i < ranges.size(); ++i)
		{
			parse_interval(ranges[i].c_str(), ax, i);
		}
	}

	for (iter_type i = 0; i < intervals_len; ++i)
	{
		set_dimensions(i);
	}

}

void SystemConf::parse_interval(const char* value, Axis ax, size_t n)
{
	if (value[0] == CONFIG_OPTION_PREFIX_C)
	{
		iter_type count = 1;
		iter_type stop = 0;

		if (sscanf(value, STR(CONFIG_OPTION_PREFIX) " %d%n", &count, &stop) != 1)
		{
			fprintf(SYMPHAS_WARN, "the grid points given by '%s' is not in the correct "
				"format, expecting the number of grid points as `@ L`\n", value);
		}

		intervals[n][ax].set_count_from_r(count, intervals[n][ax].width(), 0.0);

		const char* split = strchr(value, CONFIG_TITLE_PREFIX_C);
		if (split != NULL)
		{
			iter_type left, right;
			if (sscanf(split + 1, "%d %d", &left, &right) != 2)
			{
				fprintf(SYMPHAS_WARN, "the subinterval specification given by '%s' is not in the correct "
					"format, expecting the format to be provided as `@ L ! A B`\n", split + 1);
			}
			intervals[n][ax].set_interval_fraction((double)left / count, (double)right / count);
		}
		else
		{
			intervals[n][ax].interval_to_domain();
		}
	}
	else
	{
		double left = 0, right = 1;
		if (sscanf(value, "%lf %lf", &left, &right) != 2)
		{
			fprintf(SYMPHAS_WARN, "the interval given by '%s' is not in the correct "
				"format, expecting `A B`\n", value);
		}
		intervals[n][ax].set_domain(left, right, intervals[n][ax].width());
		intervals[n][ax].interval_to_domain();
	}
}

void SystemConf::parse_names(const char* value)
{
	auto name_list = symphas::conf::parse_options(value, true);

	if (name_list.size() > 0)
	{
		if (names_len < name_list.size())
		{
			names = new char* [names_len];
		
			for (iter_type i = 0; i < names_len; ++i)
			{
				delete[] names[i];
			}
			delete[] names;
			names_len = name_list.size();
			names = new char* [names_len];
		}

		char** name_it = names;
		for (auto&& name : name_list)
		{
			*name_it = new char[std::strlen(name.c_str()) + 1];
			std::strcpy(*name_it, name.c_str());
			name_it++;
		}
	}
}


void SystemConf::set_name(const char* name, size_t n)
{
	if (n > names_len)
	{
		char** extend = new char* [n];
		for (iter_type i = 0; i < names_len; ++i)
		{
			extend[i] = new char[std::strlen(names[i]) + 1];
			std::strcpy(extend[i], names[i]);
		}
		for (size_t i = names_len; i < n; ++i)
		{
			extend[i] = nullptr;
		}


		delete[] names;
		names = extend;
		names_len = n;
	}


	delete[] names[n];
	names[n] = new char[std::strlen(name) + 1];
	std::strcpy(names[n], name);
}


void SystemConf::set_directory(const char* directory)
{
	delete[] result_dir;
	delete[] root_dir;

	root_dir = new char[std::strlen(directory) + 1];
	std::strcpy(root_dir, directory);

	char title_dir[BUFFER_LENGTH_R2];
	symphas::lib::to_file_name(title, title_dir, BUFFER_LENGTH_R2);

	char append_dir[BUFFER_LENGTH_R1];
	snprintf(append_dir, BUFFER_LENGTH_R1, "%s/%s", root_dir, title_dir);


	// if we're using a timestamp, then the created directory using that timestamp must
	// be unique, add an index to dir if one already exists
	if (params::use_timestamp)
	{
		/* add timestamp to the directory tree
		 */

		char ts_buffer[BUFFER_LENGTH_R2];

#ifdef USING_MPI
		MPI_Request *rqst = new MPI_Request[symphas::parallel::get_num_nodes()]{};

		if (symphas::parallel::is_host_node())
		{
			symphas::lib::write_ts_str(ts_buffer);
			for (iter_type i = 0; i < symphas::parallel::get_num_nodes(); ++i)
			{
				if (!symphas::parallel::is_host_node(i))
				{
					rqst[i] = MPI_REQUEST_NULL;
					MPI_Send(&ts_buffer[0], BUFFER_LENGTH_R2, MPI_CHAR, i, 0, MPI_COMM_WORLD);
				}
			}
		}
		else
		{
			MPI_Status status;
			MPI_Recv(&ts_buffer[0], BUFFER_LENGTH_R2, MPI_CHAR, SYMPHAS_MPI_HOST_RANK, 0, MPI_COMM_WORLD, &status);
		}
#else
		symphas::lib::write_ts_str(ts_buffer);
#endif

		size_t dlen = std::strlen(append_dir) + std::strlen(ts_buffer) + 2;

#ifdef USING_MPI
		if (symphas::parallel::is_host_node())
		{
			result_dir = new char[dlen + 1];
			sprintf(result_dir, "%s/%s", append_dir, ts_buffer);
		}
		else
		{
			int rank = symphas::parallel::get_node_rank();
			result_dir = new char[dlen + 1 + symphas::lib::num_digits(rank)];
			sprintf(result_dir, "%s/%s/%d", append_dir, ts_buffer, rank);
		}
#else
		result_dir = new char[dlen];
		sprintf(result_dir, "%s/%s", append_dir, ts_buffer);
#endif

		/*
		 * decide if index should be appended
		 *
		 * if multiple runs are initiated at the same time, then the program will detect that
		 * the directory with that timestamp name already exists and it will append an index to the end
		 */

		iter_type i = 0;
		struct stat info;
		while (
			stat(result_dir, &info) == 0 
			&& (info.st_mode & S_IFMT) == S_IFDIR 
			&& i < 10)
		{

#ifdef USING_MPI
			char* new_dir;
			if (symphas::parallel::is_host_node())
			{
				char* new_dir = new char[dlen + STR_ARR_LEN(TIMESTAMP_ID_APPEND)];
				sprintf(new_dir, "%s/%s" TIMESTAMP_ID_APPEND  "", append_dir, ts_buffer, ++i);
			}
			else
			{
				int rank = symphas::parallel::get_node_rank();
				char* new_dir = new char[dlen + symphas::lib::num_digits(rank) + STR_ARR_LEN(TIMESTAMP_ID_APPEND)];
				sprintf(new_dir, "%s/%s" TIMESTAMP_ID_APPEND  "/%d", append_dir, ts_buffer, ++i, rank);
			}
#else
			char* new_dir = new char[dlen + STR_ARR_LEN(TIMESTAMP_ID_APPEND)];
			sprintf(new_dir, "%s/%s" TIMESTAMP_ID_APPEND, append_dir, ts_buffer, ++i);
#endif

			delete[] result_dir;
			result_dir = new_dir;
		}
	}
	else
	{
#ifdef USING_MPI
		int rank = symphas::parallel::get_node_rank();
		root_dir = new char[std::strlen(append_dir) + 2 + symphas::lib::num_digits(rank)];
		sprintf(result_dir, "%s/%d", append_dir, rank);
#else
		result_dir = new char[std::strlen(append_dir) + 1];
		std::strcpy(result_dir, append_dir);
#endif
	}
}

void SystemConf::select_stencil(size_t order, const char* str)
{
	const char msg[] = "derivative order '%zd' is invalid in selection of stencil point values\n";
	if (*str == CONFIG_OPTION_PREFIX_C)
	{
		StencilParams default_stp = DefaultStencil{ dimension, stp.ord }();
		if (default_stp != StencilParams{})
		{

			switch (order)
			{
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
		}
		else
		{
			fprintf(SYMPHAS_WARN, "not able to select a default stencil with the given "
				"dimension (%zd) and order of accuracy (%hu)\n", dimension, stp.ord);
		}
	}
	else
	{
		unsigned short value = static_cast<unsigned short>(std::strtoul(str, NULL, 10));
		if (value > 0 && errno != ERANGE)
		{
			switch (order)
			{
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



void SystemConf::select_stencil_accuracy(size_t dimension, const char* str)
{
	if (*str == CONFIG_OPTION_PREFIX_C)
	{
		StencilParams default_stp = DefaultStencil{ dimension, 0 }();
		if (default_stp != StencilParams{})
		{
			stp.ord = default_stp.ord;
		}
		else
		{
			fprintf(SYMPHAS_WARN, "not able to select a default stencil "
				"accuracy with the given dimension (%zd) \n", dimension);
		}
	}
	else
	{
		unsigned short value = static_cast<unsigned short>(std::strtoul(str, NULL, 10));
		stp.ord = value;
	}
}

void SystemConf::parse_model_spec(const char* value, const char* dir)
{
	char type[BUFFER_LENGTH + 1]{};
	int model_end_pos;

	if (std::strlen(value) > 1 && value[0] == CONFIG_TITLE_PREFIX_C)
	{
		iter_type start, end;
		size_t n = sscanf(value, CONFIG_TITLE_PREFIX " %" STR(BUFFER_LENGTH) "s %d %d", type, &start, &end);

		if (n == 0)
		{
			fprintf(SYMPHAS_ERR, "expected file name after '%c' in model specification\n", CONFIG_TITLE_PREFIX_C);
		}
		else
		{

			len_type type_len = std::strlen(type);
			len_type keyword_len = sizeof(STR(VIRTUAL_MODEL_KEYWORD)) / sizeof(char);
			len_type index_len = int(n > 1) * symphas::lib::num_digits(start) + int(n > 2) * symphas::lib::num_digits(end) + (n - 1);
			char* model_spec = new char[keyword_len + type_len + 1 + index_len] {};

			sprintf(model_spec, "%s%c%s", STR(VIRTUAL_MODEL_KEYWORD), VIRTUAL_MODEL_SEP_KEY_C, type);
			if (n > 1)
			{
				sprintf(model_spec + keyword_len + type_len, "%c%d", VIRTUAL_MODEL_SEP_KEY_C, start);
			}
			if (n > 2)
			{
				sprintf(model_spec + keyword_len + type_len + symphas::lib::num_digits(start) + 1, "%c%d", VIRTUAL_MODEL_SEP_KEY_C, end);
			}

			std::swap(model, model_spec);
			delete[] model_spec;

		}
	}
	else
	{

		// set the model and list of coefficients to read them, the
		// default is given by #DEFAULT_COEFF_VALUE.
		size_t n = sscanf(value, " %" STR(BUFFER_LENGTH) "s %n", type, &model_end_pos);

		if (n == 1)
		{
			set_model_name(type);

			if (model_end_pos > std::strlen(model))
			{
				const char* iter = value + model_end_pos - 1;
				while (*++iter && *iter != CONFIG_TITLE_PREFIX_C);

				if (*iter == CONFIG_TITLE_PREFIX_C)
				{
					auto options = symphas::conf::parse_options(iter + 1, ",");
					delete[] num_fields;
					for (iter_type i = 0; i < num_fields_len; ++i)
					{
						delete[] modifiers[i];
					}
					delete[] modifiers;

					num_fields_len = options.size();
					num_fields = new len_type[num_fields_len]{};
					modifiers = new char* [num_fields_len] {};


					iter_type i = 0;
					for (auto option : options)
					{
						*type = '\0';
						sscanf(option.c_str(), "%d %s", num_fields + i, type);
						if (*type)
						{
							modifiers[i] = new char[std::strlen(type) + 1];
							std::strcpy(modifiers[i], type);
						}
						else
						{
							modifiers[i] = nullptr;
						}
						++i;
					}
				}

				std::copy(value + model_end_pos, iter, type);
				type[iter - (value + model_end_pos)] = '\0';
				symphas::internal::parse_coeff_list(type, dir, coeff, coeff_len);
			}
			else
			{
				if (coeff)
				{
					coeff_len = 0;
				}
			}
		}
		else
		{
			fprintf(SYMPHAS_ERR, "model name was not specified at the configuration "
				"specification for %s\n", symphas::internal::C_MODEL);
			model = nullptr;
		}
	}
}


void SystemConf::set_model_name(const char* str)
{
	delete[] model;
	model = new char[std::strlen(str) + 1];
	symphas::lib::to_upper(str, model);
}



void SystemConf::parse_width(const char* str)
{
	auto widths = symphas::conf::parse_options(str);
	if (widths.size() > intervals_len)
	{
		extend_intervals(widths.size());

		for (iter_type i = 0; i < widths.size(); ++i)
		{
			auto values = symphas::conf::parse_options(widths[i].c_str(), true);

			if (values.size() > 0)
			{
				for (iter_type d = 0; d < values.size(); ++d)
				{
					auto& interval = intervals[i].at(symphas::index_to_axis(d));
					char* end;
					double width = strtod(values[d].c_str(), &end);

					if (*end)
					{
						fprintf(SYMPHAS_WARN, "there was an error interpreting the width "
							"in '%s'", widths[i].c_str());
						width = 1.0;
					}
					interval.set_domain(interval.left(), interval.right(), width);
				}
				for (iter_type d = static_cast<iter_type>(values.size()); d < dimension; ++d)
				{
					auto& interval0 = intervals[i].at(symphas::index_to_axis(0));
					auto& interval = intervals[i].at(symphas::index_to_axis(d));
					interval.set_domain(interval.left(), interval.right(), interval0.width());
				}
			}
		}
	}
}


void SystemConf::parse_dt(const char* str)
{
	auto dts = symphas::conf::parse_options(str, ",");
	if (dts.size() > 0)
	{
		if (dts.size() == 1)
		{
			double dt;
			if (sscanf(dts.front().c_str(), "%lf", &dt) != 1)
			{
				fprintf(SYMPHAS_WARN, "ignoring the time step specification '%s' that is not in the correct "
					"format, expected `dt`\n", dts.front().c_str());
			}
			else
			{
				dt_list.set_time_step(dt, 0, false);
			}
		}
		else
		{
			dt_list.clear_time_steps(1.0);
			double dt, time;
			for (const auto& spec : dts)
			{
				if (sscanf(spec.c_str(), "%lf " STR(CONFIG_OPTION_PREFIX) " %lf", &dt, &time) != 2)
				{
					fprintf(SYMPHAS_WARN, "ignoring the time step specification '%s' that is not in the correct "
						"format, expected `dt`\n", spec.c_str());
				}
				else
				{
					dt_list.set_time_step(dt, time);
				}
			}
		}
	}
}


void swap(SystemConf& first, SystemConf& second)
{
	using std::swap;

	swap(first.g, second.g);
	swap(first.dt_list, second.dt_list);
	swap(first.dimension, second.dimension);
	swap(first.stp, second.stp);
	swap(first.runs, second.runs);
	swap(first.save, second.save);
	swap(first.dims, second.dims);

	swap(first.intervals, second.intervals);
	swap(first.bdata, second.bdata);
	swap(first.tdata, second.tdata);

	swap(first.root_dir, second.root_dir);
	swap(first.result_dir, second.result_dir);
	swap(first.title, second.title);
	swap(first.model, second.model);
	swap(first.names, second.names);
	swap(first.coeff, second.coeff);
	swap(first.num_fields, second.num_fields);

	swap(first.intervals_len, second.intervals_len);
	swap(first.bdata_len, second.bdata_len);
	swap(first.tdata_len, second.tdata_len);
	swap(first.names_len, second.names_len);
	swap(first.coeff_len, second.coeff_len);
	swap(first.num_fields_len, second.num_fields_len);

	swap(first.modifiers, second.modifiers);
	swap(first.init_coeff_copied, second.init_coeff_copied);
}




inline void print_model_coeff(char* out, const double* coeff, iter_type i)
{
	if (coeff[i] * 100 == int(coeff[i] * 100))
	{
		sprintf(out, "[%d]=%.2lf ", i, coeff[i]);
	}
	else if (coeff[i] / 100. >= 1. || 1. / coeff[i] >= 100.)
	{
		sprintf(out, "[%d]=%" DATA_OUTPUT_ACCURACY_STR "E ", i, coeff[i]);
	}
	else
	{
		sprintf(out, "[%d]=%" DATA_OUTPUT_ACCURACY_STR "lf ", i, coeff[i]);
	}
}

inline void print_coeff(FILE* out, const double* coeff, iter_type i)
{
	if (coeff[i] == int(coeff[i]))
	{
		fprintf(out, "%d ", (int)coeff[i]);
	}
	else if (coeff[i] / 100. >= 1. || 1. / coeff[i] >= 100.)
	{
		fprintf(out, "%" DATA_OUTPUT_ACCURACY_STR "E ", coeff[i]);
	}
	else
	{
		fprintf(out, "%" DATA_OUTPUT_ACCURACY_STR "lf ", coeff[i]);
	}
}



void SystemConf::write(const char* savedir, const char* name) const
{
	char* model_name;

	if (model != nullptr)
	{
		auto sep_it = std::strchr(model, VIRTUAL_MODEL_SEP_KEY_C);

		if (sep_it == NULL)
		{
			model_name = new char[std::strlen(model) + 1] {};
			std::strcpy(model_name, model);
		}
		else
		{
			model_name = new char[sep_it - model + 1] {};
			std::copy(model, sep_it, model_name);
			model_name[sep_it - model] = '\0';
		}
	}
	else
	{
		model_name = new char[1] {};
	}
	

	char param_file[BUFFER_LENGTH]{};
	snprintf(param_file, BUFFER_LENGTH, "%s/%s.constants", savedir, model_name);


	FILE* pf;
	if ((pf = fopen(param_file, "a")) == 0)
	{
		fprintf(SYMPHAS_ERR, "error opening write configuration parameters file '%s'\n", param_file);
		exit(ERR_CODE_FILE_OPEN);
	}


	char model_spec[BUFFER_LENGTH_L4]{};
	snprintf(model_spec, BUFFER_LENGTH_L4, "%s ", model_name);
	delete[] model_name;

	// print to the parameters file and the model specification line
	for (iter_type i = 0; i < coeff_len; ++i)
	{
		char coeff_spec[BUFFER_LENGTH_R2];
		print_model_coeff(coeff_spec, coeff, i);

		fprintf(pf, "%s\n", coeff_spec);
		if (std::strlen(model_spec) + std::strlen(coeff_spec) < BUFFER_LENGTH_L4)
		{
			std::strcat(model_spec, coeff_spec);
		}
		else
		{
			fprintf(SYMPHAS_WARN, "printing configuration could not include coefficient index '%d'\n", i);
		}

	}
	fclose(pf);


	if (num_fields_len > 0)
	{
		std::strcat(model_spec, CONFIG_TITLE_PREFIX);
		for (iter_type i = 0; i < num_fields_len; ++i)
		{
			char modifier_spec[BUFFER_LENGTH_R2];
			if (i == num_fields_len - 1)
			{
				sprintf(modifier_spec, " %d %s", num_fields[i], modifiers[i]);
			}
			else
			{
				sprintf(modifier_spec, " %d %s,", num_fields[i], modifiers[i]);
			}
			std::strcat(model_spec, modifier_spec);
		}
	}

	FILE* f;
	char bname[BUFFER_LENGTH]{};
	sprintf(bname, "%s/%s" "." CONFIG_EXTENSION, savedir, name);
	if ((f = fopen(bname, "w")) == 0)
	{
		fprintf(SYMPHAS_ERR, "error opening write configuration file '%s'\n", bname);
		exit(ERR_CODE_FILE_OPEN);
	}

	fprintf(f, "%c%s\n", CONFIG_TITLE_PREFIX_C, title);

	char open = symphas::internal::option_open_bracket();
	char close = symphas::internal::option_close_bracket();

	fprintf(f, CONFIG_NAME_FMT "%zd\n", symphas::internal::C_DIM, dimension);
	fprintf(f, CONFIG_NAME_FMT "%I32d\n", symphas::internal::C_ORDER, stp.ord);
	fprintf(f, CONFIG_NAME_FMT "%I32d\n", symphas::internal::C_PTL, stp.ptl);
	fprintf(f, CONFIG_NAME_FMT "%I32d\n", symphas::internal::C_PTB, stp.ptb);
	fprintf(f, CONFIG_NAME_FMT "%I32d\n", symphas::internal::C_PTG, stp.ptg);
	fprintf(f, CONFIG_NAME_FMT "%s\n", symphas::internal::C_FORM, "");
	fprintf(f, CONFIG_NAME_FMT "%I32d\n", symphas::internal::C_SOLVERVAR, stp.type);


	fprintf(f, CONFIG_NAME_FMT, symphas::internal::C_WIDTH);
	for (iter_type i = 0; i < intervals_len; ++i)
	{
		fprintf(f, "%c ", open);
		for (iter_type d = 0; d < dimension; ++d)
		{
			auto& interval = intervals[i].at(symphas::index_to_axis(d));
			double width =  interval.width();
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
		{Side::BACK, symphas::internal::C_BNDBK}
	};


	for (iter_type side_index = 0; side_index < dimension * 2; ++side_index)
	{
		Side side = symphas::index_to_side(side_index);

		fprintf(f, CONFIG_NAME_FMT, side_key_map[side]);
		for (iter_type i = 0; i < bdata_len; ++i)
		{
			fprintf(f, "%c %s ", open, symphas::str_from_boundary(bdata[i].at(side).type));

			for (iter_type tag_index = 0; tag_index < 2; ++tag_index)
			{
				const char* tag_name = symphas::str_from_boundary_tag(bdata[i].at(side).tag[tag_index]);
				if (tag_name)
				{
					fprintf(f, "%s ", tag_name);
				}
			}
			for (iter_type n = 0; n < bdata[i].at(side).argc; ++n)
			{
				fprintf(f, "%.8E ", bdata[i].at(side).params[n]);
			}
			fprintf(f, "%c", close);
		}
		fprintf(f, "\n");
	}


	const char interval_fmt[] = "%c %lf %lf %c";
	if (dimension > 0)
	{
		fprintf(f, CONFIG_NAME_FMT, symphas::internal::C_RNGX);
		for (iter_type i = 0; i < intervals_len; ++i)
		{
			fprintf(f, interval_fmt, open, DOMAIN_X0_AT(i), DOMAIN_Xn_AT(i), close);
		}
		fprintf(f, "\n");
	}
	if (dimension > 1)
	{
		fprintf(f, CONFIG_NAME_FMT, symphas::internal::C_RNGY);
		for (iter_type i = 0; i < intervals_len; ++i)
		{
			fprintf(f, interval_fmt, open, DOMAIN_Y0_AT(i), DOMAIN_Yn_AT(i), close);
		}
		fprintf(f, "\n");
	}
	if (dimension > 2)
	{
		fprintf(f, CONFIG_NAME_FMT, symphas::internal::C_RNGZ);
		for (iter_type i = 0; i < intervals_len; ++i)
		{
			fprintf(f, interval_fmt, open, DOMAIN_Z0_AT(i), DOMAIN_Zn_AT(i), close);
		}
		fprintf(f, "\n");
	}

	fprintf(f, CONFIG_NAME_FMT, symphas::internal::C_INSIDE);
	for (iter_type i = 0; i < tdata_len; ++i)
	{
		fprintf(f, "%c ", open);
		for (auto const& [key, entry] : tdata[i])
		{
			if (tdata[i].size() > 1 && key != Axis::NONE)
			{
				fprintf(f, "%c%s ", CONFIG_TITLE_PREFIX_C,
					(key == Axis::X) ? "X" :
					(key == Axis::Y) ? "Y" :
					(key == Axis::Z) ? "Z" :
					(key == Axis::R) ? "R" :
					(key == Axis::S) ? "S" :
					(key == Axis::T) ? "T" : "?");
			}

			if (entry.in != Inside::NONE && ((tdata[i].size() > 1 && key != Axis::NONE) || (tdata[i].size() == 1)))
			{
				if (entry.in == Inside::EXPRESSION)
				{
					fprintf(f, "%s %s ", symphas::str_from_in(entry.in), entry.expr_data.get_name());
					if (!init_coeff_copied)
					{
						for (iter_type n = 0; n < entry.expr_data.get_num_coeff(); ++n)
						{
							print_coeff(f, entry.expr_data.get_coeff(), n);
						}
					}
					else
					{
						fprintf(f, STR(CONFIG_OPTION_PREFIX));
					}
				}
				else if (entry.in == Inside::FILE || entry.in == Inside::CHECKPOINT)
				{
					fprintf(f, STR(CONFIG_OPTION_PREFIX) " %s %d", entry.file.get_name(), entry.file.get_index());
				}
				else
				{
					fprintf(f, "%s ", symphas::str_from_in(entry.in));

					size_t tag = entry.intag;
					size_t pos = 0;
					while (tag > 0)
					{
						bool is_bit_set = (tag & (1ull << pos)) > 0;
						if (is_bit_set)
						{
							tag = (tag & ~(1ull << pos));
							const char* tag_name = symphas::str_from_in_tag(static_cast<InsideTag>(pos));
							if (tag_name)
							{
								fprintf(f, "%s ", tag_name);
							}
						}
						pos += 1;

					}
					for (iter_type n = 0; n < NUM_INIT_CONSTANTS; ++n)
					{
						print_coeff(f, entry.data.gp, n);
					}
				}
			}
		}

		fprintf(f, "%c", close);
	}
	fprintf(f, "\n");

	fprintf(f, CONFIG_NAME_FMT, symphas::internal::C_DELTA);
	if (dt_list.get_num_time_steps() > 1)
	{
		for (auto [time, dt] : dt_list)
		{
			fprintf(f, "%c %lf " STR(CONFIG_OPTION_PREFIX) " %lf %c ", open, dt, time, close);
		}
		fprintf(f, "\n");
	}
	else
	{
		fprintf(f, "%lf\n", dt_list.get_time_step());
	}

	if (model != nullptr)
	{
		char* sep_it;
		if ((sep_it = std::strchr(model, VIRTUAL_MODEL_SEP_KEY_C)) == NULL)
		{
			fprintf(f, CONFIG_NAME_FMT "%s\n", symphas::internal::C_MODEL, model_spec);
		}
		else
		{
			char* model_cpy = new char[std::strlen(sep_it)] {};
			auto index_it = std::strchr(sep_it + 1, VIRTUAL_MODEL_SEP_KEY_C);
			auto end_it = (index_it == NULL) ? sep_it + std::strlen(sep_it) : index_it;

			std::copy(sep_it + 1, end_it, model_cpy);
			model_cpy[end_it - sep_it - 1] = '\0';
			fprintf(f, CONFIG_NAME_FMT "%c %s", symphas::internal::C_MODEL, CONFIG_TITLE_PREFIX_C, model_cpy);

			if (index_it != NULL)
			{
				iter_type index;
				sscanf(index_it + 1, "%d", &index);
				fprintf(f, " %d", index);
				if ((index_it = std::strchr(index_it + 1, VIRTUAL_MODEL_SEP_KEY_C)) != NULL)
				{
					sscanf(index_it + 1, "%d", &index);
					fprintf(f, " %d", index);
				}
			}
			fprintf(f, "\n");

			delete[] model_cpy;
		}
	}

	fprintf(f, CONFIG_NAME_FMT "%d\n", symphas::internal::C_FRAMES, save.get_stop());

	char save_spec[BUFFER_LENGTH];
	symphas::io::save_as_str(&save, save_spec, BUFFER_LENGTH);
	fprintf(f, CONFIG_NAME_FMT "%s\n", symphas::internal::C_SAVE, save_spec);


	fprintf(f, CONFIG_NAME_FMT "%s\n", symphas::internal::C_SAVEINIT, save.get_init_flag() ? "YES" : "NO");
	fprintf(f, CONFIG_NAME_FMT "%zd\n", symphas::internal::C_RUNS, runs);
	fprintf(f, CONFIG_NAME_FMT "%s\n", symphas::internal::C_DIR, root_dir);

	fprintf(f, "\n");

	/* also copy into the backup configuration the parameters that were given
	 * in the last run
	 */
	for (auto s : params::rawparams)
	{
		if (s.second.size() == 0)
		{
			fprintf(f, CONFIG_PARAM_PREFIX "%s\n", s.first.c_str());
		}
		else
		{
			fprintf(f, CONFIG_PARAM_PREFIX "%s=%s\n", s.first.c_str(), s.second.c_str());
		}
	}
	fprintf(f, "\n");
	fclose(f);

}
