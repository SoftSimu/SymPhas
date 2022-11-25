
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
 * MODULE:  conf
 * PURPOSE: Defines the configuration functionality for defining parameters
 * of a phase field system.
 *
 * ***************************************************************************
 */

#pragma once

#include <string>
#include <stack>
#include <functional>

#include "io.h"
#include "stencilincludes.h"
#include "systemlib.h"




//! The extension of the configuration file for SymPhas.
#define CONFIG_EXTENSION "in"






//! Key used to prefix comments in the configuration file.
#define CONFIG_COMMENT_PREFIX "#"

//! The character usage of the comment key.
#define CONFIG_COMMENT_PREFIX_C ((CONFIG_COMMENT_PREFIX)[0])


//! Key used to prefix the title in the configuration file.
/*!
 * The title is typically specified at the top of the configuration file,
 * but can appear anywhere.
 */
#define CONFIG_TITLE_PREFIX "!"

//! The character usage of the title key.
#define CONFIG_TITLE_PREFIX_C ((CONFIG_TITLE_PREFIX)[0])


//! The key which separates keys and values in the configuration.
#define CONFIG_SEP_KEY ":"

//! The character usage of the separation key.
#define CONFIG_SEP_KEY_C ((CONFIG_SEP_KEY)[0])

//! The way the keys are formatted in the configuration.
#define CONFIG_NAME_FMT "%-12s" CONFIG_SEP_KEY

//! Parameters are added to the configuration using this format.
/*!
 * Parameters are added to the configuration using this format. This is
 * typically only when a backup configuration is written that needs to persist
 * the command line parameters passed to the program.
 */
#define CONFIG_PARAM_PREFIX CONFIG_COMMENT_PREFIX "*"

//! The index at which backups are written are saved with this format.
/*!
 * When backup data is written, the backup configuration file contents are
 * appended to add the current index. The index that is written needs to be
 * preceded by a specific key.
 */
#define CONFIG_INDEX_PREFIX CONFIG_COMMENT_PREFIX "$"

//! The character used to prefix special options in the configuration.
#define CONFIG_OPTION_PREFIX @
//! Character usage of #CONFIG_OPTION_PREFIX
#define CONFIG_OPTION_PREFIX_C (STR(CONFIG_OPTION_PREFIX)[0])

//! Extension given to the file specifying the coefficients.
#define COEFF_SPEC_EXTENSION "constants"

#define COEFF_STR_PREPEND_NAME "c"

//! The format string for a default phase field name.
/*!
 * The format for the default phase field name string. It only takes
 * an index format specification.
 */
#define DEFAULT_PF_NAME_FMT "Phase Field %zd"



namespace symphas
{
	//! Specifies elements used in the management of SymPhas configuration.
	/*!
	 * Specifies elements used in the management of SymPhas configuration.
	 */
	namespace conf {}
}

namespace symphas::internal
{
	char option_open_bracket();
	char option_close_bracket();
}

namespace symphas::conf
{

	std::vector<std::string> parse_options(const char* options, bool spaces_are_delimiters = false, const char* extra_delimiters = "");
	inline std::vector<std::string> parse_options(const char* options, const char* extra_delimiters)
	{
		return parse_options(options, false, extra_delimiters);
	}



	//! Append default phase field names starting at the given pointer.
	/*!
	 * The iterator to the names list is provided, which is expected to point to the
	 * start of an array of at least length `num`. The array is then filled with 
	 * that many default names, with starting index which may be chosen to be 
	 * greater than 0.
	 * 
	 * \param names_it The iterator to the position in the names list to begin 
	 * inserting the default names.
	 * \param num The number of default names to insert, and the expected length of
	 * the array.
	 * \param start_id The starting index for the default names, which is 0 by 
	 * default.
	 */
	inline void insert_default_pf_names(char** names_it, size_t num, size_t start_id = 0)
	{
		for (size_t id = start_id; id < start_id + num; ++id)
		{
			char name_buffer[BUFFER_LENGTH_R2];
			sprintf(name_buffer, DEFAULT_PF_NAME_FMT, id);
			names_it[id] = new char[std::strlen(name_buffer) + 1];
			std::strcpy(names_it[id], name_buffer);
		}

	}

	 //! Expand memory and append default names to the given names list.
	 /*
	  * Given a list of names of phase fields, append an additional number of names
	  * given by num, in the default phase field name format. The default format
	  * is given by #DEFAULT_PF_NAME_FMT. The list is expanded to fit the additional
	  * names.
	  * 
	  * \param names The list of strings with the names of the phase fields, which
	  * will be appended with num default names. It will be reallocated to the
	  * expanded length to accommodate additional names.
	  * \param names_len The current length of the names list.
	  * \param num The number of default names to append.
	  */
	inline void append_default_pf_names(char** &names, size_t names_len, size_t num)
	{
		char** extend = new char* [names_len + num];
		for (iter_type i = 0; i < names_len; ++i)
		{
			extend[i] = names[i];
		}

		insert_default_pf_names(extend + names_len, num);

		delete[] names;
		names = extend;
	}
}


// ***************************************************************************

//! Object storing configurable properties.
/*!
 * The configuration specifies a number of properties which could be used 
 * throughout the workflow of a simulation. The properties are entirely managed 
 * by this object, so  for example, copies of a configuration can be made in 
 * which case all the properties will be correctly constructed and assigned.
 * Also allows the configuration to be persisted.
 * 
 * Most of the members of the configuration are publicly accessible, and these
 * members can safely be modified. Members which are which are accessed through
 * function interfaces include the boundaries, initial conditions and
 * system intervals (these are arrays, such that each element represents
 * the description of one field, initialized in the order of the array).
 * 
 * The configuration has a default initialization, which does not represent
 * a valid set of properties of a phase field problem.
 * 
 * The parameters that are not detailed in the function docs are mentioned here:
 * 
 * |Property            | Configuration Key             | Format |
 * |--------------------|-------------------------------|--------|
 * | Dimension          | `DIM` | Either `1`, `2` or `3`. |
 * | FD Stencil Order   | `ORDER` | Either `2` or `4`. |
 * | `N`-point Laplacian Stencil | `PTL` | Refer to available stencils. |
 * | `N`-point Biaplacian Stencil | `PTB` | Refer to available stencils. |
 * | `N`-point Gradlaplacian Stencil | `PTG` | Refer to available stencils. |
 * | Time Step (\f$\Delta t\f$)  | `DELTA` | Any positive floating point number. |
 * | Number of Solution Iterations  | `FRAMES` | Any positive integer. |
 * | Number of Concurrent Simulations  | `FRAMES` | Any positive integer. |
 * | Directory of Results | `DIR` | Any valid system path (relative or absolute). |
 * | Save Interval | `SAVE` | See symphas::io::parse_save_str(). |
 * | Save Index 0? | `SAVEINIT` | `YES` or `NO`. |
 *
 * The format and key of each possible configurable property are described
 * for each function. The general configuration formatting rules are the
 * following:
 * - Each configuration property line consists of the key, the character `:`
 * and then the value of the property. The `:` character separates the key 
 * from the value, and is called the _option delimiter_.
 * - Each item in a specification is separated by any number of spaces or tabs.
 * - Numeric values can be specified using scientific notation.
 * - There is a character used for special configuration options, `@`. This
 * character is called the _detail specifier_. It changes the behaviour
 * of some properties, which are described in their respective format
 * explanations.
 * - The title of the simulation is specified on its own line, preceded by
 * the character `!`, the _title specifier_.
 * 
 * Some configurable properties support multiple options, where
 * options are simply bracket surrounded elements. Each element of the options
 * follows the format of the property it is written for. Typically, options are
 * given to specify properties about multiple phase fields. This is
 * indicated in the format specification of the respective property.
 * 
 * Options may be surrounded in either brackets, `(` and `)`, or braces,
 * `{` and `}`. Options may also be surrounded in quotes, `"` or `'`. When
 * options are surrounded in brackets, nested options can be used.
 * 
 * Additionally, the configuration property key is case insensitive. Any
 * number of spaces and tabs may separate the key and the option delimiter, and
 * any number of spaces may separate the option delimiter and the value. 
 * Moreover, as long as the key is followed by at least one tab or space,
 * then anything may be written up to the option delimiter (for example,
 * to provide a hint as to the meaning of that key).
 *
 */
struct SystemConf
{
	size_t dimension;							//!< The dimension of the system, can be 1 or 2.
	double dt;									//!< The temporal width.
	StencilParams stp;							//!< Parameters characterizing the stencils.
	size_t runs;								//!< The number of systems that will be solved simultaneously.
	SaveParams save;							//!< Encapsulation of save parameters.
	Geometry g;									//!< The type of coordinate system of the simulated problem.


protected:

	char* root_dir;								//!< The directory to save all the simulation information.
	char* result_dir;							//!< The directory to save all the simulation information.
	char* title;								//!< The name of the simulation.
	char* model;								//!< The string representing the model which is simulated.
	char** names;								//!< A list of names of each of the phase field.

	double* coeff;								//!< The coefficients in the equations of motion.


	symphas::interval_data_type* intervals;		//!< A list of the intervals of all the axes of the system.
	symphas::b_data_type* bdata;				//!< A list of the boundary data corresponding to all sides of the system.
	symphas::init_data_type* tdata;				//!< A list of the parameters required for generating the initial data.
	len_type** dims;							//!< An array of the full dimensions of the system, instead of single variables.

	size_t intervals_len;						//!< The number of interval data elements provided.
	size_t bdata_len;							//!< The number of boundary data elements provided.
	size_t tdata_len;							//!< The number of initial condition elements provided.
	size_t names_len;							//!< The number of phase field names provided.
	size_t coeff_len;							//!< The number of coefficients provided.

private:

	bool init_coeff_copied;

public:


	SystemConf();
	SystemConf(symphas::problem_parameters_type const& parameters, const char* title, const char* dir = "");

	//! Generate the configuration from the given options.
	/*!
	 * Generate the configuration from the given options. The name of the
	 * configuration is also provided.
	 * 
	 * \param options The list of key-value pairs for all the options used
	 * in the configuration.
	 * \param title The title of the configuration.
	 * \param dir The absolute parent directory of the configuration, which is
	 * optional.
	 */
	SystemConf(std::vector<std::pair<std::string, std::string>> options, const char* title, const char* dir = "");
	SystemConf(SystemConf const& other);
	SystemConf(SystemConf&& other) noexcept;
	SystemConf& operator=(SystemConf other);

	friend void swap(SystemConf& first, SystemConf& second);


	~SystemConf()
	{
		delete[] root_dir;
		delete[] result_dir;
		delete[] title;
		delete[] model;
		delete[] coeff;

		delete[] intervals;
		delete[] bdata;
		delete[] tdata;

		for (iter_type i = 0; i < intervals_len; ++i)
		{
			delete[] dims[i];
		}
		delete[] dims;

		for (iter_type i = 0; i < names_len; ++i)
		{
			delete[] names[i];
		}
		delete[] names;
	}


	/* **********************************************************************
	 * Configuration functions to set members from strings.
	 * **********************************************************************/

	//! Sets the coordinate geometry.
	/*!
	 * Parses the given string and sets the parameter SystemConf::g based on
	 * the result. The default geometry is Geometry::CARTESIAN.
	 * 
	 * > Configuration key = `FORM`
	 *
	 * The two valid string options to this parameter are:
	 *  - `CARTESIAN`
	 *  - `POLAR`
	 *
	 * \param str The string defining the desired system geometry.
	 */
	void parse_form(const char* str);

	//! Sets the boundary specifications.
	/*!
	 * The input is parsed to set the member SystemConf::bdata. Multiple options
	 * may be defined for each of the strings given in the string list. The
	 * length of the input list is expected to be twice SystemConf::dimension
	 * (because the number of boundaries is twice the dimension).
	 * 
	 * Each element in the list represents the boundary specification of
	 * one or more systems for a single side. The side that the specification
	 * index corresponds to is given by symphas::index_to_side().
	 * The input list does not require each element to specify 
	 * an equal number of options, but the boundaries which are not
	 * specified will be copied from the first specification in that list.
	 *
	 * For each of the strings in the given list, the content is first parsed
	 * for one or more options (boundaries for more than one system). The result
	 * of the parsing is a new list containing boundary specifications for as
	 * many systems as options provided. The options are used to
	 * form a new list in order to initialize the boundary specification.
	 * 
	 * > Configuration keys = `BNDLT` (Side::LEFT), `BNDRT` (Side::RIGHT), 
	 * > `BNDTP` (Side::TOP), `BNDBT` (Side::BOTTOM), `BNDFT` (Side::FRONT), `BNDBK` (Side::BACK)
	 * 
	 * See SystemConf::parse_boundaries() for formatting details of individual
	 * options.
	 *
	 * \param str_list The list of strings that define the boundary
	 * specifications. Expected to be initialized to a list of length
	 * equal to twice the dimension given by SystemConf::dimension. The
	 * boundary with #Side side is set by the element in the string list
	 * with index given by symphas::side_to_index(Side).
	 */
	void parse_boundaries_array(const char* const* str_list);

	//! Sets the specification for the nth boundary data element.
	/*!
	 * From the list of strings, which should be of length equal to twice the
	 * dimension prescribed by SystemConf::dimension, each element is parsed to
	 * set the boundary specification index for the boundary data element
	 * corresponding to the side given by symphas::index_to_side().
	 * 
	 * > Configuration keys = `BNDLT` (Side::LEFT), `BNDRT` (Side::RIGHT), 
	 * > `BNDTP` (Side::TOP), `BNDBT` (Side::BOTTOM), `BNDFT` (Side::FRONT), `BNDBK` (Side::BACK)
	 * 
	 * The format of each item in the list is:
	 * `KEY MODIFIER... VALUE...`, where zero, one or two modifiers are
	 * provided and zero or more values are provided. The
	 * key is the boundary type, one of either `PERIODIC`, `DEFAULT` or `OPEN`.
	 * This is described in more detail in ::BoundaryType.
	 * The modifiers are tags which affect the generation algorithm, and this
	 * is described in more detail in ::BoundaryTag. The additional values
	 * affect the values generated by the algorithm, and are algorithm
	 * specific. Refer to symphas::internal::new_boundary().
	 *
	 * The `KEY` is any key in symphas::internal::boundary_key_map:
	 * - `DEFAULT` BoundaryType::DEFAULT
	 * - `PERIODIC` BoundaryType::PERIODIC
	 * 
	 * If the key is not specified, then the tag needs to be specified, in which
	 * case the boundary type is automatically BoundaryType::DEFAULT.
	 * 
	 * The tag is a modifier and is any of the keys in 
	 * symphas::internal::boundary_tag_key_map, which are:
	 * - `GAUSSIAN` BoundaryTag::GAUSSIAN
	 * - `GA` BoundaryTag::GAUSSIAN
	 * - `RANDOM` BoundaryTag::RANDOM
	 * - `RA` BoundaryTag::RANDOM
	 * - `TRIGONOMETRIC` BoundaryTag::TRIG
	 * - `TR` BoundaryTag::TRIG
	 * - `CONSTANT` BoundaryTag::CONSTANT
	 * - `CO` BoundaryTag::CONSTANT
	 * - `LINEAR` BoundaryTag::LINEAR
	 * - `LI` BoundaryTag::LINEAR
	 * - `TIME` BoundaryTag::TIME
	 * - `T` BoundaryTag::TIME
	 * 
	 * These are parsed by symphas::boundary_from_str().
	 *
	 * \param str_list A list where each element defines a boundary
	 * configuration. The boundary at #Side side has its configuration given by
	 * the element in the list at index symphas::side_to_index(Side).
	 * \param n The index of SystemConf::bdata which to set.
	 */
	void parse_boundaries(const char* const* str_list, size_t n);

	//! Sets one or more initial conditions.
	/*!
	 * The string is first parsed for one or more options. Each option is then
	 * parsed in turn by SystemConf::set_in(const char*, size_t) to initialize
	 * the initial data. This initializes one or more initial conditions,
	 * applied to the phase field problem fields in order of definition.
	 *
	 * > Configuration key = `INSIDE`.
	 *
	 * See SystemConf::parse_initial_condition() for the format specification
	 * of each option.
	 * 
	 * \param str The input which contains the specification for one or more
	 * initial conditions.
	 */
	void parse_initial_condition_array(const char* str);

	//! Sets the initial condition at the given index.
	/*!
	 * The provided string configuration is parsed according to a standard
	 * format, and the initial condition for the phase field system at given
	 * index is initialized. Also an input file can be passed.
	 * 
	 * # Passing an initial condition type
	 * 
	 * > Configuration key = `INSIDE`.
	 * 
	 * The format of the specification is `KEY MODIFIER... VALUE...`, where
	 * there are zero or more modifiers provided and
	 * zero or more values are provided. The `KEY` is any of the
	 * initial conditions specified in symphas::internal::init_key_map. These
	 * are:
	 * - `GAUSSIAN` Inside::GAUSSIAN
	 * - `GA` Inside::GAUSSIAN
	 * - `UNIFORM` Inside::UNIFORM
	 * - `UN` Inside::UNIFORM
	 * - `CAPPED` Inside::CAPPED
	 * - `CA` Inside::CAPPED
	 * - `CONSTANT` Inside::CONSTANT
	 * - `CO` Inside::CONSTANT
	 * - `CIRCLE` Inside::CIRCLE
	 * - `CI` Inside::CIRCLE
	 * - `SQUARE` Inside::SQUARE
	 * - `SQ` Inside::SQUARE
	 * - `HEXAGONAL` Inside::HEXAGONAL
	 * - `HX` Inside::HEXAGONAL
	 * - `CUBIC` Inside::CUBIC
	 * - `CU` Inside::CUBIC
	 * - `SEEDSSQUARE` Inside::SEEDSSQUARE
	 * - `SS` Inside::SEEDSSQUARE
	 * - `SEEDSCIRCLE` Inside::SEEDSCIRCLE
	 * - `CS` Inside::SEEDSCIRCLE
	 * 
	 * The `MODIFIER` is specified in symphas::internal::init_tag_key_map, and 
	 * is any of:
	 * - `DEFAULT` InsideTag::DEFAULT
	 * - `RANDOM` InsideTag::RANDOM
	 * - `INVERT` InsideTag::INVERT
	 * - `INV` InsideTag::INVERT
	 * - `VARA` InsideTag::VARA
	 * - `A` InsideTag::VARA
	 * - `VARB` InsideTag::VARB
	 * - `B` InsideTag::VARB
	 * 
	 * # Passing a file with initial conditions
	 * 
	 * When the entry starts with the special character `@`, then it is assumed
	 * a file of initial conditions is being provided. The name of the file is provided
	 * immediately after `@` with no spaces in between. If there are multiple files, then
	 * each file needs to be preceded by `@`. If the file has spaces, then quotes or brackets
	 * must surround the whole entry, including the @. If `!` character is provided  
	 * first instead, then the parameter variable #input_data_file will be read
	 * in order to get the file for the initial conditions.
	 *
	 * \param str The string detailing the initial condition specification.
	 * \param n The index in the interval list which is initialized.
	 */
	void parse_initial_condition(const char* str, size_t n);

	//! For a given axis, set one or more corresponding intervals.
	/*!
	 * The string specifies the interval of more or one systems. One or more
	 * interval data sets are given as options following the option format.
	 * 
	 * > Configuration keys = `RNGX` (Axis::X), `RNGY` (Axis::Y), `RNGZ` (Axis::Z)
	 * 
	 * The format of each option is given by SystemConf::parse_interval().
	 * 
	 * Once all the intervals are specified, 
	 *
	 * \param str The specification string giving the interval information
	 * across multiple systems.
	 */
	void parse_interval_array(const char* const* str);

	//! Initializes the interval over the given axis.
	/*!
	 * Given the specification string, the interval endpoints are initialized
	 * for the given axis.
	 * 
	 * > Configuration keys = `RNGX` (Axis::X), `RNGY` (Axis::Y), `RNGZ` (Axis::Z)
	 * 
	 * Interval follows the format `START END`, where both values are real
	 * numbers, separated by a space. The interval can also be given
	 * by `@ N`, where `N` is the number of discrete elements to put in
	 * the interval. The interval `START` value will then always be zero, and
	 * the `END` will be computed based on the width.
	 *
	 * \param str The specification string giving the interval information.
	 * \param ax The axis of the interval.
	 * \param n The index in the interval information array.
	 */
	void parse_interval(const char* str, Axis ax, size_t n);

	//! Initialize the spacing along the axes
	/*!
	 * The string may define a number of values up to the size of the dimension.
	 * Each of the values then corresponds with the spacing along that axis,
	 * unless there are fewer values than the size of the dimension provided,
	 * in which case, the rest of the values will be initialized based
	 * on the first value.
	 * 
	 * The format follows: `VALUE...`, where there is at least one value, and
	 * at most, the number of values is the size of the dimension. When multiple 
	 * options are specified, the spatial widths will be set individually
	 * by system.
	 * 
	 * \param str The specification of the spatial width along the axes.
	 */
	void parse_width(const char* str);


	//! Initializes the names of the phase field systems.
	/*!
	 * Initializes the names of the phase field systems. Can specify any number
	 * of names, and the systems are named in the order they are specified in
	 * the model definition. Names are individual options separated by
	 * brackets or quotes.
	 * 
	 * Note that this does not name the variables, as this is done through the 
	 * macro definition of the model. 
	 * 
	 * The phase field names are used in the plotting configuration, in
	 * order to assign names to the plots.
	 *
	 * \param str The list of phase field system names.
	 */
	void parse_names(const char* str);


	//! Sets the coefficients and model string.
	/*!
	 * The provided input defines a model name and the coefficients of that 
	 * model. The string name of the model is
	 * associated with the linked name given in the model definition, see
	 * #LINK_WITH_NAME.
	 *
	 * The format of the input is `MODEL VALUES...`, where `MODEL` is the
	 * key name of the model, and `VALUES...` are zero or more numbers given
	 * as coefficients to the model.
	 * 
	 * The format of each value is not required to be numeric, but can also
	 * follow the format used in specifying coefficients in a file, described
	 * next.
	 *
	 * A variation of the format is `MODEL @ FILE`, where `MODEL` is the same,
	 * but `FILE` is the name of a file containing a list of coefficients in
	 * a specific format. Each line in the file follows one of the following
	 * formats:
	 * 1. `[FORMAT]=VALUE`, the coefficients specified by `FORMAT` 
	 * are initialized to `VALUE`. The `FORMAT` specifier is one or more
	 * indices separated by either a comma (`,`) or two dots (`..`). The index
	 * is the position of the coefficient in a 0-indexed list.
	 *    - If `FORMAT` includes indices separated by `..` as in `N..M`, then 
	 * the  coefficients with indices from `N` to `M`, including coefficient `M` 
	 * itself, are initialized to `VALUE`. 
	 *    - If `FORMAT` includes indices separated by `,` as in `N,M`, then the
	 * two different coefficients at indices `N` and `M` are initialized to 
	 * `VALUE`. 
	 * 2. `LETTER=VALUE`, where `LETTER` is any identifier of the form `c<N>`,
	 * where `<N>` indicates a substituted integer value, which must be greater
	 * than 0. This is using the equivalent form as in the
	 * model definition. The corresponding coefficient is then initialized.
	 * The `LETTER` identifier can also be any of the capital letters from
	 * `A` to `H` or `AA` to `HH`, where the sequence refers to the coefficient
	 * index from 0 to 7 (i.e. the alphanumeric index represented by the
	 * character) plus the offset equal to the length of the sequence (8), 
	 * times the number of repeated characters. Mixing characters is not
	 * a valid format.
	 * 
	 * As an example of using the previous format, `B` is the same as `c2`, is the 
	 * same as `[1]`. Also, `EE` is the same as `c13`, is the same as `[12]`. Two
	 * repetitions of the `E` indicates that `EE` is in the second coefficient
	 * index sequence, so 8 will be added to the alphanumeric index of `E`.
	 *
	 * If this option is used, no other parameters may be
	 * specified (i.e. no coefficients may be additionally specified).
	 * 
	 * In general, not all coefficients need to be specified. Any coefficient
	 * which is not specified is initialized to a default value,
	 * typically `1`. This is given by #DEFAULT_COEFF_VALUE.
	 * 
	 * \param str The written parameter specification.
	 * \param dir The directory to search for additional information,
	 * specifically coefficients files.
	 */
	void parse_model_spec(const char* str, const char* dir);

	//! Set the type of model.
	/*!
	 * The name of the model is set, which is typically associated with the
	 * particular phase field problem.
	 */
	void set_model_name(const char* str);


	//! Set the name at the given index.
	/*!
	 * This is the name of the phase field corresponding to the `n`th
	 * index in the problem.
	 * 
	 * \param name The given name of the phase field.
	 * \param n The index of the phase field that is named.
	 */
	void set_name(const char* name, size_t n);

	//! Set the directory in which to save the results. 
	/*!
	 * Set the directory in which to save the results.
	 */
	void set_directory(const char* directory);

	void select_stencil(size_t order, const char* str);

	//! Set the boundary on the given side for the system at the given index.
	void set_boundary(symphas::b_element_type boundary, Side side, int n)
	{
		if (n < bdata_len)
		{
			bdata[n][side] = boundary;
		}
		else
		{
			auto extend = new symphas::b_data_type[n + 1]{};
			for (iter_type i = 0; i < bdata_len; ++i)
			{
				extend[i] = bdata[i];
			}
			delete[] bdata;
			bdata = extend;
			bdata_len = n + 1;

			set_boundary(boundary, side, n);
		}
	}


	//! Set the initial condition for the system at the given index.
	void set_initial_condition(symphas::init_data_type condition, int n)
	{
		if (n < tdata_len)
		{
			tdata[n] = condition;
		}
		else
		{
			auto extend = new symphas::init_data_type[n + 1]{};
			for (iter_type i = 0; i < tdata_len; ++i)
			{
				extend[i] = tdata[i];
			}
			delete[] tdata;
			tdata = extend;
			tdata_len = n + 1;

			set_initial_condition(condition, n);
		}
	}



	//! Set the interval on the given axis for the system at the given index.
	void set_intervals(symphas::interval_data_type interval, int n)
	{
		if (n < intervals_len)
		{
			intervals[n] = interval;
			set_dimensions(n);
		}
		else
		{
			auto extendv = new symphas::interval_data_type[n + 1]{};
			auto extendm = new len_type * [n] {};
			for (iter_type i = 0; i < intervals_len; ++i)
			{
				extendv[i] = intervals[i];
				extendm[i] = dims[i];
			}

			delete[] intervals;
			delete[] dims;

			intervals = extendv;
			dims = extendm;
			intervals_len = n + 1;

			set_intervals(interval, n);
		}
	}

	//! Set the interval on the given axis for the system at the given index.
	void set_interval(symphas::interval_element_type interval, Axis ax, int n)
	{
		if (n < intervals_len)
		{
			intervals[n][ax] = interval;
			set_dimensions(n);
		}
		else
		{
			set_intervals(intervals[0], n);
			set_interval(interval, ax, n);
		}
	}

	void set_width(double h)
	{
		for (auto i = 0; i < intervals_len; ++i)
		{
			for (auto& [_, interval] : intervals[i])
			{
				interval.set_interval(interval.left(), interval.right(), h);
			}
		}
	}

	template<Axis ax>
	void set_width(double h)
	{
		for (auto i = 0; i < intervals_len; ++i)
		{
            auto& interval = intervals[i].at(ax);
			interval.set_interval(interval.left(), interval.right(), h);
		}
	}


	//! Clears all the coefficients by making the coefficient list zero.
	/*!
	 * Clears all the coefficients by making the coefficient list zero.
	 * Deallocates the coefficient array.
	 */
	void reset_coeff()
	{
		coeff_len = 0;
		delete[] coeff;
		coeff = nullptr;
	}

	//! Initializes the coefficient at the given index.
	/*!
	 * The coefficient at the given index is initialized to the value. If the
	 * index is greater than the number of coefficients, then the list is
	 * expanded.
	 */
	void set_coeff(double c, iter_type n)
	{
		if (n < coeff_len)
		{
			coeff[n] = c;
		}
		else
		{
			double* extend = new double[n + 1];
			for (iter_type i = 0; i < coeff_len; ++i)
			{
				extend[i] = coeff[i];
			}
			for (size_t i = coeff_len; i < n; ++i)
			{
				extend[i] = DEFAULT_COEFF_VALUE;
			}
			extend[n] = c;

			delete[] coeff;
			coeff = extend;
			coeff_len = n + 1;
		}
	}


	//! Get the value of the coefficient.
	double get_coeff(iter_type n) const
	{
		if (n < coeff_len)
		{
			return coeff[n];
		}
		else
		{
			return DEFAULT_COEFF_VALUE;
		}
	}

	//! Get the array of coefficients.
	const double* get_coeff_list() const
	{
		return coeff;
	}

	//! Get the number of defined coefficients.
	size_t coeff_count() const
	{
		return coeff_len;
	}

	//! Get the number of defined phase field names.
	size_t name_count() const
	{
		return names_len;
	}

	//! Get the number of systems this configuration represents.
	/*!
	 * The number of systems that the configuration represents is the
	 * greatest number between the number of boundaries and the number of
	 * intervals.
	 */
	size_t system_count() const
	{
		return std::max(bdata_len, intervals_len);
	}

	//! Return a reference to the boundary.
	/*!
	 * Obtain the boundary condition associated with the field of the given
	 * index, for the given side.
	 * 
	 * \param side The side to which the desired boundary corresponds.
	 * \param n The index of the field the boundary is taken from.
	 */
	const symphas::b_element_type& get_boundary(Side side, int n) const
	{
		if (n < bdata_len)
		{
			if (bdata[n].count(side))
			{
				return bdata[n].at(side);
			}
		}

		throw "there is no such boundary condition";
	}


	//! Return a reference to the initial condition.
	/*!
	 * Obtain the initial condition associated with the field of the given
	 * index.
	 *
	 * \param n The index of the field the boundary is taken from.
	 */
	const symphas::init_data_type& get_initial_condition(int n) const
	{
		if (n < tdata_len)
		{
			return tdata[n];
		}

		throw "there is no such initial condition";
	}

	//! Return a reference to the interval.
	/*!
	 * Obtain the interval associated with the field of the given
	 * index, for the given axis.
	 *
	 * \param ax The axis to which the desired interval corresponds.
	 * \param n The index of the field the boundary is taken from.
	 */
	const symphas::interval_element_type& get_interval(Axis ax, int n) const
	{
		if (n < intervals_len)
		{
			if (intervals[n].count(ax))
			{
				return intervals[n].at(ax);
			}
		}

		throw "there is no such axis";
	}

	const len_type* get_dimensions(int n)
	{
		return dims[n];
	}


	const char* get_name(int n) const
	{
		if (n < names_len)
		{
			return names[n];
		}
		else
		{
			return "";
		}
	}


	const char* get_result_dir() const
	{
		return result_dir;
	}

	const char* get_dir() const
	{
		return root_dir;
	}

	const char* get_title() const
	{
		return title;
	}

	const char* get_model_name() const
	{
		return model;
	}

	//! Return the problem parameters object.
	/*!
	 * Constructs the problem parameters object in order to be used in model
	 * initialization. The problem parameters are defined so that the maximum
	 * amount of information is included.
	 */
	symphas::problem_parameters_type get_problem_parameters() const
	{
		size_t sys_len = std::max(intervals_len, std::max(bdata_len, tdata_len));
		symphas::problem_parameters_type pp{ sys_len };

		pp.set_boundary_data(bdata, bdata_len);
		pp.set_initial_data(tdata, tdata_len);
		pp.set_interval_data(intervals, intervals_len);
		pp.set_time_step(dt);


		return pp;
	}




	void write(const char* savedir, const char* name = BACKUP_CONFIG_NAME) const;




protected:

	//! Extend the number of intervals.
	/*!
	 * Information from the original intervals is copied in the longer
	 * array.
	 * 
	 * \param new_intervals_len The length of the new intervals array.
	 */
	void extend_intervals(size_t new_intervals_len)
	{
		if (new_intervals_len > intervals_len)
		{
			auto extendv = new symphas::interval_data_type[new_intervals_len];
			auto extendm = new len_type*[new_intervals_len];
			for (iter_type i = 0; i < intervals_len; ++i)
			{
				extendv[i] = intervals[i];
				extendm[i] = new len_type[dimension];
				for (iter_type d = 0; d < dimension; ++d)
				{
					extendm[i][d] = dims[i][d];
				}
			}
			for (size_t i = intervals_len; i < new_intervals_len; ++i)
			{
				extendv[i] = symphas::interval_data_type{};
				extendm[i] = nullptr;

				for (iter_type d = 0; d < dimension; ++d)
				{
					extendv[i][symphas::index_to_axis(d)] = symphas::interval_element_type{};
				}
			}

			delete[] intervals;
			delete[] dims;

			intervals = extendv;
			dims = extendm;
			intervals_len = new_intervals_len;

		}
	}


	//! Set the dimensions using the nth interval set.
	/*!
		* Set the SystemConf::dims member to be an array of length prescribed by
		* the member SystemConf::dimension. Then fill the array by computing the
		* length (number of grid points) along each of the intervals from the data
		* given by the nth element from the interval list SystemConf::intervals.
		*
		* \param n The element index from SystemConf::intervals to choose the
		* interval data from.
		*/
	void set_dimensions(size_t n);

};



