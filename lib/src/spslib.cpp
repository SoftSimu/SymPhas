
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


#include "spslib.h"




len_type symphas::lib::radial_length(const len_type(&dim)[1])
{
	return dim[0];
}

len_type symphas::lib::radial_length(const len_type(&dim)[2])
{
	len_type len = std::min(dim[0], dim[1]);
	return (dim[0] * dim[1]) - ((len - 1) * len) / 2;
}

len_type symphas::lib::radial_length(const len_type(&dim)[3])
{
	len_type _dim[]{ dim[0], dim[1], dim[2] };
	std::sort(_dim, _dim + 3);

	len_type
		& l = _dim[2],
		& m = _dim[1],
		& n = _dim[0];
	return l * m * n
		- (m * (m - 1) * n / 2)
		- (l * (n - 1) * n / 2)
		+ ((n - 2) * (n - 1) * n / 2) / 3;
}




// ********************************************************************************

std::tuple<double*, size_t*> symphas::lib::make_radial_arrays(const len_type(&dims)[1],
	double dx)
{
	len_type len = radial_length(dims);
	double* rmap = new double[len];
	size_t* cmap = new size_t[len];

	double dx2 = dx * dx;

	iter_type n = 0;
	for (iter_type x = 0; x < dims[0]; ++x, ++n)
	{
		cmap[n] = 1;
		rmap[n] = sqrt(x * x * dx2);
	}
	symphas::lib::sort_data(rmap, cmap, len);
	return { rmap, cmap };
}

std::tuple<double*, size_t*> symphas::lib::make_radial_arrays(const len_type(&dims)[2],
	double dx, double dy)
{
	len_type rlen = radial_length(dims);
	double* rmap = new double[rlen];
	size_t* cmap = new size_t[rlen];

	double dx2 = dx * dx;
	double dy2 = dy * dy;
	/*
	 * the assumption is the list is sorted
	 */
	len_type
		len = std::max(dims[0], dims[1]),
		wid = std::min(dims[0], dims[1]);

	iter_type n = 0;
	for (iter_type y = 0; y < wid; ++y)
	{
		cmap[n] = 1;
		rmap[n] = sqrt(2 * y * y * dy2);

		++n;
		for (iter_type x = y + 1; x < len; ++x, ++n)
		{
			cmap[n] = (x < wid) ? 2 : 1;
			rmap[n] = sqrt(x * x * dx2 + y * y * dy2);
		}
	}
	symphas::lib::sort_data(rmap, cmap, rlen);
	return { rmap, cmap };
}

std::tuple<double*, size_t*> symphas::lib::make_radial_arrays(const len_type(&dims)[3],
	double dx, double dy, double dz)
{
	len_type rlen = radial_length(dims);
	double* rmap = new double[rlen];
	size_t* cmap = new size_t[rlen];

	double dx2 = dx * dx;
	double dy2 = dy * dy;
	double dz2 = dz * dz;

	len_type _dim[]{ dims[0], dims[1], dims[2] };
	std::sort(_dim, _dim + 3);
	size_t
		len = _dim[2],
		wid = _dim[1],
		dep = _dim[0];


	for (iter_type z = 0, n = 0; z < dep; ++z)
	{
		cmap[n] = 1;
		rmap[n] = sqrt(3 * z * z * dz2);

		++n;
		for (iter_type x = z + 1; x < len; ++x, ++n)
		{
			cmap[n] = (x < dep) ? 3 : (x < wid) ? 2 : 1;
			rmap[n] = sqrt(x * x * dx2 + 2 * z * z * dz2);
		}

		for (iter_type y = z + 1; y < wid; ++y)
		{
			cmap[n] = (y < dep) ? 3 : 2;
			rmap[n] = sqrt(2 * y * y * dy2 + z * z * dz2);

			++n;
			for (iter_type x = y + 1; x < len; ++x, ++n)
			{
				cmap[n] = (x < dep) ? 6 : (x < wid) ? 3 : 1;
				rmap[n] = sqrt(x * x * dx2 + y * y * dy2 + z * z * dz2);
			}
		}
	}
	symphas::lib::sort_data(rmap, cmap, rlen);
	return { rmap, cmap };
}



// ********************************************************************************



iter_type* symphas::lib::make_radial_index_map(const len_type(&dims)[1])
{
	len_type len = grid::length<1>(dims);
	iter_type* dmap = new iter_type[len];

	std::vector<std::tuple<double, size_t>> vecs;
	vecs.reserve(radial_length(dims));

	for (size_t i = 0; i < dims[0]; ++i)
	{
		vecs.emplace_back(sqrt(i * i), i);
	}
	std::sort(vecs.begin(), vecs.end(), [](auto a, auto b) { return std::get<0>(a) < std::get<0>(b); });

	iter_type n = 0;
	for (auto const& [_, i] : vecs)
	{
		dmap[i] = n++;
	}

	return dmap;
}

iter_type* symphas::lib::make_radial_index_map(const len_type(&dims)[2])
{
	len_type rlen = grid::length<2>(dims);
	iter_type* dmap = new iter_type[rlen]{ 0 };

	size_t
		len = std::max(dims[0], dims[1]),
		wid = std::min(dims[0], dims[1]);

	std::vector<std::tuple<double, size_t, size_t>> vecs;
	vecs.reserve(radial_length(dims));

	for (size_t j = 0, e = wid; j < e; ++j)
	{
		for (size_t i = j; i < len; ++i)
		{
			vecs.emplace_back(sqrt(i * i + j * j), i, j);
		}
	}
	std::sort(vecs.begin(), vecs.end(), [](auto a, auto b) { return std::get<0>(a) < std::get<0>(b); });

	iter_type n = 0;
	for (auto const& [_, i, j] : vecs)
	{
		for (auto ind : { i + j * dims[0], j + i * dims[0] })
		{
			if (ind < rlen && dmap[ind] <= 0)
			{
				dmap[ind] = n;
			}	
		}
		++n;
	}

	return dmap;
}


iter_type* symphas::lib::make_radial_index_map(const len_type(&dims)[3])
{
	len_type rlen = grid::length<3>(dims);
	iter_type* dmap = new iter_type[rlen]{ 0 };

	size_t
		len = std::max(dims[0], std::max(dims[1], dims[2])),
		wid = std::min(dims[0], std::max(dims[1], dims[2])),
		dep = std::min(dims[0], std::min(dims[1], dims[2]));

	std::vector<std::tuple<double, size_t, size_t, size_t>> vecs;
	vecs.reserve(radial_length(dims));

	for (size_t k = 0, n = 0; k < dep; ++k)
	{
		for (size_t j = k; j < wid; ++j)
		{
			for (size_t i = j; i < len; ++i, ++n)
			{
				vecs.emplace_back(sqrt(i * i + j * j + k * k), i, j, k);
			}
		}
	}
	std::sort(vecs.begin(), vecs.end(), [](auto a, auto b) { return std::get<0>(a) < std::get<0>(b); });

	iter_type n = 0;
	for (auto const& [_, i, j, k] : vecs)
	{
		for (auto ind : { 
			i + j * dims[0] + k * dims[0] * dims[1],
			i + k * dims[0] + j * dims[0] * dims[1],
			j + i * dims[0] + k * dims[0] * dims[1],
			j + k * dims[0] + i * dims[0] * dims[1],
			k + i * dims[0] + j * dims[0] * dims[1],
			k + j * dims[0] + i * dims[0] * dims[1]
			})
		{
			if (ind < rlen && dmap[ind] <= 0)
			{
				dmap[ind] = n;
			}
		}
		++n;
	}

	return dmap;
}






void symphas::lib::to_file_name(char const *in, char *out, size_t count)
{
	char const *inc = in;
	char *outc = out;

	while (*inc)
	{
		if (*inc == ' ')
		{
			*outc++ = '-';
		}
		else if (
			(*inc >= 'A' && *inc <= 'Z') ||
			(*inc >= '0' && *inc <= '9'))
		{
			*outc++ = static_cast<char>(tolower(*inc));
		}
		else if (
			(*inc >= 'a' && *inc <= 'z'))
		{
			*outc++ = *inc;
		}
		else if (*inc == '_' || *inc == '.' || *inc == ',')
		{
			*outc++ = '_';
		}
		++inc;

		if (count > 0 && static_cast<size_t>(outc - out) == count)
		{
			break;
		}
	}
	*outc = '\0';
}

// in place modification of the char array into lower case
void symphas::lib::to_lower(char *str)
{
	for (char *c = str; *c; ++c)
	{
		*c = static_cast<char>(tolower(*c));
	}
}

// doesn't modify the string, copies the lower case into another
void symphas::lib::to_lower(const char *str, char *out)
{
	const char *c = str;
	for (char *o = out; *c; ++c, ++o)
	{
		*o = static_cast<char>(tolower(*c));
	}
}

// in place modification of the char array into caps
void symphas::lib::to_upper(char *str)
{
	for (char *c = str; *c; ++c)
	{
		*c = static_cast<char>(toupper(*c));
	}
}

// doesn't modify the string, copies the upper case into another
void symphas::lib::to_upper(const char *str, char *out)
{
	const char *c = str;
	char *it;
	for (it = out; *c; ++c, ++it)
	{
		*it = static_cast<char>(toupper(*c));
	}
	*it = '\0';
}

// in place trimming of the string
void symphas::lib::str_trim(char *str)
{
	size_t len = std::strlen(str);
	char ws[] = { '\n', ' ', '\t', '\r' }; // whitespace chars
	char *nw = str, *mw = str + len, *it;
	bool r; // flag for empty region

	r = true;
	for (it = str; r; ++it)
	{
		r = false;
		if (it < str + len)
		{
			for (auto w : ws)
			{
				if (*it == w)
				{
					r = true;
				}
			}
		}
	}
	nw = it - 1;

	r = true;
	for (it = str + len; r; --it)
	{
		r = false;
		if (it > str)
		{
			for (auto w : ws)
			{
				if (*(it - 1) == w)
				{
					r = true;
				}
			}
		}
	}
	mw = it + 1;

	for (it = nw; it < mw; ++it, ++str)
	{
		*str = *it;
	}
	*str = '\0';
}

// in place trimming of the string
void symphas::lib::str_trim(const char* str, char* out)
{
	std::strcpy(out, str);
	str_trim(out);
}




iter_type symphas::lib::pos_after_digit(const char* content)
{
	const char* it = content;
	while (isdigit(it[0]))
	{
		it += 1;
	}
	return static_cast<iter_type>(it - content);
}

iter_type symphas::lib::pos_after_token(const char* content, const char* token)
{
	iter_type n = static_cast<iter_type>(std::strlen(token));
	for (iter_type i = 0; i < n; ++i)
	{
		if (content[i] != token[i])
		{
			return 0;
		}
	}
	return n;
}

iter_type symphas::lib::pos_after_token(const char* content, std::initializer_list<const char*> tokens)
{
	for (auto token : tokens)
	{
		iter_type n = pos_after_token(content, token);
		if (n != 0)
		{
			return n;
		}
	}
	return 0;
}




size_t symphas::lib::num_digits(size_t n)
{
	return num_digits(static_cast<int>(n));
}

size_t symphas::lib::num_digits(int n)
{
	if (((n > 0) ? n : -n) <= 9)
	{
		return 1;
	}

	return 1 + num_digits(n / 10);
}







void symphas::lib::write_ts_str(char* buffer)
{
	std::time_t t = std::time(0);
	std::tm* now = std::localtime(&t);

	sprintf(buffer, TIMESTAMP_FMT,
		now->tm_year + 1900,
		now->tm_mon + 1,
		now->tm_mday,
		now->tm_hour,
		now->tm_min,
		now->tm_sec);
}



#ifdef FILESYSTEM_HEADER_AVAILABLE
void symphas::lib::make_directory(std::filesystem::path dir, int err_no)
{
	// create parent directory for the data
	if (!std::filesystem::exists(dir))
	{
		bool created = std::filesystem::create_directories(dir);
		if (!created)
		{
			if (!std::filesystem::exists(dir))
			{
				fprintf(SYMPHAS_ERR, "%d error creating directory '%s'\n",
					err_no + 1, dir.string().c_str());
			}
		}
	}
}
#endif


void symphas::lib::make_directory(const char* dir, int err_no)
{
#ifdef _MSC_VER
	char copy_buffer[BUFFER_LENGTH]{}, build_tree[BUFFER_LENGTH]{};
	std::strncpy(copy_buffer, dir, sizeof(copy_buffer) / sizeof(char));

	char* tok = std::strtok(copy_buffer, "/\\");
	do 
	{
		std::strcat(build_tree, "/");
		std::strcat(build_tree, tok);
		if (!CreateDirectory(build_tree, NULL) && (ERROR_ALREADY_EXISTS != GetLastError()))
		{
			DWORD dwAttrib = GetFileAttributes(build_tree);
			if (dwAttrib != INVALID_FILE_ATTRIBUTES &&
				(dwAttrib & FILE_ATTRIBUTE_DIRECTORY))
			{
				fprintf(SYMPHAS_ERR, "%d error creating directory '%s'\n", err_no + 2, build_tree);
			}

		}
	} 
	while ((tok = std::strtok(NULL, "/\\")) != 0);

#else

	mode_t mode = 0755;
	struct stat sb;
	if (stat(dir, &sb) != 0)
	{
		if (mkdir(dir, mode) != 0 && errno != EEXIST)
		{
			fprintf(SYMPHAS_ERR, "%d error creating directory '%s'\n", err_no + 3, dir);
		}
	}
#endif
}


void symphas::lib::make_directory_for_file(const char* name, int err_no)
{
	char* parent = new char[std::strlen(name) + 1];
	std::strcpy(parent, name);
	get_parent_directory(parent, parent);

#ifdef FILESYSTEM_HEADER_AVAILABLE
	symphas::lib::make_directory(std::filesystem::path(parent), err_no);
#else
	symphas::lib::make_directory(copy, err_no);
#endif

	delete[] parent;
}

void symphas::lib::get_parent_directory(char* path, char* &basepath)
{
#ifdef _MSC_VER

	char drive[_MAX_DRIVE];
	char parent[_MAX_DIR];
	_splitpath(
		path,
		drive,
		parent,
		NULL,
		NULL
	);
	sprintf(basepath, "%s%s", drive, parent);

#else
	basepath = dirname(path);
#endif
}

void symphas::lib::get_parent_directory(const char* path, const char*& basepath)
{
#ifdef _MSC_VER

	char drive[_MAX_DRIVE];
	char parent[_MAX_DIR];
	_splitpath(
		path,
		drive,
		parent,
		NULL,
		NULL
	);
	sprintf(basepath, "%s%s", drive, parent);

#else
	basepath = dirname(const_cast<char*>(path));
#endif
}

#ifdef FILESYSTEM_HEADER_AVAILABLE
std::filesystem::path symphas::lib::get_parent_directory(std::filesystem::path dir)
{
	return dir.parent_path();
}
#endif

//! Applies binary sorting to axis values centered at the origin.
inline bool symphas::lib::origin_sort(axis_coord_t a, axis_coord_t b, len_type axis_len)
{
	if (a < 0)
	{
		a += axis_len;
	}
	if (b < 0)
	{
		b += axis_len;
	}

	return a < b;
}

inline bool symphas::lib::regular_sort(axis_coord_t a, axis_coord_t b, len_type)
{
	return a < b;
}

