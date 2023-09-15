
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


#include "io.h"


template<>
struct params::param_assign<symphas::IOType> : param_assign_base
{
	void assign(void* param, const char* value)
	{
		*static_cast<symphas::IOType*>(param) = extract_writer(value);
	}

	size_t print_with_name(FILE* out, void* param, const char* name)
	{
		using symphas::IOType;
		return fprintf(out, "%s=[gnu,xdr,column,movie,csv](default=%s)", name,
			(*static_cast<IOType*>(param) == IOType::GNU) ? "gnu" :
			(*static_cast<IOType*>(param) == IOType::XDR) ? "xdr" :
			(*static_cast<IOType*>(param) == IOType::COLUMN) ? "column" :
			(*static_cast<IOType*>(param) == IOType::MOVIE) ? "movie" :
			(*static_cast<IOType*>(param) == IOType::CSV) ? "csv" : "");
	}

protected:

	symphas::IOType extract_writer(const char* value)
	{
		using symphas::IOType;
		char* writer_name = new char[std::strlen(value) + 1];
		symphas::lib::str_trim(writer_name, writer_name);
		symphas::lib::to_upper(value, writer_name);
		
		// The default writer is the GNU matrix writer.
		IOType writer = IOType::GNU;

		if (std::strcmp(writer_name, "GNU") == 0)
		{
			writer = IOType::GNU;
		}
		else if (std::strcmp(writer_name, "XDR") == 0)
		{
			writer = IOType::XDR;
		}
		else if (std::strcmp(writer_name, "COLUMN") == 0)
		{
			writer = IOType::COLUMN;
		}
		else if (std::strcmp(writer_name, "MOVIE") == 0)
		{
			writer = IOType::MOVIE;
		}
		else if (std::strcmp(writer_name, "CSV") == 0)
		{
			writer = IOType::CSV;
		}
		else
		{
			fprintf(SYMPHAS_WARN, "the given writer/reader type '%s' "
				"is not recognized\n", writer_name);
		}


		delete[] writer_name;
		return writer;
	}
};


DLLIO bool params::plots_only = false;
DLLIO symphas::IOType params::writer = symphas::IOType::GNU;
DLLIO symphas::IOType params::reader = symphas::IOType::GNU;
DLLIO bool params::use_timestamp = true;
DLLIO bool params::single_output_file = true;
DLLIO bool params::single_input_file = true;
DLLIO char* params::checkpoint = nullptr;
DLLIO int params::checkpoint_count = 10;
DLLIO bool* params::single_io_file[2] { &params::single_input_file, &params::single_output_file };
DLLIO symphas::IOType* params::io_type[2] { &params::reader, &params::writer };

bool add_io_params(param_map_type& param_map)
{
	using namespace params;

	param_map["plots-only"] = { &plots_only, new param_assign<bool>, 'p', 
		"only outputs plot configuration files, without running the model" };
	param_map["writer-type"] = { &writer, new param_assign<symphas::IOType>, 'w',
		"chooses the writer type from 'gnu', 'xdr', 'column', 'movie', or 'csv'; use '-x' to set both reader and writer" };
	param_map["writer"] = { &writer, new param_assign<symphas::IOType> };
	param_map["reader-type"] = { &reader, new param_assign<symphas::IOType>, 'r',
		"chooses the writer type from 'gnu', 'xdr', 'column', 'movie', or 'csv'; use '-x' to set both reader and writer" };
	param_map["reader"] = { &reader, new param_assign<symphas::IOType> };
	param_map["use_timestamp"] = { &use_timestamp, new param_assign<bool>, 't', 
		"uses the timestamp when generating the output directory" };
	param_map["single-output_file"] = { &single_output_file, new param_assign<bool>, 'o', 
		"all data is appended to a single output file; use '-a' to set both single input and output" };
	param_map["single-input_file"] = { &single_input_file, new param_assign<bool>, 'i', 
		"all data is read from a single input file; use '-a' to set both single input and output" };
	param_map["single-output"] = { &single_output_file, new param_assign<bool> };
	param_map["single-input"] = { &single_input_file, new param_assign<bool> };
	param_map["checkpoint"] = { &checkpoint, new param_assign<char*>, 'C', 
		"load a checkpoint from this directory with a previous run, additionally appending ',N' to indicate the iteration index with N" };
	param_map["checkpoint-count"] = { &checkpoint_count, new param_assign<int>, 'c',
		"choose the number of checkpoints to save" };

	param_map["single-io-file"] = {
		params::single_io_file,
		new param_assign_multiple<bool, 2>, 'a' };

	param_map["io-type"] = {
		params::io_type,
		new param_assign_multiple<symphas::IOType, 2>, 'x' };


	return true;
}


SaveType symphas::io::get_save_type(const char* str)
{
	if (!*str)
	{
		return SaveType::DEFAULT;
	}
	else if (std::strcmp(str, "EXP") == 0)
	{
		return SaveType::EXP;
	}
	else if (std::strcmp(str, "MUL") == 0)
	{
		return SaveType::MULTIPLICATIVE;
	}
	else if (std::strcmp(str, "SPLIT") == 0)
	{
		return SaveType::SPLIT;
	}
	else
	{
		return SaveType::NONE;
	}
}

const char* symphas::io::get_save_str(SaveType type)
{
	switch (type)
	{
	case SaveType::DEFAULT:
		return "";
	case SaveType::MULTIPLICATIVE:
		return "MUL";
	case SaveType::EXP:
		return "EXP";
	case SaveType::SPLIT:
		return "SPLIT";
	default:
		return "NONE";
	}
}

void symphas::io::parse_save_str(const char* value, SaveParams* save)
{

	char const msg[] = "incorrect save configuration provided '%s'\n";
	char tag[BUFFER_LENGTH_R4]{};

	char* next0;
	double base = strtod(value, &next0);

	if (next0 > value)
	{
		char* next1;
		int start = strtol(next0, &next1, 10);

		if (*next0 == ',')
		{
			// this is a list spec

			size_t nc = std::count_if(next0, next0 + std::strlen(next0), [](char c) { return (c == ','); });
			save->set_indices(1 + nc);

			iter_type n = 1;
			save->indices[0] = static_cast<iter_type>(base);
			char* tok = std::strtok(next0, ",");
			do
			{
				save->indices[n++] = atoi(tok);
			} while ((tok = std::strtok(NULL, ",")) != NULL);
		}
		else
		{
			if (next1 > next0)
			{
				start = std::max(start, params::start_index);
			}
			else
			{
				start = params::start_index;
			}

			if (*next1)
			{
				// read a tag.
				size_t n = sscanf(next1, "%s", tag);

				if (n == 0)
				{
					save->type = SaveType::DEFAULT;
				}
				else
				{
					save->type = get_save_type(tag);
					if (save->type == SaveType::NONE)
					{
						fprintf(SYMPHAS_ERR, msg, value);
						exit(822);
					}
				}
			}

			save->base = (base <= 0) ? save->stop : base;
		}
	}

}

void symphas::io::save_as_str(SaveParams const* save, char* output, size_t buffer_size)
{
	if (save->type == SaveType::LIST)
	{
		size_t pos = 0;
		for (iter_type i = 0; i < save->count; ++i)
		{
			pos += snprintf(output + pos, buffer_size - pos, "%d", save->indices[i]);
			if (i < save->count - 1)
			{
				pos += snprintf(output + pos, buffer_size - pos, ",");
			}
		}
	}
	else
	{
		snprintf(output, buffer_size, "%lf %d %s",
			save->base, save->get_start(), get_save_str(save->type));
	}
}









std::vector<iter_type> symphas::io::get_save_indices(SaveParams const& main_save, std::vector<SaveParams> const& other_saves)
{
	std::vector<iter_type> saves;

	/* the checkpoint interval, computed by dividing the number of indices to solve
	 * by the number of checkpoints to save, then each one will be put in
	 */
	iter_type ci = (params::checkpoint_count > 0)
		? (main_save.get_stop() / params::checkpoint_count)
		: (main_save.get_stop() + 1);

	/* keep building the list until the index reaches the end
	 */
	iter_type index = params::start_index - 1;
	while (index < main_save.get_stop())
	{
		iter_type next_checkpoint = (ci > 0) ? (index / ci + 1) * ci : main_save.get_stop();
		iter_type min = std::min(
			next_checkpoint,
			main_save.next_save(index));

		// check the next save of each data parameter if they come first
		for (auto& save : other_saves)
		{
			iter_type n = save.next_save(index);
			if (n > 0)
			{
				min = std::min(n, min);
			}
		}

		index = min;
		saves.emplace_back(index);
	}

	if (saves.empty())
	{
		saves.emplace_back(params::start_index);
	}
	else
	{
		saves.back() = main_save.get_stop();
	}

	return saves;
}


iter_type symphas::io::next_save_from_list(iter_type index, SaveParams const& main_save, std::vector<SaveParams> const& other_saves)
{
	return next_save_from_list(index, symphas::io::get_save_indices(main_save, other_saves));
}

iter_type symphas::io::next_save_from_list(iter_type index, std::vector<iter_type> const& saves)
{
	auto lower = std::lower_bound(saves.begin(), saves.end(), index + 1);
	if (lower != saves.end())
	{
		return *lower;
	}
	else
	{
		return index;
	}
}


bool symphas::io::is_last_save(iter_type index, std::vector<iter_type> const& saves)
{
	return (!saves.empty()) ? index >= saves.back() : true;
}

iter_type symphas::io::next_checkpoint(iter_type index, SaveParams const& save)
{
	int ci = (params::checkpoint_count > 0)
		? (save.get_stop() / params::checkpoint_count)
		: (save.get_stop() + 1);

	return (ci > 0) ? (index / ci) * ci : save.get_stop();
}

bool symphas::io::at_checkpoint(iter_type index, SaveParams const& save)
{
	return index == symphas::io::next_checkpoint(index, save) || save.current_save(index);
}

bool symphas::io::is_checkpoint_set()
{
	return params::checkpoint != nullptr;
}





void SaveParams::set_params(SaveType type, double base, iter_type start, iter_type stop)
{
	if (type == SaveType::LIST)
	{
		if (init_flag)
		{
			set_indices(static_cast<iter_type>(base) + 1);
			indices[0] = start;
		}
		else
		{
			set_indices(static_cast<iter_type>(base));
		}

		double interval = (stop - start) / base;
		double current = start + interval;

		for (iter_type i = (init_flag) ? 1 : 0; i < count; ++i)
		{
			indices[i] = std::lround(current);
			current += interval;
		}
	}
	else
	{
		this->type = type;
		this->start = start;
		this->stop = stop;
		this->base = base;
	}
}

void SaveParams::set_params(double base, iter_type start, iter_type stop)
{
	set_params(type, base, start, stop);
}

void SaveParams::set_start(iter_type start)
{
	if (type == SaveType::LIST)
	{
		if (indices[0] > start)
		{
			iter_type* indices_extend = new iter_type[count + 1];
			std::copy(indices, indices + count, indices_extend + 1);
			indices_extend[0] = start;
			set_indices(indices_extend, count + 1);
			delete[] indices_extend;
		}
		else
		{
			indices[0] = start;
		}
	}
	else
	{
		this->start = start;
	}
}

void SaveParams::set_stop(iter_type stop)
{
	if (type == SaveType::LIST)
	{
		if (indices[count - 1] < stop)
		{
			iter_type* indices_extend = new iter_type[count + 1];
			std::copy(indices, indices + count, indices_extend);
			indices_extend[count] = stop;
			set_indices(indices_extend, count + 1);
			delete[] indices_extend;
		}
		else
		{
			indices[count - 1] = stop;
		}
	}
	else
	{
		this->stop = stop;
	}
}


bool SaveParams::index_zero_save(iter_type index) const
{
	return (index == params::start_index && init_flag) || get_stop() == params::start_index;
}

bool SaveParams::first_save(iter_type index) const
{
	return index == get_start();
}

bool SaveParams::current_save(iter_type index) const
{
	return index_zero_save(index) || (index == next_save(index - 1) && index > params::start_index);
}

bool SaveParams::is_last_save(iter_type index) const
{
	return index >= get_stop();
}


/* returns the index of the next save, based on the save type
 */
iter_type SaveParams::next_save(iter_type current) const
{

	if (type == SaveType::LIST)
	{
		for (iter_type i = 0; i < count; ++i)
		{
			if (current < indices[i])
			{
				if (indices[i] >= params::start_index || init_flag)
				{
					return indices[i];
				}
			}
		}
		return stop;
	}
	else
	{
		if (stop == start)
		{
			return stop;
		}
		else if (current < start)
		{
			if (init_flag)
			{
				return start;
			}
			else if (start < params::start_index)
			{
				return next_save(params::start_index);
			}
			else
			{
				return next_save(start);
			}
		}
		else if (base <= 0)
		{
			return stop;
		}
		else
		{
			iter_type next = current;

			switch (type)
			{
			case SaveType::SPLIT:
			{
				double step = (stop - start) / base;
				double next0 = std::floor(step * std::ceil((current + 1) / step));
				next = std::min(
					std::max(
						static_cast<iter_type>(next0), start),
					stop);
				break;
			}

			case SaveType::DEFAULT:

				next = std::min(
					std::max(
						static_cast<iter_type>(base * (std::floor(current / base) + 1)), start),
					stop);
				break;

			case SaveType::MULTIPLICATIVE:
			{
				// Sets the size of the smallest gap, relative to an even
				// gap size.
				double initial_gap_scale = 1;

				double Dn = (stop - start);
				int N = (base > stop - start) ? (stop - start) : static_cast<iter_type>(base);

				double p = 0;
				for (iter_type i = 0; i < N; ++i)
				{
					p += 1.0 / std::sqrt(i + 1);
				}
				double r = std::max(1.0, initial_gap_scale * std::pow(Dn, 1.0 / p));

				auto f = [r, Dn, N](double alph)
				{
					double poly_alph = 0;
					for (iter_type i = 0; i < N; ++i)
					{
						poly_alph += std::pow(alph, i);
					}
					return r * poly_alph - Dn;
				};

				auto df = [r, N](double alph)
				{
					double dpoly_alph = 0;
					for (iter_type i = 0; i < N; ++i)
					{
						dpoly_alph += i * std::pow(alph, i - 1);
					}
					return r * dpoly_alph;
				};

				double eps = 1e-4;
				double diff = eps;
				double alph0 = 1.0;
				while (diff >= eps)
				{
					double alph1 = alph0 - f(alph0) / df(alph0);
					diff = std::abs(alph1 - alph0);
					alph0 = alph1;
				}

				if (!std::isfinite(alph0))
				{
					double step = Dn / N;
					next = std::min(
						std::max(
							static_cast<iter_type>(step * ((current / step) + 1)), start),
						stop);
				}
				else
				{
					double end = start;
					double step = r;
					if (step > eps)
					{
						while (std::lround(end) <= current)
						{
							end += step;
							step *= alph0;
						}
						next = std::lround(end);
					}
					else
					{
						next = current;
					}
				}
				break;

			}
			case SaveType::EXP:

				if (current < 1)
				{
					next = std::max(start, static_cast<iter_type>(base));
				}
				else
				{
					double power = std::ceil(std::log1p(current) / std::log(base));
					double next0 = std::pow(base, power);
					next = std::min(std::max(static_cast<iter_type>(next0), start), stop);
				}
				break;

			default:

				next = stop;
			}

			if (next >= stop)
			{
				return stop;
			}
			else if (next <= params::start_index && !init_flag)
			{
				return next_save(params::start_index);
			}
			else
			{
				return next;
			}
		}
	}
}

/*
 * returns the number of saves that will be performed
 */
iter_type SaveParams::num_saves() const
{

	if (type == SaveType::LIST)
	{
		if (indices[0] != params::start_index && init_flag)
		{
			return static_cast<iter_type>(count) + 1;
		}
		else
		{
			return static_cast<iter_type>(count);
		}
	}
	else
	{
		if (base == 0)
		{
			return 1;
		}
		else
		{
			switch (type)
			{
			case SaveType::DEFAULT:
			{
				return static_cast<iter_type>((stop - start) / base + ((init_flag) ? 1 : 0));
			}
			case SaveType::MULTIPLICATIVE:
			{
				return static_cast<iter_type>(base + ((init_flag) ? 1 : 0));
			}
			case SaveType::SPLIT:
			{
				return static_cast<iter_type>(base + ((init_flag) ? 1 : 0));
			}
			case SaveType::EXP:
			{
				return static_cast<iter_type>(
					std::log(static_cast<double>(stop - start)) / std::log(base))
					+ ((init_flag) ? 1 : 0);
			}
			default:
				return 0;
			}
		}
	}

}





