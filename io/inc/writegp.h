
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
 * MODULE:  io
 * PURPOSE: Defines the text output functionality for plain text files.
 *
 * ***************************************************************************
 */

#pragma once

#include "writedefines.h"

//! \cond

#define MULTI_MARGIN_RATIO_KEY MARGIN
#define MULTI_ALIGN_SEPARATION_KEY ALIGN_SEP
#define COLORBOX_HEIGHT_RATIO_KEY CB_HEIGHT
#define COLORBOX_BOTTOM_PADDING_KEY CB_MARGIN
#define MULTI_OUTPUT_WIDTH_KEY OUTPUT_WIDTH

#define ALIGNMENT_FILE_NAME "alignments.txt"
#define ALIGNMENT_PAR_SEP_CHAR "="

#define MULTI_OUTPUT_WIDTH_DEFAULT 7.0
#define MULTI_MARGIN_RATIO_DEFAULT 0.03
#define MULTI_ALIGN_SEPARATION_DEFAULT 0.012
#define COLORBOX_HEIGHT_RATIO_DEFAULT 0.05
#define COLORBOX_BOTTOM_PADDING_DEFAULT 0.02

/*
 * the customizable parameters
 */

#define WIDTH_LATEX 6
#define WIDTH_WIN 500

#define MULTI_WIDTH_WIN 1000
#define MULTI_ALIGN_KEY "Right"  //!< Alignment of key (title) of subplots.

#define MULTI_COLUMN_NUM(N)  \
  ((N > 12)                  \
       ? ((int)std::sqrt(N)) \
       : (((N) % 3 == 0) ? 3 \
                         : 2))  //!< Number of columns in latex plot output.

#define MULTI_OUTPUT_WIDTH symphas::io::gp::alignment::display_width
#define MULTI_MARGIN_RATIO symphas::io::gp::alignment::mmr
#define MULTI_ALIGN_SEPARATION symphas::io::gp::alignment::mas
#define COLORBOX_HEIGHT_RATIO symphas::io::gp::alignment::chr
#define COLORBOX_BOTTOM_PADDING symphas::io::gp::alignment::cbp

#define MULTI_ALIGN_WIDTH \
  (1 - 2 * MULTI_MARGIN_RATIO)  //!< Percentage of space occupied by the plots.

/*
 * programmed parameters
 */

//! Number of leftover squares in the last row.
#define COLUMNS_LAST_ROW(N) ((N) % MULTI_COLUMN_NUM(N))
#define ADD_ROW_NUM(N) \
  ((MULTI_COLUMN_NUM(N) - COLUMNS_LAST_ROW(N)) % MULTI_COLUMN_NUM(N))

//! Number of rows.
#define MULTI_ROW_NUM(N) (((N) + ADD_ROW_NUM(N)) / MULTI_COLUMN_NUM(N))

//! Computed width of the multiplot.
#define MULTI_WIDTH_LATEX(N) \
  (MULTI_OUTPUT_WIDTH *      \
   std::pow(0.85,            \
            std::max(0.0, 0.0 + MULTI_ROW_NUM(N) - MULTI_COLUMN_NUM(N))))

#define MULTI_WIDTH_RATIO MULTI_ALIGN_WIDTH
#define MULTI_WIDTH_PLOT(N) \
  (MULTI_WIDTH_RATIO / MULTI_COLUMN_NUM(N) - 2 * MULTI_ALIGN_SEPARATION)
#define MULTI_HEIGHT_RATIO(N)                                              \
  ((MULTI_WIDTH_PLOT(N) + 2 * MULTI_ALIGN_SEPARATION) * MULTI_ROW_NUM(N) + \
   COLORBOX_HEIGHT_RATIO * MULTI_WIDTH_RATIO)
#define MULTI_HEIGHT_PLOT(N) MULTI_WIDTH_PLOT(N)

#define MULTI_HEIGHT_LATEX(N) (MULTI_WIDTH_LATEX(N) * MULTI_HEIGHT_RATIO(N))

constexpr char multi_align_fmt[] = R"~(
set tmargin at screen %.3lf;
set bmargin at screen %.3lf;
set lmargin at screen %.3lf;
set rmargin at screen %.3lf;)~";
#define MULTI_ALIGN_MARGIN_FMT multi_align_fmt

#define MULTI_ALIGN_MARGIN_VALUES(x, y, N)                                 \
  1 - (((MULTI_HEIGHT_PLOT(N) + 2 * MULTI_ALIGN_SEPARATION) /              \
        MULTI_HEIGHT_RATIO(N)) *                                           \
           (y) +                                                           \
       MULTI_ALIGN_SEPARATION / MULTI_HEIGHT_RATIO(N)),                    \
      1 - (((MULTI_HEIGHT_PLOT(N) + 2 * MULTI_ALIGN_SEPARATION) /          \
            MULTI_HEIGHT_RATIO(N)) *                                       \
               (y + 1) -                                                   \
           MULTI_ALIGN_SEPARATION / MULTI_HEIGHT_RATIO(N)),                \
      MULTI_MARGIN_RATIO +                                                 \
          ((MULTI_WIDTH_PLOT(N) + 2 * MULTI_ALIGN_SEPARATION)) * (x) +     \
          MULTI_ALIGN_SEPARATION,                                          \
      MULTI_MARGIN_RATIO +                                                 \
          ((MULTI_WIDTH_PLOT(N) + 2 * MULTI_ALIGN_SEPARATION)) * (x + 1) - \
          MULTI_ALIGN_SEPARATION

#define MULTI_COLORBOX_FMT \
  "set colorbox horizontal user origin %.3lf,%.3lf size %.3lf,%.3lf\n"

#define MULTI_COLORBOX_VALUES(N)                          \
  (MULTI_MARGIN_RATIO + MULTI_ALIGN_SEPARATION),          \
      (COLORBOX_BOTTOM_PADDING / MULTI_HEIGHT_RATIO(N)),  \
      (MULTI_WIDTH_RATIO - 2 * MULTI_ALIGN_SEPARATION),   \
      (COLORBOX_HEIGHT_RATIO - COLORBOX_BOTTOM_PADDING) / \
          MULTI_HEIGHT_RATIO(N)

#define SINGLE_WIDTH_LATEX MULTI_WIDTH_LATEX(1)
#define SINGLE_ALIGN_HEIGHT_PR \
  (0.8 * (MULTI_HEIGHT_PLOT(1) / MULTI_HEIGHT_RATIO(1)))

#define SINGLE_ALIGN_V_MARGIN                                                  \
  (((MULTI_HEIGHT_PLOT(1) / MULTI_HEIGHT_RATIO(1)) - SINGLE_ALIGN_HEIGHT_PR) / \
   2)

#define SINGLE_ALIGN_MARGIN_VALUES                           \
  1 - (MULTI_ALIGN_SEPARATION / MULTI_HEIGHT_RATIO(1)) -     \
      SINGLE_ALIGN_V_MARGIN,                                 \
      1 -                                                    \
          (SINGLE_ALIGN_HEIGHT_PR +                          \
           MULTI_ALIGN_SEPARATION / MULTI_HEIGHT_RATIO(1)) - \
          SINGLE_ALIGN_V_MARGIN,                             \
      (1 - SINGLE_ALIGN_HEIGHT_PR) / 2, (1 + SINGLE_ALIGN_HEIGHT_PR) / 2

#define SINGLE_COLORBOX_VALUES                               \
  (1 - SINGLE_ALIGN_HEIGHT_PR) / 2,                          \
      SINGLE_ALIGN_V_MARGIN / 2 +                            \
          (COLORBOX_BOTTOM_PADDING / MULTI_HEIGHT_RATIO(1)), \
      (SINGLE_ALIGN_HEIGHT_PR),                              \
      (COLORBOX_HEIGHT_RATIO - COLORBOX_BOTTOM_PADDING) /    \
          MULTI_HEIGHT_RATIO(1) / 2

#define PLOT_SQUARE_DIMENSION_LATEX "%.2lf,%.2lf"
#define PLOT_MULTI_DIMENSION_LATEX "%.2lf,%.2lf"

#define PLOT_SQUARE_DIMENSION_WIN STR(WIDTH_WIN) ",450"
#define PLOT_MULTI_DIMENSION_WIN STR(MULTI_WIDTH_WIN) ",%lld"

#define MULTI_PLOT_KEY_VPOS (1.07 + 0.01 * (7 - MULTI_OUTPUT_WIDTH))

#define LATEX_KEY_FMT "%s at index $%d$"
#define LATEX_SUBPLOT_KEY_FMT "index $%d$"

// ****************************************************************************

//! \endcond

namespace symphas::io {
//! Defines elements used in input and output for the text format.
/*!
 * The format of input/output is based on using it with the gnuplot
 * utility.
 */
namespace gp {}
}  // namespace symphas::io

namespace symphas::io::gp {

namespace alignment {

DLLIO extern double
    display_width;        //!< Base width value, in inches, of the output.
DLLIO extern double mmr;  //!< The percentage of space which is the buffer on
                          //!< each side, relative to width.
DLLIO extern double
    mas;  //!< Percentage of space between each subplot, relative to the width.
DLLIO extern double
    chr;  //!< Height of the colorbox as a percentage of the width.
DLLIO extern double cbp;  //!< Amount of space given to the colorbox legend, as
                          //!< a percentage of the width.

DLLIO extern std::pair<const char*, double*> key_par_pairs[5];

inline void load_parameters() {
  FILE* f;
  if ((f = fopen(ALIGNMENT_FILE_NAME, "r")) != 0) {
    fprintf(SYMPHAS_INFO,
            "using alignments file for latex plot "
            "configuration, %s\n",
            ALIGNMENT_FILE_NAME);

    char buffer[LINE_READ_BUFFER];
    while (fgets(buffer, LINE_READ_BUFFER, f) != 0) {
      char* key = strtok(buffer, ALIGNMENT_PAR_SEP_CHAR);
      bool key_matches = false;
      for (auto&& [match, parameter] : key_par_pairs) {
        if (strcmp(key, match) == 0) {
          key = strtok(NULL, ALIGNMENT_PAR_SEP_CHAR);
          *parameter = atof(key);
          key_matches = true;
        }
      }

      if (!key_matches) {
        fprintf(SYMPHAS_INFO,
                "the key '%s' is not a valid "
                "alignment option\n",
                key);
      }
    }
    fclose(f);
  }
}

}  // namespace alignment

constexpr char ranges_axis_2d[] = R"~(
set xrange [)~" COORD_FMT ":" COORD_FMT R"~(]
set yrange [)~" COORD_FMT ":" COORD_FMT R"~(]
)~";

constexpr char ranges_xy[] = R"~(
set xrange [)~" COORD_FMT ":" COORD_FMT R"~(]
set yrange [%.2lf:%.2lf]
)~";

constexpr char ranges_y[] = R"~(
set xrange [*:*]
set yrange [%.2lf:%.2lf]
)~";

constexpr char ranges_x[] = R"~(
set xrange [)~" COORD_FMT ":" COORD_FMT R"~(]
set yrange [*:*]
)~";

constexpr char ranges_auto[] = R"~(
set xrange [*:*]
set yrange [*:*]
)~";

constexpr char gnu_1d[] = R"~(
plot "%s" index %d using 1:2 with lines title "%s")~";

constexpr char gnu_s[] = R"~(
plot "%s" index %d matrix nonuniform with image title "%s")~";

constexpr char gnu_2v[] = R"~(
plot "%s" index %d using 1:2:3:4:5 with vectors head size 1,20,60 filled lc palette title "%s")~";

constexpr char gnu_3v[] = R"~(
plot "%s" index %d using 1:2:3:4:5:6:7 with vectors head size 1,20,60 filled lc palette title "%s")~";

constexpr char gnuset_1[] = R"~(

reset

set x2tics scale 0
set y2tics scale 0

unset x2tics
unset key
unset title


set xrange [%f:%f]
set yrange [*:*]
set xtics %f,%f,%f
set xtics axis out
set xlabel "x"
set ylabel "Order Parameter")~";

constexpr char gnuset[] = R"~(

reset

set x2tics scale 0
set y2tics scale 0

unset x2tics
unset y2tics
unset key
unset title


set xrange [%f:%f]
set yrange [%f:%f]
set xtics %f,%f,%f
set ytics %f,%f,%f
)~"

#ifdef LATEX_PLOT
                          R"~(
set cbrange [-1:1]
)~"
#else
                          R"~(
set cbrange [*:*]
set xtics axis out
set ytics axis out
set xlabel "x"
set ylabel "y"
set cblabel "Order Parameter")~"
#endif

                          R"~(

)~";

constexpr char gnu_single[] = R"~(
set title "%s"
unset key
set xtics auto
set ytics auto

)~"
#ifdef LATEX_PLOT
                              R"~(
set terminal epslatex color size )~" PLOT_SQUARE_DIMENSION_LATEX R"~(
set output "%s"
)~"
#else
                              R"~(
set terminal )~" STR(GNU_PLOTTER) R"~( size )~" PLOT_SQUARE_DIMENSION_WIN R"~(
)~"
#endif
                              R"~(
)~";

constexpr char gnu_multi[] = R"~(
unset multiplot
unset xtics
unset ytics
unset xlabel
unset ylabel
)~"
#ifdef LATEX_PLOT
                             R"~(
unset colorbox
set key at graph 1,%.3lf )~" MULTI_ALIGN_KEY R"~( reverse samplen -1
set terminal epslatex color size )~" PLOT_MULTI_DIMENSION_LATEX R"~(
set output "%s"
set multiplot layout %zd,%d
)~"
#else
                             R"~(
set colorbox
set key outside top center
set terminal )~" STR(GNU_PLOTTER) R"~( size )~" PLOT_MULTI_DIMENSION_WIN R"~(
set multiplot layout %zd,2 title "%s"
)~"
#endif
                             R"~(
)~";

constexpr char gnu_end[] = R"~(
)~"
#ifdef LATEX_PLOT
                           R"~(
unset output
)~"
#else
#endif
                           R"~(
)~";

template <size_t D>
void print_plot_ranges(FILE* gnu,
                       symphas::interval_data_type const& intervals) {
  fprintf(gnu, gnuset, DOMAIN_X0, DOMAIN_Xn, DOMAIN_Y0, DOMAIN_Yn, DOMAIN_X0,
          DOMAIN_Xn - DOMAIN_X0, DOMAIN_Xn, DOMAIN_Y0, DOMAIN_Yn - DOMAIN_Y0,
          DOMAIN_Yn);
}

template <>
inline void print_plot_ranges<1>(FILE* gnu,
                                 symphas::interval_data_type const& intervals) {
  fprintf(gnu, gnuset_1, DOMAIN_X0, DOMAIN_Xn - INTERVAL_Xh, DOMAIN_X0,
          DOMAIN_Xn - DOMAIN_X0, DOMAIN_Xn - INTERVAL_Xh);
}

template <size_t D, typename S0>
struct plot_fmt {
  inline void operator()(const char*&) {}
};

template <size_t D, typename S0>
struct plot_fmt<D, symphas::internal::field_array_t<S0>> {
  inline void operator()(const char*& gnu_set) {
    plot_fmt<D, symphas::internal::non_parameterized_type<D, S0>>{}(gnu_set);
  }
};

template <>
struct plot_fmt<1, scalar_t> {
  inline void operator()(const char*& gnu_set) {
    gnu_set = symphas::io::gp::gnu_1d;
  }
};

template <>
struct plot_fmt<2, scalar_t> {
  inline void operator()(const char*& gnu_set) {
    gnu_set = symphas::io::gp::gnu_s;
  }
};

template <>
struct plot_fmt<3, scalar_t> {
  inline void operator()(const char*& gnu_set) {
    gnu_set = symphas::io::gp::gnu_s;
  }
};

template <>
struct plot_fmt<1, complex_t> {
  inline void operator()(const char*& gnu_set) {
    gnu_set = symphas::io::gp::gnu_1d;
  }
};

template <>
struct plot_fmt<2, complex_t> {
  inline void operator()(const char*& gnu_set) {
    gnu_set = symphas::io::gp::gnu_s;
  }
};

template <>
struct plot_fmt<3, complex_t> {
  inline void operator()(const char*& gnu_set) {
    gnu_set = symphas::io::gp::gnu_s;
  }
};

template <>
struct plot_fmt<2, vector_t<2>> {
  inline void operator()(const char*& gnu_set) {
    gnu_set = symphas::io::gp::gnu_2v;
  }
};

template <>
struct plot_fmt<3, vector_t<3>> {
  inline void operator()(const char*& gnu_set) {
    gnu_set = symphas::io::gp::gnu_3v;
  }
};

//! Get the string corresponding to the gnuplot input format for the data.
template <size_t D, typename... S, size_t... Is>
void get_plot_fmt(const char*(*gnu_set), std::index_sequence<Is...>) {
  ((..., plot_fmt<D, S>{}(gnu_set[Is])));
}
//! Get the string corresponding to the gnuplot input format for the data.
template <size_t D, typename... S>
void get_plot_fmt(
    const char*(*gnu_set),
    symphas::lib::types_list<symphas::internal::field_array_t<void>, S...>) {
  get_plot_fmt<D, symphas::internal::non_parameterized_type<D, S>...>(
      gnu_set, std::make_index_sequence<sizeof...(S)>{});
}

//! Get the string corresponding to the gnuplot input format for the data.
template <size_t D, typename S0, typename... S>
void get_plot_fmt(const char*(*gnu_set), symphas::lib::types_list<S0, S...>) {
  get_plot_fmt<D, symphas::internal::non_parameterized_type<D, S0>,
               symphas::internal::non_parameterized_type<D, S>...>(
      gnu_set, std::make_index_sequence<sizeof...(S) + 1>{});
}

//! Get the string corresponding to the gnuplot input format for the data.
template <size_t D, typename... S>
void get_plot_fmt(const char*(*gnu_set)) {
  get_plot_fmt<D>(gnu_set, symphas::lib::types_list<S...>{});
}

//! Print the header, which details the information about the grid.
/*!
 * Print the header, which details the information about the grid.
 *
 * \param index The solution index of the data being saved.
 * \param id The id of the phase-field being saved.
 * \param ginfo Information about the grid.
 * \param f The file buffer to write to.
 */
void print_gp_header(int index, size_t id, symphas::grid_info const& ginfo,
                     symphas::io::write_info const& winfo, FILE* f);

template <size_t N, size_t D, typename Sp, typename S0, typename... S>
size_t get_system_offset(S0, Model<D, Sp, S...> const& model) {
  return 1;
}

template <size_t N, size_t D, typename Sp, typename S0, typename... S>
size_t get_system_offset(symphas::internal::field_array_t<S0>,
                         Model<D, Sp, S...> const& model) {
  return model.template num_fields<N>();
}

template <size_t D, typename Sp, typename... S, size_t... Is>
size_t get_system_id(size_t index, Model<D, Sp, S...> const& model,
                     std::index_sequence<Is...>) {
  size_t offsets[sizeof...(Is)];
  ((offsets[Is] = get_system_offset<Is>(S{}, model)), ...);

  size_t id = 0;
  for (size_t i = 0; i < index; ++i) {
    id += offsets[i];
  }
  return id;
}

template <size_t D, typename Sp, typename... S>
size_t get_system_id(size_t index, Model<D, Sp, S...> const& model) {
  return get_system_id(index, model, std::make_index_sequence<sizeof...(S)>{});
}

template <size_t N, size_t D, typename Sp, typename S0, typename... S>
size_t get_system_offset(S0, ArrayModel<D, Sp, S...> const& model) {
  return 1;
}

template <size_t N, size_t D, typename Sp, typename S0, typename... S>
size_t get_system_offset(symphas::internal::field_array_t<S0>,
                         ArrayModel<D, Sp, S...> const& model) {
  return model.template num_fields<N>();
}

template <size_t D, typename Sp, typename... S, size_t... Is>
size_t get_system_id(size_t index, ArrayModel<D, Sp, S...> const& model,
                     std::index_sequence<Is...>) {
  size_t offsets[sizeof...(Is)];
  ((offsets[Is] = get_system_offset<Is>(S{}, model)), ...);

  size_t id = 0;
  for (size_t i = 0; i < index; ++i) {
    id += offsets[i];
  }
  return id;
}

template <size_t D, typename Sp, typename... S>
size_t get_system_id(size_t index, ArrayModel<D, Sp, S...> const& model) {
  return get_system_id(index, model, std::make_index_sequence<sizeof...(S)>{});
}

/*
 * prints the plotting data file
 * takes as parameters the configuration and the number of systems that
 * will be plotted
 */

template <typename M>
constexpr size_t model_type_len = 0;

template <size_t D, typename Sp, typename... S>
constexpr size_t model_type_len<Model<D, Sp, S...>> = sizeof...(S);

template <size_t D, typename Sp, typename... S>
constexpr size_t model_type_len<ArrayModel<D, Sp, S...>> = sizeof...(S);

#ifdef LATEX_PLOT
template <size_t D, typename Sp, typename... S>
void write_plot_config(Model<D, Sp, S...> const& model, const char* directory,
                       const char* const* names, SaveParams const& save) {
  for (size_t id = 0; id < model_type_len<Model<D, Sp, S...>>; ++id) {
    FILE* gnu;
    char plot_name[BUFFER_LENGTH]{},     // name of the plot file
        data_name[BUFFER_LENGTH_R2]{},   // the data name of the plot file,
                                         // composed of GRID and id
        title_name[BUFFER_LENGTH_R2]{},  // the title of the simulation
                                         // formatted for filenames
        data_loc[BUFFER_LENGTH]{},       // the data file location
        latex_name[BUFFER_LENGTH]{},     // the name of the latex output file
        fid[BUFFER_LENGTH_R2]{};  // the printed numeric system identifier

    symphas::lib::to_file_name(params::title, title_name, BUFFER_LENGTH_R2);
    snprintf(latex_name, sizeof(latex_name) / sizeof(char),
             OUTPUT_LATEX_FILE_FMT, title_name, PHASEFIELD_DATA_NAME, id,
             save.get_stop());
    snprintf(fid, sizeof(fid) / sizeof(char), POSTFIX_ID_FMT, id);
    snprintf(data_name, sizeof(data_name) / sizeof(char), "%s%zd",
             PHASEFIELD_DATA_NAME, id);
    snprintf(plot_name, sizeof(plot_name) / sizeof(char),
             PHASEFIELD_PLOT_LOC_FMT, directory, data_name);

    if ((gnu = fopen(plot_name, "w")) == 0) {
      symphas::lib::make_directory_for_file(plot_name);
      if ((gnu = fopen(plot_name, "w")) == 0) {
        fprintf(SYMPHAS_ERR, SYMPHAS_MSG_ERR_FILE_OPEN, plot_name);
        exit(ERR_CODE_FILE_OPEN);
      }
    }

    const auto intervals = model.template system<0>().get_info().intervals;
    print_plot_ranges<D>(gnu, intervals);

    size_t saves = save.num_saves();
    const char* sets[model_type_len<Model<D, Sp, S...>>];
    get_plot_fmt<D, S...>(sets);

    alignment::load_parameters();
    if (saves > 1) {
      fprintf(gnu, gnu_multi, MULTI_PLOT_KEY_VPOS, MULTI_WIDTH_LATEX(saves),
              MULTI_HEIGHT_LATEX(saves), latex_name, MULTI_ROW_NUM(saves),
              MULTI_COLUMN_NUM(saves));

      for (iter_type index = save.get_start(), i = 0; i < saves;
           index = save.next_save(index), ++i) {
        char display_title[BUFFER_LENGTH];
        snprintf(display_title, sizeof(display_title) / sizeof(char),
                 LATEX_SUBPLOT_KEY_FMT, index);
        symphas::io::copy_data_file_name(DATA_DIR_RELATIVE_PLOT, index,
                                         get_system_id(id, model),
                                         DataFileType::SOLUTION_DATA, data_loc);

        /* if this is the final save, print additional plot configuration to
         * load the colorbar as well as correctly align the colorbar
         */
        if (i == saves - 1) {
          fprintf(gnu, MULTI_COLORBOX_FMT, MULTI_COLORBOX_VALUES(saves));
        }

        /* plot the grid in the multiplot with the correct width/height,
         * and provide the data file location and display title to the format
         */
        iter_type x = i % MULTI_COLUMN_NUM(saves),
                  y = i / MULTI_COLUMN_NUM(saves);
        fprintf(gnu, MULTI_ALIGN_MARGIN_FMT,
                MULTI_ALIGN_MARGIN_VALUES(x, y, saves));

        iter_type i0 = (params::single_output_file) ? std::max(0, i) : 0;

        fprintf(gnu, sets[id], data_loc, i0, display_title);
        fprintf(gnu, "\n\n");
      }
      fprintf(gnu, "\nunset multiplot");
    } else {
      char display_title[BUFFER_LENGTH];
      snprintf(display_title, sizeof(display_title) / sizeof(char),
               "%s---" LATEX_KEY_FMT, params::title, names[id],
               save.get_start());

      fprintf(gnu, gnu_single, display_title, SINGLE_WIDTH_LATEX,
              SINGLE_WIDTH_LATEX, latex_name);

      symphas::io::copy_data_file_name(
          DATA_DIR_RELATIVE_PLOT, save.next_save(save.get_start()),
          get_system_id(id, model), DataFileType::SOLUTION_DATA, data_loc);

      fprintf(gnu, MULTI_COLORBOX_FMT, SINGLE_COLORBOX_VALUES);

      iter_type index =
          (params::single_output_file) ? std::max(0, save.num_saves() - 1) : 0;

      fprintf(gnu, MULTI_ALIGN_MARGIN_FMT, SINGLE_ALIGN_MARGIN_VALUES);
      fprintf(gnu, sets[id], data_loc, index, display_title);
      fprintf(gnu, "\n");
    }

    fprintf(gnu, gnu_end);
    fclose(gnu);
  }
}

#else

//! Arrange the grid in either single or multiplot arrangement.
template <size_t N>
void print_plot_arrangement(char* gnuplot, size_t print_length) {
  snprintf(gnuplot, print_length, gnu_multi, 450 * (N + 1) / 2, (N + 1) / 2,
           params::title);
}

//! Specialization of print_plot_arrangement(char*, size_t).
template <>
inline void print_plot_arrangement<1>(char* gnuplot, size_t print_length) {
  snprintf(gnuplot, print_length, gnu_single, params::title);
}

//! Clear multiplot arrangement.
template <size_t N>
void finalize_plot_arrangement(char* gnuplot) {
  std::strcat(gnuplot, "\nunset multiplot");
}

//! Specialization of finalize_plot_arrangement(char).
template <>
inline void finalize_plot_arrangement<1>(char*) {}

template <size_t D, typename Sp, typename... S>
void write_plot_config(Model<D, Sp, S...> const& model, const char* directory,
                       const char* const* names, SaveParams const& save) {
  FILE* gnu;
  char plot_name[BUFFER_LENGTH_L2],
      data_loc[BUFFER_LENGTH];  // the data file location
  snprintf(plot_name, BUFFER_LENGTH_L2, PHASEFIELD_PLOT_LOC_FMT, directory,
           PHASEFIELD_DATA_NAME);

  if ((gnu = fopen(plot_name, "w")) == 0) {
    symphas::lib::make_directory_for_file(plot_name);
    if ((gnu = fopen(plot_name, "w")) == 0) {
      fprintf(SYMPHAS_ERR, SYMPHAS_MSG_ERR_FILE_OPEN, plot_name);
      exit(ERR_CODE_FILE_OPEN);
    }
  }
  const auto intervals = model.template system<0>().get_info().intervals;
  print_plot_ranges<D>(gnu, intervals);

  char gnuplot[BUFFER_LENGTH_L3];
  const char* sets[model_type_len<Model<D, Sp, S...>>];
  get_plot_fmt<D, S...>(sets);
  print_plot_arrangement<model_type_len<Model<D, Sp, S...>>>(gnuplot,
                                                             BUFFER_LENGTH_L3);

  char setcat[BUFFER_LENGTH_L2];

  // copy over the names of the systems for given system number
  for (size_t id = 0; id < model_type_len<Model<D, Sp, S...>>; ++id) {
    char fid[BUFFER_LENGTH_R2];
    snprintf(fid, BUFFER_LENGTH_R2, POSTFIX_ID_FMT, id);
    symphas::io::copy_data_file_name(DATA_DIR_RELATIVE_PLOT, save.get_stop(),
                                     get_system_id(id, model),
                                     DataFileType::SOLUTION_DATA, data_loc);

    /* in the case of a single output file, the correct index from the
     * list of saved data needs to be loaded
     */
    iter_type index =
        (params::single_output_file) ? std::max(0, save.num_saves() - 1) : 0;

    snprintf(setcat, STR_ARR_LEN(setcat), sets[id], data_loc, index, names[id]);
    std::strncat(
        gnuplot, setcat,
        std::min(sizeof(gnuplot) / sizeof(char) - 1 - std::strlen(gnuplot),
                 sizeof(setcat) / sizeof(char) - 1));
  }
  std::strcat(gnuplot, gnu_end);

  finalize_plot_arrangement<model_type_len<Model<D, Sp, S...>>>(gnuplot);

  fprintf(gnu, "%s", gnuplot);
  fclose(gnu);
}

#endif

//! \cond
DECLARE_SAVE_GRID_ALL_FUNCTIONS
//! \endcond

}  // namespace symphas::io::gp

//! \cond
#undef WIDTH_LATEX
#undef WIDTH_WIN
#undef MULTI_WIDTH_WIN
#undef MULTI_ALIGN_KEY
#undef MULTI_COLUMN_NUM
#undef MULTI_OUTPUT_WIDTH
#undef MULTI_MARGIN_RATIO
#undef MULTI_ALIGN_SEPARATION
#undef COLORBOX_HEIGHT_RATIO
#undef COLORBOX_BOTTOM_PADDING
#undef MULTI_ALIGN_WIDTH
#undef COLUMNS_LAST_ROW
#undef ADD_ROW_NUM
#undef MULTI_ROW_NUM
#undef MULTI_WIDTH_LATEX
#undef MULTI_WIDTH_RATIO
#undef MULTI_WIDTH_PLOT
#undef MULTI_HEIGHT_RATIO
#undef MULTI_HEIGHT_PLOT
#undef MULTI_HEIGHT_LATEX
#undef MULTI_ALIGN_MARGIN_FMT
#undef MULTI_ALIGN_MARGIN_VALUES
#undef MULTI_COLORBOX_FMT
#undef MULTI_COLORBOX_VALUES
#undef PLOT_SQUARE_DIMENSION_LATEX
#undef PLOT_SQUARE_DIMENSION_WIN
#undef PLOT_MULTI_DIMENSION_LATEX
#undef PLOT_MULTI_DIMENSION_WIN
#undef MULTI_PLOT_KEY_VPOS
#undef LATEX_KEY_FMT
#undef LATEX_SUBPLOT_KEY_FMT
#undef ALIGNMENT_FILE_NAME
#undef ALIGNMENT_PAR_SEP_CHAR
#undef SINGLE_WIDTH_LATEX
#undef SINGLE_ALIGN_MARGIN_VALUES
//! \endcond
