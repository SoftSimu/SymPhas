
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


#include "savedefines.h"


DLLIO const char* symphas::io::adefset[] = { R"~(
set xtics border autofreq
set ytics border autofreq)~"

#ifdef LATEX_PLOT 
R"~(
set term epslatex color size )~" PLOT_DISPLAY_SIZE_LATEX
#else
R"~(
set terminal windows 1 size )~" PLOT_DISPLAY_SIZE_WIN
#endif

R"~(
set ylabel "Y"
set xlabel "X"
unset key

)~",

#ifdef LATEX_PLOT 
R"~(
set title "%s Average Field Values at $t = %.2lf$"
set output %s"
plot "%s" index %d u 1:2 w lines
unset output
)~"
#else
R"~(
set title "Average Field Values"
plot "%s" index %d u 1:2 w lines
)~"
#endif
};

DLLIO const char* symphas::io::defset[] = { R"~(
set xtics border autofreq
set ytics border autofreq)~"

#ifdef LATEX_PLOT 
R"~(
set term epslatex color size )~" PLOT_DISPLAY_SQSIZE_LATEX R"~(
set ylabel "$r_y$"
set xlabel "$r_x$")~"
#else
R"~(
set terminal windows 1 size )~" PLOT_DISPLAY_SQSIZE_WIN R"~(
set ylabel "r_y"
set xlabel "r_x")~"
#endif

R"~(
set cblabel "correlation"

unset logscale
unset key

)~",

#ifdef LATEX_PLOT 
R"~(
set title "%s Data Field Values at $t = %.2lf$"
set output %s"
plot "%s" matrix nonuniform w image
unset output
)~"
#else
R"~(
set title "Data Field Values"
plot "%s" matrix nonuniform w image
)~"
#endif
};


