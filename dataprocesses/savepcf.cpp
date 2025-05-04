
#include "savepcf.h"

DLLPROC char const* apcfset[] = { R"~(
set xtics border autofreq
set ytics border autofreq)~"

#ifdef LATEX_PLOT 
R"~(
set term epslatex color size )~" PLOT_DISPLAY_SIZE_LATEX
#else
R"~(
set terminal )~" STR(GNU_PLOTTER) R"~( 1 size )~" PLOT_DISPLAY_SIZE_WIN
#endif

R"~(
set ylabel "correlation"
set xlabel "distance"

unset logscale
unset key

)~",

#ifdef LATEX_PLOT 
R"~(
set title "%s 2-Point Correlation Function Spherical Average"
set output "%s"
plot "%s" index %d w lines lc "black"
unset output
)~"
#else
R"~(
set title "2-Point Correlation Function Spherical Average"
plot "%s" index %d w lines
)~"
#endif
};

DLLPROC char const* pcfset[] = { R"~(
set xtics border autofreq
set ytics border autofreq)~"

#ifdef LATEX_PLOT 
R"~(
set term epslatex color size )~" PLOT_DISPLAY_SQSIZE_LATEX R"~(
set ylabel "$r_y$"
set xlabel "$r_x$")~"
#else
R"~(
set terminal )~" STR(GNU_PLOTTER) R"~( 1 size )~" PLOT_DISPLAY_SQSIZE_WIN R"~(
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
set title "%s 2-Point Correlation Function"
set output "%s"
plot "%s" matrix nonuniform index %d w image
unset output
)~"
#else
R"~(
set title "2-Point Correlation Function"
plot "%s" matrix nonuniform index %d w image
)~"
#endif
};


