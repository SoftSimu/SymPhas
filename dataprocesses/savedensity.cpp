
#include "savedensity.h"

DLLPROC char const* deset[] = {

R"~(
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
set ylabel "density"
set xlabel "time index"

unset logscale
unset key

)~",

#ifdef LATEX_PLOT 
R"~(
set title "%s Density Curve Over Coarsening"
set output "%s"
plot "%s" index %d w lines
unset output
)~"
#else
R"~(
set title "Density Curve Over Coarsening"
plot "%s" index %d w lines
)~"
#endif
};



DLLPROC char const* adeset[] = { R"~(
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
set ylabel "Density"
set xlabel "Time"
unset key

)~",

#ifdef LATEX_PLOT 
R"~(
set title "%s Density Curve Over Coarsening"
set output %s"
plot "%s" index %d w points
unset output
)~"
#else
R"~(
set title "Density Curve Over Coarsening"
plot "%s" index %d w points
)~"
#endif
};

