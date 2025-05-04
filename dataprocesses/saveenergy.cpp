
#include "saveenergy.h"


DLLPROC char const* enset[] = { R"~(
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
set ylabel "energy"
set xlabel "time index"

unset logscale
unset key

)~",

#ifdef LATEX_PLOT 
R"~(
set title "%s Energy Curve Over the System Time Evolution"
set output "%s"
plot "%s" index %d w lines
unset output
)~"
#else
R"~(
set title "Energy Curve Over the System Time Evolution"
plot "%s" index %d w lines
)~"
#endif
};

