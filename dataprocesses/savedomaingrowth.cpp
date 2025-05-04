
#include "savedomaingrowth.h"


DLLPROC char const* dgset[] = { R"~(
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
set ylabel "Domain Length"
set xlabel "Time (t)"

unset logscale
unset key

)~",

#ifdef LATEX_PLOT 
R"~(
set title "%s Domain Growth"
set output "%s"
plot "%s" index %d w lines
unset output
)~"
#else
R"~(
set title "Domain Length"
plot "%s" index %d w lines
)~"
#endif
};


DLLPROC char const* adgset[] = { R"~(
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
set ylabel "Domain Length"
set xlabel "Time (t)"

set logscale
set key

)~",

#ifdef LATEX_PLOT 
R"~(
set title "%s Domain Length"
set output "%s"
plot "%s" index %d u 1:2 w points notitle, f(x) title fit_title
unset output
)~"
#else
R"~(

set title "Domain Length"
plot "%s" index %d u 1:2 w linespoints notitle, f(x) w lines lt 0 lc 7 title fit_title
)~",

R"~(

#ranges
%s

d=1./3
a=1
f(x)=a *x**d
fit f(x) "%s" via a,d

#fit title
%s

)~"
#endif
};

DLLPROC const char* fit_title_str_dg = "fit_title = sprintf(\"x^{%.3f}\", d)";

