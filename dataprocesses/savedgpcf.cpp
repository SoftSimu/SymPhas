
#include "savedgpcf.h"


DLLPROC char const* dpset[] = { R"~(
set xtics border autofreq
set ytics border autofreq)~"

#ifdef LATEX_PLOT 
R"~(
set term epslatex color size )~" PLOT_DISPLAY_SIZE_LATEX
#else
R"~(
set terminal wi)~" STR(GNU_PLOTTER) R"~(ndows 1 size )~" PLOT_DISPLAY_SIZE_WIN
#endif

R"~(
set ylabel "Domain Length"
set xlabel "Time (t)"

set logscale
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
plot "%s" w lines
)~"
#endif
};

DLLPROC char const* adpset[] = { R"~(
set xtics border autofreq
set ytics border autofreq)~"

#ifdef LATEX_PLOT 
R"~(
set term epslatex color size )~" PLOT_DISPLAY_SIZE_LATEX
R"~(
set ylabel "Correlation"
set xlabel "Distance ($|r|$)"
)~"
#else
R"~(
set terminal )~" STR(GNU_PLOTTER) R"~( 1 size )~" PLOT_DISPLAY_SIZE_WIN
R"~(
set ylabel "Correlation"
set xlabel "Distance (|r|)"
)~"
#endif

R"~(

unset logscale
unset key

)~",

#ifdef LATEX_PLOT 
R"~(
set title "%s Correlation Function for Different Times"
set output "%s"
plot "%s" index %d u 1:2 w lines, 0 w lines lt 0
unset output
)~"
#else
R"~(
set title "Correlation Function for Different Times"
plot "%s" u 1:2 w linespoints pt 1 ps 0.2, 0 w lines lt 0
)~"
#endif
};

DLLPROC char const* pdpset[] = { R"~(
set xtics border autofreq
set ytics border autofreq)~"

#ifdef LATEX_PLOT 
R"~(
set term epslatex color size )~" PLOT_DISPLAY_SIZE_LATEX
R"~(
set ylabel "Domain Length"
set xlabel "Time Index"
)~"
#else
R"~(
set terminal )~" STR(GNU_PLOTTER) R"~( 1 size 700,500)~"// PLOT_DISPLAY_SIZE_WIN
R"~(
set ylabel "Domain Length"
set xlabel "Time Index"
)~"
#endif

R"~(

set logscale
set key top left



d=0.33
a=1
a_quarter=1
a_third=1
a_half=1

f_fit(x)=a *x**d
f_quarter(x)=a_quarter*x**(1./4)
f_third(x)=a_third*x**(1./3)
f_half(x)=a_half*x**(1./2)

)~",

#ifdef LATEX_PLOT 
R"~(
set title "%s Domain Length"
set output "%s"
plot "%s" index %d u 1:2 w lines, f(x) title fit_title
unset output
)~"
#else
R"~(

set title "Domain Length"
data_file_plot="%s"
data_file_index=%d
plot data_file_plot index data_file_index u 1:2 w linespoints title "2-pt", \
data_file_plot index data_file_index u 1:3 w linespoints title "full", \
f_fit(x) w lines lt 0 lc 7 title fit_title
#, f_quarter(x) w lines dt 3 lc -1 notitle, f_third(x) w lines dt 3 lc -1 notitle, f_half(x) w lines dt 3 lc -1 notitle
)~",

R"~(


#ranges
%s

source="%s"
fpt=%lf

fit [fpt:*] f_fit(x) source u 1:2 via a,d
fit [fpt:*] f_quarter(x) source u 1:2 via a_quarter
fit [fpt:*] f_third(x) source u 1:2 via a_third
fit [fpt:*] f_half(x) source u 1:2 via a_half


#fit title
%s
)~"
#endif
};

DLLPROC const char* fit_title_str_dp = "fit_title = sprintf(\"x^{%.3f}\", d)";




