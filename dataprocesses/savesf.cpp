
#include "savesf.h"


DLLPROC char const* sfset[] = { R"~(
set xtics border autofreq
set ytics border autofreq)~"

#ifdef LATEX_PLOT 
R"~(
set term epslatex color size )~" PLOT_DISPLAY_SQSIZE_LATEX R"~(
set ylabel "$q_y$"
set xlabel "$q_x$"
set cblabel "$S(\\bar{q})$")~"
#else
R"~(
set terminal )~" STR(GNU_PLOTTER) R"~( 1 size )~" PLOT_DISPLAY_SQSIZE_WIN R"~(
set ylabel "q_y"
set xlabel "q_x"
set cblabel "S(q)")~"
#endif

R"~(
unset logscale
unset key

)~",

#ifdef LATEX_PLOT 
R"~(
set title "%s Structure Factor"
set output "%s"
plot "%s" matrix nonuniform index %d w image
unset output)~"
#else
R"~(
set title "Structure Factor"
plot "%s" matrix nonuniform index %d w image
)~"
#endif
};

DLLPROC char const* asfset[] = { R"~(

C1="#d73027"
C2="#f46d43"
C3="#fdae61"
C4="#fee08b"
C5="#d9ef8b"
C6="#a6d96a"
C7="#66bd63"
C8="#1a9850"

set xtics border autofreq
set ytics border autofreq)~"

#ifdef LATEX_PLOT 
R"~(
set linetype 1 pt 7 lc rgb C1
set linetype 2 pt 7 lc rgb C2
set linetype 3 pt 7 lc rgb C3
set linetype 4 pt 7 lc rgb C4
set linetype 5 pt 7 lc rgb C5
set linetype 6 pt 7 lc rgb C6
set linetype 7 pt 7 lc rgb C7
set linetype 8 pt 7 lc rgb C8
)~"
#else
R"~(
set linetype 1 pt 7 lc rgb C1 ps 0.25
set linetype 2 pt 7 lc rgb C2 ps 0.25
set linetype 3 pt 7 lc rgb C3 ps 0.25
set linetype 4 pt 7 lc rgb C4 ps 0.25
set linetype 5 pt 7 lc rgb C5 ps 0.25
set linetype 6 pt 7 lc rgb C6 ps 0.25
set linetype 7 pt 7 lc rgb C7 ps 0.25
set linetype 8 pt 7 lc rgb C8 ps 0.25
)~"
#endif


#ifdef LATEX_PLOT 
R"~(
set term epslatex color size )~" PLOT_DISPLAY_SIZE_LATEX R"~(
set ylabel "$S(q)$"
set xlabel "$|\\bar{q}|$")~"
#else
R"~(
set terminal )~" STR(GNU_PLOTTER) R"~( 1 size )~" PLOT_DISPLAY_SIZE_WIN R"~(
set ylabel "S(q)"
set xlabel "|q|)~"
#endif

R"~(
set logscale
set key inside right top Right


d=-3
a=30
f(x)=x**d / a

)~",

#ifdef LATEX_PLOT 
R"~(
set title "%s Structure Factor Spherical Average"
set output "%s"
plot "%s" index %d w points notitle, [0.1:5] f(x) title fx_title
unset output
)~"
#else
R"~(
set title "Structure Factor Spherical Average"
plot "%s" index %d w points notitle, [0.1:5] f(x) title fx_title
)~"
#endif
,
R"~(

fit [0.1:*] f(x) "%s" via a
fx_title = sprintf("S(q) \sim q^{-(%d + 1)}")

#ranges%s

)~"

};


