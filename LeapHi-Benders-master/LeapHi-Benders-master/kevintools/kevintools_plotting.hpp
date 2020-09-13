#pragma once

#include "kevintools_typedefs.hpp"

#include "gnuplot.h"

namespace kevintools{

    std::string rgb(int r, int g, int b);
    std::string argb(int alpha, int r, int g, int b);
    std::string argb(const Vint& argb);
    void scatter(GnuplotPipe& gp, const VVdouble& data, bool replot=false, int pointtype=1, std::string color="black");
    void lines(GnuplotPipe& gp, const VVdouble& from, const VVdouble& to, bool replot, int pointtype, const std::string& color, int linewidth=1);
    void arrows(GnuplotPipe& gp, const VVdouble& from, const VVdouble& to, bool replot, const std::string& color, int linewidth);

}