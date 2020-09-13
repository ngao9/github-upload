find_path(GNUPLOT_PATH
          NAMES gnuplot.exe
          PATHS "C:/Program Files/gnuplot/bin"
         )

if(GNUPLOT_PATH)
    set(GNUPLOT_FOUND TRUE)
    MESSAGE(${GNUPLOT_PATH})
else()
    set(GNUPLOT_FOUND FALSE)
endif()
