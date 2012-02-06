set( AnalyseMCMC_SRC_FILES
  src/analyseMCMC_1dpdfs.f90
  src/analyseMCMC_2dpdf_routines.f90
  src/analyseMCMC_2dpdf_skymap.f90
  src/analyseMCMC_2dpdfs.f90
  src/analyseMCMC_animation.f90
  src/analyseMCMC_chains.f90
  src/analyseMCMC_functions.f90
  src/analyseMCMC_main.f90
  src/analyseMCMC_modules.f90
  src/analyseMCMC_plot.f90
  src/analyseMCMC_stats.f90
  src/analyseMCMC_tailored.f90
  src/analyseMCMC_textroutines.f90
  )

set( plotSignal_SRC_FILES
  src/plotSignal.f90
  )

# Source files specific to PGPlot or PLplot:
if( PLplot_FOUND )
  set( Plot_SRC_FILES
    src/PG2PLplot.f90
    src/code_version_plplot.f90
    )
else( PLplot_FOUND )
  set( Plot_SRC_FILES
    src/code_version_pgplot.f90
    )
endif( PLplot_FOUND )

