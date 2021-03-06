## SetCompileOptions.cmake: 
## Various compile/optimisation options that we may want to enable:


option( CMAKE_VERBOSE_MAKEFILE  "Verbose makefile" off )
option( CREATE_VERSION          "Create a code-version file with git version and date (for developers)" off )

option( WANT_CHECKS             "Activate runtime checks (array bounds, NaNs)" off )
option( WANT_WARNINGS           "Activate warnings" on )
option( STOP_ON_WARNING         "Stop compilation on warnings" off )

option( WANT_OPENMP             "Use OpenMP parallelisation (experimental)" off )
option( WANT_SSE42              "Enable generation of SSE4.2 code" off )
option( WANT_HOST_OPT           "Enable host-specific optimisation. Choose only when compiling and running on the same machine! Overrides WANT_SSE42" off )
option( WANT_IPO                "Inter-procedural optimisation" off )
option( WANT_STATIC             "Generate statically linked executable" off )

option( WANT_LIBRARY            "Use compiler options to create libraries" off )
option( CREATE_SHAREDLIB        "Create shared libraries" off )
option( CREATE_STATICLIB        "Create static libraries" off )
option( COMPILER_SPECIFIC_LIBS  "Create compiler-specific libraries" off )

option( FORCE_PLPLOT            "Force the use of PLplot instead of PGPLOT, even if PGPLOT is found" off )

