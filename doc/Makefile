# Compile AnalyseMCMC
#
# Dependencies:
#
# At compile/link time:
# - PGPlot:   http://www.astro.caltech.edu/~tjp/pgplot/
# or:
# - PLplot:   http://plplot.sourceforge.net/
#
# - libSUFR:  http://libsufr.sourceforge.net/
#
# At runtime:
# - ImageMagick:  http://www.imagemagick.org
# - a2ps:         http://www.inf.enst.fr/~demaille/a2ps/
# - epsropdf:     http://www.tex.ac.uk/tex-archive/help/Catalogue/entries/epstopdf.html, included in TexLive (http://www.tug.org/texlive/)


.PHONY: clean cleanwarnings code_version exit 


#No text/space after statements below! (hence #)
COMP = gfortran#           # Use GNU Fortran compiler (reads ~3x slower from disc than ifort!)
#COMP = gfortran-4.6#       # Use GNU Fortran compiler.  Need >= 4.5
#COMP = g95#               # USE g95 compiler
#COMP = ifort#             # Use Intel Fortran compiler
#COMP = nagfor#            # Use NAG Fortran compiler

WARN = yes#               # to compile with compile-time warnings
#RUNCHECK = yes#           # to compile with runtime-check options
#DEBUG = yes#             # compile with debug options
#OPTIMISE = yes#           # to compile with some extra run-time optimisation (sets RUNCHECK=no)
F2003 = yes#              # compile with Fortran 2003 standard
#STATIC = yes#              # link statically
#PLPLOT = yes#             # use PLplot instead of PGPlot



# General:
OPT = -O0
OPT = -O2
ifeq ($(PLPLOT),yes)
  LFLAGS += -lplplotf95d
  CFLAGS += -I/usr/lib64/fortran/modules/plplot/
else
  LFLAGS += -lpgplot
endif





# Directories for the source, object, module, libary and binary files:
SRCDIR = src
OBJDIR = obj
MODDIR = mod
BINDIR = .




# ifort:
ifeq ($(COMP),ifort)
   #ifeq ($(PLPLOT),yes)
   #  LFLAGS += -lplplotf95d
   #  CFLAGS += -I/home/sluys/usr/include/plplot/ifort
   #endif
   # Link against libSUFR:
   CFLAGS += -I/home/sluys/usr/include/libSUFR/ifort
   LFLAGS += -lSUFR_ifort

   #To develop code  (don't mix -check bounds with -O>0!):
   ifeq ($(RUNCHECK),yes)
      OPT = -O0
      CFLAGS += -check all -check noarg_temp_created
   endif

   ifeq ($(OPTIMISE),yes)
      OPT = -O2

      #OPT += -arch SSE2  #SSE2 needs -lsvml  (?)

      #OPT += -xsse4.1
      #OPT += -xhost
      #OPT += -static
      #OPT += -prof-gen -prof-dir./  #Don't use -ipo. Run the code thus compiled (a few times) and use -prof-use below for the final version
      #OPT += -prof-use -prof-dir./

      #OPT += -ipo
      #LFLAGS += -ipo -vec-report0
   endif
   ifeq ($(WARN),yes)
      CFLAGS += -warn all -warn errors -WB
   endif
   ifeq ($(F2003),yes)
      CFLAGS += -stand f03 -diag-disable 6894,8290,8291
   endif
   CFLAGS += -vec-report0 -traceback -error-limit 10
   CFLAGS += -module $(MODDIR) -I$(MODDIR)

   # Link statically
   ifeq ($(STATIC),yes)
      CFLAGS += -static
   endif

   # To debug code:
   ifeq ($(DEBUG),yes)
      OPT = -O0
      CFLAGS += -vec-report0 -g
   endif

   CFLAGS += -nogen-interfaces  # For ifort v.11

endif


# gfortran:
#ifeq ($(COMP),gfortran)
ifneq (,$(filter $(COMP),gfortran gfortran-4.6))
   # Link against libSUFR:
   CFLAGS += -I/home/sluys/usr/include/libSUFR/gfortran
   LFLAGS += -lSUFR_gfortran

   CFLAGS += -fwhole-file #-fwhole-program
   ifeq ($(COMP),gfortran-4.6)
      LFLAGS += -Wl,-rpath=/usr/lib/gcc/x86_64-pc-linux-gnu/4.6.2-pre9999/ #for 4.6
   endif

   ifeq ($(RUNCHECK),yes)
      OPT = -O0
      #CFLAGS += -fcheck=all
      CFLAGS += -fcheck=bounds,do,mem,pointer,recursion
      CFLAGS += -fbacktrace -g  # -g for line numbers
   endif


   ifeq ($(OPTIMISE),yes)
      OPT = -O2
   endif

   ifeq ($(WARN),yes)
      CFLAGS += -Wall -Wextra
      CFLAGS += -Werror -fmax-errors=10
      CFLAGS += -Wcharacter-truncation -Wunderflow
#      CFLAGS += -Wno-array-temporaries
#      CFLAGS += -Wconversion
#      CFLAGS += -Wimplicit-interface
# Use implicit none everywhere:
#      CFLAGS += -fimplicit-none
   endif

   #To debug code:
   ifeq ($(DEBUG),yes)
      OPT = -O0
      CFLAGS += -g
   endif

   ifeq ($(F2003),yes)
      CFLAGS += -std=f2008 -fall-intrinsics -pedantic
   else
      CFLAGS += -ffree-line-length-none
   endif

   # Link statically
   ifeq ($(STATIC),yes)
      CFLAGS += -static
   endif

   CFLAGS += -I$(MODDIR) -J$(MODDIR)
endif



# g95:
ifeq ($(COMP),g95)
   # Link against libSUFR:
   CFLAGS += -I/home/sluys/usr/include/libSUFR/g95
   LFLAGS += -lSUFR_g95

   ifeq ($(RUNCHECK),yes)
      OPT = -O0
      CFLAGS += -fbounds-check
      CFLAGS += -ftrace=full
   endif

   ifeq ($(WARN),yes)
      CFLAGS += -Wall -Wextra
      CFLAGS += -Wno=112,140,165
      #CFLAGS += -Wno=112,140,165,136
      #CFLAGS += -Wno=112,140,165,136,163
# Turn warnings into errors:
      CFLAGS += -Werror
   endif

   #To debug code:
   ifeq ($(DEBUG),yes)
      OPT = -O0
      CFLAGS += -g
   endif

   ifeq ($(F2003),yes)
      CFLAGS += -std=f2003
#   else
#      CFLAGS += 
   endif

   ifeq ($(OPENMP),yes)
      #CFLAGS += -fopenmp
      #LFLAGS += -lgomp
   endif

   # Link statically
   ifeq ($(STATIC),yes)
      CFLAGS += -static
   endif

   CFLAGS += -I$(MODDIR) -fmod=$(MODDIR)
endif


# NAGfort:
ifeq ($(COMP),nagfor)
   ifeq ($(RUNCHECK),yes)
      OPT = -O0
      CFLAGS = -C=all -gline
   endif
   ifeq ($(OPTIMISE),yes)
      OPT = -O4
   endif
   ifeq ($(WARN),yes)
      #CFLAGS += 
   endif
   CFLAGS += -colour
   CFLAGS += -I$(MODDIR) -mdir $(MODDIR)
   #CFLAGS += -I/opt/NAG/lib
endif










STDOBJ = $(OBJDIR)/analyseMCMC_modules.o $(OBJDIR)/code_version.o
ifeq ($(PLPLOT),yes)
  STDOBJ += $(OBJDIR)/PG2PLplot.o
endif

ANALYSEMCMC_OBJ =  $(OBJDIR)/analyseMCMC_functions.o $(OBJDIR)/analyseMCMC_textroutines.o 
ANALYSEMCMC_OBJ += $(OBJDIR)/analyseMCMC_stats.o $(OBJDIR)/analyseMCMC_tailored.o 
ANALYSEMCMC_OBJ += $(OBJDIR)/analyseMCMC_chains.o $(OBJDIR)/analyseMCMC_1dpdfs.o 
ANALYSEMCMC_OBJ += $(OBJDIR)/analyseMCMC_2dpdf_skymap.o $(OBJDIR)/analyseMCMC_2dpdf_routines.o $(OBJDIR)/analyseMCMC_2dpdfs.o 
ANALYSEMCMC_OBJ += $(OBJDIR)/analyseMCMC_animation.o $(OBJDIR)/analyseMCMC_plot.o $(OBJDIR)/analyseMCMC_main.o

COMPPDFS_OBJ =  $(OBJDIR)/compPDFs_functions.o


ifeq ($(COMP),g95)
  STDOBJ += $(OBJDIR)/nagfor.o
endif
ifeq ($(COMP),nagfor)
  STDOBJ += $(OBJDIR)/nagfor.o
endif




$(OBJDIR)/%.o: $(SRCDIR)/%.f90  #Makefile #$(SRCDIR)/code_version.f90
	$(COMP) $(OPT) $(CFLAGS) -c $< -o $@



all: analyseMCMC compPDFs MCMCstats plotSignal plotSpectrum #plotdata 

exit:



analyseMCMC_modules: $(SRCDIR)/analyseMCMC_modules.f90
	$(COMP) $(OPT) $(CFLAGS) -c $(SRCDIR)/analyseMCMC_modules.f90  -o $(OBJDIR)/analyseMCMC_modules.o

compPDFs_functions: $(SRCDIR)/compPDFs_functions.f90
	$(COMP) $(OPT) $(CFLAGS) -c $(SRCDIR)/compPDFs_functions.f90  -o $(OBJDIR)/compPDFs_functions.o




$(OBJDIR)/$(ANALYSEMCMC_OBJ): $(SRCDIR)/analyseMCMC_modules.f90  # All analyseMCMC*.o files depend on analyseMCMC_modules.f90
$(OBJDIR)/$(COMPPDFS_OBJ): $(SRCDIR)/compPDFs_functions.f90        # All compPDF*.o files depend on compPDFs_functions.f90

analyseMCMC: $(STDOBJ) ${ANALYSEMCMC_OBJ}
	$(COMP) -o $(BINDIR)/analyseMCMC $(LFLAGS) $(STDOBJ) ${ANALYSEMCMC_OBJ}

compPDFs: $(STDOBJ) $(OBJDIR)/compPDFs_functions.o $(OBJDIR)/compPDFs.o
	$(COMP) $(LFLAGS) -o $(BINDIR)/compPDFs $(STDOBJ) $(OBJDIR)/compPDFs_functions.o $(OBJDIR)/compPDFs.o

MCMCstats: $(STDOBJ) $(OBJDIR)/MCMCstats.o
	$(COMP) $(LFLAGS) -o $(BINDIR)/MCMCstats $(STDOBJ) $(OBJDIR)/MCMCstats.o

plotSignal: $(STDOBJ) $(OBJDIR)/plotSignal.o
	$(COMP) $(LFLAGS) -o $(BINDIR)/plotSignal $(STDOBJ) $(OBJDIR)/plotSignal.o

plotSpectrum: $(STDOBJ) $(OBJDIR)/plotSpectrum.o
	$(COMP) $(LFLAGS) -o $(BINDIR)/plotSpectrum $(STDOBJ) $(OBJDIR)/plotSpectrum.o


#$(SRCDIR)/code_version.f90: .svn/entries   #.svn/entries changes when committing/updating
$(SRCDIR)/code_version.f90: code_version    # regenerate code_version.f90; code_version is PHONY and always out-of-date
	@. ./code_version.sh src/code_version.f90 $(COMP) $(OPT) $(PLPLOT)


dirs:
	@if [ ! -e $(OBJDIR) ]; then mkdir ${OBJDIR}; fi 
	@if [ ! -e $(MODDIR) ]; then mkdir ${MODDIR}; fi
	@if [ ! -e $(BINDIR) ]; then mkdir ${BINDIR}; fi

clean:
	rm -f $(OBJDIR)/*.o $(MODDIR)/*.mod $(SRCDIR)/code_version.f90

cleanwarnings:
	rm -f $(SRCDIR)/*__genmod.f90 $(MODDIR)/*__genmod.mod









#Compile/run benchmark on Gentoo, 20/12/2008 (short run, spread in run times for each binary may be comparable to that between binaries):
#Compiler			Comp.time:	Run time:
#gfortran -O0			3.7		13.3
#gfortran -O3			27.0		12.0
#ifort -O0			2.0		12.4
#ifort -O3			14.0		11.9
#ifort -O3 -arch SSE2		69.0		11.7
#ifort -O3 -ipo			12.9		11.2	<--			
#ifort  -arch SSE2 -O3 -ipo	36.7		11.5



