# Compile AnalyseMCMC


.PHONY: clean cleanwarnings exit


#No text/space after statements below! (hence #)
COMP = gfortran#  Use GNU Fortran compiler
#COMP = gfortran-4.4.3#
#COMP = ifort#     Use Intel Fortran compiler
#COMP = nagfor#     Use NAG Fortran compiler

WARN = no
#WARN = yes#  'yes#' to compile with compile-time warnings
RUNCHECK = no
#RUNCHECK = yes#  'yes#' to compile with runtime-check options

#General:
CFLAGS = -O2
#CFLAGS = -O3
LFLAGS = -lpgplot 



# Directories for the source, object, module, libary and executable files:
SRCDIR = src
OBJDIR = obj
MODDIR = mod
EXEDIR = .




#ifort:
ifeq ($(COMP),ifort)
   #To develop code  (don't mix -check bounds with -O>0!):
   ifeq ($(WARN),yes)
      CFLAGS = -O0 -warn all
   endif
   ifeq ($(RUNCHECK),yes)
      CFLAGS += -check all -check noarg_temp_created
   endif
   CFLAGS += -vec-report0 -traceback 
   CFLAGS += -module $(MODDIR) -I$(MODDIR)

   #To debug code:
   #CFLAGS += -free -O0  -vec-report0 -g

   #To use optimum run speed
   #CFLAGS += -free -O3 -vec-report0 -traceback
   #CFLAGS += -free -O3 -arch SSE2 -vec-report0 -traceback           #SSE2 needs -lsvml
   ##CFLAGS += -free -O3 -arch SSE2 -ipo -vec-report0 -traceback

   #CFLAGS += -xsse4.1
   #CFLAGS += -xhost
   #CFLAGS += -static
   #CFLAGS += -prof-gen -prof-dir./  #Don't use -ipo. Run the code thus compiled (a few times) and use -prof-use below for the final version
   #CFLAGS += -prof-use -prof-dir./

   CFLAGS += -nogen-interfaces #For ifort v.11

   #CFLAGS += -ipo
   #LFLAGS += -ipo -vec-report0
endif


#gfortran: similar in speed to ifort, but file reading is ~2.5x slower!
#ifeq ($(COMP),gfortran)
ifneq (,$(filter $(COMP),gfortran gfortran-4.4.3))
   ifeq ($(WARN),yes)
      CFLAGS = -O0 -Wall
   endif
   ifeq ($(RUNCHECK),yes)
      CFLAGS += -fbounds-check
      #CFLAGS += -fcheck-array-temporaries  #Doesn't seem to work...
   endif
   CFLAGS += -I$(MODDIR) -J$(MODDIR)
   CFLAGS += -fmax-errors=10
   #CFLAGS += -ffree-line-length-none
endif
ifeq ($(COMP),gfortran-4.4.3)
   CFLAGS += -fbacktrace
endif



#NAGfort:
ifeq ($(COMP),nagfor)
  CFLAGS = -colour
  CFLAGS += -gline -C=all
  CFLAGS += -I$(MODDIR) -mdir $(MODDIR)
  CFLAGS += -I/opt/NAG/lib
endif










STDOBJ = $(OBJDIR)/svn_revision.o
ANALYSEMCMC_OBJ = $(OBJDIR)/analyseMCMC_modules.o $(OBJDIR)/analyseMCMC_functions.o $(OBJDIR)/analyseMCMC_textroutines.o 
ANALYSEMCMC_OBJ += $(OBJDIR)/analyseMCMC_stats.o $(OBJDIR)/analyseMCMC_tailored.o 
ANALYSEMCMC_OBJ += $(OBJDIR)/analyseMCMC_chains.o $(OBJDIR)/analyseMCMC_1dpdfs.o 
ANALYSEMCMC_OBJ += $(OBJDIR)/analyseMCMC_2dpdf_skymap.o $(OBJDIR)/analyseMCMC_2dpdf_routines.o $(OBJDIR)/analyseMCMC_2dpdfs.o 
ANALYSEMCMC_OBJ += $(OBJDIR)/analyseMCMC_animation.o $(OBJDIR)/analyseMCMC_plot.o $(OBJDIR)/analyseMCMC.o

ifeq ($(COMP),nagfor)
  ANALYSEMCMC_OBJ += $(OBJDIR)/nagfor.o
endif


$(OBJDIR)/%.o: $(SRCDIR)/%.f90 Makefile $(SRCDIR)/svn_revision.f90
	$(COMP) $(CFLAGS) -c $< -o $@


analyseMCMC:

exit:

all: compPDFs analyseMCMC MCMCstats #plotdata plotsignal 


$(OBJDIR)/$(ANALYSEMCMC_OBJ): $(SRCDIR)/analyseMCMC_modules.f90  #All analyseMCMC*.o files depend on analyseMCMC_modules.f90

analyseMCMC: $(STDOBJ) ${ANALYSEMCMC_OBJ}
	$(COMP) -o $(EXEDIR)/analyseMCMC $(LFLAGS) $(STDOBJ) ${ANALYSEMCMC_OBJ}

compPDFs: $(STDOBJ) $(OBJDIR)/compPDFs_functions.o $(OBJDIR)/compPDFs.o
	$(COMP) $(LFLAGS) -o $(EXEDIR)/compPDFs $(STDOBJ) $(OBJDIR)/compPDFs_functions.o $(OBJDIR)/compPDFs.o

plotspectrum: $(STDOBJ) $(OBJDIR)/plotspectrum.o
	$(COMP) $(LFLAGS) -o $(EXEDIR)/$(OBJDIR)/plotspectrum $(STDOBJ) $(OBJDIR)/plotspectrum.o

plotsignal: $(STDOBJ) $(OBJDIR)/plotsignal.o
	$(COMP) $(LFLAGS) -o $(EXEDIR)/$(OBJDIR)/plotsignal $(STDOBJ) $(OBJDIR)/plotsignal.o

MCMCstats: $(STDOBJ) $(OBJDIR)/MCMCstats.o
	$(COMP) $(LFLAGS) -o $(EXEDIR)/MCMCstats $(STDOBJ) $(OBJDIR)/MCMCstats.o

$(SRCDIR)/svn_revision.f90: .svn/entries   #.svn/entries changes when committing/updating
	@echo 'Generating '$(SRCDIR)'/svn_revision.f90'
	@echo '!svn_revision.f90: temporary source file created by Makefile to report the svn revision' > $(SRCDIR)/svn_revision.f90
	@echo 'module svn_revision' >> $(SRCDIR)/svn_revision.f90
	@echo "   character, parameter :: svnrevision*99 = '"`svnversion`"'" >> $(SRCDIR)/svn_revision.f90
	@echo 'end module svn_revision' >> $(SRCDIR)/svn_revision.f90

dirs:
	@if [ ! -e $(OBJDIR) ]; then mkdir ${OBJDIR}; fi 
	@if [ ! -e $(MODDIR) ]; then mkdir ${MODDIR}; fi
	@if [ ! -e $(EXEDIR) ]; then mkdir ${EXEDIR}; fi

clean:
	rm -f $(OBJDIR)/*.o $(MODDIR)/*.mod $(SRCDIR)/svn_revision.f90

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



