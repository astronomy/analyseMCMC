# Compile AnalyseMCMC


.PHONY: clean cleanwarnings exit


#No text/space after statements below! (hence #)
COMP = ifort#              # Use Intel Fortran compiler
COMP = gfortran#           # Use GNU Fortran compiler (reads ~3x slower from disc than ifort!)
#COMP = gfortran-4.4.3#
#COMP = nagfor#            # Use NAG Fortran compiler

#WARN = yes#               # to compile with compile-time warnings
#RUNCHECK = yes#           # to compile with runtime-check options
#OPTIMISE = yes#           # to compile with some extra run-time optimisation (sets RUNCHECK=no)



#General:
CFLAGS = -O0
CFLAGS = -O2
LFLAGS = -lpgplot 



# Directories for the source, object, module, libary and executable files:
SRCDIR = src
OBJDIR = obj
MODDIR = mod
EXEDIR = .




#ifort:
ifeq ($(COMP),ifort)
   #To develop code  (don't mix -check bounds with -O>0!):
   ifeq ($(RUNCHECK),yes)
      CFLAGS = -O0 -check all -check noarg_temp_created
   endif
   CFLAGS += -traceback 
   ifeq ($(OPTIMISE),yes)
      CFLAGS = -O3

      #CFLAGS += -arch SSE2  #SSE2 needs -lsvml  (?)

      #CFLAGS += -xsse4.1
      #CFLAGS += -xhost
      #CFLAGS += -static
      #CFLAGS += -prof-gen -prof-dir./  #Don't use -ipo. Run the code thus compiled (a few times) and use -prof-use below for the final version
      #CFLAGS += -prof-use -prof-dir./

      #CFLAGS += -ipo
      #LFLAGS += -ipo -vec-report0
   endif
   ifeq ($(WARN),yes)
      CFLAGS += -warn all
   endif
   CFLAGS += -vec-report0
   CFLAGS += -module $(MODDIR) -I$(MODDIR)

   #To debug code:
   #CFLAGS += -O0  -vec-report0 -g

   CFLAGS += -nogen-interfaces #For ifort v.11

endif


#gfortran:
#ifeq ($(COMP),gfortran)
ifneq (,$(filter $(COMP),gfortran gfortran-4.4.3))
   ifeq ($(RUNCHECK),yes)
      CFLAGS = -O0 -fbounds-check
      #CFLAGS += -fcheck-array-temporaries  #Doesn't seem to work...
   endif
   ifeq ($(COMP),gfortran-4.4.3)
      CFLAGS += -fbacktrace
   endif
   ifeq ($(OPTIMISE),yes)
      CFLAGS = -O3 
   endif
   ifeq ($(WARN),yes)
      CFLAGS += -Wall
   endif
   CFLAGS += -I$(MODDIR) -J$(MODDIR)
   CFLAGS += -fmax-errors=10
endif



#NAGfort:
ifeq ($(COMP),nagfor)
   ifeq ($(RUNCHECK),yes)
      CFLAGS = -O0  -C=all -gline
   endif
   ifeq ($(OPTIMISE),yes)
      CFLAGS = -O4
   endif
   ifeq ($(WARN),yes)
      #CFLAGS += 
   endif
   CFLAGS += -colour
   CFLAGS += -I$(MODDIR) -mdir $(MODDIR)
   #CFLAGS += -I/opt/NAG/lib
endif










STDOBJ = $(OBJDIR)/svn_revision.o
ANALYSEMCMC_OBJ = $(OBJDIR)/analyseMCMC_modules.o $(OBJDIR)/analyseMCMC_functions.o $(OBJDIR)/analyseMCMC_textroutines.o 
ANALYSEMCMC_OBJ += $(OBJDIR)/analyseMCMC_stats.o $(OBJDIR)/analyseMCMC_tailored.o 
ANALYSEMCMC_OBJ += $(OBJDIR)/analyseMCMC_chains.o $(OBJDIR)/analyseMCMC_1dpdfs.o 
ANALYSEMCMC_OBJ += $(OBJDIR)/analyseMCMC_2dpdf_skymap.o $(OBJDIR)/analyseMCMC_2dpdf_routines.o $(OBJDIR)/analyseMCMC_2dpdfs.o 
ANALYSEMCMC_OBJ += $(OBJDIR)/analyseMCMC_animation.o $(OBJDIR)/analyseMCMC_plot.o $(OBJDIR)/analyseMCMC_main.o

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



