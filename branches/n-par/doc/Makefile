
COMP = ifort

#To develop code  (don't mix -check bounds with -O>0!):
#CFLAGS = -free -O0  -vec-report0 -warn alignments -warn declarations -warn uncalled -warn usage -check bounds -traceback #-check format -check output_conversion
#CFLAGS = -free -O0  -traceback -warn all #Lots of compile-time checks
CFLAGS = -free -O0  -traceback -warn all -check all -check noarg_temp_created  #Lots of compile-time and run-time checks

#To debug code:
#CFLAGS = -free -O0  -vec-report0 -g

#To optimise code by hand (memory):
#CFLAGS = -free -O0  -vec-report0 -warn alignments -warn declarations -warn uncalled -warn usage -check bounds -check format -check output_conversion -traceback -warn unused
#CFLAGS = -free -O0  -vec-report0 -warn alignments -warn declarations -warn uncalled -warn usage -check bounds -check format -traceback -warn unused

#To use optimum run speed (~12s compile/link time):
CFLAGS = -free -O3 -vec-report0 -traceback
#CFLAGS = -free -O3 -ipo -vec-report0 -traceback
#CFLAGS = -free -arch SSE2 -O3 -vec-report0 -traceback
##CFLAGS = -free -arch SSE2 -O3 -ipo -vec-report0 -traceback

#Gentoo:
LFLAGS = -lpgplot_ifort -L/usr/X11R6/lib -lX11 -lm -lifcore -limf -lirc -lsvml -lifport -Xlinker 
#LFLAGS = -ipo -vec-report0 -lpgplot_ifort -L/usr/X11R6/lib -lX11 -lm -lifcore -limf -lirc -lsvml -lifport -Xlinker 


#gfortran: similar in speed to ifort, but file reading is ~2.5x slower!
#COMP = gfortran
#CFLAGS = -O0 -ffree-form -ffree-line-length-none
#CFLAGS = -O0 -ffree-form -ffree-line-length-none   -Wall 
#CFLAGS = -O -ffree-form -ffree-line-length-none
#CFLAGS = -O3 -ffree-form -ffree-line-length-none
#LFLAGS = -lpgplot_ifort -L/usr/X11R6/lib -lX11 -lm -L/opt/intel/fc/10.0.026/lib -lifcore -limf -lirc -lsvml -lifport -Xlinker 





#STDOBJ = functions.o
ANALYSEMCMC_OBJ = analysemcmc_modules.o analysemcmc_functions.o analysemcmc_stats.o analysemcmc_chains.o analysemcmc_1dpdfs.o analysemcmc_2dpdfs.o analysemcmc_animation.o analysemcmc_plot.o analysemcmc.o

#.f.o: 
#	$(COMP) $(CFLAGS) -c $?

%.o: %.f Makefile analysemcmc_modules.f
	$(COMP) $(CFLAGS) -c  $<

exit:

all: comp_pdfs analysemcmc mcmcstats #plotdata plotsignal 



analysemcmc: Makefile ${ANALYSEMCMC_OBJ}
	$(COMP) -o analysemcmc $(LFLAGS) $(STDOBJ) ${ANALYSEMCMC_OBJ}

comp_pdfs: Makefile $(STDOBJ) comp_pdfs_functions.o comp_pdfs.o
	$(COMP) $(LFLAGS) -o comp_pdfs $(STDOBJ) comp_pdfs_functions.o comp_pdfs.o

plotspectrum: Makefile $(STDOBJ) plotspectrum.o
	$(COMP) $(LFLAGS) -o plotspectrum $(STDOBJ) plotspectrum.o

plotsignal: Makefile $(STDOBJ) plotsignal.o
	$(COMP) $(LFLAGS) -o plotsignal $(STDOBJ) plotsignal.o

mcmcstats: Makefile $(STDOBJ) mcmcstats.o
	$(COMP) $(LFLAGS) -o mcmcstats $(STDOBJ) mcmcstats.o


clean:
	rm -f *.o *.mod









#Compile/run benchmark on Gentoo, 20/12/2008 (short run, spread in run times for each binary may be comparable to that between binaries):
#Compiler			Comp.time:	Run time:
#gfortran -O0			3.7		13.3
#gfortran -O3			27.0		12.0
#ifort -O0			2.0		12.4
#ifort -O3			14.0		11.9
#ifort -O3 -arch SSE2		69.0		11.7
#ifort -O3 -ipo			12.9		11.2	<--			
#ifort  -arch SSE2 -O3 -ipo	36.7		11.5



