
COMP = ifort

#To develop code  (don't mix -check bounds with -O>0!):
CFLAGS = -free -O0  -vec-report0 -warn alignments -warn declarations -warn uncalled -warn usage -check bounds #-check format -check output_conversion -traceback

#To debug code:
#CFLAGS = -free -O0  -vec-report0 -g

#To optimise code by hand (memory):
#CFLAGS = -free -O0  -vec-report0 -warn alignments -warn declarations -warn uncalled -warn usage -check bounds -check format -check output_conversion -traceback -warn unused

#To use optimum run speed (~70s compile time):
#CFLAGS = -free -O3  -vec-report0 -traceback



#Gentoo:
LFLAGS = -lpgplot_ifort -L/usr/X11R6/lib -lX11 -lm -lifcore -limf -lirc -lsvml -lifport -Xlinker 


#gfortran needs ~60% more time (mostly for reading data in plotspins, it seems)
#COMP = gfortran
#CFLAGS = -O0 -ffree-form -ffree-line-length-none
#CFLAGS = -O3 -ffree-form -ffree-line-length-none
#LFLAGS = -lpgplot_ifort -L/usr/X11R6/lib -lX11 -lm -L/opt/intel/fc/10.0.026/lib -lifcore -limf -lirc -lsvml -lifport -Xlinker



PLOTSPINS_OBJ = plotspins_functions.o plotspins.o
ANALYSEMCMC_OBJ = analysemcmc_modules.o analysemcmc_functions.o analysemcmc_stats.o analysemcmc_chains.o analysemcmc_1dpdfs.o analysemcmc_2dpdfs.o analysemcmc_animation.o analysemcmc_plot.o analysemcmc.o


.f.o: 
	$(COMP) $(CFLAGS) -c $?

exit:

all: comp_pdfs analysemcmc #mcmcstats plotdata plotsignal 


analysemcmc: $(STDOBJ) ${ANALYSEMCMC_OBJ}
	$(COMP) -o analysemcmc $(LFLAGS) $(STDOBJ) ${ANALYSEMCMC_OBJ}

comp_pdfs: $(STDOBJ) comp_pdfs_functions.o comp_pdfs.o
	$(COMP) $(LFLAGS) -o comp_pdfs $(STDOBJ) comp_pdfs_functions.o comp_pdfs.o

plotsignal: $(STDOBJ) plotsignal.o
	$(COMP) $(LFLAGS) -o plotsignal $(STDOBJ) plotsignal.o

plotspins: $(STDOBJ) ${PLOTSPINS_OBJ}
	$(COMP) -o plotspins $(LFLAGS) $(STDOBJ) ${PLOTSPINS_OBJ}

#plotspins: $(STDOBJ) plotspins_functions.o plotspins.o
#	$(COMP) -o plotspins $(LFLAGS) $(STDOBJ) plotspins_functions.o plotspins.o

mcmcstats: $(STDOBJ) mcmcstats.o
	$(COMP) $(LFLAGS) -o mcmcstats $(STDOBJ) mcmcstats.o

clean:
	rm -f *.o *.mod







