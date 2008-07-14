
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



#Comp time (plotspins, ifort, no -check bounds):
# O0: 1.5s
# O1: 8.0s
# O2: 44.4s
# O3: 68.5s





#STDOBJ = functions.o


.f.o: 
	$(COMP) $(CFLAGS) -c $?

exit:


plotsignal: $(STDOBJ) plotsignal.o
	$(COMP) $(LFLAGS) -o plotsignal $(STDOBJ) plotsignal.o

plotspins: $(STDOBJ) plotspins_functions.o plotspins.o
	$(COMP) -o plotspins $(LFLAGS) $(STDOBJ) plotspins_functions.o plotspins.o

#mcmcstats: $(STDOBJ) mcmcstats.o
#	$(COMP) $(LFLAGS) -o mcmcstats $(STDOBJ) mcmcstats.o

clean:
	rm -f *.o *.mod







