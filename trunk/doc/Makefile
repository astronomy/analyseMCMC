
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





STDOBJ = functions.o


.f.o: 
	$(COMP) $(CFLAGS) -c $?

exit:


#convergence: $(STDOBJ) convergence.o
#	$(COMP) $(LFLAGS) -o convergence $(STDOBJ) convergence.o

#plotdata: $(STDOBJ) plotdata.o
#	$(COMP) $(LFLAGS) -o plotdata $(STDOBJ) plotdata.o

#plotspectrum: $(STDOBJ) plotspectrum.o
#	$(COMP) $(LFLAGS) -o plotspectrum $(STDOBJ) plotspectrum.o

plotsignal: $(STDOBJ) plotsignal.o
	$(COMP) $(LFLAGS) -o plotsignal $(STDOBJ) plotsignal.o

plotspins: $(STDOBJ) plotspins.o
	$(COMP) -o plotspins $(LFLAGS) $(STDOBJ) plotspins.o

#plotspins_old: $(STDOBJ) plotspins_old.o
#	$(COMP) $(LFLAGS) -o plotspins_old $(STDOBJ) plotspins_old.o

#plotstats: $(STDOBJ) plotstats.o
#	$(COMP) $(LFLAGS) -o plotstats $(STDOBJ) plotstats.o

clean:
	rm -f *.o *.mod







