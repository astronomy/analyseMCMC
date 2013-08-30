
  Plot routines for SPINspiral output
  
  http://spinspiral.sourceforge.net
  
  
  To compile the plot code, you'll need	pgplot, available in the better Linux
  distributions, or from http://www.astro.caltech.edu/~tjp/pgplot/
  
  
  Compilation:
  
  1) Using a static Makefile:
    This directory (./doc/) contains an example Makefile and an example input file,
    which should reside in the directory where you run analysemcmc (and typically
    where your mcmc files are).
    These files are supposed to be changed as the code changes, you can use them as 
    reference to maintain your own files.
    
  2) Use CMake (experimental):
    The main directory of this package should contain a file called CMakeLists.txt.
    From that directory, do:
      $ mkdir build; cd build
      $ cmake ..
      $ make
    CMake is supposed to be more cross-platform compatible.

    To compile the code with your favourite compiler, prepend the cmake line with e.g. FC=gfortran:
      $ FC=gfortran cmake ..

    
  In both cases you should get a binary executable called "analyseMCMC" in the main
  directory.
  
  
  getanalysemcmcdat:
  
  The script getanalysemcmcdat does exactly that; it checks the default version of
  analysemcmc.dat (make sure you set the variable SRCFILE to wherever that file is
  on your system).  If you don't have a copy of the file, it will get you one.
  If you do, it will compare the two, and, if different, ask whether to keep or
  replace the current file, or to interactively merge it (line by line).
  You can put a copy of getanalysemcmcdat in your local executable directory,
  e.g. ~/usr/bin/
  

