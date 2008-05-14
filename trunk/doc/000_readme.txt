
  Plot routines for MCMC output
  
  
  To compile the plot code, you'll need	pgplot.
  
  This directory (./doc/) contains an example Makefile and an example input file,
  which should reside in the directory where you run plotspins (and typically
  where your mcmc files are).
  These are supposed to be changed as the code changes, you can use them as reference
  to maintain your own files.
  
  The script getplotspinsdat does exactly that; it check the default version of
  plotspins.dat (make sure you set the variable SRCFILE to wherever that file is
  on your system).  If you don't have a copy of the file, it will get you one.
  If you do, it will compare the two, and, if different, ask whether to keep or
  replace the current file, or to interactively merge it (line by line).
  Put a copy of getplotspindat in ~/bin/
  

