!> \file analyseMCMC_main.f90  AnalyseMCMC main routine

! 
! LICENCE:
! 
! Copyright (c) 2007-2013  Marc van der Sluys, Vivien Raymond, Ben Farr, Chris Chambers
!  
! This file is part of the AnalyseMCMC package.
!  
! This is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published
! by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! 
! This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License along with this code (LICENCE).  If not, see 
! <http://www.gnu.org/licenses/>.
! 


 
!> \mainpage Documentation <a href="http://analysemcmc.sourceforge.net/">AnalyseMCMC</a>
!! <a href="http://analysemcmc.sourceforge.net/">AnalyseMCMC</a> is a Fortran code that can be used to analyse and present the
!! output of <a href="http://spinspiral.sourceforge.net/">SPINspiral</a> and 
!! <a href="https://www.lsc-group.phys.uwm.edu/daswg/projects/lalsuite.html">lalinference_mcmc</a>.
!!
!! 
!! \section input  Input options in settings file analysemcmc.dat
!!
!! Place a file called analysemcmc.dat in the directory where the MCMC output files are, and where you run analyseMCMC.
!! If no settings file is present, the code will dump default output in HTML and png format in a new subdirectory html/
!! The input options in the file, with their default values, are:
!!
!! \subsection basic  Basic options
!!  Basic options:
!! -  \b THIN               =  10,         If >1, "thin" the output; read every thin-th line
!! -  \b NBURNMAX           =  500000,     If >=0: override length of the burn-in phase, for all chains! This is now the ITERATION 
!!                                         number (it becomes the line number later on).  Nburn > Nchain sets Nburn = 0.1*Nchain
!! -  \b NBURNFRAC          =  0.5,        If !=0: override length of the burn-in phase, as a fraction of the length of each 
!!                                         chain. This overrides Nburn above
!! -  \b AUTOBURNIN         =  -1.0,       If !=0: Determine burn-in automatically as the first iteration where log(L_chain) > 
!!                                         max(log(L_allchains)) - autoBurnin. If autoBurnin<0, use Npar/2. Overrides Nburn and 
!!                                         NburnFrac above
!! -  \b MAXCHLEN           =  1000000000, Maximum chain length to read in (number of iterations, not number of lines)
!! 
!! -  \b FILE               =  1,          Plot output to file:  0-no; screen,  >0-yes; 1-png, 2-eps, 3-pdf.  Give an output path 
!!                                         for files in the parameter "outputdir" below
!! -  \b COLOUR             =  1,          Use colours: 0-no (grey scales), 1-yes
!! -  \b QUALITY            =  2,          "Quality" of plot, depending on purpose: 0: draft, 1: paper, 2: talk, 3: poster
!!   
!! -  \b REVERSEREAD        =  0,          Read files reversely (anti-alphabetically), to plot coolest chain last so that it 
!!                                         becomes better visible: 0-no, 1-yes, 2-use colours in reverse order too
!! -  \b UPDATE             =  0,          Update screen plot every 10 seconds: 0-no, 1-yes
!! -  \b MERGECHAINS        =  1,          Merge the data from different files into one chain: 0-no (treat separately), 1-yes 
!!                                         (default)
!! -  \b WRAPDATA           =  1,          Wrap the data for the parameters that are in [0,2pi]: 0-no, 1-yes (useful if the peak 
!!                                         is around 0), 2-mirror inclination from 0-pi into 0-pi/2
!! -  \b CHANGEVAR          =  1,          Change MCMC parameters (e.g. logd->d, kappa->theta_SL, rad->deg)
!! -  \b OUTPUTDIR          =  ".",        Save output files/plots in this directory - this may be a relative or absolute path
!! 
!! 
!! \subsection print  Print options
!!   Select what output to print to screen and write to file:
!! 
!! -  \b PRSTDOUT           =  1,          Print standard output to 1: screen, 2: text file
!! -  \b PRPROGRESS         =  2,          Print general messages about the progress of the program: 0-no, 1-some, 2-more, 3-debug 
!!                                         output
!! -  \b PRRUNINFO          =  1,          Print run info (# iterations, seed, # detectors, SNRs, data length, etc.): 0-no, 1-only 
!!                                         for one file (eg. if all files similar), 2-for all files
!! -  \b PRCHAININFO        =  1,          Print chain info: 1-summary (tot # data points, # contributing chains),  2-details per 
!!                                         chain (file name, plot colour, # iterations, burn-in, Lmax, # data points)
!! -  \b PRINITIAL          =  3,          Print starting values for each chain: 0-no, 1-yes, 2-add injection value, 3-add Lmax 
!!                                         value, 4-add differences (~triples number of output lines for this part)
!!  
!! -  \b PRSTAT             =  1,          Print statistics: 0-no, 1-yes, for default probability interval, 2-yes, for all 
!!                                         probability intervals
!! -  \b PRCORR             =  1,          Print correlations: 0-no, 1-yes
!! -  \b PRACORR            =  0,          Print autocorrelations: 0-no, 1-some, 2-more
!! -  \b NACORR             =  100,        Compute nAcorr steps of autocorrelations if prAcorr>0 or plAcorr>0 (default: 100)
!! -  \b PRIVAL             =  1,          Print interval info: 0-no, 1-for run with injected signal, 2-for run without injection, 
!!                                         3-both
!! -  \b PRCONV             =  1,          Print convergence information for multiple chains to screen and chains plot: 0-no, 1-
!!                                         one summary line, 2-add total chain stdevs, 3-add medians, stdevs for each chain
!!  
!! -  \b SAVESTATS          =  1,          Save statistics (statistics, correlations, intervals) to file: 0-no, 1-yes, 2-yes + 
!!                                         copy in PS
!! -  \b SAVEPDF            =  0,          Save the binned data for 1d and/or 2d pdfs (depending on plPDF1D and plPDF2D).  This 
!!                                         causes all 12 parameters + m1,m2 to be saved and plotted(!), which is slighty annoying
!! -  \b WIKIOUTPUT         =  0,          Save output for the CBC wiki
!! -  \b TAILOREDOUTPUT     =  0,          Save (ascii) output for a specific purpose, e.g. table in a paper:  0-no, tO>0: use the 
!!                                         hardcoded block tO in the subroutine tailoredOutput() in analysemcmc_stats.f90
!! -  \b HTMLOUTPUT         =  0,          Save HTML output - experimental, partly implemented.  Useful when plotting all 2D PDFs 
!!                                         as png
!! 
!! \subsection plot_select  Plot selection
!!   Choose which plots to make:
!! 
!! -  \b PLOT               =  1,          0: plot nothing at all, 1: plot the items selected below
!! 
!! -  \b PLLOGL             =  1,          Plot log L chains: 0-no, 1-yes
!! -  \b PLCHAIN            =  1,          Plot parameter chains: 0-no, 1-yes
!! -  \b PLPARL             =  0,          Plot L vs. parameter value: 0-no, 1-yes
!! -  \b PLJUMP             =  0,          Plot actual jump sizes: 0-no, 1-yes: lin, 2-yes: log
!! -  \b PLACORR            =  0,          Plot autocorrelations: 0-no, 1-yes
!! -  \b PLRHAT             =  0,          Plot Rhat vs. iteration number: 0-no, 1-yes
!! 
!! -  \b PLPDF1D            =  1,          Plot 1d posterior distributions: 0-no, 1-yes: smoothed curve, 2-yes: actual histogram. 
!!                                         If plot=0 and savePDF=1, this determines whether to write the pdfs to file or not.
!! -  \b PLPDF2D            =  2,          Plot 2d posterior distributions: 0-no, 1-yes: gray + contours, 2:gray only, 3: contours 
!!                                         only. If plot=0 and savePDF=1, this determines whether to write the pdfs to file (>0) 
!!                                         or not (=0).
!! -  \b PLOTSKY            =  2,          Plot 2d pdf with stars, implies plPDF2D>0:  0-no, 1-yes, 2-full sky w/o stars, 3-full 
!!                                         sky with stars
!! -  \b MAPPROJECTION      =  1,          Choose map projection: 1-Mollweide
!! 
!! -  \b PLANIM             =  0,          Create movie frames
!! 
!! \subsection plot_options  Plot options
!!  Detailed plot settings:
!! 
!! -  \b CHAINPLI           =  0,          Plot every chainPlI-th point in chains, logL, jump plots:  chainPlI=0: autodetermine, 
!!                                         chainPlI>0: use this chainPlI.  All states in between *are* used for statistics, pdf 
!!                                         generation, etc.
!! -  \b SCLOGLPL           =  1,          Scale logL plot ranges: 0: take everything into account, including burn-in and starting 
!!                                         values;  1: take only post-burn-in and injection values into account
!! -  \b SCCHAINSPL         =  1,          Scale chains plot ranges: 0: take everything into account, including burn-in;  1: take 
!!                                         only post-burn-in and injection values into account
!!  
!! -  \b PLINJECT           =  1,          Plot injection values in the chains and pdfs: 0: no,  1: yes (all pars),  2: yes (
!!                                         selected pars), 3-4: as 1-2 + print value in PDF panel
!! -  \b PLSTART            =  1,          Plot starting values in the chains
!! -  \b PLMEDIAN           =  1,          Plot median values in the pdfs: 1-1D PDFs, 2-2D PDFs, 3-both. 4-6: as 1-3 + write value 
!!                                         in PDF panel
!! -  \b PLRANGE            =  4,          Plot the probability range in the pdfs: 1-1D PDFs, 2-2D PDFs, 3-both. 4-6: as 1-3 + 
!!                                         write value in PDF panel
!! -  \b PLBURN             =  1,          Plot the burn-in in logL, the chains, etc.
!! -  \b PLLMAX             =  0,          Plot the position of the max logL, in the chains and pdfs
!! 
!! -  \b PRVALUES           =  1,          Print values (injection, median, range) in pdfs
!! -  \b SMOOTH             =  3,          Smooth the pdfs: 0 - no, >1: smooth over smooth bins (use ~10 (3-15)?).   This is 1D 
!!                                         only for now, and can introduce artefacts on narrow peaks!
!! -  \b FILLPDF            =  1,          Fillstyle for the pdfs (pgsfs): 1-solid, 2-outline, 3-hatched, 4-cross-hatched
!! -  \b NORMPDF1D          =  1,          Normalise 1D pdfs:  0-no,  1-normalise surface area (default, a must for different bin 
!!                                         sizes),  2-normalise to height,  3-normalise to sqrt(height), nice to compare par.temp. 
!!                                         chains
!! -  \b NORMPDF2D          =  4,          'Normalise' 2D pdfs; greyscale value depends on bin height:  0-linearly,  1-
!!                                         logarithmically,  2-sqrt,  3-weigted with likelihood value,  4-2D probability intervals
!! 
!! -  \b NANIMFRAMES        =  1,          Number of frames for the movie
!! -  \b ANIMSCHEME         =  3,          AnimScheme (1-3): determines what panels to show in a movie frame; see source 
!!                                         code
!! 
!! -  \b NIVAL              =  3,          Number of probability intervals
!! -  \b IVAL0              =  2,          Number of the default probability interval (ival0<=Nival)
!! -  \b IVALS              =  0.68269, 0.9545, 0.9973, 2*0.0,  Probability intervals (ivals()). Values > 0.9999 will be treated 
!!                                         as 100%.  Maximum 5 intervals
!! 
!! \subsection format  Output format
!! Output format for plots:
!! 
!! -  \b SCRSZ              =  10.8,       Screen size for X11 windows (PGPlot units):  MacOS: 16.4, Gentoo: 10.8
!! -  \b SCRRAT             =  0.57,       Screen ratio for X11 windows (PGPlot units), MacBook: 0.57
!! -  \b BMPXSZ             =  1000,       X-size for bitmap (pixels):  1000  !! Too large values give incomplete 2D PDFs somehow !!
!! -  \b BMPYSZ             =   700,       Y-size for bitmap (pixels):  700
!! -  \b PSSZ               =  10.5,       Size for PS/PDF (PGPlot units).  Default: 10.5   \__ Gives same result as without pgpap
!! -  \b PSRAT              =  0.742,      Ratio for PS/PDF (PGPlot units). Default: 0.742  
!! -  \b SCFAC              =  1.2,        !!!Not fully implemented yet!!!  Scale .png plots up by this factor, then down to the x,
!!                                         y size indicated above to interpolate and smoothen the plot
!! -  \b UNSHARP            =  10,         Apply unsharp mask when creating .png plots. Default: 10.
!! 
!! \subsection fonts_symbols  Fonts and symbols
!!   Fonts, symbols, et cetera:
!! 
!! -  \b ORIENTATION        =  1,          Use portrait (1) or landscape (2) for eps/pdf; mainly useful when sending a plot to a 
!!                                         printer
!! -  \b FONTTYPE           =  1,          Font type used for eps/pdf: 1-simple, 2-roman, 3-italic, 4-script
!! -  \b FONTSIZE1D         =  1.0,        Scale the font size for 1D plots (chains, 1D PDFs) with respect to default. Default: 1.0
!! -  \b FONTSIZE2D         =  1.0,        Scale the font size for 2D plots (2D PDFs) with respect to default.  Default: 1.0
!! -  \b CHAINSYMBOL        =  1,          Plot symbol for the chains: 0-plot lines, !=0: plot symbols: eg: 1: dot (default), 2: 
!!                                         plus, etc.  -4: filled diamond, 16,17: filled square,circle 20: small open circle; -10/-
!!                                         11: use a selection of open/filled symbols
!! 
!! \subsection plot_parameters_binning  Plot, bin parameters
!!   Select parameters to plot, bin, et cetera:
!!
!! -  \b NPLPAR             =  14,         Number of plot parameters for 1D PDFs (and chain, jump plots, max 22).  Set to -1 to 
!!                                         plot all available parameters.  This is ignored if savePDF=1.
!! -  \b PLPARS             =  61,  62,  22,  11,  41, 31,  32,  51,  52,  13*0,    MCMC parameters to plot: 11-2:tc/t40, 21-2:d^3/
!!                                         logd, 31-2:RA/dec, 41:phase, 51-4:i/psi/thJo/phJo, 61-4:Mc/eta/M1/M2, 71-3:a/th/phi for 
!!                                         spin1, 81-3: spin2
!! -  \b PANELS             =  2*0,        Number of for 1D plots in x,y direction:  0: autodetermine
!! 
!! -  \b NBIN1D             =  0,          Number of bins for 1D PDFs:  0: autodetermine
!! -  \b NBIN2DX            =  0,          Number of bins in x-direction for 2D PDFs and 2D probability ranges:  0: autodetermine 
!!                                         (for both x and y)
!! -  \b NBIN2DY            =  0,          Number of bins in y-direction for 2D PDFs and 2D probability ranges:  0: use Nbin2Dx, 
!!                                         -1: use Nbin2Dx*(scr/bmp/ps)rat
!! 
!! -  \b NPDF2D             =  3,          Number of 2D-PDF plots to make:  -1: all plots (91 for 12+2 parameters),  >0: read 
!!                                         parameters from the lines below
!! -  \b PDF2DPAIRS         =  61,  31,  52, 47*0,  62,  32,  51, 47*0,   Parameter pairs to plot a 2D PDF for.  Maximum 50 pairs.
!!
!! 
!!
!! \section IO I/O units used
!! 
!! - \b  0: StdErr
!! - \b  6: StdOut
!! 
!! - \b 10: input (MCMC) file
!! - \b 16: temp file in getos
!! - \b 17: temp file in timestamp
!! - \b 19: StdOut redirection (__output.txt)
!! 
!! - \b 20: Output: __statistics.dat, __wiki.dat, __bayes.dat 
!! 
!! - \b 21: bsc.dat (bright star catalogue)
!! - \b 22: bsc_const.dat (BSC constellation data)
!! - \b 23: bsc_names.dat (BSC names)
!! - \b 24: milkyway*.dat (Milky Way data)
!!  
!! - \b 30: Output: __pdf1/2d.dat
!! - \b 40: Output: tailored output
!!  
!! - \b 51: Output: HTML PDF2D matrix
!! 




!***********************************************************************************************************************************
!> \brief Main routine

program analyseMCMC
  !use SUFR_version, only: print_libSUFR_version
  use SUFR_kinds, only: double,dbl
  use SUFR_constants, only: stdOut,stdErr, set_SUFR_constants
  use SUFR_random_numbers, only: get_ran_seed
  use SUFR_system, only: quit_program_error
  
  use aM_constants, only: stdOutFile, use_PLplot
  use analysemcmc_settings, only: settingsfile, panels,htmlOutput,prProgress,file,colour,prStdOut,prChainInfo
  use analysemcmc_settings, only: prCorr,saveStats,plot,plLogL,plChain,plPDF1D,plPDF2D
  use analysemcmc_settings, only: plAnim,bmpXSz,bmpYSz,Npdf2D,reverseRead
  use analysemcmc_settings, only: whiteBG,scFac,scrSz,scrRat,PSsz,PSrat,unSharp,orientation,chainSymbol,quality,plJump,savePDF
  use analysemcmc_settings, only: wrapData,update,plPars,nPlPar,mergeChains,tailoredOutput,plACorr,plRhat, maxChs
  !use analysemcmc_settings, only: phi_q_sorting
  use general_data, only: infiles,allDat,selDat,post,prior,outputDir,nchains0,nchains,ntot,outputname
  use mcmcrun_data, only: nMCMCpar,parID
  use plot_data, only: colours,symbols,colournames,maxdots,bmpsz,bmprat,bmpxpix,pltsz,pltrat,unSharplogl,unSharpchain,unSharppdf1d
  use plot_data, only: unSharppdf2d,psclr,ncolours,nsymbols,defcolour
  
  implicit none
  integer :: i,ic,io,exitcode,tempintarray(99),status,system, tmpStdOut
  real(double) :: timestamp,timestamps(9)  ! Time the progress of the code.
  character :: infile*(99), finalOutputName*(199)
  logical :: ex,timing
  
  
  call set_SUFR_constants()   ! Define constants in libSUFR
  call setconstants()         ! Define mathematical constants
  
  timestamps = 0.0_dbl
  timestamps(1) = timestamp()
  timing = .false.
  if(abs(timestamps(1)).gt.1.e-6_dbl .and. abs(timestamps(1)).lt.1.e20_dbl) timing = .true.
  
  
  
  ! Get command-line arguments:
  nchains0 = command_argument_count()
  
  settingsfile = 'analysemcmc.dat'
  if(nchains0.eq.1) then  ! Check whether this is the name of the input file
     call get_command_argument(1,infile)  ! Read file name from the command-line arguments
     if(index(trim(infile),'nalyse').ne.0 .and. infile(len_trim(infile)-3:len_trim(infile)).eq.'.dat') then
        settingsfile = trim(infile)
        nchains0 = 0
     end if
  end if
  
  ! Set and read settings:
  call set_plotsettings()     ! Set plot settings to 'default' values
  call read_settingsfile()    ! Read the plot settings (overwrite the defaults)
  if(prProgress.ge.3) call write_settingsfile()   ! Write the input file back to disc as analysemcmc.new
  
  ! New parameters that should go into the settings file(?):
  !phi_q_sorting = 0  ! Do phase/mass-ratio sorting (if phi>pi, q -> 1/q; m1 <-> m2): 0-no, 1-yes - not implemented yet
  
  
  if(nchains0.lt.1) then  ! No command-line arguments - select all files PTMCMC.output.*.00 in the current dir
     call findFiles('PTMCMC.output.*.00',maxChs,1,infiles,nchains0)
     if(nchains0.eq.0) then
        if(prProgress.ge.2) then
           write(stdErr,'(A)')'  No files matching  PTMCMC.output.*.00  were found in the current directory.'
           write(stdErr,'(A)')'  I will try SPINspiral output file names  SPINspiral.output.*.00  instead.'
        end if
        call findFiles('SPINspiral.output.*.00',maxChs,1,infiles,nchains0)
     end if
     if(nchains0.eq.0) then
        if(prProgress.ge.2) then
           write(stdErr,'(A)')'  No files matching  SPINspiral.output.*.00  were found either.'
           write(stdErr,'(A)')'  I will try the ancient file names  mcmc.output.*.00  before I give up.'
        end if
        call findFiles('mcmc.output.*.00',maxChs,1,infiles,nchains0)
        if(nchains0.eq.0) call quit_program_error('No valid input files were found in the current directory.'// &
             '  Please specify input files manually.',1)
     end if
  else
     do ic = 1,nchains0
        if(reverseRead.eq.0) then
           call get_command_argument(ic,infile)  ! Read file name from the command-line arguments
        else
           call get_command_argument(nchains0-ic+1,infile)  ! Read file name from the command-line arguments in reverse order
        end if
        infiles(ic) = infile
     end do
  end if
  
  
  if(nchains0.gt.maxChs) write(stdErr,'(A,I3,A)')'  *** WARNING:  Too many input files (chains),'// &
       ' please increase maxChs in analyseMCMC_modules.f90. Only',maxChs,' files can be read.'
  nchains0 = min(nchains0,maxChs)
  
  
  ! Print code version and set use_PLplot:
  if(prProgress.ge.1) then
     call print_code_version(stdOut, use_PLplot)
     write(StdOut,'(A)') '  Using settings file:  '//trim(settingsfile)
  end if
  
  
  if(htmlOutput.ge.1) then
     outputdir = 'html'       ! Directory where output is saved (either relative or absolute path)
     
     prStdOut=2
     file = 1

     prCorr = 1
     saveStats = 1
     savePDF = 0  ! Since this prevents the folding of the data (e.g. around 2pi) for the PDFs
     
     plot = 1
     plLogL = 1
     plChain = 1
     plPDF1D = 1
     plPDF2D = 2
     nPlPar = -1  ! Plot all
     Npdf2D = -1  ! Plot all
     plAnim = 0
     
     bmpXSz = 1000
     bmpYSz =  700
  end if
  
  
  
  inquire(file=trim(outputdir), exist=ex)  ! Check whether the directory already exists 
  if(.not.ex) then
     status = system('mkdir -p '//trim(outputdir))
     if(status.ne.0) call quit_program_error('Could not create output directory: '//trim(outputdir)//'/',1)
  end if
  
  ! Write at least the code version to file:
  tmpStdOut = 19
  write(stdOutFile,'(A,I6.6,A)') trim(outputdir)//'/analysemcmc_tempstdout_',abs(get_ran_seed(0)),'.txt'
  open(unit=tmpStdOut,action='write',form='formatted',status='replace',file=trim(stdOutFile),iostat=io)
  if(io.ne.0) call quit_program_error('Error opening output file '//trim(stdOutFile),1)
  
  if(htmlOutput.ge.1) write(tmpStdOut,'(A)') '<html><body><pre><b>'
  call print_code_version(tmpStdOut, use_PLplot)
  if(htmlOutput.ge.1) write(tmpStdOut,'(A)') '</b><br>'
  if(prProgress.ge.1) write(tmpStdOut,'(A)') '  Using settings file:  '//trim(settingsfile)
  
  
  if(prStdOut.ge.2) then
     stdOut = tmpStdOut  ! Keep writing standard output to file rather than screen
  !else
  !   close(tmpStdOut)    ! Close file and return to screen output
  end if
  write(stdOut,*)
  
  
  
  
  
  
  !Some of the stuff below will have to go to the input file
  
  ! ~Maximum number of dots to plot in e.g. chains plot, to prevent dots from being overplotted too much 
  !  and eps/pdf files from becoming huge.  Use this to autoset chainPlI
  maxdots = 25000  
  
  ! Use a white background in screen and bitmap plots: 0-no (black), 1-yes.  Used to be in input file.
  whiteBG = 1                 
  
  ! Determine plot sizes and ratios:   (ratio ~ y/x and usually < 1 ('landscape'))
  call compBitmapSize(bmpXSz,bmpYSz, scFac, bmpsz,bmprat)  
  
  write(bmpxpix,'(I4)')bmpXSz  !Used as a text string by convert
  if(file.eq.0) pltsz = scrSz
  if(file.eq.0) pltrat = scrRat
  if(file.eq.1) pltsz = bmpsz
  if(file.eq.1) pltrat = bmprat
  if(file.ge.2) pltsz = PSsz
  if(file.ge.2) pltrat = PSrat
  
  
  ! Use full unsharp-mask strength for plots with many panels and dots, weaker for those with fewer panels and/or no dots:
  write(unSharplogl,'(I4)')  max(nint(real(unSharp)/2.),1)  ! Only one panel with dots
  write(unSharpchain,'(I4)') unSharp                        ! ~15 panels with dots
  write(unSharppdf1d,'(I4)') max(nint(real(unSharp)/2.),1)  ! ~15 panels, no dots
  write(unSharppdf2d,'(I4)') max(nint(real(unSharp)/4.),1)  ! 1 panel, no dots
  
  
  
  
  
  
  
  
  ! Sort out implicit options:
  ! eps/pdf: colour and orientation
  psclr = '/cps'
  if(colour.eq.0) psclr = '/ps'
  if(orientation.eq.1) then
     psclr = '/vcps'
     if(colour.eq.0) psclr = '/vps'
  end if
  
  ncolours =  5;  colours(1:ncolours) = (/4,2,3,6,5/)              ! Paper
  ncolours = 10;  colours(1:ncolours) = (/2,3,4,5,6,7,8,9,10,11/)
  if(colour.eq.1.and.quality.eq.2) then  ! Beamer
     ncolours = 5
     colours(1:ncolours)=(/4,2,5,11,15/)
  end if
  if(colour.ne.1) then
     ncolours=3
     colours(1:ncolours)=(/1,14,15/)
     !ncolours=6
     !colours(1:ncolours)=(/1,1,1,15,1,15/)
     !print*,chainSymbol,nsymbols
  end if
  if(colour.eq.1.and.quality.eq.0.and.maxChs.gt.5) then
     ncolours = 10
     colours(1:ncolours)=(/2,3,4,5,6,7,8,9,10,11/)
  end if
  !Overrule
  !ncolours = 1
  !colours(1:ncolours)=(/6/)
  !defcolour = 2 !Red e.g. in case of 1 chain
  defcolour = colours(1)
  
  nsymbols =  1;  symbols(1:nsymbols) = (/chainSymbol/)
  if(chainSymbol.eq.-10) then
     nsymbols = 8
     symbols(1:nsymbols) = (/2,4,5,6,7,11,12,15/) !Thin/open symbols
  end if
  if(chainSymbol.eq.-11) then
     nsymbols = 6
     symbols(1:nsymbols) = (/-3,-4,16,17,18,-6/) !Filled symbols
  end if
  
  if(reverseRead.ge.2) then !Reverse colours too
     do i=1,ncolours
        tempintarray(i) = colours(i)
     end do
     do i=1,nchains0
        colours(i) = tempintarray(nchains0-i+1) !Reverse colours too, but use the same first nchains0 from the set
     end do
  end if
  
  if(plot.eq.0) then
     plLogL = 0
     plChain = 0
     plJump = 0
     if(savePDF.eq.0) then
        plPDF1D = 0
        plPDF2D = 0
     end if
     plAnim = 0
  end if
  if(savePDF.eq.1) then
     !if(nPlPar.ne.15) write(stdErr,'(/,A)')'*** WARNING:  I changed nPlPar to 15, since savePDF is selected ***'
     !nPlPar = 15; plPars(1:nPlPar) = (/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15/) !All 12 + m1,m2
     if(wrapData.ne.0) write(stdErr,'(A)')'  * Warning:  I found that savePDF = 1, so I set wrapData to 0.'
     wrapData = 0
  end if
  if(file.ge.1) update = 0
  if(plAnim.ge.1) update = 0
  
  colournames(1:15) = (/ &
       'white               ','red                 ','dark green          ','dark blue           ','cyan                ', &
       'magenta             ','yellow              ','orange              ','light green         ','brown               ', &
       'dark red            ','purple              ','red-purple          ','dark grey           ','light grey          '/)
  if(file.ge.2) colournames(1) = 'black'
  
  call set_originalParameterNames()  ! Set the names and symbols of the original MCMC parameters in the database
  
  
  if(prChainInfo.ge.1) then
     call print_rundata(stdOut)
     if(prStdOut.lt.2) call print_rundata(tmpStdOut)
  end if
  nchains = nchains0
  
  
  
  
  
  
  
  !*********************************************************************************************************************************
  !***   READ INPUT FILE(S)   ******************************************************************************************************
  !*********************************************************************************************************************************
  
101 continue
  ! Read the input files:
  exitcode = 0
  call read_mcmcfiles(exitcode)
  if(exitcode.ne.0) goto 9998
  
  ! Create HTML file:
  if(htmlOutput.ge.1) call create_html_index_file(stdOut)
  
  ! Get and print some basic chain statistics:
  if(timing) timestamps(2) = timestamp()
  exitcode = 0
  call mcmcruninfo(exitcode)
  
  ! Plot all parameters if nPlPar = -1:
  if(nPlPar.eq.-1) then
     nPlPar = nMCMCpar
     PlPars(1:nPlPar) = parID(1:nMCMCpar)
  end if
  
  if(savePDF.ge.2) then
     write(stdOut,'(A)')'  Writing after-burnin data points to file'
     call save_data(exitcode)  ! save after-burnin combined data to file
  end if
  
  
  ! More implicit options:
  if(panels(1)*panels(2).lt.min(nPlPar,nMCMCpar)) panels = 0
  if(panels(1)*panels(2).lt.1) then
     if(min(nPlPar,nMCMCpar).eq.1) panels = (/1,1/)
     if(min(nPlPar,nMCMCpar).eq.2) panels = (/2,1/)
     if(min(nPlPar,nMCMCpar).eq.3) panels = (/3,1/)
     if(min(nPlPar,nMCMCpar).eq.4) panels = (/2,2/)
     if(min(nPlPar,nMCMCpar).eq.5) panels = (/5,1/)
     if(min(nPlPar,nMCMCpar).eq.6) panels = (/3,2/)
     if(min(nPlPar,nMCMCpar).eq.7) panels = (/4,2/)
     if(min(nPlPar,nMCMCpar).eq.8) panels = (/4,2/)
     if(min(nPlPar,nMCMCpar).eq.9) panels = (/3,3/)
     if(min(nPlPar,nMCMCpar).eq.10) panels = (/5,2/)
     if(min(nPlPar,nMCMCpar).eq.11) panels = (/4,3/)
     if(min(nPlPar,nMCMCpar).eq.12) panels = (/4,3/)
     if(min(nPlPar,nMCMCpar).eq.12.and.quality.eq.3) panels = (/3,4/)
     if(min(nPlPar,nMCMCpar).eq.13) panels = (/5,3/)
     if(min(nPlPar,nMCMCpar).eq.14) panels = (/5,3/)
     if(min(nPlPar,nMCMCpar).eq.15) panels = (/5,3/)
     if(min(nPlPar,nMCMCpar).eq.16) panels = (/4,4/)
     if(min(nPlPar,nMCMCpar).eq.17) panels = (/6,3/)
     if(min(nPlPar,nMCMCpar).eq.18) panels = (/6,3/)
     if(min(nPlPar,nMCMCpar).eq.19) panels = (/5,4/)
     if(min(nPlPar,nMCMCpar).eq.20) panels = (/5,4/)
  end if
  
  
  
  
  ! ********************************************************************************************************************************
  ! ***  DO STATISTICS   ***********************************************************************************************************
  ! ********************************************************************************************************************************
  
  if(timing) timestamps(3) = timestamp()
  
  exitcode = 0
  call statistics(exitcode)
  if(exitcode.ne.0) goto 9999
  
  
  
  
  
  ! ********************************************************************************************************************************
  ! ***  CREATE PLOTS   ************************************************************************************************************
  ! ********************************************************************************************************************************
  
  if(timing) timestamps(4) = timestamp()
  
  if(prProgress.ge.2) write(stdOut,*)
  if(htmlOutput.ge.1) then
     write(stdOut,'(A)') '<br><hr><a name="plots"></a><font size="1"><a href="#top" title="Go to the top of the page">top</a>'// &
          '</font><h2>Plots</h2>'
  else
     if(plot.eq.1.and.prProgress.ge.1.and.update.eq.0) then
        write(stdOut,'(/,A)',advance="no")'  Plotting '
        if(file.eq.0) write(stdOut,'(A)',advance="no")'to screen: '
        if(file.eq.1) write(stdOut,'(A)',advance="no")'to png: '
        if(file.eq.2) write(stdOut,'(A)',advance="no")'to eps: '
        if(file.eq.3) write(stdOut,'(A)',advance="no")'to pdf: '
     end if
  end if
  
  
  !*********************************************************************************************************************************
  ! Plot (1d) chains: logL, parameter chains, jumps, etc.
  if(plot.eq.1) then
     exitcode = 0
     call chains(exitcode)
     if(exitcode.ne.0) goto 9999
  end if
  if(timing) timestamps(5) = timestamp()
  
  
  
  
  !*********************************************************************************************************************************
  ! Plot pdfs (1d)
  if(plPDF1D.ge.1) then
     exitcode = 0
     call pdfs1d(exitcode)
     if(exitcode.ne.0) goto 9999
  end if !if(plPDF1D.ge.1)
  if(timing) timestamps(6) = timestamp()
  
  
  
  
  !*********************************************************************************************************************************
  if(plPDF2D.ge.1.and.mergeChains.eq.0) then
     write(stdOut,'(A)',advance="no")', (skipping 2D PDFs since mergeChains=0), '
     plPDF2D = 0
  end if
  
  if(plPDF2D.ge.1) then
     exitcode = 0
     call pdfs2d(exitcode)
     if(exitcode.ne.0) goto 9999
  end if !if(plPDF2D.eq.1)
  
  if(Npdf2D.lt.0) then !Then we just plotted all 2D PDFs
     write(stdOut,*)
  else
     if(htmlOutput.eq.0.and.prProgress.ge.1.and.update.eq.0.and.plot.gt.0) write(stdOut,'(A,/)') ' done.  '
  end if
  
  if(timing) timestamps(7) = timestamp()  
  
  
  
  if(htmlOutput.ge.1) write(stdOut,'(A)') '<br><hr><br>'
  
  
  !*********************************************************************************************************************************
  
  if(saveStats.ge.1.and.nchains.gt.1) then
     write(stdErr,'(A)')' ******   Cannot write statistics if the number of chains is greater than one   ******'
     
     ! Write Bayes factors to file:
     exitcode = 0
     call save_bayes(exitcode)
     if(exitcode.ne.0) goto 9999
  end if
  
  ! Write statistics to file:
  if(saveStats.ge.1.and.nchains.eq.1) then
     exitcode = 0
     call save_stats(exitcode)
     if(exitcode.ne.0) goto 9999
     write(stdOut,*)
  end if !if(saveStats.ge.1.and.nchains.eq.1) then
  
  
  ! Save tailored output to a file:
  if(tailoredOutput.gt.0) then
     exitcode = 0
     call tailored_output(exitcode)
     !if(exitcode.ne.0) goto 9999
  end if
  
  !*********************************************************************************************************************************
  
  if(timing) timestamps(8) = timestamp()
  
  if(plAnim.ge.1) then
     exitcode = 0
     call animation(exitcode)
     if(exitcode.ne.0) goto 9999
  end if
  
  
  
  if(update.eq.1) then
     deallocate(allDat,selDat,post,prior)
     call sleep(5)
     if(sum(ntot).gt.nint(1.e4)) call sleep(5)
     if(sum(ntot).gt.nint(1.e5)) call sleep(10)
     if(sum(ntot).gt.nint(1.e6)) call sleep(20)
     goto 101
  end if
  
9999 continue
  deallocate(selDat)
  
9998 continue
  deallocate(allDat,post,prior)
  
  
  ! Print HTML footer:
  if(htmlOutput.ge.1) then
     write(stdOut,'(A)', advance='no') '<b>'
     call print_code_version(stdOut, use_PLplot)
     write(stdOut,'(A)') '</b>'
     write(stdOut,'(A)') '  Used settings file:  '//trim(settingsfile)
     
     !write(stdOut,'(A)') '<br><br>'
     call print_rundata(stdOut)
  end if
  
  
  ! Print run times:
  if(timing) then
     timestamps(9) = timestamp()
     
     if(prProgress.ge.1) then
        if(htmlOutput.ge.1) write(stdOut,'(A)') '<b>'
        write(stdOut,'(A)',advance="no")'  Run time: '
        write(stdOut,'(A,F5.1,A)',advance="no")'   input:',min(abs(timestamps(2)-timestamps(1)),999.9_dbl),'s,'
        write(stdOut,'(A,F5.1,A)',advance="no")'   stats:',min(abs(timestamps(4)-timestamps(2)),999.9_dbl),'s,'
        if(plot.eq.1.and.plLogL+plChain+plJump+plACorr+plRhat.gt.0) then
           write(stdOut,'(A,F5.1,A)',advance="no")'   chains:',min(abs(timestamps(5)-timestamps(4)),999.9_dbl),'s,'
        end if
        if(plot.eq.1.or.savePDF.ge.1) then
           if(plPDF1D.ge.1) write(stdOut,'(A,F5.1,A)',advance="no") &
                '   1d pdfs:',min(abs(timestamps(6)-timestamps(5)),999.9_dbl),'s,'
           if(plPDF2D.ge.1) write(stdOut,'(A,F6.1,A)',advance="no") &
                '   2d pdfs:',min(abs(timestamps(7)-timestamps(6)),999.9_dbl),'s,'
        end if
        if(plAnim.ge.1) write(stdOut,'(A,F5.1,A)',advance="no")'   movie:',min(abs(timestamps(9)-timestamps(8)),999.9_dbl),'s,'
        write(stdOut,'(A,F6.1,A)')'   total:',min(abs(timestamps(9)-timestamps(1)),999.9_dbl),'s.'
        if(htmlOutput.ge.1) write(stdOut,'(A)') '</b>'
     end if
  end if
  
  
  if(htmlOutput.ge.1) then
     write(stdOut,'(A)') '    </pre>'
     write(stdOut,'(A)') '  </body>'
     write(stdOut,'(A)') '</html>'
     close(51)
  end if
  
  write(stdOut,*)
  
  
  ! Save standard-output file under final name:
  if(prStdOut.ge.2) then
     close(stdOut)
  else
     close(tmpStdOut)
  end if
  
  if(htmlOutput.ge.1) then
     write(finalOutputName,'(A)') trim(outputdir)//'/index.html'
  else
     write(finalOutputName,'(A)') trim(outputdir)//'/'//trim(outputname)//'__output.txt'
  end if
  
  status = system('mv -f '//trim(stdOutFile)//' '//trim(finalOutputName))
  if(status.eq.0) then
     if(prStdOut.ge.2) &  ! Should be 6, not stdOut:
          write(6,'(/,A,/)')'  AnalyseMCMC:  saved standard output to '//trim(finalOutputName)
  else
     write(stdErr,'(/,A)')'  AnalyseMCMC:  Error saving standard output to '//trim(finalOutputName)
     status = system('rm -f '//trim(stdOutFile))
     write(stdErr,*)
  end if
  
end program analyseMCMC
!***********************************************************************************************************************************



