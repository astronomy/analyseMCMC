!> \file analyseMCMC_functions.f90  General routines and functions for analyseMCMC

! 
! LICENCE:
! 
! Copyright 2007-2012 Marc van der Sluys
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



!***********************************************************************************************************************************
!> \brief  Define the constants

subroutine setconstants()
  use SUFR_constants, only: stdErr,stdOut
  use aM_constants, only: detabbrs,waveforms
  
  implicit none
  
  stdErr = 0  ! Standard error unit
  stdOut = 6  ! Standard output unit - screen
  
  ! Define detector abbreviations here (don't forget to change detabbrs() in the module constants in analyseMCMC_modules.f90):
  detabbrs = (/'H1','L1','V ','H2'/)
  
  ! Define waveforms here (don't forget to change waveforms() in the module constants in analyseMCMC_modules.f90):
  waveforms    = 'Unknown'             ! All of them, except those specified below
  waveforms(1) = 'Apostolatos'
  waveforms(2) = 'SpinTaylor12'
  waveforms(3) = 'SpinTaylor15'
  waveforms(4) = 'PPN'
  waveforms(5) = 'PhenSpinInspiralRD'
  waveforms(9) = 'Ana.L'
  
end subroutine setconstants
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Read the settings file (called analysemcmc.dat by default)

subroutine read_settingsfile()
  use SUFR_kinds, only: double
  use SUFR_constants, only: stdErr
  use SUFR_system, only: quit_program
  
  use analysemcmc_settings, only: Nburn,ivals,plPars,panels,PDF2Dpairs,thin,NburnFrac,autoBurnin,maxChs,maxChLen,file,colour
  use analysemcmc_settings, only: quality,reverseRead,update,mergeChains,wrapData,changeVar,prStdOut,prProgress,prRunInfo
  use analysemcmc_settings, only: prChainInfo,prInitial,prStat,prCorr,prAcorr,nAcorr,prIval,prConv,saveStats,savePDF,tailoredOutput
  use analysemcmc_settings, only: plot,plLogL,plChain,plParL,plJump,plPDF1D,plPDF2D,plAcorr,plotSky,plAnim,chainPlI,scLogLpl
  use analysemcmc_settings, only: scChainsPl,plInject,plStart,plMedian,plRange,plBurn,plLmax,prValues,smooth,fillPDF,normPDF1D
  use analysemcmc_settings, only: normPDF2D,nAnimFrames,animScheme,Nival,ival0,scrSz,scrRat,bmpXSz,bmpYSz,PSsz,PSrat,scFac,unSharp
  use analysemcmc_settings, only: orientation,fontType,fontSize1D,fontSize2D,chainSymbol,nPlPar,Nbin1D,Nbin2Dx,Nbin2Dy,Npdf2D
  
  implicit none
  integer :: i,u,io,io1
  character :: bla,filename*(99)
  real(double) :: dblvar
  filename = 'analysemcmc.dat'
  
  ! dblvar is used when a (possibly) large integer is expected; read it as double, then convert to integer
  
  u = 15
  open(unit=u,form='formatted',status='old',action='read',file=trim(filename),iostat=io)
  if(io.ne.0) call quit_program('Error opening input file '//trim(filename))
  
  io = 0
  io1 = 0
  
  read(u,*,iostat=io) bla
  bla = bla  ! Remove 'set but never used' warning
  
  
  read(u,*,iostat=io) bla
  read(u,*,iostat=io) thin
  read(u,*,iostat=io) dblvar
  Nburn(1) = nint(dblvar)
  do i=2,maxChs
     Nburn(i) = Nburn(1)
  end do
  read(u,*,iostat=io) NburnFrac
  read(u,*,iostat=io) autoBurnin
  read(u,*,iostat=io) dblvar
  maxChLen = nint(dblvar)
  read(u,*,iostat=io) file
  read(u,*,iostat=io) colour
  read(u,*,iostat=io) quality
  read(u,*,iostat=io) reverseRead
  read(u,*,iostat=io) update
  read(u,*,iostat=io) mergeChains
  read(u,*,iostat=io) wrapData
  read(u,*,iostat=io) changeVar
  
  read(u,*,iostat=io) bla
  read(u,*,iostat=io) prStdOut
  read(u,*,iostat=io) prProgress
  read(u,*,iostat=io) prRunInfo
  read(u,*,iostat=io) prChainInfo
  read(u,*,iostat=io) prInitial
  read(u,*,iostat=io) prStat
  read(u,*,iostat=io) prCorr
  read(u,*,iostat=io) prAcorr
  read(u,*,iostat=io) dblvar
  nAcorr = nint(dblvar)
  read(u,*,iostat=io) prIval
  read(u,*,iostat=io) prConv
  read(u,*,iostat=io) saveStats
  read(u,*,iostat=io) savePDF
  read(u,*,iostat=io) tailoredOutput
  
  
  read(u,*,iostat=io) bla
  read(u,*,iostat=io) plot
  read(u,*,iostat=io) plLogL
  read(u,*,iostat=io) plChain
  read(u,*,iostat=io) plParL
  read(u,*,iostat=io) plJump
  read(u,*,iostat=io) plPDF1D
  read(u,*,iostat=io) plPDF2D
  read(u,*,iostat=io) plAcorr
  read(u,*,iostat=io) plotSky
  read(u,*,iostat=io) plAnim
  
  
  read(u,*,iostat=io) bla
  read(u,*,iostat=io) chainPlI
  read(u,*,iostat=io) scLogLpl
  read(u,*,iostat=io) scChainsPl
  read(u,*,iostat=io) plInject
  read(u,*,iostat=io) plStart
  read(u,*,iostat=io) plMedian
  read(u,*,iostat=io) plRange
  read(u,*,iostat=io) plBurn
  read(u,*,iostat=io) plLmax
  read(u,*,iostat=io) prValues
  read(u,*,iostat=io) smooth
  read(u,*,iostat=io) fillPDF
  read(u,*,iostat=io) normPDF1D
  read(u,*,iostat=io) normPDF2D
  read(u,*,iostat=io) nAnimFrames
  read(u,*,iostat=io) animScheme
  read(u,*,iostat=io) Nival,ival0
  read(u,*,iostat=io1)(ivals(i),i=1,Nival)
  
  read(u,*,iostat=io) bla
  read(u,*,iostat=io) scrSz
  read(u,*,iostat=io) scrRat
  read(u,*,iostat=io) bmpXSz
  read(u,*,iostat=io) bmpYSz
  read(u,*,iostat=io) PSsz
  read(u,*,iostat=io) PSrat
  read(u,*,iostat=io) scFac
  read(u,*,iostat=io) unSharp
  
  read(u,*,iostat=io) bla
  read(u,*,iostat=io) orientation
  read(u,*,iostat=io) fontType
  read(u,*,iostat=io) fontSize1D
  read(u,*,iostat=io) fontSize2D
  read(u,*,iostat=io) chainSymbol
  
  read(u,*,iostat=io) bla
  read(u,*,iostat=io) bla
  read(u,*,iostat=io) nPlPar
  read(u,*,iostat=io1)(plPars(i),i=1,nPlPar)
  if(io1.ne.0) nPlPar = i-1
  io1 = 0
  read(u,*,iostat=io) panels(1:2)
  read(u,*,iostat=io) Nbin1D
  read(u,*,iostat=io) Nbin2Dx
  read(u,*,iostat=io) Nbin2Dy
  read(u,*,iostat=io) Npdf2D
  do i=1,Npdf2D
     read(u,*,iostat=io1)PDF2Dpairs(i,1:2)
     if(io1.ne.0) exit
  end do
  if(io1.ne.0) Npdf2D = i-1
  close(u)
  
  if(io.ne.0) then
     write(stdErr,'(/,A,I2,A,I3,A,/)')'  Error reading input file '//trim(filename)//', aborting...'
     stop
  end if
  
end subroutine read_settingsfile
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Write a copy of the settings file (called analysemcmc.used by default)

subroutine write_settingsfile()
  use analysemcmc_settings, only: Nburn,ivals,plPars,panels,PDF2Dpairs,thin,NburnFrac,autoBurnin,maxChs,maxChLen,file,colour
  use analysemcmc_settings, only: quality,reverseRead,update,mergeChains,wrapData,changeVar,prStdOut,prProgress,prRunInfo
  use analysemcmc_settings, only: prChainInfo,prInitial,prStat,prCorr,prAcorr,nAcorr,prIval,prConv,saveStats,savePDF,tailoredOutput
  use analysemcmc_settings, only: plot,plLogL,plChain,plParL,plJump,plPDF1D,plPDF2D,plAcorr,plotSky,plAnim,chainPlI,scLogLpl
  use analysemcmc_settings, only: scChainsPl,plInject,plStart,plMedian,plRange,plBurn,plLmax,prValues,smooth,fillPDF,normPDF1D
  use analysemcmc_settings, only: normPDF2D,nAnimFrames,animScheme,Nival,ival0,scrSz,scrRat,bmpXSz,bmpYSz,PSsz,PSrat,scFac,unSharp
  use analysemcmc_settings, only: orientation,fontType,fontSize1D,fontSize2D,chainSymbol,nPlPar,Nbin1D,Nbin2Dx,Nbin2Dy,Npdf2D
  
  implicit none
  integer :: u,i
  
  u = 14
  open(unit=u,form='formatted',status='replace',action='write',file='analysemcmc.used')
  
11 format(I10,1x,A19,5x,A)
12 format(2I5,1x,A19,5x,A)
21 format(F10.5,1x,A19,5x,A)
31 format(ES10.2,1x,A19,5x,A)
  
  write(u,'(A,/)')' Input file for AnalyseMCMC'
  
  
  write(u,'(/,A)')' Basic options:'
  write(u,11)thin, 'thin',   'If >1, "thin" the output; read every thin-th line '
  write(u,11)maxval(Nburn), 'Nburn',   'If >=0: override length of the burn-in phase, for all chains!'// &
       ' This is now the ITERATION number (it becomes the line number later on).  Nburn > Nchain sets Nburn = 0.1*Nchain'
  write(u,21)NburnFrac, 'NburnFrac',   'If !=0: override length of the burn-in phase, as a fraction of the length of each chain.'//&
       ' This overrides Nburn above'
  write(u,21)autoBurnin, 'autoBurnin',   'If !=0: Determine burn-in automatically as the first iteration where'// &
       ' log(L_chain) > max(log(L_allchains)) - autoBurnin. If autoBurnin<0, use Npar/2. Overrides Nburn and NburnFrac above'
  write(u,31)dble(maxChLen), 'maxChLen',   'Maximum chain length to read in (number of iterations, not number of lines)'
  write(u,11)file, 'file',   'Plot output to file:  0-no; screen,  >0-yes; 1-png, 2-eps, 3-pdf.  Give an output path for files'// &
       ' in the parameter "outputdir" below'
  write(u,11)colour, 'colour',   'Use colours: 0-no (grey scales), 1-yes'
  write(u,11)quality, 'quality',   '"Quality" of plot, depending on purpose: 0: draft, 1: paper, 2: talk, 3: poster'
  write(u,11)reverseRead, 'reverseRead',   'Read files reversely (anti-alphabetically), to plot coolest chain last so that it'// &
       ' becomes better visible: 0-no, 1-yes, 2-use colours in reverse order too'
  write(u,11)update, 'update',   'Update screen plot every 10 seconds: 0-no, 1-yes'
  write(u,11)mergeChains, 'mergeChains',   'Merge the data from different files into one chain: 0-no (treat separately), 1-yes'// &
       ' (default)'
  write(u,11)wrapData, 'wrapData',   'Wrap the data for the parameters that are in [0,2pi]: 0-no, 1-yes (useful if the peak is'// &
       ' around 0)'
  write(u,11)changeVar, 'changeVar',   'Change MCMC parameters (e.g. logd->d, kappa->theta_SL, rad->deg): 0=no, 1=yes, 2=q->1/q'
  
  
  write(u,'(/,A)')' Select what output to print to screen and write to file:'
  write(u,11)prStdOut, 'prStdOut',   'Print standard output to 1: screen, 2: text file'
  write(u,11)prProgress, 'prProgress',   'Print general messages about the progress of the program: 0-no, 1-some, 2-more,'// &
       ' 3-debug output'
  write(u,11)prRunInfo, 'prRunInfo',   'Print run info (# iterations, seed, # detectors, SNRs, data length, etc.): 0-no, 1-only'// &
       ' for one file (eg. if all files similar), 2-for all files'
  write(u,11)prChainInfo, 'prChainInfo',   'Print chain info: 1-summary (tot # data points, # contributing chains),'// &
       '  2-details per chain (file name, plot colour, # iterations, burn-in, Lmax, # data points)'
  write(u,11)prInitial, 'prInitial',   'Print starting values for each chain: 0-no, 1-yes, 2-add injection value,'// &
       ' 3-add Lmax value, 4-add differences (~triples number of output lines for this part)'
  write(u,11)prStat, 'prStat',   'Print statistics: 0-no, 1-yes, for default probability interval,'// &
       ' 2-yes, for all probability intervals'
  write(u,11)prCorr, 'prCorr',   'Print correlations: 0-no, 1-yes'
  write(u,11)prAcorr, 'prAcorr',   'Print autocorrelations: 0-no, 1-some, 2-more'
  write(u,11)nAcorr, 'nAcorr',   'Compute nAcorr steps of autocorrelations if prAcorr>0 or plAcorr>0 (default: 100)'
  write(u,11)prIval, 'prIval',   'Print interval info: 0-no, 1-for run with injected signal, 2-for run without injection, 3-both'
  write(u,11)prConv, 'prConv',   'Print convergence information for multiple chains to screen and chains plot: 0-no,'// &
       ' 1-one summary line, 2-add total chain stdevs, 3-add medians, stdevs for each chain'
  write(u,11)saveStats, 'saveStats',   'Save statistics (statistics, correlations, intervals) to file: 0-no, 1-yes,'// &
       ' 2-yes + copy in PS'
  write(u,11)savePDF, 'savePDF',   'Save the binned data for 1d and/or 2d pdfs (depending on plPDF1D and plPDF2D).'// &
       '  This causes all 12 parameters + m1,m2 to be saved and plotted(!), which is slighty annoying'
  write(u,11)tailoredOutput, 'tailoredOutput',   'Save (ascii) output for a specific purpose, e.g. table in a paper:'// &
       '  0-no, tO>0: use the hardcoded block tO in the subroutine tailored_output() in analyseMCMC_tailored.f90'
  
  
  write(u,'(/,A)')' Select which plots to make:'
  write(u,11)plot, 'plot',   '0: plot nothing at all, 1: plot the items selected below'
  write(u,11)plLogL, 'plLogL',   'Plot log L chains: 0-no, 1-yes'
  write(u,11)plChain, 'plChain',   'Plot parameter chains: 0-no, 1-yes'
  write(u,11)plParL, 'plParL',   'Plot L vs. parameter value: 0-no, 1-yes'
  write(u,11)plJump, 'plJump',   'Plot actual jump sizes: 0-no, 1-yes: lin, 2-yes: log'
  write(u,11)plPDF1D, 'plPDF1D',   'Plot 1d posterior distributions: 0-no, 1-yes: smoothed curve, 2-yes: actual histogram.'// &
       ' If plot=0 and savePDF=1, this determines whether to write the pdfs to file or not.'
  write(u,11)plPDF2D, 'plPDF2D',   'Plot 2d posterior distributions: 0-no, 1-yes: gray + contours, 2:gray only,'// &
       ' 3: contours only. If plot=0 and savePDF=1, this determines whether to write the pdfs to file (>0) or not (=0).'
  write(u,11)plAcorr, 'plAcorr',   'Plot autocorrelations: 0-no, 1-yes'
  write(u,11)plotSky, 'plotSky',   'Plot 2d pdf with stars, implies plPDF2D>0:  0-no, 1-yes, 2-full sky w/o stars,'// &
       ' 3-full sky with stars'
  write(u,11)plAnim, 'plAnim',   'Create movie frames'
  
  
  write(u,'(/,A)')' Detailed plot settings:'
  write(u,11)chainPlI, 'chainPlI',   'Plot every chainPlI-th point in chains, logL, jump plots:  chainPlI=0: autodetermine,'// &
       ' chainPlI>0: use this chainPlI.  All states in between *are* used for statistics, pdf generation, etc.'
  write(u,11)scLogLpl, 'scLogLpl',   'Scale logL plot ranges: 0: take everything into account, including burn-in'// &
       ' and starting values;  1: take only post-burn-in and injection values into account'
  write(u,11)scChainsPl, 'scChainsPl',   'Scale chains plot ranges: 0: take everything into account, including burn-in;'// &
       '  1: take only post-burn-in and injection values into account'
  write(u,11)plInject, 'plInject',   'Plot injection values in the chains and pdfs: 0: no,  1: yes (all pars),'// &
       '  2: yes (selected pars), 3-4: as 1-2 + print value in PDF panel'
  write(u,11)plStart, 'plStart',   'Plot starting values in the chains and pdfs'
  write(u,11)plMedian, 'plMedian',   'Plot median values in the pdfs: 1-1D PDFs, 2-2D PDFs,'// &
       ' 3-both. 4-6: as 1-3 + write value in PDF panel'
  write(u,11)plRange, 'plRange',   'Plot the probability range in the pdfs: 1-1D PDFs, 2-2D PDFs,'// &
       ' 3-both. 4-6: as 1-3 + write value in PDF panel'
  write(u,11)plBurn, 'plBurn',   'Plot the burn-in in logL, the chains, etc.: 0-no, 1-vertical line, 2-colour shade, 3-both'
  write(u,11)plLmax, 'plLmax',   'Plot the position of the max logL, in the chains and pdfs: 0-no, 1-yes, 2-yes + Lmax-Npar/2'
  write(u,11)prValues, 'prValues',   'Print values (injection, median, range) in pdfs'
  write(u,11)smooth, 'smooth',   'Smooth the pdfs: 0 - no, >1: smooth over smooth bins (use ~10 (3-15)?).'// &
       '   This is 1D only for now, and can introduce artefacts on narrow peaks!'
  write(u,11)fillPDF, 'fillPDF',   'Fillstyle for the pdfs (pgsfs): 1-solid, 2-outline, 3-hatched, 4-cross-hatched'
  write(u,11)normPDF1D, 'normPDF1D',   'Normalise 1D pdfs:  0-no,  1-normalise surface area (default, a must for different'// &
       ' bin sizes),  2-normalise to height,  3-normalise to sqrt(height), nice to compare par.temp. chains'
  write(u,11)normPDF2D, 'normPDF2D',   "'Normalise' 2D pdfs; greyscale value depends on bin height:  0-linearly,"// &
       "  1-logarithmically,  2-sqrt,  3-weigted with likelihood value,  4-2D probability intervals"
  write(u,11)nAnimFrames, 'nAnimFrames',   'Number of frames for the movie'
  write(u,11)animScheme, 'animScheme',   'AnimScheme (1-3): determines what panels to show in a movie frame; see source code'
  write(u,12)Nival,ival0, 'Nival ival0',   'Number of probability intervals,  number of the default probability interval'// &
       ' (ival0<=Nival)'
  do i=1,Nival
     write(u,'(F9.5)',advance="no")ivals(i)
  end do
  write(u,*)'       Probability intervals (ivals()). Values > 0.9999 will be treated as 100%'
  
  write(u,'(/,A)')' Output format:'
  write(u,21)scrSz, 'scrSz',   'Screen size for X11 windows (PGPlot units):  MacOS: 16.4, Gentoo: 10.8'
  write(u,21)scrRat, 'scrRat',   'Screen ratio for X11 windows (PGPlot units), MacBook: 0.57'
  write(u,11)bmpXSz, 'bmpXSz',   'X-size for bitmap (pixels):  1000  !!! Too large values give incomplete 2D PDFs somehow !!!'
  write(u,11)bmpYSz, 'bmpYSz',   'Y-size for bitmap (pixels):  700'
  write(u,21)PSsz,  'PSsz',    'Size for PS/PDF (PGPlot units).  Default: 10.5   \__ Gives same result as without pgpap'
  write(u,21)PSrat, 'PSrat',   'Ratio for PS/PDF (PGPlot units). Default: 0.742  /   '
  write(u,21)scFac, 'scFac',   '!!!Not fully implemented yet!!!  Scale .png plots up by this factor, then down to the x,y'// &
       ' size indicated above to interpolate and smoothen the plot'
  write(u,11)unSharp, 'unSharp',   'Apply unsharp mask when creating .png plots. Default: 10.'
  
  write(u,'(/,A)')' Fonts, symbols, etc.:'
  write(u,11)orientation, 'orientation',   'Use portrait (1) or landscape (2) for eps/pdf; mainly useful'// &
       ' when sending a plot to a printer'
  write(u,11)fontType, 'fontType',   'Font type used for eps/pdf: 1-simple, 2-roman, 3-italic, 4-script'
  write(u,21)fontSize1D, 'fontSize1D',   'Scale the font size for 1D plots (chains, 1D PDFs) with respect to default.'// &
       '  Default: 1.0'
  write(u,21)fontSize2D, 'fontSize2D',   'Scale the font size for 2D plots (2D PDFs) with respect to default.  Default: 1.0'
  write(u,11)chainSymbol, 'chainSymbol',   'Plot symbol for the chains: 0-plot lines, !=0: plot symbols: eg: 1: dot (default),'// &
       ' 2: plus, etc.  -4: filled diamond, 16,17: filled square,circle 20: small open circle;'// &
       ' -10/-11: use a selection of open/filled symbols'
  
  write(u,'(/,A)')' Data settings:'
  write(u,'(A)')' Plot MCMC parameters (plPar()):  11-2:tc/t40, 21-2:d^3/logd, 31-2:RA/dec, 41:phase, 51-4:i/psi/thJo/phJo,'// &
       ' 61-4:Mc/eta/M1/M2, 71-3:a/th/phi for spin1, 81-3: spin2'
  write(u,11)nPlPar, 'nPlPar',   'Number of plot parameters for 1D PDFs (and chain, jump plots, max 15).'// &
       '  This is ignored when savePDF=1. Put the MCMC parameters in the line below (plPars()):'
  do i=1,nPlPar
     write(u,'(I3)',advance="no")plPars(i)
  end do
  write(u,*)
  write(u,12)panels(1:2), 'panels',   'Number of for 1D plots in x,y direction:  0: autodetermine'
  write(u,11)Nbin1D, 'Nbin1D',   'Number of bins for 1D PDFs:  0: autodetermine'
  write(u,11)Nbin2Dx, 'Nbin2Dx',   'Number of bins in x-direction for 2D PDFs and 2D probability ranges:'// &
       '  0: autodetermine (for both x and y)'
  write(u,11)Nbin2Dy, 'Nbin2Dy',   'Number of bins in y-direction for 2D PDFs and 2D probability ranges:'// &
       '  0: use Nbin2Dx, -1: use Nbin2Dx*(scr/bmp/ps)rat'
  write(u,11)Npdf2D, 'Npdf2D',     'Number of 2D-PDF plots to make:  -1: all plots (91 for 12+2 parameters),'// &
       '  >0: read parameters from the lines below'
  do i=1,Npdf2D
     write(u,12)PDF2Dpairs(i,1:2), 'PDF2Dpairs', 'Pairs of parameters to plot a 2D PDF for'
  end do
  close(u)
  
end subroutine write_settingsfile
!***********************************************************************************************************************************




!***********************************************************************************************************************************
!> \brief  Set plot settings to 'default' values

subroutine set_plotsettings()
  use analysemcmc_settings, only: Nburn,ivals,plPars,panels,PDF2Dpairs,thin,NburnFrac,autoBurnin,maxChs,maxChLen,file,colour
  use analysemcmc_settings, only: quality,reverseRead,update,mergeChains,wrapData,changeVar,prStdOut,prProgress,prRunInfo
  use analysemcmc_settings, only: prChainInfo,prInitial,prStat,prCorr,prAcorr,nAcorr,prIval,prConv,saveStats,savePDF,tailoredOutput
  use analysemcmc_settings, only: plot,plLogL,plChain,plParL,plJump,plPDF1D,plPDF2D,plAcorr,plotSky,plAnim,chainPlI,scLogLpl
  use analysemcmc_settings, only: scChainsPl,plInject,plStart,plMedian,plRange,plBurn,plLmax,prValues,smooth,fillPDF,normPDF1D
  use analysemcmc_settings, only: normPDF2D,nAnimFrames,animScheme,Nival,ival0,scrSz,scrRat,bmpXSz,bmpYSz,PSsz,PSrat,scFac,unSharp
  use analysemcmc_settings, only: orientation,fontType,fontSize1D,fontSize2D,chainSymbol,nPlPar,Nbin1D,Nbin2Dx,Nbin2Dy,Npdf2D
  
  implicit none
  
  thin = 10         ! If >1, 'thin' the output; read every thin-th line 
  Nburn = nint(1.e5) ! If >=0: override length of the burn-in phase, for all chains
  NburnFrac = 0.5   ! If !=0: override length of the burn-in phase, as a fraction of the length of each chain.
  autoBurnin = -1.  ! Determine burn-in automatically as the first iteration where log(L_chain) > max(log(L_allchains)) - autoBurnin
  maxChLen = nint(1.e8) ! Maximum chain length
  file = 1          ! Plot output to file:  0-no; screen,  >0-yes; 1-png, 2-eps, 3-pdf.
  colour = 1        ! Use colours: 0-no (grey scales), 1-yes
  quality = 0       ! 'Quality' of plot, depending on purpose: 0: draft, 1: paper, 2: talk, 3: poster
  reverseRead = 0   ! Read files reversely (anti-alphabetically), to plot coolest chain last so that it becomes better visible
  update = 0        ! Update screen plot every 10 seconds: 0-no, 1-yes
  mergeChains = 1   ! Merge the data from different files into one chain: 0-no (treat separately), 1-yes (default)
  wrapData = 1      ! Wrap the data for the parameters that are in [0,2pi]: 0-no, 1-yes (useful if the peak is around 0)
  changeVar = 1     ! Change MCMC parameters (e.g. logd->d, kappa->theta_SL, rad->deg), 2=yes + q->1/q, phi->phi-pi, m1<->m2
  
  prStdOut = 1      ! Print standard output to 1: screen, 2: text file
  prProgress = 2    ! Print general messages about the progress of the program: 0-no, 1-some, 2-more
  prRunInfo = 0     ! Print run info at read (# iterations, seed, # detectors, SNRs, data length, etc.): 0-no, 1-only for one file
  prChainInfo = 1   ! Print chain info: 1-summary (#datpts, #contr.chns),  2-detls/chain (f.name, clr, #itr, b.in, Lmax, #datpts)
  prInitial = 0     ! Print injection values, starting values and their difference
  prStat = 1        ! Print statistics: 0-no, 1-yes
  prCorr = 0        ! Print correlations: 0-no, 1-yes
  prAcorr = 0       ! Plot autocorrelations: 0-no, 1-some, 2: more
  nAcorr = 100      ! Compute prAcorr steps of autocorrelations if prAcorr>0 or plAcor>0 (default: 100)
  prIval = 0        ! Print interval info: 0-no, 1-yes
  prConv = 1        ! Print convergence information for multiple chains to screen and chains plot
  saveStats = 0     ! Save statistics (statistics, correlations, intervals) to file: 0-no, 1-yes, 2-yes + copy in PS
  savePDF = 0       ! Save the binned data for 1d and/or 2d pdfs (depending on plPDF1D and plPDF2D).
  tailoredOutput=0  ! Save output for a specific purpose, e.g. table in a paper
  
  plot = 1          ! 0: plot nothing at all, 1: plot the items selected below
  scLogLpl = 1      ! Scale logL plot ranges: 0
  scChainsPl = 1    ! Scale chains plot ranges
  plLogL = 1        ! Plot log L chains: 0-no, 1-yes
  plChain = 1       ! Plot parameter chains: 0-no, 1-yes
  plParL = 1        ! Plot L vs. parameter value: 0-no, 1-yes
  plJump = 1        ! Plot actual jump sizes
  plPDF1D = 1       ! Plot 1d posterior distributions
  plPDF2D = 2       ! Plot 2d posterior distributions
  plAcorr = 0       ! Plot autocorrelations: 0-no, 1-yes
  plotSky = 0       ! Plot 2d pdf with stars, implies plPDF2D>0
  plAnim = 0        ! Plot movie frames
  
  chainPlI = 0      ! Plot every chainPlI-th point in chains, logL, jump plots
  plInject = 1      ! Plot injection values in the chains and pdfs
  plStart = 1       ! Plot starting values in the chains and pdfs
  plMedian = 1      ! Plot median values in the pdfs: 1-1D PDFs, 2-2D PDFs, 3-both
  plRange = 1       ! Plot the probability range in the pdfs: 1-1D PDFs, 2-2D PDFs, 3-both
  plBurn = 3        ! Plot the burn-in in logL, the chains, etc.: 0-no, 1-vertical line, 2-colour shade, 3-both
  plLmax = 0        ! Plot the position of the max of logL in chains and pdfs
  prValues = 1      ! Print values (injection, median, range) in pdfs
  smooth = 3        ! Smooth the pdfs: 0 - no, >1: smooth over smooth bins (use ~10 (3-15)?)
  fillPDF = 1       ! Fillstyle for the pdfs (pgsfs): 1-solid, 2-outline, 3-hatched, 4-cross-hatched
  normPDF1D = 1     ! Normalise 1D pdfs
  normPDF2D = 0     ! 'Normalise' 2D pdfs; greyscale value depends on bin height
  nAnimFrames = 1   ! Number of frames for the movie
  animScheme = 3    ! Movie scheme: determines what panels to show in a movie frame 
  Nival = 3         ! Number of probability intervals
  ival0 = 1         ! Standard probability interval, e.g. 1 or 2, < Nival
  ivals(1:3) = (/0.6827,0.9545,0.9973/)  ! Probability intervals
  
  scrSz  = 10.8     ! Screen size for X11 windows (PGPlot units):  MacOS: 16.4, Gentoo: 10.8
  scrRat = 0.57     ! Screen ratio for X11 windows (PGPlot units), MacBook: 0.57
  bmpXSz = 1000     ! X-size for bitmap (pixels):  1000
  bmpYSz = 700      ! Y-size for bitmap (pixels):  700
  PSsz   = 10.5     ! Size for PS/PDF (PGPlot units).  Default: 10.5   \__ Gives same result as without pgpap
  PSrat  = 0.742    ! Ratio for PS/PDF (PGPlot units). Default: 0.742  /
  scFac = 1.2       ! Scale .png plots up by this factor, then down to the x,y size indicated above
  unSharp = 10      ! Apply unsharp mask when creating .png plots. Default: 10
  
  orientation = 1   ! Use portrait (1) or landscape (2) for eps/pdf; mainly useful when sending a plot to a printer
  fonttype = 1      ! Font type used for eps/pdf: 1-simple, 2-roman, 3-italic, 4-script
  fontsize1d = 1.   ! Set plot scale for 1D plots, needs to be implemented fully
  fontsize2d = 1.   ! Set plot scale for 2D plots, needs to be implemented fully. Typically, take ~1.2*fontsize1d
  chainSymbol = 1   ! Plot symbol for the chains
  
  nPlPar = 15       ! Number of plot parameters for 1D plots
  plPars(1:nPlPar) = (/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15/) ! The nPlPar plot parameters
  panels(1:2) = (/0,0/) ! Number of panels for 1D plots in x,y direction
  Nbin1D = 100      ! Number of bins for 1D PDFs
  Nbin2Dx = 60      ! Number of bins in horizontal direction for 2D PDFs
  Nbin2Dy = 40      ! Number of bins in vertical direction for 2D PDFs
  Npdf2D  = 1       ! Number of 2D PDFs to make
  PDF2Dpairs(1,1:2) = (/8,9/)  ! 2D PDFs to plot: RA,Dec
  
end subroutine set_plotsettings
!***********************************************************************************************************************************



!***********************************************************************************************************************************
!> \brief  Read the SPINspiral output files (SPINspiral.output.*)
!!
!! \retval exitcode  Exit status code (0=ok)

subroutine read_mcmcfiles(exitcode)
  use SUFR_kinds, only: double
  use SUFR_constants, only: stdErr,stdOut
  use SUFR_system, only: quit_program_error
  
  use analysemcmc_settings, only: thin,maxChLen, maxMCMCpar, prProgress
  use general_data, only: allDat,post,prior,ntot,n,nchains,nchains0,infiles,maxIter, parNames
  use mcmcrun_data, only: niter,Nburn0,detnames,detnr,parID,seed,snr,revID,ndet,flow,fhigh,t_before,nCorr,nTemps,Tmax,Tchain
  use mcmcrun_data, only: networkSNR,waveform,pnOrder,nMCMCpar,t_after,FTstart,deltaFT,samplerate,samplesize,FTsize,outputVersion
  use mcmcrun_data, only: nMCMCpar0,t0,GPStime
  use chain_data, only: is,DoverD
  
  implicit none
  integer, intent(out) :: exitcode
  integer :: i,tmpInt,io,ic,j,readerror,p
  character :: tmpStr*(99),detname*(14),firstLine*(999),infile*(99),commandline*(999),parNameStr*(999)
  real(double) :: tmpDat(maxMCMCpar),dtmpDat(maxMCMCpar),lon2ra,ra2lon
  
  exitcode = 0
  readerror = 0
  allocate(allDat(nchains,maxMCMCpar,maxIter))
  allocate(post(nchains,maxIter),prior(nchains,maxIter))
  
  
  do ic = 1,nchains0
     infile = infiles(ic)
     open(unit=10,form='formatted',status='old',file=trim(infile),iostat=io) 
     if(io.ne.0) then
        write(stdErr,'(A)')'   Error:  File not found: '//trim(infile)//', aborting.'
        exitcode = 1
        return
     end if
     rewind(10)
     
     !if(prProgress.ge.2) write(stdOut,'(A,I3,A,I3,A20)',advance="no")'    File',ic,':  '//trim(infile)// &
     !     '    Using colour',colours(mod(ic-1,ncolours)+1),': '//colournames(colours(mod(ic-1,ncolours)+1))
     
     !Read the headers
     !Determine from the length of the first line whether this is output from before (~114 characters) 
     ! or after July 2009 (~29 characters), or LALInference (~180 characters):
     !
     !  before 7/09: first line is >80 characters long header (     Niter       Nburn    seed       <d|d>...
     !  after  7/09: first line is <80 characters long version number (  SPINspiral version:    1.00)
     ! LALInference format has long first line (  LALInference version:d4cd156ea4b0174d3fbd8b67ade5584981b34aed,2010-12-15 
     !05:58:56 +0000,cbc_bayesian_devel,Vivien Raymond <vivien.raymond@ligo.org>,CLEAN: All modifications committed)
     
     read(10,'(A999)',end=199,err=199) firstLine
     outputVersion = 1.0  ! SPINspiral output, after July 2009, keep 1<=oV<2
     if(len_trim(firstLine).gt.80 .and. len_trim(firstLine).lt.140)  outputVersion = 0.0  ! SPINspiral output, before July 2009
     if(len_trim(firstLine).ge.140)  outputVersion = 2.0  ! LALInference output (after December 2010), keep 2.0<=oV<3.0
     
     if(floor(outputVersion).eq.1) read(firstLine,'(A21,F8.2)') tmpStr,outputVersion
     
     if(outputVersion > 0.5) then
        ! Read command line between version number and first header. Read nothing if no command line:
        read(10,'(A999)',end=199,err=199) commandline
        if(commandline(1:14).eq.'  Command line') outputVersion = 2.1  ! LALinference, without parameter IDs
        ! Read empty line between version number and first header:
        read(10,*,end=199,err=199) tmpStr
     end if
     
     if(outputVersion.lt.0.5) then
        read(10,*) &
             niter(ic),Nburn0(ic),seed(ic),DoverD,ndet(ic), nCorr(ic),nTemps(ic),Tmax(ic),Tchain(ic),networkSNR(ic)
     else
        !read(10,'(I10,I12,I8,F22.10,I8,  2I9,I10,F12.1,F14.6,I11,F11.1,I10)') &
        read(10,*) &
             niter(ic),Nburn0(ic),seed(ic),DoverD,ndet(ic), nCorr(ic),nTemps(ic),Tmax(ic),Tchain(ic),networkSNR(ic),waveform, &
             pnOrder,nMCMCpar
     end if
     
     nMCMCpar0 = nMCMCpar  ! nMCMCpar may change when secondary parameters are computed
     
     read(10,*,end=199,err=199) tmpStr  ! Read (empty) line above detector info
     do i=1,ndet(ic)
        read(10,*)detnames(ic,i),snr(ic,i),flow(ic,i),fhigh(ic,i),t_before(ic,i),t_after(ic,i),FTstart(ic,i),deltaFT(ic,i), &
             samplerate(ic,i),samplesize(ic,i),FTsize(ic,i)
        detname = detnames(ic,i)
        j = len_trim(detname)
        if(trim(detname).eq.'Hanford')     detnr(ic,i) = 1
        if(trim(detname).eq.'LHO_4k')      detnr(ic,i) = 1
        if(trim(detname).eq.'Livingston')  detnr(ic,i) = 2
        if(trim(detname).eq.'LLO_4k')      detnr(ic,i) = 2
        if(trim(detname).eq.'Pisa')        detnr(ic,i) = 3
        if(trim(detname).eq.'VIRGO')       detnr(ic,i) = 3
        if(trim(detname).eq.'Aigo')        detnr(ic,i) = 4
     end do
     
     parID = 0
     revID = 0
     if(outputVersion > 0.5) then
        
        if(outputVersion.lt.2.1) then  ! Then integer parameter IDs are still used
           read(10,*,iostat=io) parID(1:nMCMCpar)  ! Read parameter IDs
           if(io.ne.0) then
              write(stdErr,'(//,A,//)')'  Error reading MCMC parameter IDs, aborting...'
              stop
           end if
        end if
        
     else  !Set the parameter IDs of an old MCMC output file manually:
        
        !Assume 12-parameter Apostolatos:
        if(ic.eq.1) write(stdOut,'(A)')"  *** Note:  I'm using a hardcoded definition of the MCMC parameters,"// &
             " assuming that an Apostolatos, 1.5-pN, 12-parameter waveform was used..."
        nMCMCpar = 12  !Mc et tc ld a1 th ra de ph tJ pJ,ph_spin
        parID(1:12) = (/61,62,11,22,71,72,31,32,41,53,54,73/)
        waveform = 1
        pNorder = 1.5
        
        !Assume 12-parameter Apostolatos: lq=log(q)
        !if(ic.eq.1) write(stdOut,'(A)')"  *** Note:  I'm using a hardcoded definition of the MCMC parameters,"// &
        !     " assuming that an Apostolatos, 1.5-pN, 12-parameter waveform was used..."
        !nMCMCpar = 12  !Mc lq tc ld a1 th ra de ph tJ pJ,ph_spin
        !parID(1:12) = (/61,68,11,22,71,72,31,32,41,53,54,73/)
        !waveform = 1
        !pNorder = 1.5


        !Assume 12-parameter SpinTaylor:
        !if(ic.eq.1) write(stdOut,'(A)')"  *** Note:  I'm using a hardcoded definition of the MCMC parameters,"// &
        !" assuming that a SpinTaylor, 3.5-pN, 12-parameter waveform was used..."
        !nMCMCpar = 12  !
        !parID(1:12) = (//)  !Needs to be filled in before use
        !waveform = 2
        !pNorder = 3.5
        
        !Assume 15-parameter SpinTaylor:
        !if(ic.eq.1) write(stdOut,'(A)')"  *** Note:  I'm using a hardcoded definition of the MCMC parameters,"// &
        !" assuming that a SpinTaylor, 3.5-pN, 15-parameter waveform was used..."
        !nMCMCpar = 15  !
        !parID(1:15) = (//)  !Needs to be filled in before use
        !waveform = 3
        !pNorder = 3.5
        
        nMCMCpar0 = nMCMCpar  ! nMCMCpar may change when secondary parameters are computed
     end if
     
     read(10,'(A)',end=199,err=199) parNameStr  ! Read line with column headers (Cycle, log Post., Prior, etc)
     if(outputVersion.ge.2.1) then
        read(10,'(A)',end=199,err=199) parNameStr  ! Read empty line
        read(10,'(A)',end=199,err=199) parNameStr  ! Read empty line
        read(10,'(A)',end=199,err=199) parNameStr  ! Read empty line
        call parNames2IDs(trim(parNameStr),nMCMCpar, parID)  ! Convert parameter names to integer IDs
     end if
     
     do i=1,nMCMCpar
        revID(parID(i)) = i  ! Reverse ID
     end do
     if(prProgress.ge.3) write(stdOut,'(99A)') 'Parameters identified: ',parNames(parID(1:nMCMCpar))
     
     i=1
     do while(i.le.maxIter)
        if(outputVersion < 0.5) then
           read(10,*,iostat=io) tmpInt,post(ic,i),tmpDat(1:nMCMCpar)
        else
           read(10,*,iostat=io) tmpInt,post(ic,i),prior(ic,i),tmpDat(1:nMCMCpar)
        end if
        
        if(io.lt.0) exit  ! EOF
        if(io.gt.0) then  ! Read error
           if(readerror.eq.1) then  ! Read error in previous line as well
              if(i.lt.25) then
                 write(stdErr,'(A,I7)',advance="no")'  Read error in file '//trim(infile)//', line',i
                 write(stdErr,'(A,/)')'  Aborting program...'
                 stop
              else
                 write(stdOut,'(A,I7)',advance="no")'  Read error in file '//trim(infile)//', line',i
                 i = i-1
                 write(stdOut,'(2(A,I8))',advance='no')'  iteration:',tmpInt,', last iteration read successfully:',tmpInt
                 write(stdOut,*)
                 exit
              end if
           end if
           readerror = 1
           i = i-1
           cycle
        end if
        readerror = 0
        
        is(ic,i) = real(tmpInt)
        !if(i.lt.10) print*,tmpInt
        
        ! GPS time doesn't fit in single-precision variable
        if(ic.eq.1.and.i.eq.1) then
           dtmpDat = 0.d0
           do p=1,nMCMCpar
              if(parID(p).ge.11.and.parID(p).le.19) then
                 dtmpDat(p) = dble(floor(tmpDat(p)/10.d0)*10)  !'GPS base time', rounded off to the nearest 10s, 
                 t0 = dtmpDat(p)                               !  to allow x.xx labels on the plots for GPS-t0
                 GPStime = floor(tmpDat(p)+0.05)               !Always floor, unless >x.95s. Use e.g. in file names
              end if
           end do
        end if
        allDat(ic,1:nMCMCpar,i) = real(tmpDat(1:nMCMCpar) - dtmpDat(1:nMCMCpar))
        
        ! 'Thin' the output by reading every thin-th line, after you've read the injection and starting values:
        if(thin.gt.1.and.i.gt.1) then
           do j=1,thin-1
              read(10,*,iostat=io) tmpStr
              if(io.lt.0) exit  ! EOF
           end do
        end if
        
        tmpStr = tmpStr  ! Remove 'set but never used' warning
        if(1.eq.2) then
           !In case you ran with lon rather than RA:
           allDat(ic,revID(31),i) = real(lon2ra(dble(allDat(ic,revID(31),i)),t0))
           !allDat(ic,revID(31),i) = real(ra2lon(dble(allDat(ic,revID(31),i)),t0))  
           
           !In case only the injection value is lon rather than RA:
           if(i.eq.1) allDat(ic,revID(31),i) = real(lon2ra(dble(allDat(ic,revID(31),i)),t0))
           
           !In case only the injection value is lon rather than RA:
           if(i.ne.1) allDat(ic,revID(31),i) = real(ra2lon(dble(allDat(ic,revID(31),i)),t0))
           
           !In case all but the injection value is lon rather than RA:
           if(i.ne.1) allDat(ic,revID(31),i) = real(lon2ra(dble(allDat(ic,revID(31),i)),t0))
        end if
        
        
        i = i+1
        if(tmpInt.ge.maxChLen) exit
     end do !i
     
     if(i.ge.maxIter-2) write(stdErr,'(A)',advance="no")'   *** WARNING ***   Not all lines in this file were read    '
199  close(10)
     ntot(ic) = i-1
     n(ic) = ntot(ic)  ! n can be changed in rearranging chains, ntot wont be changed
     !if(prProgress.ge.2.and.update.ne.1) write(stdOut,'(1x,3(A,I9),A1)')' Lines:',ntot(ic),', iterations:', &
     !nint(is(ic,ntot(ic))),', burn-in:',Nburn(ic),'.'
  end do !do ic = 1,nchains0
  
  if(sum(ntot).lt.10) &
       call quit_program_error('read_mcmcfiles(): Fewer than 10 data points were read - I must have read the input file(s)'// &
       ' incorrectly somehow', stdErr)
  
end subroutine read_mcmcfiles
!***********************************************************************************************************************************



!***********************************************************************************************************************************
!> \brief  Convert parameter names to integer IDs
!!
!! \param  parNameStr  String containing parameter names
!! \param  nMCMCpar    Number of (primary and secondary) MCMC parameters
!!
!! \retval parID       Array of parameter IDs

subroutine parNames2IDs(parNameStr,nMCMCpar, parID)
  use SUFR_constants, only: stdErr
  use SUFR_system, only: quit_program_error
  use analysemcmc_settings, only: maxMCMCpar
  
  implicit none
  character, intent(in) :: parNameStr*(*)
  integer, intent(in) :: nMCMCpar
  integer, intent(out) :: parID(maxMCMCpar)
  
  integer, parameter :: npIDs = 17
  integer :: pr1,pr2,pIDs(npIDs)
  character :: pnames(npIDs)*(19),pars(npIDs)*(19)
  
  parID = 0
  pnames(1:npIDs) = [character(len=19) :: 'iota','psi','dec','ra','dist','phi_orb','time','q','mc',  &
       'a1','theta1','phi1','a2','theta2','phi2','eta','logq']  ! CHECK: time = t40? tc?
  pIDs(1:npIDs) = (/                   51,    52,   32,   31,  22,    41,       11,    67, 61,  &
       71,72,73, 81,82,83,62,68/)
  
  read(parNameStr,*) (pars(pr1),pr1=1,nMCMCpar+3)
  do pr1=1,nMCMCpar+3
     do pr2=1,npIDs
        if(trim(pars(pr1)).eq.trim(pnames(pr2))) parID(pr1-3) = pIDs(pr2)
     end do
     if(pr1.gt.3) then
        if(parID(pr1-3).eq.0) call quit_program_error('parNames2IDs: Parameter name not recognised: '//trim(pars(pr1)),StdErr)
     end if
  end do
  
end subroutine parNames2IDs
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief Extract info from the chains and print some of it to screen
!!
!! \retval exitcode  Exit status code (0=ok)
!! 
!!  \par
!!  Print MCMC run info:
!!  - determine Lmax, burn-in,  
!!  - print chain info,  
!!  - determine thinning for chains plots,  
!!  - change/add parameters, 
!!  - determine injection, start, Lmax values of chains,  
!!  - compute jumps,  
!!  - construct output file name,  
!!  - store data in selDat (from allDat)

subroutine mcmcruninfo(exitcode)  
  use SUFR_kinds, only: double
  use SUFR_constants, only: stdOut,stdErr, rpi
  use SUFR_statistics, only: compute_median_sp
  use SUFR_system, only: swapreal, warn, quit_program_error
  use aM_constants, only: waveforms,detabbrs
  
  use analysemcmc_settings, only: Nburn,update,prRunInfo,NburnFrac,thin,autoBurnin,prChainInfo,chainPlI,changeVar,prProgress
  use analysemcmc_settings, only: prInitial,mergeChains,maxMCMCpar,plInject,plStart
  
  use general_data, only: allDat,post,ntot,n,nchains,nchains0,infiles,contrChain,startval,fixedpar,selDat,iloglmax,icloglmax
  use general_data, only: contrChains,parNames,nfixedpar,outputname,maxIter
  
  use mcmcrun_data, only: niter,Nburn0,detnames,detnr,parID,seed,snr,revID,ndet,flow,fhigh,t_before,nCorr,nTemps,Tmax,Tchain
  use mcmcrun_data, only: networkSNR,waveform,pnOrder,nMCMCpar,t_after,FTstart,deltaFT,samplerate,samplesize,FTsize,outputVersion
  use mcmcrun_data, only: nMCMCpar0,t0,GPStime,totthin,loglmaxs,totiter,loglmax,totpts,totlines,offsetrun,spinningRun
  use chain_data, only: is,isburn,DoverD,jumps
  use plot_data, only: ncolours,colours,colournames,maxdots
  
  implicit none
  integer, intent(out) :: exitcode
  integer :: i,ic,j,p,maxLine
  real :: avgtotthin
  real(double) :: lon2ra,gmst
  character :: infile*(99)
  
  exitcode = 0
  totiter = 0
  do ic = 1,nchains0
     totiter = totiter + nint(is(ic,ntot(ic)))
  end do
  !if(prProgress.ge.2.and.update.ne.1) write(stdOut,'(A10,65x,2(A,I9),A1)')'Total:',' Lines:',sum(ntot(1:nchains0)),',  &
  !iterations:',totiter
  
  
  
  ! Print run info (detectors, SNR, amount of data, FFT, etc):
  if(prRunInfo.gt.0.and.update.eq.0) then
     if(prRunInfo.eq.1) write(stdOut,'(/,A)')'  Run information for chain 1:'
     if(prRunInfo.eq.2) write(stdOut,'(/,A)')'  Run information:'
     do ic = 1,nchains0
        if((prRunInfo.eq.1.and.ic.eq.1) .or. prRunInfo.eq.2) then
           infile = infiles(ic)
           if(outputVersion.lt.0.5) then
              write(stdOut,'(4x,A7,A12,A13,A10,A12,A8,A8)')'Chain','file name','colour','Niter','Nburn','seed','Ndet'
              write(stdOut,'(4x,I7,A12,A13,I10,I12,I8,I8)')ic,trim(infile(13:99)), &
                   trim(colournames(colours(mod(ic-1,ncolours)+1))),niter(ic),Nburn0(ic),seed(ic),ndet(ic)
           else
              write(stdOut,'(4x,A7,A12,A13,A10,A12,A11,A8,  2A8,2A8,  A8,  A10,  A8,  A8, 3x,A8)') 'Chain','file name','colour', &
                   'Niter','Nburn','seed','Ndet',  'Ncorr','Ntemp','Tmax','Tchain','NetSNR','<d|d>','pN','Npar','WaveForm'
              write(stdOut,'(4x,I7,A12,  A13,I10,I12,I11,I8,  2I8,2F8.1,F8.3,F10.2,F8.1,I8, 3x,A)')ic,trim(infile(19:99)), &
                   trim(colournames(colours(mod(ic-1,ncolours)+1))),niter(ic),Nburn0(ic),seed(ic),ndet(ic), &
                   nCorr(ic),nTemps(ic),real(Tmax(ic)),Tchain(ic),networkSNR(ic),DoverD,pnOrder,nMCMCpar,trim(waveforms(waveform))
           end if
           
           write(stdOut,*)
           write(stdOut,'(A14,A3,A18,4A12,A22,A17,3A14)')'Detector','Nr','SNR','f_low','f_high','before tc','after tc', &
                'Sample start (GPS)','Sample length','Sample rate','Sample size','FT size'
           do i=1,ndet(ic)
              write(stdOut,'(A14,I3,F18.8,4F12.2,F22.8,F17.7,3I14)')trim(detnames(ic,i)),detnr(ic,i),snr(ic,i),flow(ic,i), &
                   fhigh(ic,i),t_before(ic,i),t_after(ic,i),FTstart(ic,i),deltaFT(ic,i),nint(samplerate(ic,i)),samplesize(ic,i), &
                   FTsize(ic,i)
           end do
           write(stdOut,*)
        end if
     end do !do ic = 1,nchains0
  end if  !prRunInfo.gt.0
  
  maxLine = maxval(n(1:nchains0))
  
  
  
  !*** Until now, Nburn is the iteration number.
  !*** From here on, Nburn is the line number, while isburn is the iteration number
  do ic=1,nchains0
     if(Nburn(ic).le.0) Nburn(ic) = Nburn0(ic)
     if(abs(NburnFrac).gt.1.e-4.and.abs(NburnFrac).lt.1.) then
        Nburn(ic) = nint(is(ic,n(ic)) * abs(NburnFrac))
     else
        if(Nburn(ic).ge.nint(is(ic,n(ic)))) then
           !print*,Nburn(ic),nint(is(ic,n(ic)))
           if(Nburn0(ic).ge.nint(is(ic,n(ic)))) then
              write(stdErr,'(A,I3)')'   *** WARNING ***  Nburn larger than Nchain, setting Nburn to 10% for chain',ic
              Nburn(ic) = nint(is(ic,n(ic))*0.1)
           else
              Nburn(ic) = Nburn0(ic)
           end if
        end if
     end if
  end do
  
  do ic=1,nchains0
     isburn(ic) = real(Nburn(ic))
     do i=1,ntot(ic)
        if(is(ic,i).le.isburn(ic)) Nburn(ic) = i   !isburn is the true iteration number at which the burn-in ends
        totthin(ic) = nint(isburn(ic)/real(Nburn(ic)))
     end do
  end do
  avgtotthin = sum(isburn(1:nchains0))/real(sum(Nburn(1:nchains0))) !Total thinning, averaged over all chains
  
  
  
  
  ! Get point with absolute maximum likelihood over all chains:
  loglmax = -1.d99
  loglmaxs = -1.d99
  iloglmax = 0
  icloglmax = 0
  do ic=1,nchains0
     do i=3,ntot(ic)  !3: exclude injection and starting values
        if(post(ic,i).gt.loglmaxs(ic)) then
           loglmaxs(ic) = post(ic,i)
           if(post(ic,i).gt.loglmax) then
              loglmax = post(ic,i)
              iloglmax = i
              icloglmax = ic
           end if
        end if
     end do
  end do
  
  if(ilogLmax*iclogLmax.eq.0) call quit_program_error('No max logL found - perhaps I read the input files incorrectly?', &
       stdErr)
  
  if(prRunInfo.ge.1) then
     ic = icloglmax
     i = iloglmax
     infile = infiles(ic)
     write(stdOut,'(A,I4,2(A,I9),A,F10.3,A,F7.2,A)')'    Maximum likelihood point:   chain:',ic,' ('//trim(infile(19:99))// &
          '),   line:',i*thin+7+ndet(ic), &
          ',   iteration:',nint(is(ic,i)),',   max log(L):',loglmax,'  -> SNR:',sqrt(2*loglmax),'.'
     
     !Test: get parameter values for L=Lmax
     !if(prProgress.ge.3) then
     !   write(stdOut,'(I10,F14.6,1x,2F12.7,F20.8,9F12.7)')nint(is(ic,i)),loglmax,allDat(ic,2:3,i),allDat(ic,4,i)+t0, &
     !     allDat(ic,5:13,i)
     !   call mc_eta_2_m1_m2r(allDat(ic,2,i),allDat(ic,3,i), m1,m2)
     !   write(stdOut,'(F9.4,F10.4,F21.6,F11.4,F10.4,F12.4,2F11.4,F12.4,3F11.4)')m1,m2,allDat(ic,4,i)+t0,exp(allDat(ic,5,i)), &
     !     allDat(ic,6,i),acos(allDat(ic,7,i))*r2d,allDat(ic,8,i)*r2h,asin(allDat(ic,9,i))*r2d,allDat(ic,10,i)*r2d, &
     !     asin(allDat(ic,11,i))*r2d,allDat(ic,12,i)*r2d,allDat(ic,13,i)*r2d
     !end if
     if(prRunInfo.ge.1) write(stdOut,*)
  end if
  
  
  !*** AutoBurnin: for each chain, get the first point where log(L) > log(L_max)-autoBurnin:
  if(abs(autoBurnin).gt.1.e-10) then
     if(autoBurnin.lt.-1.e-10) autoBurnin = real(nMCMCpar0)/2.  ! The default value for autoBurnin = Npar/2
     loop1: do ic=1,nchains0
        isburn(ic) = is(ic,ntot(ic)) !Set burn-in to last iteration, so that chain is completely excluded if condition never met
        Nburn(ic) = ntot(ic)
        do i=2,ntot(ic) !i=1 is injection value?
           if(post(ic,i).gt.real(loglmax)-autoBurnin) then
              isburn(ic) = is(ic,i)
              Nburn(ic) = i
              cycle loop1
           end if
        end do
     end do loop1
  end if
  
  
  !Determine the total number of iterations and lines in the input and data points for statistics;
  !determine how many and which chains contribute
  totiter = 0
  totpts  = 0
  contrChains = 0
  contrChain = 0
  totlines = sum(ntot(1:nchains0))  !Total number of lines in the input files (including burn-in)
  do ic=1,nchains0
     totiter = totiter + nint(is(ic,ntot(ic)))  !Total number of iterations (including the burn-in)
     totpts = totpts + n(ic)-Nburn(ic)          !Total number of data points available for statistics, after removing the burn-in
     if(n(ic).gt.Nburn(ic)) then
        contrChains = contrChains + 1
        contrChain(ic) = 1
     end if
  end do
  
  
  !*** Print chain info to screen:
  !Print info on number of iterations, burn-in, thinning, etc.
  do ic=1,nchains0
     infile = infiles(ic)
     if(prChainInfo.ge.2.and.update.ne.1) then
        if(contrChain(ic).eq.1) then
           write(stdOut,'(A6)',advance="no")'    * '  !Flag contributing chains
        else
           write(stdOut,'(A6)',advance="no")'      '
        end if
        write(stdOut,'(A2,I2,A1,A10,A12)',advance="no") 'Ch',ic,':',trim(infile(19:99)), &
             ', '//colournames(colours(mod(ic-1,ncolours)+1))//'.'
        write(stdOut,'(A,ES7.1,A,ES7.1,A1)',advance="no") '  Lines/iter: ',real(n(ic)),'/',is(ic,n(ic)),'.'
        write(stdOut,'(A,ES7.1,A,ES7.1,A1)',advance="no") '  Burn-in: ',real(Nburn(ic)),'/',isburn(ic),'.'
        write(stdOut,'(A,F8.2,A,F6.1,A,F4.1,A1)',advance="no") '  Lmx:',loglmaxs(ic),', dLmx:',abs(loglmax-loglmaxs(ic)), &
             '/',autoBurnin,'.'
        write(stdOut,'(A,I3,A,I4,A1)',advance="no") ' Thin: file:',nint(is(ic,n(ic))/real(n(ic)*max(thin,1))), &
             ', tot:',totthin(ic),'.'
        write(stdOut,'(A,ES7.1,A1)') '  Data pts: ',abs(real(n(ic)-Nburn(ic))),'.'
     end if
  end do
  if(prChainInfo.ge.1.and.update.ne.1) &
       write(stdOut,'(4x,A, A,ES10.3, A,ES10.3, A,I4, A,ES9.2,   A,ES10.3,  A2,F5.1, A,I3,A1,I2,A1)') &
       'All chains:','  # lines:',real(totlines), ',  # iterations:',real(totiter), ',  thinning:',nint(avgtotthin), &
       'x,  med.burnin:',compute_median_sp(real(isburn(1:nChains0))),   ',  # dat.pts after burnin:',real(totpts), &
       ' (',real(totpts)/real(totlines)*100,'%), contrib.chains:',contrChains,'/',nchains0,'.'
  
  
  
  
  !*** Determine extra thinning for logL, chain, jump plots
  if(chainPlI.le.0) then
     !if(sum(ntot(1:nchains0)).gt.maxdots) then  !Change the number of points plotted in chains,logL, etc. (For all output formats)
     !Use ntot and nchains0, since n is low if many points are in the burn-in:
     chainPlI = max(1,nint(real(sum(ntot(1:nchains0)))/real(maxdots)))  
     if(prChainInfo.ge.1.and.update.eq.0) then
        if(chainPlI.gt.1) then  ! Change the number of points plotted in chains,logL, etc. (For all output formats)
           if(chainPlI.eq.1) then
              write(stdOut,'(A)', advance='no')'    Plotting every'
           else if(chainPlI.eq.2) then
              write(stdOut,'(A,I2,A)', advance='no')'    Plotting every',chainPlI,'-nd'
           else if(chainPlI.eq.2) then
              write(stdOut,'(A,I2,A)', advance='no')'    Plotting every',chainPlI,'-rd'
           else if(chainPlI.eq.2) then
              write(stdOut,'(A,I4,A)', advance='no')'    Plotting every',chainPlI,'-th'
           end if
           write(stdOut,'(A,I5,A,I5,A)')' state in likelihood, chains, jumps, etc. plots.  Average total thinning is', &
                nint(avgtotthin),'x, for these plots it is',nint(avgtotthin*real(chainPlI)),'x.'
        else
           write(stdOut,'(A,I4,A)')'    Plotting *every* state in likelihood, chains, jumps, etc. plots.'// &
                '  Average total thinning remains',nint(avgtotthin),'x for these plots.'
        end if
     end if
     write(stdOut,*)
  end if
  !if(prRunInfo.gt.0) write(stdOut,*)
  
  
  
  
  !*** Change some MCMC parameters:
  ! use SUFR_constants, only: pi
  
  if(changeVar.ge.1) then
     if(prProgress.ge.2.and.update.eq.0) write(stdOut,'(A)',advance="no")'  Changing some parameters...   '
     
     
     if(revID(61).eq.0 .and. revID(65).ne.0) then  ! Calculate Mc from Mc_16 (Mc^(1/6)):
        if(prProgress.ge.2.and.update.eq.0) write(stdOut,'(A)')'  Computing Mc from Mc^(1/6)'
        parID(nMCMCpar+1) = 61    ! Mc
        revID(61) = nMCMCpar + 1  ! Mc
        nMCMCpar = nMCMCpar + 1
        if(nMCMCpar.gt.maxMCMCpar) then
           write(stdErr,'(//,A,I4,A,I4,A,//)')'  Error:  maxMCMCpar too small.  You must increase maxMCMCpar from',maxMCMCpar, &
                ' to at least',nMCMCpar,' in order to continue.  Aborting...'
           stop
        end if
        do ic=1,nchains0
           allDat(ic,revID(61),1:ntot(ic)) = allDat(ic,revID(65),1:ntot(ic))**6   ! Mc = [Mc_16]^6
        end do
     end if
     
     
     ! Calculate q from log(q):
     if(revID(67).eq.0 .and. revID(68).ne.0) then
        if(prProgress.ge.2.and.update.eq.0) write(stdOut,'(A)')'  Computing q from log(q)'
        parID(nMCMCpar+1) = 67    ! q
        revID(67) = nMCMCpar + 1  ! q
        nMCMCpar = nMCMCpar + 1
        if(nMCMCpar.gt.maxMCMCpar) then
           write(stdErr,'(//,A,I4,A,I4,A,//)')'  Error:  maxMCMCpar too small.  You must increase maxMCMCpar from',maxMCMCpar, &
                ' to at least',nMCMCpar,' in order to continue.  Aborting...'
           stop
        end if
        do ic=1,nchains0
           if(changeVar.eq.2) then  ! for phi > pi -> logq = -logq & phi = phi - pi
              do j=1,ntot(ic)
                 if(allDat(ic,revID(41),j).gt.rpi) then
                    allDat(ic,revID(41),j) =  allDat(ic,revID(41),j) - rpi                                     ! phi = phi - pi
                    allDat(ic,revID(68),j) = -allDat(ic,revID(68),j)                                           ! log_q = -log_q
                    if(revID(63)*revID(64).ne.0) call swapreal(allDat(ic,revID(63),j),allDat(ic,revID(64),j))  ! swap m1 <-> m2
                 end if
              end do
           end if
           allDat(ic,revID(67),1:ntot(ic)) = 10.0 ** (allDat(ic,revID(68),1:ntot(ic)))                         ! q = 10**log_q
        end do
     end if
     
     
     ! Calculate eta and log(q) from q:
     if(revID(62).eq.0 .and. revID(68).eq.0 .and. revID(67).ne.0) then
        if(prProgress.ge.2.and.update.eq.0) write(stdOut,'(A)')'  Computing eta and log(q) from q'
        parID(nMCMCpar+1) = 62    ! Eta
        parID(nMCMCpar+2) = 63    ! M1
        parID(nMCMCpar+3) = 64    ! M2
        parID(nMCMCpar+4) = 66    ! Mtot
        parID(nMCMCpar+5) = 68    ! logq 
        revID(62) = nMCMCpar + 1  ! Eta
        revID(63) = nMCMCpar + 2  ! M1
        revID(64) = nMCMCpar + 3  ! M2
        revID(66) = nMCMCpar + 4  ! Mtot
        revID(68) = nMCMCpar + 5  ! logq
        nMCMCpar = nMCMCpar + 5
        if(nMCMCpar.gt.maxMCMCpar) then
           write(stdErr,'(//,A,I4,A,I4,A,//)')'  Error:  maxMCMCpar too small.  You must increase maxMCMCpar from',maxMCMCpar, &
                ' to at least',nMCMCpar,' in order to continue.  Aborting...'
           stop
        end if
        do ic=1,nchains0
           if(changeVar.eq.2) then  ! for phi > pi -> q = 1/q & phi = phi - pi
              do j=1,ntot(ic)
                 if(allDat(ic,revID(41),j).gt.rpi) then
                    allDat(ic,revID(41),j) = allDat(ic,revID(41),j) - rpi                                      ! phi = phi - pi
                    allDat(ic,revID(67),j) = 1.0 / allDat(ic,revID(67),j)                                      ! q = 1/q
                    if(revID(63)*revID(64).ne.0) call swapreal(allDat(ic,revID(63),j),allDat(ic,revID(64),j))  ! swap m1 <-> m2
                 end if
              end do
           end if
           allDat(ic,revID(62),1:ntot(ic)) =  &
                allDat(ic,revID(67),1:ntot(ic)) / (allDat(ic,revID(67),1:ntot(ic)) + 1.0)**2                ! eta = q/(1+q)^2
           allDat(ic,revID(68),1:ntot(ic)) = log10(allDat(ic,revID(67),1:ntot(ic)))                         ! log_q = log(q)
           do i=1,ntot(ic)
              call mc_q_2_m1_m2r(allDat(ic,revID(61),i),allDat(ic,revID(67),i), allDat(ic,revID(63),i),allDat(ic,revID(64),i))
           end do
           allDat(ic,revID(66),1:ntot(ic)) = allDat(ic,revID(63),1:ntot(ic)) + allDat(ic,revID(64),1:ntot(ic))     ! Mtot = m1 + m2
        end do
     end if
     
     
     ! Calculate the individual masses from Mch and eta:
     if(revID(61)*revID(62).ne.0 .and. revID(63)+revID(64).eq.0) then
        if(prProgress.ge.2.and.update.eq.0) write(stdOut,'(A)')'  Computing M1, M2, Mtot from Mc, eta'
        parID(nMCMCpar+1) = 63    ! M1
        parID(nMCMCpar+2) = 64    ! M2
        parID(nMCMCpar+3) = 66    ! Mtot
        revID(63) = nMCMCpar + 1  ! M1
        revID(64) = nMCMCpar + 2  ! M2
        revID(66) = nMCMCpar + 3  ! Mtot
        nMCMCpar = nMCMCpar + 3
        if(nMCMCpar.gt.maxMCMCpar) then
           write(stdErr,'(//,A,I4,A,I4,A,//)')'  Error:  maxMCMCpar too small.  You must increase maxMCMCpar from',maxMCMCpar, &
                ' to at least',nMCMCpar,' in order to continue.  Aborting...'
           stop
        end if
        do ic=1,nchains0
           do i=1,ntot(ic)
              call mc_eta_2_m1_m2r(allDat(ic,revID(61),i),allDat(ic,revID(62),i), allDat(ic,revID(63),i),allDat(ic,revID(64),i))
           end do
           allDat(ic,revID(66),1:ntot(ic)) = allDat(ic,revID(63),1:ntot(ic)) + allDat(ic,revID(64),1:ntot(ic))     ! Mtot = m1 + m2
        end do
        
        ! Calculate Mc, eta from the individual masses:
     else if(revID(61)+revID(62).eq.0 .and. revID(63)*revID(64).ne.0) then
        if(prProgress.ge.2.and.update.eq.0) write(stdOut,'(A)')'  Computing Mc, eta from M1, M2'
        parID(nMCMCpar+1) = 61    ! Mc
        parID(nMCMCpar+2) = 62    ! eta
        revID(61) = nMCMCpar + 1  ! Mc
        revID(62) = nMCMCpar + 2  ! eta
        nMCMCpar = nMCMCpar + 2
        if(nMCMCpar.gt.maxMCMCpar) then
           write(stdErr,'(//,A,I4,A,I4,A,//)')'  Error:  maxMCMCpar too small.  You must increase maxMCMCpar from',maxMCMCpar, &
                ' to at least',nMCMCpar,' in order to continue.  Aborting...'
           stop
        end if
        do ic=1,nchains0
           do i=1,ntot(ic)
              call m1_m2_2_mc_etar(allDat(ic,revID(63),i),allDat(ic,revID(64),i), allDat(ic,revID(61),i),allDat(ic,revID(62),i))
           end do
        end do
        
        ! Compute total mass (var 66) and mass ratio (var 67) (q=M2/M1, not \eta) from the individual masses:
        ! (var 65 is reserved for Mc^(1/6)) & convert q -> 1/q, logq -> -logq and phi -> phi -pi for phi > pi
        write(stdOut,'(A)')'  Computing Mtot, q, log(q) from masses'
        parID(nMCMCpar+1) = 66    ! Mtot
        parID(nMCMCpar+2) = 67    ! q
        parID(nMCMCpar+3) = 68    ! log(q)
        revID(66) = nMCMCpar + 1  ! Mtot
        revID(67) = nMCMCpar + 2  ! q
        revID(68) = nMCMCpar + 3  ! log(q)
        nMCMCpar = nMCMCpar + 3
        if(nMCMCpar.gt.maxMCMCpar) then
           write(stdErr,'(//,A,I4,A,I4,A,//)')'  Error:  maxMCMCpar too small.  You must increase maxMCMCpar from',maxMCMCpar, &
                ' to at least',nMCMCpar,' in order to continue.  Aborting...'
           stop
        end if
        do ic=1,nchains0
           allDat(ic,revID(67),1:ntot(ic)) = allDat(ic,revID(64),1:ntot(ic)) / allDat(ic,revID(63),1:ntot(ic))     ! q = m2 / m1
           if(changeVar.eq.2) then  ! m2/m1 for q<1, & phi<pi and m1/m2 for q>1 & phi >pi
              do j=1,ntot(ic)
                 if(allDat(ic,revID(41),j).gt.rpi) then
                    allDat(ic,revID(41),j) = allDat(ic,revID(41),j) - rpi                                          ! phi = phi - pi
                    allDat(ic,revID(67),j) = 1.0 / allDat(ic,revID(67),j)                                          ! q = 1/q = m1/m2
                    if(revID(63)*revID(64).ne.0) call swapreal(allDat(ic,revID(63),j),allDat(ic,revID(64),j))      ! swap m1 <-> m2
                 end if
              end do
           end if
           allDat(ic,revID(66),1:ntot(ic)) = allDat(ic,revID(63),1:ntot(ic)) + allDat(ic,revID(64),1:ntot(ic))     ! Mtot = m1 + m2
           allDat(ic,revID(68),1:ntot(ic)) = log10(allDat(ic,revID(67),1:ntot(ic)))                                ! log_q = log(q)
        end do
     end if !if(revID(61)+revID(62).eq.0 .and. revID(63)*revID(64).ne.0) 
     
     ! Compute inclination and polarisation angle from RA, Dec, theta_J0, phi_J0:
     if(revID(11)*revID(31)*revID(32)*revID(53)*revID(54).ne.0) then  ! Then all of these parameters are defined
        do ic=1,nchains0
           do i=1,ntot(ic)
              ! Input: RA, Dec, phi_Jo (hh->RA), theta_Jo (in rad), output: inclination, polarisation angle (rad):
              call compute_incli_polangr(allDat(ic,revID(31),i), asin(allDat(ic,revID(32),i)), &
                   real(lon2ra(dble(allDat(ic,revID(54),i)), dble(allDat(ic,revID(11),i)) + t0)), asin(allDat(ic,revID(53),i)), &
                   allDat(ic,revID(53),i),allDat(ic,revID(54),i))  
              allDat(ic,revID(53),i) = cos(allDat(ic,revID(53),i))    !i -> cos(i)
           end do
        end do !ic
        parID(revID(53)) = 51  ! Was sin(thJ0), now cos(i)
        parID(revID(54)) = 52  ! Was phi_J0, now psi
        revID(51) = revID(53)  ! Now cos(i)
        revID(52) = revID(54)  ! Now psi
        revID(53) = 0          ! No longer defined
        revID(54) = 0          ! No longer defined
     end if
     
  end if !if(changeVar.ge.1)
  
  
  
  !*** Put plot data in startval and jumps.  Print initial and starting values to screen.
  !Startval: 1: injection value, 2: starting value, 3: Lmax value
  jumps = 0.
  offsetrun = 0
  if(prInitial.ne.0) then
     if(prInitial.eq.1) write(stdOut,'(/,A)')'  Starting values for the chains:'
     if(prInitial.eq.2) write(stdOut,'(/,A)')'  Injection and starting values for the chains:'
     if(prInitial.ge.3) write(stdOut,'(/,A)')'  Injection, starting and Lmax values for the chains:'
     write(stdOut,'(15x,A10)',advance="no")'log Post'
     do p=1,nMCMCpar
        write(stdOut,'(A8)',advance="no") trim(parNames(parID(p)))
     end do
     write(stdOut,*)
  end if
  
  startval = 0.
  do ic=1,nchains
     do i=1,2
        if(nint(is(ic,i)).eq.-1) startval(ic,1:nMCMCpar,1)  = allDat(ic,1:nMCMCpar,i)  ! Injection value
        if(nint(is(ic,i)).eq.0)  startval(ic,1:nMCMCpar,2)  = allDat(ic,1:nMCMCpar,i)  ! Starting value
     end do
     
     startval(ic,1:nMCMCpar,3)    = allDat(icloglmax,1:nMCMCpar,iloglmax)  ! Lmax value
     jumps(ic,1:nMCMCpar,2:n(ic)) = allDat(ic,1:nMCMCpar,2:n(ic)) -  allDat(ic,1:nMCMCpar,1:n(ic)-1)
     if(prInitial.ne.0) then 
        if(ic.eq.1.and.prInitial.ge.2) then
           write(stdOut,'(4x,A11)',advance="no")'Injection:  '
           write(stdOut,'(F10.3)',advance="no")post(ic,1)
           do p=1,nMCMCpar
              write(stdOut,'(F8.3)',advance="no")startval(1,p,1)
           end do
           write(stdOut,*)
           if(prInitial.ge.4) write(stdOut,*)
        end if
        if(abs((sum(startval(ic,1:nMCMCpar,1))-sum(startval(ic,1:nMCMCpar,2)))/sum(startval(ic,1:nMCMCpar,1))).gt.1.e-10) then
           offsetrun = 1
           write(stdOut,'(I4,A1,A10)',advance="no")ic,':','  Start: '
           write(stdOut,'(F10.3)',advance="no")post(ic,2)
           do p=1,nMCMCpar
              write(stdOut,'(F8.3)',advance="no")startval(ic,p,2)
           end do
           write(stdOut,*)
           if(prInitial.ge.4) then
              write(stdOut,'(5x,A10)',advance="no")'Diff:  '
              write(stdOut,'(F10.3)',advance="no")abs(post(ic,1)-post(ic,2))
              do p=1,nMCMCpar
                 write(stdOut,'(F8.3)',advance="no")abs(startval(ic,p,1)-startval(ic,p,2))
              end do
              write(stdOut,'(/)')
           end if
        end if
     end if
  end do
  
  ! If the injection or starting values are not found, don't plot them:
  if(abs(sum(startval(:,1:nMCMCpar,1))).lt.1.e-10) then
     if(plInject.ne.0.and.prProgress.ge.2) write(*,'(A)') 'No injection values found, setting plInject to zero'
     plInject = 0
  end if
  if(abs(sum(startval(:,1:nMCMCpar,2))).lt.1.e-10) then
     if(plStart.ne.0.and.prProgress.ge.2) write(*,'(A)') 'No starting values found, setting plStart to zero'
     plStart  = 0
  end if
  
  if(prInitial.ge.3) then
     write(stdOut,'(5x,A10)',advance="no")'Lmax:  '
     write(stdOut,'(F10.3)',advance="no")post(icloglmax,iloglmax)
     do p=1,nMCMCpar
        write(stdOut,'(F8.3)',advance="no")startval(1,p,3)
     end do
     write(stdOut,*)
     if(prInitial.ge.4) then
        do ic=1,1 !nchains
           write(stdOut,'(I4,A1,A10)',advance="no")ic,':','Diff:  '
           write(stdOut,'(F10.3)',advance="no")abs(post(ic,1)-post(icloglmax,iloglmax))
           do p=1,nMCMCpar
              write(stdOut,'(F8.3)',advance="no")abs(startval(ic,p,1)-startval(ic,p,3))
           end do
           write(stdOut,*)
        end do
     end if
     write(stdOut,*)
  end if
  
  
  !if(prProgress.ge.1.and.update.eq.0) write(stdOut,'(A)')'  Done.'
  if(prProgress.ge.2.and.update.eq.0) write(stdOut,'(A,I12,A,F7.4)')'  t0:',nint(t0), '  GMST:',gmst(t0)
  
  
  !Check which parameters were fixed during the MCMC run:
  fixedpar = 0
  do ic=1,nchains
     do p=1,nMCMCpar
        !Doesn't matter in which chain this happens:
        if( abs(minval(allDat(ic,p,5:n(ic))) - maxval(allDat(ic,p,5:n(ic))) ) .lt. 1.d-6) fixedpar(p) = 1
     end do
  end do
  nfixedpar = sum(fixedpar)
  
  !Is this a 'spinning run' or not?  -  spinningRun: 0-no, 1-one spin, 2-two spins
  spinningRun = 0
  if(revID(71).ne.0) then
     if(fixedpar(revID(71)).eq.0) spinningRun = spinningRun + 1
  end if
  if(revID(81).ne.0) then
     if(fixedpar(revID(81)).eq.0) spinningRun = spinningRun + 1
  end if
  
  
  !*** Construct output file name:  GPS0929052945_H1L1V__Apostolatos_1.5pN_SP  for GPS time, detectors, waveform-pN Spinning
  ic = 1
  write(outputname,'(A3,I10.10,A1)')'GPS',GPStime,'_'
  do i=1,ndet(ic)
     write(outputname,'(A)')trim(outputname)//trim(detabbrs(detnr(ic,i)))
  end do
  write(outputname,'(A,F3.1,A)')trim(outputname)//'__'//trim(waveforms(waveform))//'_',pnOrder,'pN'
  if(spinningRun.gt.0) then
     write(outputname,'(A,I1,A)')trim(outputname)//'_',spinningRun,'sp'
  else
     write(outputname,'(A)')trim(outputname)//'_ns'
  end if
  
  
  
  
  !*** Put data in selDat
  if(mergeChains.eq.1) then  !Merge chains, leave out burn-in (then nchains = 1)
     allocate(selDat(1,maxMCMCpar,nchains*maxLine))
     j = 1
     do ic=1,nchains
        do i=Nburn(ic)+1,n(ic)
           !selDat has the same structure as allDat, but contains only data AFTER the burn-in:
           selDat(1,1:nMCMCpar,j) = allDat(ic,1:nMCMCpar,i)  
           j = j+1
        end do
     end do
     nchains = 1
     n(1) = j-1
     !if(prProgress.ge.1) write(stdOut,'(A,I8,A,ES7.1,A)')'  Data points in combined chains: ',n(1),'  (',real(n(1)),')'
  else
     allocate(selDat(nchains,maxMCMCpar,maxLine))
     do ic=1,nchains
        !SelDat has the same structure as allDat, but contains only info AFTER the burn-in:
        selDat(ic,1:nMCMCpar,1:n(ic)-Nburn(ic)) = allDat(ic,1:nMCMCpar,Nburn(ic)+1:n(ic))
        n(ic) = n(ic)-Nburn(ic) !n(ic)=0 if a chain is not contributing (in which case contrChain(ic)=0)!
     end do
     !if(prProgress.ge.1) write(stdOut,'(A,I8)')' Datapoints in combined chains: ',sum(n(1:nchains))
  end if
  
end subroutine mcmcruninfo
!***********************************************************************************************************************************



!***********************************************************************************************************************************
!> \brief  Save post-burnin data to file (__data.txt)
!!
!! \retval exitcode  Exit status code (0=ok)

subroutine save_data(exitcode)
  use general_data, only: outputname,outputdir,n, selDat
  use mcmcrun_data, only: nMCMCpar,parID
  
  implicit none
  integer, intent(out) :: exitcode
  integer :: i,o,p
  
  exitcode = 0
  o = 20 !Output port
  open(unit=o, form='formatted', status='replace',file=trim(outputdir)//'/'//trim(outputname)//'__data.txt')
  
  do p=1,nMCMCpar
     if(parID(p).eq.61) write(o,'(A10)',advance="no")'mchirp'
     if(parID(p).eq.62) write(o,'(A10)',advance="no")'eta'
     if(parID(p).eq.67) write(o,'(A10)',advance="no")'q'
     if(parID(p).eq.68) write(o,'(A10)',advance="no")'log(q)'
     if(parID(p).eq.11.or.parID(p).eq.12) write(o,'(A10)',advance="no")'time'
     if(parID(p).eq.22) write(o,'(A10)',advance="no")'log(dist)'
     if(parID(p).eq.31) write(o,'(A10)',advance="no")'RA'
     if(parID(p).eq.32) write(o,'(A10)',advance="no")'sin(dec)'
     if(parID(p).eq.51) write(o,'(A10)',advance="no")'cos(iota)'
     if(parID(p).eq.41) write(o,'(A10)',advance="no")'phi0'
     if(parID(p).eq.52) write(o,'(A10)',advance="no")'psi'
     !write(o,'(A10)',advance="no")trim(parNames(parID(p)))
  end do
  write(o,*)
  !  do p=1,nMCMCpar
  !    write(o,'(F10.5)',advance="no")startval(1,p,1)
  !  end do
  !  write(o,*)
  !  write(o,*)
  
  !  do ic=1,nChains0
  do i=1,n(1)
     do p=1,nMCMCpar
        write(o,'(F10.5)',advance="no")selDat(1,p,i)
     end do
     write(o,*)
  end do
  !  end do
  
  close(o)  ! data output file
  
end subroutine save_data
!***********************************************************************************************************************************



!***********************************************************************************************************************************
!> \brief  Compute right ascension (in radians) from longitude (radians) and GPS time (seconds)
!!
!! \param lon     Longitude
!! \param GPSsec  GPS time in seconds
!!
!! - Declination == latitude for equatorial coordinates

function lon2ra(lon, GPSsec)
  use SUFR_kinds, only: double
  use SUFR_constants, only: pi2
  
  implicit none
  real(double), intent(in) :: lon,GPSsec
  real(double) :: lon2ra,gmst
  
  lon2ra = mod(lon + gmst(GPSsec) + 10*pi2,pi2)
  
end function lon2ra
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Compute longitude (in radians) from right ascension (radians) and GPS time (seconds)
!!
!! \param ra      Right Ascension
!! \param GPSsec  GPS time in seconds
!!
!! - Declination == latitude for equatorial coordinates.

function ra2lon(ra, GPSsec)
  use SUFR_kinds, only: double
  use SUFR_constants, only: pi2
  
  implicit none
  real(double), intent(in) :: ra,GPSsec
  real(double) :: ra2lon,gmst
  
  ra2lon = mod(ra - gmst(GPSsec) + 10*pi2,pi2)
  
end function ra2lon
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Compute the 'Greenwich Mean Sidereal Time' (in radians) from GPS time (in seconds)
!!
!! \param GPSsec  GPS time in seconds
!!
!! \see K.R. Lang (1999), p.80sqq.

function gmst(GPSsec)
  use SUFR_kinds, only: double
  use SUFR_constants, only: stdErr
  use SUFR_constants, only: pi2  ! 2*pi
  use SUFR_system, only: warn
  
  implicit none
  real(double), intent(in) :: GPSsec
  real(double) :: gmst,seconds,days,centuries,secCurrentDay
  real(double) :: gps0,leapseconds
  
  gps0 = 630720013.d0  ! GPS time at 1/1/2000 at midnight
  leapseconds = 32.d0  ! At Jan 1st 2000
  if(GPSsec.gt.820108813.d0) leapseconds = leapseconds + 1.d0  ! leapsecond after 1/1/2006
  if(GPSsec.gt.914803214.d0) leapseconds = leapseconds + 1.d0  ! leapsecond after 1/1/2009
  if(GPSsec.lt.630720013.d0) call warn('GMSTs before 01/01/2000 are inaccurate!', stdErr)
  
  ! Time since 1/1/2000 midnight:
  seconds       = (GPSsec - gps0) + (leapseconds - 32.d0)
  days          = floor(seconds/86400.d0) - 0.5d0
  secCurrentDay = mod(seconds, 86400.d0)
  centuries     = days/36525.d0
  gmst = 24110.54841d0 + (centuries*(8640184.812866d0 + centuries*(0.093104d0 + centuries*6.2d-6)))
  gmst = gmst + secCurrentDay * 1.002737909350795d0   ! UTC day is 1.002 * MST day
  gmst = mod(gmst/86400.d0,1.d0)
  gmst = gmst * pi2
  
end function gmst
!***********************************************************************************************************************************



!***********************************************************************************************************************************
!> \brief
!!
!! Uses lubksb,ludcmp

subroutine savgol(c, np,nl,nr,ld,m)
  implicit none
  integer, intent(in) :: np,nl,nr,ld,m
  real, intent(out) :: c(np)
  
  integer, parameter :: mmax=6
  
  integer :: imj,ipj,j,k,kk,mm,indx(mmax+1)
  real :: d,fac,sum,a(mmax+1,mmax+1),b(mmax+1)
  
  !if(np.lt.nl+nr+1.or.nl.lt.0.or.nr.lt.0.or.ld.gt.m.or.m.gt.mmax.or.nl+nr.lt.m) pause 'bad args in savgol'
  if(np.lt.nl+nr+1.or.nl.lt.0.or.nr.lt.0.or.ld.gt.m.or.m.gt.mmax.or.nl+nr.lt.m) write(0,'(A)')' Bad args in savgol'
  
  do ipj=0,2*m
     sum = 0.
     if(ipj.eq.0)sum = 1.
     do k=1,nr
        sum = sum+float(k)**ipj
     end do
     do k=1,nl
        sum = sum+float(-k)**ipj
     end do
     mm = min(ipj,2*m-ipj)
     do imj=-mm,mm,2
        a(1+(ipj+imj)/2,1+(ipj-imj)/2) = sum
     end do
  end do
  
  call ludcmp(a,m+1,mmax+1,indx,d)
  
  do j=1,m+1
     b(j) = 0.
  end do
  
  b(ld+1) = 1.
  
  call lubksb(a,m+1,mmax+1,indx,b)
  
  do kk=1,np
     c(kk) = 0.
  end do
  
  do k=-nl,nr
     sum = b(1)
     fac = 1.
     do mm=1,m
        fac = fac*real(k)
        sum = sum + b(mm+1)*fac
     end do
     kk = mod(np-k,np) + 1
     c(kk) = sum
  end do
  
end subroutine savgol
!***********************************************************************************************************************************


!***********************************************************************************************************************************
subroutine lubksb(a,n,np,indx,b)
  implicit none
  integer, intent(in) :: n,np,indx(n)
  real, intent(in) :: a(np,np)
  real, intent(out) :: b(n)
  
  integer :: i,ii,j,ll
  real :: sum
  
  ii = 0
  do i=1,n
     ll = indx(i)
     sum = b(ll)
     b(ll) = b(i)
     if (ii.ne.0) then
        do j=ii,i-1
           sum = sum-a(i,j)*b(j)
        end do
     else if (sum.ne.0.) then
        ii = i
     end if
     b(i) = sum
  end do
  
  do i=n,1,-1
     sum = b(i)
     do j=i+1,n
        sum = sum - a(i,j)*b(j)
     end do
     b(i) = sum/a(i,i)
  end do
  
end subroutine lubksb
!***********************************************************************************************************************************



!***********************************************************************************************************************************
subroutine ludcmp(a,n,np,indx,d)
  implicit none
  integer, intent(in) :: n,np
  integer, intent(out) :: indx(n)
  real, intent(inout) :: a(np,np)
  real, intent(out) :: d
  
  integer, parameter :: nmax=500
  real, parameter :: tiny=1.0e-20
  integer :: i,imax,j,k
  real :: aamax,dum,sum,vv(nmax)
  
  imax = 0
  
  d = 1.
  do i=1,n
     aamax = 0.
     do j=1,n
        if (abs(a(i,j)).gt.aamax) aamax = abs(a(i,j))
     end do
     !if (aamax.eq.0.) pause 'singular matrix in ludcmp'
     if(aamax.eq.0.) write(0,'(A)')' Singular matrix in ludcmp'
     vv(i) = 1./aamax
  end do
  
  do j=1,n
     do i=1,j-1
        sum = a(i,j)
        do k=1,i-1
           sum = sum-a(i,k)*a(k,j)
        end do
        a(i,j) = sum
     end do
     
     aamax = 0.
     
     do i=j,n
        sum = a(i,j)
        do k=1,j-1
           sum = sum-a(i,k)*a(k,j)
        end do
        a(i,j) = sum
        dum = vv(i)*abs(sum)
        if (dum.ge.aamax) then
           imax = i
           aamax = dum
        end if
     end do
     
     if (j.ne.imax) then
        do k=1,n
           dum = a(imax,k)
           a(imax,k) = a(j,k)
           a(j,k) = dum
        end do
        d = -d
        vv(imax) = vv(j)
     end if
     
     indx(j) = imax
     if(a(j,j).eq.0.) a(j,j) = tiny
     
     if(j.ne.n) then
        dum = 1./a(j,j)
        do i=j+1,n
           a(i,j) = a(i,j)*dum
        end do
     end if
  end do
  
end subroutine ludcmp
!***********************************************************************************************************************************






!***********************************************************************************************************************************
!> \brief  Returns angle in radians between 0 and 2pi (double precision)
!!
!! \param x  Angle (rad)

function drev2pi(x)
  use SUFR_kinds, only: double
  use SUFR_constants, only: pi
  
  implicit none
  real(double), intent(in) :: x
  real(double) :: drev2pi
  
  drev2pi = x - floor(x/(2*pi))*2*pi
  
end function drev2pi
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Returns periodic value x between 0 and per
!!
!! \param x    Input value
!! \param per  Period of cycle

function revper(x,per)
  implicit none
  real, intent(in) :: x,per
  real :: revper
  
  revper = x - real(floor(x/per)) * per
  
end function revper
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Returns angle in radians between -pi and pi
!!
!! \param x  Angle (rad)

function revpipi(x)
  use SUFR_constants, only: rpi,rpi2
  
  implicit none
  real, intent(in) :: x
  real :: revpipi
  
  revpipi = x - real(floor(x/rpi2)) * rpi2
  if(revpipi.gt.rpi) revpipi = revpipi - rpi2
  
end function revpipi
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Returns angle in degrees between 0 and 360
!!
!! \param x  Angle (deg)

function rev360(x)
  implicit none
  real, intent(in) :: x
  real :: rev360
  
  rev360 = x - real(floor(x/360.)) * 360.
  
end function rev360
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Returns angle in degrees between 0 and 180
!!
!! \param x  Angle (deg)

function rev180(x)
  implicit none
  real, intent(in) :: x
  real :: rev180
  
  rev180 = x - real(floor(x/180.)) * 180.
  
end function rev180
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Returns angle in hours between 0 and 24
!!
!! \param x  Angle (hours)

function rev24(x)
  implicit none
  real, intent(in) :: x
  real :: rev24
  
  rev24 = x - real(floor(x/24.)) * 24.
  
end function rev24
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Returns angle in radians between 0 and 2pi
!!
!! \param x  Angle (rad)

function rev2pi(x)
  implicit none
  real, intent(in) :: x
  real :: rev2pi,pi
  
  pi = 4*atan(1.)
  rev2pi = x - real(floor(x/(2.0*pi))) * 2.0*pi
  
end function rev2pi
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Returns angle in radians between 0 and pi - double
!!
!! \param x  Angle (rad)

function drevpi(x)
  use SUFR_kinds, only: double
  use SUFR_constants, only: pi
  
  implicit none
  real(double), intent(in) :: x
  real(double) :: drevpi
  
  drevpi = x - real(floor(x/pi)) * pi
  
end function drevpi
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Returns angle in radians between 0 and pi - real
!!
!! \param x  Angle (rad)

function rrevpi(x)
  use SUFR_constants, only: rpi
  implicit none
  real, intent(in) :: x
  real :: rrevpi
  
  rrevpi = x - real(floor(x/rpi)) * rpi
  
end function rrevpi
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Print angle as mm:ss.s string, input in hours
!!
!! a1  Angle (hours)

function tms(a1)
  use SUFR_kinds, only: double
  
  implicit none
  real(double), intent(in) :: a1
  
  real(double) :: a,s
  integer :: m
  character :: tms*(8),mm*(2),ss*(4)
  
  a = a1
  m = int((a)*60.d0)
  s = (a-m/60.d0)*3600.d0
  
  write(mm,'(i2.2)') m
  write(ss,'(f4.1)') s
  if(nint(s*10).lt.100) write(ss,'(a1,f3.1)') '0',s
  write(tms,'(a2,a1,a4,a1)') mm,'m',ss,'s'
  
end function tms
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Determine the operating system type: 1-Linux, 2-MacOSX

function getos()
  use SUFR_constants, only: stdErr, homedir
  use SUFR_system, only: warn
  
  implicit none
  integer :: status,system,getos,io
  character :: ostype*(25),filename*(99)
  
  filename = trim(homedir)//'/.analysemcmc.uname.temp'
  status = system('uname &> '//trim(filename))  ! This should return "Linux" or "Darwin"
  status = status  ! Remove 'set but never used' warning
  open(unit=16,file=trim(filename), status='old', form='formatted',iostat=io)
  if(io.ne.0) then  ! Something went wrong - guess Linux
     call warn('getOS(): cannot determine OS - guessing Linux...', stdErr)
     getos = 1
     return
  end if
  read(16,'(A)')ostype
  close(16, status = 'delete')
  
  !write(stdOut,*)ostype
  getos = 1 !Linux
  if(ostype(1:5).eq.'Darwi') getos = 2 !MacOSX
  
end function getos
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Get time stamp in seconds since 1970-01-01 00:00:00 UTC, mod countmax

function timestamp()
  use SUFR_kinds, only: double
  
  implicit none
  real(double) :: timestamp
  integer :: count,countmax,countrate
  
  call system_clock(count,countrate,countmax)
  timestamp = dble(mod(count+countmax,countmax))/dble(countrate)
  
end function timestamp
!***********************************************************************************************************************************



!***********************************************************************************************************************************
!> \brief  Transforms longitude l, latitude b and radius r into a vector with length r.  Use r=1 for a unit vector
!!
!! \param l     Longitude (rad)
!! \param b     Latitude (rad)
!! \param r     Radius
!! \retval vec  3D vector with the same units as r

subroutine lbr2vec(l,b,r,vec)
  use SUFR_kinds, only: double
  implicit none
  real(double), intent(in) :: l,b,r
  real(double), intent(out) :: vec(3)
  
  real(double) :: sinb,cosb
  
  sinb = sin(b)
  cosb = sqrt(1.d0-sinb*sinb)
  vec(1) = cos(l) * cosb         ! 'Greenwich'
  vec(2) = sin(l) * cosb         ! 'Ganges'
  vec(3) = sinb;                 ! 'North Pole'
  vec = vec*r
  
end subroutine lbr2vec
!***********************************************************************************************************************************



!***********************************************************************************************************************************
!> \brief  Compute the length of a 3D cartesian vector
!!
!! \param vec  3D vector

function veclen(vec)
  use SUFR_kinds, only: double
  implicit none
  real(double), intent(in) :: vec(3)
  real(double) :: veclen
  
  veclen = sqrt(vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3))
  
end function veclen
!***********************************************************************************************************************************



!***********************************************************************************************************************************
!> \brief  Create a unit vector from a 3D cartesian vector
!!
!! \param vec  3D vector (I/O)

subroutine normvec(vec)
  use SUFR_kinds, only: double
  implicit none
  real(double), intent(inout) :: vec(3)
  real(double) :: veclen
  
  vec = vec/veclen(vec)
  
end subroutine normvec
!***********************************************************************************************************************************



!***********************************************************************************************************************************
!> \brief  Convert chirp mass and q to m1 and m2 - double precision
!!
!! \param  mc   Chirp mass (Mo)
!! \param  q    Mass ratio q
!! 
!! \retval m1   M1 (Mo)
!! \retval m2   M2 (Mo)

subroutine mc_q_2_m1_m2(mc,q, m1,m2)
  use SUFR_kinds, only: double
  implicit none
  real(double), intent(in) :: mc,q
  real(double), intent(out) :: m1,m2
  real(double) :: factor

  factor = mc*(1.d0 + q)**(0.2d0)
  m1 = factor*q**(-0.6d0)
  m2 = factor*q**(0.4d0)
  
end subroutine mc_q_2_m1_m2
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Convert chirp mass and eta to m1 and m2 - single precision
!! 
!! \param  mcr   Chirp mass (Mo)
!! \param  qr    Mass ratio q
!! 
!! \retval m1r   M1 (Mo)
!! \retval m2r   M2 (Mo)

subroutine mc_q_2_m1_m2r(mcr,qr,m1r,m2r)
  use SUFR_kinds, only: double
  implicit none
  real, intent(in) :: mcr,qr
  real, intent(out) :: m1r,m2r

  real(double) :: mc,q,m1,m2

  mc = dble(mcr)
  q = dble(qr)
  call mc_q_2_m1_m2(mc,q,m1,m2)
  m1r = real(m1)
  m2r = real(m2)

end subroutine mc_q_2_m1_m2r
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Convert chirp mass and eta to m1 and m2 - double precision
!!
!! \param  mc   Chirp mass (Mo)
!! \param  eta  Eta
!!
!! \retval m1   M1 (Mo)
!! \retval m2   M2 (Mo)

subroutine mc_eta_2_m1_m2(mc,eta, m1,m2)
  use SUFR_kinds, only: double
  implicit none
  real(double), intent(in) :: mc,eta
  real(double), intent(out) :: m1,m2
  real(double) :: dvar,mtot
  
  mtot = mc*eta**(-0.6d0)
  if(eta.le.0.25d0) then
     dvar = sqrt(1.d0-4*eta)
     m1 = mtot/2.d0 * (1.0 + dvar);
     m2 = mtot/2.d0 * (1.0 - dvar);
  else                                 ! Allow 0.25<eta<0.50
     dvar = sqrt(4*eta-1.d0)
     m1 = mtot/2.d0 * (1.0 - dvar);
     m2 = mtot/2.d0 * (1.0 + dvar);
  end if
  
end subroutine mc_eta_2_m1_m2
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Convert chirp mass and eta to m1 and m2 - single precision
!!
!! \param  mcr   Chirp mass (Mo)
!! \param  etar  Eta
!!
!! \retval m1r   M1 (Mo)
!! \retval m2r   M2 (Mo)

subroutine mc_eta_2_m1_m2r(mcr,etar,m1r,m2r)
  use SUFR_kinds, only: double
  implicit none
  real, intent(in) :: mcr,etar
  real, intent(out) :: m1r,m2r
  
  real(double) :: mc,eta,m1,m2
  
  mc = dble(mcr)
  eta = dble(etar)
  call mc_eta_2_m1_m2(mc,eta,m1,m2)
  m1r = real(m1)
  m2r = real(m2)
  
end subroutine mc_eta_2_m1_m2r
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Convert M1,M2 to Mchirp, eta  (double precision)
!!
!! \param  m1   M1 (Mo)
!! \param  m2   M2 (Mo)
!!
!! \retval mc   Chirp mass (Mo)
!! \retval eta  Eta

subroutine m1_m2_2_mc_eta(m1,m2, mc,eta)
  use SUFR_kinds, only: double
  implicit none
  real(double), intent(in) :: m1,m2
  real(double), intent(out) :: mc,eta
  real(double) :: mtot
  
  mtot = m1+m2
  eta = m1*m2/(mtot*mtot)
  mc = mtot*eta**0.6d0
  
end subroutine m1_m2_2_mc_eta
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Convert M1,M2 to Mchirp, eta  (single precision)
!!
!! \param  m1r   M1 (Mo)
!! \param  m2r   M2 (Mo)
!!
!! \retval mcr   Chirp mass (Mo)
!! \retval etar  Eta

subroutine m1_m2_2_mc_etar(m1r,m2r, mcr,etar)
  use SUFR_kinds, only: double
  implicit none
  real, intent(in) :: m1r,m2r
  real, intent(out) :: mcr,etar
  real(double) :: m1,m2,mc,eta
  
  m1 = dble(m1r)
  m2 = dble(m2r)
  call m1_m2_2_mc_eta(m1,m2,mc,eta)
  mcr = real(mc)
  etar = real(eta)
  
end subroutine m1_m2_2_mc_etar
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Convert longitude, latitude (rad) to a 3D unit vector
!!
!! \param l     Longitude, in [0,2pi[
!! \param b     Latitude, in [-pi,pi]
!!
!! \retval vec  3D unit vector

subroutine ang2vec(l,b, vec)
  use SUFR_kinds, only: double
  
  implicit none
  real(double), intent(in) :: l,b
  real(double), intent(out) :: vec(3)
  real(double) :: cosb
  
  cosb = cos(b)
  vec(1) = cos(l) * cosb
  vec(2) = sin(l) * cosb
  vec(3) = sin(b)
  
end subroutine  ang2vec
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Convert a 3D vector to longitude, latitude (rad)
!!
!! \param  vec  3D vector
!!
!! \retval l    Longitude, in [0,2pi[
!! \retval b    Latitude, in [-pi,pi]

subroutine vec2ang(vec, l,b)
  use SUFR_kinds, only: double
  implicit none
  real(double), intent(in) :: vec(3)
  real(double), intent(out) :: l,b
  real(double) :: vec1(3)
  
  vec1 = vec
  call normvec(vec1) !Make sure vec1 is normalised
  l = atan2(vec1(2),vec1(1))
  b = asin(vec1(3))
  
end subroutine  vec2ang
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Compute the dot product of two 3D cartesian vectors
!!
!! \param vec1  3D vector 1
!! \param vec2  3D vector 2

function dotproduct(vec1,vec2)
  use SUFR_kinds, only: double
  implicit none
  real(double), intent(in) :: vec1(3),vec2(3)
  real(double) :: dotproduct
  
  dotproduct = vec1(1)*vec2(1) + vec1(2)*vec2(2) + vec1(3)*vec2(3)
  
end function dotproduct
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Compute the cross (outer) product of two cartesian vectors
!!
!! \param vec1  3D vector 1
!! \param vec2  3D vector 2
!!
!! \retval crpr  Cross/outer product (vec1 x vec2)

subroutine crossproduct(vec1,vec2, crpr)
  use SUFR_kinds, only: double
  implicit none
  real(double), intent(in) :: vec1(3),vec2(3)
  real(double), intent(out) :: crpr(3)
  
  crpr(1) = vec1(2)*vec2(3) - vec1(3)*vec2(2)
  crpr(2) = vec1(3)*vec2(1) - vec1(1)*vec2(3)
  crpr(3) = vec1(1)*vec2(2) - vec1(2)*vec2(1)
  
end subroutine crossproduct
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief Compute the polarisation angle of a source with position unit vector p and orientation normal vector o
!!
!! \param p  3D position unit vector
!! \param o  3D orientation unit vector
!!
!! \see Apostolatos et al. 1994, Eq.5

function polangle(p,o)  
  use SUFR_kinds, only: double
  
  implicit none
  real(double), intent(in) :: p(3),o(3)
  real(double) :: polangle
  real(double) :: z(3),denom,ocz(3),numer,dotproduct!,datan2
  
  z = (/0.d0,0.d0,1.d0/)                                     ! Vertical unit vector
  denom = dotproduct(o,z) - dotproduct(o,p)*dotproduct(z,p)  ! Denominator
  call crossproduct(o,z,ocz)
  numer = dotproduct(p,ocz)                                  ! Numerator
  
  polangle = atan(denom/(numer+1.d-30))                      ! Take into account the degeneracy in psi, hence no atan2
  
end function polangle
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Compute the position angle of a source with position unit vector p and orientation unit vector o
!!
!! \param p  3D position unit vector
!! \param o  3D orientation unit vector

function posangle(p,o)
  use SUFR_kinds, only: double
  
  implicit none
  real(double), intent(in) :: p(3),o(3)
  real(double) :: posangle
  real(double) :: x1(3),o1(3),z(3),z1(3),dotproduct
  
  call crossproduct(p,o,x1)
  call crossproduct(x1,p,o1) !o1: projection of o in the plane of the sky
  
  z = (/0.d0,0.d0,1.d0/) !Vertical unit vector
  call crossproduct(p,z,x1)
  call crossproduct(x1,p,z1) !z1: projection of z in the plane of the sky
  
  call normvec(o1)
  call normvec(z1)
  
  posangle = acos(dotproduct(o1,z1))
  
end function posangle
!***********************************************************************************************************************************



!***********************************************************************************************************************************
!> \brief  Compute the inclination and polarisation angle for a source with position (pl,pb) and orientation (ol,ob)
!! 
!! \param  pl   Position: longtitude, in [0,2pi[  (rad)
!! \param  pb   Position: latitude, in [0,pi]  (rad)
!! \param  ol   Orientation: longtitude, in [0,2pi[  (rad)
!! \param  ob   Orientation: latitude, in [0,pi]  (rad)  
!!
!! \retval i    Inclination angle (rad)
!! \retval psi  Polarisation angle (rad)
!! 
!! \note 
!! - all variables are angles (no cos, sin)
!! - pb,ob used to be in ([-pi/2,pi/2]) now [0,pi], conf John V. & Christian R.

subroutine compute_incli_polang(pl,pb,ol,ob, i,psi) 
  use SUFR_kinds, only: double
  
  implicit none
  
  real(double), intent(in) :: pl,pb,ol,ob
  real(double), intent(out) :: i,psi
  real(double) :: p(3),o(3),dotproduct,polangle,drevpi
  
  call ang2vec(pl,pb,p)       ! Position unit vector
  call ang2vec(ol,ob,o)       ! Orientation unit vector
  
  ! Definition 1:
  ! Compute inclination angle: <0: points towards us, >0 points away from us:
  !i = pi2 - acos(dotproduct(p,o))
  ! Compute polarisation angle [-pi/2,pi/2]:
  !psi = polangle(p,o)
  
  ! Definition 2:
  ! Compute inclination angle: 0: points exactly away from us, 180 points exactly towards us, 90: in the plane of the sky:
  i = acos(dotproduct(p,o))  
  ! Compute polarisation angle [0,pi]:
  psi = drevpi(polangle(p,o))
  
end subroutine compute_incli_polang
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Compute the inclination and polarisation angle for a source with position (pl,pb) and orientation (ol,ob) - single prec.
!! 
!! \param  plr   Position: longtitude, in [0,2pi[  (rad)
!! \param  pbr   Position: latitude, in [0,pi]  (rad)
!! \param  olr   Orientation: longtitude, in [0,2pi[  (rad)
!! \param  obr   Orientation: latitude, in [0,pi]  (rad)  
!! \retval ir    Inclination angle (rad)
!! \retval psir  Polarisation angle (rad)
!! 
!! \note
!! - all variables are angles (no cos, sin)
!! - single-precision wrapper for compute_incli_polang()

subroutine compute_incli_polangr(plr,pbr,olr,obr, ir,psir)
  use SUFR_kinds, only: double
  
  implicit none
  real, intent(in) :: plr,pbr,olr,obr
  real, intent(out) :: ir,psir
  real(double) :: pl,pb,ol,ob,i,psi
  
  pl = dble(plr)
  pb = dble(pbr)
  ol = dble(olr)
  ob = dble(obr)
  call compute_incli_polang(pl,pb,ol,ob, i,psi)
  ir = real(i)
  psir = real(psi)
  
end subroutine compute_incli_polangr
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Compute the inclination and position angle for a source with position (pl,pb) and orientation (ol,ob)
!!
!! \param  pl   Position: longtitude, in [0,2pi[  (rad)
!! \param  pb   Position: latitude, in  (rad)
!! \param  ol   Orientation: longtitude, in [0,2pi[  (rad)
!! \param  ob   Orientation: latitude, in  (rad)  
!! \retval i    Inclination angle (rad)
!! \retval pa   Polarisation angle (rad)
!!
!! \note  Position angle, not polarisation angle!

subroutine compute_incli_posang(pl,pb,ol,ob, i,pa) 
  use SUFR_kinds, only: double
  
  implicit none
  real(double), intent(in) :: pl,pb,ol,ob
  real(double), intent(out) :: i,pa
  real(double) :: p(3),o(3),dotproduct,posangle
  
  call ang2vec(pl,pb,p)       ! Position unit vector
  call ang2vec(ol,ob,o)       ! Orientation unit vector
  
  ! Definition 1:
  ! Compute inclination angle: <0: points towards us, >0 points away from us
  !i = pi2 - acos(dotproduct(p,o))
  ! Compute position angle:
  !pa = drevpi(posangle(p,o))
  
  ! Definition 2:
  ! Compute inclination angle: 0: points exactly away from us, 180 points exactly towards us, 90: in the plane of the sky:
  i = acos(dotproduct(p,o))  
  ! Compute position angle:
  pa = posangle(p,o)
  
end subroutine compute_incli_posang
!***********************************************************************************************************************************



!***********************************************************************************************************************************
!> \brief  Determine the sky position pointed at by the vector that connects two detectors
!!
!! \param d1  ID Detector 1
!! \param d2  ID Detector 2
!! \param jd  Julian day
!!
!! \todo  Finish and use

subroutine detectorvector(d1,d2,jd)
  use SUFR_kinds, only: double
  implicit none
  integer, intent(in) :: d1,d2
  real(double), intent(inout) :: jd  ! should become (in) once used
  real(double) :: detcoords(3,2),vec1(3),vec2(3),dvec(3),l,b
  
  jd = 0 !get rid of warnings
  detcoords(1,:) = (/-119.41,46.45/)  !H1; l,b
  detcoords(2,:) = (/-90.77,30.56/)   !L1
  detcoords(3,:) = (/10.50,43.63/)    !V
  
  call ang2vec(detcoords(d1,1),detcoords(d1,2),vec1)
  call ang2vec(detcoords(d2,1),detcoords(d2,2),vec2)
  
  dvec = vec2 - vec1
  
  call vec2ang(dvec,l,b)  !Searched point is in zenith/nadir for an observer on this location on the globe
  
end subroutine detectorvector
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Find files in the current working directory
!!
!! \param match    Search string to match
!! \param nff      Maximum number of files to return
!! \param all      All files?  0-select manually from list, 1-always return all files in list
!!
!! \retval fnames  Array that contains the files found; make sure it has the same length as the array in the calling programme
!! \retval nf      The actual number of files returned in fnames ( = min(number found, nff))

subroutine findFiles(match,nff,all, fnames,nf)  
  use SUFR_constants, only: stdErr, homedir
  use SUFR_system, only: quit_program_error, warn
  
  implicit none
  character, intent(in) :: match*(*)
  integer, intent(in) :: nff,all
  character, intent(out) :: fnames(nff)*(99)
  integer, intent(out) :: nf
  
  integer :: i,j,k,fnum,status,system,io
  character :: names(nff)*(99),tempfile*(99)
  
  if(len_trim(homedir).ge.99) call quit_program_error('FindFiles: variable homedir not defined (forgot to call setconstants?)',1)
  
  tempfile = trim(homedir)//'/.findFile.tmp'
  ! Shell command to list all the files with the search string and pipe them to a temporary file:
  status = system('ls '//trim(match)//' 1> '//trim(tempfile)//' 2> /dev/null') 
  status = status  ! Remove 'set but never used' warning
  
  do i=1,nff
     names(i)=''
  end do
  
  k=0
  open(10,file=trim(tempfile), status='old', form='formatted',iostat=io)  ! Read the temp file and delete it when closing
  if(io.ne.0) then
     call warn('findFiles(): cannot list files in current directory...', stdErr)
     fnames = ''
     nf = 0
     return
  end if
  rewind(10)
  do i=1,nff
     read(10,'(A99)',end=100) names(i)
     k=k+1
  end do
100 continue
  close(10, status='delete')
  fnames(1) = names(1)
  nf = 1
  j = 0
  
  if(k.gt.1) then
     if(all.eq.0) then  ! Select files manually
        write(6,'(A)')'  Files found:'  ! Don't use stdOut here!
        do i=1,k
           write(6,'(I5,A3,A)')i,':  ',trim(names(i))
        end do
        write(6,*)
        write(6,'(A,I3)')'  Enter the number of the file you want to select: 1 -',k
        write(6,'(A,I3,A1)')'    (max',nff,')'
        write(6,'(A)')'      or:   0 - to select all files in the list'
        write(6,'(A)')'           -1 - when done'
        do j=1,nff
           read*,fnum
           if(fnum.lt.0) then
              nf = j-1
              return
           end if
           if(fnum.eq.0) then
              nf = min(k,nff)
              fnames(1:nf) = names(1:nf)
              return
           end if !if(fnum.eq.0)
           fnames(j) = names(fnum)
           nf = j
        end do !j 
     else  ! Select all files (all=1)
        nf = min(k,nff)
        fnames(1:nf) = names(1:nf)
        return
     end if
  end if
  
  if(k.eq.0) then
     fnames(1)='                                                                                                   '
     !write(stdErr,'(A)')'  No file found in this directory'
     nf = 0
  end if
  
end subroutine findFiles
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Compute the size needed for PGPlot to get the desired bitmap size in pixels
!!
!! \param bmpXSz   Desired x-size in pixels
!! \param bmpYSz   Desired y-size in pixels
!! \param scFac    Scale factor; produce a larger bitmap, the shrink to get smoother graphics
!!
!! \retval bmpsz   Size of the bitmap (x)
!! \retval bmprat  Aspect ration of the bitmap

subroutine compBitmapSize(bmpXSz,bmpYSz, scFac, bmpsz,bmprat)
  use aM_constants, only: use_PLplot
  implicit none
  integer, intent(in) :: bmpXSz,bmpYSz
  real, intent(in) :: scFac
  real, intent(out) :: bmpsz,bmprat
  
  if(use_PLplot) then
     bmpsz = real(bmpXSz)/300.               ! PLplot (300 dpi?), no resizing afterwards
  else
     bmpsz = real(bmpXSz-1)/85. * scFac      ! PGPlot: Make png larger, so that convert interpolates and makes the plot smoother
  end if
  bmprat = real(bmpYSz-1)/real(bmpXSz-1)
  
end subroutine compBitmapSize
!***********************************************************************************************************************************





!***********************************************************************************************************************************
!> \brief  Print a single output line to specify when and were AnalyseMCMC was run
!!
!! \param op  Output unit

subroutine print_rundata(op)
  use SUFR_constants, only: workdir,hostname,username,currenttimezonestr,currenttimestr,currentdatestr
  use general_data, only: infiles,nchains0
  implicit none
  integer, intent(in) :: op
  
  if(nchains0.eq.1) then
     write(op,'(A)', advance="no")'  Analysing 1 chain from'
  else
     write(op,'(A,I3,A)', advance="no")'  Analysing',nchains0,' chains from'
  end if
  if(index(infiles(1),'SPINspiral.output').ne.0 .or. index(infiles(1),'mcmc.output').ne.0) then
     write(op,'(A)', advance="no")' SPINspiral'
  else
     write(op,'(A)', advance="no")' LALInference'
  end if
  write(op,'(A)', advance='no')',  in '//trim(username)//'@'//trim(hostname)//':'//trim(workdir)//'/'
  write(op,'(A)')',  on '//trim(currentdatestr)//', '//trim(currenttimestr)//' ('//trim(currenttimezonestr)//').'
  
end subroutine print_rundata
!***********************************************************************************************************************************


