!Functions for analysemcmc.f


!***************************************************************************************************
subroutine setconstants
  use constants
  implicit none
  pi = 4*datan(1.d0)
  tpi = 2*pi
  pi2 = 0.5d0*pi
  r2d = 180.d0/pi
  d2r = pi/180.d0
  r2h = 12.d0/pi
  h2r = pi/12.d0
  c3rd = 1.d0/3.d0
  
  rpi = 4*atan(1.)
  rtpi = 2*rpi
  rpi2 = 0.5*rpi
  rr2d = 180./rpi
  rd2r = rpi/180.
  rr2h = 12./rpi
  rh2r = rpi/12.
  rc3rd = 1./3.
  
  
  detabbrs = (/'H1','L1','V ','H2'/)
  waveforms = (/'Apostolatos','SpinTaylor12','SpinTaylor15','PPN'/)
  
  upline = char(27)//'[2A'  !Printing this makes the cursor move up one line (actually two lines, since a hard return is included)
end subroutine setconstants
!***************************************************************************************************


!***************************************************************************************************
subroutine read_settingsfile
  use analysemcmc_settings
  implicit none
  integer :: i,u,io,io1
  character :: bla,filename*99
  real*8 :: dblvar
  filename = 'analysemcmc.dat'
  
  u = 15
  open(unit=u,form='formatted',status='old',action='read',file=trim(filename),iostat=io)
  if(io.ne.0) then
     write(0,'(A)')'  Error opening input file '//trim(filename)//', aborting...'
     stop
  end if
  
  io = 0
  io1 = 0
  
  read(u,*,iostat=io)bla
  

  read(u,*,iostat=io)bla
  read(u,*,iostat=io)thin
  read(u,*,iostat=io)nburn(1)
  do i=2,maxChs
     nburn(i) = nburn(1)
  end do
  read(u,*,iostat=io)nburnfrac
  read(u,*,iostat=io)autoburnin
  read(u,*,iostat=io)dblvar
  maxchlen = nint(dblvar)
  read(u,*,iostat=io)file
  read(u,*,iostat=io)colour
  read(u,*,iostat=io)quality
  read(u,*,iostat=io)reverseread
  read(u,*,iostat=io)update
  read(u,*,iostat=io)mergechains
  read(u,*,iostat=io)wrapdata
  read(u,*,iostat=io)changevar
  
  read(u,*,iostat=io)bla
  read(u,*,iostat=io)prprogress
  read(u,*,iostat=io)prruninfo
  read(u,*,iostat=io)prchaininfo
  read(u,*,iostat=io)prinitial
  read(u,*,iostat=io)prstat
  read(u,*,iostat=io)prcorr
  read(u,*,iostat=io)prival
  read(u,*,iostat=io)prconv
  read(u,*,iostat=io)savestats
  read(u,*,iostat=io)savepdf
  
  
  read(u,*,iostat=io)bla
  read(u,*,iostat=io)plot
  read(u,*,iostat=io)combinechainplots
  read(u,*,iostat=io)pllogl
  read(u,*,iostat=io)plchain
  read(u,*,iostat=io)plparl
  read(u,*,iostat=io)pljump
  read(u,*,iostat=io)plpdf1d
  read(u,*,iostat=io)plpdf2d
  read(u,*,iostat=io)placorr
  read(u,*,iostat=io)plotsky
  read(u,*,iostat=io)plmovie
  
  
  read(u,*,iostat=io)bla
  read(u,*,iostat=io)chainsymbol
  read(u,*,iostat=io)chainpli
  read(u,*,iostat=io)scloglpl
  read(u,*,iostat=io)scchainspl
  read(u,*,iostat=io)pltrue
  read(u,*,iostat=io)plstart
  read(u,*,iostat=io)plmedian
  read(u,*,iostat=io)plrange
  read(u,*,iostat=io)plburn
  read(u,*,iostat=io)pllmax
  read(u,*,iostat=io)prvalues
  read(u,*,iostat=io)smooth
  read(u,*,iostat=io)fillpdf
  read(u,*,iostat=io)normpdf1d
  read(u,*,iostat=io)normpdf2d
  read(u,*,iostat=io)nmovframes
  read(u,*,iostat=io)moviescheme
  read(u,*,iostat=io)nival,ival0
  read(u,*,iostat=io1)(ivals(i),i=1,nival)

  read(u,*,iostat=io)bla
  read(u,*,iostat=io)scrsz
  read(u,*,iostat=io)scrrat
  read(u,*,iostat=io)bmpxsz
  read(u,*,iostat=io)bmpysz
  read(u,*,iostat=io)pssz
  read(u,*,iostat=io)psrat
  read(u,*,iostat=io)whitebg
  read(u,*,iostat=io)scfac
  read(u,*,iostat=io)unsharp
  
  read(u,*,iostat=io)bla
  read(u,*,iostat=io)bla
  read(u,*,iostat=io)nplvar
  read(u,*,iostat=io1)(plvars(i),i=1,nplvar)
  if(io1.ne.0) nplvar = i-1
  io1 = 0
  read(u,*,iostat=io)panels(1:2)
  read(u,*,iostat=io)nbin1d
  read(u,*,iostat=io)nbin2dx
  read(u,*,iostat=io)nbin2dy
  read(u,*,iostat=io)npdf2d
  do i=1,npdf2d
     read(u,*,iostat=io1)pdf2dpairs(i,1:2)
     if(io1.ne.0) exit
  end do
  if(io1.ne.0) npdf2d = i-1
  close(u)
  
  if(io.ne.0) then
     write(0,'(/,A,I2,A,I3,A,/)')'  Error reading input file '//trim(filename)//', aborting...'
     stop
  end if
  
end subroutine read_settingsfile
!***************************************************************************************************


!***************************************************************************************************
subroutine write_settingsfile
  use analysemcmc_settings
  implicit none
  integer :: u,i
  
  u = 14
  open(unit=u,form='formatted',status='replace',file='analysemcmc.used')
  
11 format(I10,1x,A19,5x,A)
12 format(2I5,1x,A19,5x,A)
21 format(F10.5,1x,A19,5x,A)
31 format(ES10.2,1x,A19,5x,A)
  
  write(u,'(A,/)')' Input file for analysemcmc.f'
  
  
  write(u,'(/,A)')' Basic options:'
  write(u,11)thin, 'thin',   'If >1, "thin" the output; read every thin-th line '
  write(u,11)maxval(nburn), 'nburn',   'If >=0: override length of the burn-in phase, for all chains! This is now the ITERATION number (it becomes the line number later on).  Nburn > Nchain sets Nburn = 0.1*Nchain'
  write(u,21)nburnfrac, 'nburnfrac',   'If !=0: override length of the burn-in phase, as a fraction of the length of each chain. This overrides nburn above'
  write(u,21)autoburnin, 'autoburnin',   'If >0: Determine burnin automatically as the first iteration where log(L_chain) > max(log(L_allchains)) - autoburnin. Overrides nburn and nburnfrac above'
  write(u,31)dble(maxchlen), 'maxchlen',   'Maximum chain length to read in (number of iterations, not number of lines)'
  write(u,11)file, 'file',   'Plot output to file:  0-no; screen,  >0-yes; 1-png, 2-eps, 3-pdf.  Give an output path for files in the parameter "outputdir" below'
  write(u,11)colour, 'colour',   'Use colours: 0-no (grey scales), 1-yes'
  write(u,11)quality, 'quality',   '"Quality" of plot, depending on purpose: 0: draft, 1: paper, 2: talk, 3: poster'
  write(u,11)reverseread, 'reverseread',   'Read files reversely (anti-alphabetically), to plot coolest chain last so that it becomes better visible: 0-no, 1-yes, 2-use colours in reverse order too'
  write(u,11)update, 'update',   'Update screen plot every 10 seconds: 0-no, 1-yes'
  write(u,11)mergechains, 'mergechains',   'Merge the data from different files into one chain: 0-no (treat separately), 1-yes'
  write(u,11)wrapdata, 'wrapdata',   'Wrap the data for the parameters that are in [0,2pi]: 0-no, 1-yes (useful if the peak is around 0)'
  write(u,11)changevar, 'changevar',   'Change variables (e.g. logd->d, kappa->theta_SL, rad->deg)'
  
  
  write(u,'(/,A)')' Select what output to print to screen and write to file:'
  write(u,11)prprogress, 'prprogress',   'Print general messages about the progress of the program: 0-no, 1-some, 2-more, 3-debug output'
  write(u,11)prruninfo, 'prruninfo',   'Print run info (# iterations, seed, # detectors, SNRs, data length, etc.): 0-no, 1-only for one file (eg. if all files similar), 2-for all files'
  write(u,11)prchaininfo, 'prchaininfo',   'Print chain info: 1-summary (tot # data points, # contributing chains),  2-details per chain (file name, plot colour, # iterations, burnin, Lmax, # data points)'
  write(u,11)prinitial, 'prinitial',   'Print true values, starting values and their difference'
  write(u,11)prstat, 'prstat',   'Print statistics: 0-no, 1-yes, for default probability interval, 2-yes, for all probability intervals'
  write(u,11)prcorr, 'prcorr',   'Print correlations: 0-no, 1-yes'
  write(u,11)prival, 'prival',   'Print interval info: 0-no, 1-for run with injected signal, 2-for run without injection, 3-both'
  write(u,11)prconv, 'prconv',   'Print convergence information for multiple chains to screen and chains plot: 0-no, 1-one summary line, 2-add total chain stdevs, 3-add medians, stdevs for each chain'
  write(u,11)savestats, 'savestats',   'Save statistics (statistics, correlations, intervals) to file: 0-no, 1-yes, 2-yes + copy in PS'
  write(u,11)savepdf, 'savepdf',   'Save the binned data for 1d and/or 2d pdfs (depending on plpdf1d and plpdf2d).  This causes all 12 parameters + m1,m2 to be saved and plotted(!), which is slighty annoying'
  
  
  write(u,'(/,A)')' Select which plots to make:'
  write(u,11)plot, 'plot',   '0: plot nothing at all, 1: plot the items selected below'
  write(u,11)combinechainplots, 'combinechainplots',   'Combine logL and chain plots into one multipage file'
  write(u,11)pllogl, 'pllogl',   'Plot log L chains: 0-no, 1-yes'
  write(u,11)plchain, 'plchain',   'Plot parameter chains: 0-no, 1-yes'
  write(u,11)plparl, 'plparl',   'Plot L vs. parameter value: 0-no, 1-yes'
  write(u,11)pljump, 'pljump',   'Plot actual jump sizes: 0-no, 1-yes: lin, 2-yes: log'
  write(u,11)plpdf1d, 'plpdf1d',   'Plot 1d posterior distributions: 0-no, 1-yes: smoothed curve, 2-yes: actual histogram. If plot=0 and savepdf=1, this determines whether to write the pdfs to file or not.'
  write(u,11)plpdf2d, 'plpdf2d',   'Plot 2d posterior distributions: 0-no, 1-yes: gray + contours, 2:gray only, 3: contours only. If plot=0 and savepdf=1, this determines whether to write the pdfs to file (>0) or not (=0).'
  write(u,11)placorr, 'placorr',   'Plot autocorrelations: 0-no, >0-yes: plot placorr steps'
  write(u,11)plotsky, 'plotsky',   'Plot 2d pdf with stars, implies plpdf2d>0:  0-no, 1-yes, 2-full sky w/o stars, 3-full sky with stars'
  write(u,11)plmovie, 'plmovie',   'Create movie frames'
  
  
  write(u,'(/,A)')' Detailed plot settings:'
  write(u,11)chainsymbol, 'chainsymbol',   'Plot symbol for the chains: 0-plot lines, !=0: plot symbols: eg: 1: dot (default), 2: plus, etc.  -4: filled diamond, 16,17: filled square,circle 20: small open circle; -10/-11: use a selection of open/filled symbols'
  write(u,11)chainpli, 'chainpli',   'Plot every chainpli-th point in chains, logL, jump plots:  chainpli=0: autodetermine, chainpli>0: use this chainpli.  All states in between *are* used for statistics, pdf generation, etc.'
  write(u,11)scloglpl, 'scloglpl',   'Scale logL plot ranges: 0: take everything into account, including burnin and starting values;  1: take only post-burnin and true values into account'
  write(u,11)scchainspl, 'scchainspl',   'Scale chains plot ranges: 0: take everything into account, including burnin;  1: take only post-burnin and true values into account'
  write(u,11)pltrue, 'pltrue',   'Plot true values in the chains and pdfs: 0: no,  1: yes (all pars),  2: yes (selected pars), 3-4: as 1-2 + print value in PDF panel'
  write(u,11)plstart, 'plstart',   'Plot starting values in the chains and pdfs'
  write(u,11)plmedian, 'plmedian',   'Plot median values in the pdfs: 1-1D PDFs, 2-2D PDFs, 3-both. 4-6: as 1-3 + write value in PDF panel'
  write(u,11)plrange, 'plrange',   'Plot the probability range in the pdfs: 1-1D PDFs, 2-2D PDFs, 3-both. 4-6: as 1-3 + write value in PDF panel'
  write(u,11)plburn, 'plburn',   'Plot the burnin in logL, the chains, etc.'
  write(u,11)pllmax, 'pllmax',   'Plot the position of the max logL, in the chains and pdfs'
  write(u,11)prvalues, 'prvalues',   'Print values (true, median, range) in pdfs'
  write(u,11)smooth, 'smooth',   'Smooth the pdfs: 0 - no, >1: smooth over smooth bins (use ~10 (3-15)?).   This is 1D only for now, and can introduce artefacts on narrow peaks!'
  write(u,11)fillpdf, 'fillpdf',   'Fillstyle for the pdfs (pgsfs): 1-solid, 2-outline, 3-hatched, 4-cross-hatched'
  write(u,11)normpdf1d, 'normpdf1d',   'Normalise 1D pdfs:  0-no,  1-normalise surface area (default, a must for different bin sizes),  2-normalise to height,  3-normalise to sqrt(height), nice to compare par.temp. chains'
  write(u,11)normpdf2d, 'normpdf2d',   "'Normalise' 2D pdfs; greyscale value depends on bin height:  0-linearly,  1-logarithmically,  2-sqrt,  3-weigted with likelihood value,  4-2D probability intervals"
  write(u,11)nmovframes, 'nmovframes',   'Number of frames for the movie'
  write(u,11)moviescheme, 'moviescheme',   'Moviescheme (1-3): determines what panels to show in a movie frame; see source code'
  write(u,12)nival,ival0, 'nival ival0',   'Number of probability intervals,  number of the default probability interval (ival0<=nival)'
  do i=1,nival
     write(u,'(F9.5,$)')ivals(i)
  end do
  write(u,*)'       Probability intervals (ivals()). Values > 0.9999 will be treated as 100%'
  
  write(u,'(/,A)')' Output format:'
  write(u,21)scrsz, 'scrsz',   'Screen size for X11 windows (PGPlot units):  MacOS: 16.4, Gentoo: 10.8'
  write(u,21)scrrat, 'scrrat',   'Screen ratio for X11 windows (PGPlot units), MacBook: 0.57'
  write(u,11)bmpxsz, 'bmpxsz',   'X-size for bitmap (pixels):  1000  !!! Too large values give incomplete 2D PDFs somehow !!!'
  write(u,11)bmpysz, 'bmpysz',   'Y-size for bitmap (pixels):  700'
  write(u,21)pssz, 'pssz',   'Size for PS/PDF (PGPlot units).  Default: 10.5   \__ Gives same result as without pgpap'
  write(u,21)psrat, 'psrat',   'Ratio for PS/PDF (PGPlot units). Default: 0.742  /'
  write(u,11)whitebg, 'whitebg',   'Create white background for screen and .png plots: 0-no (black, default), 1-yes'
  write(u,21)scfac, 'scfac',   '!!!Not fully implemented yet!!!  Scale .png plots up by this factor, then down to the x,y size indicated above to interpolate and smoothen the plot'
  write(u,11)unsharp, 'unsharp',   'Apply unsharp mask when creating .png plots. Default: 10.'
  
  write(u,'(/,A)')' Data settings:'
  write(u,'(A)')' Plot variables:  1:logL, 2:Mc, 3:eta, 4:tc, 5:dL, 6:a, 7:th, 8:RA, 9:dec, 10:phi, 11:thJ, 12:phiJ, 13:alpha, 14:M1, 15:M2'
  write(u,11)nplvar, 'nplvar',   'Number of plot variables for 1D PDFs (and chain, jump plots, max 15).  This is ignored when savepdf=1. Put the variables in the line below:'
  do i=1,nplvar
     write(u,'(I3,$)')plvars(i)
  end do
  write(u,*)''
  !write(u,'(5x,A,5x,A)')'plvars','The actual plot variables (1-nplvar)'
  write(u,12)panels(1:2), 'panels',   'Number of for 1D plots in x,y direction:  0: autodetermine'
  write(u,11)nbin1d, 'nbin1d',   'Number of bins for 1D PDFs:  0: autodetermine'
  write(u,11)nbin2dx, 'nbin2dx',   'Number of bins in x-direction for 2D PDFs and 2D probability ranges:  0: autodetermine (for both x and y)'
  write(u,11)nbin2dy, 'nbin2dy',   'Number of bins in y-direction for 2D PDFs and 2D probability ranges:  0: use nbin2dx, -1: use nbin2dx*(scr/bmp/ps)rat'
  write(u,11)npdf2d, 'npdf2d',     'Number of 2D-PDF plots to make:  -1: all plots (91 for 12+2 parameters),  >0: read parameters from the lines below'
  do i=1,npdf2d
     write(u,12)pdf2dpairs(i,1:2), 'pdf2dpairs', 'Pairs of parameters to plot a 2D PDF for'
  end do
  close(u)
end subroutine write_settingsfile
!***************************************************************************************************




!***************************************************************************************************
subroutine set_plotsettings  !Set plot settings to 'default' values
  use analysemcmc_settings
  implicit none
  
  thin = 10         !If >1, 'thin' the output; read every thin-th line 
  nburn = 1e5       !If >=0: override length of the burn-in phase, for all chains! This is now the ITERATION number, but it becomes the line number later on in the code.  Nburn > Nchain sets Nburn = 0.1*Nchain
  nburnfrac = 0.5   !If !=0: override length of the burn-in phase, as a fraction of the length of each chain.
  autoburnin = 1.   !Determine burnin automatically as the first iteration where log(L_chain) > max(log(L_allchains)) - autoburnin
  maxchlen = 1e8    !Maximum chain length
  file = 1          !Plot output to file:  0-no; screen,  >0-yes; 1-png, 2-eps, 3-pdf.  Give an output path for files in the parameter 'outputdir' below.
  colour = 1        !Use colours: 0-no (grey scales), 1-yes
  quality = 0       !'Quality' of plot, depending on purpose: 0: draft, 1: paper, 2: talk, 3: poster
  reverseread = 0   !Read files reversely (anti-alphabetically), to plot coolest chain last so that it becomes better visible: 0-no, 1-yes, 2-use colours in reverse order too
  update = 0        !Update screen plot every 10 seconds: 0-no, 1-yes
  mergechains = 1   !Merge the data from different files into one chain: 0-no (treat separately), 1-yes
  wrapdata = 1      !Wrap the data for the parameters that are in [0,2pi]: 0-no, 1-yes (useful if the peak is around 0)
  changevar = 1     !Change variables (e.g. logd->d, kappa->theta_SL, rad->deg)
  
  prprogress = 2    !Print general messages about the progress of the program: 0-no, 1-some, 2-more
  prruninfo = 0     !Print run info at read (# iterations, seed, # detectors, SNRs, data length, etc.): 0-no, 1-only for one file (eg. if all files similar), 2-for all files
  prinitial = 0     !Print true values, starting values and their difference
  prstat = 1        !Print statistics: 0-no, 1-yes
  prcorr = 0        !Print correlations: 0-no, 1-yes
  prival = 0        !Print interval info: 0-no, 1-yes
  prconv = 1        !Print convergence information for multiple chains to screen and chains plot: 0-no, 1-one summary line, 2-medians, stdevs, etc. too.
  savestats = 0     !Save statistics (statistics, correlations, intervals) to file: 0-no, 1-yes, 2-yes + copy in PS
  savepdf = 0       !Save the binned data for 1d and/or 2d pdfs (depending on plpdf1d and plpdf2d).  This causes all 12 parameters + m1,m2 to be saved and plotted(!), which is slighty annoying
  
  plot = 1          !0: plot nothing at all, 1: plot the items selected below
  combinechainplots = 0  !Combine logL and chain plots into one multipage file
  autoburnin = 1.   !Determine burnin automatically as the first iteration where log(L_chain) > max(log(L_allchains)) - autoburnin
  scloglpl = 1      !Scale logL plot ranges: 0: take everything into account, including burnin and starting values;  1: take only post-burnin and true values into account
  scchainspl = 1    !Scale chains plot ranges: 0: take everything into account, including burnin;  1: take only post-burnin and true values into account
  pllogl = 1        !Plot log L chains: 0-no, 1-yes
  plchain = 1       !Plot parameter chains: 0-no, 1-yes
  plparl = 1        !Plot L vs. parameter value: 0-no, 1-yes
  pljump = 1        !Plot actual jump sizes
  plpdf1d = 1       !Plot 1d posterior distributions: 0-no, 1-yes: smoothed curve, 2-yes: actual histogram. If plot=0 and savepdf=1, this determines whether to write the pdfs to file or not.
  plpdf2d = 2       !Plot 2d posterior distributions: 0-no, 1-yes: gray + contours, 2:gray only, 3: contours only. If plot=0 and savepdf=1, this determines whether to write the pdfs to file (>0) or not (=0).
  placorr = 0e4     !Plot autocorrelations: 0-no, >0-yes: plot placorr steps
  plotsky = 0       !Plot 2d pdf with stars, implies plpdf2d>0:  0-no, 1-yes, 2-full sky w/o stars, 3-full sky with stars'
  plmovie = 0       !Plot movie frames
  
  chainsymbol = 1   !Plot symbol for the chains: 0-plot lines, !=0: plot symbols: eg: 1: dot (default), 2: plus, etc.  -4: filled diamond, 16,17: filled square,circle 20: small open circle
  chainpli = 0      !Plot every chainpli-th point in chains, logL, jump plots:  chainpli=0: autodetermine, chainpli>0: use this chainpli.  All states in between *are* used for statistics, pdf generation, etc.
  pltrue = 1        !Plot true values in the chains and pdfs
  plstart = 1       !Plot starting values in the chains and pdfs
  plmedian = 1      !Plot median values in the pdfs: 1-1D PDFs, 2-2D PDFs, 3-both
  plrange = 1       !Plot the probability range in the pdfs: 1-1D PDFs, 2-2D PDFs, 3-both
  plburn = 1        !Plot the burnin in logL, the chains, etc.
  pllmax = 0        !Plot the position of the max of logL in chains and pdfs
  prvalues = 1      !Print values (true, median, range) in pdfs
  smooth = 3        !Smooth the pdfs: 0 - no, >1: smooth over smooth bins (use ~10 (3-15)?).   This is 1D only for now, and can introduce artefacts on narrow peaks!
  fillpdf = 1       !Fillstyle for the pdfs (pgsfs): 1-solid, 2-outline, 3-hatched, 4-cross-hatched
  normpdf1d = 1     !Normalise 1D pdfs:  0-no,  1-normalise surface area (default, a must for different bin sizes),  2-normalise to height,  3-normalise to sqrt(height), nice to compare par.temp. chains
  normpdf2d = 0     !'Normalise' 2D pdfs; greyscale value depends on bin height:  0-linearly,  1-logarithmically,  2-sqrt,  3-weigted with likelihood value
  nmovframes = 1    !Number of frames for the movie
  moviescheme = 3   !Movie scheme: determines what panels to show in a movie frame 
  nival = 3         !Number of probability intervals
  ival0 = 1         !Standard probability interval, e.g. 1 or 2, < nival
  ivals(1:3) = (/0.6827,0.9545,0.9973/)  !Probability intervals
  
  scrsz  = 10.8     !Screen size for X11 windows (PGPlot units):  MacOS: 16.4, Gentoo: 10.8
  scrrat = 0.57     !Screen ratio for X11 windows (PGPlot units), MacBook: 0.57
  bmpxsz = 1000     !X-size for bitmap (pixels):  1000
  bmpysz = 700      !Y-size for bitmap (pixels):  700
  pssz   = 10.5     !Size for PS/PDF (PGPlot units).  Default: 10.5   \__ Gives same result as without pgpap
  psrat  = 0.742    !Ratio for PS/PDF (PGPlot units). Default: 0.742  /
  whitebg = 0       !White background for screen and .png plots: 0-no, 1-yes
  scfac = 1.2       !Scale .png plots up by this factor, then down to the x,y size indicated above to interpolate and smoothen the plot
  unsharp = 10      !Apply unsharp mask when creating .png plots. Default: 10
  
  nplvar = 15       !Number of plot variables for 1D plots
  plvars(1:nplvar) = (/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15/) !The nplvar plot variables
  panels(1:2) = (/0,0/) !Number of panels for 1D plots in x,y direction
  nbin1d = 100      !Number of bins for 1D PDFs
  nbin2dx = 60      !Number of bins in horizontal direction for 2D PDFs
  nbin2dy = 40      !Number of bins in vertical direction for 2D PDFs
  npdf2d  = 1       !Number of 2D PDFs to make
  pdf2dpairs(1,1:2) = (/8,9/)  !2D PDFs to plot: RA,Dec
end subroutine set_plotsettings
!***************************************************************************************************



!************************************************************************************************************************************
subroutine read_mcmcfiles(exitcode)  !Read the SPINspiral output files (SPINspiral.output.*)
  use analysemcmc_settings
  use general_data
  use mcmcrun_data
  use chain_data
  implicit none
  integer :: i,tmpInt,io,ic,j,exitcode,readerror,p
  character :: tmpStr*99,detname*14,firstLine*999
  real*8 :: tmpDat(maxMCMCpar),dtmpDat(maxMCMCpar)
  real :: outputVersion
  !real*8 :: lon2ra
  
  exitcode = 0
  readerror = 0
  allocate(allDat(nchains,maxMCMCpar,maxIter))
  allocate(post(nchains,maxIter),prior(nchains,maxIter))
  
  
  do ic = 1,nchains0
     if(reverseread.eq.0) then
        call getarg(ic,infile) !Read file name from the command-line arguments
     else
        call getarg(nchains0-ic+1,infile) !Read file name from the command-line arguments in reverse order
     end if
     infiles(ic) = infile
     
     open(unit=10,form='formatted',status='old',file=trim(infile),iostat=io)
     if(io.ne.0) then
        write(0,'(A)')'   Error:  File not found: '//trim(infile)//', aborting.'
        exitcode = 1
        return
     end if
     rewind(10)
     
     !if(prprogress.ge.2) write(6,'(A,I3,A,I3,A20,$)')'    File',ic,':  '//trim(infile)//'    Using colour',colours(mod(ic-1,ncolours)+1),': '//colournames(colours(mod(ic-1,ncolours)+1))
     
     !Read the headers
     !Determine from the length of the first line whether this is output from before of after July 2009
     !  before: first line is >80 characters long header (     Niter       Nburn    seed       null likelihood    Ndet    Ncorr   Ntemps      Tmax      Tchain   Network SNR)
     !  after:  first line is <80 characters long version number (  SPINspiral version:    1.00)
     
     outputVersion = 0.0  !Use old format
     read(10,'(A999)',end=199,err=199)firstLine
     if(len_trim(firstLine).lt.80) read(firstLine,'(A21,F8.2)')tmpStr,outputVersion  !Use new format

     if(outputVersion > 0.5) read(10,*,end=199,err=199)tmpStr  !Read empty line between version number and first header
     read(10,'(I10,I12,I8,F22.10,I8,  2I9,I10,F12.1,F14.6,I11,F11.1,I10)') niter(ic),nburn0(ic),seed(ic),nullh,ndet(ic), nCorr(ic),nTemps(ic),Tmax(ic),Tchain(ic),networkSNR(ic),waveform,pnOrder,nMCMCpar
     
     read(10,*,end=199,err=199)tmpStr !Read empty line above detector info
     do i=1,ndet(ic)
        !read(10,'(2x,A14,F18.8,4F12.2,F22.8,F17.7,3I14)') detnames(ic,i),snr(ic,i),flow(ic,i),fhigh(ic,i),t_before(ic,i),t_after(ic,i),FTstart(ic,i),deltaFT(ic,i),samplerate(ic,i),samplesize(ic,i),FTsize(ic,i)
        read(10,*)detnames(ic,i),snr(ic,i),flow(ic,i),fhigh(ic,i),t_before(ic,i),t_after(ic,i),FTstart(ic,i),deltaFT(ic,i),samplerate(ic,i),samplesize(ic,i),FTsize(ic,i)
        detname = detnames(ic,i)
        !write(0,'(A)')trim(detname)
        j = len_trim(detname)
        if(detname(j-3:j).eq.'ford') detnr(ic,i) = 1
        if(detname(j-3:j).eq.'ston') detnr(ic,i) = 2
        if(detname(j-3:j).eq.'Pisa') detnr(ic,i) = 3
     end do
     
     parID = 0
     revID = 0
     if(outputVersion > 0.5) then
        revID = 0
        read(10,*,iostat=io)parID(1:nMCMCpar) !Read parameter IDs
        if(io.ne.0) then
           write(0,'(//,A,//)')'  Error reading MCMC parameter IDs, aborting...'
           stop
        end if
        do i=1,nMCMCpar
           revID(parID(i)) = i  !Reverse ID
        end do
     end if
     read(10,*,end=199,err=199)tmpStr  !Read line with column headers (Cycle, log Post., Prior, etc)
     
     i=1
     do while(i.le.maxIter)
        if(outputVersion < 0.5) then
           read(10,*,iostat=io)tmpInt,post(ic,i),tmpDat(1:nMCMCpar)
        else
           read(10,*,iostat=io)tmpInt,post(ic,i),prior(ic,i),tmpDat(1:nMCMCpar)
        end if
        is(ic,i) = real(tmpInt)
        
        if(io.lt.0) exit !EOF
        if(io.gt.0) then !Read error
           if(readerror.eq.1) then !Read error in previous line as well
              if(i.lt.25) then
                 write(0,'(A,I7,$)')'  Read error in file '//trim(infile)//', line',i
                 write(0,'(A,/)')'  Aborting program...'
                 stop
              else
                 write(6,'(A,I7,$)')'  Read error in file '//trim(infile)//', line',i
                 write(6,*)
                 i = i-1
                 exit
              end if
           end if
           readerror = 1
           i = i-1
           cycle
        end if
        readerror = 0
        
        
        !GPS time doesn't fit in single-precision variable
        if(ic.eq.1.and.i.eq.1) then
           dtmpDat = 0.d0
           do p=1,nMCMCpar
              if(parID(p).ge.11.and.parID(p).le.19) then
                 dtmpDat(p) = dble(floor(tmpDat(p)/10.d0)*10)  !'GPS base time', rounded off to the nearest 10s, 
                 t0 = dtmpDat(p)                               !  to allow x.xx labels on the plots for GPS-t0
                 GPStime = floor(tmpDat(p)+0.05)               !Always floor, unless >x.95s. Use as 'name' to refer to this event, e.g. in file names
              end if
           end do
        end if
        allDat(ic,1:nMCMCpar,i) = real(tmpDat(1:nMCMCpar) - dtmpDat(1:nMCMCpar))
        
        if(thin.gt.1.and.i.gt.2) then !'Thin' the output by reading every thin-th line
           do j=1,thin-1
              read(10,*,iostat=io)tmpStr
              if(io.lt.0) exit !EOF
           end do
        end if
        
        i = i+1
        if(tmpInt.ge.maxchlen) exit
     end do !i
     
     if(i.ge.maxIter-2) write(0,'(A,$)')'   *** WARNING ***   Not all lines in this file were read    '
     goto 199
199  close(10)
     ntot(ic) = i-1
     n(ic) = ntot(ic) !n can be changed in rearranging chains, ntot won't be changed
     !if(prprogress.ge.2.and.update.ne.1) write(6,'(1x,3(A,I9),A1)')' Lines:',ntot(ic),', iterations:',nint(is(ic,ntot(ic))),', burn-in:',nburn(ic),'.'
  end do !do ic = 1,nchains0
  

end subroutine read_mcmcfiles
!************************************************************************************************************************************






!************************************************************************************************************************************
subroutine mcmcruninfo(exitcode)  !Extract info from the chains and print some of it to screen:
  !  print MCMC run info,  determine Lmax, burnin,  print chain info,  determine thinning for chains plots,  change/add MCMC parameters, 
  !  determine true, start, Lmax values of chains,  compute jumps,  construct output file name,  store data in selDat (from dat)
  use constants
  use analysemcmc_settings
  use general_data
  use mcmcrun_data
  use chain_data
  use plot_data
  implicit none
  integer :: i,ic,j,p,exitcode,maxLine
  real :: avgtotthin
  real*8 :: lon2ra,gmst
  
  exitcode = 0
  totiter = 0
  do ic = 1,nchains0
     totiter = totiter + nint(is(ic,ntot(ic)))
  end do
  !if(prprogress.ge.2.and.update.ne.1) write(6,'(A10,65x,2(A,I9),A1)')'Total:',' Lines:',sum(ntot(1:nchains0)),', iterations:',totiter
  
  
  
  !Print run info (detectors, SNR, amount of data, FFT, etc)
  if(prruninfo.gt.0.and.update.eq.0) then
     if(prruninfo.eq.1) write(6,'(/,A)')'  Run information for chain 1:'
     if(prruninfo.eq.2) write(6,'(/,A)')'  Run information:'
     do ic = 1,nchains0
        if((prruninfo.eq.1.and.ic.eq.1) .or. prruninfo.eq.2) then
           write(6,'(4x,A7,A29,A13,A10,A12,A8,A8)')'Chain','file name','colour','Niter','Nburn','seed','Ndet'
           write(6,'(4x,I7,A29,A13,I10,I12,I8,I8)')ic,trim(infiles(ic)),trim(colournames(colours(mod(ic-1,ncolours)+1))),niter(ic),nburn0(ic),seed(ic),ndet(ic)
           write(6,'(A14,A3,A18,4A12,A22,A17,3A14)')'Detector','Nr','SNR','f_low','f_high','before tc','after tc','Sample start (GPS)','Sample length','Sample rate','Sample size','FT size'
           
           do i=1,ndet(ic)
              write(6,'(A14,I3,F18.8,4F12.2,F22.8,F17.7,3I14)')trim(detnames(ic,i)),detnr(ic,i),snr(ic,i),flow(ic,i),fhigh(ic,i),t_before(ic,i),t_after(ic,i),FTstart(ic,i),deltaFT(ic,i),samplerate(ic,i),samplesize(ic,i),FTsize(ic,i)
           end do
           write(6,*)
        end if
     end do !do ic = 1,nchains0
  end if  !prruninfo.gt.0
  
  maxLine = maxval(n(1:nchains0))
  
  
  
  !*** Until now, nburn is the iteration number.
  !*** From here on, nburn is the line number, while isburn is the iteration number
  do ic=1,nchains0
     if(nburn(ic).le.0) nburn(ic) = nburn0(ic)
     if(abs(nburnfrac).gt.1.e-4.and.abs(nburnfrac).lt.1.) then
        nburn(ic) = is(ic,n(ic)) * abs(nburnfrac)
     else
        if(nburn(ic).ge.nint(is(ic,n(ic)))) then
           !print*,nburn(ic),nint(is(ic,n(ic)))
           if(nburn0(ic).ge.nint(is(ic,n(ic)))) then
              write(0,'(A,I3)')'   *** WARNING ***  Nburn larger than Nchain, setting nburn to 10% for chain',ic
              nburn(ic) = nint(is(ic,n(ic))*0.1)
           else
              nburn(ic) = nburn0(ic)
           end if
        end if
     end if
  end do
  
  do ic=1,nchains0
     isburn(ic) = real(nburn(ic))
     do i=1,ntot(ic)
        if(is(ic,i).le.isburn(ic)) nburn(ic) = i   !isburn is the true iteration number at which the burnin ends
        totthin(ic) = nint(isburn(ic)/real(nburn(ic)))
     end do
  end do
  avgtotthin = sum(isburn(1:nchains0))/real(sum(nburn(1:nchains0))) !Total thinning, averaged over all chains
  
  
  
  
  !Get point with absolute maximum likelihood over all chains
  loglmax = -1.d99
  loglmaxs = -1.d99
  do ic=1,nchains0
     do i=3,ntot(ic)  !3: exclude true and starting values
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
  
  if(prruninfo.ge.1) then
     ic = icloglmax
     i = iloglmax
     write(6,'(A,I4,2(A,I9),A,F10.3,A,F7.2,A)')'    Maximum likelihood point:   chain:',ic,' ('//trim(infiles(ic))//'),   line:',i*thin+7+ndet(ic), &
          ',   iteration:',nint(is(ic,i)),',   max log(L):',loglmax,'  -> SNR:',sqrt(2*loglmax),'.'
     
     !Test: get parameter values for L=Lmax
     !if(prprogress.ge.3) then
     !   write(6,'(I10,F14.6,1x,2F12.7,F20.8,9F12.7)')nint(is(ic,i)),loglmax,allDat(ic,2:3,i),allDat(ic,4,i)+t0,allDat(ic,5:13,i)
     !   call mc_eta_2_m1_m2r(allDat(ic,2,i),allDat(ic,3,i), m1,m2)
     !   write(6,'(F9.4,F10.4,F21.6,F11.4,F10.4,F12.4,2F11.4,F12.4,3F11.4)')m1,m2,allDat(ic,4,i)+t0,exp(allDat(ic,5,i)),allDat(ic,6,i),acos(allDat(ic,7,i))*r2d,allDat(ic,8,i)*r2h,asin(allDat(ic,9,i))*r2d,allDat(ic,10,i)*r2d,asin(allDat(ic,11,i))*r2d,allDat(ic,12,i)*r2d,allDat(ic,13,i)*r2d
     !end if
     if(prruninfo.ge.1) write(6,*)
  end if
  
  
  !***Autoburnin: for each chain, get the first point where log(L) > log(L_max)-autoburnin
  if(autoburnin.gt.1.e-10) then
     loop1: do ic=1,nchains0
        isburn(ic) = is(ic,ntot(ic)) !Set burnin to last iteration, so that chain is completely excluded if condition is never fulfilled
        nburn(ic) = ntot(ic)
        do i=2,ntot(ic) !i=1 is true value?
           if(post(ic,i).gt.real(loglmax)-autoburnin) then
              isburn(ic) = is(ic,i)
              nburn(ic) = i
              cycle loop1
           end if
        end do
     end do loop1
  end if
  
  
  
  
  !*** Print chain info to screen
  !Print info on number of iterations, burnin, thinning, etc.
  do ic=1,nchains0
     if(prchaininfo.ge.2.and.update.ne.1) write(6,'(A9,I3,A2,A23,A12,A,ES7.1,A,ES7.1,A,ES7.1,A,ES7.1, A,F8.2, A,I3,A,I4,A,ES8.2,A1)') &
          'Chain',ic,': ',trim(infiles(ic)),', '//colournames(colours(mod(ic-1,ncolours)+1))//'.', '  Lines/iter: ',real(n(ic)),'/',is(ic,n(ic)), &
          '.  Burn-in: ',real(nburn(ic)),'/',isburn(ic),'.  Lmax:',loglmaxs(ic),'.  Thin: file:',nint(is(ic,n(ic))/real(n(ic)*max(thin,1))),', tot:',totthin(ic),'.  Data pts: ',real(n(ic)-nburn(ic)),'.'
  end do
  totiter = 0
  totpts  = 0
  contrchains = 0
  contrchain = 0
  do ic=1,nchains0
     totiter = totiter + nint(is(ic,ntot(ic)))
     totpts = totpts + n(ic)-nburn(ic)
     if(n(ic).gt.nburn(ic)) then
        contrchains = contrchains + 1
        contrchain(ic) = 1
     end if
  end do
  totlines = sum(ntot(1:nchains0))
  if(prchaininfo.ge.1.and.update.ne.1) write(6,'(4x,A, A,ES10.4, A,ES10.4, A,I4, A,ES10.4,  A2,F5.1, A,I3,A1,I2,A1)') 'In all chains:','  # lines: ',real(totlines), &
       ',  # iterations: ',real(totiter),',  total thinning:',nint(avgtotthin),'x, # data points after burnin: ',real(totpts),' (',real(totpts)/real(totlines)*100,'%), contributing chains:',contrchains,'/',nchains0,'.'
  
  
  
  
  !*** Determine extra thinning for logL, chain, jump plots
  if(chainpli.le.0) then
     !if(sum(ntot(1:nchains0)).gt.maxdots) then  !Change the number of points plotted in chains,logL, etc. (For all output formats)
     chainpli = max(1,nint(real(sum(ntot(1:nchains0)))/real(maxdots)))  !Use ntot and nchains0, since n is low if many points are in the burnin
     if(prchaininfo.ge.1.and.update.eq.0) then
        if(chainpli.gt.1) then  !Change the number of points plotted in chains,logL, etc. (For all output formats)
           write(6,'(A,I4,A,I5,A,I5,A)')'    Plotting every',chainpli,'-th state in likelihood, chains, jumps, etc. plots.  Average total thinning is',nint(avgtotthin),'x, for these plots it is',nint(avgtotthin*chainpli),'x.'
        else
           write(6,'(A,I4,A)')'    Plotting *every* state in likelihood, chains, jumps, etc. plots.  Average total thinning remains',nint(avgtotthin),'x for these plots.'
        end if
     end if
     write(6,*)
  end if
  !if(prruninfo.gt.0) write(6,*)
  
  
  
  
  !*** Change some MCMC parameters:
  if(changevar.ge.1) then
     if(prprogress.ge.2.and.update.eq.0) write(6,'(A,$)')'  Changing some variables...   '
     
     !If Mc and eta are know, compute M1,M2:
     do ic=1,nchains0
        if(revID(61).ne.0 .and. revID(62).ne.0 .and. revID(63).eq.0 .and. revID(64).eq.0) then  !Calculate the individual masses from Mch and eta:
           if(ic.eq.1.and.prprogress.ge.2.and.update.eq.0) write(6,'(A)')'  Computing M1, M2 from Mc, eta'
           parID(nMCMCpar+1) = 63  !M1
           parID(nMCMCpar+2) = 64  !M2
           revID(63) = nMCMCpar + 1  !M1
           revID(64) = nMCMCpar + 2  !M2
           nMCMCpar = nMCMCpar + 2
           if(nMCMCpar.gt.maxMCMCpar) then
              write(0,'(//,A,I4,A,I4,A,//)')'  Error:  maxMCMCpar too small.  You must increase maxMCMCpar from',maxMCMCpar,' to at least',nMCMCpar,' in order to continue.  Aborting...'
              stop
           end if
           do i=1,ntot(ic)
              call mc_eta_2_m1_m2r(allDat(ic,revID(61),i),allDat(ic,revID(62),i), allDat(ic,revID(63),i),allDat(ic,revID(64),i))
           end do
        else if(revID(61).eq.0 .and. revID(62).eq.0 .and. revID(63).ne.0 .and. revID(64).ne.0) then  !Calculate Mc, eta from the individual masses:
           !Calculate Mch and eta from the individual masses:
           if(ic.eq.1.and.prprogress.ge.2.and.update.eq.0) write(6,'(A)')'  Computing Mc, eta from M1, M2'
           parID(nMCMCpar+1) = 61  !Mc
           parID(nMCMCpar+2) = 62  !eta
           revID(61) = nMCMCpar + 1  !Mc
           revID(62) = nMCMCpar + 2  !eta
           nMCMCpar = nMCMCpar + 2
           if(nMCMCpar.gt.maxMCMCpar) then
              write(0,'(//,A,I4,A,I4,A,//)')'  Error:  maxMCMCpar too small.  You must increase maxMCMCpar from',maxMCMCpar,' to at least',nMCMCpar,' in order to continue.  Aborting...'
              stop
           end if
           do i=1,ntot(ic)
              call m1_m2_2_mc_etar(allDat(ic,revID(63),i),allDat(ic,revID(64),i), allDat(ic,revID(61),i),allDat(ic,revID(62),i))
           end do
        end if
        
        !Compute inclination and polarisation angle from RA, Dec, theta_J0, phi_J0:
        do i=1,ntot(ic)
           !call compute_incli_polangr(allDat(ic,8,i), asin(allDat(ic,9,i)), real(lon2ra(dble(allDat(ic,12,i)), dble(allDat(ic,4,i)) + t0)), asin(allDat(ic,11,i)),   allDat(ic,11,i),allDat(ic,12,i))  !Input: RA, Dec, phi_Jo (hh->RA), theta_Jo (in rad), output: inclination, polarisation angle (rad)
           call compute_incli_polangr(allDat(ic,revID(31),i), asin(allDat(ic,revID(32),i)), real(lon2ra(dble(allDat(ic,revID(54),i)), dble(allDat(ic,revID(11),i)) + t0)), asin(allDat(ic,revID(53),i)),   allDat(ic,revID(53),i),allDat(ic,revID(54),i))  !Input: RA, Dec, phi_Jo (hh->RA), theta_Jo (in rad), output: inclination, polarisation angle (rad)
           allDat(ic,revID(53),i) = cos(allDat(ic,revID(53),i))    !i -> cos(i)
        end do
        parID(revID(53)) = 51  !Was sin(thJ0), now cos(i)
        parID(revID(54)) = 52  !Was phi_J0, now psi
        revID(51) = revID(53)  !Now cos(i)
        revID(52) = revID(54)  !Now psi
        revID(53) = 0          !No longer exists
        revID(54) = 0          !No longer exists
     end do
  end if
  
  
  !*** Put plot data in startval and jumps.  Print initial and starting values to screen.  Startval: 1: true value, 2: starting value, 3: Lmax value
  jumps = 0.
  offsetrun = 0
  if(prinitial.ne.0) then
     write(6,'(/,A)')'  True, starting and Lmax values for the chains:'
     write(6,'(5x,A10,$)')''
     do p=1,nMCMCpar
        write(6,'(A10,$)')trim(varnames(parID(p)))
     end do
     write(6,*)
  end if
  
  do ic=1,nchains
     startval(ic,1:nMCMCpar,1:2)  = allDat(ic,1:nMCMCpar,1:2) !True value and starting value
     startval(ic,1:nMCMCpar,3)    = allDat(icloglmax,1:nMCMCpar,iloglmax) !Lmax value
     jumps(ic,1:nMCMCpar,2:n(ic)) = allDat(ic,1:nMCMCpar,2:n(ic)) -  allDat(ic,1:nMCMCpar,1:n(ic)-1)
     if(prinitial.ne.0) then 
        if(ic.eq.1) then
           write(6,'(5x,A10,$)')'True:  '
           do p=1,nMCMCpar
              write(6,'(F10.5,$)')startval(1,p,1)
           end do
           write(6,'(/)')
        end if
        if(abs((sum(startval(ic,1:nMCMCpar,1))-sum(startval(ic,1:nMCMCpar,2)))/sum(startval(ic,1:nMCMCpar,1))).gt.1.e-10) then
           offsetrun = 1
           write(6,'(I4,A1,A10,$)')ic,':','  Start: '
           do p=1,nMCMCpar
              write(6,'(F10.5,$)')startval(ic,p,2)
           end do
           write(6,*)
           write(6,'(5x,A10,$)')'Diff:  '
           do p=1,nMCMCpar
              write(6,'(F10.5,$)')abs(startval(ic,p,1)-startval(ic,p,2))
           end do
           write(6,'(/)')
        end if
     end if
  end do
  if(prinitial.ne.0) then
     write(6,'(5x,A10,$)')'Lmax:  '
     do p=1,nMCMCpar
        write(6,'(F10.5,$)')startval(1,p,3)
     end do
     write(6,*)
     write(6,'(5x,A10,$)')'Diff:  '
     do p=1,nMCMCpar
        write(6,'(F10.5,$)')abs(startval(1,p,1)-startval(1,p,3))
     end do
     write(6,'(/)')
  end if
  
  
  !if(prprogress.ge.1.and.update.eq.0) write(6,'(A)')'  Done.'
  if(prprogress.ge.2.and.update.eq.0) write(6,'(A,I12,A,F7.4)')'  t0:',nint(t0), '  GMST:',gmst(t0)
  
  
  !Check which parameters were fixed during the MCMC run
  fixedpar = 0
  do ic=1,nchains
     do p=1,nMCMCpar
        if( abs(minval(allDat(ic,p,5:n(ic))) - maxval(allDat(ic,p,5:n(ic))) ) .lt. 1.d-6) fixedpar(p) = 1  !Doesn't matter in which chain this happens
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
  if(mergechains.eq.1) then  !Merge chains, leave out burnin (then nchains = 1)
     allocate(selDat(1,maxMCMCpar,nchains*maxLine))
     j = 1
     do ic=1,nchains
        do i=nburn(ic)+1,n(ic)
           selDat(1,1:nMCMCpar,j) = allDat(ic,1:nMCMCpar,i)  !SelDat has the same structure as allDat, but contains only info AFTER the burnin.
           j = j+1
        end do
     end do
     nchains = 1
     n(1) = j-1
     !if(prprogress.ge.1) write(6,'(A,I8,A,ES7.1,A)')'  Data points in combined chains: ',n(1),'  (',real(n(1)),')'
  else
     allocate(selDat(nchains,maxMCMCpar,maxLine))
     do ic=1,nchains
        selDat(ic,1:nMCMCpar,1:n(ic)-nburn(ic)) = allDat(ic,1:nMCMCpar,nburn(ic)+1:n(ic))  !SelDat has the same structure as allDat, but contains only info AFTER the burnin.
        n(ic) = n(ic)-nburn(ic) !n(ic)=0 if a chain is not contributing (in which case contrchain(ic)=0)!
     end do
     !if(prprogress.ge.1) write(6,'(A,I8)')' Datapoints in combined chains: ',sum(n(1:nchains))
  end if
  
  
end subroutine mcmcruninfo
!************************************************************************************************************************************





!> \brief Set the names and symbols of the original MCMC parameters
!<
!************************************************************************************************************************************
subroutine set_originalParameterNames()
  use analysemcmc_settings
  use general_data
  implicit none
  
  varnames = ''
  pgvarns = ''
  pgvarnss = ''
  pgunits = ''
  
  !Short ASCII names for text output:
  varnames(11:19) = (/'tc','t40','','','','','','',''/)
  varnames(21:29) = (/'dl^3','log_dl','','','','','','',''/)
  varnames(31:39) = (/'RA','sin_dec','','','','','','',''/)
  varnames(41:49) = (/'phase','','','','','','','',''/)
  varnames(51:59) = (/'cos_i','psi','sin_thJo','ph_Jo','','','','',''/)
  varnames(61:69) = (/'Mc','eta','M1','M2','','','','',''/)
  varnames(71:79) = (/'spin1','cos_th1','phi1','','','','','',''/)
  varnames(81:89) = (/'spin2','cos_th2','phi2','','','','','',''/)
  !varnames(1:9) = (/'','','','','','','','',''/)
  
  
  if(fonttype.eq.2) then  !Use 'roman-like' Greek font in PGPlot
     
     !Long PGPlot names (symbol + unit)
     pgvarns(11:19) = (/'t\dc\u (s)','t\d40\u (s)','','','','','','',''/)
     pgvarns(21:29) = (/'d\dL\u\u3\d (Mpc)','logd\dL\u (Mpc)','','','','','','',''/)
     pgvarns(31:39) = (/'\(2127) (rad)','sin \(2130)','','','','','','',''/)
     pgvarns(41:49) = (/'\(2147)\dc\u (rad)','','','','','','','',''/)
     pgvarns(51:59) = (/'cos \(2135)','\(2149) (rad)','sin \(2185)\dJ0\u','\(2147)\dJ0\u (rad)','','','','',''/)
     pgvarns(61:69) = (/'\(2563) (M\d\(2281)\u)','\(2133)','M\d1\u (M\d\(2281)\u)','M\d2\u (M\d\(2281)\u)','','','','',''/)
     pgvarns(71:79) = (/'a\dspin1\u','cos \(2185)\dspin1\u','\(2147)\dspin1\u (rad)','','','','','',''/)
     pgvarns(81:89) = (/'a\dspin2\u','cos \(2185)\dspin2\u','\(2147)\dspin2\u (rad)','','','','','',''/)
     !pgvarns(1:9) = (/'','','','','','','','',''/)
     
     !Short PGPlot symbols (no unit)
     pgvarnss(11:19) = (/'t\dc\u','t\d40\u','','','','','','',''/)
     pgvarnss(21:29) = (/'d\dL\u\u3\d','logd\dL\u','','','','','','',''/)
     pgvarnss(31:39) = (/'\(2127)','sin \(2130)','','','','','','',''/)
     pgvarnss(41:49) = (/'\(2147)\dc\u','','','','','','','',''/)
     pgvarnss(51:59) = (/'cos \(2135)','\(2149)','sin \(2185)\dJ0\u','\(2147)\dJ0\u','','','','',''/)
     pgvarnss(61:69) = (/'\(2563)','\(2133)','M\d1\u','M\d2\u','','','','',''/)
     pgvarnss(71:79) = (/'a\dspin1\u','cos \(2185)\dspin1\u','\(2147)\dspin1\u','','','','','',''/)
     pgvarnss(81:89) = (/'a\dspin2\u','cos \(2185)\dspin2\u','\(2147)\dspin2\u','','','','','',''/)
     !pgvarnss(1:9) = (/'','','','','','','','',''/)
     
  else  !Same, but replace '\(21' with \(06' for arial-like Greek font
     
     !Long PGPlot names (symbol + unit)
     pgvarns(11:19) = (/'t\dc\u (s)','t\d40\u (s)','','','','','','',''/)
     pgvarns(21:29) = (/'d\dL\u\u3\d (Mpc)','logd\dL\u (Mpc)','','','','','','',''/)
     pgvarns(31:39) = (/'\(0627) (rad)','sin \(0630)','','','','','','',''/)
     pgvarns(41:49) = (/'\(0647)\dc\u (rad)','','','','','','','',''/)
     pgvarns(51:59) = (/'cos \(0635)','\(0649) (rad)','sin \(0685)\dJ0\u','\(0647)\dJ0\u (rad)','','','','',''/)
     pgvarns(61:69) = (/'\(2563) (M\d\(2281)\u)','\(0633)','M\d1\u (M\d\(2281)\u)','M\d2\u (M\d\(2281)\u)','','','','',''/)
     pgvarns(71:79) = (/'a\dspin1\u','cos \(0685)\dspin1\u','\(0647)\dspin1\u (rad)','','','','','',''/)
     pgvarns(81:89) = (/'a\dspin2\u','cos \(0685)\dspin2\u','\(0647)\dspin2\u (rad)','','','','','',''/)
     !pgvarns(1:9) = (/'','','','','','','','',''/)
     
     !Short PGPlot symbols (no unit)
     pgvarnss(11:19) = (/'t\dc\u','t\d40\u','','','','','','',''/)
     pgvarnss(21:29) = (/'d\dL\u\u3\d','logd\dL\u','','','','','','',''/)
     pgvarnss(31:39) = (/'\(0627)','sin \(0630)','','','','','','',''/)
     pgvarnss(41:49) = (/'\(0647)\dc\u','','','','','','','',''/)
     pgvarnss(51:59) = (/'cos \(0635)','\(0649)','sin \(0685)\dJ0\u','\(0647)\dJ0\u','','','','',''/)
     pgvarnss(61:69) = (/'\(2563)','\(0633)','M\d1\u','M\d2\u','','','','',''/)
     pgvarnss(71:79) = (/'a\dspin1\u','cos \(0685)\dspin1\u','\(0647)\dspin1\u','','','','','',''/)
     pgvarnss(81:89) = (/'a\dspin2\u','cos \(0685)\dspin2\u','\(0647)\dspin2\u','','','','','',''/)
     !pgvarnss(1:9) = (/'','','','','','','','',''/)
     
  end if
  
  !PGPlot units (no names)
  pgunits(11:19) = (/'s','s','','','','','','',''/)
  pgunits(21:29) = (/'Mpc','Mpc','','','','','','',''/)
  pgunits(31:39) = (/'rad','','','','','','','',''/)
  pgunits(41:49) = (/'rad','','','','','','','',''/)
  pgunits(51:59) = (/'','rad','','rad','','','','',''/)
  pgunits(61:69) = (/'M\d\(2281)\u','','M\d\(2281)\u','M\d\(2281)\u','','','','',''/)
  pgunits(71:79) = (/'','','rad','','','','','',''/)
  pgunits(81:89) = (/'','','rad','','','','','',''/)
  !pgunits(1:9) = (/'','','','','','','','',''/)
  
     
  !Save the original parameter names for use after they get changed
  pgorigvarns = pgvarns
  
  
  
  
end subroutine set_originalParameterNames
!************************************************************************************************************************************






!> \brief Set the names and symbols of the derived MCMC parameters
!! e.g. d_L rather than d_L^3 or log(d_L), i rather than cos(i), etc.
!<
!************************************************************************************************************************************
subroutine set_derivedParameterNames()
  use analysemcmc_settings
  use general_data
  implicit none
  
  varnames = ''
  pgvarns = ''
  pgvarnss = ''
  pgunits = ''
  
  !Short ASCII names for text output:
  varnames(11:19) = (/'tc','t40','','','','','','',''/)
  varnames(21:29) = (/'dl','dl','','','','','','',''/)
  varnames(31:39) = (/'RA','dec','','','','','','',''/)
  varnames(41:49) = (/'phase','','','','','','','',''/)
  varnames(51:59) = (/'incl','psi','th_Jo','ph_Jo','','','','',''/)
  varnames(61:69) = (/'Mc','eta','M1','M2','','','','',''/)
  varnames(71:79) = (/'spin1','th1','phi1','','','','','',''/)
  varnames(81:89) = (/'spin2','th2','phi2','','','','','',''/)
  !varnames(1:9) = (/'','','','','','','','',''/)
  
  
  if(fonttype.eq.2) then  !Use 'roman-like' Greek font in PGPlot
     
     !Long PGPlot names (symbol + unit)
     pgvarns(11:19) = (/'t\dc\u (s)','t\d40\u (s)','','','','','','',''/)
     pgvarns(21:29) = (/'d\dL\u (Mpc)','d\dL\u (Mpc)','','','','','','',''/)
     pgvarns(31:39) = (/'\(2127) (h)','\(2130) (\(2218))','','','','','','',''/)
     pgvarns(41:49) = (/'\(2147)\dc\u (\(2218))','','','','','','','',''/)
     pgvarns(51:59) = (/'\(2135) (\(2218))','\(2149) (\(2218))','\(2185)\dJ0\u (\(2218))','\(2147)\dJ0\u (\(2218))','','','','',''/)
     pgvarns(61:69) = (/'\(2563) (M\d\(2281)\u)','\(2133)','M\d1\u (M\d\(2281)\u)','M\d2\u (M\d\(2281)\u)','','','','',''/)
     pgvarns(71:79) = (/'a\dspin1\u','\(2185)\dspin1\u (\(2218))','\(2147)\dspin1\u (\(2218))','','','','','',''/)
     pgvarns(81:89) = (/'a\dspin2\u','\(2185)\dspin2\u (\(2218))','\(2147)\dspin2\u (\(2218))','','','','','',''/)
     !pgvarns(1:9) = (/'','','','','','','','',''/)
     
     !Short PGPlot symbols (no unit)
     pgvarnss(11:19) = (/'t\dc\u','t\d40\u','','','','','','',''/)
     pgvarnss(21:29) = (/'d\dL\u\u3\d','logd\dL\u','','','','','','',''/)
     pgvarnss(31:39) = (/'\(2127)','\(2130)','','','','','','',''/)
     pgvarnss(41:49) = (/'\(2147)\dc\u','','','','','','','',''/)
     pgvarnss(51:59) = (/'\(2135)','\(2149)','\(2185)\dJ0\u','\(2147)\dJ0\u','','','','',''/)
     pgvarnss(61:69) = (/'\(2563)','\(2133)','M\d1\u','M\d2\u','','','','',''/)
     pgvarnss(71:79) = (/'a\dspin1\u','\(2185)\dspin1\u','\(2147)\dspin1\u','','','','','',''/)
     pgvarnss(81:89) = (/'a\dspin2\u','\(2185)\dspin2\u','\(2147)\dspin2\u','','','','','',''/)
     !pgvarnss(1:9) = (/'','','','','','','','',''/)
     
  else  !Same, but replace '\(21' with \(06' for arial-like Greek font
     
     !Long PGPlot names (symbol + unit)
     pgvarns(11:19) = (/'t\dc\u (s)','t\d40\u (s)','','','','','','',''/)
     pgvarns(21:29) = (/'d\dL\u (Mpc)','d\dL\u (Mpc)','','','','','','',''/)
     pgvarns(31:39) = (/'\(0627) (h)','\(0630) (\(2218))','','','','','','',''/)
     pgvarns(41:49) = (/'\(0647)\dc\u (\(2218))','','','','','','','',''/)
     pgvarns(51:59) = (/'\(0635) (\(2218))','\(0649) (\(2218))','\(0685)\dJ0\u (\(2218))','\(0647)\dJ0\u (\(2218))','','','','',''/)
     pgvarns(61:69) = (/'\(2563) (M\d\(2281)\u)','\(0633)','M\d1\u (M\d\(2281)\u)','M\d2\u (M\d\(2281)\u)','','','','',''/)
     pgvarns(71:79) = (/'a\dspin1\u','\(0685)\dspin1\u (\(2218))','\(0647)\dspin1\u (\(2218))','','','','','',''/)
     pgvarns(81:89) = (/'a\dspin2\u','\(0685)\dspin2\u (\(2218))','\(0647)\dspin2\u (\(2218))','','','','','',''/)
     !pgvarns(1:9) = (/'','','','','','','','',''/)
     
     !Short PGPlot symbols (no unit)
     pgvarnss(11:19) = (/'t\dc\u','t\d40\u','','','','','','',''/)
     pgvarnss(21:29) = (/'d\dL\u','d\dL\u','','','','','','',''/)
     pgvarnss(31:39) = (/'\(0627)','\(0630)','','','','','','',''/)
     pgvarnss(41:49) = (/'\(0647)\dc\u','','','','','','','',''/)
     pgvarnss(51:59) = (/'\(0635)','\(0649)','\(0685)\dJ0\u','\(0647)\dJ0\u','','','','',''/)
     pgvarnss(61:69) = (/'\(2563)','\(0633)','M\d1\u','M\d2\u','','','','',''/)
     pgvarnss(71:79) = (/'a\dspin1\u','\(0685)\dspin1\u','\(0647)\dspin1\u','','','','','',''/)
     pgvarnss(81:89) = (/'a\dspin2\u','\(0685)\dspin2\u','\(0647)\dspin2\u','','','','','',''/)
     !pgvarnss(1:9) = (/'','','','','','','','',''/)
     
  end if
  
  !PGPlot units (no names)
  pgunits(11:19) = (/'s','s','','','','','','',''/)
  pgunits(21:29) = (/'Mpc','Mpc','','','','','','',''/)
  pgunits(31:39) = (/'\uh\d','\(2218)','','','','','','',''/)
  pgunits(41:49) = (/'\(2218)','','','','','','','',''/)
  pgunits(51:59) = (/'\(2218)','\(2218)','\(2218)','\(2218)','','','','',''/)
  pgunits(61:69) = (/'M\d\(2281)\u','','M\d\(2281)\u','M\d\(2281)\u','','','','',''/)
  pgunits(71:79) = (/'','\(2218)','\(2218)','','','','','',''/)
  pgunits(81:89) = (/'','\(2218)','\(2218)','','','','','',''/)
  !pgunits(1:9) = (/'','','','','','','','',''/)
  
     
end subroutine set_derivedParameterNames
!************************************************************************************************************************************






!************************************************************************************************************************************
function lon2ra(lon, GPSsec)
  ! Compute right ascension (in radians) from longitude (radians) and GPS time (seconds). 
  ! Declination == latitude for equatorial coordinates.
  
  use constants
  implicit none
  real*8 :: lon2ra,lon,GPSsec,gmst
  
  lon2ra = mod(lon + gmst(GPSsec) + 10*tpi,tpi)
end function lon2ra
!************************************************************************************************************************************


!************************************************************************************************************************************
function gmst(GPSsec)
  ! Compute the 'Greenwich Mean Sidereal Time' (in radians) from GPS time (in seconds).                              
  ! See K.R. Lang (1999), p.80sqq.
  
  use constants
  implicit none
  real*8 :: gmst,seconds,days,centuries,secCurrentDay
  real*8 :: gps0,leapseconds,GPSsec
  
  gps0 = 630720013.d0 !GPS time at 1/1/2000 at midnight
  leapseconds = 32.d0 !At Jan 1st 2000
  if(GPSsec.gt.820108813.d0) leapseconds = leapseconds + 1.d0 !One more leapsecond after 1/1/2006
  if(GPSsec.gt.914803214.d0) leapseconds = leapseconds + 1.d0 !One more leapsecond after 1/1/2009
  if(GPSsec.lt.630720013.d0) write(0,'(A)')'   WARNING: GMSTs before 1.1.2000 are inaccurate!'
  !Time since 1/1/2000 midnight
  seconds       = (GPSsec - gps0) + (leapseconds - 32.d0)
  days          = floor(seconds/86400.d0) - 0.5d0
  secCurrentDay = mod(seconds, 86400.d0)
  centuries     = days/36525.d0
  gmst = 24110.54841d0 + (centuries*(8640184.812866d0 + centuries*(0.093104d0 + centuries*6.2d-6)))
  gmst = gmst + secCurrentDay * 1.002737909350795d0   !UTC day is 1.002 * MST day
  gmst = mod(gmst/86400.d0,1.d0)
  gmst = gmst * tpi
end function gmst
!************************************************************************************************************************************



!************************************************************************************************************************************
subroutine dindexx(n,arr,indx)
  implicit none
  integer, parameter :: m=7,nstack=50
  integer :: n,indx(n)
  real*8 :: arr(n),a
  integer :: i,indxt,ir,itemp,j,jstack,k,l,istack(nstack)
  
  do j=1,n
     indx(j)=j
  end do
  jstack=0
  l=1
  ir=n
1 if(ir-l.lt.m) then
     do j=l+1,ir
        indxt=indx(j)
        a=arr(indxt)
        do i=j-1,l,-1
           if(arr(indx(i)).le.a) goto 2
           indx(i+1)=indx(i)
        end do
        i=l-1
2       indx(i+1)=indxt
     end do
     if(jstack.eq.0) return
     ir=istack(jstack)
     l=istack(jstack-1)
     jstack=jstack-2
  else
     k=(l+ir)/2
     itemp=indx(k)
     indx(k)=indx(l+1)
     indx(l+1)=itemp
     if(arr(indx(l)).gt.arr(indx(ir))) then
        itemp=indx(l)
        indx(l)=indx(ir)
        indx(ir)=itemp
     end if
     if(arr(indx(l+1)).gt.arr(indx(ir))) then
        itemp=indx(l+1)
        indx(l+1)=indx(ir)
        indx(ir)=itemp
     end if
     if(arr(indx(l)).gt.arr(indx(l+1))) then
        itemp=indx(l)
        indx(l)=indx(l+1)
        indx(l+1)=itemp
     end if
     i=l+1
     j=ir
     indxt=indx(l+1)
     a=arr(indxt)
3    continue
     i=i+1
     if(arr(indx(i)).lt.a) goto 3
4    continue
     j=j-1
     if(arr(indx(j)).gt.a) goto 4
     if(j.lt.i) goto 5
     itemp=indx(i)
     indx(i)=indx(j)
     indx(j)=itemp
     goto 3
5    indx(l+1)=indx(j)
     indx(j)=indxt
     jstack=jstack+2
     !if(jstack.gt.nstack)pause 'nstack too small in indexx'
     if(jstack.gt.nstack) write(0,'(A)')' nstack too small in dindexx'
     if(ir-i+1.ge.j-l) then
        istack(jstack)=ir
        istack(jstack-1)=i
        ir=j-1
     else
        istack(jstack)=j-1
        istack(jstack-1)=l
        l=i
     end if
  end if
  goto 1
end subroutine dindexx
!************************************************************************************************************************************


!************************************************************************************************************************************
subroutine rindexx(n,arr,indx)
  integer :: n,indx(n),m,nstack
  real :: arr(n),a
  parameter (m=7,nstack=50)
  integer :: i,indxt,ir,itemp,j,jstack,k,l,istack(nstack)
  do j=1,n
     indx(j)=j
  end do
  jstack=0
  l=1
  ir=n
1 if(ir-l.lt.m) then
     do j=l+1,ir
        indxt=indx(j)
        a=arr(indxt)
        do i=j-1,l,-1
           if(arr(indx(i)).le.a) goto 2
           indx(i+1)=indx(i)
        end do
        i=l-1
2       indx(i+1)=indxt
     end do
     if(jstack.eq.0) return
     ir=istack(jstack)
     l=istack(jstack-1)
     jstack=jstack-2
  else
     k=(l+ir)/2
     itemp=indx(k)
     indx(k)=indx(l+1)
     indx(l+1)=itemp
     if(arr(indx(l)).gt.arr(indx(ir))) then
        itemp=indx(l)
        indx(l)=indx(ir)
        indx(ir)=itemp
     end if
     if(arr(indx(l+1)).gt.arr(indx(ir))) then
        itemp=indx(l+1)
        indx(l+1)=indx(ir)
        indx(ir)=itemp
     end if
     if(arr(indx(l)).gt.arr(indx(l+1))) then
        itemp=indx(l)
        indx(l)=indx(l+1)
        indx(l+1)=itemp
     end if
     i=l+1
     j=ir
     indxt=indx(l+1)
     a=arr(indxt)
3    continue
     i=i+1
     if(arr(indx(i)).lt.a) goto 3
4    continue
     j=j-1
     if(arr(indx(j)).gt.a) goto 4
     if(j.lt.i) goto 5
     itemp=indx(i)
     indx(i)=indx(j)
     indx(j)=itemp
     goto 3
5    indx(l+1)=indx(j)
     indx(j)=indxt
     jstack=jstack+2
     !if(jstack.gt.nstack)pause 'nstack too small in indexx'
     if(jstack.gt.nstack) write(0,'(A)')' nstack too small in rindexx'
     if(ir-i+1.ge.j-l) then
        istack(jstack)=ir
        istack(jstack-1)=i
        ir=j-1
     else
        istack(jstack)=j-1
        istack(jstack-1)=l
        l=i
     end if
  end if
  goto 1
end subroutine rindexx
!************************************************************************************************************************************


!************************************************************************************************************************************
subroutine savgol(c,np,nl,nr,ld,m)
  integer :: ld,m,nl,np,nr,mmax
  real :: c(np)
  parameter (mmax=6)
  !Uses lubksb,ludcmp
  integer :: imj,ipj,j,k,kk,mm,indx(mmax+1)
  real :: d,fac,sum,a(mmax+1,mmax+1),b(mmax+1)
  
  !if(np.lt.nl+nr+1.or.nl.lt.0.or.nr.lt.0.or.ld.gt.m.or.m.gt.mmax.or.nl+nr.lt.m) pause 'bad args in savgol'
  if(np.lt.nl+nr+1.or.nl.lt.0.or.nr.lt.0.or.ld.gt.m.or.m.gt.mmax.or.nl+nr.lt.m) write(0,'(A)')' Bad args in savgol'
  do ipj=0,2*m
     sum=0.
     if(ipj.eq.0)sum=1.
     do k=1,nr
        sum=sum+float(k)**ipj
     end do
     do k=1,nl
        sum=sum+float(-k)**ipj
     end do
     mm=min(ipj,2*m-ipj)
     do imj=-mm,mm,2
        a(1+(ipj+imj)/2,1+(ipj-imj)/2)=sum
     end do
  end do
  call ludcmp(a,m+1,mmax+1,indx,d)
  do j=1,m+1
     b(j)=0.
  end do
  b(ld+1)=1.
  call lubksb(a,m+1,mmax+1,indx,b)
  do kk=1,np
     c(kk)=0.
  end do
  do k=-nl,nr
     sum=b(1)
     fac=1.
     do mm=1,m
        fac=fac*k
        sum=sum+b(mm+1)*fac
     end do
     kk=mod(np-k,np)+1
     c(kk)=sum
  end do
  return
end subroutine savgol
!************************************************************************************************************************************


!************************************************************************************************************************************
subroutine lubksb(a,n,np,indx,b)
  integer :: n,np,indx(n)
  real :: a(np,np),b(n)
  integer :: i,ii,j,ll
  real :: sum
  ii=0
  do i=1,n
     ll=indx(i)
     sum=b(ll)
     b(ll)=b(i)
     if (ii.ne.0) then
        do j=ii,i-1
           sum=sum-a(i,j)*b(j)
        end do
     else if (sum.ne.0.) then
        ii=i
     end if
     b(i)=sum
  end do
  do i=n,1,-1
     sum=b(i)
     do j=i+1,n
        sum=sum-a(i,j)*b(j)
     end do
     b(i)=sum/a(i,i)
  end do
  return
end subroutine lubksb
!************************************************************************************************************************************



!************************************************************************************************************************************
subroutine ludcmp(a,n,np,indx,d)
 integer :: n,np,indx(n),nmax
 real :: d,a(np,np),tiny
 parameter (nmax=500,tiny=1.0e-20)
 integer :: i,imax,j,k
 real :: aamax,dum,sum,vv(nmax)
 d=1.
 do i=1,n
    aamax=0.
    do j=1,n
       if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
    end do
    !if (aamax.eq.0.) pause 'singular matrix in ludcmp'
    if(aamax.eq.0.) write(0,'(A)')' Singular matrix in ludcmp'
    vv(i)=1./aamax
 end do
 do j=1,n
    do i=1,j-1
       sum=a(i,j)
       do k=1,i-1
          sum=sum-a(i,k)*a(k,j)
       end do
       a(i,j)=sum
    end do
    aamax=0.
    do i=j,n
       sum=a(i,j)
       do k=1,j-1
          sum=sum-a(i,k)*a(k,j)
       end do
       a(i,j)=sum
       dum=vv(i)*abs(sum)
       if (dum.ge.aamax) then
          imax=i
          aamax=dum
       end if
    end do
    if (j.ne.imax) then
       do k=1,n
          dum=a(imax,k)
          a(imax,k)=a(j,k)
          a(j,k)=dum
       end do
       d=-d
       vv(imax)=vv(j)
    end if
    indx(j)=imax
    if(a(j,j).eq.0.)a(j,j)=tiny
    if(j.ne.n) then
       dum=1./a(j,j)
       do i=j+1,n
          a(i,j)=a(i,j)*dum
       end do
    end if
 end do
 return
end subroutine ludcmp
!************************************************************************************************************************************






!************************************************************************
function drev2pi(x)        !Returns angle in radians between 0 and 2pi (double precision)
  use constants
  real*8 :: x,drev2pi
  drev2pi = x-floor(x/(2*pi))*2*pi
end function drev2pi
!************************************************************************

!************************************************************************
!function rev(x)        !Returns angle in radians between 0 and 2pi
!  use constants
!  real :: x,rev
!  rev = x-floor(x/rpi2)*rpi2
!end function rev
!************************************************************************

!************************************************************************
function revpipi(x)      !Returns angle in radians between -pi and pi
  use constants
  real :: x,revpipi
  revpipi = x-floor(x/rpi2)*rpi2
  if(revpipi.gt.rpi) revpipi = revpipi - rpi2
end function revpipi
!************************************************************************

!************************************************************************
function rev360(x)        !Returns angle in degrees between 0 and 360
  real :: x,rev360
  rev360 = x-floor(x/(360.))*360.
end function rev360
!************************************************************************

!************************************************************************
function rev180(x)        !Returns angle in degrees between 0 and 180
  real :: x,rev180
  rev180 = x-floor(x/(180.))*180.
end function rev180
!************************************************************************

!************************************************************************
function rev24(x)        !Returns angle in hours between 0 and 24
  real :: x,rev24
  rev24 = x-floor(x/(24.))*24.
end function rev24
!************************************************************************

!************************************************************************
function rev2pi(x)        !Returns angle in radians between 0 and 2pi
  real :: x,rev2pi,pi
  pi = 4*atan(1.)
  rev2pi = x-floor(x/(2.0*pi))*2.0*pi
end function rev2pi
!************************************************************************

!************************************************************************
function drevpi(x)        !Returns angle in radians between 0 and pi - double
  use constants
  real*8 :: x,drevpi
  drevpi = x-floor(x/pi)*pi
end function drevpi
!************************************************************************

!************************************************************************
function rrevpi(x)        !Returns angle in radians between 0 and pi - real
  use constants
  real :: x,rrevpi
  rrevpi = x-floor(x/rpi)*rpi
end function rrevpi
!************************************************************************


!************************************************************************
subroutine kstwo(data1,n1,data2,n2,d,prob)  !Needs probks(), sort()
  integer :: n1,n2,j1,j2
  real*8 :: d,prob,data1(n1),data2(n2)
  real*8 :: d1,d2,dt,en1,en2,en,fn1,fn2,probks
  call sort(n1,data1)
  call sort(n2,data2)
  en1=n1
  en2=n2
  j1=1
  j2=1
  fn1=0.d0
  fn2=0.d0
  d=0.d0
1 if(j1.le.n1.and.j2.le.n2) then
     d1=data1(j1)
     d2=data2(j2)
     if(d1.le.d2) then
        fn1=j1/en1
        j1=j1+1
     end if
     if(d2.le.d1) then
        fn2=j2/en2
        j2=j2+1
     end if
     dt=dabs(fn2-fn1)
     if(dt.gt.d)d=dt
     goto 1
  end if
  en=dsqrt(en1*en2/(en1+en2))
  prob=probks((en+0.12d0+0.11d0/en)*d)
  return
end subroutine kstwo
!************************************************************************


!************************************************************************
function probks(alam)
  real*8 :: probks,alam,eps1,eps2
  parameter (eps1=1.d-3, eps2=1.d-8)
  integer :: j
  real*8 :: a2,fac,term,termbf
  a2=-2.d0*alam**2
  fac=2.d0
  probks=0.d0
  termbf=0.d0
  do j=1,100
     term=fac*dexp(a2*j**2)
     probks=probks+term
     if(dabs(term).le.eps1*termbf.or.dabs(term).le.eps2*probks) return
     fac=-fac
     termbf=dabs(term)
  end do
  probks=1.d0
  return
end function probks
!************************************************************************

!************************************************************************
subroutine sort(n,arr)
  integer :: n,m,nstack
  real*8 :: arr(n)
  parameter (m=7,nstack=50)
  integer :: i,ir,j,jstack,k,l,istack(nstack)
  real*8 :: a,temp
  jstack=0
  l=1
  ir=n
1 if(ir-l.lt.m) then
     do j=l+1,ir
        a=arr(j)
        do i=j-1,l,-1
           if(arr(i).le.a) goto 2
           arr(i+1)=arr(i)
        end do
        i=l-1
2       arr(i+1)=a
     end do
     if(jstack.eq.0) return
     ir=istack(jstack)
     l=istack(jstack-1)
     jstack=jstack-2
  else
     k=(l+ir)/2
     temp=arr(k)
     arr(k)=arr(l+1)
     arr(l+1)=temp
     if(arr(l).gt.arr(ir)) then
        temp=arr(l)
        arr(l)=arr(ir)
        arr(ir)=temp
     end if
     if(arr(l+1).gt.arr(ir)) then
        temp=arr(l+1)
        arr(l+1)=arr(ir)
        arr(ir)=temp
     end if
     if(arr(l).gt.arr(l+1)) then
        temp=arr(l)
        arr(l)=arr(l+1)
        arr(l+1)=temp
     end if
     i=l+1
     j=ir
     a=arr(l+1)
3    continue
     i=i+1
     if(arr(i).lt.a) goto 3
4    continue
     j=j-1
     if(arr(j).gt.a) goto 4
     if(j.lt.i) goto 5
     temp=arr(i)
     arr(i)=arr(j)
     arr(j)=temp
     goto 3
5    arr(l+1)=arr(j)
     arr(j)=a
     jstack=jstack+2
     !if(jstack.gt.nstack)pause 'nstack too small in sort'
     if(jstack.gt.nstack) write(0,'(A)')' nstack too small in dindexx'
     if(ir-i+1.ge.j-l) then
        istack(jstack)=ir
        istack(jstack-1)=i
        ir=j-1
     else
        istack(jstack)=j-1
        istack(jstack-1)=l
        l=i
     end if
  end if
  goto 1
end subroutine sort
!************************************************************************

!************************************************************************
function tms(a1)   !Print angle as mm:ss.s string, input in hours
  implicit none
  real*8 :: a1,a,s
  integer :: m
  character :: tms*8,mm*2,ss*4
  
  a = a1
  m = int((a)*60.d0)
  s = (a-m/60.d0)*3600.d0
  
  write(mm,'(i2.2)') m
  write(ss,'(f4.1)') s
  if(nint(s*10).lt.100) write(ss,'(a1,f3.1)') '0',s
  write(tms,'(a2,a1,a4,a1)') mm,'m',ss,'s'
  
  return
end function tms
!************************************************************************


!***********************************************************************
function getos() !Determine the operating system type: 1-Linux, 2-MacOSX
  implicit none
  integer :: i,system,getos
  character :: ostype*25,filename*99
  filename = '.analysemcmc.uname.temp'
  i = system('uname &> '//trim(filename)) !This should return "Linux" or "Darwin"
  open(16,file=trim(filename), status='old', form='formatted')
  read(16,'(A)')ostype
  close(16, status = 'delete')
  !write(6,*)ostype
  getos = 1 !Linux
  if(ostype(1:5).eq.'Darwi') getos = 2 !MacOSX
  return
end function getos
!***********************************************************************


!************************************************************************
function timestamp(os)  !Get time stamp in seconds since 1970-01-01 00:00:00 UTC
  implicit none
  real*8 :: timestamp
  integer :: os,i,system
  character :: fname*99,homedir*99
  
  call getenv('HOME',homedir)
  
  fname = trim(homedir)//'/.analysemcmc_timestamp'
  if(os.eq.1) then !Linux
     i = system('date +%s.%N >& '//trim(fname))
  else
     i = system('date +%s >& '//trim(fname)) !%N for fractional seconds doesn't work on MacOS!!! (But it does with GNU date)
  end if
  open(unit=9,status='old',file=trim(fname))
  !read(9,'(F20.9)')timestamp
  read(9,*)timestamp
  !print*,timestamp
  close(9)
  i = system('rm -f '//trim(fname))
end function timestamp
!************************************************************************



!************************************************************************
subroutine pgscidark(ci0,file,whitebg)  !Set the colour to ci, but use a darker shade if the background is black or a lighter shade if it is white
  implicit none
  integer :: ci0,ci,file,whitebg
  real :: r,g,b,weight
  ci = ci0
  call pgqcr(ci,r,g,b)
  call pgscr(99,r*0.5,g*0.5,b*0.5) !Use half the RGB value to create a darker shade
  !if(file.ge.2.or.whitebg.ge.1) call pgscr(99,(r+1)/2.,(g+1)/2.,(b+1)/2.) !Use the mean of the RGB value and 1. to create a lighter shade
  weight = 3.
  if(file.ge.2.or.whitebg.ge.1) call pgscr(99,(r+weight)/(weight+1.),(g+weight)/(weight+1.),(b+weight)/(weight+1.)) !Use the weighted mean of the RGB value and 1. to create a lighter shade
  call pgsci(99)
end subroutine pgscidark
!************************************************************************



!************************************************************************
subroutine lbr2vec(l,b,r,vec)
  !Transforms longitude l, latitude b and radius r into a vector with length r.  Use r=1 for a unit vector
  implicit none
  real*8 :: l,b,r,sinb,cosb,vec(3)
  sinb = dsin(b)
  cosb = dsqrt(1.d0-sinb*sinb)
  vec(1) = dcos(l) * cosb;  !`Greenwich'
  vec(2) = dsin(l) * cosb;  !`Ganges'
  vec(3) = sinb;            !`North Pole'
  vec = vec*r
end subroutine lbr2vec
!************************************************************************



!************************************************************************
function veclen(vec) !Compute the length of a 3D cartesian vector
  implicit none
  real*8 :: veclen,vec(3)
  veclen = dsqrt(vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3))
end function veclen
!************************************************************************

!************************************************************************
subroutine normvec(vec) !Normalise a 3D cartesian vector
  implicit none
  real*8 :: veclen,vec(3)
  vec = vec/veclen(vec)
end subroutine normvec
!************************************************************************



!************************************************************************
subroutine mc_eta_2_m1_m2(mc,eta,m1,m2)  !Convert chirp mass and eta to m1 and m2 - double precision
  implicit none
  real*8 :: mc,eta,m1,m2, dvar,mtot
  mtot = mc*eta**(-0.6d0)
  if(eta.le.0.25d0) then
     dvar = dsqrt(1.d0-4*eta)
     m1 = mtot/2.d0 * (1.0 + dvar);
     m2 = mtot/2.d0 * (1.0 - dvar);
  else                                 !Allow 0.25<eta<0.50
     dvar = dsqrt(4*eta-1.d0)
     m1 = mtot/2.d0 * (1.0 - dvar);
     m2 = mtot/2.d0 * (1.0 + dvar);
  end if
end subroutine mc_eta_2_m1_m2
!************************************************************************

!************************************************************************
subroutine mc_eta_2_m1_m2r(mcr,etar,m1r,m2r)  !Convert chirp mass and eta to m1 and m2 - single precision
  implicit none
  real*8 :: mc,eta,m1,m2
  real :: mcr,etar,m1r,m2r
  mc = dble(mcr)
  eta = dble(etar)
  call mc_eta_2_m1_m2(mc,eta,m1,m2)
  m1r = real(m1)
  m2r = real(m2)
end subroutine mc_eta_2_m1_m2r
!************************************************************************


!************************************************************************
subroutine m1_m2_2_mc_eta(m1,m2,mc,eta)
  implicit none
  real*8 :: m1,m2,mc,eta,mtot
  mtot = m1+m2
  eta = m1*m2/(mtot*mtot)
  mc = mtot*eta**0.6d0
end subroutine m1_m2_2_mc_eta
!************************************************************************


!************************************************************************
subroutine m1_m2_2_mc_etar(m1r,m2r,mcr,etar)
  implicit none
  real*8 :: m1,m2,mc,eta
  real :: m1r,m2r,mcr,etar
  m1 = dble(m1r)
  m2 = dble(m2r)
  call m1_m2_2_mc_eta(m1,m2,mc,eta)
  mcr = real(mc)
  etar = real(eta)
end subroutine m1_m2_2_mc_etar
!************************************************************************


!************************************************************************
subroutine ang2vec(l,b,vec)  !Convert longitude, latitude (rad) to a 3D normal vector
  !l in [0,2pi[; b in [-pi,pi]
  implicit none
  real*8 :: l,b,vec(3),cosb
  cosb = dcos(b)
  vec(1) = dcos(l)*cosb
  vec(2) = dsin(l)*cosb
  vec(3) = dsin(b)
end subroutine  ang2vec
!************************************************************************

!************************************************************************
subroutine vec2ang(vec,l,b)  !Convert a 3D normal vector to longitude, latitude (rad)
  !l in [0,2pi[; b in [-pi,pi]
  implicit none
  real*8 :: l,b,vec(3),vec1(3)
  vec1 = vec
  call normvec(vec1) !Make sure vec1 is normalised
  l = datan2(vec1(2),vec1(1))
  b = dasin(vec1(3))
end subroutine  vec2ang
!************************************************************************

!************************************************************************
function dotproduct(vec1,vec2) !Compute the dot product of two 3D cartesian vectors
  implicit none
  real*8 :: dotproduct,vec1(3),vec2(3)
  dotproduct = vec1(1)*vec2(1) + vec1(2)*vec2(2) + vec1(3)*vec2(3)
end function dotproduct
!************************************************************************

!************************************************************************
subroutine crossproduct(vec1,vec2,crpr) !Compute the cross (outer) product of two cartesian vectors
  implicit none
  real*8 :: vec1(3),vec2(3),crpr(3)!,veclen
  crpr(1) = vec1(2)*vec2(3) - vec1(3)*vec2(2)
  crpr(2) = vec1(3)*vec2(1) - vec1(1)*vec2(3)
  crpr(3) = vec1(1)*vec2(2) - vec1(2)*vec2(1)
  !write(6,'(3ES13.3)')veclen(vec1),veclen(vec2),veclen(crpr)
end subroutine crossproduct
!************************************************************************


!************************************************************************
function polangle(p,o)  !Compute the polarisation angle of a source with position normal vector p and orientation normal vector o, see Apostolatos et al. 1994, Eq.5
  implicit none
  real*8 :: polangle,p(3),o(3)
  real*8 :: z(3),denom,ocz(3),numer,dotproduct!,datan2
  
  z = (/0.d0,0.d0,1.d0/) !Vertical normal vector
  denom = dotproduct(o,z) - dotproduct(o,p)*dotproduct(z,p) !Denominator
  call crossproduct(o,z,ocz)
  numer = dotproduct(p,ocz) !Numerator
  
  !polangle = datan2(denom,numer)
  polangle = datan(denom/(numer+1.d-30))  !Take into account the degeneracy in psi
end function polangle
!************************************************************************

!************************************************************************
function posangle(p,o)  !Compute the position angle of a source with position normal vector p and orientation normal vector o
  implicit none
  real*8 :: posangle,p(3),o(3)
  real*8 :: x1(3),o1(3),z(3),z1(3),dotproduct
  
  !write(6,'(3ES13.3)')p
  !write(6,'(3ES13.3)')o
  
  call crossproduct(p,o,x1)
  call crossproduct(x1,p,o1) !o1: projection of o in the plane of the sky
  
  z = (/0.d0,0.d0,1.d0/) !Vertical normal vector
  call crossproduct(p,z,x1)
  call crossproduct(x1,p,z1) !z1: projection of z in the plane of the sky
  
  call normvec(o1)
  call normvec(z1)
  posangle = dacos(dotproduct(o1,z1))
  !write(6,'(3ES13.3)')p
  !write(6,'(3ES13.3)')o
  !write(6,'(3ES13.3)')o1
  !write(6,'(3ES13.3)')z1
  !write(6,'(3ES13.3)')dotproduct(o1,z1),dacos(dotproduct(o1,z1))
end function posangle
!************************************************************************


!>
!! Compute the inclination and polarisation angle for a source with position (pl,pb) and orientation (ol,ob)
!! All variables are angles (no cos, sin)
!<
!************************************************************************
subroutine compute_incli_polang(pl,pb,ol,ob, i,psi) 
  use constants
  implicit none
  !pl,ol in [0,2pi[;  pb,ob in ([-pi/2,pi/2]) now [0,pi], conf John & Christian
  real*8 :: pl,pb,ol,ob
  real*8 :: p(3),o(3),i,dotproduct,psi,polangle,drevpi
  
  call ang2vec(pl,pb,p)       !Position normal vector
  call ang2vec(ol,ob,o)       !Orientation normal vector
  
  !i = pi2 - dacos(dotproduct(p,o))  !Compute inclination angle: <0: points towards us, >0 points away from us
  !psi = polangle(p,o)         !Compute polarisation angle [-pi/2,pi/2]
  
  i = dacos(dotproduct(p,o))   !Compute inclination angle: 0: points exactly away from us, 180 points exactly towards us, 90: in the plane of the sky
  psi = drevpi(polangle(p,o))  !Compute polarisation angle [0,pi]
  
end subroutine compute_incli_polang
!************************************************************************

!>
!! Compute the inclination and polarisation angle for a source with position (pl,pb) and orientation (ol,ob)
!! All variables are angles (no cos, sin)
!! Single-precision wrapper for compute_incli_polang()
!<
!************************************************************************
subroutine compute_incli_polangr(plr,pbr,olr,obr, ir,psir)
  implicit none
  real*8 :: pl,pb,ol,ob,i,psi
  real :: plr,pbr,olr,obr,ir,psir
  
  pl = dble(plr)
  pb = dble(pbr)
  ol = dble(olr)
  ob = dble(obr)
  call compute_incli_polang(pl,pb,ol,ob, i,psi)
  ir = real(i)
  psir = real(psi)
end subroutine compute_incli_polangr
!************************************************************************

!************************************************************************
subroutine compute_incli_posang(pl,pb,ol,ob, i,psi) !Compute the inclination and position angle for a source with position (pl,pb) and orientation (ol,ob)
  use constants
  implicit none
  !pl,ol in [0,2pi[;  pb,ob in [-pi,pi]
  real*8 :: pl,pb,ol,ob
  real*8 :: p(3),o(3),i,dotproduct,psi,posangle!,drevpi
  
  call ang2vec(pl,pb,p)       !Position normal vector
  call ang2vec(ol,ob,o)       !Orientation normal vector
  !i = pi2 - dacos(dotproduct(p,o))  !Compute inclination angle: <0: points towards us, >0 points away from us
  i = dacos(dotproduct(p,o))  !Compute inclination angle: 0: points exactly away from us, 180 points exactly towards us, 90: in the plane of the sky
  psi = posangle(p,o)         !Compute position angle
  !psi = drevpi(posangle(p,o))  !Compute position angle
  
end subroutine compute_incli_posang
!************************************************************************



!************************************************************************
subroutine detectorvector(d1,d2,jd)  !Determine the sky position at which the vector that connects two detectors points
  implicit none
  integer :: d1,d2
  real*8 :: jd,detcoords(3,2),vec1(3),vec2(3),dvec(3),l,b
  jd = 0 !get rid of warnings
  detcoords(1,:) = (/-119.41,46.45/)  !H1; l,b
  detcoords(2,:) = (/-90.77,30.56/)   !L1
  detcoords(3,:) = (/10.50,43.63/)    !V
  
  call ang2vec(detcoords(d1,1),detcoords(d1,2),vec1)
  call ang2vec(detcoords(d2,1),detcoords(d2,2),vec2)
  
  dvec = vec2 - vec1
  
  call vec2ang(dvec,l,b)  !Searched point is in zenith/nadir for an observer on this location on the globe
  
end subroutine detectorvector
!************************************************************************


!************************************************************************
subroutine swapint(i1,i2)                        !Swap two integers
  implicit none
  integer :: i,i1,i2
  i = i1
  i1 = i2
  i2 = i
end subroutine swapint
!************************************************************************


!***************************************************************************************************
subroutine determine_nbin_1d(npoints,nbin)
  implicit none
  integer :: npoints,nbin
  if(npoints.le.100) then
     nbin = floor(2*sqrt(real(npoints)))
  else
     nbin = floor(10*log10(real(npoints)))
  end if
  nbin = max(nbin,5)
end subroutine determine_nbin_1d
!***************************************************************************************************
