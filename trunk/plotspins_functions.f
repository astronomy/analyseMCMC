!Functions for plotspins.f


!***************************************************************************************************
module plotspins_settings
  implicit none
  save
  integer, parameter :: nchs=10,npar1=15
  integer :: plvars(npar1),nplvar,nbin1d,nbin2dx,nbin2dy,npdf2d,pdf2dpairs(250,2),panels(2)
  integer :: thin,nburn(nchs),file,colour,quality,reverseread,update,mergechains,wrapdata,changevar,maxchlen
  integer :: prprogress,prruninfo,prinitial,prstat,prcorr,prival,prconv,savestats,savepdf       
  integer :: plot,combinechainplots,pllogl,plchain,plparl,pljump,rdsigacc,plsigacc,plpdf1d,plpdf2d,placorr,plotsky,plmovie       
  integer :: chainsymbol,chainpli,pltrue,plstart,plmedian,plrange,plburn,pllmax,prvalues,smooth,fillpdf,normpdf1d,normpdf2d
  integer :: scloglpl,scchainspl,bmpxsz,bmpysz
  integer :: nmovframes,moviescheme,whitebg,unsharp
  real :: ival0,nburnfrac,autoburnin
  real :: scrsz,scrrat,pssz,psrat,scfac
end module plotspins_settings
!***************************************************************************************************

!***************************************************************************************************
module constants
  implicit none
  save
  real*8 :: pi,tpi,pi2,r2d,d2r,r2h,h2r,c3rd
end module constants
!***************************************************************************************************

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
end subroutine setconstants
!***************************************************************************************************


!***************************************************************************************************
subroutine read_inputfile
  use plotspins_settings
  implicit none
  integer :: i,u,io,io1
  character :: bla,filename*99
  filename = 'plotspins.dat'
  
  u = 15
  open(unit=u,form='formatted',status='old',action='read',file=trim(filename),iostat=io)
  if(io.ne.0) then
     write(*,'(A)')'  Error opening input file '//trim(filename)//', aborting...'
     stop
  end if
  
  io = 0
  io1 = 0
  
  read(u,*,iostat=io)bla
  

  read(u,*,iostat=io)bla
  read(u,*,iostat=io)thin
  read(u,*,iostat=io)nburn(1)
  do i=2,nchs
     nburn(i) = nburn(1)
  end do
  read(u,*,iostat=io)nburnfrac
  read(u,*,iostat=io)autoburnin
  read(u,*,iostat=io)maxchlen
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
  read(u,*,iostat=io)rdsigacc
  read(u,*,iostat=io)plsigacc
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
  read(u,*,iostat=io)ival0
  read(u,*,iostat=io)nmovframes
  read(u,*,iostat=io)moviescheme

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
     write(*,'(/,A,I2,A,I3,A,/)')'  Error reading input file '//trim(filename)//', aborting...'
     stop
  end if
  
end subroutine read_inputfile
!***************************************************************************************************


!***************************************************************************************************
subroutine write_inputfile
  use plotspins_settings
  implicit none
  integer :: u,i
  
  u = 14
  open(unit=u,form='formatted',status='replace',file='plotspins.used')
  
11 format(I10,1x,A19,5x,A)
12 format(2I5,1x,A19,5x,A)
21 format(F10.5,1x,A19,5x,A)
31 format(ES10.2,1x,A19,5x,A)
  
  write(u,'(A,/)')' Input file for plotspins.f'
  
  
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
  write(u,11)prprogress, 'prprogress',   'Print general messages about the progress of the program: 0-no, 1-some, 2-more'
  write(u,11)prruninfo, 'prruninfo',   'Print run info at read (# iterations, seed, # detectors, SNRs, data length, etc.): 0-no, 1-only for one file (eg. if all files similar), 2-for all files'
  write(u,11)prinitial, 'prinitial',   'Print true values, starting values and their difference'
  write(u,11)prstat, 'prstat',   'Print statistics: 0-no, 1-yes'
  write(u,11)prcorr, 'prcorr',   'Print correlations: 0-no, 1-yes'
  write(u,11)prival, 'prival',   'Print interval info: 0-no, 1-yes'
  write(u,11)prconv, 'prconv',   'Print convergence information for multiple chains to screen and chains plot: 0-no, 1-one summary line, 2-medians, stdevs, etc. too.'
  write(u,11)savestats, 'savestats',   'Save statistics (statistics, correlations, intervals) to file: 0-no, 1-yes, 2-yes + copy in PS'
  write(u,11)savepdf, 'savepdf',   'Save the binned data for 1d and/or 2d pdfs (depending on plpdf1d and plpdf2d).  This causes all 12 parameters + m1,m2 to be saved and plotted(!), which is slighty annoying'
  
  
  write(u,'(/,A)')' Select which plots to make:'
  write(u,11)plot, 'plot',   '0: plot nothing at all, 1: plot the items selected below'
  write(u,11)combinechainplots, 'combinechainplots',   'Combine logL, chain, sigma and acc plots into one multipage file'
  write(u,11)pllogl, 'pllogl',   'Plot log L chains: 0-no, 1-yes'
  write(u,11)plchain, 'plchain',   'Plot parameter chains: 0-no, 1-yes'
  write(u,11)plparl, 'plparl',   'Plot L vs. parameter value: 0-no, 1-yes'
  write(u,11)pljump, 'pljump',   'Plot actual jump sizes: 0-no, 1-yes: lin, 2-yes: log'
  write(u,11)rdsigacc, 'rdsigacc',   'Read sigma and acceptance rate: 0-no, 1-yes   (0-Dont read these data, save 40% read-in time).  0 can give problems with large scale, or high-temperature chains'
  write(u,11)plsigacc, 'plsigacc',   'Plot sigma and acceptance rate: 0-no, 1-yes (lin sig), 2-yes (log sig)  (If >0, this sets rdsigacc to 1)'
  write(u,11)plpdf1d, 'plpdf1d',   'Plot 1d posterior distributions: 0-no, 1-yes. If plot=0 and savepdf=1, this determines whether to write the pdfs to file or not.'
  write(u,11)plpdf2d, 'plpdf2d',   'Plot 2d posterior distributions: 0-no, 1-yes: gray + contours, 2:gray only, 3: contours only. If plot=0 and savepdf=1, this determines whether to write the pdfs to file (>0) or not (=0).'
  write(u,11)placorr, 'placorr',   'Plot autocorrelations: 0-no, >0-yes: plot placorr steps'
  write(u,11)plotsky, 'plotsky',   'Plot 2d pdf with stars, implies plpdf2d=1'
  write(u,11)plmovie, 'plmovie',   'Plot movie frames'
  
  
  write(u,'(/,A)')' Detailed plot settings:'
  write(u,11)chainsymbol, 'chainsymbol',   'Plot symbol for the chains: 0-plot lines, !=0: plot symbols: eg: 1: dot (default), 2: plus, etc.  -4: filled diamond, 16,17: filled square,circle 20: small open circle; -10/-11: use a selection of open/filled symbols'
  write(u,11)chainpli, 'chainpli',   'Plot every chainpli-th point in chains, logL, jump plots:  chainpli=0: autodetermine, chainpli>0: use this chainpli.  All states in between *are* used for statistics, pdf generation, etc.'
  write(u,11)scloglpl, 'scloglpl',   'Scale logL plot ranges: 0: take everything into account, including burnin and starting values;  1: take only post-burnin and true values into account'
  write(u,11)scchainspl, 'scchainspl',   'Scale chains plot ranges: 0: take everything into account, including burnin;  1: take only post-burnin and true values into account'
  write(u,11)pltrue, 'pltrue',   'Plot true values in the chains and pdfs'
  write(u,11)plstart, 'plstart',   'Plot starting values in the chains and pdfs'
  write(u,11)plmedian, 'plmedian',   'Plot median values in the pdfs: 1-1D PDFs, 2-2D PDFs, 3-both'
  write(u,11)plrange, 'plrange',   'Plot the probability range in the pdfs: 1-1D PDFs, 2-2D PDFs, 3-both'
  write(u,11)plburn, 'plburn',   'Plot the burnin in logL, the chains, etc.'
  write(u,11)pllmax, 'pllmax',   'Plot the position of the max logL, in the chains and pdfs'
  write(u,11)prvalues, 'prvalues',   'Print values (true, median, range) in pdfs'
  write(u,11)smooth, 'smooth',   'Smooth the pdfs: 0 - no, >1: smooth over smooth bins (use ~10 (3-15)?).   This is 1D only for now, and can introduce artefacts on narrow peaks!'
  write(u,11)fillpdf, 'fillpdf',   'Fillstyle for the pdfs (pgsfs): 1-solid, 2-outline, 3-hatched, 4-cross-hatched'
  write(u,11)normpdf1d, 'normpdf1d',   'Normalise 1D pdfs:  0-no,  1-normalise surface area (default, a must for different bin sizes),  2-normalise to height,  3-normalise to sqrt(height), nice to compare par.temp. chains'
  write(u,11)normpdf2d, 'normpdf2d',   "'Normalise' 2D pdfs; greyscale value depends on bin height:  0-linearly,  1-logarithmically,  2-sqrt,  3-weigted with likelihood value"
  write(u,21)ival0, 'ival0',   'Standard probability interval, e.g. 0.90, 0.95'
  write(u,11)nmovframes, 'nmovframes',   'Number of frames for the movie'
  write(u,11)moviescheme, 'moviescheme',   'Moviescheme (1-3): determines what panels to show in a movie frame; see source code'
  
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
  write(u,11)nplvar, 'nplvar',   'Number of plot variables for 1D PDFs (and chain, jump plots, max 15).  Put the variables in the line below:'
  do i=1,nplvar
     write(u,'(I3,$)')plvars(i)
  end do
  write(u,*)''
  !write(u,'(5x,A,5x,A)')'plvars','The actual plot variables (1-nplvar)'
  write(u,12)panels(1:2), 'panels',   'Number of for 1D plots in x,y direction: 0: autodetermine'
  write(u,11)nbin1d, 'nbin1d',   'Number of bins for 1D PDFs'
  write(u,11)nbin2dx, 'nbin2dx',   'Number of bins in x-direction for 2D PDFs'
  write(u,11)nbin2dy, 'nbin2dy',   'Number of bins in y-direction for 2D PDFs'
  write(u,11)npdf2d, 'npdf2d',     'Number of 2D-PDF plots to make'
  do i=1,npdf2d
     write(u,12)pdf2dpairs(i,1:2), 'pdf2dpairs', 'Pairs of parameters to plot a 2D PDF for'
  end do
  close(u)
end subroutine write_inputfile
!***************************************************************************************************




!***************************************************************************************************
subroutine set_plotsettings  !Set plot settings to 'default' values
  use plotspins_settings
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
  combinechainplots = 0  !Combine logL, chain, sigma and acc plots into one multipage file
  autoburnin = 1.   !Determine burnin automatically as the first iteration where log(L_chain) > max(log(L_allchains)) - autoburnin
  scloglpl = 1      !Scale logL plot ranges: 0: take everything into account, including burnin and starting values;  1: take only post-burnin and true values into account
  scchainspl = 1    !Scale chains plot ranges: 0: take everything into account, including burnin;  1: take only post-burnin and true values into account
  pllogl = 1        !Plot log L chains: 0-no, 1-yes
  plchain = 1       !Plot parameter chains: 0-no, 1-yes
  plparl = 1        !Plot L vs. parameter value: 0-no, 1-yes
  pljump = 1        !Plot actual jump sizes
  rdsigacc = 1      !Read sigma and acceptance rate: 0-no, 1-yes   (0-Don't read these data, save 40% read-in time).  0 can give problems with large scale, or high-temperature chains
  plsigacc = 0      !Plot sigma and acceptance rate: 0-no, 1-yes   (Sets rdsigacc to 1)
  plpdf1d = 1       !Plot 1d posterior distributions: 0-no, 1-yes. If plot=0 and savepdf=1, this determines whether to write the pdfs to file or not.
  plpdf2d = 2       !Plot 2d posterior distributions: 0-no, 1-yes: gray + contours, 2:gray only, 3: contours only. If plot=0 and savepdf=1, this determines whether to write the pdfs to file (>0) or not (=0).
  placorr = 0e4     !Plot autocorrelations: 0-no, >0-yes: plot placorr steps
  plotsky = 0       !Plot 2d pdf with stars, implies plpdf2d=1
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
  ival0 = 0.90      !Standard probability interval, e.g. 0.90, 0.95
  nmovframes = 1    !Number of frames for the movie
  moviescheme = 3   !Movie scheme: determines what panels to show in a movie frame
  
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
subroutine pginitl(colour,file,whitebg)  !Initialise pgplot
  implicit none
  integer :: colour,file,i,whitebg
  if(whitebg.ge.1) then
     call pgscr(0,1.,1.,1.) !Background colour always white (also on screen, bitmap)
     call pgscr(1,0.,0.,0.) !Default foreground colour always black
     if(file.le.1) then !png: create white background
        call pgsvp(-100.,100.,-100.,100.)
        call pgswin(0.,1.,0.,1.)
        call pgsci(0)
        call pgrect(-1.,2.,-1.,2.)
        call pgsvp(0.08,0.95,0.06,0.87) !Default viewport size (?)
        call pgsci(1)
     end if
  end if
  if(colour.eq.0) then
     do i=0,99
        call pgscr(i,0.,0.,0.)
     end do
     call pgscr(0,1.,1.,1.)     !White
     call pgscr(14,0.3,0.3,0.3) !Dark grey
     call pgscr(15,0.6,0.6,0.6) !Light grey
  else
     call pgscr(2,1.,0.1,0.1)  !Make default red lighter
     call pgscr(3,0.,0.5,0.)   !Make default green darker
     call pgscr(4,0.,0.,0.8)   !Make default blue darker
     call pgscr(7,0.9,0.9,0.)  !Make default yellow darker
     call pgscr(10,0.5,0.3,0.) !10: brown
     call pgscr(11,0.6,0.,0.)  !11: dark red
  end if
end subroutine pginitl
!************************************************************************************************************************************



!************************************************************************************************************************************
subroutine bindata(n,x,norm,nbin,xmin1,xmax1,xbin,ybin)  !Count the number of points in each bin
  ! x - input: data, n points
  ! norm - input: normalise (1) or not (0)
  ! nbin - input: number of bins
  ! xmin, xmax - in/output: set xmin=xmax to auto-determine
  ! xbin, ybin - output: binned data (x, y).  The x values are the left side of the bin!

  implicit none
  integer :: i,k,n,nbin,norm
  real :: x(n),xbin(nbin+1),ybin(nbin+1),xmin,xmax,dx,xmin1,xmax1

  xmin = xmin1
  xmax = xmax1

  if(abs(xmin-xmax)/(xmax+1.e-30).lt.1.e-20) then !Autodetermine
     xmin = minval(x(1:n))
     xmax = maxval(x(1:n))
  end if
  dx = abs(xmax - xmin)/real(nbin)

  do k=1,nbin+1
     !xbin(k) = xmin + (real(k)-0.5)*dx  !x is the centre< of the bin
     xbin(k) = xmin + (k-1)*dx          !x is the left of the bin
  end do
  !ybintot=0.
  ybin = 0.
  do i=1,n
     do k=1,nbin
        if(x(i).ge.xbin(k)) then
           if(x(i).lt.xbin(k+1)) then
              ybin(k) = ybin(k) + 1.
              exit !If point i fits in this bin, don't try the others
           end if
        end if
     end do !k (bin)
     !ybintot = ybintot + ybin(k)
  end do
  !if(norm.eq.1) ybin = ybin/(ybintot+1.e-30)
  if(norm.eq.1) ybin = ybin/(sum(ybin)+1.e-30)

  if(abs(xmin1-xmax1)/(xmax1+1.e-30).lt.1.e-20) then
     xmin1 = xmin
     xmax1 = xmax
  end if

end subroutine bindata
!************************************************************************************************************************************


!************************************************************************************************************************************
subroutine bindataa(n,x,y,norm,nbin,xmin1,xmax1,xbin,ybin)  !Measure the amount of likelihood in each bin
  ! x - input: data, n points
  ! y - input: "weight" (likelihood), n points
  ! norm - input: normalise (1) or not (0)
  ! nbin - input: number of bins
  ! xmin, xmax - in/output: set xmin=xmax to auto-determine
  ! xbin, ybin - output: binned data (x, y).  The x values are the left side of the bin!

  implicit none
  integer :: i,k,n,nbin,norm
  real :: x(n),y(n),xbin(nbin+1),ybin(nbin+1),xmin,xmax,dx,ybintot,xmin1,xmax1,ymin
  
  xmin = xmin1
  xmax = xmax1
  ymin = minval(y)
  !print*,n,nbin,xmin1,xmax1
  !print*,minval(y),maxval(y)

  if(abs(xmin-xmax)/(xmax+1.e-30).lt.1.e-20) then !Autodetermine
     xmin = minval(x(1:n))
     xmax = maxval(x(1:n))
  end if
  dx = abs(xmax - xmin)/real(nbin)

  do k=1,nbin+1
     !        xbin(k) = xmin + (real(k)-0.5)*dx  !x is the centre< of the bin
     xbin(k) = xmin + (k-1)*dx          !x is the left of the bin
  end do
  ybintot=0.
  do k=1,nbin
     ybin(k) = 0.
     do i=1,n
        !if(x(i).ge.xbin(k).and.x(i).lt.xbin(k+1)) ybin(k) = ybin(k) + 1.
        if(x(i).ge.xbin(k).and.x(i).lt.xbin(k+1)) ybin(k) = ybin(k) + exp(y(i) - ymin)
     end do
     ybintot = ybintot + ybin(k)
  end do
  if(norm.eq.1) ybin = ybin/(ybintot+1.e-30)

  if(abs(xmin1-xmax1)/(xmax1+1.e-30).lt.1.e-20) then
     xmin1 = xmin
     xmax1 = xmax
  end if

end subroutine bindataa
!************************************************************************************************************************************


!************************************************************************************************************************************
subroutine bindata2d(n,x,y,norm,nxbin,nybin,xmin1,xmax1,ymin1,ymax1,z,tr)  !Count the number of points in each bin
  !x - input: data, n points
  !norm - input: normalise (1) or not (0)
  !nbin - input: number of bins
  !xmin, xmax - in/output: set xmin=xmax to auto-determine
  !xbin, ybin - output: binned data (x, y).  The x values are the left side of the bin!
  
  implicit none
  integer :: i,n,bx,by,nxbin,nybin,norm
  real :: x(n),y(n),xbin(nxbin+1),ybin(nybin+1),z(nxbin+1,nybin+1)
  real :: xmin,xmax,ymin,ymax,dx,dy,xmin1,xmax1,ymin1,ymax1,tr(6)
  
  !write(*,'(A4,5I8)')'n:',norm,nxbin,nybin
  !write(*,'(A4,2F8.3)')'x:',xmin1,xmax1
  !write(*,'(A4,2F8.3)')'y:',ymin1,ymax1
  
  xmin = xmin1
  xmax = xmax1
  ymin = ymin1
  ymax = ymax1
  
  if(abs(xmin-xmax)/(xmax+1.e-30).lt.1.e-20) then !Autodetermine
     xmin = minval(x(1:n))
     xmax = maxval(x(1:n))
  end if
  dx = abs(xmax - xmin)/real(nxbin)
  if(abs(ymin-ymax)/(ymax+1.e-30).lt.1.e-20) then !Autodetermine
     ymin = minval(y(1:n))
     ymax = maxval(y(1:n))
  end if
  dy = abs(ymax - ymin)/real(nybin)
  do bx=1,nxbin+1
     !xbin(bx) = xmin + (real(bx)-0.5)*dx  !x is the centre of the bin
     xbin(bx) = xmin + (bx-1)*dx          !x is the left of the bin
  end do
  do by=1,nybin+1
     !ybin(by) = ymin + (real(by)-0.5)*dy  !y is the centre of the bin
     ybin(by) = ymin + (by-1)*dy          !y is the left of the bin
  end do
  
  !write(*,'(50F5.2)'),x(1:50)
  !write(*,'(50F5.2)'),y(1:50)
  !write(*,'(20F8.5)'),xbin
  !write(*,'(20F8.5)'),ybin
  
  z = 0.
  !ztot = 0.
  do i=1,n
     bxl: do bx=1,nxbin
        do by=1,nybin
           if(x(i).ge.xbin(bx)) then
              if(x(i).lt.xbin(bx+1)) then
                 if(y(i).ge.ybin(by)) then
                    if(y(i).lt.ybin(by+1)) then
                       z(bx,by) = z(bx,by) + 1.
                       exit bxl !exit bx loop; if point i fits this bin, don't try other bins
                    end if
                 end if
              end if
           end if
           
        end do !by
     end do bxl !bx
  end do !i
  !if(norm.eq.1) z = z/(ztot+1.e-30)
  if(norm.eq.1) z = z/maxval(z+1.e-30)
  
  if(abs(xmin1-xmax1)/(xmax1+1.e-30).lt.1.e-20) then
     xmin1 = xmin
     xmax1 = xmax
  end if
  if(abs(ymin1-ymax1)/(ymax1+1.e-30).lt.1.e-20) then
     ymin1 = ymin
     ymax1 = ymax
  end if
  
  !Determine transformation elements for pgplot (pggray, pgcont, pgimag)
  tr(1) = xmin - dx/2.
  tr(2) = dx
  tr(3) = 0.
  tr(4) = ymin - dy/2.
  tr(5) = 0.
  tr(6) = dy
  
end subroutine bindata2d
!************************************************************************************************************************************


!************************************************************************************************************************************
subroutine bindata2da(n,x,y,z,norm,nxbin,nybin,xmin1,xmax1,ymin1,ymax1,zz,tr)  !Measure the amount of likelihood in each bin
  !x,y - input: data, n points
  !z - input: amount for each point (x,y)
  !norm - input: normalise (1) or not (0)
  !nxbin,nybin - input: number of bins in each dimension
  !xmin1,xmax1 - in/output: ranges in x dimension, set xmin=xmax as input to auto-determine
  !ymin1,ymax1 - in/output: ranges in y dimension, set ymin=ymax as input to auto-determine
  !zz - output: binned data zz(x,y).  The x,y values are the left side of the bin(?)
  !tr - output: transformation elements for pgplot (pggray, pgcont)
  
  implicit none
  integer :: i,n,bx,by,nxbin,nybin,norm
  real :: x(n),y(n),z(n),xbin(nxbin+1),ybin(nybin+1),zz(nxbin+1,nybin+1),zztot,xmin,xmax,ymin,ymax,dx,dy,xmin1,xmax1,ymin1,ymax1
  real :: tr(6),zmin
  
  !write(*,'(A4,5I8)')'n:',norm,nxbin,nybin
  !write(*,'(A4,2F8.3)')'x:',xmin1,xmax1
  !write(*,'(A4,2F8.3)')'y:',ymin1,ymax1
  
  xmin = xmin1
  xmax = xmax1
  ymin = ymin1
  ymax = ymax1
  zmin = minval(z)
  
  if(abs(xmin-xmax)/(xmax+1.e-30).lt.1.e-20) then !Autodetermine
     xmin = minval(x(1:n))
     xmax = maxval(x(1:n))
  end if
  dx = abs(xmax - xmin)/real(nxbin)
  if(abs(ymin-ymax)/(ymax+1.e-30).lt.1.e-20) then !Autodetermine
     ymin = minval(y(1:n))
     ymax = maxval(y(1:n))
  end if
  dy = abs(ymax - ymin)/real(nybin)
  do bx=1,nxbin+1
     !xbin(bx) = xmin + (real(bx)-0.5)*dx  !x is the centre of the bin
     xbin(bx) = xmin + (bx-1)*dx          !x is the left of the bin
  end do
  do by=1,nybin+1
     !ybin(by) = ymin + (real(by)-0.5)*dy  !y is the centre of the bin
     ybin(by) = ymin + (by-1)*dy          !y is the left of the bin
  end do
  
  !write(*,'(50F5.2)'),x(1:50)
  !write(*,'(50F5.2)'),y(1:50)
  !write(*,'(20F8.5)'),xbin
  !write(*,'(20F8.5)'),ybin
  
  zz = 0.
  zztot = 0.
  !print*,xmin,xmax
  !print*,ymin,ymax
  do bx=1,nxbin
     !print*,bx,xbin(bx),xbin(bx+1)
     do by=1,nybin
        zz(bx,by) = 0.
        do i=1,n
           !if(x(i).ge.xbin(bx).and.x(i).lt.xbin(bx+1) .and. y(i).ge.ybin(by).and.y(i).lt.ybin(by+1)) zz(bx,by) = zz(bx,by) + 1.
           if(x(i).ge.xbin(bx).and.x(i).lt.xbin(bx+1) .and. y(i).ge.ybin(by).and.y(i).lt.ybin(by+1)) zz(bx,by) = zz(bx,by) + exp(z(i) - zmin)
           !write(*,'(2I4,8F10.5)')bx,by,x(i),xbin(bx),xbin(bx+1),y(i),ybin(by),ybin(by+1),zz(bx,by),z(i)
        end do
        zztot = zztot + zz(bx,by) 
        !write(*,'(2I4,5x,4F6.3,5x,10I8)')bx,by,xbin(bx),xbin(bx+1),ybin(by),ybin(by+1),nint(zz(bx,by))
     end do
     !write(*,'(I4,5x,2F6.3,5x,10I8)')bx,xbin(bx),xbin(bx+1),nint(zz(bx,1:nybin))
     end do
  !if(norm.eq.1) z = z/(zztot+1.e-30)
  if(norm.eq.1) z = z/maxval(z+1.e-30)
  
  if(abs(xmin1-xmax1)/(xmax1+1.e-30).lt.1.e-20) then
     xmin1 = xmin
     xmax1 = xmax
  end if
  if(abs(ymin1-ymax1)/(ymax1+1.e-30).lt.1.e-20) then
     ymin1 = ymin
     ymax1 = ymax
  end if
  
  !Determine transformation elements for pgplot (pggray, pgcont)
  tr(1) = xmin - dx/2.
  tr(2) = dx
  tr(3) = 0.
  tr(4) = ymin - dy/2.
  tr(5) = 0.
  tr(6) = dy
  
end subroutine bindata2da
!************************************************************************************************************************************


!************************************************************************************************************************************
subroutine verthist(n,x,y)  !x is the left of the bin!
  implicit none
  integer :: j,n
  real :: x(n+1),y(n+1)

  !      call pgline(2,(/0.,x(1)/),(/y(1),y(1)/))
  !      do j=1,n
  !        call pgline(2,(/x(j),x(j)/),y(j:j+1))
  !        call pgline(2,x(j:j+1),(/y(j+1),y(j+1)/))
  !      end do

  x(n+1) = x(n) + (x(n)-x(n-1))
  y(n+1) = 0.
  call pgline(2,(/x(1),x(1)/),(/0.,y(1)/))
  do j=1,n
     call pgline(2,x(j:j+1),(/y(j),y(j)/))
     call pgline(2,(/x(j+1),x(j+1)/),y(j:j+1))
  end do

end subroutine verthist
!************************************************************************************************************************************

!************************************************************************************************************************************
subroutine horzhist(n,x,y)
  implicit none
  integer :: j,n
  real :: x(n),y(n)

  call pgline(2,(/x(1),x(1)/),(/0.,y(1)/))
  do j=1,n-2
     call pgline(2,x(j:j+1),(/y(j),y(j)/))
     call pgline(2,(/x(j+1),x(j+1)/),y(j:j+1))
  end do
  call pgline(2,x(n-1:n),(/y(n-1),y(n-1)/))
end subroutine horzhist
!************************************************************************************************************************************



!************************************************************************************************************************************
function ra(lon, GPSsec)
  ! Derives right ascension (in radians!) from longitude (radians) and GPS time. 
  ! Declination == latitude for equatorial coordinates.                        
  
  
  ! Derive the `Greenwich Mean Sidereal Time' (in radians!) 
  ! from GPS time (in seconds).                              
  ! (see K.R.Lang(1999), p.80sqq.)                           
  implicit none
  real*8 :: ra,lon,gmst,seconds,days,centuries,secCurrentDay
  real*8 :: gps0,leapseconds,GPSsec,tpi
  tpi = 8*datan(1.d0)

  gps0 = 630720013.d0 !GPS time at 1/1/2000 at midnight
  leapseconds = 32.d0 !At Jan 1st 2000
  if(GPSsec.gt.(gps0 + 189388800.d0)) leapseconds = leapseconds + 1.d0 !One more leapsecond after 1/1/2006
  if(GPSsec.lt.630720013.d0) write(*,'(A)')'WARNING: GMSTs before 1.1.2000 are inaccurate!'
  !Time since 1/1/2000 midnight
  seconds       = (GPSsec - gps0) + (leapseconds - 32.d0)
  days          = floor(seconds/86400.d0) - 0.5d0
  secCurrentDay = mod(seconds, 86400.d0)
  centuries     = days/36525.d0
  gmst = 24110.54841d0 + (centuries*(8640184.812866d0 + centuries*(0.093104d0 + centuries*6.2d-6)))
  gmst = gmst + secCurrentDay * 1.002737909350795d0   !UTC day is 1.002 * MST day
  gmst = mod(gmst/86400.d0,1.d0)
  gmst = gmst * tpi

  ra = mod(lon + gmst + 10*tpi,tpi)
end function ra
!************************************************************************************************************************************



!************************************************************************************************************************************
subroutine dindexx(n,arr,indx)
  integer :: n,indx(n),m,nstack
  real*8 :: arr(n),a
  parameter (m=7,nstack=50)
  integer :: i,indxt,ir,itemp,j,jstack,k,l,istack(nstack)
  do j=1,n
     indx(j)=j
  end do
  jstack=0
  l=1
  ir=n
1 if(ir-l.lt.m)then
     do j=l+1,ir
        indxt=indx(j)
        a=arr(indxt)
        do i=j-1,l,-1
           if(arr(indx(i)).le.a)goto 2
           indx(i+1)=indx(i)
        end do
        i=l-1
2       indx(i+1)=indxt
     end do
     if(jstack.eq.0)return
     ir=istack(jstack)
     l=istack(jstack-1)
     jstack=jstack-2
  else
     k=(l+ir)/2
     itemp=indx(k)
     indx(k)=indx(l+1)
     indx(l+1)=itemp
     if(arr(indx(l)).gt.arr(indx(ir)))then
        itemp=indx(l)
        indx(l)=indx(ir)
        indx(ir)=itemp
     end if
     if(arr(indx(l+1)).gt.arr(indx(ir)))then
        itemp=indx(l+1)
        indx(l+1)=indx(ir)
        indx(ir)=itemp
     end if
     if(arr(indx(l)).gt.arr(indx(l+1)))then
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
     if(arr(indx(i)).lt.a)goto 3
4    continue
     j=j-1
     if(arr(indx(j)).gt.a)goto 4
     if(j.lt.i)goto 5
     itemp=indx(i)
     indx(i)=indx(j)
     indx(j)=itemp
     goto 3
5    indx(l+1)=indx(j)
     indx(j)=indxt
     jstack=jstack+2
     !if(jstack.gt.nstack)pause 'nstack too small in indexx'
     if(jstack.gt.nstack) write(*,'(A)')' nstack too small in dindexx'
     if(ir-i+1.ge.j-l)then
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
1 if(ir-l.lt.m)then
     do j=l+1,ir
        indxt=indx(j)
        a=arr(indxt)
        do i=j-1,l,-1
           if(arr(indx(i)).le.a)goto 2
           indx(i+1)=indx(i)
        end do
        i=l-1
2       indx(i+1)=indxt
     end do
     if(jstack.eq.0)return
     ir=istack(jstack)
     l=istack(jstack-1)
     jstack=jstack-2
  else
     k=(l+ir)/2
     itemp=indx(k)
     indx(k)=indx(l+1)
     indx(l+1)=itemp
     if(arr(indx(l)).gt.arr(indx(ir)))then
        itemp=indx(l)
        indx(l)=indx(ir)
        indx(ir)=itemp
     end if
     if(arr(indx(l+1)).gt.arr(indx(ir)))then
        itemp=indx(l+1)
        indx(l+1)=indx(ir)
        indx(ir)=itemp
     end if
     if(arr(indx(l)).gt.arr(indx(l+1)))then
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
     if(arr(indx(i)).lt.a)goto 3
4    continue
     j=j-1
     if(arr(indx(j)).gt.a)goto 4
     if(j.lt.i)goto 5
     itemp=indx(i)
     indx(i)=indx(j)
     indx(j)=itemp
     goto 3
5    indx(l+1)=indx(j)
     indx(j)=indxt
     jstack=jstack+2
     !if(jstack.gt.nstack)pause 'nstack too small in indexx'
     if(jstack.gt.nstack) write(*,'(A)')' nstack too small in rindexx'
     if(ir-i+1.ge.j-l)then
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
  !     USES lubksb,ludcmp
  integer :: imj,ipj,j,k,kk,mm,indx(mmax+1)
  real :: d,fac,sum,a(mmax+1,mmax+1),b(mmax+1)
  !if(np.lt.nl+nr+1.or.nl.lt.0.or.nr.lt.0.or.ld.gt.m.or.m.gt.mmax.or.nl+nr.lt.m) pause 'bad args in savgol'
  if(np.lt.nl+nr+1.or.nl.lt.0.or.nr.lt.0.or.ld.gt.m.or.m.gt.mmax.or.nl+nr.lt.m) write(*,'(A)')' Bad args in savgol'
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
     if (ii.ne.0)then
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
    if(aamax.eq.0.) write(*,'(A)')' Singular matrix in ludcmp'
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
    if (j.ne.imax)then
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
    if(j.ne.n)then
       dum=1./a(j,j)
       do i=j+1,n
          a(i,j)=a(i,j)*dum
       end do
    end if
 end do
 return
end subroutine ludcmp
!************************************************************************************************************************************






!************************************************************************************************************************************
subroutine plotthesky(bx1,bx2,by1,by2)
  implicit none
  integer, parameter :: ns=9110, nsn=80
  integer :: i,j,c(100,35),nc,snr(nsn),plcst,plstar,cf,spld,n,prslbl,rv
  real*8 :: ra(ns),dec(ns),d2r,r2d,r2h,pi,dx1,dx2,dy,ra1,dec1,rev,par
  real :: pma,pmd,vm(ns),x1,y1,x2,y2,constx(99),consty(99),r1,g1,b1,r4,g4,b4
  real :: schcon,sz1,schfac,schlbl,prinf,snlim,sllim,schmag,getmag,mag,bx1,bx2,by1,by2,x,y,mlim
  character :: cn(100)*3,con(100)*20,name*10,vsopdir*99,sn(ns)*10,snam(nsn)*10,sni*10,getsname*10,mult,var*9
  
  mlim = 6.
  cf = 2
  schmag = 0.07
  schlbl = 1.
  schfac = 1.
  schcon = 1.
  plstar = 4
  plcst = 1
  
  prinf = 150.**2
  
  x = 0.
  call pgqcr(1,r1,g1,b1) !Store colours
  call pgqcr(4,r4,g4,b4)
  call pgscr(1,1.,1.,1.) !'White' (for stars)
  call pgscr(4,x,x,1.) !Blue (for constellations)
  
  pi = 4*datan(1.d0)
  d2r = pi/180.d0
  r2d = 180.d0/pi
  r2h = 12.d0/pi
  r2h = r2d
  
  
  if(bx1.gt.bx2) then
     x = bx1
     bx1 = bx2
     bx2 = x
  end if
  
  !Read BSC
  vsopdir = '/home/sluys/diverse/popular/fortran/VSOP87/'           !Linux pc
  open(unit=20,form='formatted',status='old',file=trim(vsopdir)//'data/bsc.dat')
  rewind(20)
  do i=1,ns
     read(20,320)name,ra(i),dec(i),pma,pmd,rv,vm(i),par,mult,var
320  format(A10,1x,2F10.6,1x,2F7.3,I5,F6.2,F6.3,A2,A10)
     sn(i) = getsname(name)
  end do
  close(20)


  !Read Constellation figure data
  open(unit=40,form='formatted',status='old',file=trim(vsopdir)//'data/bsc_const.dat')
  do i=1,ns
     read(40,'(I4)',end=340,advance='no')c(i,1)
     do j=1,c(i,1)
        read(40,'(I5)',advance='no')c(i,j+1)
     end do
     read(40,'(1x,A3,A20)')cn(i),con(i)
     !Get mean star position to place const. name
     dx1 = 0.d0
     dx2 = 0.d0
     dy = 0.d0
     do j=2,c(i,1)
        dx1 = dx1 + dsin(ra(c(i,j)))
        dx2 = dx2 + dcos(ra(c(i,j)))
        dy = dy + dec(c(i,j))
     end do
     dx1 = (dx1 + dsin(ra(c(i,j))))/real(c(i,1))
     dx2 = (dx2 + dcos(ra(c(i,j))))/real(c(i,1))
     ra1 = rev(datan2(dx1,dx2))
     dec1 = (dy + dec(c(i,j)))/real(c(i,1))
     !call eq2xy(ra1,dec1,l0,b0,x1,y1)
     !constx(i) = x1
     !consty(i) = y1
     constx(i) = real(ra1*r2h)
     consty(i) = real(dec1*r2d)
  end do
340 close(40)
  nc = i-1
  
  !Read Star names
  open(unit=50,form='formatted',status='old',file=trim(vsopdir)//'data/bsc_names.dat')
  do i=1,nsn
     read(50,'(I4,2x,A10)',end=350)snr(i),snam(i)
  end do
350 close(50)
  
  
  !!Read Milky Way data
  !do f=1,5
  !   write(mwfname,'(A10,I1,A4)')'milkyway_s',f,'.dat'
  !   open(unit=60,form='formatted',status='old',file=trim(vsopdir)//'data/'//mwfname)
  !   do i=1,mwn(f)
  !      read(60,'(F7.5,F9.5)')mwa(f,i),mwd(f,i)
  !      if(maptype.eq.1) call eq2az(mwa(f,i),mwd(f,i),agst)
  !      if(maptype.eq.2) call eq2ecl(mwa(f,i),mwd(f,i),eps)
  !   end do
  !end do
  !close(60)
  
  
  !Plot constellation figures
  if(plcst.gt.0) then
     !schcon = min(max(40./sz1,0.7),3.)
     call pgsch(schfac*schcon*schlbl)
     call pgscf(cf)
     call pgsci(4)
     call pgslw(2)
     do i=1,nc
        do j=2,c(i,1)
           !call eq2xy(ra(c(i,j)),dec(c(i,j)),l0,b0,x1,y1)
           !call eq2xy(ra(c(i,j+1)),dec(c(i,j+1)),l0,b0,x2,y2)
           x1 = real(ra(c(i,j))*r2h)
           y1 = real(dec(c(i,j))*r2d)
           x2 = real(ra(c(i,j+1))*r2h)
           y2 = real(dec(c(i,j+1))*r2d)
           !if((x1*x1+y1*y1.le.prinf.or.x2*x2+y2*y2.le.prinf).and.(x2-x1)**2+(y2-y1)**2.le.90.**2) & !Not too far from centre and each other 
           if((x2-x1)**2+(y2-y1)**2.le.90.**2) & !Not too far from centre and each other 
                call pgline(2,(/x1,x2/),(/y1,y2/))
	end do
        if(constx(i).lt.bx1.or.constx(i).gt.bx2.or.consty(i).lt.by1.or.consty(i).gt.by2) cycle
        if(plcst.eq.2) call pgptext(constx(i),consty(i),0.,0.5,cn(i))
        if(plcst.eq.3) call pgptext(constx(i),consty(i),0.,0.5,con(i))
     end do
     call pgsch(schfac)
     call pgscf(cf)
  end if !if(plcst.gt.0) then
  
  !Plot stars: BSC
  spld = 0
  if(plstar.gt.0) then
     n = 0
     do i=1,ns
        if(vm(i).lt.mlim.and.vm(i).ne.0.) then
           !call eq2xy(ra(i),dec(i),l0,b0,x,y)
           x = real(ra(i)*r2h)
           y = real(dec(i)*r2d)
           if(x.lt.bx1.or.x.gt.bx2.or.y.lt.by1.or.y.gt.by2) cycle
           call pgsci(1)
           mag = getmag(vm(i),mlim)*schmag
           call pgcirc(x,y,mag)
           !write(*,'(3F10.3)')x,y,mag
           call pgsch(schfac*schlbl)
           sni = sn(i)
           !if(sni(1:1).eq.'\') call pgsch(schlbl*max(1.33,schfac))  !Greek letters need larger font
           if(sni(1:1).eq.char(92)) call pgsch(schlbl*max(1.33,schfac))  !Greek letters need larger font.  Char(92) is a \, but this way it doesn't mess up emacs' parentheses count
	   call pgsci(14)
           if(vm(i).lt.sllim) then
              if((plstar.eq.2.or.plstar.eq.5)) call pgtext(x+0.02*sz1,y+0.02*sz1,sn(i))
              if(plstar.eq.4) then !Check if the name will be printed
                 prslbl = 1
                 if(vm(i).lt.snlim) then
                    do j=1,nsn
                       if(snr(j).eq.i) prslbl = 0 !Then the name will be printed, don't print the symbol
                    end do
                 end if
                 if(prslbl.eq.1) call pgtext(x+0.02*sz1,y+0.02*sz1,sn(i))
              end if
           end if
	   spld = spld+1
	end if
     end do
     if(plstar.ge.3) then !Plot star proper names
        call pgsch(schfac*schlbl)
        do i=1,nsn
           if(vm(snr(i)).lt.max(snlim,1.4)) then  !Regulus (1.35) will still be plotted, for conjunction maps
              !call eq2xy(ra(snr(i)),dec(snr(i)),l0,b0,x,y)
              x = real(ra(snr(i)))
              y = real(dec(snr(i)))
              if(x.lt.bx1.or.x.gt.bx2.or.y.lt.by1.or.y.gt.by2) cycle
              call pgtext(x+0.02*sz1,y-0.02*sz1,snam(i))
           end if
        end do
     end if !if(plstar.eq.3) then
  end if !if(plstar.gt.0) then
  
  !Restore colours
  call pgscr(1,r1,g1,b1)
  call pgscr(4,r4,g4,b4)
  
end subroutine plotthesky
!************************************************************************************************************************************

!************************************************************************
function getsname(name)               !Get star name from bsc info
  implicit none
  character :: getsname*10,name*10,num*3,grk*3,gn*1
  num = name(1:3)
  grk = name(4:6)
  gn  = name(7:7)
  !      gn = ' '
  
  getsname = '          '
  if(grk.ne.'   ') then  !Greek letter
     if(grk.eq.'Alp') getsname = '\(2127)\u'//gn
     if(grk.eq.'Bet') getsname = '\(2128)\u'//gn
     if(grk.eq.'Gam') getsname = '\(2129)\u'//gn
     if(grk.eq.'Del') getsname = '\(2130)\u'//gn
     if(grk.eq.'Eps') getsname = '\(2131)\u'//gn
     if(grk.eq.'Zet') getsname = '\(2132)\u'//gn
     if(grk.eq.'Eta') getsname = '\(2133)\u'//gn
     if(grk.eq.'The') getsname = '\(2134)\u'//gn
     if(grk.eq.'Iot') getsname = '\(2135)\u'//gn
     if(grk.eq.'Kap') getsname = '\(2136)\u'//gn
     if(grk.eq.'Lam') getsname = '\(2137)\u'//gn
     if(grk.eq.'Mu ') getsname = '\(2138)\u'//gn
     if(grk.eq.'Nu ') getsname = '\(2139)\u'//gn
     if(grk.eq.'Xi ') getsname = '\(2140)\u'//gn
     if(grk.eq.'Omi') getsname = '\(2141)\u'//gn
     if(grk.eq.'Pi ') getsname = '\(2142)\u'//gn
     if(grk.eq.'Rho') getsname = '\(2143)\u'//gn
     if(grk.eq.'Sig') getsname = '\(2144)\u'//gn
     if(grk.eq.'Tau') getsname = '\(2145)\u'//gn
     if(grk.eq.'Ups') getsname = '\(2146)\u'//gn
     if(grk.eq.'Phi') getsname = '\(2147)\u'//gn
     if(grk.eq.'Chi') getsname = '\(2148)\u'//gn
     if(grk.eq.'Psi') getsname = '\(2149)\u'//gn
     if(grk.eq.'Ome') getsname = '\(2150)\u'//gn
  else  !Then number
     if(num(1:1).eq.' ') num = num(2:3)//' '
     if(num(1:1).eq.' ') num = num(2:3)//' '
     getsname = num//'       '
  end if
  return
end function getsname
!************************************************************************

!************************************************************************
function getmag(m,mlim)  !Determine size of stellar 'disk'
  real :: getmag,m,m1,mlim
  m1 = m
  !      if(m1.lt.0.) m1 = m1*0.5  !Less excessive grow in diameter for the brightest objects
  if(m1.lt.-1.e-3) m1 = -sqrt(-m1)  !Less excessive grow in diameter for the brightest objects
  !getmag = max(mlim-m1+0.5,0.)
  getmag = max(mlim-m1+0.5,0.5) !Make sure the weakest stars are still plotted
  !getmag = max(mlim-m1+0.5,0.)+0.5
  return
end function getmag
!************************************************************************


!************************************************************************
function rev(x)        !Returns angle in radians between 0 and 2pi
  real*8 :: x,rev,pi
  pi = 4*datan(1.d0)
  rev = x-floor(x/(2*pi))*2*pi
  return
end function rev
!************************************************************************

!************************************************************************
function rev360(x)        !Returns angle in degrees between 0 and 360
  real :: x,rev360
  rev360 = x-floor(x/(360.))*360.
  return
end function rev360
!************************************************************************

!************************************************************************
function rev24(x)        !Returns angle in hours between 0 and 24
  real :: x,rev24
  rev24 = x-floor(x/(24.))*24.
  return
end function rev24
!************************************************************************

!************************************************************************
function rev2pi(x)        !Returns angle in radians between 0 and 2pi
  real :: x,rev2pi,pi
  pi = 4*atan(1.)
  rev2pi = x-floor(x/(2.0*pi))*2.0*pi
  return
end function rev2pi
!************************************************************************

!************************************************************************
function drevpi(x)        !Returns angle in radians between 0 and pi
  use constants
  real*8 :: x,drevpi!,pi
  !pi = 4*datan(1.d0)
  drevpi = x-floor(x/pi)*pi
  return
end function drevpi
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
1 if(j1.le.n1.and.j2.le.n2)then
     d1=data1(j1)
     d2=data2(j2)
     if(d1.le.d2)then
        fn1=j1/en1
        j1=j1+1
     end if
     if(d2.le.d1)then
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
     if(dabs(term).le.eps1*termbf.or.dabs(term).le.eps2*probks)return
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
1 if(ir-l.lt.m)then
     do j=l+1,ir
        a=arr(j)
        do i=j-1,l,-1
           if(arr(i).le.a)goto 2
           arr(i+1)=arr(i)
        end do
        i=l-1
2       arr(i+1)=a
     end do
     if(jstack.eq.0)return
     ir=istack(jstack)
     l=istack(jstack-1)
     jstack=jstack-2
  else
     k=(l+ir)/2
     temp=arr(k)
     arr(k)=arr(l+1)
     arr(l+1)=temp
     if(arr(l).gt.arr(ir))then
        temp=arr(l)
        arr(l)=arr(ir)
        arr(ir)=temp
     end if
     if(arr(l+1).gt.arr(ir))then
        temp=arr(l+1)
        arr(l+1)=arr(ir)
        arr(ir)=temp
     end if
     if(arr(l).gt.arr(l+1))then
        temp=arr(l)
        arr(l)=arr(l+1)
        arr(l+1)=temp
     end if
     i=l+1
     j=ir
     a=arr(l+1)
3    continue
     i=i+1
     if(arr(i).lt.a)goto 3
4    continue
     j=j-1
     if(arr(j).gt.a)goto 4
     if(j.lt.i)goto 5
     temp=arr(i)
     arr(i)=arr(j)
     arr(j)=temp
     goto 3
5    arr(l+1)=arr(j)
     arr(j)=a
     jstack=jstack+2
     !if(jstack.gt.nstack)pause 'nstack too small in sort'
     if(jstack.gt.nstack) write(*,'(A)')' nstack too small in dindexx'
     if(ir-i+1.ge.j-l)then
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


!************************************************************************
function timestamp()  !Get time stamp in seconds since 1970-01-01 00:00:00 UTC
  implicit none
  real*8 :: timestamp
  integer :: i,system

  i = system('date +%s.%N &> ~/.tmp_timestamp')
  !i = system('date +%s &> ~/.tmp_timestamp') !%N for fractional seconds doesn't work on MacOS!!! (But it does with GNU date)
  open(unit=9,status='old',file='~/.tmp_timestamp')
  read(9,'(F20.9)')timestamp
  !print*,timestamp
  close(9)
  i = system('rm -f ~/.tmp_timestamp')
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
subroutine identify_2d_ranges(ni,ivals,nx,ny,z)
  !Get the 2d probability intervals; z lies between 1 (in 100% range) and ni (in lowest-% range, e.g. 90%)
  implicit none
  integer :: ni,nx,ny,nn,indx(nx*ny),i,b,ib
  real :: ivals(ni),z(nx,ny),x1(nx*ny),x2(nx*ny),tot,np
  
  nn = nx*ny
  x1 = reshape(z,(/nn/))  !x1 is the 1D array with the same data as the 2D array z
  call rindexx(nn,-x1(1:nn),indx(1:nn)) ! -x1: sort descending
  
  np = sum(z)
  tot = 0.
  do b=1,nn !Loop over bins in 1D array
     ib = indx(b)
     x2(ib) = 0.
     if(x1(ib).eq.0.) cycle
     tot = tot + x1(ib)
     do i=ni,1,-1 !Loop over intervals
        if(tot.le.np*ivals(i)) x2(ib) = real(ni-i+1)  !e.g. x2(b) = ni if within 90%, ni-1 if within 95%, etc, and 1 if within 100%
        !write(*,'(2I4, F6.2, 3F20.5)')b,i, ivals(i), np,tot,np*ivals(i)
     end do
  end do
  
  z = reshape(x2, (/nx,ny/))  ! z lies between 1 and ni
end subroutine identify_2d_ranges
!************************************************************************



!************************************************************************
!Compute 2D probability areas
subroutine calc_2d_areas(p1,p2,changevar,ni,nx,ny,z,tr,area)
  implicit none
  integer :: p1,p2,changevar,ni,nx,ny,ix,iy,i,i1,iv
  real :: z(nx,ny),tr(6),y,dx,dy,d2r,area(ni)
  
  d2r = atan(1.)/45.
  area = 0.
  
  !print*,ni,nx,ny
  do ix = 1,nx
     do iy = 1,ny
        dx = tr(2)
        dy = tr(6)
        if(changevar.eq.1 .and. (p1.eq.8.and.p2.eq.9 .or. p1.eq.12.and.p2.eq.11) ) then !Then: RA-Dec or phi/theta_Jo plot, convert lon -> lon * 15 * cos(lat)
           !x = tr(1) + tr(2)*ix + tr(3)*iy
           !y = tr(4) + tr(5)*ix + tr(6)*iy
           y = tr(4) + tr(6)*iy
           if(p1.eq.8) then
              dx = dx*cos(y*d2r)
           else if(p1.eq.12) then
              dx = dx*abs(sin(y*d2r))  !Necessary for i-psi plot?
           end if
           !print*,p1,y,cos(y*d2r)
           if(p1.eq.8) dx = dx*15
        end if
        iv = nint(z(ix,iy))
        !if(iv.gt.0) area(iv) = area(iv) + dx*dy
        do i=1,ni
           if(iv.ge.i) then
              i1 = ni-i+1
              area(i1) = area(i1) + dx*dy
           end if
        end do
        !if(iv.eq.3) write(*,'(7F10.2)')x,y,dx,dy,dx*dy,z(ix,iy),area(iv)
     end do
  end do
end subroutine calc_2d_areas
!************************************************************************


!************************************************************************
function truerange2d(z,nx,ny,truex,truey,tr)
  !Get the smallest probability area in which the true values lie
  implicit none
  integer :: nx,ny,ix,iy,truerange2d
  real :: truex,truey,z(nx,ny),tr(6)
  
  !x = tr(1) + tr(2)*ix + tr(3)*iy
  !y = tr(4) + tr(5)*ix + tr(6)*iy
  ix = floor((truex - tr(1))/tr(2))
  iy = floor((truey - tr(4))/tr(6))
  if(ix.lt.1.or.ix.gt.nx.or.iy.lt.1.or.iy.gt.ny) then
     truerange2d = 0
  else
     truerange2d = nint(z(ix,iy))
  end if
end function truerange2d
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
  real*8 :: vec1(3),vec2(3),crpr(3),veclen
  crpr(1) = vec1(2)*vec2(3) - vec1(3)*vec2(2)
  crpr(2) = vec1(3)*vec2(1) - vec1(1)*vec2(3)
  crpr(3) = vec1(1)*vec2(2) - vec1(2)*vec2(1)
  !write(*,'(3ES13.3)')veclen(vec1),veclen(vec2),veclen(crpr)
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
  
  !write(*,'(3ES13.3)')p
  !write(*,'(3ES13.3)')o
  
  call crossproduct(p,o,x1)
  call crossproduct(x1,p,o1) !o1: projection of o in the plane of the sky
  
  z = (/0.d0,0.d0,1.d0/) !Vertical normal vector
  call crossproduct(p,z,x1)
  call crossproduct(x1,p,z1) !z1: projection of z in the plane of the sky
  
  call normvec(o1)
  call normvec(z1)
  posangle = dacos(dotproduct(o1,z1))
  !write(*,'(3ES13.3)')p
  !write(*,'(3ES13.3)')o
  !write(*,'(3ES13.3)')o1
  !write(*,'(3ES13.3)')z1
  !write(*,'(3ES13.3)')dotproduct(o1,z1),dacos(dotproduct(o1,z1))
end function posangle
!************************************************************************


!************************************************************************
subroutine compute_incli_polang(pl,pb,ol,ob, i,psi) !Compute the inclination and polarisation angle for a source with position (pl,pb) and orientation (ol,ob)
  use constants
  implicit none
  !pl,ol in [0,2pi[;  pb,ob in [-pi,pi]
  real*8 :: pl,pb,ol,ob
  real*8 :: p(3),o(3),i,dotproduct,psi,polangle,drevpi
  
  call ang2vec(pl,pb,p)       !Position normal vector
  call ang2vec(ol,ob,o)       !Orientation normal vector
  !i = pi2 - dacos(dotproduct(p,o))  !Compute inclination angle: <0: points towards us, >0 points away from us
  i = dacos(dotproduct(p,o))  !Compute inclination angle: 0: points exactly away from us, 180 points exactly towards us, 90: in the plane of the sky
  !i = dotproduct(p,o)         !Compute cos(inclination angle): 1: points exactly away from us, -1 points exactly towards us, 0: in the plane of the sky
  psi = polangle(p,o)         !Compute polarisation angle
  !psi = drevpi(polangle(p,o))  !Compute polarisation angle
  
end subroutine compute_incli_polang
!************************************************************************

!************************************************************************
subroutine compute_incli_posang(pl,pb,ol,ob, i,psi) !Compute the inclination and position angle for a source with position (pl,pb) and orientation (ol,ob)
  use constants
  implicit none
  !pl,ol in [0,2pi[;  pb,ob in [-pi,pi]
  real*8 :: pl,pb,ol,ob
  real*8 :: p(3),o(3),i,dotproduct,psi,posangle,drevpi
  
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
  
  detcoords(1,:) = (/-119.41,46.45/)  !H1; l,b
  detcoords(2,:) = (/-90.77,30.56/)   !L1
  detcoords(3,:) = (/10.50,43.63/)    !V
  
  call ang2vec(detcoords(d1,1),detcoords(d1,2),vec1)
  call ang2vec(detcoords(d2,1),detcoords(d2,2),vec2)
  
  dvec = vec2 - vec1
  
  call vec2ang(dvec,l,b)  !Searched point is in zenith/nadir for an observer on this location on the globe
  
end subroutine detectorvector
!************************************************************************

