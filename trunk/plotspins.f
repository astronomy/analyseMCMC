!Read and plot the data output from the spinning MCMC code.

program plotspins
  use plotspins_settings
  implicit none
  integer, parameter :: narr1=2.01e5+2,npar0=13,nival1=5,nr1=5,nstat1=10,ndets=3
  integer :: n(nchs),ntot(nchs),n0,n1,n2,i,j,j1,j2,nburn0(nchs),iargc,io,readerror,pgopen,system,narr,maxdots,offsetrun,imin
  integer :: niter(nchs),totiter,totpts,totlines,seed(nchs),ndet(nchs),totthin(nchs),contrchains
  integer :: index(npar1,nchs*narr1),index1(nchs*narr1),fileversion
  real :: is(nchs,narr1),isburn(nchs),jumps(nchs,npar1,narr1),startval(nchs,npar1,3)
  real :: sig(npar1,nchs,narr1),acc(npar1,nchs,narr1),avgtotthin
  real :: x(nchs,nchs*narr1),xx(nchs*narr1),yy(nchs*narr1),zz(nchs*narr1)
  real,allocatable :: xbin(:,:),ybin(:,:),xbin1(:),ybin1(:),ybin2(:),ysum(:),yconv(:),ycum(:)  !These depend on nbin1d, allocate after reading input file
  real :: a,r2d,r2h,rat,plx,ply
  real*8 :: t,t0,pi,tpi,dvar,dvar1,dvar2,nullh
  real, allocatable :: dat(:,:,:),alldat(:,:,:),pldat(:,:,:)
  character :: varnames(npar1)*8,pgunits(npar1)*99,pgvarns(npar1)*99,pgvarnss(npar1)*99,pgorigvarns(npar1)*99,infile*100,infiles(nchs)*100,str*99,str1*99,str2*99,bla*10,command*99
  
  integer :: nn,nn1,lowvar(npar1),nlowvar,highvar(npar1),nhighvar,ntotrelvar,nlogl1,nlogl2,ksn1,ksn2,iloglmax,icloglmax,ci
  real*8 :: chmean(nchs,npar1),totmean(npar1),chvar(npar1),chvar1(nchs,npar1),totvar(npar1),rhat(npar1),totrelvar,ksdat1(narr1),ksdat2(narr1),ksd,ksprob,loglmax,loglmaxs(nchs)
  character :: ch
  
  integer :: samplerate(nchs,ndets),samplesize(nchs,ndets),FTsize(nchs,ndets),detnr(nchs,ndets)
  real :: snr(nchs,ndets),flow(nchs,ndets),fhigh(nchs,ndets),t_before(nchs,ndets),t_after(nchs,ndets),deltaFT(nchs,ndets)
  real*8 :: FTstart(nchs,ndets)
  character :: detname(nchs,ndets)*14,string*99
  
  integer :: npdf,ncont!,nfx,nfy,fx,fy
  real :: x0,x1,x2,y1,y2,dx,dy,xmin,xmax,ymin,ymax,xmin1,xmax1,xpeak!,ymin1,ymax1,ymaxs(nchs+2)
  real,allocatable :: z(:,:),zs(:,:,:)  !These depend on nbin2d, allocate after reading input file
  real :: coefs(100),coefs1(100),cont(11),tr(6)
  
  integer :: nchains,nchains0,ic,i0,i1,lw,lw2
  real :: sch
  character :: outputname*99,outputdir*99,psclr*4,colournames(15)*20
  
  integer :: o,p,p1,p2,par1,par2,nr,c,c0,nstat,wrap(nchs,npar1),nival,npar,ncolours,colours(10),nsymbols,symbols(10),symbol,defcolour,plotthis,tempintarray(99)
  integer :: truerange2d,trueranges2d(npar1,npar1)
  real :: range,minrange,range1,range2,drange,maxgap,ranges(nchs,nival1+1,npar1,nr1),ival,ivals(nival1+1),centre,shift(nchs,npar1),plshift
  real :: probarea(nival1+1),probareas(npar1,npar1,nival1+1,2)
  real :: median,medians(npar1),mean(npar1),stdev1(npar1),stdev2(npar1),var1(npar1),var2(npar1),absvar1(npar1),absvar2(npar1)
  real :: stats(nchs,npar1,nstat1),corrs(npar1,npar1),acorrs(nchs,0:npar1,0:narr1)
  real :: norm
  
  real :: bmpsz,bmprat!,scrsz,scrrat,pssz,psrat
  real :: rev360,rev24,rev2pi
  character :: bmpxpix*99,unsharplogl*99,unsharpchain*99,unsharppdf1d*99,unsharppdf2d*99
  
  integer :: iframe,nplt
  real*8 :: timestamp,ts1,ts2
  character :: framename*99,tms*8,upline*4
  
  pi = 4*datan(1.d0)
  tpi = 2*pi
  r2d = real(180.d0/pi)
  r2h = real(12.d0/pi)
  
  call set_plotsettings()  !Set plot settings to 'default' values
  call read_inputfile()    !Read the plot settings (overwrite the defaults)
  call write_inputfile()   !Write the input file back to disc
  
  
  nchains0 = iargc()
  if(nchains0.lt.1) then
     write(*,'(A,/)')'  Syntax: plotspins <file1> [file2] ...'
     stop
  end if
  
  
  !Some of the stuff below will have to go to the input file
  
  par1 = 1          !First parameter to treat (stats, plot): 0-all
  par2 = 15         !Last parameter to treat (0: use npar)
  
  !whitebg = 1       !White background for screen and png plots: 0-no, 1-yes
  !scrsz  = 10.8     !Screen size for X11 windows:  MacOS: 16.4, Gentoo: 10.8
  !scrrat = 0.57     !Screen ratio for X11 windows, MacBook: 0.57
  !bmpsz  = 12.      !Size for bitmap:  10.6
  !bmprat = 0.70     !Ratio for bitmap: 0.75
  !!bmprat = 1.25    !Ratio for bitmap: 0.75
  !!bmpxpix = '1000'  !Final size of converted bitmap output in pixels (string)
  
  !scfac = 1.2
  bmpsz = real(bmpxsz-1)/85. * scfac !Make png larger, so that convert interpolates and makes the plot smoother
  bmprat = real(bmpysz-1)/real(bmpxsz-1)
  write(bmpxpix,'(I4)')bmpxsz  !Used as a text string by convert
  !print*,bmpxsz,bmpysz,bmpsz,bmprat,trim(bmpxpix)
  
  !Use full unsharp-mask strength for plots with many panels and dots, weaker for those with fewer panels and.or no dots
  write(unsharplogl,'(I4)')max(nint(real(unsharp)/2.),1)  !Only one panel with dots
  write(unsharpchain,'(I4)')unsharp                       !~12 panels with dots
  write(unsharppdf1d,'(I4)')max(nint(real(unsharp)/2.),1) !~12 panels, no dots
  write(unsharppdf2d,'(I4)')max(nint(real(unsharp)/4.),1) !1 panel, no dots
  
  !Trying to implement this, for chains plot at the moment
  !pssz   = 10.5   !Default: 10.5   \__ Gives same result as without pgpap
  !psrat  = 0.742  !Default: 0.742  /
  
  !if(quality.eq.2) then
  !   pssz   = 10.5
  !   psrat  = 0.82     !Nice for presentation (Beamer)
  !end if
  
  outputdir = '.'  !Directory where output is saved (either relative or absolute path)
  
  !NEW columns in dat: 1:logL 2:mc, 3:eta, 4:tc, 5:logdl, 6:spin, 7:kappa, 8: RA, 9:sindec,10:phase, 11:sinthJ0, 12:phiJ0, 13:alpha
  !Choose plot variables.  Number of variables to plot: 1,2,3,4,5,6,8,9,10,12,15,16
  !nplvar = 1;  plvars(1:nplvar) = (/2/)
  !nplvar = 1;  plvars(1:nplvar) = (/14/)
  !nplvar = 2;  plvars(1:nplvar) = (/2,3/)
  !nplvar = 4;  plvars(1:nplvar) = (/2,3,14,15/) !All masses
  !nplvar = 4;  plvars(1:nplvar) = (/2,3,6,7/) !Masses, spins
  !nplvar = 4;  plvars(1:nplvar) = (/4,5,8,9/) !tc, d, position
  !nplvar = 8;  plvars(1:nplvar) = (/4,5,8,9,10,11,12,13/) !All except masses, spins
  !nplvar = 9;  plvars(1:nplvar) = (/2,3,4,6,7,5,10,8,9/)
  !nplvar = 12;  plvars(1:nplvar) = (/2,3,4,5, 6,7,8,9, 10,11,12,13/) !Default: all parameters
  !nplvar = 12;  plvars(1:nplvar) = (/2,3,4, 6,7,5, 8,9,10, 11,12,13/) !For poster (3x4 rather than 4x3)
  !nplvar = 14;  plvars(1:nplvar) = (/2,3,4,5,6,7,8,9,10,11,12,13,14,15/) !All 12 + m1,m2
  !nplvar = 15;  plvars(1:nplvar) = (/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15/)
  !nplvar = 6;  plvars(1:nplvar) = (/1,4,8,9,11,12/)
  !nplvar = 9;  plvars(1:nplvar) = (/1,2,3,4,5,6,7,8,9/) !logL + 8 most important parameters
  !nplvar = 2;  plvars(1:nplvar) = (/1,2/)
  
  panels(1) = 0 ! 0 - use default values, >0 have panels(1) panels in the horizontal direction   \
  panels(2) = 0 ! 0 - use default values, >0 have panels(2) panels in the vertical direction     / Default values are use if either value = 0
  
  
  
  
  
  
  
  !Sort out implicit options:
  if(panels(1)*panels(2).lt.nplvar) panels = 0
  
  psclr = '/cps'
  if(colour.eq.0) psclr = '/ps '
  
  ncolours = 5; colours(1:ncolours)=(/4,2,3,6,5/) !Paper
  ncolours = 10; colours(1:ncolours)=(/2,3,4,5,6,7,8,9,10,11/)
  nsymbols = 1; symbols(1:nsymbols)=(/chainsymbol/)
  if(colour.eq.1.and.quality.eq.2) then !Beamer
     ncolours = 5
     colours(1:ncolours)=(/4,2,5,11,15/)
  end if
  if(colour.ne.1) then
     ncolours=3
     colours(1:ncolours)=(/1,14,15/)
     !ncolours=6
     !colours(1:ncolours)=(/1,1,1,15,1,15/)
     if(chainsymbol.eq.-10) then
        nsymbols = 8
        symbols(1:nsymbols) = (/2,4,5,6,7,11,12,15/) !Thin/open symbols
     end if
     if(chainsymbol.eq.-11) then
        nsymbols = 6
        symbols(1:nsymbols) = (/-3,-4,16,17,18,-6/) !Filled symbols
     end if
     !print*,chainsymbol,nsymbols
  end if
  if(colour.eq.1.and.quality.eq.0.and.nchs.gt.5) then
     ncolours = 10
     colours(1:ncolours)=(/2,3,4,5,6,7,8,9,10,11/)
  end if
  !Overrule
  !ncolours = 1
  !colours(1:ncolours)=(/6/)
  !defcolour = 2 !Red e.g. in case of 1 chain
  defcolour = colours(1)
  
  
  if(reverseread.ge.2) then !Reverse colours too
     do i=1,ncolours
        tempintarray(i) = colours(i)
     end do
     !do i=1,ncolours
     !   colours(i) = tempintarray(ncolours-i+1) !Reverse colours too
     !end do
     do i=1,nchains0
        colours(i) = tempintarray(nchains0-i+1) !Reverse colours too, but use the same first nchains0 from the set
     end do
  end if
  
  if(plot.eq.0) then
     pllogl = 0
     plchain = 0
     pljump = 0
     plsigacc = 0
     if(savepdf.eq.0) then
        plpdf1d = 0
        plpdf2d = 0
     end if
     plmovie = 0
  end if
  if(savepdf.eq.1) then
     nplvar = 15; plvars(1:nplvar) = (/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15/) !All 12 + m1,m2
     wrapdata = 0
  end if
  !if(par1.lt.2) par1 = 2
  if(par1.lt.1) par1 = 1 !Include log(L)
  if(plsigacc.ge.1.or.plmovie.ge.1) rdsigacc = 1
  if(file.eq.1) combinechainplots = 0
  if(file.ge.1) update = 0
  if(plmovie.eq.1) update = 0
  if(plotsky.ge.1) plpdf2d = 1
  
  colournames(1:15) = (/'white','red','dark green','dark blue','cyan','magenta','yellow','orange','light green','brown','dark red','purple','red-purple','dark grey','light grey'/)
  if(file.ge.2) colournames(1) = 'black'
  
  
  
  !NEW columns in dat: 1:logL 2:mc, 3:eta, 4:tc, 5:logdl, 6:spin, 7:kappa, 8: RA, 9:sindec,10:phase, 11:sinthJ0, 12:phiJ0, 13:alpha
  varnames(1:15) = (/'logL','Mc','eta','tc','log_dl','spin','kappa','RA','sin_dec','phase','sin_thJo','phJo','alpha','M1','M2'/)
  pgvarns(1:15)  = (/'log Likelihood        ','M\dc\u (M\d\(2281)\u) ','\(2133)               ','t\dc\u (s)            ', &
       'logd\dL\u (Mpc)       ','a\dspin\u             ','\(2136)               ','R.A. (rad)            ', &
       'sin dec.              ','\(2147)\dc\u (rad)    ','sin \(2134)\dJ0\u     ','\(2147)\dJ0\u (rad)   ', &
       '\(2127)\dc\u (rad)    ','M\d1\u (M\d\(2281)\u) ','M\d2\u (M\d\(2281)\u) '/)
  pgvarnss(1:15)  = (/'log L    ','M\dc\u ','\(2133)','t\dc\u','log d\dL\u','a\dspin\u','\(2136)','R.A.','sin dec.','\(2147)\dc\u', &
       'sin \(2134)\dJ0\u','\(2147)\dJ0\u','\(2127)\dc\u','M\d1\u','M\d2\u'/)
  pgorigvarns(1:15)  = (/'log Likelihood        ','M\dc\u (M\d\(2281)\u) ','\(2133)               ','t\dc\u (s)            ', &
       'logd\dL\u (Mpc)       ','a\dspin\u             ','\(2136)               ','R.A. (rad)            ', &
       'sin dec.              ','\(2147)\dc\u (rad)    ','sin \(2134)\dJ0\u     ','\(2147)\dJ0\u (rad)   ', &
       '\(2127)\dc\u (rad)    ','M\d1\u (M\d\(2281)\u) ','M\d2\u (M\d\(2281)\u) '/)
  !pgorigvarns(1:15)  = (/'log L    ','M\dc\u ','\(2133)','t\dc\u','log d\dL\u','a\dspin\u','\(2136)','R.A.','sin dec.','\(2147)\dc\u', &
  !     'sin \(2134)\dJ0\u','\(2147)\dJ0\u','\(2127)\dc\u','M\d1\u','M\d2\u'/)
  pgunits(1:15)  = (/'','M\d\(2281)\u ','','s','Mpc','','rad','rad','','rad','','rad','rad','M\d\(2281)\u','M\d\(2281)\u'/)
  
  
  !nival = 1; ivals(1:nival) = (/0.9/)  !Number of ranges and values of 'interval levels', one of them should be 0.90. 1.00 is added automatically.
  nival = 2; ivals(1:nival) = (/0.9,0.99/)  !Number of ranges and values of 'interval levels', one of them should be 0.90. 1.00 is added automatically.
  !nival = 3; ivals(1:nival) = (/0.683,0.9,0.997/)  !Number of ranges and values of 'interval levels', one of them should be 0.90. 1.00 is added automatically.
  !nival = 3; ivals(1:nival) = (/0.50,0.7,0.9/)  !Number of ranges and values of 'interval levels', one of them should be 0.90. 1.00 is added automatically.
  !nival = 3; ivals(1:nival) = (/0.90,0.95,0.99/)  !Number of ranges and values of 'interval levels', one of them should be 0.90. 1.00 is added automatically.
  !nival = 4; ivals(1:nival) = (/0.6827,0.9,0.9545,0.9973/)  !Number of ranges and values of 'interval levels', one of them should be 0.90. 1.00 is added automatically.
  !nival = 4; ivals(1:nival) = (/0.8,0.9,0.95,0.99/)  !Number of ranges and values of 'interval levels', one of them should be 0.90. 1.00 is added automatically.
  ivals(nival+1) = 1. !The last interval is always 100%
  
  
  
  upline = char(27)//'[2A'  !Printing this makes the cursor move up one line (actually 2, for some reason that's needed)
  
  if(prprogress+prruninfo+prinitial.ge.1) write(*,*)
  j=0
  do i=1,nival
     if(abs(ivals(i)-ival0).lt.1.e-6) j=1
  end do
  if(j.eq.0) write(*,'(A44,F5.3,A35)')' !!! Error:  standard probability interval (',ival0,') is not present in array ivals !!!'
  
  npar = 13
  if(prprogress.ge.1) then
     if(nchains0.gt.nchs) then
        write(*,'(A,I3,A)')'  Too many chains, please increase nchs. Only',nchs,' files will be read.'
     else
        write(*,'(A,I3,A)')'  Reading',nchains0,' input files: '
     end if
  end if
  nchains0 = min(nchains0,nchs)
  nchains = nchains0
  
  
  
  
  
  
  
  !*******************************************************************************************************************************
  !***   READ INPUT FILE(S)   ****************************************************************************************************
  !*******************************************************************************************************************************
  
  !fileversion = 1
  fileversion = 2  !Start trying the new format
101 continue
  readerror = 0
  narr = narr1
  allocate(dat(npar1,nchains,narr1))
  
  do ic = 1,nchains0
     if(reverseread.eq.0) then
        call getarg(ic,infile) !Read file name from the command-line arguments
     else
        call getarg(nchains0-ic+1,infile) !Read file name from the command-line arguments in reverse order
     end if
     infiles(ic) = infile
     
     open(unit=10,form='formatted',status='old',file=trim(infile),iostat=io)
     if(io.ne.0) then
        write(*,'(A)')'File not found: '//trim(infile)//', aborting.'
        goto 9999
     end if
     rewind(10)
     
     !if(prprogress.ge.2) write(*,'(A,I3,A,I3,A20,$)')'    File',ic,':  '//trim(infile)//'    Using colour',colours(mod(ic-1,ncolours)+1),': '//colournames(colours(mod(ic-1,ncolours)+1))
     
     !NEW columns in dat: 1:logL 2:mc, 3:eta, 4:tc, 5:logdl, 6:spin, 7:kappa, 8: RA, 9:sindec,10:phase, 11:sinthJ0, 12:phiJ0, 13:alpha
     !Read the headers
     read(10,*,end=199,err=199)bla
     read(10,'(I10,I12,I8,F22.10,I8)')niter(ic),nburn0(ic),seed(ic),nullh,ndet(ic)
     read(10,*,end=199,err=199)bla
     do i=1,ndet(ic)
        read(10,'(2x,A14,F18.8,4F12.2,F22.8,F17.7,3I14)') detname(ic,i),snr(ic,i),flow(ic,i),fhigh(ic,i),t_before(ic,i),t_after(ic,i),FTstart(ic,i),deltaFT(ic,i),samplerate(ic,i),samplesize(ic,i),FTsize(ic,i)
        if(detname(ic,i).eq.'       Hanford') detnr(ic,i) = 1
        if(detname(ic,i).eq.'    Livingston') detnr(ic,i) = 2
        if(detname(ic,i).eq.'          Pisa') detnr(ic,i) = 3
     end do
     !if(prprogress.ge.1.and.update.eq.0) write(*,*)''
     read(10,*,end=199,err=199)bla
     if(fileversion.eq.1) read(10,*,end=199,err=199)bla
     
     !Read the data (old format, <04/2008)
     if(fileversion.eq.1.and.rdsigacc.eq.1) then
        i=1
        do while(i.le.narr1)
           read(10,'(I12,F21.10,2(F17.10,F10.6,F7.4),F22.10,F10.6,F7.4,9(F17.10,F10.6,F7.4))',end=199,err=198)i1,dat(1,ic,i),(dat(j,ic,i),sig(j,ic,i),acc(j,ic,i),j=2,3),  t,sig(4,ic,i),acc(4,ic,i),  (dat(j,ic,i),sig(j,ic,i),acc(j,ic,i),j=5,npar0)
           !read(10,*,end=199,err=198)i1,dat(1,ic,i),(dat(j,ic,i),sig(j,ic,i),acc(j,ic,i),j=2,3),  t,sig(4,ic,i),acc(4,ic,i),  (dat(j,ic,i),sig(j,ic,i),acc(j,ic,i),j=5,npar0)  !Only works when not using new output (otherwise, no read error)
           is(ic,i) = real(i1)
           if(ic.eq.1.and.i.eq.1) t0 = dble(floor(t/10.d0)*10)
           dat(4,ic,i) = real(t - t0)
           if(thin.gt.1.and.i.gt.2) then !'Thin' the output by reading every thin-th line
              do j=1,thin-1
                 read(10,*,end=199,err=198)bla
              end do
           end if
           i = i+1
           if(i1.ge.maxchlen) exit
        end do !i
     else if(fileversion.eq.1.and.rdsigacc.eq.0) then !Read 40% quicker
        i=1
        do while(i.le.narr1)
           read(10,'(I12,F21.10,2(F17.10,17x),F22.10,17x,9(F17.10,17x))',end=199,err=198)i1,dat(1,ic,i),(dat(j,ic,i),j=2,3),  t,  (dat(j,ic,i),j=5,npar0)
           is(ic,i) = real(i1)
           if(ic.eq.1.and.i.eq.1) t0 = dble(floor(t/10.d0)*10)
           dat(4,ic,i) = real(t - t0)
           !if(abs(dat(7,ic,i)).gt.0.995) i=i-1
           if(thin.gt.1.and.i.gt.2) then !'Thin' the output by reading every thin-th line
              do j=1,thin-1
                 read(10,*,end=199,err=198)bla
              end do
           end if
           i = i+1
           if(i1.ge.maxchlen) exit
        end do !i
        
     else if(fileversion.eq.2) then !New file format (04/2008)
        i=1
        do while(i.le.narr1)
           !read(10,'(I10,F14.6,2(F13.7),F20.8,9(F12.7))',end=199,err=198)i1,dat(1,ic,i),(dat(j,ic,i),j=2,3),  t,  (dat(j,ic,i),j=5,npar0)
           read(10,'(I10,F14.6,2(F13.7),F20.8,9(F12.7))',iostat=io)i1,dat(1,ic,i),(dat(j,ic,i),j=2,3),  t,  (dat(j,ic,i),j=5,npar0)
           if(io.lt.0) exit !EOF
           if(io.gt.0) then !Read error
              if(readerror.eq.1) goto 198 !Read error in previous line as well
              readerror = 1
              i = i-1
              cycle
           end if
           readerror = 0
           is(ic,i) = real(i1)
           if(ic.eq.1.and.i.eq.1) t0 = dble(floor(t/10.d0)*10)
           dat(4,ic,i) = real(t - t0)
           if(thin.gt.1.and.i.gt.2) then !'Thin' the output by reading every thin-th line
              do j=1,thin-1
                 read(10,*,end=199,err=198)bla
              end do
           end if
           i = i+1
           if(i1.ge.maxchlen) exit
        end do !i
     end if
     
     !if(i1.lt.maxchlen) write(*,'(A,$)')'   *** WARNING ***   Not all lines in this file were read    '
     if(i.ge.narr1-2) write(*,'(A,$)')'   *** WARNING ***   Not all lines in this file were read    '
     goto 199
198  write(*,'(A,I7,$)')'  Read error in file '//trim(infile)//', line',i
     i = i-1
     if(i.lt.25) then
        !if(ic.eq.1.and.fileversion.eq.1) then
        !   fileversion = 2
        !   deallocate(dat)
        !   write(*,'(A)')"    I'll try the new file format..."
        !   goto 101
        !end if
        if(ic.eq.1.and.fileversion.eq.2) then
           fileversion = 1
           deallocate(dat)
           write(*,'(A)')"    I'll try the old file format..."
           goto 101
        end if
        write(*,'(A,/)')'  Aborting program...'
        stop
     end if
     write(*,*)
199  close(10)
     ntot(ic) = i-1
     n(ic) = ntot(ic) !n can be changed in rearranging chains, ntot won't be changed
     !if(prprogress.ge.2.and.update.ne.1) write(*,'(1x,3(A,I9),A1)')' Lines:',ntot(ic),', iterations:',nint(is(ic,ntot(ic))),', burn-in:',nburn(ic),'.'
  end do !do ic = 1,nchains0
  
  
  totiter = 0
  do ic = 1,nchains0
     totiter = totiter + nint(is(ic,ntot(ic)))
  end do
  !if(prprogress.ge.2.and.update.ne.1) write(*,'(A10,65x,2(A,I9),A1)')'Total:',' Lines:',sum(ntot(1:nchains0)),', iterations:',totiter
  
  
  
  
  !Print run info (detectors, SNR, amount of data, FFT, etc)
  if(prruninfo.gt.0.and.update.eq.0) then
     write(*,'(/,A)')'  Run inormation:'
     do ic = 1,nchains0
        if((prruninfo.eq.1.and.ic.eq.1) .or. prruninfo.eq.2) then
           write(*,'(4x,A10,A12,A8,A22,A8)')'niter','nburn','seed','null likelihood','ndet'
           write(*,'(4x,I10,I12,I8,F22.10,I8)')niter(ic),nburn0(ic),seed(ic),nullh,ndet(ic)
           write(*,'(A14,A3,A18,4A12,A22,A17,3A14)')'Detector','Nr','SNR','f_low','f_high','before tc','after tc','Sample start (GPS)','Sample length','Sample rate','Sample size','FT size'
           
           do i=1,ndet(ic)
              write(*,'(A14,I3,F18.8,4F12.2,F22.8,F17.7,3I14)')detname(ic,i),detnr(ic,i),snr(ic,i),flow(ic,i),fhigh(ic,i),t_before(ic,i),t_after(ic,i),FTstart(ic,i),deltaFT(ic,i),samplerate(ic,i),samplesize(ic,i),FTsize(ic,i)
           end do
           write(*,*)
        end if
     end do !do ic = 1,nchains0
  end if  !prruninfo.gt.0
  
  narr = maxval(n(1:nchains0))
  
  
  
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
              write(*,'(A,I3)')'   *** WARNING ***  Nburn larger than Nchain, setting nburn to 10% for chain',ic
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
        if(dat(1,ic,i).gt.loglmaxs(ic)) then
           loglmaxs(ic) = dat(1,ic,i)
           if(dat(1,ic,i).gt.loglmax) then
              loglmax = dat(1,ic,i)
              iloglmax = i
              icloglmax = ic
           end if
        end if
     end do
  end do
  
  !Autoburnin: for each chain, get the first point where log(L) > log(L_max)-autoburnin
  if(autoburnin.gt.1.e-10) then
     loop1: do ic=1,nchains0
        isburn(ic) = is(ic,ntot(ic)) !Set burnin to last iteration, so that chain is completely excluded if condition is never fulfilled
        nburn(ic) = ntot(ic)
        do i=2,ntot(ic) !i=1 is true value?
           if(dat(1,ic,i).gt.real(loglmax)-autoburnin) then
              isburn(ic) = is(ic,i)
              nburn(ic) = i
              cycle loop1
           end if
        end do
     end do loop1
  end if
  
  
  
  
  !Print info on number of iterations, burnin, thinning, etc.
  do ic=1,nchains0
     !if(prruninfo.ge.1.and.update.ne.1) write(*,'(2x,A5,I3,A,I9,A,ES7.1,A,I9,A,ES7.1,A,I4,A,ES9.1,A,ES9.1,A,I4)')'File',ic,':  Lines:',n(ic),'  (',real(n(ic)),'),  Iterations:',nint(is(ic,n(ic))),'  (',is(ic,n(ic)),'),  thinning in file:',nint(is(ic,n(ic))/real(n(ic)*max(thin,1))),'     Burn-in:  line:',real(nburn(ic)),', iteration:',isburn(ic),', total thinning:',totthin(ic)
     !if(prprogress.ge.2.and.update.ne.1) write(*,'(A9,I3,A2,A23,A12,A,ES7.1,A,ES7.1,A,ES7.1,A,ES7.1,A,I3,A,I4,A,ES8.2,A1)') &
     !     'Chain',ic,': ',trim(infiles(ic)),', '//colournames(colours(mod(ic-1,ncolours)+1))//'.', '  Lines: ',real(n(ic)),', iter: ',is(ic,n(ic)), &
     !     '.  Burn-in: line: ',real(nburn(ic)),', iter: ',isburn(ic),'.  Thin: file:',nint(is(ic,n(ic))/real(n(ic)*max(thin,1))),', tot:',totthin(ic),'.  Data pts: ',real(n(ic)-nburn(ic)),'.'
     if(prprogress.ge.2.and.update.ne.1) write(*,'(A9,I3,A2,A23,A12,A,ES7.1,A,ES7.1,A,ES7.1,A,ES7.1, A,F8.2, A,I3,A,I4,A,ES8.2,A1)') &
          'Chain',ic,': ',trim(infiles(ic)),', '//colournames(colours(mod(ic-1,ncolours)+1))//'.', '  Lines/iter: ',real(n(ic)),'/',is(ic,n(ic)), &
          !'.  Burn-in: ',real(nburn(ic)),'/',isburn(ic),'.  Lmax:',maxval(dat(1,ic,:)),'.  Thin: file:',nint(is(ic,n(ic))/real(n(ic)*max(thin,1))),', tot:',totthin(ic),'.  Data pts: ',real(n(ic)-nburn(ic)),'.'
          '.  Burn-in: ',real(nburn(ic)),'/',isburn(ic),'.  Lmax:',loglmaxs(ic),'.  Thin: file:',nint(is(ic,n(ic))/real(n(ic)*max(thin,1))),', tot:',totthin(ic),'.  Data pts: ',real(n(ic)-nburn(ic)),'.'
  end do
  totiter = 0
  totpts  = 0
  contrchains = 0
  do ic=1,nchains0
     totiter = totiter + nint(is(ic,ntot(ic)))
     totpts = totpts + n(ic)-nburn(ic)
     if(n(ic).gt.nburn(ic)) contrchains = contrchains + 1
  end do
  totlines = sum(ntot(1:nchains0))
  !if(prprogress.ge.2.and.update.ne.1) write(*,'(A10,40x,3(A,I9),A1)')'Total:',' Lines:',sum(ntot(1:nchains0)),', iterations:',totiter,', data points: ',totpts,'.'
  !if(prprogress.ge.2.and.update.ne.1) write(*,'(A10,40x, A,ES7.1, A,ES7.1,A1, 68x, A,ES8.2,A1)')'Total:',' Lines: ',real(sum(ntot(1:nchains0))),', iter: ',real(totiter),',',' Data pts: ',real(totpts),'.'
  if(prprogress.ge.1.and.update.ne.1) write(*,'(4x,A, A,ES10.4, A,ES10.4, A,ES10.4,  A2,F5.1, A,I3,A1,I2,A1)') 'In all chains:','  number of lines: ',real(totlines), &
       ',  number of iterations: ',real(totiter),',  number of data points after burnin: ',real(totpts),' (',real(totpts)/real(totlines)*100,'%), contributing chains:',contrchains,'/',nchains0,'.'
  if(prruninfo.gt.0) write(*,*)
  
  
  
  if(prprogress.ge.2.and.update.eq.0) write(*,'(A)')'  Changing some variables...   '
  !NEW columns in dat: 1:logL 2:mc, 3:eta, 4:tc, 5:logdl, 6:spin, 7:kappa, 8: RA, 9:sindec,10:phase, 11:sinthJ0, 12:phiJ0, 13:alpha
  !Calculate the masses from Mch and eta:
  do ic=1,nchains0
     do i=1,ntot(ic)
        dvar = dsqrt(0.25d0-dat(3,ic,i))
        dvar1 = (0.5d0-dvar)/(0.5d0+dvar)
        dvar2 = dvar1**0.6d0
        dat(14,ic,i) = dat(2,ic,i) * ((1.d0+dvar1)**0.2d0 / dvar2)
        dat(15,ic,i) = dat(2,ic,i) * ((1.d0+1.d0/dvar1)**0.2d0 * dvar2)
     end do
  end do
  npar = 15
  if(par2.eq.0) par2 = npar
  
  acc = acc*0.25  !Transfom back to the actual acceptance rate
  
  
  
  !*** Put plot data in pldat, startval and jumps.  Print initial and starting values to screen.  Startval: 1: true value, 2: starting value, 3: Lmax value
  allocate(pldat(nchains,npar,narr))
  jumps = 0.
  offsetrun = 0
  if(prinitial.ne.0) then
     write(*,'(/,A)')'  True, starting and Lmax values for the chains:'
     write(*,'(11x,A20,14A10)')varnames(1:15)
  end if
  do ic=1,nchains
     pldat(ic,1:npar,1:n(ic)) = real(dat(1:npar,ic,1:n(ic))) !Note the change of order of indices!!!  Pldat has the same structure as alldat.  Pldat contains all data that's in dat (also before the burnin), but in single precision.  Use it to plot log(L), jumps, chains, etc., but not for statistics (and PDF creation).
     startval(ic,1:npar,1:2)  = real(dat(1:npar,ic,1:2)) !True value and starting value
     startval(ic,1:npar,3)    = real(dat(1:npar,icloglmax,iloglmax)) !Lmax value
     jumps(ic,1:npar,2:n(ic)) = real(dat(1:npar,ic,2:n(ic)) -  dat(1:npar,ic,1:n(ic)-1))
     if(prinitial.ne.0) then 
        if(ic.eq.1) write(*,'(5x,A11,F15.5,14F10.5,/)')'True:  ',startval(1,1:15,1)
        if(abs((sum(startval(ic,1:15,1))-sum(startval(ic,1:15,2)))/sum(startval(ic,1:15,1))).gt.1.e-10) then
           offsetrun = 1
           write(*,'(I4,A1,A11,F15.5,14F10.5)')ic,':','  Start: ',startval(ic,1:15,2)
           write(*,'(5x,A11,F15.5,14F10.5,/)')'Diff:  ',abs(startval(ic,1:15,1)-startval(ic,1:15,2))!/abs(startval(ic,1:15,1))
        end if
     end if
  end do
  if(prinitial.ne.0) then
     write(*,'(5x,A11,F15.5,14F10.5)')'Lmax:  ',startval(1,1:15,3)
     write(*,'(5x,A11,F15.5,14F10.5,/)')'Diff:  ',abs(startval(1,1:15,1)-startval(1,1:15,3))
  end if
  
  
  !if(prprogress.ge.1.and.update.eq.0) write(*,'(A)')'  Done.'
  if(prprogress.ge.2.and.update.eq.0) write(*,'(A,I12)')'  t0:',nint(t0)
  
  
  !Construct output file name
  ic = 1
  outputname = ' '
  !Number of detectors
  do ic=1,nchains0
     if(ic.eq.1) write(outputname,'(A,I1,A1)')trim(outputname),ndet(ic),'d'
     if(ic.gt.1 .and. (ndet(ic).ne.ndet(1) .or. sum(detnr(ic,1:ndet(ic))).ne.sum(detnr(1,1:ndet(1))))) write(outputname,'(A,I1,A1)')trim(outputname)//'-',ndet(ic),'d'
     if(ic.eq.1 .or. (ic.gt.1.and. (ndet(ic).ne.ndet(1) .or. sum(detnr(ic,1:ndet(ic))).ne.sum(detnr(1,1:ndet(1)))) )) then 
        do i=1,ndet(ic)
           write(outputname,'(A,I1)')trim(outputname),detnr(ic,i)
        end do
     end if
  end do
  !write(outputname,'(A,F4.2,A1,I3.3)')trim(outputname)//'_',startval(ic,6,1),'_',nint(acos(startval(ic,7,1))*r2d)
  !Spin magnitudes
  do ic=1,nchains0
     if(ic.eq.1) write(outputname,'(A,F4.2)')trim(outputname)//'_',startval(ic,6,1)
     if(ic.gt.1.and.startval(ic,6,1).ne.startval(1,6,1)) write(outputname,'(A,F4.2)')trim(outputname)//'-',startval(ic,6,1)
  end do
  !Theta_SLs
  do ic=1,nchains0
     if(ic.eq.1) write(outputname,'(A,I3.3)')trim(outputname)//'_',nint(acos(startval(ic,7,1))*r2d)
     if(ic.gt.1.and.startval(ic,7,1).ne.startval(1,7,1)) write(outputname,'(A,I3.3)')trim(outputname)//'-',nint(acos(startval(ic,7,1))*r2d)
  end do
  !print*,outputname
  
  
  
  !*** Put data in alldat
  if(mergechains.eq.1) then  !Merge chains, leave out burnin (then nchains = 1)
     allocate(alldat(1,npar1,nchains*narr))
     j = 1
     do ic=1,nchains
        do i=nburn(ic)+1,n(ic)
           alldat(1,1:npar,j) = real(dat(1:npar,ic,i))  !Note the change of order of indices!!!  Alldat has the same structure as pldat, but contains only info AFTER the burnin.
           j = j+1
        end do
     end do
     nchains = 1
     n(1) = j-1
     !if(prprogress.ge.1) write(*,'(A,I8,A,ES7.1,A)')'  Data points in combined chains: ',n(1),'  (',real(n(1)),')'
  else
     allocate(alldat(nchains,npar1,narr))
     do ic=1,nchains
        alldat(ic,1:npar,1:n(ic)-nburn(ic)) = real(dat(1:npar,ic,nburn(ic)+1:n(ic)))  !Note the change of order of indices!!!  Alldat has the same structure as pldat, but contains only info AFTER the burnin.
        n(ic) = n(ic)-nburn(ic)
     end do
     !if(prprogress.ge.1) write(*,'(A,I8)')' Datapoints in combined chains: ',sum(n(1:nchains))
  end if
  
  maxdots = 25000  !~Maximum number of dots to plot in e.g. chains plot, to prevent dots from being overplotted too much.  Use this to autoset chainpli
  !if(file.ge.2) maxdots = 10000 !Smaller max for eps,pdf to reduce file size (?)
  !if(maxval(n(1:nchains)).gt.maxdots) chainpli = nint(real(maxval(n(1:nchains)))/real(maxdots).)  !Change the number of points plotted in chains
  !if(maxval(n(1:nchains)).gt.maxdots .and. (file.ge.2.and.chainpli.eq.0)) then  !Change the number of points plotted in chains, for eps,pdf, to reduce file size
  if(chainpli.le.0) then
     !if(sum(n(1:nchains)).gt.maxdots) then  !Change the number of points plotted in chains,logL, etc. (For all output formats)
     if(sum(ntot(1:nchains0)).gt.maxdots) then  !Change the number of points plotted in chains,logL, etc. (For all output formats)
        !chainpli = nint(real(sum(n(1:nchains)))/real(maxdots))
        chainpli = nint(real(sum(ntot(1:nchains0)))/real(maxdots))  !Use ntot and nchains0, since n low if many points are in the burnin
        if(prprogress.ge.1.and.update.eq.0) write(*,'(A,I5,A,I5,A,I6,A)')'  Plotting every',chainpli,'-th state in likelihood, chains, jumps, etc. plots.  Average total thinning is',nint(avgtotthin),', for these plots it is',nint(avgtotthin*chainpli),'.'
     else
        chainpli = 1
        if(prprogress.ge.1.and.update.eq.0) write(*,'(A,I5,A)')'  Plotting *every* state in likelihood, chains, jumps, etc. plots.  Average total thinning remains',nint(avgtotthin),' for these plots.'
     end if
  end if
  
  
  
  
  ! **********************************************************************************************************************************
  ! ***  DO STATISTICS   *************************************************************************************************************
  ! **********************************************************************************************************************************
  
  
  !Sort all data and find the 90% interval limits for the wrapable parameters
  if(prprogress.ge.2) write(*,*)''
  shift = 0.
  wrap = 0
  do ic=1,nchains
     ival = ival0
     index = 0
     if(prprogress.ge.2.and.mergechains.eq.0) write(*,'(A,I2.2,A,$)')' Ch',ic,' '
     if(prprogress.ge.2.and.ic.eq.i.and.wrapdata.ge.1) write(*,'(A,$)')' Wrap data. '
     do p=par1,par2
        if(wrapdata.eq.0 .or. (p.ne.8.and.p.ne.10.and.p.ne.12.and.p.ne.13)) then
           call rindexx(n(ic),alldat(ic,p,1:n(ic)),index1(1:n(ic)))
           index(p,1:n(ic)) = index1(1:n(ic))
           cycle
        end if
        !if(p.ne.8.and.p.ne.10.and.p.ne.12.and.p.ne.13) cycle
        !NEW columns in dat: 1:logL 2:mc, 3:eta, 4:tc, 5:logdl, 6:spin, 7:kappa, 8: RA, 9:sindec,10:phase, 11:sinthJ0, 12:phiJ0, 13:alpha
        !Make sure data are between 0 and 2pi to start with:
        do i=1,n(ic)
           alldat(ic,p,i) = rev2pi(alldat(ic,p,i))
        end do
        call rindexx(n(ic),alldat(ic,p,1:n(ic)),index1(1:n(ic)))
        index(p,1:n(ic)) = index1(1:n(ic))
        minrange = 1.e30
        
        
        
        do i=1,n(ic)
           x1 = alldat(ic,p,index(p,i))
           x2 = alldat(ic,p,index(p,mod(i+nint(n(ic)*ival)-1,n(ic))+1))
           range = mod(x2 - x1 + real(20*pi),real(tpi))
           if(range.lt.minrange) then
              minrange = range
              y1 = x1
              y2 = x2
              !write(*,'(2I6,7F10.5)')i,mod(nint(i+n(ic)*ival),n(ic)),x1,x2,range,minrange,y1,y2,(y1+y2)/2.
           end if
           !write(*,'(2I6,7F10.5)')i,mod(nint(i+n(ic)*ival),n(ic)),x1,x2,range,minrange,y1,y2,(y1+y2)/2.
        end do !i
        centre = (y1+y2)/2.
        if(y1.gt.y2) then
           wrap(ic,p) = 1
           centre = mod(centre + pi, tpi) !Then distribution peaks close to 0/2pi, shift centre by pi
        end if
        !if(p.eq.8) write(*,'(3I6,7F10.5)')ic,p,wrap(ic,p),y1,y2,minrange,centre

        !See whether there's a gap in the data  WHY is necessary, should it work like this???
        if(wrap(ic,p).eq.0 .and. 1.eq.2) then
           !ymin = minval(alldat(ic,p,1:n(ic)))
           !ymax = maxval(alldat(ic,p,1:n(ic)))
           i0 = -1
           maxgap = -1.e30
           !write(*,'(2I3,I8)')ic,p,n(ic)
           do i=1,n(ic)-1
              x1 = alldat(ic,p,index(p,i))
              x2 = alldat(ic,p,index(p,i+1))
              !write(*,'(2I3,2I8,4F10.5)')ic,p,i,i0,x1,x2,x2-x2,maxgap
              if(x2-x1.gt.maxgap) then
                 maxgap = x2-x1
                 i0 = i
              end if
           end do !i
           x1 = alldat(ic,p,index(p,i0))
           x2 = alldat(ic,p,index(p,i0+1))
           !if(maxgap.gt.2*tpi/sqrt(real(n(ic)))) then
           if(maxgap.gt.0.1) then 
              x0 = (x1+x2)/2.
              !write(*,'(10F10.5)')x1,x2,(x1+x2)/2.,maxgap,ymin,ymax,centre,minrange,y1,y2
              !if(y1.lt.y2.and.(x0.lt.y1.or.x0.gt.y2)) wrap(ic,p) = 1  !If centre of max gap is outside 90% range  WHY???
              !if(y1.gt.y2.and.(x0.gt.y2.and.x0.lt.y1)) wrap(ic,p) = 1
           end if
        end if
        !if(p.eq.8) write(*,'(3I6,9F10.5)')ic,p,wrap(ic,p),y1,y2,x1/pi*12,x2/pi*12,x0/pi*12

        !Now, wrap around anticentre
        shift(ic,p) = 0.
        if(wrap(ic,p).eq.1) shift(ic,p) = tpi - mod(centre + pi, tpi)
        alldat(ic,p,1:n(ic)) = mod(alldat(ic,p,1:n(ic))+shift(ic,p),tpi)-shift(ic,p)
        pldat(ic,p,1:ntot(ic)) = mod(pldat(ic,p,1:ntot(ic))+shift(ic,p),tpi)-shift(ic,p) !Original data
        y1 = mod(y1+shift(ic,p),tpi)-shift(ic,p)
        y2 = mod(y2+shift(ic,p),tpi)-shift(ic,p)
        centre = mod(centre+shift(ic,p),tpi)-shift(ic,p)
        minrange = y2-y1
        !call rindexx(n(ic),alldat(ic,p,1:n(ic)),index(p,1:n(ic)))  !Re-sort
        call rindexx(n(ic),alldat(ic,p,1:n(ic)),index1(1:n(ic)))  !Re-sort
        index(p,1:n(ic)) = index1(1:n(ic))
        !if(p.eq.8) write(*,'(I3,A8,4x,6F10.5,I4)')ic,varnames(p),y1,y2,minrange,centre,minval(alldat(ic,p,1:n(ic))),maxval(alldat(ic,p,1:n(ic))),wrap(ic,p)
        if(abs(abs(minval(alldat(ic,p,1:n(ic)))-maxval(alldat(ic,p,1:n(ic))))-2*pi).lt.1.e-3) wrap(ic,p)=1 !If centre is around pi, still needs to be flagged 'wrap' to plot PDF
     end do !p







     !Do statistics
     !if(prprogress.ge.2) write(*,'(A)')' Calculating: statistics...'
     if(prprogress.ge.2.and.ic.eq.1) write(*,'(A,$)')'  Calc: stats, '
     do p=par1,par2
        !Determine the median
        if(mod(n(ic),2).eq.0) medians(p) = 0.5*(alldat(ic,p,index(p,n(ic)/2)) + alldat(ic,p,index(p,n(ic)/2+1)))
        if(mod(n(ic),2).eq.1) medians(p) = alldat(ic,p,index(p,(n(ic)+1)/2))
        
        !Mean:
        mean(p) = sum(alldat(ic,p,1:n(ic)))/real(n(ic))
        
        !Variances, etc:
        var1(p)=0.; var2(p)=0.; absvar1(p)=0.; absvar2(p)=0.; stdev1(p)=0.; stdev2(p)=0.
        do i=1,n(ic)
           var1(p) = var1(p) + (alldat(ic,p,i) - medians(p))**2
           var2(p) = var2(p) + (alldat(ic,p,i) - mean(p))**2
           absvar1(p) = absvar1(p) + abs(alldat(ic,p,i) - medians(p))
           absvar2(p) = absvar2(p) + abs(alldat(ic,p,i) - mean(p))
           stdev1(p) = stdev1(p) + (alldat(ic,p,i) - medians(p))*(alldat(ic,p,i) - medians(p))
           stdev2(p) = stdev2(p) + (alldat(ic,p,i) - mean(p))*(alldat(ic,p,i) - mean(p))
        end do
        
        absvar1(p) = absvar1(p)/real(n(ic))
        absvar2(p) = absvar2(p)/real(n(ic))
        stdev1(p)  = sqrt(stdev1(p)/real(n(ic)-1))
        stdev2(p)  = sqrt(stdev2(p)/real(n(ic)-1))
        
        !Save statistics:
        nstat = 6
        stats(ic,p,1) = medians(p)
        stats(ic,p,2) = mean(p)
        stats(ic,p,3) = absvar1(p)
        stats(ic,p,4) = absvar2(p)
        stats(ic,p,5) = stdev1(p)
        stats(ic,p,6) = stdev2(p)
     end do
     
     
     !Correlations:
     if(prcorr.gt.0.or.savestats.gt.0) then
        !write(*,'(A)')' Calculating correlations...   '
        write(*,'(A,$)')' corrs, '
        do p1=par1,par2
           !do p2=par1,par2
           do p2=p1,par2
              corrs(p1,p2) = 0.
              do i=1,n(ic)
                 !corrs(p1,p2) = corrs(p1,p2) + (alldat(ic,p1,i) - medians(p1))*(alldat(ic,p2,i) - medians(p2))
                 corrs(p1,p2) = corrs(p1,p2) + (alldat(ic,p1,i) - mean(p1))*(alldat(ic,p2,i) - mean(p2)) !Hardly differs from median method
              end do
              !corrs(p1,p2) = corrs(p1,p2) / (stdev1(p1)*stdev1(p2)*(n(ic)-1))
              corrs(p1,p2) = corrs(p1,p2) / (stdev2(p1)*stdev2(p2)*(n(ic)-1))
           end do !p2
        end do !p1
     end if
     
     
     !Autocorrelations:
     if(placorr.gt.0) then
        !write(*,'(A)')' Calculating autocorrelations...'
        write(*,'(A,$)')' autocorrs, '
        j1 = placorr/100 !Step size to get 100 autocorrelations per var
        do p=par1,par2
           acorrs(ic,p,:) = 0.
           !do j=1,ntot(ic)-1
           !do j=1,min(placorr,ntot(ic)-1)
           do j=0,min(100,ntot(ic)-1)
              do i=1,ntot(ic)-j*j1
                 acorrs(ic,p,j) = acorrs(ic,p,j) + (pldat(ic,p,i) - medians(p))*(pldat(ic,p,i+j*j1) - medians(p))
                 !acorrs(p,j) = acorrs(ic,p,j) + (pldat(ic,p,i) - mean(p))*(pldat(ic,p,i+j*j1) - mean(p))
              end do
              !if(j.eq.0) write(*,'(3I6,A,4F9.3)')j1,j,j*j1,'  '//varnames(p),acorrs(ic,0,j),acorrs(ic,p,j),(stdev1(p)*stdev1(p)*(ntot(ic)-j*j1)),acorrs(ic,p,0)
              acorrs(ic,0,j) = real(j*j1)
              acorrs(ic,p,j) = acorrs(ic,p,j) / (stdev1(p)*stdev1(p)*(ntot(ic)-j*j1))
              !acorrs(ic,p,j) = acorrs(ic,p,j) / (stdev2(p)*stdev2(p)*(ntot(ic)-j*j1))
              !if(j.eq.0) write(*,'(3I6,A,4F9.3)')j1,j,j*j1,'  '//varnames(p),acorrs(ic,0,j),acorrs(ic,p,j),(stdev1(p)*stdev1(p)*(ntot(ic)-j*j1)),acorrs(ic,p,0)
           end do !j
           !write(*,*)''
        end do !p
     end if
     
     
     !Determine interval ranges
     !if(prprogress.ge.2) write(*,'(A29,$)')' Determining interval levels: '
     if(prprogress.ge.2.and.ic.eq.1) write(*,'(A,$)')' prob.ivals: '
     c0 = 0
     do c=1,nival
        ival = ivals(c)
        if(abs(ival-ival0).lt.0.001) c0 = c
        if(c.ne.c0.and.prival.eq.0.and.savestats.eq.0) cycle
        
        if(prprogress.ge.2.and.ic.eq.1) write(*,'(F6.3,$)')ival
        do p=par1,par2
           minrange = 1.e30
           !write(*,'(A8,4x,4F10.5,I4)')varnames(p),y1,y2,minrange,centre,wrap(ic,p)
           do i=1,floor(n(ic)*(1.-ival))
              x1 = alldat(ic,p,index(p,i))
              x2 = alldat(ic,p,index(p,i+floor(n(ic)*ival)))
              range = abs(x2 - x1)
              !range = x2 - x1
              if(range.lt.minrange) then
                 minrange = range
                 y1 = x1
                 y2 = x2
              end if
              !write(*,'(I6,7F10.5)')i,x1,x2,range,minrange,y1,y2,(y1+y2)/2.
           end do
           centre = (y1+y2)/2.
           !write(*,'(A8,4x,4F10.5,I4)')varnames(p),y1,y2,minrange,centre,wrap(ic,p)

           !Save ranges:
           nr = 4
           ranges(ic,c,p,1) = y1
           ranges(ic,c,p,2) = y2
           ranges(ic,c,p,3) = (y1+y2)/2.
           ranges(ic,c,p,4) = y2-y1
           ranges(ic,c,p,5) = ranges(ic,c,p,4)
           if(p.eq.2.or.p.eq.3.or.p.eq.5.or.p.eq.6.or.p.eq.14.or.p.eq.15) ranges(ic,c,p,5) = ranges(ic,c,p,4)/ranges(ic,c,p,3)
        end do !p
     end do !c
     !if(prprogress.ge.2) write(*,'(A34,F8.4)')'.  Standard probability interval: ',ival0
     !if(prprogress.ge.2) write(*,'(A,F8.4,$)')', default ival:',ival0
     
     
     
     
     !Change variables
     !Columns in alldat(): 1:logL, 2:Mc, 3:eta, 4:tc, 5:logdl,   6:longi, 7:sinlati:, 8:phase, 9:spin,   10:kappa,     11:sinthJ0, 12:phiJ0, 13:alpha
     if(changevar.eq.1) then
        !if(prprogress.ge.2.and.update.eq.0) write(*,'(A,$)')' Changing some variables...   '
        if(prprogress.ge.2.and.ic.eq.i.and.update.eq.0) write(*,'(A,$)')'.  Change vars. '
        do p=par1,par2
           if(p.eq.5) then
              alldat(ic,p,1:n(ic)) = exp(alldat(ic,p,1:n(ic)))     !logD -> Distance
              if(ic.eq.1) startval(1:nchains0,p,1:3) = exp(startval(1:nchains0,p,1:3))
              stats(ic,p,1:nstat) = exp(stats(ic,p,1:nstat))
              ranges(ic,1:nival,p,1:nr) = exp(ranges(ic,1:nival,p,1:nr))
              !print*,ic,p
           end if
           if(p.eq.9.or.p.eq.11) then
              alldat(ic,p,1:n(ic)) = asin(alldat(ic,p,1:n(ic)))*r2d
              if(ic.eq.1) startval(1:nchains0,p,1:3) = asin(startval(1:nchains0,p,1:3))*r2d
              stats(ic,p,1:nstat) = asin(stats(ic,p,1:nstat))*r2d
              ranges(ic,1:nival,p,1:nr) = asin(ranges(ic,1:nival,p,1:nr))*r2d
           end if
           if(p.eq.7) then
              alldat(ic,p,1:n(ic)) = acos(alldat(ic,p,1:n(ic)))*r2d
              if(ic.eq.1) startval(1:nchains0,p,1:3) = acos(startval(1:nchains0,p,1:3))*r2d
              stats(ic,p,1:nstat) = acos(stats(ic,p,1:nstat))*r2d
              ranges(ic,1:nival,p,1:nr) = acos(ranges(ic,1:nival,p,1:nr))*r2d
              do c=1,nival
                 y1 = ranges(ic,c,p,2)
                 ranges(ic,c,p,2) = ranges(ic,c,p,1)  !acos is monotonously decreasing
                 ranges(ic,c,p,1) = y1
              end do
           end if
           if(p.eq.8) then
              alldat(ic,p,1:n(ic)) = alldat(ic,p,1:n(ic))*r2h
              if(ic.eq.1) startval(1:nchains0,p,1:3) = startval(1:nchains0,p,1:3)*r2h
              stats(ic,p,1:nstat) = stats(ic,p,1:nstat)*r2h
              ranges(ic,1:nival,p,1:nr) = ranges(ic,1:nival,p,1:nr)*r2h
           end if
           if(p.eq.10.or.p.eq.12.or.p.eq.13) then
              alldat(ic,p,1:n(ic)) = alldat(ic,p,1:n(ic))*r2d
              if(ic.eq.1) startval(1:nchains0,p,1:3) = startval(1:nchains0,p,1:3)*r2d
              stats(ic,p,1:nstat) = stats(ic,p,1:nstat)*r2d
              ranges(ic,1:nival,p,1:nr) = ranges(ic,1:nival,p,1:nr)*r2d
           end if
           ranges(ic,1:nival,p,3) = 0.5*(ranges(ic,1:nival,p,1) + ranges(ic,1:nival,p,2))
           ranges(ic,1:nival,p,4) = ranges(ic,1:nival,p,2) - ranges(ic,1:nival,p,1)
           ranges(ic,1:nival,p,5) = ranges(ic,1:nival,p,4)
           if(p.eq.2.or.p.eq.3.or.p.eq.5.or.p.eq.6.or.p.eq.14.or.p.eq.15) ranges(ic,1:nival,p,5) = ranges(ic,1:nival,p,5)/ranges(ic,1:nival,p,3)
        end do !p
        !NEW columns in dat: 1:logL 2:mc, 3:eta, 4:tc, 5:dl, 6:spin,  7:theta_SL, 8: RA,   9:dec, 10:phase, 11:thJ0, 12:phiJ0, 13:alpha
        varnames(1:15) = (/'logL','Mc','eta','tc','dl','spin','th_SL','RA','Dec','phase','thJo','phJo','alpha','M1','M2'/)
        pgvarns(1:15)  = (/'log Likelihood        ','M\dc\u (M\d\(2281)\u) ','\(2133)               ','t\dc\u (s)            ',  &
             'd\dL\u (Mpc)          ','a\dspin\u             ','\(2134)\dSL\u(\(2218))','R.A. (h)              ','Dec. (\(2218))        ', &
             '\(2147)\dc\u (\(2218))','\(2134)\dJ\u (\(2218))','\(2147)\dJ\u (\(2218))','\(2127)\dc\u (\(2218))','M\d1\u (M\d\(2281)\u) ','M\d2\u(M\d\(2281)\u)  '/)
        pgvarnss(1:15)  = (/'log L','M\dc\u','\(2133)','t\dc\u','d\dL\u','a\dspin\u','\(2134)\dSL\u','R.A.','Dec.','\(2147)\dc\u',  &
             '\(2134)\dJ0\u','\(2147)\dJ0\u','\(2127)\dc\u','M\d1\u','M\d2\u'/)
        !Include units
        pgvarnss(1:15)  = (/'log L','M\dc\u (M\d\(2281)\u)','\(2133)','t\dc\u (s)','d\dL\u (Mpc)','a\dspin\u','\(2134)\dSL\u (\(2218))','R.A. (h)','Dec. (\(2218))','\(2147)\dc\u (\(2218))',  &
             '\(2134)\dJ0\u (\(2218))','\(2147)\dJ0\u (\(2218))','\(2127)\dc\u (\(2218))','M\d1\u (M\d\(2281)\u)','M\d2\u (M\d\(2281)\u)'/)
        !Units only
        pgunits(1:15)  = (/'','M\d\(2281)\u ','','s','Mpc','','\(2218)','\uh\d','\(2218)','\(2218)','\(2218)','\(2218)','\(2218)','M\d\(2281)\u','M\d\(2281)\u'/)

        !if(prprogress.ge.2.and.update.eq.0) write(*,'(A)')'  Done.'
     end if !if(changevar.eq.1)
     
     
     !Find 100% range
     do p=par1,par2
        if(p.eq.1) cycle
        ranges(ic,nival+1,p,1) = minval(alldat(ic,p,1:n(ic)))
        ranges(ic,nival+1,p,2) = maxval(alldat(ic,p,1:n(ic)))
        ranges(ic,nival+1,p,3) = 0.5*(ranges(ic,nival+1,p,1) + ranges(ic,nival+1,p,2))
        ranges(ic,nival+1,p,4) = ranges(ic,nival+1,p,2) - ranges(ic,nival+1,p,1)
        ranges(ic,nival+1,p,5) = ranges(ic,nival+1,p,4)
        if(p.eq.2.or.p.eq.3.or.p.eq.5.or.p.eq.6.or.p.eq.14.or.p.eq.15) ranges(ic,nival+1,p,5) = ranges(ic,nival+1,p,5)/ranges(ic,nival+1,p,3)
        !write(*,'(I3,2F10.4)')p,ranges(ic,nival+1,p,1),ranges(ic,nival+1,p,2)
     end do
     
     if(prprogress.ge.2) then
        if(ic.eq.nchains) then
           write(*,*)
        else
           write(*,'(A,$)')'  '
        end if
     end if
     
     
     !**********************************************************************************************
     !******   PRINT STATISTICS   ******************************************************************
     !**********************************************************************************************
     
     !Print statistics to screen
     o=6
     if(prstat.gt.0) then
        write(o,'(/,A)')'  Main statistics:'
        c = c0
        write(o,'(A10, A12,2A10,A12, 4A8, 4A10,A8,A10, A4,A12,I3,A2)')'param.','model','median','mean','Lmax','stdev1','stdev2','abvar1','abvar2',  &
             'rng_c','rng1','rng2','drng','d/drng','delta','ok?','result (',nint(ivals(c0)*100),'%)'
        do p=par1,par2
           if(stdev1(p).lt.1.d-20) cycle !Parameter was probably not fitted
           write(o,'(A10,F12.6,2F10.4,F12.6, 4F8.4,4F10.4,F8.4,F10.4,$)')varnames(p),startval(ic,p,1),stats(ic,p,1),stats(ic,p,2),startval(ic,p,3),stdev1(p),stdev2(p),absvar1(p),  &
                absvar2(p),ranges(ic,c,p,3),ranges(ic,c,p,1),ranges(ic,c,p,2),ranges(ic,c,p,4),  &
                !abs(startval(ic,p,1)-stats(ic,p,1))/ranges(ic,c,p,4),ranges(ic,c,p,5)  !d/drange wrt median
                2*abs(startval(ic,p,1)-ranges(ic,c,p,3))/ranges(ic,c,p,4),ranges(ic,c,p,5)  !d/drange wrt centre of range
           if(startval(ic,p,1).gt.ranges(ic,c,p,1).and.startval(ic,p,1).lt.ranges(ic,c,p,2)) then
              write(o,'(A4,$)')'y '
           else
              write(o,'(A4,$)')'*N*'
           end if
           write(o,'(F10.4,A3,F9.4)')ranges(ic,c,p,3),'+-',0.5*ranges(ic,c,p,4)
        end do
     end if
     
     
     !Print correlations:
     if(prcorr.gt.0) then
        write(o,'(/,A)')'  Correlations:'
        write(o,'(A8,$)')''
        do p=par1,par2
           if(stdev1(p).gt.1.d-20) write(o,'(A7,$)')trim(varnames(p))
        end do
        write(o,*)''
        do p1=par1,par2
           if(stdev1(p1).lt.1.d-20) cycle
           write(o,'(A8,$)')trim(varnames(p1))
           do p2=par1,par2
              if(stdev1(p2).lt.1.d-20) cycle
              if(abs(corrs(p1,p2)).gt.0.5) then 
              !if(abs(corrs(p1,p2)).gt.-0.5) then 
                 write(o,'(F7.2,$)')corrs(p1,p2)
              else
                 write(o,'(A7,$)')''
              end if
           end do
           write(o,'(A)')'   '//trim(varnames(p1))
        end do
     end if
     
     
     !Print intervals:
     if(prival.gt.0) then
        write(o,'(/,A)')'  Probability intervals:'
        write(o,'(A20,A8,$)')'Interval:',''
        do c=1,nival+1
           write(o,'(F20.4,A9,$)')ivals(c),''
        end do
        write(o,*)''
        
        write(o,'(A8,2x,2A9,$)')'param.','model','median'
        do c=1,nival+1
           !write(o,'(2x,2A9,A8,$)')'rng1','rng2','in rnge'
           write(o,'(2x,3A9,$)')'centre','delta','in rnge'
        end do
        write(o,*)''
        do p=par1,par2
           if(stdev1(p).lt.1.d-20) cycle
           write(o,'(A8,2x,2F9.4,$)')varnames(p),startval(ic,p,1),stats(ic,p,1)
           do c=1,nival+1
              !write(o,'(2x,2F9.4,F6.3,$)')ranges(ic,c,p,1),ranges(ic,c,p,2),2*abs(startval(ic,p,1)-stats(ic,p,1))/ranges(ic,c,p,4) !Defined with median
              !write(o,'(2x,2F9.4,F6.3,$)')ranges(ic,c,p,1),ranges(ic,c,p,2),2*abs(startval(ic,p,1)-ranges(ic,c,p,3))/ranges(ic,c,p,4) !Defined with centre of prob. range
              write(o,'(2x,2F9.4,F6.3,$)')ranges(ic,c,p,3),ranges(ic,c,p,4),2*abs(startval(ic,p,1)-ranges(ic,c,p,3))/ranges(ic,c,p,4) !Defined with centre of prob. range
              if(startval(ic,p,1).gt.ranges(ic,c,p,1).and.startval(ic,p,1).lt.ranges(ic,c,p,2)) then
                 write(o,'(A3,$)')'y '
              else
                 write(o,'(A3,$)')'N*'
              end if
           end do
           write(o,*)''
        end do
     end if
     
     
  end do !ic
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  !Check convergence for multiple chains. this works only for fixed chain length, so take the min N
  if(nchains0.gt.1 .and. (prconv.ge.1.or.savestats.ge.1)) then
     chmean = 1.d-30
     totmean = 1.d-30
     nn = minval(ntot(1:nchains0))/2
     
     do p=par1,par2
        do ic=1,nchains0
           do i=nn+1,2*nn
              chmean(ic,p) = chmean(ic,p) + dat(p,ic,i) !Can't use pldat, because it may get wrapped above
              totmean(p) = totmean(p) + dat(p,ic,i)
           end do
        end do
     end do
     chmean = chmean/dble(nn)
     totmean = totmean/dble(nn*nchains0)
     
     chvar = 1.d-30
     chvar1 = 1.d-30
     totvar = 1.d-30
     do p=par1,par2
        do ic=1,nchains0
           do i=nn+1,2*nn
              dx = (dat(p,ic,i) - chmean(ic,p))**2 !Can't use pldat, because it may get wrapped above
              chvar(p) = chvar(p) + dx
              chvar1(ic,p) = chvar1(ic,p) + dx !Keep track of the variance per chain
           end do
           totvar(p) = totvar(p) + (chmean(ic,p) - totmean(p))**2
           chvar1(ic,p) = chvar1(ic,p)/dble(nn-1)
        end do
        chvar(p) = chvar(p)/dble(nchains0*(nn-1))
        totvar(p) = totvar(p)/dble(nchains0-1)
        
        !rhat(p) = ( dble(nn-1)/dble(nn) * chvar(p)  +  totvar(p) * (1.d0 + 1.d0/dble(nchains0)) ) / chvar(p)
        rhat(p) = min( dble(nn-1)/dble(nn)  +  totvar(p)/chvar(p) * (1.d0 + 1.d0/dble(nchains0)), 99.d0)
     end do
     
     if(prconv.ge.1) then
        write(*,*)''
        if(prconv.ge.2) write(*,'(A,I7,A)')'  Convergence parameters for',nn,' iterations:'
        write(*,'(18x,20A11)')varnames(par1:min(par2,13))
        if(prconv.ge.2) then
           write(*,'(A)')'  Means:'
           do ic=1,nchains0
              write(*,'(I16,A2,20F11.6)')ic,': ',chmean(ic,par1:min(par2,13))
           end do
           write(*,'(A18,20F11.6)')'           Total: ',totmean(par1:min(par2,13))
           
           write(*,*)''
           write(*,'(A)')'  Variances:'
        end if !if(prconv.ge.2)
     end if !if(prconv.ge.1)
     do ic=1,nchains0
        !write(*,'(I16,A2,20F11.6)')ic,': ',chvar1(ic,par1:min(par2,13))
        !if(chvar1(ic,2).lt.0.5*chvar(2).and.chvar1(ic,3).lt.0.5*chvar(3).and.chvar1(ic,2).lt.0.5*chvar(2).and.chvar1(ic,2).lt.0.5*chvar(2)) then
        lowvar = 0
        highvar = 0
        totrelvar = 1.d0
        ntotrelvar = 0
        do p=par1,min(par2,13)
           if(abs(chvar1(ic,p)).gt.1.e-30) then  !The parameters that were not fitted for have a variance of 0
              if(chvar1(ic,p).lt.0.5*chvar(p)) lowvar(p) = 1  !Too (?) low variance, mark it
              if(chvar1(ic,p).gt.2*chvar(p))  highvar(p) = 1  !Too (?) high variance, mark it
              totrelvar = totrelvar * chvar1(ic,p)/chvar(p) !Take geometric mean
              ntotrelvar = ntotrelvar + 1
           end if
        end do
        nlowvar = lowvar(2)+lowvar(3)+lowvar(6)+lowvar(7)  !Sum of 2 masses and 2 spin parameters
        nhighvar = highvar(2)+highvar(3)+highvar(6)+highvar(7)  !Sum of 2 masses and 2 spin parameters
        !totrelvar = totrelvar**(1.d0/dble(abs(min(par2,13)-par1+1))) !Take geometric mean of (the variance of each chain, relative to the total variance)
        totrelvar = totrelvar**(1.d0/dble(ntotrelvar)) !Take geometric mean of (the variance of each chain, relative to the total variance)
        if(prconv.ge.2) then
           ch = ' '
           if(nlowvar.eq.4) ch = '*'
           if(nhighvar.eq.4) ch = '#'
           write(*,'(I7,A3,$)')ic,': '//ch
           ch = ' '
           if(totrelvar.lt.0.5) ch = '*'
           if(totrelvar.gt.2.0) ch = '#'
           write(*,'(F8.3,A1,$)')totrelvar,ch
           do p=par1,min(par2,13)
              ch = ' '
              if(lowvar(p).eq.1) ch = '*'
              if(highvar(p).eq.1) ch = '#'
              write(*,'(F10.5,A1,$)')chvar1(ic,p),ch
           end do
           write(*,*)''
        end if !if(prconv.ge.2)
     end do
     if(prconv.ge.2) then
        write(*,'(A9,9x,20F11.5)')'  Total: ',chvar(par1:min(par2,13))
        
        write(*,*)''
        write(*,'(A)')'  Variances:'
        write(*,'(A18,20ES11.3)')'   Within chains: ',chvar(par1:min(par2,13))
        write(*,'(A18,20ES11.3)')'  Between chains: ',totvar(par1:min(par2,13))
     end if
     
     if(prconv.ge.1) write(*,'(A18,20F11.5)')'     Convergence: ',rhat(par1:min(par2,13)),sum(rhat(par1:min(par2,13)))/dble(min(par2,13)-par1+1)
     !write(*,*)''
  end if
  
  
  
  
  !Test: get mean and stdev for log(L)
  if(1.eq.2) then
     write(*,*)''
     
     nn = minval(ntot(1:nchains0))/2
     nlogl1 = nn+1
     nlogl2 = 2*nn
     nn = abs(nlogl2-nlogl1)
     
     write(*,'(A,I7,A)')'  Convergence criterion for',nn,' parameters:'
     write(*,'(16x,16x,20A11)')'Mean','Stddev','M-S','M+S','KS d','KS prob'
     
     do ic=1,nchains0
        nn1 = ntot(ic)/20
        ksn2 = 0
        ksd = 1.
        ksprob = 0.
        do nlogl1 = 1,ntot(ic),nn1
           nlogl2 = min(nlogl1+nn1,ntot(ic))
           nn = abs(nlogl2-nlogl1)+1
           if(nn.lt.nn1) cycle
           
           chmean = 1.d-30
           totmean = 1.d-30
           chvar = 1.d-30
           chvar1 = 1.d-30
           totvar = 1.d-30
           
           p=1
           do i=nlogl1,nlogl2
              chmean(ic,p) = chmean(ic,p) + dat(p,ic,i) !Can't use pldat, because it may get wrapped above
              totmean(p) = totmean(p) + dat(p,ic,i)
           end do
           chmean = chmean/dble(nn)
           totmean = totmean/dble(nn*nchains0)
           
           do i=nlogl1,nlogl2
              dx = (dat(p,ic,i) - chmean(ic,p))**2 !Can't use pldat, because it may get wrapped above
              chvar(p) = chvar(p) + dx
              chvar1(ic,p) = chvar1(ic,p) + dx !Keep track of the variance per chain
           end do
           totvar(p) = totvar(p) + (chmean(ic,p) - totmean(p))**2
           chvar1(ic,p) = chvar1(ic,p)/dble(nn-1)
           chvar(p) = chvar(p)/dble(nchains0*(nn-1))
           totvar(p) = totvar(p)/dble(nchains0-1)
           
           !write(*,'(I16,2I8,20F11.6)')ic,nlogl1,nlogl2,chmean(ic,1),chvar1(ic,1),chmean(ic,1)-chvar1(ic,1),chmean(ic,1)+chvar1(ic,1),ksd,ksprob
           !!print*,ic,p,nlogl1,nlogl2,nn
           ksdat1(1:nn) = dble(dat(p,ic,nlogl1:nlogl2))
           ksn1   = nn
           !!call kstwo(data1,n1,data2,n2,d,prob)
           if(ksn2.ne.0) call kstwo(ksdat1(1:ksn1),ksn1,ksdat2(1:ksn2),ksn2,ksd,ksprob)
           !
           !ksdat2 = ksdat1
           ksdat2(1:nn) = dble(dat(p,ic,nlogl1:nlogl2))
           ksn2 = ksn1
           
           write(*,'(I16,2I8,20F11.6)')ic,nlogl1,nlogl2,chmean(ic,1),chvar1(ic,1),chmean(ic,1)-chvar1(ic,1),chmean(ic,1)+chvar1(ic,1),ksd,ksprob,dlog10(ksprob+1.d-100)
        end do
        write(*,*)''
     end do
     write(*,*)''
     
     !KS test
     call pgbegin(1,'21/xs',1,1)
     call pgpap(scrsz,scrrat)
     call pgsch(1.5)
     call pgsvp(0.07,0.99,0.10,0.96)
     call pgswin(0.,real(maxval(ntot(1:nchains0))),-100.,0.)
     call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)
     call pgmtxt('B',2.2,0.5,0.5,'i')
     call pgmtxt('L',1.8,0.5,0.5,'log(d\dKS\u)')
     
     do ic=1,nchains0
        call pgsci(colours(mod(ic-1,ncolours)+1))
        nn = ntot(ic)/10
        ksd = 1.
        ksprob = 0.
        do nlogl1 = 1,ntot(ic),nn
           nlogl2 = nlogl1+nn-1
           if(nlogl2.gt.ntot(ic)) cycle
           
           ksn1   = nn
           ksdat1(1:ksn1) = dble(dat(p,ic,nlogl1:nlogl2))
           ksn2   = ntot(ic)-nlogl1+1
           ksdat2(1:ksn2) = dble(dat(p,ic,nlogl1:ntot(ic)))
           
           if(ksn2.ne.0) call kstwo(ksdat1(1:ksn1),ksn1,ksdat2(1:ksn2),ksn2,ksd,ksprob)
           !write(*,'(I16,2I8,20F11.6)')ic,nlogl1,nlogl2,ksd,ksprob,dlog10(ksprob+1.d-100)
           call pgpoint(1,real(nlogl1+nlogl2)/2.,real(dlog10(ksprob+1.d-100)),2)
        end do
        write(*,*)''
     end do
     
     call pgend
     write(*,*)''
     
  end if
  
  
  
  
  
  
  
  
  
  !Change the original chain data
  if(changevar.eq.1) then
     do ic=1,nchains0
        !NEW columns in dat: 1:logL 2:mc, 3:eta, 4:tc, 5:dl, 6:spin,  7:theta_SL, 8: RA,   9:dec, 10:phase, 11:thJ0, 12:phiJ0, 13:alpha
        !if(prprogress.ge.2.and.update.eq.0) write(*,'(A,$)')'Changing some variables...   '
        do p=par1,par2
           if(p.eq.5) pldat(ic,p,1:ntot(ic)) = exp(pldat(ic,p,1:ntot(ic)))
           if(p.eq.9.or.p.eq.11) pldat(ic,p,1:ntot(ic)) = asin(pldat(ic,p,1:ntot(ic)))*r2d
           if(p.eq.7) pldat(ic,p,1:ntot(ic)) = acos(pldat(ic,p,1:ntot(ic)))*r2d
           if(p.eq.8) pldat(ic,p,1:ntot(ic)) = pldat(ic,p,1:ntot(ic))*r2h
           if(p.eq.10.or.p.eq.12.or.p.eq.13) pldat(ic,p,1:ntot(ic)) = pldat(ic,p,1:ntot(ic))*r2d
        end do !p
     end do
     !if(prprogress.ge.2.and.update.eq.0) write(*,'(A)')'  Done.'
  end if !if(changevar.eq.1)
  
  
  
  
  deallocate(dat)
  
  
  
  ! **********************************************************************************************************************************
  ! ***  CREATE PLOTS   **************************************************************************************************************
  ! **********************************************************************************************************************************
  
  
  if(prprogress.ge.2) write(*,*)''
  if(combinechainplots.eq.1.and.(pllogl.eq.1.or.plchain.eq.1.or.plsigacc.ge.1)) then
     io = pgopen('chaininfo.eps'//trim(psclr))
     call pginitl(colour,file,whitebg)
  end if
    
  
  if(plot.eq.1.and.prprogress.ge.2.and.update.eq.0) write(*,'(A,$)')'  Plotting: '
  
  
  !***********************************************************************************************************************************      
  !Plot likelihood chain
  if(pllogl.eq.1) then
     !if(prprogress.ge.2.and.update.eq.0) write(*,'(A)')' Plotting chain likelihood...'
     if(prprogress.ge.2.and.update.eq.0) write(*,'(A,$)')' chain likelihood, '
     if(file.eq.0) then
        io = pgopen('12/xs')
        call pgsch(1.5)
     end if
     if(file.ge.1.and.combinechainplots.eq.0) then
        if(file.eq.1) io = pgopen('logL.ppm/ppm')
        if(file.ge.2) io = pgopen('logL.eps'//trim(psclr))
        call pgsch(1.2)
     end if
     if(io.le.0) then
        write(*,'(A,I4)')'Cannot open PGPlot device.  Quitting the programme',io
        goto 9999
     end if
     if(file.eq.0) call pgpap(scrsz,scrrat)
     if(file.eq.0) call pgsch(1.5)
     if(file.eq.1) call pgpap(bmpsz,bmprat)
     if(file.eq.1) call pgsch(1.5)
     !call pgscr(3,0.,0.5,0.)
     call pginitl(colour,file,whitebg)
     !call pgsubp(1,2)
     
     if(quality.eq.0) call pgsvp(0.08,0.95,0.06,0.94) !To make room for title
     
     ic = 1
     p=1
     !call pgpage
     xmax = -1.e30
     ymin =  1.e30
     ymax = -1.e30
     do ic=1,nchains0
        xmin = 0.
        xmax = max(xmax,maxval(is(ic,1:ntot(ic))))
        imin = 10                                              !Take into account burn-in
        if(scloglpl.eq.1) imin = nburn(ic)                   !Scale without taking into account burnin
        ymin = min(ymin,minval(pldat(ic,p,imin:ntot(ic)))) 
        ymax = max(ymax,maxval(pldat(ic,p,imin:ntot(ic))))
     end do
     ic = 1
     p = 1
     if(ymax.gt.0.) then !This is log(L)-log(Lo) (which we started saving later on), so that nullh=Lo=0
        !ymin = min(ymin,startval(ic,p,1),startval(ic,p,2),0.)
        !ymax = max(ymax,startval(ic,p,1),startval(ic,p,2),0.)
        if(scloglpl.eq.0) then                                 !Take into account 0, true and starting values
           ymin = min(ymin,startval(ic,p,1),startval(ic,p,2),0.)
           ymax = max(ymax,startval(ic,p,1),startval(ic,p,2),0.)
        else                                                     !Take into account true values only
           ymin = min(ymin,startval(ic,p,1))
           ymax = max(ymax,startval(ic,p,1))
        end if
     else !This is log(L), so that Lo = nullh
        ymin = min(ymin,startval(ic,p,1),startval(ic,p,2),nullh)
        ymax = max(ymax,startval(ic,p,1),startval(ic,p,2),nullh)
     end if
     dx = abs(xmax-xmin)*0.01
     dy = abs(ymax-ymin)*0.05
     
     call pgswin(xmin-dx,xmax+dx,ymin-dy,ymax+dy)
     call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)
     if(abs(startval(1,1,1)-startval(1,1,2))/abs(startval(1,1,1)).gt.1.e-10) then
        call pgsls(4)
        call pgbox('',0.0,0,'G',0.0,0)
        call pgsls(1)
     end if
     
     do ic=1,nchains0
        !call pgsci(defcolour)
        !if(nchains0.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
        !do i=ic,ntot(ic),chainpli !Start at ic to reduce overplotting
        !   call pgpoint(1,is(ic,i),pldat(ic,p,i),1)
        !end do
        
        !Give pre- and post-burnin different colour
        ci = defcolour
        if(nchains0.gt.1) ci = colours(mod(ic-1,ncolours)+1)
        call pgscidark(ci,file,whitebg)
        do i=ic,nburn(ic),chainpli !Start at ic to reduce overplotting
           call pgpoint(1,is(ic,i),pldat(ic,p,i),1)
        end do
        call pgsci(ci)
        do i=nburn(ic)+ic,ntot(ic),chainpli !Start at ic to reduce overplotting
           call pgpoint(1,is(ic,i),pldat(ic,p,i),1)
        end do
     end do
     
     !Plot max likelihood
     if(pllmax.ge.1) then
        ply = pldat(icloglmax,p,iloglmax)
        call pgsci(1)
        call pgpoint(1,is(icloglmax,iloglmax),ply,18)
        call pgsls(5)
        call pgline(2,(/-1.e20,1.e20/),(/ply,ply/))
     end if
     
     do ic=1,nchains0
        call pgsci(1)
        call pgsls(2)
        call pgline(2,(/-1.e20,1.e20/),(/startval(ic,p,1),startval(ic,p,1)/))
        call pgsci(6)
        if(nchains0.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
        !call pgline(2,(/real(nburn(ic)),real(nburn(ic))/),(/-1.e20,1.e20/))
        if(plburn.ge.1.and.isburn(ic).lt.is(ic,ntot(ic))) call pgline(2,(/isburn(ic),isburn(ic)/),(/-1.e20,1.e20/))
        call pgsci(1)
        call pgsls(4)
        if(nchains0.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
        call pgline(2,(/-1.e20,1.e20/),(/startval(ic,p,2),startval(ic,p,2)/))
        call pgsci(6)
        call pgline(2,(/-1.e20,1.e20/),real((/nullh,nullh/)))
     end do
     call pgsci(1)
     call pgsls(1)
     call pgmtxt('T',0.5,0.1,0.1,trim(pgvarns(p)))
     
     if(quality.eq.0) then
        !call pgsubp(1,1)
        !call pgsvp(0.,1.,0.,1.)
        !call pgswin(-1.,1.,-1.,1.)
     
        !call pgsch(sch*0.8)
        call pgmtxt('T',0.5,0.9,0.9,trim(outputname))  !Print title
        !call pgsch(sch)
     end if
     
     if(combinechainplots.eq.1) call pgpage
     if(combinechainplots.eq.0) then
        call pgend
        if(file.ge.2) then
           if(file.eq.3) then
              i = system('eps2pdf logL.eps  -o '//trim(outputdir)//'/'//trim(outputname)//'__logL.pdf   >& /dev/null')
              if(i.ne.0) write(*,'(A,I6)')'  Error converting plot',i
           end if
           i = system('mv -f logL.eps '//trim(outputdir)//'/'//trim(outputname)//'__logL.eps')
        end if
     end if
     if(file.eq.1) then
        i = system('convert -resize '//trim(bmpxpix)//' -depth 8 -unsharp '//trim(unsharplogl)//' logL.ppm  '//trim(outputdir)//'/'//trim(outputname)//'__logL.png')
        if(i.ne.0) write(*,'(A,I6)')'  Error converting plot',i
        i = system('rm -f logL.ppm')
        !if(i.ne.0) write(*,'(A)')'  Error removing file',i
     end if
  end if !if(pllogl.eq.1) then
  
  
  
  
  
  
  
  
  !***********************************************************************************************************************************      
  !Plot chains for each parameter
  if(plchain.eq.1) then
     !if(prprogress.ge.2.and.update.eq.0) write(*,'(A)')' Plotting parameter chains...'
     if(prprogress.ge.2.and.update.eq.0) write(*,'(A,$)')' parameter chains, '
     if(file.eq.0) then
        io = pgopen('13/xs')
        sch = 1.5
        lw = 1
     end if
     if(file.ge.1.and.combinechainplots.eq.0) then
        if(file.eq.1) io = pgopen('chains.ppm/ppm')
        if(file.ge.2) io = pgopen('chains.eps'//trim(psclr))
        lw = 3
        if(nplvar.ge.10) lw = 2
        if(quality.lt.2) lw = max(lw-1,1)  !Draft/Paper
        sch = 1.2
        if(nchains.eq.1.and.nplvar.gt.9) sch = 1.2
        if(quality.eq.0) then !Draft
           sch = sch*1.75
           lw = 2
        end if
        if(quality.eq.1) then !Paper
           if(nplvar.le.12) then
              sch = sch*1.75
              lw = 2
           else
              sch = sch*1.25
              lw = 1
           end if
        end if
        if(quality.eq.2) then !Talk
           if(nplvar.gt.12) then
              sch = sch*1.5
              lw = 1
           end if
           if(nplvar.le.12) then
              sch = sch*2
              lw = 2
           end if
           !if(nplvar.le.6) then
           !   !sch = sch*4
           !   !lw = 4
           !end if
        end if
        if(quality.eq.3) then !Poster
           if(nplvar.eq.12.and.file.ge.2) then
              sch = sch*2.7
              lw = 3
           else
              !sch = sch*1.25
              sch = sch*1.5
              lw = 2
           end if
        end if
        if(quality.eq.4) then !Vivien's thesis
           sch = sch*2.5
           lw = 2
        end if
     end if
     if(io.le.0) then
        write(*,'(A,I4)')'Cannot open PGPlot device.  Quitting the programme',io
        goto 9999
     end if
     if(file.eq.0) call pgpap(scrsz,scrrat)
     if(file.eq.1) call pgpap(bmpsz,bmprat)
     if(file.ge.2) call pgpap(pssz,psrat)
     if(file.ge.2) call pgscf(2)
     !if(file.eq.1) call pgsch(1.5)
     
     call pgsch(sch)
     call pgslw(lw)
     
     
     if(nplvar.eq.2) call pgsubp(2,1)
     if(nplvar.eq.2.and.quality.ge.2) call pgsubp(1,2)
     if(nplvar.eq.3) call pgsubp(3,1)
     if(nplvar.eq.4) call pgsubp(2,2)
     if(nplvar.eq.6) call pgsubp(3,2)
     if(nplvar.eq.8) call pgsubp(4,2)
     if(nplvar.eq.9) call pgsubp(3,3)
     if(nplvar.eq.10) call pgsubp(5,2)
     if(nplvar.eq.12) call pgsubp(4,3)
     if(nplvar.eq.12.and.quality.eq.3) call pgsubp(3,4)
     if(nplvar.eq.14.or.nplvar.eq.15) call pgsubp(5,3)
     if(nplvar.eq.16) call pgsubp(4,4)
     
     !call pgsubp(4,3)
     !call pgscr(3,0.,0.5,0.)
     
     ic = 1
     !do p=2,npar0
     do j=1,nplvar
        p = plvars(j)
        
        call pgpage
        
        
        if(j.eq.1) call pginitl(colour,file,whitebg)
        if(file.eq.0.and.scrrat.gt.1.35) call pgsvp(0.08,0.95,0.1,0.95)
        if(file.eq.1.and.bmprat.gt.1.35) call pgsvp(0.08,0.95,0.1,0.95)
        if(file.ge.2.and.psrat.gt.1.35) call pgsvp(0.08,0.95,0.1,0.95)
        if(quality.eq.0) call pgsvp(0.08,0.95,0.06,0.87) !To make room for title
        if(quality.eq.4) call pgsvp(0.13,0.95,0.1,0.95)
        
        xmin = 0.
        !xmax = real(maxval(ntot(1:nchains0)))
        xmax = -1.e30
        ymin =  1.e30
        ymax = -1.e30
        do ic=1,nchains0
           xmax = max(xmax,maxval(is(ic,1:ntot(ic))))
           imin = 1                                               !Take into account burn-in
           if(scchainspl.eq.1) imin = nburn(ic)                   !Scale without taking into account burnin
           ymin = min(ymin,minval(pldat(ic,p,imin:ntot(ic))))
           ymax = max(ymax,maxval(pldat(ic,p,imin:ntot(ic))))
        end do
        
        if(p.eq.8) then
           if(ymin.lt.0..or.ymax.gt.24.) then
              ymin = 0.
              ymax = 24.
           end if
        end if
        if(p.eq.10.or.p.eq.12.or.p.eq.13) then
           if(ymin.lt.0..or.ymax.gt.360.) then
              ymin = 0.
              ymax = 360.
           end if
        end if
        dx = abs(xmax-xmin)*0.01
        dy = abs(ymax-ymin)*0.05
        if(dx.eq.0) then
           xmin = 0.5*xmin
           xmax = 2*xmax
           if(xmin.eq.0.) then
              xmin = -1.
              xmax = 1.
           end if
        end if
        if(dy.eq.0) then
           ymin = 0.5*ymin
           ymax = 2*ymax
           if(ymin.eq.0.) then
              ymin = -1.
              ymax = 1.
           end if
        end if
        
        call pgswin(xmin-dx,xmax+dx,ymin-dy,ymax+dy)
        if(quality.eq.1) call pgswin(xmin-dx,xmax+dx,ymin-dy,ymax+dy*2)
        call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)
        
        !Plot the actual chain values
        call pgsch(1.)
        if(chainsymbol.ne.1) call pgsch(0.7)
        call pgslw(1)
        !write(*,'(15I4)'),nsymbols,symbols(1:nsymbols)
        do ic=1,nchains0
           !call pgsci(mod(ic*2,10))
           !symbol = ic+1
           symbol = chainsymbol
           if(chainsymbol.le.-10) symbol = symbols(mod(ic-1,nsymbols)+1)
           call pgsci(defcolour)
           if(nchains0.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
           if(chainsymbol.eq.0) then !Plot lines rather than symbols
              call pgline(ntot(ic),is(ic,1:ntot(ic)),pldat(ic,p,1:ntot(ic)))
           else
              !do i=ic,ntot(ic),chainpli !Start at ic to reduce overplotting
              !   ply = pldat(ic,p,i)
              !   if(p.eq.8) ply = rev24(ply)
              !   if(p.eq.10.or.p.eq.12.or.p.eq.13) ply = rev360(ply)
              !   !call pgpoint(1,is(ic,i),ply,1) !Plot small dots
              !   call pgpoint(1,is(ic,i),ply,symbol) !Plot symbols
              !end do
              
              
              !Give pre- and post-burnin different colour
              ci = defcolour
              if(nchains0.gt.1) ci = colours(mod(ic-1,ncolours)+1)
              call pgscidark(ci,file,whitebg)
              do i=ic,nburn(ic),chainpli !Start at ic to reduce overplotting
                 ply = pldat(ic,p,i)
                 if(p.eq.8) ply = rev24(ply)
                 if(p.eq.10.or.p.eq.12.or.p.eq.13) ply = rev360(ply)
                 call pgpoint(1,is(ic,i),ply,symbol)
              end do
              call pgsci(ci)
              do i=nburn(ic)+ic,ntot(ic),chainpli !Start at ic to reduce overplotting
                 ply = pldat(ic,p,i)
                 if(p.eq.8) ply = rev24(ply)
                 if(p.eq.10.or.p.eq.12.or.p.eq.13) ply = rev360(ply)
                 call pgpoint(1,is(ic,i),ply,symbol)
              end do
              
              
           end if
        end do
        call pgsch(sch)
        call pgslw(lw)
        
        !Plot max likelihood
        if(pllmax.ge.1) then
           ply = pldat(icloglmax,p,iloglmax)
           if(p.eq.8) ply = rev24(ply)
           if(p.eq.10.or.p.eq.12.or.p.eq.13) ply = rev360(ply)
           call pgsci(1)
           call pgpoint(1,is(icloglmax,iloglmax),ply,12)
           call pgsls(5)
           call pgline(2,(/-1.e20,1.e20/),(/ply,ply/))
        end if
        
        !Plot burn-in, true and starting values
        do ic=1,nchains0
           call pgsls(2)
           call pgsci(6)
           
           !Plot burn-in phase
           if(nchains0.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
           !if(plburn.ge.1) call pgline(2,(/isburn(ic),isburn(ic)/),(/-1.e20,1.e20/))
           if(plburn.ge.1.and.isburn(ic).lt.is(ic,ntot(ic))) call pgline(2,(/isburn(ic),isburn(ic)/),(/-1.e20,1.e20/))
           call pgsci(1)
           
           
           !Plot true values in chains
           if(pltrue.eq.1) then
              if(mergechains.ne.1.or.ic.eq.1) then !The units of the true values haven't changed (e.g. from rad to deg) for ic>1 (but they have for the starting values, why?)
              !if(ic.eq.1) then !The units of the true values haven't changed (e.g. from rad to deg) for ic>1 (but they have for the starting values, why?)
                 plx = startval(ic,p,1) !True value
                 plx = max(min(1.e30,startval(ic,p,1)),1.e-30)
                 if(p.eq.8) plx = rev24(plx)
                 if(p.eq.10.or.p.eq.12.or.p.eq.13) plx = rev360(plx)
                 call pgline(2,(/-1.e20,1.e20/),(/plx,plx/))
                 if(p.eq.8) then
                    call pgline(2,(/-1.e20,1.e20/),(/plx-24.,plx-24./))
                    call pgline(2,(/-1.e20,1.e20/),(/plx+24.,plx+24./))
                 end if
                 if(p.eq.10.or.p.eq.12.or.p.eq.13) then
                    call pgline(2,(/-1.e20,1.e20/),(/plx-360.,plx-360./))
                    call pgline(2,(/-1.e20,1.e20/),(/plx+360.,plx+360./))
                    !print*,p,plx,plx-360.,plx+360.
                 end if
              end if
           end if
           
           
           !Plot starting values in chains
           if(plstart.eq.1.and.abs((startval(ic,p,1)-startval(ic,p,2))/startval(ic,p,1)) .gt. 1.e-10) then
              call pgsls(4)
              if(nchains0.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
              plx = startval(ic,p,2) !Initial value
              if(p.eq.8) plx = rev24(plx)
              if(p.eq.10.or.p.eq.12.or.p.eq.13) plx = rev360(plx)
              call pgline(2,(/-1.e20,1.e20/),(/plx,plx/))
              if(p.eq.8) then
                 call pgline(2,(/-1.e20,1.e20/),(/plx-24.,plx-24./))
                 call pgline(2,(/-1.e20,1.e20/),(/plx+24.,plx+24./))
              end if
              if(p.eq.10.or.p.eq.12.or.p.eq.13) then
                 call pgline(2,(/-1.e20,1.e20/),(/plx-360.,plx-360./))
                 call pgline(2,(/-1.e20,1.e20/),(/plx+360.,plx+360./))
              end if
           end if
        end do !ic=1,nchains0
        
        call pgsci(1)
        call pgsls(1)
        write(string,'(F6.3)')rhat(p)
        !call pgmtxt('T',1.,0.,0.,'Chain: '//trim(pgvarns(p)))
        call pgmtxt('T',-1.,0.,0.,' '//trim(pgvarnss(p)))
        if(nchains0.gt.1.and.prconv.ge.1) call pgmtxt('T',1.,1.,1.,'Conv: '//trim(string))
     end do !do j=1,nplvar
     
     if(quality.eq.0) then
        call pgsubp(1,1)
        call pgsvp(0.,1.,0.,1.)
        call pgswin(-1.,1.,-1.,1.)
        
        call pgsch(sch*0.8)
        call pgmtxt('T',-0.7,0.5,0.5,trim(outputname))  !Print title
        call pgsch(sch)
     end if
     
     if(combinechainplots.eq.1) call pgpage
     if(combinechainplots.eq.0) then
        call pgend
        if(file.ge.2) then
           if(file.eq.3) then
              i = system('eps2pdf chains.eps  -o '//trim(outputdir)//'/'//trim(outputname)//'__chains.pdf   >& /dev/null')
              if(i.ne.0) write(*,'(A,I6)')'  Error converting plot',i
           end if
           i = system('mv -f chains.eps '//trim(outputdir)//'/'//trim(outputname)//'__chains.eps')
        end if
     end if
     if(file.eq.1) then
        i = system('convert -resize '//trim(bmpxpix)//' -depth 8 -unsharp '//trim(unsharpchain)//' chains.ppm  '//trim(outputdir)//'/'//trim(outputname)//'__chains.png')
        if(i.ne.0) write(*,'(A,I6)')'  Error converting plot',i
        i = system('rm -f chains.ppm')
     end if
  end if !if(plchain.eq.1)








  !***********************************************************************************************************************************      
  !Plot L vs parameter value
  if(plparl.eq.1) then
     !if(prprogress.ge.2.and.update.eq.0) write(*,'(A)')' Plotting parameter-L plot...'
     if(prprogress.ge.2.and.update.eq.0) write(*,'(A,$)')' parameter-L, '
     if(file.eq.0) then
        io = pgopen('22/xs')
        sch = 1.5
        lw = 1
     end if
     if(file.ge.1.and.combinechainplots.eq.0) then
        if(file.eq.1) io = pgopen('parlogl.ppm/ppm')
        if(file.ge.2) io = pgopen('parlogl.eps'//trim(psclr))
        lw = 3
        if(nplvar.ge.10) lw = 2
        if(quality.lt.2) lw = max(lw-1,1)  !Draft/Paper
        sch = 1.2
        if(nchains.eq.1.and.nplvar.gt.9) sch = 1.2
        if(quality.eq.0) then !Draft
           sch = sch*1.75
           lw = 2
        end if
        if(quality.eq.1) then !Paper
           if(nplvar.le.12) then
              sch = sch*1.75
              lw = 2
           else
              sch = sch*1.25
              lw = 1
           end if
        end if
        if(quality.eq.2) then !Talk
           if(nplvar.gt.12) then
              sch = sch*1.5
              lw = 1
           end if
           if(nplvar.le.12) then
              sch = sch*2
              lw = 2
           end if
           !if(nplvar.le.6) then
           !   !sch = sch*4
           !   !lw = 4
           !end if
        end if
        if(quality.eq.3) then !Poster
           if(nplvar.eq.12.and.file.ge.2) then
              sch = sch*2.7
              lw = 3
           else
              !sch = sch*1.25
              sch = sch*1.5
              lw = 2
           end if
        end if
     end if
     if(io.le.0) then
        write(*,'(A,I4)')'Cannot open PGPlot device.  Quitting the programme',io
        goto 9999
     end if
     if(file.eq.0) call pgpap(scrsz,scrrat)
     if(file.eq.1) call pgpap(bmpsz,bmprat)
     if(file.ge.2) call pgpap(pssz,psrat)
     if(file.ge.2) call pgscf(2)
     !if(file.eq.1) call pgsch(1.5)
     
     call pgsch(sch)
     call pgslw(lw)
     
     if(file.eq.0.and.scrrat.gt.1.35) call pgsvp(0.08,0.95,0.1,0.95)
     if(file.eq.1.and.bmprat.gt.1.35) call pgsvp(0.08,0.95,0.1,0.95)
     if(file.ge.2.and.psrat.gt.1.35) call pgsvp(0.08,0.95,0.1,0.95)
     if(quality.eq.0) call pgsvp(0.08,0.95,0.06,0.87) !To make room for title
     
     if(nplvar.eq.2) call pgsubp(2,1)
     if(nplvar.eq.2.and.quality.ge.2) call pgsubp(1,2)
     if(nplvar.eq.3) call pgsubp(3,1)
     if(nplvar.eq.4) call pgsubp(2,2)
     if(nplvar.eq.6) call pgsubp(3,2)
     if(nplvar.eq.8) call pgsubp(4,2)
     if(nplvar.eq.9) call pgsubp(3,3)
     if(nplvar.eq.10) call pgsubp(5,2)
     if(nplvar.eq.12) call pgsubp(4,3)
     if(nplvar.eq.12.and.quality.eq.3) call pgsubp(3,4)
     if(nplvar.eq.14.or.nplvar.eq.15) call pgsubp(5,3)
     if(nplvar.eq.16) call pgsubp(4,4)
     
     !call pgsubp(4,3)
     !call pgscr(3,0.,0.5,0.)
     
     ic = 1
     !do p=2,npar0
     do j=1,nplvar
        p = plvars(j)
        call pgpage
        if(j.eq.1) call pginitl(colour,file,whitebg)
        xmin = 1.e30
        xmax = -1.e30
        ymin =  1.e30
        ymax = -1.e30
        do ic=1,nchains0
           !xmin = min(xmin,minval(is(ic,nburn(ic):ntot(ic))))
           !xmax = max(xmax,maxval(is(ic,nburn(ic):ntot(ic))))
           xmin = min(xmin,minval(pldat(ic,p,nburn(ic):ntot(ic))))
           xmax = max(xmax,maxval(pldat(ic,p,nburn(ic):ntot(ic))))
           ymin = min(ymin,minval(pldat(ic,1,nburn(ic):ntot(ic))))
           ymax = max(ymax,maxval(pldat(ic,1,nburn(ic):ntot(ic))))
        end do
        
        if(p.eq.8) then
           !ymin = max(ymin,0.)
           !ymax = min(ymax,24.)
           !ymin = max(rev24(ymin),0.)
           !ymax = min(rev24(ymax),24.)
           if(ymin.lt.0..or.ymax.gt.24.) then
              ymin = 0.
              ymax = 24.
           end if
        end if
        if(p.eq.10.or.p.eq.12.or.p.eq.13) then
           !ymin = max(ymin,0.)
           !ymax = min(ymax,360.)
           !write(*,'(I5,2F10.5,$)')p,ymin,ymax
           !ymin = max(rev360(ymin),0.)
           !ymax = min(rev360(ymax),360.)
           !write(*,'(2F10.5)')ymin,ymax
           if(ymin.lt.0..or.ymax.gt.360.) then
              ymin = 0.
              ymax = 360.
           end if
        end if
        dx = abs(xmax-xmin)*0.1
        dy = abs(ymax-ymin)*0.1
        if(dx.eq.0) then
           xmin = 0.5*xmin
           xmax = 2*xmax
           if(xmin.eq.0.) then
              xmin = -1.
              xmax = 1.
           end if
        end if
        if(dy.eq.0) then
           ymin = 0.5*ymin
           ymax = 2*ymax
           if(ymin.eq.0.) then
              ymin = -1.
              ymax = 1.
           end if
        end if
        xmin = xmin - dx
        xmax = xmax + dx
        !ymin = ymin - dy
        !ymax = ymax + dy
        
        !Plot L iso log(L)
        ymax = exp(ymax - ymin)
        !ymin = 0.
        
        
        
        call pgswin(xmin,xmax,ymin,ymax)
        if(quality.eq.1) call pgswin(xmin,xmax,ymin,ymax+dy) !An extra dy
        call pgswin(xmin,xmax,0.,ymax*1.1) !Test
        call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)
        
        !Plot the actual chain values
        call pgsch(1.)
        if(chainsymbol.ne.1) call pgsch(0.7)
        call pgslw(1)
        !write(*,'(15I4)'),nsymbols,symbols(1:nsymbols)
        do ic=1,nchains0
           !call pgsci(mod(ic*2,10))
           !symbol = ic+1
           symbol = chainsymbol
           if(chainsymbol.le.-10) symbol = symbols(mod(ic-1,nsymbols)+1)
           call pgsci(defcolour)
           if(nchains0.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
           do i=nburn(ic),ntot(ic),chainpli
              plx = pldat(ic,p,i)
              ply = pldat(ic,1,i)
              if(p.eq.8) plx = rev24(plx)
              if(p.eq.10.or.p.eq.12.or.p.eq.13) plx = rev360(plx)
              !call pgpoint(1,is(ic,i),plx,1) !Plot small dots
              call pgpoint(1,plx,exp(ply-ymin),symbol) !Plot symbols
              !print*,i,plx,ply,exp(ply-ymin)
           end do
        end do
        call pgsch(sch)
        call pgslw(lw)
        
        !Plot max likelihood
        if(pllmax.ge.1) then
           plx = pldat(icloglmax,p,iloglmax)
           if(p.eq.8) plx = rev24(plx)
           if(p.eq.10.or.p.eq.12.or.p.eq.13) plx = rev360(plx)
           ply = exp(pldat(icloglmax,1,iloglmax)-ymin)
           call pgsci(1)
           call pgpoint(1,plx,ply,12)
           call pgsls(5)
           call pgline(2,(/plx,plx/),(/-1.e20,1.e20/))
        end if
        
        
        !Plot true values
        do ic=1,nchains0
           call pgsls(2)
           call pgsci(1)
           
           !Plot true values
           if(pltrue.eq.1) then
              if(mergechains.ne.1.or.ic.eq.1) then !The units of the true values haven't changed (e.g. from rad to deg) for ic>1 (but they have for the starting values, why?)
              !if(ic.eq.1) then !The units of the true values haven't changed (e.g. from rad to deg) for ic>1 (but they have for the starting values, why?)
                 plx = startval(ic,p,1) !True value
                 plx = max(min(1.e30,startval(ic,p,1)),1.e-30)
                 if(p.eq.8) plx = rev24(plx)
                 if(p.eq.10.or.p.eq.12.or.p.eq.13) plx = rev360(plx)
                 call pgline(2,(/plx,plx/),(/-1.e20,1.e20/))
                 if(p.eq.8) then
                    call pgline(2,(/plx-24.,plx-24./),(/-1.e20,1.e20/))
                    call pgline(2,(/plx+24.,plx+24./),(/-1.e20,1.e20/))
                 end if
                 if(p.eq.10.or.p.eq.12.or.p.eq.13) then
                    call pgline(2,(/plx-360.,plx-360./),(/-1.e20,1.e20/))
                    call pgline(2,(/plx+360.,plx+360./),(/-1.e20,1.e20/))
                 end if
              end if
           end if
        end do !ic=1,nchains0
        
        call pgsci(1)
        call pgsls(1)
     end do !do j=1,nplvar
     
     if(quality.eq.0) then
        call pgsubp(1,1)
        call pgsvp(0.,1.,0.,1.)
        call pgswin(-1.,1.,-1.,1.)
        
        call pgsch(sch*0.8)
        call pgmtxt('T',-0.7,0.5,0.5,trim(outputname))  !Print title
        call pgsch(sch)
     end if
     
     if(combinechainplots.eq.1) call pgpage
     if(combinechainplots.eq.0) then
        call pgend
        if(file.ge.2) then
           if(file.eq.3) then
              i = system('eps2pdf parlogl.eps  -o '//trim(outputdir)//'/'//trim(outputname)//'__parlogl.pdf   >& /dev/null')
              if(i.ne.0) write(*,'(A,I6)')'  Error converting plot',i
           end if
           i = system('mv -f parlogl.eps '//trim(outputdir)//'/'//trim(outputname)//'__parlogl.eps')
        end if
     end if
     if(file.eq.1) then
        i = system('convert -resize '//trim(bmpxpix)//' -depth 8 -unsharp '//trim(unsharpchain)//' parlogl.ppm  '//trim(outputdir)//'/'//trim(outputname)//'__parlogl.png')
        if(i.ne.0) write(*,'(A,I6)')'  Error converting plot',i
        i = system('rm -f parlogl.ppm')
     end if
  end if !if(plparl.eq.1)








  !***********************************************************************************************************************************            
  !Plot jump sizes
  if(pljump.ge.1) then
     !if(prprogress.ge.2.and.update.eq.0) write(*,'(A)')' Plotting jump sizes...'
     if(prprogress.ge.2.and.update.eq.0) write(*,'(A,$)')' jump sizes, '
     if(file.eq.0) then
        io = pgopen('18/xs')
        sch=1.5
     end if
     if(file.ge.1.and.combinechainplots.eq.0) then
        if(file.eq.1) io = pgopen('jumps.ppm/ppm')
        if(file.ge.2) io = pgopen('jumps.eps'//trim(psclr))
        sch=1.2
     end if
     if(io.le.0) then
        write(*,'(A,I4)')'Cannot open PGPlot device.  Quitting the programme',io
        goto 9999
     end if
     if(file.eq.0) call pgpap(scrsz,scrrat)
     if(file.eq.0) sch=1.5
     if(file.eq.1) call pgpap(bmpsz,bmprat)
     if(file.eq.1) sch=1.5
     if(quality.eq.0) then !Draft
        !sch = sch*1.75
        sch = sch*1.4
        lw = 2
        call pgsvp(0.1,0.95,0.06,0.87) !To make room for title
     end if
     
     call pgsch(sch)
     
     !call pgsubp(4,3)
     !call pgscr(3,0.,0.5,0.)
     !call pginitl(colour,file,whitebg)
     
     if(nplvar.eq.2) call pgsubp(2,1)
     if(nplvar.eq.3) call pgsubp(3,1)
     if(nplvar.eq.4) call pgsubp(2,2)
     if(nplvar.eq.6) call pgsubp(3,2)
     if(nplvar.eq.8) call pgsubp(4,2)
     if(nplvar.eq.9) call pgsubp(3,3)
     if(nplvar.eq.10) call pgsubp(5,2)
     if(nplvar.eq.12) call pgsubp(4,3)
     if(nplvar.eq.12.and.quality.eq.3) call pgsubp(3,4)
     if(nplvar.eq.14.or.nplvar.eq.15) call pgsubp(5,3)
     if(nplvar.eq.16) call pgsubp(4,4)
     
     ic = 1
     !do p=2,npar0
     do j=1,nplvar
        p = plvars(j)
        call pgpage
        if(j.eq.1) call pginitl(colour,file,whitebg)
        xmax = -1.e30
        ymin =  1.e30
        ymax = -1.e30
        do ic=1,nchains0
           xmin = 0.
           xmax = max(xmax,maxval(is(ic,1:ntot(ic))))
           dx = abs(xmax-xmin)*0.01
           if(pljump.eq.1) then
              ymin = min(ymin,minval(jumps(ic,p,2:ntot(ic))))
              ymax = max(ymax,maxval(jumps(ic,p,2:ntot(ic))))
           end if
           do i=10,ntot(ic)
              if(pljump.eq.2.and.jumps(ic,p,i).gt.1.e-20) then
                 ymin = min(ymin,log10(abs(jumps(ic,p,i))))
                 ymax = max(ymax,log10(abs(jumps(ic,p,i))))
              end if
              ymin = -6.
              ymax = 1.
           end do
           !print*,p-1,ymin,ymax,dy
           dy = abs(ymax-ymin)*0.05
           if(dy.lt.1.e-10) dy = ymin*0.1
        end do
        
        call pgswin(xmin-dx,xmax+dx,ymin-dy,ymax+dy)
        if(pljump.eq.1) call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0) !lin
        if(pljump.eq.2) call pgbox('BCNTS',0.0,0,'BCLNTS',0.0,0) !log
        
        do ic=1,nchains0
           !call pgsci(mod(ic*2,10))
           call pgsci(defcolour)
           if(nchains0.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
           !if(thin.le.1) then
           !   if(pljump.eq.1) then
           !      do i=1,ntot(ic),chainpli
           !         call pgpoint(1,is(ic,i),jumps(ic,p,i),1)
           !      end do
           !   else
           !      do i=1,ntot(ic),chainpli
           !         call pgpoint(1,is(ic,i),log10(abs(jumps(ic,p,i))+1.e-30),1)
           !      end do
           !   end if
           !else
           !   if(pljump.eq.1) then
           !      call pgpoint(ntot(ic),is(ic,1:ntot(ic)),jumps(ic,p,1:ntot(ic)),1)
           !   else
           !      call pgpoint(ntot(ic),is(ic,1:ntot(ic)),log10(abs(jumps(ic,p,1:ntot(ic)))+1.e-30),1)
           !   end if
           !end if
           if(pljump.eq.1) then
              !do i=1,ntot(ic),chainpli
              do i=ic,ntot(ic),chainpli !Start at ic to reduce overplotting
                 call pgpoint(1,is(ic,i),jumps(ic,p,i),1)
              end do
           else
              !do i=1,ntot(ic),chainpli
              do i=ic,ntot(ic),chainpli !Start at ic to reduce overplotting
                 call pgpoint(1,is(ic,i),log10(abs(jumps(ic,p,i))+1.e-30),1)
              end do
           end if
        end do

        call pgsls(2)
        call pgsci(6)
        do ic=1,nchains0
           !call pgline(2,(/real(nburn(ic)),real(nburn(ic))/),(/-1.e20,1.e20/))
           if(nchains0.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
           !if(plburn.ge.1) call pgline(2,(/isburn(ic),isburn(ic)/),(/-1.e20,1.e20/))
           if(plburn.ge.1.and.isburn(ic).lt.is(ic,ntot(ic))) call pgline(2,(/isburn(ic),isburn(ic)/),(/-1.e20,1.e20/))
        end do
        call pgsci(1)
        call pgsls(1)
        !call pgmtxt('T',1.,0.5,0.5,'Jumps: '//trim(pgorigvarns(p)))
        call pgmtxt('T',-1.2,0.05,0.0,trim(pgorigvarns(p)))
     end do
     
     if(quality.eq.0) then
        call pgsubp(1,1)
        call pgsvp(0.,1.,0.,1.)
        call pgswin(-1.,1.,-1.,1.)
     
        call pgsch(sch*0.8)
        call pgmtxt('T',-0.7,0.5,0.5,trim(outputname))  !Print title
        call pgsch(sch)
     end if
     
     if(combinechainplots.eq.1) call pgpage
     if(combinechainplots.eq.0) then
        call pgend
        if(file.ge.2) then
           if(file.eq.3) then
              i = system('eps2pdf jumps.eps  -o '//trim(outputdir)//'/'//trim(outputname)//'__jumps.pdf   >& /dev/null')
              if(i.ne.0) write(*,'(A,I6)')'  Error converting plot',i
           end if
           i = system('mv -f jumps.eps '//trim(outputdir)//'/'//trim(outputname)//'__jumps.eps')
        end if
     end if
     if(file.eq.1) then
        i = system('convert -resize '//trim(bmpxpix)//' -depth 8 -unsharp '//trim(unsharpchain)//' jumps.ppm  '//trim(outputdir)//'/'//trim(outputname)//'__jumps.png')
        if(i.ne.0) write(*,'(A,I6)')'  Error converting plot',i
        i = system('rm -f jumps.ppm')
     end if
  end if !if(pljump.ge.1)









  !***********************************************************************************************************************************            
  !Plot sigma values ('jump proposal width')
  if(plsigacc.ge.1) then
     !if(prprogress.ge.2.and.update.eq.0) write(*,'(A)')' Plotting sigma...'
     if(prprogress.ge.2.and.update.eq.0) write(*,'(A,$)')' sigma, '
     if(file.eq.0) then
        io = pgopen('16/xs')
        call pgsch(1.5)
     end if
     if(file.ge.1.and.combinechainplots.eq.0) then
        if(file.eq.1) io = pgopen('sigs.ppm/ppm')
        if(file.ge.2) io = pgopen('sigs.eps'//trim(psclr))
        call pgsch(1.2)
     end if
     if(io.le.0) then
        write(*,'(A,I4)')'Cannot open PGPlot device.  Quitting the programme',io
        goto 9999
     end if
     if(file.eq.0) call pgpap(scrsz,scrrat)
     if(file.eq.0) call pgsch(1.5)
     if(file.eq.1) call pgpap(bmpsz,bmprat)
     if(file.eq.1) call pgsch(1.5)

     if(quality.eq.0) call pgsvp(0.08,0.95,0.06,0.87) !To make room for title

     !call pgsubp(4,3)
     !call pgscr(3,0.,0.5,0.)
     !call pginitl(colour,file,whitebg)

     if(nplvar.eq.2) call pgsubp(2,1)
     if(nplvar.eq.3) call pgsubp(3,1)
     if(nplvar.eq.4) call pgsubp(2,2)
     if(nplvar.eq.6) call pgsubp(3,2)
     if(nplvar.eq.8) call pgsubp(4,2)
     if(nplvar.eq.9) call pgsubp(3,3)
     if(nplvar.eq.10) call pgsubp(5,2)
     if(nplvar.eq.12) call pgsubp(4,3)
     if(nplvar.eq.12.and.quality.eq.3) call pgsubp(3,4)
     if(nplvar.eq.14.or.nplvar.eq.15) call pgsubp(5,3)
     if(nplvar.eq.16) call pgsubp(4,4)
     
     ic = 1
     !do p=2,npar0
     do j=1,nplvar
        p = plvars(j)
        call pgpage
        if(j.eq.1) call pginitl(colour,file,whitebg)
        xmax = -1.e30
        ymin =  1.e30
        ymax = -1.e30
        do ic=1,nchains0
           xmin = 0.
           xmax = max(xmax,maxval(is(ic,1:ntot(ic))))
           dx = abs(xmax-xmin)*0.01
           if(plsigacc.eq.1) then
              ymin = min(ymin,minval(sig(p,ic,10:ntot(ic))))
              ymax = max(ymax,maxval(sig(p,ic,10:ntot(ic))))
           end if
           do i=10,ntot(ic)
              if(plsigacc.eq.2.and.sig(p,ic,i).gt.1.e-20) then
                 ymin = min(ymin,log10(sig(p,ic,i)))
                 ymax = max(ymax,log10(sig(p,ic,i)))
              end if
           end do
           !print*,p-1,ymin,ymax,dy
           dy = abs(ymax-ymin)*0.05
           if(dy.lt.1.e-10) dy = ymin*0.1
        end do

        call pgswin(xmin-dx,xmax+dx,ymin-dy,ymax+dy)
        if(plsigacc.eq.1) call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0) !lin
        if(plsigacc.eq.2) call pgbox('BCNTS',0.0,0,'BCLNTS',0.0,0) !log

        do ic=1,nchains0
           !call pgsci(mod(ic*2,10))
           call pgsci(defcolour)
           if(nchains0.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
           !if(thin.le.1) then
           !   if(plsigacc.eq.1) then
           !      do i=1,ntot(ic),chainpli
           !         call pgpoint(1,is(ic,i),sig(p,ic,i),1)
           !      end do
           !   else
           !      do i=1,ntot(ic),chainpli
           !         call pgpoint(1,is(ic,i),log10(sig(p,ic,i)+1.e-30),1)
           !      end do
           !   end if
           !else
           !   if(plsigacc.eq.1) then
           !      call pgpoint(ntot(ic),is(ic,1:ntot(ic)),sig(p,ic,1:ntot(ic)),1)
           !   else
           !      call pgpoint(ntot(ic),is(ic,1:ntot(ic)),log10(sig(p,ic,1:ntot(ic))+1.e-30),1)
           !   end if
           !end if
           if(plsigacc.eq.1) then
              !do i=1,ntot(ic),chainpli
              do i=ic,ntot(ic),chainpli !Start at ic to reduce overplotting
                 call pgpoint(1,is(ic,i),sig(p,ic,i),1)
              end do
           else
              !do i=1,ntot(ic),chainpli
              do i=ic,ntot(ic),chainpli !Start at ic to reduce overplotting
                 call pgpoint(1,is(ic,i),log10(sig(p,ic,i)+1.e-30),1)
              end do
           end if
        end do

        call pgsls(2)
        call pgsci(6)
        do ic=1,nchains0
           !call pgline(2,(/real(nburn(ic)),real(nburn(ic))/),(/-1.e20,1.e20/))
           if(nchains0.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
           !if(plburn.ge.1) call pgline(2,(/isburn(ic),isburn(ic)/),(/-1.e20,1.e20/))
           if(plburn.ge.1.and.isburn(ic).lt.is(ic,ntot(ic))) call pgline(2,(/isburn(ic),isburn(ic)/),(/-1.e20,1.e20/))
        end do
        call pgsci(1)
        call pgsls(1)
        call pgmtxt('T',1.,0.5,0.5,'Sigma: '//trim(pgvarns(p)))
        !call pgmtxt('T',1.,0.5,0.5,'log Sigma: '//trim(pgvarns(p)))
     end do
     
     if(quality.eq.0) then
        call pgsubp(1,1)
        call pgsvp(0.,1.,0.,1.)
        call pgswin(-1.,1.,-1.,1.)
     
        call pgsch(sch*0.8)
        call pgmtxt('T',-0.7,0.5,0.5,trim(outputname))  !Print title
        call pgsch(sch)
     end if
     
     if(combinechainplots.eq.1) call pgpage
     if(combinechainplots.eq.0) then
        call pgend
        if(file.ge.2) then
           if(file.eq.3) then
              i = system('eps2pdf sigs.eps  -o '//trim(outputdir)//'/'//trim(outputname)//'__sigs.pdf   >& /dev/null')
              if(i.ne.0) write(*,'(A,I6)')'  Error converting plot',i
           end if
           i = system('mv -f sigs.eps '//trim(outputdir)//'/'//trim(outputname)//'__sigs.eps')
        end if
     end if
     if(file.eq.1) then
        i = system('convert -resize '//trim(bmpxpix)//' -depth 8 -unsharp '//trim(unsharpchain)//' sigs.ppm  '//trim(outputdir)//'/'//trim(outputname)//'__sigs.png')
        if(i.ne.0) write(*,'(A,I6)')'  Error converting plot',i
        i = system('rm -f sigs.ppm')
     end if
  end if !if(plsigacc.ge.1)









  !***********************************************************************************************************************************      
  !Plot acceptance rates
  if(plsigacc.ge.1) then
     !if(prprogress.ge.2.and.update.eq.0) write(*,'(A)')' Plotting acceptance rates...'
     if(prprogress.ge.2.and.update.eq.0) write(*,'(A,$)')' acceptance rates, '
     if(file.eq.0) then
        io = pgopen('17/xs')
        call pgsch(1.5)
     end if
     if(file.ge.1.and.combinechainplots.eq.0) then
        if(file.eq.1) io = pgopen('accs.ppm/ppm')
        if(file.ge.2) io = pgopen('accs.eps'//trim(psclr))
        call pgsch(1.2)
     end if
     if(io.le.0) then
        write(*,'(A,I4)')'Cannot open PGPlot device.  Quitting the programme',io
        goto 9999
     end if
     if(file.eq.0) call pgpap(scrsz,scrrat)
     if(file.eq.0) call pgsch(1.5)
     if(file.eq.1) call pgpap(bmpsz,bmprat)
     if(file.eq.1) call pgsch(1.5)

     if(quality.eq.0) call pgsvp(0.08,0.95,0.06,0.87) !To make room for title

     !call pgsubp(4,3)
     !call pgscr(3,0.,0.5,0.)
     !call pginitl(colour,file,whitebg)

     if(nplvar.eq.2) call pgsubp(2,1)
     if(nplvar.eq.3) call pgsubp(3,1)
     if(nplvar.eq.4) call pgsubp(2,2)
     if(nplvar.eq.6) call pgsubp(3,2)
     if(nplvar.eq.8) call pgsubp(4,2)
     if(nplvar.eq.9) call pgsubp(3,3)
     if(nplvar.eq.10) call pgsubp(5,2)
     if(nplvar.eq.12) call pgsubp(4,3)
     if(nplvar.eq.12.and.quality.eq.3) call pgsubp(3,4)
     if(nplvar.eq.14.or.nplvar.eq.15) call pgsubp(5,3)
     if(nplvar.eq.16) call pgsubp(4,4)
     
     ic = 1
     !do p=2,npar0
     do j=1,nplvar
        p = plvars(j)
        call pgpage
        if(j.eq.1) call pginitl(colour,file,whitebg)
        xmax = -1.e30
        ymin =  1.e30
        ymax = -1.e30
        do ic=1,nchains0
           xmin = 0.
           !xmax = max(xmax,real(ntot(ic)))
           xmax = max(xmax,maxval(is(ic,1:ntot(ic))))
           dx = abs(xmax-xmin)*0.01
           do i=1,ntot(ic)
              if(acc(p,ic,i).gt.1.e-10 .and. acc(p,ic,i).lt.1.-1.e-10) then
                 n0 = i
                 exit
              end if
           end do
           n0 = n0+10
           ymin = min(ymin,minval(acc(p,ic,n0:ntot(ic))))
           ymax = max(ymax,maxval(acc(p,ic,n0:ntot(ic))))
           dy = abs(ymax-ymin)*0.05
        end do

        call pgsci(1)
        call pgsls(1)
        call pgswin(xmin-dx,xmax+dx,ymin-dy,ymax+dy)
        call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)


        call pgsci(3)
        call pgsls(2)
        call pgline(2,(/-1.e20,1.e20/),(/0.25,0.25/))
        call pgsci(6)
        do ic=1,nchains0
           !call pgline(2,(/real(nburn(ic)),real(nburn(ic))/),(/-1.e20,1.e20/))
           if(nchains0.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
           !if(plburn.ge.1) call pgline(2,(/isburn(ic),isburn(ic)/),(/-1.e20,1.e20/))
           if(plburn.ge.1.and.isburn(ic).lt.is(ic,ntot(ic))) call pgline(2,(/isburn(ic),isburn(ic)/),(/-1.e20,1.e20/))
        end do
        call pgsci(1)
        call pgsls(1)
        call pgmtxt('T',1.,0.5,0.5,'Acceptance: '//trim(pgvarns(p)))

        do ic=1,nchains0
           !call pgsci(mod(ic*2,10))
           call pgsci(defcolour)
           if(nchains0.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
           !if(thin.le.1) then
           !   do i=1,ntot(ic),chainpli
           !      call pgpoint(1,is(ic,i),acc(p,ic,i),1)
           !   end do
           !else
           !   call pgpoint(ntot(ic),is(ic,1:ntot(ic)),acc(p,ic,1:ntot(ic)),1)
           !end if
           !do i=1,ntot(ic),chainpli
           do i=ic,ntot(ic),chainpli !Start at ic to reduce overplotting
              call pgpoint(1,is(ic,i),acc(p,ic,i),1)
           end do
        end do
     end do
     
     if(quality.eq.0) then
        call pgsubp(1,1)
        call pgsvp(0.,1.,0.,1.)
        call pgswin(-1.,1.,-1.,1.)
     
        call pgsch(sch*0.8)
        call pgmtxt('T',-0.7,0.5,0.5,trim(outputname))  !Print title
        call pgsch(sch)
     end if
     
     !if(combinechainplots.eq.1) call pgpage
     if(combinechainplots.eq.0) then
        call pgend
        if(file.ge.2) then
           if(file.eq.3) then
              i = system('eps2pdf accs.eps  -o '//trim(outputdir)//'/'//trim(outputname)//'__accs.pdf   >& /dev/null')
              if(i.ne.0) write(*,'(A,I6)')'  Error converting plot',i
           end if
           i = system('mv -f accs.eps '//trim(outputdir)//'/'//trim(outputname)//'__accs.eps')
        end if
     end if
     if(file.eq.1) then
        i = system('convert -resize '//trim(bmpxpix)//' -depth 8 -unsharp '//trim(unsharpchain)//' accs.ppm  '//trim(outputdir)//'/'//trim(outputname)//'__accs.png')
        if(i.ne.0) write(*,'(A,I6)')'  Error converting plot',i
        i = system('rm -f accs.ppm')
     end if
  end if !if(plsigacc.ge.1)
  
  
  
  if(file.ge.1.and.combinechainplots.eq.1.and.(pllogl.eq.1.or.plchain.eq.1.or.plsigacc.ge.1)) then
     call pgend
     if(file.eq.3) i = system('eps2pdf chaininfo.eps -o '//trim(outputdir)//'/'//trim(outputname)//'__chaininfo.pdf  >& /dev/null')
     i = system('mv -f chaininfo.eps '//trim(outputdir)//'/'//trim(outputname)//'__chaininfo.eps')
  end if
  !***********************************************************************************************************************************        
  
  
  
  
  
  
  
  
  
  
  !***********************************************************************************************************************************      
  !Plot autocorrelations for each parameter
  if(placorr.gt.0) then
     !if(prprogress.ge.2.and.update.eq.0) write(*,'(A)')' Plotting autocorrelations...'
     if(prprogress.ge.2.and.update.eq.0) write(*,'(A,$)')' autocorrelations, '
     if(file.eq.0) then
        io = pgopen('19/xs')
        call pgsch(1.5)
     end if
     if(file.ge.1) then
        if(file.eq.1) io = pgopen('acorrs.ppm/ppm')
        if(file.ge.2) io = pgopen('acorrs.eps'//trim(psclr))
        call pgsch(1.2)
     end if
     if(io.le.0) then
        write(*,'(A,I4)')'Cannot open PGPlot device.  Quitting the programme',io
        goto 9999
     end if
     if(file.eq.0) call pgpap(scrsz,scrrat)
     if(file.eq.0) call pgsch(1.5)
     if(file.eq.1) call pgpap(bmpsz,bmprat)
     if(file.eq.1) call pgsch(1.5)
     
     if(quality.eq.0) call pgsvp(0.08,0.95,0.06,0.87) !To make room for title
     
     !call pgsubp(4,3)
     !call pgscr(3,0.,0.5,0.)
     !call pginitl(colour,file,whitebg)
     
     if(nplvar.eq.2) call pgsubp(2,1)
     if(nplvar.eq.3) call pgsubp(3,1)
     if(nplvar.eq.4) call pgsubp(2,2)
     if(nplvar.eq.6) call pgsubp(3,2)
     if(nplvar.eq.8) call pgsubp(4,2)
     if(nplvar.eq.9) call pgsubp(3,3)
     if(nplvar.eq.10) call pgsubp(5,2)
     if(nplvar.eq.12) call pgsubp(4,3)
     if(nplvar.eq.12.and.quality.eq.3) call pgsubp(3,4)
     if(nplvar.eq.14.or.nplvar.eq.15) call pgsubp(5,3)
     if(nplvar.eq.16) call pgsubp(4,4)
     
     ic = 1
     !do p=2,npar0
     do j=1,nplvar
        p = plvars(j)
        call pgpage
        if(j.eq.1) call pginitl(colour,file,whitebg)
        xmin = 0.
        xmin = minval(acorrs(1,0,0:100))
        xmax = maxval(acorrs(1,0,0:100))
        dx = abs(xmax-xmin)*0.01
        ymin =  1.e30
        ymax = -1.e30
        do ic=1,nchains0
           ymin = min(ymin,minval(acorrs(ic,p,0:100)))
           ymax = max(ymax,maxval(acorrs(ic,p,0:100)))
        end do
        dy = abs(ymax-ymin)*0.05
        !write(*,'(I3,5F10.2)')p,xmin,xmax,ymin,ymax,acorrs(1,0,100)
        
        call pgswin(xmin-dx,xmax+dx,ymin-dy,ymax+dy)
        call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)
        
        do ic=1,nchains0
           call pgsci(defcolour)
           if(nchains.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
           do i=1,100!placorr!ntot(ic),chainpli
              call pgpoint(1,acorrs(ic,0,i),acorrs(ic,p,i),1)
           end do
        end do
        
        call pgsci(1)
        call pgsls(2)
        call pgline(2,(/-1.e20,1.e20/),(/0.,0./))
        call pgsci(1)
        call pgsls(1)
        call pgmtxt('T',1.,0.5,0.5,'Autocorrelation: '//trim(pgvarns(p)))
     end do
     
     if(quality.eq.0) then
        call pgsubp(1,1)
        call pgsvp(0.,1.,0.,1.)
        call pgswin(-1.,1.,-1.,1.)
     
        call pgsch(sch*0.8)
        call pgmtxt('T',-0.7,0.5,0.5,trim(outputname))  !Print title
        call pgsch(sch)
     end if
     
     if(combinechainplots.eq.1) call pgpage
     if(combinechainplots.eq.0) then
        call pgend
        if(file.ge.2) then
           if(file.eq.3) then
              i = system('eps2pdf acorrs.eps  -o '//trim(outputdir)//'/'//trim(outputname)//'__acorrs.pdf   >& /dev/null')
              if(i.ne.0) write(*,'(A,I6)')'  Error converting plot',i
           end if
           i = system('mv -f acorrs.eps '//trim(outputdir)//'/'//trim(outputname)//'__acorrs.eps')
        end if
     end if
     if(file.eq.1) then
        i = system('convert -resize '//trim(bmpxpix)//' -depth 8 -unsharp '//trim(unsharpchain)//' acorrs.ppm  '//trim(outputdir)//'/'//trim(outputname)//'__acorrs.png')
        if(i.ne.0) write(*,'(A,I6)')'  Error converting plot',i
        i = system('rm -f acorrs.ppm')
     end if
  end if !if(placorrs.gt.0)
  !***********************************************************************************************************************************      







  
  
  
  
  
  
  
  
  
  !***********************************************************************************************************************************      
  !Plot pdfs (1d)
  if(plpdf1d.eq.1) then
     if(plot.eq.0.and.savepdf.eq.1) write(*,'(A,$)')' Saving 1D pdfs...   '
     allocate(xbin(nchs,nbin1d+1),ybin(nchs,nbin1d+1),xbin1(nbin1d+1),ybin1(nbin1d+1),ybin2(nbin1d+1),ysum(nbin1d+1),yconv(nbin1d+1),ycum(nbin1d+1))
     if(plot.eq.1) then
        !if(prprogress.ge.2.and.update.eq.0) write(*,'(A,$)')' Plotting 1D pdfs...   '
        if(prprogress.ge.2.and.update.eq.0) write(*,'(A,$)')' 1D pdfs, '
        if(file.eq.0) then
           io = pgopen('14/xs')
           sch = 1.5
           lw = 1
        end if
        if(file.ge.1) then
           if(file.eq.1) io = pgopen('pdfs.ppm/ppm')
           if(file.ge.2) io = pgopen('pdfs.eps'//trim(psclr))
           lw = 3
           if(nplvar.ge.10) lw = 2
           if(quality.lt.2) lw = max(lw-1,1)  !Draft/Paper
           sch = 1.2!2.5
           if(nchains.eq.1.and.nplvar.gt.9) sch = 1.2
           if(quality.eq.0) then !Draft
              sch = sch*1.75
              lw = 2
           end if
           if(quality.eq.1) then !Paper
              if(nplvar.eq.12) then
                 sch = sch*1.75
                 lw = 2
              else
                 sch = sch*1.25
                 lw = 1
              end if
           end if
           if(quality.eq.2) then !Talk
              if(nplvar.le.12) then
                 sch = sch*2
                 lw = 2
              else
                 sch = sch*1.5
                 lw = 1
              end if
           end if
           if(quality.eq.3) then !Poster
              if(nplvar.eq.12.and.file.ge.2) then
                 sch = sch*2.7
                 lw = 3
              else
                 !sch = sch*1.25
                 !lw = 1
                 sch = sch*1.5
                 lw = 2
              end if
           end if
           if(quality.eq.4) then !Vivien's thesis
              sch = sch*2.5
              lw = 2
           end if
        end if !if(file.ge.1)
        if(io.le.0) then
           write(*,'(A,I4)')'Cannot open PGPlot device.  Quitting the programme',io
           goto 9999
        end if
        if(file.eq.0) call pgpap(scrsz,scrrat)
        if(file.eq.1) call pgpap(bmpsz,bmprat)
        if(file.ge.2) call pgpap(pssz,psrat)
        if(file.ge.2.and.quality.eq.3.and.nplvar.eq.12) call pgpap(10.6,0.925)
        if(file.ge.2) call pgscf(2)
        !call pgscr(3,0.,0.5,0.)
        !call pginitl(colour,file,whitebg)
        call pgslw(lw)
        call pgsch(sch)
        call pgsfs(fillpdf)
        
        
        !if(quality.eq.0) call pgsvp(0.08,0.95,0.06,0.87) !To make room for title
        
        if(panels(1)*panels(2).lt.1) then
           if(nplvar.eq.2) call pgsubp(2,1)
           if(nplvar.eq.3) call pgsubp(3,1)
           if(nplvar.eq.4) call pgsubp(2,2)
           if(nplvar.eq.6) call pgsubp(3,2)
           if(nplvar.eq.8) call pgsubp(4,2)
           if(nplvar.eq.9) call pgsubp(3,3)
           if(nplvar.eq.10) call pgsubp(5,2)
           if(nplvar.eq.12) call pgsubp(4,3)
           if(nplvar.eq.12.and.quality.eq.3) call pgsubp(3,4)
           if(nplvar.eq.14.or.nplvar.eq.15) call pgsubp(5,3)
           if(nplvar.eq.16) call pgsubp(4,4)
        else
           call pgsubp(panels(1),panels(1))
        end if
     end if !if(plot.eq.1)
     
     !Save 1D PDF data
     if(savepdf.eq.1) then
        open(unit=30,action='write',form='formatted',status='replace',file=trim(outputdir)//'/'//trim(outputname)//'__pdf1d.dat')
        write(30,'(3I6,T100,A)')nplvar,nchains,nbin1d,'Total number of plot variables, total number of chains, number of bins'
     end if
     
     !do p=par1,par2
     do j=1,nplvar
        p = plvars(j)
        if(plot.eq.1) then
           call pgpage
           if(j.eq.1) call pginitl(colour,file,whitebg)
        end if

        !Set x-ranges, bin the data and get y-ranges
        xmin = 1.e30
        xmax = -1.e30
        do ic=1,nchains
           xmin = min(xmin,minval(alldat(ic,p,1:n(ic))))
           xmax = max(xmax,maxval(alldat(ic,p,1:n(ic))))
        end do
        dx = xmax - xmin
        !dx = max(xmax - xmin,1.e-30)
        
        do ic=1,nchains
           x(ic,1:n(ic)) = alldat(ic,p,1:n(ic))
           xmin1 = minval(alldat(ic,p,1:n(ic)))
           xmax1 = maxval(alldat(ic,p,1:n(ic)))
           
           call bindata(n(ic),x(ic,1:n(ic)),1,nbin1d,xmin1,xmax1,xbin1,ybin1)
           
           !Weigh with likelihood.  I should probably do something like this at the start, to get updated ranges etc.
           !y(ic,1:n(ic)) = alldat(ic,1,1:n(ic))
           !call bindataa(n(ic),x(ic,1:n(ic)),y(ic,1:n(ic)),1,nbin1d,xmin1,xmax1,xbin1,ybin1) !Measure the amount of likelihood in each bin
           
           !NEW columns in dat: 1:logL 2:mc, 3:eta, 4:tc, 5:d_l, 6:spin, 7:theta_SL, 8: RA, 9:dec,10:phase, 11:thetaJ0, 12:phiJ0, 13:alpha, 14:M1, 15:M2
           if(p.eq.5.or.p.eq.7.or.p.eq.9.or.p.eq.11) then  !Do something about chains 'sticking to the wall'
              if(ybin1(1).gt.ybin1(2)) ybin1(1)=0.
              if(ybin1(nbin1d).gt.ybin1(nbin1d-1)) ybin1(nbin1d)=0.
           end if
           
           !Normalise 1D PDF
           if(normpdf1d.gt.0) then
              if(normpdf1d.eq.1) then !Normalise the SURFACE, not the height (because of different bin size).  This is the default
                 norm = 0.
                 do i=1,nbin1d+1
                    norm = norm + ybin1(i)
                 end do
                 norm = norm*(xmax1-xmin1)
                 ybin1 = ybin1/norm
              else !Normalise to the height of the PDF
                 if(normpdf1d.eq.2) ybin1 = ybin1/maxval(ybin1)  !Normalise to the height of the PDF
                 if(normpdf1d.eq.3) ybin1 = ybin1/(maxval(ybin1)**0.5)  !Normalise to the sqrt of the height of the PDF; Works nicely for comparing parallel-tempering chains
                 if(ic*j.eq.1)write(*,'(//,A,/)')'  *** WARNING: using non-default normalisation for PDFs ***'
              end if
           end if
           
           !Smoothen 1D PDF
           ybin2 = ybin1
           if(smooth.gt.1) then
              !i0 = nbin1d/10
              !print*,nbin1d/10,nint(min(max(real(nbin1d)/real(smooth),1.0),real(nbin1d)/2.))
              !i0 = nint(min(max(real(nbin1d)/real(smooth),1.0),real(nbin1d)/2.))
              i0 = min(max(smooth,1),floor(real(nbin1d)/2.))
              do i=1+i0,nbin1d+1-i0
                 coefs1(1:2*i0+1) = ybin1(i-i0:i+i0)
                 call savgol(coefs1(1:2*i0+1),2*i0+1,i0,i0,0,4)
                 do i1=1,i0+1
                    coefs(i0-i1+2) = coefs1(i1)
                 end do
                 do i1 = i0+2,2*i0+1
                    coefs(3*i0+3-i1) = coefs1(i1)
                 end do
                 ybin2(i) = 0.
                 do i1=1,2*i0+1
                    ybin2(i) = ybin2(i) + coefs(i1) * ybin1(i+i1-i0-1)
                 end do
              end do
              ybin1 = ybin2
           end if !if(smooth.gt.1)
           xbin(ic,1:nbin1d+1) = xbin1(1:nbin1d+1)
           ybin(ic,1:nbin1d+1) = ybin1(1:nbin1d+1)
           
           !Save binned data
           if(savepdf.eq.1) then
              write(30,'(3I6,T100,A)')ic,p,wrap(ic,p),'Chain number, variable number, and wrap'
              write(30,'(2ES15.7,T100,A)')startval(ic,p,1:2),'True and starting value'
              write(30,'(6ES15.7,T100,A)')stats(ic,p,1:6),'Stats: median, mean, absvar1, absvar2, stdev1, stdev2'
              write(30,'(5ES15.7,T100,A)')ranges(ic,c0,p,1:5),'Ranges: lower,upper limit, centre, width, relative width'
              write(30,'(2ES15.7,T100,A)')xmin1,xmax1,'Xmin and Xmax of PDF'
              do i=1,nbin1d+1
                 write(30,'(2ES15.7)')xbin1(i),ybin1(i)
              end do
           end if
        end do !ic
        
        if(plot.eq.1) then
           xmin = xmin - 0.1*dx
           xmax = xmax + 0.1*dx
           ymin = 0.
           ymax = 1.e-20
           !print*,xmin,xmax,ymin,ymax
           do ic=1,nchains
              !ymax = max(ymax,maxval(ybin(ic,1:nbin1d+1)))
              do i=1,nbin1d+1
                 !print*,ybin(ic,i),ymax
                 if(ybin(ic,i).gt.ymax) then
                    ymax = ybin(ic,i)
                    xpeak = xbin(ic,i)
                 end if
              end do
           end do
           ymax = ymax*1.1
           if(dx.eq.0) then
              xmin = 0.5*xmin
              xmax = 2*xmax
              if(xmin.eq.0.) then
                 xmin = -1.
                 xmax = 1.
              end if
           end if
           !print*,xmin,xmax,ymin,ymax
           if(ymax.lt.1.e-19) ymax = 1.
           
           
           if(file.eq.0.and.scrrat.gt.1.35) call pgsvp(0.08,0.95,0.1,0.95)
           if(file.eq.1.and.bmprat.gt.1.35) call pgsvp(0.08,0.95,0.1,0.95)
           if(file.ge.2.and.psrat.gt.1.35) call pgsvp(0.08,0.95,0.1,0.95)
           if(quality.eq.0) call pgsvp(0.08,0.95,0.06,0.87) !To make room for title
           if(quality.eq.4) call pgsvp(0.13,0.95,0.1,0.95)
           
           call pgsch(sch)
           call pgswin(xmin,xmax,ymin,ymax)
           if(abs(dx).lt.1.e-30) then !So that the program doesn't hang if a parameter is kept constant
              xbin = 0.
              ybin = 0.
           end if
           
           !Plot 1D PDF
           call pgsci(1)
           if(file.ge.2) call pgslw(lw)
           do ic=1,nchains
              if(fillpdf.ge.3) call pgshs(45.0*(-1)**ic,2.0,real(ic)/real(nchains0)) !Set hatch style: angle = +-45deg, phase between 0 and 1 (1/nchains0, 2/nchains0, ...)
              if(nchains.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
              xbin1(1:nbin1d+1) = xbin(ic,1:nbin1d+1)
              ybin1(1:nbin1d+1) = ybin(ic,1:nbin1d+1)
              if(wrap(ic,p).eq.0) then
                 if(nchains.eq.1) call pgsci(15)
                 call pgpoly(nbin1d+2,(/xbin1(1),xbin1(1:nbin1d+1)/),(/0.,ybin1(1:nbin1d+1)/))
                 !Plot pdf contour
                 !if(nchains.eq.1) call pgsci(1)
                 !call pgsci(1)
                 if(fillpdf.eq.1) call pgsci(1)
                 if(nchains.eq.1) call pgsci(2)
                 call pgline(nbin1d+1,xbin1(1:nbin1d+1),ybin1(1:nbin1d+1)) !:nbin1d) ?
                 
                 !Fix the loose ends
                 call pgline(2,(/xbin1(1),xbin1(1)/),(/0.,ybin1(1)/))
                 call pgline(2,(/xbin1(nbin1d+1),xbin1(nbin1d+1)/),(/ybin1(nbin1d+1),0./))
              else !If parameter is wrapped
                 plshift = real(2*pi)
                 if(changevar.eq.1) plshift = 360.
                 if(changevar.eq.1.and.p.eq.8) plshift = 24.  !RA in hours
                 if(nchains.eq.1) call pgsci(15)
                 call pgpoly(nbin1d+3,(/xbin1(1),xbin1(1:nbin1d),xbin1(1)+plshift,xbin1(1)+plshift/),(/0.,ybin1(1:nbin1d),ybin1(1),0./))
                 
                 !Plot pdf contour
                 !call pgsci(1)
                 !if(fillpdf.ne.1) call pgsci(colours(mod(ic-1,ncolours)+1))
                 if(fillpdf.eq.1) call pgsci(1)
                 if(nchains.eq.1) call pgsci(2)
                 call pgline(nbin1d,xbin1(1:nbin1d),ybin1(1:nbin1d))
                 
                 !Plot dotted lines outside the pdf for wrapped periodic variables
                 call pgsls(4)
                 !if(file.ge.2) call pgslw(2)
                 call pgline(nbin1d+1,(/xbin1(1:nbin1d)-plshift,xbin1(1)/),(/ybin1(1:nbin1d),ybin1(1)/))
                 call pgline(nbin1d,xbin1+plshift,ybin1)
                 
                 !Fix the loose end
                 call pgsls(1)
                 if(file.ge.2) call pgslw(lw)
                 call pgline(2,(/xbin1(nbin1d),xbin1(1)+plshift/),(/ybin1(nbin1d),ybin1(1)/))
              end if
           end do !ic
           !Plot lines again over surface of overlapping distributions
           if(nchains.gt.1.and.fillpdf.eq.1) then
              call pgsls(4)
              do ic=1,nchains
                 call pgsci(1)
                 !call pgsci(colours(mod(ic-1,ncolours)+1))
                 xbin1(1:nbin1d+1) = xbin(ic,1:nbin1d+1)
                 ybin1(1:nbin1d+1) = ybin(ic,1:nbin1d+1)
                 if(wrap(ic,p).eq.0) then
                    call pgline(nbin1d+1,xbin1(1:nbin1d+1),ybin1(1:nbin1d+1))
                 else
                    call pgline(nbin1d,xbin1(1:nbin1d),ybin1(1:nbin1d))
                 end if
              end do
              call pgsls(4)
           end if
           
           
           !Plot max likelihood
           if(pllmax.ge.1) then
              ply = pldat(icloglmax,p,iloglmax)
              if(p.eq.8) ply = rev24(ply)
              if(p.eq.10.or.p.eq.12.or.p.eq.13) ply = rev360(ply)
              call pgsci(1)
              call pgsls(5)
              call pgline(2,(/ply,ply/),(/-1.e20,1.e20/))
           end if
           
           
           !Plot median and model value
           call pgsch(sch)
           
           do ic=1,nchains
              !Draw white lines
              if(nchains.gt.1) then
                 call pgslw(lw)
                 call pgsls(1); call pgsci(0)
                 if(pltrue.ge.1) call pgline(2,(/startval(ic,p,1),startval(ic,p,1)/),(/-1.e20,1.e20/))                    !True value
                 !if(plstart.ge.1) call pgline(2,(/startval(ic,p,2),startval(ic,p,2)/),(/-1.e20,1.e20/))                   !Starting value
                 if(p.ne.1) then !Not if plotting log(L)
                    if(plmedian.eq.1.or.plmedian.eq.3) call pgline(2,(/stats(ic,p,1),stats(ic,p,1)/),(/-1.e20,1.e20/))                          !Median
                    if(plrange.ge.1) then
                       if(nchains.lt.2) call pgline(2,(/ranges(ic,c0,p,1),ranges(ic,c0,p,1)/),(/-1.e20,1.e20/)) !Left limit of 90% interval
                       if(nchains.lt.2) call pgline(2,(/ranges(ic,c0,p,2),ranges(ic,c0,p,2)/),(/-1.e20,1.e20/)) !Right limit of 90% interval
                       if(nchains.eq.1) call pgline(2,(/ranges(ic,c0,p,3),ranges(ic,c0,p,3)/),(/-1.e20,1.e20/)) !Centre of 90% interval
                    end if
                 end if
              end if
              
              call pgslw(lw+1)
              !Draw coloured lines over the white ones
              !Median
              if((plmedian.eq.1.or.plmedian.eq.3) .and. p.ne.1) then
                 call pgsls(2); call pgsci(2); if(nchains.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
                 call pgline(2,(/stats(ic,p,1),stats(ic,p,1)/),(/-1.e20,1.e20/))
              end if
              
              !Plot ranges in PDF
              if(plrange.eq.1.and.p.ne.1) then
                 call pgsls(4); call pgsci(2); if(nchains.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
                 if(nchains.lt.2) call pgline(2,(/ranges(ic,c0,p,1),ranges(ic,c0,p,1)/),(/-1.e20,1.e20/)) !Left limit of 90% interval
                 if(nchains.lt.2) call pgline(2,(/ranges(ic,c0,p,2),ranges(ic,c0,p,2)/),(/-1.e20,1.e20/)) !Right limit of 90% interval
                 !if(nchains.eq.1) call pgline(2,(/ranges(ic,c0,p,3),ranges(ic,c0,p,3)/),(/-1.e20,1.e20/)) !Centre of 90% interval
              end if
              
              !Plot true value in PDF
              if(pltrue.eq.1) then !Plot true values
                 if(mergechains.ne.1.or.ic.le.1) then !The units of the true values haven't changed (e.g. from rad to deg) for ic>1 (but they have for the starting values, why?)
                    call pgsls(2); call pgsci(1)
                    plx = startval(ic,p,1)
                    if(p.eq.8) plx = rev24(plx)
                    if(p.eq.10.or.p.eq.12.or.p.eq.13) plx = rev360(plx)
                    call pgline(2,(/plx,plx/),(/-1.e20,1.e20/)) !True value
                    if(p.eq.8) then
                       call pgline(2,(/plx-24.,plx-24./),(/-1.e20,1.e20/)) !True value
                       call pgline(2,(/plx+24.,plx+24./),(/-1.e20,1.e20/)) !True value
                    end if
                    if(p.eq.10.or.p.eq.12.or.p.eq.13) then
                       call pgline(2,(/plx-360.,plx-360./),(/-1.e20,1.e20/)) !True value
                       call pgline(2,(/plx+360.,plx+360./),(/-1.e20,1.e20/)) !True value
                    end if
                 end if
              end if
              
              !Plot starting value in PDF
              !if(plstart.eq.1.and.abs((startval(ic,p,1)-startval(ic,p,2))/startval(ic,p,1)).gt.1.e-10) then
              !   call pgsls(4); call pgsci(1); if(nchains.gt.1) call pgsci(1)
              !   call pgline(2,(/startval(ic,p,2),startval(ic,p,2)/),(/-1.e20,1.e20/))
              !end if
              
              call pgsls(1)
              call pgsci(1)
           end do !ic
           
           
           
           
           !Print median, model value and range widths in panel title
           call pgslw(lw)
           call pgsci(1)
           ic = 1
           !if(nplvar.lt.7.or.nplvar.eq.9) then  !Three or less columns
           if(quality.ne.2.and.quality.ne.3.and.quality.ne.4) then  !Not a talk/poster/thesis
              if(nplvar.lt.5) then
                 if(p.eq.2.or.p.eq.3.or.p.eq.5.or.p.eq.6.or.p.eq.14.or.p.eq.15) then
                    write(str,'(A,F7.3,A5,F7.3,A9,F6.2,A1)')trim(pgvarns(p))//': mdl:',startval(ic,p,1),' med:',stats(ic,p,1),  &
                         !' \(2030):',abs(stats(ic,p,1)-startval(ic,p,1))/startval(ic,p,1)*100,'%'
                         ' \(2030):',ranges(ic,c0,p,5)*100,'%'
                 else
                    write(str,'(A,F7.3,A5,F7.3,A9,F7.3)')trim(pgvarns(p))//': mdl:',startval(ic,p,1),' med:',stats(ic,p,1),  &
                         !' \(2030):',abs(stats(ic,p,1)-startval(ic,p,1))  
                         ' \(2030):',ranges(ic,c0,p,5)
                 end if
              else  !if nplvar>=5
                 if(p.eq.2.or.p.eq.3.or.p.eq.5.or.p.eq.6.or.p.eq.14.or.p.eq.15) then
                    !write(str,'(A,F7.3,A9,F6.2,A1)')trim(pgvarns(p))//': mdl:',startval(ic,p,1),
                    write(str,'(A4,F7.3,A9,F6.2,A1)')'mdl:',startval(ic,p,1),  &
                         !' \(2030):',abs(stats(ic,p,1)-startval(ic,p,1))/startval(ic,p,1)*100,'%'  
                         ' \(2030):',ranges(ic,c0,p,5)*100,'%'
                 else
                    !write(str,'(A,F7.3,A9,F7.3)')trim(pgvarns(p))//': mdl:',startval(ic,p,1),
                    write(str,'(A4,F7.3,A9,F7.3)')'mdl:',startval(ic,p,1),  &
                         !' \(2030):',abs(stats(ic,p,1)-startval(ic,p,1))  
                         ' \(2030):',ranges(ic,c0,p,5)
                 end if
                 call pgsch(sch*1.2)
                 call pgptxt(xmin+0.05*dx,ymax*0.9,0.,0.,trim(pgvarnss(p)))
              end if
           end if
           
           
           if(quality.eq.2.or.quality.eq.3.or.quality.eq.4) then  !Talk/poster/thesis
              !if(p.eq.2.or.p.eq.3.or.p.eq.5.or.p.eq.6.or.p.eq.14.or.p.eq.15) then
              !   write(str,'(A9,F6.2,A1)')' \(2030):',ranges(ic,c0,p,5)*100,'%'
              !else
              !   write(str,'(A9,F7.3)')' \(2030):',ranges(ic,c0,p,5)
              !end if
              x0 = ranges(ic,c0,p,5)
              if(p.eq.2.or.p.eq.3.or.p.eq.5.or.p.eq.6.or.p.eq.14.or.p.eq.15) x0 = x0*100
              !print*,p,x0,nint(x0)
              if(x0.lt.0.01) write(str,'(F6.4)')x0
              if(x0.ge.0.01.and.x0.lt.0.1) write(str,'(F5.3)')x0
              if(x0.ge.0.1.and.x0.lt.1.) write(str,'(F4.2)')x0
              if(x0.ge.1.and.x0.lt.9.95) write(str,'(F3.1)')x0
              if(x0.ge.9.95.and.x0.lt.99.5) write(str,'(I2)')nint(x0)
              if(x0.ge.99.5) write(str,'(I3)')nint(x0)
              write(str,'(A)')'\(2030): '//trim(str)
              if(p.eq.2.or.p.eq.3.or.p.eq.5.or.p.eq.6.or.p.eq.14.or.p.eq.15) then
                 write(str,'(A)')trim(str)//'%'
              else
                 write(str,'(A)')trim(str)//trim(pgunits(p))
              end if
              call pgsch(sch*1.2)
              !call pgptxt(xmin+0.05*dx,ymax*0.9,0.,0.,trim(pgvarnss(p)))
              if(abs(xmin-xpeak).lt.abs(xmax-xpeak)) then !peak is at left, put varname at right
                 call pgptxt(xmax-0.05*dx,ymax*0.9,0.,1.,trim(pgvarnss(p)))
              else
                 call pgptxt(xmin+0.05*dx,ymax*0.9,0.,0.,trim(pgvarnss(p)))
              end if
           end if
           
           
           !if(nchains.gt.1) write(str,'(A,F7.3,A)')trim(pgvarns(p))//'  mdl: ',startval(ic,p,1),''
           !if(nchains.gt.1) write(str,'(A4,F7.3,A)')'mdl:',startval(ic,p,1),''
           !if(nchains.eq.2) write(str,'(A,3(A6,F7.3))')trim(pgvarns(p)),' mdl:',startval(ic,p,1),
           !                       ' med1:',stats(1,p,1),' med2:',stats(2,p,1)
           !if(nchains.eq.2.and.abs((startval(1,p,1)-startval(2,p,1))/startval(1,p,1)).gt.1.e-10)  &
           !                        !         write(str,'(A,A5,F7.3,A2,F7.3)')trim(pgvarns(p)),'mdl:',startval(1,p,1),', ',startval(2,p,1)  
           !     write(str,'(A4,F6.2,A2,F6.2)')'mdl:',startval(1,p,1),', ',startval(2,p,1)
           
           !Write the deltas of the two pdfs
           if(nchains.eq.2.) then
              write(str,'(A8)')'\(2030)'
              if(p.eq.2.or.p.eq.3.or.p.eq.5.or.p.eq.6.or.p.eq.14.or.p.eq.15) then
                 write(str1,'(A8,F6.2,A1)')'\(2030):',ranges(1,c0,p,5)*100.,'%'
                 write(str2,'(A8,F6.2,A1)')'\(2030):',ranges(2,c0,p,5)*100.,'%'
              else
                 write(str1,'(A8,F7.3)')'\(2030):',ranges(1,c0,p,5)
                 write(str2,'(A8,F7.3)')'\(2030):',ranges(2,c0,p,5)
              end if
           end if
           
           !if(p.eq.1) then
           !   str1 = ''
           !   str2 = ''
           !   !write(str,'(A)')trim(pgvarns(p))
           !   str = ''
           !end if
           
           !call pgsch(sch*0.9)
           call pgsch(sch*1.1)
           if(prvalues.eq.1.and.p.ne.1) then  !If not plotting log(L)
              if(nchains.eq.2) then
                 call pgsci(colours(mod(0,ncolours)+1))
                 call pgmtxt('T',0.5,0.25,0.5,trim(str1))
                 call pgsci(colours(mod(1,ncolours)+1))
                 call pgmtxt('T',0.5,0.75,0.5,trim(str2))
              else
                 if(quality.eq.2.or.quality.eq.3) call pgsci(2)
                 !call pgmtxt('T',0.2,0.5,0.5,trim(str))
                 !call pgptxt(ranges(ic,c0,p,3),ymax,0.,0.5,trim(str)) !Align with centre of 90%-probability range
                 call pgptxt((xmin+xmax)/2.,ymax,0.,0.5,trim(str)) !Centre
                 !print*,trim(str)
                 !call pgarro(ranges(ic,c0,p,3),0.95*ymax,ranges(ic,c0,p,1),0.95*ymax)
                 !call pgarro(ranges(ic,c0,p,3),0.95*ymax,ranges(ic,c0,p,2),0.95*ymax)
                 call pgsci(2)
                 call pgline(2,(/ranges(ic,c0,p,1),ranges(ic,c0,p,2)/),(/0.99*ymax,0.99*ymax/))  !Plot line at top over 90%-probability width
                 call pgsci(1)
              end if
           end if
           
           
           call pgsci(1)
           call pgsch(sch)
           !call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)
           !call pgbox('BNTS',0.0,0,'BNTS',0.0,0)
           call pgbox('BNTS',0.0,0,'',0.0,0)
           !print*,sch,lw,trim(str)
        end if !if(plot.eq.1) 
     end do !p
     
     if(savepdf.eq.1) close(30)
     
     if(plot.eq.1) then
        call pgsubp(1,1)
        call pgsvp(0.,1.,0.,1.)
        call pgswin(-1.,1.,-1.,1.)
        
        if(quality.eq.0) then
           !Remove also the pgsvp at the beginning of the plot
           string=' '
           do ic=1,1!nchains !Can't do this 10x for 10 chains
              !write(string,'(A,I7,A,I6)')trim(string)//'n:',ntot(ic),', nburn:',nburn(ic)
              write(string,'(A,I8)')trim(string)//'n:',sum(ntot(1:nchains0))
           end do
           call pgsch(sch*0.7)
           call pgmtxt('T',-0.7,0.5,0.5,trim(outputname)//'  '//trim(string))  !Print title
           call pgsch(sch)
        end if
        
        
        !Make sure gv auto-reloads on change
        !      if(file.ge.2) then
        !        call pgpage
        !        call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)
        !      end if
        
        call pgend
        
        if(file.ge.2) then
           if(file.eq.3) then
              i = system('eps2pdf pdfs.eps -o '//trim(outputdir)//'/'//trim(outputname)//'__pdfs.pdf >& /dev/null')
              if(i.ne.0) write(*,'(A,I6)')'  Error converting plot',i
           end if
           i = system('mv -f pdfs.eps '//trim(outputdir)//'/'//trim(outputname)//'__pdfs.eps')
        else if(file.eq.1) then
           i = system('convert -resize '//trim(bmpxpix)//' -depth 8 -unsharp '//trim(unsharppdf1d)//' pdfs.ppm '//trim(outputdir)//'/'//trim(outputname)//'__pdfs.png')
           if(i.ne.0) write(*,'(A,I6)')'  Error converting plot',i
           i = system('rm -f pdfs.ppm')
        end if
     end if !if(plot.eq.1)
     !write(*,*)''   
  end if !if(plpdf1d.eq.1)
  
  
  
  
  
  
  
  
  
  
  
  
  !***********************************************************************************************************************************      
  if(plpdf2d.ge.1) then
     ic = 1 !Can only do one chain
     if(plot.eq.0.and.savepdf.eq.1) write(*,'(A,$)')' Saving 2D pdfs...    '
     allocate(z(nbin2dx+1,nbin2dy+1),zs(nchs,nbin2dx+1,nbin2dy+1))
     if(plot.eq.1) then
        !if(prprogress.ge.2.and.update.eq.0) write(*,'(A)')' Plotting 2D pdfs...    '
        if(prprogress.ge.2.and.update.eq.0) write(*,'(A,$)')' 2D pdfs: '
        if(file.eq.0) then
           lw = 1
           lw2 = 1
           sch = 1.5
        end if
        if(file.ge.1) then
           if(file.ge.2) io = pgopen('pdf2d.eps'//trim(psclr))
           lw = 3
           lw2 = 2 !Font lw
           sch = 1.5
           if(quality.eq.3) then !Poster
              lw = 4
              lw2 = 3 !Font lw
              sch = 2.
           end if
        end if
        if(file.ge.2.and.io.le.0) then
           write(*,'(A,I4)')'Cannot open PGPlot device.  Quitting the programme',io
           goto 9999
        end if
        if(file.ge.2) call pgpap(pssz,psrat)
        if(file.ge.2) call pgscf(2)
        if(file.gt.1) call pginitl(colour,file,whitebg)
     end if !if(plot.eq.1)
     
     !NEW columns in dat: 1:logL 2:mc, 3:eta, 4:tc, 5:logdl, 6:spin, 7:kappa, 8: RA, 9:sindec,10:phase, 11:sinthJ0, 12:phiJ0, 13:alpha, 14:M1, 15:M2
     j1 = 2
     j2 = npar
     !j2 = 13 !Don't do M1,M2
     !j1 = 6
     !j2 = 7
     if(plotsky.eq.1) then
        j1 = 8
        j2 = 9
     end if
     
     if(savepdf.eq.1) then
        open(unit=30,action='write',form='formatted',status='replace',file=trim(outputdir)//'/'//trim(outputname)//'__pdf2d.dat')
        write(30,'(5I6,T100,A)')j1,j2,1,nbin2dx,nbin2dy,'Plot variable 1,2, total number of chains, number of bins x,y'
     end if
     
     npdf=0 !Count iterations to open windows with different numbers
     !do p1=j1,j2-1
     !   do p2=p1+1,j2
     do p1=j1,j2
        do p2=j1,j2
     !do p1=3,3
        !do p2=6,6
           
           !Skip some 2d pdfs to save time:
           !if(p1.eq.4.or.p1.eq.5.or.p1.eq.10.or.p1.eq.11.or.p1.eq.12.or.p1.eq.13) cycle
           !if(p2.eq.4.or.p2.eq.5.or.p2.eq.10.or.p2.eq.11.or.p2.eq.12.or.p2.eq.13) cycle
           
           !if(p1.eq.2.and.p2.ne.3 .or. p2.eq.2.or.p1.eq.3) cycle !Mc only with eta
           !if(p1.eq.6.and.p2.ne.7 .or. p2.eq.6.or.p1.eq.7) cycle !a_spin only with theta
           !if(p1.eq.8.and.p2.ne.9 .or. p2.eq.8.or.p1.eq.9) cycle !a_spin only with theta
           !if(p1.eq.14.and.p2.ne.15 .or. p2.eq.14.or.p1.eq.15) cycle !M1 only with M2
           
           !if(p1.ne.5.and.p1.ne.8.and.p1.ne.9) cycle  !Only position and distance
           !if(p2.ne.5.and.p2.ne.8.and.p2.ne.9) cycle  !Only position and distance
           
           plotthis = 0  !Determine to plot or save this combination of j1/j2 or p1/p2
           !if(p1.eq.2.and.p2.eq.3) plotthis = 1    !Mc-eta
           !if(p1.eq.6.and.p2.eq.7) plotthis = 1    !a-theta
           !if(p1.eq.8.and.p2.eq.9) plotthis = 1    !RA-dec
           !if(p1.eq.11.and.p2.eq.12) plotthis = 1  !theta/phi_Jo
           !if(p1.eq.12.and.p2.eq.11) plotthis = 1  !phi/theta_Jo
           !if(p1.eq.14.and.p2.eq.15) plotthis = 1  !M1-M2
           
           do i=1,npdf2d
              if(p1.eq.pdf2dpairs(i,1).and.p2.eq.pdf2dpairs(i,2)) plotthis = 1  !Use the data from the input file
           end do
           
           
           if(plotthis.eq.0) cycle
           
           !print*,p1,p2
           !write(*,'(A)')'   PDF 2D:  '//trim(varnames(p1))//'-'//trim(varnames(p2))
           if(prprogress.ge.2.and.update.eq.0) write(*,'(A,$)')trim(varnames(p1))//'-'//trim(varnames(p2))//' '
           
           if(plot.eq.1) then
              if(file.eq.0) then
                 npdf=npdf+1
                 write(str,'(I3,A3)')200+npdf,'/xs'
                 io = pgopen(trim(str))
                 call pgpap(scrsz,scrrat)
                 call pginitl(colour,file,whitebg)
              end if
              if(file.eq.1) then
                 io = pgopen('pdf2d.ppm/ppm')
                 call pgpap(bmpsz,bmprat)
                 call pginitl(colour,file,whitebg)
              end if
              if(file.lt.2.and.io.le.0) then
                 write(*,'(A,I4)')'Cannot open PGPlot device.  Quitting the programme',io
                 goto 9999
              end if

              !call pgscr(3,0.,0.5,0.)
              call pgsch(sch)
           end if
           
           xmin = minval(alldat(ic,p1,1:n(ic)))
           xmax = maxval(alldat(ic,p1,1:n(ic)))
           ymin = minval(alldat(ic,p2,1:n(ic)))
           ymax = maxval(alldat(ic,p2,1:n(ic)))
           dx = xmax - xmin
           dy = ymax - ymin
           !write(*,'(A,2F10.5)')'  Xmin,Xmax: ',xmin,xmax
           !write(*,'(A,2F10.5)')'  Ymin,Ymax: ',ymin,ymax
           
           xx(1:n(ic)) = alldat(ic,p1,1:n(ic)) !Parameter 1
           yy(1:n(ic)) = alldat(ic,p2,1:n(ic)) !Parameter 2
           zz(1:n(ic)) = alldat(ic,1,1:n(ic))   !Likelihood

           xmin = xmin - 0.05*dx
           xmax = xmax + 0.05*dx
           ymin = ymin - 0.05*dy
           ymax = ymax + 0.05*dy
           
           
           !Plot a cute sky map in 2D PDF
           if(plot.eq.1.and.plotsky.eq.1) then
              xmax = 18.2!14.
              xmin = 14.8!20.
              ymin = 20.!0.
              ymax = 50.!70.
              rat = 0.75
              call pgpap(11.,rat)
              dx = xmax - xmin
              dy = ymax - ymin
              if(abs(dx)*15.lt.dy/rat) then !Expand x
                 dx = dy/(15*rat)
                 a = (xmin+xmax)*0.5
                 xmin = a - 0.5*dx
                 xmax = a + 0.5*dx
              end if
              if(abs(dx)*15.gt.dy/rat) then !Expand y
                 dy = abs(dx)*rat*15
                 a = (ymin+ymax)*0.5
                 ymin = a - 0.5*dy
                 ymax = a + 0.5*dy
              end if
           end if
           
           
           
           !'Normalise' 2D PDF
           if(normpdf2d.le.2.or.normpdf2d.eq.4) then
              call bindata2d(n(ic),xx(1:n(ic)),yy(1:n(ic)),0,nbin2dx,nbin2dy,xmin,xmax,ymin,ymax,z,tr)  !Count number of chain elements in each bin
              if(normpdf2d.eq.1) z = max(0.,log10(z + 1.e-30))
              if(normpdf2d.eq.2) z = max(0.,sqrt(z + 1.e-30))
              if(normpdf2d.eq.4) then
                 call identify_2d_ranges(nival+1,ivals,nbin2dx+1,nbin2dy+1,z) !Get 2D probability ranges; identify to which range each bin belongs
                 call calc_2d_areas(p1,p2,changevar,nival+1,nbin2dx+1,nbin2dy+1,z,tr,probarea) !Compute 2D probability areas; sum the areas of all bins
                 trueranges2d(p1,p2) = truerange2d(z,nbin2dx+1,nbin2dy+1,startval(1,p1,1),startval(1,p2,1),tr)
                 !write(*,'(/,A23,2(2x,A21))')'Probability interval:','Equivalent diameter:','Fraction of a sphere:'
                 do i=1,nival+1
                    !write(*,'(I10,F13.2,2(2x,F21.5))')i,ivals(i),sqrt(probarea(i)/pi)*2,probarea(i)*(pi/180.)**2/(4*pi)  !4pi*(180/pi)^2 = 41252.961 sq. degrees in a sphere
                    probareas(p1,p2,i,1) = probarea(i)*(pi/180.)**2/(4*pi)  !Fraction of the sky
                    probareas(p1,p2,i,2) = sqrt(probarea(i)/pi)*2           !Equivalent diameter
                 end do
                 !write(*,'(A2,$)')'  '
              end if
           end if
           if(normpdf2d.eq.3) then
              call bindata2da(n(ic),xx(1:n(ic)),yy(1:n(ic)),zz(1:n(ic)),0,nbin2dx,nbin2dy,xmin,xmax,ymin,ymax,z,tr)  !Measure amount of likelihood in each bin
           end if
           
           
           !Swap RA boundaries for RA-Dec plot in 2D PDF
           if(p1.eq.8.and.p2.eq.9) then
              a = xmin
              xmin = xmax
              xmax = a
              dx = -dx
           end if
           
           z = z/(maxval(z)+1.e-30)
           
           if(plot.eq.1.and.plotsky.eq.1.and.file.ge.2) z = 1. - z !Invert grey scales
           
           
           !Plot 2D PDF
           if(plot.eq.1) then
              
              !Force boundaries
              if(1.eq.2.and.p1.eq.8.and.p2.eq.9) then
                 !xmin = 24.
                 !xmax = 0.
                 !ymin = -90.
                 !ymax = 90.
                 
                 xmin = 14.83440
                 xmax = 10.35767
                 ymin = -32.63011
                 ymax = 31.75267
              end if
              
              call pgsch(sch)
              !call pgsvp(0.12,0.95,0.12,0.95)
              call pgsvp(0.08*sch,0.95,0.08*sch,0.95)
              call pgswin(xmin,xmax,ymin,ymax)
              if(plotsky.eq.1.and.file.ge.2) then !Need dark background
                 !call pgsvp(0.,1.,0.,1.)
                 !call pgswin(0.,1.,0.,1.)
                 call pgsci(1)
                 call pgrect(xmin,xmax,ymin,ymax)
                 !call pgsci(0)
              end if
              
              !Plot the actual 2D PDF (grey scales or colour)
              if(plpdf2d.eq.1.or.plpdf2d.eq.2) then
                 if(normpdf2d.lt.4) call pggray(z,nbin2dx+1,nbin2dy+1,1,nbin2dx+1,1,nbin2dy+1,1.,0.,tr)
                 if(normpdf2d.eq.4) then
                    call pgscr(30,1.,1.,1.) !BG colour
                    if(nival+1.eq.2) then
                       call pgscr(31,1.,1.,0.) !Yellow
                       call pgscr(32,1.,0.,0.) !Red
                    end if
                    if(nival+1.eq.3) then
                       call pgscr(31,0.,0.,1.) !Blue
                       call pgscr(32,1.,1.,0.) !Yellow
                       call pgscr(33,1.,0.,0.) !Red
                    end if
                    if(nival+1.eq.4) then
                       call pgscr(31,0.,0.,1.) !Blue
                       call pgscr(32,0.,1.,0.) !Green
                       call pgscr(33,1.,1.,0.) !Yellow
                       !call pgscr(34,1.,0.5,0.) !Orange
                       call pgscr(34,1.,0.,0.) !Red
                    end if
                    if(nival+1.eq.5) then
                       call pgscr(31,0.,0.,1.) !Blue
                       call pgscr(32,0.,1.,0.) !Green
                       call pgscr(33,1.,1.,0.) !Yellow
                       call pgscr(34,1.,0.5,0.) !Orange
                       call pgscr(35,1.,0.,0.) !Red
                    end if
                    call pgscir(30,30+nival+1)
                    call pgimag(z,nbin2dx+1,nbin2dy+1,1,nbin2dx+1,1,nbin2dy+1,0.,1.,tr)
                 end if
              end if
              
              !Plot stars in 2D PDF (over the grey scales, but underneath contours, lines, etc)
              if(plotsky.eq.1) then
                 call pgswin(xmin*15,xmax*15,ymin,ymax) !Map works in degrees
                 call plotthesky(xmin*15,xmax*15,ymin,ymax)
                 call pgswin(xmin,xmax,ymin,ymax)
              end if
              call pgsci(1)
           end if !if(plot.eq.1)
           
           
           !Plot contours in 2D PDF
           if((plpdf2d.eq.1.or.plpdf2d.eq.3) .and. plot.eq.1) then
              if(normpdf2d.lt.4) then
                 ncont = 11
                 do i=1,ncont
                    cont(i) = 0.01 + 2*real(i-1)/real(ncont-1)
                    if(plotsky.eq.1) cont(i) = 1.-cont(i)
                 end do
                 ncont = 4 !Only use the first 4
              end if
              if(normpdf2d.eq.4) then
                 ncont = nival+1
                 do i=1,ncont
                    cont(i) = max(1. - real(i-1)/real(ncont-1),0.001)
                    !if(plotsky.eq.1) cont(i) = 1.-cont(i)
                 end do
              end if
              
              call pgsls(1)
              if(plotsky.eq.0 .and. normpdf2d.ne.4) then !First in bg colour
                 call pgslw(2*lw)
                 call pgsci(0)
                 !call pgcont(z,nbin2dx+1,nbin2dy+1,1,nbin2dx+1,1,nbin2dy+1,cont,4,tr)
                 call pgcont(z,nbin2dx+1,nbin2dy+1,1,nbin2dx+1,1,nbin2dy+1,cont(1:ncont),ncont,tr)
              end if
              call pgslw(lw)
              call pgsci(1)
              if(plotsky.eq.1) call pgsci(7)
              call pgcont(z,nbin2dx+1,nbin2dy+1,1,nbin2dx+1,1,nbin2dy+1,cont(1:ncont),ncont,tr)
           end if
           
           
           !Save binned 2D PDF data
           if(savepdf.eq.1) then
              write(30,'(3I6,T100,A)')ic,p1,p2,'Chain number and variable number 1,2'
              write(30,'(2ES15.7,T100,A)')startval(ic,p1,1:2),'True and starting value p1'
              write(30,'(2ES15.7,T100,A)')startval(ic,p2,1:2),'True and starting value p2'
              write(30,'(6ES15.7,T100,A)')stats(ic,p1,1:6),'Stats: median, mean, absvar1, absvar2, stdev1, stdev2 for p1'
              write(30,'(6ES15.7,T100,A)')stats(ic,p2,1:6),'Stats: median, mean, absvar1, absvar2, stdev1, stdev2 for p2'
              write(30,'(5ES15.7,T100,A)')ranges(ic,c0,p1,1:5),'Ranges: lower,upper limit, centre, width, relative width for p1'
              write(30,'(5ES15.7,T100,A)')ranges(ic,c0,p2,1:5),'Ranges: lower,upper limit, centre, width, relative width for p2'
              write(30,'(4ES15.7,T100,A)')xmin,xmax,ymin,ymax,'Xmin,Xmax,Ymin,Ymax of PDF'
              write(30,'(6ES15.7,T100,A)')tr,'Tr'              
              do i=1,nbin2dx+1
                 do j=1,nbin2dy+1
                    write(30,'(ES15.7,$)')z(i,j)
                 end do
                 write(30,'(1x)')
              end do
           end if
           
           
           
           !Plot true value, median, ranges, etc. in 2D PDF
           if(plot.eq.1) then
              call pgsci(1)
              
              !Plot max likelihood in 2D PDF
              if(pllmax.ge.1) then
                 call pgsci(1); call pgsls(5)
                 
                 plx = pldat(icloglmax,p1,iloglmax)
                 if(p1.eq.8) plx = rev24(plx)
                 if(p1.eq.10.or.p1.eq.12.or.p1.eq.13) plx = rev360(plx)
                 call pgline(2,(/plx,plx/),(/-1.e20,1.e20/)) !Max logL
                 if(p1.eq.8) then
                    call pgline(2,(/plx-24.,plx-24./),(/-1.e20,1.e20/)) !Max logL
                    call pgline(2,(/plx+24.,plx+24./),(/-1.e20,1.e20/)) !Max logL
                 end if
                 if(p1.eq.10.or.p1.eq.12.or.p1.eq.13) then
                    call pgline(2,(/plx-360.,plx-360./),(/-1.e20,1.e20/)) !Max logL
                    call pgline(2,(/plx+360.,plx+360./),(/-1.e20,1.e20/)) !Max logL
                 end if
                 
                 ply = pldat(icloglmax,p2,iloglmax)
                 if(p2.eq.8) ply = rev24(ply)
                 if(p2.eq.10.or.p2.eq.12.or.p2.eq.13) ply = rev360(ply)
                 call pgline(2,(/-1.e20,1.e20/),(/ply,ply/)) !Max logL
                 if(p2.eq.8) then
                    call pgline(2,(/-1.e20,1.e20/),(/ply-24.,ply-24./)) !Max logL
                    call pgline(2,(/-1.e20,1.e20/),(/ply+24.,ply+24./)) !Max logL
                 end if
                 if(p2.eq.10.or.p2.eq.12.or.p2.eq.13) then
                    call pgline(2,(/-1.e20,1.e20/),(/ply-360.,ply-360./)) !Max logL
                    call pgline(2,(/-1.e20,1.e20/),(/ply+360.,ply+360./)) !Max logL
                 end if
                 
                 call pgpoint(1,plx,ply,12)
              end if
              
              
              if(plotsky.eq.1) call pgsci(0)
              call pgsls(2)
              
              !Plot true value in 2D PDF
              if(plotsky.eq.0) then
                 !call pgline(2,(/startval(ic,p1,1),startval(ic,p1,1)/),(/-1.e20,1.e20/))
                 !call pgline(2,(/-1.e20,1.e20/),(/startval(ic,p2,1),startval(ic,p2,1)/))
                 
                 if(mergechains.ne.1.or.ic.le.1) then !The units of the true values haven't changed (e.g. from rad to deg) for ic>1 (but they have for the starting values, why?)
                    !x
                    call pgsls(2); call pgsci(1)
                    plx = startval(ic,p1,1)
                    if(p1.eq.8) plx = rev24(plx)
                    if(p1.eq.10.or.p1.eq.12.or.p1.eq.13) plx = rev360(plx)
                    call pgline(2,(/plx,plx/),(/-1.e20,1.e20/)) !True value
                    if(p1.eq.8) then
                       call pgline(2,(/plx-24.,plx-24./),(/-1.e20,1.e20/)) !True value
                       call pgline(2,(/plx+24.,plx+24./),(/-1.e20,1.e20/)) !True value
                    end if
                    if(p1.eq.10.or.p1.eq.12.or.p1.eq.13) then
                       call pgline(2,(/plx-360.,plx-360./),(/-1.e20,1.e20/)) !True value
                       call pgline(2,(/plx+360.,plx+360./),(/-1.e20,1.e20/)) !True value
                    end if
                    
                    !y
                    call pgsls(2); call pgsci(1)
                    ply = startval(ic,p2,1)
                    if(p2.eq.8) ply = rev24(ply)
                    if(p2.eq.10.or.p2.eq.12.or.p2.eq.13) ply = rev360(ply)
                    call pgline(2,(/-1.e20,1.e20/),(/ply,ply/)) !True value
                    if(p2.eq.8) then
                       call pgline(2,(/-1.e20,1.e20/),(/ply-24.,ply-24./)) !True value
                       call pgline(2,(/-1.e20,1.e20/),(/ply+24.,ply+24./)) !True value
                    end if
                    if(p2.eq.10.or.p2.eq.12.or.p2.eq.13) then
                       call pgline(2,(/-1.e20,1.e20/),(/ply-360.,ply-360./)) !True value
                       call pgline(2,(/-1.e20,1.e20/),(/ply+360.,ply+360./)) !True value
                    end if
                    
                    call pgpoint(1,plx,ply,18)
                 end if
              end if
              call pgsci(1)
              call pgsls(4)
              
              
              !Plot starting values in 2D PDF
              !call pgline(2,(/startval(ic,p1,2),startval(ic,p1,2)/),(/-1.e20,1.e20/))
              !call pgline(2,(/-1.e20,1.e20/),(/startval(ic,p2,2),startval(ic,p2,2)/))
              
              call pgsci(2)
              
              !Plot interval ranges in 2D PDF
              if(plrange.eq.2.or.plrange.eq.3) then
                 call pgsls(1)
                 call pgsch(sch*0.6)
                 call pgsah(1,45.,0.1)
                 a = 0.0166667*sch
                 call pgarro(ranges(ic,c0,p1,3),ymin+dy*a,ranges(ic,c0,p1,1),ymin+dy*a)
                 call pgarro(ranges(ic,c0,p1,3),ymin+dy*a,ranges(ic,c0,p1,2),ymin+dy*a)
                 a = 0.0333333*sch
                 call pgptxt(ranges(ic,c0,p1,3),ymin+dy*a,0.,0.5,'\(2030)\d90%\u')
                 a = 0.0233333*sch
                 call pgarro(xmin+dx*a,ranges(ic,c0,p2,3),xmin+dx*a,ranges(ic,c0,p2,1))
                 call pgarro(xmin+dx*a,ranges(ic,c0,p2,3),xmin+dx*a,ranges(ic,c0,p2,2))
                 a = 0.01*sch
                 call pgptxt(xmin+dx*a,ranges(ic,c0,p2,3),90.,0.5,'\(2030)\d90%\u')
              end if
              
              call pgsch(sch)
              call pgsls(2)
              
              
              !Plot medians in 2D PDF
              if(plmedian.eq.2.or.plmedian.eq.3) then
                 call pgline(2,(/stats(ic,p1,1),stats(ic,p1,1)/),(/-1.e20,1.e20/))
                 call pgline(2,(/-1.e20,1.e20/),(/stats(ic,p2,1),stats(ic,p2,1)/))
                 call pgpoint(1,stats(ic,p1,1),stats(ic,p2,1),18)
              end if
              
              call pgsls(1)
              
              
              !Big star at true position in 2D PDF
              if(plotsky.eq.1) then
                 call pgsch(sch*2)
                 call pgsci(9)
                 call pgpoint(1,startval(ic,p1,1),startval(ic,p2,1),18)
                 call pgsch(sch)
                 call pgsci(1)
              end if
              
              
              
              
              
              !Plot coordinate axes and axis labels in 2D PDF
              call pgsls(1)
              call pgslw(lw2)
              if(plotsky.eq.1) then
                 call pgsci(0)
                 call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0) !Box, ticks, etc in white
                 call pgsci(1)
                 call pgbox('N',0.0,0,'N',0.0,0) !Number labels in black
              else
                 call pgsci(1)
                 call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)
              end if
              call pgmtxt('B',2.2,0.5,0.5,trim(pgvarns(p1)))
              call pgmtxt('L',1.7,0.5,0.5,trim(pgvarns(p2)))
              
              
              !Print 2D probability ranges in title
              if(normpdf2d.eq.4) then
                 string = ' '
                 do c = 1,nival+1
                    a = probareas(p1,p2,c,2)
                    if(a.lt.1.) then
                       write(string,'(I3,A3,F5.2,A7)')nint(ivals(c)*100),'%:',a,'\(2218)'
                    else if(a.lt.10.) then
                       write(string,'(I3,A3,F4.1,A7)')nint(ivals(c)*100),'%:',a,'\(2218)'
                    else if(a.lt.100.) then
                       write(string,'(I3,A3,F5.1,A7)')nint(ivals(c)*100),'%:',a,'\(2218)'
                    else
                       write(string,'(I3,A3,I4,A7)')nint(ivals(c)*100),'%:',nint(a),'\(2218)'
                    end if
                    a = (real(c-1)/real(nival) - 0.5)*0.7 + 0.5
                    call pgsci(30+nival+2-c)
                    call pgmtxt('T',0.5,a,0.5,trim(string))  !Print title
                 end do
                 call pgsci(1)
              end if
              
              
              !Convert plot
              if(file.eq.1) then
                 call pgend
                 i = system('convert -resize '//trim(bmpxpix)//' -depth 8 -unsharp '//trim(unsharppdf2d)//' pdf2d.ppm '//trim(outputdir)//'/'//trim(outputname)//'__pdf2d__'//trim(varnames(p1))//'-'//trim(varnames(p2))//'.png')
                 if(i.ne.0) write(*,'(A,I6)')'  Error converting plot',i
                 i = system('rm -f pdf2d.ppm')
              end if
              if(file.ge.2) call pgpage
           end if !if(plot.eq.1)
           
        end do !p2
     end do !p1
        
     
     if(savepdf.eq.1) close(30)
     
     if(plot.eq.1) then
        if(file.ne.1) call pgend
        if(file.ge.2) then
           if(abs(j2-j1).le.1) then
              if(file.eq.3) i = system('eps2pdf pdf2d.eps  -o '//trim(outputdir)//'/'//trim(outputname)//'__pdf2d_'//trim(varnames(j1))//'-'//trim(varnames(j2))//'.pdf  >& /dev/null')
              i = system('mv -f pdf2d.eps '//trim(outputdir)//'/'//trim(outputname)//'__pdf2d_'//trim(varnames(j1))//'-'//trim(varnames(j2))//'.eps')
           else
              if(file.eq.3) i = system('eps2pdf pdf2d.eps  -o '//trim(outputdir)//'/'//trim(outputname)//'__pdf2d.pdf  >& /dev/null')
              i = system('mv -f pdf2d.eps '//trim(outputdir)//'/'//trim(outputname)//'__pdf2d.eps')
           end if
        end if
     end if !plot.eq.1
     
  end if !if(plpdf2d.eq.1)
  
  
  if(prprogress.ge.2.and.update.eq.0) write(*,'(A,$)')'done.  '
  if(plot.eq.1) write(*,*)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  !***********************************************************************************************************************************      
  !***********************************************************************************************************************************      
  !***********************************************************************************************************************************      
  
  
  
  
  
  !Write statistics to file
  if(savestats.ge.1) write(*,*)''
  if(savestats.ge.1.and.nchains.gt.1) write(*,'(A)')' ******   Cannot write statistics if the number of chains is greater than one   ******'
  if(savestats.ge.1.and.nchains.eq.1) then
     ic = 1 !Use chain 1
     o = 20 !Output port
     open(unit=o, form='formatted', status='replace',file=trim(outputdir)//'/'//trim(outputname)//'__statistics.dat')
     write(o,'(A)')trim(outputname)
     
     !Print general run and detector info:
     write(o,'(//,A,/)')'GENERAL INFORMATION:'
     !write(o,'(3(A,I))')'Npoints: ',n(ic),'  Nburn: ',nburn(ic),'  Ndet: ',ndet(ic)
     !write(o,*)''
     !write(o,'(6x,A10,A12,A8,A22,A8,A10)')'niter','nburn','seed','null likelihood','ndet','nchains'
     !write(o,'(6x,I10,I12,I8,F22.10,I8,I5.2,A2,I3.2)')n(ic),nint(isburn(ic)),seed(ic),nullh,ndet(ic),contrchains,'/',nchains0
     write(o,'(6x,4A12,A12,A5  A8,A22,A8)')'totiter','totlines','totpts','totburn','nchains','used','seed','null likelihood','ndet'
     write(o,'(6x,4I12,I12,I5, I8,F22.10,I8)')totiter,totlines,totpts,totlines-totpts,contrchains,nchains0,seed(ic),nullh,ndet(ic)
     write(o,*)''
     write(o,'(A14,A3,A18,4A12,A22,A17,3A14)')'Detector','Nr','SNR','f_low','f_high','before tc','after tc','Sample start (GPS)','Sample length','Sample rate','Sample size','FT size'
     do i=1,ndet(ic)
        write(o,'(A14,I3,F18.8,4F12.2,F22.8,F17.7,3I14)')detname(ic,i),detnr(ic,i),snr(ic,i),flow(ic,i),fhigh(ic,i),t_before(ic,i),t_after(ic,i),FTstart(ic,i),deltaFT(ic,i),samplerate(ic,i),samplesize(ic,i),FTsize(ic,i)
     end do
     write(o,*)''
     
     write(o,'(A,I11)')' t0:',nint(t0)
     
     !Print statistics
     write(o,'(///,A,/)')'BASIC STATISTICS:'
     write(o,'(A,2I3)')'Npar,ncol:',par2-par1+1,7
     write(o,'(A8,7A12)')'param.','model','median','mean','stdev1','stdev2','abvar1','abvar2'
     
     do p=par1,par2
        write(o,'(A8,7F12.6)')varnames(p),startval(ic,p,1),stats(ic,p,1),stats(ic,p,2),stdev1(p),stdev2(p),absvar1(p),absvar2(p)
     end do
     write(o,*)''
     
     
     !Print correlations:
     write(o,'(//,A,/)')'CORRELATIONS:'
     write(o,'(A,I3)')'Npar:',par2-par1+1
     write(o,'(A9,$)')''
     do p=par1,par2
        write(o,'(A10,$)')trim(varnames(p))
     end do
     write(o,*)''
     do p1=par1,par2
        write(o,'(A9,$)')trim(varnames(p1))
        do p2=par1,par2
           write(o,'(F10.5,$)')corrs(p1,p2)
        end do
        write(o,'(A)')'   '//trim(varnames(p1))
     end do
     
     
     !Print probability intervals:
     write(o,'(///,A,/)')'1D PROBABILITY INTERVALS:'
     write(o,'(A,I3)')'Nival:',nival
     write(o,'(A22,$)')'Interval:'
     do c=1,nival+1
        write(o,'(F21.5,A14,$)')ivals(c),''
     end do
     write(o,*)''
     
     write(o,'(A8,2x,$)')'param.'
     do c=1,nival+1
        !write(o,'(2x,2A9,A8,$)')'rng1','rng2','in rnge'
        write(o,'(2x,2A12,A9,$)')'centre','delta','in rnge'
     end do
     write(o,*)''
     do p=par1,par2
        write(o,'(A8,2x,$)')trim(varnames(p))
        do c=1,nival+1
           !write(o,'(2x,2F11.6,F6.3,$)')ranges(ic,c,p,1),ranges(ic,c,p,2),2*abs(startval(ic,p,1)-ranges(ic,c,p,3))/ranges(ic,c,p,4) !Defined with centre of prob. range
           write(o,'(2x,2F12.6,F7.3,$)')ranges(ic,c,p,3),ranges(ic,c,p,4),2*abs(startval(ic,p,1)-ranges(ic,c,p,3))/ranges(ic,c,p,4) !Defined with centre of prob. range
           if(startval(ic,p,1).gt.ranges(ic,c,p,1).and.startval(ic,p,1).lt.ranges(ic,c,p,2)) then
              write(o,'(A2,$)')'y'
           else
              write(o,'(A2,$)')'N'
           end if
        end do
        write(o,*)''
     end do
     
     
     
     !Print 2D intervals
     write(o,'(///,A,/)')'2D PROBABILITY INTERVALS:'
     write(o,'(A,I5)')'Npdf2d: ',npdf2d
     write(o,'(A,2I5)')'Nbin2dx,nbin2dy: ',nbin2dx,nbin2dy
     !write(o,*)''
     write(o,'(A28,$)')'Interval:'
     do c=1,nival+1
        write(o,'(F19.5,$)')ivals(c)
     end do
     write(o,*)''
     
     write(o,'(A9,A19,$)')'params.',''
     do c=1,nival+1
        write(o,'(A16,A3,$)')'delta','in'
     end do
     write(o,*)''
     do p=1,npdf2d
        p1 = pdf2dpairs(p,1)
        p2 = pdf2dpairs(p,2)
        write(o,'(2I4,2(2x,A8),2x,$)')p1,p2,trim(varnames(p1)),trim(varnames(p2))
        do c=1,nival+1
           write(o,'(2x,F14.8,$)')probareas(p1,p2,c,1)
           !print*,c,trueranges2d(p1,p2),nival+2-trueranges2d(p1,p2)
           if(c.ge.nival+2-trueranges2d(p1,p2) .and. trueranges2d(p1,p2).ne.0) then
              write(o,'(A3,$)')'y'
           else
              write(o,'(A3,$)')'n'
           end if
        end do
        write(o,*)''
     end do
     
     
     
     close(o) !Statistics output file
     if(savestats.eq.2) i = system('a2ps -1rf7 '//trim(outputdir)//'/'//trim(outputname)//'__statistics.dat -o '//trim(outputdir)//'/'//trim(outputname)//'__statistics.ps')
     write(*,*)''
     if(savestats.eq.1) write(*,'(A)')' Statistics saved in '//trim(outputname)//'__statistics.dat'
     if(savestats.eq.2) write(*,'(A)')' Statistics saved in '//trim(outputname)//'__statistics.dat,ps'
  end if !if(savestats.ge.1.and.nchains.eq.1) then
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  !***********************************************************************************************************************************      
  !***********************************************************************************************************************************      
  !***********************************************************************************************************************************      
  
  
  if(plmovie.eq.1) then
     !i = system('cpu_time_type=total_alltime')
     !call cpu_time(cputime)
     ts1 = timestamp()
     write(*,*)
     p = 2 !Mc
     
     !moviescheme = 1  !Left column: Chain, sigma and acceptance, right column: numbers and PDF
     !moviescheme = 2  !Upper panel: numbers and PDF, lower panel: chain
     moviescheme = 3  !Upper panel: numbers and PDF, lower panel: chain and logL
     
     do iframe = 0,nmovframes
        nplt = nint(real(iframe)/real(nmovframes)*maxval(ntot(1:nchains)))-1  !This is the line number, not the iteration number
        !nplt = nint(real(iframe)/real(nmovframes)*maxval(is))-1
        !if(prprogress.ge.2.and.update.eq.0) write(*,'(A)')' Plotting movie frame...'
        if(prprogress.ge.2.and.update.eq.0) then
           write(6,*)upline !Move cursor up 1 line
           write(*,'(A,I5,A1,I5,A,I7,A,$)')' Plotting movie frame',iframe,'/',nmovframes,'  (',nplt,' points)'
           !Print remaining time
           !call cpu_time(cputime)
           !if(file.eq.1) write(6,'(A,A9)')'   Est.time left:',tms(dble(cputime-cputime0)*(nmovframes-iframe)/3600.d0)   !This does not take into account child processes, like convert
           !cputime0 = cputime
           ts2 = timestamp()
           if(file.eq.1) write(6,'(A,A9)')'   Est.time left:',tms((ts2-ts1)*(nmovframes-iframe+1)/3600.d0)                 !Use the system clock
           ts1 = ts2
        end if
        
        !print*,iframe,nplt
        !write(framename,'(A20,I4.4,A4)')'movies/frames/frame_',iframe,'.ppm'
        write(framename,'(A,I4.4,A4)')'frame_',iframe,'.ppm'
        
        ic = 1
        
        !print*,iframe,nplt,nburn(ic),maxval(is(ic,:))
        
        !Determine median and ranges for this selection of the chain (nburn:nplt)
        if(nplt.gt.nburn(ic)) then
           
           !Determine the median
           !call rindexx(nplt-nburn(ic),pldat(ic,p,1:nplt-nburn(ic)),index(p,1:nplt-nburn(ic)))  !Sort
           x(ic,1:nplt-nburn(ic)) = pldat(ic,p,nburn(ic)+1:nplt)
           call rindexx(nplt-nburn(ic),x(ic,1:nplt-nburn(ic)),index(p,1:nplt-nburn(ic)))  !Sort
           !print*,x(ic,index(p,1:nplt-nburn(ic)))
           
           if(mod(nplt-nburn(ic),2).eq.0) then
              !median = 0.5*(pldat(ic,p,index(p,(nplt-nburn(ic))/2)) + pldat(ic,p,index(p,(nplt-nburn(ic))/2+1)))  !Centre = nb + (n-nb)/2 = (n+nb)/2
              median = 0.5*(x(ic,index(p,(nplt-nburn(ic))/2)) + x(ic,index(p,(nplt-nburn(ic))/2+1)))  !Centre = nb + (n-nb)/2 = (n+nb)/2
           else
              !median = pldat(ic,p,index(p,(nplt-nburn(ic)+1)/2))
              median = x(ic,index(p,(nplt-nburn(ic)+1)/2))
           end if
           !print*,'median:',median
           
           
           !Determine interval ranges
           c = c0
           ival = ivals(c)
           
           minrange = 1.e30
           !do i=1,floor(nplt*(1.-ival))
           do i=1,floor((nplt-nburn(ic))*(1.-ival))
              !x1 = pldat(ic,p,index(p,i))
              !x2 = pldat(ic,p,index(p,i+floor((nplt-nburn(ic))*ival)))
              x1 = x(ic,index(p,i))
              x2 = x(ic,index(p,i+floor((nplt-nburn(ic))*ival)))
              range = abs(x2 - x1)
              if(range.lt.minrange) then
                 minrange = range
                 y1 = x1
                 y2 = x2
              end if
              !write(*,'(2I6,7F12.8)')i,i+floor((nplt-nburn(ic))*ival),x1,x2,range,minrange,y1,y2,(y1+y2)/2.
           end do !i
           centre = (y1+y2)/2.
           !write(*,'(A8,4x,4F10.5,I4)')varnames(p),y1,y2,minrange,centre,wrap(ic,p)
           
           !Save ranges:
           range1 = y1
           range2 = y2
           drange = y2-y1
           if(p.eq.2.or.p.eq.3.or.p.eq.5.or.p.eq.6.or.p.eq.14.or.p.eq.15) drange = drange/((y1+y2)/2.)
           
           !print*,'ranges:',range1,range2,drange
        end if
        
        
        
        
        
        if(file.eq.0) io = pgopen('17/xs')
        if(file.eq.1) io = pgopen(trim(framename)//'/ppm')
        if(io.le.0) then
           write(*,'(A,I4)')'Cannot open PGPlot device.  Quitting the programme',io
           goto 9999
        end if
        !if(file.eq.0) call pgpap(12.,0.72)
        if(file.eq.0) call pgpap(scrsz*0.75,scrrat)
        !if(file.eq.1) call pgpap(20.,0.72)   !for 850x612, change convert below.  1:1.388 for mpeg, make the output image twice as big, rescale in the end
        if(file.ge.1) call pgpap(24.08,0.72) !for 1024x738, change convert below. 1:1.388 for mpeg, make the output image twice as big, rescale in the end
        call pgsch(1.)
        !call pgscr(3,0.,0.5,0.)
        call pginitl(colour,file,whitebg)
        lw = 2
        if(file.ge.1) lw = 5
        
        
        
        
        
        !***********************************************************************************************************************************      
        !Plot chain for this parameter
        
        
        !if(prprogress.ge.2.and.update.eq.0) write(*,'(A)')'  - parameter chains'
        call pgslw(lw)
        if(moviescheme.eq.1) call pgsvp(0.05,0.35,0.65,0.95)
        if(moviescheme.eq.2) call pgsvp(0.05,0.95,0.05,0.25)
        if(moviescheme.eq.3) call pgsvp(0.20,0.95,0.20,0.35)
        ic = 1
        xmin = 0.
        !xmax = real(maxval(ntot(1:nchains0)))
        xmax = maxval(is)
        dx = abs(xmax-xmin)*0.01
        if(moviescheme.eq.2) dx = 0.
        !if(moviescheme.eq.3) dx = abs(xmax-xmin)*0.05
        ymin =  1.e30
        ymax = -1.e30
        do ic=1,nchains0
           ymin = min(ymin,minval(pldat(ic,p,10:ntot(ic))))
           ymax = max(ymax,maxval(pldat(ic,p,10:ntot(ic))))
           !if(moviescheme.eq.3) then
           !   ymin = min(ymin,minval(pldat(ic,p,1:ntot(ic))))
           !   ymax = max(ymax,maxval(pldat(ic,p,1:ntot(ic))))
           !end if
           if(pltrue.eq.1) then
              ymin = minval((/startval(ic,p,1),ymin/))
              ymax = maxval((/startval(ic,p,1),ymax/))
           end if
           if(plstart.eq.1) then
              ymin = minval((/startval(ic,p,2),ymin/))
              ymax = maxval((/startval(ic,p,2),ymax/))
           end if
        end do
        dy = abs(ymax-ymin)*0.05
        !if(moviescheme.eq.3) dy = abs(ymax-ymin)*0.01
        !if(offsetrun.eq.1) dy = abs(ymax-ymin)*0.02
        !if(offsetrun.eq.1) dy = -1*dy
        
        call pgswin(xmin-dx,xmax+dx,ymin-dy,ymax+dy)
        !call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)
        if(moviescheme.eq.1) call pgbox('BCTS',0.0,0,'BCNTS',0.0,0)
        if(moviescheme.eq.2) call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)
        !if(moviescheme.eq.3) call pgbox('BCTS',0.0,0,'BCNTS',0.0,0)
        if(moviescheme.eq.3) call pgbox('BCTS',0.0,0,'BCNTS',10.**floor(log10(abs(ymax-ymin))),2)
        
        !call pgslw(1)
        do ic=1,nchains0
           !call pgsci(mod(ic*2,10))
           call pgsci(5)
           if(nchains0.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
           !do i=1,ntot(ic),chainpli
           if(chainsymbol.eq.0) then !Plot lines rather than symbols
              call pgline(nplt-1,is(ic,2:nplt),pldat(ic,p,2:nplt))
           else
              call pgpoint(1,0.,startval(ic,p,2),chainsymbol)       !Starting value
              do i=1,nplt,chainpli
                 call pgpoint(1,is(ic,i),pldat(ic,p,i),chainsymbol)
              end do
           end if
        end do
        !call pgslw(lw)
        
        do ic=1,nchains0
           call pgsci(2)
           call pgsls(2)
           if(nchains0.gt.1) call pgsci(1)
           call pgline(2,(/-1.e20,1.e20/),(/startval(ic,p,1),startval(ic,p,1)/)) !True value
           call pgsci(6)
           !if(plburn.ge.1) call pgline(2,(/isburn(ic),isburn(ic)/),(/-1.e20,1.e20/))
           if(plburn.ge.1.and.isburn(ic).lt.is(ic,ntot(ic))) call pgline(2,(/isburn(ic),isburn(ic)/),(/-1.e20,1.e20/))
           call pgsci(2)
           call pgsls(4)
           if(nchains0.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
           call pgline(2,(/-1.e20,1.e20/),(/startval(ic,p,2),startval(ic,p,2)/))  !Starting value
        end do
        call pgsci(1)
        call pgsls(1)
        if(moviescheme.eq.3) then
           call pgsvp(0.00,0.20,0.20,0.35)
           call pgswin(0.,1.,0.,1.)
           call pgsch(1.5)
           call pgslw(2*lw)
           call pgptxt(0.2,0.5,0.0,0.0,'Chain:')
           call pgsch(1.)
           call pgslw(lw)
        else
           !call pgmtxt('T',1.,0.5,0.5,'Chain')
           call pgmtxt('T',-1.5,0.05,0.0,'Chain')
        end if
        
        
        
        
        !***********************************************************************************************************************************            
        !Plot log(L)
        
        if(moviescheme.eq.3) then
           
           p = 1
           !if(prprogress.ge.2.and.update.eq.0) write(*,'(A)')'  - logL'
           call pgslw(lw)
           if(moviescheme.eq.3) call pgsvp(0.20,0.95,0.05,0.20)
           ic = 1
           xmin = 0.
           !xmax = real(maxval(ntot(1:nchains0)))
           xmax = maxval(is)
           dx = abs(xmax-xmin)*0.01
           !if(moviescheme.eq.3) dx = abs(xmax-xmin)*0.05
           ymin =  1.e30
           ymax = -1.e30
           do ic=1,nchains0
              ymin = min(ymin,minval(pldat(ic,p,1:ntot(ic))))
              ymax = max(ymax,maxval(pldat(ic,p,1:ntot(ic))))
           end do
           dy = abs(ymax-ymin)*0.05
           !if(moviescheme.eq.3) dy = abs(ymax-ymin)*0.1
           
           call pgswin(xmin-dx,xmax+dx,ymin-dy,ymax+dy)
           if(offsetrun.eq.1) call pgswin(xmin-dx,xmax+dx,ymin,ymax+2*dy)
           !if(moviescheme.eq.3) call pgbox('BCNTS',0.0,0,'BCNTS',nint((ymax/2.)/100)*100.,2)
           if(moviescheme.eq.3) call pgbox('BCNTS',0.0,0,'BCNTS',10.**floor(log10(abs(ymax-ymin))),2)
           
           do ic=1,nchains0
              !call pgsci(mod(ic*2,10))
              call pgsci(5)
              if(nchains0.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
              !do i=1,ntot(ic),chainpli
              if(chainsymbol.eq.0) then !Plot lines rather than symbols
                 call pgline(nplt-1,is(ic,2:nplt),pldat(ic,p,2:nplt))
              else
                 call pgpoint(1,0.,startval(ic,p,2),chainsymbol)       !Starting value
                 do i=1,nplt,chainpli
                    call pgpoint(1,is(ic,i),pldat(ic,p,i),chainsymbol)
                 end do
              end if
           end do
           
           do ic=1,nchains0
              call pgsls(2)
              if(pltrue.eq.1) then
                 call pgsci(2)
                 if(nchains0.gt.1) call pgsci(1)
                 call pgline(2,(/-1.e20,1.e20/),(/startval(ic,p,1),startval(ic,p,1)/))  !True value
              end if
              call pgsci(6)
              !call pgline(2,(/real(nburn(ic)),real(nburn(ic))/),(/-1.e20,1.e20/))
              !if(plburn.ge.1) call pgline(2,(/isburn(ic),isburn(ic)/),(/-1.e20,1.e20/))              !Burn-in
              if(plburn.ge.1.and.isburn(ic).lt.is(ic,ntot(ic))) call pgline(2,(/isburn(ic),isburn(ic)/),(/-1.e20,1.e20/))
           !   call pgsci(2)
           !   call pgsls(4)
           !   if(nchains0.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
           !   call pgline(2,(/-1.e20,1.e20/),(/startval(ic,p,2),startval(ic,p,2)/)) !Starting value
           end do
           call pgsci(1)
           call pgsls(1)
           !call pgmtxt('T',-1.5,0.05,0.0,'log(L)')
           if(moviescheme.eq.3) then
              call pgsvp(0.00,0.20,0.05,0.20)
              call pgswin(0.,1.,0.,1.)
              call pgsch(1.5)
              call pgslw(2*lw)
              call pgptxt(0.2,0.5,0.0,0.0,'log(L):')
              call pgsch(1.)
              call pgslw(lw)
           end if
        end if
        
        
        
        
        !***********************************************************************************************************************************            
        !Plot sigma values ('jump size')
        !if(prprogress.ge.2.and.update.eq.0) write(*,'(A)')'  - sigma'
        
        if(moviescheme.eq.1) then
           call pgsvp(0.05,0.35,0.35,0.65)
           
           xmax = -1.e30
           ymin =  1.e30
           ymax = -1.e30
           do ic=1,nchains0
              xmin = 0.
              xmax = max(xmax,real(ntot(ic)))
              dx = abs(xmax-xmin)*0.01
              ymin = min(ymin,minval(sig(p,ic,10:ntot(ic))))
              ymax = max(ymax,maxval(sig(p,ic,10:ntot(ic))))
              dy = abs(ymax-ymin)*0.05
           end do
           call pgswin(xmin-dx,xmax+dx,ymin-dy,ymax+dy)
           !call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)
           call pgbox('BCTS',0.0,0,'BCNTS',0.0,0)
           
           do ic=1,nchains0
              !call pgsci(mod(ic*2,10))
              if(nchains.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
              !        !do i=1,ntot(ic),chainpli
              !do i=1,nplt,chainpli
              do i=ic,nplt,chainpli !Start at ic to reduce overplotting
                 call pgpoint(1,is(ic,i),sig(p,ic,i),chainsymbol)
              end do
           end do
           
           call pgsls(2)
           call pgsci(6)
           do ic=1,nchains0
              !call pgline(2,(/real(nburn(ic)),real(nburn(ic))/),(/-1.e20,1.e20/))
              !if(plburn.ge.1) call pgline(2,(/isburn(ic),isburn(ic)/),(/-1.e20,1.e20/))
              if(plburn.ge.1.and.isburn(ic).lt.is(ic,ntot(ic))) call pgline(2,(/isburn(ic),isburn(ic)/),(/-1.e20,1.e20/))
           end do
           call pgsci(1)
           call pgsls(1)
           !call pgmtxt('T',1.,0.5,0.5,'Sigma')
           call pgmtxt('T',-1.5,0.05,0.0,'\(2144)')
           
        end if
        
        
        !***********************************************************************************************************************************            
        !Plot acceptance rate
        !if(prprogress.ge.2.and.update.eq.0) write(*,'(A)')'  - acceptance rate'
        
        if(moviescheme.eq.1) then
           call pgsvp(0.05,0.35,0.05,0.35)
           
           xmax = -1.e30
           ymin =  1.e30
           ymax = -1.e30
           do ic=1,nchains0
              xmin = 0.
              xmax = max(xmax,real(ntot(ic)))
              dx = abs(xmax-xmin)*0.01
              do i=1,ntot(ic)
                 if(acc(p,ic,i).gt.1.e-10 .and. acc(p,ic,i).lt.1.-1.e-10) then
                    n0 = i
                    exit
                 end if
              end do
              n0 = n0+10
              ymin = min(ymin,minval(acc(p,ic,n0:ntot(ic))))
              ymax = max(ymax,maxval(acc(p,ic,n0:ntot(ic))))
              dy = abs(ymax-ymin)*0.05
           end do
           
           call pgsci(1)
           call pgsls(1)
           call pgswin(xmin-dx,xmax+dx,ymin-dy,ymax+dy)
           call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)
           
           
           call pgsci(2)
           call pgsls(2)
           call pgline(2,(/-1.e20,1.e20/),(/0.25,0.25/))
           call pgsci(6)
           do ic=1,nchains0
              !call pgline(2,(/real(nburn(ic)),real(nburn(ic))/),(/-1.e20,1.e20/))
              !if(plburn.ge.1) call pgline(2,(/isburn(ic),isburn(ic)/),(/-1.e20,1.e20/))
              if(plburn.ge.1.and.isburn(ic).lt.is(ic,ntot(ic))) call pgline(2,(/isburn(ic),isburn(ic)/),(/-1.e20,1.e20/))
           end do
           
           do ic=1,nchains0
              !call pgsci(mod(ic*2,10))
              if(nchains.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
              !do i=1,ntot(ic),chainpli
              !do i=1,nplt,chainpli
              do i=ic,nplt,chainpli !Start at ic to reduce overplotting
                 call pgpoint(1,is(ic,i),acc(p,ic,i),chainsymbol)
              end do
           end do
           
           call pgsci(1)
           call pgsls(1)
           !call pgmtxt('T',1.,0.5,0.5,'Acceptance')
           call pgmtxt('T',-1.5,0.05,0.0,'Acceptance')
           
        end if
        
        
        
        
        
        
        !***********************************************************************************************************************************            
        !Plot 1D pdf
        !if(prprogress.ge.2.and.update.eq.0) write(*,'(A)')'  - pdf'
        
        if(moviescheme.eq.1) call pgsvp(0.45,0.95,0.05,0.8)
        if(moviescheme.eq.2) call pgsvp(0.25,0.95,0.32,0.999)
        if(moviescheme.eq.3) call pgsvp(0.27,0.95,0.42,0.999)
        
        call pgsfs(fillpdf)
        
        !Set x-ranges, bin the data and get y-ranges
        xmin = 1.e30
        xmax = -1.e30
        p = 2 !Mc
        do ic=1,nchains0
           xmin = min(xmin,minval(pldat(ic,p,1:ntot(ic))))
           xmax = max(xmax,maxval(pldat(ic,p,1:ntot(ic))))
        end do
        dx = xmax - xmin
        if(offsetrun.eq.0) then
           xmin = xmin - 0.1*dx
           xmax = xmax + 0.1*dx
        end if
        
        ic = 1
        if(mergechains.eq.1) then
           n1 = 1
           n2 = 1
           do ic=1,nchains0
              if(nplt.gt.nburn(ic)) then
                 n2 = n1 + min(nplt,ntot(ic))-nburn(ic) - 1
                 x(1,n1:n2) = pldat(ic,p,nburn(ic)+1:min(nplt,ntot(ic)))
                 n1 = n2 + 1
              end if
           end do
        end if
        
        ic = 1
        if(nplt.gt.nburn(ic)) then
           do ic=1,nchains
              n1 = 1
              if(mergechains.ne.1) then
                 n2 = nplt-nburn(ic)
                 x(ic,1:n2) = pldat(ic,p,nburn(ic)+1:nplt)
              end if
              xmin1 = minval(x(ic,1:n2))
              xmax1 = maxval(x(ic,1:n2))
              call bindata(n2,x(ic,1:n2),1,nbin1d,xmin1,xmax1,xbin1,ybin1) !Count the number of points in each bin
              
              
              
              if(normpdf1d.eq.2) then !Normalise the height of the PDF
                 ybin1 = ybin1/maxval(ybin1)
              else !Normalise the SURFACE, not the height (because of different bin size).  This is the default
                 norm = 0.
                 do i=1,nbin1d+1
                    norm = norm + ybin1(i)
                 end do
                 norm = norm*(xmax1-xmin1)
                 ybin1 = ybin1/norm
              end if
              
              !Smoothen
              ybin2 = ybin1
              if(smooth.gt.1) then
                 !i0 = nbin1d/10
                 !i0 = nint(min(max(real(nbin1d)/real(smooth),1.0),real(nbin1d)/2.))
                 i0 = min(max(smooth,1),floor(real(nbin1d)/2.))
                 do i=1+i0,nbin1d+1-i0
                    coefs1(1:2*i0+1) = ybin1(i-i0:i+i0)
                    call savgol(coefs1(1:2*i0+1),2*i0+1,i0,i0,0,4)
                    do i1=1,i0+1
                       coefs(i0-i1+2) = coefs1(i1)
                    end do
                    do i1 = i0+2,2*i0+1
                       coefs(3*i0+3-i1) = coefs1(i1)
                    end do
                    ybin2(i) = 0.
                    do i1=1,2*i0+1
                       ybin2(i) = ybin2(i) + coefs(i1) * ybin1(i+i1-i0-1)
                    end do
                 end do
                 ybin1 = ybin2
              end if
              xbin(ic,1:nbin1d+1) = xbin1(1:nbin1d+1)
              ybin(ic,1:nbin1d+1) = ybin1(1:nbin1d+1)
           end do
           
           ymin = 0.
           ymax = -1.e30
           do ic=1,nchains
              ymax = max(ymax,maxval(ybin(ic,1:nbin1d+1)))
           end do
           ymax = ymax*1.05
           
           call pgsch(sch)
           call pgswin(xmin,xmax,ymin,ymax)
           !call pgswin(xmin1,xmax1,ymin,ymax)
           
           
           call pgsci(1)
           call pgslw(lw)
           
           
           
           !print*,'Plot PDF'
           !Plot 1D PDF
           do ic=1,nchains
              if(fillpdf.ge.3) call pgshs(45.0*(-1)**ic,1.0,real(ic)/real(nchains0)) !Set hatch style: angle = +-45deg, phase between 0 and 1 (1/nchains0, 2/nchains0, ...)
              if(nchains.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
              xbin1(1:nbin1d+1) = xbin(ic,1:nbin1d+1)
              ybin1(1:nbin1d+1) = ybin(ic,1:nbin1d+1)
              if(wrap(ic,p).eq.0) then
                 if(nchains.eq.1) call pgsci(15)
                 call pgpoly(nbin1d+2,(/xbin1(1),xbin1(1:nbin1d+1)/),(/0.,ybin1(1:nbin1d+1)/))
                 call pgsci(1)
                 if(nchains.eq.1) call pgsci(2)
                 call pgline(nbin1d+1,xbin1(1:nbin1d+1),ybin1(1:nbin1d+1)) !:nbin1d) ?
                 
                 call pgline(2,(/xbin1(1),xbin1(1)/),(/0.,ybin1(1)/)) !Fix the loose ends
                 call pgline(2,(/xbin1(nbin1d+1),xbin1(nbin1d+1)/),(/ybin1(nbin1d+1),0./))
              else
                 plshift = real(2*pi)
                 if(changevar.eq.1) plshift = 360.
                 if(changevar.eq.1.and.p.eq.8) plshift = 24. !RA in hours
                 if(nchains.eq.1) call pgsci(15)
                 call pgpoly(nbin1d+3,(/xbin1(1),xbin1(1:nbin1d),xbin1(1)+plshift,xbin1(1)+plshift/),(/0.,ybin1(1:nbin1d),ybin1(1),0./))
                 call pgsci(1)
                 if(nchains.eq.1) call pgsci(2)
                 call pgline(nbin1d,xbin1(1:nbin1d),ybin1(1:nbin1d))
                 
                 call pgsls(4)
                 call pgline(nbin1d+1,(/xbin1(1:nbin1d)-plshift,xbin1(1)/),(/ybin1(1:nbin1d),ybin1(1)/))
                 call pgline(nbin1d,xbin1+plshift,ybin1)
                 call pgsls(1)
                 call pgline(2,(/xbin1(nbin1d),xbin1(1)+plshift/),(/ybin1(nbin1d),ybin1(1)/))
              end if
           end do !ic
           
           !print*,'Plot PDF lines'
           !Plot lines again over surface of overlapping distributions
           if(nchains.gt.1) then
              call pgsls(4)
              do ic=1,nchains
                 call pgsci(1)
                 xbin1(1:nbin1d+1) = xbin(ic,1:nbin1d+1)
                 ybin1(1:nbin1d+1) = ybin(ic,1:nbin1d+1)
                 if(wrap(ic,p).eq.0) then
                    call pgline(nbin1d+1,xbin1(1:nbin1d+1),ybin1(1:nbin1d+1))
                 else
                    call pgline(nbin1d,xbin1(1:nbin1d),ybin1(1:nbin1d))
                 end if
              end do
              call pgsls(4)
           end if
           
           
           
           
           !print*,'Plot median'
           !Plot median and model value
           call pgsch(sch)
           
           do ic=1,nchains
              !Draw white lines
              if(nchains.gt.1) then
                 call pgslw(lw)
                 call pgsls(1); call pgsci(0)
                 call pgline(2,(/startval(ic,p,1),startval(ic,p,1)/),(/-1.e20,1.e20/))
                 call pgline(2,(/startval(ic,p,2),startval(ic,p,2)/),(/-1.e20,1.e20/))
                 call pgline(2,(/median,median/),(/-1.e20,1.e20/))
                 call pgline(2,(/range1,range1/),(/-1.e20,1.e20/)) !Left limit of 90% interval
                 call pgline(2,(/range2,range2/),(/-1.e20,1.e20/)) !Right limit of 90% interval
              end if
              
              call pgslw(lw+1)
              !Draw coloured lines over the white ones
              !Median
              call pgsls(2); call pgsci(defcolour); if(nchains.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
              call pgline(2,(/median,median/),(/-1.e20,1.e20/))
              
              !Ranges
              call pgsls(4); call pgsci(defcolour); if(nchains.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
              call pgline(2,(/range1,range1/),(/-1.e20,1.e20/)) !Left limit of 90% interval
              call pgline(2,(/range2,range2/),(/-1.e20,1.e20/)) !Right limit of 90% interval
              
              !True value
              call pgsls(2); call pgsci(1); if(nchains.gt.1) call pgsci(1)
              call pgline(2,(/startval(ic,p,1),startval(ic,p,1)/),(/-1.e20,1.e20/))
              
              !Starting value
              !if(abs((startval(ic,p,1)-startval(ic,p,2))/startval(ic,p,1)).gt.1.e-10) then
              !   call pgsls(4); call pgsci(1); if(nchains.gt.1) call pgsci(1)
              !   call pgline(2,(/startval(ic,p,2),startval(ic,p,2)/),(/-1.e20,1.e20/))
              !end if
              
              call pgsls(1)
              call pgsci(1)
           end do
           
           
           !print*,'Print median'
           !Print median, model value and range widths
           call pgslw(lw)
           call pgsci(1)
           ic = 1
           if(moviescheme.eq.1) then
              if(p.eq.2.or.p.eq.3.or.p.eq.5.or.p.eq.6.or.p.eq.14.or.p.eq.15) then
                 write(str,'(A,F7.3,A5,F7.3,A9,F6.2,A1)')'mdl:',startval(ic,p,1),' med:',median,' \(2030):',drange*100,'%'
              else
                 write(str,'(A,F7.3,A5,F7.3,A9,F7.3)')'mdl:',startval(ic,p,1),' med:',median,' \(2030):',drange
              end if
           end if
           if(moviescheme.eq.2.or.moviescheme.eq.3) then
              call pgslw(lw)
              call pgsch(1.)
              call pgbox('BNTS',0.0,0,'',0.0,0) !Box for the PDF
              
              call pgsvp(0.05,0.25,0.35,0.999)
              call pgswin(0.,1.,1.,0.)
              call pgsch(2.)
              call pgslw(lw*2)
              
              call pgptxt(0.5,0.15,0.0,0.5,trim(pgvarns(p)))
              call pgslw(lw*2)
              call pgsch(1.4)
              
              call pgptxt(0.04,0.35,0.,0.,'Model:') 
              write(str,'(F7.3)')startval(ic,p,1)
              call pgptxt(0.55,0.35,0.,0.,trim(str)) 
              
              call pgptxt(0.04,0.5,0.,0.,'Median:') 
              write(str,'(F7.3)')median
              call pgptxt(0.55,0.5,0.,0.,trim(str)) 
              
              call pgptxt(0.04,0.65,0.,0.,'\(2030)\d90%\u:')
              write(str,'(F7.3)')max(drange,0.001)
              if(p.eq.2.or.p.eq.3.or.p.eq.5.or.p.eq.6.or.p.eq.14.or.p.eq.15) write(str,'(F6.2,A1)')max(drange*100,0.01),'%' 
              call pgptxt(0.55,0.65,0.,0.,trim(str)) 

              call pgsch(1.)
              call pgslw(lw)
              
              write(str,'(A,ES9.2)')' Iteration:',max(is(ic,max(nplt,1)),0.)
              call pgptxt(0.03,0.9,0.,0.,trim(str)) 
           end if
           
        else
           call pgsch(sch)
           call pgswin(xmin,xmax,ymin,ymax)
        !end if  !if(nplt.gt.nburn(ic)) 
           
           if(moviescheme.eq.1) then
              call pgsch(1.5)
              call pgslw(lw+2)
              call pgmtxt('T',2.5,0.5,0.5,trim(pgvarns(p)))
              call pgslw(lw)
              call pgsch(1.)
              !if(iframe.gt.0) call pgmtxt('T',1.,0.5,0.5,trim(str)) 
              if(nplt.gt.nburn(ic)) call pgmtxt('T',1.,0.5,0.5,trim(str)) 
              call pgbox('BNTS',0.0,0,'',0.0,0)
           end if
           
           if(moviescheme.eq.2.or.moviescheme.eq.3) then
              call pgslw(lw)
              call pgsch(1.)
              call pgbox('BNTS',0.0,0,'',0.0,0) !Box for the PDF
              call pgsls(2)
              call pgline(2,(/startval(ic,p,1),startval(ic,p,1)/),(/-1.e20,1.e20/))
              call pgsls(1)
              
              call pgsvp(0.05,0.25,0.35,0.999)
              call pgswin(0.,1.,1.,0.)
              call pgsch(2.)
              call pgslw(lw*2)
              
              call pgptxt(0.5,0.15,0.0,0.5,trim(pgvarns(p)))
              call pgslw(lw*2)
              call pgsch(1.4)
              
              call pgptxt(0.04,0.35,0.,0.,'Model:') 
              write(str,'(F7.3)')startval(ic,p,1)
              call pgptxt(0.55,0.35,0.,0.,trim(str)) 
              
              
              call pgsch(1.)
              call pgslw(lw)
              write(str,'(A,ES9.2)')' Iteration:',max(is(ic,max(nplt,1)),0.)
              call pgptxt(0.03,0.9,0.,0.,trim(str)) 
           end if
        end if  !if(nplt.gt.nburn(ic)) 

        !************************************************************************************************************************************
        
        
        
        
        
        
        
        call pgend
        
        !if(file.ge.1) i = system('convert -resize 850x612 -depth 8 -unsharp '//trim(unsharppdf1d)//' '//trim(framename)//' '//trim(framename))  !Rescale the output frame
        if(file.ge.1) then
           call pgend
           write(command,'(A)')'convert -resize 1024x738 -depth 8 -unsharp '//trim(unsharppdf1d)//' '//trim(framename)//' '//trim(framename)
           i = system(trim(command))  !Rescale the output frame
           if(i.ne.0) write(*,'(A,I6)')'  Error converting plot',i
        end if
        
     end do
     
  end if
  
  
  
  
  
  
  
  
  
  if(update.eq.1) then
     deallocate(pldat,alldat)
     call sleep(5)
     if(sum(ntot).gt.1.e4) call sleep(5)
     if(sum(ntot).gt.1.e5) call sleep(10)
     if(sum(ntot).gt.1.e6) call sleep(20)
     goto 101
  end if
  
  !write(*,'(A)')'  Waiting for you to finish me off...'
  !pause
  
9999 continue
  deallocate(pldat,alldat)
  if(prprogress.ge.2) write(*,*)''
end program plotspins
!************************************************************************************************************************************





