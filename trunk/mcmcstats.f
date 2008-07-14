!mcmcstats.f: Read MCMC statistics output file created by plotspins, and reduce data.

program mcmcstats
  implicit none
  integer, parameter :: nf1=100,nifo1=3,npar1=15,nival1=5
  integer :: i,j,iv,iv1,iv2,f,nf,o,p,io,pgopen,system
  integer :: prinput,plfile
  character :: infile*99,bla,str*99,output*1000
  
  integer :: n(nf1),nburn(nf1),ndet(nf1),seed(nf1),detnr(nf1,nifo1),samplerate(nf1,nifo1),samplesize(nf1,nifo1),FTsize(nf1,nifo1),npar(nf1),ncol(nf1),nival(nf1),tbase(nf1)
  real :: nullh(nf1),snr(nf1,nifo1),totsnr(nf1),flow(nf1,nifo1),fhigh(nf1,nifo1),t_before(nf1,nifo1),t_after(nf1,nifo1),FTstart(nf1,nifo1),deltaFT(nf1,nifo1)
  real :: model(nf1,npar1),median(nf1,npar1),mean(nf1,npar1),stdev1(nf1,npar1),stdev2(nf1,npar1),absvar1(nf1,npar1),absvar2(nf1,npar1),corrs(nf1,npar1,1:npar1)
  real :: ivals(nf1,1:nival1),ivlcntr(nf1,npar1,nival1),ivldelta(nf1,npar1,nival1),ivlinrnge(nf1,npar1,nival1)
  character :: detname(nf1,nifo1)*25,varnames(nf1,npar1)*25,outputnames(nf1)*99,ivlok(nf1,npar1,nival1)*3,pgvarns(1:npar1)*99,pgvarnss(1:npar1)*99
  
  integer :: sym,ci,ci0,ls,ls0,p0,p1,p2,p3,p10,p20,p30,p11,p22
  real :: xmin,xmax,dx,ymin,ymax,dy,x0,y0,x1,y1,clr
  real :: par1,par2,par3,par1s(10),par2s(10),par3s(10)
  
  integer :: rel,nplpar,plpar1,plpar2,plpars(20)
  real :: x,y,pi,d2r
  real :: papsize,paprat
  
  integer :: plotdeltas,plotsnrs,plotcorrelations,plotcorrmatrix,printdeltastable
  
  
  prinput = 0  !Print input to screen: 0-no, 1-yes
  plfile  = 2  !Plot to file: 0-no (screen), 1-png, 2-eps, 3-pdf, 4-eps & pdf
  
  plotdeltas = 0  !0 or 1x
  plotsnrs = 0  !0 or 1
  plotcorrelations = 0  !0 or 1
  plotcorrmatrix = 0  !0-2; 2: swap rows/columns
  printdeltastable = 1  !0 or 1,2
  
  
  if(plfile.eq.0) then
     papsize = 10.81 !Screen size (Gentoo: 10.81, Fink: 16.4)
     paprat = 0.575  !Screen ratio (Gentoo: 0.575, Fink: 0.57)
  else
     papsize = 10.6
     paprat = 0.75
  end if
  
  nf = iargc()
  if(nf.eq.0) then
     write(*,'(/,A,/)')'  Syntax:  plotstats <file1 file2 ...>'
     stop
  end if
  if(nf.gt.nf1) then
     write(6,'(A,I3)')" Too many input files, I'll use the first",nf1
     nf = nf1
  end if
  
  pi = 4*atan(1.)
  d2r = pi/180.
  
  ivlok = ' y '
  
  pgvarns(1:14)  = (/'M\dc\u (M\d\(2281)\u) ','\(2133)               ','t\dc\u (s)            ',  &
       'd\dL\u (Mpc)          ','a\dspin\u             ','\(2134)\dSL\u(\(2218))','R.A. (h)              ','Dec. (\(2218))        ', &
       '\(2147)\dc\u (\(2218))','\(2134)\dJ\u (\(2218))','\(2147)\dJ\u (\(2218))','\(2127)\dc\u (\(2218))     ','M\d1\u (M\d\(2281)\u) ','M\d2\u(M\d\(2281)\u)  '/)
  pgvarnss(1:14)  = (/'M\dc\u','\(2133)','t\dc\u','d\dL\u','a\dspin\u','\(2134)\dSL\u','R.A.','Dec.','\(2147)\dc\u',  &
       '\(2134)\dJ0\u','\(2147)\dJ0\u','\(2127)\dc\u','M\d1\u','M\d2\u'/)
  !Include units
  !pgvarnss(1:14)  = (/'M\dc\u (M\d\(2281)\u)','\(2133)','t\dc\u (s)','d\dL\u (Mpc)','a\dspin\u','\(2134)\dSL\u (\(2218))','R.A. (h)','Dec. (\(2218))','\(2147)\dc\u (\(2218))',  &
  !     '\(2134)\dJ0\u (\(2218))','\(2147)\dJ0\u (\(2218))','\(2127) (\(2218))','M\d1\u (M\d\(2281)\u)','M\d2\u (M\d\(2281)\u)'/)
  
  
  !Read input files
  o = 20
  write(6,*)''
  do f=1,nf
     call getarg(f,infile)
     if(prinput.eq.0) write(6,'(A,$)')' Reading file: '//trim(infile)
     if(prinput.eq.1) write(6,'(A)')' Reading file: '//trim(infile)
     open(unit=o, form='formatted', status='old',file=trim(infile))
     
     read(o,'(A)')outputnames(f)
     if(prinput.eq.1) write(6,*)''
     if(prinput.eq.1) write(6,'(A)')outputnames(f)
     
     read(o,*)bla
     read(o,'(6x,I10,I12,I8,F22.10,I8)')n(f),nburn(f),seed(f),nullh(f),ndet(f)
     if(prinput.eq.1) write(6,*)''
     if(prinput.eq.1) write(6,'(A)')'           niter       nburn    seed       null likelihood    ndet'
     if(prinput.eq.1) write(6,'(6x,I10,I12,I8,F22.10,I8)')n(f),nburn(f),seed(f),nullh(f),ndet(f)
     if(prinput.eq.0) write(6,'(6x,I10,I12,2I8)')n(f),nburn(f),ndet(f),seed(f)
     
     
     read(o,*)bla
     if(prinput.eq.1) write(6,*)''
     if(prinput.eq.1) write(6,'(A)')'        Detector Nr               SNR       f_low      f_high   before tc    after tc    Sample start (GPS)    Sample length   Sample rate   Sample size       FT size'
     totsnr(f) = 0.
     do i=1,ndet(f)
        read(o,'(A16,I3,F18.8,4F12.2,F22.8,F17.7,3I14)')detname(f,i),detnr(f,i),snr(f,i),flow(f,i),fhigh(f,i),t_before(f,i),t_after(f,i),FTstart(f,i),deltaFT(f,i),samplerate(f,i),samplesize(f,i),FTsize(f,i)
        if(prinput.eq.1) write(6,'(A16,I3,F18.8,4F12.2,F22.8,F17.7,3I14)')detname(f,i),detnr(f,i),snr(f,i),flow(f,i),fhigh(f,i),t_before(f,i),t_after(f,i),FTstart(f,i),deltaFT(f,i),samplerate(f,i),samplesize(f,i),FTsize(f,i)
        totsnr(f) = totsnr(f) + snr(f,i)*snr(f,i)
     end do
     totsnr(f) = sqrt(totsnr(f))
     read(o,*)bla,tbase(f)
     if(prinput.eq.1) write(6,*)''
     if(prinput.eq.1) write(6,'(A,I)')' t0: ',tbase(f)

     
     !Read correlations:
     read(o,*)bla,bla,npar(f),ncol(f)
     if(prinput.eq.1) write(6,*)''
     if(prinput.eq.1) write(6,'(A,2I3)')' Npar,ncol: ',npar(f),ncol(f)
     read(o,*)bla !Statistics headers
     if(prinput.eq.1) write(6,'(A)')'  param.       model      median        mean      stdev1      stdev2      abvar1      abvar2'
     do p=1,npar(f)
        read(o,'(A8,7F12.6)')varnames(f,p),model(f,p),median(f,p),mean(f,p),stdev1(f,p),stdev2(f,p),absvar1(f,p),absvar2(f,p)
        if(prinput.eq.1) write(6,'(A8,7F12.6)')trim(varnames(f,p)),model(f,p),median(f,p),mean(f,p),stdev1(f,p),stdev2(f,p),absvar1(f,p),absvar2(f,p)
     end do
     
     
     !Read correlations:
     if(prinput.eq.1) write(6,'(A)')''
     read(o,*)bla
     read(o,*)bla  !Correlation headers
     if(prinput.eq.1) write(6,'(A,2I3)')' Npar: ',npar(f)
     if(prinput.eq.1) write(6,'(A)')'                Mc       eta        tc        dl      spin     th_SL        RA       Dec     phase      thJo      phJo     alpha        M1        M2 '
     do p=1,npar(f)
        read(o,*)varnames(f,p),corrs(f,p,1:npar(f))
        if(prinput.eq.1) write(6,'(A8,20F10.5)')trim(varnames(f,p)),corrs(f,p,1:npar(f))
     end do
     
     !Read intervals:
     if(prinput.eq.1) write(6,'(A)')''
     read(o,*)bla,nival(f)
     if(prinput.eq.1) write(6,'(A,I3)')' Nival: ',nival(f)
     nival(f) = nival(f) + 1  !Since 100% interval is not counted in plotspins
     read(o,*),bla,ivals(f,1:nival(f))
     if(prinput.eq.1) write(6,'(A22,10(F20.5,14x))')'Interval:',ivals(f,1:nival(f))
     
     read(o,*)bla !Interval headers
     if(prinput.eq.1) write(6,'(A)')'  param.        centre       delta in rnge        centre       delta in rnge        centre       delta in rnge        centre       delta in rnge '
     do p=1,npar(f)
        read(o,*)varnames(f,p),(ivlcntr(f,p,iv),ivldelta(f,p,iv),ivlinrnge(f,p,iv),bla,iv=1,nival(f))
        do iv=1,nival(f)
           if(ivlinrnge(f,p,iv).gt.1.) ivlok(f,p,iv) = '*N*'
        end do
        if(prinput.eq.1) write(6,'(A8,2x,5(2F12.6,F6.3,A3,1x))')trim(varnames(f,p)),(ivlcntr(f,p,iv),ivldelta(f,p,iv),ivlinrnge(f,p,iv),ivlok(f,p,iv),iv=1,nival(f))
        iv = nival(f) !100%
        iv = 3 !99%
        if(p.gt.1.and.ivlinrnge(f,p,iv).gt.1.) then !Don't print LogL (p=1)
           write(6,'(A16,A8,5x,2F12.6,5x,2F12.6,F6.3,A3,A3,F5.1,A2)')trim(outputnames(f)),trim(varnames(f,p)), &
                model(f,p),median(f,p),ivlcntr(f,p,iv),ivldelta(f,p,iv),ivlinrnge(f,p,iv),ivlok(f,p,iv),'  (',ivals(f,iv)*100,'%)'
        end if
     end do !p
     if(prinput.eq.1) write(6,*)''
     
     close(o) !Statistics file
  end do !f
  !End reading input files
  
  
  
  
  
  
  
  
  
  
  
  !***********************************************************************************************************************************      
  !Plot deltas
  if(plotdeltas.eq.1.and.nf.gt.2) then
     !write(6,*)''
     write(6,'(A)')' Plotting probability range deltas...'
     if(plfile.eq.0) then
        io = pgopen('12/xs')
     end if
     if(plfile.eq.1) io = pgopen('deltas.ppm/ppm')
     if(plfile.ge.2) io = pgopen('deltas.eps/cps')
     
     if(io.le.0) then
        write(6,'(A,I,/)')'Cannot open PGPlot device.  Quitting the programme ',io
        stop
     end if
     call pgsch(1.5)
     call pgpap(papsize,paprat)
     
     call pgsubp(4,3)
     call pgscr(3,0.,0.5,0.)
     
     iv = 1
     write(6,'(A,I2,A1,F8.3,A1)')' Plotting deltas for probability range',iv,':',ivals(1,iv)*100,'%'
     
     
     do p=1,12
        call pgpage
        call pgsch(2.)
        ci = 1
        ci0 = ci+1
        xmin = 0.
        xmax = 1.
        dx = abs(xmax-xmin)*0.1
        ymin =  1.e30
        ymax = -1.e30
        do f=1,nf
           ymin = min(ymin,ivldelta(f,p,iv))
           ymax = max(ymax,ivldelta(f,p,iv))
           !ymin = min(ymin,log10(ivldelta(f,p,iv)+1.e-30))
           !ymax = max(ymax,log10(ivldelta(f,p,iv)+1.e-30))
        end do
        ymin = 0. !Only in linear
        dy = abs(ymax-ymin)*0.1
        
        
        call pgswin(xmin-dx,xmax+dx,ymin-dy,ymax+dy)
        call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)
        
        call pgsch(3.)
        do p1=1,3        !Ndet, or detnr
           do p2=0,180,5 !Theta_SL
              do p3=0,10 !Spin magnitude * 10
                 par1 = real(p1)
                 par2 = real(p2)
                 par3 = real(p3)*0.1
                 
                 do f=1,nf
                    !if(nint(par1).ne.ndet(f)) cycle       !Multiple numbers of detectors
                    !if(nint(par1).ne.detnr(f,1)) cycle     !Multiple 1-detector cases
                    !if(nint(par2).ne.nint(model(f,6))) cycle
                    !if(nint(par3*10).ne.nint(model(f,5)*10)) cycle
                    
                    !ci = ci+1
                    
                    if(p1.ne.p10) then
                       ci = ci+1
                       if(p.eq.1) write(6,'(A,I3,A,I2)')' ci:',ci,' Ndet:',ndet(f)
                    end if
                    ci = mod(f+1,10)
                    call pgsci(ci)
                    sym = 2
                    if(ivlinrnge(f,p,iv).gt.1.0) sym = 18
                    ls = 1
                    if(p2.eq.55) ls = 2
                    x1 = model(f,5)
                    y1 = ivldelta(f,p,iv)
                    !y1 = log10(ivldelta(f,p,iv)+1.e-30)
                    call pgpoint(1,x1,y1,sym)
                    if(ci.eq.ci0.and.ls.eq.ls0) then
                       call pgsls(ls)
                       !call pgline(2,(/x0,x1/),(/y0,y1/))
                    end if
                    !print*,p,f,model(f,5),ivldelta(f,p,iv)
                    x0  = x1
                    y0  = y1
                    ci0 = ci
                    ls0 = ls
                    p10 = p1
                    p20 = p2
                 end do
              end do
           end do
        end do
        
        call pgsci(1)
        call pgsls(1)
        call pgsch(2.)
        !call pgmtxt('T',1.,0.5,0.5,trim(varnames(1,p)) )
        call pgmtxt('T',1.,0.5,0.5,trim(pgvarns(p)) )
        !write(6,*)''
     end do
     
     
     call pgend
     
     if(plfile.eq.1) then
        !i = system('convert -depth 8 deltas.ppm  '//trim(outputname)//'__deltas.png')
        i = system('convert -depth 8 deltas.ppm deltas.png')
        i = system('rm -f deltas.ppm')
     end if
     if(plfile.gt.2) then
        !i = system('eps2pdf deltas.eps  -o '//trim(outputname)//'__deltas.pdf   >& /dev/null')
        !i = system('mv -f deltas.eps '//trim(outputname)//'__deltas.eps')
        i = system('eps2pdf deltas.eps  -o deltas.pdf   >& /dev/null')
        if(plfile.eq.3) i = system('rm -f detas.eps')
     end if
  end if  !if(plotdeltas.eq.1.and.nf.gt.2) then
  
  
  
  
  
  
  
  
  
  !***********************************************************************************************************************************      
  !Plot SNRs
  if(plotsnrs.eq.1.and.nf.gt.2) then
     write(6,*)''
     write(6,'(A)')' Plotting SNRs...'
     if(plfile.eq.0) then
        io = pgopen('13/xs')
     end if
     if(plfile.eq.1) io = pgopen('snrs.ppm/ppm')
     if(plfile.ge.2) io = pgopen('snrs.eps/cps')
     
     if(io.le.0) then
        write(6,'(A,I,/)')'Cannot open PGPlot device.  Quitting the programme ',io
        stop
     end if
     call pgsch(1.5)
     call pgpap(papsize,paprat)
     
     !call pgsubp(4,3)
     call pgscr(3,0.,0.5,0.)
     
     
     !call pgpage
     call pgsch(2.)
     ci = 1
     ci0 = ci+1
     xmin = 0.
     xmax = 1.
     dx = abs(xmax-xmin)*0.1
     ymin =  1.e30
     ymax = -1.e30
     do f=1,nf
        ymin = min(ymin,totsnr(f))
        ymax = max(ymax,totsnr(f))
     end do
     dy = abs(ymax-ymin)*0.1
     
     
     call pgswin(xmin-dx,xmax+dx,ymin-dy,ymax+dy)
     call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)
     
     call pgsch(3.)
     do p1=1,3        !Ndet
        do p2=0,180,5 !Theta_SL
           do p3=0,10 !Spin magnitude * 10
              par1 = real(p1)
              par2 = real(p2)
              par3 = real(p3)*0.1
              
              do f=1,nf
                 if(nint(par1).ne.ndet(f)) cycle       !Multiple numbers of detectors
                 !if(nint(par1).ne.detnr(f,1)) cycle     !Multiple 1-detector cases
                 if(nint(par2).ne.nint(model(f,6))) cycle
                 if(nint(par3*10).ne.nint(model(f,5)*10)) cycle
                 
                 !ci = ci+1
                 
                 if(p1.ne.p10) then
                    ci = ci+1
                    if(p.eq.1) write(6,'(A,I3,A,I2)')' ci:',ci,' Ndet:',ndet(f)
                 end if
                 
                 call pgsci(ci)
                 sym = 2
                 !if(ivlinrnge(f,p,iv).gt.1.0) sym = 18
                 ls = 1
                 if(p2.eq.55) ls = 2
                 x1 = model(f,5)
                 y1 = totsnr(f)
                 call pgpoint(1,x1,y1,sym)
                 if(ci.eq.ci0.and.ls.eq.ls0) then
                    call pgsls(ls)
                    call pgline(2,(/x0,x1/),(/y0,y1/))
                 end if
                 !print*,p,f,model(f,5),ivldelta(f,p,iv)
                 x0  = x1
                 y0  = y1
                 ci0 = ci
                 ls0 = ls
                 p10 = p1
                 p20 = p2
              end do
           end do
        end do
     end do
     
     call pgsci(1)
     call pgsls(1)
     call pgsch(2.)
     call pgmtxt('T',1.,0.5,0.5,'SNR' )
     !write(6,*)''
     
     
     call pgend
     
     if(plfile.eq.1) then
        !i = system('convert -depth 8 snrs.ppm  '//trim(outputname)//'__snrs.png')
        i = system('convert -depth 8 snrs.ppm snrs.png')
        i = system('rm -f snrs.ppm')
     end if
     if(plfile.gt.2) then
        !i = system('eps2pdf snrs.eps  -o '//trim(outputname)//'__snrs.pdf   >& /dev/null')
        !i = system('mv -f snrs.eps '//trim(outputname)//'__snrs.eps')
        i = system('eps2pdf snrs.eps  -o snrs.pdf   >& /dev/null')
        if(plfile.eq.3) i = system('rm -f snrs.eps')
     end if
  end if  !if(plotsnrs.eq.1.and.nf.gt.2) then
  
  
  
  
  
  
  
  
  
  !***********************************************************************************************************************************      
  !Plot correlations
  if(plotcorrelations.eq.1.and.nf.gt.2) then
     write(6,*)''
     write(6,'(A)')' Plotting correlations...'
     
     do p0 = 1,12
        write(6,'(A)')' Plotting correlations with '//trim(varnames(1,p0))
        
        if(plfile.eq.0) then
           io = pgopen('12/xs')
        end if
        if(plfile.eq.1) io = pgopen('corrs.ppm/ppm')
        if(plfile.ge.2) io = pgopen('corrs.eps/cps')
        
        if(io.le.0) then
           write(6,'(A,I,/)')'Cannot open PGPlot device.  Quitting the programme ',io
           stop
        end if
        call pgsch(1.5)
        call pgpap(papsize,paprat)

        call pgsubp(4,3)
        call pgscr(3,0.,0.5,0.)
        
        
        do p=1,12
           call pgpage
           call pgsch(2.)
           ci = 1
           ci0 = ci+1
           xmin = 0.
           xmax = 1.
           dx = abs(xmax-xmin)*0.1
           ymin = -1.
           !ymin = 0. !If abs()
           ymax = 1.
           dy = abs(ymax-ymin)*0.1
           !ymin = ymin+dy !Force it to be zero, if abs()
           
           
           call pgswin(xmin-dx,xmax+dx,ymin-dy,ymax+dy)
           call pgbox('ABCNTS',0.0,0,'BCNTS',0.0,0)
           
           call pgsch(3.)
           do p1=1,3        !Ndet
              do p2=0,180,5 !Theta_SL
                 do p3=0,10 !Spin magnitude * 10
                    par1 = real(p1)
                    par2 = real(p2)
                    par3 = real(p3)*0.1
                    
                    do f=1,nf
                       !if(nint(par1).ne.ndet(f)) cycle       !Multiple numbers of detectors
                       !if(nint(par1).ne.detnr(f,1)) cycle     !Multiple 1-detector cases
                       !if(nint(par2).ne.nint(model(f,6))) cycle
                       !if(nint(par3*10).ne.nint(model(f,5)*10)) cycle
                       
                       !ci = ci+1
                       
                       if(p1.ne.p10) then
                          ci = ci+1
                          if(p.eq.1) write(6,'(A,I3,A,I2)')' ci:',ci,' Ndet:',ndet(f)
                       end if
                       ci = mod(f+1,10)
                       call pgsci(ci)
                       sym = 2
                       if(ivlinrnge(f,p,iv).gt.1.0) sym = 18
                       ls = 1
                       if(p2.eq.55) ls = 2
                       x1 = model(f,5)
                       y1 = corrs(f,p0,p)
                       !y1 = abs(corrs(f,p0,p))
                       !call pgpoint(1,x1,y1,sym)
                       call pgpoint(1,x1,y1,2)
                       if(ci.eq.ci0.and.ls.eq.ls0) then
                          call pgsls(ls)
                          !call pgline(2,(/x0,x1/),(/y0,y1/))
                       end if
                       x0  = x1
                       y0  = y1
                       ci0 = ci
                       ls0 = ls
                       p10 = p1
                       p20 = p2
                    end do
                 end do
              end do
           end do
           
           call pgsci(1)
           call pgsls(1)
           call pgsch(2.)
           !call pgmtxt('T',1.,0.5,0.5,trim(varnames(1,p)) )
           call pgmtxt('T',1.,0.5,0.5,trim(pgvarns(p)) )
           !write(6,*)''
        end do
        
        
        call pgend
        
        if(plfile.eq.1) then
           !i = system('convert -depth 8 corrs.ppm  '//trim(outputname)//'__corrs.png')
           i = system('convert -depth 8 corrs.ppm corrs.png')
           i = system('rm -f corrs.ppm')
        end if
        if(plfile.ge.2) then
           !i = system('eps2pdf corrs.eps  -o '//trim(outputname)//'__corrs.pdf   >& /dev/null')
           !i = system('mv -f corrs.eps '//trim(outputname)//'__corrs.eps')
           write(str,'(I2.2)')p0
           if(plfile.gt.2) i = system('eps2pdf corrs.eps  -o corrs_'//trim(str)//'.pdf   >& /dev/null')
           if(plfile.eq.2.or.plfile.eq.4) i = system('mv -f corrs.eps  corrs_'//trim(str)//'.eps')
           if(plfile.eq.3) i = system('rm -f corrs.eps')
        end if
     end do !p0
  end if !  if(plotcorrelations.eq.1.and.nf.gt.2) then
  
  
  
  
  
  
  
  !***********************************************************************************************************************************      
  !Plot correlation matrix
  if(plotcorrmatrix.eq.1.and.nf.le.2) then
     write(6,*)''
     write(6,'(A)')' Plotting correlation matrix...'
     
     if(plfile.eq.0) then
        io = pgopen('12/xs')
     end if
     if(plfile.eq.1) io = pgopen('corr_matrix.ppm/ppm')
     if(plfile.ge.2) io = pgopen('corr_matrix.eps/cps')
     
     if(io.le.0) then
        write(6,'(A,I,/)')'Cannot open PGPlot device.  Quitting the programme ',io
        stop
     end if
     call pgpap(papsize,paprat)
     
     plpars = 0
     plpar1 = 1
     plpar2 = 12
     plpars(plpar1:plpar2) = (/1,1,1,1,1,1,1,1,1,1,1,1/)
     plpars(plpar1:plpar2) = (/1,1,0,0,1,1,1,1,0,0,0,0/)
     nplpar = sum(plpars)
     !nplpar = plpar2 - plpar1 + 1
     
     
     !call pgsubp(4,3)
     call pgscr(3,0.,0.5,0.)
     !call pgswin(-1.,13.,13.,-1.)
     call pgswin(-1.,real(nplpar+1),real(nplpar+1),-1.)
     call pgslw(3)
     call pgsch(sqrt(12./real(nplpar)))
     call pgscf(2)
     !if(plfile.ge.2) then
     !   call pgsch(1.)
     !end if
     
     ci = 20
     
     call pgsci(1)
     
     dx = 6./real(nplpar)
     p11 = 0
     do p=plpar1,plpar2
        if(plpars(p).eq.1) then
           p11 = p11+1
        else
           cycle
        end if
        !call pgptxt(real(p)-0.5,-0.5,0.,0.5,trim(pgvarnss(p)))
        !if(plpar2.eq.12) call pgptxt(real(p)-0.5,12.5,0.,0.5,trim(pgvarnss(p)))
        !call pgptxt(-0.5,real(p)-0.3,0.,0.5,trim(pgvarnss(p)))
        !if(plpar2.eq.12) call pgptxt(12.5,real(p)-0.3,0.,0.5,trim(pgvarnss(p)))
        
        call pgptxt(real(p11)-0.5,-dx*0.5-0.0167*nplpar,0.,0.5,trim(pgvarnss(p))) !At top
        call pgptxt(-dx*0.5-0.0167*nplpar,real(p11)-0.5+0.0167*nplpar,0.,0.5,trim(pgvarnss(p)))  !At left
        !if(plpar2.eq.12) call pgptxt(real(p11)-dx,real(nplpar)+dx,0.,0.5,trim(pgvarnss(p)))
        !if(plpar2.eq.12) call pgptxt(real(nplpar)+dx,real(p11)-dx*0.6,0.,0.5,trim(pgvarnss(p)))
     end do
     
     do f=1,nf
        p11 = 0
        do p1=plpar1,plpar2
           if(plpars(p1).eq.0) then
              cycle
           else
              p11 = p11+1
           end if
           
           p22 = 0
           do p2=plpar1,plpar2
              if(plpars(p2).eq.0) then
                 cycle
              else
                 p22 = p22+1
              end if
              
              if(f.eq.1.and.p2.ge.p1) cycle  !Use upper triangle
              if(f.eq.2.and.p2.le.p1) cycle  !Use lower triangle
              !if(p1.gt.2.and.p1.lt.5.or.p1.gt.6) cycle !Just do masses and spins
              !if(p2.gt.2.and.p2.lt.5.or.p2.gt.6) cycle
              
              y1 = corrs(f,p2,p1)
              if(f.eq.2) y1 = corrs(f,p1,p2)
              !if(abs(y1).lt.0.8) cycle
              
              !call pgscr(ci,abs(y1),abs(y1),abs(y1)) !Grey
              if(plfile.lt.2) then                          !Screen/bitmap have black backgrounds
                 if(f.eq.1) call pgscr(ci,0.,0.,abs(y1)**2) !Blue
                 if(f.eq.2) call pgscr(ci,abs(y1)**2,0.,0.) !Red
              else                                          !eps/pdf have white backgrounds
                 if(f.eq.1) call pgscr(ci,1.-abs(y1)**2,1.-abs(y1)**2,1.) !Blue
                 if(f.eq.2) call pgscr(ci,1.,1.-abs(y1)**2,1.-abs(y1)**2) !Red 
              end if
              call pgsci(ci)
              !call pgrect(real(p1-1),real(p1),real(p2-1),real(p2))
              call pgrect(real(p11-1),real(p11),real(p22-1),real(p22))
              
              !if(plfile.lt.2) then
              !   call pgsci(0)
              !   if(abs(y1).lt.0.5) call pgsci(1) !For grey values only
              !else
              !   call pgsci(1)
              !   if(abs(y1).lt.0.5) call pgsci(0) !For grey values only
              !end if
              
              !if(abs(y1).lt.0.1) cycle
              write(str,'(F5.2)')y1
              clr = 0.
              if(plfile.lt.2) clr = 1.
              if(plfile.lt.2.and.abs(y1).lt.0.25) clr = abs(y1)*4.      !Fade out the numbers
              if(plfile.ge.2.and.abs(y1).lt.0.25) clr = 1. - abs(y1)*4.
              call pgscr(ci,clr,clr,clr)
              !call pgptxt(real(p1)-0.5,real(p2)-0.3,0.,0.5,trim(str))
              call pgptxt(real(p11)-0.5,real(p22)-0.5+0.0167*nplpar,0.,0.5,trim(str))
           end do
        end do
     end do
     
     call pgsci(1)
     call pgslw(1)
     !call pgsls(1)
     !call pgsch(2.)
     !!call pgmtxt('T',1.,0.5,0.5,trim(varnames(1,p)) )
     !call pgmtxt('T',1.,0.5,0.5,trim(pgvarns(p)) )
     !!write(6,*)''
     
     !call pgscr(14,0.7,0.7,0.7)
     !call pgsci(14)
     call pgslw(5)
     call pgline(2,(/0.,real(nplpar)/),(/0.,real(nplpar)/))
     call pgline(2,(/0.,0./),(/0.,real(nplpar)/))
     call pgline(2,(/0.,real(nplpar)/),(/0.,0.0/))
     call pgline(2,(/0.0,real(nplpar)/),(/real(nplpar),real(nplpar)/))
     call pgline(2,(/real(nplpar),real(nplpar)/),(/0.,real(nplpar)/))
     do p = 1,nplpar
        !call pgrect(real(p-1),real(p),real(p-1),real(p))
        !call pgcirc(real(p-0.5),real(p-0.5),0.2)
     end do
     call pgend
     
     if(plfile.eq.1) then
        !i = system('convert -depth 8 corr_matrix.ppm  '//trim(outputname)//'__corr_matrix.png')
        i = system('convert -depth 8 corr_matrix.ppm corr_matrix.png')
        i = system('rm -f corr_matrix.ppm')
     end if
     if(plfile.gt.2) then
        i = system('eps2pdf corr_matrix.eps  -o corr_matrix.pdf   >& /dev/null')
        if(plfile.eq.3) i = system('rm -f corr_matrix.eps')
     end if
  end if !if(plotcorrmatrix.eq.1.and.nf.le.2) then
  
  
  
  
  
  
  
  !*******************************************************************************************
  ! Create a table of Delta's
  !*******************************************************************************************
  
  if(printdeltastable.eq.1) then
     !open(unit=30, form='formatted', status='replace', file='table.tex')
     open(unit=30, form='formatted', status='unknown', position='append', file='table.tex')
     write(6,*)''
     do f=1,nf
        iv = 0
        do i=1,nival(f)
           if(abs(ivals(f,i)-0.90).lt.1.e-5) iv = i
        end do
        if(iv.eq.0) then
           write(6,'(A)')'  Interval not found, file '//trim(outputnames(f))
        else
           !Write the number of detectors, a_spin and theta_sl
           !write(output,'(I3,A,I6,A,F6.1,A)')ndet(f),'  &  ',nint(model(f,6)),'$^\circ$  &  ',model(f,5),'  &  '
           
           !Write the number of detectors, a_spin, theta_sl, and SNR:
           !write(output,'(I3,A,F6.1,A,I6,A,F6.1,A)')ndet(f),'  $\!\!\!\!$ &  ',model(f,6),'  $\!\!\!\!$ &  ',nint(model(f,7)),'  $\!\!\!\!$ &  ',totsnr(f),'  $\!\!\!\!$ &  ' 
           
           !Write the number of detectors, a_spin, theta_sl, distance and SNR:
           !write(output,'(I3,A,F6.1,A,I6,A,F6.1,A,F6.1,A)')ndet(f),'  $\!\!\!\!$ &  ',model(f,6),'  $\!\!\!\!$ &  ',nint(model(f,7)),'  $\!\!\!\!$ &  ',model(f,5),'  $\!\!\!\!$ &  ',totsnr(f),'  $\!\!\!\!$ &  ' 
           
           !Write the number of detectors and SNR:
           write(output,'(I3,A,F6.1,A)')ndet(f),'  $\!\!\!\!$ &  ',totsnr(f),'  $\!\!\!\!$ &  ' 
           
           !do p1=2,npar(f)  !Leave out logL
           do p1=4,npar(f)  !Leave out logL, M1, M2
              !print*,p,npar(f)
              
              !Display M1,M2 iso Mc,eta
              !p = p1
              !if(p.gt.12) cycle
              !if(1.eq.2) then !replace Mc,eta by M1,M2
              !   if(p1.eq.1) p=13
              !   if(p1.eq.2) p=14
              !end if
              
              !Display M1,M2 iso Mc,eta
              p = p1-2
              !if(p.gt.13) cycle
              if(1.eq.1) then !replace Mc,eta by M1,M2
                 !if(p1.eq.1) p=14
                 !if(p1.eq.2) p=15
                 if(p.eq.0) p=14
                 if(p.eq.1) p=15
              end if
              
              
              iv1 = 0
              iv2 = 1
              !do iv=iv1,iv2 !90 and 95%
              do iv=iv2,iv2 !90% only
                 x = ivldelta(f,p,iv)
                 rel = 0
                 if(p.eq.2.or.p.eq.3.or.p.eq.5.or.p.eq.6.or.p.eq.14.or.p.eq.15) rel=1  !Use relative deltas (%)
                 if(p.eq.4) x = x*1000                            !Time in ms
                 !if(p.eq.8) x = x*cos(40*d2r)*15                 !Convert RA to RA*cos(decl) and from hrs to deg
                 !if(p.eq.8) x = x*15                              !Convert RA from hrs to deg
                 if(p.eq.8) x = x*15 * cos(model(f,p)*d2r)            !Convert RA from hrs to deg and take into account decl
                 if(rel.eq.1) x = x/ivlcntr(f,p,iv)*100
                 if(x.gt.10.) then
                    write(output,'(A,I6)')trim(output),nint(x)
                 else if(x.gt.1.) then
                    write(output,'(A,F6.1)')trim(output),x
                 else  if(x.gt.0.1) then
                    write(output,'(A,F6.2)')trim(output),x
                 else
                    write(output,'(A,F6.3)')trim(output),x
                 end if
                 !print*,p1,p,x,ivldelta(f,p,iv),ivlcntr(f,p,iv)
                 
                 !if(iv.eq.iv2.and.rel.eq.1) write(output,'(A)')trim(output)//'\%'
                 if(iv.eq.iv1) write(output,'(A)')trim(output)//';'
                 !if(ivlinrnge(f,p,iv).gt.1.) write(output,'(A)')trim(output)//'*'
                 if(ivlinrnge(f,p,1).gt.1.) write(output,'(A)')trim(output)//'*'
                 if(ivlinrnge(f,p,2).gt.1.) write(output,'(A)')trim(output)//'*'
                 if(ivlinrnge(f,p,3).gt.1.) write(output,'(A)')trim(output)//'*'
                 if(ivlinrnge(f,p,4).gt.1.) write(output,'(A)')trim(output)//'*'
                 
                 !if(iv.eq.iv2.and.p.lt.npar(f)) write(output,'(A)')trim(output)//'  $\!\!\!\!$  &' !Add latex codes
                 !if(iv.eq.iv2.and.p1.lt.12) write(output,'(A)')trim(output)//'  $\!\!\!\!$  &' !Add latex codes (p1.lt.12 iso p.lt.npar, in case npar=14)
                 !if(iv.eq.iv2.and.p1.ne.npar(f)) write(output,'(A)')trim(output)//'  $\!\!\!\!$  &' !Add latex codes (p1.lt.12 iso p.lt.npar, in case npar=14)
                 !if(iv.eq.iv2.and.p1.ne.npar(f)) write(output,'(A)')trim(output)//'  &' !Add latex codes (p1.lt.12 iso p.lt.npar, in case npar=14)
                 if(iv.eq.iv2.and.p1.lt.npar(f)) write(output,'(A)')trim(output)//'   $\!\!\!$ &' !Add latex codes
              end do !iv
           end do !p1
           !write(6,*)''
           !if(f.ne.nf) write(output,'(A)')trim(output)//'  \\' !Add latex codes
           write(output,'(A)')trim(output)//'  \\' !Add latex codes
           
           !Remove ', ' from output
           do j=1,10
              do i=1,len_trim(output)-1
                 if(output(i:i+1).eq.'; ') write(output(i+1:len_trim(output)),'(A)')output(i+2:len_trim(output)+1)
              end do
           end do
           
           write(30,'(A)')trim(output)
        end if
     end do !f
     close(30)
  end if !if(printdeltastable.eq.1)
  
  
  !Swap columns and rows
  if(printdeltastable.eq.2) then
     write(6,*)''
     do p=1,npar(1)
        do f=1,nf
           iv = 0
           do i=1,nival(f)
              if(abs(ivals(f,i)-0.90).lt.1.e-5) iv = i
           end do
           if(iv.eq.0) then
              write(6,'(A)')'  Interval not found, file '//trim(outputnames(f))
           else
              !write(6,'(I3,2F10.2,$)')ndet(f),model(f,5),model(f,6) !Write the number of detectors, a_spin and theta_sl
              !write(6,'(I3,A5,F6.1,A5,I6,A5,$)')ndet(f),'  &  ',model(f,5),'  &  ',nint(model(f,6)),'  &  ' !Write the number of detectors, a_spin and theta_sl
              x = ivldelta(f,p,iv)
              rel = 0
              if(p.eq.1.or.p.eq.2.or.p.eq.4.or.p.eq.5) rel=1  !Use relative deltas (%)
              if(rel.eq.1) x = x/ivlcntr(f,p,iv)*100
              if(x.gt.10.) then
                 write(output,'(I6)')nint(x)
              else if(x.gt.1.) then
                 write(output,'(A,F6.1)')trim(output),x
              else  if(x.gt.0.1) then
                 write(output,'(A,F6.2)')trim(output),x
              else
                 write(output,'(A,F6.3)')trim(output),x
              end if
              
              if(rel.eq.1) then
                 write(output,'(A)')trim(output)//'%'
              else
                 write(output,'(A)')trim(output)//' '
              end if
              
              write(output,'(A)')trim(output)//'  &  ' !Add latex codes
              
           end if
        end do
        write(output,'(A)')trim(output)//'  \\' !Add latex codes
        
        write(6,'(A)')trim(output)
     end do
  end if !if(printdeltastable.eq.2)
  
  
  
  
  
  
  
  
  write(6,*)''
end program mcmcstats
