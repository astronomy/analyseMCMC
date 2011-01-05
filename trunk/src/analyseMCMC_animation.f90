!> \file analyseMCMC_animation.f90  Create animations for analyseMCMC


!***********************************************************************************************************************************
!> \brief  Create an animation from MCMC output
!!
!! \retval exitcode  Exit status code (0=ok)

subroutine animation(exitcode)
  use basic
  use constants
  use analysemcmc_settings
  use general_data
  use mcmcrun_data
  use chain_data
  use plot_data
  
  implicit none
  integer :: c,i,ic,io,p,iframe,nplt,pgopen,lw,n1,n2,exitcode,status,system
  integer :: index(maxMCMCpar,maxChs*maxIter),small_anim
  real :: range,range1,range2,drange,minrange,centre,median,plshift,ival,norm
  real :: x(maxChs,maxChs*maxIter),x1,x2,xmin,xmax,xmin1,xmax1,dx,y1,y2,ymin,ymax,dy,sch
  real,allocatable :: xbin(:,:),ybin(:,:),xbin1(:),ybin1(:)    !These depend on Nbin1D, allocate after reading input file
  real(double) :: ts1,ts2,timestamp
  character :: framename*(99),tms*(8),str*(99)
  
  exitcode = 0
  
  ts1 = timestamp()
  write(stdOut,*)
  !p = 1 !Parameter to plot: 1-Mc
  p = plAnim
  small_anim = 1  !Make small animation (e.g. gif) rather than ~screen size
  
  !Autodetermine number of bins:
  if(Nbin1D.le.0) then
     call determine_nbin_1d(totpts,Nbin1D)
     if(prProgress.ge.2.and.plot.eq.1.and.update.eq.0) then
        if(Nbin1D.lt.100) write(stdOut,'(A2,I2,A8)',advance="no")' (',Nbin1D,' bins), '
        if(Nbin1D.ge.100) write(stdOut,'(A2,I3,A8)',advance="no")' (',Nbin1D,' bins), '
     end if
  else
     Nbin1D = max(Nbin1D,5)
     if(prProgress.ge.1.and.plot.eq.1.and.update.eq.0) write(stdOut,'(A2)',advance="no")', '
  end if
  
  
  !animScheme = 1  !Left column: Chain, sigma and acceptance, right column: numbers and PDF
  !animScheme = 2  !Upper panel: numbers and PDF, lower panel: chain
  !animScheme = 3  !Upper panel: numbers and PDF, lower panel: chain and logL
  
  !Allocate memory:
  allocate(xbin(maxChs,Nbin1D+1),ybin(maxChs,Nbin1D+1),xbin1(Nbin1D+1),ybin1(Nbin1D+1))
  
  !nAnimFrames = nAnimFrames-1
  if(prProgress.ge.1.and.update.eq.0) write(stdOut,'(A,I5,A,I6,A,/)')'  Creating animation using',nAnimFrames,' frames and', &
       maxval(ntot(1:nchains)),' points..'
  do iframe = 0,nAnimFrames
     nplt = nint(real(iframe)/real(nAnimFrames)*maxval(ntot(1:nchains)))  !This is the line number, not the iteration number
     
     if(prProgress.ge.1.and.update.eq.0) then
        write(stdOut,*)upline !Move cursor up 1 line
        write(stdOut,'(A,I5,A1,I5,A,I6,A,I6,A)',advance="no")'  Plotting movie frame',iframe,'/',nAnimFrames, &
             '  (',nplt,'/',ntot(1:nchains),' points)'
        
        !Print remaining time
        ts2 = timestamp()
        !Use the system clock:
        if(prProgress.ge.1.and.file.eq.1) write(stdOut,'(A,A9)')'   Est.time left:',tms((ts2-ts1)*(nAnimFrames-iframe)/3600.d0)
        ts1 = ts2
     end if
     
     write(framename,'(A,I4.4,A4)')'frame_',iframe,'.ppm'
     if(small_anim.eq.1) write(framename,'(A,I4.4,A4)')'frame_',iframe,'.png'
     
     ic = 1
     
     !print*,iframe,nplt,Nburn(ic),maxval(is(ic,:))
     
     !Determine interval ranges
     c = c0
     ival = ivals(c)
     
     range1 = 0.
     range2 = 0.
        
     !Determine median and ranges for this selection of the chain (Nburn:nplt)
     if(nplt.gt.Nburn(ic)) then
        
        !Determine the median
        !call rindexx(nplt-Nburn(ic),allDat(ic,p,1:nplt-Nburn(ic)),index(p,1:nplt-Nburn(ic)))  !Sort
        x(ic,1:nplt-Nburn(ic)) = allDat(ic,p,Nburn(ic)+1:nplt)
        call rindexx(nplt-Nburn(ic),x(ic,1:nplt-Nburn(ic)),index(p,1:nplt-Nburn(ic)))  !Sort
        !print*,x(ic,index(p,1:nplt-Nburn(ic)))
        
        if(mod(nplt-Nburn(ic),2).eq.0) then
           !Centre = nb + (n-nb)/2 = (n+nb)/2:
           !median = 0.5*(allDat(ic,p,index(p,(nplt-Nburn(ic))/2)) + allDat(ic,p,index(p,(nplt-Nburn(ic))/2+1)))  
           median = 0.5*(x(ic,index(p,(nplt-Nburn(ic))/2)) + x(ic,index(p,(nplt-Nburn(ic))/2+1)))
        else
           !median = allDat(ic,p,index(p,(nplt-Nburn(ic)+1)/2))
           median = x(ic,index(p,(nplt-Nburn(ic)+1)/2))
        end if
        !print*,'median:',median
        
        
        minrange = huge(minrange)
        y1 = 0.
        y2 = 0.
        !do i=1,floor(nplt*(1.-ival))
        do i=1,floor((nplt-Nburn(ic))*(1.-ival))
           !x1 = allDat(ic,p,index(p,i))
           !x2 = allDat(ic,p,index(p,i+floor((nplt-Nburn(ic))*ival)))
           x1 = x(ic,index(p,i))
           x2 = x(ic,index(p,i+floor((nplt-Nburn(ic))*ival)))
           range = abs(x2 - x1)
           if(range.lt.minrange) then
              minrange = range
              y1 = x1
              y2 = x2
           end if
           !write(stdOut,'(2I6,7F12.8)')i,i+floor((nplt-Nburn(ic))*ival),x1,x2,range,minrange,y1,y2,(y1+y2)/2.
        end do !i
        centre = (y1+y2)/2.
        !write(stdOut,'(A8,4x,4F10.5,I4)')parNames(p),y1,y2,minrange,centre,wrap(ic,p)
        
        !Save ranges:
        range1 = y1
        range2 = y2
        drange = y2-y1
        if(parID(p).eq.21.or.parID(p).eq.22.or.parID(p).eq.61.or.parID(p).eq.62.or.parID(p).eq.63.or.parID(p).eq.64 &
             .or.parID(p).eq.71.or.parID(p).eq.81) drange = drange/((y1+y2)/2.)
        
        !print*,'ranges:',range1,range2,drange
     end if
     
     
     
     
     
     if(file.eq.0) io = pgopen('17/xs')
     if(file.eq.1) io = pgopen('analysemcmc_frame.ppm/ppm')
     if(io.le.0) then
        write(stdErr,'(A,I4)')'   ERROR:  Cannot open PGPlot device.  Quitting the programme',io
        exitcode = 1
        return
     end if
     !if(file.eq.0) call pgpap(12.,0.72)
     if(file.eq.0) call pgpap(scrSz*0.75,scrRat)
     
     if(small_anim.eq.1) then
        !for ~500x350, change convert below.  Make the output image twice as big, rescale in the end:
        if(file.eq.1) call pgpap(12.,0.75)   
     else
        !for 850x612, change convert below.  1:1.388 for mpeg, make the output image twice as big, rescale in the end:
        !if(file.eq.1) call pgpap(20.,0.72)  
        
        !for 1024x738, change convert below. 1:1.388 for mpeg, make the output image twice as big, rescale in the end:
        if(file.ge.1) call pgpap(24.08,0.72) 
     end if
     
     call pgsch(1.)
     !call pgscr(3,0.,0.5,0.)
     call pginitl(colour,file,whiteBG)
     lw = 2
     if(file.ge.1) lw = 5
     if(file.ge.1.and.small_anim.eq.1) lw = 2
     
     
     
     
     
     !******************************************************************************************************************************
     !Plot chain for this parameter
     
     
     !if(prProgress.ge.1.and.update.eq.0) write(stdOut,'(A)')'  - parameter chains'
     call pgslw(lw)
     if(animScheme.eq.1) call pgsvp(0.05,0.35,0.65,0.95)
     if(animScheme.eq.2) call pgsvp(0.05,0.95,0.05,0.25)
     if(animScheme.eq.3) call pgsvp(0.20,0.95,0.20,0.35)
     ic = 1
     xmin = 0.
     !xmax = real(maxval(ntot(1:nchains0)))
     xmax = maxval(is)
     dx = abs(xmax-xmin)*0.01
     if(animScheme.eq.2) dx = 0.
     !if(animScheme.eq.3) dx = abs(xmax-xmin)*0.05
     ymin =  1.e30
     ymax = -1.e30
     do ic=1,nchains0
        ymin = min(ymin,minval(allDat(ic,p,10:ntot(ic))))
        ymax = max(ymax,maxval(allDat(ic,p,10:ntot(ic))))
        !if(animScheme.eq.3) then
        !   ymin = min(ymin,minval(allDat(ic,p,1:ntot(ic))))
        !   ymax = max(ymax,maxval(allDat(ic,p,1:ntot(ic))))
        !end if
        if(plInject.ge.1) then
           ymin = minval((/startval(ic,p,1),ymin/))
           ymax = maxval((/startval(ic,p,1),ymax/))
        end if
        if(plStart.ge.1) then
           ymin = minval((/startval(ic,p,2),ymin/))
           ymax = maxval((/startval(ic,p,2),ymax/))
        end if
     end do
     dy = abs(ymax-ymin)*0.05
     !if(animScheme.eq.3) dy = abs(ymax-ymin)*0.01
     !if(offsetrun.eq.1) dy = abs(ymax-ymin)*0.02
     !if(offsetrun.eq.1) dy = -1*dy
     
     call pgswin(xmin-dx,xmax+dx,ymin-dy,ymax+dy)
     !call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)
     if(animScheme.eq.1) call pgbox('BCTS',0.0,0,'BCNTS',0.0,0)
     if(animScheme.eq.2) call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)
     !if(animScheme.eq.3) call pgbox('BCTS',0.0,0,'BCNTS',0.0,0)
     if(animScheme.eq.3) call pgbox('BCTS',0.0,0,'BCNTS',10.**floor(log10(abs(ymax-ymin))),2)
     
     !call pgslw(1)
     do ic=1,nchains0
        !call pgsci(mod(ic*2,10))
        call pgsci(5)
        if(nchains0.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
        !do i=1,ntot(ic),chainPlI
        if(chainSymbol.eq.0) then !Plot lines rather than symbols
           call pgline(nplt-1,is(ic,2:nplt),allDat(ic,p,2:nplt))
        else
           call pgpoint(1,0.,startval(ic,p,2),chainSymbol)       !Starting value
           do i=1,nplt,chainPlI
              call pgpoint(1,is(ic,i),allDat(ic,p,i),chainSymbol)
           end do
        end if
     end do
     !call pgslw(lw)
     
     do ic=1,nchains0
        call pgsci(2)
        call pgsls(2)
        if(nchains0.gt.1) call pgsci(1)
        call pgline(2,(/-1.e20,1.e20/),(/startval(ic,p,1),startval(ic,p,1)/)) !Injection value
        call pgsci(6)
        if(nchains0.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
        !if(plBurn.ge.1.and.isburn(ic).lt.is(ic,ntot(ic))) call pgline(2,(/isburn(ic),isburn(ic)/),(/-1.e20,1.e20/))
        if(plBurn.ge.1.and.isburn(ic).lt.is(ic,max(min(nplt,ntot(ic)),1) ))  &
             call pgline(2,(/isburn(ic),isburn(ic)/),(/-1.e20,1.e20/))
        call pgsci(2)
        call pgsls(4)
        if(nchains0.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
        call pgline(2,(/-1.e20,1.e20/),(/startval(ic,p,2),startval(ic,p,2)/))  !Starting value
     end do
     call pgsci(1)
     call pgsls(1)
     if(animScheme.eq.3) then
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
     
     
     
     
     !******************************************************************************************************************************
     !Plot log(L)
     
     if(animScheme.eq.3) then
        
        !if(prProgress.ge.1.and.update.eq.0) write(stdOut,'(A)')'  - logL'
        call pgslw(lw)
        if(animScheme.eq.3) call pgsvp(0.20,0.95,0.05,0.20)
        ic = 1
        xmin = 0.
        !xmax = real(maxval(ntot(1:nchains0)))
        xmax = maxval(is)
        dx = abs(xmax-xmin)*0.01
        !if(animScheme.eq.3) dx = abs(xmax-xmin)*0.05
        ymin =  1.e30
        ymax = -1.e30
        do ic=1,nchains0
           ymin = min(ymin,minval(post(ic,1:ntot(ic))))
           ymax = max(ymax,maxval(post(ic,1:ntot(ic))))
        end do
        ymin = 0.
        dy = abs(ymax-ymin)*0.05
        !if(animScheme.eq.3) dy = abs(ymax-ymin)*0.1
        
        call pgswin(xmin-dx,xmax+dx,ymin-dy,ymax+dy)
        if(offsetrun.eq.1) call pgswin(xmin-dx,xmax+dx,ymin,ymax+2*dy)
        !if(animScheme.eq.3) call pgbox('BCNTS',0.0,0,'BCNTS',nint((ymax/2.)/100)*100.,2)
        if(animScheme.eq.3) call pgbox('BCNTS',0.0,0,'BCNTS',10.**floor(log10(abs(ymax-ymin))),2)
        
        do ic=1,nchains0
           !call pgsci(mod(ic*2,10))
           call pgsci(5)
           if(nchains0.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
           !do i=1,ntot(ic),chainPlI
           if(chainSymbol.eq.0) then !Plot lines rather than symbols
              call pgsls(1)
              call pgline(nplt-1,is(ic,2:nplt),post(ic,2:nplt))
           else
              !call pgpoint(1,0.,startval(ic,p,2),chainSymbol)       !Starting value
              call pgsls(4)
              call pgline(2,(/-is(ic,ntot(ic)),2*is(ic,ntot(ic))/),(/startval(ic,p,2),startval(ic,p,2)/))       !Starting value
              do i=1,nplt,chainPlI
                 call pgpoint(1,is(ic,i),post(ic,i),chainSymbol)
              end do
           end if
           call pgsls(1)
        end do
        
        do ic=1,nchains0
           call pgsls(2)
           if(plInject.ge.1) then
              call pgsci(2)
              if(nchains0.gt.1) call pgsci(1)
              call pgline(2,(/-1.e20,1.e20/),(/startval(ic,p,1),startval(ic,p,1)/))  !Injection value
           end if
           call pgsci(6)
           if(nchains0.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
           !call pgline(2,(/real(Nburn(ic)),real(Nburn(ic))/),(/-1.e20,1.e20/))
           !if(plBurn.ge.1) call pgline(2,(/isburn(ic),isburn(ic)/),(/-1.e20,1.e20/))              !Burn-in
           !if(plBurn.ge.1.and.isburn(ic).lt.is(ic,ntot(ic))) call pgline(2,(/isburn(ic),isburn(ic)/),(/-1.e20,1.e20/))
           if(plBurn.ge.1.and.isburn(ic).lt.is(ic,max(min(nplt,ntot(ic)),1) ))  &
                call pgline(2,(/isburn(ic),isburn(ic)/),(/-1.e20,1.e20/))
           
           !print*,ic,nplt,isburn(ic),is(ic,nplt),is(ic,max(min(nplt,ntot(ic)),1) )
           
           !   call pgsci(2)
           !   call pgsls(4)
           !   if(nchains0.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
           !   call pgline(2,(/-1.e20,1.e20/),(/startval(ic,p,2),startval(ic,p,2)/)) !Starting value
        end do
        call pgsci(1)
        call pgsls(1)
        !call pgmtxt('T',-1.5,0.05,0.0,'log(L)')
        if(animScheme.eq.3) then
           call pgsvp(0.00,0.20,0.05,0.20)
           call pgswin(0.,1.,0.,1.)
           call pgsch(1.5)
           call pgslw(2*lw)
           call pgptxt(0.2,0.5,0.0,0.0,'log(L):')
           call pgsch(1.)
           call pgslw(lw)
        end if
     end if
     
     
     
     
     !******************************************************************************************************************************
     !Plot 1D pdf
     !if(prProgress.ge.1.and.update.eq.0) write(stdOut,'(A)')'  - pdf'
     
     if(animScheme.eq.1) call pgsvp(0.45,0.95,0.05,0.8)
     if(animScheme.eq.2) call pgsvp(0.25,0.95,0.32,0.999)
     if(animScheme.eq.3) call pgsvp(0.27,0.95,0.42,0.999)
     
     call pgsfs(fillPDF)
     
     !Set x-ranges, bin the data and get y-ranges
     xmin = 1.e30
     xmax = -1.e30
     do ic=1,nchains0
        xmin = min(xmin,minval(allDat(ic,p,1:ntot(ic))))
        xmax = max(xmax,maxval(allDat(ic,p,1:ntot(ic))))
     end do
     dx = xmax - xmin
     if(offsetrun.eq.0) then
        xmin = xmin - 0.1*dx
        xmax = xmax + 0.1*dx
     end if
     
     ic = 1
     if(mergeChains.eq.1) then
        n1 = 1
        n2 = 1
        do ic=1,nchains0
           if(nplt.gt.Nburn(ic)) then
              n2 = n1 + min(nplt,ntot(ic))-Nburn(ic) - 1
              x(1,n1:n2) = allDat(ic,p,Nburn(ic)+1:min(nplt,ntot(ic)))
              n1 = n2 + 1
           end if
        end do
     end if
     
     ic = 1
     if(nplt.gt.Nburn(ic)) then
        do ic=1,nchains
           n1 = 1
           if(mergeChains.ne.1) then
              n2 = nplt-Nburn(ic)
              x(ic,1:n2) = allDat(ic,p,Nburn(ic)+1:nplt)
           end if
           xmin1 = minval(x(ic,1:n2))
           xmax1 = maxval(x(ic,1:n2))
           call bindata1d(n2,x(ic,1:n2),1,Nbin1D,xmin1,xmax1,xbin1,ybin1) !Count the number of points in each bin
           
           if(normPDF1D.ge.1) then
              if(normPDF1D.eq.2) then !Normalise the height of the PDF
                 ybin1 = ybin1/maxval(ybin1)
              else !Normalise the SURFACE, not the height (because of different bin size).  This is the default
                 norm = 0.
                 do i=1,Nbin1D+1
                    norm = norm + ybin1(i)
                 end do
                 norm = norm*(xmax1-xmin1)
                 ybin1 = ybin1/norm
              end if
           end if
           
           !Smoothen
           if(smooth.gt.1) call smoothpdf1d(ybin1,Nbin1D+1,smooth)           
           xbin(ic,1:Nbin1D+1) = xbin1(1:Nbin1D+1)
           ybin(ic,1:Nbin1D+1) = ybin1(1:Nbin1D+1)
           
        end do !ic
        
        ymin = 0.
        ymax = -1.e30
        do ic=1,nchains
           ymax = max(ymax,maxval(ybin(ic,1:Nbin1D+1)))
        end do
        ymax = ymax*1.05
        
        call pgsch(sch)
        call pgswin(xmin,xmax,ymin,ymax)
        
        
        call pgsci(1)
        call pgslw(lw)
        
        
        
        !print*,'Plot PDF'
        !Plot 1D PDF
        do ic=1,nchains
           !Set hatch style: angle = +-45deg, phase between 0 and 1 (1/nchains0, 2/nchains0, ...)
           if(fillPDF.ge.3) call pgshs(45.0*(-1)**ic,1.0,real(ic)/real(nchains0)) 
           if(nchains.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
           xbin1(1:Nbin1D+1) = xbin(ic,1:Nbin1D+1)
           ybin1(1:Nbin1D+1) = ybin(ic,1:Nbin1D+1)
           if(wrap(ic,p).eq.0) then
              if(nchains.eq.1) call pgsci(15)
              call pgpoly(Nbin1D+2,(/xbin1(1),xbin1(1:Nbin1D+1)/),(/0.,ybin1(1:Nbin1D+1)/))
              call pgsci(1)
              if(nchains.eq.1) call pgsci(2)
              call pgline(Nbin1D+1,xbin1(1:Nbin1D+1),ybin1(1:Nbin1D+1)) !:Nbin1D) ?
              
              call pgline(2,(/xbin1(1),xbin1(1)/),(/0.,ybin1(1)/)) !Fix the loose ends
              call pgline(2,(/xbin1(Nbin1D+1),xbin1(Nbin1D+1)/),(/ybin1(Nbin1D+1),0./))
           else
              plshift = real(2*pi)
              if(changeVar.ge.1) then
                 plshift = 360.
                 if(parID(p).eq.31) plshift = 24. !RA in hours
              end if
              if(nchains.eq.1) call pgsci(15)
              call pgpoly(Nbin1D+3,(/xbin1(1),xbin1(1:Nbin1D),xbin1(1)+plshift,xbin1(1)+plshift/), &
                   (/0.,ybin1(1:Nbin1D),ybin1(1),0./))
              call pgsci(1)
              if(nchains.eq.1) call pgsci(2)
              call pgline(Nbin1D,xbin1(1:Nbin1D),ybin1(1:Nbin1D))
              
              call pgsls(4)
              call pgline(Nbin1D+1,(/xbin1(1:Nbin1D)-plshift,xbin1(1)/),(/ybin1(1:Nbin1D),ybin1(1)/))
              call pgline(Nbin1D,xbin1+plshift,ybin1)
              call pgsls(1)
              call pgline(2,(/xbin1(Nbin1D),xbin1(1)+plshift/),(/ybin1(Nbin1D),ybin1(1)/))
           end if
        end do !ic
        
        
        !Plot lines again over surface of overlapping distributions
        if(nchains.gt.1) then
           call pgsls(4)
           do ic=1,nchains
              call pgsci(1)
              xbin1(1:Nbin1D+1) = xbin(ic,1:Nbin1D+1)
              ybin1(1:Nbin1D+1) = ybin(ic,1:Nbin1D+1)
              if(wrap(ic,p).eq.0) then
                 call pgline(Nbin1D+1,xbin1(1:Nbin1D+1),ybin1(1:Nbin1D+1))
              else
                 call pgline(Nbin1D,xbin1(1:Nbin1D),ybin1(1:Nbin1D))
              end if
           end do
           call pgsls(4)
        end if
        
        
        
        
        !Plot median and signal value
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
           
           !Injection value
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
        
        
        !Print median, signal value and range widths
        call pgslw(lw)
        call pgsci(1)
        ic = 1
        if(animScheme.eq.1) then
           if(parID(p).eq.21.or.parID(p).eq.22.or.parID(p).eq.61.or.parID(p).eq.62.or.parID(p).eq.63.or.parID(p).eq.64.or. &
                parID(p).eq.71.or.parID(p).eq.81) then
              write(str,'(A,F7.3,A5,F7.3,A9,F6.2,A1)')'mdl:',startval(ic,p,1),' med:',median,' \(2030):',drange*100,'%'
           else
              write(str,'(A,F7.3,A5,F7.3,A9,F7.3)')'mdl:',startval(ic,p,1),' med:',median,' \(2030):',drange
           end if
        end if
        if(animScheme.eq.2.or.animScheme.eq.3) then
           call pgslw(lw)
           call pgsch(1.)
           call pgbox('BNTS',0.0,0,'',0.0,0) !Box for the PDF
           
           
           !Box for parameter values:
           call pgsvp(0.05,0.25,0.35,0.999)
           call pgswin(0.,1.,1.,0.)
           call pgsch(2.)
           call pgslw(lw*2)
           
           call pgptxt(0.5,0.15,0.0,0.5,trim(pgParNs(parID(p))))
           call pgslw(lw*2)
           call pgsch(1.4)
           
           call pgptxt(0.04,0.35,0.,0.,'Signal:') 
           write(str,'(F7.3)')startval(ic,p,1)
           call pgptxt(0.55,0.35,0.,0.,trim(str)) 
           
           call pgptxt(0.04,0.5,0.,0.,'Median:') 
           write(str,'(F7.3)')median
           call pgptxt(0.55,0.5,0.,0.,trim(str)) 
           
           write(str,'(A,I3,A)')'\(2030)\d',nint(ival*100.),'%\u:'
           if(nint(ival*100.).lt.100) write(str,'(A,I2,A)')'\(2030)\d',nint(ival*100.),'%\u:'
           call pgptxt(0.04,0.65,0.,0.,trim(str))
           write(str,'(F7.3)')max(drange,0.001)
           if(parID(p).eq.21.or.parID(p).eq.22.or.parID(p).eq.61.or.parID(p).eq.62.or.parID(p).eq.63.or.parID(p).eq.64.or. &
                parID(p).eq.71.or.parID(p).eq.81) write(str,'(F6.2,A1)')max(drange*100,0.01),'%' 
           call pgptxt(0.55,0.65,0.,0.,trim(str)) 
           
           call pgsch(1.)
           call pgslw(lw)
           
           write(str,'(A,ES9.2)')'Iteration:',max(is(ic,max(nplt,1)),0.)
           call pgptxt(0.15,0.84,0.,0.,trim(str)) 
           write(str,'(A,ES9.2)')'Data points:',real(n2)
           call pgptxt(0.0,0.9,0.,0.,trim(str)) 
        end if
        
     else
        
        !Before burn-in:
        call pgsch(sch)
        call pgswin(xmin,xmax,ymin,ymax)
        !end if  !if(nplt.gt.Nburn(ic)) 
        
        if(animScheme.eq.1) then
           call pgsch(1.5)
           call pgslw(lw+2)
           call pgmtxt('T',2.5,0.5,0.5,trim(pgParNs(parID(p))))
           call pgslw(lw)
           call pgsch(1.)
           !if(iframe.gt.0) call pgmtxt('T',1.,0.5,0.5,trim(str)) 
           if(nplt.gt.Nburn(ic)) call pgmtxt('T',1.,0.5,0.5,trim(str)) 
           call pgbox('BNTS',0.0,0,'',0.0,0)
        end if
        
        if(animScheme.eq.2.or.animScheme.eq.3) then
           call pgslw(lw)
           call pgsch(1.)
           call pgbox('BNTS',0.0,0,'',0.0,0) !Box for the PDF
           call pgsls(2)
           call pgline(2,(/startval(ic,p,1),startval(ic,p,1)/),(/-1.e20,1.e20/))
           call pgsls(1)
           
           
           !Box for parameter values:
           call pgsvp(0.05,0.25,0.35,0.999)
           call pgswin(0.,1.,1.,0.)
           call pgsch(2.)
           call pgslw(lw*2)
           
           call pgptxt(0.5,0.15,0.0,0.5,trim(pgParNs(parID(p))))
           call pgslw(lw*2)
           call pgsch(1.4)
           
           call pgptxt(0.04,0.35,0.,0.,'Signal:') 
           write(str,'(F7.3)')startval(ic,p,1)
           call pgptxt(0.55,0.35,0.,0.,trim(str)) 
           
           
           call pgsch(1.)
           call pgslw(lw)
           write(str,'(A,ES9.2)')'Iteration:',max(is(ic,max(nplt,1)),0.)
           call pgptxt(0.15,0.84,0.,0.,trim(str)) 
           write(str,'(A,ES9.2)')'Data points:',real(nplt)
           call pgptxt(0.0,0.9,0.,0.,trim(str)) 
        end if
     end if  !if(nplt.gt.Nburn(ic)) 
     
     !******************************************************************************************************************************
     
     
     
     
     
     
     
     call pgend
     
     if(file.ge.1) then
        call pgend
        if(small_anim.eq.1) then
           status = system('convert -resize 500 -depth 8 -unsharp '//trim(unSharppdf1d)// &
                ' analysemcmc_frame.ppm '//trim(framename))  !Rescale the output frame
        else
           status = system('convert -resize 1024x738 -depth 8 -unsharp '//trim(unSharppdf1d)// &
                ' analysemcmc_frame.ppm '//trim(framename))  !Rescale the output frame
        end if
        if(status.ne.0) write(stdErr,'(A,I6)')'  Error converting plot',status
        status = system('rm -f analysemcmc_frame.ppm')
     end if
     
  end do  !do(iframe=0,nAnimFrames)
  
  
  
  !Deallocate memory:
  deallocate(xbin,ybin,xbin1,ybin1)
  
  
end subroutine animation
!***********************************************************************************************************************************
