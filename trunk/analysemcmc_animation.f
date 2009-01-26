! Create animations for analysemcmc


subroutine animation(exitcode)
  use constants
  use analysemcmc_settings
  use general_data
  use mcmcrun_data
  use chain_data
  use plot_data
  implicit none
  integer :: c,i,i0,i1,ic,io,p,iframe,nplt,pgopen,lw,n0,n1,n2,norm,exitcode,system
  integer :: index(npar1,nchs*narr1)
  real :: range,range1,range2,drange,minrange,centre,median,plshift,ival
  real :: x(nchs,nchs*narr1),x1,x2,xmin,xmax,xmin1,xmax1,dx,y1,y2,ymin,ymax,dy,sch
  real :: coefs(100),coefs1(100)
  real,allocatable :: xbin(:,:),ybin(:,:),xbin1(:),ybin1(:),ybin2(:)    !These depend on nbin1d, allocate after reading input file
  real*8 :: ts1,ts2,timestamp
  character :: framename*99,tms*8,str*99,command*99
  
  exitcode = 0
  
  !i = system('cpu_time_type=total_alltime')
  !call cpu_time(cputime)
  ts1 = timestamp()
  write(*,*)
  p = 2 !Mc

  !moviescheme = 1  !Left column: Chain, sigma and acceptance, right column: numbers and PDF
  !moviescheme = 2  !Upper panel: numbers and PDF, lower panel: chain
  moviescheme = 3  !Upper panel: numbers and PDF, lower panel: chain and logL
  
  !Allocate memory:
  allocate(xbin(nchs,nbin1d+1),ybin(nchs,nbin1d+1),xbin1(nbin1d+1),ybin1(nbin1d+1),ybin2(nbin1d+1))

  
  do iframe = 0,nmovframes
     nplt = nint(real(iframe)/real(nmovframes)*maxval(ntot(1:nchains)))-1  !This is the line number, not the iteration number
     !nplt = nint(real(iframe)/real(nmovframes)*maxval(is))-1
     !if(prprogress.ge.1.and.update.eq.0) write(*,'(A)')' Plotting movie frame...'
     if(prprogress.ge.1.and.update.eq.0) then
        if(prprogress.ge.1) then
           write(6,*)upline !Move cursor up 1 line
           write(*,'(A,I5,A1,I5,A,I7,A,$)')' Plotting movie frame',iframe,'/',nmovframes,'  (',nplt,' points)'
        end if
        !Print remaining time
        !call cpu_time(cputime)
        !if(file.eq.1) write(6,'(A,A9)')'   Est.time left:',tms(dble(cputime-cputime0)*(nmovframes-iframe)/3600.d0)   !This does not take into account child processes, like convert
        !cputime0 = cputime
        ts2 = timestamp()
        if(prprogress.ge.1.and.file.eq.1) write(6,'(A,A9)')'   Est.time left:',tms((ts2-ts1)*(nmovframes-iframe+1)/3600.d0)                 !Use the system clock
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
        exitcode = 1
        return
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


     !if(prprogress.ge.1.and.update.eq.0) write(*,'(A)')'  - parameter chains'
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
        !if(prprogress.ge.1.and.update.eq.0) write(*,'(A)')'  - logL'
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
     !if(prprogress.ge.1.and.update.eq.0) write(*,'(A)')'  - sigma'

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
        if(fonttype.eq.2) then
           call pgmtxt('T',-1.5,0.05,0.0,'\(2144)')
        else
           call pgmtxt('T',-1.5,0.05,0.0,'\(0644)')
        end if

     end if


     !***********************************************************************************************************************************            
     !Plot acceptance rate
     !if(prprogress.ge.1.and.update.eq.0) write(*,'(A)')'  - acceptance rate'

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
     !if(prprogress.ge.1.and.update.eq.0) write(*,'(A)')'  - pdf'

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
           call bindata1d(n2,x(ic,1:n2),1,nbin1d,xmin1,xmax1,xbin1,ybin1) !Count the number of points in each bin



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

  end do  !do(iframe=0,nmovframes)
  
  
  
  !Deallocate memory:
  deallocate(xbin,ybin,xbin1,ybin1,ybin2)
  

end subroutine animation
