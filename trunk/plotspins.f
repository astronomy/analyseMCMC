

      !Read and plot the data output from Christians code
      implicit none
      integer, parameter :: narr=100000,npar=15,nchs=2,nbin=40
      integer :: n(nchs),n1,nd,i,j,j1,j2,k,nburn,iargc,io,pgopen
      real :: as(narr),fs(narr),ss(narr),is(narr),temp(narr),pldat(npar,nchs,narr),pldat0(npar),pldat1(npar)
      real :: states(narr),accepted(narr)
      real :: sig,sn,dlp,xbin(nchs,nbin+1),ybin(nchs,nbin+1),xbin1(nbin+1),ybin1(nbin+1),x(nchs,narr),y(nchs,narr),ybintot
      real :: ysum(nbin+1),yconv(nbin+1),ycum(nbin+1),a,b,medians(npar)
      real*8 :: dat(npar,nchs,narr),tbase,r2d,r2h,pi
      character :: js*4,varnames(npar)*5,pgvarns(npar)*22,infile*50,str*99,fmt*99
      
      integer :: nfx,nfy,fx,fy,file,update,ps
      real :: x1,x2,y1,y2,dx,dy,xmin,xmax,ymin,ymax,ymaxs(nchs+2),z(nbin,nbin),zs(nchs,nbin,nbin),tr(6),cont(11)
      
      integer :: nchains,ic,ip,acc(nchs,narr),i0,i1,i2,plchain,plpdf,plpdf2d
      character :: header*1000
      
      pi = 4*datan(1.d0)
      r2d = 180.d0/pi
      r2h = 12.d0/pi
      
      file = 1 !0-screen 1-file
      ps = 1 !0: create bitmap, 1: create eps
      update = 0 !Update screen plot every 10 seconds: 0-no, 1-yes
      plchain = 0
      plpdf = 1
      plpdf2d = 1
      nburn = 50000
      
      write(6,*)''
      !Columns in dat(): 1:logL, 2:eta, 3:mc, 4:tc, 5:logdl, 6:sinlati, 7longi:, 8:phase, 9:spin, 10:kappa, 11:thJ0, 12:phJ0, 13:alpha
      varnames(1:14) = (/'logL','eta','Mc','tc','dl','lat','lon','phase','spin','kappa','thJo','phJo','alpha','accept'/)
      pgvarns(1:14)  = (/'log Likelihood        ','\(2133)               ','M\dc\u (M\d\(2281)\u) ','t\dc\u (s)            ',
     &                   'd\dL\u (Mpc)          ','lat. (\(2218))        ','lon.    (\(2218))     ','\(2147)\dc\u (\(2218))',
     &                   'spin                  ','acos(\(2136))         ','\(2134)\dJ\u (\(2218))','\(2147)\dJ\u (\(2218))',
     &                   '\(2127) (\(2218))     ','Acceptance'/)
     
      infile='output.dat'
      nchains = iargc()
      if(nchains.gt.nchs) write(6,'(A)')'Too many chains, please increase nchs'
      nchains = min(nchains,nchs)
      
 101  do ic = 1,nchains
        call getarg(ic,infile)
        
        open(unit=10,form='formatted',status='old',file=trim(infile),iostat=io)
        if(io.ne.0) then
          print*,'File not found: '//trim(infile)//'. Quitting the programme.'
	  goto 9999
        endif
        rewind(10)
        
        !if(update.ne.1) 
        write(6,'(A19,A,A6,$)')'Reading input file ',trim(infile),'...   '
        
        !Columns in dat(): 1:logL, 2:eta, 3:mc, 4:tc, 5:logdl, 6:sinlati, 7longi:, 8:phase, 9:spin, 10:kappa, 11:thJ0, 12:phJ0, 13:alpha
!        nchains = 1
        do i=1,narr
          read(10,*,end=199,err=199)i1,i2,dat(1:13,ic,i)
	  is(i)=real(i)
          states(i) = i1
          dat(14,ic,i) = real(i2)/real(i1)
        enddo
  199   close(10)
        n(ic) = i-1
        write(6,'(5x,I,A12)')n(ic),'lines read.'
      enddo !do ic = 1,nchains
      
      tbase = 1.e30
      do ic=1,nchains
        tbase = floor(min(tbase,minval(dat(4,ic,1:n(ic)))))
      enddo
      
      !Columns in dat(): 1:logL, 2:eta, 3:mc, 4:tc, 5:logdl, 6:sinlati, 7longi:, 8:phase, 9:spin, 10:kappa, 11:thJ0, 12:phJ0, 13:alpha
      dat(4,:,:) = dat(4,:,:) - tbase	!To express tc as a float
      dat(5,:,:) = dexp(dat(5,:,:))     !logD -> Distance 
      dat(6,:,:) = dasin(dat(6,:,:))*r2d
      dat(7,:,:) = dat(7,:,:)*r2d
      dat(8,:,:) = mod(dat(8,:,:)+200*pi,2*pi)*r2d
      dat(10,:,:) = dacos(dat(10,:,:))*r2d
      dat(11,:,:) = dat(11,:,:)*r2d
      dat(12,:,:) = dat(12,:,:)*r2d
      dat(13,:,:) = dat(13,:,:)*r2d
      
!      do ic=1,nchains
!        do i=1,n(ic)
!	  pldat(1:n(ic),ic,i) = real(dat(1:n(ic),ic,i))
!	enddo
!      enddo
      pldat = real(dat)
      
      pldat0 = pldat(:,1,1)
      pldat1 = pldat(:,1,2)
      
!      write(6,'(13F20.10)')pldat(:,1,n-1)
!      write(6,'(13F20.10)')pldat(:,1,n)
      
!      goto 9999
      
      
      
      
      if(update.ne.1) then
         write(6,*)''
         write(6,'(A11,I6)')'Burn-in: ',nburn
      endif
      
      
      
      
      
      
      if(plchain.eq.1) then
      if(file.eq.0) then
        io = pgopen('12/xs')
        call pgsch(1.5)
      endif
      if(file.eq.1) then
        if(ps.eq.0) io = pgopen('chain.ppm/ppm')
        if(ps.eq.1) io = pgopen('chain.eps/cps')
        call pgsch(1.2)
      endif
      if(io.le.0) then
        print*,'Cannot open PGPlot device.  Quitting the programme ',io
	goto 9999
      endif
      if(file.eq.0) call pgpap(16.4,0.57)
      if(file.eq.0) call pgsch(1.5)
      if(file.eq.1.and.ps.eq.0) call pgpap(20.4,0.75)
      if(file.eq.1) call pgsch(1.5)
      call pgsubp(1,2)
      
      ic = 1
!      do j=1,12
      do j=1,14,13
	call pgpage
	ymin =  1.e30
        ymax = -1.e30
        do ic=1,nchains
          xmin = 0.
	  xmax = real(n(ic))
	  dx = abs(xmax-xmin)*0.01
	  ymin = min(ymin,minval(pldat(j,ic,1:n(ic))))
	  ymax = max(ymax,maxval(pldat(j,ic,1:n(ic))))
	  dy = abs(ymax-ymin)*0.05
	enddo
        
	call pgswin(xmin-dx,xmax+dx,ymin-dy,ymax+dy)
        call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)
        
        if(j.eq.14) then
          call pgsci(2)
          call pgsls(2)
          call pgline(2,(/-1.e20,1.e20/),(/0.25,0.25/))
          call pgsci(1)
          call pgsls(1)
        endif
        
	do ic=1,nchains
	  n1 = 1
!	  if(n(ic).gt.1000000) n1 = n(ic)/1000000
          call pgsci(mod(ic+2,10))
          do i=0,n1-1
	    call pgpoint(n(ic)/n1,is(n(ic)/n1*i+1:n(ic)/n1*(i+1)),pldat(j,ic,n(ic)/n1*i+1:n(ic)/n1*(i+1)),1)
	  enddo
	enddo
        
        call pgsci(2)
        call pgsls(2)
        call pgline(2,(/-1.e20,1.e20/),(/pldat0(j),pldat0(j)/))
        call pgsci(6)
        call pgline(2,(/real(nburn),real(nburn)/),(/-1.e20,1.e20/))
        call pgsci(2)
        call pgsls(4)
        call pgline(2,(/-1.e20,1.e20/),(/pldat1(j),pldat1(j)/))
        call pgsci(1)
        call pgsls(1)
!        call pgmtxt('L',2.2,0.5,0.5,pgvarns(j))
        call pgmtxt('T',1.,0.5,0.5,pgvarns(j))
!        call pgmtxt('B',3.,0.5,0.5,'i')
      enddo
      call pgend
      endif !if(1.eq.2) then
      
      
      
      
      
      
      
      
      
      if(plchain.eq.1) then
      if(file.eq.0) then
        io = pgopen('13/xs')
        call pgsch(1.5)
      endif
      if(file.eq.1) then
        if(ps.eq.0) io = pgopen('chains.ppm/ppm')
        if(ps.eq.1) io = pgopen('chains.eps/cps')
        call pgsch(1.2)
      endif
      if(io.le.0) then
        print*,'Cannot open PGPlot device.  Quitting the programme ',io
	goto 9999
      endif
      if(file.eq.0) call pgpap(16.4,0.57)
      if(file.eq.0) call pgsch(1.5)
      if(file.eq.1.and.ps.eq.0) call pgpap(20.4,0.75)
      if(file.eq.1) call pgsch(1.5)
      call pgsubp(4,3)
      
      ic = 1
!      do j=1,12
      do j=2,13
	call pgpage
	ymin =  1.e30
        ymax = -1.e30
        do ic=1,nchains
          xmin = 0.
	  xmax = real(n(ic))
	  dx = abs(xmax-xmin)*0.01
	  ymin = min(ymin,minval(pldat(j,ic,1:n(ic))))
	  ymax = max(ymax,maxval(pldat(j,ic,1:n(ic))))
	  dy = abs(ymax-ymin)*0.05
	enddo
        
	call pgswin(xmin-dx,xmax+dx,ymin-dy,ymax+dy)
        call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)
        
	do ic=1,nchains
	  n1 = 1
!	  if(n(ic).gt.1000000) n1 = n(ic)/1000000
          call pgsci(mod(ic+2,10))
          do i=0,n1-1
	    call pgpoint(n(ic)/n1,is(n(ic)/n1*i+1:n(ic)/n1*(i+1)),pldat(j,ic,n(ic)/n1*i+1:n(ic)/n1*(i+1)),1)
	  enddo
	enddo
        
        call pgsci(2)
        call pgsls(2)
        call pgline(2,(/-1.e20,1.e20/),(/pldat0(j),pldat0(j)/))
        call pgsci(6)
        call pgline(2,(/real(nburn),real(nburn)/),(/-1.e20,1.e20/))
        call pgsci(2)
        call pgsls(4)
        call pgline(2,(/-1.e20,1.e20/),(/pldat1(j),pldat1(j)/))
        call pgsci(1)
        call pgsls(1)
!        call pgmtxt('L',2.2,0.5,0.5,pgvarns(j))
        call pgmtxt('T',1.,0.5,0.5,pgvarns(j))
!        call pgmtxt('B',3.,0.5,0.5,'i')
      enddo
      call pgend
      endif !if(1.eq.2) then
      
      
      
      
      
      
      
      
      
      
      if(plpdf.eq.1) then
      if(file.eq.0) then
        io = pgopen('14/xs')
        call pgsch(1.5)
      endif
      if(file.eq.1) then
        if(ps.eq.0) io = pgopen('pdfs.ppm/ppm')
        if(ps.eq.1) io = pgopen('pdfs.eps/cps')
        call pgsch(1.2)
      endif
      if(io.le.0) then
        print*,'Cannot open PGPlot device.  Quitting the programme ',io
	goto 9999
      endif
      if(file.eq.0) call pgpap(16.4,0.57)
      if(file.eq.0) call pgsch(1.5)
      if(file.eq.1.and.ps.eq.0) call pgpap(20.4,0.75)
      if(file.eq.1) call pgsch(1.5)
      if(file.eq.1.and.ps.eq.1) call pgscf(2)
      call pgsubp(4,3)
      
      ic=1
      if(update.ne.1) write(6,'(7A10)')'var','y_max','x_max','median','model','diff','diff (%)'
      do j=2,13
	call pgpage
        
	xmin =  1.e30
        xmax = -1.e30
	ymin =  1.e30
        ymax = -1.e30
	do ic=1,nchains
	  xmin = min(xmin,minval(pldat(j,ic,1:n(ic))))
	  xmax = max(xmax,maxval(pldat(j,ic,1:n(ic))))
	enddo
        dx = xmax - xmin
        
        ysum = 0.
        yconv = 1.
	do ic=1,nchains
          x(ic,1:n(ic)) = pldat(j,ic,1:n(ic))
	  call bindata(n(ic)-nburn,x(ic,nburn+1:n(ic)),1,nbin,xmin,xmax,xbin1,ybin1)
          xbin(ic,1:nbin) = xbin1(1:nbin)
          ybin(ic,1:nbin) = ybin1(1:nbin)
          ysum = ysum + ybin1
          yconv = yconv * ybin1
	  ymin = min(ymin,minval(ybin(ic,1:nbin)))
	  ymax = max(ymax,maxval(ybin(ic,1:nbin)))
!          ymax = max(ymax,maxval(ysum(1:nbin)))
          ymaxs(ic) = ymax
	enddo
        
	if(file.eq.1) call pgsch(1.5)
	call pgswin(xmin,xmax,0.,ymax*1.1)
        call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)
!        write(6,'(4F7.4)')xmin,xmax,ymin,ymax
        
        ysum = ysum/real(nchains)
        yconv = yconv**(1./real(nchains))
        
        do ic=1,nchains
          if(file.eq.1.and.ps.eq.1) call pgslw(2)
          call pgsci(mod(ic+2,10))
          call pgsci(mod(ic*2,10))
          xbin1(1:nbin) = xbin(ic,1:nbin)
          ybin1(1:nbin) = ybin(ic,1:nbin)
	  call verthist(nbin,xbin1,ybin1)
	  call pgsci(1)
          if(file.eq.1.and.ps.eq.1) call pgslw(4)
          call verthist(nbin,xbin1,ysum)
!          call pgline(nbin,xbin1,ysum)
!          call pgsci(14)
!	  call verthist(nbin,xbin1,yconv)
	enddo
        if(file.eq.1.and.ps.eq.1) call pgslw(1)
        
	!Plot median and model value
        call pgsci(6)
	if(file.eq.1) call pgsch(1.5)
        
        ycum(1) = ysum(1)
        medians(j) = 0.
        do i=2,nbin
          ycum(i) = ycum(i-1) + ysum(i)
          if(medians(j).eq.0.and.ycum(i).gt.0.5) then
             a = (ycum(i)-ycum(i-1))/(xbin1(i)-xbin1(i-1))
             b = ycum(i) - a*xbin1(i)
             medians(j) = (0.5 - b)/a + dx/real(2*nbin) !xbin1 is the left side of the bin
!             write(6,'(7F10.4)')xbin1(i-1),xbin1(i),ycum(i-1),ycum(i),a,b,medians(j)
          endif
        enddo
        
        ymax = -1.e30
        do i=1,nbin
          if(ysum(i).gt.ymax) then
            i1 = i
            ymax = ysum(i)
          endif
	enddo
!	call pgpoint(1,(xbin1(i1)+xbin1(i1+1))/2.,ysum(i1),2)
        
!	if(update.ne.1) write(6,'(A10,3F10.4,F9.4,A1)')trim(varnames(j)),ysum(i1),xbin1(i1),pldat0(j),
!     &    abs(xbin1(i1)-pldat0(j))/pldat0(j)*100,'%'
	if(update.ne.1) write(6,'(A10,5F10.4,F9.4,A1)')trim(varnames(j)),ysum(i1),xbin1(i1),medians(j),pldat0(j),
     &    abs(medians(j)-pldat0(j)),abs(medians(j)-pldat0(j))/pldat0(j)*100,'%'
!	print*,''
        if(j.eq.2.or.j.eq.3.or.j.eq.5.or.j.eq.9) then
          write(str,'(A,F7.3,A7,F7.3,A11,F6.2,A1)')trim(pgvarns(j))//':  mdl: ',pldat0(j),'  med: ',medians(j),
     &                                          '  \(2030): ',abs(medians(j)-pldat0(j))/pldat0(j)*100,'%'
        else
          write(str,'(A,F7.3,A7,F7.3,A11,F7.3)')trim(pgvarns(j))//':  mdl: ',pldat0(j),'  med: ',medians(j),
     &                                          '  \(2030): ',abs(medians(j)-pldat0(j))
        endif
        if(file.eq.1.and.ps.eq.1) call pgslw(2)
        call pgsci(6)
        call pgsls(2)
        call pgline(2,(/pldat0(j),pldat0(j)/),(/-1.e20,1.e20/))
        call pgsls(4)
        call pgline(2,(/pldat1(j),pldat1(j)/),(/-1.e20,1.e20/))
        call pgsci(1)
        call pgsls(2)
        call pgline(2,(/medians(j),medians(j)/),(/-1.e20,1.e20/))
        call pgsls(1)
        if(file.eq.1.and.ps.eq.1) call pgslw(1)
	if(file.eq.1) call pgsch(1.2)
!        call pgmtxt('B',3.,0.5,0.5,varnames(j))
        call pgmtxt('T',1.,0.5,0.5,trim(str))
      enddo !j
      
      call pgsubp(1,1)
      call pgswin(-1.,1.,-1.,1.)
      call pgslw(3)

      
      !Make sure gv auto-reloads on change
      if(file.eq.1.and.ps.eq.1) then
        call pgpage
        call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)
      endif

      call pgend
      endif !if(plpdf.eq.1) then
      
      
      
      
      
      
      
      
      if(plpdf2d.eq.1) then
      if(file.eq.0) then
        io = pgopen('14/xs')
        call pgsch(1.5)
      endif
      if(file.eq.1) then
        if(ps.eq.0) io = pgopen('pdf2d.ppm/ppm')
        if(ps.eq.1) io = pgopen('pdf2d.eps/cps')
        call pgsch(1.2)
      endif
      if(io.le.0) then
        print*,'Cannot open PGPlot device.  Quitting the programme ',io
	goto 9999
      endif
      if(file.eq.0) call pgpap(16.4,0.57)
      if(file.eq.0) call pgsch(1.5)
      if(file.eq.1.and.ps.eq.0) call pgpap(20.4,0.75)
      if(file.eq.1) call pgsch(1.5)
      if(file.eq.1.and.ps.eq.1) call pgscf(2)
      
        !Columns in dat(): 1:logL, 2:eta, 3:mc, 4:tc, 5:dl, 6:sinlati, 7longi:, 8:phase, 9:spin, 10:kappa, 11:thJ0, 12:phJ0, 13:alpha
      j1 = 3
      j2 = 2
	xmin =  1.e30
        xmax = -1.e30
	ymin =  1.e30
        ymax = -1.e30
	do ic=1,nchains
	  xmin = min(xmin,minval(pldat(j1,ic,1:n(ic))))
	  xmax = max(xmax,maxval(pldat(j1,ic,1:n(ic))))
	  ymin = min(ymin,minval(pldat(j2,ic,1:n(ic))))
	  ymax = max(ymax,maxval(pldat(j2,ic,1:n(ic))))
	enddo
        dx = xmax - xmin
        dy = ymax - ymin
        
        xmin = xmin - 0.05*dx
        xmax = xmax + 0.05*dx
        ymin = ymin - 0.05*dy
        ymax = ymax + 0.05*dy
        
        ysum = 0.
        yconv = 1.
	do ic=1,nchains
          x(ic,1:n(ic)) = pldat(j1,ic,1:n(ic))
          y(ic,1:n(ic)) = pldat(j2,ic,1:n(ic))
          
	  call bindata2d(n(ic)-nburn,x(ic,nburn+1:n(ic)),y(ic,nburn+1:n(ic)),0,nbin,nbin,xmin,xmax,ymin,ymax,z,tr)
          zs(ic,:,:) = z(:,:)
        enddo
        
        z = 0.
        do ic=1,nchains
          z(:,:) = z(:,:) + zs(ic,:,:)
        end do
!        write(6,'(A5,6F8.4)')'tr: ',tr
!        write(6,'(A5,2F10.1)')'Z: ',maxval(z),sum(z)
        z = z/(maxval(z)+1.e-30)
!        do i=1,10
!          write(6,'(10F8.4)')z(i,1:10)
!        enddo
        
        
        
        
	if(file.eq.1) call pgsch(1.5)
        call pgsvp(0.12,0.95,0.12,0.95)
	call pgswin(xmin,xmax,ymin,ymax)
!        call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)  !Do this after pggray
!        write(6,'(4F7.4)')xmin,xmax,ymin,ymax
        
        
        
        
        call pggray(z,nbin,nbin,1,nbin,1,nbin,1.,0.,tr)
        
        do i=1,11
          cont(i) = 0.01 + 2*real(i-1)/10.
        enddo
        call pgsls(1)
        call pgslw(6)
        call pgsci(0)
        call pgcont(z,nbin,nbin,1,nbin,1,nbin,cont,4,tr)
        call pgslw(2)
        call pgsci(1)
        call pgcont(z,nbin,nbin,1,nbin,1,nbin,cont,4,tr)
        
        
        
        if(file.eq.1.and.ps.eq.1) call pgslw(2)
        call pgsci(6)
        call pgsls(2)
        call pgline(2,(/pldat0(j1),pldat0(j1)/),(/-1.e20,1.e20/))
        call pgline(2,(/-1.e20,1.e20/),(/pldat0(j2),pldat0(j2)/))
        call pgsls(4)
        call pgline(2,(/pldat1(j1),pldat1(j1)/),(/-1.e20,1.e20/))
        call pgline(2,(/-1.e20,1.e20/),(/pldat1(j2),pldat1(j2)/))
        call pgsci(1)
        call pgsls(2)
        call pgline(2,(/medians(j1),medians(j1)/),(/-1.e20,1.e20/))
        call pgline(2,(/-1.e20,1.e20/),(/medians(j2),medians(j2)/))
        call pgsls(1)
        if(file.eq.1.and.ps.eq.1) call pgslw(1)
	if(file.eq.1) call pgsch(1.5)
        
        
        call pgsls(1)
        call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)
        
        call pgmtxt('B',2.5,0.5,0.5,pgvarns(j1))
        call pgmtxt('L',2.5,0.5,0.5,pgvarns(j2))
        
        
        
        
        
      
      
      
      !Make sure gv auto-reloads on change
      if(1.eq.2.and.file.eq.1.and.ps.eq.1) then
        call pgpage
        call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)
      endif

      call pgend
      endif !if(plpdf2d.eq.1) then
      
      
      
      
      
      
      
      
      
      
      if(update.eq.1) then
        call sleep(10)
        goto 101
      endif
      
      
 9999 write(6,*)''
      end
************************************************************************************************************************************
      
      
      
      
************************************************************************************************************************************
      subroutine bindata(n,x,norm,nbin,xmin1,xmax1,xbin,ybin)  
************************************************************************************************************************************
      ! x - input: data, n points
      ! norm - input: normalise (1) or not (0)
      ! nbin - input: number of bins
      ! xmin, xmax - in/output: set xmin=xmax to auto-determine
      ! xbin, ybin - output: binned data (x, y).  The x values are the left of the bin!
      
      implicit none
      integer :: i,k,n,nbin,norm
      real :: x(n),xbin(nbin),ybin(nbin),xmin,xmax,dx,ybintot,xmin1,xmax1
      
      xmin = xmin1
      xmax = xmax1
      
      if(abs(xmin-xmax)/(xmax+1.e-30).lt.1.e-20) then !Autodetermine
        xmin = minval(x(1:n))
        xmax = maxval(x(1:n))
      endif
      dx = abs(xmax - xmin)/real(nbin)
        
      do k=1,nbin+1
!        xbin(k) = xmin + (real(k)-0.5)*dx  !x is the centre of the bin
        xbin(k) = xmin + (k-1)*dx          !x is the left of the bin
      enddo
      ybintot=0.
      do k=1,nbin
        ybin(k) = 0.
        do i=1,n
          if(x(i).ge.xbin(k).and.x(i).lt.xbin(k+1)) ybin(k) = ybin(k) + 1.
        enddo
        ybintot = ybintot + ybin(k)
      enddo
      if(norm.eq.1) ybin = ybin/(ybintot+1.e-30)
      
      if(abs(xmin1-xmax1)/(xmax1+1.e-30).lt.1.e-20) then
        xmin1 = xmin
	xmax1 = xmax
      endif
      
      end subroutine bindata
************************************************************************************************************************************
      
      
************************************************************************************************************************************
      subroutine bindata2d(n,x,y,norm,nxbin,nybin,xmin1,xmax1,ymin1,ymax1,z,tr)  
************************************************************************************************************************************
      ! x - input: data, n points
      ! norm - input: normalise (1) or not (0)
      ! nbin - input: number of bins
      ! xmin, xmax - in/output: set xmin=xmax to auto-determine
      ! xbin, ybin - output: binned data (x, y).  The x values are the left of the bin!
      
      implicit none
      integer :: i,k,n,ix,iy,bx,by,nxbin,nybin,norm
      real :: x(n),y(n),xbin(nxbin),ybin(nybin),z(nxbin,nybin),ztot,xmin,xmax,ymin,ymax,dx,dy,xmin1,xmax1,ymin1,ymax1
      real :: tr(6),zmax
      
!      write(6,'(A4,5I8)')'n:',nx,ny,norm,nxbin,nybin
!      write(6,'(A4,2F8.3)')'x:',xmin1,xmax1
!      write(6,'(A4,2F8.3)')'y:',ymin1,ymax1
      
      xmin = xmin1
      xmax = xmax1
      ymin = ymin1
      ymax = ymax1
      
      if(abs(xmin-xmax)/(xmax+1.e-30).lt.1.e-20) then !Autodetermine
        xmin = minval(x(1:n))
        xmax = maxval(x(1:n))
      endif
      dx = abs(xmax - xmin)/real(nxbin)
      if(abs(ymin-ymax)/(ymax+1.e-30).lt.1.e-20) then !Autodetermine
        ymin = minval(y(1:n))
        ymax = maxval(y(1:n))
      endif
      dy = abs(ymax - ymin)/real(nybin)
      do bx=1,nxbin+1
!        xbin(bx) = xmin + (real(bx)-0.5)*dx  !x is the centre of the bin
        xbin(bx) = xmin + (bx-1)*dx          !x is the left of the bin
      enddo
      do by=1,nybin+1
!        ybin(by) = ymin + (real(by)-0.5)*dy  !y is the centre of the bin
        ybin(by) = ymin + (by-1)*dy          !y is the left of the bin
      enddo
      
!      write(6,'(50F5.2)'),x(1:50)
!      write(6,'(50F5.2)'),y(1:50)
!      write(6,'(20F8.5)'),xbin
!      write(6,'(20F8.5)'),ybin
      
      z = 0.
      ztot=0.
      do bx=1,nxbin
      do by=1,nybin
        z(bx,by) = 0.
        do i=1,n
!        do ix=1,nx,10000
!        do iy=1,ny,10000
          if(x(i).ge.xbin(bx).and.x(i).lt.xbin(bx+1) .and. y(i).ge.ybin(by).and.y(i).lt.ybin(by+1)) z(bx,by) = z(bx,by) + 1.
!          write(6,'(2I4,2I9,7F10.5)')bx,by,ix,iy,x(ix),xbin(bx),xbin(bx+1),y(iy),ybin(by),ybin(by+1),z(bx,by)
!        enddo
        enddo
        ztot = ztot + z(bx,by)
      enddo
      enddo
!      if(norm.eq.1) z = z/(ztot+1.e-30)
      if(norm.eq.1) z = z/maxval(z+1.e-30)
      
      if(abs(xmin1-xmax1)/(xmax1+1.e-30).lt.1.e-20) then
        xmin1 = xmin
	xmax1 = xmax
      endif
      if(abs(ymin1-ymax1)/(ymax1+1.e-30).lt.1.e-20) then
        ymin1 = ymin
	ymax1 = ymax
      endif
      
      !Determine transformation elements for pgplot
      tr(1) = xmin
      tr(2) = dx
      tr(3) = 0.
      tr(4) = ymin
      tr(5) = 0.
      tr(6) = dy
      
      end subroutine bindata2d
************************************************************************************************************************************
      
      
************************************************************************************************************************************
      subroutine verthist(n,x,y)  !x is the left of the bin!
************************************************************************************************************************************
      implicit none
      integer :: j,n
      real :: x(n+1),y(n+1)
      
!      call pgline(2,(/0.,x(1)/),(/y(1),y(1)/))
!      do j=1,n
!        call pgline(2,(/x(j),x(j)/),y(j:j+1))
!        call pgline(2,x(j:j+1),(/y(j+1),y(j+1)/))
!      enddo
      
      x(n+1) = x(n) + (x(n)-x(n-1))
      y(n+1) = 0.
      call pgline(2,(/x(1),x(1)/),(/0.,y(1)/))
      do j=1,n
        call pgline(2,x(j:j+1),(/y(j),y(j)/))
        call pgline(2,(/x(j+1),x(j+1)/),y(j:j+1))
      enddo
      
      end
************************************************************************************************************************************
      
************************************************************************************************************************************
      subroutine horzhist(n,x,y)
************************************************************************************************************************************
      implicit none
      integer :: j,n
      real :: x(n),y(n)
      
      call pgline(2,(/x(1),x(1)/),(/0.,y(1)/))
      do j=1,n-2
        call pgline(2,x(j:j+1),(/y(j),y(j)/))
        call pgline(2,(/x(j+1),x(j+1)/),y(j:j+1))
      enddo
      call pgline(2,x(n-1:n),(/y(n-1),y(n-1)/))
      end
************************************************************************************************************************************
      
      
      
