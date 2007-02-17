

      !Read and plot the data output from Christians code
      implicit none
      integer, parameter :: narr=100000,npar=13,nchs=1,nbin=200
      integer :: n,n1,nd,i,j,k,iargc,io,pgopen
      real :: as(narr),fs(narr),ss(narr),ps(narr),is(narr),temp(narr),pldat(npar,narr*nchs),pldat0(npar),pldat1(npar)
      real :: states(narr),accepted(narr)
      real :: Ao,fo,PHIo,sig,sn,dlp,xbin(nbin+1),ybin(nbin+1),x(narr*nchs),ybintot
      real*8 :: dat(npar,nchs,narr),tbase,r2d,r2h,pi
      character :: js*4,varnames(npar)*5,pgvarns(npar)*22,infile*50,fncore*25,str*99,fmt*99
      
      integer :: nfx,nfy,fx,fy,file,update
      real :: x1,x2,y1,y2,dx,dy,xmin,xmax,ymin,ymax
      
      integer :: nchains,ic,ip,acc(nchs,narr),i0,i1,i2
      character :: header*1000
      
      pi = 4*datan(1.d0)
      r2d = 180.d0/pi
      r2h = 12.d0/pi
      
      file = 0 !0-screen 1-file
      update = 1 !Update screen plot every 10 seconds: 0-no, 1-yes
      write(6,*)''
      !Columns in dat(): 1:logL, 2:eta, 3:mc, 4:tc, 5:logdl, 6:sinlati, 7longi:, 8:phase, 9:spin, 10:kappa, 11:thJ0, 12:phJ0, 13:alpha
      varnames = (/'logL','eta','Mc','tc','dl','lat','lon','phase','spin','kappa','thJo','phJo','alpha'/)
      pgvarns  = (/'log Likelihood        ','\(2133)               ','M\dc\u (M\d\(2281)\u) ','t\dc\u (s)            ',
     &             'log d\dL\u (Mpc)      ','lat. (\(2218))        ','lon.    (\(2218))     ','\(2149) (\(2218))     ',
     &             'spin                  ','kappa                 ','th\dJ\do\u\u (\(2218))','ph\dJ\do\u\u (\(2218))',
     &             'alpha (\(2218))       '/)
     
      infile='output.dat'
      if(iargc().eq.1) call getarg(1,infile)
      
 101  open(unit=10,form='formatted',status='old',file=trim(infile),iostat=io)
      if(io.ne.0) then
        print*,'File not found: '//trim(infile)//'. Quitting the programme.'
	goto 9999
      endif
      rewind(10)
      
      if(update.ne.1) write(6,'(A19,A,A6)')'Reading input file ',trim(infile),'...   '
!      read(10,'(A)')header
!      nchains = nint((len_trim(header)-21)/88.)
!      write(6,'(I4,A25,$)')nchains,' parallel chains found,  '
!      if(nchains.gt.nchs) then
!        write(6,'(A,I2,A)')'Arrays too small for the number of chains. Only the first ',nchs,' will be read.'
!        nchains = nchs
!      endif
      
      !Number of columns is 2 + 13*nchains
      !First 2 columns: iteration, temperature
      !Columns in dat(): 1:logL, 2:eta, 3:mc, 4:tc, 5:logdl, 6:sinlati, 7longi:, 8:phase, 9:spin, 10:kappa, 11:thJ0, 12:phJ0, 13:alpha
      
      nchains = 1
      do i=1,narr
!        read(10,191,err=195,end=199,advance='no')i0
!	is(i)=i0
!	do ic=1,nchains
!          read(10,192,err=195,end=199,advance='no')dat(:,ic,i)
!	enddo
!        read(10,192,err=195,end=199)dat(:,nchains,i)
        read(10,*,end=199,err=199)i1,i2,dat(:,1,i)
	is(i)=real(i)
        states(i) = i1
        accepted(i) = real(i2)/real(i1)
      enddo
!  191 format(I10,F6.3)
!  192 format(F13.1,E17.10,F15.1,F20.9,F14.10,4F13.10,E17.10,2E28.20,E27.20)
!  195 write(6,'(A,I)')'error reading file '//trim(infile)//' line ',i+1
  199 close(10)
      n = i-1
      write(6,'(I,A12)')n,'lines read.'
      
      tbase = floor(minval(dat(4,1:nchains,1:n)))

      !Columns in dat(): 1:logL, 2:eta, 3:mc, 4:tc, 5:logdl, 6:sinlati, 7longi:, 8:phase, 9:spin, 10:kappa, 11:thJ0, 12:phJ0, 13:alpha
      dat(4,:,:) = dat(4,:,:) - tbase	!To express tc as a float
      dat(5,:,:) = dexp(dat(5,:,:))     !logD -> Distance 
      dat(6,:,:) = dasin(dat(6,:,:))*r2d
      dat(7,:,:) = dat(7,:,:)*r2d
      dat(8,:,:) = mod(dat(8,:,:)+20.d0,1.d0)
      dat(10,:,:) = dacos(dat(10,:,:))*r2d
      dat(11,:,:) = dat(11,:,:)*r2d
      dat(12,:,:) = dat(12,:,:)*r2d
      dat(13,:,:) = dat(13,:,:)*r2d
      
      j = 1
      do i=1,n
        do ic=1,nchains
	  pldat(:,j) = real(dat(:,ic,i))
	  j = j+1
	enddo
      enddo
      n = j-1
      
      pldat0 = pldat(:,1)
      pldat1 = pldat(:,2)
      
!      write(6,'(13F20.10)')pldat(:,n-1)
!      write(6,'(13F20.10)')pldat(:,n)
      
!      goto 9999
      
      
!      !Create the base of the filename
!      write(fncore,'(A7,I5.5,A1,I2.2,A1,I2)')'output_',nd,'_',nint(log10(real(n))),'_',nint(log10(sn))
!      do i=1,len_trim(fncore)
!        if(fncore(i:i).eq.' ') write(fncore(i:i),'(A1)')'0'
!      enddo
      
      fncore(1:6) = 'output'
      
      
      if(update.ne.1) write(6,*)''
      
      if(1.eq.1) then
      if(file.eq.0) then
        io = pgopen('13/xs')
        call pgsch(1.5)
      endif
      if(file.eq.1) then
!        io = pgopen(trim(fncore)//'.ppm/ppm')
        io = pgopen('output2.eps/cps')
        call pgsch(1.2)
      endif
      if(io.le.0) then
        print*,'Cannot open PGPlot device.  Quitting the programme ',io
	goto 9999
      endif
      if(file.eq.0) call pgpap(16.4,0.57)
      if(file.eq.0) call pgsch(1.5)
      call pgsubp(4,3)
      
      do j=1,12
	call pgpage
	xmin = 0.
	xmax = real(n)
	dx = abs(xmax-xmin)*0.01
	ymin = minval(pldat(j,1:n))
	ymax = maxval(pldat(j,1:n))
	dy = abs(ymax-ymin)*0.05
	
	call pgswin(xmin-dx,xmax+dx,ymin-dy,ymax+dy)
        call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)
        
!	call pgline(n,is(1:n),pldat(j,1:n))
!	call pgpoint(n,is(1:n),pldat(j,1:n),1)
!	do i=1,n
!	  call pgpoint(1,is(i),pldat(j,i),1)
!	enddo
	n1 = 1
	if(n.gt.1000000) n1 = n/1000000
	do i=0,n1-1
	  call pgpoint(n/n1,is(n/n1*i+1:n/n1*(i+1)),pldat(j,n/n1*i+1:n/n1*(i+1)),1)
	enddo
	
        call pgsci(2)
        call pgsls(2)
        call pgline(2,(/-1.e20,1.e20/),(/pldat0(j),pldat0(j)/))
        call pgsls(4)
        call pgline(2,(/-1.e20,1.e20/),(/pldat1(j),pldat1(j)/))
        call pgsci(1)
        call pgsls(1)
        call pgmtxt('L',2.2,0.5,0.5,pgvarns(j))
        call pgmtxt('B',3.,0.5,0.5,'i')
      enddo
      call pgend
      endif !if(1.eq.2) then
      
      
      
      
      
      if(1.eq.1) then
      if(file.eq.0) then
        io = pgopen('14/xs')
        call pgsch(1.5)
      endif
!      if(file.eq.1) io = pgopen(trim(fncore)//'.ppm/ppm')
!      if(file.eq.1) io = pgopen(trim(fncore)//'.eps/cps')
      if(file.eq.1) then
        io = pgopen('output.eps/cps')
        call pgsch(1.2)
      endif
      if(io.le.0) then
        print*,'Cannot open PGPlot device.  Quitting the programme ',io
	goto 9999
      endif
      if(file.eq.0) call pgpap(16.4,0.57)
      if(file.eq.0) call pgsch(1.5)
      call pgsubp(4,3)

      if(update.ne.1) write(6,'(5A10)')'var','y_max','x_max','model','diff'
      do j=2,13
	call pgpage
	x(1:n) = pldat(j,1:n)
	xmin = 0.
	xmax = 0.
	call bindata(n,x(1:n),1,nbin,xmin,xmax,xbin,ybin) 
	ymin = minval(ybin)
	ymax = maxval(ybin)
	
	call pgswin(xmin,xmax,0.,ymax*1.1)
        call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)
        
	call verthist(nbin,xbin,ybin)
	
	!Plot optimum value and model value
        call pgsci(2)
	if(file.eq.1) call pgsch(1.5)
	do i=1,nbin
	  if(ybin(i).eq.ymax) then
	    call pgpoint(1,(xbin(i)+xbin(i+1))/2.,ybin(i),2)
	    if(update.ne.1) write(6,'(A10,3F10.4,F9.4,A1)')trim(varnames(j)),ybin(i),xbin(i),pldat0(j),
     &        abs(xbin(i)-pldat0(j))/pldat0(j)*100,'%'
	    write(str,'(A,F6.2,A7,F6.2,A11,F6.2,A1)')trim(pgvarns(j))//':    mdl: ',pldat0(j),'  opt: ',xbin(i),
     &                                              '  \(2030): ',abs(xbin(i)-pldat0(j))/pldat0(j)*100,'%'
	  endif
	enddo
	
        call pgsci(2)
        call pgsls(2)
        call pgline(2,(/pldat0(j),pldat0(j)/),(/-1.e20,1.e20/))
        call pgsls(4)
        call pgline(2,(/pldat1(j),pldat1(j)/),(/-1.e20,1.e20/))
        call pgsci(1)
        call pgsls(1)
	if(file.eq.1) call pgsch(1.2)
!        call pgmtxt('B',3.,0.5,0.5,varnames(j))
        call pgmtxt('B',3.,0.5,0.5,trim(str))
      enddo
      
      call pgsubp(1,1)
      call pgswin(-1.,1.,-1.,1.)
!!      call pgsch(1.2)
      call pgslw(3)
!      write(fmt,'(A6,I2.2,A6,I2.2,A10)')'(A10,I',ceiling(log10(real(nd)))+1,',A12,I',ceiling(log10(real(n)))+1,',A14,F6.2)'
!!      print*,trim(fmt)
!!      write(str,'(A10,I,A12,I,A15,F6.2)')'n\ddat\u: ',nd,'   n\dmc\u: ',n,'   A/\(2144)): ',sn
!      write(str,fmt)'n\ddat\u: ',nd,'   n\dmc\u: ',n,'   A/\(2144): ',sn
!      call pgptxt(0.,-0.08,0.,0.5,trim(str))
      


      
      !Make sure gv auto-reloads on change
      if(file.eq.1) then
        call pgpage
        call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)
      endif

      call pgend
      endif !if(1.eq.2) then
      
      if(update.eq.1) then
        call sleep(10)
        goto 101
      endif
      
      if(file.eq.1) then
        write(6,*)''
        write(6,'(A)')'Plot saved as: '//trim(fncore)//'.eps'
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
      ! xbin, ybin - output: binned data (x, y)
      
      implicit none
      integer :: i,k,n,nbin,norm
      real :: x(n),xbin(nbin),ybin(nbin),xmin,xmax,dx,ybintot,xmin1,xmax1
      
      xmin = xmin1
      xmax = xmax1
      
      if(abs(xmin-xmax)/(xmax+1.e-30).lt.1.e-20) then
        xmin = minval(x(1:n))
        xmax = maxval(x(1:n))
      endif
      dx = abs(xmax - xmin)/real(nbin)
        
      do k=1,nbin+1
        xbin(k) = xmin + (real(k)-0.5)*dx
      enddo
      ybintot=0.
      do k=1,nbin
        ybin(k) = 0.
        do i=1,n
          if(x(i).ge.xbin(k).and.x(i).lt.xbin(k+1)) ybin(k) = ybin(k)+1.
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
      subroutine verthist(n,x,y)  !x is the left of the column 
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
      
      
      
