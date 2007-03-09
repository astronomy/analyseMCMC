      !Read and plot the data output from the spinning MCMC code
      
      implicit none
      integer, parameter :: narr=1000010,npar1=15,nchs=2,nbin=50,nconf1=5,nr1=5
      integer :: n(nchs),n0,n1,n2,nd,i,j,j1,j2,k,nburn,iargc,io,pgopen,system,cindex(npar1,nchs,narr),index(npar1,nchs*narr),npar
      real :: is(narr),pldat(npar1,nchs,narr),alldat(npar1,nchs*narr),pldat0(npar1),pldat1(npar1)
      real :: states(narr),accepted(narr),sig(npar1,nchs,narr),acc(npar1,nchs,narr)
      real :: sn,dlp,xbin(nchs,nbin+1),ybin(nchs,nbin+1),xbin1(nbin+1),ybin1(nbin+1),ybin2(nbin+1),x(nchs,narr),y(nchs,narr),ybintot
      real :: ysum(nbin+1),yconv(nbin+1),ycum(nbin+1),a,b
      real*8 :: dat(npar1,nchs,narr),tbase,r2d,r2h,pi,tpi,dvar
      character :: js*4,varnames(npar1)*5,pgvarns(npar1)*22,infile*50,str*99,fmt*99,bla
      
      integer :: nfx,nfy,fx,fy,file,update,ps,ps2pdf
      real :: x1,x2,y1,y2,dx,dy,xmin,xmax,ymin,ymax,ymaxs(nchs+2),z(nbin,nbin),zs(nchs,nbin,nbin),tr(6),cont(11)
      real :: coefs(100),coefs1(100)
      
      integer :: nchains,ic,ip,i0,i1,i2,plchain,plpdf,plpdf2d,chpli
      character :: header*1000
      
      integer :: nt,nr,c,wrap(npar1),nconf,changevar
      real :: range,minrange,ranges(nconf1,npar1,nr1),conf,confs(nconf1),centre,shift
      real :: median(npar1),mean(npar1),vari1(npar1),vari2(npar1),absvar1(npar1),absvar2(npar1)
      real :: medians(nchs,npar1)
      
      
      pi = 4*datan(1.d0)
      tpi = 2*pi
      r2d = 180.d0/pi
      r2h = 12.d0/pi
      
      file = 1     !0-screen 1-file
      ps = 1       !0: create bitmap, 1: create eps
      ps2pdf = 0   !Convert eps to pdf (and remove eps)
      update = 0   !Update screen plot every 10 seconds: 0-no, 1-yes
      plchain = 1
      plpdf = 1
      plpdf2d = 0 !Fix pldat -> alldat
      nburn = 1!00000
!      chpli = 10   !Plot every chpli'th member of a chain (to reduce file size and overlap between plotted chains)
!      if(file.eq.0.or.ps.eq.0) chpli = 1
      chpli = 1  !Change for large n
      
      
      write(6,*)''
      !Columns in dat(): 1:logL, 2:Mc, 3:eta, 4:tc, 5:logdl, 6:sinlati, 7longi:, 8:phase, 9:spin, 10:kappa, 11:thJ0, 12:phJ0, 13:alpha
      varnames(1:14) = (/'logL','Mc','eta','tc','log dl','sin lat','lon','phase','spin','kappa','sin thJo','phJo','alpha','accept'/)
      pgvarns(1:14)  = (/'log Likelihood        ','M\dc\u (M\d\(2281)\u) ','\(2133)               ','t\dc\u (s)            ',
     &                   'logd\dL\u (Mpc)       ','sin lat.              ','lon.    (rad)         ','\(2147)\dc\u (rad)    ',
     &                   'spin                  ','\(2136)               ','sin \(2134)\dJ\u      ','\(2147)\dJ\u (rad)    ',
     &                   '\(2127) (rad)         ','Acceptance'/)
      
      nconf = 1  !Number of ranges of 'confidence levels', start with 0.90!
      confs(1:nconf) = (/0.9/)
!      confs(1:nconf) = (/0.9,0.95,0.997/)
      
      
      npar = 13
      infile='output.dat'
      nchains = iargc()
      if(nchains.lt.1) then
         write(6,'(A)')'  Syntax: plotspins <file1> [file2] ...'
         goto 9999
      end if
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
        
        !Columns in dat(): 1:logL, 2:Mc, 3:eta, 4:tc, 5:logdl, 6:sinlati, 7longi:, 8:phase, 9:spin, 10:kappa, 11:thJ0, 12:phJ0, 13:alpha
!        nchains = 1
        do i=1,4
           read(10,*,end=199,err=199)bla
        end do
        do i=1,narr
!          read(10,*,end=199,err=199)i1,i2,dat(1:13,ic,i)
          read(10,*,end=199,err=199)i1,dat(1,ic,i),(dat(j,ic,i),sig(j,ic,i),acc(j,ic,i),j=2,npar)
	  is(i)=real(i)
          states(i) = i1
!          dat(14,ic,i) = real(i2)/real(i1)
        enddo
  199   close(10)
        n(ic) = i-1
        write(6,'(5x,I,A12)')n(ic),'lines read.'
      enddo !do ic = 1,nchains
      
!      n = 1000
      nburn = nburn+2
      if(update.ne.1) then
         write(6,*)''
         write(6,'(A11,I6)')'Burn-in: ',nburn
         do ic=1,nchains
            if(nburn.ge.n(ic)) then
               write(6,'(A)')'  *** WARNING ***  Nburn larger than Nchain, cannot plot PDFs'
               plpdf = 0
               plpdf2d = 0
               exit
            end if
         end do
      endif
      
      
      
      
      
      
      tbase = 1.e30
      do ic=1,nchains
        tbase = floor(min(tbase,minval(dat(4,ic,1:n(ic)))))
      enddo
      
      
      if(update.eq.0) write(6,'(A,$)')'Changing some variables...   '
      !Columns in dat(): 1:logL, 2:Mc, 3:eta, 4:tc, 5:logdl, 6:sinlati, 7longi:, 8:phase, 9:spin, 10:kappa, 11:thJ0, 12:phJ0, 13:alpha
      do ic=1,nchains
         dat(4,ic,1:n(ic)) = dat(4,ic,1:n(ic)) - tbase	!To be able to express tc as a float
!         dat(5,ic,1:n(ic)) = dexp(dat(5,ic,1:n(ic)))     !logD -> Distance 
!         dat(6,ic,1:n(ic)) = dasin(dat(6,ic,1:n(ic)))*r2d
!         dat(7,ic,1:n(ic)) = dat(7,ic,1:n(ic))*r2d
!         dat(8,ic,1:n(ic)) = dat(8,ic,1:n(ic))*r2d
!         dat(10,ic,1:n(ic)) = dacos(dat(10,ic,1:n(ic)))*r2d
!         dat(11,ic,1:n(ic)) = dasin(dat(11,ic,1:n(ic)))*r2d
!         dat(12,ic,1:n(ic)) = dat(12,ic,1:n(ic))*r2d
!         dat(13,ic,1:n(ic)) = dat(13,ic,1:n(ic))*r2d
      end do
      
      acc = acc*0.25  !Transfom back to the actual acceptance rate
      
      do ic=1,nchains
         pldat(:,ic,1:n(ic)) = real(dat(:,ic,1:n(ic)))
      enddo
      
      pldat0 = pldat(:,1,1)
      pldat1 = pldat(:,1,2)
      if(update.eq.0) write(6,'(A)')'  Done.'
      
      
      
      
      
      
      
      
      
      
      
! **********************************************************************************************************************************
! ***  DO STATISTICS   *************************************************************************************************************
! **********************************************************************************************************************************
      
      
      !Do statistics per chain
      if(1.eq.2) then
      do c=1,1!nconf
      conf = confs(c)
      write(6,'(A,F10.5)')'Confidence level: ',conf
      do i=2,npar
      minrange = -1.e30
      write(6,'(A8,$)')varnames(i)
      do ic=1,nchains
         call dindexx(n(ic)-nburn,dat(i,ic,nburn+1:n(ic)),cindex(i,ic,1:n(ic)-nburn))
         medians(ic,i) = dat(i,ic,cindex(i,ic,(n(ic)-nburn)/2))
         write(6,'(F12.4,$)')medians(ic,i)
         
         
!         do j=nburn,(n(ic)-nburn)*(1.-conf)+nburn
!            
!            b1 = dat(i,ic,j)
!            b2 = dat(i,ic,j+n(ic)*conf)
!            range = b2 - b1
!            print*,j,b1,b2,range,minrange
!         end do
         
      end do !ic
      write(6,*)''
      end do !i
      write(6,*)''
      end do !c
      end if
      
      
      
      
      
      !Add chains, leave out burnin
      j = 1
      do ic=1,nchains
         do i=nburn,n(ic)
            alldat(:,j) = real(dat(:,ic,i))
            j = j+1
         end do
      end do
      nt = j-1
      write(6,'(I8,A)') nt,' datapoints in combined chains'
      if(nt.gt.10000) chpli = nint(real(nt)/10000.)  !Change the number of points plotted in chains
      
      
      
      
      
      !Sort all data and find the 90% confidence level for the wrapping parameters
      shift = 0.
      wrap = 0
      conf = 0.90
      write(6,'(A,F10.5)')'Confidence level: ',conf
      write(6,'(A8,12A8,A4)')'param.','model','median','mean','vari1','vari2','abvar1','abvar2',
     &                          'rng_c','rng1','rng2','drng','d/drng','ok?'
      do i=2,npar
         call rindexx(nt,alldat(i,1:nt),index(i,1:nt))
         if(i.ne.7.and.i.ne.8.and.i.ne.12.and.i.ne.13) cycle
         minrange = 1.e30
         
         do j=1,nt
            x1 = alldat(i,index(i,j))
            x2 = alldat(i,index(i,mod(j+nint(nt*conf)-1,nt)+1))
            range = mod(x2 - x1 + real(20*pi),real(tpi))
            if(range.lt.minrange) then
               minrange = range
               y1 = x1
               y2 = x2
               !write(6,'(2I6,7F10.5)')j,mod(nint(j+nt*conf),nt),x1,x2,range,minrange,y1,y2,(y1+y2)/2.
            end if
            !write(6,'(2I6,7F10.5)')j,mod(nint(j+nt*conf),nt),x1,x2,range,minrange,y1,y2,(y1+y2)/2.
         end do
         centre = (y1+y2)/2.
         if(y1.gt.y2) then
            wrap(i) = 1
            centre = mod(centre + pi, tpi) !Then distribution peaks close to 0/2pi, shift centre by pi
         end if
         !write(6,'(12x,7F10.5)')y1,y2,minrange,centre
         
         !now, wrap around anticentre
         shift = tpi - mod(centre + pi, tpi)
         alldat(i,1:nt) = mod(alldat(i,1:nt)+shift,tpi)-shift
         y1 = mod(y1+shift,tpi)-shift
         y2 = mod(y2+shift,tpi)-shift
         centre = mod(centre+shift,tpi)-shift
         minrange = y2-y1
         call rindexx(nt,alldat(i,1:nt),index(i,1:nt))  !Re-sort
         write(6,'(A8,4x,4F10.5,I4)')varnames(i),y1,y2,minrange,centre,wrap(i)
      end do
      
      
!      goto 9999
      
      
      
      
      
      !Do statistics for the whole data set
      do c=1,nconf
         conf = confs(c)
         !      conf = 0.90
         write(6,'(A,F10.5)')'Confidence level: ',conf
         write(6,'(A8,12A8,A4)')'param.','model','median','mean','vari1','vari2','abvar1','abvar2',
     &                          'rng_c','rng1','rng2','drng','d/drng','ok?'
         do i=2,npar
!            minrange = 1.e30
!            if(1.eq.2) then
!            if(c.eq.1) call rindexx(nt,alldat(i,1:nt),index(i,1:nt))
!            
!            
! !           if((i.eq.7.or.i.eq.8.or.i.eq.12.or.i.eq.13).and.c.eq.1) then
!               do j=1,nt
!                  x1 = alldat(i,index(i,j))
!                  x2 = alldat(i,index(i,mod(j+nint(nt*conf)-1,nt)+1))
!                  range = mod(x2 - x1 + real(20*pi),real(tpi))
!                  if(range.lt.minrange) then
!                     minrange = range
!                     y1 = x1
!                     y2 = x2
!                     !                  write(6,'(2I6,7F10.5)')j,mod(nint(j+nt*conf),nt),x1,x2,range,minrange,y1,y2,(y1+y2)/2.
!                  end if
!                  !            write(6,'(2I6,7F10.5)')j,mod(nint(j+nt*conf),nt),x1,x2,range,minrange,y1,y2,(y1+y2)/2.
!               end do
!               centre = (y1+y2)/2.
!               if(c.eq.1.and.y1.gt.y2) then
!                  wrap(i) = 1
!                  centre = mod(centre + pi, tpi) !Then distribution peaks close to 0/2pi, shift centre by pi
!               end if
!               !          write(6,'(12x,7F10.5)')y1,y2,minrange,centre
!               
!               !now, wrap around anticentre
!               if(c.eq.1) then
!                  shift = tpi - mod(centre + pi, tpi)
!                  alldat(i,1:nt) = mod(alldat(i,1:nt)+shift,tpi)-shift
!                  y1 = mod(y1+shift,tpi)-shift
!                  y2 = mod(y2+shift,tpi)-shift
!                  centre = mod(centre+shift,tpi)-shift
!                  minrange = y2-y1
!                  call rindexx(nt,alldat(i,1:nt),index(i,1:nt))  !Re-sort
!               end if
!               !        write(6,'(12x,7F10.5)')y1,y2,minrange,centre
!               
!!            else
!            endif
            
            
            
            
            
            
            
            
           minrange = 1.e30
!           write(6,'(A8,4x,4F10.5,I4)')varnames(i),y1,y2,minrange,centre,wrap(i)
           do j=1,floor(nt*(1.-conf))
              x1 = alldat(i,index(i,j))
              x2 = alldat(i,index(i,j+floor(nt*conf)))
              range = abs(x2 - x1)
              !               range = x2 - x1
              if(range.lt.minrange) then
                 minrange = range
                 y1 = x1
                 y2 = x2
              end if
              !            write(6,'(I6,7F10.5)')j,x1,x2,range,minrange,y1,y2,(y1+y2)/2.
           end do
           centre = (y1+y2)/2.
!           write(6,'(A8,4x,4F10.5,I4)')varnames(i),y1,y2,minrange,centre,wrap(i)
           
           
               
               
 !           end if
            
            !Determine the median
            if(mod(nt,2).eq.0) median(i) = 0.5*(alldat(i,index(i,nt/2)) + alldat(i,index(i,nt/2+1)))
            if(mod(nt,2).eq.1) median(i) = alldat(i,index(i,(nt+1)/2))
            
            !Mean:
            mean(i) = sum(alldat(i,1:nt))/real(nt)
            
            !Variances, etc:
            vari1(i)=0. ; vari2(i)=0. ; absvar1(i)=0. ; absvar2(i)=0.
            do j=1,nt
               vari1(i) = vari1(i) + (alldat(i,j) - median(i))*(alldat(i,j) - median(i))
               vari2(i) = vari2(i) + (alldat(i,j) - mean(i))*(alldat(i,j) - mean(i))
               absvar1(i) = absvar1(i) + abs(alldat(i,j) - median(i))
               absvar2(i) = absvar2(i) + abs(alldat(i,j) - mean(i))
            end do
            vari1(i)   = vari1(i)/real(nt-1)
            vari2(i)   = vari2(i)/real(nt-1)
            absvar1(i) = absvar1(i)/real(nt)
            absvar2(i) = absvar2(i)/real(nt)
            
            nr = 4
            ranges(c,i,1) = y1
            ranges(c,i,2) = y2
            ranges(c,i,3) = (y1+y2)/2.
            ranges(c,i,4) = y2-y1


            
            write(6,'(A8,12F8.4,$)')varnames(i),pldat0(i),median(i),mean(i),vari1(i),vari2(i),absvar1(i),absvar2(i),
     &                              centre,y1,y2,minrange,abs(pldat0(i)-median(i))/minrange
            if(pldat0(i).gt.y1.and.pldat0(i).lt.y2) then
               write(6,'(A4)')'y '
            else
               write(6,'(A4)')'*N*'
            end if
            
         end do !i
         write(6,*)''
      end do !c
      
      
      
      !Change variables
      if(update.eq.0) write(6,'(A,$)')'Changing some variables...   '
      !Columns in alldat(): 1:logL, 2:Mc, 3:eta, 4:tc, 5:logdl, 6:sinlati, 7longi:, 8:phase, 9:spin, 10:kappa, 11:thJ0, 12:phJ0, 13:alpha
      
      changevar = 0
      if(1.eq.1) then
      changevar = 1
      do i=2,npar
         if(i.eq.5) then
            alldat(i,1:nt) = exp(alldat(i,1:nt))     !logD -> Distance
            pldat0(i) = exp(pldat0(i))
            pldat1(i) = exp(pldat1(i))
            median(i) = exp(median(i))
            ranges(1:nconf,i,1:nr) = exp(ranges(1:nconf,i,1:nr))
         end if
         if(i.eq.6.or.i.eq.11) then
            alldat(i,1:nt) = asin(alldat(i,1:nt))*real(r2d)
            pldat0(i) = asin(pldat0(i))*real(r2d)
            pldat1(i) = asin(pldat1(i))*real(r2d)
            median(i) = asin(median(i))*real(r2d)
            ranges(1:nconf,i,1:nr) = asin(ranges(1:nconf,i,1:nr))*real(r2d)
         end if
         if(i.eq.10) then
            alldat(i,1:nt) = acos(alldat(i,1:nt))*real(r2d)
            pldat0(i) = acos(pldat0(i))*real(r2d)
            pldat1(i) = acos(pldat1(i))*real(r2d)
            median(i) = acos(median(i))*real(r2d)
            ranges(1:nconf,i,1:nr) = acos(ranges(1:nconf,i,1:nr))*real(r2d)
         end if
         if(i.eq.7.or.i.eq.8.or.i.eq.12.or.i.eq.13) then
            alldat(i,1:nt) = alldat(i,1:nt)*real(r2d)
            pldat0(i) = pldat0(i)*real(r2d)
            pldat1(i) = pldat1(i)*real(r2d)
            median(i) = median(i)*real(r2d)
            ranges(1:nconf,i,1:nr) = ranges(1:nconf,i,1:nr)*real(r2d)
         end if
      end do
      varnames(1:14) = (/'logL','Mc','eta','tc','dl','lat','lon','phase','spin','kappa','thJo','phJo','alpha','accept'/)
      pgvarns(1:14)  = (/'log Likelihood        ','M\dc\u (M\d\(2281)\u) ','\(2133)               ','t\dc\u (s)            ',
     &                  'd\dL\u (Mpc)          ','lat. (\(2218))        ','lon.    (\(2218))     ','\(2147)\dc\u (\(2218))',
     &                  'spin                  ','acos(\(2136))         ','\(2134)\dJ\u (\(2218))','\(2147)\dJ\u (\(2218))',
     &                  '\(2127) (\(2218))     ','Acceptance'/)
      if(update.eq.0) write(6,'(A)')'  Done.'
      endif
      
      
      
      
      
      
      
      
      
      
      
      
! **********************************************************************************************************************************
! ***  CREATE PLOTS   **************************************************************************************************************
! **********************************************************************************************************************************
      
      
      
      
      !Plot likelihood chain
      if(plchain.eq.1) then
      if(update.eq.0) write(6,'(A)')'Plotting chain likelihood...'
      if(file.eq.0) then
        io = pgopen('12/xs')
        call pgsch(1.5)
      endif
      if(file.eq.1) then
        if(ps.eq.0) io = pgopen('logL.ppm/ppm')
        if(ps.eq.1) io = pgopen('logL.eps/cps')
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
      call pgscr(3,0.,0.5,0.)
!      call pgsubp(1,2)
      
      ic = 1
      j=1
!      do j=1,12
!      do j=1,14,13
	call pgpage
        xmax = -1.e30
	ymin =  1.e30
        ymax = -1.e30
        do ic=1,nchains
          xmin = 0.
	  xmax = max(xmax,real(n(ic)))
	  dx = abs(xmax-xmin)*0.01
	  ymin = min(ymin,minval(pldat(j,ic,1:n(ic))))
	  ymax = max(ymax,maxval(pldat(j,ic,1:n(ic))))
	  dy = abs(ymax-ymin)*0.05
	enddo
        
	call pgswin(xmin-dx,xmax+dx,ymin-dy,ymax+dy)
        call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)
        
	do ic=1,nchains
          call pgsci(mod(ic*2,10))
!	  n1 = 1
!	  if(n(ic).gt.1000000) n1 = n(ic)/1000000
!          do i=0,n1-1
!	    call pgpoint(n(ic)/n1,is(n(ic)/n1*i+1:n(ic)/n1*(i+1)),pldat(j,ic,n(ic)/n1*i+1:n(ic)/n1*(i+1)),1)
!	  enddo
          do i=1,n(ic),chpli
            call pgpoint(1,is(i),pldat(j,ic,i),1)
	  enddo
	enddo
        
        call pgsci(3)
        call pgsls(2)
        call pgline(2,(/-1.e20,1.e20/),(/pldat0(j),pldat0(j)/))
        call pgsci(6)
        call pgline(2,(/real(nburn),real(nburn)/),(/-1.e20,1.e20/))
        call pgsci(3)
        call pgsls(4)
        call pgline(2,(/-1.e20,1.e20/),(/pldat1(j),pldat1(j)/))
        call pgsci(1)
        call pgsls(1)
!        call pgmtxt('L',2.2,0.5,0.5,pgvarns(j))
        call pgmtxt('T',1.,0.5,0.5,pgvarns(j))
!        call pgmtxt('B',3.,0.5,0.5,'i')
!      enddo
      call pgend
      if(ps2pdf.eq.1) then
         i = system('eps2pdf logL.eps >& /dev/null')
         i = system('rm -f logL.eps')
      endif
      endif !if(1.eq.2) then
      
      
      
      
      
      
      
      
      !Plot chain for each parameter
      if(plchain.eq.1) then
      if(update.eq.0) write(6,'(A)')'Plotting parameter chains...'
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
      call pgscr(3,0.,0.5,0.)

      ic = 1
!      do j=1,12
      do j=2,npar
	call pgpage
        xmax = -1.e30
	ymin =  1.e30
        ymax = -1.e30
        do ic=1,nchains
          xmin = 0.
	  xmax = max(xmax,real(n(ic)))
	  dx = abs(xmax-xmin)*0.01
	  ymin = min(ymin,minval(pldat(j,ic,1:n(ic))))
	  ymax = max(ymax,maxval(pldat(j,ic,1:n(ic))))
	  dy = abs(ymax-ymin)*0.05
	enddo
        
	call pgswin(xmin-dx,xmax+dx,ymin-dy,ymax+dy)
        call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)
        
	do ic=1,nchains
          call pgsci(mod(ic*2,10))
!	  n1 = 1
!	  if(n(ic).gt.1000000) n1 = n(ic)/1000000
!          do i=0,n1-1
!	    call pgpoint(n(ic)/n1,is(n(ic)/n1*i+1:n(ic)/n1*(i+1)),pldat(j,ic,n(ic)/n1*i+1:n(ic)/n1*(i+1)),1)
!	  enddo
          do i=1,n(ic),chpli
            call pgpoint(1,is(i),pldat(j,ic,i),1)
	  enddo
	enddo
        
        call pgsci(3)
        call pgsls(2)
        call pgline(2,(/-1.e20,1.e20/),(/pldat0(j),pldat0(j)/))
        call pgsci(6)
        call pgline(2,(/real(nburn),real(nburn)/),(/-1.e20,1.e20/))
        call pgsci(3)
        call pgsls(4)
        call pgline(2,(/-1.e20,1.e20/),(/pldat1(j),pldat1(j)/))
        call pgsci(1)
        call pgsls(1)
!        call pgmtxt('L',2.2,0.5,0.5,pgvarns(j))
        call pgmtxt('T',1.,0.5,0.5,'Chain: '//pgvarns(j))
!        call pgmtxt('B',3.,0.5,0.5,'i')
      enddo
      call pgend
      if(ps2pdf.eq.1) then
         i = system('eps2pdf chains.eps >& /dev/null')
         i = system('rm -f chains.eps')
      endif
      endif !if(1.eq.2) then
      
      
      
      
      
      
      
      
      
      !Plot sigma values ('jump size')
      if(plchain.eq.1) then
      if(update.eq.0) write(6,'(A)')'Plotting sigma...'
      if(file.eq.0) then
        io = pgopen('16/xs')
        call pgsch(1.5)
      endif
      if(file.eq.1) then
        if(ps.eq.0) io = pgopen('sigs.ppm/ppm')
        if(ps.eq.1) io = pgopen('sigs.eps/cps')
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
      call pgscr(3,0.,0.5,0.)
      
      ic = 1
!      do j=1,12
      do j=2,npar
	call pgpage
        xmax = -1.e30
	ymin =  1.e30
        ymax = -1.e30
        do ic=1,nchains
          xmin = 0.
	  xmax = max(xmax,real(n(ic)))
	  dx = abs(xmax-xmin)*0.01
	  ymin = min(ymin,minval(sig(j,ic,2:n(ic))))
	  ymax = max(ymax,maxval(sig(j,ic,2:n(ic))))
	  dy = abs(ymax-ymin)*0.05
	enddo
        
	call pgswin(xmin-dx,xmax+dx,ymin-dy,ymax+dy)
        call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)
        
	do ic=1,nchains
          call pgsci(mod(ic*2,10))
!	  n1 = 1
!	  if(n(ic).gt.1000000) n1 = n(ic)/1000000
!          do i=0,n1-1
!	    call pgpoint(n(ic)/n1,is(n(ic)/n1*i+1:n(ic)/n1*(i+1)),sig(j,ic,n(ic)/n1*i+1:n(ic)/n1*(i+1)),1)
!	  enddo
!          call pgpoint(n(ic),is(1:n(ic)),sig(j,ic,1:n(ic)),1)
          do i=1,n(ic),chpli
            call pgpoint(1,is(i),sig(j,ic,i),1)
	  enddo
	enddo
        
        call pgsls(2)
        call pgsci(6)
        call pgline(2,(/real(nburn),real(nburn)/),(/-1.e20,1.e20/))
        call pgsci(1)
        call pgsls(1)
        call pgmtxt('T',1.,0.5,0.5,'Sigma: '//pgvarns(j))
!        call pgmtxt('B',3.,0.5,0.5,'i')
      enddo
      call pgend
      if(ps2pdf.eq.1) then
         i = system('eps2pdf sigs.eps >& /dev/null')
         i = system('rm -f sigs.eps')
      endif
      endif !if(1.eq.2) then
      
      
      
      
      
      
      
      
      
      !Plot acceptance rates
      if(plchain.eq.1) then
      if(update.eq.0) write(6,'(A)')'Plotting acceptance rates...'
      if(file.eq.0) then
        io = pgopen('17/xs')
        call pgsch(1.5)
      endif
      if(file.eq.1) then
        if(ps.eq.0) io = pgopen('accs.ppm/ppm')
        if(ps.eq.1) io = pgopen('accs.eps/cps')
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
      call pgscr(3,0.,0.5,0.)
      
      ic = 1
!      do j=1,12
      do j=2,npar
	call pgpage
        xmax = -1.e30
	ymin =  1.e30
        ymax = -1.e30
        do ic=1,nchains
          xmin = 0.
	  xmax = max(xmax,real(n(ic)))
	  dx = abs(xmax-xmin)*0.01
          do i=1,n(ic)
             if(acc(j,ic,i).gt.1.e-10 .and. acc(j,ic,i).lt.1.-1.e-10) then
                n0 = i
                exit
             end if
          end do
	  ymin = min(ymin,minval(acc(j,ic,n0:n(ic))))
	  ymax = max(ymax,maxval(acc(j,ic,n0:n(ic))))
	  dy = abs(ymax-ymin)*0.05
	enddo
        
        call pgsci(1)
        call pgsls(1)
	call pgswin(xmin-dx,xmax+dx,ymin-dy,ymax+dy)
        call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)
        
        
        call pgsci(3)
        call pgsls(2)
        call pgline(2,(/-1.e20,1.e20/),(/0.25,0.25/))
        call pgsci(6)
        call pgline(2,(/real(nburn),real(nburn)/),(/-1.e20,1.e20/))
        call pgsci(1)
        call pgsls(1)
        call pgmtxt('T',1.,0.5,0.5,'Acceptance: '//pgvarns(j))
        
	do ic=1,nchains
          call pgsci(mod(ic*2,10))
!	  n1 = 1
!	  if(n(ic).gt.1000000) n1 = n(ic)/1000000
!          do i=0,n1-1
!	    call pgpoint(n(ic)/n1,is(n(ic)/n1*i+1:n(ic)/n1*(i+1)),acc(j,ic,n(ic)/n1*i+1:n(ic)/n1*(i+1)),1)
!	  enddo
          do i=1,n(ic),chpli
            call pgpoint(1,is(i),acc(j,ic,i),1)
	  enddo
	enddo
      enddo
      call pgend
      if(ps2pdf.eq.1) then
         i = system('eps2pdf accs.eps >& /dev/null')
         i = system('rm -f accs.eps')
      endif
      endif !if(1.eq.2) then
      
      
      
      
      
      
      
      
      
      
      if(plpdf.eq.1) then
      if(update.eq.0) write(6,'(A)')'Plotting 1D pdfs...'
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
      call pgscr(3,0.,0.5,0.)
      
      ic=1
!      if(update.ne.1) write(6,'(7A10)')'var','y_max','x_max','median','model','diff','diff (%)'
      do j=2,npar
	call pgpage
        
        !Set x-ranges, bin the data and get y-ranges
        xmin = minval(alldat(j,1:nt))
        xmax = maxval(alldat(j,1:nt))
        dx = xmax - xmin
        x(1,1:nt) = alldat(j,1:nt)
        call bindata(nt,x(1,1:nt),1,nbin,xmin,xmax,xbin1,ybin1)
        
        
        ybin2 = ybin1
        if(1.eq.1) then
           i0 = nbin/10
           do i=1+i0,nbin+1-i0
              coefs1(1:2*i0+1) = ybin1(i-i0:i+i0)
              call savgol(coefs1(1:2*i0+1),2*i0+1,i0,i0,0,4)
              do i1=1,i0+1
                 !              print*,i1,i0-i1+2
                 coefs(i0-i1+2) = coefs1(i1)
              end do
              do i1 = i0+2,2*i0+1
                 !              print*,i1,3*i0+3-i1
                 coefs(3*i0+3-i1) = coefs1(i1)
              end do
              ybin2(i) = 0.
              do i1=1,2*i0+1
                 ybin2(i) = ybin2(i) + coefs(i1) * ybin1(i+i1-i0-1)
              end do
           end do
           ybin1 = ybin2
        end if
        
        
        xmin = xmin - 0.1*dx
        xmax = xmax + 0.1*dx
!        ymin = minval(ybin1(1:nbin))
        ymax = maxval(ybin1(1:nbin))*1.1
        
	if(file.eq.1) call pgsch(1.5)
	call pgswin(xmin,xmax,0.,ymax)
!        call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)  !After the plotting
!        write(6,'(4F7.4)')xmin,xmax,ymin,ymax
        
        
        
        call pgsci(1)
        if(file.eq.1.and.ps.eq.1) call pgslw(4)
        
        
        
!        call verthist(nbin,xbin1,ybin1)
        
        
        if(wrap(j).eq.0) then
           call pgsci(15)
           call pgpoly(nbin+2,(/xbin1(1),xbin1(1:nbin+1)/),(/0.,ybin1(1:nbin+1)/))
           call pgsci(1)
           call pgline(nbin+1,xbin1(1:nbin+1),ybin1(1:nbin+1))
           
           call pgline(2,(/xbin1(1),xbin1(1)/),(/0.,ybin1(1)/)) !Fix the loose ends
           call pgline(2,(/xbin1(nbin+1),xbin1(nbin+1)/),(/ybin1(nbin+1),0./))
        else
           call pgsci(15)
           shift = real(2*pi)
           if(changevar.eq.1) shift = 360.
           call pgpoly(nbin+3,(/xbin1(1),xbin1(1:nbin),xbin1(1)+shift,xbin1(1)+shift/),(/0.,ybin1(1:nbin),ybin1(1),0./))
           call pgsci(1)
           call pgline(nbin,xbin1(1:nbin),ybin1(1:nbin))
           
           call pgsls(4)
           if(file.eq.1.and.ps.eq.1) call pgslw(2)
           call pgline(nbin+1,(/xbin1(1:nbin)-shift,xbin1(1)/),(/ybin1(1:nbin),ybin1(1)/))
           call pgline(nbin,xbin1+shift,ybin1)
           call pgsls(1)
           if(file.eq.1.and.ps.eq.1) call pgslw(4)
           call pgline(2,(/xbin1(nbin),xbin1(1)+shift/),(/ybin1(nbin),ybin1(1)/))
        end if
        
        if(file.eq.1.and.ps.eq.1) call pgslw(1)
        
	!Plot median and model value
        call pgsci(6)
	if(file.eq.1) call pgsch(1.5)
        
!	if(update.ne.1) write(6,'(A10,5F10.4,F9.4,A1)')trim(varnames(j)),ysum(i1),xbin1(i1),median(j),pldat0(j),
!     &    abs(median(j)-pldat0(j)),abs(median(j)-pldat0(j))/pldat0(j)*100,'%'
!	print*,''
        if(j.eq.2.or.j.eq.3.or.j.eq.5.or.j.eq.9) then
          write(str,'(A,F7.3,A7,F7.3,A11,F6.2,A1)')trim(pgvarns(j))//':  mdl: ',pldat0(j),'  med: ',median(j),
     &                                          '  \(2030): ',abs(median(j)-pldat0(j))/pldat0(j)*100,'%'
        else
          write(str,'(A,F7.3,A7,F7.3,A11,F7.3)')trim(pgvarns(j))//':  mdl: ',pldat0(j),'  med: ',median(j),
     &                                          '  \(2030): ',abs(median(j)-pldat0(j))
        endif
        if(file.eq.1.and.ps.eq.1) call pgslw(2)
        call pgsci(6)
        call pgsls(2)
        call pgline(2,(/pldat0(j),pldat0(j)/),(/-1.e20,1.e20/))
        call pgsls(4)
        call pgline(2,(/pldat1(j),pldat1(j)/),(/-1.e20,1.e20/))
        call pgsci(1)
        call pgsls(2)
        call pgline(2,(/median(j),median(j)/),(/-1.e20,1.e20/))
        call pgsls(4)
        call pgline(2,(/ranges(1,j,1),ranges(1,j,1)/),(/-1.e20,1.e20/))
        call pgline(2,(/ranges(1,j,2),ranges(1,j,2)/),(/-1.e20,1.e20/))
        call pgline(2,(/ranges(1,j,3),ranges(1,j,3)/),(/-1.e20,1.e20/))
        call pgsls(1)
        if(file.eq.1.and.ps.eq.1) call pgslw(1)
	if(file.eq.1) call pgsch(1.2)
!        call pgmtxt('B',3.,0.5,0.5,varnames(j))
        call pgmtxt('T',1.,0.5,0.5,trim(str))
	if(file.eq.1) call pgsch(1.5)
        call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)
      enddo !j
      
      call pgsubp(1,1)
      call pgswin(-1.,1.,-1.,1.)
      call pgslw(3)

      
      !Make sure gv auto-reloads on change
!      if(file.eq.1.and.ps.eq.1) then
!        call pgpage
!        call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)
!      endif

      call pgend
      if(ps2pdf.eq.1) then
         i = system('eps2pdf pdfs.eps >& /dev/null')
         i = system('rm -f pdfs.eps')
      endif
      endif !if(plpdf.eq.1) then
      
      
      
      
      
      
      
      
      
      
      
      
      if(plpdf2d.eq.1) then
      if(update.eq.0) write(6,'(A)')'Plotting 2D pdfs...'
      if(file.eq.0) then
        io = pgopen('15/xs')
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
      call pgscr(3,0.,0.5,0.)
      
        !Columns in dat(): 1:logL, 2:Mc, 3:eta, 4:tc, 5:dl, 6:sinlati, 7longi:, 8:phase, 9:spin, 10:kappa, 11:thJ0, 12:phJ0, 13:alpha
      do j1=2,12
      do j2=j1+1,npar
!      do j1=2,13
!      do j2=2,13
!      if(j1.eq.j2.cycle)
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
        call pgline(2,(/median(j1),median(j1)/),(/-1.e20,1.e20/))
        call pgline(2,(/-1.e20,1.e20/),(/median(j2),median(j2)/))
        call pgsls(1)
        if(file.eq.1.and.ps.eq.1) call pgslw(1)
	if(file.eq.1) call pgsch(1.5)
        
        
        call pgsls(1)
        call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)
        
        call pgmtxt('B',2.5,0.5,0.5,pgvarns(j1))
        call pgmtxt('L',2.5,0.5,0.5,pgvarns(j2))
        
        call pgpage
      enddo !j2
      enddo !j1
        
        
      
      
      
!      !Make sure gv auto-reloads on change
!      if(1.eq.2.and.file.eq.1.and.ps.eq.1) then
!        call pgpage
!        call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)
!      endif

      call pgend
      endif !if(plpdf.eq.1) then
      
      
      
      
      
      
      
      
      
      
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
      ! xbin, ybin - output: binned data (x, y).  The x values are the left side of the bin!
      
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
!        xbin(k) = xmin + (real(k)-0.5)*dx  !x is the centre< of the bin
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
      
      
      
      
      
      
************************************************************************************************************************************
      SUBROUTINE dindexx(n,arr,indx)
************************************************************************************************************************************
      INTEGER :: n,indx(n),M,NSTACK
      REAL*8 :: arr(n),a
      PARAMETER (M=7,NSTACK=50)
      INTEGER :: i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      do 11 j=1,n
        indx(j)=j
11    continue
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do 12 i=j-1,l,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=l-1
2         indx(i+1)=indxt
13      continue
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
        endif
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(l+1)))then
          itemp=indx(l)
          indx(l)=indx(l+1)
          indx(l+1)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l+1)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l+1)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(jstack.gt.NSTACK)pause 'NSTACK too small in indexx'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END
************************************************************************************************************************************
      
      
************************************************************************************************************************************
      SUBROUTINE rindexx(n,arr,indx)
************************************************************************************************************************************
      INTEGER :: n,indx(n),M,NSTACK
      REAL :: arr(n),a
      PARAMETER (M=7,NSTACK=50)
      INTEGER :: i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      do 11 j=1,n
        indx(j)=j
11    continue
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do 12 i=j-1,l,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=l-1
2         indx(i+1)=indxt
13      continue
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
        endif
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(l+1)))then
          itemp=indx(l)
          indx(l)=indx(l+1)
          indx(l+1)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l+1)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l+1)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(jstack.gt.NSTACK)pause 'NSTACK too small in indexx'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END
************************************************************************************************************************************


************************************************************************************************************************************
      SUBROUTINE savgol(c,np,nl,nr,ld,m)
************************************************************************************************************************************
      INTEGER :: ld,m,nl,np,nr,MMAX
      REAL :: c(np)
      PARAMETER (MMAX=6)
CU    USES lubksb,ludcmp
      INTEGER :: imj,ipj,j,k,kk,mm,indx(MMAX+1)
      REAL :: d,fac,sum,a(MMAX+1,MMAX+1),b(MMAX+1)
      if(np.lt.nl+nr+
     *1.or.nl.lt.0.or.nr.lt.0.or.ld.gt.m.or.m.gt.MMAX.or.nl+nr.lt.m) 
     *pause 'bad args in savgol'
      do 14 ipj=0,2*m
        sum=0.
        if(ipj.eq.0)sum=1.
        do 11 k=1,nr
          sum=sum+float(k)**ipj
11      continue
        do 12 k=1,nl
          sum=sum+float(-k)**ipj
12      continue
        mm=min(ipj,2*m-ipj)
        do 13 imj=-mm,mm,2
          a(1+(ipj+imj)/2,1+(ipj-imj)/2)=sum
13      continue
14    continue
      call ludcmp(a,m+1,MMAX+1,indx,d)
      do 15 j=1,m+1
        b(j)=0.
15    continue
      b(ld+1)=1.
      call lubksb(a,m+1,MMAX+1,indx,b)
      do 16 kk=1,np
        c(kk)=0.
16    continue
      do 18 k=-nl,nr
        sum=b(1)
        fac=1.
        do 17 mm=1,m
          fac=fac*k
          sum=sum+b(mm+1)*fac
17      continue
        kk=mod(np-k,np)+1
        c(kk)=sum
18    continue
      return
      END
************************************************************************************************************************************


************************************************************************************************************************************
      SUBROUTINE lubksb(a,n,np,indx,b)
************************************************************************************************************************************
      INTEGER :: n,np,indx(n)
      REAL :: a(np,np),b(n)
      INTEGER :: i,ii,j,ll
      REAL :: sum
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return
      END
************************************************************************************************************************************
      
      
      
************************************************************************************************************************************
      SUBROUTINE ludcmp(a,n,np,indx,d)
************************************************************************************************************************************
      INTEGER :: n,np,indx(n),NMAX
      REAL :: d,a(np,np),TINY
      PARAMETER (NMAX=500,TINY=1.0e-20)
      INTEGER :: i,imax,j,k
      REAL :: aamax,dum,sum,vv(NMAX)
      d=1.
      do 12 i=1,n
        aamax=0.
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
        if (aamax.eq.0.) pause 'singular matrix in ludcmp'
        vv(i)=1./aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.
        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.)a(j,j)=TINY
        if(j.ne.n)then
          dum=1./a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      return
      END
************************************************************************************************************************************
