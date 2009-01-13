! Compute statistics for analysemcmc

subroutine statistics(exitcode)
  use constants
  use analysemcmc_settings
  use general_data
  use stats_data
  use chain_data
  use mcmcrun_data
  implicit none
  integer :: c,i,ic,i0,j,j1,o,p,p1,p2,nr,nstat,exitcode,io
  integer :: indexx(npar1,nchs*narr1),index1(nchs*narr1),parr(npar1)
  integer :: nn,lowvar(npar1),nlowvar,highvar(npar1),nhighvar,ntotrelvar,nrhat
  real :: rev2pi,x0,x1,x2,y1,y2,dx,x,rrevpi
  real :: range1,minrange,maxgap,ival,wrapival,centre,maxlogl,minlogl,shival,shival2
  real :: medians(npar1),mean(npar1),var1(npar1),var2(npar1),corr,corr1,corr2
  real*8 :: chmean(nchs,npar1),totmean(npar1),chvar(npar1),chvar1(nchs,npar1),totvar(npar1),totrhat,totrelvar
  real*8 :: var,total
  !real*16 :: var,total  !gfortran doesn't support this
  character :: ch,url*99,gps*19,xs10*10,xs20*20
  
  !K-S test:
  !integer :: nn1,nlogl1,nlogl2,ksn1,ksn2
  !real*8 :: ksdat1(narr1),ksdat2(narr1),ksd,ksprob
  
  exitcode = 0
  
  
  !Check which parameters were fixed during the MCMC run
  fixedpar = 0
  do ic=1,nchains
     do p=par1,par2
        if( abs(minval(alldat(ic,p,5:n(ic))) - maxval(alldat(ic,p,5:n(ic))) ) .lt. 1.d-6) fixedpar(p) = 1  !Doesn't matter in which chain this happens
     end do
  end do
  nfixedpar = sum(fixedpar)
  
  
  !Sort all data and find the interval limits for the default probability interval for the wrapable parameters
  if(prprogress.ge.2) write(*,*)''
  !ivals(nival+1) = 1. !The last probability interval is always 100% !Is this necessary?
  shift = 0.
  wrap = 0
  rashift = 0.
  do ic=1,nchains
     if(mergechains.eq.0.and.contrchain(ic).eq.0) cycle
     !wrapival = ivals(nival) !Use the largest range
     wrapival = 0.999 !Always use a very large range (?)
     indexx = 0
     if(prprogress.ge.2.and.mergechains.eq.0) write(*,'(A,I2.2,A,$)')' Ch',ic,' '
     if(prprogress.ge.2.and.ic.eq.1.and.wrapdata.ge.1) write(*,'(A,$)')'  Wrap data. '
     do p=par1,par2
        if(wrapdata.eq.0 .or. (p.ne.8.and.p.ne.10.and.p.ne.12.and.p.ne.13)) then
           call rindexx(n(ic),alldat(ic,p,1:n(ic)),index1(1:n(ic)))
           indexx(p,1:n(ic)) = index1(1:n(ic))
           cycle
        end if
        
        !Columns in dat(): 1:logL 2:mc, 3:eta, 4:tc, 5:logdl, 6:spin, 7:kappa, 8: RA, 9:sindec,10:phase, 11:sinthJ0, 12:phiJ0, 13:alpha
        !Make sure data are between 0 and 2pi or between 0 and pi to start with:
        do i=1,n(ic)
           if(p.eq.12) then
              alldat(ic,p,i) = rrevpi(alldat(ic,p,i)) !Pol.angle
           else
              alldat(ic,p,i) = rev2pi(alldat(ic,p,i)) 
           end if
        end do
        call rindexx(n(ic),alldat(ic,p,1:n(ic)),index1(1:n(ic)))
        indexx(p,1:n(ic)) = index1(1:n(ic))
        
        minrange = 1.e30
        do i=1,n(ic)
           x1 = alldat(ic,p,indexx(p,i))
           x2 = alldat(ic,p,indexx(p,mod(i+nint(n(ic)*wrapival)-1,n(ic))+1))
           if(p.eq.12) then
              range1 = mod(x2 - x1 + real(10*pi),rpi)
           else
              range1 = mod(x2 - x1 + real(20*pi),rtpi)
           end if
           if(range1.lt.minrange) then
              minrange = range1
              y1 = x1
              y2 = x2
              !write(*,'(2I6,7F10.5)')i,mod(nint(i+n(ic)*wrapival),n(ic)),x1,x2,range1,minrange,y1,y2,(y1+y2)/2.
           end if
           !write(*,'(2I6,7F10.5)')i,mod(nint(i+n(ic)*wrapival),n(ic)),x1,x2,range1,minrange,y1,y2,(y1+y2)/2.
        end do !i
        centre = (y1+y2)/2.
        
        !Define shift interval:
        if(p.eq.12) then
           shival = rpi   !Shift interval
           shival2 = rpi2 !Shift interval/2
        else
           shival = rtpi   !Shift interval
           shival2 = rpi   !Shift interval/2
        end if
        
        if(y1.gt.y2) then
           if(p.eq.12) then
              wrap(ic,p) = 2  !Wrap over pi
              !centre = mod(centre + rpi2, rpi) !Then distribution peaks close to 0/pi, shift centre by pi/2
           else
              wrap(ic,p) = 1  !Wrap over 2pi
              !centre = mod(centre + rpi, rtpi) !Then distribution peaks close to 0/2pi, shift centre by pi
           end if
           centre = mod(centre + shival2, shival) !Then distribution peaks close to 0/shival, shift centre by shival/2
        end if
        !if(p.eq.8) write(*,'(3I6,8F10.5)')ic,p,wrap(ic,p), wrapival,y1*r2h,y2*r2h,minrange*r2h,centre*r2h
        !if(p.eq.12) write(*,'(3I6,8F10.5)')ic,p,wrap(ic,p), wrapival,y1,y2,minrange,centre
        
        
        
        !See whether there's a gap in the data.  WHY is this necessary, should it work like this???
        if(wrap(ic,p).eq.0 .and. 1.eq.2) then
           !ymin = minval(alldat(ic,p,1:n(ic)))
           !ymax = maxval(alldat(ic,p,1:n(ic)))
           i0 = -1
           maxgap = -1.e30
           !write(*,'(2I3,I8)')ic,p,n(ic)
           do i=1,n(ic)-1
              x1 = alldat(ic,p,indexx(p,i))
              x2 = alldat(ic,p,indexx(p,i+1))
              !write(*,'(2I3,2I8,4F10.5)')ic,p,i,i0,x1,x2,x2-x2,maxgap
              if(x2-x1.gt.maxgap) then
                 maxgap = x2-x1
                 i0 = i
              end if
           end do !i
           x1 = alldat(ic,p,indexx(p,i0))
           x2 = alldat(ic,p,indexx(p,i0+1))
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
        
        !For the special case of 2pi:
        !if(wrap(ic,p).eq.1) shift(ic,p) = rtpi - mod(centre + rpi, rtpi)
        !if(p.eq.8.and.ic.eq.1) rashift = shift(ic,p)                         !Save RA shift to plot sky map
        !alldat(ic,p,1:n(ic))   = mod(alldat(ic,p,1:n(ic))   + shift(ic,p), rtpi) - shift(ic,p)
        !pldat(ic,p,1:ntot(ic)) = mod(pldat(ic,p,1:ntot(ic)) + shift(ic,p), rtpi) - shift(ic,p) !Original data
        !y1 = mod(y1 + shift(ic,p), rtpi) - shift(ic,p)
        !y2 = mod(y2 + shift(ic,p), rtpi) - shift(ic,p)
        !centre = mod(centre + shift(ic,p), rtpi) - shift(ic,p)
        
        !For the general case of shival (= 2pi or pi)
        if(wrap(ic,p).gt.0) shift(ic,p) = shival - mod(centre + shival2, shival)
        if(p.eq.8.and.ic.eq.1) rashift = shift(ic,p)                         !Save RA shift to plot sky map
        alldat(ic,p,1:n(ic))   = mod(alldat(ic,p,1:n(ic))   + shift(ic,p), shival) - shift(ic,p)
        pldat(ic,p,1:ntot(ic)) = mod(pldat(ic,p,1:ntot(ic)) + shift(ic,p), shival) - shift(ic,p)   !Original data
        y1 = mod(y1 + shift(ic,p), shival) - shift(ic,p)
        y2 = mod(y2 + shift(ic,p), shival) - shift(ic,p)
        centre = mod(centre + shift(ic,p), shival) - shift(ic,p)
        
        minrange = y2-y1
        !call rindexx(n(ic),alldat(ic,p,1:n(ic)),indexx(p,1:n(ic)))  !Re-sort
        call rindexx(n(ic),alldat(ic,p,1:n(ic)),index1(1:n(ic)))  !Re-sort
        indexx(p,1:n(ic)) = index1(1:n(ic))
        !if(p.eq.8) write(*,'(I3,A8,4x,6F10.5,I4)')ic,varnames(p),y1,y2,minrange,centre,minval(alldat(ic,p,1:n(ic))),maxval(alldat(ic,p,1:n(ic))),wrap(ic,p)
        
        if(abs( abs( minval(alldat(ic,p,1:n(ic))) - maxval(alldat(ic,p,1:n(ic))) ) - shival) .lt. 1.e-3)  wrap(ic,p) = 1   !If centre is around shival2, still needs to be flagged 'wrap' to plot PDF
     end do !p
     
     
     
     !Do statistics
     !if(prprogress.ge.2) write(*,'(A)')' Calculating: statistics...'
     if(prprogress.ge.1.and.ic.eq.1) write(*,'(A,$)')'  Calc: stats, '
     do p=par1,par2
        !Determine the median
        if(mod(n(ic),2).eq.0) medians(p) = 0.5*(alldat(ic,p,indexx(p,n(ic)/2)) + alldat(ic,p,indexx(p,n(ic)/2+1)))
        if(mod(n(ic),2).eq.1) medians(p) = alldat(ic,p,indexx(p,(n(ic)+1)/2))
        
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
        if(prprogress.ge.1) write(*,'(A,$)')' corrs, '
        do p1=par1,par2
           do p2=par1,par2
           !do p2=p1,par2
              corrs(p1,p2) = 0.
              if(fixedpar(p1)+fixedpar(p2).eq.0) then
                 do i=1,n(ic)
                    !corrs(p1,p2) = corrs(p1,p2) + (alldat(ic,p1,i) - medians(p1))*(alldat(ic,p2,i) - medians(p2))  !Use median
                    corrs(p1,p2) = corrs(p1,p2) + (alldat(ic,p1,i) - mean(p1))*(alldat(ic,p2,i) - mean(p2)) !Use mean; hardly differs from median method
                 end do
                 !corrs(p1,p2) = corrs(p1,p2) / (stdev1(p1)*stdev1(p2)*(n(ic)-1))  !Use median
                 corrs(p1,p2) = corrs(p1,p2) / (stdev2(p1)*stdev2(p2)*(n(ic)-1))  !Use mean
              end if
           end do !p2
        end do !p1
     end if
     
     
     !Autocorrelations:
     if(placorr.gt.0) then
        !write(*,'(A)')' Calculating autocorrelations...'
        if(prprogress.ge.1) write(*,'(A,$)')' autocorrs, '
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
     if(prprogress.ge.1.and.ic.eq.1) write(*,'(A,$)')' prob.ivals: '
     c0 = 0
     do c=1,nival
        ival = ivals(c)
        c0 = ival0
        !if(c.ne.c0.and.prival.lt.2.and.savestats.eq.0) cycle
        if(c.ne.c0 .and. prival.eq.0 .and. prstat.lt.2 .and. savestats.eq.0) cycle
        
        if(prprogress.ge.1.and.ic.eq.1) write(*,'(F6.3,$)')ival
        !print*,par1,par2
        do p=par1,par2
        !do p=2,2
           !print*,p,minval(alldat(ic,p,1:n(ic))),maxval(alldat(ic,p,1:n(ic)))
           minrange = 1.e30
           !write(*,'(A8,4x,4F10.5,I4)')varnames(p),y1,y2,minrange,centre,wrap(ic,p)
           do i=1,floor(n(ic)*(1.-ival))
              x1 = alldat(ic,p,indexx(p,i))
              x2 = alldat(ic,p,indexx(p,i+floor(n(ic)*ival)))
              range1 = abs(x2 - x1)
              !range1 = x2 - x1
              if(range1.lt.minrange) then
                 minrange = range1
                 y1 = x1
                 y2 = x2
              end if
              !write(*,'(I6,7F10.5)')i,x1,x2,range1,minrange,y1,y2,(y1+y2)/2.
           end do
           centre = (y1+y2)/2.
           !write(*,'(A8,4x,4F10.5,I4)')varnames(p),y1,y2,minrange,centre,wrap(ic,p)
           
           !Save ranges:
           nr = 4                  !Only ranges(:,:,:,1:nr) get converted
           ranges(ic,c,p,1) = y1
           ranges(ic,c,p,2) = y2
           ranges(ic,c,p,3) = centre
           ranges(ic,c,p,4) = y2-y1
           ranges(ic,c,p,5) = ranges(ic,c,p,4)
           if(p.eq.2.or.p.eq.3.or.p.eq.5.or.p.eq.6.or.p.eq.14.or.p.eq.15) ranges(ic,c,p,5) = ranges(ic,c,p,4)/ranges(ic,c,p,3)
        end do !p
     end do !c
     !if(prprogress.ge.2) write(*,'(A34,F8.4)')'.  Standard probability interval: ',ivals(ival0)
     !if(prprogress.ge.2) write(*,'(A,F8.4,$)')', default ival:',ivals(ival0)
     
     
     
     
     
     
     
     
     
     !Compute Bayes factor
     !if(prprogress.ge.1.and.ic.eq.1) write(*,'(A,$)')'  Bayes factor,'
     p=1  !Likelihood/Posterior
     total = 0
     maxlogl = -1.e30
     minlogl =  1.e30
     do i=1,n(ic)
        var = alldat(ic,p,i)          !Use quadruple precision
        total = total + exp(-var)
        !write(*,'(2ES15.5)')var,exp(-var)
        maxlogl = max(alldat(ic,p,i),maxlogl)
        minlogl = min(alldat(ic,p,i),minlogl)
     end do
     var = dble(n(ic))/total
     log10bayesfactor(ic) = real(log10(var))
     logebayesfactor(ic) = real(log(var))
     !write(*,'(2F10.3,I9)')maxlogl,minlogl,n(ic)
     
          
     
     
     
     
     
     
     
     
     
     
     !**********************************************************************************************
     !******   CHANGE VARIABLES   ******************************************************************
     !**********************************************************************************************

     
     
     
     !Change variables
     !Columns in alldat(): 1:logL, 2:Mc, 3:eta, 4:tc, 5:logdl,   6:longi, 7:sinlati:, 8:phase, 9:spin,   10:kappa,     11:sinthJ0, 12:phiJ0, 13:alpha
     if(changevar.eq.1) then
        if(prprogress.ge.1.and.ic.eq.i.and.update.eq.0) write(*,'(A,$)')'.  Change vars. '
        do p=par1,par2
           if(p.eq.5) then !Distance
              stdev1(p) = stdev1(p)*exp(stats(ic,p,1))  !Median  For exponential function y = exp(x), sig_y = exp(x) sig_x
              stdev2(p) = stdev2(p)*exp(stats(ic,p,2))  !Mean
              alldat(ic,p,1:n(ic)) = exp(alldat(ic,p,1:n(ic)))     !logD -> Distance
              if(ic.eq.1) startval(1:nchains0,p,1:3) = exp(startval(1:nchains0,p,1:3))
              stats(ic,p,1:nstat) = exp(stats(ic,p,1:nstat))
              ranges(ic,1:nival,p,1:nr) = exp(ranges(ic,1:nival,p,1:nr))
           end if
           if(p.eq.7) then !Kappa -> theta_SL
              stdev1(p) = abs(-1./(sqrt(max(1.-stats(ic,p,1)**2,1.e-30))) * stdev1(p))*rr2d  !Based on median
              stdev2(p) = abs(-1./(sqrt(max(1.-stats(ic,p,2)**2,1.e-30))) * stdev2(p))*rr2d  !Based on mean
              alldat(ic,p,1:n(ic)) = acos(alldat(ic,p,1:n(ic)))*rr2d
              if(ic.eq.1) startval(1:nchains0,p,1:3) = acos(startval(1:nchains0,p,1:3))*rr2d
              stats(ic,p,1:nstat) = acos(stats(ic,p,1:nstat))*rr2d
              ranges(ic,1:nival,p,1:nr) = acos(ranges(ic,1:nival,p,1:nr))*rr2d
              do c=1,nival
                 y1 = ranges(ic,c,p,2)
                 ranges(ic,c,p,2) = ranges(ic,c,p,1)  !acos is monotonously decreasing
                 ranges(ic,c,p,1) = y1
              end do
           end if
           if(p.eq.8) then !RA
              stdev1(p) = stdev1(p)*rr2h
              stdev2(p) = stdev2(p)*rr2h
              alldat(ic,p,1:n(ic)) = alldat(ic,p,1:n(ic))*rr2h
              if(ic.eq.1) startval(1:nchains0,p,1:3) = startval(1:nchains0,p,1:3)*rr2h
              stats(ic,p,1:nstat) = stats(ic,p,1:nstat)*rr2h
              ranges(ic,1:nival,p,1:nr) = ranges(ic,1:nival,p,1:nr)*rr2h
           end if
           !if(p.eq.9.or.p.eq.11) then  !Declination or theta_Jo
           if(p.eq.9) then  !Declination                                                               
              stdev1(p) = abs(1./(sqrt(max(1.-stats(ic,p,1)**2,1.e-30))) * stdev1(p))*rr2d  !Based on median
              stdev2(p) = abs(1./(sqrt(max(1.-stats(ic,p,2)**2,1.e-30))) * stdev2(p))*rr2d  !Based on mean
              alldat(ic,p,1:n(ic)) = asin(alldat(ic,p,1:n(ic)))*rr2d
              if(ic.eq.1) startval(1:nchains0,p,1:3) = asin(startval(1:nchains0,p,1:3))*rr2d
              stats(ic,p,1:nstat) = asin(stats(ic,p,1:nstat))*rr2d
              ranges(ic,1:nival,p,1:nr) = asin(ranges(ic,1:nival,p,1:nr))*rr2d
           end if
           !if(p.eq.10.or.p.eq.12.or.p.eq.13) then  !hi_c, phi_Jo, alpha_c
           if(p.ge.10.and.p.le.13) then  !Phi_c, incl, pol.ang, alpha_c
              stdev1(p) = stdev1(p)*rr2d
              stdev2(p) = stdev2(p)*rr2d
              alldat(ic,p,1:n(ic)) = alldat(ic,p,1:n(ic))*rr2d
              if(ic.eq.1) startval(1:nchains0,p,1:3) = startval(1:nchains0,p,1:3)*rr2d
              stats(ic,p,1:nstat) = stats(ic,p,1:nstat)*rr2d
              ranges(ic,1:nival,p,1:nr) = ranges(ic,1:nival,p,1:nr)*rr2d
           end if
           ranges(ic,1:nival,p,3) = 0.5*(ranges(ic,1:nival,p,1) + ranges(ic,1:nival,p,2))
           ranges(ic,1:nival,p,4) = ranges(ic,1:nival,p,2) - ranges(ic,1:nival,p,1)
           ranges(ic,1:nival,p,5) = ranges(ic,1:nival,p,4)
           if(p.eq.2.or.p.eq.3.or.p.eq.5.or.p.eq.6.or.p.eq.14.or.p.eq.15) ranges(ic,1:nival,p,5) = ranges(ic,1:nival,p,5)/ranges(ic,1:nival,p,3)
        end do !p
        
        !Columns in dat(): 1:logL 2:mc, 3:eta, 4:tc, 5:dl, 6:spin,  7:theta_SL, 8: RA,   9:dec, 10:phase, 11:incl, 12:pol.ang., 13:alpha
        varnames(1:15) = (/'logL','Mc','eta','tc','dl','spin','th_SL','RA','Dec','phase','incl','polang','alpha','M1','M2'/)
        pgvarns(1:15)  = (/'log Likelihood        ','M\dc\u (M\d\(2281)\u) ','\(2133)               ','t\dc\u (s)            ',  &
             'd\dL\u (Mpc)          ','a\dspin\u             ','\(2134)\dSL\u(\(2218))','R.A. (h)              ','Dec. (\(2218))        ', &
             '\(2147)\dc\u (\(2218))','acos(J\(2236)N)','\(2149)\dJ\u (\(2218))','\(2127)\dc\u (\(2218))','M\d1\u (M\d\(2281)\u) ','M\d2\u(M\d\(2281)\u)  '/)
        pgvarnss(1:15)  = (/'log L','M\dc\u','\(2133)','t\dc\u','d\dL\u','a\dspin\u','\(2134)\dSL\u','R.A.','Dec.','\(2147)\dc\u',  &
             'acos(J\(2236)N)','\(2149)\dJ0\u','\(2127)\dc\u','M\d1\u','M\d2\u'/)
        !Include units
        pgvarnss(1:15)  = (/'log L','M\dc\u (M\d\(2281)\u)','\(2133)','t\dc\u (s)','d\dL\u (Mpc)','a\dspin\u','\(2134)\dSL\u (\(2218))','R.A. (h)','Dec. (\(2218))','\(2147)\dc\u (\(2218))',  &
             'acos(J\(2236)N)','\(2149)\dJ0\u (\(2218))','\(2127)\dc\u (\(2218))','M\d1\u (M\d\(2281)\u)','M\d2\u (M\d\(2281)\u)'/)
        !Units only
        pgunits(1:15)  = (/'','M\d\(2281)\u ','','s','Mpc','','\(2218)','\uh\d','\(2218)','\(2218)','','\(2218)','\(2218)','M\d\(2281)\u','M\d\(2281)\u'/)
     end if !if(changevar.eq.1)
     
     
     
     
     
     
     
     
     !Find 100% probability range
     do c = 1,nival
        if(abs(ivals(c)-1.).lt.1.e-4) then !Then treat it as a 100% interval to prevent numerical problems
           if(prprogress.ge.1) write(*,'(A,F9.4,A)')'  Treating probability interval',ivals(c)*100,'% as 100%'
           do p=par1,par2
              if(p.eq.1) cycle
              ranges(ic,c,p,1) = minval(alldat(ic,p,1:n(ic)))
              ranges(ic,c,p,2) = maxval(alldat(ic,p,1:n(ic)))
              ranges(ic,c,p,3) = 0.5*(ranges(ic,c,p,1) + ranges(ic,c,p,2))
              ranges(ic,c,p,4) = ranges(ic,c,p,2) - ranges(ic,c,p,1)
              ranges(ic,c,p,5) = ranges(ic,c,p,4)
              if(p.eq.2.or.p.eq.3.or.p.eq.5.or.p.eq.6.or.p.eq.14.or.p.eq.15) ranges(ic,c,p,5) = ranges(ic,c,p,5)/ranges(ic,c,p,3)
           end do
        end if
     end do
     
     if(prprogress.ge.1) then
        if(ic.eq.nchains) then
           write(*,*)
        else
           write(*,'(A,$)')'  '
        end if
     end if
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     !**********************************************************************************************
     !******   PRINT STATISTICS   ******************************************************************
     !**********************************************************************************************

     o=6 !Print to this unit - 6=screen
     
     
     if(prprogress+prstat+prival+prconv.gt.0.and.ic.eq.1) then
        write(o,'(/,A,2(A,F6.2),$)')'  Bayes factor:   ','log_e(B_SN) =',logebayesfactor(ic),',  log_10(B_SN) =',log10bayesfactor(ic)
        write(o,'(8x,A,2(A,F6.2))')'  Maximum likelihood:   ','log_e(Lmax) =',startval(ic,1,3),',  log_10(Lmax) =',startval(ic,1,3)/log(10.)
     end if
     
     
     !Print statistics to screen
     if(prstat.gt.0) then
        write(o,'(/,A)')'  Main statistics:'
        do c=1,nival
           if(c.ne.c0.and.prstat.lt.2) cycle
           if(c.gt.1.and.prstat.ge.2) write(o,*)
           write(o,'(A10, A12,2A10,A12, 4A8, 4A10,A8,A10, A4,A12,F7.3,A2)')'param.','model','median','mean','Lmax','stdev1','stdev2','abvar1','abvar2',  &
                'rng_c','rng1','rng2','drng','d/drng','delta','ok?','result (',ivals(c)*100,'%)'
           do p=par1,par2
              if(fixedpar(p).eq.1) cycle !Parameter was not fitted
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
        end do
     end if
     
     
     !Print intervals as: centre, delta, in range:
     if(prival.eq.1.or.prival.eq.3) then
        write(o,'(/,A)')'  Probability intervals:'
        write(o,'(A22,A8,$)')'Interval:',''
        do c=1,nival
           write(o,'(F20.4,A9,$)')ivals(c),''
        end do
        write(o,*)''
        
        write(o,'(A10,2x,2A9,$)')'param.','model','median'
        do c=1,nival
           !write(o,'(2x,2A9,A8,$)')'rng1','rng2','in rnge'
           write(o,'(2x,3A9,$)')'centre','delta','in rnge'
        end do
        write(o,*)''
        do p=par1,par2
           !if(stdev1(p).lt.1.d-20) cycle
           if(fixedpar(p).eq.1) cycle
           write(o,'(A10,2x,2F9.4,$)')varnames(p),startval(ic,p,1),stats(ic,p,1)
           do c=1,nival
              if(mergechains.eq.0) then
                 write(o,'(2x,2F9.4,F6.3,$)')ranges(ic,c,p,3),ranges(ic,c,p,4),min(2*abs(startval(ic,p,1)-ranges(ic,c,p,3))/ranges(ic,c,p,4),9.999) !Defined with centre of prob. range, need some extra security to print correctly
              else
                 write(o,'(2x,2F9.4,F6.3,$)')ranges(ic,c,p,3),ranges(ic,c,p,4),2*abs(startval(ic,p,1)-ranges(ic,c,p,3))/ranges(ic,c,p,4) !Defined with centre of prob. range
              end if
              if(startval(ic,p,1).gt.ranges(ic,c,p,1).and.startval(ic,p,1).lt.ranges(ic,c,p,2)) then
                 write(o,'(A3,$)')'y '
              else
                 write(o,'(A3,$)')'N*'
              end if
           end do
           write(o,*)''
        end do
        write(o,*)''
     end if
     
     
     if(prival.ge.2) then
        write(o,'(/,A)')'  Statistics and probability intervals:'
        
        !Print intervals as x +- dx:
        write(o,'(A65,$)')'Interval:'
        do c=1,nival
           write(o,'(F16.3,A1,A9,$)')ivals(c)*100,'%',''
        end do
        write(o,*)''
        
        write(o,'(A11,1x,4A11,$)')'Parameter','median','mean','Lmax','stdev'
        do c=1,nival
           write(o,'(5x,A9,4x,A8,$)')'x','dx'
        end do
        write(o,*)''
        do p=max(par1,2),par2  !Leave out logL
           !if(stdev1(p).lt.1.d-20) cycle
           if(fixedpar(p).eq.1) cycle
           write(o,'(A10,2x,4F11.4,$)')varnames(p),stats(ic,p,1),stats(ic,p,2),startval(ic,p,3),stdev2(p)
           do c=1,nival
              write(o,'(5x,F9.4,A4,F8.4$)')ranges(ic,c,p,3),' +- ',0.5d0*ranges(ic,c,p,4)
           end do
           write(o,*)''
        end do
        write(o,*)''
        
        !Print intervals as low-high:
        write(o,'(A65,$)')'Interval:'
        do c=1,nival
           write(o,'(F16.3,A1,A9,$)')ivals(c)*100,'%',''
        end do
        write(o,*)''
        
        write(o,'(A11,1x,4A11,$)')'Parameter','median','mean','Lmax','stdev'
        do c=1,nival
           write(o,'(6x,A9,3x,A8,$)')'min','max'
        end do
        write(o,*)''
        do p=max(par1,2),par2  !Leave out logL
           !if(stdev1(p).lt.1.d-20) cycle
           if(fixedpar(p).eq.1) cycle
           write(o,'(A10,2x,4F11.4,$)')varnames(p),stats(ic,p,1),stats(ic,p,2),startval(ic,p,3),stdev2(p)
           do c=1,nival
              write(o,'(6x,F9.4,A3,F8.4$)')ranges(ic,c,p,1),' - ',ranges(ic,c,p,2)
           end do
           write(o,*)''
        end do
        write(o,*)''
     end if
     
     
     
     
     
     
     
     
     !Print correlations:
     if(prcorr.gt.0) then
        corr1 = 0.1
        corr2 = 0.5
        write(o,'(/,A,$)')'  Correlations  '
        write(o,'(A,3(F4.2,A))')'  (weak [',corr1,'<abs(cor)<',corr2,']: in lower triangle,  strong [abs(cor)>',corr2,']: in upper triangle):'
        write(o,'(A8,$)')''
        do p=par1,par2
           !if(stdev1(p).gt.1.d-20) write(o,'(A7,$)')trim(varnames(p))
           if(fixedpar(p).eq.0) write(o,'(A7,$)')trim(varnames(p))
        end do
        write(o,*)''
        do p1=par1,par2
           !if(stdev1(p1).lt.1.d-20) cycle
           if(fixedpar(p1).eq.1) cycle
           write(o,'(A8,$)')trim(varnames(p1))
           do p2=par1,par2
              corr = corrs(p1,p2)
              !if(stdev1(p2).lt.1.d-20) cycle
              if(fixedpar(p2).eq.1) cycle
              if( (abs(corr).ge.corr2.and.p2.gt.p1) ) then  !Print in the upper triangle
                 write(o,'(F7.2,$)')corr
              else if( (abs(corr).ge.corr1.and.abs(corr).lt.corr2.and.p1.gt.p2) ) then   !Print in the lower triangle
                 write(o,'(F7.2,$)')corr
              else if(p1.eq.p2) then !Print on the diagonal
                 !write(o,'(F7.2,$)')corr !=1
                 !write(o,'(A7,$)')'\-----\'
                 !write(o,'(A7,$)')'\______' !'
                 write(o,'(A7,$)')' ######' !'
              else
                 write(o,'(A7,$)')''
              end if
           end do
           write(o,'(A)')'   '//trim(varnames(p1))
        end do
     end if
     
     
     
     
     
     !Print output for CBC Wiki:
     if(ic.eq.1) then
        o = 10
        open(unit=o,form='formatted',status='replace',action='write',position='rewind',file='wiki.txt',iostat=io)
        if(io.ne.0) then
           write(*,'(A)')'  Error opening wiki.txt, aborting...'
           stop
        end if
        write(gps,'(I9)')nint(t0)+floor(startval(ic,4,1))
        write(url,'(A)')'http://www.astro.northwestern.edu/~sluys/CBC/gps'//trim(gps)//'/'
        write(o,'(A)')'= GPS'//trim(gps)//' - description ='
        write(o,'(/,A)')'Back to [:JointS5/InterestingEventsForBayesianFollowUp:Interesting events for Bayesian follow-up]'
        
        
        
        !True values:
        write(o,'(///,A)')'== True values =='
        write(o,'(A)')"|| '''Detectors'''  || '''M1'''    || '''M2'''    || '''Mc'''    || '''&eta;''' || '''time'''      || '''spin'''  ||'''&theta;'''|| '''Dist'''  || '''R.A.'''  || '''Dec.'''  || '''incl'''  || '''pol.'''  || '''details'''                                                                                     ||"
        write(o,'(A)')"||                  ||  (Mo)       || (Mo)        || (Mo)        ||             ||  (s)            ||             || (rad)       || (Mpc)       || (rad)       || (rad)       || (rad)       || (rad)       ||                                                                                                   ||"
        
        write(o,'(A4,$)')'|| !'
        do i=1,4
           if(i.le.ndet(ic)) then
              write(o,'(A2,$)')detabbrs(detnr(ic,i))
           else
              write(o,'(A2,$)')'  '
           end if
        end do
        write(o,'(A6,$)')'    '
        parr(1:12) = (/14,15,2,3,4,6,7,5,8,9,11,12/)
        do p=1,12
           p1 = parr(p)
           x = stats(ic,p1,1)
           if(p1.eq.8) x = rev2pi(x*rh2r)
           if(p1.eq.12.or.p1.eq.13) x = rev2pi(x*rd2r)
           if(p1.eq.7.or.p1.eq.9.or.p1.eq.11) x = x*rd2r
           if(p1.eq.4) then
              write(o,'(A5,F14.4,$)')'  || ',x+t0
           else
              write(o,'(A5,F10.4,$)')'  || ',x
           end if
        end do
        write(o,'(A)')'  || [ injection info]  ||'
        
        
        
        
        !Bayes factor:
        write(o,'(///,A)')'== Bayes Factors =='
        write(o,'(A)')"|| '''Model'''                                      || '''Detectors'''        || '''log_e Bayes Factor'''    || '''log_10 Bayes Factor'''    || '''Details'''                                                           ||"

        if(fixedpar(6).eq.0) write(o,'(A,$)')'|| 1.5pN simple precession vs. Gaussian noise       || !'
        if(fixedpar(6).eq.1) write(o,'(A,$)')'|| 1.5pN non-spinning vs. Gaussian noise            || !'
        do i=1,4
           if(i.le.ndet(ic)) then
              write(o,'(A2,$)')detabbrs(detnr(ic,i))
           else
              write(o,'(A2,$)')'  '
           end if
        end do
        write(o,'(A17,$)')'              || '
        write(o,'(F10.1,A,$)')logebayesfactor(ic),'                  || '
        write(o,'(F10.1,A,$)')log10bayesfactor(ic),'                   || '
        write(o,'(A)')'['//trim(url)//' link]       ||'
        
        
        
        
        !Medians:
        write(o,'(///,A)')'== Medians =='
        write(o,'(A)')"|| '''Code'''                    || '''Detectors'''  || '''Mc'''    || '''&eta;''' || '''time'''  || '''spin'''  ||'''&theta;'''|| '''Dist'''  || '''R.A.'''  || '''Dec.'''  || '''incl'''  || '''pol.'''  || '''details'''                                                                                     ||"
        write(o,'(A)')"||                               ||                  || (Mo)        ||             ||  (s)        ||             || (rad)       || (Mpc)       || (rad)       || (rad)       || (rad)       || (rad)       ||                                                                                                   ||"
        if(fixedpar(6).eq.0) write(o,'(A,$)')'|| 1.5pN, spinning MCMC          || !'
        if(fixedpar(6).eq.1) write(o,'(A,$)')'|| 1.5pN, non-spinning MCMC      || !'
        do i=1,4
           if(i.le.ndet(ic)) then
              write(o,'(A2,$)')detabbrs(detnr(ic,i))
           else
              write(o,'(A2,$)')'  '
           end if
        end do
        write(o,'(A10,$)')'       ||'
        parr(1:10) = (/2,3,4,6,7,5,8,9,11,12/)
        do p=1,10
           p1 = parr(p)
           x = stats(ic,p1,1)
           if(p1.eq.8) x = rev2pi(x*rh2r)
           if(p1.eq.12.or.p1.eq.13) x = rev2pi(x*rd2r)
           if(p1.eq.7.or.p1.eq.9.or.p1.eq.11) x = x*rd2r
           write(xs10,'(F10.4)')x
           if(fixedpar(6).eq.1 .and. (p1.eq.6.or.p1.eq.7.or.p1.eq.13)) xs10 = ' - '  !If non-spinning
           write(o,'(A10,A5,$)')xs10,'   ||'
        end do
        write(o,'(A)')' ['//trim(url)//' link]                                 ||'
        
        
        
        
        !Means:
        write(o,'(///,A)')'== Means =='
        write(o,'(A)')"|| '''Code'''                    || '''Detectors'''  || '''Mc'''    || '''&eta;''' || '''time'''  || '''spin'''  ||'''&theta;'''|| '''Dist'''  || '''R.A.'''  || '''Dec.'''  || '''incl'''  || '''pol.'''  || '''details'''                                                                                     ||"
        write(o,'(A)')"||                               ||                  || (Mo)        ||             ||  (s)        ||             || (rad)       || (Mpc)       || (rad)       || (rad)       || (rad)       || (rad)       ||                                                                                                   ||"
        if(fixedpar(6).eq.0) write(o,'(A,$)')'|| 1.5pN, spinning MCMC          || !'
        if(fixedpar(6).eq.1) write(o,'(A,$)')'|| 1.5pN, non-spinning MCMC      || !'
        do i=1,4
           if(i.le.ndet(ic)) then
              write(o,'(A2,$)')detabbrs(detnr(ic,i))
           else
              write(o,'(A2,$)')'  '
           end if
        end do
        write(o,'(A10,$)')'       ||'
        parr(1:10) = (/2,3,4,6,7,5,8,9,11,12/)
        do p=1,10
           p1 = parr(p)
           x = stats(ic,p1,2)
           if(p1.eq.8) x = rev2pi(x*rh2r)
           if(p1.eq.12.or.p1.eq.13) x = rev2pi(x*rd2r)
           if(p1.eq.7.or.p1.eq.9.or.p1.eq.11) x = x*rd2r
           write(xs10,'(F10.4)')x
           if(fixedpar(6).eq.1 .and. (p1.eq.6.or.p1.eq.7.or.p1.eq.13)) xs10 = ' - '  !If non-spinning
           write(o,'(A10,A5,$)')xs10,'   ||'
        end do
        write(o,'(A)')' ['//trim(url)//' link]                                 ||'
        
        
        
        
        !Lmax:
        write(o,'(///,A)')'== Maximum-likelihood points =='
        write(o,'(A)')"|| '''Code'''                    || '''Detectors'''  ||'''log(L)''' || '''Mc'''    || '''&eta;''' || '''time'''  || '''spin'''  ||'''&theta;'''|| '''Dist'''  || '''R.A.'''  || '''Dec.'''  || '''incl'''  || '''pol.'''  || '''details'''                                                                                     ||"
        write(o,'(A)')"||                               ||                  ||             || (Mo)        ||             ||  (s)        ||             || (rad)       || (Mpc)       || (rad)       || (rad)       || (rad)       || (rad)       ||                                                                                                   ||"
        if(fixedpar(6).eq.0) write(o,'(A,$)')'|| 1.5pN, spinning MCMC          || !'
        if(fixedpar(6).eq.1) write(o,'(A,$)')'|| 1.5pN, non-spinning MCMC      || !'
        do i=1,4
           if(i.le.ndet(ic)) then
              write(o,'(A2,$)')detabbrs(detnr(ic,i))
           else
              write(o,'(A2,$)')'  '
           end if
        end do
        write(o,'(A10,$)')'       ||'
        parr(1:11) = (/1,2,3,4,6,7,5,8,9,11,12/) !Include logL!
        do p=1,11
           p1 = parr(p)
           x = startval(ic,p1,3)
           if(p1.eq.8) x = rev2pi(x*rh2r)
           if(p1.eq.10.or.p1.eq.12.or.p1.eq.13) x = rev2pi(x*rd2r)
           if(p1.eq.7.or.p1.eq.9.or.p1.eq.11) x = x*rd2r
           write(xs10,'(F10.4)')x
           if(fixedpar(6).eq.1 .and. (p1.eq.6.or.p1.eq.7.or.p1.eq.13)) xs10 = ' - '  !If non-spinning
           write(o,'(A10,A5,$)')xs10,'   ||'
        end do
        write(o,'(A)')' ['//trim(url)//' link]                                 ||'
        
        
        
        
        !Stdev:
        write(o,'(///,A)')'== Standard deviations =='
        write(o,'(A)')"|| '''Code'''                    || '''Detectors'''  || '''Mc'''    || '''&eta;''' || '''time'''  || '''spin'''  ||'''&theta;'''|| '''Dist'''  || '''R.A.'''  || '''Dec.'''  || '''incl'''  || '''pol.'''  || '''details'''                                                                                     ||"
        write(o,'(A)')"||                               ||                  || (Mo)        ||             ||  (s)        ||             || (rad)       || (Mpc)       || (rad)       || (rad)       || (rad)       || (rad)       ||                                                                                                   ||"
        if(fixedpar(6).eq.0) write(o,'(A,$)')'|| 1.5pN, spinning MCMC          || !'
        if(fixedpar(6).eq.1) write(o,'(A,$)')'|| 1.5pN, non-spinning MCMC      || !'
        do i=1,4
           if(i.le.ndet(ic)) then
              write(o,'(A2,$)')detabbrs(detnr(ic,i))
           else
              write(o,'(A2,$)')'  '
           end if
        end do
        write(o,'(A10,$)')'       ||'
        parr(1:10) = (/2,3,4,6,7,5,8,9,11,12/)
        do p=1,10
           p1 = parr(p)
           x = stdev2(p1)
           if(p1.eq.8) x = x*rh2r
           if(p1.eq.10.or.p1.eq.12.or.p1.eq.13) x = x*rd2r
           if(p1.eq.7.or.p1.eq.9.or.p1.eq.11) x = x*rd2r
           write(xs10,'(F10.4)')x
           if(fixedpar(6).eq.1 .and. (p1.eq.6.or.p1.eq.7.or.p1.eq.13)) xs10 = ' - '  !If non-spinning
           write(o,'(A10,A5,$)')xs10,'   ||'
        end do
        write(o,'(A)')' ['//trim(url)//' link]                                 ||'
        
        
        
        !2-sigma range:
        write(o,'(///,A)')'== 2-sigma probability ranges =='
        write(o,'(A)')"|| '''Code'''                    || '''Detectors'''  || '''Mc'''              || '''&eta;'''           || '''time'''            || '''spin'''            || '''&theta;'''         || '''Distance'''        || '''R.A.'''            || '''Dec.'''            || '''incl.'''           || '''pol.'''            || '''details'''                                                                                     ||"
        write(o,'(A)')"||                               ||                  || (Mo)                  ||                       ||  (s)                  ||                       || (rad)                 || (Mpc)                 || (rad)                 || (rad)                 || (rad)                 || (rad)                 ||                                                                                                   ||"
        c = 0
        do i=1,nival
           if(abs(ivals(i)-0.9545).lt.0.0001) c = i
        end do
        if(c.eq.0) then
           write(0,'(A)')'  Error: 2-sigma range not found, needed for Wiki output!'
           stop
        end if
        if(fixedpar(6).eq.0) write(o,'(A,$)')'|| 1.5pN, spinning MCMC          || !'
        if(fixedpar(6).eq.1) write(o,'(A,$)')'|| 1.5pN, non-spinning MCMC      || !'
        do i=1,4
           if(i.le.ndet(ic)) then
              write(o,'(A2,$)')detabbrs(detnr(ic,i))
           else
              write(o,'(A2,$)')'  '
           end if
        end do
        write(o,'(A10,$)')'       ||'
        parr(1:10) = (/2,3,4,6,7,5,8,9,11,12/)
        do p=1,10
           p1 = parr(p)
           x1 = ranges(ic,c,p1,1)
           x2 = ranges(ic,c,p1,2)
           if(p1.eq.8) then
              x1 = x1*rh2r
              x2 = x2*rh2r
           end if
           if(p1.eq.12.or.p1.eq.13) then
              x1 = x1*rd2r
              x2 = x2*rd2r
           end if
           if(p1.eq.7.or.p1.eq.9.or.p1.eq.11) then
              x1 = x1*rd2r
              x2 = x2*rd2r
           end if
           write(xs20,'(F9.4,A2,F9.4)')x1,' -',x2
           if(fixedpar(6).eq.1 .and. (p1.eq.6.or.p1.eq.7.or.p1.eq.13)) xs20 = ' - '  !If non-spinning
           write(o,'(A20,A5,$)')xs20,'   ||'
           
        end do
        write(o,'(A)')' ['//trim(url)//' link]                                 ||'


        
        !True values:
        write(o,'(///,A)')'== True values =='
        write(o,'(A)')"|| '''Detectors'''  || '''M1'''    || '''M2'''    || '''Mc'''    || '''&eta;''' || '''time'''      || '''spin'''  ||'''&theta;'''|| '''Dist'''  || '''R.A.'''  || '''Dec.'''  || '''incl'''  || '''pol.'''  || '''details'''                                                                                     ||"
        write(o,'(A)')"||                  ||  (Mo)       || (Mo)        || (Mo)        ||             ||  (s)            ||             || (rad)       || (Mpc)       || (rad)       || (rad)       || (rad)       || (rad)       ||                                                                                                   ||"
        
        write(o,'(A4,$)')'|| !'
        do i=1,4
           if(i.le.ndet(ic)) then
              write(o,'(A2,$)')detabbrs(detnr(ic,i))
           else
              write(o,'(A2,$)')'  '
           end if
        end do
        write(o,'(A6,$)')'    '
        parr(1:12) = (/14,15,2,3,4,6,7,5,8,9,11,12/)
        do p=1,12
           p1 = parr(p)
           x = stats(ic,p1,1)
           if(p1.eq.8) x = rev2pi(x*rh2r)
           if(p1.eq.12.or.p1.eq.13) x = rev2pi(x*rd2r)
           if(p1.eq.7.or.p1.eq.9.or.p1.eq.11) x = x*rd2r
           if(p1.eq.4) then
              write(o,'(A5,F14.4,$)')'  || ',x+t0
           else
              write(o,'(A5,F10.4,$)')'  || ',x
           end if
        end do
        write(o,'(A)')'  || [ injection info]  ||'
        
        
        write(o,'(///,A)')'----'
        write(o,'(A)')'Back to [:JointS5/InterestingEventsForBayesianFollowUp:Interesting events for Bayesian follow-up]'

        
     end if
     close(o)
     







     
     
     
     
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
        write(*,'(A18,$)')''
        do p=par1,min(par2,npar0)
           if(fixedpar(p).eq.1) cycle
           write(*,'(A11,$)')trim(varnames(p))
        end do
        write(*,'(A11)')'total'
        
        if(prconv.ge.3) then
           write(*,'(A)')'  Means:'
           do ic=1,nchains0
              write(*,'(I16,A2,$)')ic,': '
              do p=par1,min(par2,npar0)
                 if(fixedpar(p).eq.1) cycle
                 write(*,'(F11.6,$)')chmean(ic,p)
              end do
              write(*,*)
           end do
           write(*,'(A18,$)')'           Total: '
           do p=par1,min(par2,npar0)
              if(fixedpar(p).eq.1) cycle
              write(*,'(F11.6,$)')totmean(p)
           end do
           write(*,*)''
           
           write(*,*)''
           write(*,'(A)')'  Variances:'
        end if !if(prconv.ge.3)
     end if !if(prconv.ge.1)
     
     do ic=1,nchains0
        !write(*,'(I16,A2,20F11.6)')ic,': ',chvar1(ic,par1:min(par2,13))
        !if(chvar1(ic,2).lt.0.5*chvar(2).and.chvar1(ic,3).lt.0.5*chvar(3).and.chvar1(ic,2).lt.0.5*chvar(2).and.chvar1(ic,2).lt.0.5*chvar(2)) then
        lowvar = 0
        highvar = 0
        totrelvar = 1.d0
        ntotrelvar = 0
        do p=par1,min(par2,13)
           !if(abs(chvar1(ic,p)).gt.1.e-30) then  !Take only the parameters that were fitted and have a variance > 0
           if(fixedpar(p).eq.0 .and.abs(chvar1(ic,p)).gt.1.e-30) then  !Take only the parameters that were fitted and have a variance > 0
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
        if(prconv.ge.3) then
           ch = ' '
           if(nlowvar.eq.4) ch = '*'
           if(nhighvar.eq.4) ch = '#'
           write(*,'(I7,A3,$)')ic,': '//ch
           ch = ' '
           if(totrelvar.lt.0.5) ch = '*'
           if(totrelvar.gt.2.0) ch = '#'
           write(*,'(F8.3,A1,$)')totrelvar,ch
           do p=par1,min(par2,13)
              if(fixedpar(p).eq.1) cycle
              ch = ' '
              if(lowvar(p).eq.1) ch = '*'
              if(highvar(p).eq.1) ch = '#'
              write(*,'(F10.5,A1,$)')chvar1(ic,p),ch
           end do
           write(*,*)''
        end if !if(prconv.ge.3)
     end do
     if(prconv.ge.3) then
        write(*,'(A18,$)')'  Total:          '
        do p=par1,min(par2,npar0)
           if(fixedpar(p).eq.1) cycle
           write(*,'(F11.5,$)')chvar(p)
        end do
        write(*,*)''
        write(*,*)''
     end if !if(prconv.ge.3)
     if(prconv.ge.2) then
        write(*,'(A)')'  Variances:'
        write(*,'(A18,$)')'   Within chains: '
        do p=par1,min(par2,npar0)
           if(fixedpar(p).eq.1) cycle
           write(*,'(ES11.3,$)')chvar(p)
        end do
        write(*,*)
        write(*,'(A18,$)')'  Between chains: '
        do p=par1,min(par2,npar0)
           if(fixedpar(p).eq.1) cycle
           write(*,'(ES11.3,$)')totvar(p)
        end do
        write(*,*)
     end if
     
     if(prconv.ge.1) then
        write(*,'(A18,$)')'     Convergence: '
        totrhat = 1.d0
        nrhat = 0
        do p=par1,min(par2,npar0)
           if(fixedpar(p).eq.1) cycle
           write(*,'(F11.5,$)')rhat(p)
           if(p.gt.1) then !Don't include logL
              totrhat = totrhat * rhat(p)
              nrhat = nrhat + 1
           end if
        end do
        !write(*,'(F11.5)')sum(rhat(par1:min(par2,13)))/dble(min(par2,13)-par1+1)
        !write(*,'(F11.5)')totrhat/dble(nrhat)
        write(*,'(F11.5)')totrhat**(1.d0/dble(nrhat))
     end if
     !write(*,*)''
  end if
  
  
  
  
!  !Test: get mean and stdev for log(L)
!  if(1.eq.2) then
!     write(*,*)''
!     
!     nn = minval(ntot(1:nchains0))/2
!     nlogl1 = nn+1
!     nlogl2 = 2*nn
!     nn = abs(nlogl2-nlogl1)
!     
!     write(*,'(A,I7,A)')'  Convergence criterion for',nn,' parameters:'
!     write(*,'(16x,16x,20A11)')'Mean','Stddev','M-S','M+S','KS d','KS prob'
!     
!     do ic=1,nchains0
!        nn1 = ntot(ic)/20
!        ksn2 = 0
!        ksd = 1.
!        ksprob = 0.
!        do nlogl1 = 1,ntot(ic),nn1
!           nlogl2 = min(nlogl1+nn1,ntot(ic))
!           nn = abs(nlogl2-nlogl1)+1
!           if(nn.lt.nn1) cycle
!           
!           chmean = 1.d-30
!           totmean = 1.d-30
!           chvar = 1.d-30
!           chvar1 = 1.d-30
!           totvar = 1.d-30
!           
!           p=1
!           do i=nlogl1,nlogl2
!              chmean(ic,p) = chmean(ic,p) + dat(p,ic,i) !Can't use pldat, because it may get wrapped above
!              totmean(p) = totmean(p) + dat(p,ic,i)
!           end do
!           chmean = chmean/dble(nn)
!           totmean = totmean/dble(nn*nchains0)
!           
!           do i=nlogl1,nlogl2
!              dx = (dat(p,ic,i) - chmean(ic,p))**2 !Can't use pldat, because it may get wrapped above
!              chvar(p) = chvar(p) + dx
!              chvar1(ic,p) = chvar1(ic,p) + dx !Keep track of the variance per chain
!           end do
!           totvar(p) = totvar(p) + (chmean(ic,p) - totmean(p))**2
!           chvar1(ic,p) = chvar1(ic,p)/dble(nn-1)
!           chvar(p) = chvar(p)/dble(nchains0*(nn-1))
!           totvar(p) = totvar(p)/dble(nchains0-1)
!           
!           !write(*,'(I16,2I8,20F11.6)')ic,nlogl1,nlogl2,chmean(ic,1),chvar1(ic,1),chmean(ic,1)-chvar1(ic,1),chmean(ic,1)+chvar1(ic,1),ksd,ksprob
!           !!print*,ic,p,nlogl1,nlogl2,nn
!           ksdat1(1:nn) = dble(dat(p,ic,nlogl1:nlogl2))
!           ksn1   = nn
!           !!call kstwo(data1,n1,data2,n2,d,prob)
!           if(ksn2.ne.0) call kstwo(ksdat1(1:ksn1),ksn1,ksdat2(1:ksn2),ksn2,ksd,ksprob)
!           !
!           !ksdat2 = ksdat1
!           ksdat2(1:nn) = dble(dat(p,ic,nlogl1:nlogl2))
!           ksn2 = ksn1
!           
!           write(*,'(I16,2I8,20F11.6)')ic,nlogl1,nlogl2,chmean(ic,1),chvar1(ic,1),chmean(ic,1)-chvar1(ic,1),chmean(ic,1)+chvar1(ic,1),ksd,ksprob,dlog10(ksprob+1.d-100)
!        end do
!        write(*,*)''
!     end do
!     write(*,*)''
!     
!     !KS test
!     call pgbegin(1,'21/xs',1,1)
!     call pgpap(scrsz,scrrat)
!     call pgsch(1.5)
!     call pgsvp(0.07,0.99,0.10,0.96)
!     call pgswin(0.,real(maxval(ntot(1:nchains0))),-100.,0.)
!     call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)
!     call pgmtxt('B',2.2,0.5,0.5,'i')
!     call pgmtxt('L',1.8,0.5,0.5,'log(d\dKS\u)')
!     
!     do ic=1,nchains0
!        call pgsci(colours(mod(ic-1,ncolours)+1))
!        nn = ntot(ic)/10
!        ksd = 1.
!        ksprob = 0.
!        do nlogl1 = 1,ntot(ic),nn
!           nlogl2 = nlogl1+nn-1
!           if(nlogl2.gt.ntot(ic)) cycle
!           
!           ksn1   = nn
!           ksdat1(1:ksn1) = dble(dat(p,ic,nlogl1:nlogl2))
!           ksn2   = ntot(ic)-nlogl1+1
!           ksdat2(1:ksn2) = dble(dat(p,ic,nlogl1:ntot(ic)))
!           
!           if(ksn2.ne.0) call kstwo(ksdat1(1:ksn1),ksn1,ksdat2(1:ksn2),ksn2,ksd,ksprob)
!           !write(*,'(I16,2I8,20F11.6)')ic,nlogl1,nlogl2,ksd,ksprob,dlog10(ksprob+1.d-100)
!           call pgpoint(1,real(nlogl1+nlogl2)/2.,real(dlog10(ksprob+1.d-100)),2)
!        end do
!        write(*,*)''
!     end do
!     
!     call pgend
!     write(*,*)''
!     
!  end if
  
  
end subroutine statistics
!***********************************************************************************************************************************







!***********************************************************************************************************************************
subroutine printstats(exitcode)
  use constants
  use analysemcmc_settings
  use general_data
  use mcmcrun_data
  use stats_data
  use chain_data
  implicit none
  integer :: c,i,ic,o,p,p1,p2,exitcode,system
  
  exitcode = 0
  ic = 1 !Use chain 1
  o = 20 !Output port
  open(unit=o, form='formatted', status='replace',file=trim(outputdir)//'/'//trim(outputname)//'__statistics.dat')
  write(o,'(A)')trim(outputname)
  
  !Print general run and detector info:
  write(o,'(//,A,/)')'GENERAL INFORMATION:'
  write(o,'(6x,4A12,A12,A5  A8,A22,A8)')'totiter','totlines','totpts','totburn','nchains','used','seed','null likelihood','ndet'
  write(o,'(6x,4I12,I12,I5, I8,F22.10,I8)')totiter,totlines,totpts,totlines-totpts,nchains0,contrchains,seed(ic),nullh,ndet(ic)
  write(o,*)''
  write(o,'(A14,A3,A18,4A12,A22,A17,3A14)')'Detector','Nr','SNR','f_low','f_high','before tc','after tc','Sample start (GPS)','Sample length','Sample rate','Sample size','FT size'
  do i=1,ndet(ic)
     write(o,'(A14,I3,F18.8,4F12.2,F22.8,F17.7,3I14)')detnames(ic,i),detnr(ic,i),snr(ic,i),flow(ic,i),fhigh(ic,i),t_before(ic,i),t_after(ic,i),FTstart(ic,i),deltaFT(ic,i),samplerate(ic,i),samplesize(ic,i),FTsize(ic,i)
  end do
  write(o,*)''
  
  write(o,'(A,I11)')' t0:',nint(t0)
  
  !Print statistics
  write(o,'(///,A,/)')'BASIC STATISTICS:'
  write(o,'(A,2I3)')'Npar,ncol:',par2-par1+1,7
  write(o,'(A8,7A12)')'param.','model','median','mean','stdev1','stdev2','abvar1','abvar2'
  
  do p=par1,par2
     if(fixedpar(p).eq.0) then
        write(o,'(A8,7F12.6)')varnames(p),startval(ic,p,1),stats(ic,p,1),stats(ic,p,2),stdev1(p),stdev2(p),absvar1(p),absvar2(p)
     else
        write(o,'(A8,7F12.6)')varnames(p),startval(ic,p,1),stats(ic,p,1),stats(ic,p,2),0.,0.,0.,0.
     end if
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
  do c=1,nival
     write(o,'(F21.5,A14,$)')ivals(c),''
  end do
  write(o,*)''
  
  write(o,'(A8,2x,$)')'param.'
  do c=1,nival
     !write(o,'(2x,2A9,A8,$)')'rng1','rng2','in rnge'
     write(o,'(2x,2A12,A9,$)')'centre','delta','in rnge'
  end do
  write(o,*)''
  do p=par1,par2
     write(o,'(A8,2x,$)')trim(varnames(p))
     do c=1,nival
        !write(o,'(2x,2F11.6,F6.3,$)')ranges(ic,c,p,1),ranges(ic,c,p,2),2*abs(startval(ic,p,1)-ranges(ic,c,p,3))/ranges(ic,c,p,4) !Defined with centre of prob. range
        if(fixedpar(p).eq.0) then
           write(o,'(2x,2F12.6,F7.3,$)')ranges(ic,c,p,3),ranges(ic,c,p,4),2*abs(startval(ic,p,1)-ranges(ic,c,p,3))/ranges(ic,c,p,4) !Defined with centre of prob. range
        else
           write(o,'(2x,2F12.6,F7.3,$)')0.,0.,99.999
        end if
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
  do c=1,nival
     write(o,'(F19.5,$)')ivals(c)
  end do
  write(o,*)''
  
  write(o,'(A9,A19,$)')'params.',''
  do c=1,nival
     write(o,'(A16,A3,$)')'delta','in'
  end do
  write(o,*)''
  do p=1,npdf2d
     p1 = pdf2dpairs(p,1)
     p2 = pdf2dpairs(p,2)
     write(o,'(2I4,2(2x,A8),2x,$)')p1,p2,trim(varnames(p1)),trim(varnames(p2))
     do c=1,nival
        write(o,'(2x,F14.8,$)')probareas(p1,p2,c,1)
        if(c.ge.nival+1-trueranges2d(p1,p2) .and. trueranges2d(p1,p2).ne.0) then
           write(o,'(A3,$)')'y'
        else
           write(o,'(A3,$)')'n'
        end if
     end do
     write(o,*)''
  end do
  
  
  
  close(o) !Statistics output file
  if(savestats.eq.2) i = system('a2ps -1rf7 '//trim(outputdir)//'/'//trim(outputname)//'__statistics.dat -o '//trim(outputdir)//'/'//trim(outputname)//'__statistics.ps')
  !write(*,*)''
  if(prprogress.ge.1) then
     if(savestats.eq.1) write(*,'(A)')'  Statistics saved in '//trim(outputname)//'__statistics.dat'
     if(savestats.eq.2) write(*,'(A)')'  Statistics saved in '//trim(outputname)//'__statistics.dat,ps'
  end if
  
end subroutine printstats
!***********************************************************************************************************************************
