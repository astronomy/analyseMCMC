! Compute statistics for analysemcmc

subroutine statistics(exitcode)
  use constants
  use analysemcmc_settings
  use general_data
  use stats_data
  use chain_data
  use mcmcrun_data
  implicit none
  integer :: c,i,ic,i0,j,j1,o,p,p1,p2,nr,nstat,exitcode,wraptype
  integer :: indexx(maxMCMCpar,maxChs*maxIter),index1(maxChs*maxIter)
  real :: rev2pi,x0,x1,x2,y1,y2,rrevpi
  real :: range1,minrange,maxgap,ival,wrapival,centre,maxlogl,minlogl,shival,shival2
  real :: medians(maxMCMCpar),mean(maxMCMCpar),var1(maxMCMCpar),var2(maxMCMCpar),corr,corr1,corr2
  
  !Need extra accuracy to compute Bayes factor
  !real*8 :: var,total
  !real*16 :: var,total  !gfortran doesn't support this
  !real(kind=10) :: var,total  !Means different things on different compilers
  real(kind=selected_real_kind(18,4931)) :: var,total  !Precision of 18, exponential range of 4931 - max gfortran supports?
  !print*,precision(var),range(var)
  
  exitcode = 0
  
  
  !Sort all data and find the interval limits for the default probability interval for the wrapable parameters
  if(prProgress.ge.2) write(6,*)''
  shift = 0.
  wrap = 0
  rashift = 0.
  do ic=1,nchains
     if(mergeChains.eq.0.and.contrchain(ic).eq.0) cycle
     !wrapival = ivals(Nival) !Use the largest range
     wrapival = 0.999 !Always use a very large range (?)
     indexx = 0
     if(prProgress.ge.2.and.mergeChains.eq.0) write(6,'(A,I2.2,A,$)')' Ch',ic,' '
     if(prProgress.ge.2.and.ic.eq.1.and.wrapData.ge.1) write(6,'(A,$)')'  Wrap data. '
     do p=1,nMCMCpar
        if(wrapData.eq.0 .or. &
             (parID(p).ne.31.and.parID(p).ne.41.and.parID(p).ne.52.and.parID(p).ne.54.and.parID(p).ne.73.and.parID(p).ne.83) ) then  !Not RA, phi_c, psi, phi_Jo, phi_1,2
           call rindexx(n(ic),selDat(ic,p,1:n(ic)),index1(1:n(ic)))
           indexx(p,1:n(ic)) = index1(1:n(ic))
           if(parID(p).eq.31) racentre = rpi !Plot 0-24h when not wrapping -> centre = 12h = pi
           cycle !No wrapping
        end if
        
        wraptype = 1  !0-2pi (e.g. phases)
        if(parID(p).eq.52) wraptype = 2  !0-pi (polarisation angle)
        
        
        
        !Make sure data are between 0 and 2pi or between 0 and pi to start with:
        do i=1,n(ic)
           if(wraptype.eq.1) selDat(ic,p,i) = rev2pi(selDat(ic,p,i)) 
           if(wraptype.eq.2) selDat(ic,p,i) = rrevpi(selDat(ic,p,i)) !Pol.angle
        end do
        call rindexx(n(ic),selDat(ic,p,1:n(ic)),index1(1:n(ic)))
        indexx(p,1:n(ic)) = index1(1:n(ic))
        
        minrange = 1.e30
        do i=1,n(ic)
           x1 = selDat(ic,p,indexx(p,i))
           x2 = selDat(ic,p,indexx(p,mod(i+nint(n(ic)*wrapival)-1,n(ic))+1))
           if(wraptype.eq.1) range1 = mod(x2 - x1 + real(20*pi),rtpi)
           if(wraptype.eq.2) range1 = mod(x2 - x1 + real(10*pi),rpi)
           if(range1.lt.minrange) then
              minrange = range1
              y1 = x1
              y2 = x2
              !write(6,'(2I6,7F10.5)')i,mod(nint(i+n(ic)*wrapival),n(ic)),x1,x2,range1,minrange,y1,y2,(y1+y2)/2.
           end if
           !write(6,'(2I6,7F10.5)')i,mod(nint(i+n(ic)*wrapival),n(ic)),x1,x2,range1,minrange,y1,y2,(y1+y2)/2.
        end do !i
        centre = (y1+y2)/2.
        
        !Define shift interval:
        shival = rtpi   !Shift interval
        shival2 = rpi   !Shift interval/2
        if(wraptype.eq.2) then
           shival = rpi   !Shift interval
           shival2 = rpi2 !Shift interval/2
        end if
        
        if(y1.gt.y2) then
           wrap(ic,p) = wraptype
           centre = mod(centre + shival2, shival) !Then distribution peaks close to 0/shival, shift centre by shival/2
        end if
        
        
        !See whether there's a gap in the data.  WHY is this necessary, should it work like this???
        if(wrap(ic,p).eq.0 .and. 1.eq.2) then
           i0 = -1
           maxgap = -1.e30
           do i=1,n(ic)-1
              x1 = selDat(ic,p,indexx(p,i))
              x2 = selDat(ic,p,indexx(p,i+1))
              !write(6,'(2I3,2I8,4F10.5)')ic,p,i,i0,x1,x2,x2-x2,maxgap
              if(x2-x1.gt.maxgap) then
                 maxgap = x2-x1
                 i0 = i
              end if
           end do !i
           x1 = selDat(ic,p,indexx(p,i0))
           x2 = selDat(ic,p,indexx(p,i0+1))
           !if(maxgap.gt.2*tpi/sqrt(real(n(ic)))) then
           if(maxgap.gt.0.1) then 
              x0 = (x1+x2)/2.
              !write(6,'(10F10.5)')x1,x2,(x1+x2)/2.,maxgap,ymin,ymax,centre,minrange,y1,y2
              !if(y1.lt.y2.and.(x0.lt.y1.or.x0.gt.y2)) wrap(ic,p) = 1  !If centre of max gap is outside 90% range  WHY???
              !if(y1.gt.y2.and.(x0.gt.y2.and.x0.lt.y1)) wrap(ic,p) = 1
           end if
        end if
        
        
        
        !Now, wrap around anticentre
        shift(ic,p) = 0.
        
        !For the general case of shival (= pi or 2pi)
        if(wrap(ic,p).gt.0) shift(ic,p) = shival - mod(centre + shival2, shival)
        if(parID(p).eq.31.and.ic.eq.1) rashift = shift(ic,p)                         !Save RA shift to plot sky map
        
        !Do the actual wrapping:
        allDat(ic,p,1:ntot(ic))  = mod(allDat(ic,p,1:ntot(ic))  + shift(ic,p), shival) - shift(ic,p)   !Original data
        selDat(ic,p,1:n(ic))     = mod(selDat(ic,p,1:n(ic))     + shift(ic,p), shival) - shift(ic,p)
        startval(ic,p,1:3)       = mod(startval(ic,p,1:3)       + shift(ic,p), shival) - shift(ic,p)   !True, starting and Lmax values
        y1 = mod(y1 + shift(ic,p), shival) - shift(ic,p)
        y2 = mod(y2 + shift(ic,p), shival) - shift(ic,p)
        
        centre = mod(centre + shift(ic,p), shival) - shift(ic,p)
        if(parID(p).eq.31.and.ic.eq.1) racentre = centre                             !Save RA centre to plot sky map
        
        minrange = y2-y1
        call rindexx(n(ic),selDat(ic,p,1:n(ic)),index1(1:n(ic)))  !Re-sort
        indexx(p,1:n(ic)) = index1(1:n(ic))
        
        if(abs( abs( minval(selDat(ic,p,1:n(ic))) - maxval(selDat(ic,p,1:n(ic))) ) - shival) .lt. 1.e-3)  wrap(ic,p) = 1   !If centre is around shival2, still needs to be flagged 'wrap' to plot PDF
     end do !p
     ! End wrapping data
     
     
     
     
     !Do statistics
     !if(prProgress.ge.2) write(6,'(A)')' Calculating: statistics...'
     if(prProgress.ge.1.and.ic.eq.1) write(6,'(A,$)')'  Calc: stats, '
     do p=1,nMCMCpar
        !Determine the median
        if(mod(n(ic),2).eq.0) medians(p) = 0.5*(selDat(ic,p,indexx(p,n(ic)/2)) + selDat(ic,p,indexx(p,n(ic)/2+1)))
        if(mod(n(ic),2).eq.1) medians(p) = selDat(ic,p,indexx(p,(n(ic)+1)/2))
        
        !Mean:
        mean(p) = sum(selDat(ic,p,1:n(ic)))/real(n(ic))
        
        !Variances, etc:
        var1(p)=0.; var2(p)=0.; absVar1(p)=0.; absVar2(p)=0.; stdev1(p)=0.; stdev2(p)=0.
        do i=1,n(ic)
           var1(p) = var1(p) + (selDat(ic,p,i) - medians(p))**2       !Based on median
           var2(p) = var2(p) + (selDat(ic,p,i) - mean(p))**2          !Based on mean
           absVar1(p) = absVar1(p) + abs(selDat(ic,p,i) - medians(p)) !Based on median
           absVar2(p) = absVar2(p) + abs(selDat(ic,p,i) - mean(p))    !Based on mean
           stdev1(p) = stdev1(p) + (selDat(ic,p,i) - medians(p))**2   !Based on median
           stdev2(p) = stdev2(p) + (selDat(ic,p,i) - mean(p))**2      !Based on mean
        end do
        
        absVar1(p) = absVar1(p)/real(n(ic))
        absVar2(p) = absVar2(p)/real(n(ic))
        stdev1(p)  = sqrt(stdev1(p)/real(n(ic)-1))
        stdev2(p)  = sqrt(stdev2(p)/real(n(ic)-1))
        
        !Save statistics:
        nstat = 6
        stats(ic,p,1) = medians(p)
        stats(ic,p,2) = mean(p)
        stats(ic,p,3) = absVar1(p)  !Based on median
        stats(ic,p,4) = absVar2(p)  !Based on mean
        stats(ic,p,5) = stdev1(p)   !Based on median
        stats(ic,p,6) = stdev2(p)   !Based on mean
     end do
     
     
     
     !Correlations:
     if(prCorr.gt.0.or.saveStats.gt.0) then
        !write(6,'(A)')' Calculating correlations...   '
        if(prProgress.ge.1) write(6,'(A,$)')' corrs, '
        do p1=1,nMCMCpar
           do p2=1,nMCMCpar
           !do p2=p1,nMCMCpar
              corrs(p1,p2) = 0.
              if(fixedpar(p1)+fixedpar(p2).eq.0) then
                 do i=1,n(ic)
                    !corrs(p1,p2) = corrs(p1,p2) + (selDat(ic,p1,i) - medians(p1))*(selDat(ic,p2,i) - medians(p2))  !Use median
                    corrs(p1,p2) = corrs(p1,p2) + (selDat(ic,p1,i) - mean(p1))*(selDat(ic,p2,i) - mean(p2)) !Use mean; hardly differs from median method
                 end do
                 !corrs(p1,p2) = corrs(p1,p2) / (stdev1(p1)*stdev1(p2)*(n(ic)-1))  !Use median
                 corrs(p1,p2) = corrs(p1,p2) / (stdev2(p1)*stdev2(p2)*(n(ic)-1))  !Use mean
              end if
           end do !p2
        end do !p1
     end if
     
     
     !Autocorrelations:
     if(plACorr.gt.0) then
        !write(6,'(A)')' Calculating autocorrelations...'
        if(prProgress.ge.1) write(6,'(A,$)')' autocorrs, '
        j1 = plACorr/100 !Step size to get 100 autocorrelations per var
        do p=1,nMCMCpar
           acorrs(ic,p,:) = 0.
           !do j=1,ntot(ic)-1
           !do j=1,min(plACorr,ntot(ic)-1)
           do j=0,min(100,ntot(ic)-1)
              do i=1,ntot(ic)-j*j1
                 acorrs(ic,p,j) = acorrs(ic,p,j) + (allDat(ic,p,i) - medians(p))*(allDat(ic,p,i+j*j1) - medians(p))
                 !acorrs(p,j) = acorrs(ic,p,j) + (allDat(ic,p,i) - mean(p))*(allDat(ic,p,i+j*j1) - mean(p))
              end do
              acorrs(ic,0,j) = real(j*j1)
              acorrs(ic,p,j) = acorrs(ic,p,j) / (stdev1(p)*stdev1(p)*(ntot(ic)-j*j1))
              !acorrs(ic,p,j) = acorrs(ic,p,j) / (stdev2(p)*stdev2(p)*(ntot(ic)-j*j1))
           end do !j
           !write(6,*)''
        end do !p
     end if
     
     
     !Determine interval ranges
     !if(prProgress.ge.2) write(6,'(A29,$)')' Determining interval levels: '
     if(prProgress.ge.1.and.ic.eq.1) write(6,'(A,$)')' prob.ivals: '
     c0 = 0
     do c=1,Nival
        ival = ivals(c)
        c0 = ival0
        if(c.ne.c0 .and. prIval.eq.0 .and. prStat.lt.2 .and. saveStats.eq.0) cycle
        
        if(prProgress.ge.1.and.ic.eq.1) write(6,'(F6.3,$)')ival
        do p=1,nMCMCpar
           minrange = 1.e30
           !write(6,'(A8,4x,4F10.5,I4)')parNames(parID(p)),y1,y2,minrange,centre,wrap(ic,p)
           do i=1,floor(n(ic)*(1.-ival))
              x1 = selDat(ic,p,indexx(p,i))
              x2 = selDat(ic,p,indexx(p,i+floor(n(ic)*ival)))
              range1 = abs(x2 - x1)
              !range1 = x2 - x1
              if(range1.lt.minrange) then
                 minrange = range1
                 y1 = x1
                 y2 = x2
              end if
              !write(6,'(I6,7F10.5)')i,x1,x2,range1,minrange,y1,y2,(y1+y2)/2.
           end do
           centre = (y1+y2)/2.
           !write(6,'(A8,4x,4F10.5,I4)')parNames(parID(p)),y1,y2,minrange,centre,wrap(ic,p)
           
           !Save ranges:
           nr = 4                  !Only ranges(:,:,:,1:nr) get converted later on
           ranges(ic,c,p,1) = y1
           ranges(ic,c,p,2) = y2
           ranges(ic,c,p,3) = centre
           ranges(ic,c,p,4) = y2-y1
           ranges(ic,c,p,5) = ranges(ic,c,p,4)
           if(parID(p).eq.21.or.parID(p).eq.22 .or. parID(p).eq.61.or.parID(p).eq.63.or.parID(p).eq.64) ranges(ic,c,p,5) = ranges(ic,c,p,4)/ranges(ic,c,p,3)  !Distance or mass
        end do !p
     end do !c
     !if(prProgress.ge.2) write(6,'(A34,F8.4)')'.  Standard probability interval: ',ivals(ival0)
     !if(prProgress.ge.2) write(6,'(A,F8.4,$)')', default ival:',ivals(ival0)
     
     
     
     
     
     
     
     
     
     !Compute Bayes factor
     !if(prProgress.ge.1.and.ic.eq.1) write(6,'(A,$)')'  Bayes factor,'
     total = 0
     maxlogl = -1.e30
     minlogl =  1.e30
     do i=Nburn(ic),ntot(ic)
        var = post(ic,i)          !Use quadruple precision
        total = total + exp(-var)
        maxlogl = max(post(ic,i),maxlogl)
        minlogl = min(post(ic,i),minlogl)
     end do
     var = dble(n(ic))/total
     logebayesfactor(ic) = real(log(var))
     log10bayesfactor(ic) = real(log10(var))
     !write(6,'(4F10.3,I9)')logebayesfactor(ic),log10bayesfactor(ic),maxlogl,minlogl,n(ic)
     
     
     
     
     
     
     
     
     
     
     
     !**********************************************************************************************
     !******   CHANGE MCMC PARAMETERS   ******************************************************************
     !**********************************************************************************************
     
     
     
     
     !Change MCMC parameters
     if(changeVar.ge.1) then
        if(prProgress.ge.1.and.ic.eq.i.and.update.eq.0) write(6,'(A,$)')'.  Change vars. '
        do p=1,nMCMCpar
           !CHECK: need d^3!
           if(parID(p).eq.22) then !Take exp
              stdev1(p) = stdev1(p)*exp(stats(ic,p,1))  !Median  For exponential function y = exp(x), sig_y = exp(x) sig_x
              stdev2(p) = stdev2(p)*exp(stats(ic,p,2))  !Mean
              selDat(ic,p,1:n(ic)) = exp(selDat(ic,p,1:n(ic)))     !logD -> Distance
              if(ic.eq.1) startval(1:nchains0,p,1:3) = exp(startval(1:nchains0,p,1:3))
              stats(ic,p,1:nstat) = exp(stats(ic,p,1:nstat))
              ranges(ic,1:Nival,p,1:nr) = exp(ranges(ic,1:Nival,p,1:nr))
           end if
           if(parID(p).eq.51.or.parID(p).eq.72.or.parID(p).eq.82) then !cos -> deg
              stdev1(p) = abs(-1./(sqrt(max(1.-stats(ic,p,1)**2,1.e-30))) * stdev1(p))*rr2d  !Based on median
              stdev2(p) = abs(-1./(sqrt(max(1.-stats(ic,p,2)**2,1.e-30))) * stdev2(p))*rr2d  !Based on mean
              selDat(ic,p,1:n(ic)) = acos(selDat(ic,p,1:n(ic)))*rr2d
              if(ic.eq.1) startval(1:nchains0,p,1:3) = acos(startval(1:nchains0,p,1:3))*rr2d
              stats(ic,p,1:nstat) = acos(stats(ic,p,1:nstat))*rr2d
              ranges(ic,1:Nival,p,1:nr) = acos(ranges(ic,1:Nival,p,1:nr))*rr2d
              do c=1,Nival
                 y1 = ranges(ic,c,p,2)
                 ranges(ic,c,p,2) = ranges(ic,c,p,1)  !acos is monotonously decreasing
                 ranges(ic,c,p,1) = y1
              end do
           end if
           if(parID(p).eq.31) then !rad -> h
              stdev1(p) = stdev1(p)*rr2h
              stdev2(p) = stdev2(p)*rr2h
              selDat(ic,p,1:n(ic)) = selDat(ic,p,1:n(ic))*rr2h
              if(ic.eq.1) startval(1:nchains0,p,1:3) = startval(1:nchains0,p,1:3)*rr2h
              stats(ic,p,1:nstat) = stats(ic,p,1:nstat)*rr2h
              ranges(ic,1:Nival,p,1:nr) = ranges(ic,1:Nival,p,1:nr)*rr2h
           end if
           if(parID(p).eq.32.or.parID(p).eq.53) then !sin -> deg
              stdev1(p) = abs(1./(sqrt(max(1.-stats(ic,p,1)**2,1.e-30))) * stdev1(p))*rr2d  !Based on median
              stdev2(p) = abs(1./(sqrt(max(1.-stats(ic,p,2)**2,1.e-30))) * stdev2(p))*rr2d  !Based on mean
              selDat(ic,p,1:n(ic)) = asin(selDat(ic,p,1:n(ic)))*rr2d
              if(ic.eq.1) startval(1:nchains0,p,1:3) = asin(startval(1:nchains0,p,1:3))*rr2d
              stats(ic,p,1:nstat) = asin(stats(ic,p,1:nstat))*rr2d
              ranges(ic,1:Nival,p,1:nr) = asin(ranges(ic,1:Nival,p,1:nr))*rr2d
           end if
           if(parID(p).eq.41.or.parID(p).eq.52.or.parID(p).eq.54.or.parID(p).eq.73.or.parID(p).eq.83) then  !rad -> deg
              stdev1(p) = stdev1(p)*rr2d
              stdev2(p) = stdev2(p)*rr2d
              selDat(ic,p,1:n(ic)) = selDat(ic,p,1:n(ic))*rr2d
              if(ic.eq.1) startval(1:nchains0,p,1:3) = startval(1:nchains0,p,1:3)*rr2d
              stats(ic,p,1:nstat) = stats(ic,p,1:nstat)*rr2d
              ranges(ic,1:Nival,p,1:nr) = ranges(ic,1:Nival,p,1:nr)*rr2d
           end if
           
           ranges(ic,1:Nival,p,3) = 0.5*(ranges(ic,1:Nival,p,1) + ranges(ic,1:Nival,p,2))
           ranges(ic,1:Nival,p,4) = ranges(ic,1:Nival,p,2) - ranges(ic,1:Nival,p,1)
           ranges(ic,1:Nival,p,5) = ranges(ic,1:Nival,p,4)
           if(parID(p).eq.21.or.parID(p).eq.22 .or. parID(p).eq.61.or.parID(p).eq.63.or.parID(p).eq.64) ranges(ic,1:Nival,p,5) = ranges(ic,1:Nival,p,4)/ranges(ic,1:Nival,p,3)  !Distance or mass
        end do !p
        
        !Change the parameter names:
        call set_derivedParameterNames()
        
     end if !if(changeVar.ge.1)
     
     
     !Find 100% probability range
     do c = 1,Nival
        if(abs(ivals(c)-1.).lt.1.e-4) then !Then treat it as a 100% interval to prevent numerical problems
           if(prProgress.ge.1) write(6,'(A,F9.4,A)')'  Treating probability interval',ivals(c)*100,'% as 100%'
           do p=1,nMCMCpar
              ranges(ic,c,p,1) = minval(selDat(ic,p,1:n(ic)))
              ranges(ic,c,p,2) = maxval(selDat(ic,p,1:n(ic)))
              ranges(ic,c,p,3) = 0.5*(ranges(ic,c,p,1) + ranges(ic,c,p,2))
              ranges(ic,c,p,4) = ranges(ic,c,p,2) - ranges(ic,c,p,1)
              ranges(ic,c,p,5) = ranges(ic,c,p,4)
              if(parID(p).eq.21.or.parID(p).eq.22 .or. parID(p).eq.61.or.parID(p).eq.63.or.parID(p).eq.64) ranges(ic,c,p,5) = ranges(ic,c,p,4)/ranges(ic,c,p,3)  !Distance or mass
           end do
        end if
     end do
     
     if(prProgress.ge.1) then
        if(ic.eq.nchains) then
           write(6,*)
        else
           write(6,'(A,$)')'  '
        end if
     end if
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     !**********************************************************************************************
     !******   PRINT STATISTICS   ******************************************************************
     !**********************************************************************************************

     o=6 !Print to this unit - 6=screen
     
     
     if(prProgress+prStat+prIval+prConv.gt.0.and.ic.eq.1) then
        write(o,'(/,A,2(A,F7.2),$)')'  Bayes factor:   ','log_e(B_SN) =',logebayesfactor(ic),',  log_10(B_SN) =',log10bayesfactor(ic)
        !write(o,'(8x,A,3(A,F7.2),A)')'  Maximum likelihood:   ','log_e(Lmax) =',startval(ic,1,3),',  log_10(Lmax) =',startval(ic,1,3)/log(10.),',  -> SNR =',sqrt(2*startval(ic,1,3)),'.'
        write(o,'(8x,A,3(A,F7.2),A)')'  Maximum likelihood:   ','log_e(Lmax) =',loglmax,',  log_10(Lmax) =',loglmax/log(10.),',  -> SNR =',sqrt(2*loglmax),'.'
     end if
     
     
     !Print statistics to screen
     if(prStat.gt.0) then
        write(o,'(/,A)')'  Main statistics:'
        do c=1,Nival
           if(c.ne.c0.and.prStat.lt.2) cycle
           if(c.gt.1.and.prStat.ge.2) write(o,*)
           write(o,'(A10, A12,2A10,A12, 4A8, 4A10,A8,A10, A4,A12,F7.3,A2)')'Param.  ','model','median','mean','Lmax','stdev1','stdev2','abvar1','abvar2',  &
                'rng_c','rng1','rng2','drng','d/drng','delta','ok?','result (',ivals(c)*100,'%)'
           do p=1,nMCMCpar
              if(fixedpar(p).eq.1) cycle !Parameter was not fitted
              write(o,'(A10,F12.6,2F10.4,F12.6, 4F8.4,4F10.4,F8.4,F10.4,$)')parNames(parID(p)),startval(ic,p,1),stats(ic,p,1),stats(ic,p,2),startval(ic,p,3),stdev1(p),stdev2(p),absVar1(p),  &
                   absVar2(p),ranges(ic,c,p,3),ranges(ic,c,p,1),ranges(ic,c,p,2),ranges(ic,c,p,4),  &
                   2*abs(startval(ic,p,1)-ranges(ic,c,p,3))/ranges(ic,c,p,4),ranges(ic,c,p,5)  !d/drange wrt centre of range
              if(startval(ic,p,1).ge.ranges(ic,c,p,1).and.startval(ic,p,1).le.ranges(ic,c,p,2)) then
                 write(o,'(A4,$)')'y '
              else
                 write(o,'(A4,$)')'*N*'
              end if
              write(o,'(F10.4,A3,F9.4)')ranges(ic,c,p,3),'+-',0.5*ranges(ic,c,p,4)
           end do !p
        end do !c
     end if
     
     
     !Print intervals as: centre, delta, in range:
     if(prIval.eq.1.or.prIval.eq.3) then
        write(o,'(/,A)')'  Probability intervals:'
        write(o,'(A22,A8,$)')'Interval:',''
        do c=1,Nival
           write(o,'(F20.4,A9,$)')ivals(c),''
        end do
        write(o,*)''
        
        write(o,'(A10,2x,2A9,$)')'Param.  ','model','median'
        do c=1,Nival
           !write(o,'(2x,2A9,A8,$)')'rng1','rng2','in rnge'
           write(o,'(2x,3A9,$)')'centre','delta','in rnge'
        end do
        write(o,*)''
        do p=1,nMCMCpar
           !if(stdev1(p).lt.1.d-20) cycle
           if(fixedpar(p).eq.1) cycle
           write(o,'(A10,2x,2F9.4,$)')parNames(parID(p)),startval(ic,p,1),stats(ic,p,1)
           do c=1,Nival
              if(mergeChains.eq.0) then
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
     
     
     if(prIval.ge.2) then
        write(o,'(/,A)')'  Statistics and probability intervals:'
        
        !Print intervals as x +- dx:
        write(o,'(A61,$)')'Interval:'
        do c=1,Nival
           write(o,'(F14.3,A1,A9,$)')ivals(c)*100,'%',''
        end do
        write(o,*)''
        
        write(o,'(A11,1x,4A10,$)')'Parameter','median','mean','Lmax','stdev'
        do c=1,Nival
           write(o,'(5x,A8,4x,A7,$)')'x','dx'
        end do
        write(o,*)''
        do p=1,nMCMCpar
           if(fixedpar(p).eq.1) cycle
           write(o,'(A10,2x,4F10.3,$)')parNames(parID(p)),stats(ic,p,1),stats(ic,p,2),startval(ic,p,3),stdev2(p)
           do c=1,Nival
              write(o,'(5x,F8.3,A4,F7.3$)')ranges(ic,c,p,3),' +- ',0.5d0*ranges(ic,c,p,4)
           end do
           write(o,*)''
        end do
        write(o,*)''
        
        !Print intervals as low-high:
        write(o,'(A61,$)')'Interval:'
        do c=1,Nival
           write(o,'(F14.3,A1,A9,$)')ivals(c)*100,'%',''
        end do
        write(o,*)''
        
        write(o,'(A11,1x,4A10,$)')'Parameter','median','mean','Lmax','stdev'
        do c=1,Nival
           write(o,'(6x,A8,3x,A7,$)')'min','max'
        end do
        write(o,*)''
        do p=1,nMCMCpar
           if(fixedpar(p).eq.1) cycle
           write(o,'(A10,2x,4F10.3,$)')parNames(parID(p)),stats(ic,p,1),stats(ic,p,2),startval(ic,p,3),stdev2(p)
           do c=1,Nival
              write(o,'(6x,F8.3,A3,F7.3$)')ranges(ic,c,p,1),' - ',ranges(ic,c,p,2)
           end do
           write(o,*)''
        end do
        write(o,*)''
        
        
        !Print intervals as x -dx1 +dx2 in LaTeX format:
        if(1.eq.2) then
           write(o,'(A19,$)')'Interval:   '
           do c=1,Nival
              write(o,'(F27.3,A1,A9,$)')ivals(c)*100,'%',''
           end do
           write(o,*)''
           
           do p=1,nMCMCpar
              if(fixedpar(p).eq.1) cycle
              write(o,'(A10,2x,$)')parNames(parID(p))
              do c=1,Nival
                 write(o,'(5x,A1,F8.3,2(A,F7.3),A,$)')'$',stats(ic,p,1),'_{-',abs(stats(ic,p,1)-ranges(ic,c,p,1)),'}^{+',abs(stats(ic,p,1)-ranges(ic,c,p,2)),'}$'
              end do
              write(o,*)''
           end do
           write(o,*)''
        end if
        
     end if
     
     
     
     
     
     
     
     
     
     
     !Print correlations:
     if(prCorr.gt.0) then
        corr1 = 0.1
        corr2 = 0.5
        write(o,'(/,A,$)')'  Correlations  '
        write(o,'(A,3(F4.2,A))')'  (weak [',corr1,'<abs(cor)<',corr2,']: in lower triangle,  strong [abs(cor)>',corr2,']: in upper triangle):'
        write(o,'(A8,$)')''
        do p=1,nMCMCpar
           if(fixedpar(p).eq.0) write(o,'(A7,$)')trim(parNames(parID(p)))
        end do
        write(o,*)''
        do p1=1,nMCMCpar
           if(fixedpar(p1).eq.1) cycle
           write(o,'(A8,$)')trim(parNames(parID(p1)))
           do p2=1,nMCMCpar
              corr = corrs(p1,p2)
              if(fixedpar(p2).eq.1) cycle
              if( (abs(corr).ge.corr2.and.p2.gt.p1) ) then  !Print in the upper triangle
                 write(o,'(F7.2,$)')corr
              else if( (abs(corr).ge.corr1.and.abs(corr).lt.corr2.and.p1.gt.p2) ) then   !Print in the lower triangle
                 write(o,'(F7.2,$)')corr
              else if(p1.eq.p2) then !Print on the diagonal
                 write(o,'(A7,$)')' ######' !'
              else
                 write(o,'(A7,$)')''
              end if
           end do
           write(o,'(A)')'   '//trim(parNames(parID(p1)))
        end do
     end if
     
     
     
     
     
     !Print output for CBC Wiki:
     if(ic.eq.1.and.wikioutput.eq.1) call save_cbc_wiki_data(ic)
     
     
  end do !ic
  
  
  
  
  
  
  
  
  !Compute and print convergence:
  if(nchains0.gt.1 .and. (prConv.ge.1.or.saveStats.ge.1)) call compute_convergence()  !Need unwrapped data for this (?)
  
  
  !Change the original chain data:
  if(changeVar.ge.1) then
     do ic=1,nchains0
        do p=1,nMCMCpar
           if(parID(p).eq.21) allDat(ic,p,1:ntot(ic)) = allDat(ic,p,1:ntot(ic))**c3rd  !^(1/3)
           if(parID(p).eq.22) allDat(ic,p,1:ntot(ic)) = exp(allDat(ic,p,1:ntot(ic)))   !exp
           if(parID(p).eq.51.or.parID(p).eq.72.or.parID(p).eq.82) allDat(ic,p,1:ntot(ic)) = acos(allDat(ic,p,1:ntot(ic)))*r2d  !acos -> deg
           if(parID(p).eq.31) allDat(ic,p,1:ntot(ic)) = allDat(ic,p,1:ntot(ic))*r2h  !rad -> h
           if(parID(p).eq.32.or.parID(p).eq.53) allDat(ic,p,1:ntot(ic)) = asin(allDat(ic,p,1:ntot(ic)))*r2d  !asin -> deg
           if(parID(p).eq.41.or.parID(p).eq.52.or.parID(p).eq.54.or.parID(p).eq.73.or.parID(p).eq.83) allDat(ic,p,1:ntot(ic)) = allDat(ic,p,1:ntot(ic))*r2d  !rad -> deg
        end do !p
     end do
  end if !if(changeVar.ge.1)
  
  
  
  
  
end subroutine statistics
!***********************************************************************************************************************************







!***********************************************************************************************************************************
subroutine save_stats(exitcode)  !Save statistics to file  
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
  write(o,'(A,2I3)')'Npar,ncol:',nMCMCpar,7
  write(o,'(A8,7A12)')'Param.  ','model','median','mean','stdev1','stdev2','abvar1','abvar2'
  
  do p=1,nMCMCpar
     if(fixedpar(p).eq.0) then
        write(o,'(A8,7F12.6)')parNames(parID(p)),startval(ic,p,1),stats(ic,p,1),stats(ic,p,2),stdev1(p),stdev2(p),absVar1(p),absVar2(p)
     else
        write(o,'(A8,7F12.6)')parNames(parID(p)),startval(ic,p,1),stats(ic,p,1),stats(ic,p,2),0.,0.,0.,0.
     end if
  end do
  write(o,*)''
  
  
  !Print correlations:
  write(o,'(//,A,/)')'CORRELATIONS:'
  write(o,'(A,I3)')'Npar:',nMCMCpar
  write(o,'(A9,$)')''
  do p=1,nMCMCpar
     write(o,'(A10,$)')trim(parNames(parID(p)))
  end do
  write(o,*)''
  do p1=1,nMCMCpar
     write(o,'(A9,$)')trim(parNames(parID(p1)))
     do p2=1,nMCMCpar
        write(o,'(F10.5,$)')corrs(p1,p2)
     end do
     write(o,'(A)')'   '//trim(parNames(parID(p1)))
  end do
  
  
  !Print probability intervals:
  write(o,'(///,A,/)')'1D PROBABILITY INTERVALS:'
  write(o,'(A,I3)')'Nival:',Nival
  write(o,'(A22,$)')'Interval:'
  do c=1,Nival
     write(o,'(F21.5,A14,$)')ivals(c),''
  end do
  write(o,*)''
  
  write(o,'(A8,2x,$)')'param.'
  do c=1,Nival
     !write(o,'(2x,2A9,A8,$)')'rng1','rng2','in rnge'
     write(o,'(2x,2A12,A9,$)')'centre','delta','in rnge'
  end do
  write(o,*)''
  do p=1,nMCMCpar
     write(o,'(A8,2x,$)')trim(parNames(parID(p)))
     do c=1,Nival
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
  write(o,'(A,I5)')'Npdf2D: ',Npdf2D
  write(o,'(A,2I5)')'Nbin2Dx,Nbin2Dy: ',Nbin2Dx,Nbin2Dy
  !write(o,*)''
  write(o,'(A28,$)')'Interval:'
  do c=1,Nival
     write(o,'(F19.5,$)')ivals(c)
  end do
  write(o,*)''
  
  write(o,'(A9,A19,$)')'params.',''
  do c=1,Nival
     write(o,'(A16,A3,$)')'delta','in'
  end do
  write(o,*)''
  do p=1,Npdf2D
     p1 = PDF2Dpairs(p,1)
     p2 = PDF2Dpairs(p,2)
     if(parID(p1)*parID(p2).eq.0) then
        if(parID(p1).eq.0) write(0,'(/,A,I4,A,/,A,//)')'  ***  ERROR:  save_stats():  parameter',p1,' not defined, check PDF2Dpairs in the input file ***','  Aborting...'
        if(parID(p2).eq.0) write(0,'(/,A,I4,A,/,A,//)')'  ***  ERROR:  save_stats():  parameter',p2,' not defined, check PDF2Dpairs in the input file ***','  Aborting...'
        stop
     end if
     write(o,'(2I4,2(2x,A8),2x,$)')p1,p2,trim(parNames(parID(p1))),trim(parNames(parID(p2)))
     do c=1,Nival
        write(o,'(2x,F14.8,$)')probareas(p1,p2,c,1)
        if(c.ge.Nival+1-trueranges2d(p1,p2) .and. trueranges2d(p1,p2).ne.0) then
           write(o,'(A3,$)')'y'
        else
           write(o,'(A3,$)')'n'
        end if
     end do
     write(o,*)''
  end do
  
  
  
  close(o) !Statistics output file
  if(saveStats.eq.2) i = system('a2ps -1rf7 '//trim(outputdir)//'/'//trim(outputname)//'__statistics.dat -o '//trim(outputdir)//'/'//trim(outputname)//'__statistics.ps')
  !write(6,*)''
  if(prProgress.ge.1) then
     if(saveStats.eq.1) write(6,'(A)')'  Statistics saved in '//trim(outputname)//'__statistics.dat'
     if(saveStats.eq.2) write(6,'(A)')'  Statistics saved in '//trim(outputname)//'__statistics.dat,ps'
  end if
  
end subroutine save_stats
!***********************************************************************************************************************************





!***********************************************************************************************************************************
subroutine save_cbc_wiki_data(ic)
  use constants
  use analysemcmc_settings
  use general_data
  use stats_data
  use chain_data
  use mcmcrun_data
  implicit none
  
  integer :: c,i,ic,io,o,p,p1,parr(maxMCMCpar)
  real :: x,rev2pi,rrevpi,x1,x2
  character :: url*99,gps*19,xs11*11,xs20*20,ans,wikifilename*99,pnstr*3
  
  write(pnstr,'(F3.1)')pnOrder
  
  !Print output for CBC Wiki:
  o = 10
  write(wikifilename,'(A)')trim(outputname)//'__wiki.dat'
  open(unit=o,form='formatted',status='replace',action='write',position='rewind',file=trim(wikifilename),iostat=io)
  if(io.ne.0) then
     write(0,'(A)')'  Error opening '//trim(wikifilename)//', aborting...'
     stop
  end if
  write(gps,'(I10.10)')GPStime
  if(GPStime.lt.1e9) write(gps,'(I9.9)')GPStime
  
  write(url,'(A)')'http://www.astro.northwestern.edu/~lsc/E14/GPS'//trim(gps)//'/'
  write(o,'(A)')'= GPS'//trim(gps)//' - description ='
  write(o,'(/,A)')'Back to [:JointS5/BayesianFollowUpOfE14Events:Bayesian follow-up in E14]'
  
  
  
  !Injection values:
  write(o,'(///,A)')'== Injection values =='
  write(o,'(A)')"|| '''Detectors'''  || '''M1'''     || '''M2'''     || '''Mc'''     || '''&eta;'''  || '''time'''      || '''spin1'''  ||'''&theta;1'''|| '''spin2'''  ||'''&theta;2'''|| '''Dist'''   || '''R.A.'''   || '''Dec.'''   || '''incl'''   || '''pol.'''   || '''details'''                                                                                     ||"
  write(o,'(A)')"||                  ||  (Mo)        || (Mo)         || (Mo)         ||              ||  (s)            ||              || (rad)        ||              || (rad)        || (Mpc)        || (rad)        || (rad)        || (rad)        || (rad)        ||                                                                                                   ||"
  
  write(o,'(A4,$)')'||  '
  do i=1,4
     if(i.le.ndet(ic)) then
        write(o,'(A2,$)')detabbrs(detnr(ic,i))
     else
        write(o,'(A2,$)')'  '
     end if
  end do
  write(o,'(A6,$)')'    '
  parr(1:14) = (/63,64,61,62,11,71,72,81,82,22,31,32,51,52/)
  do p=1,14
     p1 = parr(p)
     if(revID(p1).eq.0) then  !Parameter not used
        write(o,'(A5,A11,$)')'  || ',' - '
     else
        x = allDat(ic,revID(p1),1)
        if(p1.eq.31) x = rev2pi(x*rh2r)
        if(p1.eq.52) x = rrevpi(x*rd2r)  !Polarisation angle
        if(p1.eq.41.or.p1.eq.54.or.p1.eq.73.or.p1.eq.83) x = rev2pi(x*rd2r)
        if(p1.eq.32.or.p1.eq.51.or.p1.eq.53.or.p1.eq.72.or.p1.eq.82) x = x*rd2r
        if(p1.ge.11.and.p1.le.19) then
           write(o,'(A5,F14.4,$)')'  || ',x+t0
        else
           write(o,'(A5,F11.4,$)')'  || ',x
        end if
     end if
  end do
  write(o,'(A,79x,A)')'  || [ injection info] ',' ||'
  
  
  
  
  !Bayes factor:
  write(o,'(///,A)')'== Bayes Factors =='
  write(o,'(A)')"|| '''Code'''                                     || '''Model'''                                                || '''Detectors'''  || '''log_e Bayes Factor'''    || '''log_10 Bayes Factor'''    || '''Details'''                                                           ||"
  
  write(o,'(A3,A47,$)')'|| ','[http://tinyurl.com/SPINspiral SPINspiral]     '
  if(spinningRun.eq.0) write(o,'(A3,A59,$)')'|| ',trim(waveforms(waveform))//' '//pnstr//'pN non-spinning vs. Gaussian noise '
  if(spinningRun.eq.1) write(o,'(A3,A59,$)')'|| ',trim(waveforms(waveform))//' '//pnstr//'pN one spin vs. Gaussian noise     '
  if(spinningRun.eq.2) write(o,'(A3,A59,$)')'|| ',trim(waveforms(waveform))//' '//pnstr//'pN two spins vs. Gaussian noise    '
  write(o,'(A,$)')'||  '
  do i=1,4
     if(i.le.ndet(ic)) then
        write(o,'(A2,$)')detabbrs(detnr(ic,i))
     else
        write(o,'(A2,$)')'  '
     end if
  end do
  write(o,'(A11,$)')'        || '
  write(o,'(F10.1,A,$)')logebayesfactor(ic),'                  || '
  write(o,'(F10.1,A,$)')log10bayesfactor(ic),'                   || '
  write(o,'(A)')'['//trim(url)//' link]         ||'
  
  
  
  
  !Medians:
  write(o,'(///,A)')'== Medians =='
  write(o,'(A)')"|| '''Code'''                                     || '''Waveform'''                                 || '''Detectors'''  || '''Mc'''     || '''&eta;'''  || '''time'''   || '''spin1'''  ||'''&theta;1'''|| '''spin2'''  ||'''&theta;2'''|| '''Dist'''   || '''R.A.'''   || '''Dec.'''   || '''incl'''   || '''pol.'''   || '''details'''                                                                                     ||"
  write(o,'(A)')"||                                                ||                                                ||                  || (Mo)         ||              ||  (s)         ||              || (rad)        ||              || (rad)        || (Mpc)        || (rad)        || (rad)        || (rad)        || (rad)        ||                                                                                                   ||"
  write(o,'(A3,A47,$)')'|| ','[http://tinyurl.com/SPINspiral SPINspiral]     '
  if(spinningRun.eq.0) write(o,'(A3,A47,$)')'|| ',trim(waveforms(waveform))//' '//pnstr//'pN non-spinning '
  if(spinningRun.eq.1) write(o,'(A3,A47,$)')'|| ',trim(waveforms(waveform))//' '//pnstr//'pN one spin     '
  if(spinningRun.eq.2) write(o,'(A3,A47,$)')'|| ',trim(waveforms(waveform))//' '//pnstr//'pN two spins    '
  write(o,'(A,$)')'||  '
  do i=1,4
     if(i.le.ndet(ic)) then
        write(o,'(A2,$)')detabbrs(detnr(ic,i))
     else
        write(o,'(A2,$)')'  '
     end if
  end do
  write(o,'(A10,$)')'       ||'
  !parr(1:10) = (/2,3,4,6,7,5,8,9,11,12/)
  parr(1:12) = (/61,62,11,71,72,81,82,22,31,32,51,52/)
  do p=1,12
     p1 = parr(p)
     if(revID(p1).eq.0) then  !Parameter not used
        xs11 = ' - '
     else
        x = stats(ic,revID(p1),1)
        if(p1.eq.31) x = rev2pi(x*rh2r)
        if(p1.eq.52) x = rrevpi(x*rd2r)  !Polarisation angle
        if(p1.eq.41.or.p1.eq.54.or.p1.eq.73.or.p1.eq.83) x = rev2pi(x*rd2r)
        if(p1.eq.32.or.p1.eq.51.or.p1.eq.53.or.p1.eq.72.or.p1.eq.82) x = x*rd2r
        write(xs11,'(F11.4)')x
     end if
     write(o,'(A11,A5,$)')xs11,'   ||'
  end do
  write(o,'(A)')' ['//trim(url)//' link]                                   ||'
  
  
  
  
  !Means:
  write(o,'(///,A)')'== Means =='
  write(o,'(A)')"|| '''Code'''                                     || '''Waveform'''                                 || '''Detectors'''  || '''Mc'''     || '''&eta;'''  || '''time'''   || '''spin1'''  ||'''&theta;1'''|| '''spin2'''  ||'''&theta;2'''|| '''Dist'''   || '''R.A.'''   || '''Dec.'''   || '''incl'''   || '''pol.'''   || '''details'''                                                                                     ||"
  write(o,'(A)')"||                                                ||                                                ||                  || (Mo)         ||              ||  (s)         ||              || (rad)        ||              || (rad)        || (Mpc)        || (rad)        || (rad)        || (rad)        || (rad)        ||                                                                                                   ||"
  write(o,'(A3,A47,$)')'|| ','[http://tinyurl.com/SPINspiral SPINspiral]     '
  if(spinningRun.eq.0) write(o,'(A3,A47,$)')'|| ',trim(waveforms(waveform))//' '//pnstr//'pN non-spinning '
  if(spinningRun.eq.1) write(o,'(A3,A47,$)')'|| ',trim(waveforms(waveform))//' '//pnstr//'pN one spin     '
  if(spinningRun.eq.2) write(o,'(A3,A47,$)')'|| ',trim(waveforms(waveform))//' '//pnstr//'pN two spins    '
  write(o,'(A,$)')'||  '
  do i=1,4
     if(i.le.ndet(ic)) then
        write(o,'(A2,$)')detabbrs(detnr(ic,i))
     else
        write(o,'(A2,$)')'  '
     end if
  end do
  write(o,'(A10,$)')'       ||'
  parr(1:12) = (/61,62,11,71,72,81,82,22,31,32,51,52/)
  do p=1,12
     p1 = parr(p)
     if(revID(p1).eq.0) then  !Parameter not used
        xs11 = ' - '
     else
        x = stats(ic,revID(p1),2)
        if(p1.eq.31) x = rev2pi(x*rh2r)
        if(p1.eq.52) x = rrevpi(x*rd2r)  !Polarisation angle
        if(p1.eq.41.or.p1.eq.54.or.p1.eq.73.or.p1.eq.83) x = rev2pi(x*rd2r)
        if(p1.eq.32.or.p1.eq.51.or.p1.eq.53.or.p1.eq.72.or.p1.eq.82) x = x*rd2r
        write(xs11,'(F11.4)')x
     end if
     write(o,'(A11,A5,$)')xs11,'   ||'
  end do
  write(o,'(A)')' ['//trim(url)//' link]                                   ||'
  
  
  
  
  !Lmax:
  write(o,'(///,A)')'== Maximum-likelihood points =='
  write(o,'(A)')"|| '''Code'''                                     || '''Waveform'''                                 || '''Detectors'''  ||'''log(L)'''  || '''Mc'''     || '''&eta;'''  || '''time'''   || '''spin1'''  ||'''&theta;1'''|| '''spin2'''  ||'''&theta;2'''|| '''Dist'''   || '''R.A.'''   || '''Dec.'''   || '''incl'''   || '''pol.'''   || '''details'''                                                                                     ||"
  write(o,'(A)')"||                                                ||                                                ||                  ||              || (Mo)         ||              ||  (s)         ||              || (rad)        ||              || (rad)        || (Mpc)        || (rad)        || (rad)        || (rad)        || (rad)        ||                                                                                                   ||"
  write(o,'(A3,A47,$)')'|| ','[http://tinyurl.com/SPINspiral SPINspiral]     '
  if(spinningRun.eq.0) write(o,'(A3,A47,$)')'|| ',trim(waveforms(waveform))//' '//pnstr//'pN non-spinning '
  if(spinningRun.eq.1) write(o,'(A3,A47,$)')'|| ',trim(waveforms(waveform))//' '//pnstr//'pN one spin     '
  if(spinningRun.eq.2) write(o,'(A3,A47,$)')'|| ',trim(waveforms(waveform))//' '//pnstr//'pN two spins    '
  write(o,'(A,$)')'||  '
  do i=1,4
     if(i.le.ndet(ic)) then
        write(o,'(A2,$)')detabbrs(detnr(ic,i))
     else
        write(o,'(A2,$)')'  '
     end if
  end do
  write(o,'(A10,$)')'       ||'
  
  !Print log L:
  write(xs11,'(F11.4)')logLmax
  write(o,'(A11,A5,$)')xs11,'   ||'
  
  parr(1:12) = (/61,62,11,71,72,81,82,22,31,32,51,52/)
  do p=1,12
     p1 = parr(p)
     if(revID(p1).eq.0) then  !Parameter not used
        xs11 = ' - '
     else
        x = startval(ic,revID(p1),3)
        if(p1.eq.31) x = rev2pi(x*rh2r)
        if(p1.eq.52) x = rrevpi(x*rd2r)  !Polarisation angle
        if(p1.eq.41.or.p1.eq.54.or.p1.eq.73.or.p1.eq.83) x = rev2pi(x*rd2r)
        if(p1.eq.32.or.p1.eq.51.or.p1.eq.53.or.p1.eq.72.or.p1.eq.82) x = x*rd2r
        write(xs11,'(F11.4)')x
     end if
     write(o,'(A11,A5,$)')xs11,'   ||'
  end do
  write(o,'(A)')' ['//trim(url)//' link]                                   ||'
  
  
  
  
  !Stdev:
  write(o,'(///,A)')'== Standard deviations =='
  write(o,'(A)')"|| '''Code'''                                     || '''Waveform'''                                 || '''Detectors'''  || '''Mc'''     || '''&eta;'''  || '''time'''   || '''spin1'''  ||'''&theta;1'''|| '''spin2'''  ||'''&theta;2'''|| '''Dist'''   || '''R.A.'''   || '''Dec.'''   || '''incl'''   || '''pol.'''   || '''details'''                                                                                     ||"
  write(o,'(A)')"||                                                ||                                                ||                  || (Mo)         ||              ||  (s)         ||              || (rad)        ||              || (rad)        || (Mpc)        || (rad)        || (rad)        || (rad)        || (rad)        ||                                                                                                   ||"
  write(o,'(A3,A47,$)')'|| ','[http://tinyurl.com/SPINspiral SPINspiral]     '
  if(spinningRun.eq.0) write(o,'(A3,A47,$)')'|| ',trim(waveforms(waveform))//' '//pnstr//'pN non-spinning '
  if(spinningRun.eq.1) write(o,'(A3,A47,$)')'|| ',trim(waveforms(waveform))//' '//pnstr//'pN one spin     '
  if(spinningRun.eq.2) write(o,'(A3,A47,$)')'|| ',trim(waveforms(waveform))//' '//pnstr//'pN two spins    '
  write(o,'(A,$)')'||  '
  do i=1,4
     if(i.le.ndet(ic)) then
        write(o,'(A2,$)')detabbrs(detnr(ic,i))
     else
        write(o,'(A2,$)')'  '
     end if
  end do
  write(o,'(A10,$)')'       ||'
  parr(1:12) = (/61,62,11,71,72,81,82,22,31,32,51,52/)
  do p=1,12
     p1 = parr(p)
     if(revID(p1).eq.0) then  !Parameter not used
        xs11 = ' - '
     else
        x = stdev2(revID(p1))
        if(p1.eq.31) x = x*rh2r
        if(p1.eq.52) x = rrevpi(x*rd2r)  !Polarisation angle
        if(p1.eq.41.or.p1.eq.54.or.p1.eq.73.or.p1.eq.83) x = x*rd2r
        if(p1.eq.32.or.p1.eq.51.or.p1.eq.53.or.p1.eq.72.or.p1.eq.82) x = x*rd2r
        write(xs11,'(F11.4)')x
     end if
     write(o,'(A11,A5,$)')xs11,'   ||'
  end do
  write(o,'(A)')' ['//trim(url)//' link]                                   ||'
  
  
  
  !2-sigma range:
  write(o,'(///,A)')'== 2-sigma probability ranges =='
  write(o,'(A)')"|| '''Code'''                                     || '''Waveform'''                                 || '''Detectors'''  || '''Mc'''              || '''&eta;'''           || '''time'''            || '''spin1'''           || '''&theta;1'''        || '''spin2'''           || '''&theta;2'''        || '''Distance'''        || '''R.A.'''            || '''Dec.'''            || '''incl.'''           || '''pol.'''            || '''details'''                                                                                     ||"
  write(o,'(A)')"||                                                ||                                                ||                  || (Mo)                  ||                       ||  (s)                  ||                       || (rad)                 ||                       || (rad)                 || (Mpc)                 || (rad)                 || (rad)                 || (rad)                 || (rad)                 ||                                                                                                   ||"
  c = 0
  do i=1,Nival
     if(abs(ivals(i)-0.9545).lt.0.0001) c = i
  end do
  if(c.eq.0) then
     write(0,'(A)')'  Error: 2-sigma range not found, needed for Wiki output!'
     write(6,'(A,$)')'  Do you want to continue?  (y/n)  '
     read(5,*)ans
     if(ans.eq.'y'.or.ans.eq.'Y') then
        c = 1
        write(6,'(A,F6.2,A)')'  Continuing with ',ivals(c)*100,"% probability interval, don't use wiki.txt!!!"
     else
        stop
     end if
  end if
  write(o,'(A3,A47,$)')'|| ','[http://tinyurl.com/SPINspiral SPINspiral]     '
  if(spinningRun.eq.0) write(o,'(A3,A47,$)')'|| ',trim(waveforms(waveform))//' '//pnstr//'pN non-spinning '
  if(spinningRun.eq.1) write(o,'(A3,A47,$)')'|| ',trim(waveforms(waveform))//' '//pnstr//'pN one spin     '
  if(spinningRun.eq.2) write(o,'(A3,A47,$)')'|| ',trim(waveforms(waveform))//' '//pnstr//'pN two spins    '
  write(o,'(A,$)')'||  '
  do i=1,4
     if(i.le.ndet(ic)) then
        write(o,'(A2,$)')detabbrs(detnr(ic,i))
     else
        write(o,'(A2,$)')'  '
     end if
  end do
  write(o,'(A10,$)')'       ||'
  parr(1:12) = (/61,62,11,71,72,81,82,22,31,32,51,52/)
  do p=1,12
     p1 = parr(p)
     if(revID(p1).eq.0) then  !Parameter not used
        xs20 = ' - '
     else
        x1 = ranges(ic,c,revID(p1),1)
        x2 = ranges(ic,c,revID(p1),2)
        if(p1.eq.31) then  !RA
           x1 = x1*rh2r
           x2 = x2*rh2r
        end if
        if(p1.eq.41.or.p1.eq.54.or.p1.eq.73.or.p1.eq.83) then
           x1 = x1*rd2r
           x2 = x2*rd2r
        end if
        if(p1.eq.32.or.p1.eq.51.or.p1.eq.52.or.p1.eq.53.or.p1.eq.72.or.p1.eq.82) then
           x1 = x1*rd2r
           x2 = x2*rd2r
        end if
        write(xs20,'(F9.4,A2,F9.4)')x1,' -',x2
     end if
     write(o,'(A20,A5,$)')xs20,'   ||'
     
  end do
  write(o,'(A)')' ['//trim(url)//' link]                                   ||'
  
  
  
  !Injection values:
  write(o,'(///,A)')'== Injection values =='
  write(o,'(A)')"|| '''Detectors'''  || '''M1'''     || '''M2'''     || '''Mc'''     || '''&eta;'''  || '''time'''      || '''spin1'''  ||'''&theta;1'''|| '''spin2'''  ||'''&theta;2'''|| '''Dist'''   || '''R.A.'''   || '''Dec.'''   || '''incl'''   || '''pol.'''   || '''details'''                                                                                     ||"
  write(o,'(A)')"||                  ||  (Mo)        || (Mo)         || (Mo)         ||              ||  (s)            ||              || (rad)        ||              || (rad)        || (Mpc)        || (rad)        || (rad)        || (rad)        || (rad)        ||                                                                                                   ||"
  
  write(o,'(A4,$)')'||  '
  do i=1,4
     if(i.le.ndet(ic)) then
        write(o,'(A2,$)')detabbrs(detnr(ic,i))
     else
        write(o,'(A2,$)')'  '
     end if
  end do
  write(o,'(A6,$)')'    '
  parr(1:14) = (/63,64,61,62,11,71,72,81,82,22,31,32,51,52/)
  do p=1,14
     p1 = parr(p)
     if(revID(p1).eq.0) then  !Parameter not used
        write(o,'(A5,A11,$)')'  || ',' - '
     else
        x = allDat(ic,revID(p1),1)
        if(p1.eq.31) x = rev2pi(x*rh2r)
        if(p1.eq.52) x = rrevpi(x*rd2r)  !Polarisation angle
        if(p1.eq.41.or.p1.eq.54.or.p1.eq.73.or.p1.eq.83) x = rev2pi(x*rd2r)
        if(p1.eq.32.or.p1.eq.51.or.p1.eq.53.or.p1.eq.72.or.p1.eq.82) x = x*rd2r
        if(p1.ge.11.and.p1.le.19) then
           write(o,'(A5,F14.4,$)')'  || ',x+t0
        else
           write(o,'(A5,F11.4,$)')'  || ',x
        end if
     end if
  end do
  write(o,'(A,79x,A)')'  || [ injection info] ',' ||'
  
  
  
  
  
  
  write(o,'(///,A)')'----'
  write(o,'(A)')'Back to [:JointS5/BayesianFollowUpOfE14Events:Bayesian follow-up in E14]'
  
  close(o)
  
end subroutine save_cbc_wiki_data
!***********************************************************************************************************************************







!***********************************************************************************************************************************
!>
!!  Check convergence for multiple chains. This works only for fixed chain length, so take the min N
!!  This is probably taken from (at least equal to):
!!    Brooks & Gelman, Journal of Computational and Graphical Statistics, Vol.7, Nr.4, p.434-455, 1998:
!!    "General Methods for Monitoring Convergence of Iterative Simulations"
!!    http://www.jstor.org/pss/1390675 (for purchase)
!!    http://www.stat.columbia.edu/~gelman/research/published/brooksgelman.pdf (Author's website)
!!    See Eq.1.1,  where:  B/n := meanVar  and  W := chVar
!!  \todo:  use only data selected after (auto)burnin
!!  \todo:  do we need unwrapped data for this?
!<
!***********************************************************************************************************************************
subroutine compute_convergence()
  use constants
  use analysemcmc_settings
  use general_data
  use stats_data
  use chain_data
  use mcmcrun_data
  
  implicit none
  integer :: i,ic,p
  integer :: nn,lowVar(maxMCMCpar),nLowVar,highVar(maxMCMCpar),nHighVar,nmeanRelVar,nRhat,IDs(maxMCMCpar),nUsedPar
  real :: dx
  real*8 :: chmean(maxChs,maxMCMCpar),avgMean(maxMCMCpar),chVar(maxMCMCpar),chVar1(maxChs,maxMCMCpar),meanVar(maxMCMCpar),totRhat,meanRelVar
  character :: ch
  
  
  !Compute the means for each chain and for all chains:
  chmean = 1.d-30
  avgMean = 1.d-30
  nn = minval(ntot(1:nchains0))/2
  do p=1,nMCMCpar
     do ic=1,nchains0
        do i=nn+1,2*nn
           chmean(ic,p) = chmean(ic,p) + allDat(ic,p,i) !We used to take unwrapped data for this...
        end do
        avgMean(p) = avgMean(p) + chmean(ic,p)
     end do
  end do
  chmean = chmean/dble(nn)
  avgMean = avgMean/dble(nn*nchains0)
  
  
  !Compute variances per chain, for all chains and Rhat:
  chVar = 1.d-30
  chVar1 = 1.d-30
  meanVar = 1.d-30
  do p=1,nMCMCpar
     do ic=1,nchains0
        do i=nn+1,2*nn
           dx = (allDat(ic,p,i) - chmean(ic,p))**2 !We used to take unwrapped data for this...
           chVar(p) = chVar(p) + dx
           chVar1(ic,p) = chVar1(ic,p) + dx !Keep track of the variance per chain
        end do
        meanVar(p) = meanVar(p) + (chmean(ic,p) - avgMean(p))**2
        chVar1(ic,p) = chVar1(ic,p)/dble(nn-1)
     end do
     chVar(p) = chVar(p)/dble(nchains0*(nn-1))
     meanVar(p) = meanVar(p)/dble(nchains0-1)
     
     Rhat(p) = min( dble(nn-1)/dble(nn)  +  meanVar(p)/chVar(p) * (1.d0 + 1.d0/dble(nchains0)), 99.d0)
  end do
  
  
  !Print means per chain:
  if(prConv.ge.1) then
     write(6,*)''
     if(prConv.ge.2) write(6,'(A,I7,A)')'  Convergence parameters for',nn,' iterations:'
     write(6,'(A14,$)')''
     do p=1,nMCMCpar
        if(fixedpar(p).eq.1) cycle
        write(6,'(A9,$)')trim(parNames(parID(p)))
     end do
     write(6,'(A9)')'Mean'
     
     if(prConv.ge.3) then
        write(6,'(A)')'  Means:'
        do ic=1,nchains0
           write(6,'(I12,A2,$)')ic,': '
           do p=1,nMCMCpar
              if(fixedpar(p).eq.1) cycle
              write(6,'(F9.5,$)')chmean(ic,p)
           end do
           write(6,*)
        end do
        write(6,'(A14,$)')'        Mean: '
        do p=1,nMCMCpar
           if(fixedpar(p).eq.1) cycle
           write(6,'(F9.5,$)')avgMean(p)
        end do
        write(6,*)''
        
     end if !if(prConv.ge.3)
  end if !if(prConv.ge.1)
  
  
  !Flag and print variances:
  if(prConv.ge.3) then
     write(6,*)''
     write(6,'(A)')'  Variances:'
  end if
  do ic=1,nchains0
     lowVar = 0
     highVar = 0
     meanRelVar = 1.d0
     nmeanRelVar = 0
     do p=1,nMCMCpar
        if(fixedpar(p).eq.0 .and.abs(chVar1(ic,p)).gt.1.e-30) then  !Take only the parameters that were fitted and have a variance > 0
           if(chVar1(ic,p).lt.0.5*chVar(p)) lowVar(p) = 1  !Too (?) low variance, mark it
           if(chVar1(ic,p).gt.2*chVar(p))  highVar(p) = 1  !Too (?) high variance, mark it
           meanRelVar = meanRelVar * chVar1(ic,p)/chVar(p) !Take geometric mean
           nmeanRelVar = nmeanRelVar + 1
        end if
     end do
     
     !Find and flag extraordinarily low and high variances:
     IDs(1:4) = (/61,62,71,81/)  !Mass and spin parameters
     nLowVar = 0
     nHighVar = 0
     nUsedPar = 0
     do p=1,4
        if(revID(IDs(p)).ne.0) then
           nLowVar  = nLowVar  + lowVar(revID(IDs(p)))
           nHighVar = nHighVar + highVar(revID(IDs(p)))
           nUsedPar = nUsedPar + 1
        end if
     end do
     meanRelVar = meanRelVar**(1.d0/dble(nmeanRelVar)) !Take geometric mean of (the variance of each chain, relative to the total variance)
     
     !Print and flag mean variance and variances for each chain:
     if(prConv.ge.3) then
        ch = ' '
        if(nLowVar.eq.nUsedPar) ch = '*'
        if(nHighVar.eq.nUsedPar) ch = '#'
        write(6,'(I12,A2,$)')ic,':'//ch
        do p=1,nMCMCpar
           if(fixedpar(p).eq.1) cycle
           ch = ' '
           if(lowVar(p).eq.1) ch = '*'
           if(highVar(p).eq.1) ch = '#'
           write(6,'(F8.4,A1,$)')chVar1(ic,p),ch
        end do
        ch = ' '
        if(meanRelVar.lt.0.5) ch = '*'
        if(meanRelVar.gt.2.0) ch = '#'
        write(6,'(F8.3,A1,$)')meanRelVar,ch
        write(6,*)''
     end if !if(prConv.ge.3)
  end do
  
  !Print mean variance for all parameters :
  if(prConv.ge.3) then
     write(6,'(A13,$)')'   Mean:'
     do p=1,nMCMCpar
        if(fixedpar(p).eq.1) cycle
        write(6,'(F9.4,$)')chVar(p)
     end do
     write(6,*)''
     write(6,*)''
  end if !if(prConv.ge.3)
  
  !Print the variances within chains and between chains:
  if(prConv.ge.2) then
     write(6,'(A)')'  Variances:'
     write(6,'(A14,$)')'      In chs: '
     do p=1,nMCMCpar
        if(fixedpar(p).eq.1) cycle
        write(6,'(ES9.1,$)')chVar(p)
     end do
     write(6,*)
     write(6,'(A14,$)')'   Betw. chs: '
     do p=1,nMCMCpar
        if(fixedpar(p).eq.1) cycle
        write(6,'(ES9.1,$)')meanVar(p)
     end do
     write(6,*)
  end if
  
  !Print R-hat:
  if(prConv.ge.1) then
     write(6,'(A14,$)')'       R-hat: '
     totRhat = 1.d0
     nRhat = 0
     do p=1,nMCMCpar
        if(fixedpar(p).eq.1) cycle
        write(6,'(F9.4,$)')Rhat(p)
        totRhat = totRhat * Rhat(p)
        nRhat = nRhat + 1
     end do
     !write(6,'(F9.4)')totRhat/dble(nRhat)
     write(6,'(F9.4)')totRhat**(1.d0/dble(nRhat))
  end if
  
end subroutine compute_convergence
!***********************************************************************************************************************************


