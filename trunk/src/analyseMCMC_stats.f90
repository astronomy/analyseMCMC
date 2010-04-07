! Compute statistics for analysemcmc

subroutine statistics(exitcode)
  use constants
  use analysemcmc_settings
  use general_data
  use stats_data
  use chain_data
  use mcmcrun_data
  implicit none
  integer :: c,i,ic,p,p1,p2,nr,nstat,exitcode,wraptype
  integer :: indexx(maxMCMCpar,maxChs*maxIter),index1(maxChs*maxIter)
  real :: revper
  real :: x1,x2,y1,y2
  real :: range1,minrange,ival,wrapival,centre,maxlogl,minlogl,shift,shIval,shIval2
  real :: medians(maxMCMCpar),mean(maxMCMCpar),var1(maxMCMCpar),var2(maxMCMCpar),corr,corr1,corr2
  
  !Need extra accuracy to compute Bayes factor
  !real*8 :: var,total,total2,total3,deltab
  real(kind=selected_real_kind(18,4931)) :: var,total,total2,total3,deltab  !Precision of 18, exponential range of 4931 - max gfortran supports?
  
  exitcode = 0
  
  
  !Convert MCMC parameters/PDFs (cos/sin->ang, rad->deg, etc):
  if(changeVar.ge.1) then
     do ic=1,nChains0  !selDat consists of nChains chains, allDat of nChains0 chains;  nChains <= nChains0
        !if(prProgress.ge.1.and.ic.eq.1.and.update.eq.0) write(stdOut,'(A)',advance="no")'.  Change vars. '
        do p=1,nMCMCpar
           select case(parID(p))
              
           case(21) !Take cube root: d^3 -> Distance:
              allDat(ic,p,1:Ntot(ic)) = allDat(ic,p,1:Ntot(ic))**c3rd
              if(ic.le.nChains) selDat(ic,p,1:n(ic)) = selDat(ic,p,1:n(ic))**c3rd
              if(ic.eq.1) startval(1:nChains0,p,1:3) = startval(1:nChains0,p,1:3)**c3rd
              
           case(22) !Take exp: logD -> Distance:
              allDat(ic,p,1:Ntot(ic)) = exp(allDat(ic,p,1:Ntot(ic)))   
              if(ic.le.nChains) selDat(ic,p,1:n(ic)) = exp(selDat(ic,p,1:n(ic)))
              if(ic.eq.1) startval(1:nChains0,p,1:3) = exp(startval(1:nChains0,p,1:3))

           case(51,72,82) !cos -> deg:
              allDat(ic,p,1:Ntot(ic)) = acos(allDat(ic,p,1:Ntot(ic)))*r2d
              if(ic.le.nChains) selDat(ic,p,1:n(ic)) = acos(selDat(ic,p,1:n(ic)))*rr2d
              if(ic.eq.1) startval(1:nChains0,p,1:3) = acos(startval(1:nChains0,p,1:3))*rr2d
              
           case(31) !rad -> h:
              allDat(ic,p,1:Ntot(ic)) = allDat(ic,p,1:Ntot(ic))*r2h  !rad -> h
              if(ic.le.nChains) selDat(ic,p,1:n(ic)) = selDat(ic,p,1:n(ic))*rr2h
              if(ic.eq.1) startval(1:nChains0,p,1:3) = startval(1:nChains0,p,1:3)*rr2h
              
           case(32,53) !sin -> deg:
              allDat(ic,p,1:Ntot(ic)) = asin(allDat(ic,p,1:Ntot(ic)))*r2d
              if(ic.le.nChains) selDat(ic,p,1:n(ic)) = asin(selDat(ic,p,1:n(ic)))*rr2d
              if(ic.eq.1) startval(1:nChains0,p,1:3) = asin(startval(1:nChains0,p,1:3))*rr2d
              
           case(41,52,54,73,83) !rad -> deg:
              allDat(ic,p,1:Ntot(ic)) = allDat(ic,p,1:Ntot(ic))*r2d
              if(ic.le.nChains) selDat(ic,p,1:n(ic)) = selDat(ic,p,1:n(ic))*rr2d
              if(ic.eq.1) startval(1:nChains0,p,1:3) = startval(1:nChains0,p,1:3)*rr2d
           end select
           
        end do !p
        
     end do !ic
     
     !Change the parameter names:
     call set_derivedParameterNames()
     
  end if !if(changeVar.ge.1)
  
  
  
  
  !Sort all data and find the interval limits for the default probability interval for the wrapable parameters
  if(prProgress.ge.2) write(stdOut,*)''
  shift = 0.
  wrap = 0
  raShift = 0.
  do ic=1,nChains
     if(mergeChains.eq.0.and.contrChain(ic).eq.0) cycle
     !wrapival = ivals(Nival) !Use the largest range
     wrapival = 0.999 !Always use a very large range (?)
     indexx = 0
     if(prProgress.ge.2.and.mergeChains.eq.0) write(stdOut,'(A,I2.2,A)',advance="no")' Ch',ic,' '
     if(prProgress.ge.2.and.ic.eq.1.and.wrapData.ge.1) write(stdOut,'(A)',advance="no")'  Wrap data. '
     do p=1,nMCMCpar
        if(wrapData.eq.0 .or. &
             (parID(p).ne.31.and.parID(p).ne.41.and.parID(p).ne.52.and.parID(p).ne.54.and.parID(p).ne.73.and.parID(p).ne.83) ) then  !Not RA, phi_c, psi, phi_Jo, phi_1,2
           call rindexx(n(ic),selDat(ic,p,1:n(ic)),index1(1:n(ic)))
           indexx(p,1:n(ic)) = index1(1:n(ic))
           if(parID(p).eq.31) raCentre = rpi                      !Plot 0-24h when not wrapping -> centre = 12h = pi
           cycle !No wrapping necessary
        end if
        
        
        !Determine 'wraptype':
        wraptype = 1                                !0-2pi (e.g. phases)
        if(parID(p).eq.52) wraptype = 2             !0-pi (polarisation angle)
        if(changeVar.ge.1) then
           if(parID(p).eq.31) wraptype = 3          !0-24h (RA)
           wraptype = wraptype+10                   !11,12,13 iso 1,2
        end if
        
        
        !Determine shift interval from wraptype:
        select case(wraptype)
        case(1)
           shIval = rtpi        ! "Phase": 0-2pi
        case(2)
           shIval = rpi         ! Pol.angle: 0-pi
        case(11)
           shIval = 360.        ! "Phase": 0-360
        case(12)
           shIval = 180.        ! Pol.angle: 0-180
        case(13)
           shIval = 24.         ! RA: 0-24h
        end select
        shIval2 = shIval/2.     ! Shift interval/2
        shIvals(ic,p) = shIval
        
        
        
        !Make sure data are between 0 and 2pi or between 0 and pi to start with:
        do i=1,n(ic)
           selDat(ic,p,i) = revper(selDat(ic,p,i),shIval)          !Bring selDat(i) between 0 and shIval
        end do
        call rindexx(n(ic),selDat(ic,p,1:n(ic)),index1(1:n(ic)))
        indexx(p,1:n(ic)) = index1(1:n(ic))
        
        minrange = 1.e30
        do i=1,n(ic)
           x1 = selDat(ic,p,indexx(p,i))
           x2 = selDat(ic,p,indexx(p,mod(i+nint(n(ic)*wrapival)-1,n(ic))+1))
           range1 = mod(x2 - x1 + real(10*shIval),shIval)    !0-shIval
           
           if(range1.lt.minrange) then
              minrange = range1
              y1 = x1
              y2 = x2
              !write(stdOut,'(2I6,7F10.5)')i,mod(nint(i+n(ic)*wrapival),n(ic)),x1,x2,range1,minrange,y1,y2,(y1+y2)/2.
           end if
           !write(stdOut,'(2I6,7F10.5)')i,mod(nint(i+n(ic)*wrapival),n(ic)),x1,x2,range1,minrange,y1,y2,(y1+y2)/2.
        end do !i
        centre = (y1+y2)/2.
        
        
        if(y1.gt.y2) then
           wrap(ic,p) = wraptype
           centre = mod(centre + shIval2, shIval)   !Distribution peaks close to 0/shIval; shift centre by shIval/2
        end if
        
        
        
        !Wrap around anticentre:
        shift = 0.
        
        !For the general case of shIval (= pi or 2pi; 180, 360 or 24)
        if(wrap(ic,p).gt.0) shift = shIval - mod(centre + shIval2, shIval)
        shifts(ic,p) = shift
        if(parID(p).eq.31.and.ic.eq.1) raShift = shift                                     !Save RA shift to plot sky map
        
        !Do the actual wrapping:
        allDat(ic,p,1:Ntot(ic))  = mod(allDat(ic,p,1:Ntot(ic))  + shift, shIval) - shift   !Original data
        selDat(ic,p,1:n(ic))     = mod(selDat(ic,p,1:n(ic))     + shift, shIval) - shift
        startval(ic,p,1:3)       = mod(startval(ic,p,1:3)       + shift, shIval) - shift   !Injection, starting and Lmax values
        y1 = mod(y1 + shift, shIval) - shift
        y2 = mod(y2 + shift, shIval) - shift
        
        centre = mod(centre + shift, shIval) - shift
        if(parID(p).eq.31.and.ic.eq.1) raCentre = centre                                   !Save RA centre to plot sky map
        
        minrange = y2-y1
        call rindexx(n(ic),selDat(ic,p,1:n(ic)),index1(1:n(ic)))  !Re-sort
        indexx(p,1:n(ic)) = index1(1:n(ic))
        
        if(abs( abs( minval(selDat(ic,p,1:n(ic))) - maxval(selDat(ic,p,1:n(ic))) ) - shIval) .lt. 1.e-3)  wrap(ic,p) = 1   !If centre is around shIval2, still needs to be flagged 'wrap' to plot PDF
     end do !p
     ! End wrapping data
     
     
     
     
     
     
     
     !Do statistics:
     if(prProgress.ge.1.and.ic.eq.1) write(stdOut,'(A)',advance="no")'  Calc: stats, '
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
     
     
     
     
     
     !Compute correlations:
     if(prCorr.gt.0.or.saveStats.gt.0) then
        !write(stdOut,'(A)')' Calculating correlations...   '
        if(prProgress.ge.1) write(stdOut,'(A)',advance="no")' corrs, '
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
     
     
     
     
     !Determine interval ranges:
     if(prProgress.ge.1.and.ic.eq.1) write(stdOut,'(A)',advance="no")' prob.ivals: '
     c0 = 0
     do c=1,Nival
        ival = ivals(c)
        c0 = ival0
        if(c.ne.c0 .and. prIval.eq.0 .and. prStat.lt.2 .and. saveStats.eq.0) cycle
        
        if(prProgress.ge.1.and.ic.eq.1) write(stdOut,'(F6.3)',advance="no")ival
        do p=1,nMCMCpar
           minrange = 1.e30
           !write(stdOut,'(A8,4x,4F10.5,I4)')parNames(parID(p)),y1,y2,minrange,centre,wrap(ic,p)
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
              !write(stdOut,'(I6,7F10.5)')i,x1,x2,range1,minrange,y1,y2,(y1+y2)/2.
           end do
           centre = (y1+y2)/2.
           !write(stdOut,'(A8,4x,4F10.5,I4)')parNames(parID(p)),y1,y2,minrange,centre,wrap(ic,p)
           
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
     !if(prProgress.ge.2) write(stdOut,'(A34,F8.4)')'.  Standard probability interval: ',ivals(ival0)
     !if(prProgress.ge.2) write(stdOut,'(A,F8.4)',advance="no")', default ival:',ivals(ival0)
     
     
     
     
     
     !Compute Bayes factor:
     !if(prProgress.ge.1.and.ic.eq.1) write(stdOut,'(A)',advance="no")'  Bayes factor,'
     total = 0
     total2 = 0
     total3 = 0
     maxlogl = -1.e30
     minlogl =  1.e30
     do i = Nburn(ic),Ntot(ic)
        var = post(ic,i)          !Use highest precision available
        total = total + exp(-var)**(1/Tchain(ic))
        total2 = total2 + exp(var)**(1-(1/Tchain(ic)))
        total3 = total3 + var
        maxlogl = max(post(ic,i),maxlogl)
        minlogl = min(post(ic,i),minlogl)
     end do
     !var = dble(n(ic))/total
     var = total2/total
     logebayesfactor(ic) = real(log(var))
     log10bayesfactor(ic) = real(log10(var))
     
     if(ic.eq.nChains) then
        deltab = 1/Tchain(ic)
     else
        deltab = 1/Tchain(ic) - 1/Tchain(ic+1)
     end if
     
     logebayestempfactor(ic) = real((total3/(n(ic)))*deltab)
     !write(stdOut,'(A,4F10.3,2I9,F10.3)')'ln Bayes',logebayesfactor(ic),log10bayesfactor(ic),maxlogl,minlogl,n(ic),ic,Tchain(ic)
     write(stdOut,'(A,F10.3,I9,2F10.3,3I9)')'ln Bayes',logebayestempfactor(ic),ic,Tchain(ic),deltab,Nburn(ic),Ntot(ic),n(ic)
     
     
     
     
     
     
     
     !Change MCMC parameters
     !(moved to beginning of routine)
     
     
     
     
     
     
     !Find 100% probability range
     do c = 1,Nival
        if(abs(ivals(c)-1.).lt.1.e-4) then !Then treat it as a 100% interval to prevent numerical problems
           if(prProgress.ge.1) write(stdOut,'(A,F9.4,A)')'  Treating probability interval',ivals(c)*100,'% as 100%'
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
        if(ic.eq.nChains) then
           write(stdOut,*)
        else
           write(stdOut,'(A)',advance="no")'  '
        end if
     end if
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
     !**********************************************************************************************
     !******   PRINT STATISTICS   ******************************************************************
     !**********************************************************************************************
     
     
     if(prProgress+prStat+prIval+prConv.gt.0.and.ic.eq.1) then
        write(stdOut,'(/,A,2(A,F7.2))',advance="no")'  Bayes factor:   ','log_e(B_SN) =',logebayesfactor(ic),',  log_10(B_SN) =',log10bayesfactor(ic)
        !write(stdOut,'(8x,A,3(A,F7.2),A)')'  Maximum likelihood:   ','log_e(Lmax) =',startval(ic,1,3),',  log_10(Lmax) =',startval(ic,1,3)/log(10.),',  -> SNR =',sqrt(2*startval(ic,1,3)),'.'
        write(stdOut,'(8x,A,3(A,F7.2),A)')'  Maximum likelihood:   ','log_e(Lmax) =',loglmax,',  log_10(Lmax) =',loglmax/log(10.),',  -> SNR =',sqrt(2*loglmax),'.'
     end if
     
     
     !Print statistics to screen
     if(prStat.gt.0) then
        write(stdOut,'(/,A)')'  Main statistics:'
        do c=1,Nival
           if(c.ne.c0.and.prStat.lt.2) cycle
           if(c.gt.1.and.prStat.ge.2) write(stdOut,*)
           write(stdOut,'(A10, A12,2A10,A12, 4A8, 4A10,A8,A10, A4,A12,F7.3,A2)')'Param.  ','model','median','mean','Lmax','stdev1','stdev2','abvar1','abvar2',  &
                'rng_c','rng1','rng2','drng','d/drng','delta','ok?','result (',ivals(c)*100,'%)'
           do p=1,nMCMCpar
              if(fixedpar(p).eq.1) cycle  !Varying parameters only
              write(stdOut,'(A10,F12.6,2F10.4,F12.6, 4F8.4,4F10.4,F8.4,F10.4)',advance="no")parNames(parID(p)),startval(ic,p,1),stats(ic,p,1),stats(ic,p,2),startval(ic,p,3),stdev1(p),stdev2(p),absVar1(p),  &
                   absVar2(p),ranges(ic,c,p,3),ranges(ic,c,p,1),ranges(ic,c,p,2),ranges(ic,c,p,4),  &
                   2*abs(startval(ic,p,1)-ranges(ic,c,p,3))/ranges(ic,c,p,4),ranges(ic,c,p,5)  !d/drange wrt centre of range
              if(startval(ic,p,1).ge.ranges(ic,c,p,1).and.startval(ic,p,1).le.ranges(ic,c,p,2)) then
                 write(stdOut,'(A4)',advance="no")'y '
              else
                 write(stdOut,'(A4)',advance="no")'*N*'
              end if
              write(stdOut,'(F10.4,A3,F9.4)')ranges(ic,c,p,3),'+-',0.5*ranges(ic,c,p,4)
           end do !p
        end do !c
     end if
     
     
     !Print intervals as: centre, delta, in range:
     if(prIval.eq.1.or.prIval.eq.3) then
        write(stdOut,'(/,A)')'  Probability intervals:'
        write(stdOut,'(A22,A8)',advance="no")'Interval:',''
        do c=1,Nival
           write(stdOut,'(F20.4,A9)',advance="no")ivals(c),''
        end do
        write(stdOut,*)''
        
        write(stdOut,'(A10,2x,2A9)',advance="no")'Param.  ','model','median'
        do c=1,Nival
           !write(stdOut,'(2x,2A9,A8)',advance="no")'rng1','rng2','in rnge'
           write(stdOut,'(2x,3A9)',advance="no")'centre','delta','in rnge'
        end do
        write(stdOut,*)''
        do p=1,nMCMCpar
           !if(stdev1(p).lt.1.d-20) cycle
           if(fixedpar(p).eq.1) cycle  !Varying parameters only
           write(stdOut,'(A10,2x,2F9.4)',advance="no")parNames(parID(p)),startval(ic,p,1),stats(ic,p,1)
           do c=1,Nival
              if(mergeChains.eq.0) then
                 write(stdOut,'(2x,2F9.4,F6.3)',advance="no")ranges(ic,c,p,3),ranges(ic,c,p,4),min(2*abs(startval(ic,p,1)-ranges(ic,c,p,3))/ranges(ic,c,p,4),9.999) !Defined with centre of prob. range, need some extra security to print correctly
              else
                 write(stdOut,'(2x,2F9.4,F6.3)',advance="no")ranges(ic,c,p,3),ranges(ic,c,p,4),2*abs(startval(ic,p,1)-ranges(ic,c,p,3))/ranges(ic,c,p,4) !Defined with centre of prob. range
              end if
              if(startval(ic,p,1).gt.ranges(ic,c,p,1).and.startval(ic,p,1).lt.ranges(ic,c,p,2)) then
                 write(stdOut,'(A3)',advance="no")'y '
              else
                 write(stdOut,'(A3)',advance="no")'N*'
              end if
           end do !c
           write(stdOut,*)
        end do !p
        write(stdOut,*)
     end if
     
     
     if(prIval.ge.2) then
        write(stdOut,'(/,A)')'  Statistics and probability intervals:'
        
        !Print intervals as x +- dx:
        write(stdOut,'(A61)',advance="no")'Interval:'
        do c=1,Nival
           write(stdOut,'(F14.3,A1,A9)',advance="no")ivals(c)*100,'%',''
        end do
        write(stdOut,*)''
        
        write(stdOut,'(A11,1x,4A10)',advance="no")'Parameter','median','mean','Lmax','stdev'
        do c=1,Nival
           write(stdOut,'(5x,A8,4x,A7)',advance="no")'x','dx'
        end do
        write(stdOut,*)''
        do p=1,nMCMCpar
           if(fixedpar(p).eq.1) cycle  !Varying parameters only
           write(stdOut,'(A10,2x,4F10.3)',advance="no")parNames(parID(p)),stats(ic,p,1),stats(ic,p,2),startval(ic,p,3),stdev2(p)
           do c=1,Nival
              write(stdOut,'(5x,F8.3,A4,F7.3)',advance="no")ranges(ic,c,p,3),' +- ',0.5d0*ranges(ic,c,p,4)
           end do
           write(stdOut,*)''
        end do
        write(stdOut,*)''
        
        !Print intervals as low-high:
        write(stdOut,'(A61)',advance="no")'Interval:'
        do c=1,Nival
           write(stdOut,'(F14.3,A1,A9)',advance="no")ivals(c)*100,'%',''
        end do
        write(stdOut,*)''
        
        write(stdOut,'(A11,1x,4A10)',advance="no")'Parameter','median','mean','Lmax','stdev'
        do c=1,Nival
           write(stdOut,'(6x,A8,3x,A7)',advance="no")'min','max'
        end do
        write(stdOut,*)''
        do p=1,nMCMCpar
           if(fixedpar(p).eq.1) cycle  !Varying parameters only
           write(stdOut,'(A10,2x,4F10.3)',advance="no")parNames(parID(p)),stats(ic,p,1),stats(ic,p,2),startval(ic,p,3),stdev2(p)
           do c=1,Nival
              write(stdOut,'(6x,F8.3,A3,F7.3)',advance="no")ranges(ic,c,p,1),' - ',ranges(ic,c,p,2)
           end do
           write(stdOut,*)''
        end do
        write(stdOut,*)''
        
        
        !Print intervals as x -dx1 +dx2 in LaTeX format:
        if(1.eq.2) then
           write(stdOut,'(A19)',advance="no")'Interval:   '
           do c=1,Nival
              write(stdOut,'(F27.3,A1,A9)',advance="no")ivals(c)*100,'%',''
           end do
           write(stdOut,*)''
           
           do p=1,nMCMCpar
              if(fixedpar(p).eq.1) cycle  !Varying parameters only
              write(stdOut,'(A10,2x)',advance="no")parNames(parID(p))
              do c=1,Nival
                 write(stdOut,'(5x,A1,F8.3,2(A,F7.3),A)',advance="no")'$',stats(ic,p,1),'_{-',abs(stats(ic,p,1)-ranges(ic,c,p,1)),'}^{+',abs(stats(ic,p,1)-ranges(ic,c,p,2)),'}$'
              end do
              write(stdOut,*)''
           end do
           write(stdOut,*)''
        end if
        
     end if
     
     
     
     
     
     
     
     
     
     
     !Print correlations:
     if(prCorr.gt.0) then
        corr1 = 0.1
        corr2 = 0.5
        write(stdOut,'(/,A)',advance="no")'  Correlations  '
        write(stdOut,'(A,3(F4.2,A))')'  (weak [',corr1,'<abs(cor)<',corr2,']: in lower triangle,  strong [abs(cor)>',corr2,']: in upper triangle):'
        write(stdOut,'(A8)',advance="no")''
        do p=1,nMCMCpar
           if(fixedpar(p).eq.0) write(stdOut,'(A7)',advance="no")trim(parNames(parID(p)))
        end do
        write(stdOut,*)''
        do p1=1,nMCMCpar
           if(fixedpar(p1).eq.1) cycle  !Varying parameters only
           write(stdOut,'(A8)',advance="no")trim(parNames(parID(p1)))
           do p2=1,nMCMCpar
              corr = corrs(p1,p2)
              if(fixedpar(p2).eq.1) cycle  !Varying parameters only
              if( (abs(corr).ge.corr2.and.p2.gt.p1) ) then  !Print in the upper triangle
                 write(stdOut,'(F7.2)',advance="no")corr
              else if( (abs(corr).ge.corr1.and.abs(corr).lt.corr2.and.p1.gt.p2) ) then   !Print in the lower triangle
                 write(stdOut,'(F7.2)',advance="no")corr
              else if(p1.eq.p2) then !Print on the diagonal
                 write(stdOut,'(A7)',advance="no")' ######' !'
              else
                 write(stdOut,'(A7)',advance="no")''
              end if
           end do
           write(stdOut,'(A)')'   '//trim(parNames(parID(p1)))
        end do
     end if
     
     
     
     !Print output for CBC Wiki:
     if(ic.eq.1.and.wikioutput.eq.1) call save_cbc_wiki_data(ic)
     
     
  end do !ic
  
  
  
  !Average the Bayes factors from all chains:
  logebayesfactortotalgeom = 0.0
  logebayesfactortotalarith = 0.0
  logebayesfactortotalharmo = 0.0
  logebayesfactortotal = 0.0
  var = 0.0
  total = 0.0
  total2 = 0.0
  do ic = 1,nchains0
     var = logebayesfactor(ic)
     !write(6,'(A,F10.3)')'ln Bayes',logebayesfactor(ic)
     logebayesfactortotalgeom = logebayesfactortotalgeom+var
     !write(6,'(A,F10.3)')'logebayesfactortotalgeom',logebayesfactortotalgeom
     total = total+exp(var)
     !write(6,'(A,F10.3)')'logebayesfactortotalarith',logebayesfactortotalarith
     total2 = total2+(exp(-var))
     !write(6,'(A,F10.3)')'logebayesfactortotalharmo',logebayesfactortotalharmo
     logebayesfactortotal = logebayesfactortotal+logebayestempfactor(ic)
  end do
  logebayesfactortotalgeom = logebayesfactortotalgeom/dble(nchains)
  total = total/dble(nchains)
  logebayesfactortotalarith = real(log(total))
  total2 = dble(nchains)/total2
  logebayesfactortotalharmo = real(log(total2))
  
  !write(6,'(A,F10.3)')'ln Bayes Total',logebayesfactortotal
  
  
  !!Compute autocorrelations:
  if(prAcorr.gt.0.or.plAcorr.gt.0) call compute_autocorrelations()
  
  !Compute and print convergence:
  if(nChains0.gt.1 .and. (prConv.ge.1.or.saveStats.ge.1)) call compute_convergence()  !Need unwrapped data for this (?)
  
  
  !Change the original chain data:
  !moved to top of the routine
  
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
  integer :: c,i,ic,o,p,p1,p2,exitcode,status,system
  
  exitcode = 0
  ic = 1 !Use chain 1
  o = 20 !Output port
  open(unit=o, form='formatted', status='replace',file=trim(outputdir)//'/'//trim(outputname)//'__statistics.dat')
  write(o,'(A)')trim(outputname)
  
  !Print general run and detector info:
  write(o,'(//,A,/)')'GENERAL INFORMATION:'
  write(o,'(6x,4A12,A12,A5, A8,A22,A8)')'totiter','totlines','totpts','totburn','nChains','used','seed','<d|d>','ndet'
  write(o,'(6x,4I12,I12,I5, I8,F22.10,I8)')totiter,totlines,totpts,totlines-totpts,nChains0,contrChains,seed(ic),DoverD,ndet(ic)
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
  write(o,'(A9)',advance="no")''
  do p=1,nMCMCpar
     write(o,'(A10)',advance="no")trim(parNames(parID(p)))
  end do
  write(o,*)''
  do p1=1,nMCMCpar
     write(o,'(A9)',advance="no")trim(parNames(parID(p1)))
     do p2=1,nMCMCpar
        write(o,'(F10.5)',advance="no")corrs(p1,p2)
     end do
     write(o,'(A)')'   '//trim(parNames(parID(p1)))
  end do
  
  
  !Print probability intervals:
  write(o,'(///,A,/)')'1D PROBABILITY INTERVALS:'
  write(o,'(A,I3)')'Nival:',Nival
  write(o,'(A22)',advance="no")'Interval:'
  do c=1,Nival
     write(o,'(F21.5,A14)',advance="no")ivals(c),''
  end do
  write(o,*)''
  
  write(o,'(A8,2x)',advance="no")'param.'
  do c=1,Nival
     !write(o,'(2x,2A9,A8)',advance="no")'rng1','rng2','in rnge'
     write(o,'(2x,2A12,A9)',advance="no")'centre','delta','in rnge'
  end do
  write(o,*)''
  do p=1,nMCMCpar
     write(o,'(A8,2x)',advance="no")trim(parNames(parID(p)))
     do c=1,Nival
        !write(o,'(2x,2F11.6,F6.3)',advance="no")ranges(ic,c,p,1),ranges(ic,c,p,2),2*abs(startval(ic,p,1)-ranges(ic,c,p,3))/ranges(ic,c,p,4) !Defined with centre of prob. range
        if(fixedpar(p).eq.0) then
           write(o,'(2x,2F12.6,F7.3)',advance="no")ranges(ic,c,p,3),ranges(ic,c,p,4),2*abs(startval(ic,p,1)-ranges(ic,c,p,3))/ranges(ic,c,p,4) !Defined with centre of prob. range
        else
           write(o,'(2x,2F12.6,F7.3)',advance="no")0.,0.,99.999
        end if
        if(startval(ic,p,1).gt.ranges(ic,c,p,1).and.startval(ic,p,1).lt.ranges(ic,c,p,2)) then
           write(o,'(A2)',advance="no")'y'
        else
           write(o,'(A2)',advance="no")'N'
        end if
     end do
     write(o,*)''
  end do
  
  
  
  !Print 2D intervals
  write(o,'(///,A,/)')'2D PROBABILITY INTERVALS:'
  write(o,'(A,I5)')'Npdf2D: ',Npdf2D
  write(o,'(A,2I5)')'Nbin2Dx,Nbin2Dy: ',Nbin2Dx,Nbin2Dy
  !write(o,*)''
  write(o,'(A28)',advance="no")'Interval:'
  do c=1,Nival
     write(o,'(F19.5)',advance="no")ivals(c)
  end do
  write(o,*)''
  
  write(o,'(A9,A19)',advance="no")'params.',''
  do c=1,Nival
     write(o,'(A16,A3)',advance="no")'delta','in'
  end do
  write(o,*)''
  do p=1,Npdf2D
     p1 = revID(PDF2Dpairs(p,1))
     p2 = revID(PDF2Dpairs(p,2))
     if(p1*p2.eq.0) cycle
     if(parID(p1)*parID(p2).eq.0) then
        if(parID(p1).eq.0) write(stdErr,'(/,A,/,A,//)')'  ***  ERROR:  save_stats():  parameter '//trim(parNames(parID(p1)))//' not defined, check PDF2Dpairs in the input file ***','  Aborting...'
        if(parID(p2).eq.0) write(stdErr,'(/,A,/,A,//)')'  ***  ERROR:  save_stats():  parameter '//trim(parNames(parID(p2)))//' not defined, check PDF2Dpairs in the input file ***','  Aborting...'
        stop
     end if
     write(o,'(2I4,2(2x,A8),2x)',advance="no")p1,p2,trim(parNames(parID(p1))),trim(parNames(parID(p2)))
     do c=1,Nival
        write(o,'(2x,F14.8)',advance="no")probareas(p1,p2,c,1)
        if(c.ge.Nival+1-injectionranges2d(p1,p2) .and. injectionranges2d(p1,p2).ne.0) then
           write(o,'(A3)',advance="no")'y'
        else
           write(o,'(A3)',advance="no")'n'
        end if
     end do
     write(o,*)''
  end do
  
  
  
  close(o) !Statistics output file
  if(saveStats.eq.2) status = system('a2ps -1rf7 '//trim(outputdir)//'/'//trim(outputname)//'__statistics.dat -o '//trim(outputdir)//'/'//trim(outputname)//'__statistics.ps')
  !write(stdOut,*)''
  if(prProgress.ge.1) then
     if(saveStats.eq.1) write(stdOut,'(A)')'  Statistics saved in '//trim(outputdir)//'/'//trim(outputname)//'__statistics.dat'
     if(saveStats.eq.2) write(stdOut,'(A)')'  Statistics saved in '//trim(outputdir)//'/'//trim(outputname)//'__statistics.dat,ps'
  end if
  
end subroutine save_stats
!***********************************************************************************************************************************


!***********************************************************************************************************************************
subroutine save_bayes(exitcode)  !Save Bayes-factor statistics to file  
  use constants
  use analysemcmc_settings
  use general_data
  use mcmcrun_data
  use chain_data
  
  implicit none
  integer :: ic,o,exitcode,status,system
  
  exitcode = 0
  !ic = 1 !Use chain 1
  o = 20 !Output port
  open(unit=o, form='formatted', status='replace',file=trim(outputdir)//'/'//trim(outputname)//'__bayes.dat')
  write(o,'(A)')trim(outputname)
  
  write(o,'(//,A,/)')'GENERAL INFORMATION:'
  write(o,'(6x,7A12,5A22)')'nChains','used','totiter','totlines','totpts','totburn','seed','<d|d>','ln(Bayes_total_arith)','ln(Bayes_total_harmo)','ln(Bayes_total_geom)','ln(Bayes_total_temp)'
  write(o,'(6x,7I12,F22.5,4F22.5)')nchains0,contrChains,totiter,totlines,totpts,totlines-totpts,seed(1),DoverD,logebayesfactortotalarith,logebayesfactortotalharmo,logebayesfactortotalgeom,logebayesfactortotal

  
  write(o,'(//,A,/)')'EVIDENCES:'
  write(o,'(6x,5A12,2A22)')'chain','totiter','totlines','totpts','totburn','temperature','ln(Bayes)'
  do ic=1,nChains0
     write(o,'(6x,5I12,F22.3,F22.5)')ic,totiter,totlines,totpts,totlines-totpts,Tchain(ic),logebayesfactor(ic)
  end do
  
  close(o) !Bayes output file
  
  if(saveStats.eq.2) status = system('a2ps -1rf7 '//trim(outputdir)//'/'//trim(outputname)//'__bayes.dat -o '//trim(outputdir)//'/'//trim(outputname)//'__bayes.ps')
  !write(stdOut,*)''
  if(prProgress.ge.1) then
     if(saveStats.eq.1) write(stdOut,'(A)')'  Bayes factors saved in '//trim(outputdir)//'/'//trim(outputname)//'__bayes.dat'
     if(saveStats.eq.2) write(stdOut,'(A)')'  Bayes factors saved in '//trim(outputdir)//'/'//trim(outputname)//'__bayes.dat,ps'
  end if
  
end subroutine save_bayes
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
  o = 20
  write(wikifilename,'(A)')trim(outputname)//'__wiki.dat'
  open(unit=o,form='formatted',status='replace',action='write',position='rewind',file=trim(outputdir)//'/'//trim(wikifilename),iostat=io)
  if(io.ne.0) then
     write(stdErr,'(A)')'  Error opening '//trim(wikifilename)//', aborting...'
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
  
  write(o,'(A4)',advance="no")'||  '
  do i=1,4
     if(i.le.ndet(ic)) then
        write(o,'(A2)',advance="no")detabbrs(detnr(ic,i))
     else
        write(o,'(A2)',advance="no")'  '
     end if
  end do
  write(o,'(A6)',advance="no")'    '
  parr(1:14) = (/63,64,61,62,11,71,72,81,82,22,31,32,51,52/)
  do p=1,14
     p1 = parr(p)
     if(revID(p1).eq.0) then  !Parameter not used
        write(o,'(A5,A11)',advance="no")'  || ',' - '
     else
        x = allDat(ic,revID(p1),1)
        if(p1.eq.31) x = rev2pi(x*rh2r)
        if(p1.eq.52) x = rrevpi(x*rd2r)  !Polarisation angle
        if(p1.eq.41.or.p1.eq.54.or.p1.eq.73.or.p1.eq.83) x = rev2pi(x*rd2r)
        if(p1.eq.32.or.p1.eq.51.or.p1.eq.53.or.p1.eq.72.or.p1.eq.82) x = x*rd2r
        if(p1.ge.11.and.p1.le.19) then
           write(o,'(A5,F14.4)',advance="no")'  || ',x+t0
        else
           write(o,'(A5,F11.4)',advance="no")'  || ',x
        end if
     end if
  end do
  write(o,'(A,79x,A)')'  || [ injection info] ',' ||'
  
  
  
  
  !Bayes factor:
  write(o,'(///,A)')'== Bayes Factors =='
  write(o,'(A)')"|| '''Code'''                                     || '''Model'''                                                || '''Detectors'''  || '''log_e Bayes Factor'''    || '''log_10 Bayes Factor'''    || '''Details'''                                                           ||"
  
  write(o,'(A3,A47)',advance="no")'|| ','[http://tinyurl.com/SPINspiral SPINspiral]     '
  if(spinningRun.eq.0) write(o,'(A3,A59)',advance="no")'|| ',trim(waveforms(waveform))//' '//pnstr//'pN non-spinning vs. Gaussian noise '
  if(spinningRun.eq.1) write(o,'(A3,A59)',advance="no")'|| ',trim(waveforms(waveform))//' '//pnstr//'pN one spin vs. Gaussian noise     '
  if(spinningRun.eq.2) write(o,'(A3,A59)',advance="no")'|| ',trim(waveforms(waveform))//' '//pnstr//'pN two spins vs. Gaussian noise    '
  write(o,'(A)',advance="no")'||  '
  do i=1,4
     if(i.le.ndet(ic)) then
        write(o,'(A2)',advance="no")detabbrs(detnr(ic,i))
     else
        write(o,'(A2)',advance="no")'  '
     end if
  end do
  write(o,'(A11)',advance="no")'        || '
  write(o,'(F10.1,A)',advance="no")logebayesfactor(ic),'                  || '
  write(o,'(F10.1,A)',advance="no")log10bayesfactor(ic),'                   || '
  write(o,'(A)')'['//trim(url)//' link]         ||'
  
  
  
  
  !Medians:
  write(o,'(///,A)')'== Medians =='
  write(o,'(A)')"|| '''Code'''                                     || '''Waveform'''                                 || '''Detectors'''  || '''Mc'''     || '''&eta;'''  || '''time'''   || '''spin1'''  ||'''&theta;1'''|| '''spin2'''  ||'''&theta;2'''|| '''Dist'''   || '''R.A.'''   || '''Dec.'''   || '''incl'''   || '''pol.'''   || '''details'''                                                                                     ||"
  write(o,'(A)')"||                                                ||                                                ||                  || (Mo)         ||              ||  (s)         ||              || (rad)        ||              || (rad)        || (Mpc)        || (rad)        || (rad)        || (rad)        || (rad)        ||                                                                                                   ||"
  write(o,'(A3,A47)',advance="no")'|| ','[http://tinyurl.com/SPINspiral SPINspiral]     '
  if(spinningRun.eq.0) write(o,'(A3,A47)',advance="no")'|| ',trim(waveforms(waveform))//' '//pnstr//'pN non-spinning '
  if(spinningRun.eq.1) write(o,'(A3,A47)',advance="no")'|| ',trim(waveforms(waveform))//' '//pnstr//'pN one spin     '
  if(spinningRun.eq.2) write(o,'(A3,A47)',advance="no")'|| ',trim(waveforms(waveform))//' '//pnstr//'pN two spins    '
  write(o,'(A)',advance="no")'||  '
  do i=1,4
     if(i.le.ndet(ic)) then
        write(o,'(A2)',advance="no")detabbrs(detnr(ic,i))
     else
        write(o,'(A2)',advance="no")'  '
     end if
  end do
  write(o,'(A10)',advance="no")'       ||'
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
     write(o,'(A11,A5)',advance="no")xs11,'   ||'
  end do
  write(o,'(A)')' ['//trim(url)//' link]                                   ||'
  
  
  
  
  !Means:
  write(o,'(///,A)')'== Means =='
  write(o,'(A)')"|| '''Code'''                                     || '''Waveform'''                                 || '''Detectors'''  || '''Mc'''     || '''&eta;'''  || '''time'''   || '''spin1'''  ||'''&theta;1'''|| '''spin2'''  ||'''&theta;2'''|| '''Dist'''   || '''R.A.'''   || '''Dec.'''   || '''incl'''   || '''pol.'''   || '''details'''                                                                                     ||"
  write(o,'(A)')"||                                                ||                                                ||                  || (Mo)         ||              ||  (s)         ||              || (rad)        ||              || (rad)        || (Mpc)        || (rad)        || (rad)        || (rad)        || (rad)        ||                                                                                                   ||"
  write(o,'(A3,A47)',advance="no")'|| ','[http://tinyurl.com/SPINspiral SPINspiral]     '
  if(spinningRun.eq.0) write(o,'(A3,A47)',advance="no")'|| ',trim(waveforms(waveform))//' '//pnstr//'pN non-spinning '
  if(spinningRun.eq.1) write(o,'(A3,A47)',advance="no")'|| ',trim(waveforms(waveform))//' '//pnstr//'pN one spin     '
  if(spinningRun.eq.2) write(o,'(A3,A47)',advance="no")'|| ',trim(waveforms(waveform))//' '//pnstr//'pN two spins    '
  write(o,'(A)',advance="no")'||  '
  do i=1,4
     if(i.le.ndet(ic)) then
        write(o,'(A2)',advance="no")detabbrs(detnr(ic,i))
     else
        write(o,'(A2)',advance="no")'  '
     end if
  end do
  write(o,'(A10)',advance="no")'       ||'
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
     write(o,'(A11,A5)',advance="no")xs11,'   ||'
  end do
  write(o,'(A)')' ['//trim(url)//' link]                                   ||'
  
  
  
  
  !Lmax:
  write(o,'(///,A)')'== Maximum-likelihood points =='
  write(o,'(A)')"|| '''Code'''                                     || '''Waveform'''                                 || '''Detectors'''  ||'''log(L)'''  || '''Mc'''     || '''&eta;'''  || '''time'''   || '''spin1'''  ||'''&theta;1'''|| '''spin2'''  ||'''&theta;2'''|| '''Dist'''   || '''R.A.'''   || '''Dec.'''   || '''incl'''   || '''pol.'''   || '''details'''                                                                                     ||"
  write(o,'(A)')"||                                                ||                                                ||                  ||              || (Mo)         ||              ||  (s)         ||              || (rad)        ||              || (rad)        || (Mpc)        || (rad)        || (rad)        || (rad)        || (rad)        ||                                                                                                   ||"
  write(o,'(A3,A47)',advance="no")'|| ','[http://tinyurl.com/SPINspiral SPINspiral]     '
  if(spinningRun.eq.0) write(o,'(A3,A47)',advance="no")'|| ',trim(waveforms(waveform))//' '//pnstr//'pN non-spinning '
  if(spinningRun.eq.1) write(o,'(A3,A47)',advance="no")'|| ',trim(waveforms(waveform))//' '//pnstr//'pN one spin     '
  if(spinningRun.eq.2) write(o,'(A3,A47)',advance="no")'|| ',trim(waveforms(waveform))//' '//pnstr//'pN two spins    '
  write(o,'(A)',advance="no")'||  '
  do i=1,4
     if(i.le.ndet(ic)) then
        write(o,'(A2)',advance="no")detabbrs(detnr(ic,i))
     else
        write(o,'(A2)',advance="no")'  '
     end if
  end do
  write(o,'(A10)',advance="no")'       ||'
  
  !Print log L:
  write(xs11,'(F11.4)')logLmax
  write(o,'(A11,A5)',advance="no")xs11,'   ||'
  
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
     write(o,'(A11,A5)',advance="no")xs11,'   ||'
  end do
  write(o,'(A)')' ['//trim(url)//' link]                                   ||'
  
  
  
  
  !Stdev:
  write(o,'(///,A)')'== Standard deviations =='
  write(o,'(A)')"|| '''Code'''                                     || '''Waveform'''                                 || '''Detectors'''  || '''Mc'''     || '''&eta;'''  || '''time'''   || '''spin1'''  ||'''&theta;1'''|| '''spin2'''  ||'''&theta;2'''|| '''Dist'''   || '''R.A.'''   || '''Dec.'''   || '''incl'''   || '''pol.'''   || '''details'''                                                                                     ||"
  write(o,'(A)')"||                                                ||                                                ||                  || (Mo)         ||              ||  (s)         ||              || (rad)        ||              || (rad)        || (Mpc)        || (rad)        || (rad)        || (rad)        || (rad)        ||                                                                                                   ||"
  write(o,'(A3,A47)',advance="no")'|| ','[http://tinyurl.com/SPINspiral SPINspiral]     '
  if(spinningRun.eq.0) write(o,'(A3,A47)',advance="no")'|| ',trim(waveforms(waveform))//' '//pnstr//'pN non-spinning '
  if(spinningRun.eq.1) write(o,'(A3,A47)',advance="no")'|| ',trim(waveforms(waveform))//' '//pnstr//'pN one spin     '
  if(spinningRun.eq.2) write(o,'(A3,A47)',advance="no")'|| ',trim(waveforms(waveform))//' '//pnstr//'pN two spins    '
  write(o,'(A)',advance="no")'||  '
  do i=1,4
     if(i.le.ndet(ic)) then
        write(o,'(A2)',advance="no")detabbrs(detnr(ic,i))
     else
        write(o,'(A2)',advance="no")'  '
     end if
  end do
  write(o,'(A10)',advance="no")'       ||'
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
     write(o,'(A11,A5)',advance="no")xs11,'   ||'
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
     write(stdErr,'(A)')'  Error: 2-sigma range not found, needed for Wiki output!'
     write(stdOut,'(A)',advance="no")'  Do you want to continue?  (y/n)  '
     read(5,*)ans
     if(ans.eq.'y'.or.ans.eq.'Y') then
        c = 1
        write(stdOut,'(A,F6.2,A)')'  Continuing with ',ivals(c)*100,"% probability interval, don't use wiki.txt!!!"
     else
        stop
     end if
  end if
  write(o,'(A3,A47)',advance="no")'|| ','[http://tinyurl.com/SPINspiral SPINspiral]     '
  if(spinningRun.eq.0) write(o,'(A3,A47)',advance="no")'|| ',trim(waveforms(waveform))//' '//pnstr//'pN non-spinning '
  if(spinningRun.eq.1) write(o,'(A3,A47)',advance="no")'|| ',trim(waveforms(waveform))//' '//pnstr//'pN one spin     '
  if(spinningRun.eq.2) write(o,'(A3,A47)',advance="no")'|| ',trim(waveforms(waveform))//' '//pnstr//'pN two spins    '
  write(o,'(A)',advance="no")'||  '
  do i=1,4
     if(i.le.ndet(ic)) then
        write(o,'(A2)',advance="no")detabbrs(detnr(ic,i))
     else
        write(o,'(A2)',advance="no")'  '
     end if
  end do
  write(o,'(A10)',advance="no")'       ||'
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
     write(o,'(A20,A5)',advance="no")xs20,'   ||'
     
  end do
  write(o,'(A)')' ['//trim(url)//' link]                                   ||'
  
  
  
  !Injection values:
  write(o,'(///,A)')'== Injection values =='
  write(o,'(A)')"|| '''Detectors'''  || '''M1'''     || '''M2'''     || '''Mc'''     || '''&eta;'''  || '''time'''      || '''spin1'''  ||'''&theta;1'''|| '''spin2'''  ||'''&theta;2'''|| '''Dist'''   || '''R.A.'''   || '''Dec.'''   || '''incl'''   || '''pol.'''   || '''details'''                                                                                     ||"
  write(o,'(A)')"||                  ||  (Mo)        || (Mo)         || (Mo)         ||              ||  (s)            ||              || (rad)        ||              || (rad)        || (Mpc)        || (rad)        || (rad)        || (rad)        || (rad)        ||                                                                                                   ||"
  
  write(o,'(A4)',advance="no")'||  '
  do i=1,4
     if(i.le.ndet(ic)) then
        write(o,'(A2)',advance="no")detabbrs(detnr(ic,i))
     else
        write(o,'(A2)',advance="no")'  '
     end if
  end do
  write(o,'(A6)',advance="no")'    '
  parr(1:14) = (/63,64,61,62,11,71,72,81,82,22,31,32,51,52/)
  do p=1,14
     p1 = parr(p)
     if(revID(p1).eq.0) then  !Parameter not used
        write(o,'(A5,A11)',advance="no")'  || ',' - '
     else
        x = allDat(ic,revID(p1),1)
        if(p1.eq.31) x = rev2pi(x*rh2r)
        if(p1.eq.52) x = rrevpi(x*rd2r)  !Polarisation angle
        if(p1.eq.41.or.p1.eq.54.or.p1.eq.73.or.p1.eq.83) x = rev2pi(x*rd2r)
        if(p1.eq.32.or.p1.eq.51.or.p1.eq.53.or.p1.eq.72.or.p1.eq.82) x = x*rd2r
        if(p1.ge.11.and.p1.le.19) then
           write(o,'(A5,F14.4)',advance="no")'  || ',x+t0
        else
           write(o,'(A5,F11.4)',advance="no")'  || ',x
        end if
     end if
  end do
  write(o,'(A,79x,A)')'  || [ injection info] ',' ||'
  
  
  
  
  
  
  write(o,'(///,A)')'----'
  write(o,'(A)')'Back to [:JointS5/BayesianFollowUpOfE14Events:Bayesian follow-up in E14]'
  
  close(o)  !Wiki output file
  
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
!!  \todo:  do we need unwrapped data for this? - probably not, since the different chains are wrapped in the same way
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
  integer :: i,ic,p,nn,nn1
  integer :: lowVar(maxMCMCpar),nLowVar,highVar(maxMCMCpar),nHighVar,nmeanRelVar,nRhat,IDs(maxMCMCpar),nUsedPar
  real :: dx
  real*8 :: chMean(maxChs,maxMCMCpar),avgMean(maxMCMCpar),chVar(maxMCMCpar),chVar1(maxChs,maxMCMCpar),meanVar(maxMCMCpar),totRhat,meanRelVar,compute_median
  character :: ch
  
  if(contrChains.le.1) then
     if(prConv.ge.1) write(stdOut,'(A)')'  At least two chains are needed to compute R-hat...'
     return
  end if
  
  !< Determine the number of points to use for the computation of R-hat
  !<  Since we need a fixed number of points for all chains, take nn = the minimum of (the length of each chain after the burn-in) and use the last nn data points of each chain
  nn = nint(1e9)
  do ic=1,nChains0
     if(contrChain(ic).eq.0) cycle  !Contributing chains only
     nn1 = Ntot(ic)-Nburn(ic)
     if(nn1.lt.nn) nn = nn1
  end do
  
  if(prConv.ge.2) write(stdOut,'(A,I8,A)')'  Convergence parameters for',nn,' data points in each chain:'
  
  
  !Compute the means for each chain and for all chains:
  chMean = 1.d-30
  avgMean = 1.d-30
  do p=1,nMCMCpar
     do ic=1,nChains0
        if(contrChain(ic).eq.0) cycle  !Contributing chains only
        do i=Ntot(ic)-nn+1,Ntot(ic)
           chMean(ic,p) = chMean(ic,p) + allDat(ic,p,i)
        end do
        avgMean(p) = avgMean(p) + chMean(ic,p)
     end do
  end do
  chMean = chMean/dble(nn)
  avgMean = avgMean/dble(nn*contrChains)
  
  
  !Compute variances per chain, for all chains:
  chVar = 1.d-30
  chVar1 = 1.d-30
  meanVar = 1.d-30
  do p=1,nMCMCpar
     do ic=1,nChains0
        if(contrChain(ic).eq.0) cycle  !Contributing chains only
        do i=Ntot(ic)-nn+1,Ntot(ic)
           dx = (allDat(ic,p,i) - chMean(ic,p))**2
           chVar(p) = chVar(p) + dx
           chVar1(ic,p) = chVar1(ic,p) + dx !Keep track of the variance per chain
        end do
        meanVar(p) = meanVar(p) + (chMean(ic,p) - avgMean(p))**2
        chVar1(ic,p) = chVar1(ic,p)/dble(nn-1)
     end do
     chVar(p) = chVar(p)/dble(contrChains*(nn-1))
     meanVar(p) = meanVar(p)/dble(contrChains-1)
     
     !Compute Rhat:
     Rhat(p) = min( dble(nn-1)/dble(nn)  +  meanVar(p)/chVar(p) * (1.d0 + 1.d0/dble(contrChains)), 99.d0)
  end do
  
  
  !Print means per chain:
  if(prConv.ge.1) then
     write(stdOut,*)''
     write(stdOut,'(A14)',advance="no")''
     do p=1,nMCMCpar
        if(fixedpar(p).eq.1) cycle  !Varying parameters only
        write(stdOut,'(A9)',advance="no")trim(parNames(parID(p)))
     end do
     write(stdOut,'(A9)')'Mean'
     
     if(prConv.ge.3) then
        write(stdOut,'(A)')'  Means:'
        do ic=1,nChains0
           if(contrChain(ic).eq.0) cycle  !Contributing chains only
           write(stdOut,'(I12,A2)',advance="no")ic,': '
           do p=1,nMCMCpar
              if(fixedpar(p).eq.1) cycle  !Varying parameters only
              write(stdOut,'(F9.5)',advance="no")chMean(ic,p)
           end do
           write(stdOut,*)
        end do
        write(stdOut,'(A14)',advance="no")'        Mean: '
        do p=1,nMCMCpar
           if(fixedpar(p).eq.1) cycle  !Varying parameters only
           write(stdOut,'(F9.5)',advance="no")avgMean(p)
        end do
        write(stdOut,*)''
        
     end if !if(prConv.ge.3)
  end if !if(prConv.ge.1)
  
  
  !Flag and print variances:
  if(prConv.ge.3) write(stdOut,'(/,A)')'  Variances:'
  do ic=1,nChains0
     if(contrChain(ic).eq.0) cycle  !Contributing chains only
     lowVar = 0
     highVar = 0
     meanRelVar = 1.d0
     nmeanRelVar = 0
     do p=1,nMCMCpar
        if(fixedpar(p).eq.0 .and.abs(chVar1(ic,p)).gt.1.e-30) then  !Take only the parameters that were fitted and have a variance > 0
           if(chVar1(ic,p).lt.0.5*chVar(p)) lowVar(p) = 1  !Too (?) low variance, mark it
           if(chVar1(ic,p).gt.2*chVar(p))  highVar(p) = 1  !Too (?) high variance, mark it
           meanRelVar = meanRelVar * chVar1(ic,p)/chVar(p) !Multiply in order to take the geometric mean
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
        write(stdOut,'(I12,A2)',advance="no")ic,':'//ch
        do p=1,nMCMCpar
           if(fixedpar(p).eq.1) cycle  !Varying parameters only
           ch = ' '
           if(lowVar(p).eq.1) ch = '*'
           if(highVar(p).eq.1) ch = '#'
           write(stdOut,'(F8.4,A1)',advance="no")chVar1(ic,p),ch
        end do
        ch = ' '
        if(meanRelVar.lt.0.5) ch = '*'
        if(meanRelVar.gt.2.0) ch = '#'
        write(stdOut,'(F8.3,A1)',advance="no")meanRelVar,ch
        write(stdOut,*)''
     end if !if(prConv.ge.3)
  end do
  
  
  !Print mean variance for all parameters :
  if(prConv.ge.3) then
     write(stdOut,'(A13)',advance="no")'   Mean:'
     do p=1,nMCMCpar
        if(fixedpar(p).eq.1) cycle  !Varying parameters only
        write(stdOut,'(F9.4)',advance="no")chVar(p)
     end do
     write(stdOut,*)''
     write(stdOut,*)''
  end if !if(prConv.ge.3)
  
  !Print the variances within chains and between chains:
  if(prConv.ge.2) then
     write(stdOut,'(A)')'  Variances:'
     write(stdOut,'(A14)',advance="no")'      In chs: '
     do p=1,nMCMCpar
        if(fixedpar(p).eq.1) cycle  !Varying parameters only
        write(stdOut,'(ES9.1)',advance="no")chVar(p)
     end do
     write(stdOut,*)
     write(stdOut,'(A14)',advance="no")'   Betw. chs: '
     do p=1,nMCMCpar
        if(fixedpar(p).eq.1) cycle  !Varying parameters only
        write(stdOut,'(ES9.1)',advance="no")meanVar(p)
     end do
     write(stdOut,*)
  end if
  
  !Print R-hat:
  if(prConv.ge.1) then
     write(stdOut,'(A14)',advance="no")'       R-hat: '
     totRhat = 1.d0
     nRhat = 0
     do p=1,nMCMCpar
        if(fixedpar(p).eq.1) cycle  !Varying parameters only
        write(stdOut,'(F9.4)',advance="no")Rhat(p)
        !print*,parID(p)
        if(parID(p).ne.63.and.parID(p).ne.64) then !If not one of M1,M2
           !totRhat = totRhat + Rhat(p)   !Arithmetic mean
           !totRhat = totRhat * Rhat(p)   !Geometric mean
           nRhat = nRhat + 1
        end if
     end do
     !write(stdOut,'(F9.4)')totRhat/dble(nRhat)          !Arithmetic mean
     !write(stdOut,'(F9.4)')totRhat**(1.d0/dble(nRhat))  !Geometric mean
     write(stdOut,'(F9.4,A)')compute_median(Rhat(1:nMCMCpar),nMCMCpar),' (med)'  !Median
  end if
  
end subroutine compute_convergence
!***********************************************************************************************************************************



!***********************************************************************************************************************************
subroutine compute_autocorrelations()
  !> Compute autocorrelations
  !! Use the original nChains0 chains and allDat()
  !! Results are printed to stdout if prAcorr>0
  !< 
  
  use constants
  use mcmcrun_data
  use chain_data
  implicit none
  integer :: ic,j1,p,j,i,np
  real :: median,medians(nMCMCpar),compute_median_real,stdev,compute_stdev_real
  
  if(prAcorr.eq.0) then
     write(stdOut,'(A)',advance="no")' autocorrs, '
  else
     write(stdOut,*)
     if(prAcorr.ge.2) write(stdOut,'(A)')' Autocorrelations:'
     write(stdOut,'(A14)',advance="no")''
     do p=1,nMCMCpar
        if(fixedpar(p).eq.1) cycle  !Varying parameters only
        write(stdOut,'(A9)',advance="no")trim(parNames(parID(p)))
     end do
     write(stdOut,'(A11)')'Median'
  end if
  
  
  do ic=1,nChains0  !Use original data
     acorrs(ic,:,:) = 0.
     lAcorrs(ic,:) = 0.
     
     j1 = Ntot(ic)/nAcorr   !Step size to get nAcorr autocorrelations per parameter - don't waste any CPU
     if(j1.lt.1) then
        write(stdErr,'(A)')"  *** WARNING:  nAcorr too small or nAcorr too large to compute autocorrelations properly. I can't use all data points and the results may not be what you expect ***"
        j1 = 1
     end if
     
     if(prAcorr.ge.2) write(stdOut,'(A,I3,A)',advance="no")'  Chain',ic,':   '
     do p=1,nMCMCpar
        if(fixedpar(p).eq.1) cycle  !Varying parameters only
        
        median = compute_median_real(allDat(ic,p,1:Ntot(ic)),Ntot(ic))
        !median = sum(selDat(ic,p,1:Ntot(ic)))/real(Ntot(ic))  !Replace median with mean
        stdev  = compute_stdev_real(allDat(ic,p,1:Ntot(ic)),Ntot(ic),median)
        
        do j=0,min(nAcorr,Ntot(ic)-1)
           do i=1,Ntot(ic)-j*j1
              acorrs(ic,p,j) = acorrs(ic,p,j) + (allDat(ic,p,i) - median)*(allDat(ic,p,i+j*j1) - median)
           end do
           acorrs(ic,p,j) = acorrs(ic,p,j) / (stdev*stdev*(Ntot(ic)-j*j1))
           
           if(lAcorrs(ic,p).lt.1. .and. acorrs(ic,p,j).lt.0) lAcorrs(ic,p) = real(j*j1*totthin(ic))
           if(p.eq.1) acorrs(ic,0,j) = real(j*j1*totthin(ic))  !Make sure you get the iteration number, not the data-point number
        end do !j
        
        if(prAcorr.ge.2) write(stdOut,'(ES9.1)',advance="no")lAcorrs(ic,p)
     end do !p
     
     if(prAcorr.ge.2) write(stdOut,'(ES11.1)')compute_median_real(lAcorrs(ic,1:nMCMCpar),nMCMCpar)
  end do !ic
  
  
  if(prAcorr.ge.1) then
     if(prAcorr.eq.1) write(stdOut,'(A)',advance="no")' Med.Autocor: '
     if(prAcorr.ge.2) write(stdOut,'(/,A)',advance="no")'    Median:   '
     np = 0
     do p=1,nMCMCpar
        if(fixedpar(p).eq.1) cycle  !Varying parameters only
        np = np+1
        medians(np) = compute_median_real(lAcorrs(1:nChains0,p),nChains0)
        write(stdOut,'(ES9.1)',advance="no")medians(np)
     end do
     write(stdOut,'(ES11.1)')compute_median_real(medians(1:np),np)
  end if
  
end subroutine compute_autocorrelations
!***********************************************************************************************************************************