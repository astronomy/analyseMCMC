!> \file analyseMCMC_chains.f90  Plot chains (posterior, parameters, jumps, etc.) for analyseMCMC

! 
! LICENCE:
! 
! Copyright 2007-2013 Marc van der Sluys
!  
! This file is part of the AnalyseMCMC package.
!  
! This is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published
! by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! 
! This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License along with this code (LICENCE).  If not, see 
! <http://www.gnu.org/licenses/>.
! 



!***********************************************************************************************************************************
!> \brief  Plot chains (posterior, parameters, jumps, etc.) for analyseMCMC
!!
!! \retval exitcode  Exit status code (0=ok)

subroutine chains(exitcode)
  use analysemcmc_settings, only: plLogL,plChain,plParL,plJump,plAcorr,plRhat
  
  implicit none
  integer, intent(out) :: exitcode
  
  
  exitcode = 0
  
  ! Plot posterior chain:
  if(plLogL.gt.0) call plot_posterior_chain(exitcode)
  if(exitcode.ne.0) return
  
  
  ! Plot Markov chains per parameter:
  if(plChain.gt.0) call plot_parameter_chains(exitcode)
  if(exitcode.ne.0) return
  
  
  ! Plot L vs parameter value:
  if(plParL.gt.0) call plot_par_L(exitcode)
  if(exitcode.ne.0) return
  
  
  ! Plot jump sizes:
  if(plJump.gt.0) call plot_Jump_sizes(exitcode)
  if(exitcode.ne.0) return
  
  
  ! Plot autocorrelations:
  if(plAcorr.gt.0) call plot_Acorr_chains(exitcode)
  if(exitcode.ne.0) return
  
  
  ! Plot R-hat:
  if(plRhat.gt.0) call plot_Rhat_chains(exitcode)
  if(exitcode.ne.0) return
  
end subroutine chains
!***********************************************************************************************************************************




!***********************************************************************************************************************************
!> \brief  Plot posterior chain
!!
!! \retval exitcode  Exit status code (0=ok)
!!

subroutine plot_posterior_chain(exitcode)
  use SUFR_constants, only: stdOut,stdErr
  
  use aM_constants, only: use_PLplot, rmeps
  use analysemcmc_settings, only: update,prProgress,file,scrsz,scrrat,pssz,psrat,fonttype,colour,whitebg,quality,scLogLpl
  use analysemcmc_settings, only: Nburn,plLmax,plBurn,chainPlI,fontsize1d,nPlPar, chainSymbol,plInject, autoBurnin, htmlOutput
  use general_data, only: post,outputname,outputdir,nChains0,Ntot,startval,icloglmax,iloglmax
  use plot_data, only: psclr,defcolour,bmpsz,bmprat,ncolours,colours,unSharplogl,bmpxpix,nsymbols,symbols
  use chain_data, only: is,isburn
  
  implicit none
  integer, intent(out) :: exitcode
  
  integer :: i,pgopen,imin,ci,lw,symbol,io,ic,status,system
  real :: dx,dy,xmin,xmax,ymin,ymax,ply, sch
  character :: tempfile*(199), convopts*(99)
  logical :: ex
  
  
  exitcode = 0
  
  if(htmlOutput.ge.1) then
     write(stdOut,'(A)') '<a name="postchains"></a>'
     write(stdOut,'(A)') '<font size="1"><a href="#top" title="Go to the top of the page">top</a></font>'
     write(stdOut,'(A)') '<h3>Posterior chain:</h3>'
     write(stdOut,'(A)') '<img src="'//trim(outputname)//'__posterior.png" title="log Posterior">'
  else
     if(prProgress.ge.1.and.update.eq.0) write(stdOut,'(A)',advance="no")' posterior chain, '
  end if
  
  if(file.eq.0) then
     io = pgopen('12/xs')
     sch = 1.5*fontsize1D
     lw = 1
  end if
  if(file.ge.1) then
     tempfile = trim(outputdir)//'/'//trim(outputname)//'__posterior'
     if(file.eq.1) io = pgopen(trim(tempfile)//'.ppm/ppm')
     if(file.ge.2) io = pgopen(trim(tempfile)//'.eps'//trim(psclr))
     sch = 1.2*fontsize1D
     lw = 1
     if(file.ge.2) then
        lw = lw*2
        if(nPlPar.eq.1) lw = nint(2*fontsize1d)
     end if
  end if
  
  if(io.le.0) then
     write(stdErr,'(A,I4)')'   Error:  Cannot open PGPlot device.  Quitting the programme',io
     exitcode = 1
     return
  end if
  
  if(file.eq.0) then
     call pgpap(scrSz,scrRat)
     sch = 1.5*fontsize1D
  end if
  if(file.eq.1) then
     call pgpap(bmpsz,bmprat)
     sch = 1.5*fontsize1D
  end if
  if(file.ge.2) then
     call pgpap(PSsz,PSrat)
     call pgscf(fonttype)
     sch = fontsize1D
  end if
  
  ! Custom initialisation of PGPlot/PLplot:
  call pginitl(colour,file,whiteBG)
  call pgslw(lw)
  call pgsch(sch)
  
  
  !call pgsvp(0.08,0.92,0.07,0.94)  ! Need wider margin for title on right
  call pgsvp(0.08*fontsize1D,1.-0.08*fontsize1D,0.07*fontsize1D,1.-0.06*fontsize1D)  ! Need wider margin for title on right
  if(quality.eq.0) call pgsvp(0.08,0.92,0.07,0.94)  ! To make room for title
  
  ic = 1
  xmax = -1.e30
  ymin =  1.e30
  ymax = -1.e30
  do ic=1,nChains0
     xmin = 0.
     xmax = max(xmax,maxval(is(ic,1:Ntot(ic))))
     imin = 10                                              ! Take into account burn-in
     if(scLogLpl.eq.1) imin = Nburn(ic)                     ! Scale without taking into account burn-in
     ymin = min(ymin,minval(post(ic,imin:Ntot(ic)))) 
     ymax = max(ymax,maxval(post(ic,imin:Ntot(ic))))
  end do
  ic = 1
  if(scLogLpl.eq.0) then                                    ! Take into account 0, injection and starting values
     ymin = min(ymin,post(ic,1),post(ic,2),0.)
     ymax = max(ymax,post(ic,1),post(ic,2),0.)
  else                                                      ! Take into account injection values only
     ymin = min(ymin,post(ic,1))
     ymax = max(ymax,post(ic,1))
  end if
  
  ymin = max(0.,ymin)
  dx = abs(xmax-xmin)*0.01
  dy = abs(ymax-ymin)*0.05
  xmin = xmin - dx
  xmax = xmax + dx
  ymin = ymin - dy
  ymax = ymax + dy
  
  
  call pgswin(xmin,xmax,ymin,ymax)
  !call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)
  call plot_posterior_snr_axes(xmin,xmax,ymin,ymax)
  
  if(quality.ne.1 .and. abs(startval(1,1,1)-startval(1,1,2))/abs((startval(1,1,1))+tiny(startval)).gt.1.e-10) then
     call pgsls(4)
     call pgbox('',0.0,0,'G',0.0,0)  ! Plot a grid of horizontal lines
     call pgsls(1)
  end if
  
  do ic=1,nChains0
     ! Give pre- and post-burn-in different colour in posterior chain:
     ci = defcolour
     if(nChains0.gt.1) ci = colours(mod(ic-1,ncolours)+1)
     
     symbol = chainSymbol
     if(chainSymbol.le.-10) symbol = symbols(mod(ic-1,nsymbols)+1)
     
     ! Pre-burnin:
     call pgsci(ci)
     if(plBurn.ge.2) call pgscidark(ci,file,whiteBG)
     do i=ic,Nburn(ic),chainPlI  ! Start at ic to reduce overplotting
        call pgpoint(1,is(ic,i),post(ic,i), symbol)
     end do
     
     ! Post-burnin:
     call pgsci(ci)
     do i=Nburn(ic)+ic,Ntot(ic),chainPlI  ! Start at ic to reduce overplotting
        call pgpoint(1,is(ic,i),post(ic,i), symbol)
     end do
  end do
  
  
  ! Plot max posterior:
  if(plLmax.ge.1) then
     ply = post(icloglmax,iloglmax)
     call pgsci(1)
     call pgpoint(1,is(icloglmax,iloglmax),ply,18)
     call pgsls(5)
     call pgline(2,(/xmin,xmax/),(/ply,ply/))
     if(plLmax.ge.2 .and. abs(autoBurnin).gt.1.e-10) then
        ply = ply - autoBurnin
        call pgline(2,(/xmin,xmax/),(/ply,ply/))
        print*,ply,autoBurnin
     end if
  end if
  
  do ic=1,nChains0
     call pgsci(1)
     call pgsls(2)
     
     ! Plot injection value, only if injection was done: dashed horizontal line:
     if(plInject.ge.1 .and. abs(post(ic,1)).gt.1.e-4) call pgline(2,(/xmin,xmax/),(/post(ic,1),post(ic,1)/))  
     call pgsci(colours(mod(ic-1,ncolours)+1))
     
     ! Mark the end of the burn-in phase: dashed vertical line:
     if((plBurn.eq.1.or.plBurn.ge.3).and.isburn(ic).lt.is(ic,Ntot(ic))) call pgline(2,(/isburn(ic),isburn(ic)/), &
          (/ymin,ymax/))
     
     ! Starting value: horizontal dotted line:
     call pgsls(4)
     call pgline(2,(/xmin,xmax/),(/post(ic,2),post(ic,2)/))
  end do
  
  call pgsci(1)
  call pgsls(1)
  if(use_PLplot) then
     call pgmtxt('L',5.0,0.5,0.5,'log Posterior')
     call pgmtxt('R',5.0,0.5,0.5,'\(2267)(2logP) \(0248) SNR')
  else
     call pgmtxt('L',2.3,0.5,0.5,'log Posterior')
     call pgmtxt('R',2.5,0.5,0.5,'\(2267)(2logP) \(0248) SNR')
  end if
  if(nPlPar.eq.1.and.quality.eq.1) call pgmtxt('B',2.5,0.5,0.5,'iteration')
  
  if(quality.eq.0) then
     call pgmtxt('T',0.5,0.9,0.9,trim(outputname))  !Print title
  end if
  
  call pgend()
  
  if(file.ge.2) then
     if(file.eq.3) then
        status = system('eps2pdf '//trim(tempfile)//'.eps  -o '//trim(tempfile)//'.pdf   >& /dev/null')
        if(status.ne.0) then
           write(stdErr,'(A,I6)')'  Error converting plot eps -> pdf',status
        else
           if(rmeps) status = system('rm -f '//trim(tempfile)//'.eps  >& /dev/null')
        end if
     end if
  else if(file.eq.1) then
     inquire(file=trim(tempfile)//'.ppm', exist=ex)
     if(ex) then
        convopts = '-resize '//trim(bmpxpix)//' -depth 8 -unsharp '//trim(unSharplogl)
        status = system('convert '//trim(convopts)//' '//trim(tempfile)//'.ppm '//trim(tempfile)//'.png')
        if(status.ne.0) write(stdErr,'(A,I6)')'  Error converting plot ppm -> png',status
        status = system('rm -f '//trim(tempfile)//'.ppm')
        !if(status.ne.0) write(stdErr,'(A)')'  Error removing file',status
     end if
  end if
  
end subroutine plot_posterior_chain
!***********************************************************************************************************************************


  

!***********************************************************************************************************************************
!> \brief  Plot chains for each parameter
!!
!! \retval exitcode  Exit status code (0=ok)
!!

subroutine plot_parameter_chains(exitcode)
  use SUFR_constants, only: stdOut,stdErr
  
  use aM_constants, only: use_PLplot, rmeps
  use analysemcmc_settings, only: update,prProgress,file,scrsz,scrrat,pssz,psrat,fonttype,colour,whitebg,quality
  use analysemcmc_settings, only: Nburn,plLmax,plBurn,chainPlI,fontsize1d,nPlPar,panels,plPars,scChainsPl,changeVar
  use analysemcmc_settings, only: chainSymbol,plInject,mergeChains,plStart,prConv, htmlOutput
  use general_data, only: allDat,outputname,outputdir,nChains0,Ntot,startval,icloglmax,iloglmax,parNames,pgParNs,rhat
  use mcmcrun_data, only: revID,parID
  use plot_data, only: psclr,defcolour,bmpsz,bmprat,ncolours,colours,bmpxpix,nsymbols,symbols,unSharpchain
  use chain_data, only: is,isburn
  
  implicit none
  integer, intent(out) :: exitcode
  
  integer :: i,j,pgopen,imin,ci,lw,symbol,io,ic,p,status,system
  real :: rev360,rev24,rev180
  real :: dx,dy,xmin,xmax,ymin,ymax,sch,plx,ply
  character :: title*(99), tempfile*(199), convopts*(99)
  logical :: ex
  
  
  exitcode = 0
  
  if(htmlOutput.ge.1) then
     write(stdOut,'(A)') '<a name="parchains"></a>'
     write(stdOut,'(A)') '<font size="1"><a href="#top" title="Go to the top of the page">top</a></font>'
     write(stdOut,'(A)') '<h3>Parameter chains:</h3>'
     write(stdOut,'(A)') '<img src="'//trim(outputname)//'__chains.png" title="Parameter chains">'
  else
     if(prProgress.ge.1.and.update.eq.0) write(stdOut,'(A)',advance="no")' parameter chains, '
  end if
  
  if(file.eq.0) then
     io = pgopen('13/xs')
     call pgpap(scrSz,scrRat)
     
     sch = fontsize1d*1.5
     lw = 1
  end if
  if(file.ge.1) then
     tempfile = trim(outputdir)//'/'//trim(outputname)//'__chains'
     if(file.eq.1) then
        io = pgopen(trim(tempfile)//'.ppm/ppm')
        call pgpap(bmpsz,bmprat)
     end if
     if(file.ge.2) then
        io = pgopen(trim(tempfile)//'.eps'//trim(psclr))
        call pgpap(PSsz,PSrat)
        
        call pgscf(fonttype)
     end if
     
     
     lw = 3
     if(nPlPar.ge.10) lw = 2
     if(quality.lt.2) lw = max(lw-1,1)  ! Draft/Paper
     sch = fontsize1d*1.2
     !if(nChains.eq.1.and.nPlPar.gt.9) sch = fontsize1d*1.2
     if(quality.eq.0) then  ! Draft
        sch = sch*1.75
        lw = 2
     end if
     if(quality.eq.1) then  ! Paper
        if(nPlPar.eq.1) then
           sch = fontsize1d*2
           lw = nint(2*fontsize1d)
        else if(nPlPar.eq.12) then
           sch = sch*1.75
           lw = 2
        else
           sch = sch*1.25
           lw = 1
        end if
     end if
     if(quality.eq.2) then  ! Talk
        if(nPlPar.gt.12) then
           sch = sch*1.5
           lw = 1
        end if
        if(nPlPar.le.12) then
           sch = sch*2
           lw = 2
        end if
        !if(nPlPar.le.6) then
        !   !sch = sch*4
        !   !lw = 4
        !end if
     end if
     if(quality.eq.3) then !Poster
        if(nPlPar.eq.12.and.file.ge.2) then
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
  end if  ! if(file.ge.1)
  
  if(io.le.0) then
     write(stdErr,'(A,I4)')'   Error:  Cannot open PGPlot device.  Quitting the programme',io
     exitcode = 1
     return
  end if
  
  if(file.eq.1) call pgsch(1.5)
  if(use_PLplot) then
     sch = sch*0.6
  else
     if(file.eq.1) sch = sch*0.8
  end if
  call pgsch(sch)
  call pgslw(lw)
  
  
  call pgsubp(panels(1),panels(2))
  
  ic = 1
  do j=1,nPlPar
     p = revID(plPars(j))
     if(p.eq.0) then
        call report_undefined_parameter(trim(parNames(plPars(j))), 'plot_parameter_chains')
        cycle
     end if
     
     if(nPlPar.gt.1) call pgpage()
     if(j.eq.1) call pginitl(colour,file,whiteBG)
     
     call pgsch(sch/2.*sqrt(real(nPlPar)))
     
     if(file.eq.0.and.scrRat.gt.1.35) call pgsvp(0.08,0.95,0.1,0.95)
     if(file.eq.1) then
        if(bmprat.gt.1.35) then
           call pgsvp(0.08,0.95,0.1,0.95)
        else
           if(use_PLplot) then
              call pgsvp(0.1,0.92,0.12,0.9)  ! default case?
           else
              call pgsvp(0.1,0.95,0.08,0.9)  ! default case?
           end if
        end if
     end if
     if(file.ge.2.and.PSrat.gt.1.35) & 
          !call pgsvp(0.08,0.95,0.1,0.94)
          call pgsvp(0.08*sqrt(fontsize1D),1.0-0.04-0.04*sqrt(fontsize1D),0.07*fontsize1D**1.5,1.0-0.05-0.03*fontsize1D**1.5)
     if(file.ge.2.and.PSrat.le.1.35.and.nPlPar.eq.1) &
          !call pgsvp(0.08,0.92,0.07,0.94)  ! Same as logP plot
          call pgsvp(0.08*fontsize1D,1.0-0.08*fontsize1D,0.07*fontsize1D,1.0-0.06*fontsize1D)
     
     if(quality.eq.0) call pgsvp(0.08,0.95,0.06,0.87)  ! To make room for title
     if(quality.eq.4) call pgsvp(0.13,0.95,0.1,0.95)
     
     xmin = 0.
     !xma x = real(maxval(Ntot(1:nChains0)))
     xmax = -1.e30
     ymin =  1.e30
     ymax = -1.e30
     do ic=1,nChains0
        xmax = max(xmax,maxval(is(ic,1:Ntot(ic))))
        imin = 1                                               ! Take into account burn-in
        if(scChainsPl.eq.1) imin = Nburn(ic)                   ! Scale without taking into account burn-in
        ymin = min(ymin,minval(allDat(ic,p,imin:Ntot(ic))))
        ymax = max(ymax,maxval(allDat(ic,p,imin:Ntot(ic))))
     end do
     
     if(changeVar.ge.1) then
        select case(parID(p))
        case(31) !RA
           if(ymin.lt.0..or.ymax.gt.24.) then
              ymin = 0.
              ymax = 24.
           end if
        case(41,54,73,83) !Phases
           if(ymin.lt.0..or.ymax.gt.360.) then
              ymin = 0.
              ymax = 360.
           end if
        case(52) !Psi
           if(ymin.lt.0..or.ymax.gt.180.) then
              ymin = 0.
              ymax = 180.
           end if
        end select
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
     xmin = xmin - dx
     xmax = xmax + dx
     
     call pgswin(xmin,xmax,ymin-dy,ymax+dy)
     if(quality.eq.1) call pgswin(xmin,xmax,ymin-dy,ymax+dy*2)
     call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)
     
     ! Plot the actual chain values:
     call pgsch(sch*0.5)
     if(nPlPar.ne.1.and.chainSymbol.ne.1) call pgsch(sch*0.7*0.5)
     
     call pgslw(1)
     !write(stdOut,'(15I4)'),nsymbols,symbols(1:nsymbols)
     do ic=1,nChains0
        symbol = chainSymbol
        if(chainSymbol.le.-10) symbol = symbols(mod(ic-1,nsymbols)+1)
        call pgsci(defcolour)
        if(nChains0.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
        if(chainSymbol.eq.0) then !Plot lines rather than symbols
           call pgline(Ntot(ic),is(ic,1:Ntot(ic)),allDat(ic,p,1:Ntot(ic)))
        else
           
           ! Give pre- and post-burn-in different colour in parameter-chain plots:
           ci = defcolour
           if(nChains0.gt.1) ci = colours(mod(ic-1,ncolours)+1)
           
           ! Pre-burnin:
           call pgsci(ci)
           if(plBurn.ge.2) call pgscidark(ci,file,whiteBG)
           do i=ic,Nburn(ic),chainPlI !Start at ic to reduce overplotting
              ply = allDat(ic,p,i)
              if(changeVar.ge.1) then
                 select case(parID(p))
                 case(31) 
                    ply = rev24(ply)  !RA
                 case(41,54,73,83) 
                    ply = rev360(ply)
                 case(52) 
                    ply = rev180(ply)
                 end select
              end if
              call pgpoint(1,is(ic,i),ply,symbol)
           end do
           
           ! Post-burnin:
           call pgsci(ci)
           do i=Nburn(ic)+ic,Ntot(ic),chainPlI !Start at ic to reduce overplotting
              ply = allDat(ic,p,i)
              if(changeVar.ge.1) then
                 select case(parID(p))
                 case(31) 
                    ply = rev24(ply)  !RA
                 case(41,54,73,83) 
                    ply = rev360(ply)
                 case(52) 
                    ply = rev180(ply)
                 end select
              end if
              call pgpoint(1,is(ic,i),ply,symbol)
           end do
           
           
        end if
     end do
     
     call pgsch(sch/2.*sqrt(real(nPlPar)))
     
     call pgslw(lw)
     
     ! Plot max posterior:
     if(plLmax.ge.1) then
        ply = allDat(icloglmax,p,iloglmax)
        if(changeVar.ge.1) then
           select case(parID(p))
           case(31) 
              ply = rev24(ply)  !RA
           case(41,54,73,83) 
              ply = rev360(ply)
           case(52) 
              ply = rev180(ply)
           end select
        end if
        call pgsci(1)
        call pgpoint(1,is(icloglmax,iloglmax),ply,18)
        call pgsls(5)
        call pgline(2,(/xmin,xmax/),(/ply,ply/))
     end if
     
     
     ! Plot burn-in, injection and starting values:
     do ic=1,nChains0
        call pgsls(2)
        call pgsci(6)
        
        ! Mark end of burn-in phase:
        if(nChains0.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
        if((plBurn.eq.1.or.plBurn.ge.3).and.isburn(ic).lt.is(ic,Ntot(ic))) call pgline(2,(/isburn(ic),isburn(ic)/), &
             (/ymin-dy,ymax+2*dy/))
        call pgsci(1)
        
        
        ! Plot injection values in chains:
        if(plInject.ge.1) then
           if(mergeChains.ne.1.or.ic.eq.1) then 
              !The units of the injection values haven't changed (e.g. from rad to deg) for ic>1
              !(but they have for the starting values, why?)
              
              plx = startval(ic,p,1) !Injection value
              plx = max(min(1.e30,startval(ic,p,1)),1.e-30)
              if(changeVar.ge.1) then
                 select case(parID(p))
                 case(31) 
                    plx = rev24(plx)  !RA
                 case(41,54,73,83) 
                    plx = rev360(plx)
                 case(52) 
                    plx = rev180(plx)
                 end select
              end if
              call pgline(2,(/xmin,xmax/),(/plx,plx/))
              if(changeVar.ge.1) then
                 select case(parID(p))
                 case(31)
                    call pgline(2,(/xmin,xmax/),(/plx-24.,plx-24./))
                    call pgline(2,(/xmin,xmax/),(/plx+24.,plx+24./))
                 case(41,54,73,83)
                    call pgline(2,(/xmin,xmax/),(/plx-360.,plx-360./))
                    call pgline(2,(/xmin,xmax/),(/plx+360.,plx+360./))
                 case(52)
                    call pgline(2,(/xmin,xmax/),(/plx-180.,plx-180./))
                    call pgline(2,(/xmin,xmax/),(/plx+180.,plx+180./))
                 end select
              end if
           end if
        end if
        
        ! Plot starting values in chains:
        if(plStart.ge.1.and.abs((startval(ic,p,1)-startval(ic,p,2))/(startval(ic,p,1)+sqrt(tiny(startval)))) .gt. 1.e-10) then
           call pgsls(4)
           if(nChains0.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
           plx = startval(ic,p,2) !Initial value
           if(changeVar.ge.1) then
              select case(parID(p))
              case(31) 
                 plx = rev24(plx)
              case(41,54,73,83) 
                 plx = rev360(plx)
              case(52) 
                 plx = rev180(plx)
              end select
           end if
           call pgline(2,(/xmin,xmax/),(/plx,plx/))
           if(changeVar.ge.1) then
              select case(parID(p))
              case(31)
                 call pgline(2,(/xmin,xmax/),(/plx-24.,plx-24./))
                 call pgline(2,(/xmin,xmax/),(/plx+24.,plx+24./))
              case(41,54,73,83)
                 call pgline(2,(/xmin,xmax/),(/plx-360.,plx-360./))
                 call pgline(2,(/xmin,xmax/),(/plx+360.,plx+360./))
              case(52)
                 call pgline(2,(/xmin,xmax/),(/plx-180.,plx-180./))
                 call pgline(2,(/xmin,xmax/),(/plx+180.,plx+180./))
              end select
           end if
        end if
     end do !ic=1,nChains0
     
     
     ! Print labels:
     call pgsci(1)
     call pgsls(1)
     write(title,'(F6.3)') rhat(p)
     if(nPlPar.eq.1) then
        call pgmtxt('L',2.0,0.5,0.5,' '//trim(pgParNs(parID(p))))
        if(quality.eq.1) call pgmtxt('B',2.5,0.5,0.5,'iteration')
        if(nChains0.gt.1.and.prConv.ge.1) call pgmtxt('T',2.0,1.,1.,'R-hat: '//trim(title))
     else
        if(use_PLplot) then
           call pgmtxt('T',2.4,0.,0.,' '//trim(pgParNs(parID(p))))
           if(nChains0.gt.1.and.prConv.ge.1) call pgmtxt('T',2.4,1.,1.,'R-hat: '//trim(title))
        else
           call pgmtxt('T',0.6,0.,0.,' '//trim(pgParNs(parID(p))))
           if(nChains0.gt.1.and.prConv.ge.1) call pgmtxt('T',0.6,1.,1.,'R-hat: '//trim(title))
        end if
     end if
     
     
  end do !do j=1,nPlPar
  
  
  if(quality.eq.0) then
     call pgsubp(1,1)
     call pgsvp(0.,1.,0.,1.)
     call pgswin(-1.,1.,-1.,1.)
     
     call pgsch(sch*0.8)
     call pgmtxt('T',-0.7,0.5,0.5,trim(outputname))  ! Print titel
     call pgsch(sch)
  end if
  
  call pgend()
  if(file.ge.2) then
     if(file.eq.3) then
        status = system('eps2pdf '//trim(tempfile)//'.eps  -o '//trim(tempfile)//'.pdf  >& /dev/null')
        if(status.ne.0) then
           write(stdErr,'(A,I6)')'  Error converting plot eps -> pdf',status
        else
           if(rmeps) status = system('rm -f '//trim(tempfile)//'.eps  >& /dev/null')
        end if
     end if
  else if(file.eq.1) then
     inquire(file=trim(tempfile)//'.ppm', exist=ex)
     if(ex) then
        convopts = '-resize '//trim(bmpxpix)//' -depth 8 -unsharp '//trim(unSharpchain)
        status = system('convert '//trim(convopts)//' '//trim(tempfile)//'.ppm  '//trim(tempfile)//'.png')
        if(status.ne.0) write(stdErr,'(A,I6)')'  Error converting plot ppm -> png',status
        status = system('rm -f '//trim(tempfile)//'.ppm')
     end if
  end if
  
end subroutine plot_parameter_chains
!*********************************************************************************************************************************
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  


!***********************************************************************************************************************************
!> \brief  Plot L vs parameter value
!!
!! \retval exitcode  Exit status code (0=ok)
!!

subroutine plot_par_L(exitcode)
  use SUFR_constants, only: stdOut,stdErr
  
  use aM_constants, only: rmeps
  use analysemcmc_settings, only: update,prProgress,file,scrsz,scrrat,pssz,psrat,fonttype,colour,whitebg,quality, htmlOutput
  use analysemcmc_settings, only: Nburn,plLmax,chainPlI,fontsize1d,nPlPar,panels,plPars,changeVar, chainSymbol,plInject,mergeChains
  use general_data, only: post,allDat,outputname,outputdir,nChains0,Ntot,startval,icloglmax,iloglmax,nChains,parNames
  use mcmcrun_data, only: revID,parID
  use plot_data, only: psclr,defcolour,bmpsz,bmprat,ncolours,colours,bmpxpix,nsymbols,symbols,unSharpchain
  
  implicit none
  integer, intent(out) :: exitcode
  
  integer :: i,j,pgopen,lw,symbol,io,ic,p,status,system
  real :: rev360,rev24,rev180
  real :: dx,dy,xmin,xmax,ymin,ymax,sch,plx,ply
  character :: tempfile*(199), convopts*(99)
  logical :: ex
  
  
  exitcode = 0
  
  if(htmlOutput.ge.1) then
     write(stdOut,'(A)') '<a name="par-l"></a>'
     write(stdOut,'(A)') '<font size="1"><a href="#top" title="Go to the top of the page">top</a></font>'
     write(stdOut,'(A)') '<h3>Parameter-L:</h3>'
     write(stdOut,'(A)') '<img src="'//trim(outputname)//'__logl.png" title="Posterior vs. parameter value">'
  else
     !if(prProgress.ge.1.and.update.eq.0) write(stdOut,'(A)')' Plotting parameter-L plot...'
     if(prProgress.ge.1.and.update.eq.0) write(stdOut,'(A)',advance="no")' parameter-L, '
  end if
  
  if(file.eq.0) then
     io = pgopen('22/xs')
     call pgpap(scrSz,scrRat)
     
     sch = fontsize1d*1.5
     lw = 1
  end if
  if(file.ge.1) then
     tempfile = trim(outputdir)//'/'//trim(outputname)//'__logl'
     if(file.eq.1) then
        io = pgopen(trim(tempfile)//'.ppm/ppm')
        call pgpap(bmpsz,bmprat)
     end if
     if(file.ge.2) then
        io = pgopen(trim(tempfile)//'.eps'//trim(psclr))
        call pgpap(PSsz,PSrat)
        
        call pgscf(fonttype)
     end if
     
     
     lw = 3
     if(nPlPar.ge.10) lw = 2
     if(quality.lt.2) lw = max(lw-1,1)  !Draft/Paper
     sch = fontsize1d*1.2
     if(nChains.eq.1.and.nPlPar.gt.9) sch = fontsize1d*1.2
     if(quality.eq.0) then !Draft
        sch = sch*1.75
        lw = 2
     end if
     if(quality.eq.1) then !Paper
        if(nPlPar.le.12) then
           sch = sch*1.75
           lw = 2
        else
           sch = sch*1.25
           lw = 1
        end if
     end if
     if(quality.eq.2) then !Talk
        if(nPlPar.gt.12) then
           sch = sch*1.5
           lw = 1
        end if
        if(nPlPar.le.12) then
           sch = sch*2
           lw = 2
        end if
        !if(nPlPar.le.6) then
        !   !sch = sch*4
        !   !lw = 4
        !end if
     end if
     if(quality.eq.3) then !Poster
        if(nPlPar.eq.12.and.file.ge.2) then
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
     write(stdErr,'(A,I4)')'   Error:  Cannot open PGPlot device.  Quitting the programme',io
     exitcode = 1
     return
  end if
  !if(file.eq.1) call pgsch(1.5)
  
  call pgsch(sch)
  call pgslw(lw)
  
  if(file.eq.0.and.scrRat.gt.1.35) call pgsvp(0.08,0.95,0.1,0.95)
  if(file.eq.1.and.bmprat.gt.1.35) call pgsvp(0.08,0.95,0.1,0.95)
  if(file.ge.2.and.PSrat.gt.1.35) call pgsvp(0.08,0.95,0.1,0.95)
  if(quality.eq.0) call pgsvp(0.08,0.95,0.06,0.87) !To make room for title
  
  call pgsubp(panels(1),panels(2))
  
  ic = 1
  do j=1,nPlPar
     p = revID(plPars(j))
     if(p.eq.0) then
        call report_undefined_parameter(trim(parNames(plPars(j))), 'plot_par_L')
        cycle
     end if
     
     call pgpage()
     if(j.eq.1) call pginitl(colour,file,whiteBG)
     call pgsch(sch)
     
     xmin = 1.e30
     xmax = -1.e30
     ymin =  1.e30
     ymax = -1.e30
     do ic=1,nChains0
        !xmin = min(xmin,minval(is(ic,Nburn(ic):Ntot(ic))))
        !xmax = max(xmax,maxval(is(ic,Nburn(ic):Ntot(ic))))
        xmin = min(xmin,minval(allDat(ic,p,Nburn(ic):Ntot(ic))))
        xmax = max(xmax,maxval(allDat(ic,p,Nburn(ic):Ntot(ic))))
        ymin = min(ymin,minval(post(ic,Nburn(ic):Ntot(ic))))
        ymax = max(ymax,maxval(post(ic,Nburn(ic):Ntot(ic))))
     end do
     
     if(changeVar.ge.1) then
        select case(parID(p))
        case(31)
           if(ymin.lt.0..or.ymax.gt.24.) then
              ymin = 0.
              ymax = 24.
           end if
        case(41,54,73,83)
           if(ymin.lt.0..or.ymax.gt.360.) then
              ymin = 0.
              ymax = 360.
           end if
        case(52) 
           if(ymin.lt.0..or.ymax.gt.180.) then
              ymin = 0.
              ymax = 180.
           end if
        end select
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
     if(chainSymbol.ne.1) call pgsch(0.7)
     call pgslw(1)
     !write(stdOut,'(15I4)'),nsymbols,symbols(1:nsymbols)
     do ic=1,nChains0
        !call pgsci(mod(ic*2,10))
        !symbol = ic+1
        symbol = chainSymbol
        if(chainSymbol.le.-10) symbol = symbols(mod(ic-1,nsymbols)+1)
        call pgsci(defcolour)
        if(nChains0.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
        do i=Nburn(ic),Ntot(ic),chainPlI
           plx = allDat(ic,p,i)
           ply = post(ic,i)
           if(changeVar.ge.1) then
              select case(parID(p))
              case(31) 
                 plx = rev24(plx)
              case(41,54,73,83) 
                 plx = rev360(plx)
              case(52) 
                 plx = rev180(plx)
              end select
              !call pgpoint(1,is(ic,i),plx,1) !Plot small dots
           end if
           !call pgpoint(1,plx,ply,symbol) !Plot symbols
           call pgpoint(1,plx,exp(ply-ymin),symbol) !Plot symbols
           !print*,i,plx,ply,exp(ply-ymin)
        end do
     end do
     call pgsch(sch)
     call pgslw(lw)
     
     
     ! Plot max posterior:
     if(plLmax.ge.1) then
        plx = allDat(icloglmax,p,iloglmax)
        if(changeVar.ge.1) then
           select case(parID(p))
           case(31) 
              plx = rev24(plx)
           case(41,54,73,83) 
              plx = rev360(plx)
           case(52) 
              plx = rev180(plx)
           end select
        end if
        ply = exp(post(icloglmax,iloglmax)-ymin)
        call pgsci(1)
        call pgpoint(1,plx,ply,12)
        call pgsls(5)
        call pgline(2,(/plx,plx/),(/ymin,ymax/))
     end if
     
     
     !Plot injection values
     do ic=1,nChains0
        call pgsls(2)
        call pgsci(1)
        
        !Plot injection values
        if(plInject.ge.1) then
           if(mergeChains.ne.1.or.ic.eq.1) then 
              !The units of the injection values haven't changed (e.g. from rad to deg) for ic>1 
              ! (but they have for the starting values, why?)
              plx = startval(ic,p,1) !Injection value
              plx = max(min(1.e30,startval(ic,p,1)),1.e-30)
              if(changeVar.ge.1) then
                 select case(parID(p))
                 case(31) 
                    plx = rev24(plx)
                 case(41,54,73,83) 
                    plx = rev360(plx)
                 case(52) 
                    plx = rev180(plx)
                 end select
              end if
              call pgline(2,(/plx,plx/),(/ymin,ymax/))
              if(changeVar.ge.1) then
                 select case(parID(p))
                 case(31)
                    call pgline(2,(/plx-24.,plx-24./),(/ymin,ymax/))
                    call pgline(2,(/plx+24.,plx+24./),(/ymin,ymax/))
                 case(41,54,73,83)
                    call pgline(2,(/plx-360.,plx-360./),(/ymin,ymax/))
                    call pgline(2,(/plx+360.,plx+360./),(/ymin,ymax/))
                 case(52)
                    call pgline(2,(/plx-180.,plx-180./),(/ymin,ymax/))
                    call pgline(2,(/plx+180.,plx+180./),(/ymin,ymax/))
                 end select
              end if
           end if
        end if
     end do !ic=1,nChains0
     
     call pgsci(1)
     call pgsls(1)
  end do !do j=1,nPlPar
  
  if(quality.eq.0) then
     call pgsubp(1,1)
     call pgsvp(0.,1.,0.,1.)
     call pgswin(-1.,1.,-1.,1.)
     
     call pgsch(sch*0.8)
     call pgmtxt('T',-0.7,0.5,0.5,trim(outputname))  !Print title
     call pgsch(sch)
  end if
  
  call pgend()
  if(file.ge.2) then
     if(file.eq.3) then
        status = system('eps2pdf '//trim(tempfile)//'.eps  -o '//trim(tempfile)//'.pdf   >& /dev/null')
        if(status.ne.0) then
           write(stdErr,'(A,I6)')'  Error converting plot eps -> pdf',status
        else
           if(rmeps) status = system('rm -f '//trim(tempfile)//'.eps  >& /dev/null')
        end if
     end if
  else if(file.eq.1) then
     inquire(file=trim(tempfile)//'.ppm', exist=ex)
     if(ex) then
        convopts = '-resize '//trim(bmpxpix)//' -depth 8 -unsharp '//trim(unSharpchain)
        status = system('convert '//trim(convopts)//' '//trim(tempfile)//'.ppm  '//trim(tempfile)//'.png')
        if(status.ne.0) write(stdErr,'(A,I6)')'  Error converting plot ppm -> png',status
        status = system('rm -f '//trim(tempfile)//'.ppm')
     end if
  end if
  
end subroutine plot_par_L
!***********************************************************************************************************************************
  
  


!***********************************************************************************************************************************
!> \brief  Plot jump sizes
!!
!! \retval exitcode  Exit status code (0=ok)
!!

subroutine plot_Jump_sizes(exitcode)
  use SUFR_constants, only: stdOut,stdErr
  
  use aM_constants, only: rmeps
  use analysemcmc_settings, only: update,prProgress,file,scrsz,scrrat,pssz,psrat,colour,whitebg,quality
  use analysemcmc_settings, only: plBurn,chainPlI,fontsize1d,nPlPar,panels,plPars, chainSymbol,plJump, htmlOutput
  use general_data, only: outputname,outputdir,nChains0,Ntot, parNames,pgOrigParns
  use mcmcrun_data, only: revID,parID
  use plot_data, only: psclr,defcolour,bmpsz,bmprat,ncolours,colours,bmpxpix,nsymbols,symbols,unSharpchain
  use chain_data, only: is,jumps,isburn
  
  implicit none
  integer, intent(out) :: exitcode
  
  integer :: i,j,pgopen,symbol,io,ic,p,status,system
  real :: dx,dy,xmin,xmax,ymin,ymax,sch
  character :: tempfile*(199), convopts*(99)
  logical :: ex
  
  
  exitcode = 0
  
  if(htmlOutput.ge.1) then
     write(stdOut,'(A)') '<a name="jumps"></a>'
     write(stdOut,'(A)') '<font size="1"><a href="#top" title="Go to the top of the page">top</a></font>'
     write(stdOut,'(A)') '<h3>Jump sizes:</h3>'
     write(stdOut,'(A)') '<img src="'//trim(outputname)//'__jumps.png" title="Jump sizes">'
  else
     !if(prProgress.ge.1.and.update.eq.0) write(stdOut,'(A)')' Plotting jump sizes...'
     if(prProgress.ge.1.and.update.eq.0) write(stdOut,'(A)',advance="no")' jump sizes, '
  end if
  
  if(file.eq.0) then
     io = pgopen('18/xs')
     call pgpap(scrSz,scrRat)
     
     sch = fontsize1d*1.5
  end if
  if(file.ge.1) then
     tempfile = trim(outputdir)//'/'//trim(outputname)//'__jumps'
     if(file.eq.1) then
        io = pgopen(trim(tempfile)//'.ppm/ppm')
        call pgpap(bmpsz,bmprat)
        
        sch = fontsize1d*1.5
     end if
     if(file.ge.2) then
        io = pgopen(trim(tempfile)//'.eps'//trim(psclr))
        call pgpap(PSsz,PSrat)
        
        sch = fontsize1d*1.2
     end if
  end if
  
  if(io.le.0) then
     write(stdErr,'(A,I4)')'   Error:  Cannot open PGPlot device.  Quitting the programme',io
     exitcode = 1
     return
  end if
  
  if(quality.eq.0) then !Draft
     !sch = sch*1.75
     sch = sch*1.4
     call pgsvp(0.1,0.95,0.06,0.87) !To make room for title
  end if
  
  call pgsch(sch)
  
  call pgsubp(panels(1),panels(2))
  
  ic = 1
  xmin = 0.
  do j=1,nPlPar
     p = revID(plPars(j))
     if(p.eq.0) then
        call report_undefined_parameter(trim(parNames(plPars(j))), 'plot_Jump_sizes')
        cycle
     end if
     
     call pgpage()
     if(j.eq.1) call pginitl(colour,file,whiteBG)
     call pgsch(sch)
     
     xmax = -1.e30
     ymin =  1.e30
     ymax = -1.e30
     do ic=1,nChains0
        xmin = 0.
        xmax = max(xmax,maxval(is(ic,1:Ntot(ic))))
        if(plJump.eq.1) then
           ymin = min(ymin,minval(jumps(ic,p,2:Ntot(ic))))
           ymax = max(ymax,maxval(jumps(ic,p,2:Ntot(ic))))
        end if
        do i=10,Ntot(ic)
           if(plJump.eq.2.and.jumps(ic,p,i).gt.1.e-20) then
              ymin = min(ymin,log10(abs(jumps(ic,p,i))))
              ymax = max(ymax,log10(abs(jumps(ic,p,i))))
           end if
           ymin = -6.
           ymax = 1.
        end do
     end do
     dx = abs(xmax-xmin)*0.01
     dy = abs(ymax-ymin)*0.05
     if(dy.lt.1.e-10) dy = ymin*0.1
     
     call pgswin(xmin-dx,xmax+dx,ymin-dy,ymax+dy)
     if(plJump.eq.1) call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0) !lin
     if(plJump.eq.2) call pgbox('BCNTS',0.0,0,'BCLNTS',0.0,0) !log
     
     do ic=1,nChains0
        
        !call pgsci(mod(ic*2,10))
        call pgsci(defcolour)
        if(nChains0.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
        
        symbol = chainSymbol
        if(chainSymbol.le.-10) symbol = symbols(mod(ic-1,nsymbols)+1)
        
        if(plJump.eq.1) then
           !do i=1,Ntot(ic),chainPlI
           do i=ic,Ntot(ic),chainPlI !Start at ic to reduce overplotting
              call pgpoint(1,is(ic,i),jumps(ic,p,i), symbol)
           end do
        else
           !do i=1,Ntot(ic),chainPlI
           do i=ic,Ntot(ic),chainPlI !Start at ic to reduce overplotting
              call pgpoint(1,is(ic,i),log10(abs(jumps(ic,p,i))+1.e-30), symbol)
           end do
        end if
     end do
     
     ! Mark the end of the burn-in phase:
     call pgsls(2)
     call pgsci(6)
     do ic=1,nChains0
        if(nChains0.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
        if((plBurn.eq.1.or.plBurn.ge.3).and.isburn(ic).lt.is(ic,Ntot(ic))) call pgline(2,(/isburn(ic),isburn(ic)/), &
             (/ymin-dy,ymax+dy/))
     end do
     
     call pgsci(1)
     call pgsls(1)
     !call pgmtxt('T',1.,0.5,0.5,'Jumps: '//trim(pgOrigParns(parID(p))))
     call pgmtxt('T',-1.2,0.05,0.0,trim(pgOrigParns(parID(p))))
  end do
  
  if(quality.eq.0) then
     call pgsubp(1,1)
     call pgsvp(0.,1.,0.,1.)
     call pgswin(-1.,1.,-1.,1.)
     
     call pgsch(sch*0.8)
     call pgmtxt('T',-0.7,0.5,0.5,trim(outputname))  !Print title
     call pgsch(sch)
  end if
  
  call pgend()
  if(file.ge.2) then
     if(file.eq.3) then
        status = system('eps2pdf '//trim(tempfile)//'.eps  -o '//trim(tempfile)//'.pdf   >& /dev/null')
        if(status.ne.0) then
           write(stdErr,'(A,I6)')'  Error converting plot eps -> pdf',status
        else
           if(rmeps) status = system('rm -f '//trim(tempfile)//'.eps  >& /dev/null')
        end if
     end if
  else if(file.eq.1) then
     inquire(file=trim(tempfile)//'.ppm', exist=ex)
     if(ex) then
        convopts = '-resize '//trim(bmpxpix)//' -depth 8 -unsharp '//trim(unSharpchain)
        status = system('convert '//trim(convopts)//' '//trim(tempfile)//'.ppm  '//trim(tempfile)//'.png')
        if(status.ne.0) write(stdErr,'(A,I6)')'  Error converting plot ppm -> png',status
        status = system('rm -f '//trim(tempfile)//'.ppm')
     end if
  end if
end subroutine plot_Jump_sizes
!***********************************************************************************************************************************
  



!***********************************************************************************************************************************
!> \brief  Plot autocorrelations for each parameter
!!
!! \retval exitcode  Exit status code (0=ok)
!!

subroutine plot_Acorr_chains(exitcode)
  use SUFR_constants, only: stdOut,stdErr
  use SUFR_statistics, only: compute_median_sp
  
  use aM_constants, only: rmeps
  use analysemcmc_settings, only: update,prProgress,file,scrsz,scrrat,pssz,psrat,fonttype,colour,whitebg,quality
  use analysemcmc_settings, only: fontsize1d,nPlPar,panels,plPars, chainSymbol,nAcorr, htmlOutput
  use general_data, only: outputname,outputdir,nChains0,parNames,pgParNss
  use mcmcrun_data, only: revID,parID
  use plot_data, only: psclr,defcolour,bmpsz,bmprat,ncolours,colours,bmpxpix,nsymbols,symbols,unSharpchain
  use chain_data, only: acorrs,lAcorrs
  
  implicit none
  integer, intent(out) :: exitcode
  
  integer :: i,j,pgopen,symbol,io,ic,p,status,system
  real :: dx,dy,xmin,xmax,ymin,ymax,sch
  character :: title*(99), tempfile*(199), convopts*(99)
  logical :: ex
  
  
  exitcode = 0
  
  if(htmlOutput.ge.1) then
     write(stdOut,'(A)') '<a name="acorrs"></a>'
     write(stdOut,'(A)') '<font size="1"><a href="#top" title="Go to the top of the page">top</a></font>'
     write(stdOut,'(A)') '<h3>Autocorrelations:</h3>'
     write(stdOut,'(A)') '<img src="'//trim(outputname)//'__acorrs.png" title="Autocorrelations">'
  else
     !if(prProgress.ge.1.and.update.eq.0) write(stdOut,'(A)')' Plotting autocorrelations...'
     if(prProgress.ge.1.and.update.eq.0) write(stdOut,'(A)',advance="no")' autocorrelations, '
  end if
  
  if(file.eq.0) then
     io = pgopen('19/xs')
     call pgpap(scrSz,scrRat)
     
     sch = fontsize1d*1.5
     sch = sch*1.5
  end if
  tempfile = trim(outputdir)//'/'//trim(outputname)//'__acorrs'
  if(file.eq.1) then
     io = pgopen(trim(tempfile)//'.ppm/ppm')
     call pgpap(bmpsz,bmprat)
     sch = fontsize1d*1.2
     sch = sch*1.5
  end if
  if(file.ge.2) then
     io = pgopen(trim(tempfile)//'.eps'//trim(psclr))
     call pgpap(PSsz,PSrat)
     
     sch = fontsize1d*1.2
     sch = sch*1.5
     call pgscf(fonttype)
  end if
  if(io.le.0) then
     write(stdErr,'(A,I4)')'   Error:  Cannot open PGPlot device.  Quitting the programme',io
     exitcode = 1
     return
  end if
  
  call pgsch(sch)
  
  if(quality.eq.0) call pgsvp(0.08,0.95,0.06,0.87) !To make room for title
  
  call pgsubp(panels(1),panels(2))
  
  ic = 1
  do j=1,nPlPar
     p = revID(plPars(j))
     if(p.eq.0) then
        call report_undefined_parameter(trim(parNames(plPars(j))), 'plot_Acor_chains')
        cycle
     end if
     
     call pgpage()
     if(j.eq.1) call pginitl(colour,file,whiteBG)
     call pgsch(sch)
     
     xmin = 0.
     xmin = minval(acorrs(1:nChains0,0,0:nAcorr))
     xmax = maxval(acorrs(1:nChains0,0,0:nAcorr))
     dx = abs(xmax-xmin)*0.01
     
     !ymin =  1.e30
     !ymax = -1.e30
     !do ic=1,nChains0
     !   ymin = min(ymin,minval(acorrs(ic,p,0:nAcorr)))
     !   ymax = max(ymax,maxval(acorrs(ic,p,0:nAcorr)))
     !end do
     !ymin = max(ymin,-1.)
     !ymax = min(ymax,1.)
     ymin = -1.
     ymax = 1.
     dy = abs(ymax-ymin)*0.05
     !write(stdOut,'(3I3,5F12.2)')p,nChains,nChains0,xmin,xmax,ymin,ymax,acorrs(1,0,nAcorr)
     
     call pgswin(xmin-dx,xmax+dx,ymin-dy,ymax+dy)
     call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)
     
     call pgsci(defcolour)
     do ic=1,nChains0
        
        if(nChains0.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
        symbol = chainSymbol
        if(chainSymbol.le.-10) symbol = symbols(mod(ic-1,nsymbols)+1)
        
        do i=1,nAcorr
           call pgpoint(1,acorrs(ic,0,i),acorrs(ic,p,i),symbol)
        end do
     end do
     
     ! Plot autocorrelation length:
     call pgsls(2)
     call pgsci(defcolour)
     do ic=1,nChains0
        if(nChains0.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
        call pgline(2,(/lAcorrs(ic,p),lAcorrs(ic,p)/),(/ymin-dy,ymax+dy/))
     end do
     
     ! Plot horizontal line at 0:
     call pgsci(1)
     call pgline(2,(/xmin,xmax/),(/0.,0./))
     call pgsci(1)
     call pgsls(1)
     !write(title,'(A,ES9.2)')'Autocorr.: '//trim(pgParNss(parID(p)))//', mean length:',sum(lAcorrs(1:nChains0,p))/real(nChains0)
     write(title,'(A,ES9.2)')'Autocorr.: '//trim(pgParNss(parID(p)))//', med. length:', &
          compute_median_sp(lAcorrs(1:nChains0,p))
     call pgmtxt('T',1.,0.5,0.5,trim(title))
  end do
  
  if(quality.eq.0) then
     call pgsubp(1,1)
     call pgsvp(0.,1.,0.,1.)
     call pgswin(-1.,1.,-1.,1.)
     
     call pgsch(sch*0.8)
     call pgmtxt('T',-0.7,0.5,0.5,trim(outputname))  !Print title
     call pgsch(sch)
  end if
  
  call pgend()
  if(file.ge.2) then
     if(file.eq.3) then
        status = system('eps2pdf '//trim(tempfile)//'.eps  -o '//trim(tempfile)//'.pdf   >& /dev/null')
        if(status.ne.0) then
           write(stdErr,'(A,I6)')'  Error converting plot eps -> pdf',status
        else
           if(rmeps) status = system('rm -f '//trim(tempfile)//'.eps  >& /dev/null')
        end if
     end if
  else if(file.eq.1) then
     inquire(file=trim(tempfile)//'.ppm', exist=ex)
     if(ex) then
        convopts = '-resize '//trim(bmpxpix)//' -depth 8 -unsharp '//trim(unSharpchain)
        status = system('convert '//trim(convopts)//' '//trim(tempfile)//'.ppm  '//trim(tempfile)//'.png')
        if(status.ne.0) write(stdErr,'(A,I6)')'  Error converting plot ppm -> png',status
        status = system('rm -f '//trim(tempfile)//'.ppm')
     end if
  end if
  
end subroutine plot_Acorr_chains
!***********************************************************************************************************************************
  
  







!***********************************************************************************************************************************
!> \brief  Plot R-hat vs. iteration number
!!
!! \retval exitcode  Exit status code (0=ok)
!!

subroutine plot_Rhat_chains(exitcode)
  use SUFR_constants, only: stdOut,stdErr
  
  use aM_constants, only: use_PLplot, rmeps
  use analysemcmc_settings, only: update,prProgress,file,scrsz,scrrat,pssz,psrat,fonttype,colour,whitebg,quality
  use analysemcmc_settings, only: fontsize1d,nPlPar, chainSymbol, htmlOutput
  use general_data, only: outputname,outputdir
  use chain_data, only: Rhats,RhatsN
  use plot_data, only: psclr,defcolour,bmpsz,bmprat,bmpxpix,unSharpchain
  
  implicit none
  integer, intent(out) :: exitcode
  
  integer :: pgopen,symbol,io,status,system, lw
  real :: dx,dy,xmin,xmax,ymin,ymax,sch
  character :: tempfile*(199), convopts*(99)
  logical :: ex
  
  
  exitcode = 0
  
  if(htmlOutput.ge.1) then
     write(stdOut,'(A)') '<a name="rhat"></a>'
     write(stdOut,'(A)') '<font size="1"><a href="#top" title="Go to the top of the page">top</a></font>'
     write(stdOut,'(A)') '<h3>R-hat:</h3>'
     write(stdOut,'(A)') '<img src="'//trim(outputname)//'__rhat.png" title="R-hat">'
  else
     if(prProgress.ge.1.and.update.eq.0) write(stdOut,'(A)',advance="no")' R-hat, '
  end if
  
  if(file.eq.0) then
     io = pgopen('12/xs')
     sch = 1.5*fontsize1D
     lw = 1
  else  ! file.gt.0
     tempfile = trim(outputdir)//'/'//trim(outputname)//'__rhat'
     if(file.eq.1) then
        io = pgopen(trim(tempfile)//'.ppm/ppm')
     else  ! file.ge.2
        io = pgopen(trim(tempfile)//'.eps'//trim(psclr))
     end if
     sch = 1.2*fontsize1D
     lw = 1
     if(file.ge.2) then
        lw = lw*2
        if(nPlPar.eq.1) lw = nint(2*fontsize1d)
     end if
  end if
  
  if(io.le.0) then
     write(stdErr,'(A,I4)')'   Error:  Cannot open PGPlot device.  Quitting the programme',io
     exitcode = 1
     return
  end if
  
  if(file.eq.0) then
     call pgpap(scrSz,scrRat)
     sch = 1.5*fontsize1D
  end if
  if(file.eq.1) then
     call pgpap(bmpsz,bmprat)
     sch = 1.5*fontsize1D
  end if
  if(file.ge.2) then
     call pgpap(PSsz,PSrat)
     call pgscf(fonttype)
     sch = fontsize1D
  end if
  
  
  
  
  call pginitl(colour,file,whiteBG)
  call pgsch(sch)
  call pgslw(lw)
  
  call pgsvp(0.08*fontsize1D,1.-0.03*fontsize1D,0.07*fontsize1D,1.-0.06*fontsize1D)
  if(quality.eq.0) call pgsvp(0.08,0.97,0.07,0.94)  ! To make room for title
  
  
  xmin = 0.
  xmax = maxval(Rhats(1,1:RhatsN))
  dx = abs(xmax-xmin)*0.01
  
  ymin = log10(minval(abs(Rhats(2,1:RhatsN) - 1.0)))
  ymax = min(log10(maxval(Rhats(2,1:RhatsN) - 1.0)),1.0)
  !ymin = log10(minval(Rhats(2,1:RhatsN)))
  !ymax = min(log10(maxval(Rhats(2,1:RhatsN))),1.0)
  dy = abs(ymax-ymin)*0.05
  
  call pgswin(xmin-dx,xmax+dx,ymin-dy,ymax+dy)
  call pgbox('BCNTS',0.0,0,'BCLNTS',0.0,0)
  
  
  call pgsci(defcolour)
  symbol = chainSymbol
  call pgpoint(RhatsN, Rhats(1,1:RhatsN), log10(abs(Rhats(2,1:RhatsN)-1.0)), symbol)
  !call pgpoint(RhatsN, Rhats(1,1:RhatsN), log10(abs(Rhats(2,1:RhatsN))), symbol)
  
  
  call pgsci(1)
  call pgsls(1)
  if(use_PLplot) then
     call pgmtxt('L',5.0,0.5,0.5,'R-hat - 1.0')
     !call pgmtxt('L',5.0,0.5,0.5,'R-hat')
  else
     call pgmtxt('L',2.3,0.5,0.5,'R-hat - 1.0')
     !call pgmtxt('L',2.3,0.5,0.5,'R-hat')
  end if
  if(nPlPar.eq.1.and.quality.eq.1) call pgmtxt('B',2.5,0.5,0.5,'iteration')
  
  if(quality.eq.0) call pgmtxt('T',0.5,0.9,0.9,trim(outputname))  ! Print title
  
  call pgend()
  
  if(file.ge.2) then
     if(file.eq.3) then
        status = system('eps2pdf '//trim(tempfile)//'.eps  -o '//trim(tempfile)//'.pdf   >& /dev/null')
        if(status.ne.0) then
           write(stdErr,'(A,I6)')'  Error converting plot eps -> pdf',status
        else
           if(rmeps) status = system('rm -f '//trim(tempfile)//'.eps  >& /dev/null')
        end if
     end if
  else if(file.eq.1) then
     inquire(file=trim(tempfile)//'.ppm', exist=ex)
     if(ex) then
        convopts = '-resize '//trim(bmpxpix)//' -depth 8 -unsharp '//trim(unSharpchain)
        status = system('convert '//trim(convopts)//' '//trim(tempfile)//'.ppm  '//trim(tempfile)//'.png')
        if(status.ne.0) write(stdErr,'(A,I6)')'  Error converting plot ppm -> png',status
        status = system('rm -f '//trim(tempfile)//'.ppm')
     end if
  end if
  
end subroutine plot_Rhat_chains
!***********************************************************************************************************************************
  
  







!***********************************************************************************************************************************
!> \brief  Plot the vertical axes for the log Posterior plot: logP on the left, sqrt(2logP)~SNR on the right
!!
!! \param itermin  Minimum iteration number in plot (horizontal axis)
!! \param itermax  Maximum iteration number in plot (horizontal axis)
!! \param logpmin  Minimum log Posterior in plot (vertical axis)
!! \param logpmax  Maximum log Posterior in plot (vertical axis)

subroutine plot_posterior_snr_axes(itermin,itermax,logpmin,logpmax)
  use aM_constants, only: use_PLplot
  
  implicit none
  real, intent(in) :: itermin,itermax,logpmin,logpmax
  integer :: i,imin,imax, tick_omi, n
  real :: logpmin0,snrmin,snrmax,dsnr,snr, dlogp,logp, ntick0,tick_om,dtick, tick
  real :: len, dist, ori
  character :: label*(19),fmt*(19)
  
  
  call pgbox('BCNTS',0.0,0,'BNTS',0.0,0)  ! No right border
  call pgbox('',0.0,0,'C',0.0,0)          ! Right border without ticks, labels
  
  dlogp = itermin                         ! Remove 'unused variable' compiler warning
  logpmin0 = max(logpmin,0.)
  dlogp = logpmax - logpmin0
  snrmin = sqrt(max(logpmin0*2,0.))
  snrmax = sqrt(logpmax*2)
  dsnr = snrmax-snrmin
  
  
  ntick0 = 4.                                          ! ~ desired number of ticks
  tick_om = 10.**(floor(log10(abs(dsnr/ntick0))))      ! Order of magnitude of distance between ticks
  tick_omi = nint(log10(tick_om))
  dtick = real(nint(dsnr/ntick0 / tick_om)) * tick_om  ! Distance between ticks
  
  
  ! Format for snr label:
  n = abs(tick_omi)
  if(tick_omi.lt.0) then
     write(fmt,'(A2,I1,A1,I1,A1)')'(F',n+3,'.',n,')'
  else if(tick_omi.le.6) then
     write(fmt,'(A2,I1,A1)')'(I',tick_omi,')'
  else
     write(fmt,'(A)')'(ES7.1)'
  end if
  
  
  ! Print ticks and labels:
  len  = 0.5  ! Tick length
  dist = 0.3  ! Distance axis - label
  if(use_PLplot) dist = 2.0
  ori  = 0.0  ! Orientation (0-360 degrees)
  
  imin = floor(snrmin/dtick)
  imax = ceiling(snrmax/dtick)
  do i=imin,imax
     snr = real(i)*dtick
     logp = (snr**2)/2.
     tick = (logp - logpmin0)/dlogp
     if(tick.lt.0..or.tick.gt.1.) cycle
     
     if(tick_omi.lt.0) then
        write(label,trim(fmt)) snr
     else
        write(label,trim(fmt)) nint(snr)
     end if
     
     call pgtick(itermax,logpmin0,itermax,logpmax,  tick,  len, 0.0,  dist, ori, trim(label))
  end do
  
  
  ! Print sub-ticks:
  len  = 0.25  ! Tick length
  
  imin = floor(snrmin/tick_om)
  imax = ceiling(snrmax/tick_om)
  do i=imin,imax
     snr = real(i)*tick_om
     logp = (snr**2)/2.
     tick = (logp - logpmin0)/dlogp
     if(tick.lt.0..or.tick.gt.1.) cycle
     
     call pgtick(itermax,logpmin0,itermax,logpmax,  tick,  len, 0.0,  dist, ori, '')
  end do
  
end subroutine plot_posterior_snr_axes
!***********************************************************************************************************************************

