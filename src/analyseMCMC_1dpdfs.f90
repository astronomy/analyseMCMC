!> \file analyseMCMC_1dpdfs.f90  Routines to compute and plot one-dimensional PDFs

! 
! LICENCE:
! 
! Copyright 2007-2011 Marc van der Sluys
!  
! This file is part of the AnalyseMCMC package.
!  
! This is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published
! by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! 
! This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License along with this code (LICENSE).  If not, see 
! <http://www.gnu.org/licenses/>.
! 



!***********************************************************************************************************************************
!> \brief  Plot 1D marginalised PDFs
!!
!! \retval exitcode  Exit code: 0=ok

subroutine pdfs1d(exitcode)
  use SUFR_constants, only: stdOut,stdErr, rpi2
  use aM_constants, only: use_PLplot
  use analysemcmc_settings, only: update,prProgress,file,scrsz,scrrat,pssz,psrat,fonttype,colour,whitebg,quality
  use analysemcmc_settings, only: plLmax,fontsize1d,nPlPar,panels,plPars,changeVar
  use analysemcmc_settings, only: plInject,mergeChains,maxChs
  use analysemcmc_settings, only: savePDF,plot,ivals,Nbin1D,Nival,fillPDF,normPDF1D,smooth,plPDF1D,plMedian,plRange,prIval,prValues
  use general_data, only: allDat,outputname,outputdir,nChains0,startval,icloglmax,iloglmax,nChains,parNames,pgParNs
  use general_data, only: pgParNss,selDat,stats,ranges,c0,contrchain,n,maxIter,wrap,fixedpar,shifts,shIvals,pgunits
  use mcmcrun_data, only: totpts,revID,parID,Tchain
  use plot_data, only: psclr,bmpsz,bmprat,ncolours,colours,bmpxpix,unSharppdf1d
  
  implicit none
  integer, intent(out) :: exitcode
  
  integer :: i,j,p,ic,io,pgopen,lw,status,system
  real :: rev24,rev360,rev180
  real :: x(maxChs,maxChs*maxIter),xmin,xmax,xmin1,xmax1,xpeak,dx,ymin,ymax,sch
  
  ! These depend on Nbin1D, allocate after reading input file:
  real,allocatable :: xbin(:,:),ybin(:,:),xbin1(:),ybin1(:),ysum(:),yconv(:),ycum(:)  
  real :: plshift,plx,ply,x0,norm,bindx
  character :: string*(99),str*(99),str1*(99),str2*(99),delta*(19), tempfile*(199), convopts*(99)
  logical :: ex
  
  
  exitcode=0
  if(prProgress.ge.1.and.plot.eq.0.and.savePDF.eq.1) write(stdOut,'(A)',advance="no")'  Saving 1D pdfs'
  if(prProgress.ge.1.and.plot.eq.1.and.update.eq.0) write(stdOut,'(A)',advance="no")' 1D pdfs'
  
  write(delta,'(A,I3.3,A)')'\(2030)\d',nint(ivals(c0)*100),'%\u'
  if(nint(ivals(c0)*100).lt.100) write(delta,'(A,I2.2,A)')'\(2030)\d',nint(ivals(c0)*100),'%\u'
  
  ! Autodetermine number of bins:
  if(Nbin1D.le.0) then
     call determine_nbin_1d(totpts,Nbin1D)
  else
     Nbin1D = max(Nbin1D,5)
  end if
  
  ! Report number of bins used:
  if(prProgress.ge.2.and.plot.eq.1.and.update.eq.0) then
     if(Nbin1D.lt.100) write(stdOut,'(A2,I2,A8)',advance="no")' (',Nbin1D,' bins), '
     if(Nbin1D.ge.100) write(stdOut,'(A2,I3,A8)',advance="no")' (',Nbin1D,' bins), '
  else
     if(prProgress.ge.1.and.plot.eq.1.and.update.eq.0) write(stdOut,'(A2)',advance="no")', '
  end if
  
  ! Allocate memory:
  allocate(xbin(maxChs,Nbin1D+1),ybin(maxChs,Nbin1D+1),xbin1(Nbin1D+1), &
       ybin1(Nbin1D+1),ysum(Nbin1D+1),yconv(Nbin1D+1),ycum(Nbin1D+1))
  
  if(plot.eq.1) then
     
     if(use_PLplot) then  ! then call pgpap *before* pgopen
        if(file.eq.0) call pgpap(scrSz,scrRat)
        if(file.eq.1) call pgpap(bmpsz,bmprat)
        if(file.ge.2) then
           call pgpap(PSsz,PSrat)
           !if(quality.eq.3.and.nPlPar.eq.12) call pgpap(10.6,0.925)
           if(quality.eq.3.and.nPlPar.eq.12) call pgpap(10.6,0.85)
           !call pgscf(fonttype)
        end if
     end if
     
     if(file.eq.0) then
        io = pgopen('14/xs')
        sch = 1.5*fontsize1d
        lw = nint(1*fontsize1d)
     end if
     if(file.ge.1) then
        tempfile = trim(outputdir)//'/'//trim(outputname)//'__pdfs'
        if(file.eq.1) io = pgopen(trim(tempfile)//'.ppm/ppm')
        if(file.ge.2) io = pgopen(trim(tempfile)//'.eps'//trim(psclr))
        lw = nint(3*fontsize1d)
        if(nPlPar.ge.10) lw = nint(2*fontsize1d)
        if(quality.lt.2) lw = max(lw-1,1)  ! Draft/Paper
        sch = 1.2*fontsize1d !2.5
        if(nchains.eq.1.and.nPlPar.gt.9) sch = 1.2*fontsize1d
        if(quality.eq.0) then  ! Draft
           sch = sch*1.75
           lw = nint(2*fontsize1d)
        end if
        if(quality.eq.1) then  ! Paper
           if(nPlPar.eq.12) then
              sch = sch*1.75
              lw = nint(2*fontsize1d)
           else
              sch = sch*1.25
              lw = nint(1*fontsize1d)
           end if
        end if
        if(quality.eq.2) then  ! Talk
           if(nPlPar.le.12) then
              sch = sch*2
              lw = nint(2*fontsize1d)
           else
              sch = sch*1.5
              lw = nint(1*fontsize1d)
           end if
        end if
        if(quality.eq.3) then  ! Poster
           if(nPlPar.eq.12.and.file.ge.2) then
              sch = sch*2.7
              lw = nint(3*fontsize1d)
           else
              !sch = sch*1.25
              !lw = nint(1*fontsize1d)
              sch = sch*1.5
              lw = nint(2*fontsize1d)
           end if
        end if
        if(quality.eq.4) then  ! Vivien's thesis
           sch = sch*2.5
           lw = nint(2*fontsize1d)
        end if
     end if !if(file.ge.1)
     if(io.le.0) then
        write(stdErr,'(A,I4)')'  ERROR:  Cannot open PGPlot device.  Quitting the programme',io
        exitcode = 1
        return
     end if
     
     if(.not.use_PLplot) then  ! then call pgpap *after* pgopen
        if(file.eq.0) call pgpap(scrSz,scrRat)
        if(file.eq.1) call pgpap(bmpsz,bmprat)
        if(file.ge.2) then
           call pgpap(PSsz,PSrat)
           !if(quality.eq.3.and.nPlPar.eq.12) call pgpap(10.6,0.925)
           if(quality.eq.3.and.nPlPar.eq.12) call pgpap(10.6,0.85)
           call pgscf(fonttype)
        end if
     end if
     
     !call pgscr(3,0.,0.5,0.)
     !call pginitl(colour,file,whiteBG)
     call pgslw(lw)
     call pgsch(sch)
     call pgsfs(fillPDF)
     
     
     !if(quality.eq.0) call pgsvp(0.08,0.95,0.06,0.87)  ! To make room for title
     
     call pgsubp(panels(1),panels(2))
  end if !if(plot.eq.1)
  
  !Save 1D PDF data
  if(savePDF.eq.1) then
     open(unit=30,action='write',form='formatted',status='replace',file=trim(outputdir)//'/'//trim(outputname)//'__pdf1d.dat')
     write(30,'(3I6,T101,A)')nPlPar,nchains,Nbin1D,'Total number of plot parameters, total number of chains, number of bins'
     !write(30,'(3I6,T101,A)')nPlPar-nfixedpar,nchains,Nbin1D,&
     !'Total number of plot parameters, total number of chains, number of bins'
  end if
  
  do j=1,nPlPar
     p = revID(plPars(j))
     if(p.eq.0) then
        write(stdErr,'(/,A)')'  * Warning:  pdfs1d():  parameter '//trim(parNames(plPars(j)))// &
             ' is not defined, check plPars() in the input file.  Skipping...'
        cycle
     end if
     
     if(plot.eq.1) then
        call pgpage()
        if(j.eq.1 .or. use_PLplot) call pginitl(colour,file,whiteBG)
        call pgsch(sch)
     end if
     
     ! Set x-ranges for plotting, bin the data and get y-ranges
     ! Use widest probability range (hopefully ~3-sigma) - doesn't always work well...
     if(1.eq.2) then  ! This can only be used if ranges are computed - isn't this always the case?
        xmin = 1.e30
        xmax = -1.e30
        do ic=1,nchains
           xmin = min(xmin,ranges(ic,Nival,p,1))
           xmax = max(xmax,ranges(ic,Nival,p,2))
           !write(stdOut,'(3I4,6F10.3)')ic,p,Nival,xmin,xmax,ranges(ic,Nival,p,1),ranges(ic,Nival,p,2)
        end do
     end if
     !print*,xmin,huge(xmin)
     !if(xmin.le.-huge(xmin).or.xmin.ge.huge(xmin).or.xmax.le.-huge(xmax).or.xmax.ge.huge(xmax)) then !NaN
     !if(version.eq.2.or.xmin.ne.xmin.or.xmax.ne.xmax) then
     xmin = 1.e30
     xmax = -1.e30
     do ic=1,nchains
        if(mergeChains.eq.0.and.contrchain(ic).eq.0) cycle
        xmin = min(xmin,minval(selDat(ic,p,1:n(ic))))
        xmax = max(xmax,maxval(selDat(ic,p,1:n(ic))))
     end do
     !end if
     dx = xmax - xmin
     !dx = max(xmax - xmin,1.e-30)
     
     do ic=1,nchains
        if(mergeChains.eq.0.and.contrchain(ic).eq.0) cycle
        x(ic,1:n(ic)) = selDat(ic,p,1:n(ic))
        xmin1 = minval(selDat(ic,p,1:n(ic)))
        xmax1 = maxval(selDat(ic,p,1:n(ic)))
        if(wrap(ic,p).eq.0) then !Use the plotting ranges
           !xmin1 = xmin - 0.1*dx
           !xmax1 = xmax + 0.1*dx
        end if
        
        call bindata1d(n(ic),x(ic,1:n(ic)),1,Nbin1D,xmin1,xmax1,xbin1,ybin1)
        
        ! Weigh with likelihood.  I should probably do something like this at the start, to get updated ranges etc.:
        !y(ic,1:n(ic)) = selDat(ic,1,1:n(ic))
        !call bindata1da(n(ic),x(ic,1:n(ic)),y(ic,1:n(ic)),1,Nbin1D,xmin1,xmax1,xbin1,ybin1) !Measure amount of L in each bin
        
        if(parID(p).eq.21.or.parID(p).eq.22.or.parID(p).eq.72.or.parID(p).eq.82.or.parID(p).eq.32.or.parID(p).eq.53) then
           if(ybin1(1).gt.ybin1(2)) ybin1(1)=0.
           if(ybin1(Nbin1D).gt.ybin1(Nbin1D-1)) ybin1(Nbin1D)=0.
        end if
        
        ! Normalise 1D PDF:
        if(normPDF1D.gt.0) then
           if(normPDF1D.eq.1) then ! Normalise the SURFACE, not the height (because of different bin size).  This is the default
              norm = 0.
              do i=1,Nbin1D+1
                 norm = norm + ybin1(i)
              end do
              norm = norm*(xmax1-xmin1)
              ybin1 = ybin1/norm
           else ! Normalise to the height of the PDF
              if(normPDF1D.eq.2) ybin1 = ybin1/maxval(ybin1)  ! Normalise to the height of the PDF
              if(normPDF1D.eq.3) ybin1 = ybin1/(maxval(ybin1)**0.5)  ! Normalise to the sqrt of the height of the PDF;
              !                                                       Works nicely for comparing parallel-tempering chains
              if(normPDF1D.eq.4) ybin1 = (ybin1/maxval(ybin1))**(Tchain(ic)) ! Extrapolate the true PDF from temperature chains
              if(normPDF1D.eq.5) ybin1 = (ybin1/maxval(ybin1))**(2.0)
              if(ic*j.eq.1) write(stdErr,'(//,A,/)')'  *** WARNING: using non-default normalisation for PDFs ***'
           end if
        end if
        
        ! Smoothen 1D PDF:
        if(smooth.gt.1) call smoothpdf1d(ybin1,Nbin1D+1,smooth)
        xbin(ic,1:Nbin1D+1) = xbin1(1:Nbin1D+1)
        ybin(ic,1:Nbin1D+1) = ybin1(1:Nbin1D+1)
        
        
        ! Save binned data:
        if(savePDF.eq.1) then
           !if(savePDF.eq.1.and.fixedpar(p).eq.0) then
           if(fixedpar(p).eq.1) ybin1 = 0.  ! Prevent NaNs
           write(30,'(A)')'--------------------------------------------------------------------------------------------------'// &
                '------------------------------------------------------------------------------------------------------'
           write(30,'(3I6,T31,A10,T101,A)')ic,parID(p),wrap(ic,p),parNames(parID(p)), &
                'Chain number, parameter ID, wrap (1/0 = y/n) and parameter name  (ic, parID(), wrap(), parNames())'
           write(30,'(2ES15.7,T101,A)')startval(ic,p,1:2),'Injection and starting value  (startval(1:2)'
           write(30,'(6ES15.7,T101,A)')stats(ic,p,1:6),'Stats: median, mean, absVar1, absVar2, stdev1, stdev2  (stats(1:6))'
           write(30,'(5ES15.7,T101,A)')ranges(ic,c0,p,1:5),'Ranges: lower,upper limit, centre, width, relative width  '// &
                '(ranges(1:5))'
           write(30,'(2ES15.7,T101,A)')xmin1,xmax1,'Xmin and Xmax of PDF  (xmin1,xmax1)'
           
           ! Bin contents:
           write(30,'(2ES15.7,T101,A,I4,A)')xbin1(1),ybin1(1),'The X and Y values of the',Nbin1D,' bins  (xbin1,ybin1)'
           do i=2,Nbin1D+1
              write(30,'(2ES15.7)')xbin1(i),ybin1(i)
           end do
        end if
     end do !ic
     
     if(plot.eq.1) then
        ! Ranges for plot panel:
        xmin = xmin - 0.1*dx
        xmax = xmax + 0.1*dx
        ymin = 0.
        ymax = tiny(ymax)
        xpeak = 0.
        !print*,xmin,xmax,ymin,ymax
        do ic=1,nchains
           if(mergeChains.eq.0.and.contrchain(ic).eq.0) cycle
           !ymax = max(ymax,maxval(ybin(ic,1:Nbin1D+1)))
           do i=1,Nbin1D+1
              !print*,ybin(ic,i),ymax
              if(ybin(ic,i).gt.ymax) then
                 ymax = ybin(ic,i)
                 xpeak = xbin(ic,i)
              end if
           end do
        end do
        ymax = ymax*1.1
        if(dx.eq.0) then
           xmin = 0.5*xmin
           xmax = 2*xmax
           if(xmin.eq.0.) then
              xmin = -1.
              xmax = 1.
           end if
        end if
        if(ymax.lt.1.e-19) ymax = 1.
        
        
        if(file.eq.0.and.scrRat.gt.1.35) call pgsvp(0.08,0.95,0.1,0.95)
        if(file.eq.1.and.bmprat.gt.1.35) call pgsvp(0.08,0.95,0.1,0.95)
        if(file.ge.2.and.PSrat.gt.1.35) call pgsvp(0.08,0.95,0.1,0.95)
        if(quality.eq.0) call pgsvp(0.08,0.95,0.06,0.87)  ! To make room for title
        if(quality.eq.4) call pgsvp(0.13,0.95,0.1,0.95)
        
        call pgsch(sch)
        call pgswin(xmin,xmax,ymin,ymax)
        if(abs(dx).lt.1.e-30) then  ! So that the program doesn't hang if a parameter is kept constant
           xbin = 0.
           ybin = 0.
        end if
        
        ! Plot 1D PDF:
        !if(fixedpar(p).eq.0) then
        call pgsci(1)
        bindx = 0.
        if(file.ge.2) call pgslw(lw)
        do ic=1,nchains
           if(mergeChains.eq.0.and.contrchain(ic).eq.0) cycle
           
           ! Set hatch style: angle = +-45deg, phase between 0 and 1 (1/nchains0, 2/nchains0, ...):
           if(fillPDF.ge.3) call pgshs(45.0*(-1)**ic,2.0,real(ic)/real(nchains0)) 
           if(nchains.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
           xbin1(1:Nbin1D+1) = xbin(ic,1:Nbin1D+1)
           ybin1(1:Nbin1D+1) = ybin(ic,1:Nbin1D+1)
           if(wrap(ic,p).eq.0) then
              if(nchains.eq.1) call pgsci(15)
              if(plPDF1D.eq.1) then
                 call pgpoly(Nbin1D+2,(/xbin1(1),xbin1(1:Nbin1D+1)/)+bindx/2.,(/0.,ybin1(1:Nbin1D+1)/))
              else
                 call verthist(Nbin1D+2,(/xbin1(1),xbin1(1:Nbin1D+1)/),(/0.,ybin1(1:Nbin1D+1)/),2)
              end if
              
              ! Plot pdf contour:
              !if(nchains.eq.1) call pgsci(1)
              !call pgsci(1)
              if(fillPDF.eq.1) call pgsci(1)
              if(nchains.eq.1) call pgsci(2)
              if(plPDF1D.eq.1) then
                 call pgline(Nbin1D+1,xbin1(1:Nbin1D+1)+bindx/2.,ybin1(1:Nbin1D+1)) !:Nbin1D) ?
              else
                 call verthist(Nbin1D+2,(/xbin1(1),xbin1(1:Nbin1D+1)/),(/0.,ybin1(1:Nbin1D+1)/),1)
              end if
              
              ! Fix the loose ends:
              call pgline(2,(/xbin1(1)-bindx,xbin1(1)/)+bindx/2.,(/0.,ybin1(1)/))
              call pgline(2,(/xbin1(Nbin1D+1),xbin1(Nbin1D+1)/)+bindx/2.,(/ybin1(Nbin1D+1),0./))
           else  ! If parameter is wrapped
              plshift = rpi2  !2pi
              if(changeVar.ge.1) then
                 plshift = 360.
                 if(parID(p).eq.31) plshift = 24.   ! RA in hours
                 if(parID(p).eq.52) plshift = 180.  ! Pol.angle
              end if
              if(nchains.eq.1) call pgsci(15)
              if(plPDF1D.eq.1) then
                 call pgpoly(Nbin1D+3,(/xbin1(1),xbin1(1:Nbin1D),xbin1(1)+plshift,xbin1(1)+plshift/),(/0.,ybin1(1:Nbin1D), &
                      ybin1(1),0./))
              else
                 call verthist(Nbin1D+3,(/xbin1(1),xbin1(1:Nbin1D),xbin1(1)+plshift,xbin1(1)+plshift/),(/0.,ybin1(1:Nbin1D), &
                      ybin1(1),0./),2)
              end if
              
              ! Plot pdf contour:
              !call pgsci(1)
              !if(fillPDF.ne.1) call pgsci(colours(mod(ic-1,ncolours)+1))
              if(fillPDF.eq.1) call pgsci(1)
              if(nchains.eq.1) call pgsci(2)
              if(plPDF1D.eq.1) then
                 call pgline(Nbin1D,xbin1(1:Nbin1D)+bindx/2.,ybin1(1:Nbin1D))
              else
                 call verthist(Nbin1D,xbin1(1:Nbin1D)+bindx/2.,ybin1(1:Nbin1D),0)
              end if
              
              ! Plot dotted lines outside the pdf for wrapped periodic parameters:
              call pgsls(4)
              !if(file.ge.2) call pgslw(2)
              if(plPDF1D.eq.1) then
                 call pgline(Nbin1D+1,(/xbin1(1:Nbin1D)-plshift,xbin1(1)/)+bindx/2.,(/ybin1(1:Nbin1D),ybin1(1)/))
                 call pgline(Nbin1D,xbin1+plshift+bindx/2.,ybin1)
              else
                 call verthist(Nbin1D+1,(/xbin1(1:Nbin1D)-plshift,xbin1(1)/)+bindx/2.,(/ybin1(1:Nbin1D),ybin1(1)/),0)
                 call verthist(Nbin1D,xbin1+plshift+bindx/2.,ybin1,0)
              end if
              
              ! Fix the loose end:
              call pgsls(1)
              if(file.ge.2) call pgslw(lw)
              if(plPDF1D.eq.1) then
                 call pgline(2,(/xbin1(Nbin1D),xbin1(1)+plshift/)+bindx/2.,(/ybin1(Nbin1D),ybin1(1)/))
              else
                 call verthist(2,(/xbin1(Nbin1D),xbin1(1)+plshift/)+bindx/2.,(/ybin1(Nbin1D),ybin1(1)/),0)
              end if
           end if
        end do !ic
        
        
        ! Plot lines again over surface of overlapping distributions:
        if(nchains.gt.1.and.fillPDF.eq.1) then
           call pgsls(4)
           do ic=1,nchains
              if(mergeChains.eq.0.and.contrchain(ic).eq.0) cycle
              call pgsci(1)
              !call pgsci(colours(mod(ic-1,ncolours)+1))
              xbin1(1:Nbin1D+1) = xbin(ic,1:Nbin1D+1)
              ybin1(1:Nbin1D+1) = ybin(ic,1:Nbin1D+1)
              if(wrap(ic,p).eq.0) then
                 if(plPDF1D.eq.1) then
                    call pgline(Nbin1D+1,xbin1(1:Nbin1D+1)+bindx/2.,ybin1(1:Nbin1D+1))
                 else
                    call verthist(Nbin1D+1,xbin1(1:Nbin1D+1)+bindx/2.,ybin1(1:Nbin1D+1),0)
                 end if
              else
                 if(plPDF1D.eq.1) then
                    call pgline(Nbin1D,xbin1(1:Nbin1D)+bindx/2.,ybin1(1:Nbin1D))
                 else
                    call verthist(Nbin1D,xbin1(1:Nbin1D)+bindx/2.,ybin1(1:Nbin1D),0)
                 end if
              end if
           end do
           call pgsls(4)
        end if
        
        
        ! Plot max likelihood:
        if(plLmax.ge.1) then
           ply = allDat(icloglmax,p,iloglmax)
           if(parID(p).eq.31) ply = rev24(ply)
           if(parID(p).eq.41.or.parID(p).eq.54.or.parID(p).eq.73.or.parID(p).eq.83) ply = rev360(ply)
           if(parID(p).eq.52) ply = rev180(ply)
           call pgsci(1)
           call pgsls(5)
           call pgline(2,(/ply,ply/),(/ymin,ymax/))
        end if
        
        
        ! Plot median and model value:
        call pgsch(sch)
        
        do ic=1,nchains
           if(mergeChains.eq.0.and.contrchain(ic).eq.0) cycle
           
           ! Draw white lines:
           if(nchains.gt.1) then
              call pgslw(lw)
              call pgsls(1); call pgsci(0)
              ! Injection value:
              if(plInject.eq.1.or.plInject.eq.3) call pgline(2,(/startval(ic,p,1),startval(ic,p,1)/),(/ymin,ymax/))
              
              ! Injection value - mass, t_c and spin only (?) - CHECK put in input file?
              if((plInject.eq.2.or.plInject.eq.4).and.(parID(p).eq.11.or.parID(p).eq.12 &
                   .or.parID(p).eq.61.or.parID(p).eq.62.or.parID(p).eq.63.or.parID(p).eq.64.or. &
                   parID(p).eq.71.or.parID(p).eq.81)) call pgline(2,(/startval(ic,p,1),startval(ic,p,1)/),(/ymin,ymax/))  
              
              !if(plStart.ge.1) call pgline(2,(/startval(ic,p,2),startval(ic,p,2)/),(/ymin,ymax/))              ! Starting value
              
              ! Median:
              if(plMedian.eq.1.or.plMedian.eq.3.or.plMedian.eq.4.or.plMedian.eq.6)  &
                   call pgline(2,(/stats(ic,p,1),stats(ic,p,1)/),(/ymin,ymax/))
              
              ! Ranges:
              if(plRange.eq.1.or.plRange.eq.3.or.plRange.eq.4.or.plRange.eq.6) then
                 if(nchains.lt.2) call pgline(2,(/ranges(ic,c0,p,1),ranges(ic,c0,p,1)/),(/ymin,ymax/)) ! Left limit of interval
                 if(nchains.lt.2) call pgline(2,(/ranges(ic,c0,p,2),ranges(ic,c0,p,2)/),(/ymin,ymax/)) ! Right limit of interval
                 if(nchains.eq.1) call pgline(2,(/ranges(ic,c0,p,3),ranges(ic,c0,p,3)/),(/ymin,ymax/)) ! Centre of interval
              end if
           end if
           
           call pgslw(lw+1)
           
           ! Draw coloured lines over the white ones:
           ! Median:
           if((plMedian.eq.1.or.plMedian.eq.3.or.plMedian.eq.4.or.plMedian.eq.6)) then
              call pgsls(2); call pgsci(2); if(nchains.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
              call pgline(2,(/stats(ic,p,1),stats(ic,p,1)/),(/ymin,ymax/))
           end if
           
           ! Plot ranges in 1D PDF:
           if((plRange.eq.1.or.plRange.eq.3.or.plRange.eq.4.or.plRange.eq.6)) then
              call pgsls(4); call pgsci(2); if(nchains.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
              if(nchains.lt.2) call pgline(2,(/ranges(ic,c0,p,1),ranges(ic,c0,p,1)/),(/ymin,ymax/))   ! Left limit of interval
              if(nchains.lt.2) call pgline(2,(/ranges(ic,c0,p,2),ranges(ic,c0,p,2)/),(/ymin,ymax/))   ! Right limit of interval
              !if(nchains.eq.1) call pgline(2,(/ranges(ic,c0,p,3),ranges(ic,c0,p,3)/),(/ymin,ymax/))  ! Centre of interval
           end if
           
           ! Plot injection value in PDF:
           !if(plInject.ge.1) then !Plot injection values
           if(plInject.eq.1.or.plInject.eq.3.or.((plInject.eq.2.or.plInject.eq.4).and.(parID(p).eq.11.or.parID(p).eq.12 &
                .or.parID(p).eq.61.or.parID(p).eq.62.or.parID(p).eq.63.or.parID(p).eq.64.or.parID(p).eq.71.or.parID(p).eq.81))) then
              if(mergeChains.ne.1.or.ic.le.1) then 
                 !The units of the injection values haven't changed (e.g. from rad to deg) for ic>1 
                 !(but they have for the starting values, why?)
                 call pgsls(2); call pgsci(1)
                 
                 ! Dash-dotted line for injection value when Lmax line isn't plotted (should we do this always?):
                 if(plLmax.eq.0) call pgsls(3)  
                 plx = startval(ic,p,1)
                 if(wrap(ic,p).ne.0) plx = mod(plx + shifts(ic,p), shIvals(ic,p)) - shifts(ic,p)
                 call pgline(2,(/plx,plx/),(/ymin,ymax/)) !Injection value
              end if
           end if
           
           ! Plot starting value in 1D PDF:
           !if(plStart.eq.1.and.abs((startval(ic,p,1)-startval(ic,p,2))/startval(ic,p,1)).gt.1.e-10) then
           !   call pgsls(4); call pgsci(1); if(nchains.gt.1) call pgsci(1)
           !   call pgline(2,(/startval(ic,p,2),startval(ic,p,2)/),(/ymin,ymax/))
           !end if
           
           call pgsls(1)
           call pgsci(1)
        end do !ic
        
        
        
        
        ! Print median, model value and range widths in 1D PDF panel title:
        call pgslw(lw)
        call pgsci(1)
        ic = 1
        !if(nPlPar.lt.7.or.nPlPar.eq.9) then  !Three or less columns
        if(quality.ne.2.and.quality.ne.3.and.quality.ne.4) then  !Not a talk/poster/thesis
           if(nPlPar.le.5) then
              write(str,'(A,F7.3,A5,F7.3)')trim(pgParNs(parID(p)))//': mdl:',startval(ic,p,1),' med:',stats(ic,p,1)
              if(prIval.ge.1) then
                 if(plRange.eq.4.or.plRange.eq.5.or.plRange.eq.6) then
                    
                    ! Distance, Mc, M1, M2:
                    if(parID(p).eq.21.or.parID(p).eq.22.or.parID(p).eq.61.or.parID(p).eq.63.or.parID(p).eq.64) then  
                       write(str,'(A,F6.2,A1)')trim(str)//' '//trim(delta)//':',ranges(ic,c0,p,5)*100,'%'
                    else
                       write(str,'(A,F7.3)')trim(str)//' '//trim(delta)//':',ranges(ic,c0,p,5)
                    end if
                 end if
              end if
           else  !if nPlPar>=5
              str = ' '
              if(parID(p).eq.21.or.parID(p).eq.22.or.parID(p).eq.61.or.parID(p).eq.63.or.parID(p).eq.64) then  !Distance, Mc,M1,M2
                 if(plInject.eq.3.or.plInject.eq.4) write(str,'(A,F7.3)')trim(str)//' mdl:',startval(ic,p,1)
                 if(plMedian.eq.4.or.plMedian.eq.5.or.plMedian.eq.6) write(str,'(A,F8.3)')trim(str)//' med:',stats(ic,p,1)
                 if(prIval.ge.1.and.(plRange.eq.4.or.plRange.eq.5.or.plRange.eq.6))  &
                      write(str,'(A,F6.2,A1)')trim(str)//' '//trim(delta)//':',ranges(ic,c0,p,5)*100,'%'
              else
                 if(plInject.eq.3.or.plInject.eq.4) write(str,'(A,F7.3)')trim(str)//' mdl:',startval(ic,p,1)
                 if(plMedian.eq.4.or.plMedian.eq.5.or.plMedian.eq.6) write(str,'(A,F8.3)')trim(str)//' med:',stats(ic,p,1)
                 if(prIval.ge.1.and.(plRange.eq.4.or.plRange.eq.5.or.plRange.eq.6))  &
                      write(str,'(A,F7.3)')trim(str)//' '//trim(delta)//':',ranges(ic,c0,p,5)
              end if
              call pgsch(sch*1.2)
              call pgptxt(xmin+0.05*dx,ymax*(1.0-0.1*fontsize1d),0.,0.,trim(pgParNss(parID(p))))
           end if
        end if
        
        
        if(quality.eq.2.or.quality.eq.3.or.quality.eq.4) then  ! Talk/poster/thesis quality for 1D PDF
           if(plRange.eq.1.or.plRange.eq.3.or.plRange.eq.4.or.plRange.eq.6) then
              x0 = ranges(ic,c0,p,5)
              if(parID(p).eq.21.or.parID(p).eq.22.or.parID(p).eq.61.or.parID(p).eq.63.or.parID(p).eq.64) x0 = x0*100
              !print*,p,x0,nint(x0)
              if(x0.lt.0.01) write(str,'(F6.4)')x0
              if(x0.ge.0.01.and.x0.lt.0.1) write(str,'(F5.3)')x0
              if(x0.ge.0.1.and.x0.lt.1.) write(str,'(F4.2)')x0
              if(x0.ge.1.and.x0.lt.9.95) write(str,'(F3.1)')x0
              if(x0.ge.9.95.and.x0.lt.99.5) write(str,'(I2)')nint(x0)
              if(x0.ge.99.5) write(str,'(I3)')nint(x0)
              write(str,'(A)')trim(delta)//': '//trim(str)
              if(parID(p).eq.21.or.parID(p).eq.22.or.parID(p).eq.61.or.parID(p).eq.63.or.parID(p).eq.64) then
                 write(str,'(A)')trim(str)//'%'
              else
                 write(str,'(A)')trim(str)//trim(pgunits(parID(p)))
              end if
           else if(plInject.eq.2.or.plInject.eq.4) then  ! If not plotting ranges, but do plot injection values
              x0 = startval(ic,p,1)
              !print*,p,x0,nint(x0)
              if(x0.lt.0.01) write(str,'(F7.4)')x0
              if(x0.ge.0.01.and.x0.lt.0.1) write(str,'(F6.3)')x0
              if(x0.ge.0.1.and.x0.lt.1.) write(str,'(F6.3)')x0
              if(x0.ge.1.and.x0.lt.9.995) write(str,'(F6.3)')x0
              if(x0.ge.9.995.and.x0.lt.99.9) write(str,'(F5.1)')x0
              if(x0.ge.99.9) write(str,'(F6.1)')x0
              write(str,'(A)')'inj: '//trim(str)//trim(pgunits(parID(p)))
           end if
           
           
           ! Print MCMC parameter name in top of panel:
           call pgsch(sch*1.2)
           if(abs(xmin-xpeak).lt.abs(xmax-xpeak)) then  ! Peak is at left, put varname at right
              call pgptxt(xmax-0.05*dx,ymax*(1.-0.1*fontsize1d),0.,1.,trim(pgParNs(parID(p))))
           else
              call pgptxt(xmin+0.05*dx,ymax*(1.-0.1*fontsize1d),0.,0.,trim(pgParNs(parID(p))))
           end if
        end if
        
        
        ! Write the deltas of the two pdfs:
        if(nchains.eq.2..and.(plRange.eq.4.or.plRange.eq.5.or.plRange.eq.6)) then
           write(str,'(A)')trim(delta)
           if(parID(p).eq.21.or.parID(p).eq.22.or.parID(p).eq.61.or.parID(p).eq.63.or.parID(p).eq.64) then  ! Distance, Mc, M1, M2
              write(str1,'(A,F6.2,A1)')trim(delta)//':',ranges(1,c0,p,5)*100.,'%'
              write(str2,'(A,F6.2,A1)')trim(delta)//':',ranges(2,c0,p,5)*100.,'%'
           else
              write(str1,'(A,F7.3)')trim(delta)//':',ranges(1,c0,p,5)
              write(str2,'(A,F7.3)')trim(delta)//':',ranges(2,c0,p,5)
           end if
        end if
        call pgsch(sch*1.1)
        if(prValues.eq.1) then
           if(nchains.eq.2) then
              call pgsci(colours(mod(0,ncolours)+1))
              call pgmtxt('T',0.5,0.25,0.5,trim(str1))
              call pgsci(colours(mod(1,ncolours)+1))
              call pgmtxt('T',0.5,0.75,0.5,trim(str2))
           else
              if(quality.eq.2.or.quality.eq.3) call pgsci(2)
              !call pgptxt(ranges(ic,c0,p,3),ymax,0.,0.5,trim(str))  ! Align with centre of 90%-probability range
              call pgptxt((xmin+xmax)/2.,ymax,0.,0.5,trim(str))      ! Centre
              call pgsci(2)
              if(plRange.eq.1.or.plRange.eq.3.or.plRange.eq.4.or.plRange.eq.6)  &
                   call pgline(2,(/ranges(ic,c0,p,1),ranges(ic,c0,p,2)/),(/0.99*ymax,0.99*ymax/))  ! Plot line at top over P.range
              call pgsci(1)
           end if
        end if
        !else  !If parameter was fixed, plot parameter name in empty panel
        !   call pgsch(sch*1.2)
        !   call pgptxt(xmin+0.05*dx,ymax*0.9,0.,0.,trim(pgParNss(parID(p))))
        !end if !if(fixedpar(p).eq.0)
        
        call pgsci(1)
        call pgsch(sch)
        call pgbox('BNTS',0.0,0,'',0.0,0)
     end if !if(plot.eq.1) 
  end do !p
  
  if(savePDF.eq.1) close(30)
  
  if(plot.eq.1) then
     call pgsubp(1,1)
     
     ! CHECK: PLplot complains, needed for PGPlot?
     !call pgsvp(0.,1.,0.,1.)
     !call pgswin(-1.,1.,-1.,1.)
     
     
     if(quality.eq.0) then
        ! Remove also the pgsvp at the beginning of the plot:
        string=' '
        write(string,'(A,I7,A,I4)')trim(string)//'n:',totpts,', nbin:',Nbin1D
        call pgsch(sch*0.5)
        call pgmtxt('T',-0.7,0.5,0.5,trim(outputname)//'  '//trim(string))  !Print title
        call pgsch(sch)
     end if
     
     call pgend
     
     if(file.ge.2) then
        if(file.eq.3) then
           status = system('eps2pdf '//trim(tempfile)//'.eps -o '//trim(tempfile)//'.pdf  >& /dev/null')
           if(status.ne.0) write(stdErr,'(A,I6)')'  Error converting plot eps - >pdf',status
        end if
     else if(file.eq.1) then
        inquire(file=trim(tempfile)//'.ppm', exist=ex)
        if(ex) then
           convopts = '-resize '//trim(bmpxpix)//' -depth 8 -unsharp '//trim(unSharppdf1d)
           status = system('convert '//trim(convopts)//' '//trim(tempfile)//'.ppm '//trim(tempfile)//'.png')
           if(status.ne.0) write(stdErr,'(A,I6)')'  Error converting plot ppm -> png',status
           status = system('rm -f '//trim(tempfile)//'.ppm')
        end if
     end if
  end if !if(plot.eq.1)
  
  
  ! Deallocate memory:
  deallocate(xbin,ybin,xbin1,ybin1,ysum,yconv,ycum)
  
  
end subroutine pdfs1d
!***********************************************************************************************************************************









!***********************************************************************************************************************************
!> \brief  Bin data 1D by counting the number of points in each bin
!! 
!! \param n      Number of data points
!! \param x      Data to be binned (n points)
!! \param norm   Normalise histogram (1) or not (0)
!! \param nbin   Desired number of bins
!! \param xmin1  Minimum value of the binning range.  Set xmin=xmax to auto-determine (I/O)
!! \param xmax1  Maximum value of the binning range.  Set xmin=xmax to auto-determine (I/O)
!! \retval xbin  Binned data, location of the bins.  The x values are the left side of the bin!
!! \retval ybin  Binned data, height of the bins.

subroutine bindata1d(n,x,norm,nbin,xmin1,xmax1,xbin,ybin)
  implicit none
  integer, intent(in) :: n,nbin,norm
  real, intent(in) :: x(n)
  real, intent(inout) :: xmin1,xmax1
  real, intent(out) :: xbin(nbin+1),ybin(nbin+1)
  integer :: i,k
  real :: xmin,xmax,dx
  
  xmin = xmin1
  xmax = xmax1
  
  if(abs((xmin-xmax)/(xmax+1.e-30)).lt.1.e-20) then  !Autodetermine
     xmin = minval(x(1:n))
     xmax = maxval(x(1:n))
     xmin1 = xmin                                    !And return new values
     xmax1 = xmax
  end if
  dx = abs(xmax - xmin)/real(nbin)
  
  do k=1,nbin+1
     !xbin(k) = xmin + (real(k)-0.5)*dx  !x is the centre of the bin
     xbin(k) = xmin + (k-1)*dx          !x is the left of the bin
  end do
  
  !ybintot=0.
  ybin = 0.
  do i=1,n
     do k=1,nbin
        if(x(i).ge.xbin(k)) then
           if(x(i).lt.xbin(k+1)) then
              ybin(k) = ybin(k) + 1.
              exit !If point i fits in this bin, don't try the others
           end if
        end if
     end do !k (bin)
     !ybintot = ybintot + ybin(k)
  end do
  !if(norm.eq.1) ybin = ybin/(ybintot+1.e-30)
  if(norm.eq.1) ybin = ybin/(sum(ybin)+1.e-30)
  
end subroutine bindata1d
!***********************************************************************************************************************************



!***********************************************************************************************************************************
!> \brief  "Bin data" 1D by measuring the amount of likelihood in each bin
!! 
!! \param n      Number of data points
!! \param x      Data to be binned (n points)
!! \param y      Likelihoods of the n points
!! \param norm   Normalise histogram (1) or not (0)
!! \param nbin   Desired number of bins
!! \param xmin1  Minimum value of the binning range.  Set xmin=xmax to auto-determine (I/O)
!! \param xmax1  Maximum value of the binning range.  Set xmin=xmax to auto-determine (I/O)
!!
!! \retval xbin  Binned data, location of the bins.  The x values are the left side of the bin!
!! \retval ybin  Binned data, height of the bins.

subroutine bindata1da(n,x,y, norm,nbin,xmin1,xmax1, xbin,ybin)
  implicit none
  integer, intent(in) :: n,nbin,norm
  real, intent(in) :: x(n),y(n)
  real, intent(inout) :: xmin1,xmax1
  real, intent(out) :: xbin(nbin+1),ybin(nbin+1)
  integer :: i,k
  real :: xmin,xmax,dx,ymin,ybintot
  
  xmin = xmin1
  xmax = xmax1
  ymin = minval(y)
  !print*,n,nbin,xmin1,xmax1
  !print*,minval(y),maxval(y)
  
  if(abs((xmin-xmax)/(xmax+1.e-30)).lt.1.e-20) then  ! Autodetermine
     xmin = minval(x(1:n))
     xmax = maxval(x(1:n))
  end if
  dx = abs(xmax - xmin)/real(nbin)
  
  do k=1,nbin+1
     !        xbin(k) = xmin + (real(k)-0.5)*dx  !x is the centre< of the bin
     xbin(k) = xmin + (k-1)*dx          !x is the left of the bin
  end do
  ybintot=0.
  do k=1,nbin
     ybin(k) = 0.
     do i=1,n
        !if(x(i).ge.xbin(k).and.x(i).lt.xbin(k+1)) ybin(k) = ybin(k) + 1.
        if(x(i).ge.xbin(k).and.x(i).lt.xbin(k+1)) ybin(k) = ybin(k) + exp(y(i) - ymin)
     end do
     ybintot = ybintot + ybin(k)
  end do
  if(norm.eq.1) ybin = ybin/(ybintot+1.e-30)
  
  if(abs((xmin1-xmax1)/(xmax1+1.e-30)).lt.1.e-20) then   ! Autodetermine
     xmin1 = xmin
     xmax1 = xmax
  end if
  
end subroutine bindata1da
!***********************************************************************************************************************************



!***********************************************************************************************************************************
!> \brief Plot a 1D vertical histogram.  x is the left of the bin!
!!
!! \param n      Number of bins
!! \param x1     Location of the bins
!! \param y1     Height of the bins
!! \param style  Histogram style:  0: don't drop vertical lines to 0 at each bin, 1: do this, 2: fill the histogram

subroutine verthist(n,x1,y1,style)
  use analysemcmc_settings, only: fillPDF
  
  implicit none
  integer, intent(in) :: n,style
  real, intent(in) :: x1(n+1),y1(n+1)
  
  integer :: j
  real :: dx,c, x(n+1),y(n+1)
  
  ! Make columns overlap a little:
  dx = x1(2)-x1(1)
  c = 0.001*dx
  
  x = x1
  y = y1
  
  x(n+1) = x(n) + (x(n)-x(n-1))
  y(n+1) = 0.
  
  ! Don't close line at left and right (don't drop to 0):
  if(style.eq.0) then
     do j=1,n-1
        call pgline(2,x(j:j+1),(/y(j),y(j)/))
        call pgline(2,(/x(j+1),x(j+1)/),y(j:j+1))
     end do
     call pgline(2,x(n:n+1),(/y(n),y(n)/))
  end if
  
  ! Close line at left and right (drop to 0):
  if(style.eq.1) then
     call pgline(2,(/x(1),x(1)/),(/0.,y(1)/))
     do j=1,n
        call pgline(2,x(j:j+1),(/y(j),y(j)/))
        call pgline(2,(/x(j+1),x(j+1)/),y(j:j+1))
     end do
  end if
  
  ! Fill the histogram:
  if(style.eq.2) then
     call pgsfs(fillPDF)
     do j=1,n
        call pgpoly(4,(/x(j)-c,x(j)-c,x(j+1)+c,x(j+1)+c/),(/0.,y(j),y(j),0./))
     end do
  end if
  
end subroutine verthist
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Plot a 1D horizontal histogram
!!
!! \param n  Number of bins
!! \param x  Locations of the bins
!! \param y  Heights of the bins

subroutine horzhist(n,x,y)
  implicit none
  integer, intent(in) :: n
  real, intent(in) :: x(n),y(n)
  integer :: j
  
  call pgline(2,(/x(1),x(1)/),(/0.,y(1)/))
  do j=1,n-2
     call pgline(2,x(j:j+1),(/y(j),y(j)/))
     call pgline(2,(/x(j+1),x(j+1)/),y(j:j+1))
  end do
  call pgline(2,x(n-1:n),(/y(n-1),y(n-1)/))
  
end subroutine horzhist
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Smooth a 1D PDF
!!
!! \param ybin    Height of the bins
!! \param nbin    Number of bins
!! \param smooth  Smoothing factor

subroutine smoothpdf1d(ybin,nbin,smooth)
  implicit none
  integer, intent(in) :: nbin,smooth
  integer :: i,i0,i1,i00
  real, intent(inout) :: ybin(nbin)
  real :: ybin1(nbin),coefs(100),coefs1(100)
  
  ybin1 = ybin
  i0 = min(max(smooth,1),floor(real(nbin-1)/2.))
  
  ! Do all points, except the first and last i0:
  do i=1+i0,nbin-i0
     coefs1(1:2*i0+1) = ybin(i-i0:i+i0)
     call savgol(coefs1(1:2*i0+1),2*i0+1,i0,i0,0,4)
     do i1=1,i0+1
        coefs(i0-i1+2) = coefs1(i1)
     end do
     do i1 = i0+2,2*i0+1
        coefs(3*i0+3-i1) = coefs1(i1)
     end do
     ybin1(i) = 0.
     do i1=1,2*i0+1
        ybin1(i) = ybin1(i) + coefs(i1) * ybin(i+i1-i0-1)
     end do
  end do
  
  i00 = i0-1
  ! Do the first and last i0=i00 points:
  if(i00.ge.2) then
     do i0 = i00,2,-1
        i=1+i0                ! Point at beginning
        coefs1(1:2*i0+1) = ybin(i-i0:i+i0)
        call savgol(coefs1(1:2*i0+1),2*i0+1,i0,i0,0,4)
        do i1=1,i0+1
           coefs(i0-i1+2) = coefs1(i1)
        end do
        do i1 = i0+2,2*i0+1
           coefs(3*i0+3-i1) = coefs1(i1)
        end do
        ybin1(i) = 0.
        do i1=1,2*i0+1
           ybin1(i) = ybin1(i) + coefs(i1) * ybin(i+i1-i0-1)
        end do
        
        i=nbin-i0              ! Point at end
        coefs1(1:2*i0+1) = ybin(i-i0:i+i0)
        call savgol(coefs1(1:2*i0+1),2*i0+1,i0,i0,0,4)
        do i1=1,i0+1
           coefs(i0-i1+2) = coefs1(i1)
        end do
        do i1 = i0+2,2*i0+1
           coefs(3*i0+3-i1) = coefs1(i1)
        end do
        ybin1(i) = 0.
        do i1=1,2*i0+1
           ybin1(i) = ybin1(i) + coefs(i1) * ybin(i+i1-i0-1)
        end do
     end do
  end if
  
  ybin = ybin1
  
end subroutine smoothpdf1d
!***********************************************************************************************************************************

