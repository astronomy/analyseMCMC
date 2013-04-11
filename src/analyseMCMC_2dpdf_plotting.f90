!> \file  analyseMCMC_2dpdf_plotting.f90  Routines and functions to help plot 2D marginalised PDFs

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
!> \brief  Open a screen or plot file for the 2D PDF plot
!!
!! \param  p1           ID of parameter 1
!! \param  p2           ID of parameter 2
!! \param  npdf         Number of the current PDF (I/O)
!! \param  sch          Default (character) scaling
!! \param  project_map  Use map projection?
!!
!! \retval exitcode     Exit code: 0=ok

subroutine open_2D_PDF_plot_file(p1,p2, npdf, sch, project_map, exitcode)
  use SUFR_constants, only: stdErr, stdOut
  use aM_constants, only: use_PLplot
  use analysemcmc_settings, only: outputbasefile,outputtempfile, file, Npdf2D, htmlOutput
  use analysemcmc_settings, only: scrsz,scrrat,pssz,psrat, colour,whitebg,fonttype
  use general_data, only: outputname,outputdir, parNames, htParNs
  use mcmcrun_data, only: parID
  use plot_data, only: bmpsz,bmprat,psclr
  
  implicit none
  integer, intent(in) :: p1,p2
  integer, intent(inout) :: npdf
  real, intent(in) :: sch
  logical, intent(in) :: project_map
  integer, intent(out) :: exitcode
  
  integer :: io, pgopen
  character :: str*(99)
  
  exitcode = 0
  
  write(outputbasefile,'(A)') trim(outputname)//'__pdf2d__'//trim(parNames(parID(p1)))//'-'//trim(parNames(parID(p2)))
  
  if(htmlOutput.ge.1 .and. Npdf2D.gt.0) then
     write(stdOut,'(A)') '<h4>'//trim(htParNs(parID(p1)))//'-'//trim(htParNs(parID(p2)))//':</h4>'
     write(stdOut,'(A)') '<img src="'//trim(outputbasefile)//'.png" title="2D PDF: '//trim(htParNs(parID(p1)))//'-'// &
          trim(htParNs(parID(p2)))//'">'
  end if
  
  write(outputbasefile,'(A)') trim(outputdir)//'/'//trim(outputbasefile)
  
  if(file.eq.0) then
     npdf=npdf+1
     write(str,'(I3,A3)') 200+npdf,'/xs'
     if(.not.use_PLplot) io = pgopen(trim(str))
     if(project_map) then
        call pgpap(scrSz/0.5*scrRat,0.5)
     else
        call pgpap(scrSz,scrRat)
     end if
     if(use_PLplot) io = pgopen(trim(str))
     call pginitl(colour,file,whiteBG)
  end if
  
  if(file.eq.1) then
     write(outputtempfile,'(A)') trim(outputbasefile)
     if(.not.use_PLplot) io = pgopen(trim(outputtempfile)//'.ppm/ppm')
     if(project_map) then
        call pgpap(bmpsz/0.5*bmprat,0.5)
     else
        call pgpap(bmpsz,bmprat)
     end if
     if(use_PLplot) io = pgopen(trim(outputtempfile)//'.ppm/ppm')
     call pginitl(colour,file,whiteBG)
  end if
  
  if(file.ge.2) then
     write(outputtempfile,'(A)') trim(outputbasefile)
     if(.not.use_PLplot) io = pgopen(trim(outputtempfile)//'.eps'//trim(psclr))
     if(project_map) then
        call pgpap(PSsz/0.5*PSrat,0.5)
     else
        call pgpap(PSsz,PSrat)
     end if
     if(use_PLplot) io = pgopen(trim(outputtempfile)//'.eps'//trim(psclr))
     call pginitl(colour,file,whiteBG)
     call pgscf(fonttype)
  end if
  
  if(io.le.0) then
     write(stdErr,'(A,I4)')'   Error:  Cannot open PGPlot device.  Quitting the programme',io
     exitcode = 1
     return
  end if
  
  !call pgscr(3,0.,0.5,0.)
  call pgsch(sch)
  
end subroutine open_2D_PDF_plot_file
!***********************************************************************************************************************************



!***********************************************************************************************************************************
!> \brief  Plot the actual 2D PDF (grey-scale or colour pixels)
!!
!! \param z            2D binned data
!! \param tr           Transformation elements used by PGPlot
!! \param project_map  Use map projection?

subroutine plot_2D_PDF(z, tr, project_map)
  use SUFR_constants, only: stdOut
  use analysemcmc_settings, only: Nbin2Dx,Nbin2Dy, prProgress, plotSky,mapProjection
  
  implicit none
  real, intent(in) :: z(Nbin2Dx+1,Nbin2Dy+1), tr(6)
  logical, intent(in) :: project_map
  
  integer :: clr1,clr2
  
  
  ! Set the colour scheme:
  call set_2D_probability_colours(clr1, clr2)  ! Define the grey scales or colours for the 2D probability areas
  
  
  ! Plot the PDF:
  if(project_map .and. plotSky.ge.2) then
     
     if(prProgress.ge.4) write(stdOut,'(A)',advance="no")'  plotting map projection...'
     call pgimag_project(z, Nbin2Dx+1, Nbin2Dy+1, 1,Nbin2Dx+1, 1,Nbin2Dy+1, 0.,1., clr1,clr2, tr, mapProjection)
     
  else
     
     if(prProgress.ge.4) write(stdOut,'(A)',advance="no")'  plotting 2D PDF...'
     
     ! Plot 2D image - 0: no projection:
     call pgimag_project(z, Nbin2Dx+1, Nbin2Dy+1, 1,Nbin2Dx+1, 1,Nbin2Dy+1, 0.,1., clr1,clr2, tr, 0)
     
     ! Plot 2D image - produces ~2.5x smaller plots - used to give segfaults (still does when plotting to ppm):
     !call pgimag(z,Nbin2Dx+1,Nbin2Dy+1,1,Nbin2Dx+1,1,Nbin2Dy+1,0.,1.,tr)
     
  end if
  
end subroutine plot_2D_PDF
!***********************************************************************************************************************************



!***********************************************************************************************************************************
!> \brief  Define the grey scales or colours for the 2D probability areas
!!
!! \retval clr1  First colour of colour-index range for pgimag
!! \retval clr2  Second colour of colour-index range for pgimag

subroutine set_2D_probability_colours(clr1,clr2)
  use SUFR_system, only: warn
  use analysemcmc_settings, only: colour, Nival, file, normPDF2D
  
  implicit none
  integer, intent(out) :: clr1,clr2
  integer :: ci, clr, minclr,maxclr
  real :: tmpflt
  
  
  if(normPDF2D.lt.4) then  ! Use grey scales
     
     call pgscir(0,nint(1e9))
     call pgqcir(clr,maxclr)  ! Maxclr is device-dependent
     minclr = 30
     if(maxclr.lt.minclr) call warn('Not enough colours on device for 2D plot!',0)
     
     do ci=0,maxclr-minclr  ! Colour indices typically run 0-255, but this is device-dependent.
        ! Reserve ~0-29 for other purposes -> (maxclr-minclr) for these grey scales:
        tmpflt = real((maxclr-minclr) - ci)/real(maxclr-minclr)          ! White background
        call pgscr(minclr+ci,tmpflt,tmpflt,tmpflt)
     end do
     
     call pgscir(minclr,maxclr)  ! Set colour-index range for pgimag
     
     
  else if(normPDF2D.eq.4) then  ! Use colour
     
     minclr = 30
     clr1 = minclr
     clr2 = minclr+Nival
     call pgscir(clr1,clr2)  ! Set colour-index range for pgimag
     
     
     if(colour.eq.0) then
        if(Nival.eq.2) then
           call pgscr(minclr+1,0.5,0.5,0.5)  ! Grey
           call pgscr(minclr+2,0.,0.,0.)     ! Black
        end if
        if(Nival.eq.3) then
           call pgscr(minclr+1,0.7,0.7,0.7)  ! Light grey
           call pgscr(minclr+2,0.4,0.4,0.4)  ! Dark grey
           call pgscr(minclr+3,0.0,0.0,0.0)  ! Black
        end if
        if(Nival.eq.4) then
           call pgscr(minclr+1,0.75,0.75,0.75)  ! Light grey
           call pgscr(minclr+2,0.50,0.50,0.50)  ! Medium grey
           call pgscr(minclr+3,0.25,0.25,0.25)  ! Dark grey
           call pgscr(minclr+4,0.00,0.00,0.00)  ! Black
        end if
        if(Nival.eq.5) then
           call pgscr(minclr+1,0.8,0.8,0.8)  ! Light grey
           call pgscr(minclr+2,0.6,0.6,0.6)  ! Semi-light grey
           call pgscr(minclr+3,0.4,0.4,0.4)  ! Semi-dark grey
           call pgscr(minclr+4,0.2,0.2,0.2)  ! Dark grey
           call pgscr(minclr+5,0.0,0.0,0.0)  ! Black
        end if
     end if
     
     if(colour.ge.1) then
        if(Nival.eq.2) then
           call pgscr(minclr+1,1.,1.,0.)  ! Yellow
           if(file.ge.2) call pgscr(minclr+1,0.8,0.7,0.)  ! Dark yellow
           call pgscr(minclr+2,1.,0.,0.)  ! Red
        end if
        if(Nival.eq.3) then
           call pgscr(minclr+1,0.,0.,1.)  ! Blue
           call pgscr(minclr+2,1.,1.,0.)  ! Yellow
           if(file.ge.2) call pgscr(minclr+2,0.8,0.7,0.)  ! Dark yellow
           call pgscr(minclr+3,1.,0.,0.)  ! Red
        end if
        if(Nival.eq.4) then
           call pgscr(minclr+1,0.,0.,1.)  ! Blue
           call pgscr(minclr+2,0.,1.,0.)  ! Green
           call pgscr(minclr+3,1.,1.,0.)  ! Yellow
           if(file.ge.2) call pgscr(minclr+3,0.8,0.7,0.)  ! Dark yellow
           call pgscr(minclr+4,1.,0.,0.)  ! Red
        end if
        if(Nival.eq.5) then
           call pgscr(minclr+1,0.,0.,1.)  ! Blue
           call pgscr(minclr+2,0.,1.,0.)  ! Green
           call pgscr(minclr+3,1.,1.,0.)  ! Yellow
           if(file.ge.2) call pgscr(minclr+3,0.8,0.7,0.)  ! Dark yellow
           call pgscr(minclr+4,1.,0.5,0.)  ! Orange
           call pgscr(minclr+5,1.,0.,0.)  ! Red
        end if
     end if
     
  end if  ! if(normPDF2D .eq. 4)  ! Use colour
  
end subroutine set_2D_probability_colours
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Plot contours in a 2D PDF
!!
!! \param z            2D binned data
!! \param tr           Transformation elements used by PGPlot
!! \param project_map  Use map projection?
!! \param lw           Default line width

subroutine plot_2D_contours(z, tr, project_map, lw)
  use analysemcmc_settings, only: normPDF2D, plotSky, Nival, Nbin2Dx,Nbin2Dy
  
  implicit none
  real, intent(in) :: z(Nbin2Dx+1,Nbin2Dy+1), tr(6)
  logical, intent(in) :: project_map
  integer, intent(in) :: lw
  
  integer :: i, Ncont
  real :: cont(11)
  
  
  if(normPDF2D.lt.4) then
     Ncont = 11
     do i=1,Ncont
        cont(i) = 0.01 + 2*real(i-1)/real(Ncont-1)
        if(project_map .and. (plotSky.eq.1.or.plotSky.eq.3)) cont(i) = 1.-cont(i)
     end do
     Ncont = min(4,Ncont)  ! Only use the first 4
  else if(normPDF2D.eq.4) then
     Ncont = Nival
     do i=1,Ncont
        cont(i) = max(1. - real(i-1)/real(Ncont-1),0.001)
        !if(project_map) cont(i) = 1.-cont(i)
     end do
  end if
  
  call pgsls(1)
  if((.not.project_map .or. plotSky.ne.1.or.plotSky.ne.3) .and. normPDF2D.ne.4) then  ! First in bg colour
     call pgslw(2*lw)
     call pgsci(0)
     call pgcont(z,Nbin2Dx+1,Nbin2Dy+1,1,Nbin2Dx+1,1,Nbin2Dy+1,cont(1:Ncont),Ncont,tr)
  end if
  
  call pgslw(lw)
  call pgsci(1)
  if(project_map .and. (plotSky.eq.1.or.plotSky.eq.3)) call pgsci(7)
  
  call pgcont(z,Nbin2Dx+1,Nbin2Dy+1,1,Nbin2Dx+1,1,Nbin2Dy+1,cont(1:Ncont),Ncont,tr)  ! Plot contours
  
end subroutine plot_2D_contours
!***********************************************************************************************************************************



!***********************************************************************************************************************************
!> \brief Plot injection value, median, ranges, etc. in 2D PDF
!!
!! \param ic           Chain ID
!! \param p1           ID of parameter 1
!! \param p2           ID of parameter 2
!!
!! \param xmin         Lower limit of horizontal plot range
!! \param xmax         Upper limit of horizontal plot range
!! \param ymin         Lower limit of vertical plot range
!! \param ymax         Upper limit of vertical plot range
!! \param dx           Width of horizontal plot range
!! \param dy           Width of vertical plot range
!!
!! \param sch          Default (character) scaling
!! \param lw           Default line width
!! \param project_map  Use map projection?


subroutine plot_values_in_2D_PDF(ic, p1,p2, xmin,xmax, ymin,ymax, dx,dy, sch,lw, project_map)
  use analysemcmc_settings, only: plotSky, plLmax,plInject,plRange,plMedian, mergeChains, ivals, normPDF2D, mapProjection
  use general_data, only: allDat,startval, c0, ranges,stats, icloglmax,iloglmax, wrap,shifts,shIvals,raCentre
  use mcmcrun_data, only: parID
  
  implicit none
  integer, intent(in) :: ic, p1,p2, lw
  real, intent(in) :: xmin,xmax, ymin,ymax, dx,dy, sch
  logical, intent(in) :: project_map
  
  real :: plx,ply
  character :: delta*(19)
  
  
  if(.not.project_map.or.plotSky.eq.1) then
     
     call pgsci(1)
     
     
     ! Plot max likelihood in 2D PDF:
     if(plLmax.ge.1) then
        call pgsci(1)
        call pgsls(5)
        
        plx = allDat(icloglmax,p1,iloglmax)
        if(wrap(ic,p1).ne.0) plx = mod(plx + shifts(ic,p1), shIvals(ic,p1)) - shifts(ic,p1)
        call pgline(2,(/plx,plx/),(/ymin,ymax/))  ! Max logL
        
        ply = allDat(icloglmax,p2,iloglmax)
        if(wrap(ic,p2).ne.0) ply = mod(ply + shifts(ic,p2), shIvals(ic,p2)) - shifts(ic,p2)
        call pgline(2,(/xmin,xmax/),(/ply,ply/))  ! Max logL
        
        call pgpoint(1,plx,ply,18)
     end if
     
     
     if(project_map .and. (plotSky.eq.1.or.plotSky.eq.3)) call pgsci(0)
     call pgsls(2)
     
     
     ! Plot injection value in 2D PDF:
     if((plInject.eq.1.or.plInject.eq.3).and.(.not.project_map) .or. &
          ((plInject.eq.2.or.plInject.eq.4) .and.  &
          (parID(p1).eq.61.and.parID(p2).eq.62 .or. parID(p1).eq.63.and.parID(p2).eq.64 .or. &
          parID(p2).eq.61.and.parID(p2).eq.67  .or. parID(p2).eq.61.and.parID(p2).eq.68)) ) then
        
        ! CHECK The units of the injection values haven't changed (e.g. from rad to deg) for ic>1 
        ! (but they have for the starting values, why?)
        if(mergeChains.ne.1.or.ic.le.1) then
           
           call pgsls(3)  ! Dash-dotted line for injection value
           call pgsci(1)
           
           ! x:
           plx = startval(ic,p1,1)
           if(wrap(ic,p1).ne.0) plx = mod(plx + shifts(ic,p1), shIvals(ic,p1)) - shifts(ic,p1)
           call pgline(2,(/plx,plx/),(/ymin,ymax/))  ! Injection value
           
           ! y:
           ply = startval(ic,p2,1)
           if(wrap(ic,p2).ne.0) ply = mod(ply + shifts(ic,p2), shIvals(ic,p2)) - shifts(ic,p2)
           call pgline(2,(/xmin,xmax/),(/ply,ply/))  ! Injection value
           
           call pgpoint(1,plx,ply,18)
        end if
     end if  !If plotting injection values in 2D plot
     
     call pgsci(1)
     call pgsls(4)
     
     
     ! Plot starting values in 2D PDF:
     !call pgline(2,(/startval(ic,p1,2),startval(ic,p1,2)/),(/ymin,ymax/))
     !call pgline(2,(/xmin,xmax/),(/startval(ic,p2,2),startval(ic,p2,2)/))
     
     call pgsci(2)
     
     
     ! Plot probability ranges in 2D PDF:
     if(plRange.eq.2.or.plRange.eq.3.or.plRange.eq.5.or.plRange.eq.6.or.plRange.eq.7) then
        write(delta,'(A,I3.3,A)')'\(2030)\d',nint(ivals(c0)*100),'%\u'
        if(nint(ivals(c0)*100).lt.100) write(delta,'(A,I2.2,A)')'\(2030)\d',nint(ivals(c0)*100),'%\u'
        
        call pgsls(1)
        call pgsch(sch*0.6)
        call pgsah(1,45.,0.1)
        
        call pgarro( ranges(ic,c0,p1,3), ymin+dy*0.017*sch, ranges(ic,c0,p1,1), ymin+dy*0.017*sch)
        call pgarro( ranges(ic,c0,p1,3), ymin+dy*0.017*sch, ranges(ic,c0,p1,2), ymin+dy*0.017*sch)
        call pgptxt( ranges(ic,c0,p1,3), ymin+dy*0.033*sch, 0., 0.5, trim(delta) )
        
        call pgarro( xmin+dx*0.023*sch, ranges(ic,c0,p2,3), xmin+dx*0.023*sch, ranges(ic,c0,p2,1) )
        call pgarro( xmin+dx*0.023*sch, ranges(ic,c0,p2,3), xmin+dx*0.023*sch, ranges(ic,c0,p2,2) )
        call pgptxt( xmin+dx*0.01*sch, ranges(ic,c0,p2,3), 90., 0.5, trim(delta) )
        
        if(plRange.eq.7) then  ! Plot dotted lines for 1D probability ranges in 2D plot
           call pgsls(4)
           call pgline(2, (/ranges(ic,c0,p1,1), ranges(ic,c0,p1,1)/), (/ymin, ymax/) )
           call pgline(2, (/ranges(ic,c0,p1,2), ranges(ic,c0,p1,2)/), (/ymin, ymax/) )
           
           call pgline(2, (/xmin, xmax/), (/ranges(ic,c0,p2,1), ranges(ic,c0,p2,1)/) )
           call pgline(2, (/xmin, xmax/), (/ranges(ic,c0,p2,2), ranges(ic,c0,p2,2)/) )
           call pgsls(1)
        end if
     end if
     
     call pgsch(sch)
     call pgsls(2)
     
     
     ! Plot medians in 2D PDF:
     if(plMedian.eq.2.or.plMedian.eq.3.or.plMedian.eq.5.or.plMedian.eq.6) then
        call pgline(2,(/stats(ic,p1,1),stats(ic,p1,1)/),(/ymin,ymax/))
        call pgline(2,(/xmin,xmax/),(/stats(ic,p2,1),stats(ic,p2,1)/))
        call pgpoint(1,stats(ic,p1,1),stats(ic,p2,1),18)
     end if
     
     call pgsls(1)
     
  end if  ! if(.not.project_map.or.plotSky.eq.1)
  
  
  
  ! Plot big symbol at injection position in sky map:
  if(project_map .and. (plInject.eq.1.or.plInject.eq.3)) then
     call pgsch(sch*1.5)               ! Use 1.5 for plsym=8, 2 for plsym=18
     call pgslw(lw*2)
     call pgsci(9)
     if(normPDF2D.eq.4) call pgsci(1)  ! Black
     
     plx = startval(ic,p1,1)
     ply = startval(ic,p2,1)
     
     if(plotSky.eq.2.or.plotSky.eq.4) call project_skymap(plx,ply,raCentre,mapProjection)
     call pgpoint(1,plx,ply,8)
     
     call pgsch(sch)
     call pgslw(lw)
     call pgsci(1)
  end if
  
  
end subroutine plot_values_in_2D_PDF
!***********************************************************************************************************************************



!***********************************************************************************************************************************
!> \brief  Print 2D PDF plot axes, axis labels and plot title


subroutine plot_2D_PDF_axes_labels_titles(p1,p2, sch,flw, project_map)
  use SUFR_text, only: replace_substring
  
  use aM_constants, only: use_PLplot
  use analysemcmc_settings, only: prIval, normPDF2D, nIval,ivals, fonttype, quality, plotSky, changeVar
  use general_data, only: pgParNs,pgUnits
  use mcmcrun_data, only: parID
  use stats_data, only: probAreas
  
  implicit none
  
  integer, intent(in) :: p1,p2, flw
  real, intent(in) :: sch
  logical, intent(in) :: project_map
  
  integer :: c, i
  real :: area, just
  character :: string*(99), ivalstr*(99), areaunit*(99), tmpStr1*(19),tmpStr2*(19)
  
  
  ! Plot coordinate axes and axis labels in 2D PDF:
  call pgsls(1)
  call pgslw(flw)
  call pgsci(1)
  
  if(.not.project_map) call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)
  
  !if(plotSky.eq.1) then
  !   call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0) !Box, ticks, etc in white
  !   call pgsci(1)
  !   call pgbox('N',0.0,0,'N',0.0,0) !Number labels in black
  !end if
  
  if(use_PLplot) then
     call pgmtxt('B',5.0,0.5,0.5,trim(pgParNs(parID(p1))))
     call pgmtxt('L',5.0,0.5,0.5,trim(pgParNs(parID(p2))))
  else
     call pgmtxt('B',2.2,0.5,0.5,trim(pgParNs(parID(p1))))
     call pgmtxt('L',2.2,0.5,0.5,trim(pgParNs(parID(p2))))
  end if
  
  
  ! Print 2D probability ranges in plot title:
  if(prIval.ge.1.and.normPDF2D.eq.4) then
     string = ' '
     do c = 1,Nival
        write(ivalstr,'(F5.1,A1)') ivals(c)*100,'%'
        if(fonttype.eq.2) then
           if(abs(ivals(c)-0.6827).lt.0.001) write(ivalstr,'(A)') '1\(2144)'
           if(abs(ivals(c)-0.9545).lt.0.001) write(ivalstr,'(A)') '2\(2144)'
           if(abs(ivals(c)-0.9973).lt.0.001) write(ivalstr,'(A)') '3\(2144)'
        else
           if(abs(ivals(c)-0.6827).lt.0.001) write(ivalstr,'(A)') '1\(0644)'
           if(abs(ivals(c)-0.9545).lt.0.001) write(ivalstr,'(A)') '2\(0644)'
           if(abs(ivals(c)-0.9973).lt.0.001) write(ivalstr,'(A)') '3\(0644)'
        end if
        
        areaunit = trim(pgUnits(parID(p1)))//' '//trim(pgUnits(parID(p2)))
        if(trim(pgUnits(parID(p1))) .eq. trim(pgUnits(parID(p2)))) areaunit = trim(pgUnits(parID(p1)))//'\u2\d'  ! mm->m^2
        if(trim(areaunit).eq.'\(2218)\u2\d') areaunit = 'deg\u2\d'  ! Square degrees
        if(changeVar.ge.1 .and. ((parID(p1).eq.31.and.parID(p2).eq.32) .or. (parID(p1).eq.52.and.parID(p2).eq.51)) ) &
             areaunit = 'deg\u2\d'  ! Square degrees for sky map and iota-psi 'map' - exception in calc_2d_areas()
        areaunit = ' '//trim(areaunit)  ! Add space between value and unit
        
        
        area = probAreas(p1,p2,c,3)
        
        if(area.lt.1.) then
           write(string,'(A,F5.2,A)') trim(ivalstr)//':',area,trim(areaunit)
        else if(area.lt.10.) then
           write(string,'(A,F4.1,A)') trim(ivalstr)//':',area,trim(areaunit)
        else if(area.lt.100.) then
           write(string,'(A,F5.1,A)') trim(ivalstr)//':',area,trim(areaunit)
        else if(area.lt.1000.) then
           write(string,'(A,I4,A)') trim(ivalstr)//':',nint(area),trim(areaunit)
        else if(area.lt.10000.) then
           write(string,'(A,I5,A)') trim(ivalstr)//':',nint(area),trim(areaunit)
        else
           write(string,'(A,I6,A)') trim(ivalstr)//':',nint(area),trim(areaunit)
        end if
        
        if(area.lt.0.09.or.area.ge.1.e5) write(string,'(A,1p,G8.1,A)') trim(ivalstr)//':',abs(area),trim(areaunit)
        
        ! Replace exponentials by x10^:
        do i=-99,99
           write(tmpStr1,'(A,I3.2)') 'E',i
           call replace_substring(tmpStr1, 'E ', 'E+')
           write(tmpStr2,'(A,I3.2,A)') '\(727)10\u',i,'\d'
           call replace_substring(string, trim(tmpStr1), trim(tmpStr2))
        end do
        call replace_substring(string, '\(727)10\u-0', '\(727)10\u-')
        call replace_substring(string, '\(727)10\u 0', '\(727)10\u')
        call replace_substring(string, '\(727)10\u ', '\(727)10\u')
        
        
        call pgsch(sch*0.8)  ! Needed to fit three sigma values in
        if(quality.eq.3) call pgsch(sch*0.6)  ! Poster
        if(quality.eq.91) then  ! NINJA
           call pgsch(sch)
           if(fonttype.eq.2) then
              write(string,'(I2,A7,A10)') c,'\(2144)',''
           else
              write(string,'(I2,A7,A10)') c,'\(0644)',''
           end if
        end if
        
        just = (real(c-1)/real(Nival-1) - 0.5)*0.7 + 0.5
        call pgsci(30+Nival+1-c)
        if(project_map .and. plotSky.ge.2) then
           if(use_PLplot) then
              call pgmtxt('T',3.0,just,0.5,trim(string))  ! Print title
           else
              call pgmtxt('T',1.0,just,0.5,trim(string))  ! Print title
           end if
        else
           if(use_PLplot) then
              call pgmtxt('T',2.0,just,0.5,trim(string))  ! Print title
           else
              call pgmtxt('T',0.5,just,0.5,trim(string))  ! Print title
           end if
        end if
        call pgsch(sch)
        
     end do  ! c=1,Nival
     
     call pgsci(1)
  end if  ! if(prIval.ge.1.and.normPDF2D.eq.4)
  
end subroutine plot_2D_PDF_axes_labels_titles
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Convert 2D PDF plot
!!
!! \param p1          ID of parameter 1
!! \param p2          ID of parameter 2
!! \param countplots  Count of the current plot

subroutine convert_2D_PDF_plot(p1,p2, countplots)
  use SUFR_constants, only: stdErr
  use analysemcmc_settings, only: file,outputtempfile,outputbasefile, Npdf2D
  use general_data, only: parNames
  use mcmcrun_data, only: parID
  use plot_data, only: bmpxpix, unSharppdf2d
  implicit none
  
  integer, intent(in) :: p1,p2, countplots
  
  integer :: status,system
  character :: convopts*(99)
  logical :: ex
  
  if(file.eq.1) then
     call pgend
     
     inquire(file=trim(outputtempfile)//'.ppm', exist=ex)
     if(ex) then
        convopts = '-resize '//trim(bmpxpix)//' -depth 8 -unsharp '//trim(unSharppdf2d)
        !if(p2.eq.1 .or. countplots.eq.Npdf2D) then        ! Doesn't seem to matter much which one you choose
        if(mod(p2,6).eq.0 .or. countplots.eq.Npdf2D) then  ! Doesn't seem to matter much which one you choose
           
           ! Convert the last plot in the foreground, so that the process finishes before deleting the original file:
           status = system('convert '//trim(convopts)//' '//trim(outputtempfile)//'.ppm '//trim(outputbasefile)//'.png')
           
        else  ! in the background
           
           status = system('convert '//trim(convopts)//' '//trim(outputtempfile)//'.ppm '//trim(outputbasefile)//'.png &')
           
        end if
        
        if(status.ne.0) write(stdErr,'(A)')'  Error converting plot for '//trim(parNames(parID(p1)))//'-'// &
             trim(parNames(parID(p2)))
     end if
  end if
  
  
  if(file.ge.2) then
     call pgend
     if(file.eq.3) then
        status = system('eps2pdf '//trim(outputtempfile)//'.eps &> /dev/null')
        if(status.ne.0) write(stdErr,'(A)')'  Error converting plot for '//trim(parNames(parID(p1)))//'-'// &
             trim(parNames(parID(p2)))
     end if
  end if
  
end subroutine convert_2D_PDF_plot
!***********************************************************************************************************************************



!***********************************************************************************************************************************
!> \brief Remove all the .ppm files, create thumbnails and create the index.html
!!
!! \param j1  First parameter to plot
!! \param j2  Last parameter to plot

subroutine removeppm_createthumbnails_createhtml_2D_PDF(j1,j2)
  use SUFR_constants, only: stdErr
  
  use analysemcmc_settings, only: file, htmlOutput, Npdf2D,PDF2Dpairs, bmpXSz,bmpYSz, scFac
  use aM_constants, only: use_PLplot
  use mcmcrun_data, only: revID,parID
  use general_data, only: outputname,outputdir,parNames,htParNs, fixedpar
  use plot_data, only: bmpsz,bmprat,bmpxpix,pltsz,pltrat
  
  implicit none
  
  integer, intent(in) :: j1,j2
  integer :: i, p1,p2, countplots,plotthis, status,system
  character :: basefile*(199),tempfile*(99)
  logical :: ex
  
  
  if(file.eq.1) then
     countplots = 0
     
     if(htmlOutput.ge.1) write(51,'(4x,A)')'<table>'
     
     do p2=j1,j2  ! p2 is in the outer loop to get the proper HTML matrix
        
        if(fixedpar(p2).ge.1) cycle
        
        if(htmlOutput.ge.1) then
           if(p2.eq.j1) then
              write(51,'(6x,A)')'<tr>'
              write(51,'(8x,A)')'<td></td>'
              do p1=j1,j2
                 if(fixedpar(p1).ge.1) cycle
                 write(51,'(8x,A)')'<td align="center"><h1>'//trim(htParNs(parID(p1)))//'</h1></td>'
              end do
              write(51,'(8x,A)')'<td></td>'
              write(51,'(6x,A)')'</tr>'
           end if
           
           write(51,'(6x,A)')'<tr>'
           write(51,'(8x,A)')'<td align="center"><h1>'//trim(htParNs(parID(p2)))//'</h1></td>'
        end if
        
        
        
        do p1=j1,j2
           
           if(fixedpar(p1).ge.1) cycle
           
           if(Npdf2D.ge.0) then
              plotthis = 0  ! Determine to plot or save this combination of j1/j2 or p1/p2
              do i=1,Npdf2D
                 if(p1.eq.revID(PDF2Dpairs(i,1)).and.p2.eq.revID(PDF2Dpairs(i,2))) plotthis = 1  
              end do
              if(plotthis.eq.0) cycle
           end if
           
           if(p1.eq.p2) then
              if(htmlOutput.ge.1) then
                 write(51,'(8x,A)')'<td align="center"><h1>'//trim(htParNs(parID(p1)))//'</h1></td>'
              end if
              cycle
           end if
           
           countplots = countplots + 1  ! The current plot is number countplots
           write(basefile,'(A)') trim(outputname)//'__pdf2d__'//trim(parNames(parID(p1)))//'-'//trim(parNames(parID(p2)))
           write(tempfile,'(A)') trim(outputdir)//'/'//trim(basefile)
           status = system('rm -f '//trim(tempfile)//'.ppm')
           
           
           if(htmlOutput.ge.1) then
              
              inquire(file=trim(tempfile)//'.png', exist=ex)
              if(ex) then
                 
                 ! Convert the last plot in the foreground, so that the process finishes before deleting the original file:
                 !if(p1.eq.1 .or. countplots.eq.Npdf2D) then        ! Doesn't seem to matter much which one you choose
                 if(mod(p1,6).eq.0 .or. countplots.eq.Npdf2D) then  ! Doesn't seem to matter much which one you choose
                    status = system('convert -resize 200x200 '//trim(tempfile)//'.png '//trim(tempfile)//'_thumb.png')
                 else
                    status = system('convert -resize 200x200 '//trim(tempfile)//'.png '//trim(tempfile)//'_thumb.png &')
                 end if
                 if(status.ne.0) write(stdErr,'(A)')'  Error creating thumbnail for '//trim(parNames(parID(p1)))//'-'// &
                      trim(parNames(parID(p2)))
              end if
              
              write(51,'(8x,A)')'<td>'
              write(51,'(8x,A)')'  <a href="'//trim(basefile)//'.png">'
              write(51,'(8x,A)')'    <img src="'//trim(basefile)//'_thumb.png" title="2D PDF: '//trim(ParNames(parID(p1)))//'-'// &
                   trim(ParNames(parID(p2)))//'">'
              write(51,'(8x,A)')'  </a>'
              write(51,'(8x,A)')'</td>'
           end if
           
        end do  ! p1=j1,j2
        
        
        if(htmlOutput.ge.1) then
           write(51,'(8x,A)')'<td align="center"><h1>'//trim(htParNs(parID(p2)))//'</h1></td>'
           write(51,'(6x,A)')'</tr>'
           
           if(p2.eq.j2) then
              write(51,'(6x,A)')'<tr>'
              write(51,'(8x,A)')'<td></td>'
              do p1=j1,j2
                 if(fixedpar(p1).ge.1) cycle
                 write(51,'(8x,A)')'<td align="center"><h1>'//trim(htParNs(parID(p1)))//'</h1></td>'
              end do
              write(51,'(8x,A)')'<td></td>'
              write(51,'(6x,A)')'</tr>'
           end if
        end if
        
     end do  ! p2=j1,j2
     
     
     if(htmlOutput.ge.1) then
        write(51,'(4x,A)')'</table>'
        
        
        write(51,'(A)')'  <br><hr><br>'
        
        write(51,'(A)', advance='no') '<b>'
        call print_code_version(51, use_PLplot)
        write(51,'(A)') '</b>'
        
        write(51,'(A)') '<br><br>'
        call print_rundata(51)
        
        
        write(51,'(A)') '<br><br>'
        write(51,'(4x,A)') '<font size="2">'
        write(51,'(6x,A)') '<a href="index.html#2dpdfs" title="Go back to the main page">Main page</a>'
        write(51,'(4x,A)') '</font>'
        
        
        write(51,'(A)')'  </body>'
        write(51,'(A)')'</html>'
     end if
     
  end if  ! if(file.eq.1)
  
  
  if(htmlOutput.eq.1) then
     bmpXSz = 1000
     bmpYSz =  700
     
     call compBitmapSize(bmpXSz,bmpYSz, scFac, bmpsz,bmprat)  ! Determine plot size and ratio
     write(bmpxpix,'(I4)') bmpXSz  ! Used as a text string by convert
     pltsz = bmpsz
     pltrat = bmprat
  end if
  
end subroutine removeppm_createthumbnails_createhtml_2D_PDF
!***********************************************************************************************************************************

