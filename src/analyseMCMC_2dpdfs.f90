!> \file analyseMCMC_2dpdfs.f90  Routines to compute and plot two-dimensional PDFs

! 
! LICENCE:
! 
! Copyright 2007-2012 Marc van der Sluys
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
!> \brief  Plot 2D marginalised PDFs
!!
!! \retval exitcode  Exit code: 0=ok

subroutine pdfs2d(exitcode)
  use SUFR_constants, only: stdOut,stdErr, rh2r
  use SUFR_system, only: swapreal
  use SUFR_text, only: replace_substring
  
  use analysemcmc_settings, only: update,prProgress,file,pssz,quality,fontsize2d,maxChs
  use analysemcmc_settings, only: Npdf2D,PDF2Dpairs,htmlOutput,bmpXSz,bmpYSz,scFac,Nbin2Dx,Nbin2Dy,plotSky
  use analysemcmc_settings, only: savePDF,plot,plPDF1D,plPDF2D
  use general_data, only: outputname,outputdir,parNames, nfixedpar, maxIter, raCentre,raShift
  use mcmcrun_data, only: totpts,revID, nMCMCpar
  use plot_data, only: bmpsz,bmprat,bmpxpix,pltsz,pltrat
  
  implicit none
  integer, intent(out) :: exitcode
  
  real,allocatable :: z(:,:),zs(:,:,:)  ! These depend on nbin2d, allocate after reading input file
  
  integer :: i,j,j1,j2,p1,p2,ic,lw,flw
  integer :: npdf,countplots,totplots
  real :: tr(6), sch, xmin,xmax, ymin,ymax, dx,dy
  logical :: project_map,sky_position,binary_orientation, create_this_2D_PDF
  
  
  exitcode = 0
  countplots = 0
  ic = 1  ! Can only do one chain
  
  j1 = 1
  j2 = nMCMCpar
  
  if(htmlOutput.ge.1) then
     write(stdOut,'(A)') '<a name="2dpdfs"></a><h3>2D PDFs:</h3>'
  else
     if(prProgress.ge.1.and.plot.eq.0.and.savePDF.eq.1.and.plPDF1D.eq.0) write(stdOut,'(A)',advance="no")'  Saving'
     if(prProgress.ge.1.and.update.eq.0.and.Npdf2D.ge.0) write(stdOut,'(A)',advance="no")'  2D pdfs: '
  end if
  
  if(Npdf2D.lt.0) then  ! Plot all 2D PDFs
     totplots = 0
     do j=j1,j2-nfixedpar
        totplots = totplots + j - j1
     end do
     totplots = totplots*2  ! Since we're plotting Mc-eta as well as eta
     
     if(htmlOutput.ge.1) then
        write(stdOut,'(A,I3,A,I3,A)') '<a href="2dpdfs.html"><b>*all* ',totplots, ' 2D PDFs for all combinations of the', &
             j2-j1+1-nfixedpar,' non-fixed parameters</b></a>'
     else if(prProgress.ge.1.and.update.eq.0) then
        write(stdOut,'(A,I3,A,I3,A,/)') '  *all* ',totplots,' 2D PDFs for all combinations of the', &
             j2-j1+1-nfixedpar,' non-fixed parameters:'
     end if
  end if
  
  
  ! Check consistency of PDF2Dpairs():
  do i=1,Npdf2D
     if(revID(PDF2Dpairs(i,1)).eq.0) write(stdErr,'(/,A)')'  * Warning:  pdfs2d():  parameter '// &
          trim(parNames(PDF2Dpairs(i,1)))//' is not defined, check plPars() in the input file.  Skipping...'
     if(revID(PDF2Dpairs(i,2)).eq.0) write(stdErr,'(/,A)')'  * Warning:  pdfs2d():  parameter '// &
          trim(parNames(PDF2Dpairs(i,2)))//' is not defined, check plPars() in the input file.  Skipping...'
  end do
  
  
  if(htmlOutput.eq.1) then
     bmpXSz = 700
     bmpYSz = 700
     
     call compBitmapSize(bmpXSz,bmpYSz, scFac, bmpsz,bmprat)  ! Determine plot size and ratio
     write(bmpxpix,'(I4)') bmpXSz                             ! Used as a text string by convert
     pltsz = bmpsz
     pltrat = bmprat
     
     ! Open the 2D PDF matrix HTML page using unit 51:
     call create_html_2dpdf_file(51)
  end if
  
  
  ! Autodetermine number of bins for 2D PDFs:
  if(Nbin2Dx.le.0) then
     if(totpts.le.100) then
        Nbin2Dx = floor(2*sqrt(real(totpts))/pltrat)
        Nbin2Dx = max(Nbin2Dx,5)
        Nbin2Dy = floor(2*sqrt(real(totpts)))           ! Same as for 1D case (~50)
     else
        Nbin2Dx = floor(10*log10(real(totpts))/pltrat)  
        Nbin2Dx = max(Nbin2Dx,5)
        Nbin2Dy = floor(10*log10(real(totpts)))         ! Same as for 1D case (~50)
     end if
  end if
  if(Nbin2Dy.eq.0) Nbin2Dy = Nbin2Dx
  if(Nbin2Dy.le.-1) Nbin2Dy = nint(real(Nbin2Dx)*pltrat)
  
  
  ! Report number of bins used:
  if(prProgress.ge.1.and.plot.eq.1.and.update.eq.0.and.Npdf2D.ge.0) then
     if(Nbin2Dx.lt.100) then
        write(stdOut,'(A1,I2,A1)',advance="no") '(',Nbin2Dx,'x'
     else
        write(stdOut,'(A1,I3,A1)',advance="no") '(',Nbin2Dx,'x'
     end if
     if(Nbin2Dy.lt.100) then
        write(stdOut,'(I2,A7)',advance="no") Nbin2Dy,' bins) '
     else
        write(stdOut,'(I3,A7)',advance="no") Nbin2Dy,' bins) '
     end if
  end if
  
  
  ! Allocate memory:
  allocate(z(Nbin2Dx+1,Nbin2Dy+1),zs(maxChs,Nbin2Dx+1,Nbin2Dy+1))
  
  
  sch = fontsize2d
  if(plot.eq.1) then
     if(file.eq.0) then
        lw = 1
        flw = nint(1*fontsize2d)  ! Font lw
        sch = 1.5*fontsize2d
     end if
     if(file.ge.1) then
        !if(file.ge.2.and.multipagefile) io = pgopen(trim(outputdir)//'/pdf2d.eps'//trim(psclr))
        lw = 3
        flw = ceiling(3.*fontsize2d)  ! Font lw
        if(quality.eq.91) flw = nint(3*fontsize2d)  ! NINJA
        sch = 1.5*fontsize2d
        if(quality.eq.3) then  ! Poster
           lw = 4
           flw = nint(3*fontsize2d)  ! Font lw
           sch = 2.*fontsize2d
        end if
        if(PSsz.lt.5) sch = sch * sqrt(5.0/PSsz)
     end if
     !if(file.ge.2.and.multipagefile) then
     !   if(io.le.0) then
     !      write(stdErr,'(A,I4)')'  Error:  cannot open PGPlot device.  Quitting the programme',io
     !      exitcode = 1
     !      return
     !   end if
     !   call pgscf(fonttype)
     !   call pginitl(colour,file,whiteBG)
     !end if
  end if !if(plot.eq.1)
  
  
  if(savePDF.eq.1) then
     open(unit=30,action='write',form='formatted',status='replace',file=trim(outputdir)//'/'//trim(outputname)//'__pdf2d.dat')
     write(30,'(5I6,T101,A)')j1,j2,1,Nbin2Dx,Nbin2Dy,'Plot parameter 1,2, total number of chains, number of bins x,y '// &
          ' (j1,j2,(1),Nbin2Dx,Nbin2Dy)'
  end if
  
  npdf=0  ! Count iterations to open windows with different numbers
  do p1=j1,j2
     do p2=j1,j2
        
        ! Create a 2D PDF for this combination of j1/j2 or p1/p2?
        if(.not. create_this_2D_PDF(p1,p2, countplots,totplots)) cycle
        
        ! Identify special combinations of parameters:
        call identify_special_combinations_of_parameters(p1,p2, sky_position, binary_orientation, project_map)
        
        
        ! Open the plot file:
        exitcode = 0
        if(plot.eq.1) call open_2D_PDF_plot_file(p1,p2, npdf, sch, project_map, exitcode)
        if(exitcode.ne.0) return
        
        
        ! Determine plot/binning ranges:
        call determine_2D_PDF_binning_plot_ranges(ic,p1,p2, project_map, xmin,xmax, ymin,ymax, dx,dy)
        
        ! Prepare binning for a 2D sky map:
        if(plot.eq.1 .and. project_map) call prepare_skymap_binning(xmin,xmax, ymin,ymax)
        
        
        ! Bin data and 'normalise' 2D PDF:
        call bin_and_normalise_2D_data(ic,p1,p2, xmin,xmax, ymin,ymax, z,tr, sky_position,binary_orientation)
        
        
        
        ! Swap RA boundaries for RA-Dec plot in 2D PDF, when plotting a map:
        if(sky_position .and. plotSky.ge.1) then
           call swapreal(xmin, xmax)
           dx = -dx
        end if
        
        
        
        ! Plot 2D PDF:
        if(plot.eq.1) then
           
           ! Set plot ranges for whole-sky map.  Does not affect binning:
           if(project_map .and. plotSky.ge.2) then
              raCentre = raCentre
              xmin = raCentre + 12.  ! Must be the larger of the two
              xmax = raCentre - 12.
              ymin = -90.
              ymax = 90.
           end if
           
           
           ! Set viewport on page:
           call pgsch(sch)
           if(project_map .and. plotSky.ge.2) then
              call pgsvp(0.08*sch,0.95,0.08*sch,1.0-0.05*sch)   ! Make room for title and +90deg label
           else
              call pgsvp(0.08*sch,0.98,0.07*sch,1.0-0.033*sch)  ! Make room for title.
              ! Since sch is typically ~1.5*fontsize2d: 0.95 -> 1-0.05*fontsize ~ 1-0.03*sch
           end if
           
           ! Setup plot window:
           call pgswin(xmin,xmax,ymin,ymax)
           if(project_map .and. (plotSky.eq.1.or.plotSky.eq.3).and.file.ge.2) then  ! Need dark background
              call pgsci(1)
              call pgrect(xmin,xmax,ymin,ymax)
           end if
           
           
           ! Plot the actual 2D PDF (grey-scale or colour pixels):
           if(plPDF2D.eq.1.or.plPDF2D.eq.2) call plot_2D_PDF(z, tr, project_map)
           
           ! Plot stars in 2D PDF (over the grey scales, but underneath contours, lines, etc):
           if(project_map .and. (plotSky.eq.1.or.plotSky.eq.3)) then
              call pgswin(xmin*15,xmax*15,ymin,ymax)  ! Map works in degrees
              call plotthesky(xmin*15,xmax*15,ymin,ymax,raShift*rh2r)
              call pgswin(xmin,xmax,ymin,ymax)
           end if
           call pgsci(1)
           
           ! Plot contours in 2D PDF:
           if((plPDF2D.eq.1.or.plPDF2D.eq.3) .and. (.not.project_map .or. plotSky.eq.1.or.plotSky.eq.3)) &
                call plot_2D_contours(z, tr, project_map, lw)
           
           
           ! Plot injection value, median, ranges, etc. in 2D PDF:
           call plot_values_in_2D_PDF(ic, p1,p2, xmin,xmax, ymin,ymax, dx,dy, sch,lw, project_map)
           
           
           ! Print axes, axis labels and plot title:
           call plot_2D_PDF_axes_labels_titles(p1,p2, sch,flw, project_map)
           
           
           countplots = countplots + 1  ! The current plot is number countplots
           
           ! Convert plot:
           call convert_2D_PDF_plot(p1,p2, countplots)
           
        end if  ! if(plot.eq.1)
        
        
        ! Save binned 2D PDF data:
        if(savePDF.eq.1) call save_binned_2D_PDF_data(ic,p1,p2, xmin,xmax,ymin,ymax, z,tr, 30)  ! output unit: 30
        
        
        
     end do  ! p2
  end do  ! p1
  
  
  
  
  
  
  
  if(savePDF.eq.1) close(30)
  
  if(plot.eq.1.and.file.ne.1) call pgend
  
  !if(plot.eq.1.and.file.ge.2.and.multipagefile) then
  !   if(abs(j2-j1).le.1) then
  !      if(file.eq.3) status = system('eps2pdf pdf2d.eps  -o '//trim(outputdir)//'/'//trim(outputname)//'__pdf2d_'// &
  !trim(parNames(parID(j1)))//'-'//trim(parNames(parID(j2)))//'.pdf  &> /dev/null')
  !      status = system('mv -f pdf2d.eps '//trim(outputdir)//'/'//trim(outputname)//'__pdf2d_'// &
  !trim(parNames(parID(j1)))//'-'//trim(parNames(parID(j2)))//'.eps')
  !   else
  !      if(file.eq.3) status = system('eps2pdf pdf2d.eps  -o '//trim(outputdir)//'/'//trim(outputname)// &
  !'__pdf2d.pdf  &> /dev/null')
  !      status = system('mv -f pdf2d.eps '//trim(outputdir)//'/'//trim(outputname)//'__pdf2d.eps')
  !   end if
  !end if
  
  
  ! Remove all the .ppm files, create thumbnails and create the index.html:
  if(plot.eq.1) call removeppm_createthumbnails_createhtml_2D_PDF(j1,j2)
  
  
  
end subroutine pdfs2d
!***********************************************************************************************************************************





