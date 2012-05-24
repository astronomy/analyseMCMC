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
  use SUFR_constants, only: stdOut,stdErr, cursorup, rh2r
  use SUFR_system, only: swapreal
  use SUFR_text, only: replace_substring
  
  use aM_constants, only: use_PLplot
  use analysemcmc_settings, only: update,prProgress,file,scrsz,scrrat,pssz,psrat,fonttype,colour,whitebg,quality
  use analysemcmc_settings, only: fontsize2d,maxChs
  use analysemcmc_settings, only: Npdf2D,PDF2Dpairs,html,bmpXSz,bmpYSz,scFac,Nbin2Dx,Nbin2Dy,plotSky
  use analysemcmc_settings, only: savePDF,plot,plPDF1D,plPDF2D
  use general_data, only: outputname,outputdir,startval,parNames, nfixedpar
  use general_data, only: selDat,stats,ranges,c0,n,maxIter,fixedpar, raCentre,raShift
  use mcmcrun_data, only: totpts,revID,parID, nMCMCpar
  use plot_data, only: psclr,bmpsz,bmprat,bmpxpix,unSharppdf2d,pltsz,pltrat
  
  implicit none
  integer, intent(out) :: exitcode
  
  real,allocatable :: z(:,:),zs(:,:,:)  ! These depend on nbin2d, allocate after reading input file
  
  integer :: i,j,j1,j2,p1,p2,ic,lw,flw,io,status,system,pgopen
  integer :: npdf,plotthis,countplots,totplots
  real :: tr(6), sch, xmin,xmax, ymin,ymax, dx,dy
  character :: str*(99),tempfile*(99), outputbasefile*(199), convopts*(99)
  logical :: project_map,sky_position,binary_orientation, ex
  
  
  exitcode = 0
  countplots = 0
  ic = 1  ! Can only do one chain
  
  j1 = 1
  j2 = nMCMCpar
  
  if(prProgress.ge.1.and.plot.eq.0.and.savePDF.eq.1.and.plPDF1D.eq.0) write(stdOut,'(A)',advance="no")'  Saving'
  if(prProgress.ge.1.and.update.eq.0.and.Npdf2D.ge.0) write(stdOut,'(A)',advance="no")'  2D pdfs: '
  if(Npdf2D.lt.0) then
     totplots = 0
     do j=j1,j2-nfixedpar
        totplots = totplots + j - j1
     end do
     if(prProgress.ge.1.and.update.eq.0) write(stdOut,'(A,I3,A,I3,A,/)')'  *all* ',totplots, &
          ' 2D PDFs for the all combinations of the',j2-j1+1-nfixedpar,' non-fixed parameters: '
  end if
  
  
  ! Check consistency of PDF2Dpairs():
  do i=1,Npdf2D
     if(revID(PDF2Dpairs(i,1)).eq.0) write(stdErr,'(/,A)')'  * Warning:  pdfs2d():  parameter '// &
          trim(parNames(PDF2Dpairs(i,1)))//' is not defined, check plPars() in the input file.  Skipping...'
     if(revID(PDF2Dpairs(i,2)).eq.0) write(stdErr,'(/,A)')'  * Warning:  pdfs2d():  parameter '// &
          trim(parNames(PDF2Dpairs(i,2)))//' is not defined, check plPars() in the input file.  Skipping...'
  end do
  
  
  if(html.eq.1) then
     bmpXSz = 700
     bmpYSz = 700
     
     call compBitmapSize(bmpXSz,bmpYSz, scFac, bmpsz,bmprat)  ! Determine plot size and ratio
     write(bmpxpix,'(I4)')bmpXSz                              ! Used as a text string by convert
     pltsz = bmpsz
     pltrat = bmprat
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
        
        if(Npdf2D.ge.0) then
           plotthis = 0  ! Determine to plot or save this combination of j1/j2 or p1/p2
           do i=1,Npdf2D
              if(p1.eq.revID(PDF2Dpairs(i,1)).and.p2.eq.revID(PDF2Dpairs(i,2))) plotthis = 1  ! Use PDF2Dpairs from the input file
           end do
           if(plotthis.eq.0) cycle
           if(prProgress.ge.1.and.update.eq.0) write(stdOut,'(A)',advance="no")trim(parNames(parID(p1)))//'-'// &
                trim(parNames(parID(p2)))//' '
           
        else  ! Npdf2D.lt.0 - all combinations of non-fixed parameters
           
           if(p2.le.p1) cycle
           if(fixedpar(p1)+fixedpar(p2).ge.1) cycle
           if(stdOut.lt.10) then
              write(stdOut,*)cursorup  ! Move cursor up 1 line
              if(prProgress.ge.1.and.update.eq.0) write(stdOut,'(F7.1,A)') real(countplots+1)/real(totplots)*100, &
                   '%    ('//trim(parNames(parID(p1)))//'-'//trim(parNames(parID(p2)))//')                                      '
           end if
        end if
        
        
        ! Identify special combinations of parameters:
        call identify_special_combinations_of_parameters(p1,p2, sky_position, binary_orientation, project_map)
        
        
        if(plot.eq.1) then
           
           !*** Open the plot file:
           write(outputbasefile,'(A)') trim(outputdir)//'/'//trim(outputname)//'__pdf2d__'// &
                trim(parNames(parID(p1)))//'-'//trim(parNames(parID(p2)))
           
           if(file.eq.0) then
              npdf=npdf+1
              write(str,'(I3,A3)')200+npdf,'/xs'
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
              write(tempfile,'(A)') trim(outputbasefile)
              if(.not.use_PLplot) io = pgopen(trim(tempfile)//'.ppm/ppm')
              if(project_map) then
                 call pgpap(bmpsz/0.5*bmprat,0.5)
              else
                 call pgpap(bmpsz,bmprat)
              end if
              if(use_PLplot) io = pgopen(trim(tempfile)//'.ppm/ppm')
              call pginitl(colour,file,whiteBG)
           end if
           if(file.ge.2) then
              write(tempfile,'(A)') trim(outputbasefile)
              if(.not.use_PLplot) io = pgopen(trim(tempfile)//'.eps'//trim(psclr))
              if(project_map) then
                 call pgpap(PSsz/0.5*PSrat,0.5)
              else
                 call pgpap(PSsz,PSrat)
              end if
              if(use_PLplot) io = pgopen(trim(tempfile)//'.eps'//trim(psclr))
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
        end if
        
        
        
        !*** Determine plot/binning ranges:
        xmin = minval(selDat(ic,p1,1:n(ic)))
        xmax = maxval(selDat(ic,p1,1:n(ic)))
        ymin = minval(selDat(ic,p2,1:n(ic)))
        ymax = maxval(selDat(ic,p2,1:n(ic)))
        
        dx = xmax - xmin
        dy = ymax - ymin
        !write(stdOut,'(A,2F10.5)')'  Xmin,Xmax: ',xmin,xmax
        !write(stdOut,'(A,2F10.5)')'  Ymin,Ymax: ',ymin,ymax
        
        if(.not.project_map) then
           xmin = xmin - 0.05*dx
           xmax = xmax + 0.05*dx
           ymin = ymin - 0.05*dy
           ymax = ymax + 0.05*dy
        end if
        
        
        
        
        
        
        if(plot.eq.1 .and. project_map) call prepare_skymap_binning(xmin,xmax, ymin,ymax)  ! Prepare binning for a 2D sky map
        
        
        
        !*** Bin data and 'normalise' 2D PDF:
        call bin_and_normalise_2D_data(ic,p1,p2, xmin,xmax, ymin,ymax, z,tr, sky_position,binary_orientation)
        
        
        ! Swap RA boundaries for RA-Dec plot in 2D PDF:
        if(sky_position) then
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
           
           
           call pgsch(sch)
           if(project_map .and. plotSky.ge.2) then
              call pgsvp(0.08*sch,0.95,0.08*sch,1.0-0.05*sch)   ! Make room for title and +90deg label
           else
              call pgsvp(0.08*sch,0.98,0.07*sch,1.0-0.033*sch)  ! Make room for title.
              ! Since sch is typically ~1.5*fontsize2d: 0.95 -> 1-0.05*fontsize ~ 1-0.03*sch
           end if
           
           
           call pgswin(xmin,xmax,ymin,ymax)
           if(project_map .and. (plotSky.eq.1.or.plotSky.eq.3).and.file.ge.2) then  ! Need dark background
              call pgsci(1)
              call pgrect(xmin,xmax,ymin,ymax)
           end if
           
           
           !*** Plot the actual 2D PDF (grey-scale or colour pixels):
           if(plPDF2D.eq.1.or.plPDF2D.eq.2) call plot_2D_PDF(z, tr, project_map)
           
           
           
           ! Plot stars in 2D PDF (over the grey scales, but underneath contours, lines, etc):
           if(project_map .and. (plotSky.eq.1.or.plotSky.eq.3)) then
              call pgswin(xmin*15,xmax*15,ymin,ymax)  ! Map works in degrees
              call plotthesky(xmin*15,xmax*15,ymin,ymax,raShift*rh2r)
              call pgswin(xmin,xmax,ymin,ymax)
           end if
           call pgsci(1)
           
           
           
           
           !*** Plot contours in 2D PDF:
           if((plPDF2D.eq.1.or.plPDF2D.eq.3) .and. (.not.project_map .or. plotSky.eq.1.or.plotSky.eq.3)) &
                call plot_2D_contours(z, tr, project_map, lw)
           
        end if  ! if(plot.eq.1)
        
        
        
        !*** Save binned 2D PDF data:
        if(savePDF.eq.1) then
           write(30,'(A)')'--------------------------------------------------------------------------------------------------'// &
                '------------------------------------------------------------------------------------------------------'
           write(30,'(3I6,T26,2A15,T101,A)')ic,parID(p1),parID(p2),parNames(parID(p1)),parNames(parID(p2)), &
                'Chain number, parameter ID 1,2 and parameter names  (ic,parID(1:2),parNames(1:2))'
           write(30,'(2ES15.7,T101,A)')startval(ic,p1,1:2),'Injection and starting value p1  (startval(1,1:2))'
           write(30,'(2ES15.7,T101,A)')startval(ic,p2,1:2),'Injection and starting value p2  (startval(2,1:2))'
           write(30,'(6ES15.7,T101,A)')stats(ic,p1,1:6), &
                'Stats: median, mean, absVar1, absVar2, stdev1, stdev2 for p1  (stats(1,1:6))'
           write(30,'(6ES15.7,T101,A)')stats(ic,p2,1:6), &
                'Stats: median, mean, absVar1, absVar2, stdev1, stdev2 for p2  (stats(2,1:6))'
           write(30,'(5ES15.7,T101,A)')ranges(ic,c0,p1,1:5), &
                'Ranges: lower,upper limit, centre, width, relative width for p1  (ranges(1,1:5))'
           write(30,'(5ES15.7,T101,A)')ranges(ic,c0,p2,1:5), &
                'Ranges: lower,upper limit, centre, width, relative width for p2  (ranges(1,1:5))'
           write(30,'(4ES15.7,T101,A)')xmin,xmax,ymin,ymax,'Xmin,Xmax,Ymin,Ymax of PDF  (xmin,xmax,ymin,ymax)'
           write(30,'(6ES15.7,T101,A)')tr,'Tr; transformation matrix used by PGPlot to project data  (tr)'
           write(30,'(A)')'  2D bins:'
           do i=1,Nbin2Dx+1
              do j=1,Nbin2Dy+1
                 write(30,'(ES15.7)',advance="no") z(i,j)
              end do
              write(30,'(1x)')
           end do
        end if
        
        
        
        
        if(plot.eq.1) then
           
           !*** Plot injection value, median, ranges, etc. in 2D PDF:
           call plot_values_in_2D_PDF(ic, p1,p2, xmin,xmax, ymin,ymax, dx,dy, sch,lw, project_map)
           
           
           !*** Print axes, axis labels and plot title:
           call plot_2D_PDF_axes_labels_titles(p1,p2, sch,flw, project_map)
           
           
           !*** Finish the current plot:
           countplots = countplots + 1  ! The current plot is number countplots
           
           
           !*** Convert plot:
           if(file.eq.1) then
              call pgend
              
              inquire(file=trim(tempfile)//'.ppm', exist=ex)
              if(ex) then
                 convopts = '-resize '//trim(bmpxpix)//' -depth 8 -unsharp '//trim(unSharppdf2d)
                 if(countplots.eq.Npdf2D) then
                    
                    ! Convert the last plot in the foreground, so that the process finishes before deleting the original file:
                    status = system('convert '//trim(convopts)//' '//trim(tempfile)//'.ppm '//trim(outputbasefile)//'.png')
                    
                 else  ! in the background
                    
                    status = system('convert '//trim(convopts)//' '//trim(tempfile)//'.ppm '//trim(outputbasefile)//'.png &')
                    
                 end if
                 
                 if(status.ne.0) write(stdErr,'(A)')'  Error converting plot for '//trim(parNames(parID(p1)))//'-'// &
                      trim(parNames(parID(p2)))
              end if
           end if
           
           
           if(file.ge.2) then
              call pgend
              if(file.eq.3) then
                 status = system('eps2pdf '//trim(tempfile)//'.eps &> /dev/null')
                 if(status.ne.0) write(stdErr,'(A)')'  Error converting plot for '//trim(parNames(parID(p1)))//'-'// &
                      trim(parNames(parID(p2)))
              end if
           end if
           
        end if  ! if(plot.eq.1)
        
        
     end do  ! p2
  end do  ! p1
  
  
  
  
  if(savePDF.eq.1) close(30)
  
  if(plot.eq.1) then
     if(file.ne.1) call pgend
     !if(file.ge.2.and.multipagefile) then
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
     
     ! Remove all the .ppm files:
     if(file.eq.1) then
        countplots = 0
        if(html.ge.1) write(51,'(A)')'<table>'
        do p1=j1,j2
           if(html.ge.1) write(51,'(A)')'<tr>'
           do p2=j1,j2
              
              
              if(Npdf2D.ge.0) then
                 plotthis = 0  ! Determine to plot or save this combination of j1/j2 or p1/p2
                 do i=1,Npdf2D
                    ! Use PDF2Dpairs from the input file:
                    if(p1.eq.revID(PDF2Dpairs(i,1)).and.p2.eq.revID(PDF2Dpairs(i,2))) plotthis = 1  
                 end do
                 if(plotthis.eq.0) cycle
              else
                 if(p2.le.p1) then
                    if(html.ge.1) then
                       if(p1.eq.p2) then
                          write(51,'(A)')'<td></td>'
                       else
                          write(outputbasefile,'(A)') trim(outputname)//'__pdf2d__'// &
                               trim(parNames(parID(p2)))//'-'//trim(parNames(parID(p1)))
                          write(51,'(A)')'<td>'
                          write(51,'(A)')'  <a href="'//trim(outputbasefile)//'.png">'
                          write(51,'(A)')'    <img src="'//trim(outputbasefile)//'_thumb.png">'
                          write(51,'(A)')'  </a>'
                          write(51,'(A)')'</td>'
                       end if
                    end if
                    cycle
                 end if
                 if(fixedpar(p1)+fixedpar(p2).ge.1) cycle
              end if
              
              countplots = countplots + 1  ! The current plot is number countplots
              write(outputbasefile,'(A)') trim(outputname)//'__pdf2d__'//trim(parNames(parID(p1)))//'-'//trim(parNames(parID(p2)))
              write(tempfile,'(A)') trim(outputdir)//'/'//trim(outputbasefile)
              status = system('rm -f '//trim(tempfile)//'.ppm')
              
              if(html.ge.1) then
                 
                 inquire(file=trim(tempfile)//'.png', exist=ex)
                 if(ex) then
                    
                    ! Convert the last plot in the foreground, so that the process finishes before deleting the original file:
                    if(countplots.eq.Npdf2D) then
                       status = system('convert -resize 200x200 '//trim(tempfile)//'.png '//trim(tempfile)//'_thumb.png')
                    else
                       status = system('convert -resize 200x200 '//trim(tempfile)//'.png '//trim(tempfile)//'_thumb.png &')
                    end if
                    if(status.ne.0) write(stdErr,'(A)')'  Error creating thumbnail for '//trim(parNames(parID(p1)))//'-'// &
                         trim(parNames(parID(p2)))
                 end if
                 
                 write(51,'(A)')'<td>'
                 write(51,'(A)')'  <a href="'//trim(outputbasefile)//'.png">'
                 write(51,'(A)')'    <img src="'//trim(outputbasefile)//'_thumb.png">'
                 write(51,'(A)')'  </a>'
                 write(51,'(A)')'</td>'
              end if
              
           end do
           if(html.ge.1) write(51,'(A)')'</tr>'
        end do
        if(html.ge.1) write(51,'(A)')'</table>'
     end if
     
     if(html.eq.1) then
        bmpXSz = 1000
        bmpYSz =  700
        
        call compBitmapSize(bmpXSz,bmpYSz, scFac, bmpsz,bmprat)  ! Determine plot size and ratio
        write(bmpxpix,'(I4)')bmpXSz  ! Used as a text string by convert
        pltsz = bmpsz
        pltrat = bmprat
     end if
     
  end if !plot.eq.1
  
  
  
end subroutine pdfs2d
!***********************************************************************************************************************************





