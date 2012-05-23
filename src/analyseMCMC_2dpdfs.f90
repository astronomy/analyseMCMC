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
  use SUFR_system, only: warn, swapreal
  use SUFR_text, only: replace_substring
  
  use aM_constants, only: use_PLplot
  use analysemcmc_settings, only: update,prProgress,file,scrsz,scrrat,pssz,psrat,fonttype,colour,whitebg,quality
  use analysemcmc_settings, only: fontsize2d,map_projection,maxChs
  use analysemcmc_settings, only: Npdf2D,PDF2Dpairs,html,bmpXSz,bmpYSz,scFac,Nbin2Dx,Nbin2Dy,plotSky
  use analysemcmc_settings, only: savePDF,plot,ivals,Nival,normPDF2D,plPDF1D,plPDF2D,prIval
  use general_data, only: outputname,outputdir,startval,parNames,pgParNs,pgUnits, nfixedpar
  use general_data, only: selDat,stats,ranges,c0,n,maxIter,fixedpar, raCentre,raShift
  use mcmcrun_data, only: totpts,revID,parID, nMCMCpar
  use plot_data, only: psclr,bmpsz,bmprat,bmpxpix,unSharppdf2d,pltsz,pltrat
  use stats_data, only: probAreas
  
  implicit none
  integer, intent(out) :: exitcode
  
  real,allocatable :: z(:,:),zs(:,:,:)  ! These depend on nbin2d, allocate after reading input file
  
  integer :: i,j,j1,j2,p1,p2,ic,lw,io,c,status,system,pgopen,clr,maxclr
  integer :: npdf,flw,plotthis,countplots,totplots, clr1,clr2
  real :: tr(6), sch, just
  real :: x,xmin,xmax, ymin,ymax, dx,dy, area
  character :: string*(99),str*(99),tempfile*(99),ivalstr*(99), outputbasefile*(199), convopts*(99), areaunit*(19)
  character :: tmpStr1*(19),tmpStr2*(19)
  logical :: project_map,sky_position,binary_orientation, ex
  !real :: xmin1,xmax1,ymin1,ymax1
  
  
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
        else
           if(p2.le.p1) cycle
           if(fixedpar(p1)+fixedpar(p2).ge.1) cycle
           if(stdOut.lt.10) then
              write(stdOut,*)cursorup  ! Move cursor up 1 line
              if(prProgress.ge.1.and.update.eq.0) write(stdOut,'(F7.1,A)')real(countplots+1)/real(totplots)*100, &
                   '%    ('//trim(parNames(parID(p1)))//'-'//trim(parNames(parID(p2)))//')                                      '
           end if
        end if
        
        
        ! Identify special combinations of parameters:
        sky_position = .false.
        binary_orientation = .false.
        project_map = .false.
        if(parID(p1).eq.31.and.parID(p2).eq.32) sky_position = .true.
        if(parID(p1).eq.31.and.parID(p2).eq.33) sky_position = .true.
        if(parID(p1).eq.52.and.parID(p2).eq.51) binary_orientation = .true.
        
        ! Make a special sky plot (i.e., plot stars or use projection) if plotSky>0 and RA,Dec are plotted:
        if(sky_position .and. plotSky.ge.1) project_map = .true.
        
        
        
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
        
        
        
        
        
        
        if(plot.eq.1 .and. project_map) call prepare_skymap_binning()  ! Prepare binning for a cute sky map in 2D PDF
        
        
        
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
           
           
           ! Force plotting boundaries (not binning boundaries):
           if(1.eq.2.and.sky_position) then
              xmin = 24.
              xmax = 0.
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
           
           
           !*** Plot the actual 2D PDF (grey scales or colour):
           if(plPDF2D.eq.1.or.plPDF2D.eq.2) then
              
              ! Set the colour schemes:
              if(normPDF2D.lt.4) then  ! Grey scales
                 call pgscir(0,nint(1e9))
                 call pgqcir(clr,maxclr)  ! Maxclr is device-dependent
                 if(maxclr.lt.30) call warn('Not enough colours on device for 2D plot!',0)
                 do i=0,maxclr-30  ! Colour indices typically run 0-255, but this is device-dependent. 
                    ! Reserve ~0-29 for other purposes -> (maxclr-30) for these grey scales:
                    x = real((maxclr-30) - i)/real(maxclr-30)          ! White background
                    call pgscr(30+i,x,x,x)
                 end do
                 call pgscir(30,maxclr)  ! Set colour-index range for pgimag
              else if(normPDF2D.eq.4) then
                 call set_2D_probability_colours(clr1, clr2)  ! Define the colours for the 2D probability areas
              end if
              
              
              ! Plot the PDF:
              if(project_map .and. plotSky.ge.2) then
                 if(prProgress.ge.3) write(stdOut,'(A)',advance="no")'  plotting map projection...'
                 call pgimag_project(z, Nbin2Dx+1, Nbin2Dy+1, 1,Nbin2Dx+1, 1,Nbin2Dy+1, 0.,1., clr1,clr2, tr, map_projection)
              else
                 if(prProgress.ge.3) write(stdOut,'(A)',advance="no")'  plotting 2D PDF...'
                 
                 ! Plot 2D image - 0: no projection:
                 !call pgimag_project(z, Nbin2Dx+1, Nbin2Dy+1, 1,Nbin2Dx+1, 1,Nbin2Dy+1, 0.,1., clr1,clr2, tr, 0)
                 
                 ! Plot 2D image - produces ~2.5x smaller plots - used to give segfaults:
                 call pgimag(z,Nbin2Dx+1,Nbin2Dy+1,1,Nbin2Dx+1,1,Nbin2Dy+1,0.,1.,tr)
              end if
              
           end if  !if(plPDF2D.eq.1.or.plPDF2D.eq.2)
           
           
           
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
           
           ! Plot injection value, median, ranges, etc. in 2D PDF:
           call plot_values_in_2D_PDF(ic, p1,p2, xmin,xmax, ymin,ymax, dx,dy, sch,lw, project_map)
           
           
           
           
           
           
           
           
           !*** Print axes, axis labels and plot title:
           
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
                 call replace_substring(string, '\(727)10\u-', '\(727)10\u\(240)')
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
              end do
              call pgsci(1)
           end if  ! if(prIval.ge.1.and.normPDF2D.eq.4)
           
           
           
           !*** Finish the current plot:
           countplots = countplots + 1  ! The current plot is number countplots
           
           ! Convert plot:
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
           
           !if(file.ge.2.and.multipagefile) call pgpage
           if(file.ge.2) then
              call pgend
              if(file.eq.3) then
                 status = system('eps2pdf '//trim(tempfile)//'.eps &> /dev/null')
                 if(status.ne.0) write(stdErr,'(A)')'  Error converting plot for '//trim(parNames(parID(p1)))//'-'// &
                      trim(parNames(parID(p2)))
              end if
           end if
        end if !if(plot.eq.1)
        
        
     end do !p2
  end do !p1
  
  
  
  
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





