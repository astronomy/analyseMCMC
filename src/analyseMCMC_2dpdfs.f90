!> \file analyseMCMC_2dpdfs.f90  Routines to compute and plot two-dimensional PDFs


!***********************************************************************************************************************************
!> \brief  Plot 2D marginalised PDFs
!!
!! \retval exitcode  Exit code: 0=ok

subroutine pdfs2d(exitcode)
  use SUFR_constants, only: stdOut,stdErr
  use SUFR_constants, only: cursorup, pi,rpi,rh2r
  use SUFR_system, only: warn
  
  use analysemcmc_settings, only: update,prProgress,file,scrsz,scrrat,pssz,psrat,fonttype,colour,whitebg,quality
  use analysemcmc_settings, only: plLmax,fontsize2d,map_projection,maxChs
  use analysemcmc_settings, only: plInject,mergeChains,Npdf2D,PDF2Dpairs,html,bmpXSz,bmpYSz,scFac,Nbin2Dx,Nbin2Dy,plotSky,wrapData
  use analysemcmc_settings, only: savePDF,plot,ivals,Nival,normPDF2D,plPDF1D,plPDF2D,plMedian,plRange,prIval
  use general_data, only: allDat,outputname,outputdir,startval,icloglmax,iloglmax,parNames,pgParNs,nfixedpar
  use general_data, only: selDat,stats,ranges,c0,n,maxIter,wrap,fixedpar,shifts,shIvals,raCentre,raShift
  use mcmcrun_data, only: totpts,revID,parID, nMCMCpar
  use plot_data, only: psclr,bmpsz,bmprat,bmpxpix,unSharppdf2d,pltsz,pltrat
  use stats_data, only: probArea,probAreas,injectionranges2d
  
  implicit none
  integer, intent(out) :: exitcode
  
  integer :: i,j,j1,j2,p1,p2,ic,lw,io,c,status,system,pgopen,clr,maxclr
  integer :: npdf,ncont,flw,plotthis,injectionrange2d,countplots,totplots, clr1,clr2
  real :: a,rat,cont(11),tr(6),sch,plx,ply
  real :: x,xmin,xmax,ymin,ymax,dx,dy,xx(maxChs*maxIter),yy(maxChs*maxIter),zz(maxChs*maxIter)
  real,allocatable :: z(:,:),zs(:,:,:)  !These depend on nbin2d, allocate after reading input file
  character :: string*(99),str*(99),tempfile*(99),ivalstr*(99),delta*(19),outputbasefile*(199), convopts*(99)
  logical :: project_map,sky_position,binary_orientation, ex
  !real :: xmin1,xmax1,ymin1,ymax1
  
  
  exitcode = 0
  countplots = 0
  ic = 1 !Can only do one chain
  
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
     
     bmpsz = real(bmpXSz-1)/85. * scFac  ! Make png larger, so that convert interpolates and makes the plot smoother
     bmprat = real(bmpYSz-1)/real(bmpXSz-1)
     write(bmpxpix,'(I4)')bmpXSz  ! Used as a text string by convert
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
  if(Nbin2Dy.le.-1) Nbin2Dy = nint(Nbin2Dx*pltrat)
  
  ! Report number of bins used:
  if(prProgress.ge.1.and.plot.eq.1.and.update.eq.0.and.Npdf2D.ge.0) then
     if(Nbin2Dx.lt.100) then
        write(stdOut,'(A1,I2,A1)',advance="no")'(',Nbin2Dx,'x'
     else
        write(stdOut,'(A1,I3,A1)',advance="no")'(',Nbin2Dx,'x'
     end if
     if(Nbin2Dy.lt.100) then
        write(stdOut,'(I2,A7)',advance="no")Nbin2Dy,' bins) '
     else
        write(stdOut,'(I3,A7)',advance="no")Nbin2Dy,' bins) '
     end if
  end if
  
  ! Allocate memory:
  allocate(z(Nbin2Dx+1,Nbin2Dy+1),zs(maxChs,Nbin2Dx+1,Nbin2Dy+1))
  
  
  if(plot.eq.1) then
     if(file.eq.0) then
        lw = 1
        flw = nint(1*fontsize2d)  ! Font lw
        sch = 1.5*fontsize2d
     end if
     if(file.ge.1) then
        !if(file.ge.2.and.multipagefile) io = pgopen(trim(outputdir)//'/pdf2d.eps'//trim(psclr))
        lw = 3
        flw = nint(2*fontsize2d)  ! Font lw
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
           write(outputbasefile,'(A)') trim(outputdir)//'/'//trim(outputname)//'__pdf2d__'// &
                trim(parNames(parID(p1)))//'-'//trim(parNames(parID(p2)))
           
           if(file.eq.0) then
              npdf=npdf+1
              write(str,'(I3,A3)')200+npdf,'/xs'
              io = pgopen(trim(str))
              if(project_map) then
                 call pgpap(scrSz/0.5*scrRat,0.5)
              else
                 call pgpap(scrSz,scrRat)
              end if
              call pginitl(colour,file,whiteBG)
           end if
           if(file.eq.1) then
              write(tempfile,'(A)') trim(outputbasefile)
              io = pgopen(trim(tempfile)//'.ppm/ppm')
              if(project_map) then
                 call pgpap(bmpsz/0.5*bmprat,0.5)
              else
                 call pgpap(bmpsz,bmprat)
              end if
              call pginitl(colour,file,whiteBG)
           end if
           if(file.ge.2) then
              write(tempfile,'(A)') trim(outputbasefile)
              io = pgopen(trim(tempfile)//'.eps'//trim(psclr))
              if(project_map) then
                 call pgpap(PSsz/0.5*PSrat,0.5)
              else
                 call pgpap(PSsz,PSrat)
              end if
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
        
        xmin = minval(selDat(ic,p1,1:n(ic)))
        xmax = maxval(selDat(ic,p1,1:n(ic)))
        ymin = minval(selDat(ic,p2,1:n(ic)))
        ymax = maxval(selDat(ic,p2,1:n(ic)))
        dx = xmax - xmin
        dy = ymax - ymin
        !write(stdOut,'(A,2F10.5)')'  Xmin,Xmax: ',xmin,xmax
        !write(stdOut,'(A,2F10.5)')'  Ymin,Ymax: ',ymin,ymax
        
        xx(1:n(ic)) = selDat(ic,p1,1:n(ic))  ! Parameter 1
        yy(1:n(ic)) = selDat(ic,p2,1:n(ic))  ! Parameter 2
        zz(1:n(ic)) = selDat(ic,1,1:n(ic))   ! Likelihood
        
        if(.not.project_map) then
           xmin = xmin - 0.05*dx
           xmax = xmax + 0.05*dx
           ymin = ymin - 0.05*dy
           ymax = ymax + 0.05*dy
        end if
        
        
        
        
        
        ! Prepare binning for a cute sky map in 2D PDF:
        if(plot.eq.1 .and. project_map) then
           rat = 0.5 !scrRat !0.75
           !call pgpap(11.,rat)
           !call pgpap(scrSz,scrRat)  ! This causes a 'pgpage' when pggray is called...
           dx = xmax - xmin
           dy = ymax - ymin
           if(abs(dx)*15.lt.dy/rat) then  ! Expand x
              dx = dy/(15*rat)
              a = (xmin+xmax)*0.5
              xmin = a - 0.5*dx
              xmax = a + 0.5*dx
              if(prProgress.ge.3) write(stdOut,'(A,F6.1,A3,F6.1,A)',advance="no")'  Changing RA binning range to ' &
                   ,xmin,' - ',xmax,' h.'
           end if
           if(abs(dx)*15.gt.dy/rat) then  ! Expand y
              dy = abs(dx)*rat*15
              a = (ymin+ymax)*0.5
              ymin = a - 0.5*dy
              ymax = a + 0.5*dy
              if(prProgress.ge.3) write(stdOut,'(A,F6.1,A3,F6.1,A)',advance="no")'  Changing declination binning range to ' &
                   ,ymin,' - ',ymax,' deg.'
           end if
        end if !if(plot.eq.1 .and. project_map)
        
        ! Force plotting and binning boundaries:
        ! CHECK: lose this? - yes: these are the binning ranges, not the plotting ranges; 
        ! Don't necessarily want to bin the whole sky when plotting it.
        if(1.eq.2.and.wrapData.eq.0.and.sky_position) then
           xmin = 0.
           xmax = 24.
           ymin = -90.
           ymax = 90.
        end if
        
        
        ! Bin data and 'normalise' 2D PDF:
        if(normPDF2D.le.2.or.normPDF2D.eq.4) then
           
           ! Bin data:  compute bin number rather than find it, ~10x faster:
           call bindata2d(n(ic),xx(1:n(ic)),yy(1:n(ic)),0,Nbin2Dx,Nbin2Dy,xmin,xmax,ymin,ymax,z,tr)
           
           !Test
           !call check_binned_data(Nbin2Dx,Nbin2Dy,z)
           
           !do Nbin2Dx = 10,200,10
           !   Nbin2Dy = Nbin2Dx
           !   xmin1 = xmin
           !   xmax1 = xmax
           !   ymin1 = ymin
           !   ymax1 = ymax
           !   
           !   !Bin data:  compute bin number rather than find it, ~10x faster:
           !   call bindata2d(n(ic),xx(1:n(ic)),yy(1:n(ic)),0,Nbin2Dx,Nbin2Dy,xmin1,xmax1,ymin1,ymax1,z,tr)
           !   
           !   !Test!
           !   call check_binned_data(Nbin2Dx,Nbin2Dy,z)
           !   
           !end do
           !stop
           
           
           if(normPDF2D.eq.1) z = max(0.,log10(z + 1.e-30))
           if(normPDF2D.eq.2) z = max(0.,sqrt(z + 1.e-30))
           
           if(normPDF2D.eq.4) then
              
              ! Get 2D probability ranges; identify to which range each bin belongs:
              if(prProgress.ge.3) write(stdOut,'(A)',advance="no")'  identifying 2D ranges...'
              call identify_2d_ranges(p1,p2,Nival,Nbin2Dx+1,Nbin2Dy+1,z,tr)
              
              ! Compute 2D probability areas; sum the areas of all bins:
              if(prProgress.ge.3) write(stdOut,'(A)',advance="no")'  computing 2D areas...'
              call calc_2d_areas(p1,p2,Nival,Nbin2Dx+1,Nbin2Dy+1,z,tr,probArea)
              injectionranges2d(p1,p2) = injectionrange2d(z,Nbin2Dx+1,Nbin2Dy+1,startval(1,p1,1),startval(1,p2,1),tr)
              
              do i=1,Nival
                 if(prIval.ge.1.and.prProgress.ge.2 .and. (sky_position .or. binary_orientation)) then  
                    ! For sky position and orientation only:
                    if(i.eq.1) write(stdOut,'(/,1x,A10,A13,3A23)')'Nr.','Ival frac.','Area (sq.deg) ', &
                         'Circ. area rad. (deg) ','Fraction of sky '
                    write(stdOut,'(I10,F13.2,3(2x,F21.5))')i,ivals(i),probArea(i),sqrt(probArea(i)/pi)*2, &
                         probArea(i)*(pi/180.)**2/(4*pi)  ! 4pi*(180/pi)^2 = 41252.961 sq. degrees in a sphere
                 end if
                 probAreas(p1,p2,i,1) = probArea(i)*(rpi/180.)**2/(4*rpi)  ! Fraction of the sky
                 probAreas(p1,p2,i,2) = sqrt(probArea(i)/rpi)*2            ! Equivalent diameter
                 probAreas(p1,p2,i,3) = probArea(i)                        ! Square degrees
              end do
           end if
        end if
        if(normPDF2D.eq.3) then  ! Weigh by likelihood value
           if(prProgress.ge.3) write(stdOut,'(A)',advance="no")'  binning 2D data...'
           ! Measure amount of likelihood in each bin:
           call bindata2da(n(ic),xx(1:n(ic)),yy(1:n(ic)),zz(1:n(ic)),0,Nbin2Dx,Nbin2Dy,xmin,xmax,ymin,ymax,z,tr)
        end if
        
        
        
        
        ! Swap RA boundaries for RA-Dec plot in 2D PDF:
        if(sky_position) then
           a = xmin
           xmin = xmax
           xmax = a
           dx = -dx
        end if
        
        z = z/(maxval(z)+1.e-30)
        
        
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
              call pgsvp(0.08*sch,0.95,0.08*sch,1.0-0.033*sch)  ! Make room for title.  
              ! Since sch is typically ~1.5*fontsize2d: 0.95 -> 1-0.05*fontsize ~ 1-0.03*sch
           end if
           
           
           call pgswin(xmin,xmax,ymin,ymax)
           if(project_map .and. (plotSky.eq.1.or.plotSky.eq.3).and.file.ge.2) then  ! Need dark background
              call pgsci(1)
              call pgrect(xmin,xmax,ymin,ymax)
           end if
           
           ! Plot the actual 2D PDF (grey scales or colour):
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
              end if
              
              if(normPDF2D.eq.4) then  ! Colour
                 if(colour.eq.0) then
                    call pgscr(30,1.,1.,1.)  ! BG colour
                    if(Nival.eq.2) then
                       call pgscr(31,0.5,0.5,0.5)  ! Grey
                       call pgscr(32,0.,0.,0.)     ! Black
                    end if
                    if(Nival.eq.3) then
                       call pgscr(31,0.7,0.7,0.7)  ! Light grey
                       call pgscr(32,0.4,0.4,0.4)  ! Dark grey
                       call pgscr(33,0.0,0.0,0.0)  ! Black
                    end if
                    if(Nival.eq.4) then
                       call pgscr(31,0.75,0.75,0.75)  ! Light grey
                       call pgscr(32,0.50,0.50,0.50)  ! 
                       call pgscr(33,0.25,0.25,0.25)  ! Dark grey
                       call pgscr(34,0.00,0.00,0.00)  ! Black
                    end if
                    if(Nival.eq.5) then
                       call pgscr(31,0.8,0.8,0.8)  ! Light grey
                       call pgscr(32,0.6,0.6,0.6)  ! 
                       call pgscr(33,0.4,0.4,0.4)  ! 
                       call pgscr(34,0.2,0.2,0.2)  ! Dark grey
                       call pgscr(35,0.0,0.0,0.0)  ! Black
                    end if
                 end if
                 if(colour.ge.1) then
                    call pgscr(30,1.,1.,1.)  ! BG colour
                    if(Nival.eq.2) then
                       call pgscr(31,1.,1.,0.)  ! Yellow
                       if(file.ge.2) call pgscr(31,0.8,0.7,0.)  ! Dark yellow
                       call pgscr(32,1.,0.,0.)  ! Red
                    end if
                    if(Nival.eq.3) then
                       call pgscr(31,0.,0.,1.)  ! Blue
                       call pgscr(32,1.,1.,0.)  ! Yellow
                       if(file.ge.2) call pgscr(32,0.8,0.7,0.)  ! Dark yellow
                       call pgscr(33,1.,0.,0.)  ! Red
                    end if
                    if(Nival.eq.4) then
                       call pgscr(31,0.,0.,1.)  ! Blue
                       call pgscr(32,0.,1.,0.)  ! Green
                       call pgscr(33,1.,1.,0.)  ! Yellow
                       if(file.ge.2) call pgscr(33,0.8,0.7,0.)  ! Dark yellow
                       call pgscr(34,1.,0.,0.)  ! Red
                    end if
                    if(Nival.eq.5) then
                       call pgscr(31,0.,0.,1.)  ! Blue
                       call pgscr(32,0.,1.,0.)  ! Green
                       call pgscr(33,1.,1.,0.)  ! Yellow
                       if(file.ge.2) call pgscr(33,0.8,0.7,0.)  ! Dark yellow
                       call pgscr(34,1.,0.5,0.)  ! Orange
                       call pgscr(35,1.,0.,0.)  ! Red
                    end if
                 end if
                 clr1 = 30
                 clr2 = 30+Nival
                 call pgscir(clr1,clr2)  ! set colour-index range for pgimag
              end if  !if(normPDF2D.eq.4)
              
              
              
              
              ! Plot the PDF:
              if(project_map .and. plotSky.ge.2) then
                 if(prProgress.ge.3) write(stdOut,'(A)',advance="no")'  plotting map projection...'
                 call pgimag_project(z, Nbin2Dx+1, Nbin2Dy+1, 1,Nbin2Dx+1, 1,Nbin2Dy+1, 0.,1., clr1,clr2, tr, map_projection)
              else
                 if(prProgress.ge.3) write(stdOut,'(A)',advance="no")'  plotting 2D PDF...'
                 
                 ! CHECK Gives seg.fault (on amd64) - but why? - need more memory?
                 !call pgimag(z,Nbin2Dx+1,Nbin2Dy+1,1,Nbin2Dx+1,1,Nbin2Dy+1,0.,1.,tr)
                 
                 call pgimag_project(z, Nbin2Dx+1, Nbin2Dy+1, 1,Nbin2Dx+1, 1,Nbin2Dy+1, 0.,1., clr1,clr2, tr, 0)  ! 0-no projection
              end if
              
           end if  !if(plPDF2D.eq.1.or.plPDF2D.eq.2)
           
           
           
           ! Plot stars in 2D PDF (over the grey scales, but underneath contours, lines, etc):
           if(project_map .and. (plotSky.eq.1.or.plotSky.eq.3)) then
              call pgswin(xmin*15,xmax*15,ymin,ymax)  ! Map works in degrees
              call plotthesky(xmin*15,xmax*15,ymin,ymax,raShift*rh2r)
              call pgswin(xmin,xmax,ymin,ymax)
           end if
           call pgsci(1)
           
        end if !if(plot.eq.1)
        
        
        
        ! Plot contours in 2D PDF:
        if(plot.eq.1 .and. (plPDF2D.eq.1.or.plPDF2D.eq.3) .and. (.not.project_map .or. plotSky.eq.1.or.plotSky.eq.3)) then
           if(normPDF2D.lt.4) then
              ncont = 11
              do i=1,ncont
                 cont(i) = 0.01 + 2*real(i-1)/real(ncont-1)
                 if(project_map .and. (plotSky.eq.1.or.plotSky.eq.3)) cont(i) = 1.-cont(i)
              end do
              ncont = min(4,ncont)  ! Only use the first 4
           end if
           if(normPDF2D.eq.4) then
              ncont = Nival
              do i=1,ncont
                 cont(i) = max(1. - real(i-1)/real(ncont-1),0.001)
                 !if(project_map) cont(i) = 1.-cont(i)
              end do
           end if
           
           call pgsls(1)
           if((.not.project_map .or. plotSky.ne.1.or.plotSky.ne.3) .and. normPDF2D.ne.4) then  ! First in bg colour
              call pgslw(2*lw)
              call pgsci(0)
              call pgcont(z,Nbin2Dx+1,Nbin2Dy+1,1,Nbin2Dx+1,1,Nbin2Dy+1,cont(1:ncont),ncont,tr)
           end if
           call pgslw(lw)
           call pgsci(1)
           if(project_map .and. (plotSky.eq.1.or.plotSky.eq.3)) call pgsci(7)
           call pgcont(z,Nbin2Dx+1,Nbin2Dy+1,1,Nbin2Dx+1,1,Nbin2Dy+1,cont(1:ncont),ncont,tr)
        end if
        
        
        ! Save binned 2D PDF data:
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
                 write(30,'(ES15.7)',advance="no")z(i,j)
              end do
              write(30,'(1x)')
           end do
        end if
        
        
        
        ! Plot injection value, median, ranges, etc. in 2D PDF:
        if(plot.eq.1) then
           if(.not.project_map.or.plotSky.eq.1) then
              call pgsci(1)
              
              ! Plot max likelihood in 2D PDF:
              if(plLmax.ge.1) then
                 call pgsci(1); call pgsls(5)
                 
                 plx = allDat(icloglmax,p1,iloglmax)
                 if(wrap(ic,p1).ne.0) plx = mod(plx + shifts(ic,p1), shIvals(ic,p1)) - shifts(ic,p1)
                 call pgline(2,(/plx,plx/),(/ymin,ymax/))  ! Max logL
                 
                 ply = allDat(icloglmax,p2,iloglmax)
                 if(wrap(ic,p2).ne.0) ply = mod(ply + shifts(ic,p2), shIvals(ic,p2)) - shifts(ic,p2)
                 call pgline(2,(/xmin,xmax/),(/ply,ply/))  ! Max logL
                 
                 call pgpoint(1,plx,ply,12)
              end if
              
              
              if(project_map .and. (plotSky.eq.1.or.plotSky.eq.3)) call pgsci(0)
              call pgsls(2)
              
              ! Plot injection value in 2D PDF:
              if((plInject.eq.1.or.plInject.eq.3).and.(.not.project_map) .or. &
                   ((plInject.eq.2.or.plInject.eq.4) .and.  &
                   (parID(p1).eq.61.and.parID(p2).eq.62).or.(parID(p1).eq.63.and.parID(p2).eq.64)) ) then
                 
                 ! CHECK The units of the injection values haven't changed (e.g. from rad to deg) for ic>1 
                 ! (but they have for the starting values, why?)
                 if(mergeChains.ne.1.or.ic.le.1) then 
                    
                    ! x:
                    call pgsls(2); call pgsci(1)
                    
                    ! Dash-dotted line for injection value when Lmax line isn't plotted (should we do this always?):
                    if(plLmax.eq.0) call pgsls(3)  
                    plx = startval(ic,p1,1)
                    if(wrap(ic,p1).ne.0) plx = mod(plx + shifts(ic,p1), shIvals(ic,p1)) - shifts(ic,p1)
                    call pgline(2,(/plx,plx/),(/ymin,ymax/)) !Injection value
                    
                    ! y:
                    call pgsls(2); call pgsci(1)
                    
                    ! Dash-dotted line for injection value when Lmax line isn't plotted (should we do this always?):
                    if(plLmax.eq.0) call pgsls(3)  
                    ply = startval(ic,p2,1)
                    if(wrap(ic,p2).ne.0) ply = mod(ply + shifts(ic,p2), shIvals(ic,p2)) - shifts(ic,p2)
                    call pgline(2,(/xmin,xmax/),(/ply,ply/)) !Injection value
                    
                    call pgpoint(1,plx,ply,18)
                 end if
              end if  !If plotting injection values in 2D plot
              call pgsci(1)
              call pgsls(4)
              
              
              ! Plot starting values in 2D PDF:
              !call pgline(2,(/startval(ic,p1,2),startval(ic,p1,2)/),(/ymin,ymax/))
              !call pgline(2,(/xmin,xmax/),(/startval(ic,p2,2),startval(ic,p2,2)/))
              
              call pgsci(2)
              
              
              ! Plot interval ranges in 2D PDF:
              if(plRange.eq.2.or.plRange.eq.3.or.plRange.eq.5.or.plRange.eq.6) then
                 write(delta,'(A,I3.3,A)')'\(2030)\d',nint(ivals(c0)*100),'%\u'
                 if(nint(ivals(c0)*100).lt.100) write(delta,'(A,I2.2,A)')'\(2030)\d',nint(ivals(c0)*100),'%\u'
                 call pgsls(1)
                 call pgsch(sch*0.6)
                 call pgsah(1,45.,0.1)
                 a = 0.0166667*sch
                 call pgarro(ranges(ic,c0,p1,3),ymin+dy*a,ranges(ic,c0,p1,1),ymin+dy*a)
                 call pgarro(ranges(ic,c0,p1,3),ymin+dy*a,ranges(ic,c0,p1,2),ymin+dy*a)
                 a = 0.0333333*sch
                 call pgptxt(ranges(ic,c0,p1,3),ymin+dy*a,0.,0.5,trim(delta))
                 a = 0.0233333*sch
                 call pgarro(xmin+dx*a,ranges(ic,c0,p2,3),xmin+dx*a,ranges(ic,c0,p2,1))
                 call pgarro(xmin+dx*a,ranges(ic,c0,p2,3),xmin+dx*a,ranges(ic,c0,p2,2))
                 a = 0.01*sch
                 call pgptxt(xmin+dx*a,ranges(ic,c0,p2,3),90.,0.5,trim(delta))
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
              
              
           end if  !if(.not.project_map.or.plotSky.eq.1)
           
           
           ! Plot big symbol at injection position in sky map:
           if(project_map .and. (plInject.eq.1.or.plInject.eq.3)) then
              call pgsch(sch*1.5)               ! Use 1.5 for plsym=8, 2 for plsym=18
              call pgslw(lw*2)
              call pgsci(9)
              if(normPDF2D.eq.4) call pgsci(1)  ! Black
              plx = startval(ic,p1,1)
              ply = startval(ic,p2,1)
              if(plotSky.eq.2.or.plotSky.eq.4) call project_skymap(plx,ply,raCentre,map_projection)
              call pgpoint(1,plx,ply,8)
              call pgsch(sch)
              call pgslw(lw)
              call pgsci(1)
           end if
           
        end if  !if(plot.eq.1)
        
        
        
        
        
        
        ! Print axes, axis labels and plot title:
        if(plot.eq.1) then
           
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
           call pgmtxt('B',2.2,0.5,0.5,trim(pgParNs(parID(p1))))
           call pgmtxt('L',1.7,0.5,0.5,trim(pgParNs(parID(p2))))
           
           
           ! Print 2D probability ranges in plot title:
           ! For sky position and orientation only:
           if(prIval.ge.1.and.normPDF2D.eq.4.and. (sky_position .or. binary_orientation)) then  
              string = ' '
              do c = 1,Nival
                 write(ivalstr,'(F5.1,A1)')ivals(c)*100,'%'
                 if(fonttype.eq.2) then
                    if(abs(ivals(c)-0.6827).lt.0.001) write(ivalstr,'(A)')'1\(2144)'
                    if(abs(ivals(c)-0.9545).lt.0.001) write(ivalstr,'(A)')'2\(2144)'
                    if(abs(ivals(c)-0.9973).lt.0.001) write(ivalstr,'(A)')'3\(2144)'
                 else
                    if(abs(ivals(c)-0.6827).lt.0.001) write(ivalstr,'(A)')'1\(0644)'
                    if(abs(ivals(c)-0.9545).lt.0.001) write(ivalstr,'(A)')'2\(0644)'
                    if(abs(ivals(c)-0.9973).lt.0.001) write(ivalstr,'(A)')'3\(0644)'
                 end if
                 
                 i = 3  ! 2-use degrees, 3-square degrees
                 a = probAreas(p1,p2,c,i)
                 if(i.eq.2) then
                    if(a.lt.1.) then
                       write(string,'(A,F5.2,A7)')trim(ivalstr)//':',a,'\(2218)'
                    else if(a.lt.10.) then
                       write(string,'(A,F4.1,A7)')trim(ivalstr)//':',a,'\(2218)'
                    else if(a.lt.100.) then
                       write(string,'(A,F5.1,A7)')trim(ivalstr)//':',a,'\(2218)'
                    else
                       write(string,'(A,I4,A7)')trim(ivalstr)//':',nint(a),'\(2218)'
                    end if
                 end if
                 if(i.eq.3) then
                    call pgsch(sch*0.85)  ! Needed to fit the square-degree sign in
                    if(quality.eq.3) call pgsch(sch*0.6)  ! Poster
                    if(a.lt.1.) then
                       write(string,'(A,F5.2,A9)')trim(ivalstr)//':',a,'deg\u2\d'
                    else if(a.lt.10.) then
                       write(string,'(A,F4.1,A9)')trim(ivalstr)//':',a,'deg\u2\d'
                    else if(a.lt.100.) then
                       write(string,'(A,F5.1,A9)')trim(ivalstr)//':',a,'deg\u2\d'
                    else if(a.lt.1000.) then
                       write(string,'(A,I4,A9)')trim(ivalstr)//':',nint(a),'deg\u2\d'
                    else if(a.lt.10000.) then
                       write(string,'(A,I5,A9)')trim(ivalstr)//':',nint(a),'deg\u2\d'
                    else
                       write(string,'(A,I6,A9)')trim(ivalstr)//':',nint(a),'deg\u2\d'
                    end if
                 end if
                 if(quality.eq.91) then  ! NINJA
                    call pgsch(sch)
                    if(fonttype.eq.2) then
                       write(string,'(I2,A7,A10)')c,'\(2144)',''
                    else
                       write(string,'(I2,A7,A10)')c,'\(0644)',''
                    end if
                 end if
                 a = (real(c-1)/real(Nival-1) - 0.5)*0.7 + 0.5
                 call pgsci(30+Nival+1-c)
                 if(project_map .and. plotSky.ge.2) then
                    call pgmtxt('T',1.0,a,0.5,trim(string))  ! Print title
                 else
                    call pgmtxt('T',0.5,a,0.5,trim(string))  ! Print title
                 end if
                 call pgsch(sch)
              end do
              call pgsci(1)
           end if
        end if  !if(plot.eq.1)
        
        
        
        if(plot.eq.1) then
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
        
        bmpsz = real(bmpXSz-1)/85. * scFac  ! Make png larger, so that convert interpolates and makes the plot smoother
        bmprat = real(bmpYSz-1)/real(bmpXSz-1)
        write(bmpxpix,'(I4)')bmpXSz  ! Used as a text string by convert
        pltsz = bmpsz
        pltrat = bmprat
     end if
     
  end if !plot.eq.1
  
  
  
end subroutine pdfs2d
!***********************************************************************************************************************************





