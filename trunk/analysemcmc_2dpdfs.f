!Routines to compute and plot two-dimensional PDFs


subroutine pdfs2d(exitcode)
  use constants
  use analysemcmc_settings
  use general_data
  use mcmcrun_data
  use plot_data
  use stats_data
  implicit none
  integer :: i,j,j1,j2,p1,p2,ic,lw,io,exitcode,c,system,pgopen,clr,maxclr
  integer :: npdf,ncont,lw2,plotthis,truerange2d,countplots,totplots
  real :: rev360,rev180,rev24
  real :: a,rat,cont(11),tr(6),sch,plx,ply
  real :: x,xmin,xmax,ymin,ymax,dx,dy,xx(nchs*narr1),yy(nchs*narr1),zz(nchs*narr1)
  real,allocatable :: z(:,:),zs(:,:,:)  !These depend on nbin2d, allocate after reading input file
  character :: string*99,str*99,tempfile*99,ivalstr*99
  logical :: project_map,sky_position,binary_orientation
  
  exitcode = 0
  countplots = 0
  ic = 1 !Can only do one chain
  
  !Columns in dat(): 1:logL 2:mc, 3:eta, 4:tc, 5:logdl, 6:spin, 7:kappa, 8: RA, 9:sindec,10:phase, 11:sinthJ0, 12:phiJ0, 13:alpha, 14:M1, 15:M2
  j1 = 2
  j2 = npar
  
  if(prprogress.ge.1.and.plot.eq.0.and.savepdf.eq.1.and.plpdf1d.eq.0) write(6,'(A,$)')'  Saving'
  if(prprogress.ge.1.and.update.eq.0.and.npdf2d.ge.0) write(6,'(A,$)')'  2D pdfs: '
  if(npdf2d.lt.0) then
     totplots = 0
     do i=j1,j2
        totplots = totplots + i - j1
     end do
     if(prprogress.ge.1.and.update.eq.0) then
        if(totplots.lt.100) write(6,'(A,I2,A,/)')'  *ALL* (',totplots,') 2D pdfs: '
        if(totplots.ge.100) write(6,'(A,I3,A,/)')'  *ALL* (',totplots,') 2D pdfs: '
     end if
  end if
  

  
  
  !Autodetermine number of bins for 2D PDFs:
  if(nbin2dx.le.0) then
     if(totpts.le.100) then
        nbin2dx = floor(2*sqrt(real(totpts))/pltrat)
        nbin2dx = max(nbin2dx,5)
        nbin2dy = floor(2*sqrt(real(totpts)))           !Same as for 1D case (~50)
     else
        nbin2dx = floor(10*log10(real(totpts))/pltrat)  
        nbin2dx = max(nbin2dx,5)
        nbin2dy = floor(10*log10(real(totpts)))         !Same as for 1D case (~50)
     end if
     if(prprogress.ge.2.and.plot.eq.1.and.update.eq.0) then
        if(nbin2dx.lt.100) write(6,'(A1,I2,A1,$)')'(',nbin2dx,'x'
        if(nbin2dx.ge.100) write(6,'(A1,I3,A1,$)')'(',nbin2dx,'x'
        if(nbin2dy.lt.100) write(6,'(I2,A7,$)')nbin2dy,' bins) '
        if(nbin2dy.ge.100) write(6,'(I3,A7,$)')nbin2dy,' bins) '
     end if
  end if
  if(nbin2dy.eq.0) nbin2dy = nbin2dx
  if(nbin2dy.le.-1) nbin2dy = nbin2dx*pltrat
  
  !Allocate memory:
  allocate(z(nbin2dx+1,nbin2dy+1),zs(nchs,nbin2dx+1,nbin2dy+1))
  
  if(plot.eq.1) then
     if(file.eq.0) then
        lw = 1
        lw2 = 1
        sch = 1.5*fontsize2d
     end if
     if(file.ge.1) then
        !if(file.ge.2.and.multipagefile) io = pgopen('pdf2d.eps'//trim(psclr))
        lw = 3
        lw2 = 2 !Font lw
        if(quality.eq.91) lw2 = 3 !NINJA
        sch = 1.5*fontsize2d
        if(quality.eq.3) then !Poster
           lw = 4
           lw2 = 3 !Font lw
           sch = 2.*fontsize2d
        end if
        if(pssz.lt.5) sch = sch * sqrt(5.0/pssz)
     end if
     !if(file.ge.2.and.multipagefile) then
     !   if(io.le.0) then
     !      write(0,'(A,I4)')'  Error:  cannot open PGPlot device.  Quitting the programme',io
     !      exitcode = 1
     !      return
     !   end if
     !   call pgscf(fonttype)
     !   call pginitl(colour,file,whitebg)
     !end if
  end if !if(plot.eq.1)
  
  if(savepdf.eq.1) then
     open(unit=30,action='write',form='formatted',status='replace',file=trim(outputdir)//'/'//trim(outputname)//'__pdf2d.dat')
     write(30,'(5I6,T100,A)')j1,j2,1,nbin2dx,nbin2dy,'Plot variable 1,2, total number of chains, number of bins x,y'
  end if
  
  npdf=0 !Count iterations to open windows with different numbers
  do p1=j1,j2
     do p2=j1,j2
        
        !Identify special combinations of parameters:
        sky_position = .false.
        binary_orientation = .false.
        project_map = .false.
        if(version.eq.1.and.p1.eq.8.and.p2.eq.9) sky_position = .true.
        if(version.eq.1.and.p1.eq.12.and.p2.eq.11) binary_orientation = .true.
        if(version.eq.2.and.p1.eq.6.and.p2.eq.7) sky_position = .true.
        if(version.eq.2.and.p1.eq.12.and.p2.eq.11) binary_orientation = .true.
        if(sky_position .and. plotsky.ge.1) project_map = .true.  !Make a special sky plot (i.e., plot stars or use projection) if plotsky>0 and RA,Dec are plotted
        
        
        if(npdf2d.ge.0) then
           plotthis = 0  !Determine to plot or save this combination of j1/j2 or p1/p2
           do i=1,npdf2d
              if(p1.eq.pdf2dpairs(i,1).and.p2.eq.pdf2dpairs(i,2)) plotthis = 1  !Use the data from the input file
           end do
           if(plotthis.eq.0) cycle
           if(prprogress.ge.1.and.update.eq.0) write(6,'(A,$)')trim(varnames(p1))//'-'//trim(varnames(p2))//' '
        else
           if(p2.le.p1) cycle
           write(6,*)upline !Move cursor up 1 line
           if(prprogress.ge.1.and.update.eq.0) write(6,'(F7.1,A)')real(countplots+1)/real(totplots)*100,'%    ('//trim(varnames(p1))//'-'//trim(varnames(p2))//')                                      '
        end if
        
        
        
        if(plot.eq.1) then
           if(file.eq.0) then
              npdf=npdf+1
              write(str,'(I3,A3)')200+npdf,'/xs'
              io = pgopen(trim(str))
              if(project_map) then
                 call pgpap(scrsz/0.5*scrrat,0.5)
              else
                 call pgpap(scrsz,scrrat)
              end if
              call pginitl(colour,file,whitebg)
           end if
           if(file.eq.1) then
              write(tempfile,'(A)') trim(outputname)//'__pdf2d__'//trim(varnames(p1))//'-'//trim(varnames(p2))//'.ppm'
              io = pgopen(trim(tempfile)//'/ppm')
              if(project_map) then
                 call pgpap(bmpsz/0.5*bmprat,0.5)
              else
                 call pgpap(bmpsz,bmprat)
              end if
              call pginitl(colour,file,whitebg)
           end if
           if(file.ge.2) then
              write(tempfile,'(A)') trim(outputname)//'__pdf2d__'//trim(varnames(p1))//'-'//trim(varnames(p2))//'.eps'
              io = pgopen(trim(tempfile)//trim(psclr))
              if(project_map) then
                 call pgpap(pssz/0.5*psrat,0.5)
              else
                 call pgpap(pssz,psrat)
              end if
              call pginitl(colour,file,whitebg)
              call pgscf(fonttype)
           end if
           if(io.le.0) then
              write(0,'(A,I4)')'   Error:  Cannot open PGPlot device.  Quitting the programme',io
              exitcode = 1
              return
           end if
           
           !call pgscr(3,0.,0.5,0.)
           call pgsch(sch)
        end if
        
        xmin = minval(alldat(ic,p1,1:n(ic)))
        xmax = maxval(alldat(ic,p1,1:n(ic)))
        ymin = minval(alldat(ic,p2,1:n(ic)))
        ymax = maxval(alldat(ic,p2,1:n(ic)))
        dx = xmax - xmin
        dy = ymax - ymin
        !write(6,'(A,2F10.5)')'  Xmin,Xmax: ',xmin,xmax
        !write(6,'(A,2F10.5)')'  Ymin,Ymax: ',ymin,ymax

        xx(1:n(ic)) = alldat(ic,p1,1:n(ic)) !Parameter 1
        yy(1:n(ic)) = alldat(ic,p2,1:n(ic)) !Parameter 2
        zz(1:n(ic)) = alldat(ic,1,1:n(ic))   !Likelihood
        
        if(.not.project_map) then
           xmin = xmin - 0.05*dx
           xmax = xmax + 0.05*dx
           ymin = ymin - 0.05*dy
           ymax = ymax + 0.05*dy
        end if
        
        
        
        
        
        !Prepare binning for a cute sky map in 2D PDF
        if(plot.eq.1 .and. project_map) then
           rat = 0.5 !scrrat !0.75
           !call pgpap(11.,rat)
           !call pgpap(scrsz,scrrat) !This causes a 'pgpage' when pggray is called...
           dx = xmax - xmin
           dy = ymax - ymin
           if(abs(dx)*15.lt.dy/rat) then !Expand x
              dx = dy/(15*rat)
              a = (xmin+xmax)*0.5
              xmin = a - 0.5*dx
              xmax = a + 0.5*dx
              if(prprogress.ge.1) write(6,'(A,F6.1,A3,F6.1,A,$)')'  Changing RA range to ',xmin,' - ',xmax,' h.'
           end if
           if(abs(dx)*15.gt.dy/rat) then !Expand y
              dy = abs(dx)*rat*15
              a = (ymin+ymax)*0.5
              ymin = a - 0.5*dy
              ymax = a + 0.5*dy
              if(prprogress.ge.1) write(6,'(A,F6.1,A3,F6.1,A,$)')'  Changing declination range to ',ymin,' - ',ymax,' deg.'
           end if
        end if !if(plot.eq.1 .and. project_map)
        
        !Force plotting and binning boundaries  CHECK: lose this? - yes: these are the binning ranges, not the plotting ranges; Don't necessarily want to bin the whole sky when plotting it.
        if(1.eq.2.and.wrapdata.eq.0.and.sky_position) then
           xmin = 0.
           xmax = 24.
           ymin = -90.
           ymax = 90.
        end if
        
        
        !Bin data and 'normalise' 2D PDF
        if(normpdf2d.le.2.or.normpdf2d.eq.4) then
           
           !Bin data:
           call bindata2d(n(ic),xx(1:n(ic)),yy(1:n(ic)),0,nbin2dx,nbin2dy,xmin,xmax,ymin,ymax,z,tr)  !Compute bin number rather than find it, ~10x faster
           
           if(normpdf2d.eq.1) z = max(0.,log10(z + 1.e-30))
           if(normpdf2d.eq.2) z = max(0.,sqrt(z + 1.e-30))
           if(normpdf2d.eq.4) then
              call identify_2d_ranges(p1,p2,nival,nbin2dx+1,nbin2dy+1,z,tr) !Get 2D probability ranges; identify to which range each bin belongs
              call calc_2d_areas(p1,p2,nival,nbin2dx+1,nbin2dy+1,z,tr,probarea) !Compute 2D probability areas; sum the areas of all bins
              trueranges2d(p1,p2) = truerange2d(z,nbin2dx+1,nbin2dy+1,startval(1,p1,1),startval(1,p2,1),tr)
              !write(6,'(/,A23,2(2x,A21))')'Probability interval:','Equivalent diameter:','Fraction of a sphere:'
              do i=1,nival
                 if(prival.ge.1.and.prprogress.ge.2 .and. (sky_position .or. binary_orientation)) then  !For sky position and orientation only
                    if(i.eq.1) write(6,*)
                    write(6,'(I10,F13.2,3(2x,F21.5))')i,ivals(i),probarea(i),sqrt(probarea(i)/pi)*2,probarea(i)*(pi/180.)**2/(4*pi)  !4pi*(180/pi)^2 = 41252.961 sq. degrees in a sphere
                 end if
                 probareas(p1,p2,i,1) = probarea(i)*(pi/180.)**2/(4*pi)  !Fraction of the sky
                 probareas(p1,p2,i,2) = sqrt(probarea(i)/pi)*2           !Equivalent diameter
                 probareas(p1,p2,i,3) = probarea(i)                      !Square degrees
              end do
           end if
        end if
        if(normpdf2d.eq.3) then
           call bindata2da(n(ic),xx(1:n(ic)),yy(1:n(ic)),zz(1:n(ic)),0,nbin2dx,nbin2dy,xmin,xmax,ymin,ymax,z,tr)  !Measure amount of likelihood in each bin
        end if
        
        
        
        
        !Swap RA boundaries for RA-Dec plot in 2D PDF
        if(sky_position) then
           a = xmin
           xmin = xmax
           xmax = a
           dx = -dx
        end if
        
        z = z/(maxval(z)+1.e-30)
        
        
        !Plot 2D PDF
        if(plot.eq.1) then
           
           !Set plot ranges for whole-sky map.  Does not affect binning
           if(project_map .and. plotsky.ge.2) then
              racentre = racentre*r2h
              xmin = racentre + 12.  !Must be the larger of the two
              xmax = racentre - 12.
              ymin = -90.
              ymax = 90.
           end if
           
           !Force plotting boundaries (not binning boundaries)
           if(1.eq.2.and.sky_position) then
              xmin = 24.
              xmax = 0.
              ymin = -90.
              ymax = 90.
           end if
           
           call pgsch(sch)
           if(project_map .and. plotsky.ge.2) then
              call pgsvp(0.08*sch,0.95,0.08*sch,1.0-0.05*sch)  !Make room for title and +90deg label
           else
              call pgsvp(0.08*sch,0.95,0.08*sch,1.0-0.033*sch)  !Make room for title.  Since sch is typically ~1.5*fontsize2d: 0.95 -> 1-0.05*fontsize ~ 1-0.03*sch
           end if
           
           call pgswin(xmin,xmax,ymin,ymax)
           if(project_map .and. (plotsky.eq.1.or.plotsky.eq.3).and.file.ge.2) then !Need dark background
              call pgsci(1)
              call pgrect(xmin,xmax,ymin,ymax)
           end if
           
           
           !Plot the actual 2D PDF (grey scales or colour)
           if(plpdf2d.eq.1.or.plpdf2d.eq.2) then
              
              !Set the colour schemes:
              if(normpdf2d.lt.4) then  !Grey scales
                 call pgscir(0,1e9)
                 call pgqcir(clr,maxclr)  !Maxclr is device-dependent
                 do i=0,maxclr-30  !Colour indices typically run 0-255, but this is device-dependent. Reserve ~0-29 for other purposes -> (maxclr-30) for these grey scales
                    x = real((maxclr-30) - i)/real(maxclr-30)          !White background
                    call pgscr(30+i,x,x,x)
                 end do
                 call pgscir(30,maxclr) !set colour-index range for pgimag
              end if
              
              if(normpdf2d.eq.4) then  !Colour
                 if(colour.eq.0) then
                    call pgscr(30,1.,1.,1.) !BG colour
                    if(nival.eq.2) then
                       call pgscr(31,0.5,0.5,0.5) !Grey
                       call pgscr(32,0.,0.,0.) !Black
                    end if
                    if(nival.eq.3) then
                       call pgscr(31,0.7,0.7,0.7) !
                       call pgscr(32,0.4,0.4,0.4) !
                       call pgscr(33,0.0,0.0,0.0) !
                    end if
                    if(nival.eq.4) then
                       call pgscr(31,0.75,0.75,0.75) !
                       call pgscr(32,0.50,0.50,0.50) !
                       call pgscr(33,0.25,0.25,0.25) !
                       call pgscr(34,0.00,0.00,0.00) !
                    end if
                    if(nival.eq.5) then
                       call pgscr(31,0.8,0.8,0.8) !
                       call pgscr(32,0.6,0.6,0.6) !
                       call pgscr(33,0.4,0.4,0.4) !
                       call pgscr(34,0.2,0.2,0.2) !
                       call pgscr(35,0.0,0.0,0.0) !
                    end if
                 end if
                 if(colour.ge.1) then
                    call pgscr(30,1.,1.,1.) !BG colour
                    if(nival.eq.2) then
                       call pgscr(31,1.,1.,0.) !Yellow
                       if(file.ge.2) call pgscr(31,0.8,0.7,0.) !Dark yellow
                       call pgscr(32,1.,0.,0.) !Red
                    end if
                    if(nival.eq.3) then
                       call pgscr(31,0.,0.,1.) !Blue
                       call pgscr(32,1.,1.,0.) !Yellow
                       if(file.ge.2) call pgscr(32,0.8,0.7,0.) !Dark yellow
                       call pgscr(33,1.,0.,0.) !Red
                    end if
                    if(nival.eq.4) then
                       call pgscr(31,0.,0.,1.) !Blue
                       call pgscr(32,0.,1.,0.) !Green
                       call pgscr(33,1.,1.,0.) !Yellow
                       if(file.ge.2) call pgscr(33,0.8,0.7,0.) !Dark yellow
                       call pgscr(34,1.,0.,0.) !Red
                    end if
                    if(nival.eq.5) then
                       call pgscr(31,0.,0.,1.) !Blue
                       call pgscr(32,0.,1.,0.) !Green
                       call pgscr(33,1.,1.,0.) !Yellow
                       if(file.ge.2) call pgscr(33,0.8,0.7,0.) !Dark yellow
                       call pgscr(34,1.,0.5,0.) !Orange
                       call pgscr(35,1.,0.,0.) !Red
                    end if
                 end if
                 call pgscir(30,30+nival) !set colour-index range for pgimag
              end if  !if(normpdf2d.eq.4)
              
              
              
              !Plot the PDF
              if(project_map .and. plotsky.ge.2) then
                 call pgimag_project(z,nbin2dx+1,nbin2dy+1,1,nbin2dx+1,1,nbin2dy+1,0.,1.,tr,map_projection)
              else
                 call pgimag(z,nbin2dx+1,nbin2dy+1,1,nbin2dx+1,1,nbin2dy+1,0.,1.,tr)
              end if
           end if
           
           
           !Plot stars in 2D PDF (over the grey scales, but underneath contours, lines, etc)
           if(project_map .and. (plotsky.eq.1.or.plotsky.eq.3)) then
              call pgswin(xmin*15,xmax*15,ymin,ymax) !Map works in degrees
              call plotthesky(xmin*15,xmax*15,ymin,ymax,rashift)
              call pgswin(xmin,xmax,ymin,ymax)
           end if
           call pgsci(1)
        end if !if(plot.eq.1)
        
        
        !Plot contours in 2D PDF
        if(plot.eq.1 .and. (plpdf2d.eq.1.or.plpdf2d.eq.3) .and. (.not.project_map .or. plotsky.eq.1.or.plotsky.eq.3)) then
           if(normpdf2d.lt.4) then
              ncont = 11
              do i=1,ncont
                 cont(i) = 0.01 + 2*real(i-1)/real(ncont-1)
                 if(project_map .and. (plotsky.eq.1.or.plotsky.eq.3)) cont(i) = 1.-cont(i)
              end do
              ncont = min(4,ncont) !Only use the first 4
           end if
           if(normpdf2d.eq.4) then
              ncont = nival
              do i=1,ncont
                 cont(i) = max(1. - real(i-1)/real(ncont-1),0.001)
                 !if(project_map) cont(i) = 1.-cont(i)
              end do
           end if
           
           call pgsls(1)
           if((.not.project_map .or. plotsky.ne.1.or.plotsky.ne.3) .and. normpdf2d.ne.4) then !First in bg colour
              call pgslw(2*lw)
              call pgsci(0)
              call pgcont(z,nbin2dx+1,nbin2dy+1,1,nbin2dx+1,1,nbin2dy+1,cont(1:ncont),ncont,tr)
           end if
           call pgslw(lw)
           call pgsci(1)
           if(project_map .and. (plotsky.eq.1.or.plotsky.eq.3)) call pgsci(7)
           call pgcont(z,nbin2dx+1,nbin2dy+1,1,nbin2dx+1,1,nbin2dy+1,cont(1:ncont),ncont,tr)
        end if
        
        
        !Save binned 2D PDF data
        if(savepdf.eq.1) then
           write(30,'(3I6,T100,A)')ic,p1,p2,'Chain number and variable number 1,2'
           write(30,'(2ES15.7,T100,A)')startval(ic,p1,1:2),'True and starting value p1'
           write(30,'(2ES15.7,T100,A)')startval(ic,p2,1:2),'True and starting value p2'
           write(30,'(6ES15.7,T100,A)')stats(ic,p1,1:6),'Stats: median, mean, absvar1, absvar2, stdev1, stdev2 for p1'
           write(30,'(6ES15.7,T100,A)')stats(ic,p2,1:6),'Stats: median, mean, absvar1, absvar2, stdev1, stdev2 for p2'
           write(30,'(5ES15.7,T100,A)')ranges(ic,c0,p1,1:5),'Ranges: lower,upper limit, centre, width, relative width for p1'
           write(30,'(5ES15.7,T100,A)')ranges(ic,c0,p2,1:5),'Ranges: lower,upper limit, centre, width, relative width for p2'
           write(30,'(4ES15.7,T100,A)')xmin,xmax,ymin,ymax,'Xmin,Xmax,Ymin,Ymax of PDF'
           write(30,'(6ES15.7,T100,A)')tr,'Tr'              
           do i=1,nbin2dx+1
              do j=1,nbin2dy+1
                 write(30,'(ES15.7,$)')z(i,j)
              end do
              write(30,'(1x)')
           end do
        end if
        
        
        
        !Plot true value, median, ranges, etc. in 2D PDF
        if(plot.eq.1) then
           if(.not.project_map.or.plotsky.eq.1) then
              call pgsci(1)
              
              !Plot max likelihood in 2D PDF
              if(pllmax.ge.1) then
                 call pgsci(1); call pgsls(5)
                 
                 plx = pldat(icloglmax,p1,iloglmax)
                 if(version.eq.1.and.p1.eq.8 .or. version.eq.2.and.p1.eq.6) plx = rev24(plx)
                 if(version.eq.1.and.(p1.eq.10.or.p1.eq.13) .or. version.eq.2.and.(p1.eq.9.or.p1.eq.13.or.p1.eq.16)) plx = rev360(plx)
                 if(version.eq.1.and.p1.eq.12 .or. version.eq.2.and.p1.eq.8) plx = rev180(plx)
                 call pgline(2,(/plx,plx/),(/-1.e20,1.e20/)) !Max logL
                 if(version.eq.1.and.p1.eq.8 .or. version.eq.2.and.p1.eq.6) then
                    call pgline(2,(/plx-24.,plx-24./),(/-1.e20,1.e20/)) !Max logL
                    call pgline(2,(/plx+24.,plx+24./),(/-1.e20,1.e20/)) !Max logL
                 end if
                 if(version.eq.1.and.(p1.eq.10.or.p1.eq.13) .or. version.eq.2.and.(p1.eq.9.or.p1.eq.13.or.p1.eq.16)) then
                    call pgline(2,(/plx-360.,plx-360./),(/-1.e20,1.e20/)) !Max logL
                    call pgline(2,(/plx+360.,plx+360./),(/-1.e20,1.e20/)) !Max logL
                 end if
                 if(version.eq.1.and.p1.eq.12 .or. version.eq.2.and.p1.eq.8) then
                    call pgline(2,(/plx-180.,plx-180./),(/-1.e20,1.e20/)) !Max logL
                    call pgline(2,(/plx+180.,plx+180./),(/-1.e20,1.e20/)) !Max logL
                 end if
                 
                 ply = pldat(icloglmax,p2,iloglmax)
                 if(version.eq.1.and.p2.eq.8 .or. version.eq.2.and.p2.eq.6) ply = rev24(ply)
                 if(version.eq.1.and.(p2.eq.10.or.p2.eq.13) .or. version.eq.2.and.(p2.eq.9.or.p2.eq.13.or.p2.eq.16)) ply = rev360(ply)
                 if(version.eq.1.and.p2.eq.12 .or. version.eq.2.and.p2.eq.8) ply = rev180(ply)
                 call pgline(2,(/-1.e20,1.e20/),(/ply,ply/)) !Max logL
                 if(version.eq.1.and.p2.eq.8 .or. version.eq.2.and.p2.eq.6) then
                    call pgline(2,(/-1.e20,1.e20/),(/ply-24.,ply-24./)) !Max logL
                    call pgline(2,(/-1.e20,1.e20/),(/ply+24.,ply+24./)) !Max logL
                 end if
                 if(version.eq.1.and.(p2.eq.10.or.p2.eq.13) .or. version.eq.2.and.(p2.eq.9.or.p2.eq.13.or.p2.eq.16)) then
                    call pgline(2,(/-1.e20,1.e20/),(/ply-360.,ply-360./)) !Max logL
                    call pgline(2,(/-1.e20,1.e20/),(/ply+360.,ply+360./)) !Max logL
                 end if
                 if(version.eq.1.and.p2.eq.12 .or. version.eq.2.and.p2.eq.8) then
                    call pgline(2,(/-1.e20,1.e20/),(/ply-180.,ply-180./)) !Max logL
                    call pgline(2,(/-1.e20,1.e20/),(/ply+180.,ply+180./)) !Max logL
                 end if
                 
                 call pgpoint(1,plx,ply,12)
              end if
              
              
              if(project_map .and. (plotsky.eq.1.or.plotsky.eq.3)) call pgsci(0)
              call pgsls(2)
              
              !Plot true value in 2D PDF
              if((pltrue.eq.1.or.pltrue.eq.3).and.(.not.project_map) .or. &
                   !((pltrue.eq.2.or.pltrue.eq.4) .and. (p1.eq.2.and.p2.eq.3).or.(p1.eq.6.and.p2.eq.7).or.(p1.eq.14.and.p2.eq.15)) ) then
                   ((pltrue.eq.2.or.pltrue.eq.4) .and. (p1.eq.2.and.p2.eq.3).or.(p1.eq.14.and.p2.eq.15)) ) then
                 !call pgline(2,(/startval(ic,p1,1),startval(ic,p1,1)/),(/-1.e20,1.e20/))
                 !call pgline(2,(/-1.e20,1.e20/),(/startval(ic,p2,1),startval(ic,p2,1)/))
                 
                 if(mergechains.ne.1.or.ic.le.1) then !The units of the true values haven't changed (e.g. from rad to deg) for ic>1 (but they have for the starting values, why?)
                    !x
                    call pgsls(2); call pgsci(1)
                    if(pllmax.eq.0) call pgsls(3)  !Dash-dotted line for true value when Lmax line isn't plotted (should we do this always?)
                    plx = startval(ic,p1,1)
                    if(version.eq.1.and.p1.eq.8 .or. version.eq.2.and.p1.eq.6) plx = rev24(plx)
                    if(version.eq.1.and.(p1.eq.10.or.p1.eq.13) .or. version.eq.2.and.(p1.eq.9.or.p1.eq.13.or.p1.eq.16)) plx = rev360(plx)
                    if(version.eq.1.and.p1.eq.12 .or. version.eq.2.and.p1.eq.8) plx = rev180(plx)
                    call pgline(2,(/plx,plx/),(/-1.e20,1.e20/)) !True value
                    if(version.eq.1.and.p1.eq.8 .or. version.eq.2.and.p1.eq.6) then
                       call pgline(2,(/plx-24.,plx-24./),(/-1.e20,1.e20/)) !True value
                       call pgline(2,(/plx+24.,plx+24./),(/-1.e20,1.e20/)) !True value
                    end if
                    if(version.eq.1.and.(p1.eq.10.or.p1.eq.13) .or. version.eq.2.and.(p1.eq.9.or.p1.eq.13.or.p1.eq.16)) then
                       call pgline(2,(/plx-360.,plx-360./),(/-1.e20,1.e20/)) !True value
                       call pgline(2,(/plx+360.,plx+360./),(/-1.e20,1.e20/)) !True value
                    end if
                    if(version.eq.1.and.p1.eq.12 .or. version.eq.2.and.p1.eq.8) then
                       call pgline(2,(/plx-180.,plx-180./),(/-1.e20,1.e20/)) !True value
                       call pgline(2,(/plx+180.,plx+180./),(/-1.e20,1.e20/)) !True value
                    end if
                    
                    !y
                    call pgsls(2); call pgsci(1)
                    if(pllmax.eq.0) call pgsls(3)  !Dash-dotted line for true value when Lmax line isn't plotted (should we do this always?)
                    ply = startval(ic,p2,1)
                    if(version.eq.1.and.p2.eq.8 .or. version.eq.2.and.p2.eq.6) ply = rev24(ply)
                    if(version.eq.1.and.(p2.eq.10.or.p2.eq.13) .or. version.eq.2.and.(p2.eq.9.or.p2.eq.13.or.p2.eq.16)) ply = rev360(ply)
                    if(version.eq.1.and.p2.eq.12 .or. version.eq.2.and.p2.eq.8) ply = rev180(ply)
                    call pgline(2,(/-1.e20,1.e20/),(/ply,ply/)) !True value
                    if(version.eq.1.and.p2.eq.8 .or. version.eq.2.and.p2.eq.6) then
                       call pgline(2,(/-1.e20,1.e20/),(/ply-24.,ply-24./)) !True value
                       call pgline(2,(/-1.e20,1.e20/),(/ply+24.,ply+24./)) !True value
                    end if
                    if(version.eq.1.and.(p2.eq.10.or.p2.eq.13) .or. version.eq.2.and.(p2.eq.9.or.p2.eq.13.or.p2.eq.16)) then
                       call pgline(2,(/-1.e20,1.e20/),(/ply-360.,ply-360./)) !True value
                       call pgline(2,(/-1.e20,1.e20/),(/ply+360.,ply+360./)) !True value
                    end if
                    if(version.eq.1.and.p2.eq.12 .or. version.eq.2.and.p2.eq.8) then
                       call pgline(2,(/-1.e20,1.e20/),(/ply-180.,ply-180./)) !True value
                       call pgline(2,(/-1.e20,1.e20/),(/ply+180.,ply+180./)) !True value
                    end if
                    
                    call pgpoint(1,plx,ply,18)
                 end if
              end if
              call pgsci(1)
              call pgsls(4)
              
              
              !Plot starting values in 2D PDF
              !call pgline(2,(/startval(ic,p1,2),startval(ic,p1,2)/),(/-1.e20,1.e20/))
              !call pgline(2,(/-1.e20,1.e20/),(/startval(ic,p2,2),startval(ic,p2,2)/))
              
              call pgsci(2)
              
              !Plot interval ranges in 2D PDF
              if(plrange.eq.2.or.plrange.eq.3.or.plrange.eq.5.or.plrange.eq.6) then
                 call pgsls(1)
                 call pgsch(sch*0.6)
                 call pgsah(1,45.,0.1)
                 a = 0.0166667*sch
                 call pgarro(ranges(ic,c0,p1,3),ymin+dy*a,ranges(ic,c0,p1,1),ymin+dy*a)
                 call pgarro(ranges(ic,c0,p1,3),ymin+dy*a,ranges(ic,c0,p1,2),ymin+dy*a)
                 a = 0.0333333*sch
                 call pgptxt(ranges(ic,c0,p1,3),ymin+dy*a,0.,0.5,'\(2030)\d90%\u')
                 a = 0.0233333*sch
                 call pgarro(xmin+dx*a,ranges(ic,c0,p2,3),xmin+dx*a,ranges(ic,c0,p2,1))
                 call pgarro(xmin+dx*a,ranges(ic,c0,p2,3),xmin+dx*a,ranges(ic,c0,p2,2))
                 a = 0.01*sch
                 call pgptxt(xmin+dx*a,ranges(ic,c0,p2,3),90.,0.5,'\(2030)\d90%\u')
              end if
              
              call pgsch(sch)
              call pgsls(2)
              
              
              !Plot medians in 2D PDF
              if(plmedian.eq.2.or.plmedian.eq.3.or.plmedian.eq.5.or.plmedian.eq.6) then
                 call pgline(2,(/stats(ic,p1,1),stats(ic,p1,1)/),(/-1.e20,1.e20/))
                 call pgline(2,(/-1.e20,1.e20/),(/stats(ic,p2,1),stats(ic,p2,1)/))
                 call pgpoint(1,stats(ic,p1,1),stats(ic,p2,1),18)
              end if
              
              call pgsls(1)
              
              
           end if  !if(.not.project_map.or.plotsky.eq.1)
           
           !Plot big symbol at true position in sky map
           if(project_map .and. (pltrue.eq.1.or.pltrue.eq.3)) then
              call pgsch(sch*1.5) !Use 1.5 for plsym=8, 2 for plsym=18
              call pgslw(lw*2)
              call pgsci(9)
              if(normpdf2d.eq.4) call pgsci(1)  !Black
              plx = startval(ic,p1,1)
              ply = startval(ic,p2,1)
              if(plotsky.eq.2.or.plotsky.eq.4) call project_skymap(plx,ply,racentre,map_projection)
              call pgpoint(1,plx,ply,8)
              call pgsch(sch)
              call pgslw(lw)
              call pgsci(1)
           end if
           
        end if  !if(plot.eq.1)
        
        
        
        
        
        
        !Print axes, axis labels and plot title
        if(plot.eq.1) then
           
           !Plot coordinate axes and axis labels in 2D PDF
           call pgsls(1)
           call pgslw(lw2)
           call pgsci(1)
           if(.not.project_map) call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)
           !if(plotsky.eq.1) then
           !   call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0) !Box, ticks, etc in white
           !   call pgsci(1)
           !   call pgbox('N',0.0,0,'N',0.0,0) !Number labels in black
           !end if
           call pgmtxt('B',2.2,0.5,0.5,trim(pgvarns(p1)))
           call pgmtxt('L',1.7,0.5,0.5,trim(pgvarns(p2)))
           
           
           !Print 2D probability ranges in plot title
           if(prival.ge.1.and.normpdf2d.eq.4.and. (sky_position .or. binary_orientation)) then  !For sky position and orientation only
              string = ' '
              do c = 1,nival
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
                 
                 i = 3  !2-use degrees, 3-square degrees
                 a = probareas(p1,p2,c,i)
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
                    call pgsch(sch*0.85) !Needed to fit the square-degree sign in
                    if(quality.eq.3) call pgsch(sch*0.6) !Poster
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
                 if(quality.eq.91) then !NINJA
                    call pgsch(sch)
                    if(fonttype.eq.2) then
                       write(string,'(I2,A7,A10)')c,'\(2144)',''
                    else
                       write(string,'(I2,A7,A10)')c,'\(0644)',''
                    end if
                 end if
                 a = (real(c-1)/real(nival-1) - 0.5)*0.7 + 0.5
                 call pgsci(30+nival+1-c)
                 if(project_map .and. plotsky.ge.2) then
                    call pgmtxt('T',1.0,a,0.5,trim(string))  !Print title
                 else
                    call pgmtxt('T',0.5,a,0.5,trim(string))  !Print title
                 end if
                 call pgsch(sch)
              end do
              call pgsci(1)
           end if
        end if  !if(plot.eq.1)
        
        
        
        if(plot.eq.1) then
           countplots = countplots + 1  !The current plot is number countplots
           
           !Convert plot
           if(file.eq.1) then
              call pgend
              if(countplots.eq.npdf2d) then !Convert the last plot in the foreground, so that the process finishes before deleting the original file
                 i = system('convert -resize '//trim(bmpxpix)//' -depth 8 -unsharp '//trim(unsharppdf2d)//' '//trim(tempfile)//' '// &
                      trim(outputdir)//'/'//trim(outputname)//'__pdf2d__'//trim(varnames(p1))//'-'//trim(varnames(p2))//'.png')
              else !in the background
                 i = system('convert -resize '//trim(bmpxpix)//' -depth 8 -unsharp '//trim(unsharppdf2d)//' '//trim(tempfile)//' '// &
                      trim(outputdir)//'/'//trim(outputname)//'__pdf2d__'//trim(varnames(p1))//'-'//trim(varnames(p2))//'.png &')
              end if
              if(i.ne.0) write(0,'(A,I6)')'  Error converting plot',i
           end if
           !if(file.ge.2.and.multipagefile) call pgpage
           if(file.ge.2) then
              call pgend
              if(file.eq.3) then
                 i = system('eps2pdf '//trim(tempfile))
                 if(i.ne.0) write(0,'(A,I6)')'  Error converting plot',i
              end if
           end if
        end if !if(plot.eq.1)
        
     end do !p2
  end do !p1
  
  
  if(savepdf.eq.1) close(30)
  
  if(plot.eq.1) then
     if(file.ne.1) call pgend
     !if(file.ge.2.and.multipagefile) then
     !   if(abs(j2-j1).le.1) then
     !      if(file.eq.3) i = system('eps2pdf pdf2d.eps  -o '//trim(outputdir)//'/'//trim(outputname)//'__pdf2d_'//trim(varnames(j1))//'-'//trim(varnames(j2))//'.pdf  &> /dev/null')
     !      i = system('mv -f pdf2d.eps '//trim(outputdir)//'/'//trim(outputname)//'__pdf2d_'//trim(varnames(j1))//'-'//trim(varnames(j2))//'.eps')
     !   else
     !      if(file.eq.3) i = system('eps2pdf pdf2d.eps  -o '//trim(outputdir)//'/'//trim(outputname)//'__pdf2d.pdf  &> /dev/null')
     !      i = system('mv -f pdf2d.eps '//trim(outputdir)//'/'//trim(outputname)//'__pdf2d.eps')
     !   end if
     !end if
  
     !Remove all the .ppm files
     if(file.eq.1) then
        do p1=j1,j2
           do p2=j1,j2
              if(npdf2d.ge.0) then
                 plotthis = 0  !Determine to plot or save this combination of j1/j2 or p1/p2
                 do i=1,npdf2d
                    if(p1.eq.pdf2dpairs(i,1).and.p2.eq.pdf2dpairs(i,2)) plotthis = 1  !Use the data from the input file
                 end do
                 if(plotthis.eq.0) cycle
              else
                 if(p2.le.p1) cycle
              end if
              write(tempfile,'(A)') trim(outputname)//'__pdf2d__'//trim(varnames(p1))//'-'//trim(varnames(p2))//'.ppm'
              i = system('rm -f '//trim(tempfile))
           end do
        end do
     end if
     
  end if !plot.eq.1
  
  
  
end subroutine pdfs2d
!************************************************************************************************************************************






!************************************************************************************************************************************
subroutine bindata2dold(n,x,y,norm,nxbin,nybin,xmin1,xmax1,ymin1,ymax1,z,tr)  !Count the number of points in each bin
  !x - input: data, n points
  !norm - input: normalise (1) or not (0)
  !nbin - input: number of bins
  !xmin, xmax - in/output: set xmin=xmax to auto-determine
  !xbin, ybin - output: binned data (x, y).  The x values are the left side of the bin!
  
  implicit none
  integer :: i,n,bx,by,nxbin,nybin,norm
  real :: x(n),y(n),xbin(nxbin+1),ybin(nybin+1),z(nxbin+1,nybin+1)
  real :: xmin,xmax,ymin,ymax,dx,dy,xmin1,xmax1,ymin1,ymax1,tr(6)
  
  !write(6,'(A4,5I8)')'n:',norm,nxbin,nybin
  !write(6,'(A4,2F8.3)')'x:',xmin1,xmax1
  !write(6,'(A4,2F8.3)')'y:',ymin1,ymax1
  
  xmin = xmin1
  xmax = xmax1
  ymin = ymin1
  ymax = ymax1
  
  if(abs(xmin-xmax)/(xmax+1.e-30).lt.1.e-20) then !Autodetermine
     xmin = minval(x(1:n))
     xmax = maxval(x(1:n))
  end if
  dx = abs(xmax - xmin)/real(nxbin)
  if(abs(ymin-ymax)/(ymax+1.e-30).lt.1.e-20) then !Autodetermine
     ymin = minval(y(1:n))
     ymax = maxval(y(1:n))
  end if
  dy = abs(ymax - ymin)/real(nybin)
  do bx=1,nxbin+1
     !xbin(bx) = xmin + (real(bx)-0.5)*dx  !x is the centre of the bin
     xbin(bx) = xmin + (bx-1)*dx          !x is the left of the bin
  end do
  do by=1,nybin+1
     !ybin(by) = ymin + (real(by)-0.5)*dy  !y is the centre of the bin
     ybin(by) = ymin + (by-1)*dy          !y is the left of the bin
  end do
  
  !write(6,'(50F5.2)'),x(1:50)
  !write(6,'(50F5.2)'),y(1:50)
  !write(6,'(20F8.5)'),xbin
  !write(6,'(20F8.5)'),ybin
  
  z = 0.
  !ztot = 0.
  do i=1,n
     bxl: do bx=1,nxbin
        do by=1,nybin
           if(x(i).ge.xbin(bx)) then
              if(x(i).lt.xbin(bx+1)) then
                 if(y(i).ge.ybin(by)) then
                    if(y(i).lt.ybin(by+1)) then
                       z(bx,by) = z(bx,by) + 1.
                       exit bxl !exit bx loop; if point i fits this bin, don't try other bins. Speeds things up ~2x
                    end if
                 end if
              end if
           end if
           
        end do !by
     end do bxl !bx
  end do !i
  !if(norm.eq.1) z = z/(ztot+1.e-30)
  if(norm.eq.1) z = z/maxval(z+1.e-30)
  
  if(abs(xmin1-xmax1)/(xmax1+1.e-30).lt.1.e-20) then
     xmin1 = xmin
     xmax1 = xmax
  end if
  if(abs(ymin1-ymax1)/(ymax1+1.e-30).lt.1.e-20) then
     ymin1 = ymin
     ymax1 = ymax
  end if
  
  !Determine transformation elements for pgplot (pggray, pgcont, pgimag)
  tr(1) = xmin - dx/2.
  tr(2) = dx
  tr(3) = 0.
  tr(4) = ymin - dy/2.
  tr(5) = 0.
  tr(6) = dy
  
end subroutine bindata2dold
!************************************************************************************************************************************


!************************************************************************************************************************************
subroutine bindata2d(n,x,y,norm,nxbin,nybin,xmin1,xmax1,ymin1,ymax1,z,tr)  !Compute bin number rather than search for it ~10x faster
  !x - input: data, n points
  !norm - input: normalise (1) or not (0)
  !nbin - input: number of bins
  !xmin, xmax - in/output: set xmin=xmax to auto-determine
  
  implicit none
  integer :: i,n,bx,by,nxbin,nybin,norm
  real :: x(n),y(n),z(nxbin+1,nybin+1)
  real :: xmin,xmax,ymin,ymax,dx,dy,xmin1,xmax1,ymin1,ymax1,tr(6)
  
  xmin = xmin1
  xmax = xmax1
  ymin = ymin1
  ymax = ymax1
  
  if(abs(xmin-xmax)/(xmax+1.e-30).lt.1.e-20) then !Autodetermine
     xmin = minval(x(1:n))
     xmax = maxval(x(1:n))
  end if
  dx = abs(xmax - xmin)/real(nxbin)
  if(abs(ymin-ymax)/(ymax+1.e-30).lt.1.e-20) then !Autodetermine
     ymin = minval(y(1:n))
     ymax = maxval(y(1:n))
  end if
  dy = abs(ymax - ymin)/real(nybin)
  
  
  
  !Determine transformation elements for pgplot (pggray, pgcont, pgimag)
  tr(1) = xmin - dx/2.
  tr(2) = dx
  tr(3) = 0.
  tr(4) = ymin - dy/2.
  tr(5) = 0.
  tr(6) = dy
  
  z = 0.
  do i=1,n
     bx = floor((x(i) - xmin)/dx) + 1 
     by = floor((y(i) - ymin)/dy) + 1
     if(bx.lt.1.or.bx.gt.nxbin.or.by.lt.1.or.by.gt.nybin) then
        if(bx.lt.0.or.bx.gt.nxbin+1.or.by.lt.0.or.by.gt.nybin+1) then  !Treat an error of 1 bin as round-off
           bx = max(min(bx,nxbin),1)
           by = max(min(by,nybin),1)
           z(bx,by) = z(bx,by) + 1.
        else
           if(bx.lt.0.or.bx.gt.nxbin+1) write(0,'(A,I7,A2,F8.3,A,I4,A,I4,A1)')'  Bindata2d:  error for X data point',i,' (',x(i),').  I found bin',bx,', but it should lie between 1 and',nxbin,'.'
           if(by.lt.0.or.by.gt.nybin+1) write(0,'(A,I7,A2,F8.3,A,I4,A,I4,A1)')'  Bindata2d:  error for Y data point',i,' (',y(i),').  I found bin',by,', but it should lie between 1 and',nybin,'.'
        end if
     else
        z(bx,by) = z(bx,by) + 1.
     end if
  end do
  
  !if(norm.eq.1) z = z/(ztot+1.e-30)
  if(norm.eq.1) z = z/maxval(z+1.e-30)
  
  if(abs(xmin1-xmax1)/(xmax1+1.e-30).lt.1.e-20) then
     xmin1 = xmin
     xmax1 = xmax
  end if
  if(abs(ymin1-ymax1)/(ymax1+1.e-30).lt.1.e-20) then
     ymin1 = ymin
     ymax1 = ymax
  end if
  
end subroutine bindata2d
!************************************************************************************************************************************


!************************************************************************************************************************************
subroutine bindata2da(n,x,y,z,norm,nxbin,nybin,xmin1,xmax1,ymin1,ymax1,zz,tr)  !Measure the amount of likelihood in each bin
  !x,y - input: data, n points
  !z - input: amount for each point (x,y)
  !norm - input: normalise (1) or not (0)
  !nxbin,nybin - input: number of bins in each dimension
  !xmin1,xmax1 - in/output: ranges in x dimension, set xmin=xmax as input to auto-determine
  !ymin1,ymax1 - in/output: ranges in y dimension, set ymin=ymax as input to auto-determine
  !zz - output: binned data zz(x,y).  The x,y values are the left side of the bin(?)
  !tr - output: transformation elements for pgplot (pggray, pgcont)
  
  implicit none
  integer :: i,n,bx,by,nxbin,nybin,norm
  real :: x(n),y(n),z(n),xbin(nxbin+1),ybin(nybin+1),zz(nxbin+1,nybin+1),zztot,xmin,xmax,ymin,ymax,dx,dy,xmin1,xmax1,ymin1,ymax1
  real :: tr(6),zmin
  
  !write(6,'(A4,5I8)')'n:',norm,nxbin,nybin
  !write(6,'(A4,2F8.3)')'x:',xmin1,xmax1
  !write(6,'(A4,2F8.3)')'y:',ymin1,ymax1
  
  xmin = xmin1
  xmax = xmax1
  ymin = ymin1
  ymax = ymax1
  zmin = minval(z)
  
  if(abs(xmin-xmax)/(xmax+1.e-30).lt.1.e-20) then !Autodetermine
     xmin = minval(x(1:n))
     xmax = maxval(x(1:n))
  end if
  dx = abs(xmax - xmin)/real(nxbin)
  if(abs(ymin-ymax)/(ymax+1.e-30).lt.1.e-20) then !Autodetermine
     ymin = minval(y(1:n))
     ymax = maxval(y(1:n))
  end if
  dy = abs(ymax - ymin)/real(nybin)
  do bx=1,nxbin+1
     !xbin(bx) = xmin + (real(bx)-0.5)*dx  !x is the centre of the bin
     xbin(bx) = xmin + (bx-1)*dx          !x is the left of the bin
  end do
  do by=1,nybin+1
     !ybin(by) = ymin + (real(by)-0.5)*dy  !y is the centre of the bin
     ybin(by) = ymin + (by-1)*dy          !y is the left of the bin
  end do
  
  !write(6,'(50F5.2)'),x(1:50)
  !write(6,'(50F5.2)'),y(1:50)
  !write(6,'(20F8.5)'),xbin
  !write(6,'(20F8.5)'),ybin
  
  zz = 0.
  zztot = 0.
  do bx=1,nxbin
     do by=1,nybin
        zz(bx,by) = 0.
        do i=1,n
           !if(x(i).ge.xbin(bx).and.x(i).lt.xbin(bx+1) .and. y(i).ge.ybin(by).and.y(i).lt.ybin(by+1)) zz(bx,by) = zz(bx,by) + 1.
           if(x(i).ge.xbin(bx).and.x(i).lt.xbin(bx+1) .and. y(i).ge.ybin(by).and.y(i).lt.ybin(by+1)) zz(bx,by) = zz(bx,by) + exp(z(i) - zmin)
           !write(6,'(2I4,8F10.5)')bx,by,x(i),xbin(bx),xbin(bx+1),y(i),ybin(by),ybin(by+1),zz(bx,by),z(i)
        end do
        zztot = zztot + zz(bx,by) 
        !write(6,'(2I4,5x,4F6.3,5x,10I8)')bx,by,xbin(bx),xbin(bx+1),ybin(by),ybin(by+1),nint(zz(bx,by))
     end do
     !write(6,'(I4,5x,2F6.3,5x,10I8)')bx,xbin(bx),xbin(bx+1),nint(zz(bx,1:nybin))
     end do
  !if(norm.eq.1) z = z/(zztot+1.e-30)
  if(norm.eq.1) z = z/maxval(z+1.e-30)
  
  if(abs(xmin1-xmax1)/(xmax1+1.e-30).lt.1.e-20) then
     xmin1 = xmin
     xmax1 = xmax
  end if
  if(abs(ymin1-ymax1)/(ymax1+1.e-30).lt.1.e-20) then
     ymin1 = ymin
     ymax1 = ymax
  end if
  
  !Determine transformation elements for pgplot (pggray, pgcont)
  tr(1) = xmin - dx/2.
  tr(2) = dx
  tr(3) = 0.
  tr(4) = ymin - dy/2.
  tr(5) = 0.
  tr(6) = dy
  
end subroutine bindata2da
!************************************************************************************************************************************




!************************************************************************
subroutine identify_2d_ranges(p1,p2,ni,nx,ny,z,tr)
  !Get the 2d probability intervals; z lies between 1 (in 100% range) and ni (in lowest-% range, e.g. 90%)
  use constants
  use analysemcmc_settings
  implicit none
  integer :: p1,p2,ni,nx,ny,nn,indx(nx*ny),i,b,ib,full(ni),iy
  real :: z(nx,ny),x1(nx*ny),x2(nx*ny),tot,np,tr(6),y
  
  
  !Weight number of points in each bin by bin size for position/orientation plots
  do iy = 1,ny
     if(changevar.ge.1) then
        if(version.eq.1 .and. (p1.eq.8.and.p2.eq.9 .or. p1.eq.12.and.p2.eq.11)) then  !Then: RA-Dec or (phi/theta_Jo)/(psi/i) plot, convert lon -> lon * 15 * cos(lat)
           y = tr(4) + tr(6)*iy
           if(p1.eq.8) then
              if(abs(y).le.90.) then
                 z(1:nx,iy) = z(1:nx,iy)/(cos(y*rd2r)+1.e-30)
              else  !This can happen when the PDF lies close to the pole
                 z(1:nx,iy) = 0.
              end if
           else if(p1.eq.12) then
              if(y.ge.0..and.y.lt.180.) then
                 z(1:nx,iy) = z(1:nx,iy)/(abs(sin(y*rd2r))+1.e-30)
              else  !This can happen when the PDF lies close to the pole
                 z(1:nx,iy) = 0.
                 !write(0,'(//,A,//)')'  *** identify_2d_ranges:  sin(y)<0.  Please check whether the if(y.ge.0..and.y.lt.180.) statement works properly ***'
              end if
           end if
        end if
        if(version.eq.2 .and. (p1.eq.6.and.p2.eq.7 .or. p1.eq.10.and.p2.eq.8)) then  !Then: RA-Dec or (phi/theta_Jo)/(psi/i) plot, convert lon -> lon * 15 * cos(lat)
           y = tr(4) + tr(6)*iy
           if(p1.eq.6) then
              if(abs(y).le.90.) then
                 z(1:nx,iy) = z(1:nx,iy)/(cos(y*rd2r)+1.e-30)
              else  !This can happen when the PDF lies close to the pole
                 z(1:nx,iy) = 0.
              end if
           else if(p1.eq.10) then
              if(y.ge.0..and.y.lt.180.) then
                 z(1:nx,iy) = z(1:nx,iy)/(abs(sin(y*rd2r))+1.e-30)
              else  !This can happen when the PDF lies close to the pole
                 z(1:nx,iy) = 0.
                 !write(0,'(//,A,//)')'  *** identify_2d_ranges:  sin(y)<0.  Please check whether the if(y.ge.0..and.y.lt.180.) statement works properly ***'
              end if
           end if
        end if
     end if !if(changevar.ge.1)
  end do !iy
  
  
  nn = nx*ny
  x1 = reshape(z,(/nn/))  !x1 is an 1D array with the same data as the 2D array z
  call rindexx(nn,-x1(1:nn),indx(1:nn)) ! -x1: sort the 1D array to descending value
  
  np = sum(z)
  tot = 0.
  full = 0
  do b=1,nn !Loop over bins in 1D array
     ib = indx(b)
     x2(ib) = 0.
     if(x1(ib).eq.0.) cycle
     tot = tot + x1(ib)
     do i=ni,1,-1 !Loop over intervals
        if(tot.le.np*ivals(i)) then
           x2(ib) = real(ni-i+1)  !e.g. x2(b) = ni if within 68%, ni-1 if within 95%, etc, and 1 if within 99.7%
        else
           if(prprogress.ge.3.and.full(i).eq.0) then !Report the number of points in the lastly selected bin
              if(i.eq.1) write(6,'(A,$)')'Last bin:'
              !write(6,'(F6.3,I5,$)')ivals(i),nint(x1(ib))
              write(6,'(I5,$)')nint(x1(ib))
              full(i) = 1
           end if
        end if
        !write(6,'(2I4, F6.2, 3F20.5)')b,i, ivals(i), np,tot,np*ivals(i)
     end do
  end do
  
  z = reshape(x2, (/nx,ny/))  ! z lies between 1 and ni
end subroutine identify_2d_ranges
!************************************************************************



!************************************************************************
!Compute 2D probability areas
subroutine calc_2d_areas(p1,p2,ni,nx,ny,z,tr,area)
  use constants
  use analysemcmc_settings
  implicit none
  integer :: p1,p2,ni,nx,ny,ix,iy,i,i1,iv
  real :: z(nx,ny),tr(6),y,dx,dy,area(ni)
  
  area = 0.
  
  do ix = 1,nx
     do iy = 1,ny
        dx = tr(2)
        dy = tr(6)
        if(changevar.ge.1) then
           if(version.eq.1 .and. (p1.eq.8.and.p2.eq.9 .or. p1.eq.12.and.p2.eq.11)) then  !Then: RA-Dec or (phi/theta_Jo)/(psi/i) plot, convert lon -> lon * 15 * cos(lat)
              y = tr(4) + tr(6)*iy
              if(p1.eq.8) then
                 dx = dx*cos(y*rd2r)
              else if(p1.eq.12) then
                 dx = dx*abs(sin(y*rd2r))  !Necessary for i-psi plot?
              end if
              if(p1.eq.8) dx = dx*15
           end if
           if(version.eq.2 .and. (p1.eq.6.and.p2.eq.7 .or. p1.eq.10.and.p2.eq.8)) then  !Then: RA-Dec or (phi/theta_Jo)/(psi/i) plot, convert lon -> lon * 15 * cos(lat)
              y = tr(4) + tr(6)*iy
              if(p1.eq.6) then
                 dx = dx*cos(y*rd2r)
              else if(p1.eq.10) then
                 dx = dx*abs(sin(y*rd2r))  !Necessary for i-psi plot?
              end if
              if(p1.eq.6) dx = dx*15
           end if
        end if
        iv = nint(z(ix,iy))
        do i=1,ni
           if(iv.ge.i) then
              i1 = ni-i+1
              area(i1) = area(i1) + dx*dy
           end if
        end do !i
        
     end do !iy
  end do !ix
end subroutine calc_2d_areas
!************************************************************************


!************************************************************************
function truerange2d(z,nx,ny,truex,truey,tr)
  !Get the smallest probability area in which the true values lie
  implicit none
  integer :: nx,ny,ix,iy,truerange2d
  real :: truex,truey,z(nx,ny),tr(6)
  
  !x = tr(1) + tr(2)*ix + tr(3)*iy
  !y = tr(4) + tr(5)*ix + tr(6)*iy
  ix = floor((truex - tr(1))/tr(2))
  iy = floor((truey - tr(4))/tr(6))
  if(ix.lt.1.or.ix.gt.nx.or.iy.lt.1.or.iy.gt.ny) then
     truerange2d = 0
  else
     truerange2d = nint(z(ix,iy))
  end if
end function truerange2d
!************************************************************************








!************************************************************************************************************************************
subroutine plotthesky(bx1,bx2,by1,by2,rashift)
  use plot_data
  implicit none
  integer, parameter :: ns=9110, nsn=80
  integer :: i,j,c(100,35),nc,snr(nsn),plcst,plstar,spld,n,prslbl,rv
  real*8 :: ra(ns),dec(ns),d2r,r2d,r2h,pi,tpi,dx1,dx2,dy,ra1,dec1,drev2pi,par
  real :: pma,pmd,vm(ns),x1,y1,x2,y2,constx(99),consty(99),r1,g1,b1,r4,g4,b4
  real :: schcon,sz1,schfac,schlbl,prinf,snlim,sllim,schmag,getmag,mag,bx1,bx2,by1,by2,x,y,mlim,rashift
  character :: cn(100)*3,con(100)*20,name*10,vsopdir*99,sn(ns)*10,snam(nsn)*10,sni*10,getsname*10,mult,var*9
  
  mlim = 6.            !Magnitude limit for stars
  sllim = 2.5          !Limit for labels
  snlim = 1.4          !Limit for names
  fonttype = 2
  schmag = 0.07
  schlbl = fontsize1d
  schfac = 1.
  schcon = 1.
  plstar = 1  !0-no, 1-yes no label, 2-symbol, 3-name, 4-name or symbol, 5-name and symbol
  plcst = 2   !0-no, 1-figures, 2-figures+abbreviations, 3-figures+names
  
  prinf = 150.**2
  
  x = 0.
  call pgqcr(1,r1,g1,b1) !Store colours
  call pgqcr(4,r4,g4,b4)
  call pgscr(1,1.,1.,1.) !'White' (for stars)
  call pgscr(4,x,x,1.) !Blue (for constellations)
  
  pi = 4*datan(1.d0)
  tpi = 2*pi
  d2r = pi/180.d0
  r2d = 180.d0/pi
  r2h = 12.d0/pi
  r2h = r2d
  
  
  if(bx1.gt.bx2) then
     x = bx1
     bx1 = bx2
     bx2 = x
  end if
  
  !Read BSC
  vsopdir = '/home/sluys/diverse/popular/TheSky/'           !Linux pc
  open(unit=20,form='formatted',status='old',file=trim(vsopdir)//'data/bsc.dat')
  rewind(20)
  do i=1,ns
     read(20,320)name,ra(i),dec(i),pma,pmd,rv,vm(i),par,mult,var
320  format(A10,1x,2F10.6,1x,2F7.3,I5,F6.2,F6.3,A2,A10)
     sn(i) = getsname(name)
     ra(i) = mod(ra(i)+rashift,tpi)-rashift
  end do
  close(20)


  !Read Constellation figure data
  open(unit=40,form='formatted',status='old',file=trim(vsopdir)//'data/bsc_const.dat')
  do i=1,ns
     read(40,'(I4)',end=340,advance='no')c(i,1)
     do j=1,c(i,1)
        read(40,'(I5)',advance='no')c(i,j+1)
     end do
     read(40,'(1x,A3,A20)')cn(i),con(i)
     !Get mean star position to place const. name
     dx1 = 0.d0
     dx2 = 0.d0
     dy = 0.d0
     do j=2,c(i,1)
        dx1 = dx1 + dsin(ra(c(i,j)))
        dx2 = dx2 + dcos(ra(c(i,j)))
        dy = dy + dec(c(i,j))
     end do
     dx1 = (dx1 + dsin(ra(c(i,j))))/real(c(i,1))
     dx2 = (dx2 + dcos(ra(c(i,j))))/real(c(i,1))
     ra1 = drev2pi(datan2(dx1,dx2))
     dec1 = (dy + dec(c(i,j)))/real(c(i,1))
     !call eq2xy(ra1,dec1,l0,b0,x1,y1)
     !constx(i) = x1
     !consty(i) = y1
     !constx(i) = real(ra1*r2h)
     constx(i) = real((mod(ra1+rashift,tpi)-rashift)*r2h)
     consty(i) = real(dec1*r2d)
  end do
340 close(40)
  nc = i-1
  
  !Read Star names
  open(unit=50,form='formatted',status='old',file=trim(vsopdir)//'data/bsc_names.dat')
  do i=1,nsn
     read(50,'(I4,2x,A10)',end=350)snr(i),snam(i)
  end do
350 close(50)
  
  
  !!Read Milky Way data
  !do f=1,5
  !   write(mwfname,'(A10,I1,A4)')'milkyway_s',f,'.dat'
  !   open(unit=60,form='formatted',status='old',file=trim(vsopdir)//'data/'//mwfname)
  !   do i=1,mwn(f)
  !      read(60,'(F7.5,F9.5)')mwa(f,i),mwd(f,i)
  !      if(maptype.eq.1) call eq2az(mwa(f,i),mwd(f,i),agst)
  !      if(maptype.eq.2) call eq2ecl(mwa(f,i),mwd(f,i),eps)
  !   end do
  !end do
  !close(60)
  
  
  !Plot constellation figures
  if(plcst.gt.0) then
     !schcon = min(max(40./sz1,0.7),3.)
     call pgsch(schfac*schcon*schlbl)
     call pgscf(fonttype)
     call pgsci(4)
     call pgslw(2)
     do i=1,nc
        do j=2,c(i,1)
           !call eq2xy(ra(c(i,j)),dec(c(i,j)),l0,b0,x1,y1)
           !call eq2xy(ra(c(i,j+1)),dec(c(i,j+1)),l0,b0,x2,y2)
           x1 = real(ra(c(i,j))*r2h)
           y1 = real(dec(c(i,j))*r2d)
           x2 = real(ra(c(i,j+1))*r2h)
           y2 = real(dec(c(i,j+1))*r2d)
           !if((x1*x1+y1*y1.le.prinf.or.x2*x2+y2*y2.le.prinf).and.(x2-x1)**2+(y2-y1)**2.le.90.**2) & !Not too far from centre and each other 
           if((x2-x1)**2+(y2-y1)**2.le.90.**2)  call pgline(2,(/x1,x2/),(/y1,y2/))  !Not too far from centre and each other 
	end do
        if(constx(i).lt.bx1.or.constx(i).gt.bx2.or.consty(i).lt.by1.or.consty(i).gt.by2) cycle
        if(plcst.eq.2) call pgptext(constx(i),consty(i),0.,0.5,cn(i))
        if(plcst.eq.3) call pgptext(constx(i),consty(i),0.,0.5,con(i))
     end do
     call pgsch(schfac)
     call pgscf(fonttype)
  end if !if(plcst.gt.0) then
  
  !Plot stars: BSC
  spld = 0
  if(plstar.gt.0) then
     n = 0
     do i=1,ns
        if(vm(i).lt.mlim.and.vm(i).ne.0.) then
           !call eq2xy(ra(i),dec(i),l0,b0,x,y)
           x = real(ra(i)*r2h)
           y = real(dec(i)*r2d)
           if(x.lt.bx1.or.x.gt.bx2.or.y.lt.by1.or.y.gt.by2) cycle
           call pgsci(1)
           mag = getmag(vm(i),mlim)*schmag
           call pgcirc(x,y,mag)
           !write(6,'(3F10.3)')x,y,mag
           call pgsch(schfac*schlbl)
           sni = sn(i)
           !if(sni(1:1).eq.'\') call pgsch(schlbl*max(1.33,schfac))  !Greek letters need larger font
           if(sni(1:1).eq.char(92)) call pgsch(schlbl*max(1.33,schfac))  !Greek letters need larger font.  Char(92) is a \, but this way it doesn't mess up emacs' parentheses count
	   call pgsci(14)
           if(vm(i).lt.sllim) then
              if((plstar.eq.2.or.plstar.eq.5)) call pgtext(x+0.02*sz1,y+0.02*sz1,sn(i))
              if(plstar.eq.4) then !Check if the name will be printed
                 prslbl = 1
                 if(vm(i).lt.snlim) then
                    do j=1,nsn
                       if(snr(j).eq.i) prslbl = 0 !Then the name will be printed, don't print the symbol
                    end do
                 end if
                 if(prslbl.eq.1) call pgtext(x+0.02*sz1,y+0.02*sz1,sn(i))
              end if
           end if
	   spld = spld+1
	end if
     end do
     if(plstar.ge.3) then !Plot star proper names
        call pgsch(schfac*schlbl)
        do i=1,nsn
           if(vm(snr(i)).lt.max(snlim,1.4)) then  !Regulus (1.35) will still be plotted, for conjunction maps
              !call eq2xy(ra(snr(i)),dec(snr(i)),l0,b0,x,y)
              x = real(ra(snr(i))*r2h)
              y = real(dec(snr(i))*r2d)
              if(x.lt.bx1.or.x.gt.bx2.or.y.lt.by1.or.y.gt.by2) cycle
              call pgtext(x+0.02*sz1,y-0.02*sz1,snam(i))
           end if
        end do
     end if !if(plstar.eq.3) then
  end if !if(plstar.gt.0) then
  
  !Restore colours
  call pgscr(1,r1,g1,b1)
  call pgscr(4,r4,g4,b4)
  
end subroutine plotthesky
!************************************************************************************************************************************

!************************************************************************
function getsname(name)               !Get star name from bsc info
  use analysemcmc_settings
  implicit none
  character :: getsname*10,name*10,num*3,grk*3,gn*1
  num = name(1:3)
  grk = name(4:6)
  gn  = name(7:7)
  !      gn = ' '
  
  getsname = '          '
  if(grk.ne.'   ') then  !Greek letter
     if(fonttype.eq.2) then
        if(grk.eq.'Alp') getsname = '\(2127)\u'//gn
        if(grk.eq.'Bet') getsname = '\(2128)\u'//gn
        if(grk.eq.'Gam') getsname = '\(2129)\u'//gn
        if(grk.eq.'Del') getsname = '\(2130)\u'//gn
        if(grk.eq.'Eps') getsname = '\(2131)\u'//gn
        if(grk.eq.'Zet') getsname = '\(2132)\u'//gn
        if(grk.eq.'Eta') getsname = '\(2133)\u'//gn
        if(grk.eq.'The') getsname = '\(2134)\u'//gn
        if(grk.eq.'Iot') getsname = '\(2135)\u'//gn
        if(grk.eq.'Kap') getsname = '\(2136)\u'//gn
        if(grk.eq.'Lam') getsname = '\(2137)\u'//gn
        if(grk.eq.'Mu ') getsname = '\(2138)\u'//gn
        if(grk.eq.'Nu ') getsname = '\(2139)\u'//gn
        if(grk.eq.'Xi ') getsname = '\(2140)\u'//gn
        if(grk.eq.'Omi') getsname = '\(2141)\u'//gn
        if(grk.eq.'Pi ') getsname = '\(2142)\u'//gn
        if(grk.eq.'Rho') getsname = '\(2143)\u'//gn
        if(grk.eq.'Sig') getsname = '\(2144)\u'//gn
        if(grk.eq.'Tau') getsname = '\(2145)\u'//gn
        if(grk.eq.'Ups') getsname = '\(2146)\u'//gn
        if(grk.eq.'Phi') getsname = '\(2147)\u'//gn
        if(grk.eq.'Chi') getsname = '\(2148)\u'//gn
        if(grk.eq.'Psi') getsname = '\(2149)\u'//gn
        if(grk.eq.'Ome') getsname = '\(2150)\u'//gn
     else
        if(grk.eq.'Alp') getsname = '\(0627)\u'//gn
        if(grk.eq.'Bet') getsname = '\(0628)\u'//gn
        if(grk.eq.'Gam') getsname = '\(0629)\u'//gn
        if(grk.eq.'Del') getsname = '\(0630)\u'//gn
        if(grk.eq.'Eps') getsname = '\(0631)\u'//gn
        if(grk.eq.'Zet') getsname = '\(0632)\u'//gn
        if(grk.eq.'Eta') getsname = '\(0633)\u'//gn
        if(grk.eq.'The') getsname = '\(0634)\u'//gn
        if(grk.eq.'Iot') getsname = '\(0635)\u'//gn
        if(grk.eq.'Kap') getsname = '\(0636)\u'//gn
        if(grk.eq.'Lam') getsname = '\(0637)\u'//gn
        if(grk.eq.'Mu ') getsname = '\(0638)\u'//gn
        if(grk.eq.'Nu ') getsname = '\(0639)\u'//gn
        if(grk.eq.'Xi ') getsname = '\(0640)\u'//gn
        if(grk.eq.'Omi') getsname = '\(0641)\u'//gn
        if(grk.eq.'Pi ') getsname = '\(0642)\u'//gn
        if(grk.eq.'Rho') getsname = '\(0643)\u'//gn
        if(grk.eq.'Sig') getsname = '\(0644)\u'//gn
        if(grk.eq.'Tau') getsname = '\(0645)\u'//gn
        if(grk.eq.'Ups') getsname = '\(0646)\u'//gn
        if(grk.eq.'Phi') getsname = '\(0647)\u'//gn
        if(grk.eq.'Chi') getsname = '\(0648)\u'//gn
        if(grk.eq.'Psi') getsname = '\(0649)\u'//gn
        if(grk.eq.'Ome') getsname = '\(0650)\u'//gn
     end if
  else  !Then number
     if(num(1:1).eq.' ') num = num(2:3)//' '
     if(num(1:1).eq.' ') num = num(2:3)//' '
     getsname = num//'       '
  end if
  return
end function getsname
!************************************************************************

!************************************************************************
function getmag(m,mlim)  !Determine size of stellar 'disk'
  real :: getmag,m,m1,mlim
  m1 = m
  !      if(m1.lt.0.) m1 = m1*0.5  !Less excessive grow in diameter for the brightest objects
  if(m1.lt.-1.e-3) m1 = -sqrt(-m1)  !Less excessive grow in diameter for the brightest objects
  !getmag = max(mlim-m1+0.5,0.)
  getmag = max(mlim-m1+0.5,0.5) !Make sure the weakest stars are still plotted
  !getmag = max(mlim-m1+0.5,0.)+0.5
  return
end function getmag
!************************************************************************




!*****************************************************************************************************************************************************
subroutine pgimag_project(z,nbx,nby,xb1,xb2,yb1,yb2,z1,z2,tr,projection)  !Clone of pgimag, use projection
  use constants
  use general_data
  implicit none
  integer, parameter :: nell=100
  integer :: nbx,nby,xb1,xb2,yb1,yb2
  real :: z(nbx,nby),tr(6),z1,z2,dz,dcdz
  real :: x,y,dx,dy,xs(5),ys(5),xell(nell),yell(nell),sch
  integer :: i,ix,iy,clr1,clr2,dc,ci,projection,lw
  character :: str*99
  
  call pgqcir(clr1,clr2)
  dz = z2-z1
  dc = clr2-clr1
  dcdz = real(dc)/dz
  
  call pgbbuf  !Buffer output to speed up screen plotting
  dx = tr(2)/2.*1.05  !Distance between pixel centres / 2 = half width of pixels
  dy = tr(6)/2.*1.05  !Spaces between pixels seem to go away when multiplying with a factor between 1.02 and 1.04
  
  
  !Loop over pixels (each dimension has one array row/column too many)
  do ix = xb1,xb2-1
     do iy = yb1,yb2-1
        
        !Get colour for this pixel:
        dz = z(ix,iy)-z1
        ci = min(clr1 + nint(dz*dcdz),clr2)
        if(ci.eq.clr1) cycle  !Don't draw background pixels
        
        call pgsci(ci)
        
        !Get central coordinates for this pixel:
        x = tr(1) + tr(2)*ix + tr(3)*iy
        y = tr(4) + tr(5)*ix + tr(6)*iy
        
        !Get the coordinates of the four corners (projected rectangular pixel is not necessarily rectangular!)
        xs(1) = x-dx
        ys(1) = y-dy
        xs(2) = xs(1)
        ys(2) = y+dy
        xs(3) = x+dx
        ys(3) = ys(2)
        xs(4) = xs(3)
        ys(4) = ys(1)
        
        !Do the projection:
        if(projection.ge.1) then
           do i=1,4
              call project_skymap(xs(i),ys(i),racentre,projection)
           end do
        end if
        xs(5) = xs(1)
        ys(5) = ys(1)
        
        !Plot the pixel:
        call pgpoly(5,xs,ys)
        
     end do
  end do
  
  
  !Draw lines on map:
  if(projection.eq.1) then
     call pgqch(sch) !Save current character height
     call pgsch(0.5*sch)
     
     !Get data to plot ellipses:
     do i=1,nell
        x = real(i-1)/real(nell-1)*rtpi
        xell(i) = sin(x)
        yell(i) = cos(x)
     end do
     call pgsci(1)
     
     !Plot meridians:
     do i=-24,24,3  !In hours
        call pgsci(14)
        if(i.eq.0.or.abs(i).eq.24) call pgsci(1) !Null-meridian in black
        if(real(i).gt.-racentre-12.and.real(i).lt.-racentre+12) then
           !Plot line:
           x = -(racentre+real(i))
           call pgline(nell/2+1,xell*x+racentre,yell*90.)
           
           !Print label:
           write(str,'(I2,A)')mod(48-i,24),'\uh\d'
           if(mod(48-i,24).lt.10) write(str,'(I1,A)')mod(48-i,24),'\uh\d'
           call pgptext(x+racentre-0.1,2.,0.,0.,trim(str))
        end if
     end do
     
     !Plot lines of constant declination:
     do i=-90,90,15  !In degrees
        if(abs(i).eq.90) cycle
        call pgsci(14)
        if(i.eq.0) call pgsci(1) !Equator in black
        
        !Get start and end point on line and project them:
        xs(1) = racentre-12.
        xs(2) = racentre+12.
        ys(1) = real(i)
        ys(2) = real(i)
        call project_skymap(xs(1),ys(1),racentre,projection)
        call project_skymap(xs(2),ys(2),racentre,projection)
        
        !Plot line:
        call pgline(2,xs(1:2),ys(1:2))
        
        !Print labels:
        if(i.gt.0) then
           write(str,'(A1,I2,A)')'+',i,'\(2218)'
           call pgptext(xs(2)+0.2,ys(1),0.,1.,trim(str))
        else
           write(str,'(I3,A)')i,'\(2218)'
           call pgptext(xs(2)+0.4,ys(1)-2.,0.,1.,trim(str))
        end if
     end do
     
     !Overplot main lines:
     call pgsci(1)
     call pgptext(racentre,92.,0.,0.5,'+90\(2218)')           !NP
     call pgptext(racentre,-95.,0.,0.5,'-90\(2218)')          !SP
     call pgline(2,(/racentre-12.,racentre+12./),(/0.,0./))   !Equator
     
     !Plot null-meridian:
     do i=-24,24,24
        if(real(i).gt.-racentre-12.and.real(i).lt.-racentre+12) call pgline(nell/2+1,-(racentre+real(i))*xell+racentre,yell*90.)
     end do
     
     call pgqlw(lw)     !Save current line width
     call pgslw(lw*2)
     call pgline(nell,xell*12.+racentre,yell*90.)             !Outline
     call pgslw(lw)     !Restore line width
     
     
     
  end if  !if(projection.eq.1)
  
  call pgsch(sch)  !Restore character height
  call pgebuf      !Release buffer
  
end subroutine pgimag_project
!*****************************************************************************************************************************************************




!*****************************************************************************************************************************************************
subroutine project_skymap(x,y,racentre,projection)  !Project a sky map, using projection 'projection'
  use constants
  implicit none
  integer :: projection
  real :: x,y,theta,siny,th2,dth2,delta,racentre
  
  if(projection.eq.1) then
     !Mollweide projection:
     !http://en.wikipedia.org/wiki/Mollweide_projection
     !Newton-Rapson scheme to solve equation:  2*theta + sin(2*theta) = pi*sin(y*rd2r)
     !Convergence is relatively fast, somewhat slower near poles
     
     delta = 1.e-6        !Radians
     siny  = sin(y*rd2r)
     th2 = y*rd2r
     dth2 = 1.e30
     
     do while(abs(dth2).gt.delta)
        dth2 = -(th2 + sin(th2) - rpi*siny)/(1.+cos(th2))
        th2 = th2 + dth2
     end do
     
     theta = th2/2.
     
     !Original projection:
     !x = 2*sqrt2/rpi * x * cos(theta)
     !y = sqrt(2) * sin(theta) * r2d
     
     !Map it back to a 24hx180d plot:
     x = (x-racentre) * cos(theta) + racentre
     y = sin(theta)*90.
  else
     write(0,'(A,I3)')'  ERROR:  Projection not defined:',projection
     write(0,'(A)')'  Aborting...'
     stop
  end if
  
end subroutine project_skymap
!*****************************************************************************************************************************************************
