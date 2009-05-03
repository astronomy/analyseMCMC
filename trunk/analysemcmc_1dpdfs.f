!Routines to compute and plot one-dimensional PDFs


subroutine pdfs1d(exitcode)
  use constants
  use analysemcmc_settings
  use general_data
  use mcmcrun_data
  use plot_data
  implicit none
  
  integer :: i,j,p,ic,io,pgopen,lw,exitcode,system
  real :: rev24,rev360,rev180
  real :: x(nchs,nchs*narr1),xmin,xmax,xmin1,xmax1,xpeak,dx,ymin,ymax,sch
  real,allocatable :: xbin(:,:),ybin(:,:),xbin1(:),ybin1(:),ybin2(:),ysum(:),yconv(:),ycum(:)  !These depend on nbin1d, allocate after reading input file
  real :: plshift,plx,ply,x0,norm,bindx
  character :: string*99,str*99,str1*99,str2*99
  
  exitcode=0
  if(prprogress.ge.1.and.plot.eq.0.and.savepdf.eq.1) write(*,'(A,$)')'  Saving 1D pdfs'
  if(prprogress.ge.1.and.plot.eq.1.and.update.eq.0) write(*,'(A,$)')' 1D pdfs'
  
  !Autodetermine number of bins:
  if(nbin1d.le.0) then
     if(totpts.le.100) then
        nbin1d = floor(2*sqrt(real(totpts)))
     else
        nbin1d = floor(10*log10(real(totpts)))
     end if
     nbin1d = max(nbin1d,5)
     if(prprogress.ge.2.and.plot.eq.1.and.update.eq.0) then
        if(nbin1d.lt.100) write(*,'(A2,I2,A8,$)')' (',nbin1d,' bins), '
        if(nbin1d.ge.100) write(*,'(A2,I3,A8,$)')' (',nbin1d,' bins), '
     end if
  else
     nbin1d = max(nbin1d,5)
     if(prprogress.ge.1.and.plot.eq.1.and.update.eq.0) write(*,'(A2,$)')', '
  end if

  !Allocate memory:
  allocate(xbin(nchs,nbin1d+1),ybin(nchs,nbin1d+1),xbin1(nbin1d+1),ybin1(nbin1d+1),ybin2(nbin1d+1),ysum(nbin1d+1),yconv(nbin1d+1),ycum(nbin1d+1))

  if(plot.eq.1) then
     if(file.eq.0) then
        io = pgopen('14/xs')
        sch = 1.5*fontsize1d
        lw = 1
     end if
     if(file.ge.1) then
        if(file.eq.1) io = pgopen('pdfs.ppm/ppm')
        if(file.ge.2) io = pgopen('pdfs.eps'//trim(psclr))
        lw = 3
        if(nplvar.ge.10) lw = 2
        if(quality.lt.2) lw = max(lw-1,1)  !Draft/Paper
        sch = 1.2*fontsize1d !2.5
        if(nchains.eq.1.and.nplvar.gt.9) sch = 1.2*fontsize1d
        if(quality.eq.0) then !Draft
           sch = sch*1.75
           lw = 2
        end if
        if(quality.eq.1) then !Paper
           if(nplvar.eq.12) then
              sch = sch*1.75
              lw = 2
           else
              sch = sch*1.25
              lw = 1
           end if
        end if
        if(quality.eq.2) then !Talk
           if(nplvar.le.12) then
              sch = sch*2
              lw = 2
           else
              sch = sch*1.5
              lw = 1
           end if
        end if
        if(quality.eq.3) then !Poster
           if(nplvar.eq.12.and.file.ge.2) then
              sch = sch*2.7
              lw = 3
           else
              !sch = sch*1.25
              !lw = 1
              sch = sch*1.5
              lw = 2
           end if
        end if
        if(quality.eq.4) then !Vivien's thesis
           sch = sch*2.5
           lw = 2
        end if
     end if !if(file.ge.1)
     if(io.le.0) then
        write(*,'(A,I4)')'Cannot open PGPlot device.  Quitting the programme',io
        exitcode = 1
        return
     end if
     if(file.eq.0) call pgpap(scrsz,scrrat)
     if(file.eq.1) call pgpap(bmpsz,bmprat)
     if(file.ge.2) call pgpap(pssz,psrat)
     !if(file.ge.2.and.quality.eq.3.and.nplvar.eq.12) call pgpap(10.6,0.925)
     if(file.ge.2.and.quality.eq.3.and.nplvar.eq.12) call pgpap(10.6,0.85)
     if(file.ge.2) call pgscf(fonttype)
     !call pgscr(3,0.,0.5,0.)
     !call pginitl(colour,file,whitebg)
     call pgslw(lw)
     call pgsch(sch)
     call pgsfs(fillpdf)


     !if(quality.eq.0) call pgsvp(0.08,0.95,0.06,0.87) !To make room for title

     call pgsubp(panels(1),panels(2))
  end if !if(plot.eq.1)

  !Save 1D PDF data
  if(savepdf.eq.1) then
     open(unit=30,action='write',form='formatted',status='replace',file=trim(outputdir)//'/'//trim(outputname)//'__pdf1d.dat')
     write(30,'(3I6,T100,A)')nplvar,nchains,nbin1d,'Total number of plot variables, total number of chains, number of bins'
     !write(30,'(3I6,T100,A)')nplvar-nfixedpar,nchains,nbin1d,'Total number of plot variables, total number of chains, number of bins'
  end if

  !do p=par1,par2
  do j=1,nplvar
     p = plvars(j)
     if(plot.eq.1) then
        call pgpage
        if(j.eq.1) call pginitl(colour,file,whitebg)
     end if
     
     !Set x-ranges for plotting, bin the data and get y-ranges
     !Use widest probability range (hopefully ~3-sigma) - doesn't always work well...
     if(1.eq.2.and.version.eq.1) then  !This can only be used if ranges are computed - isn't this always the case?
        xmin = 1.e30
        xmax = -1.e30
        do ic=1,nchains
           xmin = min(xmin,ranges(ic,nival,p,1))
           xmax = max(xmax,ranges(ic,nival,p,2))
           !write(*,'(3I4,6F10.3)')ic,p,nival,xmin,xmax,ranges(ic,nival,p,1),ranges(ic,nival,p,2)
        end do
     end if
     !print*,xmin,huge(xmin)
     !if(xmin.le.-huge(xmin).or.xmin.ge.huge(xmin).or.xmax.le.-huge(xmax).or.xmax.ge.huge(xmax)) then !NaN
     !if(version.eq.2.or.xmin.ne.xmin.or.xmax.ne.xmax) then
        xmin = 1.e30
        xmax = -1.e30
        do ic=1,nchains
           if(mergechains.eq.0.and.contrchain(ic).eq.0) cycle
           xmin = min(xmin,minval(alldat(ic,p,1:n(ic))))
           xmax = max(xmax,maxval(alldat(ic,p,1:n(ic))))
        end do
     !end if
     dx = xmax - xmin
     !dx = max(xmax - xmin,1.e-30)
     
     do ic=1,nchains
        if(mergechains.eq.0.and.contrchain(ic).eq.0) cycle
        x(ic,1:n(ic)) = alldat(ic,p,1:n(ic))
        xmin1 = minval(alldat(ic,p,1:n(ic)))
        xmax1 = maxval(alldat(ic,p,1:n(ic)))
        if(wrap(ic,p).eq.0) then !Use the plotting ranges
           !xmin1 = xmin - 0.1*dx
           !xmax1 = xmax + 0.1*dx
        end if
        
        call bindata1d(n(ic),x(ic,1:n(ic)),1,nbin1d,xmin1,xmax1,xbin1,ybin1)
        
        !Weight with likelihood.  I should probably do something like this at the start, to get updated ranges etc.
        !y(ic,1:n(ic)) = alldat(ic,1,1:n(ic))
        !call bindata1da(n(ic),x(ic,1:n(ic)),y(ic,1:n(ic)),1,nbin1d,xmin1,xmax1,xbin1,ybin1) !Measure the amount of likelihood in each bin
        
        !Columns in dat(): 1:logL 2:mc, 3:eta, 4:tc, 5:d_l, 6:spin, 7:theta_SL, 8: RA, 9:dec,10:phase, 11:thetaJ0, 12:phiJ0, 13:alpha, 14:M1, 15:M2
        if(p.eq.5.or.p.eq.7.or.p.eq.9.or.p.eq.11) then  !Do something about chains 'sticking to the wall'
           if(ybin1(1).gt.ybin1(2)) ybin1(1)=0.
           if(ybin1(nbin1d).gt.ybin1(nbin1d-1)) ybin1(nbin1d)=0.
        end if
        
        !Normalise 1D PDF
        if(normpdf1d.gt.0) then
           if(normpdf1d.eq.1) then !Normalise the SURFACE, not the height (because of different bin size).  This is the default
              norm = 0.
              do i=1,nbin1d+1
                 norm = norm + ybin1(i)
              end do
              norm = norm*(xmax1-xmin1)
              ybin1 = ybin1/norm
           else !Normalise to the height of the PDF
              if(normpdf1d.eq.2) ybin1 = ybin1/maxval(ybin1)  !Normalise to the height of the PDF
              if(normpdf1d.eq.3) ybin1 = ybin1/(maxval(ybin1)**0.5)  !Normalise to the sqrt of the height of the PDF; Works nicely for comparing parallel-tempering chains
              if(ic*j.eq.1)write(*,'(//,A,/)')'  *** WARNING: using non-default normalisation for PDFs ***'
           end if
        end if
        
        !Smoothen 1D PDF
        ybin2 = ybin1
        if(smooth.gt.1) call smoothpdf1d(ybin1,nbin1d+1,smooth)
        xbin(ic,1:nbin1d+1) = xbin1(1:nbin1d+1)
        ybin(ic,1:nbin1d+1) = ybin1(1:nbin1d+1)
        
        
        !Save binned data
        if(savepdf.eq.1) then
           !if(savepdf.eq.1.and.fixedpar(p).eq.0) then
           if(fixedpar(p).eq.1) ybin1 = 0.  !Prevent NaNs
           write(30,'(3I6,T100,A)')ic,p,wrap(ic,p),'Chain number, variable number, and wrap'
           write(30,'(2ES15.7,T100,A)')startval(ic,p,1:2),'True and starting value'
           write(30,'(6ES15.7,T100,A)')stats(ic,p,1:6),'Stats: median, mean, absvar1, absvar2, stdev1, stdev2'
           write(30,'(5ES15.7,T100,A)')ranges(ic,c0,p,1:5),'Ranges: lower,upper limit, centre, width, relative width'
           write(30,'(2ES15.7,T100,A)')xmin1,xmax1,'Xmin and Xmax of PDF'
           do i=1,nbin1d+1
              write(30,'(2ES15.7)')xbin1(i),ybin1(i)
           end do
        end if
     end do !ic
     
     if(plot.eq.1) then
        !Ranges for plot panel
        xmin = xmin - 0.1*dx
        xmax = xmax + 0.1*dx
        ymin = 0.
        ymax = 1.e-20
        !print*,xmin,xmax,ymin,ymax
        do ic=1,nchains
           if(mergechains.eq.0.and.contrchain(ic).eq.0) cycle
           !ymax = max(ymax,maxval(ybin(ic,1:nbin1d+1)))
           do i=1,nbin1d+1
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
        
        
        if(file.eq.0.and.scrrat.gt.1.35) call pgsvp(0.08,0.95,0.1,0.95)
        if(file.eq.1.and.bmprat.gt.1.35) call pgsvp(0.08,0.95,0.1,0.95)
        if(file.ge.2.and.psrat.gt.1.35) call pgsvp(0.08,0.95,0.1,0.95)
        if(quality.eq.0) call pgsvp(0.08,0.95,0.06,0.87) !To make room for title
        if(quality.eq.4) call pgsvp(0.13,0.95,0.1,0.95)

        call pgsch(sch)
        call pgswin(xmin,xmax,ymin,ymax)
        if(abs(dx).lt.1.e-30) then !So that the program doesn't hang if a parameter is kept constant
           xbin = 0.
           ybin = 0.
        end if
        
        !Plot 1D PDF
        !if(fixedpar(p).eq.0) then
        call pgsci(1)
        if(file.ge.2) call pgslw(lw)
        do ic=1,nchains
           if(mergechains.eq.0.and.contrchain(ic).eq.0) cycle
           if(fillpdf.ge.3) call pgshs(45.0*(-1)**ic,2.0,real(ic)/real(nchains0)) !Set hatch style: angle = +-45deg, phase between 0 and 1 (1/nchains0, 2/nchains0, ...)
           if(nchains.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
           xbin1(1:nbin1d+1) = xbin(ic,1:nbin1d+1)
           ybin1(1:nbin1d+1) = ybin(ic,1:nbin1d+1)
           if(wrap(ic,p).eq.0) then
              if(nchains.eq.1) call pgsci(15)
              if(plpdf1d.eq.1) then
                 call pgpoly(nbin1d+2,(/xbin1(1),xbin1(1:nbin1d+1)/)+bindx/2.,(/0.,ybin1(1:nbin1d+1)/))
              else
                 call verthist(nbin1d+2,(/xbin1(1),xbin1(1:nbin1d+1)/),(/0.,ybin1(1:nbin1d+1)/),2)
              end if
              !Plot pdf contour
              !if(nchains.eq.1) call pgsci(1)
              !call pgsci(1)
              if(fillpdf.eq.1) call pgsci(1)
              if(nchains.eq.1) call pgsci(2)
              if(plpdf1d.eq.1) then
                 call pgline(nbin1d+1,xbin1(1:nbin1d+1)+bindx/2.,ybin1(1:nbin1d+1)) !:nbin1d) ?
              else
                 call verthist(nbin1d+2,(/xbin1(1),xbin1(1:nbin1d+1)/),(/0.,ybin1(1:nbin1d+1)/),1)
              end if
              
              !Fix the loose ends
              call pgline(2,(/xbin1(1)-bindx,xbin1(1)/)+bindx/2.,(/0.,ybin1(1)/))
              call pgline(2,(/xbin1(nbin1d+1),xbin1(nbin1d+1)/)+bindx/2.,(/ybin1(nbin1d+1),0./))
           else !If parameter is wrapped
              plshift = real(2*pi)
              if(changevar.eq.1) plshift = 360.
              if(changevar.eq.1.and.p.eq.8) plshift = 24.  !RA in hours
              if(changevar.eq.1.and.p.eq.12) plshift = 180.  !Pol.angle
              if(nchains.eq.1) call pgsci(15)
              if(plpdf1d.eq.1) then
                 call pgpoly(nbin1d+3,(/xbin1(1),xbin1(1:nbin1d),xbin1(1)+plshift,xbin1(1)+plshift/),(/0.,ybin1(1:nbin1d),ybin1(1),0./))
              else
                 call verthist(nbin1d+3,(/xbin1(1),xbin1(1:nbin1d),xbin1(1)+plshift,xbin1(1)+plshift/),(/0.,ybin1(1:nbin1d),ybin1(1),0./),2)
              end if
              !Plot pdf contour
              !call pgsci(1)
              !if(fillpdf.ne.1) call pgsci(colours(mod(ic-1,ncolours)+1))
              if(fillpdf.eq.1) call pgsci(1)
              if(nchains.eq.1) call pgsci(2)
              if(plpdf1d.eq.1) then
                 call pgline(nbin1d,xbin1(1:nbin1d)+bindx/2.,ybin1(1:nbin1d))
              else
                 call verthist(nbin1d,xbin1(1:nbin1d)+bindx/2.,ybin1(1:nbin1d),0)
              end if
              
              !Plot dotted lines outside the pdf for wrapped periodic variables
              call pgsls(4)
              !if(file.ge.2) call pgslw(2)
              if(plpdf1d.eq.1) then
                 call pgline(nbin1d+1,(/xbin1(1:nbin1d)-plshift,xbin1(1)/)+bindx/2.,(/ybin1(1:nbin1d),ybin1(1)/))
                 call pgline(nbin1d,xbin1+plshift+bindx/2.,ybin1)
              else
                 call verthist(nbin1d+1,(/xbin1(1:nbin1d)-plshift,xbin1(1)/)+bindx/2.,(/ybin1(1:nbin1d),ybin1(1)/),0)
                 call verthist(nbin1d,xbin1+plshift+bindx/2.,ybin1,0)
              end if
              
              !Fix the loose end
              call pgsls(1)
              if(file.ge.2) call pgslw(lw)
              if(plpdf1d.eq.1) then
                 call pgline(2,(/xbin1(nbin1d),xbin1(1)+plshift/)+bindx/2.,(/ybin1(nbin1d),ybin1(1)/))
              else
                 call verthist(2,(/xbin1(nbin1d),xbin1(1)+plshift/)+bindx/2.,(/ybin1(nbin1d),ybin1(1)/),0)
              end if
           end if
        end do !ic
        
        
        !Plot lines again over surface of overlapping distributions
        if(nchains.gt.1.and.fillpdf.eq.1) then
           call pgsls(4)
           do ic=1,nchains
              if(mergechains.eq.0.and.contrchain(ic).eq.0) cycle
              call pgsci(1)
              !call pgsci(colours(mod(ic-1,ncolours)+1))
              xbin1(1:nbin1d+1) = xbin(ic,1:nbin1d+1)
              ybin1(1:nbin1d+1) = ybin(ic,1:nbin1d+1)
              if(wrap(ic,p).eq.0) then
                 if(plpdf1d.eq.1) then
                    call pgline(nbin1d+1,xbin1(1:nbin1d+1)+bindx/2.,ybin1(1:nbin1d+1))
                 else
                    call verthist(nbin1d+1,xbin1(1:nbin1d+1)+bindx/2.,ybin1(1:nbin1d+1),0)
                 end if
              else
                 if(plpdf1d.eq.1) then
                    call pgline(nbin1d,xbin1(1:nbin1d)+bindx/2.,ybin1(1:nbin1d))
                 else
                    call verthist(nbin1d,xbin1(1:nbin1d)+bindx/2.,ybin1(1:nbin1d),0)
                 end if
              end if
           end do
           call pgsls(4)
        end if
        
        
        !Plot max likelihood
        if(pllmax.ge.1) then
           ply = pldat(icloglmax,p,iloglmax)
           if(p.eq.8) ply = rev24(ply)
           if(p.eq.10.or.p.eq.13) ply = rev360(ply)
           if(p.eq.12) ply = rev180(ply)
           call pgsci(1)
           call pgsls(5)
           call pgline(2,(/ply,ply/),(/-1.e20,1.e20/))
        end if
        
        
        !Plot median and model value
        call pgsch(sch)
        
        do ic=1,nchains
           if(mergechains.eq.0.and.contrchain(ic).eq.0) cycle
           !Draw white lines
           if(nchains.gt.1) then
              call pgslw(lw)
              call pgsls(1); call pgsci(0)
              if(pltrue.eq.1.or.pltrue.eq.3) call pgline(2,(/startval(ic,p,1),startval(ic,p,1)/),(/-1.e20,1.e20/))                    !True value
              !if((pltrue.eq.2.or.pltrue.eq.4).and.version.eq.1.and.(p.eq.2.or.p.eq.3.or.p.eq.4.or.p.eq.6.or.p.eq.7.or.p.eq.14.or.p.eq.15)) call pgline(2,(/startval(ic,p,1),startval(ic,p,1)/),(/-1.e20,1.e20/))                    !True value - mass and spin only
              if((pltrue.eq.2.or.pltrue.eq.4).and.version.eq.1.and.(p.eq.2.or.p.eq.3.or.p.eq.4.or.p.eq.6.or.p.eq.14.or.p.eq.15)) call pgline(2,(/startval(ic,p,1),startval(ic,p,1)/),(/-1.e20,1.e20/))                    !True value - mass and spin only
              !if(plstart.ge.1) call pgline(2,(/startval(ic,p,2),startval(ic,p,2)/),(/-1.e20,1.e20/))                   !Starting value
              if(p.ne.1) then !Not if plotting log(L)
                 if(plmedian.eq.1.or.plmedian.eq.3.or.plmedian.eq.4.or.plmedian.eq.6) call pgline(2,(/stats(ic,p,1),stats(ic,p,1)/),(/-1.e20,1.e20/))                          !Median
                 if(plrange.eq.1.or.plrange.eq.3.or.plrange.eq.4.or.plrange.eq.6) then
                    if(nchains.lt.2) call pgline(2,(/ranges(ic,c0,p,1),ranges(ic,c0,p,1)/),(/-1.e20,1.e20/)) !Left limit of 90% interval
                    if(nchains.lt.2) call pgline(2,(/ranges(ic,c0,p,2),ranges(ic,c0,p,2)/),(/-1.e20,1.e20/)) !Right limit of 90% interval
                    if(nchains.eq.1) call pgline(2,(/ranges(ic,c0,p,3),ranges(ic,c0,p,3)/),(/-1.e20,1.e20/)) !Centre of 90% interval
                 end if
              end if
           end if

           call pgslw(lw+1)
           !Draw coloured lines over the white ones
           !Median
           if((plmedian.eq.1.or.plmedian.eq.3.or.plmedian.eq.4.or.plmedian.eq.6) .and. p.ne.1) then
              call pgsls(2); call pgsci(2); if(nchains.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
              call pgline(2,(/stats(ic,p,1),stats(ic,p,1)/),(/-1.e20,1.e20/))
           end if

           !Plot ranges in 1D PDF
           if((plrange.eq.1.or.plrange.eq.3.or.plrange.eq.4.or.plrange.eq.6) .and. p.ne.1) then
              call pgsls(4); call pgsci(2); if(nchains.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
              if(nchains.lt.2) call pgline(2,(/ranges(ic,c0,p,1),ranges(ic,c0,p,1)/),(/-1.e20,1.e20/)) !Left limit of 90% interval
              if(nchains.lt.2) call pgline(2,(/ranges(ic,c0,p,2),ranges(ic,c0,p,2)/),(/-1.e20,1.e20/)) !Right limit of 90% interval
              !if(nchains.eq.1) call pgline(2,(/ranges(ic,c0,p,3),ranges(ic,c0,p,3)/),(/-1.e20,1.e20/)) !Centre of 90% interval
           end if
           
           !Plot true value in PDF
           !if(pltrue.ge.1) then !Plot true values
           !if(pltrue.eq.1.or.pltrue.eq.3.or.((pltrue.eq.2.or.pltrue.eq.4).and.version.eq.1.and.(p.eq.2.or.p.eq.3.or.p.eq.4.or.p.eq.6.or.p.eq.7.or.p.eq.14.or.p.eq.15))) then
           if(pltrue.eq.1.or.pltrue.eq.3.or.((pltrue.eq.2.or.pltrue.eq.4).and.version.eq.1.and.(p.eq.2.or.p.eq.3.or.p.eq.4.or.p.eq.6.or.p.eq.14.or.p.eq.15))) then
              if(mergechains.ne.1.or.ic.le.1) then !The units of the true values haven't changed (e.g. from rad to deg) for ic>1 (but they have for the starting values, why?)
                 call pgsls(2); call pgsci(1)
                 if(pllmax.eq.0) call pgsls(3)  !Dash-dotted line when Lmax line isn't plotted (should we do this always?)
                 plx = startval(ic,p,1)
                 if(p.eq.8) plx = rev24(plx)
                 if(p.eq.10.or.p.eq.13) ply = rev360(ply)
                 if(p.eq.12) ply = rev180(ply)
                 call pgline(2,(/plx,plx/),(/-1.e20,1.e20/)) !True value
                 if(p.eq.8) then
                    call pgline(2,(/plx-24.,plx-24./),(/-1.e20,1.e20/)) !True value
                    call pgline(2,(/plx+24.,plx+24./),(/-1.e20,1.e20/)) !True value
                 end if
                 if(p.eq.10.or.p.eq.13) then
                    call pgline(2,(/plx-360.,plx-360./),(/-1.e20,1.e20/)) !True value
                    call pgline(2,(/plx+360.,plx+360./),(/-1.e20,1.e20/)) !True value
                 end if
                 if(p.eq.12) then
                    call pgline(2,(/plx-180.,plx-180./),(/-1.e20,1.e20/)) !True value
                    call pgline(2,(/plx+180.,plx+180./),(/-1.e20,1.e20/)) !True value
                 end if
              end if
           end if

           !Plot starting value in 1D PDF
           !if(plstart.eq.1.and.abs((startval(ic,p,1)-startval(ic,p,2))/startval(ic,p,1)).gt.1.e-10) then
           !   call pgsls(4); call pgsci(1); if(nchains.gt.1) call pgsci(1)
           !   call pgline(2,(/startval(ic,p,2),startval(ic,p,2)/),(/-1.e20,1.e20/))
           !end if

           call pgsls(1)
           call pgsci(1)
        end do !ic




        !Print median, model value and range widths in 1D PDF panel title
        call pgslw(lw)
        call pgsci(1)
        ic = 1
        !if(nplvar.lt.7.or.nplvar.eq.9) then  !Three or less columns
        if(quality.ne.2.and.quality.ne.3.and.quality.ne.4) then  !Not a talk/poster/thesis
           if(nplvar.le.5) then
              write(str,'(A,F7.3,A5,F7.3)')trim(pgvarns(p))//': mdl:',startval(ic,p,1),' med:',stats(ic,p,1)
              if(prival.ge.1) then
                 if(plrange.eq.4.or.plrange.eq.5.or.plrange.eq.6) then
                    !if(p.eq.2.or.p.eq.3.or.p.eq.5.or.p.eq.6.or.p.eq.14.or.p.eq.15) then
                    if(p.eq.2.or.p.eq.5.or.p.eq.14.or.p.eq.15) then
                       write(str,'(A,F6.2,A1)')trim(str)//' \(2030):',ranges(ic,c0,p,5)*100,'%'
                    else
                       write(str,'(A,F7.3)')trim(str)//' \(2030):',ranges(ic,c0,p,5)
                    end if
                 end if
              end if
           else  !if nplvar>=5
              str = ' '
              !if(p.eq.2.or.p.eq.3.or.p.eq.5.or.p.eq.6.or.p.eq.14.or.p.eq.15) then
              if(p.eq.2.or.p.eq.5.or.p.eq.14.or.p.eq.15) then
                 if(pltrue.eq.3.or.pltrue.eq.4) write(str,'(A,F7.3)')trim(str)//' mdl:',startval(ic,p,1)
                 if(plmedian.eq.4.or.plmedian.eq.5.or.plmedian.eq.6) write(str,'(A,F8.3)')trim(str)//' med:',stats(ic,p,1)
                 if(prival.ge.1.and.(plrange.eq.4.or.plrange.eq.5.or.plrange.eq.6)) write(str,'(A,F6.2,A1)')trim(str)//' \(2030):',ranges(ic,c0,p,5)*100,'%'
              else
                 if(pltrue.eq.3.or.pltrue.eq.4) write(str,'(A,F7.3)')trim(str)//' mdl:',startval(ic,p,1)
                 if(plmedian.eq.4.or.plmedian.eq.5.or.plmedian.eq.6) write(str,'(A,F8.3)')trim(str)//' med:',stats(ic,p,1)
                 if(prival.ge.1.and.(plrange.eq.4.or.plrange.eq.5.or.plrange.eq.6)) write(str,'(A,F7.3)')trim(str)//' \(2030):',ranges(ic,c0,p,5)
              end if
              call pgsch(sch*1.2)
              call pgptxt(xmin+0.05*dx,ymax*(1.0-0.1*fontsize1d),0.,0.,trim(pgvarnss(p)))
           end if
        end if
        
        
        if(quality.eq.2.or.quality.eq.3.or.quality.eq.4) then  !Talk/poster/thesis quality for 1D PDF
           !if(p.eq.2.or.p.eq.3.or.p.eq.5.or.p.eq.6.or.p.eq.14.or.p.eq.15) then
           !   write(str,'(A9,F6.2,A1)')' \(2030):',ranges(ic,c0,p,5)*100,'%'
           !else
           !   write(str,'(A9,F7.3)')' \(2030):',ranges(ic,c0,p,5)
           !end if
           if(plrange.eq.1.or.plrange.eq.3.or.plrange.eq.4.or.plrange.eq.6) then
              x0 = ranges(ic,c0,p,5)
              if(p.eq.2.or.p.eq.3.or.p.eq.5.or.p.eq.6.or.p.eq.14.or.p.eq.15) x0 = x0*100
              !print*,p,x0,nint(x0)
              if(x0.lt.0.01) write(str,'(F6.4)')x0
              if(x0.ge.0.01.and.x0.lt.0.1) write(str,'(F5.3)')x0
              if(x0.ge.0.1.and.x0.lt.1.) write(str,'(F4.2)')x0
              if(x0.ge.1.and.x0.lt.9.95) write(str,'(F3.1)')x0
              if(x0.ge.9.95.and.x0.lt.99.5) write(str,'(I2)')nint(x0)
              if(x0.ge.99.5) write(str,'(I3)')nint(x0)
              write(str,'(A)')'\(2030): '//trim(str)
              !if(p.eq.2.or.p.eq.3.or.p.eq.5.or.p.eq.6.or.p.eq.14.or.p.eq.15) then
              if(p.eq.2.or.p.eq.5.or.p.eq.14.or.p.eq.15) then
                 write(str,'(A)')trim(str)//'%'
              else
                 write(str,'(A)')trim(str)//trim(pgunits(p))
              end if
           else if(pltrue.eq.2.or.pltrue.eq.4) then  !If not plotting ranges, but do plot true values
              x0 = startval(ic,p,1)
              !print*,p,x0,nint(x0)
              if(x0.lt.0.01) write(str,'(F7.4)')x0
              if(x0.ge.0.01.and.x0.lt.0.1) write(str,'(F6.3)')x0
              if(x0.ge.0.1.and.x0.lt.1.) write(str,'(F6.3)')x0
              if(x0.ge.1.and.x0.lt.9.995) write(str,'(F6.3)')x0
              if(x0.ge.9.995.and.x0.lt.99.9) write(str,'(F5.1)')x0
              if(x0.ge.99.9) write(str,'(F6.1)')x0
              write(str,'(A)')'true: '//trim(str)//trim(pgunits(p))
           end if
           
           
           !Print variable name in top of panel:
           call pgsch(sch*1.2)
           if(abs(xmin-xpeak).lt.abs(xmax-xpeak)) then !peak is at left, put varname at right
              call pgptxt(xmax-0.05*dx,ymax*(1.-0.1*fontsize1d),0.,1.,trim(pgvarnss(p)))
           else
              call pgptxt(xmin+0.05*dx,ymax*(1.-0.1*fontsize1d),0.,0.,trim(pgvarnss(p)))
           end if
        end if
        
        
        !Write the deltas of the two pdfs
        if(nchains.eq.2..and.(plrange.eq.4.or.plrange.eq.5.or.plrange.eq.6)) then
           write(str,'(A8)')'\(2030)'
           if(p.eq.2.or.p.eq.5.or.p.eq.14.or.p.eq.15) then
              write(str1,'(A8,F6.2,A1)')'\(2030):',ranges(1,c0,p,5)*100.,'%'
              write(str2,'(A8,F6.2,A1)')'\(2030):',ranges(2,c0,p,5)*100.,'%'
           else
              write(str1,'(A8,F7.3)')'\(2030):',ranges(1,c0,p,5)
              write(str2,'(A8,F7.3)')'\(2030):',ranges(2,c0,p,5)
           end if
        end if
        call pgsch(sch*1.1)
        if(prvalues.eq.1.and.p.ne.1) then  !If not plotting log(L)
           if(nchains.eq.2) then
              call pgsci(colours(mod(0,ncolours)+1))
              call pgmtxt('T',0.5,0.25,0.5,trim(str1))
              call pgsci(colours(mod(1,ncolours)+1))
              call pgmtxt('T',0.5,0.75,0.5,trim(str2))
           else
              if(quality.eq.2.or.quality.eq.3) call pgsci(2)
              !call pgptxt(ranges(ic,c0,p,3),ymax,0.,0.5,trim(str)) !Align with centre of 90%-probability range
              call pgptxt((xmin+xmax)/2.,ymax,0.,0.5,trim(str)) !Centre
              call pgsci(2)
              if(plrange.eq.1.or.plrange.eq.3.or.plrange.eq.4.or.plrange.eq.6) call pgline(2,(/ranges(ic,c0,p,1),ranges(ic,c0,p,2)/),(/0.99*ymax,0.99*ymax/))  !Plot line at top over 90%-probability width
              call pgsci(1)
           end if
        end if
        !else  !If parameter was fixed, plot variable name in empty panel
        !   call pgsch(sch*1.2)
        !   call pgptxt(xmin+0.05*dx,ymax*0.9,0.,0.,trim(pgvarnss(p)))
        !end if !if(fixedpar(p).eq.0)
        
        call pgsci(1)
        call pgsch(sch)
        call pgbox('BNTS',0.0,0,'',0.0,0)
     end if !if(plot.eq.1) 
  end do !p
  
  if(savepdf.eq.1) close(30)
  
  if(plot.eq.1) then
     call pgsubp(1,1)
     call pgsvp(0.,1.,0.,1.)
     call pgswin(-1.,1.,-1.,1.)
     
     if(quality.eq.0) then
        !Remove also the pgsvp at the beginning of the plot
        string=' '
        write(string,'(A,I7,A,I4)')trim(string)//'n:',totpts,', nbin:',nbin1d
        call pgsch(sch*0.5)
        call pgmtxt('T',-0.7,0.5,0.5,trim(outputname)//'  '//trim(string))  !Print title
        call pgsch(sch)
     end if
     
     call pgend
     
     if(file.ge.2) then
        if(file.eq.3) then
           i = system('eps2pdf pdfs.eps -o '//trim(outputdir)//'/'//trim(outputname)//'__pdfs.pdf >& /dev/null')
           if(i.ne.0) write(*,'(A,I6)')'  Error converting plot',i
        end if
        i = system('mv -f pdfs.eps '//trim(outputdir)//'/'//trim(outputname)//'__pdfs.eps')
     else if(file.eq.1) then
        i = system('convert -resize '//trim(bmpxpix)//' -depth 8 -unsharp '//trim(unsharppdf1d)//' pdfs.ppm '//trim(outputdir)//'/'//trim(outputname)//'__pdfs.png')
        if(i.ne.0) write(*,'(A,I6)')'  Error converting plot',i
        i = system('rm -f pdfs.ppm')
     end if
  end if !if(plot.eq.1)
  !write(*,*)''   
  
  
  !Deallocate memory:
  deallocate(xbin,ybin,xbin1,ybin1,ybin2,ysum,yconv,ycum)
  
  
end subroutine pdfs1d
!************************************************************************************************************************************









!************************************************************************************************************************************
subroutine bindata1d(n,x,norm,nbin,xmin1,xmax1,xbin,ybin)  !Count the number of points in each bin (1D)
  ! x - input: data, n points
  ! norm - input: normalise (1) or not (0)
  ! nbin - input: number of bins
  ! xmin, xmax - in/output: set xmin=xmax to auto-determine
  ! xbin, ybin - output: binned data (x, y).  The x values are the left side of the bin!
  
  implicit none
  integer :: i,k,n,nbin,norm
  real :: x(n),xbin(nbin+1),ybin(nbin+1),xmin,xmax,dx,xmin1,xmax1
  
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
!************************************************************************************************************************************


!************************************************************************************************************************************
subroutine bindata1da(n,x,y,norm,nbin,xmin1,xmax1,xbin,ybin)  !Measure the amount of likelihood in each bin (1D)
  ! x - input: data, n points
  ! y - input: "weight" (likelihood), n points
  ! norm - input: normalise (1) or not (0)
  ! nbin - input: number of bins
  ! xmin, xmax - in/output: set xmin=xmax to auto-determine
  ! xbin, ybin - output: binned data (x, y).  The x values are the left side of the bin!

  implicit none
  integer :: i,k,n,nbin,norm
  real :: x(n),y(n),xbin(nbin+1),ybin(nbin+1),xmin,xmax,dx,ybintot,xmin1,xmax1,ymin
  
  xmin = xmin1
  xmax = xmax1
  ymin = minval(y)
  !print*,n,nbin,xmin1,xmax1
  !print*,minval(y),maxval(y)

  if(abs(xmin-xmax)/(xmax+1.e-30).lt.1.e-20) then !Autodetermine
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

  if(abs(xmin1-xmax1)/(xmax1+1.e-30).lt.1.e-20) then
     xmin1 = xmin
     xmax1 = xmax
  end if

end subroutine bindata1da
!************************************************************************************************************************************



!************************************************************************************************************************************
subroutine verthist(n,x,y,style)  !Plot a 1D vertical histogram.  x is the left of the bin!
  implicit none
  integer :: j,n,style
  real :: x(n+1),y(n+1),dx,c
  
  !Make columns overlap a little
  dx = x(2)-x(1)
  c = 0.001*dx
  
  x(n+1) = x(n) + (x(n)-x(n-1))
  y(n+1) = 0.
  
  !Don't close line at left and right (don't drop to 0)
  if(style.eq.0) then
     do j=1,n-1
        call pgline(2,x(j:j+1),(/y(j),y(j)/))
        call pgline(2,(/x(j+1),x(j+1)/),y(j:j+1))
     end do
     call pgline(2,x(n:n+1),(/y(n),y(n)/))
  end if
  
  !Close line at left and right (drop to 0)
  if(style.eq.1) then
     call pgline(2,(/x(1),x(1)/),(/0.,y(1)/))
     do j=1,n
        call pgline(2,x(j:j+1),(/y(j),y(j)/))
        call pgline(2,(/x(j+1),x(j+1)/),y(j:j+1))
     end do
  end if
  
  !Fill the histogram
  if(style.eq.2) then
     do j=1,n
        call pgpoly(4,(/x(j)-c,x(j)-c,x(j+1)+c,x(j+1)+c/),(/0.,y(j),y(j),0./))
     end do
  end if
  
end subroutine verthist
!************************************************************************************************************************************

!************************************************************************************************************************************
subroutine horzhist(n,x,y)  !Plot a 1D horizontal histogram
  implicit none
  integer :: j,n
  real :: x(n),y(n)

  call pgline(2,(/x(1),x(1)/),(/0.,y(1)/))
  do j=1,n-2
     call pgline(2,x(j:j+1),(/y(j),y(j)/))
     call pgline(2,(/x(j+1),x(j+1)/),y(j:j+1))
  end do
  call pgline(2,x(n-1:n),(/y(n-1),y(n-1)/))
end subroutine horzhist
!************************************************************************************************************************************


!************************************************************************************************************************************
subroutine smoothpdf1d(ybin,nbin,smooth)
  implicit none
  integer :: nbin,smooth,i,i0,i1,i00
  real :: ybin(nbin),ybin1(nbin),coefs(100),coefs1(100)
  
  ybin1 = ybin
  i0 = min(max(smooth,1),floor(real(nbin-1)/2.))
  
  !Do all points, except the first and last i0
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
  !Do the first and last i0=i00 points
  if(i00.ge.2) then
  !if(i00.ge.2.and.1.eq.2) then
     do i0 = i00,2,-1
        i=1+i0  !Point at beginning
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
        
        i=nbin-i0  !Point at end
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
!************************************************************************************************************************************
