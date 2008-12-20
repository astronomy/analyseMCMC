!Plot chains (log, parameters, jumps, etc.) for analysemcmc

subroutine chains(exitcode)
  use constants
  use analysemcmc_settings
  use general_data
  use plot_data
  use chain_data
  implicit none
  integer :: i,j,n0,pgopen,imin,ci,lw,symbol,io,ic,p,system,exitcode
  real :: rev360,rev24
  real :: dx,dy,xmin,xmax,ymin,ymax,sch,plx,ply
  character :: string*99
  
  exitcode = 0
  if(combinechainplots.eq.1.and.(pllogl.eq.1.or.plchain.eq.1.or.plsigacc.ge.1)) then
     io = pgopen('chaininfo.eps'//trim(psclr))
     call pginitl(colour,file,whitebg)
  end if



  !***********************************************************************************************************************************      
  !Plot likelihood chain
  if(pllogl.eq.1) then
     !if(prprogress.ge.1.and.update.eq.0) write(*,'(A)')' Plotting chain likelihood...'
     if(prprogress.ge.1.and.update.eq.0) write(*,'(A,$)')' chain likelihood, '
     if(file.eq.0) then
        io = pgopen('12/xs')
        call pgsch(1.5)
     end if
     if(file.ge.1.and.combinechainplots.eq.0) then
        if(file.eq.1) io = pgopen('logL.ppm/ppm')
        if(file.ge.2) io = pgopen('logL.eps'//trim(psclr))
        call pgsch(1.2)
     end if
     if(io.le.0) then
        write(*,'(A,I4)')'Cannot open PGPlot device.  Quitting the programme',io
        exitcode = 1
        return
     end if
     if(file.eq.0) call pgpap(scrsz,scrrat)
     if(file.eq.0) call pgsch(1.5)
     if(file.eq.1) call pgpap(bmpsz,bmprat)
     if(file.eq.1) call pgsch(1.5)
     !call pgscr(3,0.,0.5,0.)
     call pginitl(colour,file,whitebg)
     !call pgsubp(1,2)

     if(quality.eq.0) call pgsvp(0.08,0.95,0.06,0.94) !To make room for title

     ic = 1
     p=1
     !call pgpage
     xmax = -1.e30
     ymin =  1.e30
     ymax = -1.e30
     do ic=1,nchains0
        xmin = 0.
        xmax = max(xmax,maxval(is(ic,1:ntot(ic))))
        imin = 10                                              !Take into account burn-in
        if(scloglpl.eq.1) imin = nburn(ic)                   !Scale without taking into account burnin
        ymin = min(ymin,minval(pldat(ic,p,imin:ntot(ic)))) 
        ymax = max(ymax,maxval(pldat(ic,p,imin:ntot(ic))))
     end do
     ic = 1
     p = 1
     if(ymax.gt.0.) then !This is log(L)-log(Lo) (which we started saving later on), so that nullh=Lo=0
        !ymin = min(ymin,startval(ic,p,1),startval(ic,p,2),0.)
        !ymax = max(ymax,startval(ic,p,1),startval(ic,p,2),0.)
        if(scloglpl.eq.0) then                                 !Take into account 0, true and starting values
           ymin = min(ymin,startval(ic,p,1),startval(ic,p,2),0.)
           ymax = max(ymax,startval(ic,p,1),startval(ic,p,2),0.)
        else                                                     !Take into account true values only
           ymin = min(ymin,startval(ic,p,1))
           ymax = max(ymax,startval(ic,p,1))
        end if
     else !This is log(L), so that Lo = nullh
        ymin = min(ymin,startval(ic,p,1),startval(ic,p,2),nullh)
        ymax = max(ymax,startval(ic,p,1),startval(ic,p,2),nullh)
     end if
     ymin = max(0.,ymin)  !Since we're using null logL := 0 now
     dx = abs(xmax-xmin)*0.01
     dy = abs(ymax-ymin)*0.05

     call pgswin(xmin-dx,xmax+dx,ymin-dy,ymax+dy)
     call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)
     if(abs(startval(1,1,1)-startval(1,1,2))/abs(startval(1,1,1)).gt.1.e-10) then
        call pgsls(4)
        call pgbox('',0.0,0,'G',0.0,0)
        call pgsls(1)
     end if

     do ic=1,nchains0
        !call pgsci(defcolour)
        !if(nchains0.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
        !do i=ic,ntot(ic),chainpli !Start at ic to reduce overplotting
        !   call pgpoint(1,is(ic,i),pldat(ic,p,i),1)
        !end do

        !Give pre- and post-burnin different colour
        ci = defcolour
        if(nchains0.gt.1) ci = colours(mod(ic-1,ncolours)+1)
        call pgscidark(ci,file,whitebg)
        do i=ic,nburn(ic),chainpli !Start at ic to reduce overplotting
           call pgpoint(1,is(ic,i),pldat(ic,p,i),1)
        end do
        call pgsci(ci)
        do i=nburn(ic)+ic,ntot(ic),chainpli !Start at ic to reduce overplotting
           call pgpoint(1,is(ic,i),pldat(ic,p,i),1)
        end do
     end do

     !Plot max likelihood
     if(pllmax.ge.1) then
        ply = pldat(icloglmax,p,iloglmax)
        call pgsci(1)
        call pgpoint(1,is(icloglmax,iloglmax),ply,18)
        call pgsls(5)
        call pgline(2,(/-1.e20,1.e20/),(/ply,ply/))
     end if

     do ic=1,nchains0
        call pgsci(1)
        call pgsls(2)
        call pgline(2,(/-1.e20,1.e20/),(/startval(ic,p,1),startval(ic,p,1)/))
        call pgsci(6)
        if(nchains0.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
        !call pgline(2,(/real(nburn(ic)),real(nburn(ic))/),(/-1.e20,1.e20/))
        if(plburn.ge.1.and.isburn(ic).lt.is(ic,ntot(ic))) call pgline(2,(/isburn(ic),isburn(ic)/),(/-1.e20,1.e20/))
        call pgsci(1)
        call pgsls(4)
        if(nchains0.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
        call pgline(2,(/-1.e20,1.e20/),(/startval(ic,p,2),startval(ic,p,2)/))
        call pgsci(6)
        call pgline(2,(/-1.e20,1.e20/),real((/nullh,nullh/)))
     end do
     call pgsci(1)
     call pgsls(1)
     call pgmtxt('T',0.5,0.1,0.1,trim(pgvarns(p)))

     if(quality.eq.0) then
        !call pgsubp(1,1)
        !call pgsvp(0.,1.,0.,1.)
        !call pgswin(-1.,1.,-1.,1.)

        !call pgsch(sch*0.8)
        call pgmtxt('T',0.5,0.9,0.9,trim(outputname))  !Print title
        !call pgsch(sch)
     end if

     if(combinechainplots.eq.1) call pgpage
     if(combinechainplots.eq.0) then
        call pgend
        if(file.ge.2) then
           if(file.eq.3) then
              i = system('eps2pdf logL.eps  -o '//trim(outputdir)//'/'//trim(outputname)//'__logL.pdf   >& /dev/null')
              if(i.ne.0) write(*,'(A,I6)')'  Error converting plot',i
           end if
           i = system('mv -f logL.eps '//trim(outputdir)//'/'//trim(outputname)//'__logL.eps')
        end if
     end if
     if(file.eq.1) then
        i = system('convert -resize '//trim(bmpxpix)//' -depth 8 -unsharp '//trim(unsharplogl)//' logL.ppm  '//trim(outputdir)//'/'//trim(outputname)//'__logL.png')
        if(i.ne.0) write(*,'(A,I6)')'  Error converting plot',i
        i = system('rm -f logL.ppm')
        !if(i.ne.0) write(*,'(A)')'  Error removing file',i
     end if
  end if !if(pllogl.eq.1) then








  !***********************************************************************************************************************************      
  !Plot chains for each parameter
  if(plchain.eq.1) then
     !if(prprogress.ge.1.and.update.eq.0) write(*,'(A)')' Plotting parameter chains...'
     if(prprogress.ge.1.and.update.eq.0) write(*,'(A,$)')' parameter chains, '
     if(file.eq.0) then
        io = pgopen('13/xs')
        sch = 1.5
        lw = 1
     end if
     if(file.ge.1.and.combinechainplots.eq.0) then
        if(file.eq.1) io = pgopen('chains.ppm/ppm')
        if(file.ge.2) io = pgopen('chains.eps'//trim(psclr))
        lw = 3
        if(nplvar.ge.10) lw = 2
        if(quality.lt.2) lw = max(lw-1,1)  !Draft/Paper
        sch = 1.2
        if(nchains.eq.1.and.nplvar.gt.9) sch = 1.2
        if(quality.eq.0) then !Draft
           sch = sch*1.75
           lw = 2
        end if
        if(quality.eq.1) then !Paper
           if(nplvar.le.12) then
              sch = sch*1.75
              lw = 2
           else
              sch = sch*1.25
              lw = 1
           end if
        end if
        if(quality.eq.2) then !Talk
           if(nplvar.gt.12) then
              sch = sch*1.5
              lw = 1
           end if
           if(nplvar.le.12) then
              sch = sch*2
              lw = 2
           end if
           !if(nplvar.le.6) then
           !   !sch = sch*4
           !   !lw = 4
           !end if
        end if
        if(quality.eq.3) then !Poster
           if(nplvar.eq.12.and.file.ge.2) then
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
     end if
     if(io.le.0) then
        write(*,'(A,I4)')'Cannot open PGPlot device.  Quitting the programme',io
        exitcode = 1
        return
     end if
     if(file.eq.0) call pgpap(scrsz,scrrat)
     if(file.eq.1) call pgpap(bmpsz,bmprat)
     if(file.ge.2) call pgpap(pssz,psrat)
     if(file.ge.2) call pgscf(2)
     !if(file.eq.1) call pgsch(1.5)

     call pgsch(sch)
     call pgslw(lw)


     call pgsubp(panels(1),panels(2))

     ic = 1
     !do p=2,npar0
     do j=1,nplvar
        p = plvars(j)

        call pgpage


        if(j.eq.1) call pginitl(colour,file,whitebg)
        if(file.eq.0.and.scrrat.gt.1.35) call pgsvp(0.08,0.95,0.1,0.95)
        if(file.eq.1.and.bmprat.gt.1.35) call pgsvp(0.08,0.95,0.1,0.95)
        if(file.ge.2.and.psrat.gt.1.35) call pgsvp(0.08,0.95,0.1,0.95)
        if(quality.eq.0) call pgsvp(0.08,0.95,0.06,0.87) !To make room for title
        if(quality.eq.4) call pgsvp(0.13,0.95,0.1,0.95)

        xmin = 0.
        !xmax = real(maxval(ntot(1:nchains0)))
        xmax = -1.e30
        ymin =  1.e30
        ymax = -1.e30
        do ic=1,nchains0
           xmax = max(xmax,maxval(is(ic,1:ntot(ic))))
           imin = 1                                               !Take into account burn-in
           if(scchainspl.eq.1) imin = nburn(ic)                   !Scale without taking into account burnin
           ymin = min(ymin,minval(pldat(ic,p,imin:ntot(ic))))
           ymax = max(ymax,maxval(pldat(ic,p,imin:ntot(ic))))
        end do

        if(p.eq.8) then
           if(ymin.lt.0..or.ymax.gt.24.) then
              ymin = 0.
              ymax = 24.
           end if
        end if
        if(p.eq.10.or.p.eq.12.or.p.eq.13) then
           if(ymin.lt.0..or.ymax.gt.360.) then
              ymin = 0.
              ymax = 360.
           end if
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

        call pgswin(xmin-dx,xmax+dx,ymin-dy,ymax+dy)
        if(quality.eq.1) call pgswin(xmin-dx,xmax+dx,ymin-dy,ymax+dy*2)
        call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)

        !Plot the actual chain values
        call pgsch(1.)
        if(chainsymbol.ne.1) call pgsch(0.7)
        call pgslw(1)
        !write(*,'(15I4)'),nsymbols,symbols(1:nsymbols)
        do ic=1,nchains0
           !call pgsci(mod(ic*2,10))
           !symbol = ic+1
           symbol = chainsymbol
           if(chainsymbol.le.-10) symbol = symbols(mod(ic-1,nsymbols)+1)
           call pgsci(defcolour)
           if(nchains0.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
           if(chainsymbol.eq.0) then !Plot lines rather than symbols
              call pgline(ntot(ic),is(ic,1:ntot(ic)),pldat(ic,p,1:ntot(ic)))
           else
              !do i=ic,ntot(ic),chainpli !Start at ic to reduce overplotting
              !   ply = pldat(ic,p,i)
              !   if(p.eq.8) ply = rev24(ply)
              !   if(p.eq.10.or.p.eq.12.or.p.eq.13) ply = rev360(ply)
              !   !call pgpoint(1,is(ic,i),ply,1) !Plot small dots
              !   call pgpoint(1,is(ic,i),ply,symbol) !Plot symbols
              !end do


              !Give pre- and post-burnin different colour
              ci = defcolour
              if(nchains0.gt.1) ci = colours(mod(ic-1,ncolours)+1)
              call pgscidark(ci,file,whitebg)
              do i=ic,nburn(ic),chainpli !Start at ic to reduce overplotting
                 ply = pldat(ic,p,i)
                 if(p.eq.8) ply = rev24(ply)
                 if(p.eq.10.or.p.eq.12.or.p.eq.13) ply = rev360(ply)
                 call pgpoint(1,is(ic,i),ply,symbol)
              end do
              call pgsci(ci)
              do i=nburn(ic)+ic,ntot(ic),chainpli !Start at ic to reduce overplotting
                 ply = pldat(ic,p,i)
                 if(p.eq.8) ply = rev24(ply)
                 if(p.eq.10.or.p.eq.12.or.p.eq.13) ply = rev360(ply)
                 call pgpoint(1,is(ic,i),ply,symbol)
              end do


           end if
        end do
        call pgsch(sch)
        call pgslw(lw)

        !Plot max likelihood
        if(pllmax.ge.1) then
           ply = pldat(icloglmax,p,iloglmax)
           if(p.eq.8) ply = rev24(ply)
           if(p.eq.10.or.p.eq.12.or.p.eq.13) ply = rev360(ply)
           call pgsci(1)
           call pgpoint(1,is(icloglmax,iloglmax),ply,12)
           call pgsls(5)
           call pgline(2,(/-1.e20,1.e20/),(/ply,ply/))
        end if

        !Plot burn-in, true and starting values
        do ic=1,nchains0
           call pgsls(2)
           call pgsci(6)

           !Plot burn-in phase
           if(nchains0.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
           !if(plburn.ge.1) call pgline(2,(/isburn(ic),isburn(ic)/),(/-1.e20,1.e20/))
           if(plburn.ge.1.and.isburn(ic).lt.is(ic,ntot(ic))) call pgline(2,(/isburn(ic),isburn(ic)/),(/-1.e20,1.e20/))
           call pgsci(1)


           !Plot true values in chains
           if(pltrue.eq.1) then
              if(mergechains.ne.1.or.ic.eq.1) then !The units of the true values haven't changed (e.g. from rad to deg) for ic>1 (but they have for the starting values, why?)
                 !if(ic.eq.1) then !The units of the true values haven't changed (e.g. from rad to deg) for ic>1 (but they have for the starting values, why?)
                 plx = startval(ic,p,1) !True value
                 plx = max(min(1.e30,startval(ic,p,1)),1.e-30)
                 if(p.eq.8) plx = rev24(plx)
                 if(p.eq.10.or.p.eq.12.or.p.eq.13) plx = rev360(plx)
                 call pgline(2,(/-1.e20,1.e20/),(/plx,plx/))
                 if(p.eq.8) then
                    call pgline(2,(/-1.e20,1.e20/),(/plx-24.,plx-24./))
                    call pgline(2,(/-1.e20,1.e20/),(/plx+24.,plx+24./))
                 end if
                 if(p.eq.10.or.p.eq.12.or.p.eq.13) then
                    call pgline(2,(/-1.e20,1.e20/),(/plx-360.,plx-360./))
                    call pgline(2,(/-1.e20,1.e20/),(/plx+360.,plx+360./))
                    !print*,p,plx,plx-360.,plx+360.
                 end if
              end if
           end if


           !Plot starting values in chains
           if(plstart.eq.1.and.abs((startval(ic,p,1)-startval(ic,p,2))/startval(ic,p,1)) .gt. 1.e-10) then
              call pgsls(4)
              if(nchains0.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
              plx = startval(ic,p,2) !Initial value
              if(p.eq.8) plx = rev24(plx)
              if(p.eq.10.or.p.eq.12.or.p.eq.13) plx = rev360(plx)
              call pgline(2,(/-1.e20,1.e20/),(/plx,plx/))
              if(p.eq.8) then
                 call pgline(2,(/-1.e20,1.e20/),(/plx-24.,plx-24./))
                 call pgline(2,(/-1.e20,1.e20/),(/plx+24.,plx+24./))
              end if
              if(p.eq.10.or.p.eq.12.or.p.eq.13) then
                 call pgline(2,(/-1.e20,1.e20/),(/plx-360.,plx-360./))
                 call pgline(2,(/-1.e20,1.e20/),(/plx+360.,plx+360./))
              end if
           end if
        end do !ic=1,nchains0

        call pgsci(1)
        call pgsls(1)
        write(string,'(F6.3)')rhat(p)
        !call pgmtxt('T',1.,0.,0.,'Chain: '//trim(pgvarns(p)))
        call pgmtxt('T',-1.,0.,0.,' '//trim(pgvarnss(p)))
        if(nchains0.gt.1.and.prconv.ge.1) call pgmtxt('T',1.,1.,1.,'Conv: '//trim(string))
     end do !do j=1,nplvar

     if(quality.eq.0) then
        call pgsubp(1,1)
        call pgsvp(0.,1.,0.,1.)
        call pgswin(-1.,1.,-1.,1.)

        call pgsch(sch*0.8)
        call pgmtxt('T',-0.7,0.5,0.5,trim(outputname))  !Print title
        call pgsch(sch)
     end if

     if(combinechainplots.eq.1) call pgpage
     if(combinechainplots.eq.0) then
        call pgend
        if(file.ge.2) then
           if(file.eq.3) then
              i = system('eps2pdf chains.eps  -o '//trim(outputdir)//'/'//trim(outputname)//'__chains.pdf   >& /dev/null')
              if(i.ne.0) write(*,'(A,I6)')'  Error converting plot',i
           end if
           i = system('mv -f chains.eps '//trim(outputdir)//'/'//trim(outputname)//'__chains.eps')
        end if
     end if
     if(file.eq.1) then
        i = system('convert -resize '//trim(bmpxpix)//' -depth 8 -unsharp '//trim(unsharpchain)//' chains.ppm  '//trim(outputdir)//'/'//trim(outputname)//'__chains.png')
        if(i.ne.0) write(*,'(A,I6)')'  Error converting plot',i
        i = system('rm -f chains.ppm')
     end if
  end if !if(plchain.eq.1)








  !***********************************************************************************************************************************      
  !Plot L vs parameter value
  if(plparl.eq.1) then
     !if(prprogress.ge.1.and.update.eq.0) write(*,'(A)')' Plotting parameter-L plot...'
     if(prprogress.ge.1.and.update.eq.0) write(*,'(A,$)')' parameter-L, '
     if(file.eq.0) then
        io = pgopen('22/xs')
        sch = 1.5
        lw = 1
     end if
     if(file.ge.1.and.combinechainplots.eq.0) then
        if(file.eq.1) io = pgopen('parlogl.ppm/ppm')
        if(file.ge.2) io = pgopen('parlogl.eps'//trim(psclr))
        lw = 3
        if(nplvar.ge.10) lw = 2
        if(quality.lt.2) lw = max(lw-1,1)  !Draft/Paper
        sch = 1.2
        if(nchains.eq.1.and.nplvar.gt.9) sch = 1.2
        if(quality.eq.0) then !Draft
           sch = sch*1.75
           lw = 2
        end if
        if(quality.eq.1) then !Paper
           if(nplvar.le.12) then
              sch = sch*1.75
              lw = 2
           else
              sch = sch*1.25
              lw = 1
           end if
        end if
        if(quality.eq.2) then !Talk
           if(nplvar.gt.12) then
              sch = sch*1.5
              lw = 1
           end if
           if(nplvar.le.12) then
              sch = sch*2
              lw = 2
           end if
           !if(nplvar.le.6) then
           !   !sch = sch*4
           !   !lw = 4
           !end if
        end if
        if(quality.eq.3) then !Poster
           if(nplvar.eq.12.and.file.ge.2) then
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
        write(*,'(A,I4)')'Cannot open PGPlot device.  Quitting the programme',io
        exitcode = 1
        return
     end if
     if(file.eq.0) call pgpap(scrsz,scrrat)
     if(file.eq.1) call pgpap(bmpsz,bmprat)
     if(file.ge.2) call pgpap(pssz,psrat)
     if(file.ge.2) call pgscf(2)
     !if(file.eq.1) call pgsch(1.5)

     call pgsch(sch)
     call pgslw(lw)

     if(file.eq.0.and.scrrat.gt.1.35) call pgsvp(0.08,0.95,0.1,0.95)
     if(file.eq.1.and.bmprat.gt.1.35) call pgsvp(0.08,0.95,0.1,0.95)
     if(file.ge.2.and.psrat.gt.1.35) call pgsvp(0.08,0.95,0.1,0.95)
     if(quality.eq.0) call pgsvp(0.08,0.95,0.06,0.87) !To make room for title

     call pgsubp(panels(1),panels(2))

     ic = 1
     !do p=2,npar0
     do j=1,nplvar
        p = plvars(j)
        call pgpage
        if(j.eq.1) call pginitl(colour,file,whitebg)
        xmin = 1.e30
        xmax = -1.e30
        ymin =  1.e30
        ymax = -1.e30
        do ic=1,nchains0
           !xmin = min(xmin,minval(is(ic,nburn(ic):ntot(ic))))
           !xmax = max(xmax,maxval(is(ic,nburn(ic):ntot(ic))))
           xmin = min(xmin,minval(pldat(ic,p,nburn(ic):ntot(ic))))
           xmax = max(xmax,maxval(pldat(ic,p,nburn(ic):ntot(ic))))
           ymin = min(ymin,minval(pldat(ic,1,nburn(ic):ntot(ic))))
           ymax = max(ymax,maxval(pldat(ic,1,nburn(ic):ntot(ic))))
        end do

        if(p.eq.8) then
           !ymin = max(ymin,0.)
           !ymax = min(ymax,24.)
           !ymin = max(rev24(ymin),0.)
           !ymax = min(rev24(ymax),24.)
           if(ymin.lt.0..or.ymax.gt.24.) then
              ymin = 0.
              ymax = 24.
           end if
        end if
        if(p.eq.10.or.p.eq.12.or.p.eq.13) then
           !ymin = max(ymin,0.)
           !ymax = min(ymax,360.)
           !write(*,'(I5,2F10.5,$)')p,ymin,ymax
           !ymin = max(rev360(ymin),0.)
           !ymax = min(rev360(ymax),360.)
           !write(*,'(2F10.5)')ymin,ymax
           if(ymin.lt.0..or.ymax.gt.360.) then
              ymin = 0.
              ymax = 360.
           end if
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
        if(chainsymbol.ne.1) call pgsch(0.7)
        call pgslw(1)
        !write(*,'(15I4)'),nsymbols,symbols(1:nsymbols)
        do ic=1,nchains0
           !call pgsci(mod(ic*2,10))
           !symbol = ic+1
           symbol = chainsymbol
           if(chainsymbol.le.-10) symbol = symbols(mod(ic-1,nsymbols)+1)
           call pgsci(defcolour)
           if(nchains0.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
           do i=nburn(ic),ntot(ic),chainpli
              plx = pldat(ic,p,i)
              ply = pldat(ic,1,i)
              if(p.eq.8) plx = rev24(plx)
              if(p.eq.10.or.p.eq.12.or.p.eq.13) plx = rev360(plx)
              !call pgpoint(1,is(ic,i),plx,1) !Plot small dots
              call pgpoint(1,plx,exp(ply-ymin),symbol) !Plot symbols
              !print*,i,plx,ply,exp(ply-ymin)
           end do
        end do
        call pgsch(sch)
        call pgslw(lw)

        !Plot max likelihood
        if(pllmax.ge.1) then
           plx = pldat(icloglmax,p,iloglmax)
           if(p.eq.8) plx = rev24(plx)
           if(p.eq.10.or.p.eq.12.or.p.eq.13) plx = rev360(plx)
           ply = exp(pldat(icloglmax,1,iloglmax)-ymin)
           call pgsci(1)
           call pgpoint(1,plx,ply,12)
           call pgsls(5)
           call pgline(2,(/plx,plx/),(/-1.e20,1.e20/))
        end if


        !Plot true values
        do ic=1,nchains0
           call pgsls(2)
           call pgsci(1)

           !Plot true values
           if(pltrue.eq.1) then
              if(mergechains.ne.1.or.ic.eq.1) then !The units of the true values haven't changed (e.g. from rad to deg) for ic>1 (but they have for the starting values, why?)
                 !if(ic.eq.1) then !The units of the true values haven't changed (e.g. from rad to deg) for ic>1 (but they have for the starting values, why?)
                 plx = startval(ic,p,1) !True value
                 plx = max(min(1.e30,startval(ic,p,1)),1.e-30)
                 if(p.eq.8) plx = rev24(plx)
                 if(p.eq.10.or.p.eq.12.or.p.eq.13) plx = rev360(plx)
                 call pgline(2,(/plx,plx/),(/-1.e20,1.e20/))
                 if(p.eq.8) then
                    call pgline(2,(/plx-24.,plx-24./),(/-1.e20,1.e20/))
                    call pgline(2,(/plx+24.,plx+24./),(/-1.e20,1.e20/))
                 end if
                 if(p.eq.10.or.p.eq.12.or.p.eq.13) then
                    call pgline(2,(/plx-360.,plx-360./),(/-1.e20,1.e20/))
                    call pgline(2,(/plx+360.,plx+360./),(/-1.e20,1.e20/))
                 end if
              end if
           end if
        end do !ic=1,nchains0

        call pgsci(1)
        call pgsls(1)
     end do !do j=1,nplvar

     if(quality.eq.0) then
        call pgsubp(1,1)
        call pgsvp(0.,1.,0.,1.)
        call pgswin(-1.,1.,-1.,1.)

        call pgsch(sch*0.8)
        call pgmtxt('T',-0.7,0.5,0.5,trim(outputname))  !Print title
        call pgsch(sch)
     end if

     if(combinechainplots.eq.1) call pgpage
     if(combinechainplots.eq.0) then
        call pgend
        if(file.ge.2) then
           if(file.eq.3) then
              i = system('eps2pdf parlogl.eps  -o '//trim(outputdir)//'/'//trim(outputname)//'__parlogl.pdf   >& /dev/null')
              if(i.ne.0) write(*,'(A,I6)')'  Error converting plot',i
           end if
           i = system('mv -f parlogl.eps '//trim(outputdir)//'/'//trim(outputname)//'__parlogl.eps')
        end if
     end if
     if(file.eq.1) then
        i = system('convert -resize '//trim(bmpxpix)//' -depth 8 -unsharp '//trim(unsharpchain)//' parlogl.ppm  '//trim(outputdir)//'/'//trim(outputname)//'__parlogl.png')
        if(i.ne.0) write(*,'(A,I6)')'  Error converting plot',i
        i = system('rm -f parlogl.ppm')
     end if
  end if !if(plparl.eq.1)








  !***********************************************************************************************************************************            
  !Plot jump sizes
  if(pljump.ge.1) then
     !if(prprogress.ge.1.and.update.eq.0) write(*,'(A)')' Plotting jump sizes...'
     if(prprogress.ge.1.and.update.eq.0) write(*,'(A,$)')' jump sizes, '
     if(file.eq.0) then
        io = pgopen('18/xs')
        sch=1.5
     end if
     if(file.ge.1.and.combinechainplots.eq.0) then
        if(file.eq.1) io = pgopen('jumps.ppm/ppm')
        if(file.ge.2) io = pgopen('jumps.eps'//trim(psclr))
        sch=1.2
     end if
     if(io.le.0) then
        write(*,'(A,I4)')'Cannot open PGPlot device.  Quitting the programme',io
        exitcode = 1
        return
     end if
     if(file.eq.0) call pgpap(scrsz,scrrat)
     if(file.eq.0) sch=1.5
     if(file.eq.1) call pgpap(bmpsz,bmprat)
     if(file.eq.1) sch=1.5
     if(quality.eq.0) then !Draft
        !sch = sch*1.75
        sch = sch*1.4
        lw = 2
        call pgsvp(0.1,0.95,0.06,0.87) !To make room for title
     end if

     call pgsch(sch)

     call pgsubp(panels(1),panels(2))

     ic = 1
     !do p=2,npar0
     do j=1,nplvar
        p = plvars(j)
        call pgpage
        if(j.eq.1) call pginitl(colour,file,whitebg)
        xmax = -1.e30
        ymin =  1.e30
        ymax = -1.e30
        do ic=1,nchains0
           xmin = 0.
           xmax = max(xmax,maxval(is(ic,1:ntot(ic))))
           dx = abs(xmax-xmin)*0.01
           if(pljump.eq.1) then
              ymin = min(ymin,minval(jumps(ic,p,2:ntot(ic))))
              ymax = max(ymax,maxval(jumps(ic,p,2:ntot(ic))))
           end if
           do i=10,ntot(ic)
              if(pljump.eq.2.and.jumps(ic,p,i).gt.1.e-20) then
                 ymin = min(ymin,log10(abs(jumps(ic,p,i))))
                 ymax = max(ymax,log10(abs(jumps(ic,p,i))))
              end if
              ymin = -6.
              ymax = 1.
           end do
           !print*,p-1,ymin,ymax,dy
           dy = abs(ymax-ymin)*0.05
           if(dy.lt.1.e-10) dy = ymin*0.1
        end do

        call pgswin(xmin-dx,xmax+dx,ymin-dy,ymax+dy)
        if(pljump.eq.1) call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0) !lin
        if(pljump.eq.2) call pgbox('BCNTS',0.0,0,'BCLNTS',0.0,0) !log

        do ic=1,nchains0
           !call pgsci(mod(ic*2,10))
           call pgsci(defcolour)
           if(nchains0.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
           !if(thin.le.1) then
           !   if(pljump.eq.1) then
           !      do i=1,ntot(ic),chainpli
           !         call pgpoint(1,is(ic,i),jumps(ic,p,i),1)
           !      end do
           !   else
           !      do i=1,ntot(ic),chainpli
           !         call pgpoint(1,is(ic,i),log10(abs(jumps(ic,p,i))+1.e-30),1)
           !      end do
           !   end if
           !else
           !   if(pljump.eq.1) then
           !      call pgpoint(ntot(ic),is(ic,1:ntot(ic)),jumps(ic,p,1:ntot(ic)),1)
           !   else
           !      call pgpoint(ntot(ic),is(ic,1:ntot(ic)),log10(abs(jumps(ic,p,1:ntot(ic)))+1.e-30),1)
           !   end if
           !end if
           if(pljump.eq.1) then
              !do i=1,ntot(ic),chainpli
              do i=ic,ntot(ic),chainpli !Start at ic to reduce overplotting
                 call pgpoint(1,is(ic,i),jumps(ic,p,i),1)
              end do
           else
              !do i=1,ntot(ic),chainpli
              do i=ic,ntot(ic),chainpli !Start at ic to reduce overplotting
                 call pgpoint(1,is(ic,i),log10(abs(jumps(ic,p,i))+1.e-30),1)
              end do
           end if
        end do

        call pgsls(2)
        call pgsci(6)
        do ic=1,nchains0
           !call pgline(2,(/real(nburn(ic)),real(nburn(ic))/),(/-1.e20,1.e20/))
           if(nchains0.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
           !if(plburn.ge.1) call pgline(2,(/isburn(ic),isburn(ic)/),(/-1.e20,1.e20/))
           if(plburn.ge.1.and.isburn(ic).lt.is(ic,ntot(ic))) call pgline(2,(/isburn(ic),isburn(ic)/),(/-1.e20,1.e20/))
        end do
        call pgsci(1)
        call pgsls(1)
        !call pgmtxt('T',1.,0.5,0.5,'Jumps: '//trim(pgorigvarns(p)))
        call pgmtxt('T',-1.2,0.05,0.0,trim(pgorigvarns(p)))
     end do

     if(quality.eq.0) then
        call pgsubp(1,1)
        call pgsvp(0.,1.,0.,1.)
        call pgswin(-1.,1.,-1.,1.)

        call pgsch(sch*0.8)
        call pgmtxt('T',-0.7,0.5,0.5,trim(outputname))  !Print title
        call pgsch(sch)
     end if

     if(combinechainplots.eq.1) call pgpage
     if(combinechainplots.eq.0) then
        call pgend
        if(file.ge.2) then
           if(file.eq.3) then
              i = system('eps2pdf jumps.eps  -o '//trim(outputdir)//'/'//trim(outputname)//'__jumps.pdf   >& /dev/null')
              if(i.ne.0) write(*,'(A,I6)')'  Error converting plot',i
           end if
           i = system('mv -f jumps.eps '//trim(outputdir)//'/'//trim(outputname)//'__jumps.eps')
        end if
     end if
     if(file.eq.1) then
        i = system('convert -resize '//trim(bmpxpix)//' -depth 8 -unsharp '//trim(unsharpchain)//' jumps.ppm  '//trim(outputdir)//'/'//trim(outputname)//'__jumps.png')
        if(i.ne.0) write(*,'(A,I6)')'  Error converting plot',i
        i = system('rm -f jumps.ppm')
     end if
  end if !if(pljump.ge.1)









  !***********************************************************************************************************************************            
  !Plot sigma values ('jump proposal width')
  if(plsigacc.ge.1) then
     !if(prprogress.ge.1.and.update.eq.0) write(*,'(A)')' Plotting sigma...'
     if(prprogress.ge.1.and.update.eq.0) write(*,'(A,$)')' sigma, '
     if(file.eq.0) then
        io = pgopen('16/xs')
        call pgsch(1.5)
     end if
     if(file.ge.1.and.combinechainplots.eq.0) then
        if(file.eq.1) io = pgopen('sigs.ppm/ppm')
        if(file.ge.2) io = pgopen('sigs.eps'//trim(psclr))
        call pgsch(1.2)
     end if
     if(io.le.0) then
        write(*,'(A,I4)')'Cannot open PGPlot device.  Quitting the programme',io
        exitcode = 1
        return
     end if
     if(file.eq.0) call pgpap(scrsz,scrrat)
     if(file.eq.0) call pgsch(1.5)
     if(file.eq.1) call pgpap(bmpsz,bmprat)
     if(file.eq.1) call pgsch(1.5)

     if(quality.eq.0) call pgsvp(0.08,0.95,0.06,0.87) !To make room for title

     call pgsubp(panels(1),panels(2))

     ic = 1
     !do p=2,npar0
     do j=1,nplvar
        p = plvars(j)
        call pgpage
        if(j.eq.1) call pginitl(colour,file,whitebg)
        xmax = -1.e30
        ymin =  1.e30
        ymax = -1.e30
        do ic=1,nchains0
           xmin = 0.
           xmax = max(xmax,maxval(is(ic,1:ntot(ic))))
           dx = abs(xmax-xmin)*0.01
           if(plsigacc.eq.1) then
              ymin = min(ymin,minval(sig(p,ic,10:ntot(ic))))
              ymax = max(ymax,maxval(sig(p,ic,10:ntot(ic))))
           end if
           do i=10,ntot(ic)
              if(plsigacc.eq.2.and.sig(p,ic,i).gt.1.e-20) then
                 ymin = min(ymin,log10(sig(p,ic,i)))
                 ymax = max(ymax,log10(sig(p,ic,i)))
              end if
           end do
           !print*,p-1,ymin,ymax,dy
           dy = abs(ymax-ymin)*0.05
           if(dy.lt.1.e-10) dy = ymin*0.1
        end do

        call pgswin(xmin-dx,xmax+dx,ymin-dy,ymax+dy)
        if(plsigacc.eq.1) call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0) !lin
        if(plsigacc.eq.2) call pgbox('BCNTS',0.0,0,'BCLNTS',0.0,0) !log

        do ic=1,nchains0
           !call pgsci(mod(ic*2,10))
           call pgsci(defcolour)
           if(nchains0.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
           !if(thin.le.1) then
           !   if(plsigacc.eq.1) then
           !      do i=1,ntot(ic),chainpli
           !         call pgpoint(1,is(ic,i),sig(p,ic,i),1)
           !      end do
           !   else
           !      do i=1,ntot(ic),chainpli
           !         call pgpoint(1,is(ic,i),log10(sig(p,ic,i)+1.e-30),1)
           !      end do
           !   end if
           !else
           !   if(plsigacc.eq.1) then
           !      call pgpoint(ntot(ic),is(ic,1:ntot(ic)),sig(p,ic,1:ntot(ic)),1)
           !   else
           !      call pgpoint(ntot(ic),is(ic,1:ntot(ic)),log10(sig(p,ic,1:ntot(ic))+1.e-30),1)
           !   end if
           !end if
           if(plsigacc.eq.1) then
              !do i=1,ntot(ic),chainpli
              do i=ic,ntot(ic),chainpli !Start at ic to reduce overplotting
                 call pgpoint(1,is(ic,i),sig(p,ic,i),1)
              end do
           else
              !do i=1,ntot(ic),chainpli
              do i=ic,ntot(ic),chainpli !Start at ic to reduce overplotting
                 call pgpoint(1,is(ic,i),log10(sig(p,ic,i)+1.e-30),1)
              end do
           end if
        end do

        call pgsls(2)
        call pgsci(6)
        do ic=1,nchains0
           !call pgline(2,(/real(nburn(ic)),real(nburn(ic))/),(/-1.e20,1.e20/))
           if(nchains0.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
           !if(plburn.ge.1) call pgline(2,(/isburn(ic),isburn(ic)/),(/-1.e20,1.e20/))
           if(plburn.ge.1.and.isburn(ic).lt.is(ic,ntot(ic))) call pgline(2,(/isburn(ic),isburn(ic)/),(/-1.e20,1.e20/))
        end do
        call pgsci(1)
        call pgsls(1)
        call pgmtxt('T',1.,0.5,0.5,'Sigma: '//trim(pgvarns(p)))
        !call pgmtxt('T',1.,0.5,0.5,'log Sigma: '//trim(pgvarns(p)))
     end do

     if(quality.eq.0) then
        call pgsubp(1,1)
        call pgsvp(0.,1.,0.,1.)
        call pgswin(-1.,1.,-1.,1.)

        call pgsch(sch*0.8)
        call pgmtxt('T',-0.7,0.5,0.5,trim(outputname))  !Print title
        call pgsch(sch)
     end if

     if(combinechainplots.eq.1) call pgpage
     if(combinechainplots.eq.0) then
        call pgend
        if(file.ge.2) then
           if(file.eq.3) then
              i = system('eps2pdf sigs.eps  -o '//trim(outputdir)//'/'//trim(outputname)//'__sigs.pdf   >& /dev/null')
              if(i.ne.0) write(*,'(A,I6)')'  Error converting plot',i
           end if
           i = system('mv -f sigs.eps '//trim(outputdir)//'/'//trim(outputname)//'__sigs.eps')
        end if
     end if
     if(file.eq.1) then
        i = system('convert -resize '//trim(bmpxpix)//' -depth 8 -unsharp '//trim(unsharpchain)//' sigs.ppm  '//trim(outputdir)//'/'//trim(outputname)//'__sigs.png')
        if(i.ne.0) write(*,'(A,I6)')'  Error converting plot',i
        i = system('rm -f sigs.ppm')
     end if
  end if !if(plsigacc.ge.1)









  !***********************************************************************************************************************************      
  !Plot acceptance rates
  if(plsigacc.ge.1) then
     !if(prprogress.ge.1.and.update.eq.0) write(*,'(A)')' Plotting acceptance rates...'
     if(prprogress.ge.1.and.update.eq.0) write(*,'(A,$)')' acceptance rates, '
     if(file.eq.0) then
        io = pgopen('17/xs')
        call pgsch(1.5)
     end if
     if(file.ge.1.and.combinechainplots.eq.0) then
        if(file.eq.1) io = pgopen('accs.ppm/ppm')
        if(file.ge.2) io = pgopen('accs.eps'//trim(psclr))
        call pgsch(1.2)
     end if
     if(io.le.0) then
        write(*,'(A,I4)')'Cannot open PGPlot device.  Quitting the programme',io
        exitcode = 1
        return
     end if
     if(file.eq.0) call pgpap(scrsz,scrrat)
     if(file.eq.0) call pgsch(1.5)
     if(file.eq.1) call pgpap(bmpsz,bmprat)
     if(file.eq.1) call pgsch(1.5)

     if(quality.eq.0) call pgsvp(0.08,0.95,0.06,0.87) !To make room for title

     call pgsubp(panels(1),panels(2))

     ic = 1
     !do p=2,npar0
     do j=1,nplvar
        p = plvars(j)
        call pgpage
        if(j.eq.1) call pginitl(colour,file,whitebg)
        xmax = -1.e30
        ymin =  1.e30
        ymax = -1.e30
        do ic=1,nchains0
           xmin = 0.
           !xmax = max(xmax,real(ntot(ic)))
           xmax = max(xmax,maxval(is(ic,1:ntot(ic))))
           dx = abs(xmax-xmin)*0.01
           do i=1,ntot(ic)
              if(acc(p,ic,i).gt.1.e-10 .and. acc(p,ic,i).lt.1.-1.e-10) then
                 n0 = i
                 exit
              end if
           end do
           n0 = n0+10
           ymin = min(ymin,minval(acc(p,ic,n0:ntot(ic))))
           ymax = max(ymax,maxval(acc(p,ic,n0:ntot(ic))))
           dy = abs(ymax-ymin)*0.05
        end do

        call pgsci(1)
        call pgsls(1)
        call pgswin(xmin-dx,xmax+dx,ymin-dy,ymax+dy)
        call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)


        call pgsci(3)
        call pgsls(2)
        call pgline(2,(/-1.e20,1.e20/),(/0.25,0.25/))
        call pgsci(6)
        do ic=1,nchains0
           !call pgline(2,(/real(nburn(ic)),real(nburn(ic))/),(/-1.e20,1.e20/))
           if(nchains0.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
           !if(plburn.ge.1) call pgline(2,(/isburn(ic),isburn(ic)/),(/-1.e20,1.e20/))
           if(plburn.ge.1.and.isburn(ic).lt.is(ic,ntot(ic))) call pgline(2,(/isburn(ic),isburn(ic)/),(/-1.e20,1.e20/))
        end do
        call pgsci(1)
        call pgsls(1)
        call pgmtxt('T',1.,0.5,0.5,'Acceptance: '//trim(pgvarns(p)))

        do ic=1,nchains0
           !call pgsci(mod(ic*2,10))
           call pgsci(defcolour)
           if(nchains0.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
           !if(thin.le.1) then
           !   do i=1,ntot(ic),chainpli
           !      call pgpoint(1,is(ic,i),acc(p,ic,i),1)
           !   end do
           !else
           !   call pgpoint(ntot(ic),is(ic,1:ntot(ic)),acc(p,ic,1:ntot(ic)),1)
           !end if
           !do i=1,ntot(ic),chainpli
           do i=ic,ntot(ic),chainpli !Start at ic to reduce overplotting
              call pgpoint(1,is(ic,i),acc(p,ic,i),1)
           end do
        end do
     end do

     if(quality.eq.0) then
        call pgsubp(1,1)
        call pgsvp(0.,1.,0.,1.)
        call pgswin(-1.,1.,-1.,1.)

        call pgsch(sch*0.8)
        call pgmtxt('T',-0.7,0.5,0.5,trim(outputname))  !Print title
        call pgsch(sch)
     end if

     !if(combinechainplots.eq.1) call pgpage
     if(combinechainplots.eq.0) then
        call pgend
        if(file.ge.2) then
           if(file.eq.3) then
              i = system('eps2pdf accs.eps  -o '//trim(outputdir)//'/'//trim(outputname)//'__accs.pdf   >& /dev/null')
              if(i.ne.0) write(*,'(A,I6)')'  Error converting plot',i
           end if
           i = system('mv -f accs.eps '//trim(outputdir)//'/'//trim(outputname)//'__accs.eps')
        end if
     end if
     if(file.eq.1) then
        i = system('convert -resize '//trim(bmpxpix)//' -depth 8 -unsharp '//trim(unsharpchain)//' accs.ppm  '//trim(outputdir)//'/'//trim(outputname)//'__accs.png')
        if(i.ne.0) write(*,'(A,I6)')'  Error converting plot',i
        i = system('rm -f accs.ppm')
     end if
  end if !if(plsigacc.ge.1)



  if(file.ge.1.and.combinechainplots.eq.1.and.(pllogl.eq.1.or.plchain.eq.1.or.plsigacc.ge.1)) then
     call pgend
     if(file.eq.3) i = system('eps2pdf chaininfo.eps -o '//trim(outputdir)//'/'//trim(outputname)//'__chaininfo.pdf  >& /dev/null')
     i = system('mv -f chaininfo.eps '//trim(outputdir)//'/'//trim(outputname)//'__chaininfo.eps')
  end if
  !***********************************************************************************************************************************        










  !***********************************************************************************************************************************      
  !Plot autocorrelations for each parameter
  if(placorr.gt.0) then
     !if(prprogress.ge.1.and.update.eq.0) write(*,'(A)')' Plotting autocorrelations...'
     if(prprogress.ge.1.and.update.eq.0) write(*,'(A,$)')' autocorrelations, '
     if(file.eq.0) then
        io = pgopen('19/xs')
        call pgsch(1.5)
     end if
     if(file.ge.1) then
        if(file.eq.1) io = pgopen('acorrs.ppm/ppm')
        if(file.ge.2) io = pgopen('acorrs.eps'//trim(psclr))
        call pgsch(1.2)
     end if
     if(io.le.0) then
        write(*,'(A,I4)')'Cannot open PGPlot device.  Quitting the programme',io
        exitcode = 1
        return
     end if
     if(file.eq.0) call pgpap(scrsz,scrrat)
     if(file.eq.0) call pgsch(1.5)
     if(file.eq.1) call pgpap(bmpsz,bmprat)
     if(file.eq.1) call pgsch(1.5)

     if(quality.eq.0) call pgsvp(0.08,0.95,0.06,0.87) !To make room for title

     call pgsubp(panels(1),panels(2))

     ic = 1
     !do p=2,npar0
     do j=1,nplvar
        p = plvars(j)
        call pgpage
        if(j.eq.1) call pginitl(colour,file,whitebg)
        xmin = 0.
        xmin = minval(acorrs(1,0,0:100))
        xmax = maxval(acorrs(1,0,0:100))
        dx = abs(xmax-xmin)*0.01
        ymin =  1.e30
        ymax = -1.e30
        do ic=1,nchains0
           ymin = min(ymin,minval(acorrs(ic,p,0:100)))
           ymax = max(ymax,maxval(acorrs(ic,p,0:100)))
        end do
        dy = abs(ymax-ymin)*0.05
        !write(*,'(I3,5F10.2)')p,xmin,xmax,ymin,ymax,acorrs(1,0,100)

        call pgswin(xmin-dx,xmax+dx,ymin-dy,ymax+dy)
        call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)

        do ic=1,nchains0
           call pgsci(defcolour)
           if(nchains.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
           do i=1,100!placorr!ntot(ic),chainpli
              call pgpoint(1,acorrs(ic,0,i),acorrs(ic,p,i),1)
           end do
        end do

        call pgsci(1)
        call pgsls(2)
        call pgline(2,(/-1.e20,1.e20/),(/0.,0./))
        call pgsci(1)
        call pgsls(1)
        call pgmtxt('T',1.,0.5,0.5,'Autocorrelation: '//trim(pgvarns(p)))
     end do

     if(quality.eq.0) then
        call pgsubp(1,1)
        call pgsvp(0.,1.,0.,1.)
        call pgswin(-1.,1.,-1.,1.)

        call pgsch(sch*0.8)
        call pgmtxt('T',-0.7,0.5,0.5,trim(outputname))  !Print title
        call pgsch(sch)
     end if

     if(combinechainplots.eq.1) call pgpage
     if(combinechainplots.eq.0) then
        call pgend
        if(file.ge.2) then
           if(file.eq.3) then
              i = system('eps2pdf acorrs.eps  -o '//trim(outputdir)//'/'//trim(outputname)//'__acorrs.pdf   >& /dev/null')
              if(i.ne.0) write(*,'(A,I6)')'  Error converting plot',i
           end if
           i = system('mv -f acorrs.eps '//trim(outputdir)//'/'//trim(outputname)//'__acorrs.eps')
        end if
     end if
     if(file.eq.1) then
        i = system('convert -resize '//trim(bmpxpix)//' -depth 8 -unsharp '//trim(unsharpchain)//' acorrs.ppm  '//trim(outputdir)//'/'//trim(outputname)//'__acorrs.png')
        if(i.ne.0) write(*,'(A,I6)')'  Error converting plot',i
        i = system('rm -f acorrs.ppm')
     end if
  end if !if(placorrs.gt.0)
  !***********************************************************************************************************************************      


end subroutine chains
