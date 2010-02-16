!Plot chains (posterior, parameters, jumps, etc.) for analysemcmc

subroutine chains(exitcode)
  use constants
  use analysemcmc_settings
  use general_data
  use mcmcrun_data
  use plot_data
  use chain_data
  implicit none
  integer :: i,j,pgopen,imin,ci,lw,symbol,io,ic,p,system,exitcode
  real :: rev360,rev24,rev180,compute_median_real
  real :: dx,dy,xmin,xmax,ymin,ymax,sch,plx,ply
  character :: title*99
  
  exitcode = 0
  
  
  
  !***********************************************************************************************************************************      
  !Plot posterior chain
  if(plLogL.eq.1) then
     if(prProgress.ge.1.and.update.eq.0) write(stdOut,'(A,$)')' posterior chain, '
     if(file.eq.0) then
        io = pgopen('12/xs')
        call pgsch(1.5)
     end if
     if(file.ge.1) then
        if(file.eq.1) io = pgopen('posterior.ppm/ppm')
        if(file.ge.2) io = pgopen('posterior.eps'//trim(psclr))
        call pgsch(1.2)
     end if
     if(io.le.0) then
        write(stdErr,'(A,I4)')'   Error:  Cannot open PGPlot device.  Quitting the programme',io
        exitcode = 1
        return
     end if
     if(file.eq.0) call pgpap(scrSz,scrRat)
     if(file.eq.0) call pgsch(1.5)
     if(file.eq.1) call pgpap(bmpsz,bmprat)
     if(file.eq.1) call pgsch(1.5)
     if(file.ge.2) call pgpap(PSsz,PSrat)
     if(file.ge.2) call pgscf(fonttype)
     !call pgscr(3,0.,0.5,0.)
     call pginitl(colour,file,whiteBG)
     !call pgsubp(1,2)
     
     if(quality.eq.0) call pgsvp(0.08,0.95,0.06,0.94) !To make room for title
     
     ic = 1
     p=1
     !call pgpage
     xmax = -1.e30
     ymin =  1.e30
     ymax = -1.e30
     do ic=1,nChains0
        xmin = 0.
        xmax = max(xmax,maxval(is(ic,1:Ntot(ic))))
        imin = 10                                              !Take into account burn-in
        if(scLogLpl.eq.1) imin = Nburn(ic)                   !Scale without taking into account burnin
        ymin = min(ymin,minval(post(ic,imin:Ntot(ic)))) 
        ymax = max(ymax,maxval(post(ic,imin:Ntot(ic))))
     end do
     ic = 1
     p = 1
     if(scLogLpl.eq.0) then                                 !Take into account 0, injection and starting values
        ymin = min(ymin,post(ic,1),post(ic,2),0.)
        ymax = max(ymax,post(ic,1),post(ic,2),0.)
     else                                                     !Take into account injection values only
        ymin = min(ymin,post(ic,1))
        ymax = max(ymax,post(ic,1))
     end if
     
     ymin = max(0.,ymin)
     dx = abs(xmax-xmin)*0.01
     dy = abs(ymax-ymin)*0.05
     
     call pgswin(xmin-dx,xmax+dx,ymin-dy,ymax+dy)
     call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)
     if(abs(startval(1,1,1)-startval(1,1,2))/abs(startval(1,1,1)).gt.1.e-10) then
        call pgsls(4)
        call pgbox('',0.0,0,'G',0.0,0)
        call pgsls(1)
     end if
     
     do ic=1,nChains0
        !Give pre- and post-burnin different colour
        ci = defcolour
        if(nChains0.gt.1) ci = colours(mod(ic-1,ncolours)+1)
        call pgscidark(ci,file,whiteBG)
        do i=ic,Nburn(ic),chainPlI !Start at ic to reduce overplotting
           call pgpoint(1,is(ic,i),post(ic,i),1)
        end do
        call pgsci(ci)
        do i=Nburn(ic)+ic,Ntot(ic),chainPlI !Start at ic to reduce overplotting
           call pgpoint(1,is(ic,i),post(ic,i),1)
        end do
     end do
     
     !Plot max posterior
     if(plLmax.ge.1) then
        ply = post(icloglmax,iloglmax)
        call pgsci(1)
        call pgpoint(1,is(icloglmax,iloglmax),ply,18)
        call pgsls(5)
        call pgline(2,(/-1.e20,1.e20/),(/ply,ply/))
     end if
     
     
     do ic=1,nChains0
        call pgsci(1)
        call pgsls(2)
        if(abs(post(ic,1)).gt.1.e-4) call pgline(2,(/-1.e20,1.e20/),(/post(ic,1),post(ic,1)/))  !Plot injection value, only if injection was done
        call pgsci(colours(mod(ic-1,ncolours)+1))
        if(plBurn.ge.1.and.isburn(ic).lt.is(ic,Ntot(ic))) call pgline(2,(/isburn(ic),isburn(ic)/),(/-1.e20,1.e20/)) !Vertical line at burn-in
        call pgsls(4)
        call pgline(2,(/-1.e20,1.e20/),(/post(ic,2),post(ic,2)/))  !Horizontal dotted line at starting value
     end do
     call pgsci(1)
     call pgsls(1)
     call pgmtxt('T',0.5,0.1,0.1,'log Posterior')
     
     if(quality.eq.0) then
        !call pgsch(sch*0.8)
        call pgmtxt('T',0.5,0.9,0.9,trim(outputname))  !Print title
        !call pgsch(sch)
     end if
     
     call pgend
     if(file.ge.2) then
        if(file.eq.3) then
           i = system('eps2pdf posterior.eps  -o '//trim(outputdir)//'/'//trim(outputname)//'__posterior.pdf   >& /dev/null')
           if(i.ne.0) write(stdErr,'(A,I6)')'  Error converting plot',i
        end if
        i = system('mv -f posterior.eps '//trim(outputdir)//'/'//trim(outputname)//'__posterior.eps')
     end if
     if(file.eq.1) then
        i = system('convert -resize '//trim(bmpxpix)//' -depth 8 -unsharp '//trim(unSharplogl)//' posterior.ppm  '//trim(outputdir)//'/'//trim(outputname)//'__posterior.png')
        if(i.ne.0) write(stdErr,'(A,I6)')'  Error converting plot',i
        i = system('rm -f posterior.ppm')
        !if(i.ne.0) write(stdErr,'(A)')'  Error removing file',i
     end if
  end if !if(plLogL.eq.1) then
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  !***********************************************************************************************************************************      
  !Plot chains for each parameter
  if(plChain.eq.1) then
     if(prProgress.ge.1.and.update.eq.0) write(stdOut,'(A,$)')' parameter chains, '
     if(file.eq.0) then
        io = pgopen('13/xs')
        sch = fontsize1d*1.5
        lw = 1
     end if
     if(file.ge.1) then
        if(file.eq.1) io = pgopen('chains.ppm/ppm')
        if(file.ge.2) io = pgopen('chains.eps'//trim(psclr))
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
        if(quality.eq.4) then !Vivien's thesis
           sch = sch*2.5
           lw = 2
        end if
     end if
     if(io.le.0) then
        write(stdErr,'(A,I4)')'   Error:  Cannot open PGPlot device.  Quitting the programme',io
        exitcode = 1
        return
     end if
     if(file.eq.0) call pgpap(scrSz,scrRat)
     if(file.eq.1) call pgpap(bmpsz,bmprat)
     if(file.ge.2) call pgpap(PSsz,PSrat)
     if(file.ge.2) call pgscf(fonttype)
     !if(file.eq.1) call pgsch(1.5)
     
     call pgsch(sch)
     call pgslw(lw)
     
     
     call pgsubp(panels(1),panels(2))
     
     ic = 1
     do j=1,nPlPar
        p = revID(plPars(j))
        if(p.eq.0) then
           write(stdErr,'(/,A)')'  * Warning:  chains():  parameter '//trim(parNames(plPars(j)))//' is not defined, check plPars() in the input file.  Skipping...'
           cycle
        end if
        
        call pgpage
        
        
        if(j.eq.1) call pginitl(colour,file,whiteBG)
        if(file.eq.0.and.scrRat.gt.1.35) call pgsvp(0.08,0.95,0.1,0.95)
        if(file.eq.1.and.bmprat.gt.1.35) call pgsvp(0.08,0.95,0.1,0.95)
        if(file.ge.2.and.PSrat.gt.1.35) call pgsvp(0.08,0.95,0.1,0.95)
        if(quality.eq.0) call pgsvp(0.08,0.95,0.06,0.87) !To make room for title
        if(quality.eq.4) call pgsvp(0.13,0.95,0.1,0.95)
        
        xmin = 0.
        !xmax = real(maxval(Ntot(1:nChains0)))
        xmax = -1.e30
        ymin =  1.e30
        ymax = -1.e30
        do ic=1,nChains0
           xmax = max(xmax,maxval(is(ic,1:Ntot(ic))))
           imin = 1                                               !Take into account burn-in
           if(scChainsPl.eq.1) imin = Nburn(ic)                   !Scale without taking into account burnin
           ymin = min(ymin,minval(allDat(ic,p,imin:Ntot(ic))))
           ymax = max(ymax,maxval(allDat(ic,p,imin:Ntot(ic))))
        end do
        
        if(changeVar.gt.0) then
           if(parID(p).eq.31) then  !RA
              if(ymin.lt.0..or.ymax.gt.24.) then
                 ymin = 0.
                 ymax = 24.
              end if
           end if
           if(parID(p).eq.41.or.parID(p).eq.54.or.parID(p).eq.73.or.parID(p).eq.83) then  !Phases
              if(ymin.lt.0..or.ymax.gt.360.) then
                 ymin = 0.
                 ymax = 360.
              end if
           end if
           if(parID(p).eq.52) then  !Psi
              if(ymin.lt.0..or.ymax.gt.180.) then
                 ymin = 0.
                 ymax = 180.
              end if
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
           if(chainSymbol.eq.0) then !Plot lines rather than symbols
              call pgline(Ntot(ic),is(ic,1:Ntot(ic)),allDat(ic,p,1:Ntot(ic)))
           else
              !Give pre- and post-burnin different colour
              ci = defcolour
              if(nChains0.gt.1) ci = colours(mod(ic-1,ncolours)+1)
              call pgscidark(ci,file,whiteBG)
              do i=ic,Nburn(ic),chainPlI !Start at ic to reduce overplotting
                 ply = allDat(ic,p,i)
                 if(changeVar.gt.0) then
                    if(parID(p).eq.31) ply = rev24(ply)  !RA
                    if(parID(p).eq.41.or.parID(p).eq.54.or.parID(p).eq.73.or.parID(p).eq.83) ply = rev360(ply)
                    if(parID(p).eq.52) ply = rev180(ply)
                 end if
                 call pgpoint(1,is(ic,i),ply,symbol)
              end do
              call pgsci(ci)
              do i=Nburn(ic)+ic,Ntot(ic),chainPlI !Start at ic to reduce overplotting
                 ply = allDat(ic,p,i)
                 if(changeVar.gt.0) then
                    if(parID(p).eq.31) ply = rev24(ply)  !RA
                    if(parID(p).eq.41.or.parID(p).eq.54.or.parID(p).eq.73.or.parID(p).eq.83) ply = rev360(ply)
                    if(parID(p).eq.52) ply = rev180(ply)
                 end if
                 call pgpoint(1,is(ic,i),ply,symbol)
              end do
              
              
           end if
        end do
        call pgsch(sch)
        call pgslw(lw)
        
        !Plot max posterior
        if(plLmax.ge.1) then
           ply = allDat(icloglmax,p,iloglmax)
           if(changeVar.gt.0) then
              if(parID(p).eq.31) ply = rev24(ply)  !RA
              if(parID(p).eq.41.or.parID(p).eq.54.or.parID(p).eq.73.or.parID(p).eq.83) ply = rev360(ply)
              if(parID(p).eq.52) ply = rev180(ply)
           end if
           call pgsci(1)
           call pgpoint(1,is(icloglmax,iloglmax),ply,12)
           call pgsls(5)
           call pgline(2,(/-1.e20,1.e20/),(/ply,ply/))
        end if
        
        !Plot burn-in, injection and starting values
        do ic=1,nChains0
           call pgsls(2)
           call pgsci(6)
           
           !Plot burn-in phase
           if(nChains0.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
           !if(plBurn.ge.1) call pgline(2,(/isburn(ic),isburn(ic)/),(/-1.e20,1.e20/))
           if(plBurn.ge.1.and.isburn(ic).lt.is(ic,Ntot(ic))) call pgline(2,(/isburn(ic),isburn(ic)/),(/-1.e20,1.e20/))
           call pgsci(1)
           
           
           !Plot injection values in chains
           if(plInject.ge.1) then
              if(mergeChains.ne.1.or.ic.eq.1) then !The units of the injection values haven't changed (e.g. from rad to deg) for ic>1 (but they have for the starting values, why?)
                 !if(ic.eq.1) then !The units of the injection values haven't changed (e.g. from rad to deg) for ic>1 (but they have for the starting values, why?)
                 plx = startval(ic,p,1) !Injection value
                 plx = max(min(1.e30,startval(ic,p,1)),1.e-30)
                 if(changeVar.gt.0) then
                    if(parID(p).eq.31) plx = rev24(plx)  !RA
                    if(parID(p).eq.41.or.parID(p).eq.54.or.parID(p).eq.73.or.parID(p).eq.83) plx = rev360(plx)
                    if(parID(p).eq.52) plx = rev180(plx)
                 end if
                 call pgline(2,(/-1.e20,1.e20/),(/plx,plx/))
                 if(changeVar.gt.0) then
                    if(parID(p).eq.31) then
                       call pgline(2,(/-1.e20,1.e20/),(/plx-24.,plx-24./))
                       call pgline(2,(/-1.e20,1.e20/),(/plx+24.,plx+24./))
                    end if
                    if(parID(p).eq.41.or.parID(p).eq.54.or.parID(p).eq.73.or.parID(p).eq.83) then
                       call pgline(2,(/-1.e20,1.e20/),(/plx-360.,plx-360./))
                       call pgline(2,(/-1.e20,1.e20/),(/plx+360.,plx+360./))
                    end if
                    if(parID(p).eq.52) then
                       call pgline(2,(/-1.e20,1.e20/),(/plx-180.,plx-180./))
                       call pgline(2,(/-1.e20,1.e20/),(/plx+180.,plx+180./))
                    end if
                 end if
              end if
           end if
           
           
           !Plot starting values in chains
           if(plStart.ge.1.and.abs((startval(ic,p,1)-startval(ic,p,2))/startval(ic,p,1)) .gt. 1.e-10) then
              call pgsls(4)
              if(nChains0.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
              plx = startval(ic,p,2) !Initial value
              if(changeVar.gt.0) then
                 if(parID(p).eq.31) plx = rev24(plx)
                 if(parID(p).eq.41.or.parID(p).eq.54.or.parID(p).eq.73.or.parID(p).eq.83) plx = rev360(plx)
                 if(parID(p).eq.52) plx = rev180(plx)
              end if
              call pgline(2,(/-1.e20,1.e20/),(/plx,plx/))
              if(changeVar.gt.0) then
                 if(parID(p).eq.31) then
                    call pgline(2,(/-1.e20,1.e20/),(/plx-24.,plx-24./))
                    call pgline(2,(/-1.e20,1.e20/),(/plx+24.,plx+24./))
                 end if
                 if(parID(p).eq.41.or.parID(p).eq.54.or.parID(p).eq.73.or.parID(p).eq.83) then
                    call pgline(2,(/-1.e20,1.e20/),(/plx-360.,plx-360./))
                    call pgline(2,(/-1.e20,1.e20/),(/plx+360.,plx+360./))
                 end if
                 if(parID(p).eq.52) then
                    call pgline(2,(/-1.e20,1.e20/),(/plx-180.,plx-180./))
                    call pgline(2,(/-1.e20,1.e20/),(/plx+180.,plx+180./))
                 end if
              end if
           end if
        end do !ic=1,nChains0
        
        call pgsci(1)
        call pgsls(1)
        write(title,'(F6.3)')rhat(p)
        if(nChains0.gt.1.and.prConv.ge.1) then
           call pgmtxt('T',1.,1.,1.,'R-hat: '//trim(title))
        else
           call pgmtxt('T',-1.,0.,0.,' '//trim(pgParNs(parID(p))))
        end if
     end do !do j=1,nPlPar
     
     if(quality.eq.0) then
        call pgsubp(1,1)
        call pgsvp(0.,1.,0.,1.)
        call pgswin(-1.,1.,-1.,1.)
        
        call pgsch(sch*0.8)
        call pgmtxt('T',-0.7,0.5,0.5,trim(outputname))  !Print title
        call pgsch(sch)
     end if
     
     call pgend
     if(file.ge.2) then
        if(file.eq.3) then
           i = system('eps2pdf chains.eps  -o '//trim(outputdir)//'/'//trim(outputname)//'__chains.pdf   >& /dev/null')
           if(i.ne.0) write(stdErr,'(A,I6)')'  Error converting plot',i
        end if
        i = system('mv -f chains.eps '//trim(outputdir)//'/'//trim(outputname)//'__chains.eps')
     end if
     if(file.eq.1) then
        i = system('convert -resize '//trim(bmpxpix)//' -depth 8 -unsharp '//trim(unSharpchain)//' chains.ppm  '//trim(outputdir)//'/'//trim(outputname)//'__chains.png')
        if(i.ne.0) write(stdErr,'(A,I6)')'  Error converting plot',i
        i = system('rm -f chains.ppm')
     end if
  end if !if(plChain.eq.1)
  !***********************************************************************************************************************************      
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  !***********************************************************************************************************************************      
  !Plot L vs parameter value
  if(plParL.eq.1) then
     !if(prProgress.ge.1.and.update.eq.0) write(stdOut,'(A)')' Plotting parameter-L plot...'
     if(prProgress.ge.1.and.update.eq.0) write(stdOut,'(A,$)')' parameter-L, '
     if(file.eq.0) then
        io = pgopen('22/xs')
        sch = fontsize1d*1.5
        lw = 1
     end if
     if(file.ge.1) then
        if(file.eq.1) io = pgopen('parlogl.ppm/ppm')
        if(file.ge.2) io = pgopen('parlogl.eps'//trim(psclr))
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
     if(file.eq.0) call pgpap(scrSz,scrRat)
     if(file.eq.1) call pgpap(bmpsz,bmprat)
     if(file.ge.2) call pgpap(PSsz,PSrat)
     if(file.ge.2) call pgscf(fonttype)
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
           write(stdErr,'(/,A)')'  * Warning:  chains():  parameter '//trim(parNames(plPars(j)))//' is not defined, check plPars() in the input file.  Skipping...'
           cycle
        end if
        
        call pgpage
        if(j.eq.1) call pginitl(colour,file,whiteBG)
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
        
        if(changeVar.gt.0) then
           if(parID(p).eq.31) then
              if(ymin.lt.0..or.ymax.gt.24.) then
                 ymin = 0.
                 ymax = 24.
              end if
           end if
           if(parID(p).eq.41.or.parID(p).eq.54.or.parID(p).eq.73.or.parID(p).eq.83) then
              !ymin = max(ymin,0.)
              !ymax = min(ymax,360.)
              !write(stdOut,'(I5,2F10.5,$)')p,ymin,ymax
              !ymin = max(rev360(ymin),0.)
              !ymax = min(rev360(ymax),360.)
              !write(stdOut,'(2F10.5)')ymin,ymax
              if(ymin.lt.0..or.ymax.gt.360.) then
                 ymin = 0.
                 ymax = 360.
              end if
           end if
           if(parID(p).eq.52) then
              if(ymin.lt.0..or.ymax.gt.180.) then
                 ymin = 0.
                 ymax = 180.
              end if
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
              if(changeVar.gt.0) then
                 if(parID(p).eq.31) plx = rev24(plx)
                 if(parID(p).eq.41.or.parID(p).eq.54.or.parID(p).eq.73.or.parID(p).eq.83) plx = rev360(plx)
                 if(parID(p).eq.52) plx = rev180(plx)
                 !call pgpoint(1,is(ic,i),plx,1) !Plot small dots
              end if
              !call pgpoint(1,plx,ply,symbol) !Plot symbols
              call pgpoint(1,plx,exp(ply-ymin),symbol) !Plot symbols
              !print*,i,plx,ply,exp(ply-ymin)
           end do
        end do
        call pgsch(sch)
        call pgslw(lw)
        
        !Plot max posterior
        if(plLmax.ge.1) then
           plx = allDat(icloglmax,p,iloglmax)
           if(changeVar.gt.0) then
              if(parID(p).eq.31) plx = rev24(plx)
              if(parID(p).eq.41.or.parID(p).eq.54.or.parID(p).eq.73.or.parID(p).eq.83) plx = rev360(plx)
              if(parID(p).eq.52) plx = rev180(plx)
           end if
           ply = exp(post(icloglmax,iloglmax)-ymin)
           call pgsci(1)
           call pgpoint(1,plx,ply,12)
           call pgsls(5)
           call pgline(2,(/plx,plx/),(/-1.e20,1.e20/))
        end if
        
        
        !Plot injection values
        do ic=1,nChains0
           call pgsls(2)
           call pgsci(1)
           
           !Plot injection values
           if(plInject.ge.1) then
              if(mergeChains.ne.1.or.ic.eq.1) then !The units of the injection values haven't changed (e.g. from rad to deg) for ic>1 (but they have for the starting values, why?)
                 !if(ic.eq.1) then !The units of the injection values haven't changed (e.g. from rad to deg) for ic>1 (but they have for the starting values, why?)
                 plx = startval(ic,p,1) !Injection value
                 plx = max(min(1.e30,startval(ic,p,1)),1.e-30)
                 if(changeVar.gt.0) then
                    if(parID(p).eq.31) plx = rev24(plx)
                    if(parID(p).eq.41.or.parID(p).eq.54.or.parID(p).eq.73.or.parID(p).eq.83) plx = rev360(plx)
                    if(parID(p).eq.52) plx = rev180(plx)
                 end if
                 call pgline(2,(/plx,plx/),(/-1.e20,1.e20/))
                 if(changeVar.gt.0) then
                    if(parID(p).eq.31) then
                       call pgline(2,(/plx-24.,plx-24./),(/-1.e20,1.e20/))
                       call pgline(2,(/plx+24.,plx+24./),(/-1.e20,1.e20/))
                    end if
                    if(parID(p).eq.41.or.parID(p).eq.54.or.parID(p).eq.73.or.parID(p).eq.83) then
                       call pgline(2,(/plx-360.,plx-360./),(/-1.e20,1.e20/))
                       call pgline(2,(/plx+360.,plx+360./),(/-1.e20,1.e20/))
                    end if
                    if(parID(p).eq.52) then
                       call pgline(2,(/plx-180.,plx-180./),(/-1.e20,1.e20/))
                       call pgline(2,(/plx+180.,plx+180./),(/-1.e20,1.e20/))
                    end if
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
     
     call pgend
     if(file.ge.2) then
        if(file.eq.3) then
           i = system('eps2pdf parlogl.eps  -o '//trim(outputdir)//'/'//trim(outputname)//'__parlogl.pdf   >& /dev/null')
           if(i.ne.0) write(stdErr,'(A,I6)')'  Error converting plot',i
        end if
        i = system('mv -f parlogl.eps '//trim(outputdir)//'/'//trim(outputname)//'__parlogl.eps')
     end if
     if(file.eq.1) then
        i = system('convert -resize '//trim(bmpxpix)//' -depth 8 -unsharp '//trim(unSharpchain)//' parlogl.ppm  '//trim(outputdir)//'/'//trim(outputname)//'__parlogl.png')
        if(i.ne.0) write(stdErr,'(A,I6)')'  Error converting plot',i
        i = system('rm -f parlogl.ppm')
     end if
  end if !if(plParL.eq.1)
  
  
  
  
  
  
  
  
  !***********************************************************************************************************************************            
  !Plot jump sizes
  if(plJump.ge.1) then
     !if(prProgress.ge.1.and.update.eq.0) write(stdOut,'(A)')' Plotting jump sizes...'
     if(prProgress.ge.1.and.update.eq.0) write(stdOut,'(A,$)')' jump sizes, '
     if(file.eq.0) then
        io = pgopen('18/xs')
        sch = fontsize1d*1.5
     end if
     if(file.ge.1) then
        if(file.eq.1) io = pgopen('jumps.ppm/ppm')
        if(file.ge.2) io = pgopen('jumps.eps'//trim(psclr))
        sch = fontsize1d*1.2
     end if
     if(io.le.0) then
        write(stdErr,'(A,I4)')'   Error:  Cannot open PGPlot device.  Quitting the programme',io
        exitcode = 1
        return
     end if
     if(file.eq.0) call pgpap(scrSz,scrRat)
     if(file.eq.0) sch = fontsize1d*1.5
     if(file.eq.1) call pgpap(bmpsz,bmprat)
     if(file.eq.1) sch = fontsize1d*1.5
     if(quality.eq.0) then !Draft
        !sch = sch*1.75
        sch = sch*1.4
        lw = 2
        call pgsvp(0.1,0.95,0.06,0.87) !To make room for title
     end if
     
     call pgsch(sch)
     
     call pgsubp(panels(1),panels(2))
     
     ic = 1
     do j=1,nPlPar
        p = revID(plPars(j))
        if(p.eq.0) then
           write(stdErr,'(/,A)')'  * Warning:  chains():  parameter '//trim(parNames(plPars(j)))//' is not defined, check plPars() in the input file.  Skipping...'
           cycle
        end if
        
        call pgpage
        if(j.eq.1) call pginitl(colour,file,whiteBG)
        xmax = -1.e30
        ymin =  1.e30
        ymax = -1.e30
        do ic=1,nChains0
           xmin = 0.
           xmax = max(xmax,maxval(is(ic,1:Ntot(ic))))
           dx = abs(xmax-xmin)*0.01
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
           !print*,p-1,ymin,ymax,dy
           dy = abs(ymax-ymin)*0.05
           if(dy.lt.1.e-10) dy = ymin*0.1
        end do
        
        call pgswin(xmin-dx,xmax+dx,ymin-dy,ymax+dy)
        if(plJump.eq.1) call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0) !lin
        if(plJump.eq.2) call pgbox('BCNTS',0.0,0,'BCLNTS',0.0,0) !log
        
        do ic=1,nChains0
           !call pgsci(mod(ic*2,10))
           call pgsci(defcolour)
           if(nChains0.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
           if(plJump.eq.1) then
              !do i=1,Ntot(ic),chainPlI
              do i=ic,Ntot(ic),chainPlI !Start at ic to reduce overplotting
                 call pgpoint(1,is(ic,i),jumps(ic,p,i),1)
              end do
           else
              !do i=1,Ntot(ic),chainPlI
              do i=ic,Ntot(ic),chainPlI !Start at ic to reduce overplotting
                 call pgpoint(1,is(ic,i),log10(abs(jumps(ic,p,i))+1.e-30),1)
              end do
           end if
        end do
        
        call pgsls(2)
        call pgsci(6)
        do ic=1,nChains0
           !call pgline(2,(/real(Nburn(ic)),real(Nburn(ic))/),(/-1.e20,1.e20/))
           if(nChains0.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
           !if(plBurn.ge.1) call pgline(2,(/isburn(ic),isburn(ic)/),(/-1.e20,1.e20/))
           if(plBurn.ge.1.and.isburn(ic).lt.is(ic,Ntot(ic))) call pgline(2,(/isburn(ic),isburn(ic)/),(/-1.e20,1.e20/))
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
     
     call pgend
     if(file.ge.2) then
        if(file.eq.3) then
           i = system('eps2pdf jumps.eps  -o '//trim(outputdir)//'/'//trim(outputname)//'__jumps.pdf   >& /dev/null')
           if(i.ne.0) write(stdErr,'(A,I6)')'  Error converting plot',i
        end if
        i = system('mv -f jumps.eps '//trim(outputdir)//'/'//trim(outputname)//'__jumps.eps')
     end if
     if(file.eq.1) then
        i = system('convert -resize '//trim(bmpxpix)//' -depth 8 -unsharp '//trim(unSharpchain)//' jumps.ppm  '//trim(outputdir)//'/'//trim(outputname)//'__jumps.png')
        if(i.ne.0) write(stdErr,'(A,I6)')'  Error converting plot',i
        i = system('rm -f jumps.ppm')
     end if
  end if !if(plJump.ge.1)
  
  
  
  
  
  
  
  !***********************************************************************************************************************************      
  !Plot autocorrelations for each parameter
  if(plAcorr.gt.0) then
     !if(prProgress.ge.1.and.update.eq.0) write(stdOut,'(A)')' Plotting autocorrelations...'
     if(prProgress.ge.1.and.update.eq.0) write(stdOut,'(A,$)')' autocorrelations, '
     if(file.eq.0) then
        io = pgopen('19/xs')
        sch = fontsize1d*1.5
        sch = sch*1.5
        call pgpap(scrSz,scrRat)
     end if
     if(file.eq.1) then
        io = pgopen('acorrs.ppm/ppm')
        sch = fontsize1d*1.2
        sch = sch*1.5
        call pgpap(bmpsz,bmprat)
     end if
     if(file.ge.2) then
        io = pgopen('acorrs.eps'//trim(psclr))
        sch = fontsize1d*1.2
        sch = sch*1.5
        call pgpap(PSsz,PSrat)
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
           write(stdErr,'(/,A)')'  * Warning:  chains():  parameter '//trim(parNames(plPars(j)))//' is not defined, check plPars() in the input file.  Skipping...'
           cycle
        end if
        
        call pgpage
        if(j.eq.1) call pginitl(colour,file,whiteBG)
        
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
           do i=1,nAcorr
              call pgpoint(1,acorrs(ic,0,i),acorrs(ic,p,i),1)
           end do
        end do
        
        !Plot autocorrelation length:
        call pgsls(2)
        call pgsci(defcolour)
        do ic=1,nChains0
           if(nChains0.gt.1) call pgsci(colours(mod(ic-1,ncolours)+1))
           call pgline(2,(/lAcorrs(ic,p),lAcorrs(ic,p)/),(/-1.e20,1.e20/))
        end do
        
        !Plot horizontal line at 0
        call pgsci(1)
        call pgline(2,(/-1.e20,1.e20/),(/0.,0./))
        call pgsci(1)
        call pgsls(1)
        !write(title,'(A,ES9.2)')'Autocorr.: '//trim(pgParNss(parID(p)))//', mean length:',sum(lAcorrs(1:nChains0,p))/real(nChains0)
        write(title,'(A,ES9.2)')'Autocorr.: '//trim(pgParNss(parID(p)))//', med. length:',compute_median_real(lAcorrs(1:nChains0,p),nChains0)
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
     
     call pgend
     if(file.ge.2) then
        if(file.eq.3) then
           i = system('eps2pdf acorrs.eps  -o '//trim(outputdir)//'/'//trim(outputname)//'__acorrs.pdf   >& /dev/null')
           if(i.ne.0) write(stdErr,'(A,I6)')'  Error converting plot',i
        end if
        i = system('mv -f acorrs.eps '//trim(outputdir)//'/'//trim(outputname)//'__acorrs.eps')
     end if
     if(file.eq.1) then
        i = system('convert -resize '//trim(bmpxpix)//' -depth 8 -unsharp '//trim(unSharpchain)//' acorrs.ppm  '//trim(outputdir)//'/'//trim(outputname)//'__acorrs.png')
        if(i.ne.0) write(stdErr,'(A,I6)')'  Error converting plot',i
        i = system('rm -f acorrs.ppm')
     end if
  end if !if(plAcorr.gt.0)
  
  
end subroutine chains
!***********************************************************************************************************************************      
