!> \mainpage Documentation analyseMCMC
!! AnalyseMCMC is a Fortran code that can be used to analyse the output of <a href="http://www.astro.northwestern.edu/~sluys/index.php?title=SPINspiral">SPINspiral</a>.
!! 
!! \file analysemcmc.f
!! \brief Contains analyseMCMC main routine
!<

! This program replaces plotspins.



!> \brief Main routine of AnalyseMCMC
!<
!************************************************************************************************************************************
program analyseMCMC
  use constants
  use analysemcmc_settings
  use general_data
  use mcmcrun_data
  use stats_data
  use plot_data
  use chain_data
  implicit none
  integer :: i,ic,p,iargc,exitcode,tempintarray(99),getos
  real :: pltsz
  real*8 :: timestamp,timestamps(9)  !< Time the progress of the code.
  
  wikioutput = 1  !Produce output for CBC Wiki: 0-no, 1-yes (requires one of the probability intervals to be 2-sigma)
  map_projection = 1  !Choose map projection: 1-Mollweide
  
  timestamps(1) = timestamp(os)
  write(6,*)
  
  os = getos() !1-Linux, 2-MacOS
  
  call setconstants           !Define mathematical constants
  call set_plotsettings()     !Set plot settings to 'default' values
  call read_settingsfile()    !Read the plot settings (overwrite the defaults)
  call write_settingsfile()   !Write the input file back to disc
  
  fontsize1d = 1.             !Set plot scale for 1D plots, needs to be implemented fully
  fontsize2d = 1.             !Set plot scale for 2D plots, needs to be implemented fully. Typically, take ~1.2*fontsize1d
  if(quality.eq.91) then !NINJA
     fontsize1d = 1.3
     fontsize2d = 1.55
  end if
  orientation = 1             !Use portrait (1) or landscape (2) for eps/pdf
  if(quality.eq.0) orientation = 2  !Easier to print eps
  fonttype = 1                !Font type used for eps/pdf: 1-simple, 2-roman, 3-italic 4-script
  
  nchains0 = iargc()
  if(nchains0.lt.1) then
     write(0,'(A,/)')'  Syntax: analysemcmc <file1> [file2] ...'
     stop
  end if
  if(nchains0.gt.nchs) write(0,'(A,I3,A)')'  *** WARNING:  Too many input files (chains), please increase nchs in analysemcmc_functions.f. Only',nchs,' files can be read.'
  nchains0 = min(nchains0,nchs)
  
  
  !Some of the stuff below will have to go to the input file
  
  par1 = 1          !First parameter to treat (stats, plot): 0-all
  par2 = 15         !Last parameter to treat (0: use npar)
  if(version.eq.2) par2 = 16
  
  maxdots = 25000  !~Maximum number of dots to plot in e.g. chains plot, to prevent dots from being overplotted too much and eps/pdf files from becoming huge.  Use this to autoset chainpli
  
  
  !Determine plot sizes and ratios:   (ratio ~ y/x and usually < 1 ('landscape'))
  bmpsz = real(bmpxsz-1)/85. * scfac !Make png larger, so that convert interpolates and makes the plot smoother
  bmprat = real(bmpysz-1)/real(bmpxsz-1)
  write(bmpxpix,'(I4)')bmpxsz  !Used as a text string by convert
  if(file.eq.0) pltsz = scrsz
  if(file.eq.0) pltrat = scrrat
  if(file.eq.1) pltsz = bmpsz
  if(file.eq.1) pltrat = bmprat
  if(file.ge.2) pltsz = pssz
  if(file.ge.2) pltrat = psrat
  
  
  !Use full unsharp-mask strength for plots with many panels and dots, weaker for those with fewer panels and/or no dots
  write(unsharplogl,'(I4)')max(nint(real(unsharp)/2.),1)  !Only one panel with dots
  write(unsharpchain,'(I4)')unsharp                       !~12 panels with dots
  write(unsharppdf1d,'(I4)')max(nint(real(unsharp)/2.),1) !~12 panels, no dots
  write(unsharppdf2d,'(I4)')max(nint(real(unsharp)/4.),1) !1 panel, no dots
  
  outputdir = '.'  !Directory where output is saved (either relative or absolute path)
  
  
  
  
  
  
  
  !Sort out implicit options:
  if(panels(1)*panels(2).lt.nplvar) panels = 0
  if(panels(1)*panels(2).lt.1) then
     if(nplvar.eq.1) panels = (/1,1/)
     if(nplvar.eq.2) panels = (/2,1/)
     if(nplvar.eq.3) panels = (/3,1/)
     if(nplvar.eq.4) panels = (/2,2/)
     if(nplvar.eq.5) panels = (/5,1/)
     if(nplvar.eq.6) panels = (/3,2/)
     if(nplvar.eq.7) panels = (/4,2/)
     if(nplvar.eq.8) panels = (/4,2/)
     if(nplvar.eq.9) panels = (/3,3/)
     if(nplvar.eq.10) panels = (/5,2/)
     if(nplvar.eq.11) panels = (/4,3/)
     if(nplvar.eq.12) panels = (/4,3/)
     if(nplvar.eq.12.and.quality.eq.3) panels = (/3,4/)
     if(nplvar.eq.13) panels = (/5,3/)
     if(nplvar.eq.14) panels = (/5,3/)
     if(nplvar.eq.15) panels = (/5,3/)
     if(nplvar.eq.16) panels = (/4,4/)
     if(nplvar.eq.17) panels = (/6,3/)
     if(nplvar.eq.18) panels = (/6,3/)
     if(nplvar.eq.19) panels = (/5,4/)
     if(nplvar.eq.20) panels = (/5,4/)
  end if
  
  !eps/pdf: colour and orientation
  psclr = '/cps'
  if(colour.eq.0) psclr = '/ps'
  if(orientation.eq.1) then
     psclr = '/vcps'
     if(colour.eq.0) psclr = '/vps'
  end if
  
  ncolours = 5; colours(1:ncolours)=(/4,2,3,6,5/) !Paper
  ncolours = 10; colours(1:ncolours)=(/2,3,4,5,6,7,8,9,10,11/)
  nsymbols = 1; symbols(1:nsymbols)=(/chainsymbol/)
  if(colour.eq.1.and.quality.eq.2) then !Beamer
     ncolours = 5
     colours(1:ncolours)=(/4,2,5,11,15/)
  end if
  if(colour.ne.1) then
     ncolours=3
     colours(1:ncolours)=(/1,14,15/)
     !ncolours=6
     !colours(1:ncolours)=(/1,1,1,15,1,15/)
     if(chainsymbol.eq.-10) then
        nsymbols = 8
        symbols(1:nsymbols) = (/2,4,5,6,7,11,12,15/) !Thin/open symbols
     end if
     if(chainsymbol.eq.-11) then
        nsymbols = 6
        symbols(1:nsymbols) = (/-3,-4,16,17,18,-6/) !Filled symbols
     end if
     !print*,chainsymbol,nsymbols
  end if
  if(colour.eq.1.and.quality.eq.0.and.nchs.gt.5) then
     ncolours = 10
     colours(1:ncolours)=(/2,3,4,5,6,7,8,9,10,11/)
  end if
  !Overrule
  !ncolours = 1
  !colours(1:ncolours)=(/6/)
  !defcolour = 2 !Red e.g. in case of 1 chain
  defcolour = colours(1)
  
  
  if(reverseread.ge.2) then !Reverse colours too
     do i=1,ncolours
        tempintarray(i) = colours(i)
     end do
     do i=1,nchains0
        colours(i) = tempintarray(nchains0-i+1) !Reverse colours too, but use the same first nchains0 from the set
     end do
  end if
  
  if(plot.eq.0) then
     pllogl = 0
     plchain = 0
     pljump = 0
     plsigacc = 0
     if(savepdf.eq.0) then
        plpdf1d = 0
        plpdf2d = 0
     end if
     plmovie = 0
  end if
  if(savepdf.eq.1) then
     if(nplvar.ne.15) write(0,'(/,A)')'*** WARNING:  I changed nplvar to 15, since savepdf is selected ***'
     nplvar = 15; plvars(1:nplvar) = (/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15/) !All 12 + m1,m2
     wrapdata = 0
  end if
  !if(par1.lt.2) par1 = 2
  if(par1.lt.1) par1 = 1 !Include log(L)
  if(plsigacc.ge.1.or.plmovie.ge.1) rdsigacc = 1
  if(file.eq.1) combinechainplots = 0
  if(file.ge.1) update = 0
  if(plmovie.ge.1) update = 0
  
  colournames(1:15) = (/'white','red','dark green','dark blue','cyan','magenta','yellow','orange','light green','brown','dark red','purple','red-purple','dark grey','light grey'/)
  if(file.ge.2) colournames(1) = 'black'
  
  
  
  !Columns in dat() (12 par): 1:logL 2:mc, 3:eta, 4:tc, 5:logdl, 6:spin, 7:kappa, 8: RA, 9:sindec,10:phase, 11:sinthJ0, 12:phiJ0, 13:alpha
  !Columns in dat() (15 par): 1:logL 2:mc, 3:eta, 4:t0, 5:logdl, 6:RA, 7:sindec, 8: cosi, 9:phase,10:psi, 11:spin1, 12:theta1, 13:phi1, 14:spin2, 15:theta2, 16:phi2
  if(fonttype.eq.2) then  !Use 'roman-like' Greek font
     varnames(1:15) = (/'logL','Mc','eta','tc','log_dl','spin','kappa','RA','sin_dec','phase','sin_thJo','phJo','alpha','M1','M2'/)
     pgvarns(1:15)  = (/'log Likelihood        ','\(2563) (M\d\(2281)\u) ','\(2133)               ','t\dc\u (s)            ', &
          'logd\dL\u (Mpc)       ','a\dspin\u             ','\(2136)               ','\(2127) (rad)         ', &
          'sin \(2130)           ','\(2147)\dc\u (rad)    ','sin \(2185)\dJ0\u     ','\(2147)\dJ0\u (rad)   ', &
          '\(2127)\dc\u (rad)    ','M\d1\u (M\d\(2281)\u) ','M\d2\u (M\d\(2281)\u) '/)
     pgvarnss(1:15)  = (/'log L    ','\(2563) ','\(2133)','t\dc\u','log d\dL\u','a\dspin\u','\(2136)','\(2127)','sin \(0630)','\(2147)\dc\u', &
          'sin \(2185)\dJ0\u','\(2147)\dJ0\u','\(2127)\dc\u','M\d1\u','M\d2\u'/)
     pgorigvarns(1:15)  = (/'log Likelihood        ','\(2563) (M\d\(2281)\u) ','\(2133)               ','t\dc\u (s)            ', &
          'logd\dL\u (Mpc)       ','a\dspin\u             ','\(2136)               ','\(2127) (rad)         ', &
          'sin \(2130)           ','\(2147)\dc\u (rad)    ','sin \(2185)\dJ0\u     ','\(2147)\dJ0\u (rad)   ', &
          '\(2127)\dc\u (rad)    ','M\d1\u (M\d\(2281)\u) ','M\d2\u (M\d\(2281)\u) '/)
     pgunits(1:15)  = (/'','M\d\(2281)\u ','','s','Mpc','','rad','rad','','rad','','rad','rad','M\d\(2281)\u','M\d\(2281)\u'/)
     
     if(version.eq.2) then !15par, 2 spins
        varnames(1:16) = (/'logL','Mc','eta','t0','log_dl','RA','sin_dec','cosi','phase','psi','spin1','th1','phi1','spin2','th2','phi2'/)
        pgvarns(1:16)  = (/'log Likelihood        ','\(2563) (M\d\(2281)\u) ','\(2133)               ','t\d0\u (s)            ', &
             'log d\dL\u (Mpc)      ','\(2127) (rad)         ','sin \(2130)           ','cos \(2135)           ', &
             '\(2147)\dc\u (rad)    ','\(2149) (rad)         ','a\dspin1\u            ','cos \(2185)\d1\u      ', &
             '\(2147)\d1\u (rad)    ','a\dspin2\u (rad)      ','cos \(2185)\d2\u      ','\(2147)\d2\u (rad)    '/)
        pgvarnss(1:16)  = (/'log L    ','\(2563) ','\(2133)','t\dc\u','log d\dL\u','\(2127)','sin \(0630)','cos \(2135)','\(2147)\dc\u', '\(2149)', &
             'a\dspin1\u','cos \(2185)\d1\u','\(2147)\d1\u','a\dspin2\u','cos \(2185)\d2\u','\(2147)\d2\u'/)
        pgorigvarns(1:16) = pgvarns(1:16)
        pgunits(1:16)  = (/'','M\d\(2281)\u ','','s','Mpc','rad','','','rad','rad','','rad','rad','','rad','rad'/)
     end if
  else  !Same, but replace '\(21' with \(06' for arial-like Greek font
     varnames(1:15) = (/'logL','Mc','eta','tc','log_dl','spin','kappa','RA','sin_dec','phase','sin_thJo','phJo','alpha','M1','M2'/)
     pgvarns(1:15)  = (/'log Likelihood        ','\(2563) (M\d\(2281)\u) ','\(0633)               ','t\dc\u (s)            ', &
          'logd\dL\u (Mpc)       ','a\dspin\u             ','\(0636)               ','\(0627) (rad)         ', &
          'sin \(0630)           ','\(0647)\dc\u (rad)    ','sin \(0685)\dJ0\u     ','\(0647)\dJ0\u (rad)   ', &
          '\(0627)\dc\u (rad)    ','M\d1\u (M\d\(2281)\u) ','M\d2\u (M\d\(2281)\u) '/)
     pgvarnss(1:15)  = (/'log L    ','\(2563) ','\(0633)','t\dc\u','log d\dL\u','a\dspin\u','\(0636)','\(0627)','sin \(0630)','\(0647)\dc\u', &
          'sin \(0685)\dJ0\u','\(0647)\dJ0\u','\(0627)\dc\u','M\d1\u','M\d2\u'/)
     pgorigvarns(1:15)  = (/'log Likelihood        ','\(2563) (M\d\(2281)\u) ','\(0633)               ','t\dc\u (s)            ', &
          'logd\dL\u (Mpc)       ','a\dspin\u             ','\(0636)               ','\(0627) (rad)         ', &
          'sin \(0630)           ','\(0647)\dc\u (rad)    ','sin \(0685)\dJ0\u     ','\(0647)\dJ0\u (rad)   ', &
          '\(0627)\dc\u (rad)    ','M\d1\u (M\d\(2281)\u) ','M\d2\u (M\d\(2281)\u) '/)
     pgunits(1:15)  = (/'','M\d\(2281)\u ','','s','Mpc','','rad','rad','','rad','','rad','rad','M\d\(2281)\u','M\d\(2281)\u'/)
     
     if(version.eq.2) then !15par, 2 spins
        varnames(1:16) = (/'logL','Mc','eta','t0','log_dl','RA','sin_dec','cosi','phase','psi','spin1','th1','phi1','spin2','th2','phi2'/)
        pgvarns(1:16)  = (/'log Likelihood        ','\(2563) (M\d\(2281)\u) ','\(0633)               ','t\d0\u (s)            ', &
             'log d\dL\u (Mpc)      ','\(0627) (rad)         ','sin \(0630)           ','cos \(0635)           ', &
             '\(0647)\dc\u (rad)    ','\(0649) (rad)         ','a\dspin1\u            ','cos \(0685)\d1\u      ', &
             '\(0647)\d1\u (rad)    ','a\dspin2\u (rad)      ','cos \(0685)\d2\u      ','\(0647)\d2\u (rad)    '/)
        pgvarnss(1:16)  = (/'log L    ','\(2563) ','\(0633)','t\dc\u','log d\dL\u','\(0627)','sin \(0630)','cos \(0635)','\(0647)\dc\u', '\(0649)', &
             'a\dspin1\u','cos \(0685)\d1\u','\(0647)\d1\u','a\dspin2\u','cos \(0685)\d2\u','\(0647)\d2\u'/)
        pgorigvarns(1:16) = pgvarns(1:16)
        pgunits(1:16)  = (/'','M\d\(2281)\u ','','s','Mpc','rad','','','rad','rad','','rad','rad','','rad','rad'/)
     end if
  end if
  
  
  
  !if(prprogress+prruninfo+prinitial.ge.1) write(6,*)
  npar = 13
  if(prchaininfo.ge.1) then
     if(nchains0.eq.1) then
        write(6,'(A)')'  Analysing 1 chain from SPINspiral'
     else
        write(6,'(A,I3,A)')'  Analysing',nchains0,' chains from SPINspiral'
     end if
  end if
  nchains = nchains0
  
  
  
  
  
  
  
  !*******************************************************************************************************************************
  !***   READ INPUT FILE(S)   ****************************************************************************************************
  !*******************************************************************************************************************************
  
101 continue
  !Read the input files:
  call read_mcmcfiles(exitcode)
  if(exitcode.ne.0) goto 9999
  
  
  !Get and print some basic chain statistics:
  timestamps(2) = timestamp(os)
  call mcmcruninfo(exitcode)
  
  
  
  
  
  ! **********************************************************************************************************************************
  ! ***  DO STATISTICS   *************************************************************************************************************
  ! **********************************************************************************************************************************
  
  timestamps(3) = timestamp(os)
  
  call statistics(exitcode)
  if(exitcode.ne.0) goto 9999
  
  
  
  
  
  
  
  !Change the original chain data
  if(changevar.ge.1) then
     if(version.eq.1) then
        do ic=1,nchains0
           !Columns in dat() (12 par): 1:logL 2:mc, 3:eta, 4:tc, 5:dl, 6:spin,  7:theta_SL, 8: RA,   9:dec, 10:phase, 11:thJ0, 12:phiJ0, 13:alpha
           !if(prprogress.ge.2.and.update.eq.0) write(6,'(A,$)')'Changing some variables...   '
           do p=par1,par2
              if(p.eq.5) pldat(ic,p,1:ntot(ic)) = exp(pldat(ic,p,1:ntot(ic)))
              !if(p.eq.9.or.p.eq.11) pldat(ic,p,1:ntot(ic)) = asin(pldat(ic,p,1:ntot(ic)))*r2d
              if(p.eq.9) pldat(ic,p,1:ntot(ic)) = asin(pldat(ic,p,1:ntot(ic)))*r2d
              if(p.eq.7) pldat(ic,p,1:ntot(ic)) = acos(pldat(ic,p,1:ntot(ic)))*r2d
              if(p.eq.8) pldat(ic,p,1:ntot(ic)) = pldat(ic,p,1:ntot(ic))*r2h
              !if(p.eq.10.or.p.eq.12.or.p.eq.13) pldat(ic,p,1:ntot(ic)) = pldat(ic,p,1:ntot(ic))*r2d
              if(p.ge.10.and.p.le.13) pldat(ic,p,1:ntot(ic)) = pldat(ic,p,1:ntot(ic))*r2d
           end do !p
        end do
     end if
     if(version.eq.2) then
        do ic=1,nchains0
           !Columns in dat() (15 par): 1:logL 2:mc, 3:eta, 4:t0, 5:logdl, 6:RA, 7:sindec, 8: cosi, 9:phase,10:psi, 11:spin1, 12:theta1, 13:phi1, 14:spin2, 15:theta2, 16:phi2
           !if(prprogress.ge.2.and.update.eq.0) write(6,'(A,$)')'Changing some variables...   '
           do p=par1,par2
              if(p.eq.5) pldat(ic,p,1:ntot(ic)) = exp(pldat(ic,p,1:ntot(ic)))
              if(p.eq.8.or.p.eq.12.or.p.eq.15) pldat(ic,p,1:ntot(ic)) = acos(pldat(ic,p,1:ntot(ic)))*r2d
              if(p.eq.6) pldat(ic,p,1:ntot(ic)) = pldat(ic,p,1:ntot(ic))*r2h
              if(p.eq.7) pldat(ic,p,1:ntot(ic)) = asin(pldat(ic,p,1:ntot(ic)))*r2d
              if(p.eq.9.or.p.eq.10.or.p.eq.13.or.p.eq.16) pldat(ic,p,1:ntot(ic)) = pldat(ic,p,1:ntot(ic))*r2d
           end do !p
        end do
     end if
     !if(prprogress.ge.2.and.update.eq.0) write(6,'(A)')'  Done.'
  end if !if(changevar.ge.1)
  
  
  
  deallocate(dat)
  
  
  
  ! **********************************************************************************************************************************
  ! ***  CREATE PLOTS   **************************************************************************************************************
  ! **********************************************************************************************************************************
  
  timestamps(4) = timestamp(os)
  
  if(prprogress.ge.2) write(6,*)''
  if(plot.eq.1.and.prprogress.ge.1.and.update.eq.0) then
     write(6,'(/,A,$)')'  Plotting '
     if(file.eq.0) write(6,'(A,$)')'to screen: '
     if(file.eq.1) write(6,'(A,$)')'to png: '
     if(file.eq.2) write(6,'(A,$)')'to eps: '
     if(file.eq.3) write(6,'(A,$)')'to pdf: '
  end if
  
  
  
  !***********************************************************************************************************************************      
  !Plot (1d) chains: logL, parameter chains, jumps, etc.
  if(plot.eq.1) then
     call chains(exitcode)
     if(exitcode.ne.0) goto 9999
  end if
  timestamps(5) = timestamp(os)
  
  
  
  
  !***********************************************************************************************************************************      
  !Plot pdfs (1d)
  if(plpdf1d.ge.1) then
     call pdfs1d(exitcode)
     if(exitcode.ne.0) goto 9999
  end if !if(plpdf1d.ge.1)
  timestamps(6) = timestamp(os)
  
  
  
  
  !***********************************************************************************************************************************      
  if(plpdf2d.ge.1.and.mergechains.eq.0) then
     write(6,'(A,$)')', (skipping 2D PDFs since mergechains=0), '
     plpdf2d = 0
  end if
  
  if(plpdf2d.ge.1) then
     call pdfs2d(exitcode)
     if(exitcode.ne.0) goto 9999
  end if !if(plpdf2d.eq.1)
  
  if(npdf2d.lt.0) then !Then we just plotted all 2D PDFs
     write(6,*)
  else
     if(prprogress.ge.1.and.update.eq.0.and.plot.gt.0) write(6,'(A,/)')' done.  '
  end if
  
  timestamps(7) = timestamp(os)  
  
  
  
  
  !***********************************************************************************************************************************      
  
  !Write statistics to file
  if(savestats.ge.1.and.nchains.gt.1) write(0,'(A)')' ******   Cannot write statistics if the number of chains is greater than one   ******'
  if(savestats.ge.1.and.nchains.eq.1) then
     call save_stats(exitcode)
     if(exitcode.ne.0) goto 9999
     write(6,*)''
  end if !if(savestats.ge.1.and.nchains.eq.1) then
  
  
  
  
  !***********************************************************************************************************************************      
  
  timestamps(8) = timestamp(os)
  
  if(plmovie.ge.1) then
     call animation(exitcode)
     if(exitcode.ne.0) goto 9999
  end if
  
  
  
  
  if(update.eq.1) then
     deallocate(pldat,alldat,post,prior)
     call sleep(5)
     if(sum(ntot).gt.1.e4) call sleep(5)
     if(sum(ntot).gt.1.e5) call sleep(10)
     if(sum(ntot).gt.1.e6) call sleep(20)
     goto 101
  end if
  
  !write(6,'(A)')'  Waiting for you to finish me off...'
  !pause
  
9999 continue
  deallocate(pldat,alldat,post,prior)
  !if(prprogress.ge.1) write(6,*)''
  
  timestamps(9) = timestamp(os)
  
  if(prprogress.ge.1) then
     write(6,'(A,$)')'  Run time: '
     write(6,'(A,F5.1,A,$)')'   input:',min(dabs(timestamps(2)-timestamps(1)),999.9),'s,'
     !write(6,'(A,F5.1,A,$)')'   info:',min(dabs(timestamps(3)-timestamps(2)),999.9),'s,'
     !write(6,'(A,F5.1,A,$)')'   stats:',min(dabs(timestamps(4)-timestamps(3)),999.9),'s,'
     write(6,'(A,F5.1,A,$)')'   stats:',min(dabs(timestamps(4)-timestamps(2)),999.9),'s,'
     if(plot.eq.1.and.pllogl+plchain+pljump+plsigacc+placorr.gt.0) then
        write(6,'(A,F5.1,A,$)')'   chains:',min(dabs(timestamps(5)-timestamps(4)),999.9),'s,'
     end if
     if(plot.eq.1.or.savepdf.ge.1) then
        if(plpdf1d.ge.1) write(6,'(A,F5.1,A,$)')'   1d pdfs:',min(dabs(timestamps(6)-timestamps(5)),999.9),'s,'
        if(plpdf2d.ge.1) write(6,'(A,F6.1,A,$)')'   2d pdfs:',min(dabs(timestamps(7)-timestamps(6)),999.9),'s,'
     end if
     !write(6,'(A,F6.1,A,$)')'   plots:',min(dabs(timestamps(7)-timestamps(4)),999.9),'s,'
     !write(6,'(A,F5.1,A,$)')'   save stats:',min(dabs(timestamps(8)-timestamps(7)),999.9),'s,'
     if(plmovie.ge.1) write(6,'(A,F5.1,A,$)')'   movie:',min(dabs(timestamps(9)-timestamps(8)),999.9),'s,'
     write(6,'(A,F6.1,A)')'   total:',min(dabs(timestamps(9)-timestamps(1)),999.9),'s.'
  end if
  
  write(6,*)''
end program analyseMCMC
!************************************************************************************************************************************





