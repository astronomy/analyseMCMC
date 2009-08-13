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
  integer :: i,iargc,exitcode,tempintarray(99),getos
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
  if(nchains0.gt.maxChs) write(0,'(A,I3,A)')'  *** WARNING:  Too many input files (chains), please increase maxChs in analysemcmc_functions.f. Only',maxChs,' files can be read.'
  nchains0 = min(nchains0,maxChs)
  
  
  !Some of the stuff below will have to go to the input file
  
  maxdots = 25000  !~Maximum number of dots to plot in e.g. chains plot, to prevent dots from being overplotted too much and eps/pdf files from becoming huge.  Use this to autoset chainPlI
  
  
  !Determine plot sizes and ratios:   (ratio ~ y/x and usually < 1 ('landscape'))
  bmpsz = real(bmpXSz-1)/85. * scFac !Make png larger, so that convert interpolates and makes the plot smoother
  bmprat = real(bmpYSz-1)/real(bmpXSz-1)
  write(bmpxpix,'(I4)')bmpXSz  !Used as a text string by convert
  if(file.eq.0) pltsz = scrSz
  if(file.eq.0) pltrat = scrRat
  if(file.eq.1) pltsz = bmpsz
  if(file.eq.1) pltrat = bmprat
  if(file.ge.2) pltsz = PSsz
  if(file.ge.2) pltrat = PSrat
  
  
  !Use full unsharp-mask strength for plots with many panels and dots, weaker for those with fewer panels and/or no dots
  write(unSharplogl,'(I4)')max(nint(real(unSharp)/2.),1)  !Only one panel with dots
  write(unSharpchain,'(I4)')unSharp                       !~12 panels with dots
  write(unSharppdf1d,'(I4)')max(nint(real(unSharp)/2.),1) !~12 panels, no dots
  write(unSharppdf2d,'(I4)')max(nint(real(unSharp)/4.),1) !1 panel, no dots
  
  outputdir = '.'  !Directory where output is saved (either relative or absolute path)
  
  
  
  
  
  
  
  !Sort out implicit options:
  if(panels(1)*panels(2).lt.nPlPar) panels = 0
  if(panels(1)*panels(2).lt.1) then
     if(nPlPar.eq.1) panels = (/1,1/)
     if(nPlPar.eq.2) panels = (/2,1/)
     if(nPlPar.eq.3) panels = (/3,1/)
     if(nPlPar.eq.4) panels = (/2,2/)
     if(nPlPar.eq.5) panels = (/5,1/)
     if(nPlPar.eq.6) panels = (/3,2/)
     if(nPlPar.eq.7) panels = (/4,2/)
     if(nPlPar.eq.8) panels = (/4,2/)
     if(nPlPar.eq.9) panels = (/3,3/)
     if(nPlPar.eq.10) panels = (/5,2/)
     if(nPlPar.eq.11) panels = (/4,3/)
     if(nPlPar.eq.12) panels = (/4,3/)
     if(nPlPar.eq.12.and.quality.eq.3) panels = (/3,4/)
     if(nPlPar.eq.13) panels = (/5,3/)
     if(nPlPar.eq.14) panels = (/5,3/)
     if(nPlPar.eq.15) panels = (/5,3/)
     if(nPlPar.eq.16) panels = (/4,4/)
     if(nPlPar.eq.17) panels = (/6,3/)
     if(nPlPar.eq.18) panels = (/6,3/)
     if(nPlPar.eq.19) panels = (/5,4/)
     if(nPlPar.eq.20) panels = (/5,4/)
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
  nsymbols = 1; symbols(1:nsymbols)=(/chainSymbol/)
  if(colour.eq.1.and.quality.eq.2) then !Beamer
     ncolours = 5
     colours(1:ncolours)=(/4,2,5,11,15/)
  end if
  if(colour.ne.1) then
     ncolours=3
     colours(1:ncolours)=(/1,14,15/)
     !ncolours=6
     !colours(1:ncolours)=(/1,1,1,15,1,15/)
     if(chainSymbol.eq.-10) then
        nsymbols = 8
        symbols(1:nsymbols) = (/2,4,5,6,7,11,12,15/) !Thin/open symbols
     end if
     if(chainSymbol.eq.-11) then
        nsymbols = 6
        symbols(1:nsymbols) = (/-3,-4,16,17,18,-6/) !Filled symbols
     end if
     !print*,chainSymbol,nsymbols
  end if
  if(colour.eq.1.and.quality.eq.0.and.maxChs.gt.5) then
     ncolours = 10
     colours(1:ncolours)=(/2,3,4,5,6,7,8,9,10,11/)
  end if
  !Overrule
  !ncolours = 1
  !colours(1:ncolours)=(/6/)
  !defcolour = 2 !Red e.g. in case of 1 chain
  defcolour = colours(1)
  
  
  if(reverseRead.ge.2) then !Reverse colours too
     do i=1,ncolours
        tempintarray(i) = colours(i)
     end do
     do i=1,nchains0
        colours(i) = tempintarray(nchains0-i+1) !Reverse colours too, but use the same first nchains0 from the set
     end do
  end if
  
  if(plot.eq.0) then
     plLogL = 0
     plChain = 0
     plJump = 0
     if(savePDF.eq.0) then
        plPDF1D = 0
        plPDF2D = 0
     end if
     plAnim = 0
  end if
  if(savePDF.eq.1) then
     !if(nPlPar.ne.15) write(0,'(/,A)')'*** WARNING:  I changed nPlPar to 15, since savePDF is selected ***'
     !nPlPar = 15; plPars(1:nPlPar) = (/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15/) !All 12 + m1,m2
     wrapData = 0
  end if
  if(file.eq.1) combineChainPlots = 0
  if(file.ge.1) update = 0
  if(plAnim.ge.1) update = 0
  
  colournames(1:15) = (/'white','red','dark green','dark blue','cyan','magenta','yellow','orange','light green','brown','dark red','purple','red-purple','dark grey','light grey'/)
  if(file.ge.2) colournames(1) = 'black'
  
  call set_originalParameterNames()  !Set the names and symbols of the original MCMC parameters in the database
  
  
  if(prChainInfo.ge.1) then
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
  
  
  
  
  
  ! **********************************************************************************************************************************
  ! ***  CREATE PLOTS   **************************************************************************************************************
  ! **********************************************************************************************************************************
  
  timestamps(4) = timestamp(os)
  
  if(prProgress.ge.2) write(6,*)''
  if(plot.eq.1.and.prProgress.ge.1.and.update.eq.0) then
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
  if(plPDF1D.ge.1) then
     call pdfs1d(exitcode)
     if(exitcode.ne.0) goto 9999
  end if !if(plPDF1D.ge.1)
  timestamps(6) = timestamp(os)
  
  
  
  
  !***********************************************************************************************************************************      
  if(plPDF2D.ge.1.and.mergeChains.eq.0) then
     write(6,'(A,$)')', (skipping 2D PDFs since mergeChains=0), '
     plPDF2D = 0
  end if
  
  if(plPDF2D.ge.1) then
     call pdfs2d(exitcode)
     if(exitcode.ne.0) goto 9999
  end if !if(plPDF2D.eq.1)
  
  if(Npdf2D.lt.0) then !Then we just plotted all 2D PDFs
     write(6,*)
  else
     if(prProgress.ge.1.and.update.eq.0.and.plot.gt.0) write(6,'(A,/)')' done.  '
  end if
  
  timestamps(7) = timestamp(os)  
  
  
  
  
  !***********************************************************************************************************************************      
  
  !Write statistics to file
  if(saveStats.ge.1.and.nchains.gt.1) write(0,'(A)')' ******   Cannot write statistics if the number of chains is greater than one   ******'
  if(saveStats.ge.1.and.nchains.eq.1) then
     call save_stats(exitcode)
     if(exitcode.ne.0) goto 9999
     write(6,*)''
  end if !if(saveStats.ge.1.and.nchains.eq.1) then
  
  
  
  
  !***********************************************************************************************************************************      
  
  timestamps(8) = timestamp(os)
  
  if(plAnim.ge.1) then
     call animation(exitcode)
     if(exitcode.ne.0) goto 9999
  end if
  
  
  
  if(update.eq.1) then
     deallocate(allDat,selDat,post,prior)
     call sleep(5)
     if(sum(ntot).gt.1.e4) call sleep(5)
     if(sum(ntot).gt.1.e5) call sleep(10)
     if(sum(ntot).gt.1.e6) call sleep(20)
     goto 101
  end if
  
  !write(6,'(A)')'  Waiting for you to finish me off...'
  !pause
  
9999 continue
  deallocate(allDat,selDat,post,prior)
  !if(prProgress.ge.1) write(6,*)''
  
  timestamps(9) = timestamp(os)
  
  if(prProgress.ge.1) then
     write(6,'(A,$)')'  Run time: '
     write(6,'(A,F5.1,A,$)')'   input:',min(dabs(timestamps(2)-timestamps(1)),999.9),'s,'
     !write(6,'(A,F5.1,A,$)')'   info:',min(dabs(timestamps(3)-timestamps(2)),999.9),'s,'
     !write(6,'(A,F5.1,A,$)')'   stats:',min(dabs(timestamps(4)-timestamps(3)),999.9),'s,'
     write(6,'(A,F5.1,A,$)')'   stats:',min(dabs(timestamps(4)-timestamps(2)),999.9),'s,'
     if(plot.eq.1.and.plLogL+plChain+plJump+plACorr.gt.0) then
        write(6,'(A,F5.1,A,$)')'   chains:',min(dabs(timestamps(5)-timestamps(4)),999.9),'s,'
     end if
     if(plot.eq.1.or.savePDF.ge.1) then
        if(plPDF1D.ge.1) write(6,'(A,F5.1,A,$)')'   1d pdfs:',min(dabs(timestamps(6)-timestamps(5)),999.9),'s,'
        if(plPDF2D.ge.1) write(6,'(A,F6.1,A,$)')'   2d pdfs:',min(dabs(timestamps(7)-timestamps(6)),999.9),'s,'
     end if
     !write(6,'(A,F6.1,A,$)')'   plots:',min(dabs(timestamps(7)-timestamps(4)),999.9),'s,'
     !write(6,'(A,F5.1,A,$)')'   save stats:',min(dabs(timestamps(8)-timestamps(7)),999.9),'s,'
     if(plAnim.ge.1) write(6,'(A,F5.1,A,$)')'   movie:',min(dabs(timestamps(9)-timestamps(8)),999.9),'s,'
     write(6,'(A,F6.1,A)')'   total:',min(dabs(timestamps(9)-timestamps(1)),999.9),'s.'
  end if
  
  write(6,*)''
end program analyseMCMC
!************************************************************************************************************************************





