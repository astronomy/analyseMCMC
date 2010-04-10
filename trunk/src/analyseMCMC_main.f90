!> \mainpage Documentation analyseMCMC
!! AnalyseMCMC is a Fortran code that can be used to analyse the output of 
!! <a href="http://www.phys.ualberta.ca/~sluys/index.php?title=SPINspiral">SPINspiral</a>.
!!
!! \file analyseMCMC_main.f90
!! 
!! \brief Contains analyseMCMC main routine
!!
!! I/O units used:
!! 
!!  0: StdErr
!!  6: StdOut
!!
!! 10: input (MCMC) file
!! 14: write settings file (analysemcmc.dat)
!! 15: read settings file (analysemcmc.dat)
!! 16: temp file in getos
!! 17: temp file in timestamp
!! 19: StdOut redirection (__output.txt)
!!
!! 20: Output: __statistics.dat, __wiki.dat, __bayes.dat 
!!
!! 21: bsc.dat (bright star catalogue)
!! 22: bsc_const.dat (BSC constellation data)
!! 23: bsc_names.dat (BSC names)
!! 24: milkyway*.dat (Milky Way data)
!! 
!! 30: Output: __pdf1/2d.dat
!! 40: Output: tailored output
!! 
!! 51: Output: HTML main file
!! 
!! This program replaces plotspins.
!<




!***********************************************************************************************************************************
program analyseMCMC
  use svn_revision
  use basic
  use constants
  use analysemcmc_settings
  use general_data
  use mcmcrun_data
  use stats_data
  use plot_data
  use chain_data
  
  implicit none
  integer :: i,ic,io,iargc,exitcode,tempintarray(99),getos,get_ran_seed,status,system
  real(double) :: timestamp,timestamps(9)  ! Time the progress of the code.
  character :: infile*99
  logical :: ex,timing
  
  
  outputdir = '.'     !Directory where output is saved (either relative or absolute path)
  wikioutput = 1      !Produce output for CBC Wiki: 0-no, 1-yes (requires one of the probability intervals to be 2-sigma)
  map_projection = 1  !Choose map projection: 1-Mollweide
  html = 0            !Produce HTML output
  
  call setconstants()         !Define mathematical constants
  os = getos()                !1-Linux, 2-MacOS
  timestamps(1) = timestamp()
  timing = .false.
  if(abs(timestamps(1)).gt.1.e-6_dbl .and. abs(timestamps(1)).lt.1.e20_dbl) timing = .true.
  
  
  call set_plotsettings()     !Set plot settings to 'default' values
  call read_settingsfile()    !Read the plot settings (overwrite the defaults)
  if(prProgress.ge.3) call write_settingsfile()   !Write the input file back to disc
  !call write_settingsfile()   !Write the input file back to disc
  
  if(html.ge.1) then
     outputdir = 'html'     !Directory where output is saved (either relative or absolute path)
     
     file = 1
     colour = 1
     
     prStdOut = 2
     prProgress = 2
     prRunInfo = 2
     prChainInfo = 2
     prInitial = 4
     
     prStat = 2
     prCorr = 1
     prIval = 3
     prConv = 3
     saveStats = 1
     
     plot = 1
     plLogL = 1
     plChain = 1
     plPDF1D = 1
     plPDF2D = 2
     plotSky = 2
     plAnim = 0
     
     plInject = 0
     plStart = 1
     plMedian = 1
     plRange = 4
     plBurn = 1
     plLmax = 1
     normPDF2D = 4
     
     bmpXSz = 1000
     bmpYSz =  700
     
     !nPlPar = 
     !plPars = (//)
     Npdf2D = -1
  end if
  
  
  
  inquire(file=trim(outputdir), exist=ex)  !Check whether the directory already exists 
  if(.not.ex) then
     status = system('mkdir -p '//trim(outputdir))
     if(status.ne.0) call quit_program('Could not create output directory: '//trim(outputdir)//'/')
  end if
  
  if(prStdOut.ge.2) then  !Write standard output to file rather than screen
     stdOut = 19
     write(stdOutFile,'(A,I6.6,A)')trim(outputdir)//'/analysemcmc_tempstdout_',abs(get_ran_seed(0)),'.txt'
     open(unit=stdOut,action='write',form='formatted',status='replace',file=trim(stdOutFile),iostat=io)
     if(io.ne.0) call quit_program('Error opening output file '//trim(stdOutFile))
  end if
  write(stdOut,*)
  
  if(html.ge.1) then  !Write standardised HTML output
     open(unit=51,action='write',form='formatted',status='replace',file=trim(outputdir)//'/index.html',iostat=io)
     if(io.ne.0) call quit_program('Error opening HTML output file '//trim(outputdir)//'/index.html')
  end if
  
  
  
  
  
  !Get command-line arguments
  nchains0 = iargc()
  if(nchains0.lt.1) then  !No command-line arguments - select all files SPINspiral.output.*.00 in the current dir
     call findFiles('SPINspiral.output.*.00',maxChs,1,infiles,nchains0)
     if(nchains0.eq.0) then
        write(stdErr,'(A)')'  No files matching  SPINspiral.output.*.00  were found in the current directory.'
        write(stdErr,'(A)')'  I will try the old file names  mcmc.output.*.00  instead.'
        call findFiles('mcmc.output.*.00',maxChs,1,infiles,nchains0)
        if(nchains0.eq.0) call quit_program('No valid input files were found in the current directory.'// &
             '  Please specify input files manually.')
     end if
  else
     do ic = 1,nchains0
        if(reverseRead.eq.0) then
           call getarg(ic,infile) !Read file name from the command-line arguments
        else
           call getarg(nchains0-ic+1,infile) !Read file name from the command-line arguments in reverse order
        end if
        infiles(ic) = infile
     end do
  end if
  
  if(nchains0.gt.maxChs) write(stdErr,'(A,I3,A)')'  *** WARNING:  Too many input files (chains),'// &
       ' please increase maxChs in analyseMCMC_modules.f90. Only',maxChs,' files can be read.'
  nchains0 = min(nchains0,maxChs)
  
  
  
  
  !Some of the stuff below will have to go to the input file
  
  !~Maximum number of dots to plot in e.g. chains plot, to prevent dots from being overplotted too much 
  !  and eps/pdf files from becoming huge.  Use this to autoset chainPlI
  maxdots = 25000  
  
  !Use a white background in screen and bitmap plots: 0-no (black), 1-yes.  Used to be in input file.
  whiteBG = 1                 
  
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
  write(unSharpchain,'(I4)')unSharp                       !~15 panels with dots
  write(unSharppdf1d,'(I4)')max(nint(real(unSharp)/2.),1) !~15 panels, no dots
  write(unSharppdf2d,'(I4)')max(nint(real(unSharp)/4.),1) !1 panel, no dots
  
  
  
  
  
  
  
  
  !Sort out implicit options:
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
     !if(nPlPar.ne.15) write(stdErr,'(/,A)')'*** WARNING:  I changed nPlPar to 15, since savePDF is selected ***'
     !nPlPar = 15; plPars(1:nPlPar) = (/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15/) !All 12 + m1,m2
     wrapData = 0
  end if
  if(file.ge.1) update = 0
  if(plAnim.ge.1) update = 0
  
  colournames(1:15) = (/ &
       'white               ','red                 ','dark green          ','dark blue           ','cyan                ', &
       'magenta             ','yellow              ','orange              ','light green         ','brown               ', &
       'dark red            ','purple              ','red-purple          ','dark grey           ','light grey          '/)
  if(file.ge.2) colournames(1) = 'black'
  
  call set_originalParameterNames()  !Set the names and symbols of the original MCMC parameters in the database
  
  
  if(prChainInfo.ge.1) then
     write(stdOut,'(A)',advance="no")'  AnalyseMCMC, svn revision: '//trim(svnrevision)//'.'
     if(nchains0.eq.1) then
        write(stdOut,'(A)')'  Analysing 1 chain from SPINspiral,'
     else
        write(stdOut,'(A,I3,A)')'  Analysing',nchains0,' chains from SPINspiral,'
     end if
     write(stdOut,'(A)',advance='no')'  in '//trim(username)//'@'//trim(hostname)//':'//trim(workdir)//'/'
     write(stdOut,'(A)')',  at '//trim(currentdatestr)//' '//trim(currenttimestr)//' ('//trim(currenttimezonestr)//').'
  end if
  nchains = nchains0
  
  
  
  
  
  
  
  !*********************************************************************************************************************************
  !***   READ INPUT FILE(S)   ******************************************************************************************************
  !*********************************************************************************************************************************
  
101 continue
  !Read the input files:
  exitcode = 0
  call read_mcmcfiles(exitcode)
  if(exitcode.ne.0) goto 9998
  
  !Get and print some basic chain statistics:
  if(timing) timestamps(2) = timestamp()
  exitcode = 0
  call mcmcruninfo(exitcode)
  
  if(savePDF.ge.2) then
       write(stdOut,'(A)')'  Writing after-burnin data points to file'
       call save_data(exitcode) !save after-burnin combined data to file
  end if

  
  !More implicit options:
  if(panels(1)*panels(2).lt.min(nPlPar,nMCMCpar)) panels = 0
  if(panels(1)*panels(2).lt.1) then
     if(min(nPlPar,nMCMCpar).eq.1) panels = (/1,1/)
     if(min(nPlPar,nMCMCpar).eq.2) panels = (/2,1/)
     if(min(nPlPar,nMCMCpar).eq.3) panels = (/3,1/)
     if(min(nPlPar,nMCMCpar).eq.4) panels = (/2,2/)
     if(min(nPlPar,nMCMCpar).eq.5) panels = (/5,1/)
     if(min(nPlPar,nMCMCpar).eq.6) panels = (/3,2/)
     if(min(nPlPar,nMCMCpar).eq.7) panels = (/4,2/)
     if(min(nPlPar,nMCMCpar).eq.8) panels = (/4,2/)
     if(min(nPlPar,nMCMCpar).eq.9) panels = (/3,3/)
     if(min(nPlPar,nMCMCpar).eq.10) panels = (/5,2/)
     if(min(nPlPar,nMCMCpar).eq.11) panels = (/4,3/)
     if(min(nPlPar,nMCMCpar).eq.12) panels = (/4,3/)
     if(min(nPlPar,nMCMCpar).eq.12.and.quality.eq.3) panels = (/3,4/)
     if(min(nPlPar,nMCMCpar).eq.13) panels = (/5,3/)
     if(min(nPlPar,nMCMCpar).eq.14) panels = (/5,3/)
     if(min(nPlPar,nMCMCpar).eq.15) panels = (/5,3/)
     if(min(nPlPar,nMCMCpar).eq.16) panels = (/4,4/)
     if(min(nPlPar,nMCMCpar).eq.17) panels = (/6,3/)
     if(min(nPlPar,nMCMCpar).eq.18) panels = (/6,3/)
     if(min(nPlPar,nMCMCpar).eq.19) panels = (/5,4/)
     if(min(nPlPar,nMCMCpar).eq.20) panels = (/5,4/)
  end if
  
  
  
  
  ! ********************************************************************************************************************************
  ! ***  DO STATISTICS   ***********************************************************************************************************
  ! ********************************************************************************************************************************
  
  if(timing) timestamps(3) = timestamp()
  
  exitcode = 0
  call statistics(exitcode)
  if(exitcode.ne.0) goto 9999
  
  
  
  
  
  ! ********************************************************************************************************************************
  ! ***  CREATE PLOTS   ************************************************************************************************************
  ! ********************************************************************************************************************************
  
  if(timing) timestamps(4) = timestamp()
  
  if(prProgress.ge.2) write(stdOut,*)
  if(plot.eq.1.and.prProgress.ge.1.and.update.eq.0) then
     write(stdOut,'(/,A)',advance="no")'  Plotting '
     if(file.eq.0) write(stdOut,'(A)',advance="no")'to screen: '
     if(file.eq.1) write(stdOut,'(A)',advance="no")'to png: '
     if(file.eq.2) write(stdOut,'(A)',advance="no")'to eps: '
     if(file.eq.3) write(stdOut,'(A)',advance="no")'to pdf: '
  end if
  
  
  
  !*********************************************************************************************************************************
  !Plot (1d) chains: logL, parameter chains, jumps, etc.
  if(plot.eq.1) then
     exitcode = 0
     call chains(exitcode)
     if(exitcode.ne.0) goto 9999
  end if
  if(timing) timestamps(5) = timestamp()
  
  
  
  
  !*********************************************************************************************************************************
  !Plot pdfs (1d)
  if(plPDF1D.ge.1) then
     exitcode = 0
     call pdfs1d(exitcode)
     if(exitcode.ne.0) goto 9999
  end if !if(plPDF1D.ge.1)
  if(timing) timestamps(6) = timestamp()
  
  
  
  
  !*********************************************************************************************************************************
  if(plPDF2D.ge.1.and.mergeChains.eq.0) then
     write(stdOut,'(A)',advance="no")', (skipping 2D PDFs since mergeChains=0), '
     plPDF2D = 0
  end if
  
  if(plPDF2D.ge.1) then
     exitcode = 0
     call pdfs2d(exitcode)
     if(exitcode.ne.0) goto 9999
  end if !if(plPDF2D.eq.1)
  
  if(Npdf2D.lt.0) then !Then we just plotted all 2D PDFs
     write(stdOut,*)
  else
     if(prProgress.ge.1.and.update.eq.0.and.plot.gt.0) write(stdOut,'(A,/)')' done.  '
  end if
  
  if(timing) timestamps(7) = timestamp()  
  
  
  
  
  !*********************************************************************************************************************************
  
  if(saveStats.ge.1.and.nchains.gt.1) then
     write(stdErr,'(A)')' ******   Cannot write statistics if the number of chains is greater than one   ******'
     
     !Write Bayes factors to file:
     exitcode = 0
     call save_bayes(exitcode)
     if(exitcode.ne.0) goto 9999
  end if
  
  !Write statistics to file:
  if(saveStats.ge.1.and.nchains.eq.1) then
     exitcode = 0
     call save_stats(exitcode)
     if(exitcode.ne.0) goto 9999
     write(stdOut,*)
  end if !if(saveStats.ge.1.and.nchains.eq.1) then
  
  
  !Save tailored output to a file:
  if(tailoredOutput.gt.0) then
     exitcode = 0
     call tailored_output(exitcode)
     !if(exitcode.ne.0) goto 9999
  end if
  
  !*********************************************************************************************************************************
  
  if(timing) timestamps(8) = timestamp()
  
  if(plAnim.ge.1) then
     exitcode = 0
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
  
  !write(stdOut,'(A)')'  Waiting for you to finish me off...'
  !pause
  
9999 continue
  deallocate(selDat)
9998 continue
  deallocate(allDat,post,prior)
  !if(prProgress.ge.1) write(stdOut,*)
  
  if(timing) then
     timestamps(9) = timestamp()
  
     if(prProgress.ge.1) then
        write(stdOut,'(A)',advance="no")'  Run time: '
        write(stdOut,'(A,F5.1,A)',advance="no")'   input:',min(abs(timestamps(2)-timestamps(1)),999.9_dbl),'s,'
        !write(stdOut,'(A,F5.1,A)',advance="no")'   info:',min(abs(timestamps(3)-timestamps(2)),999.9_dbl),'s,'
        !write(stdOut,'(A,F5.1,A)',advance="no")'   stats:',min(abs(timestamps(4)-timestamps(3)),999.9_dbl),'s,'
        write(stdOut,'(A,F5.1,A)',advance="no")'   stats:',min(abs(timestamps(4)-timestamps(2)),999.9_dbl),'s,'
        if(plot.eq.1.and.plLogL+plChain+plJump+plACorr.gt.0) then
           write(stdOut,'(A,F5.1,A)',advance="no")'   chains:',min(abs(timestamps(5)-timestamps(4)),999.9_dbl),'s,'
        end if
        if(plot.eq.1.or.savePDF.ge.1) then
           if(plPDF1D.ge.1) write(stdOut,'(A,F5.1,A)',advance="no") &
                '   1d pdfs:',min(abs(timestamps(6)-timestamps(5)),999.9_dbl),'s,'
           if(plPDF2D.ge.1) write(stdOut,'(A,F6.1,A)',advance="no") &
                '   2d pdfs:',min(abs(timestamps(7)-timestamps(6)),999.9_dbl),'s,'
        end if
        !write(stdOut,'(A,F6.1,A)',advance="no")'   plots:',min(abs(timestamps(7)-timestamps(4)),999.9_dbl),'s,'
        !write(stdOut,'(A,F5.1,A)',advance="no")'   save stats:',min(abs(timestamps(8)-timestamps(7)),999.9_dbl),'s,'
        if(plAnim.ge.1) write(stdOut,'(A,F5.1,A)',advance="no")'   movie:',min(abs(timestamps(9)-timestamps(8)),999.9_dbl),'s,'
        write(stdOut,'(A,F6.1,A)')'   total:',min(abs(timestamps(9)-timestamps(1)),999.9_dbl),'s.'
     end if
  end if
  
  
  if(html.ge.1) close(51)
  
  write(stdOut,*)
  
  !Save standard-output file under final name:
  if(prStdOut.ge.2) then
     status = system('mv -f '//trim(stdOutFile)//' '//trim(outputdir)//'/'//trim(outputname)//'__output.txt')
     if(status.eq.0) then
        !Should be 6, not stdOut:
        write(6,'(/,A,/)')'  AnalyseMCMC:  saved standard output to '//trim(outputdir)//'/'//trim(outputname)//'__output.txt'
     else
        write(stdErr,'(/,A)')'  AnalyseMCMC:  Error saving standard output to '// &
             trim(outputdir)//'/'//trim(outputname)//'__output.txt'
        !write(stdErr,'(A,/)')'                Check the file '//trim(stdOutFile)
        status = system('rm -f '//trim(stdOutFile))
        write(stdErr,*)
     end if
  end if
  
end program analyseMCMC
!***********************************************************************************************************************************





