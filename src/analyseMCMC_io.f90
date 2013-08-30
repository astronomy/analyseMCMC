!> \file analyseMCMC_io.f90  Input/output routines for AnalyseMCMC

! 
! LICENCE:
! 
! Copyright 2007-2013 Marc van der Sluys
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
!> \brief  Define the constants

subroutine setconstants()
  use SUFR_constants, only: stdErr,stdOut
  use aM_constants, only: detabbrs,waveforms, rmeps
  
  implicit none
  
  stdErr = 0  ! Standard error unit
  stdOut = 6  ! Standard output unit - screen
  
  ! Define detector abbreviations here (don't forget to change detabbrs() in the module constants in analyseMCMC_modules.f90):
  detabbrs = (/'H1','L1','V ','H2'/)
  
  ! Define waveforms here (don't forget to change waveforms() in the module constants in analyseMCMC_modules.f90):
  ! In LAL, these numbers change so rapidly that its pointless to specify them here
  waveforms    = 'Unknown'             ! All of them, except those specified below
  waveforms(1) = 'Apostolatos'
  waveforms(2) = 'SpinTaylor12'
  waveforms(3) = 'SpinTaylor15'
  waveforms(4) = 'PPN'
  waveforms(5) = 'PhenSpinInspiralRD'
  waveforms(9) = 'Ana.L'
  
  rmeps = .true.  ! Remove eps files after conversion to pdf
  
end subroutine setconstants
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Read the settings file (called analysemcmc.dat by default) - old version (< 2012-06)
!!
!! \retval status  Status: 0-ok

subroutine read_settingsfile_old(status)
  use SUFR_kinds, only: double
  use SUFR_system, only: quit_program, find_free_io_unit
  
  use analysemcmc_settings, only: settingsfile
  use analysemcmc_settings, only: Nburn,ivals,plPars,panels,PDF2Dpairs,thin,NburnFrac,autoBurnin,maxChs,maxChLen,file,colour
  use analysemcmc_settings, only: quality,reverseRead,update,mergeChains,wrapData,changeVar,prStdOut,prProgress,prRunInfo
  use analysemcmc_settings, only: prChainInfo,prInitial,prStat,prCorr,prAcorr,nAcorr,prIval,prConv,saveStats,savePDF,tailoredOutput
  use analysemcmc_settings, only: plot,plLogL,plChain,plParL,plJump,plPDF1D,plPDF2D,plAcorr,plotSky,plAnim,chainPlI,scLogLpl
  use analysemcmc_settings, only: scChainsPl,plInject,plStart,plMedian,plRange,plBurn,plLmax,prValues,smooth,fillPDF,normPDF1D
  use analysemcmc_settings, only: normPDF2D,nAnimFrames,animScheme,Nival,ival0,scrSz,scrRat,bmpXSz,bmpYSz,PSsz,PSrat,scFac,unSharp
  use analysemcmc_settings, only: orientation,fontType,fontSize1D,fontSize2D,chainSymbol,nPlPar,Nbin1D,Nbin2Dx,Nbin2Dy,Npdf2D
  
  implicit none
  integer, intent(out) :: status
  
  integer :: i,ip,io,io1
  character :: bla
  real(double) :: dblvar
  
  status = 0
  ! dblvar is used when a (possibly) large integer is expected; read it as double, then convert to integer
  
  call find_free_io_unit(ip)
  open(unit=ip,form='formatted',status='old',action='read',file=trim(settingsfile), iostat=io)
  if(io.ne.0) call quit_program('Error opening input file '//trim(settingsfile))
  
  io = 0
  io1 = 0
  
  read(ip,*,iostat=io) bla
  bla = bla  ! Remove 'set but never used' warning
  
  
  read(ip,*,iostat=io) bla
  read(ip,*,iostat=io) thin
  read(ip,*,iostat=io) dblvar
  Nburn(1) = nint(dblvar)
  do i=2,maxChs
     Nburn(i) = Nburn(1)
  end do
  read(ip,*,iostat=io) NburnFrac
  read(ip,*,iostat=io) autoBurnin
  read(ip,*,iostat=io) dblvar
  maxChLen = nint(dblvar)
  read(ip,*,iostat=io) file
  read(ip,*,iostat=io) colour
  read(ip,*,iostat=io) quality
  read(ip,*,iostat=io) reverseRead
  read(ip,*,iostat=io) update
  read(ip,*,iostat=io) mergeChains
  read(ip,*,iostat=io) wrapData
  read(ip,*,iostat=io) changeVar
  
  read(ip,*,iostat=io) bla
  read(ip,*,iostat=io) prStdOut
  read(ip,*,iostat=io) prProgress
  read(ip,*,iostat=io) prRunInfo
  read(ip,*,iostat=io) prChainInfo
  read(ip,*,iostat=io) prInitial
  read(ip,*,iostat=io) prStat
  read(ip,*,iostat=io) prCorr
  read(ip,*,iostat=io) prAcorr
  read(ip,*,iostat=io) dblvar
  nAcorr = nint(dblvar)
  read(ip,*,iostat=io) prIval
  read(ip,*,iostat=io) prConv
  read(ip,*,iostat=io) saveStats
  read(ip,*,iostat=io) savePDF
  read(ip,*,iostat=io) tailoredOutput
  
  
  read(ip,*,iostat=io) bla
  read(ip,*,iostat=io) plot
  read(ip,*,iostat=io) plLogL
  read(ip,*,iostat=io) plChain
  read(ip,*,iostat=io) plParL
  read(ip,*,iostat=io) plJump
  read(ip,*,iostat=io) plPDF1D
  read(ip,*,iostat=io) plPDF2D
  read(ip,*,iostat=io) plAcorr
  read(ip,*,iostat=io) plotSky
  read(ip,*,iostat=io) plAnim
  
  
  read(ip,*,iostat=io) bla
  read(ip,*,iostat=io) chainPlI
  read(ip,*,iostat=io) scLogLpl
  read(ip,*,iostat=io) scChainsPl
  read(ip,*,iostat=io) plInject
  read(ip,*,iostat=io) plStart
  read(ip,*,iostat=io) plMedian
  read(ip,*,iostat=io) plRange
  read(ip,*,iostat=io) plBurn
  read(ip,*,iostat=io) plLmax
  read(ip,*,iostat=io) prValues
  read(ip,*,iostat=io) smooth
  read(ip,*,iostat=io) fillPDF
  read(ip,*,iostat=io) normPDF1D
  read(ip,*,iostat=io) normPDF2D
  read(ip,*,iostat=io) nAnimFrames
  read(ip,*,iostat=io) animScheme
  read(ip,*,iostat=io) Nival,ival0
  read(ip,*,iostat=io1)(ivals(i),i=1,Nival)
  
  read(ip,*,iostat=io) bla
  read(ip,*,iostat=io) scrSz
  read(ip,*,iostat=io) scrRat
  read(ip,*,iostat=io) bmpXSz
  read(ip,*,iostat=io) bmpYSz
  read(ip,*,iostat=io) PSsz
  read(ip,*,iostat=io) PSrat
  read(ip,*,iostat=io) scFac
  read(ip,*,iostat=io) unSharp
  
  read(ip,*,iostat=io) bla
  read(ip,*,iostat=io) orientation
  read(ip,*,iostat=io) fontType
  read(ip,*,iostat=io) fontSize1D
  read(ip,*,iostat=io) fontSize2D
  read(ip,*,iostat=io) chainSymbol
  
  read(ip,*,iostat=io) bla
  read(ip,*,iostat=io) bla
  read(ip,*,iostat=io) nPlPar
  read(ip,*,iostat=io1)(plPars(i),i=1,nPlPar)
  if(io1.ne.0) nPlPar = i-1
  io1 = 0
  read(ip,*,iostat=io) panels(1:2)
  read(ip,*,iostat=io) Nbin1D
  read(ip,*,iostat=io) Nbin2Dx
  read(ip,*,iostat=io) Nbin2Dy
  read(ip,*,iostat=io) Npdf2D
  do i=1,Npdf2D
     read(ip,*,iostat=io1)PDF2Dpairs(i,1:2)
     if(io1.ne.0) exit
  end do
  if(io1.ne.0) Npdf2D = i-1
  close(ip)
  
  if(io.ne.0) then
     status = 1
     return
  end if
  
end subroutine read_settingsfile_old
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Read the settings file (called analysemcmc.dat by default)

subroutine read_settingsfile()
  use SUFR_system, only: find_free_io_unit, quit_program_error
  
  use analysemcmc_settings, only: settingsfile
  use analysemcmc_settings, only: Nburn,ivals,plPars,panels,PDF2Dpairs,thin,NburnFrac,autoBurnin,maxChLen,file,colour
  use analysemcmc_settings, only: quality,reverseRead,update,mergeChains,wrapData,changeVar,prStdOut,prProgress,prRunInfo
  use analysemcmc_settings, only: prChainInfo,prInitial,prStat,prCorr,prAcorr,nAcorr,prIval,prConv,saveStats,savePDF,tailoredOutput
  use analysemcmc_settings, only: plot,plLogL,plChain,plParL,plJump,plPDF1D,plPDF2D,plAcorr,plRhat,plotSky,plAnim,chainPlI,scLogLpl
  use analysemcmc_settings, only: scChainsPl,plInject,plStart,plMedian,plRange,plBurn,plLmax,prValues,smooth,fillPDF,normPDF1D
  use analysemcmc_settings, only: normPDF2D,nAnimFrames,animScheme,Nival,ival0,scrSz,scrRat,bmpXSz,bmpYSz,PSsz,PSrat,scFac,unSharp
  use analysemcmc_settings, only: orientation,fontType,fontSize1D,fontSize2D,chainSymbol,nPlPar,Nbin1D,Nbin2Dx,Nbin2Dy,Npdf2D
  use analysemcmc_settings, only: mapProjection, wikiOutput, htmlOutput
  use general_data, only: outputDir
  
  implicit none
  integer :: ip, io, NburnMax
  
  ! Basic options:
  namelist /basic_options/ thin, NburnMax, NburnFrac, autoBurnin, maxChLen, file, colour, quality, reverseRead, &
       update, mergeChains, wrapData, changeVar, outputDir
  
  ! Select what output to print to screen and write to file:
  namelist /print_options/ prStdOut, prProgress, prRunInfo, prChainInfo, prInitial, prStat, prCorr, prAcorr, nAcorr, prIval, &
       prConv, saveStats, savePDF, wikiOutput, tailoredOutput, htmlOutput
  
  ! Select which plots to make:
  namelist /plot_select/ plot, plLogL, plChain, plParL, plJump, plPDF1D, plPDF2D, plAcorr, plRhat, plotSky, mapProjection, plAnim
  
  ! Detailed plot settings:
  namelist /plot_options/ chainPlI, scLogLpl, scChainsPl, plInject, plStart, plMedian, plRange, plBurn, plLmax, prValues, smooth, &
       fillPDF, normPDF1D, normPDF2D, nAnimFrames, animScheme, Nival, ival0, ivals  ! was: ivals(1:Nival)
  
  ! Output format:
  namelist /output_format/ scrSz, scrRat, bmpXSz, bmpYSz, PSsz, PSrat, scFac, unSharp
  
  ! Fonts, symbols, etc.:
  namelist /fonts_symbols/ orientation, fontType, fontSize1D, fontSize2D, chainSymbol
  
  ! Data settings:
  namelist /plot_parameters_binning/ nPlPar, plPars, panels, Nbin1D, Nbin2Dx, Nbin2Dy, Npdf2D, PDF2Dpairs
  ! was: plPars(1:nPlPar), panels(1:2), PDF2Dpairs(1:Npdf2D,1:2)
  
  
  
  call find_free_io_unit(ip)
  write(*,'(A)') '  Reading settings file '//trim(settingsfile)//'...'
  open(unit=ip, status='old', action='read', file=trim(settingsfile), iostat=io)
  if(io.ne.0) call quit_program_error('readsettingsfile(): error opening settings file '//trim(settingsfile), 0)
  
  read(ip, nml=basic_options, iostat=io)             ! Basic options
  if(io.ne.0) then
     call try_old_settings_file(ip)
     call quit_program_error('readsettingsfile(): error reading settings file '//trim(settingsfile)//', Basic options', 0)
  end if
  
  read(ip, nml=print_options, iostat=io)             ! Print options
  if(io.ne.0) then
     call try_old_settings_file(ip)
     call quit_program_error('readsettingsfile(): error reading settings file '//trim(settingsfile)//', Print options', 0)
  end if
  
  read(ip, nml=plot_select, iostat=io)               ! Select which plots to make
  if(io.ne.0) then
     call try_old_settings_file(ip)
     call quit_program_error('readsettingsfile(): error reading settings file '//trim(settingsfile)//', Plot select', 0)
  end if
  
  read(ip, nml=plot_options, iostat=io)              ! Detailed plot settings
  if(io.ne.0) then
     call try_old_settings_file(ip)
     call quit_program_error('readsettingsfile(): error reading settings file '//trim(settingsfile)//', Plot options', 0)
  end if
  
  read(ip, nml=output_format, iostat=io)             ! Output format
  if(io.ne.0) then
     call try_old_settings_file(ip)
     call quit_program_error('readsettingsfile(): error reading settings file '//trim(settingsfile)//', Output format', 0)
  end if
  
  read(ip, nml=fonts_symbols, iostat=io)             ! Fonts, symbols, etc.
  if(io.ne.0) then
     call try_old_settings_file(ip)
     call quit_program_error('readsettingsfile(): error reading settings file '//trim(settingsfile)//', Fonts, symbols, etc.', 0)
  end if
  
  read(ip, nml=plot_parameters_binning, iostat=io)   ! Select parameters to plot, binning, etc.
  if(io.ne.0) then
     call try_old_settings_file(ip)
     call quit_program_error('readsettingsfile(): error reading settings file '//trim(settingsfile)// &
          ', Plot parameters, binning, etc.', 0)
  end if
  
  close(ip)
  
  Nburn = NburnMax
  if(Npdf2D.lt.0) plotsky = 0  ! Do not plot a (full) sky for the RA-dec 2D PDF when plotting all 2D PDFs
  
end subroutine read_settingsfile
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Reading the settings file in the new (2012-06) namelist format didn't work, try the old format and save the data in the
!!         new format
!!
!! \brief old_ip  Input unit currently open with the input file

subroutine try_old_settings_file(old_ip)
  use SUFR_constants, only: stdOut
  implicit none
  integer, intent(in) :: old_ip
  integer :: status
  
  close(old_ip)
  call read_settingsfile_old(status)
  if(status.ne.0) return
  
  call write_settingsfile()
  
  write(stdOut,*) ''
  write(stdOut,'(A)') '  Your input file could not be read using the current format.  I tried to read it using the old (< 2012-06)'
  write(stdOut,'(A)') '  format and saved the input file in the new format, as analysemcmc.new.  Please consider comparing that'
  write(stdOut,'(A)') '  file to the analysemcmc.dat file in the doc/ directory of AnalyseMCMC.  The latter has comments, which'
  write(stdOut,'(A)') '  makes it much easier to use the file.  Then rename the file to analysemcmc.dat and try running the code'
  write(stdOut,'(A)') '  again.'
  write(stdOut,*) ''
  
  stop
  
end subroutine try_old_settings_file
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Write a copy of the settings file (called analysemcmc.new by default)

subroutine write_settingsfile()
  use SUFR_system, only: find_free_io_unit
  use SUFR_constants, only: stdOut
  
  use analysemcmc_settings, only: Nburn,ivals,plPars,panels,PDF2Dpairs,thin,NburnFrac,autoBurnin,maxChLen,file,colour
  use analysemcmc_settings, only: quality,reverseRead,update,mergeChains,wrapData,changeVar,prStdOut,prProgress,prRunInfo
  use analysemcmc_settings, only: prChainInfo,prInitial,prStat,prCorr,prAcorr,nAcorr,prIval,prConv,saveStats,savePDF,tailoredOutput
  use analysemcmc_settings, only: plot,plLogL,plChain,plParL,plJump,plPDF1D,plPDF2D,plAcorr,plRhat,plotSky,plAnim,chainPlI,scLogLpl
  use analysemcmc_settings, only: scChainsPl,plInject,plStart,plMedian,plRange,plBurn,plLmax,prValues,smooth,fillPDF,normPDF1D
  use analysemcmc_settings, only: normPDF2D,nAnimFrames,animScheme,Nival,ival0,scrSz,scrRat,bmpXSz,bmpYSz,PSsz,PSrat,scFac,unSharp
  use analysemcmc_settings, only: orientation,fontType,fontSize1D,fontSize2D,chainSymbol,nPlPar,Nbin1D,Nbin2Dx,Nbin2Dy,Npdf2D
  use analysemcmc_settings, only: mapProjection, wikiOutput, htmlOutput
  use general_data, only: outputDir
  
  implicit none
  integer :: op, NburnMax
  
  ! Basic options:
  namelist /basic_options/ thin, NburnMax, NburnFrac, autoBurnin, maxChLen, file, colour, quality, reverseRead, &
       update, mergeChains, wrapData, changeVar, outputDir
  
  ! Select what output to print to screen and write to file:
  namelist /print_options/ prStdOut, prProgress, prRunInfo, prChainInfo, prInitial, prStat, prCorr, prAcorr, nAcorr, prIval, &
       prConv, saveStats, savePDF, wikiOutput, tailoredOutput, htmlOutput
  
  ! Select which plots to make:
  namelist /plot_select/ plot, plLogL, plChain, plParL, plJump, plAcorr, plRhat, plPDF1D, plPDF2D, plotSky, mapProjection, plAnim
  
  ! Detailed plot settings:
  namelist /plot_options/ chainPlI, scLogLpl, scChainsPl, plInject, plStart, plMedian, plRange, plBurn, plLmax, prValues, smooth, &
       fillPDF, normPDF1D, normPDF2D, nAnimFrames, animScheme, Nival, ival0, ivals  ! was: ivals(1:Nival)
  
  ! Output format:
  namelist /output_format/ scrSz, scrRat, bmpXSz, bmpYSz, PSsz, PSrat, scFac, unSharp
  
  ! Fonts, symbols, etc.:
  namelist /fonts_symbols/ orientation, fontType, fontSize1D, fontSize2D, chainSymbol
  
  ! Data settings:
  namelist /plot_parameters_binning/ nPlPar, plPars, panels, Nbin1D, Nbin2Dx, Nbin2Dy, Npdf2D, PDF2Dpairs
  ! was: plPars(1:nPlPar), panels(1:2), PDF2Dpairs(1:Npdf2D,1:2)
  
  
  NburnMax = maxval(Nburn)
  
  
  call find_free_io_unit(op)
  write(stdOut,'(A)') '  Saving settings file as analysemcmc.new'
  open(unit=op,form='formatted',status='replace',action='write',file='analysemcmc.new')
  
  write(op,'(A,/)') '# Input file for AnalyseMCMC'
  
  write(op, nml=basic_options)             ! Basic options
  write(op, nml=print_options)             ! Print options
  write(op, nml=plot_select)               ! Select which plots to select
  write(op, nml=plot_options)              ! Detailed plot settings
  write(op, nml=output_format)             ! Output format
  write(op, nml=fonts_symbols)             ! Fonts, symbols, etc.
  write(op, nml=plot_parameters_binning)   ! Select parameters to plot, binning, etc.
  
  close(op)
  
end subroutine write_settingsfile
!***********************************************************************************************************************************




!***********************************************************************************************************************************
!> \brief  Set plot settings to 'default' values

subroutine set_plotsettings()
  use analysemcmc_settings, only: Nburn,ivals,plPars,panels,PDF2Dpairs,thin,NburnFrac,autoBurnin,maxChLen,file,colour
  use analysemcmc_settings, only: quality,reverseRead,update,mergeChains,wrapData,changeVar,prStdOut,prProgress,prRunInfo
  use analysemcmc_settings, only: prChainInfo,prInitial,prStat,prCorr,prAcorr,nAcorr,prIval,prConv,saveStats,savePDF,tailoredOutput
  use analysemcmc_settings, only: plot,plLogL,plChain,plParL,plJump,plPDF1D,plPDF2D,plAcorr,plRhat,plotSky,plAnim,chainPlI,scLogLpl
  use analysemcmc_settings, only: scChainsPl,plInject,plStart,plMedian,plRange,plBurn,plLmax,prValues,smooth,fillPDF,normPDF1D
  use analysemcmc_settings, only: normPDF2D,nAnimFrames,animScheme,Nival,ival0,scrSz,scrRat,bmpXSz,bmpYSz,PSsz,PSrat,scFac,unSharp
  use analysemcmc_settings, only: orientation,fontType,fontSize1D,fontSize2D,chainSymbol,nPlPar,Nbin1D,Nbin2Dx,Nbin2Dy,Npdf2D
  use analysemcmc_settings, only: htmlOutput, mapProjection, wikiOutput
  use general_data, only: outputDir
  
  implicit none
  
  
  ! Basic options:
  thin = 10         ! If >1, 'thin' the output; read every thin-th line 
  Nburn = nint(5.e5) ! If >=0: override length of the burn-in phase, for all chains
  NburnFrac = 0.5   ! If !=0: override length of the burn-in phase, as a fraction of the length of each chain.
  autoBurnin = -1.  ! Determine burn-in automatically as the first iteration where log(L_chain) > max(log(L_allchains)) - autoBurnin
  maxChLen = nint(1.e9) ! Maximum chain length
  
  file = 1          ! Plot output to file:  0-no; screen,  >0-yes; 1-png, 2-eps, 3-pdf.
  colour = 1        ! Use colours: 0-no (grey scales), 1-yes
  quality = 2       ! 'Quality' of plot, depending on purpose: 0: draft, 1: paper, 2: talk, 3: poster
  
  reverseRead = 0   ! Read files reversely (anti-alphabetically), to plot coolest chain last so that it becomes better visible
  update = 0        ! Update screen plot every 10 seconds: 0-no, 1-yes
  mergeChains = 1   ! Merge the data from different files into one chain: 0-no (treat separately), 1-yes (default)
  wrapData = 1      ! Wrap the data for the parameters that are in [0,2pi]: 0-no, 1-yes (useful if the peak is around 0)
  changeVar = 1     ! Change MCMC parameters (e.g. logd->d, kappa->theta_SL, rad->deg), 2=yes + q->1/q, phi->phi-pi, m1<->m2
  outputDir = "."   ! Save output files/plots in this directory - this may be a relative or absolute path
  
  
  ! Select what output to print to screen and write to file:
  prStdOut = 1      ! Print standard output to 1: screen, 2: text file
  prProgress = 2    ! Print general messages about the progress of the program: 0-no, 1-some, 2-more
  prRunInfo = 1     ! Print run info at read (# iterations, seed, # detectors, SNRs, data length, etc.): 0-no, 1-only for one file
  prChainInfo = 1   ! Print chain info: 1-summary (#datpts, #contr.chns),  2-detls/chain (f.name, clr, #itr, b.in, Lmax, #datpts)
  prInitial = 3     ! Print injection values, starting values and their difference: 0-no, 1-yes, 2-+ injection, 3-+ Lmax, 4-+ diffs
  
  prStat = 1        ! Print statistics: 0-no, 1-yes
  prCorr = 1        ! Print correlations: 0-no, 1-yes
  prAcorr = 0       ! Plot autocorrelations: 0-no, 1-some, 2: more
  nAcorr = 100      ! Compute prAcorr steps of autocorrelations if prAcorr>0 or plAcor>0 (default: 100)
  prIval = 1        ! Print interval info: 0-no, 1-yes
  prConv = 1        ! Print convergence information for multiple chains to screen and chains plot
  
  saveStats = 1     ! Save statistics (statistics, correlations, intervals) to file: 0-no, 1-yes, 2-yes + copy in PS
  savePDF = 0       ! Save the binned data for 1d and/or 2d pdfs (depending on plPDF1D and plPDF2D).
  wikiOutput = 0    ! Save output for the CBC wiki
  tailoredOutput = 0 ! Save output for a specific purpose, e.g. table in a paper
  htmlOutput = 0    ! Save HTML output - experimental, partly implemented.  Useful when plotting all 2D PDFs as png
  
  
  ! Select which plots to make:
  plot = 1          ! 0: plot nothing at all, 1: plot the items selected below
  
  plLogL = 1        ! Plot log L chains: 0-no, 1-yes
  plChain = 1       ! Plot parameter chains: 0-no, 1-yes
  plParL = 0        ! Plot L vs. parameter value: 0-no, 1-yes
  plJump = 0        ! Plot actual jump sizes 
  plAcorr = 0       ! Plot autocorrelations: 0-no, 1-yes
  plRhat = 0        ! Plot Rhat: 0-no, 1-yes
 
  plPDF1D = 1       ! Plot 1d posterior distributions
  plPDF2D = 2       ! Plot 2d posterior distributions
  plotSky = 2       ! Plot 2d pdf with stars, implies plPDF2D>0
  mapProjection = 1 ! Choose map projection: 1-Mollweide
  
  plAnim = 0        ! Plot movie frames
  
  
  ! Detailed plot settings:
  chainPlI = 0      ! Plot every chainPlI-th point in chains, logL, jump plots
  scLogLpl = 1      ! Scale logL plot ranges: 0
  scChainsPl = 1    ! Scale chains plot ranges
  
  plInject = 0      ! Plot injection values in the chains and pdfs
  plStart = 1       ! Plot starting values in the chains and pdfs
  plMedian = 1      ! Plot median values in the pdfs: 1-1D PDFs, 2-2D PDFs, 3-both
  plRange = 4       ! Plot the probability range in the pdfs: 1-1D PDFs, 2-2D PDFs, 3-both. 4-6: as 1-3 + write value in PDF panel
  plBurn = 1        ! Plot the burn-in in logL, the chains, etc.: 0-no, 1-vertical line, 2-colour shade, 3-both
  plLmax = 1        ! Plot the position of the max of logL in chains and pdfs
  
  prValues = 1      ! Print values (injection, median, range) in pdfs
  smooth = 3        ! Smooth the pdfs: 0 - no, >1: smooth over smooth bins (use ~10 (3-15)?)
  fillPDF = 1       ! Fillstyle for the pdfs (pgsfs): 1-solid, 2-outline, 3-hatched, 4-cross-hatched
  normPDF1D = 1     ! Normalise 1D pdfs
  normPDF2D = 4     ! 'Normalise' 2D pdfs: 0-linearly, 1-log, 2-sqrt, 3-weigted with likelihood value,  4-2D probability intervals
  
  nAnimFrames = 1   ! Number of frames for the movie
  animScheme = 3    ! Movie scheme: determines what panels to show in a movie frame 
  
  Nival = 3         ! Number of probability intervals
  ival0 = 2         ! Standard probability interval, e.g. 1 or 2, < Nival
  ivals = (/0.68269, 0.9545, 0.9973, 0., 0./)  ! Probability intervals - 1,2,3 sigma
  
  
  ! Output format:
  scrSz  = 10.8     ! Screen size for X11 windows (PGPlot units):  MacOS: 16.4, Gentoo: 10.8
  scrRat = 0.57     ! Screen ratio for X11 windows (PGPlot units), MacBook: 0.57
  bmpXSz = 1000     ! X-size for bitmap (pixels):  1000
  bmpYSz = 700      ! Y-size for bitmap (pixels):  700
  PSsz   = 10.5     ! Size for PS/PDF (PGPlot units).  Default: 10.5   \__ Gives same result as without pgpap
  PSrat  = 0.742    ! Ratio for PS/PDF (PGPlot units). Default: 0.742  /
  
  scFac = 1.2       ! Scale .png plots up by this factor, then down to the x,y size indicated above
  unSharp = 10      ! Apply unsharp mask when creating .png plots. Default: 10
  
  
  ! Fonts, symbols, etc.:
  orientation = 1   ! Use portrait (1) or landscape (2) for eps/pdf; mainly useful when sending a plot to a printer
  fonttype = 1      ! Font type used for eps/pdf: 1-simple, 2-roman, 3-italic, 4-script
  fontsize1d = 1.   ! Set plot scale for 1D plots, needs to be implemented fully
  fontsize2d = 1.   ! Set plot scale for 2D plots, needs to be implemented fully. Typically, take ~1.2*fontsize1d
  chainSymbol = 1   ! Plot symbol for the chains
  
  
  ! Select parameters to plot, binning, etc.:
  nPlPar = 9        ! Number of plot parameters for 1D plots
  plPars(1:nPlPar) = (/61,62, 22,11, 41, 31,32, 51,52/) ! The nPlPar plot parameters
  panels(1:2) = (/0,0/) ! Number of panels for 1D plots in x,y direction
  
  Nbin1D = 100      ! Number of bins for 1D PDFs
  Nbin2Dx = 60      ! Number of bins in horizontal direction for 2D PDFs
  Nbin2Dy = 40      ! Number of bins in vertical direction for 2D PDFs
  
  Npdf2D  = 1       ! Number of 2D PDFs to make
  PDF2Dpairs(1,1:2) = (/31,32/)  ! 2D PDFs to plot: RA,Dec
  
end subroutine set_plotsettings
!***********************************************************************************************************************************



!***********************************************************************************************************************************
!> \brief  Read the MCMC output files (SPINspiral.output.*, PTMCMC.output.*)
!!
!! \retval exitcode  Exit status code (0=ok)

subroutine read_mcmcfiles(exitcode)
  use SUFR_kinds, only: double
  use SUFR_constants, only: stdErr,stdOut
  use SUFR_system, only: quit_program_error, warn
  
  use analysemcmc_settings, only: thin,maxChLen, maxMCMCpar, prProgress
  use general_data, only: allDat, startval, post,prior, ntot,n,nchains,nchains0,infiles,maxIter, parNames
  use mcmcrun_data, only: niter,Nburn0,detnames,detnr,parID,seed,snr,revID,ndet,flow,fhigh,t_before,nCorr,nTemps,Tmax,Tchain
  use mcmcrun_data, only: networkSNR,waveform,waveformName,pnOrder,nMCMCpar,t_after,FTstart,deltaFT,samplerate,samplesize,FTsize
  use mcmcrun_data, only: nMCMCpar0,t0,GPStime,outputVersion
  use chain_data, only: is,DoverD
  use aM_constants, only: waveforms
  
  implicit none
  integer, intent(out) :: exitcode
  integer :: i,tmpInt,io,ic,j,readerror,p
  character :: tmpStr*(99),tmpStr999*(999),detname*(14),firstLine*(999),infile*(99),commandline*(999),parNameStr*(999)
  real :: injPost,injPrior
  real(double) :: tmpDat(maxMCMCpar),dtmpDat(maxMCMCpar), lon2ra,ra2lon
  
  exitcode = 0
  readerror = 0
  allocate(allDat(nchains,maxMCMCpar,maxIter))
  allocate(post(nchains,maxIter),prior(nchains,maxIter))
  
  t0 = 0.d0
  
  do ic = 1,nchains0
     infile = infiles(ic)
     open(unit=10,form='formatted',status='old',file=trim(infile),iostat=io) 
     if(io.ne.0) then
        write(stdErr,'(A)')'   Error:  File not found: '//trim(infile)//', aborting.'
        exitcode = 1
        return
     end if
     rewind(10)
     
     !if(prProgress.ge.2) write(stdOut,'(A,I3,A,I3,A20)',advance="no")'    File',ic,':  '//trim(infile)// &
     !     '    Using colour',colours(mod(ic-1,ncolours)+1),': '//colournames(colours(mod(ic-1,ncolours)+1))
     
     !Read the headers
     !Determine from the length of the first line whether this is output from before (~114 characters) 
     ! or after July 2009 (~29 characters), or LALInference (~180 characters):
     !
     !  before 7/09: first line is >80 characters long header (     Niter       Nburn    seed       <d|d>...
     !  after  7/09: first line is <80 characters long version number (  SPINspiral version:    1.00)
     ! LALInference format has long first line (  LALInference version:d4cd156ea4b0174d3fbd8b67ade5584981b34aed,2010-12-15 
     !05:58:56 +0000,cbc_bayesian_devel,Vivien Raymond <vivien.raymond@ligo.org>,CLEAN: All modifications committed)
     
     read(10,'(A999)',end=199,err=199) firstLine
     outputVersion = 1.0  ! SPINspiral output, after July 2009, keep 1<=oV<2
     if(len_trim(firstLine).gt.80 .and. len_trim(firstLine).lt.140)  outputVersion = 0.0  ! SPINspiral output, before July 2009
     if(len_trim(firstLine).ge.140)  outputVersion = 2.0  ! LALInference output (after December 2010), keep 2.0<=oV<3.0
     
     if(floor(outputVersion).eq.1) read(firstLine,'(A21,F8.2)') tmpStr,outputVersion
     
     if(outputVersion > 0.5) then
        ! Read command line between version number and first header. Read nothing if no command line:
        read(10,'(A999)',end=199,err=199) commandline
        if(commandline(1:14).eq.'  Command line') then
           outputVersion = 2.1  ! LALinference, without parameter IDs
           
           ! Determine waveformName from commandline, since IDs change rapidly
           waveformName = 'Unknown'
           i = index(trim(commandline),'--approximant', back=.false.)
           if(i.ne.0) read(commandline(i+13:),*,iostat=io) waveformName
        end if
        
        ! Read empty line between version number and first header:
        read(10,'(A)',end=199,err=199) tmpStr999
        if(index(trim(tmpStr999),'Tmax').eq.0) outputVersion = outputVersion + 0.1  ! No Tmax, LIM > Aug 2012 - v2.2
     end if
     
     
     if(outputVersion.lt.0.5) then
        read(10,*) &
             niter(ic),Nburn0(ic),seed(ic),DoverD,ndet(ic), nCorr(ic),nTemps(ic),Tmax(ic),Tchain(ic),networkSNR(ic)
     else
        
        if(outputVersion.gt.2.15) then  ! No Tmax
           read(10,*, iostat=io) niter(ic),Nburn0(ic),seed(ic),DoverD,ndet(ic), nCorr(ic), &
                nTemps(ic),Tchain(ic),networkSNR(ic),waveform, pnOrder,nMCMCpar
        else
           !read(10,'(I10,I12,I8,F22.10,I8,  2I9,I10,F12.1,F14.6,I11,F11.1,I10)') &
           read(10,*, iostat=io) niter(ic),Nburn0(ic),seed(ic),DoverD,ndet(ic), nCorr(ic), &
                nTemps(ic),Tmax(ic),Tchain(ic),networkSNR(ic),waveform, pnOrder,nMCMCpar
        end if
        
        if(io.ne.0) call quit_program_error('Error reading header line (nIter, Nburn, etc)',stdErr)
        
        if(outputVersion.lt.1.999) waveformName = waveforms(waveform)  ! SPINspiral > v0.5
     end if
     
     
     nMCMCpar0 = nMCMCpar  ! nMCMCpar may change when secondary parameters are computed
     
     read(10,*,end=199,err=199) tmpStr  ! Read (empty) line above detector info
     do i=1,ndet(ic)
        read(10,*) detnames(ic,i),snr(ic,i),flow(ic,i),fhigh(ic,i),t_before(ic,i),t_after(ic,i),FTstart(ic,i),deltaFT(ic,i), &
             samplerate(ic,i),samplesize(ic,i),FTsize(ic,i)
        detname = detnames(ic,i)
        j = len_trim(detname)
        if(trim(detname).eq.'Hanford')     detnr(ic,i) = 1
        if(trim(detname).eq.'LHO_4k')      detnr(ic,i) = 1
        if(trim(detname).eq.'Livingston')  detnr(ic,i) = 2
        if(trim(detname).eq.'LLO_4k')      detnr(ic,i) = 2
        if(trim(detname).eq.'Pisa')        detnr(ic,i) = 3
        if(trim(detname).eq.'VIRGO')       detnr(ic,i) = 3
        if(trim(detname).eq.'Aigo')        detnr(ic,i) = 4
     end do
     
     parID = 0
     revID = 0
     if(outputVersion > 0.5) then
        
        if(outputVersion.lt.2.1) then  ! Then integer parameter IDs are still used
           read(10,*,iostat=io) parID(1:nMCMCpar)  ! Read parameter IDs
           if(io.ne.0) then
              write(stdErr,'(//,A,//)')'  Error reading MCMC parameter IDs, aborting...'
              stop
           end if
        end if
        
     else  !Set the parameter IDs of an old MCMC output file manually:
        
        ! Assume 12-parameter Apostolatos:
        if(ic.eq.1) write(stdOut,'(A)')"  *** Note:  I'm using a hardcoded definition of the MCMC parameters,"// &
             " assuming that an Apostolatos, 1.5-pN, 12-parameter waveform was used..."
        nMCMCpar = 12  !Mc et tc ld a1 th ra de ph tJ pJ,ph_spin
        parID(1:12) = (/61,62,11,22,71,72,31,32,41,53,54,73/)
        waveform = 1
        pNorder = 1.5
        
        ! Assume 12-parameter Apostolatos: lq=log(q)
        !if(ic.eq.1) write(stdOut,'(A)')"  *** Note:  I'm using a hardcoded definition of the MCMC parameters,"// &
        !     " assuming that an Apostolatos, 1.5-pN, 12-parameter waveform was used..."
        !nMCMCpar = 12  !Mc lq tc ld a1 th ra de ph tJ pJ,ph_spin
        !parID(1:12) = (/61,68,11,22,71,72,31,32,41,53,54,73/)
        !waveform = 1
        !pNorder = 1.5
        
        ! Assume 12-parameter SpinTaylor:
        !if(ic.eq.1) write(stdOut,'(A)')"  *** Note:  I'm using a hardcoded definition of the MCMC parameters,"// &
        !" assuming that a SpinTaylor, 3.5-pN, 12-parameter waveform was used..."
        !nMCMCpar = 12  !
        !parID(1:12) = (//)  !Needs to be filled in before use
        !waveform = 2
        !pNorder = 3.5
        
        ! Assume 15-parameter SpinTaylor:
        !if(ic.eq.1) write(stdOut,'(A)')"  *** Note:  I'm using a hardcoded definition of the MCMC parameters,"// &
        !" assuming that a SpinTaylor, 3.5-pN, 15-parameter waveform was used..."
        !nMCMCpar = 15  !
        !parID(1:15) = (//)  !Needs to be filled in before use
        !waveform = 3
        !pNorder = 3.5
        
        nMCMCpar0 = nMCMCpar  ! nMCMCpar may change when secondary parameters are computed
     end if
     
     read(10,'(A)',end=199,err=199) parNameStr  ! Read line with column headers (Cycle, log Post., Prior, etc)
     if(outputVersion.ge.2.1) then
        read(10,'(A)',end=199,err=199) parNameStr  ! Read empty line
        read(10,'(A)',end=199,err=199) parNameStr  ! Read empty line
        read(10,'(A)',end=199,err=199) parNameStr  ! Read empty line
        call parNames2IDs(trim(parNameStr),nMCMCpar, parID)  ! Convert parameter names to integer IDs
     end if
     
     if(prProgress.ge.3) write(stdOut,'(A,I3,A)', advance='no') '  Chain',ic, ':  parameters identified:  '
     do i=1,nMCMCpar
        revID(parID(i)) = i  ! Reverse ID
        if(prProgress.ge.3) write(stdOut,'(A,1x)', advance='no') trim(parNames(parID(i)))
     end do
     if(prProgress.ge.3) write(stdOut,*)
     
     
     i = 1
     if(outputVersion.gt.1.999) then  ! LALInferenceMCMC, no injection value in file - get it from injection file:
        if(outputVersion.gt.1.999) call get_LIM_injection_values(ic,nMCMCpar, startval, injPost,injPrior)
        is(ic,i) = -1
        post(ic,i)  = injPost
        prior(ic,i) = injPrior
        allDat(ic,1:nMCMCpar,i) = startval(ic,1:nMCMCpar,1)
        
        i = 2  ! Continue with reading the "second" line of the Markov chains
     end if
     
     
     do while(i.le.maxIter)
        if(outputVersion < 0.5) then
           read(10,*,iostat=io) tmpInt,post(ic,i),tmpDat(1:nMCMCpar)
        else
           read(10,*,iostat=io) tmpInt,post(ic,i),prior(ic,i),tmpDat(1:nMCMCpar)
        end if
        
        if(io.lt.0) exit  ! EOF
        if(io.gt.0) then  ! Read error
           if(readerror.eq.1) then  ! Read error in previous line as well
              if(i.lt.25) then
                 write(stdErr,'(A,I7)',advance="no")'  Read error in file '//trim(infile)//', line',i
                 write(stdErr,'(A,/)')'  Aborting program...'
                 stop
              else
                 write(stdOut,'(A,I7)',advance="no")'  Read error in file '//trim(infile)//', line',i
                 i = i-1
                 write(stdOut,'(2(A,I8))',advance='no')'  iteration:',tmpInt,', last iteration read successfully:',tmpInt
                 write(stdOut,*)
                 exit
              end if
           end if
           readerror = 1
           i = i-1
           cycle
        end if
        readerror = 0
        
        is(ic,i) = real(tmpInt)
        !if(i.lt.10) print*,tmpInt
        
        ! GPS time doesn't fit in single-precision variable
        if(ic.eq.1 .and. i.le.2 .and. t0.lt.1.d0) then
           dtmpDat = 0.d0
           do p=1,nMCMCpar
              if(parID(p).ge.11.and.parID(p).le.19) then
                 dtmpDat(p) = dble(nint(tmpDat(p)/10.d0)*10)   !'GPS base time', rounded off to the nearest 10s,
                 t0 = dtmpDat(p)                               !  to allow x.xx labels on the plots for GPS-t0
                 GPStime = floor(tmpDat(p)+0.05)               !Always floor, unless >x.95s. Use e.g. in file names
              end if
           end do
        end if
        allDat(ic,1:nMCMCpar,i) = real(tmpDat(1:nMCMCpar) - dtmpDat(1:nMCMCpar))
        
        ! 'Thin' the output by reading every thin-th line, after you've read the injection and starting values:
        if(thin.gt.1.and.i.gt.1) then
           do j=1,thin-1
              read(10,*,iostat=io) tmpStr
              if(io.lt.0) exit  ! EOF
           end do
        end if
        
        tmpStr = tmpStr  ! Remove 'set but never used' warning
        if(1.eq.2) then
           ! In case you ran with lon rather than RA:
           allDat(ic,revID(31),i) = real(lon2ra(dble(allDat(ic,revID(31),i)),t0))
           !allDat(ic,revID(31),i) = real(ra2lon(dble(allDat(ic,revID(31),i)),t0))  
           
           ! In case only the injection value is lon rather than RA:
           if(i.eq.1) allDat(ic,revID(31),i) = real(lon2ra(dble(allDat(ic,revID(31),i)),t0))
           
           ! In case only the injection value is lon rather than RA:
           if(i.ne.1) allDat(ic,revID(31),i) = real(ra2lon(dble(allDat(ic,revID(31),i)),t0))
           
           ! In case all but the injection value is lon rather than RA:
           if(i.ne.1) allDat(ic,revID(31),i) = real(lon2ra(dble(allDat(ic,revID(31),i)),t0))
        end if
        
        
        i = i+1
        if(tmpInt.ge.maxChLen) exit
     end do !i
     
     if(i.ge.maxIter-2) write(stdErr,'(A)',advance="no")'   *** WARNING ***   Not all lines in this file were read    '
     
199  continue
     close(10)
     ntot(ic) = i-1
     n(ic) = ntot(ic)  ! n can be changed in rearranging chains, ntot wont be changed
     !if(prProgress.ge.2.and.update.ne.1) write(stdOut,'(1x,3(A,I9),A1)')' Lines:',ntot(ic),', iterations:', &
     !nint(is(ic,ntot(ic))),', burn-in:',Nburn(ic),'.'
  end do !do ic = 1,nchains0
  
  if(sum(ntot).lt.10) &
       call quit_program_error('read_mcmcfiles(): Fewer than 10 data points were read - I must have read the input file(s)'// &
       ' incorrectly somehow', stdErr)
  
end subroutine read_mcmcfiles
!***********************************************************************************************************************************



!***********************************************************************************************************************************
!> \brief  Convert parameter names to integer IDs
!!
!! \param  parNameStr  String containing parameter names
!! \param  nMCMCpar    Number of (primary and secondary) MCMC parameters
!!
!! \retval parID       Array of parameter IDs

subroutine parNames2IDs(parNameStr,nMCMCpar, parID)
  use SUFR_constants, only: stdErr
  use SUFR_system, only: quit_program_error
  use analysemcmc_settings, only: maxMCMCpar
  
  implicit none
  character, intent(in) :: parNameStr*(*)
  integer, intent(in) :: nMCMCpar
  integer, intent(out) :: parID(maxMCMCpar)
  
  integer, parameter :: npIDs = 17
  integer :: pr1,pr2,pIDs(npIDs)
  character :: pnames(npIDs)*(19),pars(nMCMCpar+3)*(19)
  
  parID = 0
  pnames(1:npIDs) = [character(len=19) :: 'iota','psi','dec','ra','dist','phi_orb','time','q','mc',  &
       'a1','theta1','phi1','a2','theta2','phi2','eta','logq']  ! CHECK: time = t40? tc?
  pIDs(1:npIDs) = (/                       51,    52,   32,   31,  22,    41,       11,    67, 61,  &
       71,   72,      73,    81,  82,      83,    62,   68/)
  
  read(parNameStr,*) (pars(pr1),pr1=1,nMCMCpar+3)
  do pr1=1,nMCMCpar+3
     do pr2=1,npIDs
        if(trim(pars(pr1)).eq.trim(pnames(pr2))) parID(pr1-3) = pIDs(pr2)
     end do
     if(pr1.gt.3) then
        if(parID(pr1-3).eq.0) call quit_program_error('parNames2IDs: Parameter name not recognised: '//trim(pars(pr1)),StdErr)
     end if
  end do
  
end subroutine parNames2IDs
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief Extract info from the chains and print some of it to screen
!!
!! \retval exitcode  Exit status code (0=ok)
!! 
!!  \par
!!  Print MCMC run info:
!!  - determine Lmax, burn-in,  
!!  - print chain info,  
!!  - determine thinning for chains plots,  
!!  - change/add parameters, 
!!  - determine injection, start, Lmax values of chains,  
!!  - compute jumps,  
!!  - construct output file name,  
!!  - store data in selDat (from allDat)

subroutine mcmcruninfo(exitcode)  
  use SUFR_kinds, only: double
  use SUFR_constants, only: stdOut,stdErr, rpi
  use SUFR_statistics, only: compute_median_sp
  use SUFR_system, only: swapreal, quit_program_error
  use aM_constants, only: detabbrs
  
  use analysemcmc_settings, only: Nburn,update,prRunInfo,NburnFrac,thin,autoBurnin,prChainInfo,chainPlI,changeVar,prProgress
  use analysemcmc_settings, only: prInitial,mergeChains,maxMCMCpar,plInject,plStart,htmlOutput
  
  use general_data, only: allDat,post,ntot,n,nchains,nchains0,infiles,contrChain,startval,fixedpar,selDat,iloglmax,icloglmax
  use general_data, only: contrChains,parNames,nfixedpar,outputname
  
  use mcmcrun_data, only: niter,Nburn0,detnames,detnr,parID,seed,snr,revID,ndet,flow,fhigh,t_before,nCorr,nTemps,Tmax,Tchain
  use mcmcrun_data, only: networkSNR,waveformName,pnOrder,nMCMCpar,t_after,FTstart,deltaFT,samplerate,samplesize,FTsize
  use mcmcrun_data, only: nMCMCpar0,t0,GPStime,totthin,loglmaxs,totiter,loglmax,totpts,totlines,offsetrun,spinningRun, avgTotThin
  use mcmcrun_data, only: outputVersion
  use chain_data, only: is,isburn,DoverD,jumps
  use plot_data, only: ncolours,colours,colournames,maxdots
  
  implicit none
  integer, intent(out) :: exitcode
  integer :: i,ic,j,p,maxLine
  real(double) :: lon2ra,gmst
  character :: infile*(99)
  
  exitcode = 0
  totiter = 0
  do ic = 1,nchains0
     totiter = totiter + nint(is(ic,ntot(ic)))
  end do
  !if(prProgress.ge.2.and.update.ne.1) write(stdOut,'(A10,65x,2(A,I9),A1)')'Total:',' Lines:',sum(ntot(1:nchains0)),',  &
  !iterations:',totiter
  
  
  
  ! Print run info (detectors, SNR, amount of data, FFT, etc):
  if(prRunInfo.gt.0.and.update.eq.0) then
     write(stdOut,*)
     if(htmlOutput.ge.1) then
        write(stdOut,'(A)')'<br><hr><a name="runinfo"></a><font size="1"><a href="#top" title="Go to the top of the page">'// &
             'top</a></font><h2>Run info</h2>'
     else
        if(prRunInfo.eq.1) write(stdOut,'(/,A)')'  Run information for chain 1:'
        if(prRunInfo.eq.2) write(stdOut,'(/,A)')'  Run information:'
     end if
     
     do ic = 1,nchains0
        if((prRunInfo.eq.1.and.ic.eq.1) .or. prRunInfo.eq.2) then
           infile = infiles(ic)
           if(outputVersion.lt.0.5) then
              if(prRunInfo.le.2.and.ic.eq.1 .or. prRunInfo.eq.3) &
                   write(stdOut,'(4x,A7,A12,A13,A10,A12,A8,A8)')'Chain','file name','colour','Niter','Nburn','seed','Ndet'
              write(stdOut,'(4x,I7,A12,A13,I10,I12,I8,I8)')ic,trim(infile(13:99)), &
                   trim(colournames(colours(mod(ic-1,ncolours)+1))),niter(ic),Nburn0(ic),seed(ic),ndet(ic)
           else
              if(prRunInfo.le.2.and.ic.eq.1 .or. prRunInfo.eq.3) then
                 if(htmlOutput.ge.1) then
                    write(stdOut,'(4x,A3,A7,A12,A13,A10,A12,A11,A8,  2A8,2A8, A8, A16, A8, A8, 3x,A8, A4)') '<b>','Chain', &
                         'file name','colour','Niter','Nburn','seed','Ndet',  'Ncorr','Ntemp','Tmax','Tchain','NetSNR', &
                         '&lt;d|d&gt;','pN','Npar','WaveForm','</b>'
                 else
                    write(stdOut,'(4x,A7,A12,A13,A10,A12,A11,A8,  2A8,2A8, A8, A10, A8, A8, 3x,A8)') 'Chain','file name','colour', &
                         'Niter','Nburn','seed','Ndet',  'Ncorr','Ntemp','Tmax','Tchain','NetSNR','<d|d>','pN','Npar','WaveForm'
                 end if
              end if
              write(stdOut,'(4x,I7,A12,  A13,I10,I12,I11,I8,  2I8,2F8.1,F8.3,F10.2,F8.1,I8, 3x,A)') ic,trim(infile(19:99)), &
                   trim(colournames(colours(mod(ic-1,ncolours)+1))),niter(ic),Nburn0(ic),seed(ic),ndet(ic), &
                   nCorr(ic),nTemps(ic),real(Tmax(ic)),Tchain(ic),networkSNR(ic),DoverD,pnOrder,nMCMCpar,trim(waveformName)
           end if
        end if
        
        if((prRunInfo.le.2.and.ic.eq.nchains0) .or. prRunInfo.eq.3) then
           write(stdOut,*)
           if(htmlOutput.ge.1) then
              write(stdOut,'(A3,A14,A3,A18,4A12,A22,A17,3A14,A4)') '<b>','Detector','Nr','SNR','f_low','f_high','before tc', &
                   'after tc','Sample start (GPS)','Sample length','Sample rate','Sample size','FT size','</b>'
           else
              write(stdOut,'(A14,A3,A18,4A12,A22,A17,3A14)') 'Detector','Nr','SNR','f_low','f_high','before tc','after tc', &
                   'Sample start (GPS)','Sample length','Sample rate','Sample size','FT size'
           end if
           
           do i=1,ndet(ic)
              write(stdOut,'(A14,I3,F18.8,4F12.2,F22.8,F17.7,3I14)')trim(detnames(ic,i)),detnr(ic,i),snr(ic,i),flow(ic,i), &
                   fhigh(ic,i),t_before(ic,i),t_after(ic,i),FTstart(ic,i),deltaFT(ic,i),nint(samplerate(ic,i)),samplesize(ic,i), &
                   FTsize(ic,i)
           end do
           write(stdOut,*)
        end if
     end do !do ic = 1,nchains0
  end if  !prRunInfo.gt.0
  
  maxLine = maxval(n(1:nchains0))
  
  
  
  !*** Until now, Nburn is the iteration number.
  !*** From here on, Nburn is the line number, while isburn is the iteration number
  do ic=1,nchains0
     if(Nburn(ic).le.0) Nburn(ic) = Nburn0(ic)
     if(abs(NburnFrac).gt.1.e-4.and.abs(NburnFrac).lt.1.) then
        Nburn(ic) = nint(is(ic,n(ic)) * abs(NburnFrac))
     else
        if(Nburn(ic).ge.nint(is(ic,n(ic)))) then
           !print*,Nburn(ic),nint(is(ic,n(ic)))
           if(Nburn0(ic).ge.nint(is(ic,n(ic)))) then
              write(stdErr,'(A,I3)')'   *** WARNING ***  Nburn larger than Nchain, setting Nburn to 10% for chain',ic
              Nburn(ic) = nint(is(ic,n(ic))*0.1)
           else
              Nburn(ic) = Nburn0(ic)
           end if
        end if
     end if
  end do
  
  do ic=1,nchains0
     isburn(ic) = real(Nburn(ic))
     do i=1,ntot(ic)
        if(is(ic,i).le.isburn(ic)) Nburn(ic) = i   !isburn is the true iteration number at which the burn-in ends
        totthin(ic) = nint(isburn(ic)/real(Nburn(ic)))
     end do
  end do
  avgTotThin = sum(isburn(1:nchains0))/real(sum(Nburn(1:nchains0)))  ! Total thinning, averaged over all chains
  
  
  
  
  ! Get point with absolute maximum likelihood over all chains:
  loglmax = -1.d99
  loglmaxs = -1.d99
  iloglmax = 0
  icloglmax = 0
  do ic=1,nchains0
     do i=3,ntot(ic)  !3: exclude injection and starting values
        if(post(ic,i).gt.loglmaxs(ic)) then
           loglmaxs(ic) = post(ic,i)
           if(post(ic,i).gt.loglmax) then
              loglmax = post(ic,i)
              iloglmax = i
              icloglmax = ic
           end if
        end if
     end do
  end do
  
  if(ilogLmax*iclogLmax.eq.0) call quit_program_error('No max logL found - perhaps I read the input files incorrectly?', &
       stdErr)
  
  if(prRunInfo.ge.1) then
     ic = icloglmax
     i = iloglmax
     infile = infiles(ic)
     write(stdOut,'(A,I4,2(A,I9))', advance='no')'    Maximum likelihood point:   chain:',ic,' ('//trim(infile(19:99))// &
          '),   line:', i*thin+7+ndet(ic), ',   iteration:',nint(is(ic,i))
     if(htmlOutput.ge.1) then
        write(stdOut,'(A,F10.3,A,F7.2,A)') ',   <b>max log(L):',loglmax,'  -> SNR:',sqrt(2*loglmax),'</b>.'
     else
        write(stdOut,'(A,F10.3,A,F7.2,A)') ',   max log(L):',loglmax,'  -> SNR:',sqrt(2*loglmax),'.'
     end if
     
     !Test: get parameter values for L=Lmax
     !if(prProgress.ge.3) then
     !   write(stdOut,'(I10,F14.6,1x,2F12.7,F20.8,9F12.7)')nint(is(ic,i)),loglmax,allDat(ic,2:3,i),allDat(ic,4,i)+t0, &
     !     allDat(ic,5:13,i)
     !   call mc_eta_2_m1_m2r(allDat(ic,2,i),allDat(ic,3,i), m1,m2)
     !   write(stdOut,'(F9.4,F10.4,F21.6,F11.4,F10.4,F12.4,2F11.4,F12.4,3F11.4)')m1,m2,allDat(ic,4,i)+t0,exp(allDat(ic,5,i)), &
     !     allDat(ic,6,i),acos(allDat(ic,7,i))*r2d,allDat(ic,8,i)*r2h,asin(allDat(ic,9,i))*r2d,allDat(ic,10,i)*r2d, &
     !     asin(allDat(ic,11,i))*r2d,allDat(ic,12,i)*r2d,allDat(ic,13,i)*r2d
     !end if
     if(prRunInfo.ge.1) write(stdOut,*)
  end if
  
  
  !*** AutoBurnin: for each chain, get the first point where log(L) > log(L_max)-autoBurnin:
  if(abs(autoBurnin).gt.1.e-10) then
     if(autoBurnin.lt.-1.e-10) autoBurnin = real(nMCMCpar0)/2.  ! The default value for autoBurnin = Npar/2
     loop1: do ic=1,nchains0
        isburn(ic) = is(ic,ntot(ic))  ! Set burn-in to last iteration, so that chain is completely excluded if condition never met
        Nburn(ic) = ntot(ic)
        do i=2,ntot(ic)  ! i=1 is injection value?
           if(post(ic,i).gt.real(loglmax)-autoBurnin) then
              isburn(ic) = is(ic,i)
              Nburn(ic) = i
              cycle loop1
           end if
        end do
     end do loop1
  end if
  
  
  ! Determine the total number of iterations and lines in the input and data points for statistics;
  ! determine how many and which chains contribute
  totiter = 0
  totpts  = 0
  contrChains = 0
  contrChain = 0
  totlines = sum(ntot(1:nchains0))  !Total number of lines in the input files (including burn-in)
  do ic=1,nchains0
     totiter = totiter + nint(is(ic,ntot(ic)))  !Total number of iterations (including the burn-in)
     totpts = totpts + n(ic)-Nburn(ic)          !Total number of data points available for statistics, after removing the burn-in
     if(n(ic).gt.Nburn(ic)) then
        contrChains = contrChains + 1
        contrChain(ic) = 1
     end if
  end do
  
  
  !*** Print chain info to screen:
  ! Print info on number of iterations, burn-in, thinning, etc.:
  if(prChainInfo.ge.2.and.update.ne.1.and.htmlOutput.ge.1) write(stdOut,'(A)')'<br><hr><a name="chaininfo"></a>'// &
       '<font size="1"><a href="#top" title="Go to the top of the page">top</a></font><h2>Chain info</h2>'
  do ic=1,nchains0
     infile = infiles(ic)
     if(prChainInfo.ge.2.and.update.ne.1) then
        if(contrChain(ic).eq.1) then
           write(stdOut,'(A6)',advance="no")'    * '  !Flag contributing chains
        else
           write(stdOut,'(A6)',advance="no")'      '
        end if
        write(stdOut,'(A2,I2,A1,A10,A12)',advance="no") 'Ch',ic,':',trim(infile(19:99)), &
             ', '//colournames(colours(mod(ic-1,ncolours)+1))//'.'
        write(stdOut,'(A,ES7.1,A,ES7.1,A1)',advance="no") '  Lines/iter: ',real(n(ic)),'/',is(ic,n(ic)),'.'
        write(stdOut,'(A,ES7.1,A,ES7.1,A1)',advance="no") '  Burn-in: ',real(Nburn(ic)),'/',isburn(ic),'.'
        write(stdOut,'(A,F8.2,A,F6.1,A,F4.1,A1)',advance="no") '  Lmx:',loglmaxs(ic),', dLmx:',abs(loglmax-loglmaxs(ic)), &
             '/',autoBurnin,'.'
        write(stdOut,'(A,I3,A,I4,A1)',advance="no") ' Thin: file:',nint(is(ic,n(ic))/real(n(ic)*max(thin,1))), &
             ', tot:',totthin(ic),'.'
        write(stdOut,'(A,ES7.1,A1)') '  Data pts: ',abs(real(n(ic)-Nburn(ic))),'.'
     end if
  end do
  if(prChainInfo.ge.1.and.update.ne.1) then
     write(stdOut,'(4x,A, A,ES10.3, A,ES10.3, A,I4, A,ES9.2)', advance='no') 'All chains:','  # lines:',real(totlines), &
          ',  # iterations:',real(totiter), ',  thinning:',nint(avgTotThin), 'x,  med.burnin:', &
          compute_median_sp(real(isburn(1:nChains0)))
     
     if(htmlOutput.ge.1) write(stdOut,'(A3)', advance='no') '<b>'
     write(stdOut,'(A,ES10.3,  A2,F5.1, A,I3,A1,I2,A1)', advance='no') ',  # dat.pts after burnin:',real(totpts), &
          ' (',real(totpts)/real(totlines)*100,'%), contrib.chains:',contrChains,'/',nchains0,'.'
     if(htmlOutput.ge.1) then
        write(stdOut,'(A)') '</b>'
     else
        write(stdOut,*) ''
     end if
  end if
  
  
  
  !*** Determine extra thinning for logL, chain, jump plots
  if(chainPlI.le.0) then
     !if(sum(ntot(1:nchains0)).gt.maxdots) then  !Change the number of points plotted in chains,logL, etc. (For all output formats)
     !Use ntot and nchains0, since n is low if many points are in the burn-in:
     chainPlI = max(1,nint(real(sum(ntot(1:nchains0)))/real(maxdots)))  
     if(prChainInfo.ge.1.and.update.eq.0) then
        if(chainPlI.gt.1) then  ! Change the number of points plotted in chains,logL, etc. (For all output formats)
           if(chainPlI.eq.1) then
              write(stdOut,'(A)', advance='no')'    Plotting every'
           else if(chainPlI.eq.2) then
              write(stdOut,'(A,I2,A)', advance='no')'    Plotting every',chainPlI,'-nd'
           else if(chainPlI.eq.3) then
              write(stdOut,'(A,I2,A)', advance='no')'    Plotting every',chainPlI,'-rd'
           else
              write(stdOut,'(A,I4,A)', advance='no')'    Plotting every',chainPlI,'-th'
           end if
           write(stdOut,'(A,I5,A,I5,A)')' state in likelihood, chains, jumps, etc. plots.  Average total thinning is', &
                nint(avgTotThin),'x, for these plots it is',nint(avgTotThin*real(chainPlI)),'x.'
        else
           write(stdOut,'(A,I5,A)')'    Plotting *every* state in likelihood, chains, jumps, etc. plots.'// &
                '  Average total thinning remains',nint(avgTotThin),'x for these plots.'
        end if
     end if
     write(stdOut,*)
  end if
  !if(prRunInfo.gt.0) write(stdOut,*)
  
  
  
  
  ! *** Change some MCMC parameters:
  
  if(changeVar.ge.1) then
     if(htmlOutput.eq.0.and.prProgress.ge.2.and.update.eq.0) write(stdOut,'(A)',advance="no")'  Changing some parameters...   '
     
     
     ! Calculate Mc from Mc_16 (Mc^(1/6)):
     if(revID(61).eq.0 .and. revID(65).ne.0) then
        if(htmlOutput.eq.0.and.prProgress.ge.2.and.update.eq.0) write(stdOut,'(A)')'  Computing Mc from Mc^(1/6)'
        parID(nMCMCpar+1) = 61    ! Mc
        revID(61) = nMCMCpar + 1  ! Mc
        nMCMCpar = nMCMCpar + 1
        
        if(nMCMCpar.gt.maxMCMCpar) call error_maxMCMCpar_quit(StdErr, maxMCMCpar, nMCMCpar)  ! Stop program with error
        
        do ic=1,nchains0
           allDat(ic,revID(61),1:ntot(ic)) = allDat(ic,revID(65),1:ntot(ic))**6   ! Mc = [Mc_16]^6
        end do
     end if
     
     
     ! Calculate q from log(q):
     if(revID(67).eq.0 .and. revID(68).ne.0) then
        if(htmlOutput.eq.0.and.prProgress.ge.2.and.update.eq.0) write(stdOut,'(A)')'  Computing q from log(q)'
        parID(nMCMCpar+1) = 67    ! q
        revID(67) = nMCMCpar + 1  ! q
        nMCMCpar = nMCMCpar + 1
        
        if(nMCMCpar.gt.maxMCMCpar) call error_maxMCMCpar_quit(StdErr, maxMCMCpar, nMCMCpar)  ! Stop program with error
        
        do ic=1,nchains0
           
           if(changeVar.eq.4) then ! Folding log(q) for comparison
              do j=1,ntot(ic)
                 if(allDat(ic,revID(68),j).gt.0.0) then
                 allDat(ic,revID(68),j) = -allDat(ic,revID(68),j)
                 end if
              end do
           end if
           if(changeVar.eq.3) then  ! for phi > pi -> logq = -logq & phi = phi - pi
              do j=1,ntot(ic)
                 if(allDat(ic,revID(41),j).gt.rpi) then
                    allDat(ic,revID(41),j) =  allDat(ic,revID(41),j) - rpi                                     ! phi = phi - pi
                    allDat(ic,revID(68),j) = -allDat(ic,revID(68),j)                                           ! log_q = -log_q
                    if(revID(63)*revID(64).ne.0) call swapreal(allDat(ic,revID(63),j),allDat(ic,revID(64),j))  ! swap m1 <-> m2
                 end if
              end do
           end if
           
           allDat(ic,revID(67),1:ntot(ic)) = 10.0 ** (allDat(ic,revID(68),1:ntot(ic)))                         ! q = 10**log_q
        end do
     end if
     
     
     ! Calculate eta and log(q) from q:
     if(revID(62).eq.0 .and. revID(67).ne.0) then
        if(htmlOutput.eq.0.and.prProgress.ge.2.and.update.eq.0) write(stdOut,'(A)')'  Computing eta and log(q) from q'
        parID(nMCMCpar+1) = 62    ! Eta
        parID(nMCMCpar+2) = 63    ! M1
        parID(nMCMCpar+3) = 64    ! M2
        parID(nMCMCpar+4) = 66    ! Mtot
        parID(nMCMCpar+5) = 68    ! logq 
        revID(62) = nMCMCpar + 1  ! Eta
        revID(63) = nMCMCpar + 2  ! M1
        revID(64) = nMCMCpar + 3  ! M2
        revID(66) = nMCMCpar + 4  ! Mtot
        revID(68) = nMCMCpar + 5  ! logq
        nMCMCpar = nMCMCpar + 5
        
        if(nMCMCpar.gt.maxMCMCpar) call error_maxMCMCpar_quit(StdErr, maxMCMCpar, nMCMCpar)  ! Stop program with error
        
        do ic=1,nchains0
           if(changeVar.eq.4) then ! Folding log(q) for comparison
              if(htmlOutput.eq.0.and.prProgress.ge.2.and.update.eq.0.and.ic.eq.1) write(stdOut,'(A)') ' Folding log(q)'
              do j=1,ntot(ic)
                 if(allDat(ic,revID(67),j).gt.1.0) then
                    allDat(ic,revID(67),j) = 1.0 / allDat(ic,revID(67),j)
                 end if
              end do
           end if
           if(changeVar.eq.3) then  ! for phi > pi -> q = 1/q & phi = phi - pi
              if(htmlOutput.eq.0.and.prProgress.ge.2.and.update.eq.0.and.ic.eq.1) write(stdOut,'(A)') ' Swapping q, phi'
              do j=1,ntot(ic)
                 if(allDat(ic,revID(41),j).gt.rpi) then
                    allDat(ic,revID(41),j) = allDat(ic,revID(41),j) - rpi                                      ! phi = phi - pi
                    allDat(ic,revID(67),j) = 1.0 / allDat(ic,revID(67),j)                                      ! q = 1/q
                    if(revID(63)*revID(64).ne.0) call swapreal(allDat(ic,revID(63),j),allDat(ic,revID(64),j))  ! swap m1 <-> m2
                 end if
              end do
           end if
           
           allDat(ic,revID(62),1:ntot(ic)) =  &
                allDat(ic,revID(67),1:ntot(ic)) / (allDat(ic,revID(67),1:ntot(ic)) + 1.0)**2                ! eta = q/(1+q)^2
           allDat(ic,revID(68),1:ntot(ic)) = log10(allDat(ic,revID(67),1:ntot(ic)))                         ! log_q = log(q)
           
           
           if(htmlOutput.eq.0.and.prProgress.ge.2.and.update.eq.0.and.ic.eq.1) write(stdOut,'(A)') &
                '  Computing M1, M2 and Mtot from Mc and q'
           do i=1,ntot(ic)
              call mc_q_2_m1_m2r(allDat(ic,revID(61),i),allDat(ic,revID(67),i), allDat(ic,revID(63),i),allDat(ic,revID(64),i))
           end do
           allDat(ic,revID(66),1:ntot(ic)) = allDat(ic,revID(63),1:ntot(ic)) + allDat(ic,revID(64),1:ntot(ic))     ! Mtot = m1 + m2
        end do
     end if
     
     
     ! Calculate the individual masses from Mch and eta:
     if(revID(61)*revID(62).ne.0 .and. revID(63)+revID(64).eq.0) then
        if(htmlOutput.eq.0.and.prProgress.ge.2.and.update.eq.0) write(stdOut,'(A)')'  Computing M1, M2, Mtot from Mc, eta'
        parID(nMCMCpar+1) = 63    ! M1
        parID(nMCMCpar+2) = 64    ! M2
        parID(nMCMCpar+3) = 66    ! Mtot
        revID(63) = nMCMCpar + 1  ! M1
        revID(64) = nMCMCpar + 2  ! M2
        revID(66) = nMCMCpar + 3  ! Mtot
        nMCMCpar = nMCMCpar + 3
        
        if(nMCMCpar.gt.maxMCMCpar) call error_maxMCMCpar_quit(StdErr, maxMCMCpar, nMCMCpar)  ! Stop program with error
        
        do ic=1,nchains0
           do i=1,ntot(ic)
              call mc_eta_2_m1_m2r(allDat(ic,revID(61),i),allDat(ic,revID(62),i), allDat(ic,revID(63),i),allDat(ic,revID(64),i))
           end do
           allDat(ic,revID(66),1:ntot(ic)) = allDat(ic,revID(63),1:ntot(ic)) + allDat(ic,revID(64),1:ntot(ic))     ! Mtot = m1 + m2
        end do
        
        
        ! Calculate Mc, eta from the individual masses:
     else if(revID(61)+revID(62).eq.0 .and. revID(63)*revID(64).ne.0) then
        
        if(htmlOutput.eq.0.and.prProgress.ge.2.and.update.eq.0) write(stdOut,'(A)')'  Computing Mc, eta from M1, M2'
        parID(nMCMCpar+1) = 61    ! Mc
        parID(nMCMCpar+2) = 62    ! eta
        revID(61) = nMCMCpar + 1  ! Mc
        revID(62) = nMCMCpar + 2  ! eta
        nMCMCpar = nMCMCpar + 2
        
        if(nMCMCpar.gt.maxMCMCpar) call error_maxMCMCpar_quit(StdErr, maxMCMCpar, nMCMCpar)  ! Stop program with error
        
        do ic=1,nchains0
           do i=1,ntot(ic)
              call m1_m2_2_mc_etar(allDat(ic,revID(63),i),allDat(ic,revID(64),i), allDat(ic,revID(61),i),allDat(ic,revID(62),i))
           end do
        end do
        
        ! Compute total mass (var 66) and mass ratio (var 67) (q=M2/M1, not \eta) from the individual masses:
        ! (var 65 is reserved for Mc^(1/6)) & convert q -> 1/q, logq -> -logq and phi -> phi -pi for phi > pi
        if(htmlOutput.eq.0.and.prProgress.ge.2.and.update.eq.0) write(stdOut,'(A)')'  Computing Mtot, q, log(q) from masses'
        parID(nMCMCpar+1) = 66    ! Mtot
        parID(nMCMCpar+2) = 67    ! q
        parID(nMCMCpar+3) = 68    ! log(q)
        revID(66) = nMCMCpar + 1  ! Mtot
        revID(67) = nMCMCpar + 2  ! q
        revID(68) = nMCMCpar + 3  ! log(q)
        nMCMCpar = nMCMCpar + 3
        
        if(nMCMCpar.gt.maxMCMCpar) call error_maxMCMCpar_quit(StdErr, maxMCMCpar, nMCMCpar)  ! Stop program with error
        
        do ic=1,nchains0
           allDat(ic,revID(67),1:ntot(ic)) = allDat(ic,revID(64),1:ntot(ic)) / allDat(ic,revID(63),1:ntot(ic))     ! q = m2 / m1
           
           if(changeVar.eq.4) then ! Folding log(q) for comparison
              do j=1,ntot(ic)
                 if(allDat(ic,revID(67),j).gt.1.0) then
                    allDat(ic,revID(67),j) = 1.0 / allDat(ic,revID(67),j)
                 end if
              end do
           end if
           if(changeVar.eq.3) then  ! m2/m1 for q<1, & phi<pi and m1/m2 for q>1 & phi >pi
              do j=1,ntot(ic)
                 if(allDat(ic,revID(41),j).gt.rpi) then
                    allDat(ic,revID(41),j) = allDat(ic,revID(41),j) - rpi                                          ! phi = phi - pi
                    allDat(ic,revID(67),j) = 1.0 / allDat(ic,revID(67),j)                                          ! q = 1/q = m1/m2
                    if(revID(63)*revID(64).ne.0) call swapreal(allDat(ic,revID(63),j),allDat(ic,revID(64),j))      ! swap m1 <-> m2
                 end if
              end do
           end if
           
           allDat(ic,revID(66),1:ntot(ic)) = allDat(ic,revID(63),1:ntot(ic)) + allDat(ic,revID(64),1:ntot(ic))     ! Mtot = m1 + m2
           allDat(ic,revID(68),1:ntot(ic)) = log10(allDat(ic,revID(67),1:ntot(ic)))                                ! log_q = log(q)
        end do
     end if !if(revID(61)+revID(62).eq.0 .and. revID(63)*revID(64).ne.0) 
     
     
     ! Compute inclination and polarisation angle from RA, Dec, theta_J0, phi_J0:
     if(revID(11)*revID(31)*revID(32)*revID(53)*revID(54).ne.0) then  ! Then all of these parameters are defined
        do ic=1,nchains0
           do i=1,ntot(ic)
              ! Input: RA, Dec, phi_Jo (hh->RA), theta_Jo (in rad), output: inclination, polarisation angle (rad):
              call compute_incli_polangr(allDat(ic,revID(31),i), asin(allDat(ic,revID(32),i)), &
                   real(lon2ra(dble(allDat(ic,revID(54),i)), dble(allDat(ic,revID(11),i)) + t0)), asin(allDat(ic,revID(53),i)), &
                   allDat(ic,revID(53),i),allDat(ic,revID(54),i))  
              allDat(ic,revID(53),i) = cos(allDat(ic,revID(53),i))    !i -> cos(i)
           end do
        end do !ic
        parID(revID(53)) = 51  ! Was sin(thJ0), now cos(i)
        parID(revID(54)) = 52  ! Was phi_J0, now psi
        revID(51) = revID(53)  ! Now cos(i)
        revID(52) = revID(54)  ! Now psi
        revID(53) = 0          ! No longer defined
        revID(54) = 0          ! No longer defined
     end if
     
  end if !if(changeVar.ge.1)
  
  
  
  ! *** Put plot data in startval and jumps.  Print initial and starting values to screen.
  ! Startval: 1: injection value, 2: starting value, 3: Lmax value
  jumps = 0.
  offsetrun = 0
  if(prInitial.ne.0) then
     if(htmlOutput.ge.1) then
        if(prInitial.eq.1) write(stdOut,'(/,A)')'  <b>Starting values for the chains:</b>'
        if(prInitial.eq.2) write(stdOut,'(/,A)')'  <b>Injection and starting values for the chains:</b>'
        if(prInitial.ge.3) write(stdOut,'(/,A)')'  <b>Injection, starting and Lmax values for the chains:</b>'
        write(stdOut,'(15x,A3,A10)',advance="no") '<b>','log Post'
     else
        if(prInitial.eq.1) write(stdOut,'(/,A)')'  Starting values for the chains:'
        if(prInitial.eq.2) write(stdOut,'(/,A)')'  Injection and starting values for the chains:'
        if(prInitial.ge.3) write(stdOut,'(/,A)')'  Injection, starting and Lmax values for the chains:'
        write(stdOut,'(15x,A10)',advance="no") 'log Post'
     end if
     do p=1,nMCMCpar
        write(stdOut,'(A8)',advance="no") trim(parNames(parID(p)))
     end do
     if(htmlOutput.ge.1) then
        write(stdOut,'(A4)') '</b>'
     else
        write(stdOut,*)
     end if
  end if
  
  startval = 0.
  do ic=1,nchains
     do i=1,2
        if(nint(is(ic,i)).eq.-1) startval(ic,1:nMCMCpar,1)  = allDat(ic,1:nMCMCpar,i)  ! Injection value
        if(nint(is(ic,i)).eq.0)  startval(ic,1:nMCMCpar,2)  = allDat(ic,1:nMCMCpar,i)  ! Starting value
     end do
     
     startval(ic,1:nMCMCpar,3)    = allDat(icloglmax,1:nMCMCpar,iloglmax)  ! Lmax value
     jumps(ic,1:nMCMCpar,2:n(ic)) = allDat(ic,1:nMCMCpar,2:n(ic)) -  allDat(ic,1:nMCMCpar,1:n(ic)-1)
     if(prInitial.ne.0) then 
        if(ic.eq.1.and.prInitial.ge.2) then
           write(stdOut,'(4x,A11)',advance="no") 'Injection:  '
           write(stdOut,'(F10.3)',advance="no") post(ic,1)
           do p=1,nMCMCpar
              write(stdOut,'(F8.3)',advance="no") startval(1,p,1)
           end do
           write(stdOut,*)
           if(prInitial.ge.4) write(stdOut,*)
        end if
        if( abs( (sum(startval(ic,1:nMCMCpar,1))-sum(startval(ic,1:nMCMCpar,2)) ) / &
             (sum(startval(ic,1:nMCMCpar,1))+sqrt(tiny(startval)) ) ) .gt. 1.e-10) then
           offsetrun = 1
           write(stdOut,'(I4,A1,A10)',advance="no")ic,':','  Start: '
           write(stdOut,'(F10.3)',advance="no")post(ic,2)
           do p=1,nMCMCpar
              write(stdOut,'(F8.3)',advance="no")startval(ic,p,2)
           end do
           write(stdOut,*)
           if(prInitial.ge.4) then
              write(stdOut,'(5x,A10)',advance="no")'Diff:  '
              write(stdOut,'(F10.3)',advance="no")abs(post(ic,1)-post(ic,2))
              do p=1,nMCMCpar
                 write(stdOut,'(F8.3)',advance="no")abs(startval(ic,p,1)-startval(ic,p,2))
              end do
              write(stdOut,'(/)')
           end if
        end if
     end if
  end do
  
  ! If the injection or starting values are not found, don't plot them:
  if(abs(sum(startval(:,1:nMCMCpar,1))).lt.1.e-10) then
     if(plInject.ne.0.and.prProgress.ge.2) write(*,'(A)') 'No injection values found, setting plInject to zero'
     plInject = 0
  end if
  if(abs(sum(startval(:,1:nMCMCpar,2))).lt.1.e-10) then
     if(plStart.ne.0.and.prProgress.ge.2) write(*,'(A)') 'No starting values found, setting plStart to zero'
     plStart  = 0
  end if
  
  if(prInitial.ge.3) then
     write(stdOut,'(5x,A10)',advance="no")'Lmax:  '
     write(stdOut,'(F10.3)',advance="no")post(icloglmax,iloglmax)
     do p=1,nMCMCpar
        write(stdOut,'(F8.3)',advance="no")startval(1,p,3)
     end do
     write(stdOut,*)
     if(prInitial.ge.4) then
        do ic=1,1 !nchains
           write(stdOut,'(I4,A1,A10)',advance="no")ic,':','Diff:  '
           write(stdOut,'(F10.3)',advance="no")abs(post(ic,1)-post(icloglmax,iloglmax))
           do p=1,nMCMCpar
              write(stdOut,'(F8.3)',advance="no")abs(startval(ic,p,1)-startval(ic,p,3))
           end do
           write(stdOut,*)
        end do
     end if
     write(stdOut,*)
  end if
  
  
  !if(prProgress.ge.1.and.update.eq.0) write(stdOut,'(A)')'  Done.'
  if(prProgress.ge.2.and.update.eq.0) write(stdOut,'(A,I12,A,F7.4)')'  t0:',nint(t0), '  GMST:',gmst(t0)
  
  
  !Check which parameters were fixed during the MCMC run:
  fixedpar = 0
  do ic=1,nchains
     do p=1,nMCMCpar
        !Doesn't matter in which chain this happens:
        if( abs(minval(allDat(ic,p,5:n(ic))) - maxval(allDat(ic,p,5:n(ic))) ) .lt. 1.d-6) fixedpar(p) = 1
     end do
  end do
  nfixedpar = sum(fixedpar)
  
  !Is this a 'spinning run' or not?  -  spinningRun: 0-no, 1-one spin, 2-two spins
  spinningRun = 0
  if(revID(71).ne.0) then
     if(fixedpar(revID(71)).eq.0) spinningRun = spinningRun + 1
  end if
  if(revID(81).ne.0) then
     if(fixedpar(revID(81)).eq.0) spinningRun = spinningRun + 1
  end if
  
  
  !*** Construct output file name:  GPS0929052945_H1L1V__Apostolatos_1.5pN_SP  for GPS time, detectors, waveform-pN Spinning
  ic = 1
  write(outputname,'(A3,I10.10,A1)') 'GPS',GPStime,'_'
  do i=1,ndet(ic)
     write(outputname,'(A)') trim(outputname)//trim(detabbrs(detnr(ic,i)))
  end do
  write(outputname,'(A,F3.1,A)') trim(outputname)//'__'//trim(waveformName)//'_',pnOrder,'pN'
  if(spinningRun.gt.0) then
     write(outputname,'(A,I1,A)') trim(outputname)//'_',spinningRun,'sp'
  else
     write(outputname,'(A)') trim(outputname)//'_ns'
  end if
  
  
  
  
  !*** Put data in selDat
  if(mergeChains.eq.1) then  !Merge chains, leave out burn-in (then nchains = 1)
     allocate(selDat(1,maxMCMCpar,nchains*maxLine))
     j = 1
     do ic=1,nchains
        do i=Nburn(ic)+1,n(ic)
           !selDat has the same structure as allDat, but contains only data AFTER the burn-in:
           selDat(1,1:nMCMCpar,j) = allDat(ic,1:nMCMCpar,i)  
           j = j+1
        end do
     end do
     nchains = 1
     n(1) = j-1
     !if(prProgress.ge.1) write(stdOut,'(A,I8,A,ES7.1,A)')'  Data points in combined chains: ',n(1),'  (',real(n(1)),')'
  else
     allocate(selDat(nchains,maxMCMCpar,maxLine))
     do ic=1,nchains
        !SelDat has the same structure as allDat, but contains only info AFTER the burn-in:
        selDat(ic,1:nMCMCpar,1:n(ic)-Nburn(ic)) = allDat(ic,1:nMCMCpar,Nburn(ic)+1:n(ic))
        n(ic) = n(ic)-Nburn(ic) !n(ic)=0 if a chain is not contributing (in which case contrChain(ic)=0)!
     end do
     !if(prProgress.ge.1) write(stdOut,'(A,I8)')' Datapoints in combined chains: ',sum(n(1:nchains))
  end if
  
end subroutine mcmcruninfo
!***********************************************************************************************************************************



!***********************************************************************************************************************************
!> \brief  Save post-burnin data to file (__data.txt)
!!
!! \retval exitcode  Exit status code (0=ok)

subroutine save_data(exitcode)
  use general_data, only: outputname,outputdir,n, selDat
  use mcmcrun_data, only: nMCMCpar,parID
  
  implicit none
  integer, intent(out) :: exitcode
  integer :: i,o,p
  
  exitcode = 0
  o = 20 !Output port
  open(unit=o, form='formatted', status='replace',file=trim(outputdir)//'/'//trim(outputname)//'__data.txt')
  
  do p=1,nMCMCpar
     if(parID(p).eq.61) write(o,'(A10)',advance="no")'mchirp'
     if(parID(p).eq.62) write(o,'(A10)',advance="no")'eta'
     if(parID(p).eq.67) write(o,'(A10)',advance="no")'q'
     if(parID(p).eq.68) write(o,'(A10)',advance="no")'log(q)'
     if(parID(p).eq.11.or.parID(p).eq.12) write(o,'(A10)',advance="no")'time'
     if(parID(p).eq.22) write(o,'(A10)',advance="no")'log(dist)'
     if(parID(p).eq.31) write(o,'(A10)',advance="no")'RA'
     if(parID(p).eq.32) write(o,'(A10)',advance="no")'sin(dec)'
     if(parID(p).eq.51) write(o,'(A10)',advance="no")'cos(iota)'
     if(parID(p).eq.41) write(o,'(A10)',advance="no")'phi0'
     if(parID(p).eq.52) write(o,'(A10)',advance="no")'psi'
     !write(o,'(A10)',advance="no")trim(parNames(parID(p)))
  end do
  write(o,*)
  !  do p=1,nMCMCpar
  !    write(o,'(F10.5)',advance="no")startval(1,p,1)
  !  end do
  !  write(o,*)
  !  write(o,*)
  
  !  do ic=1,nChains0
  do i=1,n(1)
     do p=1,nMCMCpar
        write(o,'(F10.5)',advance="no")selDat(1,p,i)
     end do
     write(o,*)
  end do
  !  end do
  
  close(o)  ! data output file
  
end subroutine save_data
!***********************************************************************************************************************************



!***********************************************************************************************************************************
!> \brief  Get injection values for LALInferenceMCMC (LIM) output
!!
!! \param ic  Current chain/input file
!! \param startval  Injection/Staring/Lmax values

subroutine get_LIM_injection_values(ic, nMCMCpar, startval, post,prior)
  use SUFR_kinds, only: double
  use SUFR_system, only: warn,find_free_io_unit
  use analysemcmc_settings, only: maxChs,maxMCMCpar
  use general_data, only: infiles
  use mcmcrun_data, only: ndet,parID
  
  implicit none
  integer, intent(in) :: ic, nMCMCpar
  real, intent(inout) :: startval(maxChs,maxMCMCpar,3)
  real, intent(out) :: post,prior
  
  integer :: ip,in,status,line, tmpInt, p
  real(double) :: tmpDat(maxMCMCpar),dtmpDat(maxMCMCpar)
  character :: infile*(99), tmpStr*(99)
  logical :: ex
  
  startval(ic,:,1) = 0.
  
  infile = infiles(ic)
  in = index(trim(infile),'.', back=.true.)
  if(len_trim(infile)-in.ne.2) return        ! Third-last character in file name is not a period
  
  infile = infile(1:in)//'injection'         ! Compose injection file name from MCMC output file name
  inquire(file=trim(infile), exist=ex)       ! Check whether the file exists
  if(.not.ex) return
  
  call find_free_io_unit(ip)
  open(unit=ip,form='formatted',status='old',position='rewind',file=trim(infile),iostat=status)
  if(status.ne.0) then
     call warn('File not found: '//trim(infile), 0)
     return
  end if
  
  do line=1,ndet(ic)+6
     read(ip,*) tmpStr
  end do
  
  read(ip,*,iostat=status) tmpInt,post,prior,tmpDat(1:nMCMCpar)        ! Don't have injection posterior (first tmpReal)...
  close(ip)
  
  ! Get rid of first zillion decimals of GPS time, so that the number can be saved in single-precision:
  dtmpDat = 0.d0
  do p=1,nMCMCpar
     if(parID(p).ge.11.and.parID(p).le.19) &
          dtmpDat(p) = dble(nint(tmpDat(p)/10.d0)*10)  ! 'GPS base time', rounded off to the nearest 10s
  end do
  
  startval(ic,1:nMCMCpar,1) = real(tmpDat(1:nMCMCpar) - dtmpDat(1:nMCMCpar))
  
  ! Set but never used:
  tmpInt = tmpInt
  tmpStr = tmpStr
  
end subroutine get_LIM_injection_values
!***********************************************************************************************************************************



!***********************************************************************************************************************************
subroutine error_maxMCMCpar_quit(op, maxMCMCpar, nMCMCpar)
  implicit none
  integer, intent(in) :: op, maxMCMCpar, nMCMCpar
  
  write(op,'(//,A,I4,A,I4,A,//)')'  Error:  maxMCMCpar too small.  You must increase maxMCMCpar from',maxMCMCpar, &
       ' to at least',nMCMCpar,' in order to continue.  Aborting...'
  
  stop 1
end subroutine error_maxMCMCpar_quit
!***********************************************************************************************************************************
