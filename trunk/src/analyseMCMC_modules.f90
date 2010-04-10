!> \file analyseMCMC_modules.f90
!! \brief Modules used by analyseMCMC
!<



!***********************************************************************************************************************************
!> \brief Module with data types, etc.
!<
module basic
   implicit none
   
   ! Determine real double-precision kind for the current compiler/OS:
   integer, parameter :: double = selected_real_kind(15,307),  dbl = double
   
   ! Determine integer double-precision kind for the current compiler/OS:
   integer, parameter :: long = selected_int_kind(18),  lng = long
   
   ! Determine maximum-precision kinds for the current compiler/OS:
   !   selected_int/real_kind return negative integers for unsupported precision, and supposedly larger values for supported kinds
   !   with larger precision.  Int 9,18 are standard int,long; real 6,15 are standard single,double precision.  Larger precisions
   !   used are found in some compilers, except 99, which I made up.
   integer, parameter :: maxint = max(selected_int_kind(9),selected_int_kind(18),selected_int_kind(38),selected_int_kind(99)) 
   integer, parameter :: maxreal = max(selected_real_kind(6),selected_real_kind(15),selected_real_kind(18), &
        selected_real_kind(31),selected_real_kind(33),selected_real_kind(99))
   
   integer :: stdOut,stdErr
end module basic
!***********************************************************************************************************************************

!***********************************************************************************************************************************
!> \brief Module with settings from the input file (e.g. analysemcmc.dat)
!<
module analysemcmc_settings
  implicit none
  save
  integer, parameter :: maxNival=5    !< maxNival: Maximum number of probability intervals that can be used
  integer, parameter :: maxChs=25     !< macChs: Maximum number of chains that can be read
  integer, parameter :: maxMCMCpar=17 !< MaxMCMCpar: MCMCpar+secondary parameters, e.g. 12+2(M1M2) = 14 for 12 par; 17 for 15 par
  integer, parameter :: nParDB=199     !< nParDB: size of the parameter database
  integer :: plPars(maxMCMCpar),nPlPar,Nbin1D,Nbin2Dx,Nbin2Dy,Npdf2D,PDF2Dpairs(250,2),panels(2)
  integer :: thin,Nburn(maxChs),reverseRead,update,mergeChains,wrapData,changeVar,maxChLen
  integer :: file,colour,orientation,quality,fonttype
  integer :: prStdOut,prProgress,prRunInfo,prChainInfo,prInitial,prStat,prCorr,prAcorr,nAcorr,prIval,prConv
  integer :: saveStats,savePDF,tailoredOutput
  integer :: plot,plLogL,plChain,plParL,plJump,plPDF1D,plPDF2D,plACorr,plotSky,plAnim       
  integer :: chainSymbol,chainPlI,plInject,plStart,plMedian,plRange,plBurn,plLmax,prValues,smooth,fillPDF,normPDF1D,normPDF2D
  integer :: scLogLpl,scChainsPl,bmpXSz,bmpYSz,map_projection
  integer :: nAnimFrames,animScheme,whiteBG,unSharp,Nival,ival0,wikioutput,html
  real :: NburnFrac,autoBurnin,ivals(maxNival)
  real :: scrSz,scrRat,PSsz,PSrat,scFac,fontsize1d,fontsize2d
end module analysemcmc_settings
!***********************************************************************************************************************************

!***********************************************************************************************************************************
!> \brief Module with (currently) mathematical and string constants
!<
module constants
  use basic
  
  implicit none
  save
  integer :: os
  real(double) :: pi,tpi,pi2,r2d,d2r,r2h,h2r,c3rd
  real :: rpi,rtpi,rpi2,rr2d,rd2r,rr2h,rh2r,rc3rd
  character :: upline*4,detabbrs(4)*2,waveforms(0:9)*99,currentdatestr*19,currenttimestr*19,currenttimezonestr*19
  character(len=99) :: homedir,workdir,hostname,username,stdOutFile
end module constants
!***********************************************************************************************************************************

!***********************************************************************************************************************************
!> \brief Module with Markov-chain data from the SPINspiral output files
!< 
module general_data
  use basic
  use analysemcmc_settings
  
  implicit none
  save
  integer, parameter :: maxIter=150000                    !< maxIter: Maximum number of iterations (output lines) that can be stored
  integer, parameter :: nr1=5,nstat1=10,ndets=4
  integer :: n(maxChs),ntot(maxChs),iloglmax,icloglmax,c0,nchains,nchains0
  integer :: fixedpar(maxMCMCpar),nfixedpar,contrchains,contrchain(maxChs)
  real, allocatable :: selDat(:,:,:),allDat(:,:,:),post(:,:),prior(:,:)
  real :: startval(maxChs,maxMCMCpar,3)
  real :: ranges(maxChs,maxNival,maxMCMCpar,nr1),stats(maxChs,maxMCMCpar,nstat1)
  real :: log10bayesfactor(maxChs),logebayesfactor(maxChs),logebayesfactortotalharmo,logebayesfactortotalgeom
  real :: logebayesfactortotalarith,logebayesfactortotal,logebayestempfactor(maxChs)
  real(double) :: rhat(maxMCMCpar)
  
  character :: parNames(nParDB)*8,infiles(maxChs)*99,outputname*99,outputdir*99
  character :: pgunits(nParDB)*99,pgParNs(nParDB)*99,pgParNss(nParDB)*99,pgOrigParns(nParDB)*99
  
  integer :: wrap(maxChs,maxMCMCpar)
  real :: raShift,raCentre,shifts(maxChs,maxMCMCpar),shIvals(maxChs,maxMCMCpar)
end module general_data
!***********************************************************************************************************************************

!***********************************************************************************************************************************
!> \brief Module with MCMC run data from the SPINspiral output files
!< 
module mcmcrun_data
  use basic
  use analysemcmc_settings
  use general_data
  
  implicit none
  save
  integer :: niter(maxChs),totiter,totlines,totpts,Nburn0(maxChs),seed(maxChs),ndet(maxChs),totthin(maxChs)
  integer :: nCorr(maxChs),nTemps(maxChs),waveform,nMCMCpar,nMCMCpar0,Tmax(maxChs)
  integer :: samplerate(maxChs,ndets),samplesize(maxChs,ndets),FTsize(maxChs,ndets),detnr(maxChs,ndets),offsetrun
  integer :: parID(maxMCMCpar),revID(nParDB),spinningRun
  integer(long) :: GPStime
  real :: snr(maxChs,ndets),flow(maxChs,ndets),fhigh(maxChs,ndets),t_before(maxChs,ndets),t_after(maxChs,ndets)
  real :: deltaFT(maxChs,ndets)
  real :: Tchain(maxChs),networkSNR(maxChs),pnOrder,outputVersion
  real(double) :: FTstart(maxChs,ndets),t0,loglmax,loglmaxs(maxChs)
  character :: detnames(maxChs,ndets)*14
end module mcmcrun_data
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief Module with generated statistics
!< 
module stats_data
  use analysemcmc_settings
  implicit none
  save
  real :: stdev1(maxMCMCpar),stdev2(maxMCMCpar),absVar1(maxMCMCpar),absVar2(maxMCMCpar)
  integer :: injectionranges2d(maxMCMCpar,maxMCMCpar)
  real :: probarea(maxNival),probareas(maxMCMCpar,maxMCMCpar,maxNival,3)
end module stats_data
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief Module with plot settings
!< 
module plot_data
  use analysemcmc_settings
  implicit none
  save
  integer :: ncolours,colours(10),defcolour,nsymbols,symbols(10),maxdots
  real :: pltsz,pltrat,bmpsz,bmprat
  character :: bmpxpix*99,unSharplogl*99,unSharpchain*99,unSharppdf1d*99,unSharppdf2d*99
  character :: psclr*9,colournames(15)*20
end module plot_data
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief Module with secondary Markov-chain data
!< 
module chain_data
  use general_data
  implicit none
  save
  real :: is(maxChs,maxIter),isburn(maxChs)
  real :: jumps(maxChs,maxMCMCpar,maxIter)
  real :: corrs(maxMCMCpar,maxMCMCpar),acorrs(maxChs,0:maxMCMCpar,0:maxIter),lAcorrs(maxChs,0:maxMCMCpar)
  real(double) :: DoverD
end module chain_data
!***********************************************************************************************************************************

