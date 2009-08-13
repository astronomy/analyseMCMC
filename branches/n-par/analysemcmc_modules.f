!Modules for analysemcmc


!***************************************************************************************************
!> Module with settings from the input file (e.g. analysemcmc.dat)
!<
module analysemcmc_settings
  implicit none
  save
  integer, parameter :: nival1=9
  integer, parameter :: maxChs=20     !< Maximum number of chains that can be read
  integer, parameter :: maxMCMCpar=18 !< MaxMCMCpar: logL+MCMCpar+secondary variables, e.g. 1+12+2(M1M2) = 15 for 12 par; 18 for 15 par
  integer, parameter :: nParDB=99     !< nParDB: size of the parameter database
  integer :: plvars(maxMCMCpar),nplvar,nbin1d,nbin2dx,nbin2dy,npdf2d,pdf2dpairs(250,2),panels(2)
  integer :: thin,nburn(maxChs),reverseread,update,mergechains,wrapdata,changevar,maxchlen
  integer :: file,colour,orientation,quality,fonttype
  integer :: prprogress,prruninfo,prchaininfo,prinitial,prstat,prcorr,prival,prconv,savestats,savepdf       
  integer :: plot,combinechainplots,pllogl,plchain,plparl,pljump,plpdf1d,plpdf2d,placorr,plotsky,plmovie       
  integer :: chainsymbol,chainpli,pltrue,plstart,plmedian,plrange,plburn,pllmax,prvalues,smooth,fillpdf,normpdf1d,normpdf2d
  integer :: scloglpl,scchainspl,bmpxsz,bmpysz,map_projection
  integer :: nmovframes,moviescheme,whitebg,unsharp,nival,ival0,wikioutput
  real :: nburnfrac,autoburnin,ivals(nival1)
  real :: scrsz,scrrat,pssz,psrat,scfac,fontsize1d,fontsize2d
end module analysemcmc_settings
!***************************************************************************************************

!***************************************************************************************************
!> Module with (currently) mathematical and string constants
!<
module constants
  implicit none
  save
  integer :: os
  real*8 :: pi,tpi,pi2,r2d,d2r,r2h,h2r,c3rd
  real :: rpi,rtpi,rpi2,rr2d,rd2r,rr2h,rh2r,rc3rd
  character :: upline*4,detabbrs(4)*2,waveforms(4)*99
end module constants
!***************************************************************************************************

!***************************************************************************************************
!> Module with Markov-chain data from the SPINspiral output files
!< 
module general_data
  use analysemcmc_settings
  implicit none
  save
  integer, parameter :: maxIter=1.01e5+2,nr1=5,nstat1=10,ndets=3
  integer :: n(maxChs),ntot(maxChs),iloglmax,icloglmax,c0,nchains,nchains0
  integer :: fixedpar(maxMCMCpar),nfixedpar,contrchains,contrchain(maxChs)
  real, allocatable :: selDat(:,:,:),allDat(:,:,:),post(:,:),prior(:,:)
  real :: startval(maxChs,maxMCMCpar,3)
  real :: ranges(maxChs,nival1,maxMCMCpar,nr1),stats(maxChs,maxMCMCpar,nstat1),log10bayesfactor(maxChs),logebayesfactor(maxChs)
  real*8 :: rhat(maxMCMCpar)
  
  character :: varnames(nParDB)*8,infiles(maxChs)*99,outputname*99,outputdir*99
  character :: pgunits(nParDB)*99,pgvarns(nParDB)*99,pgvarnss(nParDB)*99,pgorigvarns(nParDB)*99
  
  integer :: wrap(maxChs,maxMCMCpar)
  real :: rashift,racentre,shift(maxChs,maxMCMCpar)
end module general_data
!***************************************************************************************************

!***************************************************************************************************
!> Module with MCMC run data from the SPINspiral output files
!< 
module mcmcrun_data
  use analysemcmc_settings
  use general_data
  implicit none
  save
  integer :: niter(maxChs),totiter,totlines,totpts,nburn0(maxChs),seed(maxChs),ndet(maxChs),totthin(maxChs)
  integer :: nCorr(maxChs),nTemps(maxChs),waveform,nMCMCpar,Tmax(maxChs)
  integer :: samplerate(maxChs,ndets),samplesize(maxChs,ndets),FTsize(maxChs,ndets),detnr(maxChs,ndets),offsetrun
  integer :: parID(maxMCMCpar),revID(nParDB),spinningRun
  integer*8 :: GPStime
  real :: snr(maxChs,ndets),flow(maxChs,ndets),fhigh(maxChs,ndets),t_before(maxChs,ndets),t_after(maxChs,ndets),deltaFT(maxChs,ndets)
  real :: Tchain(maxChs),networkSNR(maxChs),pnOrder
  real*8 :: FTstart(maxChs,ndets),t0,loglmax,loglmaxs(maxChs)
  character :: detnames(maxChs,ndets)*14
end module mcmcrun_data
!***************************************************************************************************


!***************************************************************************************************
!> Module with generated statistics
!< 
module stats_data
  use analysemcmc_settings
  implicit none
  save
  real :: stdev1(maxMCMCpar),stdev2(maxMCMCpar),absvar1(maxMCMCpar),absvar2(maxMCMCpar)
  integer :: trueranges2d(maxMCMCpar,maxMCMCpar)
  real :: probarea(nival1),probareas(maxMCMCpar,maxMCMCpar,nival1,3)
end module stats_data
!***************************************************************************************************


!***************************************************************************************************
!> Module with plot settings
!< 
module plot_data
  use analysemcmc_settings
  implicit none
  save
  integer :: ncolours,colours(10),defcolour,nsymbols,symbols(10),maxdots
  real :: pltrat,bmpsz,bmprat
  character :: bmpxpix*99,unsharplogl*99,unsharpchain*99,unsharppdf1d*99,unsharppdf2d*99
  character :: psclr*9,colournames(15)*20
end module plot_data
!***************************************************************************************************


!***************************************************************************************************
!> Module with secondary Markov-chain data
!< 
module chain_data
  use general_data
  implicit none
  save
  real :: is(maxChs,maxIter),isburn(maxChs)
  real :: jumps(maxChs,maxMCMCpar,maxIter)
  real :: corrs(maxMCMCpar,maxMCMCpar),acorrs(maxChs,0:maxMCMCpar,0:maxIter)
  real*8 :: nullh
end module chain_data
!***************************************************************************************************

