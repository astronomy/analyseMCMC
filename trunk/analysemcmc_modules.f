!Modules for analysemcmc


!***************************************************************************************************
module analysemcmc_settings
  implicit none
  save
  integer, parameter :: nchs=20,npar1=18,nival1=9  !Npar1: logL+MCMCpar+secundary vars, e.g. 1+12+2(M1M2) = 15. For 15-par: 18?  
  integer :: plvars(npar1),nplvar,nbin1d,nbin2dx,nbin2dy,npdf2d,pdf2dpairs(250,2),panels(2)
  integer :: thin,nburn(nchs),file,colour,quality,reverseread,update,mergechains,wrapdata,changevar,maxchlen
  integer :: prprogress,prruninfo,prchaininfo,prinitial,prstat,prcorr,prival,prconv,savestats,savepdf       
  integer :: plot,combinechainplots,pllogl,plchain,plparl,pljump,rdsigacc,plsigacc,plpdf1d,plpdf2d,placorr,plotsky,plmovie       
  integer :: chainsymbol,chainpli,pltrue,plstart,plmedian,plrange,plburn,pllmax,prvalues,smooth,fillpdf,normpdf1d,normpdf2d
  integer :: scloglpl,scchainspl,bmpxsz,bmpysz
  integer :: nmovframes,moviescheme,whitebg,unsharp,nival,ival0
  real :: nburnfrac,autoburnin,ivals(nival1)
  real :: scrsz,scrrat,pssz,psrat,scfac
end module analysemcmc_settings
!***************************************************************************************************

!***************************************************************************************************
module constants
  implicit none
  save
  integer :: version
  real*8 :: pi,tpi,pi2,r2d,d2r,r2h,h2r,c3rd
  real :: rpi,rtpi,rpi2,rr2d,rd2r,rr2h,rh2r,rc3rd
  character :: upline*4,detabbrs(4)*2
end module constants
!***************************************************************************************************

!***************************************************************************************************
module general_data
  use analysemcmc_settings
  implicit none
  save
  !integer, parameter :: narr1=1.01e5+2,npar0=13,nr1=5,nstat1=10,ndets=3  !npar0: logL+MCMCpar; for 12-par
  integer, parameter :: narr1=1.01e5+2,npar0=16,nr1=5,nstat1=10,ndets=3  !npar0: logL+MCMCpar; for 15-par
  integer :: n(nchs),ntot(nchs),npar,iloglmax,icloglmax,c0,nchains,nchains0
  integer :: fixedpar(npar1),nfixedpar,contrchains,contrchain(nchs)
  integer :: par1,par2
  real, allocatable :: dat(:,:,:),alldat(:,:,:),pldat(:,:,:)
  real :: startval(nchs,npar1,3)
  real :: ranges(nchs,nival1,npar1,nr1),stats(nchs,npar1,nstat1),log10bayesfactor(nchs),logebayesfactor(nchs)
  real*8 :: rhat(npar1)
  
  character :: varnames(npar1)*8,infile*99,infiles(nchs)*99,outputname*99,outputdir*99
  character :: pgunits(npar1)*99,pgvarns(npar1)*99,pgvarnss(npar1)*99,pgorigvarns(npar1)*99
  
  integer :: wrap(nchs,npar1)
  real :: rashift,shift(nchs,npar1)
end module general_data
!***************************************************************************************************

!***************************************************************************************************
module mcmcrun_data
  use general_data
  implicit none
  save
  integer :: niter(nchs),totiter,totlines,totpts,nburn0(nchs),seed(nchs),ndet(nchs),totthin(nchs)
  integer :: samplerate(nchs,ndets),samplesize(nchs,ndets),FTsize(nchs,ndets),detnr(nchs,ndets),offsetrun
  real :: snr(nchs,ndets),flow(nchs,ndets),fhigh(nchs,ndets),t_before(nchs,ndets),t_after(nchs,ndets),deltaFT(nchs,ndets)
  real*8 :: FTstart(nchs,ndets),t0,loglmax,loglmaxs(nchs)
  character :: detnames(nchs,ndets)*14
end module mcmcrun_data
!***************************************************************************************************


!***************************************************************************************************
module stats_data
  use analysemcmc_settings
  implicit none
  save
  real :: stdev1(npar1),stdev2(npar1),absvar1(npar1),absvar2(npar1)
  integer :: trueranges2d(npar1,npar1)
  real :: probarea(nival1),probareas(npar1,npar1,nival1,3)
end module stats_data
!***************************************************************************************************


!***************************************************************************************************
module plot_data
  use analysemcmc_settings
  implicit none
  save
  integer :: ncolours,colours(10),defcolour,nsymbols,symbols(10),maxdots
  real :: pltrat
  real :: bmpsz,bmprat
  character :: bmpxpix*99,unsharplogl*99,unsharpchain*99,unsharppdf1d*99,unsharppdf2d*99
  character :: psclr*9,colournames(15)*20
end module plot_data
!***************************************************************************************************


!***************************************************************************************************
module chain_data
  use general_data
  implicit none
  save
  real :: is(nchs,narr1),isburn(nchs)
  real :: sig(npar1,nchs,narr1),acc(npar1,nchs,narr1),jumps(nchs,npar1,narr1)
  real :: corrs(npar1,npar1),acorrs(nchs,0:npar1,0:narr1)
  real*8 :: nullh
end module chain_data
!***************************************************************************************************

