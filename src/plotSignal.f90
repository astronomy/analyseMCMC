!> \file plotSignal.f90  Read and plot the signal output from SPINspiral

! 
! LICENCE:
! 
! Copyright (c) 2007-2014  Marc van der Sluys
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
program plotsignal
  use SUFR_kinds, only: double
  use SUFR_numerics, only: seq,sne
  
  implicit none
  integer, parameter :: n1=1000000, nf0=5
  
  integer :: n(nf0),i,j,io,pgopen,file,f,prname,prtitle,system,thin
  real(double) :: t0,t1(nf0,n1)
  real :: t(nf0,n1),h(nf0,n1),dx,dy,xmin,xmax,ymin,ymax
  character :: title*(1000),dirname*(99),subdir*(99),fnames(nf0)*(99),fname*(199),bla*(199),t0s*(10),detnames(nf0)*(99)
  
  integer :: nf, nfrx,nfry,frx,fry,fr,colours(nf0), lw
  real :: xwinmin,xwinmax,ywinmin,ywinmax,dxwin,dywin,xfrmin,xfrmax,yfrmin,yfrmax,sch
  real(double) :: m1,m2, dl,spin !,mc,eta,tc,lat,lon,phase,kappa,thJ0,phJ0,alpha
  character :: m1s*(20),m2s*(20),dls*(20),spins*(20)
  logical :: inOnePanel
  
  file = 1 !0-screen, 1-file: eps, 2-file: pdf, 3-file: png
  thin = 1 !Thin the number of points by this factor at read
  prname = 0 !Print detector name in frame: 0-no, 1-yes
  prtitle = 0 !Print plot title: 0-no, 1-yes
  inOnePanel = .true.  ! Make all plots in the same panel
  
  colours = 2
  !colours = (/1,1/)
  !colours = (/4,2/)
  !colours = (/4,3,2/)
  colours(1:4) = (/2,4,6,3/)
  !colours = (/4,2,4,2/)
  
  ! Directory that contains the signals
  dirname = '/home/sluys/work/GW/programs/MCMC/SPINspiral/signals/'
  
  !fnames = (/'Hanford-signal.dat','Livingston-signal.dat','Pisa-signal.dat'/)
  !detnames = (/'LIGO/Hanford','LIGO/Livingston','VIRGO/Pisa'/)
  !fnames = (/'Hanford-signal.dat','Pisa-signal.dat'/)
  !fnames = (/'ChrisB/Hanford_a0.0_th35.dat','ChrisB/Pisa_a0.0_th35.dat'/)
  !detnames = (/'LIGO','VIRGO'/)
  
  !fnames = (/'Hanford-signal.dat','Livingston-signal.dat'/)
  !detnames = (/'LIGO/Hanford','LIGO/Livingston'/)
  !fnames = (/'Hanford-signal_run003.dat','Hanford-signal_run001.dat','bla'/)
  !fnames = (/'Hanford-signal_run003.dat','Hanford-signal_run001.dat'/)
  !detnames = (/'\(2134)\dSL\u=55\(2218)','\(2134)\dSL\u=25\(2218)','bla'/)
  !fnames = (/'Hanford-run002.dat','Hanford-run004.dat','bla'/)
  !detnames = (/'a\dspin\u = 0.1','a\dspin\u = 0.8','bla'/)
  !fnames = (/'Hanford_0.10_055.dat','Hanford_0.80_055.dat'/)
  !fnames = (/'../signal11H.dat','../signal12H.dat'/)
  !fnames = (/'../Hanford-signal.dat'/)
  !fnames = (/'../signal_00.dat','../signal_de.dat','../signal_ra.dat'/)
  !fnames = (/'../signal11L.dat','../signal12L.dat'/)
  !fnames = (/'Hanford-signal0.dat','Hanford-signal5.dat','Hanford-signal8.dat','Hanford-signal9.dat'/)
  !fnames = (/'../Hanford-signal.dat'/)
  !fnames = (/'../Hanford-signal0.dat','../Hanford-signal3.dat','../Hanford-signal0.dat','../Hanford-signal3.dat'/)
  !fnames = (/'Hanford-signal.dat','Hanford-data.dat','Hanford-injection.dat'/)
  !fnames = (/'Hanford-signal_run001.dat','Hanford-signal_run003.dat'/)
  !fnames = (/'../Hanford-data.dat','../Livingston-data.dat','../Pisa-data.dat'/)
  !fnames = (/'../Hanford-signal.dat','../Livingston-signal.dat','../Pisa-signal.dat'/)
  !fnames = (/'../Hanford-signal.dat','../Hanford-template.dat','../Pisa-signal.dat'/)
  !fnames = (/'Hanford_a0.0_th20.dat','Hanford_a0.1_th20.dat','Hanford_a0.5_th20.dat'/)
  !fnames = (/'Hanford_a0.0_th55.dat','Hanford_a0.1_th55.dat','Hanford_a0.5_th55.dat'/)
  !fnames = (/'Hanford_a0.5_th20.dat','Hanford_a0.5_th55.dat'/)
  !detnames = (/'\(0634)=20\(2218)','\(0634)=55\(2218)'/)
  !fnames = (/'~/work/GW/programs/MCMC/mcmc_code/trunk/Hanford-data.dat', &
  !'~/work/GW/programs/MCMC/backup/spinning.backup/Hanford-data.dat'/)
  !detnames = (/'new','old'/)
  
  !For GWDAW poster:
  !fnames = (/'0.00_020-signal.dat', '0.10_020-signal.dat', '0.50_020-signal.dat'/)
  !detnames = (/'','',''/)
  
  !! Methods paper 1:
  !nf = 3
  !subdir = 'Apo03/'
  !!fnames = (/'0.00_020-signal.dat', '0.50_020-signal.dat', '1.00_020-signal.dat'/)  ! Apo01
  !fnames = (/'0.00_050-signal.dat', '0.50_050-signal.dat', '1.00_050-signal.dat'/)  ! Apo02-03
  !detnames = (/'','',''/)
  
  ! Methods paper 2:
  !nf = 3
  !subdir = 'ST01/'
  !fnames(1:nf) = (/'NS_NS-signal.dat', 'HL_NS-signal.dat', 'HL_HL-signal.dat'/)  ! ST01
  !detnames(1:nf = (/'0 spins','1 spin ','2 spins'/)
  !nf = 4
  !fnames(1:nf) = (/'NS_NS-signal.dat', 'HL_NS-signal.dat', 'NS_HL-signal.dat', 'HL_HL-signal.dat'/)  ! ST0
  !detnames(1:nf) = (/'','','',''/)

  ! Methods paper 3:
  !nf = 2
  !subdir = 'ST01/'
  !!fnames(1:nf) = (/'1.5pN-signal.dat', '3.5pN-signal.dat'/)  ! ST01
  !fnames(1:nf) = (/'1.5pN_ns-signal.dat', '3.5pN_ns-signal.dat'/)  ! ST01
  !detnames(1:nf) = (/'',''/)
  nf = 1
  subdir = 'ST01/'
  !fnames(1:nf) = (/'1.5pN-signal.dat', '3.5pN-signal.dat'/)  ! ST01
  fnames(1:nf) = (/'3.5pN_ns1-signal.dat'/)  ! ST01 - eta=0.04
  detnames(1:nf) = (/''/)
  
  
  !fnames = (/'0.00_020-signal.dat'/)
  !detnames = (/''/)
  
  !fnames = (/'Hanford-signal-NS.dat', 'Hanford-signal-Sp.dat'/)
  !detnames = (/'GeneratePPN','SpinTaylor '/)
  
  !fnames = (/'Hanford-signal-NS.dat'/)
  !detnames = (/''/)
  
  !detnames = (/'Ref','Dec','RA'/)
  !detnames = (/'a\dspin\u = 0.0','a\dspin\u = 0.1','a\dspin\u = 0.5'/)
  !fnames = (/'../Pisa-data.dat'/)
  !detnames = (/'a\dspin\u=0.8, \(2134)\dSL\u=25\(2218)','a\dspin\u=0.1, \(2134)\dSL\u=25\(2218)'/)
  !detnames = (/'','','',''/)
  
  
  write(6,*)''
  do f=1,nf
     fname = trim(dirname)//trim(subdir)//trim(fnames(f))
     print*,trim(dirname),' - ',trim(subdir), ' - ', trim(fnames(f))
     open(unit=10,form='formatted',status='old',file=trim(fname),iostat=io)
     if(io.ne.0) then
        write(6,'(A,/)')'File not found: '//trim(fname)//'. Quitting the programme.'  
        stop
     end if
     rewind(10)
     
     write(6,'(A)') 'Reading input file '//trim(fname)//'...     '
     read(10,'(A)')bla
     print*,trim(bla)
     
     !read(10,*)m1,m2,mc,eta,tc,dl,lat,lon,phase,spin,kappa,thJ0,phJ0,alpha
     m1 = 0.
     m2 = 0.
     dl = 0.
     spin = 0.
     
     read(10,'(A)')bla
     print*,trim(bla)
     read(10,*)bla
     do i=1,n1
        read(10,*,err=195,end=199)t1(f,i),h(f,i)
        if(thin.gt.1) then
           do j=1,thin-1
              read(10,*,end=199)bla
           end do
        end if
     end do
     
195  write(6,'(A,I5)')'error reading file '//trim(fname)//' line ',i+1
199  close(10)
     n(f) = i-1
     write(6,'(I5,A12)')n(f),'lines read.'
  end do !f
  
  !write(6,'(6(A10,F8.3))')'m1:',m1,'m2:',m2,'mc:',mc,'eta:',eta,'dl:',dl,'spin:',spin
  !write(title,'(A,F5.1,A,F5.1,A,F4.1,A,F5.1,A)') 'M\d1\u =',m1,'M\d\(2281)\u,  M\d2\u =',m2,'M\d\(2281)\u,  a\dspin\u =',spin,', &
  !d\dL\u=',dl,'Mpc'
  write(m1s,'(F5.1)')m1
  if(m1.lt.10.) write(m1s,'(F4.1)')m1
  write(m2s,'(F5.1)')m2
  if(m2.lt.10.) write(m2s,'(F4.1)')m2
  write(dls,'(F5.1)')dl
  if(dl.lt.10.) write(dls,'(F4.1)')dl
  write(spins,'(F5.1)')spin
  if(spin.lt.10.) write(spins,'(F4.1)')spin
  !write(title,'(A)')'M\d1\u ='//trim(m1s)//'M\d\(2281)\u,  M\d2\u ='//trim(m2s)//'M\d\(2281)\u,  a\dspin\u ='//trim(spins)//', &
  !d\dL\u='//trim(dls)//'Mpc'
  write(title,'(A)')'M\d1\u ='//trim(m1s)//'M\d\(2281)\u,  M\d2\u ='//trim(m2s)//'M\d\(2281)\u,  d\dL\u='//trim(dls)//'Mpc'
  !write(title,'(A)')'M\d1\u ='//trim(m1s)//'M\d\(2281)\u,  M\d2\u ='//trim(m2s)//'M\d\(2281)\u'
  !print*,trim(title)
  !write(6,'(8(A10,F15.5))')'lat:',lat,'lon:',lon,'phase:',phase,'kappa:',kappa,'theta_Jo:',thJ0,'phi_Jo:',phJ0,'alpha:',alpha, &
  !'tc:',tc
  
  do f=1,nf
     do i=n(f),2,-1
        if(seq(h(f,i),0.) .and. sne(h(f,i-1),0.)) then
           !n(f) = i-1
           write(6,'(A20,F15.5)')'  '//trim(detnames(f))//':  ',t1(f,i-1)
           exit
        end if
     end do
  end do
  
  t0 = 1.d99
  do f=1,nf
     t0 = min(minval(t1(f,1:n(f))),t0)
  end do
  
  t0 = (nint(t0*1.d-3)-1)*1.d3
  write(t0s,'(I10)') nint(t0)
  write(6,*)''
  write(6,'(A)')'  t0: '//t0s
  do f=1,nf
     t(f,1:n(f)) = real(t1(f,1:n(f))-t0)
     !if(f.eq.2) t(f,1:n(f)) = real(t1(f,1:n(f))-t0 - 1.3935709e8) !Waveform 2 comes from a different run, with a different time
  end do
  
  
  
  do f=1,nf
     write(detnames(f),'(A,F3.1,A)')'a\dspin\u = ',spin,', '//trim(detnames(f))
  end do
  
  
  
  
  
  
  
  !*** PLOT ***
  
  !h = h*1.e22
  
  write(6,*)''
  if(1.eq.1) then
     select case(file)
        case(0)
           io = pgopen('12/xs')
           lw = 1
        case(1)
           io = pgopen('detectorsignal.eps/cps')
           lw = 1
        case(2)
           io = pgopen('detectorsignal.eps/cps')
           lw = 2
        case(3)
           io = pgopen('detectorsignal.ppm/ppm')
           lw = 1
     end select
     
     if(io.le.0) then
        write(0,'(A,I6,/)')'Cannot open PGPlot device.  Quitting the programme ',io
        stop
     end if
     call pgscf(1)
     !call pgpap(10.,min(max(0.25*real(nf),0.3),0.8))
     !call pgpap(10.,0.75)
     !call pgpap(15.,0.40)  !Poster plot with 3 waveforms
     if(inOnePanel) then
        call pgpap(15.,0.1333*8)     ! All in one panel
     else
        call pgpap(15.,0.1333*real(nf))  ! Poster plot with nf waveforms
     end if
     !call pgpap(20.,0.20)  !
     !sch = max(3./real(nf),1.5)
     sch = 6./real(nf)
     if(inOnePanel) sch = 6./8.
     !sch = 3.  !Poster plot with 3 waveforms
     call pgsch(sch)
     call pgscr(3,0.,0.6,0.)
     call pgslw(lw*2)
     
     xmin =  1.e30
     xmax = -1.e30
     ymin =  1.e30
     ymax = -1.e30
     do f=1,nf
        !xmin = min(minval(t(f,1:n(f))),xmin)
        !xmax = max(maxval(t(f,1:n(f))),xmax)
        do i=1,n(f)
           if(abs(h(f,i)).gt.1.e-33.and.t(f,i).lt.xmin) xmin = t(f,i)
           if(abs(h(f,i)).gt.1.e-33.and.t(f,i).gt.xmax) xmax = t(f,i)
        end do
        
        !ymin = min(minval(h(f,1:n(f))),ymin)
        !ymax = max(maxval(h(f,1:n(f))),ymax)
        ymax = maxval(abs(h))
        ymin = -ymax
     end do
     print*,'hmax:',maxval(abs(h))
     
     if(1.eq.1) then
        t = t - xmax
        xmin = xmin - xmax
        xmax = xmax - xmax
     end if
     
     !xmin  = 994.6
     !xmin  = 999.
     !xmin  = 1000.
     xmin = -1.
     !xmin = -0.2
     !xmax = 997.
     !xmax = 1000.1
     
     xmin = -5.0
     xmax = -3.0
     
     dx = abs(xmax-xmin)*0.01
     !dx = abs(xmax-xmin)*0.05
     dy = abs(ymax-ymin)*0.05
     !dy = abs(ymax-ymin)*0.01
     
     
     !do f=1,nf
     nfrx = 1      !Number of frames
     nfry = nf      
     !if(inOnePanel) nfry = 1
     xwinmin = 0.10    !Absolute limits
     xwinmax = 0.95
     !ywinmin = 0.10
     ywinmin = 0.20
     ywinmax = 0.95
     if(nf.eq.1) ywinmin = 0.25    
     if(nf.eq.1) ywinmax = 0.8
     if(nf.eq.2) ywinmax = 0.9
     if(nf.eq.3) ywinmax = 0.95
     !ywinmax = 0.5
     
     dxwin = (xwinmax-xwinmin)/real(nfrx)
     dywin = (ywinmax-ywinmin)/real(nfry)
     do fry = 1,nfry
        do frx = 1,nfrx
           fr = (fry-1)*nfrx + frx
           
           !  write(6,'(A)')"Changing ranges per frame !!!'
           !  if(fr.le.2) then
           !     xmin = 1006.8
           !     xmax = 1007.1
           !  end if
           !  if(fr.gt.2) then
           !     xmin = 1012.3
           !     xmax = 1012.35
           !  end if
           !  dx = abs(xmax-xmin)*0.01
           !  dy = abs(ymax-ymin)*0.05
           
           
           xfrmin = xwinmin+dxwin*real((frx-1))
           xfrmax = xwinmin+dxwin*real(frx)
           yfrmin = ywinmax-dywin*real(fry)
           yfrmax = ywinmax-dywin*real((fry-1))
           !print*,xfrmin,xfrmax,yfrmin,yfrmax
           
           call pgsci(1)
           if(fr.eq.1 .or. .not.inOnePanel) then
              call pgsvp(xfrmin,xfrmax,yfrmin,yfrmax)
              call pgswin(xmin-dx,xmax+dx,ymin-dy,ymax+dy)
              !call pgswin(xmin-dx,xmax+dx,ymin,ymax+2*dy)
           end if
           
           if(frx.eq.1) then
              !call pgmtxt('L',2.,0.5,0.5,'h (10\u-22\d)')
              !call pgmtxt('L',2.,0.5,0.5,'h/10\u-22\d')
              !call pgbox('',0.0,0,'BCNTSV',0.0,0)
              !call pgbox('',0.0,0,'BCNTS',0.0,0)
              call pgbox('',0.0,0,'BC',0.0,0)
           else
              call pgbox('',0.0,0,'BC',0.0,0)
           end if
           
           !Print binary parameters
           !  if(fry.eq.1) call pgmtxt('T',1.,0.5,0.5,'M\dbh\u=10M\d\(2281)\u, M\dns\u=1.4M\d\(2281)\u, a\dspin\u=0, d\dL\u=30Mpc')
           if(fry.eq.1.and.prtitle.eq.1) call pgmtxt('T',1.,0.5,0.5,trim(title))
           if(fry.eq.nfry) then
              !call pgmtxt('B',2.2,0.5,0.5,'GPStime - '//t0s//' (s)')
              call pgmtxt('B',2.5,0.5,0.5,'time (s)')
              !call pgmtxt('B',2.5,0.5,0.5,'tijd (s)')
              call pgbox('BCNTS',0.0,0,'',0.0,0)
           else
              call pgbox('BCTS',0.0,0,'',0.0,0)
           end if
           
           !Plot signal
           call pgsci(colours(fr))
           call pgslw(lw)
           call pgline(n(fr),t(fr,1:n(fr)),h(fr,1:n(fr)))
           !print*,frx,fry,fr
           call pgsch(1.)
           !  call pgpoint(n(fr),t(fr,1:n(fr)),h(fr,1:n(fr)),1)
           call pgsch(sch)
           call pgsci(1)
           call pgslw(lw*2)
           
           !Print detector name
           call pgsci(colours(fr))
           if(prname.eq.1) call pgmtxt('T',-2.0*ywinmax,0.04,0.,trim(detnames(fr)))
           call pgsci(1)
           
        end do  !frx
     end do  !fry
     !end do !f=1,nf
  end if !if(1.eq.2) then
  
  
  !Make gv autoreload 
  if(1.eq.2.and.file.eq.1) then
     call pgpage()
     call pgbox('BCTS',0.0,0,'BCNTS',0.0,0)
  end if
  call pgend()
  
  !Convert eps to pdf
  if(file.eq.2) then
     i = system('eps2pdf detectorsignal.eps >& /dev/null')
     i = system('rm -f detectorsignal.eps')
  end if
  !Convert ppm to png
  if(file.eq.3) then
     i = system('convert -depth 8 detectorsignal.ppm detectorsignal.png')
     i = system('rm -f detectorsignal.ppm')
  end if
  
  
  
  
  
  
  write(6,*)''
end program plotsignal
!***********************************************************************************************************************************



