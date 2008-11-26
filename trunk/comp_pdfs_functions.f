!Plot functions for comp_pdfs.f (1d, 2d PDFs and waveforms)


!***************************************************************************************************
module comp_pdfs_settings
  implicit none
  save
  integer, parameter :: nf1=9
  integer :: nf,file,type,dim,clr,fs,frames(2),plpars(20),plpars2d(2),clrs(nf1)
  integer :: pltrue,plmedian,plrange
  real :: sch
  character :: fnames(nf1)*99,dirnames(nf1)*99,outnamebase*99,settingsfile*99
end module comp_pdfs_settings
!***************************************************************************************************





!***************************************************************************************************
subroutine plotpdf1d(pp,lbl)
  use comp_pdfs_settings
  implicit none
  integer, parameter :: np=15,nbin1=500
  integer :: pp,b,p1,io,f,nplvar,nchains,nbin,p(np),pp1,ic,wrap(nf,np),lw,detnan,identical
  real :: x,startval(nf,np,2),stats(nf,np,6),ranges(nf,np,5),xmin1(nf,np),xmax1(nf,np),plshift
  real :: xbin1(nf,np,nbin1),ybin1(nf,np,nbin1),xbin(nbin1),ybin(nbin1),xmin,xmax,ymin,ymax,dx,yrange(2),xpeak
  character :: fname1*99,fname2*99,fname*99,varnss(np)*99,pgvarnss(np)*99,pgunits(np)*99,lbl*99,str*99
  
  
  !varnss(1:15)  = (/'log L','Mc','eta','t_c','d_L','a_spin','theta_SL','R.A.','Dec.','phi_c','theta_J0','phi_J0','alpha_c','M1','M2'/)
  varnss(1:15)  = (/'log L','Mc','eta','t_c','d_L','a_spin','theta_SL','R.A.','Dec.','phi_c','incl','polang','alpha_c','M1','M2'/)
  pgvarnss(1:15)  = (/'log L','M\dc\u','\(2133)','t\dc\u','d\dL\u','a\dspin\u','\(2134)\dSL\u','R.A.','Dec.','\(2147)\dc\u','acos(J\(2236)N)','\(2149)\dJ0\u','\(2127)\dc\u','M\d1\u','M\d2\u'/)
  !Include units
  pgvarnss(1:15)  = (/'log L','M\dc\u (M\d\(2281)\u)','\(2133)','t\dc\u (s)','d\dL\u (Mpc)','a\dspin\u','\(2134)\dSL\u (\(2218))','R.A. (h)','Dec. (\(2218))','\(2147)\dc\u (\(2218))',  &
       'acos(J\(2236)N) (\(2218))','\(2149)\dJ0\u (\(2218))','\(2127)\dc\u (\(2218))','M\d1\u (M\d\(2281)\u)','M\d2\u (M\d\(2281)\u)'/)
  !Units only
  pgunits(1:15)  = (/'','M\d\(2281)\u ','','s','Mpc','','\(2218)','\uh\d','\(2218)','\(2218)','\(2218)','\(2218)','\(2218)','M\d\(2281)\u','M\d\(2281)\u'/)
  
  
  detnan = 0 !Used to detect NaNs
  do f=1,nf
     !if(f.eq.1) fname = fname1
     !if(f.eq.2) fname = fname2
     fname = fnames(f)
     if(fname(1:3).eq.'   ') cycle
     !print*,f,trim(fname)
     open(unit=10,action='read',form='formatted',status='old',position='rewind',file=trim(fname),iostat=io)
     !print*,trim(fname),io
     if(io.ne.0) then
        write(*,'(A,I5)')'  Error opening file '//trim(fname)//', error:',io
        !cycle
        write(*,'(A,/)')'  Aborting.'
        stop
     end if
     read(10,'(3I6)')nplvar,nchains,nbin !'Total number of plot variables, total number of chains, number of bins'
     !print*,nplvar,nchains,nbin
     do p1=1,nplvar
        read(10,'(3I6)')ic,p(p1),wrap(f,p1) !'Chain number, variable number, and wrap'
        read(10,'(2E15.7)')startval(f,p1,1:2) !'True and starting value'
        read(10,'(6E15.7)')stats(f,p1,1:6) !'Stats: median, mean, absvar1, absvar2, stdev1, stdev2'
        read(10,'(5E15.7)')ranges(f,p1,1:5) !'Ranges: lower,upper limit, centre, width, relative width'
        read(10,'(2E15.7)')xmin1(f,p1),xmax1(f,p1) !'Xmin and Xmax of PDF'
        do b=1,nbin+1
           !read(10,'(2E15.7)',iostat=io)xbin1(f,p1,b),ybin1(f,p1,b)
           read(10,*,iostat=io)xbin1(f,p1,b),ybin1(f,p1,b)
           if(io.ne.0) write(*,'(A,I4,A1,I4,A1)')'  Error reading file '//trim(fname)//', variable '//trim(varnss(p1))//', bin', &
                b,'/',nbin,'.'
           if(p(p1).eq.pp .and. (.not.ybin1(f,p1,b).gt.-1.e30) .and. (.not.ybin1(f,p1,b).lt.1.e30)) then  !Then it's probably a NaN
              detnan = 1
              ybin1(f,p1,b) = 0.
           end if
        end do
        if(p(p1).eq.pp) pp1 = p1
     end do !p1
     
     close(10)
  end do !f
  
  !print*,pp1
  if(detnan.gt.0) write(*,'(A)')'  Warning:  I think I detected NaNs and set them to zero in the PDF for '//trim(varnss(pp))//'.'
  
  xmin = minval(xmin1(1:nf,pp1))
  xmax = maxval(xmax1(1:nf,pp1))
  dx = abs(xmax-xmin)
  xmin = xmin - 0.1*dx
  xmax = xmax + 0.1*dx
  ymin = 0.
  ymax = -1.e30
  do f=1,nf
     !ymax = max(ymax,maxval(ybin1(f,pp1,1:nbin+1)))
     do b=1,nbin+1
        if(ybin1(f,pp1,b).gt.ymax) then
           ymax = ybin1(f,pp1,b)
           xpeak = xbin1(f,pp1,b)
        end if
     end do
  end do
  !print*,f,pp1,xpeak
  
  ymax = ymax*1.1
  if(nf.eq.1) ymax = ymax*1.2 !Make room for numbers
  
  yrange = (/-1.e30,1.e30/) !Used to plot true value, probability range
  if(nf.eq.1) yrange = (/-1.e30,ymax*0.9/)
  
  !print*,xmin,xmax,ymin,ymax
  
  lw = 3 !Line width for pdf contours, not hatches
  call pgslw(lw)
  call pgsch(sch*0.5)
  call pgswin(xmin,xmax,ymin,ymax)
  
  !if(nf.gt.1) call pgsfs(3) !Hatches
  call pgsfs(fs)
  do f=1,nf
     !if(nf.gt.1.and.(clr.eq.0.or.clr.eq.2).and.fs.eq.2) call pgsls(mod(f,4))
     xbin = xbin1(f,pp1,:)
     ybin = ybin1(f,pp1,:)
     call pgshs(45.,0.7,0.) !Hatches slanted up
     if(f.eq.2) call pgshs(-45.,0.7,0.) !Hatches slanted down
     if(f.eq.3) call pgshs(0.,0.7,0.) !Hatches horizontal
     if(f.eq.4) call pgshs(90.,0.7,0.) !Hatches vertical
     
     call pgsci(clrs(f))
     if(nf.eq.1) call pgsci(15)
     
     if(wrap(f,pp1).eq.0) then
        !if(clr.eq.0) call pgsci(15)
        call pgslw(1)
        call pgpoly(nbin+2,(/xbin(1),xbin(1:nbin+1)/),(/0.,ybin(1:nbin+1)/))
        !Plot pdf contour
        if(nf.eq.1.and.clr.eq.1) call pgsci(2)
        if(clr.eq.0.or.clr.eq.2) call pgsci(1)
        call pgslw(lw)
        call pgline(nbin+1,xbin(1:nbin+1),ybin(1:nbin+1)) !:nbin) ?
        
        !Fix the loose ends
        call pgline(2,(/xbin(1),xbin(1)/),(/0.,ybin(1)/))
        call pgline(2,(/xbin(nbin+1),xbin(nbin+1)/),(/ybin(nbin+1),0./))
     else
        !plshift = real(2*pi)
        plshift = 360.
        !if(clr.eq.0) call pgsci(15)
        call pgslw(1)
        call pgpoly(nbin+3,(/xbin(1),xbin(1:nbin),xbin(1)+plshift,xbin(1)+plshift/),(/0.,ybin(1:nbin),ybin(1),0./))
        
        !Plot pdf contour
        if(nf.eq.1.and.clr.eq.1) call pgsci(2)
        if(clr.eq.0.or.clr.eq.2) call pgsci(1)
        !if(clr.eq.0) call pgsci(1)
        call pgslw(lw)
        call pgline(nbin,xbin(1:nbin),ybin(1:nbin))
        
        !Plot dotted lines outside the pdf for wrapped periodic variables
        call pgsls(4)
        call pgline(nbin+1,(/xbin(1:nbin)-plshift,xbin(1)/),(/ybin(1:nbin),ybin(1)/))
        call pgline(nbin,xbin+plshift,ybin)
        
        !Fix the loose end
        call pgsls(1)
        call pgslw(lw)
        call pgline(2,(/xbin(nbin),xbin(1)+plshift/),(/ybin(nbin),ybin(1)/))
     end if
     
     !plot probability ranges
     if(plrange.ge.1.and.nf.eq.1) then
        !Plot limits
        call pgsls(4)
        if(clr.eq.1) call pgsci(2)
        if(clr.eq.0.or.clr.eq.2) call pgsci(14)
        !if(clr.eq.0) call pgsci(1)
        call pgline(2,(/ranges(f,pp1,1),ranges(f,pp1,1)/),yrange)
        call pgline(2,(/ranges(f,pp1,2),ranges(f,pp1,2)/),yrange)
        
        !Print number
        if(plrange.ge.2) then
           x = ranges(f,pp1,5)
           if(pp1.eq.2.or.pp1.eq.3.or.pp1.eq.5.or.pp1.eq.6.or.pp1.eq.14.or.pp1.eq.15) x = x*100
           write(str,'(F10.3)')x
           if(pp1.eq.2.or.pp1.eq.3.or.pp1.eq.5.or.pp1.eq.6.or.pp1.eq.14.or.pp1.eq.15) write(str,'(A)')trim(str)//'%'
           
           if(x.lt.0.01) write(str,'(F6.4)')x
           if(x.ge.0.01.and.x.lt.0.1) write(str,'(F5.3)')x
           if(x.ge.0.1.and.x.lt.1) write(str,'(F4.2)')x
           if(x.ge.1.and.x.lt.10) write(str,'(F3.1)')x
           if(x.ge.10.and.x.lt.100) write(str,'(I2)')nint(x)
           if(x.ge.100) write(str,'(I3)')nint(x)
           write(str,'(A)')'\(2030): '//trim(str)
           if(pp1.eq.2.or.pp1.eq.3.or.pp1.eq.5.or.pp1.eq.6.or.pp1.eq.14.or.pp1.eq.15) then
              write(str,'(A)')trim(str)//'%'
           else
              write(str,'(A)')trim(str)//trim(pgunits(pp1))
           end if
           !call pgsch(sch*1.2)
           !call pgptxt(xmin+0.05*dx,ymax*0.9,0.,0.,trim(pgvarnss(pp1)))
           call pgptxt(ranges(f,pp1,3),yrange(2)*1.05,0.,0.5,trim(str))
           call pgsls(1)
           call pgline(2,(/ranges(f,pp1,1),ranges(f,pp1,2)/),(/yrange(2),yrange(2)/))
        end if
     end if
     call pgsls(1)
  end do !f
  
  !Plot lines again over surface of overlapping distributions
  if(1.eq.2) then
     call pgsls(4)
     !if(clr.eq.0) call pgsci(1)
     do f=1,nf
        xbin = xbin1(f,pp1,:)
        ybin = ybin1(f,pp1,:)
        if(wrap(f,pp1).eq.0) then
           call pgline(nbin+1,xbin(1:nbin+1),ybin(1:nbin+1))
        else
           call pgline(nbin,xbin(1:nbin),ybin(1:nbin))
        end if
     end do
     call pgsls(1)
  end if
  
  
  
  !Plot true values
  if(pltrue.eq.1) then
     identical = 1
     if(nf.gt.1) then
        do f=2,nf
           if(abs(startval(f,pp1,1)-startval(1,pp1,1)).gt.1.e-10) identical = 0
        end do
     end if
     
     if(identical.eq.1) then
        call pgsci(1); call pgslw(lw);call pgsls(2) !; if(clr.eq.0) call pgsci(1)
        !if(nf.eq.1) call pgsci(2)
        if(.not.(pp1.eq.6.and.abs(startval(1,5,1)).lt.0.001)) call pgline(2,(/startval(1,pp1,1),startval(1,pp1,1)/),yrange)   !If a=0, don't plot theta_SL (pp1=6)
     else
        do f=1,nf
           if(pp1.eq.6.and.abs(startval(f,5,1)).lt.0.001) cycle  !If a=0, don't plot theta_SL (pp1=6)
           call pgslw(lw); call pgsls(2); call pgsci(1)
           !call pgsci(clrs(f))  !it seems clearer when all true values are white
           call pgline(2,(/startval(f,pp1,1),startval(f,pp1,1)/),yrange)
        end do
     end if
  end if !if(pltrue.eq.1)
  
  
  
  
  !Plot median (if one file)
  if(plmedian.eq.1.and.nf.eq.1) then
     if(clr.eq.1) call pgsci(2); call pgslw(lw);call pgsls(2)  !; if(clr.eq.0) call pgsci(1)
     call pgline(2,(/stats(1,pp1,1),stats(1,pp1,1)/),yrange)
  end if !if(plmedian.eq.1.and.nf.eq.1)
  

  
  
  call pgsls(1)
  call pgsci(1)
  !if(nf.eq.1) then
  call pgbox('BNTS',0.0,0,'',0.0,0)
  !else
  !   call pgbox('BCNTS',0.0,0,'BC',0.0,0)
  !end if
  
  !call pgmtxt('B',2.4,0.5,0.5,trim(pgvarnss(pp1)))  !Plot label under x axis
  if(abs(xpeak-xmin).gt.abs(xpeak-xmax)) then  !Peak is right, plot varname left
     call pgmtxt('T',-1.,0.05,0.,trim(pgvarnss(pp1)))  !Plot label in upper-left corner
     call pgmtxt('T',-1.5,0.95,1.,trim(lbl))
  else   !Peak is left, plot varname right
     call pgmtxt('T',-1.,0.95,1.,trim(pgvarnss(pp1)))  !Plot label in upper-left corner
     call pgmtxt('T',-1.5,0.05,0.,trim(lbl))
  end if
  
  
end subroutine plotpdf1d
!***************************************************************************************************











!***************************************************************************************************
subroutine plotpdf2d(pp1,pp2,lbl)
  use comp_pdfs_settings
  implicit none
  integer, parameter :: np=15,nbinx1=500,nbiny1=500
  integer :: pp,pp1,pp2,bx,by,p1,p2,p11,p22,pp11,pp22,pp12,p12,io,f,nplvar,nplvar1,nplvar2,nchains,nbinx,nbiny,p(np),ic,lw,c,foundit
  integer :: identical
  real :: startval(nf,np,2),stats(nf,np,6),ranges(nf,np,5),xmin1(nf,np),xmax1(nf,np),ymin1(nf,np),ymax1(nf,np),x
  real :: xmin,xmax,ymin,ymax,dx,dy,z(nf,nbinx1,nbiny1),z1(nbinx1,nbiny1),tr(nf,np*np,6),cont(11)
  character :: fname1*99,fname2*99,fname*99,pgvarnss(np)*99,pgunits(np)*99,lbl*99,str*99
  
  
  pgvarnss(1:15)  = (/'log L','M\dc\u','\(2133)','t\dc\u','d\dL\u','a\dspin\u','\(2134)\dSL\u','R.A.','Dec.','\(2147)\dc\u','acos(J\(2236)N)','\(2149)\dJ0\u','\(2127)\dc\u','M\d1\u','M\d2\u'/)
  !Include units
  pgvarnss(1:15)  = (/'log L','M\dc\u (M\d\(2281)\u)','\(2133)','t\dc\u (s)','d\dL\u (Mpc)','a\dspin\u','\(2134)\dSL\u (\(2218))','R.A. (h)','Dec. (\(2218))','\(2147)\dc\u (\(2218))',  &
       'acos(J\(2236)N) (\(2218))','\(2149)\dJ0\u (\(2218))','\(2127)\dc\u (\(2218))','M\d1\u (M\d\(2281)\u)','M\d2\u (M\d\(2281)\u)'/)
  
  !Units only
  pgunits(1:15)  = (/'','M\d\(2281)\u ','','s','Mpc','','\(2218)','\uh\d','\(2218)','\(2218)','\(2218)','\(2218)','\(2218)','M\d\(2281)\u','M\d\(2281)\u'/)

  
  
  do f=1,nf
     foundit = 0
     !if(f.eq.1) fname = fname1
     !if(f.eq.2) fname = fname2
     fname = fnames(f)
     if(fname(1:3).eq.'   ') cycle
     !print*,f,trim(fname)
     open(unit=10,action='read',form='formatted',status='old',position='rewind',file=trim(fname),iostat=io)
     !print*,io
     if(io.ne.0) then
        write(*,'(A,I5)')'  Error reading file '//trim(fname)//', error:',io
        cycle
     end if
     read(10,'(5I6)')nplvar1,nplvar2,nchains,nbinx,nbiny !Plot variable 1,2, total number of chains, number of bins x,y
     nplvar = nplvar2-nplvar1+1
     !print*,nplvar1,nplvar2,nchains,nbinx,nbiny
     p12 = 0
     !print*,nplvar1,nplvar2-1,p11,nplvar2
     do p11=nplvar1,nplvar2-1
        do p22=p11+1,nplvar2
           if(foundit.eq.1) cycle
           p12 = p12+1
           !print*,p11,p22
           read(10,'(3I6)')ic,p1,p2 !'Chain number and variable number 1,2'
           !print*,p1,p2
           !read(10,'(2E15.7)')startval(f,p11,1:2) !'True and starting value p1'
           !read(10,'(2E15.7)')startval(f,p22,1:2) !True and starting value p2'
           !read(10,'(6E15.7)')stats(f,p11,1:6) !Stats: median, mean, absvar1, absvar2, stdev1, stdev2 for p1'
           !read(10,'(6E15.7)')stats(f,p22,1:6) !Stats: median, mean, absvar1, absvar2, stdev1, stdev2 for p2'
           !read(10,'(5E15.7)')ranges(f,p11,1:5) !Ranges: lower,upper limit, centre, width, relative width for p1'
           !read(10,'(5E15.7)')ranges(f,p22,1:5) !Ranges: lower,upper limit, centre, width, relative width for p2'
           !read(10,'(4E15.7)')xmin1(f,p11),xmax1(f,p11),ymin1(f,p22),ymax1(f,p22) !'Xmin,Xmax,Ymin,Ymax of PDF'
           read(10,'(2E15.7)')startval(f,p1,1:2) !'True and starting value p1'
           read(10,'(2E15.7)')startval(f,p2,1:2) !True and starting value p2'
           read(10,'(6E15.7)')stats(f,p1,1:6) !Stats: median, mean, absvar1, absvar2, stdev1, stdev2 for p1'
           read(10,'(6E15.7)')stats(f,p2,1:6) !Stats: median, mean, absvar1, absvar2, stdev1, stdev2 for p2'
           read(10,'(5E15.7)')ranges(f,p1,1:5) !Ranges: lower,upper limit, centre, width, relative width for p1' The limits are NOT used as plot limits!!!
           read(10,'(5E15.7)')ranges(f,p2,1:5) !Ranges: lower,upper limit, centre, width, relative width for p2'
           read(10,'(4E15.7)')xmin1(f,p1),xmax1(f,p1),ymin1(f,p2),ymax1(f,p2) !'Xmin,Xmax,Ymin,Ymax of PDF'
           read(10,'(6E15.7)')tr(f,p12,1:6) !'Tr'              
           !print*,tr
           do bx=1,nbinx+1
              do by=1,nbiny
                 read(10,'(E15.7)',advance='no')z1(bx,by)
              end do
              read(10,'(E15.7)')z1(bx,by)
           end do
           
           !print*,p1,p2
           !print*,pp1,pp2
           
           !if(p11.eq.pp1.and.p22.eq.pp2) then !This only works if you print the whole set of pdfs to file
           if(p1.eq.pp1.and.p2.eq.pp2) then
              pp12 = p12
              !pp11 = p11
              !pp22 = p22
              pp11 = p1
              pp22 = p2
              z(f,:,:) = z1(:,:)
              !print*,pp1,pp2,p11,p22
              !print*,f,xmin1(f,p11),xmax1(f,p11),ymin1(f,p22),ymax1(f,p22)
              !write(6,'(6F10.4)')tr(f,p12,1:6)
              foundit = 1 !Found it!
              !print*,'  Found it !',f
           end if
           !print*,foundit
           
           !end do !while foundit.eq.0
        end do !p11
     end do !p22
     
     close(10)
  end do !f
  
  !print*,p12
  
  
  xmin = minval(xmin1(1:nf,pp11))
  xmax = maxval(xmax1(1:nf,pp11))
  ymin = minval(ymin1(1:nf,pp22))
  ymax = maxval(ymax1(1:nf,pp22))
  
  !print*,xmin,xmin1(1:nf,pp11)
  !print*,xmax,xmax1(1:nf,pp11)
  !print*,ymin,ymin1(1:nf,pp22)
  !print*,ymax,ymax1(1:nf,pp22)
  !print*,''
  
  if(pp1.eq.8) then
     xmin = maxval(xmin1(1:nf,pp11))
     xmax = minval(xmax1(1:nf,pp11))
  end if
  
  if(1.eq.2.and.pp1.eq.8.and.pp2.eq.9) then !Get whole sky for sky plot
     xmin = 24.
     xmax = 0.
     ymin = -90.
     ymax = 90.
  end if
  
  dx = xmax - xmin
  dy = ymax - ymin
  
  lw = 3
  call pgsls(1)
  call pgsci(1)
  call pgslw(lw)
  call pgsch(sch)
  call pgswin(xmin,xmax,ymin,ymax)
  call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)
  call pgmtxt('B',2.4,0.5,0.5,trim(pgvarnss(pp1)))  !Cheat a bit
  call pgmtxt('L',2.2,0.5,0.5,trim(pgvarnss(pp2)))  !Cheat a bit
  
  !select contours
  do c=1,11
     cont(c) = 0.01 + 2*real(c-1)/10.
  end do
  
  do f=1,nf
     !Fill contours with hatches
     call pgsci(1)
     !call pgsfs(3) !Hatches
     call pgsfs(fs)
     call pgshs(45.,0.7,0.) !Hatches slanted up
     if(f.eq.2) call pgshs(-45.,0.7,0.) !Hatches slanted down
     if(f.eq.3) call pgshs(0.,0.7,0.) !Hatches horizontal
     if(f.eq.4) call pgshs(90.,0.7,0.) !Hatches vertical
     
     !Fill contours with solid colour, override hatches
     !call pgsci(15) !Light grey
     !call pgsfs(1)
     
     call pgslw(1)
     if(f.eq.2) call pgsci(14) !Dark grey
     !if(f.eq.2) call pgshs(-45.,0.7,0.) !Hatches down
     call pgsci(clrs(f))
     call pgconf(z(f,1:nbinx+1,1:nbiny+1),nbinx+1,nbiny+1,1,nbinx+1,1,nbiny,cont(1),cont(11),tr(f,pp12,1:6))
     !call pgsci(1) !When not using hatches
     call pgslw(lw)
     call pgcont(z(f,1:nbinx+1,1:nbiny+1),nbinx+1,nbiny+1,1,nbinx+1,1,nbiny,cont(1:1),1,tr(f,pp12,1:6))
     
     
     !plot probability ranges
     !if(plrange.eq.1) then
     if(plrange.ge.1.and.nf.eq.1) then
        !Plot limits
        call pgsls(1)
        if(clr.eq.1) call pgsci(2)
        call pgsci(clrs(f))
        x = 0.05+f*0.02
        call pgarro(ranges(f,pp1,3),ymin+dy*x,ranges(f,pp1,1),ymin+dy*x)
        call pgarro(ranges(f,pp1,3),ymin+dy*x,ranges(f,pp1,2),ymin+dy*x)
        call pgarro(xmin+dx*x*0.7,ranges(f,pp2,3),xmin+dx*x*0.7,ranges(f,pp2,1))
        call pgarro(xmin+dx*x*0.7,ranges(f,pp2,3),xmin+dx*x*0.7,ranges(f,pp2,2))
        
        !Print number
        !if(nf.eq.1) then
        if(plrange.ge.2.and.nf.eq.1) then
           x = ranges(f,pp1,5)
           if(pp1.eq.2.or.pp1.eq.3.or.pp1.eq.5.or.pp1.eq.6.or.pp1.eq.14.or.pp1.eq.15) x = x*100
           write(str,'(F10.3)')x
           if(pp1.eq.2.or.pp1.eq.3.or.pp1.eq.5.or.pp1.eq.6.or.pp1.eq.14.or.pp1.eq.15) write(str,'(A)')trim(str)//'%'
           
           if(x.lt.0.01) write(str,'(F6.4)')x
           if(x.ge.0.01.and.x.lt.0.1) write(str,'(F5.3)')x
           if(x.ge.0.1.and.x.lt.1) write(str,'(F4.2)')x
           if(x.ge.1.and.x.lt.10) write(str,'(F3.1)')x
           if(x.ge.10.and.x.lt.100) write(str,'(I2)')nint(x)
           if(x.ge.100) write(str,'(I3)')nint(x)
           write(str,'(A)')'\(2030): '//trim(str)
           if(pp1.eq.2.or.pp1.eq.3.or.pp1.eq.5.or.pp1.eq.6.or.pp1.eq.14.or.pp1.eq.15) then
              write(str,'(A)')trim(str)//'%'
           else
              write(str,'(A)')trim(str)//trim(pgunits(pp1))
           end if
           !call pgptxt(ranges(f,pp1,3),yrange(2)*1.05,0.,0.5,trim(str))
           call pgsls(1)
           !call pgline(2,(/ranges(f,pp1,1),ranges(f,pp1,2)/),(/yrange(2),yrange(2)/))
        end if
     end if !if(plrange.ge.1.and.nf.eq.1)
     call pgsls(1)
     
  end do !do f=1,nf
  
  
  !Plot true values 
  if(pltrue.eq.1) then
     !call pgline(2,(/startval(1,pp11,1),startval(1,pp11,1)/),(/-1.e20,1.e20/))
     !call pgline(2,(/-1.e20,1.e20/),(/startval(1,pp22,1),startval(1,pp22,1)/))
     
     identical = 1
     if(nf.gt.1) then
        do f=2,nf
           if(abs(startval(f,pp11,1)-startval(1,pp11,1)).gt.1.e-10) identical = 0
        end do
     end if
     !if(abs((startval(1,pp11,1) - sum(startval(1:nf,pp11,1))/(real(nf)))/startval(1,pp11,1)).lt.0.001) then !Then all have the same value
     if(identical.eq.1) then
        call pgslw(lw);call pgsls(2); call pgsci(1)
        call pgline(2,(/startval(1,pp11,1),startval(1,pp11,1)/),(/-1.e20,1.e20/))
     else
        do f=1,nf
           call pgslw(lw);call pgsls(2)
           call pgsci(1)
           !call pgsci(clrs(f))  !it seems clearer when all true values are white
           call pgline(2,(/startval(f,pp11,1),startval(f,pp11,1)/),(/-1.e20,1.e20/))
        end do
     end if
     
     identical = 1
     do f=2,nf
        if(abs(startval(f,pp22,1)-startval(1,pp22,1)).gt.1.e-10) identical = 0
     end do
     !if(abs((startval(1,pp22,1) - sum(startval(1:nf,pp22,1))/(real(nf)))/startval(1,pp22,1)).lt.0.001) then !Then all have the same value
     if(identical.eq.1) then
        call pgslw(lw);call pgsls(2); call pgsci(1)
        call pgline(2,(/-1.e20,1.e20/),(/startval(1,pp22,1),startval(1,pp22,1)/))
     else
        do f=1,nf
           call pgslw(lw);call pgsls(2)
           call pgsci(1)
           !call pgsci(clrs(f))  !it seems clearer when all true values are white
           call pgline(2,(/-1.e20,1.e20/),(/startval(f,pp22,1),startval(f,pp22,1)/))
        end do
     end if
  end if !if(pltrue.eq.1)
  
  
  !Plot median (if one file)
  if(plmedian.eq.1.and.nf.eq.1) then
     if(clr.eq.1) call pgsci(2); call pgslw(lw);call pgsls(2)  !; if(clr.eq.0) call pgsci(1)
     call pgline(2,(/stats(1,pp1,1),stats(1,pp1,1)/),(/-1.e20,1.e20/))
     call pgline(2,(/(/-1.e20,1.e20/),stats(1,pp2,1),stats(1,pp2,1)/))
  end if !if(plmedian.eq.1.and.nf.eq.1)
  
  
  call pgmtxt('T',-1.5,0.05,0.,trim(lbl))
  
end subroutine plotpdf2d
!***************************************************************************************************










!***************************************************************************************************
!subroutine plotwave(fname1,fname2)
subroutine plotwave(fname1,thingy,lbl)
  implicit none
  integer, parameter :: nf=1,n1=1e6
  integer :: i,f,n(nf),io,thingy,lw
  integer :: nfrx,nfry,frx,fry,fr,colours(nf)
  real :: xwinmin,xwinmax,ywinmin,ywinmax,dxwin,dywin,xfrmin,xfrmax,yfrmin,yfrmax,sch
  real :: t(nf,n1),h(nf,n1),dx,dy,xmin,xmax,ymin,ymax
  real*8 :: t1,t0,m1,m2,mc,eta,tc,dl,lat,lon,phase,spin,kappa,thJ0,phJ0,alpha
  character :: fname*99,fname1*99,fname2*99,bla,lbl*99
  
  do f=1,nf
     fname = fname1
     if(f.eq.2) fname = fname2
     open(unit=10,form='formatted',status='old',file=trim(fname),iostat=io)
     rewind(10)
     
     !write(6,'(A,$)')'Reading input file '//trim(fname)//'...     '
     read(10,*)bla
     read(10,*)m1,m2,mc,eta,tc,dl,lat,lon,phase,spin,kappa,thJ0,phJ0,alpha
     read(10,*)bla
     do i=1,n1
        read(10,*,err=195,end=199)t1,h(f,i)
        if(i.eq.1) t0 = (nint(t1*1.d-3)-1)*1.d3
        t(f,i) = real(t1-t0)
     end do
195  write(6,'(A,I8)')'error reading file '//trim(fname)//' line',i+1
199  close(10)
     n(f) = i-1
     !write(6,'(I,A12)')n(f),'lines read.'
  end do
  
  h = h*1.e21
  
  xmin =  1.e30
  xmax = -1.e30
  ymin =  1.e30
  ymax = -1.e30
  do f=1,nf
     xmin = min(minval(t(f,1:n(f))),xmin)
     xmax = max(maxval(t(f,1:n(f))),xmax)
     ymin = min(minval(h(f,1:n(f))),ymin)
     ymax = max(maxval(h(f,1:n(f))),ymax)
  end do
  
  ymax = maxval(abs(h))
  ymin = -ymax
  
  xmin = 1006.65
  xmin = 1010.55
  xmin = 1011.6
  xmax = 1012.32
  !xmin = 1012.3
  xmax = 1012.36
  
  
  
  dx = abs(xmax-xmin)
  dy = abs(ymax-ymin)
  
  lw = 3
  call pgslw(lw)
  nfrx = 1      !Number of frames
  nfry = nf      
  xwinmin = 0.10    !Absolute limits
  xwinmax = 0.95
  ywinmin = 0.20
  ywinmax = 0.95
  if(nf.eq.1) ywinmin = 0.25    
  if(nf.eq.1) ywinmax = 0.8
  if(nf.eq.2) ywinmax = 0.9
  if(nf.eq.3) ywinmax = 0.95
  !call pgswin(xmin,xmax,2*ymin-dy*0.05,2*ymax+dy*0.05)
  call pgswin(xmin,xmax,ymin-dy*0.05,ymax+dy*0.05)
  if(thingy.eq.1) call pgbox('BCTS',0.0,0,'BCNTS',0.0,0)
  if(thingy.eq.2) call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)
  !call pgbox('BCNTS',0.0,0,'BC',0.0,0)
  if(thingy.eq.2) call pgmtxt('B',2.4,0.5,0.5,'time (s)')
  call pgmtxt('L',2.2,0.5,0.5,'h (10\u-21\d)')
  do f=1,nf
     call pgslw(1)
     !call pgline(n(f),t(f,1:n(f)),h(f,1:n(f))-(f-1.5)*dy)
     call pgline(n(f),t(f,1:n(f)),h(f,1:n(f)))
     !call pgsch(sch)
     call pgsci(1)
     call pgslw(lw)
  end do  !f
  call pgsch(1.)
  call pgmtxt('T',-1.5,0.03,0.,trim(lbl))
  
  
end subroutine plotwave
!***************************************************************************************************



!***************************************************************************************************
subroutine read_inputfile
  use comp_pdfs_settings
  implicit none
  integer :: u,i,io
  character :: bla
  
  u = 14
  open(unit=u,form='formatted',status='old',action='read',file=trim(settingsfile),iostat=io)
  if(io.ne.0) write(*,'(A)')'  Error opening input file.'
  
  
    read(u,*)bla
    
    read(u,*)bla
    read(u,*)nf
    read(u,*)file
    read(u,*)type
    read(u,*)dim
    read(u,*)pltrue
    read(u,*)plmedian
    read(u,*)plrange
    read(u,*)clr
    read(u,*)fs
    read(u,*)sch

    read(u,*)bla
    read(u,*)bla
    read(u,*)frames
    read(u,*)bla
    read(u,*)plpars

    read(u,*)bla
    read(u,*)bla
    read(u,*)plpars2d
    
    read(u,*)bla
    read(u,*)bla
    read(u,*)fnames(1:nf)
    read(u,*)bla
    read(u,*)dirnames(1:nf)
    read(u,*)bla
    read(u,*)outnamebase
    
    close(u)
end subroutine read_inputfile
!***************************************************************************************************



!***************************************************************************************************
subroutine write_inputfile
  use comp_pdfs_settings
  implicit none
  integer :: u,i,io
  
  u = 14
  open(unit=u,form='formatted',status='replace',file='comp_pdfs.used',iostat=io)
  if(io.ne.0) write(*,'(A)')'  Error writing input file.'
  
11 format(I10,1x,A9,5x,A)
12 format(1x,2I5)
13 format(3x,20I3)

21 format(F10.5,1x,A9,5x,A)

  
  write(u,'(A,/)')' Input file for comp_pdfs.f'
  
  
  write(u,'(/,A)')' General plot options:'
  write(u,11)nf, 'nf',   'Number of input files'
  write(u,11)file, 'file',   'Output file type: 0:screen, 1: png, 2:eps, 3:pdf'
  write(u,11)type, 'type',   'Type of plot:  1: default, 2: poster, 3: talk'
  write(u,11)dim, 'dim',   'Dimension of PDF: 1 or 2'
  write(u,11)pltrue, 'pltrue',   'Plot true values for 1D or 2D plots'
  write(u,11)plmedian, 'plmedian',   'Plot medians for 1D or 2D plots, only if reading 1 file'
  write(u,11)plrange, 'plrange',   'Plot ranges  for 1D or 2D plots, only if reading 1 file: 0-no, 1:plot lines, 2:plot numbers too'
  write(u,11)clr, 'clr',   'Use colour: 0-no (grey), 1-yes'
  write(u,11)fs, 'fs',   'Fill style for PDFs: 1-solid, 2-outline, 3-hatched, 4-cross-hatched'
  write(u,21)sch, 'sch',   'Scale of font, etc.  Default: 2.'
  
  
  write(u,'(//,A)')' 1D plot options:'
  write(u,'(A)')'   Number of frames in 1D plot: (w, h  or  x, y)  (frames):'
  write(u,12)frames
  write(u,'(A)')'   Labels of the plot parameters for the 1d frames (plpars):'
  write(u,13)plpars !Array of size 20
  
  
  write(u,'(//,A)')' 2D plot options:'
  write(u,'(A)')'   Labels of the plot parameters for the 2d plot (plpars2d):'
  write(u,12)plpars2d
  
  
  write(u,'(//,A)')' Input files:'
  write(u,'(A)')'   Base of input files names (without extension) (fnames):'
  write(u,'(A,$)')'     '
  do i=1,nf
     write(u,'(A,$)')trim(fnames(i))//'   '
  end do
  write(u,*)
  
  write(u,'(/,A)')'   Directories of the input files (dirnames):'
  write(u,'(A,$)')'     '
  do i=1,nf
     write(u,'(A,$)')trim(dirnames(i))//'   '
  end do
  write(u,*)
  
  write(u,'(/,A)')'   Base of the output file name (without extention, can start with dir) (outnamebase)'
  write(u,'(A,$)')'     '//trim(outnamebase)
  
  
  write(u,*)
  close(u)
end subroutine write_inputfile
!***************************************************************************************************


