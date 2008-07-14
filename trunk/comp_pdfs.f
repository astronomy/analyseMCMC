!Compare 1d pdfs, use pdf output produced by plotspins


program comp_pdfs
  use comp_pdfs_settings
  integer :: nfrx,nfry,frx,fry,fr,fr1,f,i,system
  real :: size,rat,xwinmin,xwinmax,ywinmin,ywinmax,dxwin,dywin,xfrmin,xfrmax,yfrmin,yfrmax,dxfr,dyfr,space
  character :: fname1*99,fname2*99,lbl*99,varnss(15)*99,outname*99,exts(0:3)*4
  
  write(*,*)
  if(iargc().eq.1) then
     call getarg(1,settingsfile)
  else
     write(*,'(A,/)')'  Syntax:  comp_pdfs <input_file>'
     stop
  end if
  
  file = 1 !0:screen, 1: png, 2:eps, 3:pdf
  type = 1 !1: default, 2: poster, 3: talk
  dim = 1  !Dimension: 1 or 2 for 1d or 2d pdfs
  clr = 1  !Use colour: 0-no (grey), 1-yes
  fs  = 3  !Fill style for PDFs: 1-solid, 2-outline, 3-hatched, 4-cross-hatched
  
  sch = 2. !Scale of font etc.
  
  call read_inputfile()  !Read plot settings from file
  call write_inputfile() !Write plot settings back to file
  
  clrs(1:9) = (/2,4,6,3,5,8,9,15,7/)
  !if(type.eq.2) clrs(1:9) = (/2,4,5,3,6,8,9,15,7/) !Avoid red and magenta in posters
  if(type.eq.2) clrs(1:9) = (/2,4,3,5,6,8,9,15,7/) !Screw the colour blind and use red and green?
  !clrs(1:9) = (/2,4,6,3,2,4,6,3,5/) !Repeat the first 4
  !clrs(1:9) = (/2,4,2,4,2,4,2,4,5/) !Repeat the first 2
  !clrs(1:9) = (/2,2,4,4,5,5,5,5,5/) !Repeat 2x
  
  if(clr.eq.0) clrs = 1 !All white or black
  if(clr.eq.2) clrs(1:9) = (/15,14,16,15,14,1,15,14,1/) !Grey values
  
  size = 10.8
  rat  = 0.57
  if(file.eq.1) then
     size = 9.99
     rat  = 0.75
  end if
  if(file.ge.2) then
     size = 10.6
     rat  = 0.75 !3x4
     if(dim.eq.1.and.type.eq.2) rat  = 0.925 !4x3 for poster
     !if(dim.eq.1.and.type.eq.3) rat  = 0.82 !for talk (beamer)
     if(type.eq.3) rat  = 0.82 !for talk (beamer), 1&2D
  end if
  
  
  
  
  
  varnss(1:15)  = (/'logL','Mc','eta','t_c','d_L','a_spin','theta_SL','RA','Dec','phi_c','theta_J0','phi_J0','alpha_c','M1','M2'/)
  do f=1,nf
     if(dim.eq.2) then
        write(fnames(f),'(A)')trim(dirnames(f))//'/'//trim(fnames(f))//'__pdf2d.dat'
     else
        write(fnames(f),'(A)')trim(dirnames(f))//'/'//trim(fnames(f))//'__pdf1d.dat'
     end if
     print*,f,trim(fnames(f))
  end do
  
  if(file.eq.0) then
     call pgbegin(1,'/xs',1,1)
     call pgpap(size,rat)
     call pgsch(1.3)
  end if
  if(file.eq.1) then
     call pgbegin(1,'plot.ppm/ppm',1,1)
     call pgpap(size,rat)
     call pgsch(1.)
     !call pgscf(2)
  end if
  if(file.ge.2) then
     call pgbegin(1,'plot.eps/cps',1,1)
     call pgpap(size,rat)
     call pgscf(2)
  end if
  
  call pgscr(15,0.8,0.8,0.8)  !Default: 0.7
  call pgscr(14,0.45,0.45,0.45) !Default: 0.33
  
  nfrx = frames(1)      !Number of frames
  nfry = frames(2)
  if(dim.eq.2) then
     nfrx = 1
     nfry = 1
  end if
  space = 0.3 !Fraction of space between the panels (>=0.0)
  xwinmin = 0.0    !Absolute limits
  xwinmax = 0.99
  ywinmin = 0.0
  ywinmax = 0.99
  dxwin = (xwinmax-xwinmin)/real(nfrx)
  dywin = (ywinmax-ywinmin)/real(nfry)
  dxfr = dxwin/real(nfrx)*real(nfrx-1) * space
  dyfr = dxfr/rat
  do fry = 1,nfry
     do frx = 1,nfrx
        fr = (fry-1)*nfrx + frx
        xfrmin = xwinmin + dxwin*(frx-1) + dxfr
        xfrmax = xwinmin + dxwin*frx
        yfrmin = ywinmax - dywin*fry     + dyfr
        yfrmax = ywinmax - dywin*(fry-1)
        !if(fr.eq.1) xfrmax = xfrmax + dxwin !Combine frame 1 and 2
        !if(fr.eq.2) cycle
        if(nfrx*nfry.eq.1) then
           xfrmin = 0.1
           xfrmax = 0.9
           yfrmin = 0.1
           yfrmax = 0.9
           if(sch.gt.1.5) then
              xfrmin = 0.15
              yfrmin = 0.15
           end if
        end if
        call pgsvp(xfrmin,xfrmax,yfrmin,yfrmax)
        
        !fname1 = '1d1_0.50_055_pdf1d.dat'
        !fname2 = '2d12_0.50_055_pdf1d.dat'
        lbl = ' '
        !fr1 = fr+1 !The parameter that gets plotted: 2-13
        fr1 = plpars(fr) !The parameter that gets plotted
        
        !if(fr1.lt.4) fr1 = fr1 + 12  !Use m1,m2 iso Mc,eta
        !fr1 = fr1+2  !Skip the masses
        !if(fr1.gt.5) fr1 = fr1+2 !Skip the spin parameters
        !print*,fr,fr1

        
        if(dim.eq.1) then
           !call plotpdf1d(nf,fnames,fr1,lbl,clr,sch)
           call plotpdf1d(fr1,lbl)
        else
           !call plotpdf2d(nf,fnames,plpars2d(1),plpars2d(2),lbl,clr,sch)
           call plotpdf2d(plpars2d(1),plpars2d(2),lbl)
        end if
        
        !fname1 = '1d1_0.50_055_pdf2d.dat'
        !fname2 = '2d12_0.50_055_pdf2d.dat'
        !lbl = 'd)'
        !if(fr.eq.4) call plotpdf2d(fname1,fname2,14,15,lbl)    !M1,M2
        
     end do  !frx
  end do  !fry
  call pgend
  
  exts = (/'    ','.png','.eps','.pdf'/) !Extensions for the different file types
  if(dim.eq.1) write(outname,'(A,I1,A)')trim(outnamebase)//'__',dim,'d'//exts(file)
  if(dim.eq.2) write(outname,'(A,I1,A)')trim(outnamebase)//'__',dim,'d__'//trim(varnss(plpars2d(1)))//'-'//trim(varnss(plpars2d(2)))//exts(file)
  if(file.eq.1) then
     !if(dim.eq.1) i = system('convert -depth 8 plot.ppm comp_pdfs1d.png')
     !if(dim.eq.2) i = system('convert -depth 8 plot.ppm comp_pdfs2d__'//trim(varnss(plpars2d(1)))//'-'//trim(varnss(plpars2d(2)))//'.png')
     i = system('convert -depth 8 plot.ppm '//trim(outname))
     i = system('rm -f plot.ppm')
  end if
  if(file.eq.2) then
     !i = system('mv -f plot.eps comp_pdfs.eps')
     i = system('mv -f plot.eps '//trim(outname))
  end if
  if(file.eq.3) then
     i = system('eps2pdf plot.eps >& /dev/null')
     i = system('mv -f plot.pdf '//trim(outname))
     i = system('rm -f plot.eps')
  end if
  if(file.ge.1) write(*,'(/,A,/)')'  Plot saved as '//trim(outname)
  
end program comp_pdfs
!***************************************************************************************************

