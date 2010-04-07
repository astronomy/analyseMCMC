!> AnalyseMCMC_tailored.f:
!! Produce tailored ASCII output, e.g. for (LaTeX) tables, etc.
!! 
!<





!***********************************************************************************************************************************
subroutine tailored_output(exitcode)
   use constants
   !use analysemcmc_settings
   use mcmcrun_data
   
   implicit none
   integer :: exitcode
   integer :: out
   character :: outname*99
   
   exitcode = 1
   out = 40  !Output unit
   outname = trim(outputdir)//'/'//trim(outputname)//'__tailoredOutput.txt'
   open(unit=out,action='write',status='replace',form='formatted',file=trim(outname))
   
   
   
   !Output format for ApJL 2008
   if(tailoredOutput.eq.1) then
      
      !call tailored_output_0002(out)
      !close(out)
      !exitcode = 0
   end if
   
   
   
   !Output format for methods paper 2010:
   if(tailoredOutput.eq.2) then
      exitcode = 0
      call tailored_output_0002(out,exitcode)
      close(out)
      if(exitcode.ne.0) return
   end if
   
   
   
   close(out)
   
   if(prProgress.ge.1) then
      if(exitcode.eq.0) then
         write(stdOut,'(A,/)')'  Tailored output was saved in '//trim(outName)//'.'
      else
         write(stdOut,'(A,I4,A,/)')'  tailoredOutput =',tailoredOutput,' has not been implemented yet...'
      end if
   end if
   
end subroutine tailored_output
!***********************************************************************************************************************************




!***********************************************************************************************************************************
subroutine tailored_output_0002(out,exitcode)
   !< Output format for methods paper 2010
   use constants
   use analysemcmc_settings
   use general_data
   use mcmcrun_data
   use stats_data
   use chain_data
   
   implicit none
   integer :: out,exitcode
   integer :: par,par1,par2,ic
   real*8 :: x
   character :: runID*99,col*99,row*99
   logical :: sky_position,binary_orientation
   
   ic = 1
   col = '  &'
   row = '  \\'
   
   
   !Do we have all the necessary settings?
   if(waveform.ne.3 .or. nMCMCpar0.ne.15) then
      write(stdOut,'(A,/)')'  I need output from the 15-parameter SpinTaylor waveform template for tailoredOutput=2'
      exitcode = 1
      return
   end if
   
   if(prIval.eq.0.or.plPDF2D.eq.0.or.normPDF2D.ne.4) then
      write(stdOut,'(A,/)')'  You need to set prIval>0, plPDF2D>0 and normPDF2D=4 for tailoredOutput=2'
      exitcode = 1
      return
   end if
   
   sky_position = .false.
   binary_orientation = .false.
   do par=1,Npdf2D
      if(PDF2Dpairs(par,1).eq.31.and.PDF2Dpairs(par,2).eq.32) sky_position = .true.
      if(PDF2Dpairs(par,1).eq.52.and.PDF2Dpairs(par,2).eq.51) binary_orientation = .true.
   end do
   if(.not.(sky_position.and.binary_orientation)) then
      write(stdOut,'(A,/)')'  You need to enable 2D PDFs for both sky position and binary orientation for tailoredOutput=2'
      exitcode = 1
      return
   end if
   
   
   
   
   
   !Assemble runID:
   !injection:
   runID = ''
   if(startval(ic,revID(71),1).lt.0.01d0) then  !Non-spinning
      write(runID,'(A)')'N-'
   else 
      if(startval(ic,revID(71),1).lt.0.5d0) then
         write(runID,'(A)')'L'           !Low spins
      else
         write(runID,'(A)')'H'           !High spins
      end if
      if(startval(ic,revID(72),1).lt.45.d0) then
         write(runID,'(A)')trim(runID)//'S-'  !Small angles
      else
         write(runID,'(A)')trim(runID)//'L-'  !Large angles
      end if
   end if
   
   
   !recovery template:
   write(runID,'(A,I1)')trim(runID),spinningRun  !spinningRun = 0,1,2 for the number of spins allowed for in the recovery template
   write(out,'(A6)',advance="no")trim(runID)
   write(out,'(A)',advance="no")trim(col)
   
   
   !Injection parameters:
   par = revID(22)  !Distance
   x = startval(ic,par,1)
   write(out,'(F6.1)',advance="no") x
   write(out,'(A)',advance="no")trim(col)
   !write(out,'(A)',advance="no")trim(col)
   
   
   
   
   !Accuracies:
   write(out,'(A)',advance="no")'     '
   
   !Chirp mass:
   par = revID(61)
   x = ranges(ic,ival0,par,4)/ranges(ic,ival0,par,3)*100 
   if(x.lt.9.95) then
      write(out,'(F6.1)',advance="no") x
   else
      write(out,'(I6)',advance="no") nint(x)
   end if
   write(out,'(A)',advance="no")trim(col)
   
   !Eta:
   par = revID(62)
   x = ranges(ic,ival0,par,4)
   write(out,'(F7.3)',advance="no") x
   write(out,'(A)',advance="no")trim(col)
   
   !M1:
   par = revID(63)
   x = ranges(ic,ival0,par,4)/ranges(ic,ival0,par,3)*100 
   if(x.lt.9.95) then
      write(out,'(F6.1)',advance="no") x
   else
      write(out,'(I6)',advance="no") nint(x)
   end if
   write(out,'(A)',advance="no")trim(col)
   
   !M2:
   par = revID(64)
   x = ranges(ic,ival0,par,4)/ranges(ic,ival0,par,3)*100 
   if(x.lt.9.95) then
      write(out,'(F6.1)',advance="no") x
   else
      write(out,'(I6)',advance="no") nint(x)
   end if
   write(out,'(A)',advance="no")trim(col)
   
   
   
   
   !Distance:
   par = revID(22)
   x = ranges(ic,ival0,par,4)/ranges(ic,ival0,par,3)*100
   write(out,'(I5)',advance="no")nint(x)
   write(out,'(A)',advance="no")trim(col)
   
   !t_c:
   par = revID(11)
   x = ranges(ic,ival0,par,4)*1000 !s->ms
   if(x.lt.9.95) then
      write(out,'(F6.1)',advance="no") x
   else
      write(out,'(I6)',advance="no") nint(x)
   end if
   write(out,'(A)',advance="no")trim(col)
   
   
   
   !Spins:
   if(spinningRun.ge.1) then
      
      !a_spin1
      par = revID(71)
      x = ranges(ic,ival0,par,4)
      if(x.lt.0.0995) then
         write(out,'(F7.3)',advance="no") x
      else
         write(out,'(F7.2)',advance="no") x
      end if
      write(out,'(A)',advance="no")trim(col)
      
      !theta_spin1
      par = revID(72)
      x = ranges(ic,ival0,par,4)
      if(x.lt.9.95) then
         write(out,'(F6.1)',advance="no") x
      else
         write(out,'(I6)',advance="no") nint(x)
      end if
      write(out,'(A)',advance="no")trim(col)
      
      !phi_spin1
      par = revID(73)
      x = ranges(ic,ival0,par,4)
      if(x.lt.9.95) then
         write(out,'(F6.1)',advance="no") x
      else
         write(out,'(I6)',advance="no") nint(x)
      end if
      write(out,'(A)',advance="no")trim(col)
      
      
      
      if(spinningRun.ge.2) then
         
         !a_spin2
         par = revID(81)
         x = ranges(ic,ival0,par,4)
         if(x.lt.0.0995) then
            write(out,'(F7.3)',advance="no") x
         else
            write(out,'(F7.2)',advance="no") x
         end if
         write(out,'(A)',advance="no")trim(col)
         
         !theta_spin2
         par = revID(82)
         x = ranges(ic,ival0,par,4)
         if(x.lt.9.95) then
            write(out,'(F6.1)',advance="no") x
         else
            write(out,'(I6)',advance="no") nint(x)
         end if
         write(out,'(A)',advance="no")trim(col)
         
         !phi_spin2
         par = revID(83)
         x = ranges(ic,ival0,par,4)
         if(x.lt.9.95) then
            write(out,'(F6.1)',advance="no") x
         else
            write(out,'(I6)',advance="no") nint(x)
         end if
         write(out,'(A)',advance="no")trim(col)
         
      else
         
         write(out,'(A7,A,2(A6,A))',advance="no")' -- ',trim(col),' -- ',trim(col),' -- ',trim(col)
         
      end if ! 2 spins
      
   else
      write(out,'(A7,A,2(A6,A))',advance="no")' -- ',trim(col),' -- ',trim(col),' -- ',trim(col)
      write(out,'(A7,A,2(A6,A))',advance="no")' -- ',trim(col),' -- ',trim(col),' -- ',trim(col)
   end if ! Spins
   
   
   
   !RA:
   if(1.eq.2) then
      par = revID(31)
      x = ranges(ic,ival0,par,4)
      if(x.lt.9.95) then
         write(out,'(F6.1)',advance="no") x
      else
         write(out,'(I6)',advance="no") nint(x)
      end if
      write(out,'(A)',advance="no")trim(col)
   end if
   
   !Dec:
   if(1.eq.2) then
      par = revID(32)
      x = ranges(ic,ival0,par,4)
      if(x.lt.9.95) then
         write(out,'(F6.1)',advance="no") x
      else
         write(out,'(I6)',advance="no") nint(x)
      end if
      write(out,'(A)',advance="no")trim(col)
   end if
   
   !Sky position:
   par1 = revID(31)  !RA
   par2 = revID(32)  !Dec
   x = probAreas(par1,par2,ival0,3)  !3: area in sqare degrees
   if(x.lt.99.95d0) then
      write(out,'(F6.1)',advance="no") x
   else
      write(out,'(I6)',advance="no") nint(x)
   end if
   write(out,'(A)',advance="no")trim(col)
   
   
   !psi:
   par = revID(52)
   x = ranges(ic,ival0,par,4)
   if(x.lt.9.95) then
      write(out,'(F6.1)',advance="no") x
   else
      write(out,'(I6)',advance="no") nint(x)
   end if
   write(out,'(A)',advance="no")trim(col)
   
   !i:
   par = revID(51)
   x = ranges(ic,ival0,par,4)
   if(x.lt.9.95) then
      write(out,'(F6.1)',advance="no") x
   else
      write(out,'(I6)',advance="no") nint(x)
   end if
   write(out,'(A)',advance="no")trim(col)
   
   
   !Binary orientation:
   if(1.eq.2) then
      par1 = revID(52)  !Polarisation angle
      par2 = revID(51)  !Inclination
      x = probAreas(par1,par2,ival0,3)  !3: area in sqare degrees
      if(x.lt.99.95d0) then
         write(out,'(F6.1)',advance="no") x
      else
         write(out,'(I6)',advance="no") nint(x)
      end if
      write(out,'(A)',advance="no")trim(col)
   end if
   
   
   !Orbital phase:
   par = revID(41)
   x = ranges(ic,ival0,par,4)
   if(x.lt.9.95) then
      write(out,'(F6.1)',advance="no") x
   else
      write(out,'(I6)',advance="no") nint(x)
   end if
   !write(out,'(A)',advance="no")trim(col)
   
   
   
   write(out,'(A)')trim(row)
   
end subroutine tailored_output_0002
!***********************************************************************************************************************************
      
      
