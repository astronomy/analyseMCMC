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
      if(waveform.ne.3 .or. nMCMCpar0.ne.15) then
         write(stdOut,'(A,/)')'  I need output from the 15-parameter SpinTaylor waveform template for tailoredOutput=2'
         return
      end if
      
      call tailored_output_0002(out)
      close(out)
      exitcode = 0
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
subroutine tailored_output_0002(out)
   !< Output format for methods paper 2010
   use constants
   use analysemcmc_settings
   use general_data
   use mcmcrun_data
   use stats_data
   use chain_data
   
   implicit none
   integer :: out
   integer :: par,ic
   real*8 :: x
   character :: runID*99,col*99,row*99
   
   ic = 1
   col = '  &'
   row = '  \\'
   
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
      if(startval(ic,revID(72),1).gt.0.7d0) then
         write(runID,'(A)')trim(runID)//'S-'  !Small angles
      else
         write(runID,'(A)')trim(runID)//'L-'  !Large angles
      end if
   end if
   
   !recovery template:
   write(runID,'(A,I1)')trim(runID),spinningRun  !spinningRun = 0,1,2 for the number of spins allowed for in the recovery template
   write(out,'(A6,$)')trim(runID)
   write(out,'(A,$)')trim(col)
   
   !Injection parameters:
   par = revID(22)  !Distance
   x = startval(ic,par,1)
   write(out,'(F6.1,$)') x
   write(out,'(A,$)')trim(col)
   write(out,'(A,$)')trim(col)
   
   
   !Accuracies:
   write(out,'(A,$)')'     '
   
   !Chirp mass:
   par = revID(61)
   x = ranges(ic,ival0,par,4)/ranges(ic,ival0,par,3)*100 
   write(out,'(F6.2,$)') x
   write(out,'(A,$)')trim(col)
   
   !Eta:
   par = revID(62)
   x = ranges(ic,ival0,par,4)
   write(out,'(F7.3,$)') x
   write(out,'(A,$)')trim(col)
   
   !Distance:
   par = revID(22)
   x = ranges(ic,ival0,par,4)/ranges(ic,ival0,par,3)*100
   write(out,'(I5,$)')nint(x)
   write(out,'(A,$)')trim(col)
   
   !t_c (ms):
   par = revID(11)
   x = ranges(ic,ival0,par,4)*1000
   if(x.lt.9.95) then
      write(out,'(F6.1,$)') x
   else
      write(out,'(I6,$)') nint(x)
   end if
   write(out,'(A,$)')trim(col)
   
   
   write(out,'(A)')trim(row)
   
end subroutine tailored_output_0002
!***********************************************************************************************************************************
      
      
