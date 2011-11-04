!> \file analyseMCMC_textroutines.f90  Text functions and routines for analyseMCMC

! 
! LICENCE:
! 
! Copyright 2007-2011 Marc van der Sluys
!  
! This file is part of the AnalyseMCMC package.
!  
! This is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published
! by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! 
! This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License along with this code (LICENSE).  If not, see 
! <http://www.gnu.org/licenses/>.
! 



!***********************************************************************************************************************************
!> \brief Define the names and symbols of the original MCMC parameters

subroutine set_originalParameterNames()
  use analysemcmc_settings, only: fonttype
  use general_data, only: parNames,pgParNs,pgParNss,pgUnits, pgOrigParns
  implicit none
  
  parNames = ''
  pgParNs = ''
  pgParNss = ''
  pgUnits = ''
  
  ! Short ASCII names for text output:
  parNames(11:19)   = (/'tc      ','t40     ','        ','        ','        ','        ','        ','        ','        '/)
  parNames(21:29)   = (/'dl^3    ','log_dl  ','dl      ','        ','        ','        ','        ','        ','        '/)
  parNames(31:39)   = (/'RA      ','sin_dec ','dec     ','        ','        ','        ','        ','        ','        '/)
  parNames(41:49)   = (/'phase   ','        ','        ','        ','        ','        ','        ','        ','        '/)
  parNames(51:59)   = (/'cos_i   ','psi     ','sin_thJo','ph_Jo   ','i       ','        ','        ','        ','        '/)
  parNames(61:69)   = (/'Mc      ','eta     ','M1      ','M2      ','Mc_16   ','Mtot    ','q       ','        ','        '/)
  parNames(71:79)   = (/'spin1   ','cos_th1 ','phi1    ','theta1  ','        ','        ','        ','        ','        '/)
  parNames(81:89)   = (/'spin2   ','cos_th2 ','phi2    ','theta2  ','        ','        ','        ','        ','        '/)
  parNames(185:190) = (/'x1      ','x2      ','x3      ','x4      ','x5      ','x6      '/)
  parNames(191:199) = (/'x7      ','x8      ','x9      ','x10     ','x11     ','x12     ','x13     ','x14     ','x15     '/)
  !parNames(1:9) = (/'','','','','','','','',''/)
  
  
  if(fonttype.eq.2) then  ! Use 'roman-like' Greek font in PGPlot
     
     ! Long PGPlot names (symbol + unit):
     pgParNs(11) = 't\dc\u (s)'
     pgParNs(12) = 't\d40\u (s)'
     
     pgParNs(21) = 'd\dL\u\u3\d (Mpc)'
     pgParNs(22) = 'logd\dL\u (Mpc)'
     pgParNs(23) = 'd\dL\u (Mpc)'
     
     pgParNs(31) = '\(2127) (rad)'
     pgParNs(32) = 'sin \(2130)'
     pgParNs(33) = '\(2130) (rad)'
     
     pgParNs(41) = '\(2147)\dc\u (rad)'
     
     pgParNs(51) = 'cos \(2135)'
     pgParNs(52) = '\(2149) (rad)'
     pgParNs(53) = 'sin \(2185)\dJ0\u'
     pgParNs(54) = '\(2147)\dJ0\u (rad)'
     pgParNs(55) = '\(2135) (rad)'
     
     pgParNs(61) = '\(2563) (M\d\(2281)\u)'
     pgParNs(62) = '\(2133)'
     pgParNs(63) = 'M\d1\u (M\d\(2281)\u)'
     pgParNs(64) = 'M\d2\u (M\d\(2281)\u)'
     pgParNs(65) = '\(2563)\u1/6\d (M\d\(2281)\u\u1/6\d)'
     pgParNs(66) = 'M\dtot\u (M\d\(2281)\u)'
     pgParNs(67) = 'q'
     pgParNs(68) = 'log(q)'

     pgParNs(71) = 'a\dspin1\u'
     pgParNs(72) = 'cos \(2185)\dspin1\u'
     pgParNs(73) = '\(2147)\dspin1\u (rad)'
     pgParNs(74) = '\(2185)\dspin1\u (rad)'
     
     pgParNs(81) = 'a\dspin2\u'
     pgParNs(82) = 'cos \(2185)\dspin2\u'
     pgParNs(83) = '\(2147)\dspin2\u (rad)'
     pgParNs(84) = '\(2185)\dspin2\u (rad)'
     
     pgParNs(185) = 'x\d1\u'
     pgParNs(186) = 'x\d2\u'
     pgParNs(187) = 'x\d3\u'
     pgParNs(188) = 'x\d4\u'
     pgParNs(189) = 'x\d5\u'
     pgParNs(190) = 'x\d6\u'
     pgParNs(191) = 'x\d7\u'
     pgParNs(192) = 'x\d8\u'
     pgParNs(193) = 'x\d9\u'
     pgParNs(194) = 'x\d10\u'
     pgParNs(195) = 'x\d11\u'
     pgParNs(196) = 'x\d12\u'
     pgParNs(197) = 'x\d13\u'
     pgParNs(198) = 'x\d14\u'
     pgParNs(199) = 'x\d15\u'
     
     
     
     ! Short PGPlot symbols (no unit):
     pgParNss(11) = 't\dc\u'
     pgParNss(12) = 't\d40\u'
     
     pgParNss(21) = 'd\dL\u\u3\d'
     pgParNss(22) = 'logd\dL\u'
     pgParNss(23) = 'd\dL\u'
     
     pgParNss(31) = '\(2127)'
     pgParNss(32) = 'sin \(2130)'
     pgParNss(33) = '\(2130)'
     
     pgParNss(41) = '\(2147)\dc\u'
     
     pgParNss(51) = 'cos \(2135)'
     pgParNss(52) = '\(2149)'
     pgParNss(53) = 'sin \(2185)\dJ0\u'
     pgParNss(54) = '\(2147)\dJ0\u'
     pgParNss(55) = '\(2135)'
     
     pgParNss(61) = '\(2563)'
     pgParNss(62) = '\(2133)'
     pgParNss(63) = 'M\d1\u'
     pgParNss(64) = 'M\d2\u'
     pgParNss(65) = '\(2563)\u1/6\d'
     pgParNss(66) = 'M\dt\u'
     pgParNss(67) = 'q'
     pgParNss(68) = 'log(q)'

     pgParNss(71) = 'a\dspin1\u'
     pgParNss(72) = 'cos \(2185)\dspin1\u'
     pgParNss(73) = '\(2147)\dspin1\u'
     pgParNss(74) = '\(2185)\dspin1\u'
     
     pgParNss(81) = 'a\dspin2\u'
     pgParNss(82) = 'cos \(2185)\dspin2\u'
     pgParNss(83) = '\(2147)\dspin2\u'
     pgParNss(84) = '\(2185)\dspin2\u'
     
     pgParNss(185) = 'x\d1\u'
     pgParNss(186) = 'x\d2\u'
     pgParNss(187) = 'x\d3\u'
     pgParNss(188) = 'x\d4\u'
     pgParNss(189) = 'x\d5\u'
     pgParNss(190) = 'x\d6\u'
     pgParNss(191) = 'x\d7\u'
     pgParNss(192) = 'x\d8\u'
     pgParNss(193) = 'x\d9\u'
     pgParNss(194) = 'x\d10\u'
     pgParNss(195) = 'x\d11\u'
     pgParNss(196) = 'x\d12\u'
     pgParNss(197) = 'x\d13\u'
     pgParNss(198) = 'x\d14\u'
     pgParNss(199) = 'x\d15\u'
     
     
  else  ! Same, but replace '\(21' with \(06' for arial-like Greek font
     
     
     ! Long PGPlot names (symbol + unit):
     pgParNs(11) = 't\dc\u (s)'
     pgParNs(12) = 't\d40\u (s)'
     
     pgParNs(21) = 'd\dL\u\u3\d (Mpc)'
     pgParNs(22) = 'logd\dL\u (Mpc)'
     pgParNs(23) = 'd\dL\u (Mpc)'
     
     pgParNs(31) = '\(0627) (rad)'
     pgParNs(32) = 'sin \(0630)'
     pgParNs(33) = '\(0630) (rad)'
     
     pgParNs(41) = '\(0647)\dc\u (rad)'
     
     pgParNs(51) = 'cos \(0635)'
     pgParNs(52) = '\(0649) (rad)'
     pgParNs(53) = 'sin \(0685)\dJ0\u'
     pgParNs(54) = '\(0647)\dJ0\u (rad)'
     pgParNs(55) = '\(0635) (rad)'
     
     pgParNs(61) = '\(2563) (M\d\(2281)\u)'
     pgParNs(62) = '\(0633)'
     pgParNs(63) = 'M\d1\u (M\d\(2281)\u)'
     pgParNs(64) = 'M\d2\u (M\d\(2281)\u)'
     pgParNs(65) = '\(2563)\u1/6\d (M\d\(2281)\u\u1/6\d)'
     pgParNs(66) = 'M\dtot\u (M\d\(2281)\u)'
     pgParNs(67) = 'q'
     pgParNs(68) = 'log(q)'

     pgParNs(71) = 'a\dspin1\u'
     pgParNs(72) = 'cos \(0685)\dspin1\u'
     pgParNs(73) = '\(0647)\dspin1\u (rad)'
     pgParNs(74) = '\(0685)\dspin1\u (rad)'
     
     pgParNs(81) = 'a\dspin2\u'
     pgParNs(82) = 'cos \(0685)\dspin2\u'
     pgParNs(83) = '\(0647)\dspin2\u (rad)'
     pgParNs(84) = '\(0685)\dspin2\u (rad)'
     
     pgParNs(185) = 'x\d1\u'
     pgParNs(186) = 'x\d2\u'
     pgParNs(187) = 'x\d3\u'
     pgParNs(188) = 'x\d4\u'
     pgParNs(189) = 'x\d5\u'
     pgParNs(190) = 'x\d6\u'
     pgParNs(191) = 'x\d7\u'
     pgParNs(192) = 'x\d8\u'
     pgParNs(193) = 'x\d9\u'
     pgParNs(194) = 'x\d10\u'
     pgParNs(195) = 'x\d11\u'
     pgParNs(196) = 'x\d12\u'
     pgParNs(197) = 'x\d13\u'
     pgParNs(198) = 'x\d14\u'
     pgParNs(199) = 'x\d15\u'
     
     
     
     ! Short PGPlot symbols (no unit):
     pgParNss(11) = 't\dc\u'
     pgParNss(12) = 't\d40\u'
     
     pgParNss(21) = 'd\dL\u\u3\d'
     pgParNss(22) = 'logd\dL\u'
     pgParNss(23) = 'd\dL\u'
     
     pgParNss(31) = '\(0627)'
     pgParNss(32) = 'sin \(0630)'
     pgParNss(33) = '\(0630)'
     
     pgParNss(41) = '\(0647)\dc\u'
     
     pgParNss(51) = 'cos \(0635)'
     pgParNss(52) = '\(0649)'
     pgParNss(53) = 'sin \(0685)\dJ0\u'
     pgParNss(54) = '\(0647)\dJ0\u'
     pgParNss(55) = '\(0635)'
     
     pgParNss(61) = '\(2563)'
     pgParNss(62) = '\(0633)'
     pgParNss(63) = 'M\d1\u'
     pgParNss(64) = 'M\d2\u'
     pgParNss(61) = '\(2563)\u1/6\d'
     pgParNss(66) = 'M\dt\u'
     pgParNss(67) = 'q'
     pgParNss(68) = 'log(q)'

     pgParNss(71) = 'a\dspin1\u'
     pgParNss(72) = 'cos \(0685)\dspin1\u'
     pgParNss(73) = '\(0647)\dspin1\u'
     pgParNss(74) = '\(0685)\dspin1\u'
     
     pgParNss(81) = 'a\dspin2\u'
     pgParNss(82) = 'cos \(0685)\dspin2\u'
     pgParNss(83) = '\(0647)\dspin2\u'
     pgParNss(84) = '\(0685)\dspin2\u'
     
     pgParNss(185) = 'x\d1\u'
     pgParNss(186) = 'x\d2\u'
     pgParNss(187) = 'x\d3\u'
     pgParNss(188) = 'x\d4\u'
     pgParNss(189) = 'x\d5\u'
     pgParNss(190) = 'x\d6\u'
     pgParNss(191) = 'x\d7\u'
     pgParNss(192) = 'x\d8\u'
     pgParNss(193) = 'x\d9\u'
     pgParNss(194) = 'x\d10\u'
     pgParNss(195) = 'x\d11\u'
     pgParNss(196) = 'x\d12\u'
     pgParNss(197) = 'x\d13\u'
     pgParNss(198) = 'x\d14\u'
     pgParNss(199) = 'x\d15\u'
     
  end if
  
  
  ! PGPlot units (no names)
  pgUnits(11) = 's'
  pgUnits(12) = 's'
  
  pgUnits(21) = 'Mpc'
  pgUnits(22) = 'Mpc'
  pgUnits(23) = 'Mpc'
  
  pgUnits(31) = 'rad'
  pgUnits(33) = 'rad'
  
  pgUnits(41) = 'rad'
  
  pgUnits(52) = 'rad'
  pgUnits(54) = 'rad'
  pgUnits(55) = 'rad'
  
  pgUnits(61) = 'M\d\(2281)\u'
  pgUnits(63) = 'M\d\(2281)\u'
  pgUnits(64) = 'M\d\(2281)\u'
  pgUnits(65) = 'M\d\(2281)\u\u1/6\d'
  pgUnits(66) = 'M\d\(2281)\u'
  
  pgUnits(73) = 'rad'
  pgUnits(74) = 'rad'
  
  pgUnits(83) = 'rad'
  pgUnits(84) = 'rad'
  
  ! Save the original parameter names for use after they get changed
  pgOrigParns = pgParNs
  
  
  
  
end subroutine set_originalParameterNames
!***********************************************************************************************************************************






!***********************************************************************************************************************************
!> \brief Define the names and symbols of the derived MCMC parameters
!!
!! - e.g. d_L rather than d_L^3 or log(d_L), i rather than cos(i), etc.

subroutine set_derivedParameterNames()
  use analysemcmc_settings, only: fonttype
  use general_data, only: parNames,pgParNs,pgParNss,pgUnits
  
  implicit none
  
  
  parNames = ''
  pgParNs = ''
  pgParNss = ''
  pgUnits = ''
  
  ! Short ASCII names for text output:
  parNames(11:19)   = (/'tc      ','t40     ','        ','        ','        ','        ','        ','        ','        '/)
  parNames(21:29)   = (/'dl      ','dl      ','dl      ','        ','        ','        ','        ','        ','        '/)
  parNames(31:39)   = (/'RA      ','dec     ','dec     ','        ','        ','        ','        ','        ','        '/)
  parNames(41:49)   = (/'phase   ','        ','        ','        ','        ','        ','        ','        ','        '/)
  parNames(51:59)   = (/'incl    ','psi     ','th_Jo   ','ph_Jo   ','incl    ','        ','        ','        ','        '/)
  parNames(61:69)   = (/'Mc      ','eta     ','M1      ','M2      ','Mc      ','Mtot    ','q       ','        ','        '/)
  parNames(71:79)   = (/'spin1   ','th1     ','phi1    ','theta1  ','        ','        ','        ','        ','        '/)
  parNames(81:89)   = (/'spin2   ','th2     ','phi2    ','theta2  ','        ','        ','        ','        ','        '/)
  parNames(185:190) = (/'x1      ','x2      ','x3      ','x4      ','x5      ','x6      '/)
  parNames(191:199) = (/'x7      ','x8      ','x9      ','x10     ','x11     ','x12     ','x13     ','x14     ','x15     '/)
  !parNames(1:9) = (/'','','','','','','','',''/)
  
  
  if(fonttype.eq.2) then  ! Use 'roman-like' Greek font in PGPlot
     
     ! Long PGPlot names (symbol + unit):
     pgParNs(11) = 't\dc\u (s)'
     pgParNs(12) = 't\d40\u (s)'
     
     pgParNs(21) = 'd\dL\u (Mpc)'
     pgParNs(22) = 'd\dL\u (Mpc)'
     pgParNs(23) = 'd\dL\u (Mpc)'
     
     pgParNs(31) = '\(2127) (h)'
     pgParNs(32) = '\(2130) (\(2218))'
     pgParNs(33) = '\(2130) (\(2218))'
     
     pgParNs(41) = '\(2147)\dc\u (\(2218))'
     
     pgParNs(51) = '\(2135) (\(2218))'
     pgParNs(52) = '\(2149) (\(2218))'
     pgParNs(53) = '\(2185)\dJ0\u (\(2218))'
     pgParNs(54) = '\(2147)\dJ0\u (\(2218))'
     pgParNs(55) = '\(2135) (\(2218))'
     
     pgParNs(61) = '\(2563) (M\d\(2281)\u)'
     pgParNs(62) = '\(2133)'
     pgParNs(63) = 'M\d1\u (M\d\(2281)\u)'
     pgParNs(64) = 'M\d2\u (M\d\(2281)\u)'
     pgParNs(65) = '\(2563) (M\d\(2281)\u)'
     pgParNs(66) = 'M\dtot\u (M\d\(2281)\u)'
     pgParNs(67) = 'q'
     pgParNs(68) = 'log(q)'

     pgParNs(71) = 'a\dspin1\u'
     pgParNs(72) = '\(2185)\dspin1\u (\(2218))'
     pgParNs(73) = '\(2147)\dspin1\u (\(2218))'
     pgParNs(74) = '\(2185)\dspin1\u (\(2218))'
     
     pgParNs(81) = 'a\dspin2\u'
     pgParNs(82) = '\(2185)\dspin2\u (\(2218))'
     pgParNs(83) = '\(2147)\dspin2\u (\(2218))'
     pgParNs(84) = '\(2185)\dspin2\u (\(2218))'
     
     pgParNs(185) = 'x\d1\u'
     pgParNs(186) = 'x\d2\u'
     pgParNs(187) = 'x\d3\u'
     pgParNs(188) = 'x\d4\u'
     pgParNs(189) = 'x\d5\u'
     pgParNs(190) = 'x\d6\u'
     pgParNs(191) = 'x\d7\u'
     pgParNs(192) = 'x\d8\u'
     pgParNs(193) = 'x\d9\u'
     pgParNs(194) = 'x\d10\u'
     pgParNs(195) = 'x\d11\u'
     pgParNs(196) = 'x\d12\u'
     pgParNs(197) = 'x\d13\u'
     pgParNs(198) = 'x\d14\u'
     pgParNs(199) = 'x\d15\u'
     
     
     
     ! Short PGPlot symbols (no unit):
     pgParNss(11) = 't\dc\u'
     pgParNss(12) = 't\d40\u'
     
     pgParNss(21) = 'd\dL\u\u3\d'
     pgParNss(22) = 'logd\dL\u'
     pgParNss(21) = 'd\dL\u'
     
     pgParNss(31) = '\(2127)'
     pgParNss(32) = '\(2130)'
     pgParNss(33) = '\(2130)'
     
     pgParNss(41) = '\(2147)\dc\u'
     
     pgParNss(51) = '\(2135)'
     pgParNss(52) = '\(2149)'
     pgParNss(53) = '\(2185)\dJ0\u'
     pgParNss(54) = '\(2147)\dJ0\u'
     pgParNss(51) = '\(2135)'
     
     pgParNss(61) = '\(2563)'
     pgParNss(62) = '\(2133)'
     pgParNss(63) = 'M\d1\u'
     pgParNss(64) = 'M\d2\u'
     pgParNss(65) = '\(2563)'
     pgParNss(66) = 'M\dt\u'
     pgParNss(67) = 'q'
     pgParNss(68) = 'log(q)'

     pgParNss(71) = 'a\dspin1\u'
     pgParNss(72) = '\(2185)\dspin1\u'
     pgParNss(73) = '\(2147)\dspin1\u'
     pgParNss(74) = '\(2185)\dspin1\u'
     
     pgParNss(81) = 'a\dspin2\u'
     pgParNss(82) = '\(2185)\dspin2\u'
     pgParNss(83) = '\(2147)\dspin2\u'
     pgParNss(84) = '\(2185)\dspin2\u'
     
     pgParNss(185) = 'x\d1\u'
     pgParNss(186) = 'x\d2\u'
     pgParNss(187) = 'x\d3\u'
     pgParNss(188) = 'x\d4\u'
     pgParNss(189) = 'x\d5\u'
     pgParNss(190) = 'x\d6\u'
     pgParNss(191) = 'x\d7\u'
     pgParNss(192) = 'x\d8\u'
     pgParNss(193) = 'x\d9\u'
     pgParNss(194) = 'x\d10\u'
     pgParNss(195) = 'x\d11\u'
     pgParNss(196) = 'x\d12\u'
     pgParNss(197) = 'x\d13\u'
     pgParNss(198) = 'x\d14\u'
     pgParNss(199) = 'x\d15\u'
     
     
     
  else  ! Same, but replace '\(21' with \(06' for arial-like Greek font
     
     
     ! Long PGPlot names (symbol + unit):
     pgParNs(11) = 't\dc\u (s)'
     pgParNs(12) = 't\d40\u (s)'
     
     pgParNs(21) = 'd\dL\u (Mpc)'
     pgParNs(22) = 'd\dL\u (Mpc)'
     pgParNs(23) = 'd\dL\u (Mpc)'
     
     pgParNs(31) = '\(0627) (h)'
     pgParNs(32) = '\(0630) (\(2218))'
     pgParNs(33) = '\(0630) (\(2218))'
     
     pgParNs(41) = '\(0647)\dc\u (\(2218))'
     
     pgParNs(51) = '\(0635) (\(2218))'
     pgParNs(52) = '\(0649) (\(2218))'
     pgParNs(53) = '\(0685)\dJ0\u (\(2218))'
     pgParNs(54) = '\(0647)\dJ0\u (\(2218))'
     pgParNs(55) = '\(0635) (\(2218))'
     
     pgParNs(61) = '\(2563) (M\d\(2281)\u)'
     pgParNs(62) = '\(0633)'
     pgParNs(63) = 'M\d1\u (M\d\(2281)\u)'
     pgParNs(64) = 'M\d2\u (M\d\(2281)\u)'
     pgParNs(65) = '\(2563) (M\d\(2281)\u)'
     pgParNs(66) = 'M\dtot\u (M\d\(2281)\u)'
     pgParNs(67) = 'q'
     pgParNs(68) = 'log(q)'

     pgParNs(71) = 'a\dspin1\u'
     pgParNs(72) = '\(0685)\dspin1\u (\(2218))'
     pgParNs(73) = '\(0647)\dspin1\u (\(2218))'
     pgParNs(74) = '\(0685)\dspin1\u (\(2218))'
     
     pgParNs(81) = 'a\dspin2\u'
     pgParNs(82) = '\(0685)\dspin2\u (\(2218))'
     pgParNs(83) = '\(0647)\dspin2\u (\(2218))'
     pgParNs(84) = '\(0685)\dspin2\u (\(2218))'
     
     pgParNs(185) = 'x\d1\u'
     pgParNs(186) = 'x\d2\u'
     pgParNs(187) = 'x\d3\u'
     pgParNs(188) = 'x\d4\u'
     pgParNs(189) = 'x\d5\u'
     pgParNs(190) = 'x\d6\u'
     pgParNs(191) = 'x\d7\u'
     pgParNs(192) = 'x\d8\u'
     pgParNs(193) = 'x\d9\u'
     pgParNs(194) = 'x\d10\u'
     pgParNs(195) = 'x\d11\u'
     pgParNs(196) = 'x\d12\u'
     pgParNs(197) = 'x\d13\u'
     pgParNs(198) = 'x\d14\u'
     pgParNs(199) = 'x\d15\u'
     
     
     
     ! Short PGPlot symbols (no unit):
     pgParNss(11) = 't\dc\u'
     pgParNss(12) = 't\d40\u'
     
     pgParNss(21) = 'd\dL\u\u3\d'
     pgParNss(22) = 'logd\dL\u'
     pgParNss(23) = 'd\dL\u'
     
     pgParNss(31) = '\(0627)'
     pgParNss(32) = '\(0630)'
     pgParNss(33) = '\(0630)'
     
     pgParNss(41) = '\(0647)\dc\u'
     
     pgParNss(51) = '\(0635)'
     pgParNss(52) = '\(0649)'
     pgParNss(53) = '\(0685)\dJ0\u'
     pgParNss(54) = '\(0647)\dJ0\u'
     pgParNss(55) = '\(0635)'
     
     pgParNss(61) = '\(2563)'
     pgParNss(62) = '\(0633)'
     pgParNss(63) = 'M\d1\u'
     pgParNss(64) = 'M\d2\u'
     pgParNss(65) = '\(2563)'
     pgParNss(66) = 'M\dt\u'
     pgParNss(67) = 'q'
     pgParNss(68) = 'log(q)'

     pgParNss(71) = 'a\dspin1\u'
     pgParNss(72) = '\(0685)\dspin1\u'
     pgParNss(73) = '\(0647)\dspin1\u'
     pgParNss(74) = '\(0685)\dspin1\u'
     
     pgParNss(81) = 'a\dspin2\u'
     pgParNss(82) = '\(0685)\dspin2\u'
     pgParNss(83) = '\(0647)\dspin2\u'
     pgParNss(84) = '\(0685)\dspin2\u'
     
     pgParNss(185) = 'x\d1\u'
     pgParNss(186) = 'x\d2\u'
     pgParNss(187) = 'x\d3\u'
     pgParNss(188) = 'x\d4\u'
     pgParNss(189) = 'x\d5\u'
     pgParNss(190) = 'x\d6\u'
     pgParNss(191) = 'x\d7\u'
     pgParNss(192) = 'x\d8\u'
     pgParNss(193) = 'x\d9\u'
     pgParNss(194) = 'x\d10\u'
     pgParNss(195) = 'x\d11\u'
     pgParNss(196) = 'x\d12\u'
     pgParNss(197) = 'x\d13\u'
     pgParNss(198) = 'x\d14\u'
     pgParNss(199) = 'x\d15\u'
     
     
  end if
  
  ! PGPlot units (no names)
  pgUnits(11:19) = 's'
  pgUnits(12) = 's'
  
  pgUnits(21) = 'Mpc'
  pgUnits(22) = 'Mpc'
  pgUnits(23) = 'Mpc'
  
  pgUnits(31) = '\uh\d'
  pgUnits(32) = '\(2218)'
  pgUnits(33) = '\(2218)'
  
  pgUnits(41) = '\(2218)'
  
  pgUnits(51) = '\(2218)'
  pgUnits(52) = '\(2218)'
  pgUnits(53) = '\(2218)'
  pgUnits(54) = '\(2218)'
  pgUnits(55) = '\(2218)'
  
  pgUnits(61) = 'M\d\(2281)\u'
  pgUnits(63) = 'M\d\(2281)\u'
  pgUnits(64) = 'M\d\(2281)\u'
  pgUnits(65) = 'M\d\(2281)\u'
  pgUnits(66) = 'M\d\(2281)\u'
  
  pgUnits(72) = '\(2218)'
  pgUnits(73) = '\(2218)'
  pgUnits(74) = '\(2218)'
  
  pgUnits(82) = '\(2218)'
  pgUnits(83) = '\(2218)'
  pgUnits(84) = '\(2218)'
  
  
  
end subroutine set_derivedParameterNames
!***********************************************************************************************************************************



