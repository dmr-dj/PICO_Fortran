!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
!>      VUA and IPSL/LSCE by the iLOVECLIM / LUDUS coding group / Within the LUDUS code environement
!
!       LICENSING TERMS:
!>      \copyright
!! ---|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----|
!!   Copyright [2021] [Didier M. Roche]
!!
!!   Licensed under the Apache License, Version 2.0 (the "License");
!!   you may not use this file except in compliance with the License.
!!   You may obtain a copy of the License at
!!
!!       http://www.apache.org/licenses/LICENSE-2.0
!!
!!   Unless required by applicable law or agreed to in writing, software
!!   distributed under the License is distributed on an "AS IS" BASIS,
!!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!!   See the License for the specific language governing permissions and
!!   limitations under the License.
!! ---|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----| 
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      MODULE: global_constants_mod
!
!>     @author  Didier M. Roche (dmr)
!
!
!>     @brief This module global_constants_mod is storing the global physical constants needed for the computation within
!>              the different modules in LUDUS
!
!>     @date Creation date: March, 1st, 2016
!>     @date Last modification: $LastChangedDate$
!>     @author Last modified by : dmr
!
!>     @version This is svn version: $LastChangedRevision$
!
!>     Here add the long_description of the module ...
!!     The constants in the different modules of the LOVECLIM and thus iLOVECLIM and LUDUS framework have been used without internal
!!       consistency between the different components. While this is probably fine in common applications, it may complicated the
!!       checkings of mass, energy etc. conservations between the different components. It also hampers the possibility of updating
!!       the value of the different constants easily.
!!     Thus the current module is developped with explicit relation to the physical constants, their units and the reference that
!!        have been used. The module will be slowly expanded to be included in as much components of the framework as possible.
!!     It will not be included in external svn resources to avoid conflict at external sub-component level.
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module global_constants_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   stdin, stdout definition ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|


#define f2003
#ifdef f2003
       use, intrinsic :: iso_fortran_env, only: stdin=>input_unit, stdout=>output_unit, stderr=>error_unit
#else
#define stdin  5
#define stdout 6
#define stderr 0
#endif

       implicit none
       public

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!   Selected double precision as taken from:
!   Metcalf, M., J. Reid, and M. Cohen (2004). Fortran 95/2003 Explained. Oxford University Press.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       integer, parameter  :: sp = 4      ! kind(1.0) does not work if compilation is forced to double as in e.g. -fdefault-real-8
       integer, parameter  :: dp = selected_real_kind(2*precision(1.0_sp))
       integer, parameter  :: qp = selected_real_kind(2*precision(1.0_dp))

       integer, parameter  :: sip = 2, ip = 4, dip = 8

       integer, parameter  :: str_len =256

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   Physical constants from http://arxiv.org/abs/1507.07956v1
!       CODATA group for physical constants
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!       speed_lightspeed of light in vacuum c, c0 299 792 458 m s−1 exact
       real(dp), parameter :: ligth_speed = 299792458.0_dp

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   Other Physical constants
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

! dmr  Temperature of 0°C in K                   (in K)
       real(dp), parameter :: tK_zero_C  = 273.15_dp

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Time constants (taken in part from https://www.cfa.harvard.edu/~dfabricant/huchra/ay145/constants.html)
!      Given in seconds (SI)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       real(dp), parameter :: second_s     = 1.0_dp
       real(dp), parameter :: minute_s     = 60.0_dp * second_s
       real(dp), parameter :: hour_s       = 60.0_dp * minute_s
       real(dp), parameter :: solar_day    = 24.0_dp * hour_s
       real(dp), parameter :: sidereal_day = 23.0_dp * hour_s + 56.0d0 * minute_s + 04.09054d0 * second_s

       real(dp), parameter :: months_year  = 12.0_dp    ! in days
       real(dp), parameter :: days_year360d= 360.0_dp   ! in days
       real(dp), parameter :: days_year365d= 365.0_dp   ! in days

       integer(ip), parameter :: months_year_i = nint(months_year)
       integer(ip), parameter :: days_year360d_i = nint(days_year360d)
       integer(ip), parameter :: days_month360d_i = nint(days_year360d/months_year)
       integer(ip), parameter :: days_year365d_i = nint(days_year365d)
       integer(ip), parameter :: seconds_month30d_i = solar_day * nint(days_year360d/months_year)
       integer(ip), parameter :: seconds_solar_day_i = nint(solar_day)

       real(dp), parameter :: sidereal_yr  = 3.155815E7_dp                 ! in seconds
       real(dp), parameter :: sidereal_yr_d= 3.155815E7_dp / sidereal_day  ! in sidereal days

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Physical constants for Earth (planet)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

! dmr  Radius of the Earth ...              (in meters)
       real(dp), parameter :: radius_earth = 6.37E6_dp

! dmr  Solar Constant (as seen from the Earth s upper atmosphere) ...    (in W.m^-2)
       real(dp), parameter :: solar_constant = 1365.0_dp ! should be 1368 from ECBilt value / afq, back to pre-rev 638 value

!dmr   Gravitational acceleration for the surface of the earth (World Geodetic System 1984)

       real(dp), parameter :: gravit_constant = 9.80665_dp  ! [m/s^2]

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Physical constants for water (H2O)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

! dmr  Volumetric mass density of water ... (in kg.m-3)
       real(dp), parameter :: vol_mass_dens_wat= 1000.0_dp
! dmr  Latent heat of vaporization ...      (in J.kg-1) ! ams value = 2.501+06
       real(dp), parameter :: latheat_vap_wat  = 2.5E6_dp
! dmr  Latent heat of sublimation ...       (in J.kg-1) ! ams value = 2.834+06
       real(dp), parameter :: latheat_sub_wat  =2.8E6_dp
! dmr  Latent heat of fusion ...            (in J.kg-1) ! ams value = 0.337+06
       real(dp), parameter :: latheat_fus_wat  =0.3E6_dp

       
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Mathematical constants
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

! dmr  PI, the number ...                              (no unit)
       real(dp), parameter :: pi_dp      = 4.0_dp*datan(1.0_dp)
! dmr  radian
       real(dp), parameter :: rad_to_deg = 180.0_dp/pi_dp
! dmr  degrees
       real(dp), parameter :: deg_to_rad = pi_dp/180.0_dp

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   Other constants
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

! dmr  Altitude of Olympus Mons, on Planet Mars. Can be used as a big constant ...
       real(dp), parameter :: alt_olympus_mons = 21229.0_dp
! afq  Tiny numbers
       real(sp), parameter :: eps_sp = epsilon(eps_sp)
       real(dp), parameter :: eps_dp = epsilon(eps_dp)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the module here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       end module global_constants_mod

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
