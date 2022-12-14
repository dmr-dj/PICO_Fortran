!> file pico_functions.f90
!>      @author: Maxence Menthon, Didier M. Roche (dmr)
!>      @date creation: May 2021
!>      @date last modification: 31st May 2021
!>      @brief: This script gather all functions for the PICO model of sub-ice-shelves melting

!       LICENSING TERMS:
!>      \copyright
!! ---|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----|
!!   Copyright [2021] [Didier M. Roche, Pepijn Bakker, Maxence Menthon]
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

! Inputs : T0, S0, Ak, pk
        
!=========================================================================================================== 
! Excecution order of the functions with their inputs and output
!1.	Input : So 					Function : T_star_fun		Output : T_s
!2.	Input : Ak 					Function : g1_fun		Output : g1_pico
!3.	Input : g1_pico					Function : g2_fun 		Output : g2_pico
!4.	Input : T_s, S0, g1_pico			Function : x_b1_fun		Output : x_b1_pico       /!\  Only for box 1
!       Input : T_s, gi_pico, g2_pico, q_pico, Sk_m1 	Function : x_bk_fun		Output : x_bk_pico	 /!\  For all other boxes, after calculation of q_pico
!5.	Input : Tk_m1 and (x_b1_pico OR x_bk_pico)	Function : Tk_fun		Output : Tk_pico 	 /!\  In case of box 1 : Tk_m1 = T0
!6.	Input : x_b1_pico OR x_bk_pico 			Function : y_fun		Output : y_pico		
!7.	Input : y_pico, Sk_m1				Function : Sk_fun 		Output : Sk_pico
!8.	Input : Sk_pico, Tk_pico, pk		 	Function : melting_pico_fun	Output : mk_pico
!9.	Input : T0, S0, Tk_pico, Sk_pico		Function : overturn_fun 	Output : q_pico		 /!\  Only once, when we have Tk_pico and Sk_pico for the box1 	
!=========================================================================================================== 

module pico_functions

use global_constants_mod, only : dp, g=>gravit_constant, sip

implicit none

public :: pico_output, pico_go_array

private 

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----|
! shared constants

! dmr [NOTA]: latheat_fus_wat from global constant could be used, but value is 3., AMS value is 3.337 ...
real(kind=dp), parameter         :: L            = 3.35e5_dp      ! [J/kg]      Latent heat of fusion for ice

real(kind=dp), parameter         :: cp           = 3.974e3_dp     ! [J/kg/degC] specific heat capacity of ocean
real(kind=dp), parameter         :: spy          = 86400.*365._dp ! [s/yr]      seconds per yer

! PICO model
real(kind=dp), parameter         :: a_pico       = -0.0572_dp   ! [degC /psu]
real(kind=dp), parameter         :: b_pico       =  0.0788_dp   ! [degC]
real(kind=dp), parameter         :: c_pico       =  7.77e-8_dp  ! [degC/Pa]
real(kind=dp), parameter         :: alpha_pico   =  7.5e-5_dp   ! [1/degC]
real(kind=dp), parameter         :: beta_pico    =  7.7e-4_dp   ! [1/psu]
real(kind=dp), parameter         :: rho0         =  1033.0_dp   ! [kg/m^3]    ! = rho_star
real(kind=dp), parameter         :: rhoi         =   910.0_dp   ! [kg/m^3]
real(kind=dp), parameter         :: rhow         =  1028.0_dp   ! [kg/m^3]
real(kind=dp), parameter         :: gammaS       =  2.e-6_dp    ! [m/s]       salinity mixing coefficient
real(kind=dp), parameter         :: gammaT       =  5e-5_dp     ! [m/s]       temperature mixing coefficient $$
real(kind=dp), parameter         :: gammae       =  2e-5_dp     ! [m/s]       effective mixing coefficient $\gamma^\star_T$
real(kind=dp), parameter         :: C_b1_pico    =  1e6_dp      ! [m^6/s/kg]
real(kind=dp), parameter         :: T_rho        =  0_dp        ! [degC]
real(kind=dp), parameter         :: S_rho        =  34_dp       ! [psu]

real(kind=dp), parameter         :: nulambda     = (rhoi/rhow)*(L/cp) ! [degC]

!real(kind=dp), parameter        :: gT_Ml        =  2.02e-5_dp
!real(kind=dp), parameter        :: gT_Mq        =  99.32e-5_dp
!real(kind=dp), parameter        :: gT_Mp        =  132.9e-5_dp

type :: pico_output
    real(dp) :: melt
    real(dp) :: q
    real(dp) :: T_k_out
    real(dp) :: S_k_out
end type pico_output

contains

function pico_go_array(S_0, T_0, p_k, A_k) result (pico_out)

        real(kind=dp), dimension(:),   INTENT(IN) :: S_0             ! DIM == [nregions,]
        real(kind=dp), dimension(:),   INTENT(IN) :: T_0             ! DIM == [nregions,]
        real(kind=dp), dimension(:,:), INTENT(IN) :: p_k             ! DIM == [nregions,n_k]
        real(kind=dp), dimension(:,:), INTENT(IN) :: A_k             ! DIM == [nregions,n_k]
        !real(kind=dp), dimension(:,:,:,:), allocatable  :: dim_outputs     ! To have good shape of the outputs
        
        type(pico_output), dimension(:,:), allocatable :: pico_out    ! DIM == [nregions,n_k]


                                                                     ! DIM == [nregions]
        real(kind=dp), dimension(:), allocatable  :: T_ess, g_1, g_2, x_box, y_box, T_k,  S_k, T_km1, S_km1, q_out
        
        
        integer(kind=sip), parameter :: grounding_line=1                ! assumes that grounding line is at index = 1, can be changed to n_k if necessary
        integer(kind=sip)            :: indx_k
        
        ! dmr --- Allocate working arrays
        
        allocate(pico_out(LBOUND(p_k,dim=1):UBOUND(p_k,dim=1),& ! melt
                          LBOUND(p_k,dim=2):UBOUND(p_k,dim=2)))!,& ! q 
                         ! LBOUND(dim_outputs,dim=3):UBOUND(dim_outputs,dim=3),& ! T_k_out
                         ! LBOUND(dim_outputs,dim=4):UBOUND(dim_outputs,dim=4))) ! S_k_out
        allocate(T_ess(LBOUND(S_0,dim=1):UBOUND(S_0,dim=1)))
        allocate(g_1(LBOUND(S_0,dim=1):UBOUND(S_0,dim=1)))
        allocate(g_2(LBOUND(S_0,dim=1):UBOUND(S_0,dim=1)))
        allocate(x_box(LBOUND(S_0,dim=1):UBOUND(S_0,dim=1)))
        allocate(y_box(LBOUND(S_0,dim=1):UBOUND(S_0,dim=1)))
        allocate(T_k(LBOUND(S_0,dim=1):UBOUND(S_0,dim=1)))
        allocate(S_k(LBOUND(S_0,dim=1):UBOUND(S_0,dim=1)))
        allocate(T_km1(LBOUND(S_0,dim=1):UBOUND(S_0,dim=1)))
        allocate(S_km1(LBOUND(S_0,dim=1):UBOUND(S_0,dim=1)))
        allocate(q_out(LBOUND(S_0,dim=1):UBOUND(S_0,dim=1)))                             
        !dmr - from here I assume that the input variables are n_k, nregions
        
        ! initial work, on the first element == grounding_line
        
        !print *, '--- BOX NB : 1'

        T_km1 = T_0
        S_km1 = S_0
        
        T_ess = T_star_fun(S_km1(:), p_k(:,grounding_line), T_km1(:))  
        call gS_fun(A_k(:,grounding_line), g_1, g_2)

        x_box = x_b1_fun(T_ess, S_0, g_1)
        T_k = Tk_fun(T_km1, x_box) 

        y_box = y_pico_fun(x_box, S_km1)
        S_k = Sk_fun(y_box, S_km1)

        q_out = overturn_fun(T_0, S_0, T_k, S_k)
        
        pico_out(:, grounding_line)%melt        = melting_pico_fun(S_k, T_k, p_k(:,grounding_line))
        pico_out(:, grounding_line)%q           = q_out
        pico_out(:, grounding_line)%T_k_out     = T_k
        pico_out(:, grounding_line)%S_k_out     = S_k

        T_km1 = T_k
        S_km1 = S_k
        
        ! rest of the elements:
        do indx_k = grounding_line+1, UBOUND(p_k,dim=2) ! assumes that shelf front is at index = 1 or any low value, not maximum
       
           !print *, '--- BOX NB : ',indx_k
  
           T_ess = T_star_fun(S_km1(:), p_k(:,indx_k), T_km1(:))
           call gS_fun(A_k(:,indx_k), g_1, g_2)

           x_box = x_bk_fun(T_ess, g_1, g_2, q_out, S_km1)
           T_k = Tk_fun(T_km1, x_box)

           y_box = y_pico_fun(x_box, S_km1)
           S_k = Sk_fun(y_box, S_km1)
           
           !pico_out(:, indx_k)%melt       = melting_pico_fun(S_k, T_k, p_k(:,indx_k))
           pico_out(:, indx_k)%q          = q_out
           pico_out(:, indx_k)%T_k_out    = T_k
           pico_out(:, indx_k)%S_k_out    = S_k

           T_km1 = T_k
           S_km1 = S_k
           
        enddo
       
        do indx_k = grounding_line+1, UBOUND(p_k,dim=2) ! assumes that shelf front is at index = 1 or any low value, not maximum
           pico_out(:, indx_k)%melt       = melting_pico_fun(pico_out(:, indx_k)%S_k_out, pico_out(:,indx_k)%T_k_out, p_k(:,indx_k))
        enddo


end function pico_go_array


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----|
! Temperature star :
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----|
pure elemental function T_star_fun(Sk_m1, pk, Tk_m1) result(T_s)

        real(kind=dp), INTENT(IN) :: Sk_m1
        real(kind=dp), INTENT(IN) :: pk
        real(kind=dp), INTENT(IN) :: Tk_m1

        real(kind=dp)             :: T_s

        T_s = a_pico * Sk_m1 + b_pico - c_pico * pk - Tk_m1

        ! output : T_s
        return
end function T_star_fun

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----|
! g1_pico & g2_pico analysis
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----|

pure elemental subroutine gS_fun(Ak, g1_pico, g2_pico)


        real(kind=dp), INTENT(IN) :: Ak
        
        real(kind=dp), INTENT(OUT):: g1_pico
        real(kind=dp), INTENT(OUT):: g2_pico
                         
        g1_pico = Ak*gammae             
        g2_pico = g1_pico/nulambda

        return
end subroutine gS_fun

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----|
! x in box 1 :
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----|

pure elemental function x_b1_fun(T_s, S_0bis, g1_pico) result(x_b1_pico)
        real(kind=dp), INTENT(IN) :: T_s
        real(kind=dp), INTENT(IN) :: S_0bis
        real(kind=dp), INTENT(IN) :: g1_pico

        real(kind=dp)             :: x_b1_pico

        real(kind=dp)             :: crbsa
        real(kind=dp)             :: s_pico

        s_pico = S_0bis/nulambda
        crbsa = C_b1_pico*rho0*(beta_pico*s_pico-alpha_pico)
        x_b1_pico = - (g1_pico/(2*crbsa))+ sqrt((g1_pico/(2*crbsa))**2 - (g1_pico*T_s)/crbsa)

        ! outputs : x_b1_pico
        return
end function x_b1_fun

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----|
! temperature in box k (including the first one) :
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----|

pure elemental function Tk_fun(Tk_m1, x_bk_pico) result(Tk_pico)
        real(kind=dp), INTENT(IN) :: x_bk_pico
        real(kind=dp), INTENT(IN) :: Tk_m1

        !--- dmr output element
        real(kind=dp)             :: Tk_pico
        
        Tk_pico = Tk_m1 - x_bk_pico

        ! output : Tk_pico
        return
end function Tk_fun

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----|
! y to find salinity in all boxes :
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----|

pure elemental function y_pico_fun(x_pico, Sk_m1) result(y_pico)

        real(kind=dp), INTENT(IN) :: x_pico
        real(kind=dp), INTENT(IN) :: Sk_m1

        real(kind=dp)             :: y_pico
                        
        y_pico = (Sk_m1 * x_pico)/nulambda
   
        ! output : y_pico
        return
end function y_pico_fun

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----|
! Salinity in all boxes with y :
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----|

pure elemental function Sk_fun(y_pico, Sk_m1) result(Sk_pico)

        real(kind=dp), INTENT(IN) :: y_pico
        real(kind=dp), INTENT(IN) :: Sk_m1
                        
        real(kind=dp)             :: Sk_pico

        Sk_pico = Sk_m1 - y_pico

        ! output : Sk_pico
        return
end function Sk_fun

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----|
! mk melt in all boxes : (Put this function at the top because that is the one your are interested in)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----|

pure elemental function  melting_pico_fun(Sk_pico, Tk_pico, pk) result(mk_pico)
        real(kind=dp), INTENT(IN) :: Sk_pico
        real(kind=dp), INTENT(IN) :: Tk_pico
        real(kind=dp), INTENT(IN) :: pk

        real(kind=dp)             :: temp
        real(kind=dp)             :: mk_pico

        temp = a_pico*Sk_pico + b_pico - c_pico*pk - Tk_pico
        mk_pico = - (gammae/nulambda)*temp

        ! output : mk_pico
        return
end function melting_pico_fun

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----|
! q the overturning circulation :
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----|

pure elemental function overturn_fun(T_0bis, S_0bis, Tk_pico, Sk_pico) result(q_pico)

        real(kind=dp), INTENT(IN) :: T_0bis
        real(kind=dp), INTENT(IN) :: S_0bis
        real(kind=dp), INTENT(IN) :: Tk_pico
        real(kind=dp), INTENT(IN) :: Sk_pico
 
        real(kind=dp)             :: rho_box1
        real(kind=dp)             :: rho_box0
        
        !--- dmr output element
        real(kind=dp)             :: q_pico

        rho_box1 = rho0*(1 - alpha_pico*(Tk_pico-T_rho) + beta_pico*(Sk_pico-S_rho)) ! Here Tk_pico and Sk_pico must be T1 and S1
        rho_box0 = rho0*(1 - alpha_pico*(T_0bis-T_rho) + beta_pico*(S_0bis-S_rho))
        q_pico = C_b1_pico * (rho_box0 - rho_box1)

        ! output : q_pico
        return
end function overturn_fun

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----|
! x_k to find temperature in box k (except the first one) :
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----|

pure elemental function x_bk_fun(T_s, g1_pico, g2_pico, q_pico, Sk_m1) result(x_bk_pico)
        real(kind=dp), INTENT(IN) :: T_s
        real(kind=dp), INTENT(IN) :: g1_pico
        real(kind=dp), INTENT(IN) :: g2_pico
        real(kind=dp), INTENT(IN) :: q_pico
        real(kind=dp), INTENT(IN) :: Sk_m1


        !--- dmr output element
        real(kind=dp)             :: x_bk_pico
        
        x_bk_pico = (-g1_pico*T_s)/(q_pico+g1_pico-g2_pico*a_pico*Sk_m1)

        ! output : x_bk_pico
        return
end function x_bk_fun

end module pico_functions

