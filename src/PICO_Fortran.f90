!! ---|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----|
!   Copyright [2021] [Didier M. Roche, Pepijn Bakker, Maxence Menthon]
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

! create the variables T0, S0 etc. with the proper dimensions
! check go_array the inputs 
! test_output pico_go_array

module PICO_Fortran
  
 ! use pico_physics,             only: pico_output, pico_go_array
  use pico_functions,          only: pico_output, pico_go_array
  use pico_geometry_functions, only: pico_output_geo, pico_geometry_run, boxes_mask, var_x, var_y
  use global_constants_mod,    only: dp, ip!, g ! use g from here
  use rea_input_geom,          only: read_input_geom ! to read netcdf files
  use ncio,                    only: nc_write_attr, nc_create, nc_read, nc_read_attr, nc_write, nc_size, nc_write_dim

  implicit none


  contains


  subroutine run_PICO_geom_ISOMIP

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----|
  
 ! INPUTS GEOMETRY : 
  !real(kind=dp), dimension(:,:), allocatable :: low_surf
  real(kind=dp), dimension(:,:), allocatable :: draft_ice
  !real(kind=dp), dimension(:,:), allocatable :: flo_mask
  real(kind=dp), dimension(:,:), allocatable :: grd_mask
  real(kind=dp), dimension(:,:), allocatable :: oce_mask
  integer(kind=ip), parameter                :: nregions = 1
  integer(kind=ip), parameter                :: nboxes = 5
  real(kind=dp), parameter                   :: rho0 = 1033.0_dp ! value taken in the Python code
  real(kind=dp), parameter                   :: g = 9.81_dp ! value taken in the Python code >> use from global constants instead

  real(kind=dp), dimension(:,:), allocatable :: melt_boxes

  
  integer                                    :: i

  !real(kind=dp), dimension(:), allocatable   :: var_x
  !real(kind=dp), dimension(:), allocatable   :: var_y
  !real(kind=dp), dimension(:,:), allocatable :: size_var 
  
  !real(kind=dp), dimension(:), allocatable   :: nb_max_box
  !real(kind=dp), dimension(:), allocatable   :: nb_boxes
  !real(kind=dp), dimension(:,:), allocatable   :: avg_p_box
  !real(kind=dp), dimension(:,:), allocatable   :: area_box

  !integer(kind=ip) :: size_x, size_y
 
  ! INPUTS PHYSICS : 
  type(pico_output),     dimension(nregions,nboxes)  :: pico_phy_output
  type(pico_output_geo), dimension(nregions,nboxes)  :: pico_geo_main ! nb_boxes not supposed to be know here 

  real(kind=dp), dimension(nregions)           :: T_0in = 0.691885613013199_dp 
  real(kind=dp), dimension(nregions)           :: S_0in = 34.6043782936937_dp 
  !real(kind=dp), dimension(nregions,nboxes)    :: p_kin
  !real(kind=dp), dimension(nregions,nboxes)    :: A_kin 
  !p_kin(nregions,1:nboxes)      =  (/ 4022924.228457864_dp, 3262169.7641403656_dp, 2604339.5094784987_dp, 1922704.4817063196_dp, 1748411.9167280812_dp/) ! ISOMIP values 
  !p_kin(nregions,1:nboxes)      = (/9696931.2931_dp, 8648614.3965_dp, 7600297.5_dp, 6551980.6034_dp, 5503663.7068_dp/)   ! DIM == [nregions,nboxes] 
  !A_kin(nregions,1:nboxes)      = (/ 1.353e+09_dp, 1.505e+09_dp , 1.786e+09_dp, 2.5595e+09_dp, 5.6775e+09_dp/) ! ISOMIP values
  !A_kin(nregions,1:nboxes)      = (/1.6646849e+09_dp, 1.99762188e+09_dp, 1.99762188e+09_dp, 1.99762188e+09_dp, 1.6646849e+09_dp/)     ! DIM == [nregions,nboxes]


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----|
  ! READ ISOMIP+ DATA :
  ! IN : netcdf adress ; variable name
  ! OUT : variables read
  !up_surf = read_input_geom("data/Ocean1_input_geom_v1.01.nc", "upperSurface")
  draft_ice = read_input_geom("data/Ocean1_input_geom_v1.01.nc", "lowerSurface")
  !bed_topo = read_input_geom("data/Ocean1_input_geom_v1.01.nc", "bedrockTopography")
  !flo_mask = read_input_geom("data/Ocean1_input_geom_v1.01.nc", "floatingMask")
  grd_mask = read_input_geom("data/Ocean1_input_geom_v1.01.nc", "groundedMask")
  oce_mask = read_input_geom("data/Ocean1_input_geom_v1.01.nc", "openOceanMask")

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----|
  ! RUN PICO GEOMETRY:
  pico_geo_main = pico_geometry_run(grd_mask, oce_mask, draft_ice)

  ! print outputs geometry : 
  print*, 'outputs geometry : ', pico_geo_main
  print*, 'outputs geometry 1 element : ', pico_geo_main(:,:)%avg_p_box_out

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----|

  print *, '--- INPUTS PHYSICS : '
  print *, '[T_0in] : ', T_0in 
  print *, '[S_0in] : ', S_0in
  !!print *, '[p_kin] : ', p_kin
  !!print *, '[A_kin] : ', A_kin
  !print *, '[avg_p_box] : ', avg_p_box
  !print *, '[area_box] : ', area_box
  !print *, '[nb_boxes] : ', nb_boxes
  
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----|
  ! RUN PICO PHYSICS :
  !test_pico_output=pico_go_array(S_0in, T_0in, p_kin, A_kin) ! to run with given inputs p and A 
  !pico_phy_output=pico_go_array(S_0in, T_0in, avg_p_box, area_box) ! run with computed inputs p and A
  pico_phy_output=pico_go_array(S_0in, T_0in, pico_geo_main(:,:)%avg_p_box_out, pico_geo_main(:,:)%area_box_out) ! with input from new geometry

  print *, '--- OUTPUTS PHYSICS : '
  do i = 1,nboxes
        print *, 'pico_phy_output values [mk, q, T_k, S_k] box nb. =',i,':', pico_phy_output(1,i)
  end do

  ! save outputs : 

  allocate(melt_boxes(LBOUND(grd_mask,dim=1):UBOUND(grd_mask,dim=1),& 
                          LBOUND(grd_mask,dim=2):UBOUND(grd_mask,dim=2)))
  where(boxes_mask == 0)
          melt_boxes = 0
  endwhere
  do i = 1, nboxes
  print*, i
     where(boxes_mask == i)
          melt_boxes(:,:) = pico_phy_output(1,i)%melt * 3600*24*365
     endwhere
  enddo

  call nc_create("data/grl_isomip.nc") ! This create the netcdf file where the data are stored

  call nc_write_dim("data/grl_isomip.nc","x",x=var_x)
  call nc_write_dim("data/grl_isomip.nc","y",x=var_y)

  call nc_write("data/grl_isomip.nc", "boxes_mask", boxes_mask, dim1="x", dim2="y")
  call nc_write("data/grl_isomip.nc", "melt_boxes", melt_boxes, dim1="x", dim2="y")

  end subroutine run_PICO_geom_ISOMIP

end module PICO_Fortran
