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


module pico_geometry_functions

use global_constants_mod, only : dp, ip ! use of sip? ! ADD G LATER ONCE OKAY WITH PYTHON <<<<<<
use ncio,                    only: nc_write_attr, nc_create, nc_read, nc_read_attr, nc_write, nc_size, nc_write_dim

implicit none 

public :: pico_output_geo, pico_geometry_run, boxes_mask, var_x, var_y

real(kind=dp), dimension(:,:), allocatable :: boxes_mask
real(kind=dp), dimension(:),   allocatable :: var_x, var_y


private 

! Shared variables between the subroutines 
! Put the variables passed from one subroutine to another here 
integer(kind=ip), parameter                :: nb_max_box  = 5.0
real(kind=dp), parameter                   :: rho0 = 1033.0_dp ! value taken in the Python code
real(kind=dp), parameter                   :: g = 9.81_dp ! value taken in the Python code >> use from global constants instead


type :: pico_output_geo
        real(dp) :: avg_p_box_out
        real(dp) :: area_box_out
        real(dp) :: nb_boxes_out
end type pico_output_geo

contains 

function pico_geometry_run(grd_mask, oce_mask, draft_ice) result(pico_out_geo)

        ! Variables again 
        real(kind=dp), dimension(:,:), INTENT(IN) :: grd_mask, oce_mask, draft_ice

        ! Put the variables passed from one subroutine to another here 
        real(kind=dp), dimension(:,:), allocatable :: grl_mask, isf_mask, is_mask, area_cells, p_grid, rel_dist_ndim! , boxes_mask
        real(kind=dp), dimension(:),   allocatable :: val_grl, val_isf, val_dist 
        real(kind=dp), dimension(:),   allocatable :: avg_p_box, area_box !, var_x, var_y
        integer(kind=ip)                           :: nb_boxes


        type(pico_output_geo), dimension(:,:), allocatable :: pico_out_geo


        allocate(grl_mask(LBOUND(grd_mask,dim=1):UBOUND(grd_mask,dim=1),& 
                          LBOUND(grd_mask,dim=2):UBOUND(grd_mask,dim=2)))
        allocate(isf_mask(LBOUND(grd_mask,dim=1):UBOUND(grd_mask,dim=1),& 
                          LBOUND(grd_mask,dim=2):UBOUND(grd_mask,dim=2)))
        allocate(is_mask(LBOUND(grd_mask,dim=1):UBOUND(grd_mask,dim=1),& 
                          LBOUND(grd_mask,dim=2):UBOUND(grd_mask,dim=2)))
        allocate(area_cells(LBOUND(grd_mask,dim=1):UBOUND(grd_mask,dim=1),& 
                          LBOUND(grd_mask,dim=2):UBOUND(grd_mask,dim=2)))  
        allocate(p_grid(LBOUND(grd_mask,dim=1):UBOUND(grd_mask,dim=1),& 
                          LBOUND(grd_mask,dim=2):UBOUND(grd_mask,dim=2)))
        allocate(boxes_mask(LBOUND(grd_mask,dim=1):UBOUND(grd_mask,dim=1),& 
                          LBOUND(grd_mask,dim=2):UBOUND(grd_mask,dim=2)))
        allocate(rel_dist_ndim(LBOUND(grd_mask,dim=1):UBOUND(grd_mask,dim=1),& 
                          LBOUND(grd_mask,dim=2):UBOUND(grd_mask,dim=2)))

        ! changed to dim=2

        allocate(val_grl(LBOUND(grd_mask,dim=2):UBOUND(grd_mask,dim=2)))! might be dim=2 instead of dim=1
        allocate(val_isf(LBOUND(grd_mask,dim=2):UBOUND(grd_mask,dim=2)))
        allocate(val_dist(LBOUND(grd_mask,dim=2):UBOUND(grd_mask,dim=2)))
        
        allocate(var_x(LBOUND(grd_mask,dim=1):UBOUND(grd_mask,dim=1)))

        allocate(var_y(LBOUND(grd_mask,dim=2):UBOUND(grd_mask,dim=2)))

        call nc_read("data/Ocean1_input_geom_v1.01.nc", "x", var_x)
        call nc_read("data/Ocean1_input_geom_v1.01.nc", "y", var_y)
 
        call nc_create("data/grl_isomip.nc") ! This create the netcdf file where the data are stored


        ! call subroutines 
        call get_grl_sbr(grd_mask, grl_mask, val_grl) ! should I use different variable names in the call?
        call get_isf_sbr(oce_mask, isf_mask, val_isf) 
        call dist_grl_isf_sbr(val_grl, val_isf, val_dist)
        call rel_dist_sbr(grd_mask, oce_mask, val_grl, val_isf, rel_dist_ndim, is_mask)
        call pico_boxes_sbr(val_dist, rel_dist_ndim, is_mask, oce_mask, nb_boxes, boxes_mask, nb_max_box)

        allocate(avg_p_box(1:nb_boxes)) ! Now nb_boxes 
        allocate(area_box(1:nb_boxes))
        allocate(pico_out_geo(1,1:nb_boxes)) ! 1 should be nregions <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        call area_cells_sbr(oce_mask, var_x, var_y, area_cells)!, nb_boxes, boxes_mask, area_box)
        call pressure_cells_sbr(draft_ice, p_grid)
        call p_area_boxes_sbr(oce_mask, nb_boxes, boxes_mask, p_grid, area_cells, avg_p_box, area_box)
        ! OUTPUTS ARE : avg_p_box and area_box (2 vectors with corresponding values for each PICO boxes)

        pico_out_geo(1,:)%avg_p_box_out = avg_p_box ! ! 1 should be nregions <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        pico_out_geo(1,:)%area_box_out = area_box
        pico_out_geo(1,:)%nb_boxes_out = nb_boxes

end function pico_geometry_run       


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----|
! Retrieve grounding line
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----|
subroutine get_grl_sbr(grd_mask, grl_mask, val_grl)
        real(kind=dp), dimension(:,:), INTENT(IN)    :: grd_mask
        real(kind=dp), dimension(:,:), INTENT(INOUT) :: grl_mask
        real(kind=dp), dimension(:),   INTENT(INOUT) :: val_grl
        !real(kind=dp), dimension(:)                  :: val_temp ! size_x, size_y
        integer(kind=ip)                             :: i, j, val ! KIND IP TO INCLUDE <<<<<<<<<<<<<<<<<<<<<<<<

        grl_mask = grd_mask
        grl_mask(:,:) = 0
        !size_x = size(grd_mask, dim=1)
        !size_y = size(grd_mask, dim=2)
        !size_y = nc_size("data/Ocean1_input_geom_v1.01.nc", "y") ! <<< IS IT OKAY TO NOT USE THIS LINE? 
        !allocate(val_dist(1:size_y))
        !allocate(val_temp(1:size_y))
        !allocate(val_grl(1:size_y))

        !write(*,*) "Loop for the Grounding line"
        i = 0
        j = 0
        val = 1
        do j = 1,80
           !write(*,*) "j = ", j 
           do while (val == 1)
              i = i+1
              !write(*,*) "i = ", i
              val = grd_mask(i,j)
              !write(*,*) "-- Value val : ", val
           end do
           grl_mask(i-1,j) = 1
           !grl_isf_mask(i-1,j) = 1
           val_grl(j) = i-1
           !val_temp(j) = i-1
           val = 1
           i = 0
        end do
        
        return 
end subroutine get_grl_sbr


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----|
! Retrieve ice shelf front
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----|
subroutine get_isf_sbr(oce_mask, isf_mask, val_isf)
        real(kind=dp), dimension(:,:), INTENT(IN)  :: oce_mask
        real(kind=dp), dimension(:,:), INTENT(OUT) :: isf_mask
        real(kind=dp), dimension(:),   INTENT(OUT) :: val_isf
        integer(kind=ip)                           :: val, i, j

        isf_mask = oce_mask
        !isf_mask = read_input_geom("data/Ocean1_input_geom_v1.01.nc", "openOceanMask") 

        isf_mask(:,:) = 0
        !size_x = size(isf_mask, dim=1)
        !size_y = size(isf_mask, dim=2)

        !size_y = nc_size("data/Ocean1_input_geom_v1.01.nc", "y")
        !allocate(val_isf(1:size_y))

        !write(*,*) "Loop for the Ice shelf front"
        i = 0
        j = 0
        val = 0
        do j = 1,80
           !write(*,*) "j = ", j 
           do while (val == 0)
              i = i+1
              !write(*,*) "i = ", i
              val = oce_mask(i,j)
              !write(*,*) "-- Value val : ", val
           end do
           isf_mask(i,j) = 1 ! i-1 ?
           !grl_isf_mask(i,j) = 1 ! i-1?
           val_isf(j) = i !i-1?
           !val_dist(j) = i-val_temp(j) !i-1? ! disctance btw grl and isf
           val = 0
           i = 0
        end do
end subroutine get_isf_sbr


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----|
! Linear distance grl - isf
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----|
subroutine dist_grl_isf_sbr(val_grl, val_isf, val_dist)
        real(kind=dp), dimension(:), INTENT(IN)  :: val_grl, val_isf
        real(kind=dp), dimension(:), INTENT(OUT) :: val_dist
        
        val_dist(:) = val_isf(:) - val_grl(:)
        return 
end subroutine dist_grl_isf_sbr


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----|
! Relative distance (eq.10 Reese et al. 2018) + distances to grl and isf + ice-shelf mask
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----|
subroutine rel_dist_sbr(grd_mask, oce_mask, val_grl, val_isf, rel_dist_ndim, is_mask)
        real(kind=dp), dimension(:,:), INTENT(IN)  :: grd_mask, oce_mask
        real(kind=dp), dimension(:),   INTENT(IN)  :: val_grl, val_isf
        real(kind=dp), dimension(:,:), INTENT(OUT) :: is_mask
        real(kind=dp), dimension(:,:), INTENT(OUT) :: rel_dist_ndim 
        real(kind=dp), dimension(:,:), allocatable :: dist_grl, dist_isf
        integer(kind=ip)                           :: i, j
        
        allocate(dist_grl(LBOUND(grd_mask,dim=1):UBOUND(grd_mask,dim=1),& 
                          LBOUND(grd_mask,dim=2):UBOUND(grd_mask,dim=2)))
        allocate(dist_isf(LBOUND(grd_mask,dim=1):UBOUND(grd_mask,dim=1),& 
                          LBOUND(grd_mask,dim=2):UBOUND(grd_mask,dim=2)))

        dist_grl = grd_mask
        dist_isf = oce_mask
        is_mask = oce_mask
        rel_dist_ndim = oce_mask
        !boxes_mask = oce_mask ! SHOULD BE IN THE GEOMETRY PART

        is_mask(:,:) = 0
        !boxes_mask(:,:) = 0 ! SHOULD BE IN THE GEO
        rel_dist_ndim(:,:) = 0 

        do i = 1,480 ! UBOUNDS
           do j = 1,80 
              dist_grl(i,j) = abs(i-val_grl(j))
              !dist_grl_temp(i,j) = (i-val_grl(j))
              dist_isf(i,j) = abs(i-val_isf(j))
              !dist_isf_temp(i,j) = (i-val_isf(j))
              !rel_dist_ndim(i,j) = dist_grl(i,j)/(dist_grl(i,j) + dist_isf(i,j)) ! (eq.10 Reese et al. 2018)

              if ((i-val_grl(j))>0 .AND. (i-val_isf(j))<0) then
                 rel_dist_ndim(i,j) = dist_grl(i,j)/(dist_grl(i,j) + dist_isf(i,j)) ! (eq.10 Reese et al. 2018)
                 is_mask(i,j) = 1
              end if
           end do 
        end do

        ! SAVING INTO NETCDF ? (distances to grl and isf?)

        return 
end subroutine rel_dist_sbr


!!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----|
!! WRITE DATA IN NETCDF : to be added in the corresponding functions higher 
!  !write_output_geom("data/grl_isomip.nc", "grl_var_name", grl_mask)
!  call nc_create("data/grl_isomip.nc") ! This create the netcdf file where the data are stored
!
!  ! Define the dims
!  !write(*,*) "dim grd_mask : ", size(grd_mask, dim=1), size(grd_mask, dim=2)
!  allocate(var_x(size_x))
!  allocate(var_y(size_y))
!  call nc_read("data/Ocean1_input_geom_v1.01.nc", "x", var_x)
!  call nc_read("data/Ocean1_input_geom_v1.01.nc", "y", var_y)
!  !write(*,*) "var_x(1,1) : ", var_x(1)
!  call nc_write_dim("data/grl_isomip.nc","x",x=var_x)
!  call nc_write_dim("data/grl_isomip.nc","y",x=var_y)
!
!  call nc_write("data/grl_isomip.nc", "grl_mask", grl_mask, dim1="x", dim2="y")
!  call nc_write("data/grl_isomip.nc", "isf_mask", isf_mask, dim1="x", dim2="y")
!  call nc_write("data/grl_isomip.nc", "grl_isf_mask", grl_isf_mask, dim1="x", dim2="y")
!  call nc_write("data/grl_isomip.nc", "val_dist", val_dist, dim1="y")
!  call nc_write("data/grl_isomip.nc", "dist_grl", dist_grl, dim1="x", dim2="y")
!  call nc_write("data/grl_isomip.nc", "dist_isf", dist_isf, dim1="x", dim2="y")
!  call nc_write("data/grl_isomip.nc", "rel_dist_ndim", rel_dist_ndim, dim1="x", dim2="y")
!  call nc_write("data/grl_isomip.nc", "is_mask", is_mask, dim1="x", dim2="y")


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----|
! PICO boxes (eq.9 and 11 Reese et al. 2018)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----|
subroutine pico_boxes_sbr(val_dist, rel_dist_ndim, is_mask, oce_mask, nb_boxes, boxes_mask, nb_max_box)
        real(kind=dp), dimension(:,:), INTENT(IN)  :: is_mask, oce_mask, rel_dist_ndim
        real(kind=dp), dimension(:),   INTENT(IN)  :: val_dist ! val_dist_global
        integer(kind=ip),              INTENT(IN)  :: nb_max_box
        real(kind=dp), dimension(:,:), INTENT(OUT) :: boxes_mask
        integer(kind=ip),              INTENT(OUT) :: nb_boxes
        integer(kind=ip)                           :: k, i, j
        real(kind=dp)                              :: temp1, temp2, dist_max_ish, dist_max_global

        ! IN : nb_max_box ; dist_max_reg ; dist_max_ish ; dims(replace hard nb) ; rel_dist_ndim ; is_mask
        ! OUT : boxes_mask ; nb_boxes
        !allocate(dist_max_ish(1))
        !allocate(dist_max_global(1))
        !allocate(nb_max_box(1))
        !allocate(nb_boxes(1))
        !allocate(temp1(1))
        !allocate(temp2(1))
        !dist_max_global = 1.0  ! <<< previous value 
        !nb_max_box(1) = 5.0 ! CHECK THIS NUMBER ? 
        nb_boxes = 1
        dist_max_ish = MAXVAL(val_dist)
        dist_max_global = dist_max_ish !! <<<<<<<<<<< TEMPORARY FIX
        !write(*,*) "dist_max_reg (global) : ", dist_max_global ! <<<<<<< NEED TO COME FROM A FIRST CALL 
        !write(*,*) "dist_max_ish (ice-shelf) : ", dist_max_ish
        ! ERROR IF MAX_DISTANCE < LOCAL DISTANCES    <<<<<<<<<<<< TO FIX !!! -- !! ..<<>><<><><><><>??
        
        boxes_mask = oce_mask ! SHOULD BE IN THE GEOMETRY PART
        boxes_mask(:,:) = 0 ! SHOULD BE IN THE GEO

        ! nb boxes (eq.9 Reese et al. 2018)
        nb_boxes = 1 + nint(sqrt(dist_max_ish/dist_max_global)*(nb_max_box-1))
        write(*,*) "nb of boxes for this ice-shelf : ", nb_boxes
              if (nb_boxes>nb_max_box) then ! <<<<<< IF TO CHANGE/REMOVE ?
           nb_boxes = nb_max_box
        end if

        ! box assignements (eq.11 Reese et al. 2018)
        do k = 1,nb_boxes
           temp1 = 1 - sqrt((real(nb_boxes)-k+1)/real(nb_boxes))
           temp2 = 1 - sqrt((real(nb_boxes)-k)/real(nb_boxes))
           do i = 1,480 ! <<<<<<<<<<<<<<<< Change hard number to variable
                    do j = 1,80
                       if ((temp1 <= rel_dist_ndim(i,j)) .AND. (rel_dist_ndim(i,j) <= temp2) &
                       .AND. (is_mask(i,j)==1)) then
                    boxes_mask(i,j) = k
                 end if
              end do
           end do
        end do

        !call nc_write("data/grl_isomip.nc", "boxes_mask", boxes_mask, dim1="x", dim2="y")

        return
end subroutine pico_boxes_sbr


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----|
! Areas (per cell and per box)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----|
! IN : one data set (size) ; dims (var_x and var_y)
! OUT : area_cells
! Need to be done in different ways for regular and not regulat grids 
! Here it is regular grid, with only x and y in 1D as inputs 

subroutine area_cells_sbr(oce_mask, var_x, var_y, area_cells)
        real(kind=dp), dimension(:,:), INTENT(IN)  :: oce_mask
        real(kind=dp), dimension(:),   INTENT(IN)  :: var_x, var_y
        real(kind=dp), dimension(:,:), INTENT(OUT) :: area_cells
        integer(kind=ip)                           :: i,j

        !----- Start area per cells
        ! main part: 
        area_cells = oce_mask
        area_cells(:,:) = 0
        do i = 1+1,480-1 ! <<<<<<<<<<<<<< CHANGE HARD NUMBERS !!!
           do j = 1+1,80-1 
              area_cells(i,j)=((var_x(i+1)-var_x(i))/2+(var_x(i)-var_x(i-1))/2)*((var_y(j+1)-var_y(j))/2+(var_y(j)-var_y(j-1))/2)
           end do
        end do
        
        ! edges:
        i = 1
        do j = 1+1,80-1 
           area_cells(i,j) = (var_x(i+1)-var_x(i)) * ((var_y(j+1)-var_y(j))/2 + (var_y(j)-var_y(j-1))/2)
        end do
        i = 480
        do j = 1+1,80-1 
           area_cells(i,j) = (var_x(i)-var_x(i-1)) * ((var_y(j+1)-var_y(j))/2 + (var_y(j)-var_y(j-1))/2)
        end do
        j = 1
        do i = 1+1,480-1 
           area_cells(i,j) = ((var_x(i+1)-var_x(i))/2 + (var_x(i)-var_x(i-1))/2) * (var_y(j+1)-var_y(j))
        end do
        j = 80
        do i = 1+1,480-1 
           area_cells(i,j) = ((var_x(i+1)-var_x(i))/2 + (var_x(i)-var_x(i-1))/2) * (var_y(j)-var_y(j-1))
        end do
        
        ! corners:
        area_cells(1,1) = area_cells(2,1)
        area_cells(1,80) = area_cells(2,80)
        area_cells(480,1) = area_cells(479,1)
        area_cells(480,80) = area_cells(479,1)

        !call nc_write("data/grl_isomip.nc", "area_cells", area_cells, dim1="x", dim2="y")
        !----- End area per cells
        
!        !----- Start area per box
!        one_box = oce_mask
!        one_box_area = oce_mask
!        allocate(area_box(1,nint(nb_boxes(1))))
!        do i = 1,nint(nb_boxes(1))
!           one_box(:,:) = 0.0_dp
!           one_box_area(:,:) = 0.0_dp
!
!           one_box = (boxes_mask==i)
!           one_box = one_box * (-1)
!           
!           !print*, "SUM ONE_BOX : ", i,  SUM(one_box)
!           one_box_area = one_box * area_cells
!           area_box(1,i) = SUM(one_box_area)
!        end do
!        !print*, "area per box : ", area_box
!        !----- End area per box

!        call nc_write("data/grl_isomip.nc", "one_box", one_box, dim1="x", dim2="y")
!        call nc_write("data/grl_isomip.nc", "one_box_area", one_box_area, dim1="x", dim2="y")
!        call nc_write("data/grl_isomip.nc", "area_box", area_box, dim1="x", dim2="y")

         return
end subroutine area_cells_sbr


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----|
! Pressure per cells :
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----|
  ! IN : draft (ISMOMIP or GRISLI) ; rho ; g ; area_cells
  ! OUT : p_grid ; avg_p_box 
 
subroutine pressure_cells_sbr(draft_ice, p_grid)
        real(kind=dp), dimension(:,:),  INTENT(IN)  :: draft_ice
        real(kind=dp), dimension(:,:),  INTENT(OUT) :: p_grid

        !draft_ice = read_input_geom("data/Ocean1_input_geom_v1.01.nc", "lowerSurface") ! THIS IS DRAFT_ICE (to put out of this function)

        ! ADD HERE A CHECKING IF THE GIVEN VALUES ARE POSITIVES OR NEGATIVES <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        p_grid = (-draft_ice) * rho0 * g  
        !call nc_write("data/grl_isomip.nc", "p_grid", p_grid, dim1="x", dim2="y")

        return
end subroutine pressure_cells_sbr


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----|
! Pressure and area per box :
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----|
subroutine p_area_boxes_sbr(oce_mask, nb_boxes, boxes_mask, p_grid, area_cells, avg_p_box, area_box)
        real(kind=dp), dimension(:,:), INTENT(IN)  :: oce_mask, boxes_mask, p_grid, area_cells
        integer(kind=ip),              INTENT(IN)  :: nb_boxes
        real(kind=dp), dimension(:),   INTENT(OUT) :: avg_p_box, area_box
        real(kind=dp), dimension(:,:), allocatable :: one_box, one_box_p, one_box_area
        integer(kind=ip)                           :: i

        allocate(one_box(LBOUND(oce_mask,dim=1):UBOUND(oce_mask,dim=1),& 
                          LBOUND(oce_mask,dim=2):UBOUND(oce_mask,dim=2)))
        allocate(one_box_p(LBOUND(oce_mask,dim=1):UBOUND(oce_mask,dim=1),& 
                          LBOUND(oce_mask,dim=2):UBOUND(oce_mask,dim=2)))
        allocate(one_box_area(LBOUND(oce_mask,dim=1):UBOUND(oce_mask,dim=1),& 
                          LBOUND(oce_mask,dim=2):UBOUND(oce_mask,dim=2)))

        one_box = oce_mask
        one_box_p = oce_mask
        one_box_area = oce_mask
        !allocate(avg_p_box(1,nb_boxes))
        !allocate(area_box(1,nb_boxes))

        do i = 1,nb_boxes
           one_box(:,:) = 0.0_dp
           one_box_p(:,:) = 0.0_dp
           one_box_area(:,:) = 0.0_dp

           one_box = MERGE(1.0,0.0,(int(boxes_mask)==i))
!~            one_box = one_box * (-1)

           !print*, "SUM ONE_BOX : ", i,  SUM(one_box)
           one_box_p = one_box * p_grid
           one_box_area = one_box * area_cells
   
           !print*, SUM(p_grid)
           !print*, SUM(one_box)
           !print*, SUM(area_cells)
           !print*, SUM(one_box * area_cells)
           !print*, SUM(one_box_p * area_cells)
           !print*, 'ok'


           avg_p_box(i) = SUM(one_box_p * area_cells) / SUM(one_box * area_cells)
           area_box(i) = SUM(one_box_area)

           !print*, i
           !print*, avg_p_box
           !print*, "Area one box : ", area_box
        end do

        !print*, i
        print*, "average pressure per box : ", avg_p_box
        print*, "area per box : ", area_box

        !call nc_write("data/grl_isomip.nc", "one_box", one_box, dim1="x", dim2="y")
        !call nc_write("data/grl_isomip.nc", "one_box_p", one_box_p, dim1="x", dim2="y")
        !call nc_write("data/grl_isomip.nc", "one_box_area", one_box_area, dim1="x", dim2="y")
        !call nc_write("data/grl_isomip.nc", "avg_p_box", avg_p_box, dim1="x", dim2="y") ! NOT 2 DIM NOW
        !call nc_write("data/grl_isomip.nc", "area_box", area_box, dim1="x", dim2="y") ! NOT 2D NOW
        
        return
end subroutine p_area_boxes_sbr


!!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----|
!  ! RUN PICO :
!  !pico_phy_output=pico_go_array(S_0in, T_0in, p_kin, A_kin) ! to run with given inputs p and A 
!  pico_phy_output=pico_go_array(S_0in, T_0in, avg_p_box, area_box) ! run with computed inputs p and A
!
!
!!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----|
!
!  ! PRINTS RESULTS PHYSICS :
!  print *, '--- INPUTS PHYSICS : '
!  print *, '[T_0in] : ', T_0in 
!  print *, '[S_0in] : ', S_0in
!  !print *, '[p_kin] : ', p_kin
!  !print *, '[A_kin] : ', A_kin
!  print *, '[avg_p_box] : ', avg_p_box
!  print *, '[area_box] : ', area_box
!  
!  print *, '--- OUTPUTS PHYSICS : '
!  do i = 1,nboxes
!        print *, 'test_pico_output values [mk, q, T_k, S_k] box nb. =',i,':', test_pico_output(1,i)
!  end do
!
!
!end program test_program



end module pico_geometry_functions
