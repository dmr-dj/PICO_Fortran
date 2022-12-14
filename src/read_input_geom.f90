module rea_input_geom

!use global_constants_mod, only : dp

implicit none

!real(kind=dp), dimension(:,:), allocatable :: var
!
!var = read_input_geom("Ocean1_input_geom_v1.01.nc", "bedrockTopography")
!
!write(*,*) "MAX = ", MAXVAL(var)
!write(*,*) "MIN = ", MINVAL(var)
private
public :: read_input_geom

contains

function read_input_geom(file_name, var_name) result (var_to_read)

use ncio, only: nc_read, nc_size
use global_constants_mod, only : dp, ip

character(len=*), intent(in)               :: file_name, var_name
character(len=1), parameter                :: nameX="x", nameY="y"
real(kind=dp), dimension(:,:), allocatable :: var_to_read

integer(kind=ip) :: size_x, size_y

size_y = nc_size(file_name, nameY)
size_x = nc_size(file_name, nameX)

allocate(var_to_read(1:size_x,1:size_y))
!print*, "x :", size_x
!print*, "y :", size_y

call nc_read(file_name,var_name,var_to_read)


end function read_input_geom

end module rea_input_geom




