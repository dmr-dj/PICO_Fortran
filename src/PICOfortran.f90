module PICOfortran
  implicit none
  private

  public :: say_hello
contains
  subroutine say_hello
    print *, "Hello, PICOfortran!"
  end subroutine say_hello
end module PICOfortran
