module varGlobales
  
  implicit none

  integer :: typsoli, typsolf, nsolref, nsolcomp, nsolf
  real(kind=8), dimension(:), pointer :: solf
  real(kind=8), dimension(:), pointer :: solref
  real(kind=8), dimension(:), pointer :: solcomp
  
  double precision, dimension(:,:), pointer :: node
  integer,dimension(:,:), pointer :: elt
  integer :: nbn,nbe

end module varGlobales
