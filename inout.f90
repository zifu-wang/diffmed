module inout
  use varGlobales
  use lire_med

  implicit none

  character(len=32) :: mesh, medsolref, medsolcomp, medsolf
  character(len=20) :: name_champ, name_solf

contains


subroutine lire_maillage_utilise

    call Lire_Maillage(mesh, node, elt, nbn, nbe)
    
end subroutine lire_maillage_utilise
  
 
subroutine lire_soli
    print*,"solution ref=",medsolref,"champ=",name_champ
    call Lire_Solution(medsolref,name_champ, solref, nsolref)
    print*,"solution a comparer=",medsolcomp,"champ=",name_champ
    call Lire_Solution(medsolcomp,name_champ, solcomp, nsolcomp)   

end subroutine lire_soli


subroutine ecrire_solf
    call ecrire_solution(mesh,medsolf,name_solf, solf, nsolf, typsolf, 0)
    
end subroutine ecrire_solf

subroutine lire_param

    open(unit=99,file='param_comp', action='read', position = 'rewind')

    read(99,*) mesh
    print*, "maillage = ", mesh
    read(99,*) medsolref
    read(99,*) medsolcomp
    read(99,*) name_champ
    read(99,*) medsolf
    read(99,*) name_solf
   
end subroutine lire_param

end module inout
