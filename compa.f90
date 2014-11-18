
program comparaison

  use varGlobales
  use inout
  implicit none
  
  integer:: i,j
 

  print*, " "
  print*, "**********************************"
  print*, "*    Comparaison des 2 champs    *"
  print*, "**********************************"

 
  call lire_param
  call lire_maillage_utilise
  print*, " nbr noeuds = ", nbn,", nbr elements = ", nbe 

  call lire_soli 
    if (nsolref.NE.nsolcomp)then
            print*,"nombre de valeurs a comparer: NON IDENTIQUE"
            stop
    end if
print*,"nsolref:",nsolref,",nbr valeurs ref: ",size(solref)
print*,"nsolcomp:",nsolcomp,",nbr valeur a comparer: ",size(solcomp)

    if (nsolref==nbe)then
            typsoli=0
    elseif(nsolref==nbn)then
            typsoli=1
    else
            print*,"type de valeurs INCONNU"
            stop
    end if

!do i=1,nsolref
!   j=3*(i-1)+1
!write(1111,*) solref(j),solref(j+1),solref(j+2) 
!write(2222,*) solcomp(j),solcomp(j+1),solcomp(j+2) 
!end do  

  typsolf=typsoli
  nsolf=nsolref
  allocate(solf(nsolf*3)) 
  solf=0.0d0

  do i=1,nsolf*3
    solf(i)=solcomp(i)-solref(i)
  end do



  call ecrire_solf

  deallocate(solf)
end program comparaison


