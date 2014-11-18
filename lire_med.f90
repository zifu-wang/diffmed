!Module d'interfa�age MED
!

module lire_med
  !
  implicit none 
  !
  integer :: ierr=0


  !
  include 'med.hf'
  integer, parameter :: MED_TAILLE_NOM = 80
  integer, parameter :: MED_TAILLE_LNOM = 200
  integer, parameter :: MED_TAILLE_DESC = 200
  integer, parameter :: MED_TAILLE_PNOM = 16

contains 

  !! Lecture d'un maillage au format MED, 
  !! et assemblage d'un tableau de noeuds et un tableau d'elements.

  subroutine Lire_Maillage(mesh_file, vert_tab, elt_tab, nvert, nelt)
    ! Dummy arguments 
    character(len=*), intent(in) :: mesh_file 
    integer, intent(out) :: nelt, nvert
    double precision, dimension(:,:), pointer :: vert_tab
    integer, dimension(:,:), pointer :: elt_tab 
    ! Local variables
    integer :: ionum 
    ! Number of meshes
    integer :: nmaa
    !     ! Name of mesh
    character(len=MED_TAILLE_NOM) :: maa
    !     ! dimension du maillage (3D -> mdim=3) 
    integer :: mdim
    !     ! Mesh type : structure/non structure 
    integer :: typ
    !     ! Description associee au maillage
    character(len=MED_TAILLE_DESC) :: descr



    !impr! write(6,*) '   ******************************* '  
    !impr! write(6,*) '   *     LECTURE DU MAILLAGE     * '
    !impr! write(6,*) '   ******************************* '
    !impr! write(6,*) '   '


    ! Verifier que le fichier maillage est au bon format 
    call effoco(mesh_file, ierr)
    call trap_error(ierr, &
         "  Erreur effocco: mauvais format du fichier maillage ")
    call efveco(mesh_file, ierr)
    call trap_error(ierr, &
         " ERROR efveco: mesh MED incompatibility")
    
    ! Open mesh file (action = read only)
    call efouvr(ionum, mesh_file, MED_LECTURE, ierr)
   
    if (ierr/=0) then 
       stop " Error while opening mesh_file"
    endif
    ! Read the number of meshes
    call efnmaa (ionum, nmaa, ierr)
    if (ierr/=0) then 
       stop  " Erreur de lecture du nombre de mailage"
    else if (nmaa/=1) then 
       write(6,*) " Number of meshes to read: ", nmaa
       stop " Le fichier med contient plusieurs maillages: ambiguite!"
    endif
    !
    ! Read mesh informations 
    !
    call efmaai(ionum,nmaa,maa,mdim,typ,descr,ierr)
 
    if (mdim/=3) then 
       write(6,*)  " Mesh dimension: ", mdim
       stop " Maillage volumique necessaire"
    endif
    !
    if (ierr/=0) then 
       stop "Erreur dans la lecture du maillage"
    endif
    ! Read vertices
    call read_vertices(ionum, maa, vert_tab, nvert)
    ! Read Elements (an element is given by the indices of its vertices)
    call read_elements(ionum, maa, elt_tab, nelt)
    ! Read Groups
    !    call read_groups(ionum, maa, group_name, elt2group, vert2group)
    ! Close med file 
    call efferm(ionum, ierr)
    !
  end subroutine Lire_Maillage




  subroutine read_vertices(ionum, maa, vert_tab, nvert) 
    ! Dummy arguments 
    integer, intent(in) :: ionum
    character(len=MED_TAILLE_NOM), intent(in) :: maa
    double precision, dimension(:,:), pointer :: vert_tab
    integer, intent(out) :: nvert
    ! Local variables 
    integer :: ivert
    integer :: rep 
    double precision, dimension(:), allocatable :: coo
    character(len=16),dimension(3) :: nomcoo,unicoo 
    integer :: mdim
    integer :: profilsize
    integer :: typegeo, typecon  
    integer, dimension(2) :: profil  ! not used 
    !
    ! Read the total number of vertices in nvert 
    typegeo = 0  ! vertex => typegeo =0 
    typecon = 0 
    !impr! print*, ionum,maa,MED_COOR,MED_NOEUD,typegeo,typecon,nvert,ierr

    call efnema(ionum,maa,MED_COOR,MED_NOEUD,typegeo,typecon,nvert,ierr)

    !impr! print*, ionum,maa,MED_COOR,MED_NOEUD,typegeo,typecon,nvert,ierr
    if (ierr/=0) then 
       stop " Error while reading the number of vertices"
    endif
    !impr! write(6,*) '      - Number of vertices  : ', nvert
    !
    ! Allocate the array of coordinates 
    allocate(coo(nvert*3),STAT = ierr)
    if (ierr/=0) then 
       stop "Pb allocation coo"
    endif
    !
    ! Get coordinates of the vertices 
    profil(:)=0
    profilsize = 0 ! no profil is used 
    mdim = 3 ! 3D
    call efcool(ionum,maa,mdim,coo,MED_FULL_INTERLACE,   &
         &        MED_ALL,profil,profilsize,rep,nomcoo,unicoo,ierr)
       ! Dans le mode MED_FULL_INTERLACE, les coordonnees d'un sommet
       ! sont stockees à la suite : x1,y1,z1,x2,y2,z2... 
    !
    ! Allocate the array of vertices
    allocate(vert_tab(nvert,3),STAT=ierr)
    if (ierr/=0) then 
       stop " Allocation failed : vertices in read_vertices"
    endif
    !
    ! Loop over vertices
    do ivert=1, nvert

       vert_tab(ivert,:)=coo(3*(ivert-1)+1:3*ivert)
    enddo
    ! Free memory
    deallocate(coo)
    !
  end subroutine read_vertices
!
  subroutine read_elements(ionum, maa, elt_tab, nelt) 
    ! Dummy arguments
    integer, intent(in) :: ionum
    character(len=*), intent(in) :: maa
    integer, dimension(:,:), pointer :: elt_tab 
    integer, intent(out) :: nelt 
    ! Local variables
    integer :: ielt, i
    integer :: nvert
    integer :: typegeo, typecon 
    integer, parameter :: mdim=3 
    integer, dimension(:), allocatable :: connectivity
    integer, dimension(:), allocatable :: numelem,numvert,vertnum
    integer, dimension(1) :: profil  ! not used
   
    ! Initialize ierr to 0 
    ierr=0
    ! Read the total number of elements (tetrahedra)
    call efnema(ionum,maa,MED_CONN,MED_MAILLE,MED_TETRA4,&
         MED_NOD,nelt,ierr)
    call trap_error(ierr, " Reading the number of tetrahedra in med_interface_module")
    !impr! write(6,*) '      - Number of tetrahedra: ', nelt

    !
    ! Read the total number of vertices in nvert 
    !

    typegeo = 0  ! vertex => typegeo =0 
    typecon = 0 
    call efnema(ionum,maa,MED_COOR,MED_NOEUD,typegeo,typecon,nvert,ierr)
    call trap_error(ierr, " Reading the number of vertices in read_elements")
    ! 
    ! Allocation 
    allocate(connectivity(nelt*4),STAT=ierr)
    call trap_error(ierr,'Allocation of connectivity in med_interface_module')
    allocate(numelem(nelt), STAT=ierr)
    call trap_error(ierr,'Allocation of numelem in med_interface_module')
    allocate(numvert(nvert), STAT=ierr)
    call trap_error(ierr,'Allocation of numvert in med_interface_module')
    !
    ! Read vertices labels 
    ! Type geometrique des entites dont on veut lire les labels
    typegeo=0 ! vertex => typegeo =0 
    ! Lecture des labels des sommets du maillage 
    call efnuml(ionum,maa,numvert,nvert,MED_NOEUD,typegeo,ierr)
    call trap_error(ierr,"Problem when reading vertices labels")
    !
    ! Read a list of four vertices labels per element in "connectivity" 
    ! (connectivite nodale)
    call efconl(ionum, maa, mdim, connectivity,MED_FULL_INTERLACE,&
         profil,MED_ALL,MED_MAILLE,MED_TETRA4,MED_NOD,ierr)
    ! Translate vertices labels into indices. The translation is in-place 
!!$    allocate(vertnum(nvert), STAT=ierr)
    call trap_error(ierr,'Allocation of numvert in med_interface_module')
!!$    do i=1,nvert
!!$       write(6,*) ' ', i,'>', numvert(i) ,','
!!$    enddo
!!$    do i=1,nelt*4
!!$       connectivity(i)=vertnum(connectivity(i))
!!$    enddo
    !
    ! Allocate the array of elements 
    nullify(elt_tab)
    allocate(elt_tab(nelt,4),STAT=ierr)
    ! Loop over elements 
    do ielt=1, nelt
      elt_tab(ielt,:)=connectivity((ielt-1)*4+1:ielt*4)
    enddo
    ! Free memory 
    deallocate(connectivity,numelem, numvert)
    !
  end subroutine read_elements


!_______________________________________________________________________________________
! Lecture d'un champs 


  subroutine Lire_Solution(mesh_file, field_name, vfield, nval)

    character(len=*), intent(in) :: mesh_file
    character(len=*), intent(in) :: field_name
    real(kind=8), dimension(:), pointer :: vfield
    !nombre de valeurs 
    integer, intent(out) :: nval 
    ! Local variables
    integer :: ionum, i, typ2, ncomp
    character(len=32) :: cha
    !!character(len=MED_TAILLE_NOM) :: cha
    character(len=MED_TAILLE_PNOM) :: comp, unit
    !     ! Name of mesh
    character(len=80) :: maa
    !     ! dimension du maillage (3D -> mdim=3) 
    integer :: mdim
    !     ! Mesh type : structure/non structure 
    integer :: typ
    !     ! Description associee au maillage
    character(len=MED_TAILLE_DESC) :: descr
    !pour debuggage:     double precision, dimension(:), allocatable   :: vfield2
    !     ! Nb de maillages et de champs dans le fichier,
    integer :: nmaa, ncha,nelt
    character(len=MED_TAILLE_NOM) :: profilname
    character(len=MED_TAILLE_NOM) :: locname 

    !impr! write(6,*) '   ************************************************************* '  
    !impr! write(6,*) "Lecture du champ ",field_name ,"a partir du fichier ",mesh_file
    !impr! write(6,*) '   ************************************************************* '  


    ! Verifier que le fichier maillage est au bon format !!!!!!! 
    call effoco(mesh_file, ierr)
    call trap_error(ierr, &
         " ERROR effoco: Mauvais format du fichier de champ")
    call efveco(mesh_file, ierr)
    call trap_error(ierr, &
         "ERROR efveco: field MED incompatibility ")
    !
    ! Open mesh file (action = read only)
    call efouvr(ionum, mesh_file, MED_LECTURE, ierr)
    !
    if (ierr/=0) then 
       stop " Error while opening mesh_file"
    endif
    
    ! Read the number of meshes
    call efnmaa (ionum, nmaa, ierr)
    if (ierr/=0) then 
       stop  " Erreur de lecture du nombre de mailage"
    else if (nmaa/=1) then 
       write(6,*) " Number of meshes to read: ", nmaa
       stop " Le fichier med contient plusieurs maillages: ambiguite!"
    endif
    !!!!!!!!!!! fin verification maillage au bon format   
    !
    ! Read mesh informations 
    !
    !impr! print*, "nombre de maillages =", nmaa
    call efmaai(ionum,nmaa,maa,mdim,typ,descr,ierr)
    ! EFMAAI nous renseigne 
    
    ! Lire le nombre de champs dans un fichier 
    i = 0
    call efncha(ionum,i,ncha,ierr)
    !on pourra faire un test et exiger un seul champ
    print*, "nombre de champs =", ncha
    if (ierr/=0) then 
       stop  " Erreur de lecture du nombre de mailage"
    endif

    !! lire les informations sur le champ
    do i = 1, ncha
       !impr! print*, i
       call efchai(ionum,i, cha, typ2, comp, unit, ncomp, ierr)

       !impr! print*, "info sur champ numero", i
       !impr! print*, cha 
       !impr! print*, comp
       !impr! print*, unit
       !impr! print*, ncomp
       !impr! print*, "type =", typ2, "(", MED_FLOAT64, MED_INT32, MED_INT64,")"

    end do

    ! Read the total number of elements (tetrahedra)
    call efnema(ionum,maa,MED_CONN,MED_MAILLE,MED_TETRA4,&
         MED_NOD,nelt,ierr)
    call trap_error(ierr, " Reading the number of tetrahedra in Lire_Solution")
    !impr! write(6,*) '      - Number of tetrahedra: ', nelt

    !Verification que le nombre de valeurs est egal au nombre d'elements
    call efnval(ionum,field_name,MED_MAILLE,MED_TETRA4,MED_NOPDT,MED_NONOR,maa,&
         MED_GLOBAL,nval,ierr)
!    print*,'nb val',nval

    if(nval > 0)then

        print*, "champ aux elements"
       call trap_error(ierr, " Reading the number of values of field in Lire_Solution")
       !impr! write(6,*) '      - Number of values in the field ' ,field_name,' : ', nval
       !Allocaton de memoire pour vfield
       allocate(vfield(nval),STAT = ierr)
       if (ierr/=0) then 
          stop "Pb allocation coo" 
       endif

       ! Lire les valeurs du champ vfield
       
       
       call efchal(ionum,maa,field_name,vfield,MED_FULL_INTERLACE,MED_ALL,locname,&
            profilname,MED_GLOBAL,MED_MAILLE,MED_TETRA4,MED_NOPDT,MED_NONOR,ierr)
       

    else

       
        print*, "champ aux noeuds"
       !! champ aux noeuds
        call efnval(ionum,field_name,MED_NOEUD,MED_POINT1,MED_NOPDT,MED_NONOR,maa,&
         MED_GLOBAL,nval,ierr)

        call trap_error(ierr, " Reading the number of values of field in Lire_Solution")
       !impr! write(6,*) '      - Number of values in the field ' ,field_name,' : ', nval
       !Allocaton de memoire pour vfield
       allocate(vfield(3*nval),STAT = ierr)
       if (ierr/=0) then 
          stop "Pb allocation coo" 
       endif

       ! Lire les valeurs du champ vfield
 !     print*,'nb val',nval 
       
       call efchal(ionum,maa,field_name,vfield,MED_FULL_INTERLACE,MED_ALL,locname,&
            profilname,MED_GLOBAL,MED_NOEUD,MED_POINT1,MED_NOPDT,MED_NONOR,ierr)

    end if

       
!      subroutine efchal(fid,maa,cha,val,
!     1                  modswt,numco,locname,profil,pflmod,
!     1                  typent,typgeo,
!     1                  numdt, numo,cret)
!      character *(*) cha,maa,locname,profil
!      integer fid,val(*),n,typent,typgeo,cret
!      integer modswt,numco,pflmod,numdt,numo
! #  Paramètres :
! 
!     * fid               (IN) : descripteur du fichier.
!     * maa            (IN/OUT) :
!           o IN : Nom du maillage sur lequel porte les resultats.
!           o IN/OUT :
!                1. IN : MED_NOREF
!                2. OUT : nom du maillage par defaut. 
!     * cha             (IN) : nom du champ où lire les valeurs.
!     * val          (OUT) : tableau des valeurs.
!     * modswt     (IN) : mode de stockage des donnees dans "val".
!     * numco       (IN) : numero de composante à selectionner pour la lecture s'il y a lieu.
!     * profil     (OUT) : nom du profil utilise.
!     * pflmod      (IN) : Indique comment lire les informations en memoire.
!     * typent       (IN) : type d'entite sur lequel porte les resultats.
!     * typgeo       (IN) : type geometrique de l'entite s'il y a lieu.
!     * numdt       (IN) : numero du pas de temps s'il y a lieu.
!     * numo         (IN) : numero d'ordre s'il y a lieu.
! 
! # Code retourne : 0 si reussite, -1 sinon. 
    if (ierr/=0) then 
       stop  " Erreur de lecture du nombre de mailage"
    else 
       !impr! write(6,*) " Lecture du champ ",field_name ," reussie "
    endif
    call efferm(ionum, ierr)
  end subroutine Lire_Solution


 subroutine Ecrire_Solution(mesh_file,resu_file,field_name, vfield, nval, typsolf, typgauss)
     character(len=*), intent(in) :: mesh_file
     character(len=*), intent(in) :: resu_file
     character(len=*), intent(in) :: field_name
     real(kind=8), dimension(:), pointer :: vfield
     integer, intent(in) :: typsolf, typgauss
     integer, intent(in) :: nval !nombre de valeurs 
      ! Local variables
     integer :: ionum
!     ! Name of mesh
     character(len=MED_TAILLE_NOM) :: maa
!     ! dimension du maillage (3D -> mdim=3)
     integer :: mdim
!     ! dimension de la solution: scalaire:
     integer, parameter :: ncomp = 1
!     ! Mesh type : structure/non structure 
     integer :: typ
!     ! Description associee au maillage
     character(len=MED_TAILLE_DESC) :: descr
!     character(len = MED_TAILLE_PNOM) :: comp
!     character(len = MED_TAILLE_PNOM) :: unit
     character(len = MED_TAILLE_PNOM),dimension(3) :: comp
     character(len = MED_TAILLE_PNOM),dimension(3) :: unit
     character(len= MED_TAILLE_PNOM) :: dtunit
!     ! 
     integer :: dt 

     ! Pour les points de gauss
     integer :: ngauss1, i, j, nval2
     real(kind=16), dimension(:), allocatable :: refcoo1, gscoo1, wg1
     real(kind=16), dimension(:), allocatable :: vfield2
     character*32 :: gauss1
     real :: aa, bb, pp

!     !
   !impr! write(6,*) '   ******************************* '  
   !impr! write(6,*) "Ecriture du champ ", field_name ," dans le fichier fichier ",resu_file
   !impr! write(6,*) '   ******************************* '  


    ! Preparation du fichier resultat
    call open_visu_file_from_mesh( &
         mesh_file, resu_file, maa, ionum)

    call trap_error(ierr, &
         " Allocation problem in module visu_field_projection_module")
    dtunit = ''
    dt = 0
    comp(1)='DX';comp(2)='DY';comp(3)='DZ'
    unit= ''
    call efchac(ionum,field_name , MED_FLOAT64,comp,unit,3,ierr)
    if (ierr/=0) then 
       stop " Erreur au retour d'efchac dans losses_into_med"
    endif

    if(typsolf == 1) then


       !impr! print*, "Ecriture d'un champ au noeud"
       call efchae(ionum,maa,field_name,vfield,&
            MED_FULL_INTERLACE,nval,MED_NOGAUSS, MED_ALL,&
            MED_NOPFL,MED_NO_PFLMOD,MED_NOEUD,MED_POINT1,MED_NOPDT,&
            dtunit,0.,MED_NONOR,ierr)

      



       if (ierr/=0) then 
          stop " Erreur au retour d'efchae dans losses_into_med"
       else
          !impr! write(6,*) " Ecriture du champ ",field_name ," dans le fichier ",resu_file," reussie "
       endif

    else
 print*, "typgauss =", typgauss
       !impr! print*, "Ecriture d'un champ aux elems"
!!! SI champ par elem
       if(typgauss == 0) then

          call efchae(ionum,maa,field_name,vfield,&
               MED_FULL_INTERLACE,nval,MED_NOGAUSS, MED_ALL,&
               MED_NOPFL,MED_NO_PFLMOD,MED_MAILLE,MED_TETRA4,MED_NOPDT,&
               dtunit,0.,MED_NONOR,ierr)
          if (ierr/=0) then 
             stop " Erreur au retour d'efchae dans losses_into_med"
          else
             !impr! write(6,*) " Ecriture du champ ",field_name ," dans le fichier ",resu_file," reussie "
          endif

       else


!!! SI champ par points de Gauss

          allocate(refcoo1(12), gscoo1(12), wg1(4))
          ngauss1 = 4
          aa = (5.-sqrt(5.))/20.
          bb = (5. + 3*sqrt(5.))/20.
          pp = 1/24.
          gauss1 = '4points'

          refcoo1(1) = 0.
          refcoo1(2) = 1. 
          refcoo1(3) = 0.
          refcoo1(4) = 0.
          refcoo1(5) = 0.
          refcoo1(6) = 1.
          refcoo1(7) = 0.
          refcoo1(8) = 0.
          refcoo1(9) = 0.
          refcoo1(10) = 1.
          refcoo1(11) = 0.
          refcoo1(12) = 0.

          gscoo1(1) = aa
          gscoo1(2) = aa
          gscoo1(3) = aa
          gscoo1(4) = aa
          gscoo1(5) = aa
          gscoo1(6) = bb
          gscoo1(7) = aa
          gscoo1(8) = bb
          gscoo1(9) = aa
          gscoo1(10) = bb
          gscoo1(11) = aa
          gscoo1(12) = aa

          wg1(1) = pp
          wg1(2) = pp
          wg1(3) = pp
          wg1(4) = pp


          nval2 = nval*4
          allocate(vfield2(nval2))

          vfield2 = 0.

          do i = 1, nval 
             do j=1, 4
                vfield2(j+(i-1)*4) = vfield(i)
             end do
          end do


          call efgaue(ionum, MED_TETRA4, refcoo1, MED_FULL_INTERLACE, &
               ngauss1, gscoo1, wg1, gauss1, ierr)

          if (ierr/=0) then 
             stop " Erreur dans la definition des points de gauss"
          end if

          call efchae(ionum,maa,field_name,vfield2,&
               MED_FULL_INTERLACE,nval2, gauss1, MED_ALL,&
               MED_NOPFL,MED_NO_PFLMOD,MED_MAILLE,MED_TETRA4,MED_NOPDT,&
               dtunit,0.,MED_NONOR,ierr)

          if (ierr/=0) then 
             stop " Erreur au retour d'efchae dans losses_into_med"
          else
             !impr! write(6,*) " Ecriture du champ ",field_name ," dans le fichier ",resu_file," reussie "
          endif
       end if
!!! fin si
    endif

    ! Fin de la visu des champs sources 
    call efferm(ionum, ierr)
    call trap_error(ierr, &
       " Problem with efferm ")

    !!deallocate(vfield)
  end subroutine Ecrire_Solution


! La routine open_visu_file_from_mesh prepare un fichier med 
!! dans lequel on recopie le maillage initial du problème.
!! Ensuite on ecrira dans ce fichier les champs de resultats 
!! Developpee par Carmel 
subroutine open_visu_file_from_mesh(mesh_file, visu_file, maa, visu_ionum)
  ! Dummy arguments
  character(len=*), intent(in) :: mesh_file
  character(len=*), intent(in) :: visu_file
  ! Name of mesh
  character(len=80), intent(out) :: maa
  integer, intent(out) :: visu_ionum
  ! Local variables
  integer :: mesh_ionum
  integer :: nmaa
  ! dimension du maillage (3D -> mdim=3) 
  integer :: mdim
  ! Mesh type : structure/non structure 
  integer :: typ
  ! Description associee au maillage
  character(len=200) :: descr, descr_out
  integer :: rep 
  double precision, dimension(:), allocatable :: coo  
  character(len=16),dimension(3) :: nomcoo,unicoo 
  ! Nom (noeuds/elements) 
  character(len=16), dimension(:), allocatable :: nom
  ! Numeros (noeuds/elements) 
  integer, dimension(:), allocatable :: num
  ! Numero des familles (noeuds/elements) 
  integer, dimension(:), allocatable :: numfam
  ! Connectivite des elements
  integer, dimension(:), allocatable :: connectivity 
  ! 
  ! Indicateur de presence des noms et des numeros des sommets 
  logical :: inom, inum
  integer :: nvert, nelt 
  integer :: typegeo, typecon
  !
  ! Ouverture du fichier contenant le maillage 
  call efouvr(mesh_ionum, mesh_file, MED_LECTURE, ierr)
  call trap_error(ierr, &
       " Pb ouverture du fichier de maillage")
  ! Lecture du nombre de maillages dans ce fichier
  call efnmaa(mesh_ionum, nmaa, ierr)
  call trap_error(ierr, &
       " Error in efnmaa in module med_interface_module")
  if (nmaa/=1) then 
       stop " Only one mesh in a med file allowed !"
  endif
  ! Lecture des informations du maillage 
  call efmaai(mesh_ionum, nmaa , maa , mdim , typ, descr, ierr)
  call trap_error(ierr, &
       " Error in efmaai in module med_interface_module")
  !
  ! Ouverture du fichier de resultats pour la visualisation 
  call efouvr(visu_ionum, visu_file, MED_CREATION, ierr)
  call trap_error(ierr, &
       " Pb ouverture du fichier de maillage")
  ! Creation du maillage dans le fichier de resultats 
  descr_out='Maillage Visu'
  call efmaac(visu_ionum,maa,mdim,MED_NON_STRUCTURE, &
       descr_out,ierr)
   call trap_error(ierr, &
       " Error in efmaac in module med_interface_module")
   !!*********************
   ! Sommets du maillage
   !!*********************
   ! Lecture du nombre de sommets du maillage
   typegeo=0  ! vertex => typegeo =0 
   typecon=0 
   call efnema(mesh_ionum,maa,MED_COOR,MED_NOEUD,typegeo,typecon,nvert,ierr)
   call trap_error(ierr, &
        " Error while reading the number of vertices : efnema")
   ! Allocation des tableaux  pour la lecture des sommets du maillage
   allocate(coo(nvert*3),STAT = ierr)
   call trap_error(ierr, " Allocation of coo")
   allocate(nom(nvert), STAT = ierr)
   call trap_error(ierr, " Allocation of nom")
   allocate(num(nvert), STAT = ierr)
   call trap_error(ierr, " Allocation of num")
   allocate(numfam(nvert), STAT = ierr)
   call trap_error(ierr, " Allocation of numfam")
  
   ! Appel de la routine med de lecture des sommets 
   call efnoel(mesh_ionum,maa,mdim,coo,MED_FULL_INTERLACE,rep,nomcoo,&
        unicoo,nom,inom,num,inum,&
        numfam,nvert,ierr)
   call trap_error(ierr, &
        " Error while reading vertices : efnoel")
   ! Ecriture des sommets du maillage dans le fichier de resultat 
   call efnoee(visu_ionum,maa,mdim,coo,MED_FULL_INTERLACE,rep,nomcoo,&
        unicoo,nom,inom,num,inum,&
        numfam,nvert,ierr)
   call trap_error(ierr, &
        " Error while writing vertices : efnoee")
   ! Liberation des tableaux de travail 
   deallocate(coo, nom, num, numfam)
   !!**********************
   ! Elements du maillage
   !!**********************
   ! Read the total number of elements (tetrahedra)
    call efnema(mesh_ionum,maa,MED_CONN,MED_MAILLE,MED_TETRA4,&
         MED_NOD,nelt,ierr)
    call trap_error(ierr, &
         " Reading the number of tetrahedra in med_interface_module")
    ! Allocation des tableaux de travail 
    allocate(connectivity(nelt*4),STAT=ierr)
    call trap_error(ierr,'Allocation of connectivity in med_interface_module')
    allocate(nom(nelt), STAT = ierr)
    call trap_error(ierr, " Allocation of nom")
   allocate(num(nelt), STAT = ierr)
   call trap_error(ierr, " Allocation of num")
   allocate(numfam(nelt), STAT = ierr)
   call trap_error(ierr, " Allocation of numfam")
   ! Lecture des elements (tetraèdres à 4 noeuds) du maillage 
   call efelel(mesh_ionum,maa,mdim,connectivity,MED_FULL_INTERLACE,&
        nom,inom,num,inum,numfam,nelt,MED_MAILLE,MED_TETRA4,&
        MED_NOD, ierr)
   call trap_error(ierr, " Lecture des elements : efelel")
   ! Ecriture des elements dans le fichier de resultats 
   call efelee(visu_ionum,maa,mdim,connectivity,MED_FULL_INTERLACE,&
        nom,inom,num,inum,numfam,nelt,MED_MAILLE,MED_TETRA4,&
        MED_NOD, ierr)
   call trap_error(ierr, " Ecriture  des elements : efelee")
   !!************************
   !Familles du maillage 
   !!************************
   
   ! Liberation des tableaux de travail 
   deallocate(connectivity, nom, num, numfam)
   ! Fermeture du fichier de maillage 
   call efferm(mesh_ionum, ierr)
   call trap_error(ierr, &
       " Pb fermeture du fichier de maillage")
   !
end subroutine open_visu_file_from_mesh

  subroutine trap_error(err,message)
    ! Dummy arguments 
    integer, intent(in)          :: err
    character(len=*), intent(in) :: message
    !
    if (err /= 0) then
       write(6,*) "error code : ",err
       write(6,'(a)') message
       stop " "
    endif
    !
  end subroutine trap_error
end module Lire_med






