module lib

use common_parameters

implicit none

contains

subroutine get_N_bands(N_bands)
  integer(i_p) :: N_bands
!local variable
  integer(i_p) :: k_local,p_local
  logical :: local_pk_dble_counted(1:N_of_particles,5:8)

  local_pk_dble_counted=.false.
  N_bands=0
  do p_local=1,N_of_particles
    do k_local=5,8
      if (bonds(p_local,k_local) .ne. 0 .and. .not. local_pk_dble_counted(p_local,k_local)) then
        N_bands=N_bands+1
        local_pk_dble_counted(bonds(p_local,k_local),K_opposite(k_local))=.true.
        local_pk_dble_counted(p_local,k_local)=.true.
      end if
    end do
  end do


end subroutine

function  N_of_particles_with_different_deH(A_number)
  integer(i_p),intent(in) :: A_number
  integer(i_p) :: N_of_particles_with_different_deH,i_local,k_local,CorresN
  if(A_number .lt. 0 .or. A_number .gt. 4) stop 'error N_of_particles_with_different_deH'
  N_of_particles_with_different_deH=0
  do i_local=1,N_of_particles
    CorresN=0
    do k_local=5,8
      if(dehydrated(i_local,k_local))CorresN=CorresN+1
    end do
    if(CorresN .eq. A_number)N_of_particles_with_different_deH=N_of_particles_with_different_deH+1
  end do
end function

function N_of_hydrogen_on_surface()
  integer(i_p) :: N_of_hydrogen_on_surface,i_local
  N_of_hydrogen_on_surface=0
  do i_local=1,N_of_hydrogen
      if(HYDROGEN_EXISTED(i_local))then
        N_of_hydrogen_on_surface=N_of_hydrogen_on_surface+1
      end if
  end do
end function

function N_of_hydrogen_in_molecules()
  integer(i_p) :: N_of_hydrogen_in_molecules,i_local,k_local
  N_of_hydrogen_in_molecules=0
  do i_local=1,N_of_particles
    do k_local =5,8
      if(.not. dehydrated(i_local,k_local))N_of_hydrogen_in_molecules=N_of_hydrogen_in_molecules+1
    end do
  end do
end function

subroutine Numb_clusters(N_clusters,average,to_print_on_screen)
  integer(i_p),intent(out) :: N_clusters
  logical,intent(in) :: to_print_on_screen
  integer(i_p) :: p,k,p_prime,loop_index
  real(r_p),intent(out):: average
  integer(i_p),allocatable :: N_connected(:),N_elements(:)
  logical(i_p) :: double_counted(1:N_of_particles,1:N_of_particles)
  logical(i_p) :: clusters_p(1:N_of_particles,1:N_of_particles)
  logical(i_p) :: p_list(1:N_of_particles),pf_check_passed
  logical :: double_counted_center(1:N_of_particles)
  
  N_clusters=0
  p_list=.false.;loop_index=0
  clusters_p=.false.;double_counted=.false.
  double_counted_center=.false.
  
  do p=1,N_of_particles
    if(DIMER_P(p) .ne.  0 .and. .not. double_counted_center(p))then
      N_clusters=N_clusters+1
      p_list(p)=.true.
      p_prime=p
      clusters_p(N_clusters,p)=.true.
      do while(.not. p_list_empty(p_list))
        pf_check_passed=.true.
        do k=5,8
          if(bonds(p_prime,k) .ne. 0) then
            if(.not. p_list((bonds(p_prime,k))))then
              pf_check_passed=.false.
            end if
          end if
        end do

        if(pf_check_passed) then
          p_list(p_prime)=.false.
          double_counted_center(p_prime)=.true.
          clusters_p(N_clusters,p_prime)=.true.
        else if(.not. pf_check_passed) then
          ! find p',p'' are connected to p_prime
          do k=5,8
            if(bonds(p_prime,k) .ne. 0 ) then
              if(.not. p_list(bonds(p_prime,k)))then
                if(.not. double_counted_center(bonds(p_prime,k))) then
                  p_list(bonds(p_prime,k))=.true.
                end if
              end if
            end if
          end do
          double_counted_center(p_prime)=.true.
          p_list(p_prime)=.false.
          clusters_p(N_clusters,p_prime)=.true.
        end if
        !find next p_prime
        if(.not. p_list_empty(p_list)) then
          do loop_index=1,N_of_particles
            if (p_list(loop_index)) then
              p_prime=loop_index
              exit
            end if
          end do
        end if
        
      end do
    endif
  enddo
  allocate(N_connected(1:N_clusters));N_connected=0
  allocate(N_elements(1:N_clusters));N_elements=0
  
  
    !also find how many particles in there.
    do loop_index=1,N_clusters
      do p=1,N_of_particles
        if(clusters_p(loop_index,p)) then
          do k=5,8
            if(bonds(P,K) .ne. 0) then
              if(.not. double_counted(P,bonds(P,K))) then
                double_counted(P,bonds(P,K))=.true.
                double_counted(bonds(P,K),P)=.true.
                N_connected(loop_index)=N_connected(loop_index)+1
              end if
            end if
          end do
          N_elements(loop_index)=N_elements(loop_index)+1
        end if
      end do
    end do
    !find elements in of each cluster
    average=0.0
    do loop_index=1,N_clusters
      if(to_print_on_screen) write(*,'("Cluster ",i0,":",2x)',advance='no')loop_index
      if(to_print_on_screen) write(*,'("Number of Bonds:",i0,2x)',advance='no')N_connected(loop_index)
      average=average+dble(N_connected(loop_index))/dble(N_elements(loop_index))
      if(to_print_on_screen) write(*,'("Compact index :",f8.3,2x)',advance='no')dble(N_connected(loop_index))/dble(N_elements(loop_index))
      do p=1,N_of_particles
        if(clusters_p(loop_index,p)) then
          if(to_print_on_screen) write(*,'(i0,",")',advance='no')p
        end if
      end do
      if(to_print_on_screen) write(*,*)
    end do
    if (N_clusters .ne. 0) then
      average=average/dble(N_clusters)
      if(to_print_on_screen) write(*,'("Average Compact index :",f8.3,2x)')average
    else if (N_clusters .eq. 0) then
      if(to_print_on_screen) write(*,'("Average Compact index :",f8.3,2x)')0.000
    else
      stop 'error on Numb_clusters subroutine'
    end if
  
end subroutine Numb_clusters


function p_list_empty(p_list)
  logical :: p_list_empty
  integer :: list_index
  logical(i_p),intent(in) :: p_list(1:N_of_particles)
  p_list_empty=.true.
  do list_index=1,N_of_particles
    if(p_list(list_index)) p_list_empty=.false.
  end do
end function

FUNCTION N_DIMER()
integer :: N_DIMER
integer :: TEMP_P
N_DIMER=0
do TEMP_P=1,N_of_particles
    if(DIMER_P(TEMP_P) .ne. 0) then !unconnected particles
        N_DIMER=N_DIMER+1
    end if
end do
END FUNCTION N_DIMER

subroutine check_dehyd_or_creat_h(P_LOCAL,K_LOCAL,space_left,ACTION_LOCAL)
![in]P,K
![out]HYDROGEN_X(H),HYDROGEN_Y(H),HYDROGEN_X(H'),HYDROGEN_Y(H')
![out]OCCUPATION of them equal to 1

integer,intent(in) :: P_LOCAL,K_LOCAL,ACTION_LOCAL
logical,intent(out) :: space_left
integer :: AVILIBLE_SPACE
integer :: POINT_X_LOCAL,POINT_Y_LOCAL
integer :: LOOP_Y_LOCAL
integer :: LOOP_X_LOCAL
integer :: LOCAL_RANDOM_UP
integer :: H_LOCAL
integer :: i_local
integer,dimension(1:3*(l_particle+1)**2,1:2) :: LOOKUP_LOCAL

!000000000^ !000000| X-1 to X-1-L, Y to Y+L and X to X+L, Y-1 to Y-1-L
!000000000| !005160|
!000617000| !001110|
!000111000| !007180|
!000518000| !000000X
!000000000| !----->Y
!000000000Y
!x------->
AVILIBLE_SPACE=0
LOOKUP_LOCAL=0
LOOP_Y_LOCAL=0
space_left=.false.

POINT_X_LOCAL=CORNERS(X(P_LOCAL),K_LOCAL,1)
POINT_Y_LOCAL=CORNERS(Y(P_LOCAL),K_LOCAL,2)

do i_local=0,2
  do LOOP_X_LOCAL=C_H_BOUNDARY(K_LOCAL,POINT_X_LOCAL,1,1,1+i_local), &
        C_H_BOUNDARY(K_LOCAL,POINT_X_LOCAL,1,2,1+i_local)
      
    do LOOP_Y_LOCAL=C_H_BOUNDARY(K_LOCAL,POINT_Y_LOCAL,2,1,1+i_local), &
           C_H_BOUNDARY(K_LOCAL,POINT_Y_LOCAL,2,2,1+i_local)
      if(OCC(LOOP_X_LOCAL,LOOP_Y_LOCAL) .eq. 0) THEN
        AVILIBLE_SPACE=AVILIBLE_SPACE+1
        if(ACTION_LOCAL .eq. 0)then
          if(AVILIBLE_SPACE .ne. 0)then
            space_left=.true.
            return
          end if
        end if
        LOOKUP_LOCAL(AVILIBLE_SPACE,1) = LOOP_X_LOCAL
        LOOKUP_LOCAL(AVILIBLE_SPACE,2) = LOOP_Y_LOCAL
      end if
    end do
  end do
end do

!=====================================!
!========chose one_of_the_space=======!
if(ACTION_LOCAL .eq. 1)then
   if(AVILIBLE_SPACE .ne. 0) then
!==========upper area=========
      H_LOCAL=1
      do
        if(.not. HYDROGEN_EXISTED(H_LOCAL))then
          HYDROGEN_EXISTED(H_LOCAL)=.true.
          exit
        end if
          H_LOCAL=H_LOCAL+1
      end do
      
      LOCAL_RANDOM_UP = 1+nint((AVILIBLE_SPACE-1)*ran1(idum))
      HYDROGEN_X(H_LOCAL)=LOOKUP_LOCAL(LOCAL_RANDOM_UP,1)
      HYDROGEN_Y(H_LOCAL)=LOOKUP_LOCAL(LOCAL_RANDOM_UP,2)
      HYDROGEN_X(H_LOCAL)=mod1(HYDROGEN_X(H_LOCAL))
      HYDROGEN_Y(H_LOCAL)=mod1(HYDROGEN_Y(H_LOCAL))

!=========lower area===========
        !====================checking point===================
        !write(*,*)H_LOCAL
        !write(*,*)"HYDROGEN", HYDROGEN_X(H_LOCAL),HYDROGEN_Y(H_LOCAL)
        !write(*,*)mod1(HYDROGEN_X(H_LOCAL)),mod1(HYDROGEN_Y(H_LOCAL))
        !====================end checking point================
      if(OCC(HYDROGEN_X(H_LOCAL),HYDROGEN_Y(H_LOCAL)) .ne. 0) then
         STOP 'error on CTREAT HYDROGEN' !... no need for this case
      else if(OCC(HYDROGEN_X(H_LOCAL),HYDROGEN_Y(H_LOCAL)) .eq. 0) then
         OCCUPATION(mod1(HYDROGEN_X(H_LOCAL)),mod1(HYDROGEN_Y(H_LOCAL)))=1 !.... no need to use mod1 again
      end if
        !========================checking point===========================
        !write(*,'("direction at",X,I0)')K_LOCAL
        !write(*,'(A,X,2(I0,X))')'HYDROGEN AT:',HYDROGEN_X(H_LOCAL),HYDROGEN_Y(H_LOCAL)
        !write(*,'((I0,X))')OCC(HYDROGEN_X(H_LOCAL),HYDROGEN_Y(H_LOCAL))
        !=========================end checking point=======================
      end if
   else
    
end if
end subroutine check_dehyd_or_creat_h


integer function mod1(i_local)
  integer,intent(in)::i_local
  if (i_local .lt. 0)then
      mod1=N+i_local+1
  else if (i_local .gt. N) then
      mod1=i_local-N-1
  else
      mod1=i_local
  end if
end function

integer function OCC(X_LOCAL,Y_LOCAL)
  integer,intent(in)::X_LOCAL,Y_LOCAL
  OCC=OCCUPATION(MOD1(X_LOCAL),MOD1(Y_LOCAL))
end function


subroutine UPDATE_PARTICLE_SHAPE(P_LOCAL,ACTION_LOCAL)
![in] P
![out] OCCUPATION
![do]creat or annilate 3by3 particles
integer,intent(in) :: P_LOCAL,ACTION_LOCAL
integer :: LOCAL_X,LOCAL_Y,LOOP_X,loop_y
LOCAL_X=X(P_LOCAL)
LOCAL_Y=Y(P_LOCAL)
!without doing this could cause memory issues

if(ACTION_LOCAL .eq. 1) then !creat particle
   do loop_x=0,l_particle
      do loop_y=0,l_particle
         OCCUPATION(mod1(LOCAL_X+loop_x),mod1(LOCAL_Y+loop_y)) = 1
      end do
   end do
   
else if(ACTION_LOCAL .eq. 0)then !annilate particle
   do loop_x=0,l_particle
      do loop_y=0,l_particle
         OCCUPATION(mod1(LOCAL_X+loop_x),mod1(LOCAL_Y+loop_y)) = 0
      end do
   end do
else if(ACTION_LOCAL .eq. 2)then
   do loop_x=0,l_particle
      do loop_y=0,l_particle
         OCCUPATION(mod1(LOCAL_X+loop_x),mod1(LOCAL_Y+loop_y)) = 2
      end do
   end do
!==================check point================

!==================end check point============

end if
end subroutine


subroutine CHECK_OCCUPANCY(P_LOCAL,ACTION_LOCAL,P_or_PK_OCCUPIED)
integer(i_p),intent(in) :: P_LOCAL,ACTION_LOCAL
logical,intent(out)::P_or_PK_OCCUPIED
integer :: LOOP_X,LOOP_Y,LOCAL_X,LOCAL_Y
![in]X(P),Y(P),OCCUPATION
![out] P_or_PK_OCCUPIED
![do] CHECK_OCCUPANCY
!ACTION_LOCAL = 0 left corner type check space from X to X+size and Y to Y+size
!ACTION_LOCAL = 1-4 direction at 4 sides of the particle
!ACTION_LOCAL = 5-8 direction at 4 corners
LOCAL_X=X(P_LOCAL)
LOCAL_Y=Y(P_LOCAL)
P_or_PK_OCCUPIED=.false.

if (ACTION_LOCAL .eq. 0) then
    do loop_x = 0,l_particle
        do loop_y= 0,l_particle
            if (OCC(LOCAL_X+loop_x,LOCAL_Y+loop_y) .eq. 1) then
                P_or_PK_OCCUPIED=.true.
                exit
            end if
        end do
    end do
!==============sides check=============================

   else if(ACTION_LOCAL .eq. 1) then
        do loop_y=0,l_particle
            if(OCC(LOCAL_X+l_particle+1,LOCAL_Y+loop_y) .eq. 1) then
                P_or_PK_OCCUPIED=.true.
                exit
            end if
        end do

    else if(ACTION_LOCAL .eq. 2) then
        do loop_y=0,l_particle
            if(OCC(LOCAL_X-1,LOCAL_Y+loop_y) .eq. 1) then
                P_or_PK_OCCUPIED=.true.
                exit
            end if
        end do

    else if(ACTION_LOCAL .eq. 3) then
        do loop_x=0,l_particle
            if(OCC(LOCAL_X+loop_x,LOCAL_Y+1+l_particle) .eq. 1) then
                P_or_PK_OCCUPIED=.true.
                exit
            end if
        end do

    else if(ACTION_LOCAL .eq. 4) then
        do loop_x=0,l_particle
            if(OCC(LOCAL_X+loop_x,LOCAL_Y-1) .eq. 1) then
                P_or_PK_OCCUPIED=.true.
                exit
            end if
        end do
!===========corner direction check===============
    else if(ACTION_LOCAL .eq. 5) then
        !X-1, Y-1 left down
        if(OCC(LOCAL_X-1,LOCAL_Y-1) .eq. 1)then
            P_or_PK_OCCUPIED=.true.
        end if

    else if(ACTION_LOCAL .eq. 6) then
        !X-1, Y'+1 left up
        if(OCC(LOCAL_X-1,LOCAL_Y+l_particle+1) .eq. 1)then
            P_or_PK_OCCUPIED=.true.
        end if

    else if(ACTION_LOCAL .eq. 7) then
        !X'+1, Y'+1 right top
        if(OCC(LOCAL_X+l_particle+1,LOCAL_Y+l_particle+1) .eq. 1)then
            P_or_PK_OCCUPIED=.true.
        end if

    else if(ACTION_LOCAL .eq. 8) then
        !X'+1, Y-1 right down
        if(OCC(LOCAL_X+l_particle+1,LOCAL_Y-1) .eq. 1)then
            P_or_PK_OCCUPIED=.true.
        end if
    else
        stop "error on input of CHECK_OCCUPANCY"
    end if
end subroutine



integer function CORNER_POINTED(VAR,K_LOCAL,X_or_Y)
integer,intent(in) :: K_LOCAL,VAR,X_or_Y
!X_or_Y 1 for X, 2 for Y
!VAR for X or Y
!============
!00000 left down(K=5) X-1 and Y-1
!01110 left up(K=6) X-1 and Y+l_particle+1
!01110 right up(K=7) X+l_particle+1 and Y+l_particle+1
!0X110 right down(K=8) x+l_particle+1 and Y-1
!00000

if (X_or_Y .eq. 1 ) then
    if(K_LOCAL .eq. 5 .or. K_LOCAL .eq. 6) then
        CORNER_POINTED=mod1(VAR-1)
    else if(K_LOCAL .eq. 7 .or. K_LOCAL .eq. 8) then
        CORNER_POINTED=mod1(VAR+l_particle+1)
    end if
else if(X_or_Y .eq. 2) then
    if (K_LOCAL .eq. 5 .or. K_LOCAL .eq. 8) then
        CORNER_POINTED=mod1(VAR-1)
    else if(K_LOCAL .eq. 6 .or. K_LOCAL .eq. 7)then
        CORNER_POINTED=mod1(VAR+l_particle+1)
    end if

end if

end function

integer function CORNERS(VAR,K_LOCAL,X_or_Y)
integer,intent(in) :: K_LOCAL,VAR,X_or_Y
!X_or_Y 1 for X, 2 for Y
!VAR for X or Y

!00000 left down(K=5) X and Y
!01110 left up(K=6) X and Y+l_particle
!01110 right up(K=7) X+l_particle and Y+l_particle
!0X110 right down(K=8) x+l_particle+1 and Y
!00000
!================X================
if (X_or_Y .eq. 1 ) then
    if(K_LOCAL .eq. 5 .or. K_LOCAL .eq. 6) then
        CORNERS=mod1(VAR)
    else if(K_LOCAL .eq. 7 .or. K_LOCAL .eq. 8) then
        CORNERS=mod1(VAR+l_particle)
    end if
!================Y=================
else if(X_or_Y .eq. 2) then

    if (K_LOCAL .eq. 5 .or. K_LOCAL .eq. 8) then
        CORNERS=mod1(VAR)
    else if(K_LOCAL .eq. 6 .or. K_LOCAL .eq. 7)then
        CORNERS=mod1(VAR+l_particle)
    end if

end if

end function


integer function K_opposite(K_LOCAL)
integer,intent(in) :: K_LOCAL
!K=1 X+1
!K=2 X-1
!k=3 Y+1
!k=4 Y-1
!k=5 X-1 Y-1 should attach to X+l Y+l K_oppo=7
!k=6 X-1 Y+l+1 should attach to X+l Y K_oppo=8
if (K_LOCAL .eq. 1)then
    K_opposite=2
else if (K_LOCAL .eq. 2)then
    K_opposite=1
else if (K_LOCAL .eq. 3)then
    K_opposite=4
else if (K_LOCAL .eq. 4)then
    K_opposite=3
else if (K_LOCAL .eq. 5)then
    K_opposite=7
else if (K_LOCAL .eq. 6)then
    K_opposite=8
else if (K_LOCAL .eq. 7)then
    K_opposite=5
else if (K_LOCAL .eq. 8)then
    K_opposite=6
end if

end function


subroutine CHECK_CORNER(P_LOCAL,K_LOCAL,PKP_LOCAL,local_corner_corner_connectable)
![in] P,K
![out] local_corner_corner_connectable,PKP_LOCAL
![do] check whether local_corner_corner_connectable, do this check if P there before running this subroutine
integer,intent(in) :: P_LOCAL,K_LOCAL
integer,intent(out) :: PKP_LOCAL
logical,intent(out) :: local_corner_corner_connectable
integer :: P_LOOP_LOCAL,MY_LOCAL_CORNER,MY_LOCAL_CORNER_POINTED,MY_LOCAL_CORNER_Y,MY_LOCAL_CORNER_POINTED_Y


!check if CORNER AND CORNER CAN TOUCH
PKP_LOCAL=-1
local_corner_corner_connectable=.false.
do P_LOOP_LOCAL=1,N_of_particles
   if (P_LOOP_LOCAL .eq. P_LOCAL) then
      cycle
   else
      MY_LOCAL_CORNER_POINTED = CORNER_POINTED(X(P_LOCAL),K_LOCAL,1)
      MY_LOCAL_CORNER =    CORNERS(X(P_LOOP_LOCAL),K_opposite(K_LOCAL),1)
      
      MY_LOCAL_CORNER_POINTED_Y   =  CORNER_POINTED(Y(P_LOCAL),K_LOCAL,2)
      MY_LOCAL_CORNER_Y   =  CORNERS(Y(P_LOOP_LOCAL),K_opposite(K_LOCAL),2)
      if(MY_LOCAL_CORNER_POINTED .eq. MY_LOCAL_CORNER) then
         if( MY_LOCAL_CORNER_Y .eq. MY_LOCAL_CORNER_POINTED_Y )then
            local_corner_corner_connectable=.true.
            PKP_LOCAL=P_LOOP_LOCAL
            if(PKP_LOCAL .lt. 1 .or. PKP_LOCAL .gt. N_of_particles) then
               STOP 'error on CHECK CORNER'
            end if
            exit
         else
            cycle
         end if
      else
         cycle
      end if
   end if
end do
 
!........ this could have been done much simplier!
!.. you dimply need to check if there is a molecule at a particular site
!   that can connect to P,K;
!.... no need for cycle; 
end subroutine CHECK_CORNER


subroutine randomly_choose_a_sequence(sequence_of_events_local,R_c_local)
![in]R_c
![out]sequence_of_events
![do]as it's name
integer,intent(out)::sequence_of_events_local
real(r_p),intent(in) :: R_c_local

!local variable
real(r_p) :: KMC_random
real(r_p):: LOCAL_R_C

LOCAL_R_C=0
KMC_random=ran1(idum)*R_c_local
sequence_of_events_local=1
LOCAL_R_C=R(1)

do while(KMC_random .ge. LOCAL_R_C)
    sequence_of_events_local=sequence_of_events_local+1
    LOCAL_R_C=LOCAL_R_C+R(sequence_of_events_local)
end do

end subroutine

integer function C_H_BOUNDARY(K_LOCAL,POINT_LOCAL,X_or_Y,LEFT_OR_RIGHT,UP_OR_DOWN)
integer,intent(in) :: K_LOCAL,POINT_LOCAL,X_or_Y,LEFT_OR_RIGHT,UP_OR_DOWN

if(K_LOCAL .eq. 5) then
   ! for K_LOCAL =5, X+0 to X+L, Y-1 to Y-1-L and ! X-1 to X-1-L, Y to Y+L
   ! X-1 to X-1-L, Y to Y+L and X to X+L, Y-1 to Y-1-L
   if(UP_OR_DOWN .eq. 1) then
      if(X_or_Y .eq. 1) then ! X
         if(LEFT_OR_RIGHT .eq. 1) then ! left
            C_H_BOUNDARY=POINT_LOCAL
         else if(LEFT_OR_RIGHT .eq. 2)then !right
            C_H_BOUNDARY=POINT_LOCAL+l_particle
         end if
      else if(X_or_Y .eq. 2) then ! Y
         if(LEFT_OR_RIGHT .eq. 1) then !left
            C_H_BOUNDARY=POINT_LOCAL-1-l_particle
         else if(LEFT_OR_RIGHT .eq. 2) then !right
            C_H_BOUNDARY=POINT_LOCAL-1
         end if!left or right
      end if!x or y
!! X-1 to X-1-L, Y to Y+L
   else if(UP_OR_DOWN .eq. 2) then !DOWN
      if(X_or_Y .eq. 1) then ! X
         if(LEFT_OR_RIGHT .eq. 1) then ! left
            C_H_BOUNDARY=POINT_LOCAL-1-l_particle
         else if(LEFT_OR_RIGHT .eq. 2)then !right
            C_H_BOUNDARY=POINT_LOCAL-1
         end if
         
      else if(X_or_Y .eq. 2) then ! Y
         if(LEFT_OR_RIGHT .eq. 1) then !left
            C_H_BOUNDARY=POINT_LOCAL
         else if(LEFT_OR_RIGHT .eq. 2) then !right
            C_H_BOUNDARY=POINT_LOCAL+l_particle
         end if !left or right
      end if!x or y
!and X-1 to X-1-L, Y-1 to Y-1-L
   else if(UP_OR_DOWN .eq. 3) then !MIDDLE
      
      if(X_or_Y .eq. 1) then ! X
         if(LEFT_OR_RIGHT .eq. 1) then ! left
            C_H_BOUNDARY=POINT_LOCAL-1-l_particle
         else if(LEFT_OR_RIGHT .eq. 2)then !right
            C_H_BOUNDARY=POINT_LOCAL-1
         end if
      else if(X_or_Y .eq. 2) then ! Y
         if(LEFT_OR_RIGHT .eq. 1) then !left
                C_H_BOUNDARY=POINT_LOCAL-1-l_particle
             else if(LEFT_OR_RIGHT .eq. 2) then !right
                C_H_BOUNDARY=POINT_LOCAL-1
            end if
         end if
      end if

   else if(K_LOCAL .eq. 6) then
!K_LOCAL=6,X-1:X-1-l, Y:Y-l and X:X+l, Y+1:Y+1+L
! for K_LOCAL =6, X+0 to X+L, Y+1 to Y+1+L and X-1 to X-1-l, Y to Y-L
      if(UP_OR_DOWN .eq. 1) then
         if(X_or_Y .eq. 1) then ! X
            if(LEFT_OR_RIGHT .eq. 1) then ! left
                C_H_BOUNDARY=POINT_LOCAL
             else if(LEFT_OR_RIGHT .eq. 2)then !right
                C_H_BOUNDARY=POINT_LOCAL+l_particle
             end if
          else if(X_or_Y .eq. 2) then ! Y
             if(LEFT_OR_RIGHT .eq. 1) then !left
                C_H_BOUNDARY=POINT_LOCAL+1
             else if(LEFT_OR_RIGHT .eq. 2) then !down
                C_H_BOUNDARY=POINT_LOCAL+1+l_particle
             end if
          end if
! X-1 to X-1-l, Y to Y-L
       else if(UP_OR_DOWN .eq. 2) then !DOWN
          if(X_or_Y .eq. 1) then ! X
             if(LEFT_OR_RIGHT .eq. 1) then ! left
                C_H_BOUNDARY=POINT_LOCAL-1-l_particle
            else if(LEFT_OR_RIGHT .eq. 2)then !right
               C_H_BOUNDARY=POINT_LOCAL-1
            end if
            
         else if(X_or_Y .eq. 2) then ! Y
            if(LEFT_OR_RIGHT .eq. 1) then !left
               C_H_BOUNDARY=POINT_LOCAL-l_particle
            else if(LEFT_OR_RIGHT .eq. 2) then !right
               C_H_BOUNDARY=POINT_LOCAL
            end if
         end if
!and Y+1 to Y+1+L, X-1 to X-1-l
      else if(UP_OR_DOWN .eq. 3) then !MIDDLE
         if(X_or_Y .eq. 1) then ! X
            if(LEFT_OR_RIGHT .eq. 1) then ! left
               C_H_BOUNDARY=POINT_LOCAL-1-l_particle
            else if(LEFT_OR_RIGHT .eq. 2)then !right
               C_H_BOUNDARY=POINT_LOCAL-1
            end if
            
         else if(X_or_Y .eq. 2) then ! Y
            if(LEFT_OR_RIGHT .eq. 1) then !left
               C_H_BOUNDARY=POINT_LOCAL+1
            else if(LEFT_OR_RIGHT .eq. 2) then !down
               C_H_BOUNDARY=POINT_LOCAL+1+l_particle
            end if
         end if
      end if
   else if(K_LOCAL .eq. 7) then

! for K_LOCAL =7, X+1 to X+1+L, Y to Y-L and X to X-L , Y+1 to Y+1+L
      if(UP_OR_DOWN .eq. 1) then
         if(X_or_Y .eq. 1) then ! X
            if(LEFT_OR_RIGHT .eq. 1) then ! left
               C_H_BOUNDARY=POINT_LOCAL+1
            else if(LEFT_OR_RIGHT .eq. 2)then !right
               C_H_BOUNDARY=POINT_LOCAL+1+l_particle
            end if
         else if(X_or_Y .eq. 2) then ! Y
            if(LEFT_OR_RIGHT .eq. 1) then !left
               C_H_BOUNDARY=POINT_LOCAL-l_particle
            else if(LEFT_OR_RIGHT .eq. 2) then !right
               C_H_BOUNDARY=POINT_LOCAL
            end if
         end if
! X to X-L , Y+1 to Y+1+L
      else if(UP_OR_DOWN .eq. 2) then !DOWN
         if(X_or_Y .eq. 1) then ! X
            if(LEFT_OR_RIGHT .eq. 1) then ! up
               C_H_BOUNDARY=POINT_LOCAL-l_particle
            else if(LEFT_OR_RIGHT .eq. 2)then !right
               C_H_BOUNDARY=POINT_LOCAL
            end if

         else if(X_or_Y .eq. 2) then ! Y
            if(LEFT_OR_RIGHT .eq. 1) then !left
               C_H_BOUNDARY=POINT_LOCAL+1
            else if(LEFT_OR_RIGHT .eq. 2) then !right
               C_H_BOUNDARY=POINT_LOCAL+1+l_particle
            end if
         end if
!and X+1 to X+1+L, Y+1 to Y+1+L
      else if(UP_OR_DOWN .eq. 3) then
         if(X_or_Y .eq. 1) then ! X
            if(LEFT_OR_RIGHT .eq. 1) then ! left
               C_H_BOUNDARY=POINT_LOCAL+1
            else if(LEFT_OR_RIGHT .eq. 2)then !right
               C_H_BOUNDARY=POINT_LOCAL+1+l_particle
            end if
         else if(X_or_Y .eq. 2) then ! Y
            if(LEFT_OR_RIGHT .eq. 1) then !left
               C_H_BOUNDARY=POINT_LOCAL+1
            else if(LEFT_OR_RIGHT .eq. 2) then !right
               C_H_BOUNDARY=POINT_LOCAL+1+l_particle
            end if
         end if
      end if

   else if(K_LOCAL .eq. 8) then

! for K_LOCAL =8, X+1 to X+1+L, Y to Y+L and X to X-L , Y-1 to Y-1-L
      if(UP_OR_DOWN .eq. 1) then
         if(X_or_Y .eq. 1) then ! X
            if(LEFT_OR_RIGHT .eq. 1) then ! left
               C_H_BOUNDARY=POINT_LOCAL+1
            else if(LEFT_OR_RIGHT .eq. 2)then !right
               C_H_BOUNDARY=POINT_LOCAL+1+l_particle
            end if
         else if(X_or_Y .eq. 2) then ! Y
            if(LEFT_OR_RIGHT .eq. 1) then !left
               C_H_BOUNDARY=POINT_LOCAL
            else if(LEFT_OR_RIGHT .eq. 2) then !down
               C_H_BOUNDARY=POINT_LOCAL+l_particle
            end if
         end if
!X to X-L , Y-1 to Y-1-L
      else if(UP_OR_DOWN .eq. 2) then !DOWN
         if(X_or_Y .eq. 1) then ! X
            if(LEFT_OR_RIGHT .eq. 1) then ! up
               C_H_BOUNDARY=POINT_LOCAL-l_particle
            else if(LEFT_OR_RIGHT .eq. 2)then !right
               C_H_BOUNDARY=POINT_LOCAL
            end if
         else if(X_or_Y .eq. 2) then ! Y
            if(LEFT_OR_RIGHT .eq. 1) then !left
               C_H_BOUNDARY=POINT_LOCAL-1-l_particle
            else if(LEFT_OR_RIGHT .eq. 2) then !right
               C_H_BOUNDARY=POINT_LOCAL-1
            end if
         end if
!and X+1 to X+1+L,Y-1 to Y-1-L
      else if(UP_OR_DOWN .eq. 3) then
         if(X_or_Y .eq. 1) then ! X
            if(LEFT_OR_RIGHT .eq. 1) then ! left
               C_H_BOUNDARY=POINT_LOCAL+1
            else if(LEFT_OR_RIGHT .eq. 2)then !right
               C_H_BOUNDARY=POINT_LOCAL+1+l_particle
            end if
         else if(X_or_Y .eq. 2) then ! Y
            if(LEFT_OR_RIGHT .eq. 1) then !left
               C_H_BOUNDARY=POINT_LOCAL-1-l_particle
            else if(LEFT_OR_RIGHT .eq. 2) then !right
               C_H_BOUNDARY=POINT_LOCAL-1
            end if!left or right
         end if !x or y
      end if! 1,2 or 3
   end if !K_LOCAL

 end function C_H_BOUNDARY


subroutine UPDATE_OCCUPATION(P_LOCAL,K_LOCAL)
    integer,intent(in) :: P_LOCAL,K_LOCAL
    ![in]P,K
    ![out]OCCUPATION

!=====================check point=========================
!write(*,'("at kmc steps", I0)')KMC_STEPS
!write(*,'("P_",I0, "K_",I0,X,"events",I0)')P_LOCAL,K_LOCAL,LOOKUP(sequence_of_events,3)
!110 format(20(i0,X))
!
    call UPDATE_PARTICLE_SHAPE(P_LOCAL,0)

!=====================check point=========================
!    do loop_lattice=0,N
!        write(*,110) OCCUPATION(loop_lattice,:)
!    end do


    if (K_LOCAL .eq. 1) then
        X(P_LOCAL)=mod1(X(P_LOCAL)+1)
        CALL UPDATE_PARTICLE_SHAPE(P_LOCAL,0) !.... not needed
    else if (K_LOCAL .eq. 2) then
        X(P_LOCAL)=mod1(X(P_LOCAL)-1)
        CALL UPDATE_PARTICLE_SHAPE(P_LOCAL,0) !.... not needed
    else if (K_LOCAL .eq. 3) then
        Y(P_LOCAL)=mod1(Y(P_LOCAL)+1)
        CALL UPDATE_PARTICLE_SHAPE(P_LOCAL,0) !.... not needed
    else if (K_LOCAL .eq. 4) then
        Y(P_LOCAL)=mod1(Y(P_LOCAL)-1)
        CALL UPDATE_PARTICLE_SHAPE(P_LOCAL,0) !.... not needed
    else
        STOP 'error on the input of UPDATE_OCCUPATION subrountine'
    end if

    call UPDATE_PARTICLE_SHAPE(P_LOCAL,1)


!====================check point====================
!write(*,*)KMC_STEPS,'after updated'
!do loop_lattice=0,N
!write(*,110) OCCUPATION(loop_lattice,:)
!end do
!====================end check point================

end subroutine



subroutine HYDROGEN_OCCUPATION_UPDATE(H_LOCAL,K_LOCAL)
![in] series_#_of_HYDROGEN, DIRECTION TO MOVE
![out] HYDROGEN_X,HYDROGEN_Y
![do] as the name
integer,intent(in) :: H_LOCAL,K_LOCAL

OCCUPATION(HYDROGEN_X(H_LOCAL),HYDROGEN_Y(H_LOCAL))=0

if (K_LOCAL .eq. 1) then
    HYDROGEN_X(H_LOCAL)=mod1(HYDROGEN_X(H_LOCAL)+1)
else if (K_LOCAL .eq. 2) then
    HYDROGEN_X(H_LOCAL)=mod1(HYDROGEN_X(H_LOCAL)-1)
else if (K_LOCAL .eq. 3) then
    HYDROGEN_Y(H_LOCAL)=mod1(HYDROGEN_Y(H_LOCAL)+1)
else if (K_LOCAL .eq. 4) then
    HYDROGEN_Y(H_LOCAL)=mod1(HYDROGEN_Y(H_LOCAL)-1)
end if

OCCUPATION(HYDROGEN_X(H_LOCAL),HYDROGEN_Y(H_LOCAL))=1

end subroutine


subroutine HYDROGEN_CHECK_HYDROGEN(H_LOCAL,K_LOCAL,HKH,HYDROGEN_TOUCH_HYDROGEN)
![in] hydrogen series number
![out] OCCUPATION AND
![do] check number if there is any Hydrogen contact with H_LOCAL
integer,intent(in) :: H_LOCAL,K_LOCAL
integer,intent(out) :: HKH
logical,intent(out) :: HYDROGEN_TOUCH_HYDROGEN

!LOCALVARIABLE
integer :: H_LOOP_LOCAL


HYDROGEN_TOUCH_HYDROGEN=.false.

if(occ(HYDROGEN_X(H_LOCAL)+K_ARRAY(K_LOCAL,1), &
    HYDROGEN_Y(H_LOCAL)+K_ARRAY(K_LOCAL,2)) .eq. 1) then
    do H_LOOP_LOCAL=1,N_of_hydrogen!_created
        if(HYDROGEN_EXISTED(H_LOOP_LOCAL))then
            !check if there is another hydrogen contact with it
            if(HYDROGEN_X(H_LOOP_LOCAL) .eq. &
            HYDROGEN_X(H_LOCAL)+K_ARRAY(K_LOCAL,1))then
                if(HYDROGEN_Y(H_LOOP_LOCAL) .eq. &
                HYDROGEN_Y(H_LOCAL)+K_ARRAY(K_LOCAL,2)) then
                    HYDROGEN_TOUCH_HYDROGEN=.true.
                    HKH=H_LOOP_LOCAL
                end if
            end if
        end if
    end do
end if

end subroutine


subroutine check_h_hydrated_p_corner(H_LOCAL,K_LOCAL,HKP,HKPK,h_deP_touched_LOCAL)
integer,intent(in) :: K_LOCAL,H_LOCAL
integer,intent(out) :: HKP,HKPK
logical,intent(out) :: h_deP_touched_LOCAL
integer :: i_local,j_local,X_LOCAL,Y_LOCAL,H_X_LOCAL,H_Y_LOCAL
HKP=-1
h_deP_touched_LOCAL=.false.
HKPK=-1
do  i_local=1,N_of_particles
    do j_local=5,8
        X_LOCAL=CORNERS(X(i_local),j_local,1)
        Y_LOCAL=CORNERS(Y(i_local),j_local,2)
        H_X_LOCAL=mod1(HYDROGEN_X(H_LOCAL)+K_ARRAY(K_LOCAL,1))
        H_Y_LOCAL=mod1(HYDROGEN_Y(H_LOCAL)+K_ARRAY(K_LOCAL,2))
        if(X_LOCAL .eq. H_X_LOCAL)then
            if(Y_LOCAL .eq. H_Y_LOCAL) then
                h_deP_touched_LOCAL=.true.
                HKP=i_local
                HKPK=j_local
                return
            end if
        end if
    end do
end do


end subroutine


subroutine hydrogen_gone(H_LOCAL,HKH_LOCAL)
integer,intent(in) :: H_LOCAL,HKH_LOCAL
if (HKH_LOCAL .ne. 0)then

    HYDROGEN_EXISTED(H_LOCAL)=.false.
    OCCUPATION(HYDROGEN_X(H_LOCAL),HYDROGEN_Y(H_LOCAL))=0
    HYDROGEN_X(H_LOCAL)=-1
    HYDROGEN_Y(H_LOCAL)=-1
    
    HYDROGEN_EXISTED(HKH_LOCAL)=.false.
    OCCUPATION(HYDROGEN_X(HKH_LOCAL),HYDROGEN_Y(HKH_LOCAL))=0
    HYDROGEN_X(HKH_LOCAL)=-1
    HYDROGEN_Y(HKH_LOCAL)=-1
else if(HKH_LOCAL .eq. 0)then

    HYDROGEN_EXISTED(H_LOCAL)=.false.
    OCCUPATION(HYDROGEN_X(H_LOCAL),HYDROGEN_Y(H_LOCAL))=0
    HYDROGEN_X(H_LOCAL)=-1
    HYDROGEN_Y(H_LOCAL)=-1
end if

end subroutine

real function ran1(idum)
INTEGER IA,IM,IQ,IR,NTAB,NDIV,idum
REAL AM,EPS,RNMX
PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
INTEGER j,k,iv(NTAB),iy
SAVE iv,iy
DATA iv /NTAB*0/, iy /0/
if (idum.le.0.or.iy.eq.0)then
    idum=max(-idum,1)
    do j=NTAB+8,1,-1
        k=idum/IQ
        idum=IA*(idum-k*IQ)-IR*k
        if (idum.lt.0) idum=idum+IM
        if (j.le.NTAB) iv(j)=idum
    enddo
    iy=iv(1)
endif
k=idum/IQ
idum=IA*(idum-k*IQ)-IR*k
if (idum.lt.0) idum=idum+IM
j=1+iy/NDIV
iy=iv(j)
iv(j)=idum
ran1=min(AM*iy,RNMX)
return
END FUNCTION ran1


end module




