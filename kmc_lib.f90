module kmc_lib

  use common_parameters !initialise
  use lib
  use events
  
  implicit none

  contains
  
  
subroutine choose_kmc_events(MAX_SEQ_LOCAL,R_c_local,SEQUENCE_OF_EVENTS_LOCAL)

  real(r_p),intent(in) :: R_c_local !Cumulative transition rate
  integer(i_p),intent(in) :: MAX_SEQ_LOCAL
  integer(i_p),intent(out) :: SEQUENCE_OF_EVENTS_LOCAL !what to happen

  if(random_or_manual)then
     call randomly_choose_a_sequence(SEQUENCE_OF_EVENTS_LOCAL,R_c_local)!sequence_events as output
  else if(.not. random_or_manual)then
     call manual_seq(MAX_SEQ_LOCAL,KMC_STEPS,SEQUENCE_OF_EVENTS_LOCAL)
  end if
  
  delta_time=delta_time-log(ran1(idum))/R_c_local
  
end subroutine

subroutine KMC_EVENTS(KMC_STEPS_LOCAL,R_c_LOCAL,MAX_SEQ_LOCAL)

!GLOBAL VARIABLE
integer(i_p),intent(inout) :: KMC_STEPS_LOCAL !LOOP NUMBER
real(r_p),intent(out) :: R_c_LOCAL
integer(i_p),intent(out) :: MAX_SEQ_LOCAL


!LOCAL VARIABLE
logical dehydratable,h_deP_touched
logical,allocatable,dimension(:,:) :: DOUBEL_COUNTED,DOUBEL_H_H_COUNTED
integer(i_p) :: sequence_of_events_local
integer(i_p) :: P,H,K,HKH,HKP,HKPK,PKP
logical :: HYDROGEN_TOUCH_HYDROGEN,corner_corner_connectable,PK_OCCUPIED


allocate(DOUBEL_COUNTED(1:N_of_particles,1:N_of_particles)); DOUBEL_COUNTED=.false.
allocate(DOUBEL_H_H_COUNTED(1:N_of_hydrogen,1:N_of_hydrogen)); DOUBEL_H_H_COUNTED=.false.

KMC_STEPS_LOCAL=KMC_STEPS_LOCAL+1
if(endless_kmc_step) MAX_MAX_KMC_STEPS=MAX_MAX_KMC_STEPS+1

   LOOKUP=0;LOOKUP_PKP=0;LOOKUP_HKH=0
   sequence_of_events_local=0;R_c_LOCAL=0.0_r_p;R=0.0_r_p
   
!================random walk hydrogen========== !
   if(H_random) then
      do H = 1,N_of_hydrogen
         if(HYDROGEN_EXISTED(H))then
         OCCUPATION(HYDROGEN_X(H),HYDROGEN_Y(H))=0
            do
               Hydrogen_X(H)=nint(N*ran1(idum))
               Hydrogen_Y(H)=nint(N*ran1(idum))
               if(OCC(HYDROGEN_X(H),Hydrogen_Y(H)) .eq. 0)then
                  OCCUPATION(HYDROGEN_X(H),HYDROGEN_Y(H))=1
                  exit
               end if
            end do
         end if
      end do
   end if
!===============end random walk================ !
   
   do P=1,N_of_particles
!=====================creat events for free particles=================================
      if (DIMER_P(P) .eq. 0) then !not connected
         !............... diffusion of a free molecule
         do k=1,4
            call CHECK_OCCUPANCY(P,K,PK_OCCUPIED)
            if (.not. PK_OCCUPIED) then
               sequence_of_events_local=sequence_of_events_local+1
               LOOKUP(sequence_of_events_local,3) = 1
               if(K .eq. 1 .or. K .eq. 2)then
                  R(sequence_of_events_local)=TST(1)
               else if(K .eq. 3 .or. K.eq. 4)then
                  R(sequence_of_events_local)=TST(0)
               end if
               R_c_LOCAL=R_c_LOCAL+R(sequence_of_events_local)
               LOOKUP(sequence_of_events_local,1)=P
               LOOKUP(sequence_of_events_local,2)=K
            end if
         end do
         !............... corner events
         do k=5,8 !free particle
            !check dehydration events_type 3
            if(.not. dehydrated(P,K))then
               !......... add dehydration (if possible) for this corner K
               dehydratable=.false.
               call check_dehyd_or_creat_h(P,K,dehydratable,0)
               if(dehydratable)then
                  sequence_of_events_local=sequence_of_events_local+1
                  LOOKUP(sequence_of_events_local,3) = 3
                  R(sequence_of_events_local)=TST(3)
                  R_c_LOCAL=R_c_LOCAL+R(sequence_of_events_local)
                  LOOKUP(sequence_of_events_local,1)=P
                  LOOKUP(sequence_of_events_local,2)=K
               end if
            else if(dehydrated(P,K))then
               !........ add association at this corner K between two de-H molecules
               call CHECK_OCCUPANCY(P,K,PK_OCCUPIED) !???.... not clear if it could also be occ by H?
               if (PK_OCCUPIED) then
                  call CHECK_CORNER(P,K,PKP,corner_corner_connectable)
                  !this will gives PKP
                  if(corner_corner_connectable) then
                     if(dehydrated(PKP,K_opposite(K)))then
                        if(.not. DOUBEL_COUNTED(P,PKP)) then
!======================creat events for association====================================
                           sequence_of_events_local=sequence_of_events_local+1
                           LOOKUP(sequence_of_events_local,3) = 2 !ASSOCIATION
                           R(sequence_of_events_local)=TST(2)
                           R_c_LOCAL=R_c_LOCAL+R(sequence_of_events_local)
                           LOOKUP_PKP(sequence_of_events_local)=PKP
                           LOOKUP(sequence_of_events_local,1)=P
                           LOOKUP(sequence_of_events_local,2)=K
                           DOUBEL_COUNTED(P,PKP)=.true.
                           DOUBEL_COUNTED(PKP,P)=.true.
                        end if!double counting
                     end if !dehydrated
                  end if !corner_corner_connectable
               end if !PK OCCUPIED
            end if !dehydrated
         end do !direction from 5 to 8

!===================================creat events for frozen particle======================
      else if (DIMER_P(P) .gt. 0) then !consider frozen particles
         do k=5,8
            !two kind of events dehydration and association
            if(.not. dehydrated(P,K)) then
               !........ add dehydration (if possible) at this corner K
               dehydratable=.false.
               call check_dehyd_or_creat_h(P,K,dehydratable,0)
               if(dehydratable)then
                  sequence_of_events_local=sequence_of_events_local+1
                  LOOKUP(sequence_of_events_local,3) = 3
                  R(sequence_of_events_local)=TST(3)
                  R_c_LOCAL=R_c_LOCAL+R(sequence_of_events_local)
                  LOOKUP(sequence_of_events_local,1)=P
                  LOOKUP(sequence_of_events_local,2)=K
               end if
            ELSE if(dehydrated(P,K)) then
               !......... add association betwween two de-H molecules at this corner K
               call CHECK_OCCUPANCY(P,K,PK_OCCUPIED)
               if (PK_OCCUPIED) then
                  call CHECK_CORNER(P,K,PKP,corner_corner_connectable)
                  if(.not. corner_corner_connectable) then!means contact with a side of another particle
                     cycle
                  else if(corner_corner_connectable) then!prevent connected already, so check if it is a dimer
                     if(dehydrated(PKP,K_opposite(K)))then
                        if(.not. DOUBEL_COUNTED(P,PKP)) then
                           if(bonds(P,K) .eq. PKP) then
                              cycle
                           else if(bonds(P,K) .eq. 0)then !frozen connected with unbonded particle
!======================creat events for association====================================
                              sequence_of_events_local=sequence_of_events_local+1
                              LOOKUP(sequence_of_events_local,3) = 2
                              R(sequence_of_events_local)=TST(2)
                              R_c_LOCAL=R_c_LOCAL+R(sequence_of_events_local)
                              LOOKUP(sequence_of_events_local,1)=P
                              LOOKUP(sequence_of_events_local,2)=K
                              LOOKUP_PKP(sequence_of_events_local) = PKP
                              DOUBEL_COUNTED(P,PKP)=.true.
                              DOUBEL_COUNTED(PKP,P)=.true.
                           end if!bonds
                        end if !DOUBLE_COUNTING
                     end if!dehydrated
                  end if!corner_corner_connectable
               end if !PK_OCCUPIED
!           end if!!dehydratable
            end if!dehydrated
         end do !K
      end if !DIMER_P
   end do !P

   

!============================HYDROGEN EVENTS================================
!  if(N_of_hydrogen!_created .ne. 0)then
      do H = 1,N_of_hydrogen
!===================MOVEMENT================
! first check if this atom disappered
         if(HYDROGEN_EXISTED(H))then
            do K = 1,4
               !....... add diffusion of H atom
               if(OCC(HYDROGEN_X(H)+K_ARRAY(K,1), &
                    HYDROGEN_Y(H)+K_ARRAY(K,2)) .eq. 0) then
                  if(.not.H_random) then
                     sequence_of_events_local=sequence_of_events_local+1
                     LOOKUP(sequence_of_events_local,3) = 4
                     R(sequence_of_events_local)=TST(4)
                     R_c_LOCAL=R_c_LOCAL+R(sequence_of_events_local)
                     LOOKUP(sequence_of_events_local,1)=H
                     LOOKUP(sequence_of_events_local,2)=K
                  end if
               else if(OCC(HYDROGEN_X(H)+K_ARRAY(K,1), &
                    HYDROGEN_Y(H)+K_ARRAY(K,2)) .ne. 0) then
                  !hydrogen leave cu
                  !...... add formation of H2 and desorption
                  call HYDROGEN_CHECK_HYDROGEN(H,K,HKH,HYDROGEN_TOUCH_HYDROGEN)
                  if(HYDROGEN_TOUCH_HYDROGEN)then
                     if(.not. DOUBEL_H_H_COUNTED(H,HKH))then
                        sequence_of_events_local=sequence_of_events_local+1
                        R(sequence_of_events_local)=TST(5)
                        R_c_LOCAL=R_c_LOCAL+R(sequence_of_events_local)
                        LOOKUP(sequence_of_events_local,3) = 5
                        LOOKUP(sequence_of_events_local,1) = H
                        LOOKUP(sequence_of_events_local,2) = K
                        LOOKUP_HKH(sequence_of_events_local,1) = HKH
                        DOUBEL_H_H_COUNTED(H,HKH)=.true.
                        DOUBEL_H_H_COUNTED(HKH,H)=.true.
                     end if
                  else if(.not. HYDROGEN_TOUCH_HYDROGEN) then !it must touched a TMTPP
                     !if it touched a corner that is dehydrated
                     !...... add hydration back to the de-H molecule
                     h_deP_touched=.false.
                     call check_h_hydrated_p_corner(H,K,HKP,HKPK,h_deP_touched) !HKP
                     !hydrate
                     if(h_deP_touched)then
                        if(bonds(HKP,HKPK) .eq. 0)then
                           if(dehydrated(HKP,HKPK))then
                              sequence_of_events_local=sequence_of_events_local+1
                              R(sequence_of_events_local)=TST(6)
                              R_c_LOCAL=R_c_LOCAL+R(sequence_of_events_local)
                              LOOKUP(sequence_of_events_local,3) = 6
                              LOOKUP(sequence_of_events_local,1) = H
                              LOOKUP(sequence_of_events_local,2) = K
                              LOOKUP_HKH(sequence_of_events_local,1) = HKP
                              LOOKUP_HKH(sequence_of_events_local,2) = HKPK
                           end if!dehydrated
                        end if!bonds
                     end if!touched
                  end if !h touch h
               end if !if OCC
            end do !K=1,4
         end if !HYDROGEN_EXISTED
      end do !H = 1,N_of_hydrogen!_created
!   end if !N_of_hydrogen .ne. 0

   MAX_SEQ_LOCAL=sequence_of_events_local

end subroutine

subroutine update_kmc_events(SEQUENCE_OF_EVENTS_LOCAL)
integer(i_p),intent(out) :: SEQUENCE_OF_EVENTS_LOCAL

!LOCAL VARIABLE
logical :: dehydratable
integer(i_p) :: K !direction from 1 to 8
integer(i_p) :: P,H,HKH,HKPK,HKP,PKP

   if(LOOKUP(sequence_of_events_local,3) .ge. 1 .and. LOOKUP(sequence_of_events_local,3) .le. 3) then
      P=LOOKUP(sequence_of_events_local,1)
      K=LOOKUP(sequence_of_events_local,2)
    !===================free movement================
      if(LOOKUP(sequence_of_events_local,3) .eq. 1) then !free movement
         if(K .eq. 1 .or. K .eq. 2)then
            counter_type(1)=counter_type(1)+1
         else if(K .eq. 3 .or. K.eq. 4)then
            counter_type(0)=counter_type(0)+1
         end if
         counter_K(K)=counter_K(K)+1
         call UPDATE_OCCUPATION(P,K)
    !===================association==================
      else if(LOOKUP(sequence_of_events_local,3) .eq. 2)then !association
         PKP=LOOKUP_PKP(sequence_of_events_local)
         DIMER_P(P)=DIMER_P(P)+1
         DIMER_P(PKP)=DIMER_P(PKP)+1
         bonds(P,K)=PKP
         bonds(PKP,K_opposite(K))=P
         counter_type(2)=counter_type(2)+1
      else if(LOOKUP(sequence_of_events_local,3) .eq. 3)then !dehydration
         dehydrated(P,K)=.true.
         counter_type(3)=counter_type(3)+1
!======================CREAT HYDROGEN ATOM===================================
         call check_dehyd_or_creat_h(P,K,dehydratable,1)
!======================END CREAT HYDROGEN ATOM===============================
      end if

   else if(LOOKUP(sequence_of_events_local,3) .eq. 4) then !HYDROGEN MOVEMENT
      H = LOOKUP(sequence_of_events_local,1)
      K = LOOKUP(sequence_of_events_local,2)
      counter_type(4)=counter_type(4)+1
      counter_hydrogen_move=counter_hydrogen_move+1
      counter_not_hydrogen_move=0
      call HYDROGEN_OCCUPATION_UPDATE(H,K)!move HYDROGEN IN THE OCCUPATION ARRAY
   else if(LOOKUP(sequence_of_events_local,3) .eq. 5) then !HYDROGEN GONE
      H = LOOKUP(sequence_of_events_local,1)
      HKH = LOOKUP_HKH(sequence_of_events_local,1)
      counter_type(5)=counter_type(5)+1
      call hydrogen_gone(H,HKH)
   else if(LOOKUP(sequence_of_events_local,3) .eq. 6)then !hydrate
      H = LOOKUP(sequence_of_events_local,1)
      K = LOOKUP(sequence_of_events_local,2)
      HKP = LOOKUP_HKH(sequence_of_events_local,1)
      HKPK = LOOKUP_HKH(sequence_of_events_local,2)
      call hydrogen_gone(H,0)!only one hydrogen gone
      dehydrated(HKP,HKPK)=.false.
      counter_type(6)=counter_type(6)+1
   end if
   
   if(record_H) then
      if(LOOKUP(sequence_of_events_local,3) .ne. 4)then
         counter_not_hydrogen_move=counter_not_hydrogen_move+1
         call FRQ_of_HYDROGEN(delta_time,counter_hydrogen_move,counter_not_hydrogen_move)
         counter_hydrogen_move=0
      end if
   end if
end subroutine


end module
