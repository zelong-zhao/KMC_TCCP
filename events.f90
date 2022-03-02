module events
    use lib
    use common_parameters
    implicit none
    
contains

subroutine how_to_print(SEQUENCE_OF_EVENTS_LOCAL,force_to_print)

  integer(i_p),intent(in) :: SEQUENCE_OF_EVENTS_LOCAL
  logical,intent(in)::force_to_print
  logical :: yes_to_print
  yes_to_print=.false.
  if(force_to_print)yes_to_print=.true.
  if(random_or_manual)then
  !==============print 100 data===============
    if(LOOKUP(SEQUENCE_OF_EVENTS_LOCAL,3) .eq. 5 .or. LOOKUP(SEQUENCE_OF_EVENTS_LOCAL,3) .eq. 2) then
      yes_to_print=.true.
    else if (KMC_STEPS .eq. MAX_MAX_KMC_STEPS) then
      yes_to_print=.true.
    else if(type_of_simulation .eq. 'TPD' .and. temp .ge. temp_plot) then
      yes_to_print=.true.
      write(22,'(3(e12.6,X),(i10))') delta_time,temp,P_STAY,counter_type(5)
      temp_plot=temp+temp_per_plot
    end if
    
      if(yes_to_print)then
        call discontinue_data()
        write(20,'(e12.6,8(x,i10))') delta_time,KMC_STEPS,counter_type(0:6)
        write(90,'(e12.6,4(x,i10))')delta_time,KMC_STEPS,N_of_hydrogen_on_surface(),counter_type(5),N_of_hydrogen_in_molecules()
        write(89,'(e12.6,6(x,i10))')delta_time,KMC_STEPS,N_of_particles_with_different_deH(0),N_of_particles_with_different_deH(1),N_of_particles_with_different_deH(2),N_of_particles_with_different_deH(3),N_of_particles_with_different_deH(4)
        
        call save_frames(KMC_STEPS)
        if(LOOKUP(SEQUENCE_OF_EVENTS_LOCAL,3) .eq. 2) then
          call plot()
          call Numb_clusters(N_clusters,Compactness_index,.false.)
          write(77,'(e12.6,2(5x,i0),5x,f8.3)')delta_time,KMC_STEPS,N_clusters,Compactness_index
        end if
        call print_()
      end if
    !==============end print====================
  end if
   
end subroutine


subroutine open_kmc_files()
integer(i_p)::i
if(state .eq. 'N') then
   OPEN(UNIT=4,FILE='CONECTIVITY.DATA')
   OPEN(20,file='stat_moves.DATA')
   write(4,'(7a)')'# KMC_STEPS ','delta_time ','zero ','one ','two ',&
   'three ','four '
   write(20,'(9a)')'#  TIME    ','    KMC_STEP  ',(' '//counter_type_name(i)//'  ',i=0,6)
   OPEN(77,file='COMPACTNESS_INDEX.DATA')
   write(77,'(4a)')'#  TIME    ','    KMC_STEP  ',' N_clusters ',' COMPACTNESS_INDEX  '
   OPEN(90,file='COMPLETE_HYDROGEN.DATA')
   write(90,'(5a)')'#  TIME    ','    KMC_STEP   ',' FREE ELECTRON ','    H_GONE    ','   H_BONDED   '
   OPEN(89,file='CompleteDeH.DATA')
   write(89,'(7a)')'#  TIME    ','    KMC_STEP   ','     Zero    ','     One    ','    Two    ','   Three   ','   four   '
   if(record_H) then
     open(UNIT=10,FILE='FEQ_H.DATA')
     write(10,'(3a)')'#  TIME    ','    H_MOVE    ','   H_NOT_MOVE   '
   end if
else if(state .eq. 'C')then
   OPEN(UNIT=4,FILE='CONECTIVITY.DATA',STATUS='OLD',ACCESS='APPEND')
   OPEN(20,file='stat_moves.DATA',STATUS='OLD',ACCESS='APPEND')
   OPEN(77,file='COMPACTNESS_INDEX.DATA',STATUS='OLD',ACCESS='APPEND')
   OPEN(90,file='COMPLETE_HYDROGEN.DATA',STATUS='OLD',ACCESS='APPEND')
   OPEN(89,file='CompleteDeH.DATA',STATUS='OLD',ACCESS='APPEND')
   if(record_H) open(UNIT=10,FILE='FEQ_H.DATA',STATUS='OLD',ACCESS='APPEND')
end if


end subroutine

subroutine print_()
![in] bonds
!global varible
character,allocatable,dimension(:) :: character_print*3
!local variable
integer(i_p) :: i_local,j_local,k_local
integer(i_p) :: P_LOCAL,H_LOCAL,p_loc

if(.not. print_map_in_terminal .or. N.gt.99) return


allocate(character_print(0:N+1)); character_print='-'

!=================================================!
do i_local=1,N_of_particles
    do k_local=5,8
        if(bonds(i_local,k_local) .ne. 0) then
            CALL UPDATE_PARTICLE_SHAPE(i_local,2)
        end if
    end do
end do

do i_local=1,N_of_hydrogen
    if(HYDROGEN_EXISTED(i_local))then
        OCCUPATION(HYDROGEN_X(i_local),HYDROGEN_Y(i_local))=3
    end if
end do

do i_local=1,N_of_particles
    do k_local=5,8
        if(dehydrated(i_local,k_local))then
            OCCUPATION(CORNERS(X(i_local),k_local,1),CORNERS(Y(i_local),k_local,2))=4
        end if
    end do
end do
!=================================================!
!..... printing the map

character_print=' - '

do i_local=0,N
  !character_print(i+1)=char(i_local)
  if (i_local .le. 9) then
    write(character_print(i_local+1),'(i1,2X)')i_local
  else if(i_local .le. 99) then
    write(character_print(i_local+1),'(i2,X)')i_local
  else if(i_local .le. 999) then
    write(character_print(i_local+1),'(i3)')i_local
  else
    stop 'N too large'
  end if
end do
write(*,'('//len_grid//'(A3))')character_print(0:N+1)

do  i_local=0,N
   if (i_local .le. 9) then
      write(character_print(0),'(i1,2X)')i_local
   else if(i_local .le. 99) then
      write(character_print(0),'(i2,X)')i_local
   end if
   do j_local=0,N
      if(OCCUPATION(i_local,j_local) .eq. 0)then
         character_print(j_local+1)='-  '
      else if(OCCUPATION(i_local,j_local) .eq. 1)then
         character_print(j_local+1)='1  '
      else if(OCCUPATION(i_local,j_local) .eq. 2)then
         character_print(j_local+1)='C  '
      else if(OCCUPATION(i_local,j_local) .eq. 3)then
         character_print(j_local+1)='H  '
      else if(OCCUPATION(i_local,j_local) .eq. 4)then
         character_print(j_local+1)='D  '
      end if
   end do
!....... put the paeticvle number in the centre of the square of its image
   do p_loc=1,N_of_particles
      do j_local=0,N
         if(j_local.eq.mod1(Y(p_loc)+1) .and. i_local.eq.mod1(X(p_loc)+1)) then
            write(character_print(j_local+1),'(i1,2x)') p_loc
         end if
      end do
   end do
   write(*,'('//len_grid//'(A3))')character_print(0:N+1)
end do
!...... finished printing the map



do i_local=1,N_of_hydrogen!_created
   if(HYDROGEN_EXISTED(i_local))then
      write(*,*)'hydrogen at',HYDROGEN_X(i_local),HYDROGEN_Y(i_local)
   end if
end do

!=================================================!
do P_LOCAL=1,N_of_particles
   do k_local=5,8
      if(bonds(P_LOCAL,k_local) .ne. 0) then
         CALL UPDATE_PARTICLE_SHAPE(P_LOCAL,1)
      end if
   end do
end do

do H_LOCAL=1,N_of_hydrogen!_created
   if(HYDROGEN_EXISTED(H_LOCAL))then
      OCCUPATION(HYDROGEN_X(H_LOCAL),HYDROGEN_Y(H_LOCAL))=1
   end if
end do

do i_local=1,N_of_particles
   do k_local=5,8
      if(dehydrated(i_local,k_local))then
         OCCUPATION(CORNERS(X(i_local),k_local,1),CORNERS(Y(i_local),k_local,2))=1
      end if
   end do
end do
!==================================================!
deallocate(character_print)
write(*,*)'KMC',KMC_STEPS
write(*,*)'<------------------------------------------------->'
end subroutine print_

subroutine manual_seq(LOCAL_SEQUENCE_OF_EVENTS,KMC_STEPS_LOCAL,out_seq)
integer,intent(in) :: local_SEQUENCE_OF_EVENTS,KMC_STEPS_LOCAL
integer,intent(out) :: out_seq

integer :: i_local,j_local,read_seq
character :: CHAR_LOCAL*8,dir(8)*2
integer :: counter_type_local(1:6)
dir(1)='R ';dir(2)='L ';dir(3)='U ';dir(4)='D '
dir(5)='LD';dir(6)='LU';dir(7)='RU';dir(8)='RD'
counter_type_local=0
write(*,'(A,I0)')'KMC_STEPS:',KMC_STEPS_LOCAL
write(*,'("total step= ",I0)')LOCAL_SEQUENCE_OF_EVENTS
j_local=LOCAL_SEQUENCE_OF_EVENTS

do i_local=1,j_local
    CHAR_LOCAL=counter_type_name(LOOKUP(i_local,3))
!    write(*,'(i0," sequence_of_events")')i_local
    if(LOOKUP(i_local,3) .eq. 2) then
        counter_type_local(LOOKUP(i_local,3))=counter_type_local(LOOKUP(i_local,3))+1
!1010        format(i3,x,2(A,I0),X,2(A,X),2(A,I0),A1)
!1010        format(i3,x,A,I0,a,a,X,2(A,X),2(A,I0),A1)
!        write(*,1010) i_local,'P=',LOOKUP(i_local,1),' K=',LOOKUP(i_local,2),' type=',CHAR_LOCAL,&
        write(*,'(i2,a,i2,5a,2(i2,a),6(a,i2))') &
             i_local,' P=',LOOKUP(i_local,1),' dir=',dir(LOOKUP(i_local,2)),' type=',CHAR_LOCAL,&
             ' (X,Y)=(',X(LOOKUP(i_local,1)),',',Y(LOOKUP(i_local,1)),')',&
             ' PKP=',LOOKUP_PKP(i_local),' K=',K_opposite(LOOKUP(i_local,2)),&
             ' DIMER_P(P)=',DIMER_P(LOOKUP(i_local,1)),' DIMER_P(PKP)=',DIMER_P(LOOKUP_PKP(i_local)),&
             ' BONDS(P,K)=',BONDS(LOOKUP(i_local,1),LOOKUP(i_local,2)),&
             ' BONDS(PKP,K)=',BONDS(LOOKUP_PKP(i_local),K_opposite(LOOKUP(i_local,2)))
    else if(LOOKUP(i_local,3) .eq. 1 .or. LOOKUP(i_local,3) .eq. 3)then
        counter_type_local(LOOKUP(i_local,3))=counter_type_local(LOOKUP(i_local,3))+1
        write(*,'(i2,a,i2,5a,i2,a,i2,a)') i_local,' P=',LOOKUP(i_local,1),' dir=',dir(LOOKUP(i_local,2)),' type=',CHAR_LOCAL,&
        ' (X,Y)=(',X(LOOKUP(i_local,1)),',',Y(LOOKUP(i_local,1)),')'
    else
        counter_type_local(LOOKUP(i_local,3))=counter_type_local(LOOKUP(i_local,3))+1
        if(LOOKUP(i_local,3) .ne. 6)then
           write(*,'(i2,a,i2,4a,2(a,i2),a)')i_local,' H=',LOOKUP(i_local,1),' dir=',dir(LOOKUP(i_local,2)),' type=',CHAR_LOCAL,&
                '(X,Y)=(',HYDROGEN_X(LOOKUP(i_local,1)),',',HYDROGEN_Y(LOOKUP(i_local,1)),')'
        else
           write(*,'(i2,a,i2,4a,2(a,i2),2a,i2)')i_local,' H=',LOOKUP(i_local,1),' dir=',dir(LOOKUP(i_local,2)),' type=',CHAR_LOCAL,&
                ' (X,Y)=(',HYDROGEN_X(LOOKUP(i_local,1)),',',HYDROGEN_Y(LOOKUP(i_local,1)),')',&
                'hydrate =',LOOKUP_HKH(i_local,1)
        end if
    end if
end do

do i_local=1,6
    write(*,'(A8,": ",i0,X)')counter_type_name(i_local),counter_type_local(i_local)
end do
call print_()

12 write(*,*)'CHOOSE AN EVENT_________> '

read(*,*)read_seq
!call randomly_choose_a_sequence(read_seq)
!write(*,*)read_seq
!read(*,*)

if(read_seq .gt. j_local .or. read_seq .lt. 1) go to 12
    write(*,'(A,X,I0,X,I0,X,A)')'AT KMC STEPS', read_seq,LOOKUP(read_seq,3),'HAPPENED'
out_seq=read_seq
end subroutine

subroutine FRQ_of_HYDROGEN(delta_time_local,counter_hydrogen_move_local,counter_not_hydrogen_move_local)
integer,intent(in) :: counter_hydrogen_move_local,counter_not_hydrogen_move_local
real(r_p),intent(in) :: delta_time_local
    write(10,'(e12.6,2(x,i10))') delta_time_local,counter_hydrogen_move_local,&
    counter_not_hydrogen_move_local
end subroutine

subroutine plot()

integer :: TEMP_P
integer :: TEMP_UNCONNECTED,TEMP_CONNECTED_ONE,TEMP_CONNECTED_TWO,TEMP_CONNECTED_THREE,TEMP_CONNECTED_FOUR

TEMP_UNCONNECTED=0
TEMP_CONNECTED_ONE=0
TEMP_CONNECTED_TWO=0
TEMP_CONNECTED_THREE=0
TEMP_CONNECTED_FOUR=0
do TEMP_P=1,N_of_particles
   if(DIMER_P(TEMP_P) .eq. 0) then !unconnected particles
      TEMP_UNCONNECTED=TEMP_UNCONNECTED+1
   else if (DIMER_P(TEMP_P) .eq. 1) then
      TEMP_CONNECTED_ONE=TEMP_CONNECTED_ONE+1
   else if (DIMER_P(TEMP_P) .eq. 2) then
      TEMP_CONNECTED_TWO=TEMP_CONNECTED_TWO+1
   else if (DIMER_P(TEMP_P) .eq. 3) then
      TEMP_CONNECTED_THREE=TEMP_CONNECTED_THREE+1
   else if (DIMER_P(TEMP_P) .eq. 4) then
      TEMP_CONNECTED_FOUR=TEMP_CONNECTED_FOUR+1
   end if
end do
120 format((I0,A1)(e10.5,A1)4(I0,A1)(I0)) !CSV
write(4,120)KMC_STEPS,','&
,delta_time,',',TEMP_UNCONNECTED,',',TEMP_CONNECTED_ONE,',',TEMP_CONNECTED_TWO,','&
,TEMP_CONNECTED_THREE,',',TEMP_CONNECTED_FOUR

end subroutine

subroutine discontinue_data()
    integer(i_p) :: i_local
    
    OPEN(UNIT=7,FILE='DATA_LAST.dat')
    write(7,'(3(A,","))')'# KMC_STEPS','delta_time','temp'
    write(7,'((i0,",")(e10.5,",")(f9.5,","))')KMC_STEPS,delta_time,temp
    write(7,'(("#",X)7(A,","))')counter_type_name(0:6)
    write(7,'(7(i0,","))')counter_type(0:6)
    write(7,'(("#",X)11(A,","))')'X','Y','DIMER_P','bonds_5','bonds_6','bonds_7','bonds_8', &
    'dehydrated_5','dehydrated_6','dehydrated_7','dehydrated_8'
    do  i_local=1,N_of_particles
        write(7,'(7(i0,",")4(L1,","))')X(i_local),Y(i_local),DIMER_P(i_local),bonds(i_local,5:8), &
            dehydrated(i_local,5:8)
    end do
    write(7,'(("#",X),(A,","))')'N_of_hydrogen_on_surface'
    write(7,'(i0,",")')N_of_hydrogen_on_surface()
    write(7,'(("#",X),3(A,","))')'HYDROGEN_X','HYDROGEN_Y','HYDROGEN_EXISTED'
    do i_local=1,N_of_hydrogen
        if(HYDROGEN_EXISTED(i_local))then
        write(7,'(2(i0,",")(L1,","))')HYDROGEN_X(i_local),HYDROGEN_Y(i_local),HYDROGEN_EXISTED(i_local)
        end if
    end do
    write(7,'("#",X,A)')'random_seed'
    write(7,'(i0)') idum
    close(7)

end subroutine discontinue_data

subroutine save_frames(frame)
  integer(i_p),intent(in)::frame
  character :: this_used_to_print*100
  integer(i_p) :: i_local

  write(this_used_to_print,'(I0,".frame")')frame
  
  OPEN(UNIT=88,file=this_used_to_print)
  write(88,'(3(A,","))')'# KMC_STEPS','delta_time','temp'
  write(88,'((i0,",")(e10.5,",")(f9.5,","))')KMC_STEPS,delta_time,temp
  write(88,'(("#",X)7(A,","))')counter_type_name(0:6)
  write(88,'(7(i0,","))')counter_type(0:6)
  write(88,'(("#",X)11(A,","))')'X','Y','DIMER_P','bonds_5','bonds_6','bonds_7','bonds_8', &
  'dehydrated_5','dehydrated_6','dehydrated_7','dehydrated_8'
  do  i_local=1,N_of_particles
      write(88,'(7(i0,",")4(L1,","))')X(i_local),Y(i_local),DIMER_P(i_local),bonds(i_local,5:8), &
          dehydrated(i_local,5:8)
  end do
  write(88,'(("#",X),(A,","))')'N_of_hydrogen_on_surface'
  write(88,'(i0,",")')N_of_hydrogen_on_surface()
  write(88,'(("#",X),3(A,","))')'HYDROGEN_X','HYDROGEN_Y','HYDROGEN_EXISTED'
  do i_local=1,N_of_hydrogen
      if(HYDROGEN_EXISTED(i_local))then
      write(88,'(2(i0,",")(L1,","))')HYDROGEN_X(i_local),HYDROGEN_Y(i_local),HYDROGEN_EXISTED(i_local)
      end if
  end do
  write(88,'("#",X,A)')'random_seed'
  write(88,'(i0)') idum
  close(88)
  
end subroutine

subroutine export_to_povray()
  integer(i_p) :: p_loc,k_local,x_local,y_local,bx,by
  open(UNIT=99,file='POV.DATA')
358 format(2(I0,","))
  write(99,358)N*100,N*100
  do p_loc=1,N_of_particles
    write(99,358)X(p_loc),Y(p_loc)
  end do
  close(99)
  
  open(UNIT=102,file='HYDROGEN.DATA')
  write(102,358)N*100,N*100
    do p_loc=1,N_of_hydrogen
        if(HYDROGEN_EXISTED(p_loc))write(102,358)HYDROGEN_X(p_loc),HYDROGEN_Y(p_loc)
    end do
  close(102)
  
  open(UNIT=100,file='fH.DATA')
  write(100,358)N*100,N*100
  do p_loc=1,N_of_particles
    do k_local=5,8
      if(.not. dehydrated(p_loc,k_local))then
        x_local=CORNERS(X(p_loc),k_local,1)
        y_local=CORNERS(Y(p_loc),k_local,2)
        if(abs(x_local-CORNERS(X(p_loc),5,1)) .gt. 2)x_local=CORNERS(X(p_loc),5,1)+2
        if(abs(y_local-CORNERS(Y(p_loc),5,2)) .gt. 2)y_local=CORNERS(Y(p_loc),5,2)+2
        write(100,358)x_local,y_local
      end if
    end do
  end do
  close(100)
  
386 format(4(I0,","))
  open(UNIT=101,file='bond.DATA')
    write(101,386)N*100,N*100,2*N*100,2*N*100
    do p_loc=1,N_of_particles
      do k_local=5,8
          bx=-1
          by=-1
        if(bonds(p_loc,k_local).ne. 0 .and. bonds(p_loc,k_local).ne. -1)then
          x_local=CORNERS(X(p_loc),k_local,1)
          y_local=CORNERS(Y(p_loc),k_local,2)
          if(abs(x_local-CORNERS(X(p_loc),5,1)) .gt. 2)x_local=CORNERS(X(p_loc),5,1)+2
          if(abs(y_local-CORNERS(Y(p_loc),5,2)) .gt. 2)y_local=CORNERS(Y(p_loc),5,2)+2
          if(k_local .eq. 5 .or. k_local .eq. 6) bx=x_local-1_i_p
          if(k_local .eq. 7 .or. k_local .eq. 8) bx=x_local+1_i_p
          if(k_local .eq. 5 .or. k_local .eq. 8) by=y_local-1_i_p
          if(k_local .eq. 6 .or. k_local .eq. 7) by=y_local+1_i_p
          write(101,386)x_local,y_local,bx,by
        end if
      end do
    end do
  close(101)
  
end subroutine

end module
