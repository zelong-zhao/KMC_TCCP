module memory
  use common_parameters
  use events
  use lib

implicit none
  contains
 
subroutine common_input()

  integer :: i_local
  
  open(3,file='input.dat')
  read(3,*) type_of_simulation
  read(3,*) state
  read(3,*) random_or_manual
  read(3,*) N
  read(3,*) N_of_particles
  read(3,*) temp
  read(3,*) MAX_KMC_STEPS
  read(3,*) l_particle
  read(3,*) energy_barriers(0)
  read(3,*) energy_barriers(1)
  read(3,*) energy_barriers(2)
  read(3,*) energy_barriers(3)
  read(3,*) energy_barriers(4)!hydrogen move
  read(3,*) energy_barriers(5)
  read(3,*) energy_barriers(6)
  read(3,*) prefactor
  read(3,*) print_pov_freq
  read(3,*) print_gnuplot_freq
  read(3,*) idum
  read(3,*) H_random
  read(3,*) record_H
  read(3,*) stop_if_all_de_H
  read(3,*) print_map_in_terminal
  read(3,*) endless_kmc_step

  write(*,*)'......> '//type_of_simulation//'               <......'
  if(random_or_manual) then
    write(*,*)'......> Random run        <......'
  else
    write(*,*)'......> Manual run        <......'
  end if
  if(H_random) then
    write(*,*)'......> H are at random   <......'
  else
    write(*,*)'......> H moves included  <......'
  end if
  if(record_H) then
    write(*,*)'......> Record H - TRUE   <......'
  else
    write(*,*)'......> Record H - FALSE  <......'
  end if
  if(stop_if_all_de_H) then
    write(*,*)'......> Stop if ALL de-H  <......'
  else
    write(*,*)'......> Cont. if ALL de-H <......'
  end if
  
  if(endless_kmc_step) then
    if(type_of_simulation .ne. 'TPD' .and. .not. stop_if_all_de_H) stop 'turn off endless kmc step'
  end if
  if(record_H .and. H_random) stop 'turn off record_H'
 

 N_of_hydrogen=4*N_of_particles
 delta_time=0.0
 !=========TMTPP============
 if(type_of_simulation .eq. 'NOR') then
   call transition(temp)
   write(*,'("beta",f15.10)')beta
   beta=1.0/(temp*Kb_eV)
   
   do i_local=0,6
       write(*,'(a,e10.5)') counter_type_name(i_local)//' rate= ',TST(i_local)
   end do
 end if
 close(3)
 
 call set_print_format()
 end subroutine
 
  subroutine set_print_format()
   !write(*,*) N+2
   if(N+2 .le. 9) then
      write(len_grid,'(2x,i1)') N+2
   else if(N+2 .le. 99) then
      write(len_grid,'(x,i2)') N+2
   else if(N+2 .le. 999) then
      write(len_grid,'(i3)') N+2
   end if
   if(N .gt. 99 .and. print_map_in_terminal) stop 'turn off print map in terminal'
  end subroutine set_print_format

  subroutine transition(temp_local)
  real(r_p),intent(in) :: temp_local
  integer :: i_local
  beta=1.0/(temp_local*Kb_eV)
  do i_local=0,6
    TST(i_local)=exp(-beta*energy_barriers(i_local)+log(prefactor))
  end do
    !write(*,*)'=============================='
  end subroutine

  subroutine allocation(state_local)
  !global variable
  character,intent(in) :: state_local*1 !N or C !N means New, C means continue
  
  !local variable
  integer(i_p) :: P,i_local,N_of_hydrogen_on_surface_local
  logical :: P_OCCUPIED

  !MAP
  allocate(OCCUPATION(0:N,0:N));  OCCUPATION=0
  allocate(X(1:N_of_particles));  X=0
  allocate(Y(1:N_of_particles));  Y=0
  allocate(HYDROGEN_X(1:N_of_hydrogen)) ;HYDROGEN_X=-1
  allocate(HYDROGEN_Y(1:N_of_hydrogen)) ;HYDROGEN_Y=-1
  allocate(HYDROGEN_EXISTED(1:N_of_hydrogen)) ;HYDROGEN_EXISTED=.false.

  !LOOKUP
  allocate(LOOKUP_PKP(1:8*N_of_particles)); LOOKUP_PKP=0
  allocate(LOOKUP_HKH(1:8*N_of_hydrogen,1:2));LOOKUP_HKH=0
  allocate(LOOKUP(1:8*N_of_particles+8*N_of_hydrogen,1:3)); LOOKUP=0
  allocate(R(1:N_of_particles*8+8*N_of_hydrogen)); R=0.0d0
  !TPD
  allocate(Prob(1:N_of_particles*8+8*N_of_hydrogen)); Prob=0.0d0

  !TMTPP bonds & de
  allocate(DIMER_P(1:N_of_particles));  DIMER_P=0
  allocate(bonds(1:N_of_particles,5:8));  bonds=0
  allocate(dehydrated(1:N_of_particles,5:8)); dehydrated=.false.




  if(state_local .eq. 'N')then

  !=============check point==========

    delta_time=0.0_r_p
    counter_type=0_i_p
    KMC_STEPS=0_i_p;
    P=1
    do while (P .le. N_of_particles)
       X(P) = nint(N*ran1(idum))
       Y(P) = nint(N*ran1(idum))
       call CHECK_OCCUPANCY(P,0,P_OCCUPIED)!check whether P_OCCUPIED
       if(.not. P_OCCUPIED) then
          call UPDATE_PARTICLE_SHAPE(P,1)
          P=P+1
       end if
    end do
    START_KMC_STEPS=0
    KMC_STEPS=START_KMC_STEPS
    MAX_MAX_KMC_STEPS=MAX_KMC_STEPS
    
    if(.not. random_or_manual)then
       call print_()
       write(*,'(a)')'-----------------------------------------------------------'
    end if

  else if(state_local .eq. 'C')then
      OPEN(UNIT=7,FILE='DATA_LAST.dat')
      read(7,*)
      read(7,*)START_KMC_STEPS,delta_time,temp
      read(7,*)
      read(7,*)counter_type(0:6)
      read(7,*)
      do  i_local=1,N_of_particles
          read(7,*)X(i_local),Y(i_local),DIMER_P(i_local),bonds(i_local,5:8), &
              dehydrated(i_local,5:8)
          call UPDATE_PARTICLE_SHAPE(i_local,1)
      end do
      read(7,*)
      read(7,*)N_of_hydrogen_on_surface_local
      read(7,*)
      if(N_of_hydrogen_on_surface_local .ne. 0) then
        do i_local=1,N_of_hydrogen_on_surface_local
            read(7,*)HYDROGEN_X(i_local),HYDROGEN_Y(i_local),HYDROGEN_EXISTED(i_local)
        end do
      end if
      read(7,*)
      read(7,*) idum
      close(7)

      KMC_STEPS=START_KMC_STEPS
      MAX_MAX_KMC_STEPS=MAX_KMC_STEPS+KMC_STEPS

      call print_()
      write(*,'(a)')'-----------------------------------------------------------'
      write(*,*)'   KMC step:',KMC_STEPS
      write(*,*)' delta_time:',delta_time
      write(*,*)'temperature:',temp
  else
    stop 'error on input'
  end if
  end subroutine
  
  subroutine load_frame(frame_to_print)
  
  integer(i_p),intent(in) :: frame_to_print
  character :: what_to_load*100
  
  !local variable
  integer(i_p) :: P,i_local,N_of_hydrogen_on_surface_local
  logical :: P_OCCUPIED

  !MAP
  allocate(OCCUPATION(0:N,0:N));  OCCUPATION=0
  allocate(X(1:N_of_particles));  X=0
  allocate(Y(1:N_of_particles));  Y=0
  allocate(HYDROGEN_X(1:N_of_hydrogen)) ;HYDROGEN_X=-1
  allocate(HYDROGEN_Y(1:N_of_hydrogen)) ;HYDROGEN_Y=-1
  allocate(HYDROGEN_EXISTED(1:N_of_hydrogen)) ;HYDROGEN_EXISTED=.false.

  !LOOKUP
  allocate(LOOKUP_PKP(1:8*N_of_particles)); LOOKUP_PKP=0
  allocate(LOOKUP_HKH(1:8*N_of_hydrogen,1:2));LOOKUP_HKH=0
  allocate(LOOKUP(1:8*N_of_particles+8*N_of_hydrogen,1:3)); LOOKUP=0
  allocate(R(1:N_of_particles*8+8*N_of_hydrogen)); R=0.0d0
  !TPD
  allocate(Prob(1:N_of_particles*8+8*N_of_hydrogen)); Prob=0.0d0

  !TMTPP bonds & de
  allocate(DIMER_P(1:N_of_particles));  DIMER_P=0
  allocate(bonds(1:N_of_particles,5:8));  bonds=0
  allocate(dehydrated(1:N_of_particles,5:8)); dehydrated=.false.

  
  write(what_to_load,'("./frames/",I0,".frame")')frame_to_print
  
  OPEN(UNIT=7,FILE=what_to_load)
  read(7,*)
  read(7,*)START_KMC_STEPS,delta_time,temp
  read(7,*)
  read(7,*)counter_type(0:6)
  read(7,*)
  do  i_local=1,N_of_particles
      read(7,*)X(i_local),Y(i_local),DIMER_P(i_local),bonds(i_local,5:8), &
          dehydrated(i_local,5:8)
      call UPDATE_PARTICLE_SHAPE(i_local,1)
  end do
  read(7,*)
  read(7,*)N_of_hydrogen_on_surface_local
  read(7,*)
  if(N_of_hydrogen_on_surface_local .ne. 0) then
    do i_local=1,N_of_hydrogen_on_surface_local
        read(7,*)HYDROGEN_X(i_local),HYDROGEN_Y(i_local),HYDROGEN_EXISTED(i_local)
    end do
  end if
  read(7,*)
  read(7,*) idum
  close(7)

  KMC_STEPS=START_KMC_STEPS
  MAX_MAX_KMC_STEPS=MAX_KMC_STEPS+KMC_STEPS

!  call print_()
!  write(*,'(a)')'-----------------------------------------------------------'
!  write(*,*)'   KMC step:',KMC_STEPS
!  write(*,*)' delta_time:',delta_time
!  write(*,*)'temperature:',temp

  end subroutine
  
  subroutine deallocate_everything()
   deallocate(OCCUPATION,X,Y,HYDROGEN_X,HYDROGEN_Y,HYDROGEN_EXISTED,LOOKUP_PKP,LOOKUP_HKH,LOOKUP,R,Prob,DIMER_P,bonds,dehydrated)
  end subroutine
  
  
 
 end module
