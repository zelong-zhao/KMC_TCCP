module common_parameters

  implicit none

integer(4), parameter :: i_p=4
integer(4), parameter :: r_p=8

!Variables from old KMC code

integer(i_p) :: counter_K(1:4)
integer(i_p) :: counter_hydrogen_move,counter_not_hydrogen_move
integer(i_p) :: START_KMC_STEPS,MAX_MAX_KMC_STEPS,MAX_SEQ
character*3 :: type_of_simulation

integer(i_p) :: N_clusters
real(r_p) :: Compactness_index

!Variables from events.f90
character :: len_grid*3,state*1
logical :: random_or_manual,H_random,record_H,stop_if_all_de_H,print_map_in_terminal,endless_kmc_step
integer(i_p) :: MAX_KMC_STEPS
integer(i_p) :: counter_type(0:7)
character,dimension(0:7) :: counter_type_name*8
integer(i_p) :: print_pov_freq, print_gnuplot_freq
integer(i_p) :: N_of_hydrogen
real(r_p) :: Kb_eV=8.617385d-5
real(r_p) :: prefactor,U_ev=1.0d-3 !prefactor_of_transition_probability=1.0
real(r_p) :: beta,temp
real(r_p) :: delta_time
real(r_p) :: energy_barriers(0:6)!1 for escape  2 for for association, 3 for dissociation
real(r_p) :: TST(0:6)
integer(i_p),allocatable,dimension(:) :: DIMER_P
integer(i_p),allocatable,dimension(:,:) :: bonds
integer(i_p),allocatable,dimension(:,:) :: LOOKUP
integer(i_p),allocatable,dimension(:) :: LOOKUP_PKP
integer(i_p),allocatable,dimension(:,:) :: LOOKUP_HKH

!variables from mymodule
integer(i_p) :: N, N_of_particles !N size of the lattice, P which particle
integer(i_p) :: l_particle,KMC_STEPS
integer(i_p) :: sequence_of_events !max Number of particlr * 4.
integer(i_p),allocatable,dimension(:) :: X
integer(i_p),allocatable,dimension(:) :: Y

integer(i_p) idum

logical,allocatable,dimension(:,:) :: dehydrated
integer(i_p),allocatable,dimension(:,:) :: OCCUPATION
integer(i_p),allocatable,dimension(:) :: HYDROGEN_X, HYDROGEN_Y
logical,allocatable,dimension(:) :: HYDROGEN_EXISTED
integer(i_p),dimension(1:4,1:2) :: K_ARRAY
real(r_p),allocatable,dimension(:) :: R
real(r_p) :: R_c

!from TPD
real(r_p) :: time_increa,temp_min,temp_max,temp_increa,temp_per_plot,P_STAY
real(r_p) :: temp_plot
real(r_p),allocatable,dimension(:) :: Prob

contains

subroutine initialise()
  counter_K=0
  counter_type_name(0)='diff 1-2'
  counter_type_name(1)='diff 3-4'
  counter_type_name(2)='associat'
  counter_type_name(3)='dehydrat'
  counter_type_name(4)='h_moveme'
  counter_type_name(5)=' h_gone '
  counter_type_name(6)=' hydrate'
  counter_type_name(7)='  stayed'

  !K_ARRAY(K,X or Y)
  K_ARRAY=0; K_ARRAY(1,1)=1 !X+1,K=1
  K_ARRAY(2,1)=-1 !X-1,
  K_ARRAY(3,2)=1 !Y+1
  K_ARRAY(4,2)=-1 !Y-1
  counter_hydrogen_move=0
  counter_not_hydrogen_move=0
  
  temp_plot=0
    
end subroutine


end module
