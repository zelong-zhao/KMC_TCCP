module TPD_lib
  use common_parameters !parameters
  use events
  use memory
  
  implicit none
  
  
  contains
  
subroutine choose_kmc_TPD_events(MAX_SEQ_LOCAL,R_c_local,SEQUENCE_OF_EVENTS_LOCAL)
!GLOBAL VARIABLE
  real(r_p),intent(in) :: R_c_local !Cumulative transition rate
  integer(i_p),intent(in) :: MAX_SEQ_LOCAL
  integer(i_p),intent(out) :: SEQUENCE_OF_EVENTS_LOCAL !what to happen
!LOCAL VARIABLE
  logical :: happened
  integer(i_p) :: H
    if(random_or_manual)then
        happened=.false.
        !do while(.not. happened .and. temp .le. temp_max)
!============================calculate probability==========================
            P_STAY=exp(-R_c_local*time_increa)
!make r then transition rate to be probability
            Prob=R/R_c_local*(1-P_STAY)
!============================CHOSSING EVENTS================================
            call randomly_choose_a_sequence_TPD(SEQUENCE_OF_EVENTS,happened)!this gives me a sequence_events back
            delta_time=delta_time+time_increa !delta_time
            temp=temp+temp_increa*time_increa
            call transition(temp) !update the transition rate array
!            if(.not. happened) then
!                counter_type(7)=counter_type(7)+1
!                call updating_R_and_R_c(MAX_SEQ,TST)
!            end if
        !end do
    else if(.not. random_or_manual)then
        stop 'no such function'
    end if
end subroutine

subroutine randomly_choose_a_sequence_TPD(sequence_of_events_local,happened_local)
![out]sequence_of_events
integer(i_p),intent(out)::sequence_of_events_local
logical,intent(out) :: happened_local

real(r_p) :: KMC_random
real(r_p):: PROB_LOCAL

PROB_LOCAL=P_STAY
sequence_of_events_local=1
KMC_random=ran1(idum)

if(KMC_random .le. P_STAY) then
    happened_local = .false.
    return
else
    do while(KMC_random .gt. PROB_LOCAL)
        sequence_of_events_local=sequence_of_events_local+1
        PROB_LOCAL=PROB_LOCAL+PROB(sequence_of_events_local)
        happened_local = .true.
    end do
end if

end subroutine

subroutine updating_R_and_R_c(MAX_SEQ_LOCAL,TST_LOCAL)
  integer(i_p),intent(in) :: MAX_SEQ_LOCAL
  real(r_p),intent(in) :: TST_LOCAL(0:6)
! local variable
  integer(i_p) :: i_local
  
  R_c=0.0
  R=0.0
  do i_local=1,MAX_SEQ_LOCAL
      R(i_local)=TST_LOCAL(LOOKUP(i_local,3))
      R_c=R_c+R(i_local)
  end do

end subroutine

subroutine input_TPD()

integer :: i_local

open(33,file='input_TPD.dat')

read(33,*) time_increa
read(33,*) temp_min
read(33,*) temp_max
read(33,*) temp_increa
read(33,*) temp_per_plot

!initialize

temp=temp_min
!=========TMTPP============
write(*,*)'......>   temp,temp_max   <......'
write(*,*)temp,temp_max
write(*,*)'......> Time Increa(s/kmc)<......'
write(*,*)time_increa
write(*,*)'......> Temp Increa(k/s)  <......'
write(*,*)temp_increa
write(*,*)'......> Temp Increa(k/kmc)<......'
write(*,*)time_increa*temp_increa
write(*,*)'......> MAX KMC STEPS     <......'
write(*,*)(temp_max-temp_min)/(time_increa*temp_increa)

call transition(temp)
do i_local=0,6
    write(*,'(a,e12.5)') counter_type_name(i_local)//' rate= ',TST(i_local)
end do

close(33)
end subroutine


subroutine open_TPD_files()

if(state .eq. 'N') then
  open(22,file='hydrogen_evolution.DATA')
  write(22,'(4a)')'#  TIME    ',' TEMPERATURE','    P_Stay  ',(' '//counter_type_name(5)//'  ')
else if(state .eq. 'C')then
  open(22,file='hydrogen_evolution.DATA',STATUS='OLD',ACCESS='APPEND')
end if

end subroutine


end module
