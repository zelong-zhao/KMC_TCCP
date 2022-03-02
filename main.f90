program main

    use common_parameters !initialise
    use events !open_kmc_files
    use memory !allocation,common_input,set_print_format
    use lib !all the functions and subrountines to analysis
    use kmc_lib !KMC_EVENTS,choose_kmc_events,!update_kmc_events,how_to_print
    use TPD_lib !tpd_lib
    
implicit none

call initialise() !predefined value !call before input
call common_input() ! parameters read from input.f90

if (type_of_simulation .eq. 'NOR') then
  call kmc_main()
else if(type_of_simulation .eq. 'TPD') then
  call tpd_main()
else if(type_of_simulation .eq. 'CLU') then
  call NCLU_main()
else if(type_of_simulation .eq. 'POV') then
  call pov_main()
else if(type_of_simulation .eq. 'REP') then
  call replot_CLU()
end if

contains

subroutine pov_main()
integer(i_p) :: frame
logical :: if_call_pov
  write(*,*)'------------------------------------------'
  write(*,*)'post processing'
  write(*,*)'frames to print: '
  read(*,*)frame
  call load_frame(frame)
  call export_to_povray()
  write(*,*)'T call povray F -not'
  read(*,*)if_call_pov
  if(if_call_pov) then
    call system('cp /Users/zelongzhao/Google_Drive/KMC_project/code/code_ZZ/tools/{povray.pov,GRID.dat,RES.ini} ./')
    call system('povray RES.ini povray.pov')
    call system('rm RES.ini bond.DATA povray.pov GRID.dat fH.DATA POV.DATA HYDROGEN.DATA')
  end if
end subroutine

subroutine kmc_main()

call open_kmc_files() ! CONECTIVITY.DATA,stat_moves.DATA
call allocation(state) ! allocate arrays

  do !conditions to stop
    call KMC_EVENTS(KMC_STEPS,R_c,MAX_SEQ)
    call choose_kmc_events(MAX_SEQ,R_c,SEQUENCE_OF_EVENTS)
    call update_kmc_events(SEQUENCE_OF_EVENTS)
    call how_to_print(SEQUENCE_OF_EVENTS,.false.)
    call condition_to_stop()
  end do !KMC_LOOP

end subroutine


subroutine tpd_main()

call input_TPD() ! tpd parameters input_TPD.dat
call open_kmc_files() ! open some files
call open_TPD_files() ! open hydrogen_evolution file
call allocation(state) ! allocate arrays New simulations or old


do
  call KMC_EVENTS(KMC_STEPS,R_c,MAX_SEQ)
  call choose_kmc_TPD_events(MAX_SEQ,R_c,SEQUENCE_OF_EVENTS)
  call update_kmc_events(SEQUENCE_OF_EVENTS)
  call how_to_print(SEQUENCE_OF_EVENTS,.false.)
  call condition_to_stop()
end do !KMC_LOOP

end subroutine

subroutine NCLU_main()
  integer(i_p) :: frame
  write(*,*)'post processing'
  write(*,*)'frames to calculate compactness index: '
  read(*,*)frame
  call load_frame(frame)
  call Numb_clusters(N_clusters,Compactness_index,.true.)
  write(*,*)'------------------------------------------'
  write(*,*)'N_of_clusters',N_clusters
  write(*,*)'------------------------------------------'
  write(*,*)'------------------------------------------'
  write(*,*)'      N_dimer',N_dimer()
  write(*,*)'------------------------------------------'
end subroutine

subroutine Connectivity_index()
  

end subroutine

subroutine replot_CLU()
  real(r_p) :: delta_time_local,Compactness_index_LOCAL
  integer(i_p) :: KMC_STEPS_LOCAL,N_clusters_LOCAL,Reason,i,nb_lines,N_bands
  write(*,*)'post processing'
  
  OPEN(77,file='COMPACTNESS_INDEX.DATA')
  open(88,file='NEW_COMPACTNESS_INDEX.DATA')
  write(88,'(4a)')'#  TIME    ','    KMC_STEP  ',' N_bands ',' COMPACTNESS_INDEX  '
  read(77,*)
  
  i=0
  do
    i=i+1
    read(77,*,IOSTAT=Reason)delta_time_local,KMC_STEPS_LOCAL,N_clusters_LOCAL,Compactness_index_LOCAL
    call load_frame(KMC_STEPS_LOCAL)
    call get_N_bands(N_bands)
    
    write(88,'(e12.6,2(5x,i0),5x,f8.3)')delta_time_local,KMC_STEPS_LOCAL,N_bands,2.0_r_p*dble(N_bands)/dble(N_of_Particles)
    !write(*,*)KMC_STEPS_LOCAL,N_bands,2.0_r_p*dble(N_bands)/dble(N_of_Particles)
    if(Reason .lt. 0) then
      nb_lines=i
      write(*,*)"NUMBER OF LINES:",nb_lines
      exit
    end if
    
    call deallocate_everything()
    
  end do
    write(*,*)'NEW_COMPACTNESS_INDEX.DATA is produced'
  
  close(88)
  close(77)
 !
end subroutine


subroutine condition_to_stop()
  logical :: simulation_stop_satisfied
  simulation_stop_satisfied=.false.

  if(KMC_STEPS .ge. MAX_MAX_KMC_STEPS)simulation_stop_satisfied=.true.
  
  if(stop_if_all_de_H .and. (counter_type(5) .eq. N_of_hydrogen/2)) then
    simulation_stop_satisfied=.true.
    write(*,*)'STOP,ALL MOLECULES ARE DEHYDROGENATED'
  end if
  
  if(N_dimer() .eq. N_of_Particles) then
    simulation_stop_satisfied=.true.
    write(*,*)'STOP,ALL DIMER'
  end if
  
  if(type_of_simulation .eq. 'TPD') then
    if(temp .ge. temp_max)simulation_stop_satisfied=.true.
  end if
  

  if(simulation_stop_satisfied)then
    call how_to_print(SEQUENCE_OF_EVENTS,.true.)
    close(4) !CONECTIVITY.DATA
    close(20) !stat_moves.DATA
    close(77) !compactness index
    close(90) !COMPLETE_HYDROGEN.DATA
    close(89) !CompleteDeH.DATA
    if(record_H) close(10)
    write(*,*)'==================================================='
    if(type_of_simulation .eq. 'TPD') then
      close(22)
      write(*,*)'TEMP',temp,'TEMP_MAX',temp_max
      write(*,*)'STOP,TPD'
    end if
    write(*,*)'KMC_STEPS',KMC_STEPS
    write(*,*)'TIME',delta_time
    STOP
  end if
end subroutine

end program


