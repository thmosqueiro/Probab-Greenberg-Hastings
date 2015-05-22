!
! Probab-Greenberg-Hastings simulation
! 

program avalanches

  implicit real*8(a-h,o-z)

  character(len=200)     :: bugger, log_file, date, data_filename
  integer                :: Nneurons, Nstates, sizer, Tsimul, pr, Nrs, av, lifetime, maior_av, maior_lf
  real(8)                :: ext_stimulus, pmax, sigma, fraction, r, r0, h, dr, Avg_Conn, mbr
  real(8), allocatable   :: mat_adj(:,:), rho(:)
  integer, allocatable   :: netstate(:), seed(:)
  
! Funcoes
  real(8)                :: getActivity, getMeanConnectivity, getMeanBR, ran1, mean_activity
  integer                :: getMeanNofConnections, get_avalanche



! Lendo o seed como uma flag
  call getarg(1, bugger)
  read(bugger,*) sizer

  pr = 8
  allocate(seed(pr))
  CALL RANDOM_SEED (SIZE=pr)
  seed = sizer + 37 * (/ (i - 1, i = 1, pr) /)
  CALL RANDOM_SEED (PUT=seed)

! Lendo os outros argumentos
  call getarg(2, bugger)
  read(bugger,*) Nneurons
  call getarg(3, bugger)
  read(bugger,*) Ntrials
  call getarg(4, bugger)
  read(bugger,*) sigma
  call getarg(5, bugger)
  read(bugger,*) Avg_Conn
  call getarg(6, bugger)
  read(bugger,*) Tsimul
  call getarg(7, bugger)
  read(bugger,*) data_filename
  call getarg(8, bugger)
  read(bugger,*) date
  call getarg(9, bugger)
  read(bugger,*) log_file


! Definicao de parametros -- parte de implementacao
  Nstates = 10
  ext_stimulus = 0.0d0 !1d0 - dexp(0.5d0)

! Calculando grandezas derivadas
  pmax = 2d0*sigma/Avg_Conn

  allocate(netstate(Nneurons), mat_adj(Nneurons, 0:5*int(Avg_Conn)), rho(Nneurons))

!!! Setting the initial values
  CALL setup_network(Nneurons, netstate, fraction)


!!! Construindo o grafo
!! Imprimindo LOG
  write(*,*) "- Construindo topologia com pedgree..."
  write(*,"(A)",advance='no') " - Construindo a topologia..."
  mbr = 0.d0

  
  do while ( dabs(mbr - sigma) > 1d-3)
     CALL connect_network_selection(Nneurons, Avg_Conn, sigma, mat_adj, M, 1, 0.d0)
     mbr = getMeanBR(Nneurons, M, mat_adj)
  end do
  
  cont_PASS = 0
  do j = 1, Nneurons
     rho(j) = 0.d0
  enddo
  maior_av = 0
  maior_lf = 0

  write(1,*) '=============================================================='  
  write(1,*) 
  write(1,*) '          Log da simulacao!' 
  write(1,*) 
  write(1,*) ' Todos os parametros estao arquivados!' 
  write(1,*) 
  write(1,*) '==============================================================' 
  write(1,*) 
  write(1,*) 
  write(1,*) 'N          = ', Nneurons
  write(1,*) 'K          = ', Avg_Conn
  write(1,*) 'm          = ', Nstates
  write(1,*) 'Tsimul     = ', Tsimul
  write(1,*) 'Ntrials    = ', ntrials
  write(1,*) 'Madj       = ', 5*int(Avg_Conn)
  write(1,*) 'Seed       = ', seed
  write(1,*)
  write(1,103) getMeanConnectivity(Nneurons, 5*int(Avg_Conn), mat_adj) 
  write(1,106) getMeanBR(Nneurons, 5*int(Avg_Conn), mat_adj)
  
!!! Comecando a simulacao
  write(*,"(A)",advance='no') " - Simulacao correndo... "
  
  open(27, file=data_filename)
  
  write(*,*)
  write(*,*)
  
  do j = 1, Ntrials
     
     CALL exciteOneNeuron(Nneurons, netstate)
     
     av = get_avalanche(Nneurons, M, netstate, mat_adj, Nstates, ext_stimulus, Tsimul, lifetime, rho)

     ! If you want more performance: commend this next line!!
     write(*,*) 'Iteracao:', j, 'Resultado:', av, 'Lifetime:', lifetime
     
     write(27,'(A)') " =========================================== "
     write(27,'(A, I8)') "  Size:   ", av
     write(27,'(A, I8)') "  Tempo:  ", lifetime
     write(27,'(A)') " =========================================== "
     write(27,*) 
     write(27,*) 
     
     if ( maior_av .lt. av ) maior_av = av
     if ( maior_lf .lt. lifetime ) maior_lf = av
     
  enddo
  

  close(27)


  write(1,*)
  write(1,*) 'Largest Avch = ', maior_av
  write(1,*) 'Longest Avch = ', maior_lf
  write(1,*)
  write(1,*)
  write(1,*) '=============================================================='  


!!! Desalocando as matrizes e vetores  
  deallocate(netstate, mat_adj)

!!! Formatacoes
100 format("(A)")
102 format(" - Atividade inicial: ", F4.3)
105 format(" - Conectividade media: ", F6.2)
104 format(" - Peso maximo de cada noh: ", F4.3)
103 format(" - Pronto. Conectividade media empirica: ", F6.3)
106 format(" - Average branching ratio: ", F8.5)

  stop
end program avalanches
