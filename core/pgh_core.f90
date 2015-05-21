!
! Update
!
! Salva a atividade em um arquivo externo
!
subroutine update(N, Madj, netstate, mat_adj, m, h, Nsteps)

  integer, intent(in)    :: N, Nsteps, m, Madj
  real(8), intent(in)    :: h, mat_adj(N, 0:Madj)
  integer                :: contador, aux, NN, j, jj, p, progressbar, netstate_new(N)
!  real(8)                :: prob, auxp
  integer, intent(inout) :: netstate(N)

  ! External functions
  real(8)                :: getActivity, ran1

  open(unit=1, file="saving", access="append") !action

! Como o algoritmo eh sincrono, entao eh necessario
! guardar os estados atualizados em um novo vetor.
  netstate_new = netstate

! Primeira marcacao de atividade
  write(1,*) 0, getActivity(N, netstate)

! Setando o progressbar
  contador = 0
  progressbar = Nsteps/30

!  write(*,"(A)",advance='no') " >> ["

  do j = 1, Nsteps

     do p = 1, N


        ! Neuronio pode ser excitado
        if ( ( netstate(p) .eq. 0 ) .and. ( netstate_new(p) .eq. 0 ) ) then
           
           ! O estimulo externo pode excitar o neuronio
           if ( ran1() .lt. h ) then
              netstate_new(p) = 1
           end if
           
        else if ( netstate(p) .eq.  1 ) then
           
           netstate_new(p) = 2
                      
           ! Para nao ficar acessando desnecessariamente mat_adj(p, 0)
           aux = int( mat_adj(p, 0) )
           
           do jj = 1, aux

              idx = int( mat_adj(p,jj) )
              
              if ( ( netstate( idx ) + netstate_new( idx ) .eq. 0 ) .and. ( ran1() .lt. mod( mat_adj(p,jj), 1d0) ) ) then
                 netstate_new( idx ) = 1
              end if

!              write(*,*) jj, aux,  mod( mat_adj(p,jj), 1d0)
!              read(*,*)

           end do

        else
           
           netstate_new(p) = modulo(netstate(p) + 1, m)
           
        end if
        
        
     end do

     ! Atualizando o estado da rede
     netstate = netstate_new

     ! Imprimindo em arquivo a historia
     write(1,*) j, getActivity(N, netstate)

!     CALL progressbar_update(contador, progressbar)

  end do

! Finalizando o contador
!  write(*,"(A)") "]"

  close(1)

end subroutine update




!
! mean_activity
!
! Funcao que retorna a atividade media apos um periodo dado
! de iteracoes
!
real(8) function mean_activity(N, Madj, netstate, mat_adj, m, h, Nsteps) result(media)

  integer, intent(in)    :: N, Nsteps, m, Madj
  real(8), intent(in)    :: h, mat_adj(N, 0:Madj)
  integer, intent(inout) :: netstate(N)
  integer                :: netstate_new(N)

! Indices
  integer                :: j, jj, NN, p, aux, idx
  integer                :: contador, progressbar, chave

! Funcoes externas
  real(8)                :: getActivity, ran1
  logical                :: isDead


! Como o algoritmo eh sincrono, entao eh necessario
! guardar os estados atualizados em um novo vetor.
  netstate_new = netstate

! Inicializando a media
  media = 0d0
  
  do j = 1, Nsteps
     
     do p = 1, N
        
        ! Neuronio pode ser excitado
        if ( ( netstate(p) .eq. 0 ) .and. ( netstate_new(p) .eq. 0 ) ) then
           
           ! O estimulo externo pode excitar o neuronio
           if ( ran1() .lt. h ) then
              netstate_new(p) = 1
           end if
           
        else if ( netstate(p) .eq.  1 ) then
           
           ! This was wrong: it don't catch the case where m = 2
           ! netstate_new(p) = 2
           netstate_new(p) = modulo(netstate(p) + 1, m)
                      
           ! Para nao ficar acessando desnecessariamente mat_adj(p, 0)
           aux = int( mat_adj(p, 0) )
           
           
           do jj = 1, aux

              idx = int( mat_adj(p,jj) )
              
              if ( ( netstate( idx ) + netstate_new( idx ) .eq. 0 ) .and. &
                   ( ran1() .lt. mod( mat_adj(p,jj), 1d0) ) ) then
                 
                 netstate_new( idx ) = 1
                 
              end if
              
           end do

        else
           
           netstate_new(p) = modulo(netstate(p) + 1, m)
           
        end if
        
     end do

     ! Se a rede estiver morta e não houver estimulo externo, não adianta prosseguir.
     if ( isDead(N, netstate_new) .and. ( h .eq. 0d0 ) ) then
        
        media = media/float(Nsteps)

        return
     end if

     ! Atualizando o estado da rede
     netstate = netstate_new

     ! Imprimindo em arquivo a historia
     media = media + getActivity(N, netstate)
     
     if ( .FALSE. ) then
        open(unit=27, file="activity.data", access="append")
        write(27,*) 'Iteracao: ', j
        do jj = 1, N
           if ( netstate(jj) .eq. 1 ) then
              write(27,*) jj
           endif
        enddo
        write(27,*) 'fim'
        write(27,*) 

        close(27)
     end if
     
     
  end do

  media = media/float(Nsteps)

  return
end function mean_activity








!
! mean_activity -- OLD
!
! Funcao que retorna a atividade media apos um periodo dado
! de iteracoes
!
real(8) function mean_activity_old(N, Madj, netstate, mat_adj, m, h, Nsteps) result(media)

  integer, intent(in)    :: N, Nsteps, m, Madj
  real(8), intent(in)    :: h, mat_adj(N, 0:Madj)
  integer                :: contador, aux, NN, j, jj, p, progressbar, netstate_new(N)
!  real(8)                :: prob, auxp
  integer, intent(inout) :: netstate(N)
  real(8)                :: getActivity, ran1


! Como o algoritmo eh sincrono, entao eh necessario
! guardar os estados atualizados em um novo vetor.
  netstate_new = netstate

! Inicializando a media
  media = 0d0

  do j = 1, Nsteps

     do p = 1, N
        
        ! Neuronio pode ser excitado
        if ( netstate(p) .eq. 0 ) then

           ! Para nao ficar acessando desnecessariamente mat_adj(p, 0)
           aux = int( mat_adj(p, 0) )
           
           ! O estimulo externo pode excitar o neuronio
           if ( ran1() .lt. h ) then

              netstate_new(p) = 1

           ! Ou seus vizinhos podem excita-lo...
           else

              do jj = 1, aux
                 if ( ( netstate(int( mat_adj(p,jj) )) .eq. 1 ) .and. ( ran1() .lt. mod( mat_adj(p,jj), 1d0 ) ) ) then
                    netstate_new(p) = 1
                 end if
              end do

           end if

        else
           
           netstate_new(p) = modulo(netstate(p) + 1, m)
           
        end if
        
     end do

     ! Atualizando o estado da rede
     netstate = netstate_new

     ! Imprimindo em arquivo a historia
     media = media + getActivity(N, netstate)

  end do

  media = media/float(Nsteps)

end function mean_activity_old




!
! probabilities
!
! Retorna as probabilidades marginais de neuronios estarem excitados
!
subroutine probabilities(N, Madj, netstate, mat_adj, m, h, Nsteps, probab)

  integer, intent(in)    :: N, Nsteps, m, Madj
  real(8), intent(in)    :: h, mat_adj(N, 0:Madj)
  integer                :: contador, aux, NN, j, jj, p, progressbar, netstate_new(N)
  integer, intent(inout) :: netstate(N)
  real(8), external      :: getActivity, ran1
  real(8), intent(inout) :: probab(N)


! Como o algoritmo eh sincrono, entao eh necessario
! guardar os estados atualizados em um novo vetor.
  netstate_new = netstate
  
! Inicializando o vetor de probabilidades
  do j = 1,N
     probab(j) = 0d0
  end do


! Setando o progressbar
  contador = 0
  progressbar = Nsteps/30

  write(*,"(A)",advance='no') " >> ["



! Realizando o experimento
  do j = 1, Nsteps

     do p = 1, N
        
        ! Neuronio pode ser excitado
        if ( netstate(p) .eq. 0 ) then

           ! Para nao ficar acessando desnecessariamente mat_adj(p, 0)
           aux = int( mat_adj(p, 0) )
           
           do jj = 1, aux
              if ( ( netstate(int( mat_adj(p,jj) )) .eq. 1 ) .and. ( ran1() .lt. mod( mat_adj(p,jj), 1d0 ) ) ) then
                 netstate_new(p) = 1
              end if
           end do

        else
           
           netstate_new(p) = modulo(netstate(p) + 1, m)
           
        end if

     end do

     ! Atualizando o estado da rede
     netstate = netstate_new

     CALL progressbar_update(contador, progressbar)

  end do

! Normalizando as probabilidades calculadas
  do j = 1,N
     probab(j) = probab(j)/dfloat(Nsteps)
  end do

! Finalizando o contador
  write(*,"(A)") "]"
  
  return
end subroutine probabilities



!
! death_time
!
! Retorna o tempo em que a rede parou de funcionar sem um
! estimulo externo
!
integer function death_time(N, Madj, netstate, mat_adj, m, Nsteps_max) result(dt)

  integer, intent(in)    :: N, Nsteps_max, m, Madj
  real(8), intent(in)    :: mat_adj(N, 0:Madj)
  integer                :: contador, aux, NN, j, jj, p, progressbar, netstate_new(N), idx
  integer, intent(inout) :: netstate(N)
  logical                :: condicao

  ! funcoes
  real(8)                :: ran1, getActivity
  logical                :: isDead


! Como o algoritmo eh sincrono, entao eh necessario
! guardar os estados atualizados em um novo vetor.
  netstate_new = netstate

  j = 0

  do while (j .lt. Nsteps_max)

     j = j + 1

     do p = 1, N
        
        if ( netstate(p) .eq. 1 ) then
           
           ! Atualizando o proprio neuronio
           netstate_new(p) = 2
           
           ! Para nao ficar acessando desnecessariamente mat_adj(p, 0)
           aux = int( mat_adj(p, 0) )
           
           do jj = 1, aux
           
              ! Indice para nao ficar acessando desnecessariamente mat_adj
              idx = int( mat_adj(p,jj) )
              
              if ( ( netstate( idx ) + netstate_new( idx ) .eq. 0 ) .and. ( ran1() .lt. mod( mat_adj(p,jj), 1d0 ) ) ) then
                 
                 netstate_new(idx) = 1
                 
              end if
           end do
           
        else if ( netstate(p) .gt. 1 ) then
           
           netstate_new(p) = modulo(netstate(p) + 1, m)
           
        end if
        
     end do

     ! Testa se a rede jah morreu
     if ( isDead(N, netstate_new) ) then
        dt = j
        return
     end if

     ! Atualizando o estado da rede
     netstate = netstate_new

  end do

  dt = Nsteps_max

  return
end function death_time




!
! death_time
!
! Retorna o tempo em que a rede parou de funcionar sem um
! estimulo externo
!
subroutine death_time_subsampling(N, Madj, netstate, mat_adj, m, Nsteps_max, nume, ElectrodeIndex, outputfile)

  character(len=100)     :: outputfile
  integer, intent(in)    :: N, Nsteps_max, m, Madj, nume
  real(8), intent(in)    :: mat_adj(N, 0:Madj)
  integer                :: contador, aux, NN, j, jj, p, progressbar, netstate_new(N), idx
  integer, intent(inout) :: netstate(N), ElectrodeIndex(nume)
  logical                :: condicao, Teste_Are_ElecDead, new_av_started

  ! funcoes
  real(8)                :: ran1, getActivity
  logical                :: isDead, isDead_subsampling


! Como o algoritmo eh sincrono, entao eh necessario
! guardar os estados atualizados em um novo vetor.
  netstate_new = netstate

! Arquivo para escrever os resulstados
  open(unit=27, file=outputfile, access="append")

! Iniciando as variaveis que contam o numero de geracoes
! passadas depois do primeiro seed
  j = 0

! jrel conta o numero de novas geracoes desde a ultima avalanche
  jrel = j

! O seed eh colocado em um dos neuronios avistados por um eletrodo,
! caso contrario nao teria se iniciado avalanche nenhuma.
  new_av_started = .true.

  do while (j .lt. Nsteps_max)

     j = j + 1

     do p = 1, N
        
        if ( netstate(p) .eq. 1 ) then
           
           ! Atualizando o proprio neuronio
           netstate_new(p) = 2
           
           ! Para nao ficar acessando desnecessariamente mat_adj(p, 0)
           aux = int( mat_adj(p, 0) ) 
           
           if ( aux .gt. 0 ) then
              
              do jj = 1, aux
           
                 ! Indice para nao ficar acessando desnecessariamente mat_adj
                 idx = int( mat_adj(p,jj) )
                 
                 if ( ( netstate( idx ) + netstate_new( idx ) .eq. 0 ) .and. ( ran1() .lt. mod( mat_adj(p,jj), 1d0 ) ) ) then
                 
                    netstate_new(idx) = 1
                 
                 end if
              end do
           end if
              
        else if ( netstate(p) .gt. 1 ) then
           
           netstate_new(p) = modulo(netstate(p) + 1, m)
           
        end if
        
     end do


     ! Teste nos eletrodos

     ! Pode acontecer de uma avalanche terminar (para os eletrodos) e haver um
     ! espaço de tempo em que nenhum eletrodo sente potencial. Por isso, eh
     ! necessario testar se uma nova avalanche, para os eletrodos, jah se iniciou
     ! antes de marcar mais tempos. Isso eh feito via chaveamento da variavel
     ! logica new_av_started.

     ! Para nao realizar a mesma conta duas vezes:
     Teste_Are_ElecDead = isDead_subsampling(N, netstate_new, nume, ElectrodeIndex)

     ! Testa se hah algum eletrodo ativo
     if ( Teste_Are_ElecDead .and. new_av_started ) then

        write(27,*) j - jrel
        jrel = j

        new_av_started = .false.

     end if

     ! Se nao houver nenhum ativo, testa se uma nova avalanche (para eletrodos)
     ! se iniciou.
     if ( .not. new_av_started ) then
        new_av_started = .not. Teste_Are_ElecDead
     end if


     ! Testa se a rede jah morreu
     ! Testamos se a rede morre apenas no final porque pode ser que a rede morra sem
     ! que os eletrodos sejam excitados em suas ultimas excitacoes.
     if ( isDead(N, netstate_new) ) then

        return

     end if

     ! Em caso de continuacao do algoritmo, atualizando o estado da rede
     netstate = netstate_new

  end do

  close(27)

  return
end subroutine death_time_subsampling






!
! death_time
!
! Retorna o tempo em que a rede parou de funcionar sem um
! estimulo externo
!
subroutine death_time_loose(N, Madj, netstate, mat_adj, m, Nsteps_max, outputfile, pl)

  character(len=100)     :: outputfile
  integer, intent(in)    :: N, Nsteps_max, m, Madj
  real(8), intent(in)    :: mat_adj(N, 0:Madj), pl
  integer                :: contador, aux, NN, j, jj, p, progressbar, netstate_new(N), idx
  integer, intent(inout) :: netstate(N)
  logical                :: condicao, Teste_Are_ElecDead, new_av_started

  ! funcoes
  real(8)                :: ran1, getActivity
  logical                :: isDead, isDead_loose


! Como o algoritmo eh sincrono, entao eh necessario
! guardar os estados atualizados em um novo vetor.
  netstate_new = netstate

! Arquivo para escrever os resulstados
  open(unit=27, file=outputfile, access="append")

! Iniciando as variaveis que contam o numero de geracoes
! passadas depois do primeiro seed
  j = 0

! jrel conta o numero de novas geracoes desde a ultima avalanche
  jrel = j

! O seed eh colocado em um dos neuronios avistados por um eletrodo,
! caso contrario nao teria se iniciado avalanche nenhuma.
  new_av_started = .true.

  do while (j .lt. Nsteps_max)

     j = j + 1

     do p = 1, N
        
        if ( netstate(p) .eq. 1 ) then
           
           ! Atualizando o proprio neuronio
           netstate_new(p) = 2
           
           ! Para nao ficar acessando desnecessariamente mat_adj(p, 0)
           aux = int( mat_adj(p, 0) ) 
           
           if ( aux .gt. 0 ) then
              
              do jj = 1, aux
           
                 ! Indice para nao ficar acessando desnecessariamente mat_adj
                 idx = int( mat_adj(p,jj) )
                 
                 if ( ( netstate( idx ) + netstate_new( idx ) .eq. 0 ) .and. ( ran1() .lt. mod( mat_adj(p,jj), 1d0 ) ) ) then
                 
                    netstate_new(idx) = 1
                 
                 end if
              end do
           end if
              
        else if ( netstate(p) .gt. 1 ) then
           
           netstate_new(p) = modulo(netstate(p) + 1, m)
           
        end if
        
     end do


     ! Teste nos eletrodos

     ! Pode acontecer de uma avalanche terminar (para os eletrodos) e haver um
     ! espaço de tempo em que nenhum eletrodo sente potencial. Por isso, eh
     ! necessario testar se uma nova avalanche, para os eletrodos, jah se iniciou
     ! antes de marcar mais tempos. Isso eh feito via chaveamento da variavel
     ! logica new_av_started.

     ! Para nao realizar a mesma conta duas vezes:
     Teste_Are_ElecDead = isDead_loose(N, netstate_new, pl)

     ! Testa se hah algum eletrodo ativo
     if ( Teste_Are_ElecDead .and. new_av_started ) then

        write(27,*) j - jrel
        jrel = j

        new_av_started = .false.

     end if

     ! Se nao houver nenhum ativo, testa se uma nova avalanche (para eletrodos)
     ! se iniciou.
     if ( .not. new_av_started ) then
        new_av_started = .not. Teste_Are_ElecDead
     end if


     ! Testa se a rede jah morreu
     ! Testamos se a rede morre apenas no final porque pode ser que a rede morra sem
     ! que os eletrodos sejam excitados em suas ultimas excitacoes.
     if ( isDead(N, netstate_new) ) then

        return

     end if

     ! Em caso de continuacao do algoritmo, atualizando o estado da rede
     netstate = netstate_new

  end do

  close(27)

  return
end subroutine death_time_loose






!
! Avalanche duration usando mesma tecnica que Beggs e Plenz
!
! Esta rotina foi completamente baseada na versão final de death_time_subsampling, mas
! com detalhes um pouco mais realisticos se comparados ao experimento de Beggs e Plenz.
! Em vez de ficar observando quando a avalanche morre segundo o eletrodo, nos integraremos
! a atividade em um intervalo dt. De dt em dt, verificaremos a morte da avalanche.
!
! O valor original de dt nos trabalhos de Beggs e Plenz eh 4 ms.
!
subroutine death_time_BeggsPlenz(N, Madj, netstate, mat_adj, m, Nsteps_max, dt, nume, ElectrodeIndex, outputfile)

  character(len=100)     :: outputfile
  integer, intent(in)    :: N, Nsteps_max, m, Madj, nume, dt
  real(8), intent(in)    :: mat_adj(N, 0:Madj)
  integer                :: NN, progressbar, netstate_new(N)
  integer, intent(inout) :: netstate(N), ElectrodeIndex(nume)
  logical                :: condicao, Teste_Are_ElecDead, new_av_started,  atLeastOneBinAlive

  ! Indices
  integer                :: j, jj, p, jrel, idx, aux, jaux, j0
  integer                :: contador

  ! funcoes
  real(8)                :: ran1, getActivity, getActivity_subsampling
  logical                :: isDead, isDead_subsampling


! Como o algoritmo eh sincrono, entao eh necessario
! guardar os estados atualizados em um novo vetor.
  netstate_new = netstate

! Arquivo para escrever os resulstados
  open(unit=27, file=outputfile, access="append")

! Iniciando as variaveis que contam o numero de geracoes
! passadas depois do primeiro seed
  j = 0

! jrel conta o numero de novas geracoes desde a ultima avalanche
  jrel = 2

! jrel comeca com este valor para que a contagem do numero de dts 
! em que houve uma avalanche (nos eletrodos) estar correta.

! O seed eh colocado em um dos neuronios avistados por um eletrodo,
! caso contrario nao teria se iniciado avalanche nenhuma.
  new_av_started = .true.

! Dah conta de incluir a primeira excitacao no primeiro bin
  j0 = 2

  do while (j .lt. Nsteps_max)

     ! Atualizando o tempo da simulacao
     j = j + dt

     ! A menos que ao menos um eletrodo se ascenda durante as dt
     ! proximas evolucoes temporais, serah dado como a nenhuma atividade
     ! Por isso, setamos inicialmente a variavel como falsa.
     atLeastOneBinAlive = .false.

     ! Realiza dt evolucoes temporais
     ! sem se preocupar se a rede morre ou nao.
     do jaux = j0, dt
        
        do p = 1, N
        
           if ( netstate(p) .eq. 1 ) then
              
              ! Atualizando o proprio neuronio
              netstate_new(p) = 2
              
              ! Para nao ficar acessando desnecessariamente mat_adj(p, 0)
              aux = int( mat_adj(p, 0) ) 
              
              if ( aux .gt. 0 ) then
                 
                 do jj = 1, aux
                    
                    ! Indice para nao ficar acessando desnecessariamente mat_adj
                    idx = int( mat_adj(p,jj) )
                    
                    if ( ( netstate( idx ) + netstate_new( idx ) .eq. 0 ) .and. ( ran1() .lt. mod( mat_adj(p,jj), 1d0 ) ) ) then
                       
                       netstate_new(idx) = 1
                       
                    end if
                 end do
              end if
              
           else if ( netstate(p) .gt. 1 ) then
              
              netstate_new(p) = modulo(netstate(p) + 1, m)
              
           end if
           
        end do

        ! Testando se os eletrodos estao ativos nesta rodada
        if ( .not. isDead_subsampling(N, netstate_new, nume, ElectrodeIndex) ) then
           atLeastOneBinAlive = .true.
        end if

        ! Atualizando o estado da rede
        netstate = netstate_new

!        write(*,*) j, jrel, &
!             ' | Electrodes alive: ', .not. isDead_subsampling(N, netstate_new, nume, ElectrodeIndex), &
!             'Network alive: ', .not. isDead(N, netstate_new)
        
     end do
!     read(*,*)
     
     
     ! Teste nos eletrodos
     
     ! Pode acontecer de uma avalanche terminar (para os eletrodos) e haver um
     ! espaço de tempo em que nenhum eletrodo sente potencial. Por isso, eh
     ! necessario testar se uma nova avalanche, para os eletrodos, jah se iniciou
     ! antes de marcar mais tempos. Isso eh feito via chaveamento da variavel
     ! logica new_av_started.
     
     ! Se uma avalanche nao estiver correndo, marca como uma nova avalanche encontrada
     if ( .not. new_av_started ) then
!        new_av_started = .not. isDead_subsampling(N, netstate_new, nume, ElectrodeIndex)
        new_av_started = atLeastOneBinAlive
        jrel = j 
     end if

!    write(*,*) 'Avalanche ocurring: ', new_av_started

     ! Testa durante os ultimos dt tempos nem mesmo um eletrodo foi ativado
     if ( ( .not. atLeastOneBinAlive ) .and. new_av_started ) then

        ! Imprimindo o tamanho da ultima avalanche
        if ( j .eq. 2 ) then
           write(27,*) 2
!           write(*,*) 'Escreveu resultado: ', 2
        else
           write(27,*) j - jrel
!           write(*,*) 'Escreveu resultado: ', j - jrel
        end if

        ! Marcando o termino da ultima avalanche
        new_av_started = .false.

     end if

     ! Testa se a rede jah morreu
     ! Testamos se a rede morre apenas no final porque pode ser que a rede morra sem
     ! que os eletrodos sejam excitados em suas ultimas excitacoes.
     if ( isDead(N, netstate_new) ) then

        ! Pode ocorrer de a rede toda estar morta no ultimo bin temporal,
        ! no entando nos ultimos dts bins ao menos um eletrodo ter se
        ! excitado. Sendo assim, precisamos marcar mais essa avalache.
        if ( atLeastOneBinAlive ) then
           write(27,*) j - jrel
!           write(*,*) 'Ultimo caso: ', j - jrel
        end if

        ! Agora sim basta só retornar
        return
     end if

     ! Nao mais precisa contar como primeiro bin
     j0 = 1

  end do

  close(27)

  return
end subroutine death_time_BeggsPlenz








!
! Calcula estimador para as probabilidades de 
! criação de novos nós ativos
!
subroutine pr_estimate(N, Madj, netstate, mat_adj, m, Nsteps_max, VP, x)

  integer, intent(in)    :: N, Nsteps_max, m, Madj
  real(8), intent(in)    :: mat_adj(N, 0:Madj)

  integer, intent(inout) :: netstate(N)
  real(8), intent(inout) :: VP(0:N), x

  integer                :: contador, contador_old, aux, NN, j, jj, p, netstate_new(N) , dt
  real(8)                :: norma
  integer                :: num_exc(N), estado, nbindex
  logical                :: condicao

  real(8)                :: getActivity, ran1, delta


! Como o algoritmo eh sincrono, entao eh necessario
! guardar os estados atualizados em um novo vetor.
  netstate_new = netstate

! Zerando o vetor de probabilidades e zerando a contagem
! de excitacoes novas
  do j = 1, N
     VP(j) = 0.d0
  enddo

! Como o run jah inicia-se com um cara excitado...
  Vp(1) = 1
  contador = 1
  x = 0

! Indice principal que vamos percorrer
  j = 0

  dt = Nsteps_max

  do while (j .lt. Nsteps_max)

     ! Incremento do indice
     j = j + 1

     ! Zerando a contagem de excitacoes
     do j = 1, N
        num_exc(j) = -1
     enddo

     ! contadores em igualdade
     contador_old = contador

     ! Percorrendo a rede
     do p = 1, N

        estado = netstate(p)

        ! Neuronio excitado
        if ( estado .eq. 1 ) then

           ! Para nao ficar acessando desnecessariamente mat_adj(p, 0)
           aux = int( mat_adj(p, 0) )

           ! Marcando esse noh como participante da avalanche
           num_exc( p ) = 0

           do jj = 1, aux

              ! indice do vizinho
              nbindex = int( mat_adj(p,jj) )

              ! Excita ou nao excita?
              condicao = ( netstate( nbindex ) .eq. 0 ) .and. ( ran1() .lt. mod( mat_adj(p,jj), 1d0 )  )

              if ( condicao ) then

                 netstate_new( nbindex ) = 1
                 num_exc( p ) = num_exc( p ) + 1

              end if

           end do

        else if ( estado .gt. 1 ) then
           
           netstate_new(p) = modulo(netstate(p) + 1, m)
           
        end if
        
     end do
     
     ! Atualizando a probabilidade
     do jj = 1, N
        if ( num_exc(jj) .ge. 0 ) then
           VP( num_exc(jj) ) = VP( num_exc(jj) ) + 1
        end if

        ! Contando numero de caras jah excitados ateh o momento
        contador = contador + delta(1, netstate_new(jj))
        
     enddo

     x = dfloat( contador_old - 1 )/dfloat(contador)

     ! Atualizando o estado da rede
     netstate = netstate_new

     ! Testa se a rede jah morreu
     if ( getActivity(N, netstate_new) .eq. 0.d0 ) then
        dt = j
        j = Nsteps_max + 10
     end if

  end do
  
  ! Finalizando o estimador para a probabilidade
  norma = 1.d0/dfloat(contador)

  do j = 1, N
     VP( num_exc(j) ) = VP( num_exc(j) )*norma
  enddo

  x = x/dt

  return
end subroutine pr_estimate



!
! Routine optimized to get avalanches sizes without repetition
!
!
integer function get_avalanche(N, Madj, netstate, mat_adj, m, h, Nsteps, lifetime, rho) result(soma)

  integer, intent(in)    :: N, Nsteps, m, Madj
  real(8), intent(in)    :: h, mat_adj(N, 0:Madj)
  integer                :: aux, NN, j, jj, p, netstate_new(N), norep(N)
  real(8)                :: prob, auxp, active_neurons
  real(8), intent(inout) :: rho(N)
  integer, intent(inout) :: netstate(N), lifetime
  
! Coisas auxiliares
  integer                :: contador, progressbar, idx
  integer                :: k
  
! Funcoes externas
  real(8)                :: getActivity, ran1
  integer                :: getActiveNeuronsNoRepeat
    
  
  
! Replacing garbage with zeros.
  do j = 1, N
     netstate(j) = 0
     norep(j) = 0
  end do
  
! Exciting one randomly chosen hub neuron
  aux = int( N*ran1() + 1)
  netstate(aux) = 1

! Initial condition
  rho(1) = getActivity(N, netstate)
  
! Contagem inicial dos neuronios 
  soma = getActiveNeuronsNoRepeat(N, netstate, norep)

! Prova dos 9: soma tem que ser, necessariamente, 1 neste
! momento do programa!
!    write(*,*) soma
  
! Como o algoritmo eh sincrono, entao eh necessario
! guardar os estados atualizados em um novo vetor.
  netstate_new = netstate  
  
  do j = 1, Nsteps
     
     do p = 1, N
        
        if ( netstate(p) .eq. 1 ) then
           
           ! Atualizando o proprio neuronio
           netstate_new(p) = modulo(netstate(p) + 1, m)
           !netstate_new(p) = 2 -> este caso nao compreende quando nao ha estados refratarios
           
           ! Para nao ficar acessando desnecessariamente mat_adj(p, 0)
           aux = int( mat_adj(p, 0) )
           
           do jj = 1, aux
              
              ! Indice para nao ficar acessando desnecessariamente mat_adj
              idx = int( mat_adj(p,jj) )
              
              if ( idx .eq. 0 ) then
                 write(*,*) mat_adj(p,0), mat_adj(p,jj), jj
                 read(*,*) 
              end if
              
              if ( ( netstate( idx ) + netstate_new( idx ) .eq. 0 ) .and. ( ran1() .lt. mod( mat_adj(p,jj), 1d0 ) ) ) then
                 
                 netstate_new(idx) = 1
                 
              end if
           end do
           
        else if ( netstate(p) .gt. 1 ) then
           
           netstate_new(p) = modulo(netstate(p) + 1, m)
           
        end if
        
     end do
     
     
     ! Atualizando o estado da rede
     netstate = netstate_new
     
     ! Calculando quantos neuronios "novos" participaram na
     ! ultima iteracao!
     soma = soma + getActiveNeuronsNoRepeat(N, netstate, norep)
     
     ! Verificando se nao ha nenhum neuronio ativo, contando 
     ! as repeticoes.
     active_neurons = getActivity(N, netstate)

     rho(j + 1) = active_neurons

     if ( active_neurons .eq. 0 ) then
        lifetime = j
        return
     end if
     
  end do

  lifetime = Nsteps
  return
end function get_avalanche







!
! Routine optimized to get avalanches sizes with repetition
!
!
integer function get_avalanche_wrep(N, Madj, netstate, mat_adj, m, h, Nsteps) result(soma)

  integer, intent(in)    :: N, Nsteps, m, Madj
  real(8), intent(in)    :: h, mat_adj(N, 0:Madj)
  integer                :: aux, NN, j, jj, p, netstate_new(N)
  real(8)                :: prob, auxp, active_neurons
  integer, intent(inout) :: netstate(N)
  
! Coisas auxiliares
  integer                :: contador, progressbar, idx
  integer                :: k
  
! Funcoes externas
  real(8)                :: getActivity, ran1
  integer                :: getActiveNeurons
  
  
! Replacing garbage with zeros.
  do j = 1, N
     netstate(j) = 0
  end do
  
! Exciting one randomly chosen hub neuron
  aux = int( N*ran1() + 1)
  netstate(aux) = 1
  
! Contagem inicial dos neuronios 
  soma = getActiveNeurons(N, netstate)

! Prova dos 9: soma tem que ser, necessariamente, 1 neste
! momento do programa!
!    write(*,*) soma
  
! Como o algoritmo eh sincrono, entao eh necessario
! guardar os estados atualizados em um novo vetor.
  netstate_new = netstate
  
  do j = 1, Nsteps
     
     do p = 1, N
        
        if ( netstate(p) .eq. 1 ) then
           
           ! Atualizando o proprio neuronio
           netstate_new(p) = 2
           
           ! Para nao ficar acessando desnecessariamente mat_adj(p, 0)
           aux = int( mat_adj(p, 0) )
           
           do jj = 1, aux
              
              ! Indice para nao ficar acessando desnecessariamente mat_adj
              idx = int( mat_adj(p,jj) )
              
              if ( ( netstate( idx ) + netstate_new( idx ) .eq. 0 ) .and. ( ran1() .lt. mod( mat_adj(p,jj), 1d0 ) ) ) then
                 
                 netstate_new(idx) = 1
                 
              end if
           end do
           
        else if ( netstate(p) .gt. 1 ) then
           
           netstate_new(p) = modulo(netstate(p) + 1, m)
           
        end if
        
     end do
     
     
     ! Atualizando o estado da rede
     netstate = netstate_new
     
     ! Calculando quantos neuronios "novos" participaram na
     ! ultima iteracao!
     soma = soma + getActiveNeurons(N, netstate)
     
     ! Verificando se nao ha nenhum neuronio ativo, contando 
     ! as repeticoes.
     active_neurons = getActivity(N, netstate)
     if ( active_neurons .eq. 0 ) then
        return
     end if
     
  end do
  
  return
end function get_avalanche_wrep




!
! death_time
!
! Retorna o tempo em que a rede parou de funcionar sem um
! estimulo externo
!
subroutine sizeavalanche_loose(N, Madj, netstate, mat_adj, m, Nsteps_max, outputfile, pl)

  character(len=100)     :: outputfile
  integer, intent(in)    :: N, Nsteps_max, m, Madj
  real(8), intent(in)    :: mat_adj(N, 0:Madj), pl
  integer                :: contador, aux, NN, j, jj, p, progressbar, netstate_new(N), idx, soma
  integer                :: norep(N)
  integer, intent(inout) :: netstate(N)
  logical                :: condicao, Teste_Are_ElecDead, new_av_started

  ! funcoes
  real(8)                :: ran1, getActivity
  logical                :: isDead, isDead_loose


! Como o algoritmo eh sincrono, entao eh necessario
! guardar os estados atualizados em um novo vetor.
  netstate_new = netstate


! Defining the norep vector
  do j = 1, N
     norep(j) = 0
  end do

! Arquivo para escrever os resulstados
  open(unit=27, file=outputfile, access="append")

! Iniciando as variaveis que contam o numero de geracoes
! passadas depois do primeiro seed
  j = 0

! Starting the avalanche size properly
  soma = 1

! O seed eh colocado em um dos neuronios avistados por um eletrodo,
! caso contrario nao teria se iniciado avalanche nenhuma.
  new_av_started = .true.

  do while (j .lt. Nsteps_max)

     j = j + 1

     do p = 1, N
        
        if ( netstate(p) .eq. 1 ) then
           
           ! Atualizando o proprio neuronio
           netstate_new(p) = 2
           
           ! Para nao ficar acessando desnecessariamente mat_adj(p, 0)
           aux = int( mat_adj(p, 0) ) 
           
           if ( aux .gt. 0 ) then
              
              do jj = 1, aux
           
                 ! Indice para nao ficar acessando desnecessariamente mat_adj
                 idx = int( mat_adj(p,jj) )
                 
                 if ( ( netstate( idx ) + netstate_new( idx ) .eq. 0 ) .and. ( ran1() .lt. mod( mat_adj(p,jj), 1d0 ) ) ) then
                 
                    netstate_new(idx) = 1
                 
                 end if
              end do
           end if
              
        else if ( netstate(p) .gt. 1 ) then
           
           netstate_new(p) = modulo(netstate(p) + 1, m)
           
        end if
        
     end do


     ! Teste nos eletrodos

     ! Pode acontecer de uma avalanche terminar (para os eletrodos) e haver um
     ! espaço de tempo em que nenhum eletrodo sente potencial. Por isso, eh
     ! necessario testar se uma nova avalanche, para os eletrodos, jah se iniciou
     ! antes de marcar mais tempos. Isso eh feito via chaveamento da variavel
     ! logica new_av_started.

     ! Calculando quantos neuronios "novos" participaram na
     ! ultima iteracao!
     call getActiveNeuronsNoRepeat_loose(N, netstate_new, norep, pl, soma, Teste_Are_ElecDead)

     ! Testa se hah algum eletrodo ativo
     if ( Teste_Are_ElecDead .and. new_av_started ) then

        write(27,*) soma
        soma = 0

        ! RE-defining the norep vector
        do j = 1, N
           norep(j) = 0
        end do
        
        new_av_started = .false.

     end if

     ! Se nao houver nenhum ativo, testa se uma nova avalanche (para eletrodos)
     ! se iniciou.
     if ( .not. new_av_started ) then
        new_av_started = .not. Teste_Are_ElecDead
     end if


     ! Testa se a rede jah morreu
     ! Testamos se a rede morre apenas no final porque pode ser que a rede morra sem
     ! que os eletrodos sejam excitados em suas ultimas excitacoes.
     if ( isDead(N, netstate_new) ) then

        return

     end if

     ! Em caso de continuacao do algoritmo, atualizando o estado da rede
     netstate = netstate_new

  end do

  close(27)

  return
end subroutine sizeavalanche_loose




!
! Update
!
! Salva a atividade em um arquivo externo
!
subroutine returnmap(N, Madj, netstate, mat_adj, m, h, Nsteps, RM, jindx)

  integer, intent(in)    :: N, Nsteps, m, Madj, jindx
  real(8), intent(in)    :: h, mat_adj(N, 0:Madj)
  integer                :: contador, aux, NN, j, jj, p, progressbar, netstate_new(N), last_spktm
!  real(8)                :: prob, auxp
  integer, intent(inout) :: netstate(N), RM(0:Nsteps)
  real(8)                :: getActivity, ran1

! Como o algoritmo eh sincrono, entao eh necessario
! guardar os estados atualizados em um novo vetor.
  netstate_new = netstate

  last_spktm = 0
  if ( netstate(jindx) .ne. 1 ) then
     netstate(jindx) = 1
  end if
  
  nisis = 0

! Primeira marcacao de atividade
  write(1,*) 0, getActivity(N, netstate)

! Setando o progressbar
  contador = 0
  progressbar = Nsteps/30

  write(*,"(A)",advance='no') " >> ["

  do j = 1, Nsteps

     do p = 1, N


        ! Neuronio pode ser excitado
        if ( ( netstate(p) .eq. 0 ) .and. ( netstate_new(p) .eq. 0 ) ) then
           
           ! O estimulo externo pode excitar o neuronio
           if ( ran1() .lt. h ) then
              netstate_new(p) = 1
           end if
           
        else if ( netstate(p) .eq.  1 ) then
           
           netstate_new(p) = 2
                      
           ! Para nao ficar acessando desnecessariamente mat_adj(p, 0)
           aux = int( mat_adj(p, 0) )
           
           do jj = 1, aux

              idx = int( mat_adj(p,jj) )
              
              if ( ( netstate( idx ) + netstate_new( idx ) .eq. 0 ) .and. ( ran1() .lt. mod( mat_adj(p,jj), 1d0) ) ) then
                 netstate_new( idx ) = 1
              end if

!              write(*,*) jj, aux,  mod( mat_adj(p,jj), 1d0)
!              read(*,*)

           end do

        else
           
           netstate_new(p) = modulo(netstate(p) + 1, m)
           
        end if
        
        
     end do

     ! Atualizando o estado da rede
     netstate = netstate_new

     CALL progressbar_update(contador, progressbar)

     if ( netstate(jindx) .eq. 1 ) then

        RM(nisis) = j - last_spktm
        nisis = nisis + 1
        
        last_spktm = j

     end if

  end do

! Finalizando o contador
  write(*,"(A)") "]"

end subroutine returnmap






!
! phasecode
!
! Funcao que retorna a atividade media apos um periodo dado
! de iteracoes
!
real(8) function phasecode(N, Madj, netstate, mat_adj, phi, m, h, Nsteps) result(media)

  integer, intent(in)    :: N, Nsteps, m, Madj
  real(8), intent(in)    :: h, mat_adj(N, 0:Madj)  ! Matriz de adjacencias
  integer, intent(inout) :: netstate(N)            ! Vetor de estados
  integer, intent(inout) :: phi(N)                 ! Vetor de fases
  integer                :: netstate_new(N)

! Indices
  integer                :: j, jj, NN, p, aux, idx
  integer                :: contador, progressbar, chave

! Funcoes externas
  real(8)                :: getActivity, ran1, syncOverlap
  logical                :: isDead


! Como o algoritmo eh sincrono, entao eh necessario
! guardar os estados atualizados em um novo vetor.
  netstate_new = netstate

! Inicializando a media
  media = 0d0
  
  do j = 1, Nsteps
     
     do p = 1, N
        
        ! Neuronio pode ser excitado
        if ( ( netstate(p) .eq. 0 ) .and. ( netstate_new(p) .eq. 0 ) .and. p .lt. 11 ) then
           
           ! O estimulo externo pode excitar o neuronio
           if ( ran1() .lt. h ) then
              netstate_new(p) = 1
           end if
           
        else if ( netstate(p) .eq.  1 ) then
           
           ! This was wrong: it don't catch the case where m = 2
           ! netstate_new(p) = 2
           netstate_new(p) = modulo(netstate(p) + 1, m)
                      
           ! Para nao ficar acessando desnecessariamente mat_adj(p, 0)
           aux = int( mat_adj(p, 0) )
           
           do jj = 1, aux

              idx = int( mat_adj(p,jj) )
              
              if ( ( netstate( idx ) + netstate_new( idx ) .eq. 0 ) .and. ( ran1() .lt. mod( mat_adj(p,jj), 1d0) ) ) then
                 netstate_new( idx ) = 1
              end if

!              write(*,*) jj, aux,  mod( mat_adj(p,jj), 1d0)
!              read(*,*)
              
           end do

        else
           
           netstate_new(p) = modulo(netstate(p) + 1, m)
           
        end if
        
     end do

     ! Se a rede estiver morta e não houver estimulo externo, não adianta prosseguir.
     if ( isDead(N, netstate_new) .and. ( h .eq. 0.d0 ) ) then
        
        media = media/float(Nsteps)

        return
     end if

     ! Atualizando o estado da rede
     netstate = netstate_new

     ! Imprimindo em arquivo a historia
     media = media + syncOverlap(N, netstate, phi)
     !write(73,*) syncOverlap(N, netstate, phi)
     
     if ( .true. ) then
        open(unit=27, file="activity.data", access="append")
        write(27,*) 'Iteracao: ', j
        do jj = 1, N
           if ( netstate(jj) .eq. 1 ) then
              write(27,*) jj
           endif
        enddo
        write(27,*) 'fim'
        write(27,*) 

        close(27)
     end if
     
     
  end do
  
  media = media/dfloat(Nsteps)
  
  return
end function phasecode
