subroutine setup_network(Nneurons, netstate, frac)

  integer, intent(in)    :: Nneurons
  real(8), intent(in)    :: frac
  integer                :: netstate(Nneurons)
  integer                :: contador, fraction, aux

! Funcoes
  real(8)                :: ran1

  contador = 0

  fraction = Nneurons*frac

  do j = 1, Nneurons
     netstate(j) = 0
  end do

  do while (contador .lt. fraction)

     aux = int( Nneurons*ran1() + 1)

     if ( netstate( aux ) .eq. 0 ) then

        netstate(aux) = 1
        contador = contador + 1

     end if
  end do

  return
end subroutine setup_network





subroutine setup_network_hubs(Nneurons, netstate, frac_to_excite, frac_of_network)

  integer, intent(in)    :: Nneurons, frac_of_network
  real(8), intent(in)    :: frac_to_excite
  integer, intent(inout) :: netstate(Nneurons)
  integer                :: contador, fte, fon, aux

! Funcoes
  real(8)                :: ran1

  contador = 0

  fte = int(Nneurons*frac_to_excite)

  do j = 1, Nneurons
     netstate(j) = 0
  end do

  do while (contador .lt. fte)

     aux = int( frac_of_network*ran1() + 1 )

     if ( netstate( aux ) .eq. 0 ) then

        netstate(aux) = 1
        contador = contador + 1

     end if
  end do

  return
end subroutine setup_network_hubs




!
! Verifica se a parte da rede sendo provada está morta.
!
! Mais rapido que getactivity em caso medio, compensa se nao for necessario
! o valor da atividade (c.c., usar getActivity).
! *  V( EI(j) ) significa a leitura do j-esimon eletrodo.
!
subroutine SelectNeuronsToProbe(N, M, EI)

  integer, intent(in)    :: N, M
  integer, intent(inout) :: EI(M)
  integer                :: j, idx, jj
  logical                :: noh_novo
  real(8)                :: ran1

  if ( M .eq. N ) then
     do j = 1, M
        EI(j) = j
     end do     
  else 
     do j = 1, M
        EI(j) = 0
     end do
     
     j = 1
     
     do while ( j .le. M )
        
        idx = int( ( N  )*ran1() + 1)
        
        noh_novo = .true.
        
        do jj = 1, j
           if ( EI(jj) .eq. idx ) then
              noh_novo = .false.
           end if
        end do
        
        if ( noh_novo ) then
           
           EI(j) = idx
           
           j = j + 1
           
        end if
        
     end do
  end if

  return
end subroutine SelectNeuronsToProbe





subroutine exciteOneNeuron(N, v)

  integer, intent(in)    :: N
  integer, intent(inout) :: v(N)
  integer                :: j, aux

! Funcoes
  real(8)                :: ran1

  do j = 1, N
     v(j) = 0
  end do

  aux = int( N*ran1() + 1)

  v(aux) = 1

  return
end subroutine exciteOneNeuron



subroutine exciteOneNeuron_subsampling(N, V, M, EI)

  integer, intent(in)    :: N, M, EI(M)
  integer, intent(inout) :: V(N)
  integer                :: j, aux

! Funcoes
  real(8)                :: ran1

  do j = 1, N
     V(j) = 0
  end do

  aux = int( M*ran1() + 1)

  V( EI(aux) ) = 1

  return
end subroutine exciteOneNeuron_subsampling




subroutine exciteOneNeuronNoRepeat(Nneurons, netstate, norep)

  integer, intent(in)    :: Nneurons
  integer, intent(inout) :: netstate(Nneurons), norep(Nneurons)
  integer                :: contador, fraction, aux

! Funcoes
  real(8)                :: ran1

  do j = 1, Nneurons
     netstate(j) = 0
     norep(j) = 0
  end do

  aux = int( Nneurons*ran1() + 1)

  netstate(aux) = 1

  return
end subroutine exciteOneNeuronNoRepeat



!
! Selection of the Network
!
! Type: integer selecting the random graph to be constructed.
!
! -1: Mean-field Erdos-Renyi
! 1 : Erdos-Renyi
! 2 : Cayley w free group
! 
!
subroutine connect_network_selection(N, K, sigma, mat_adj, M, type, delta)
  
  integer, intent(in)        :: N, type, M
  real(8), intent(in)        :: sigma, K, delta
  real(8), intent(inout)     :: mat_adj(N, 0:M)
  
  if ( type .eq. 1 ) then
     
     ! Erdos-Renyi with fixed number of connections
     write(*,*) 'Erdos-Renyi'
     call connect_network(N, K, sigma, mat_adj)
     
  else if ( type .eq. 2 ) then
     
     ! Cayley Tree of free group with K+1 generators
     call Cayley_Tree(N, K, sigma, mat_adj)
     
  else if (type .eq. 3) then
     
     ! BA random graph - power law with exponent = 3
     write(*,*) '    |->> Construindo Barabasi-Albert '
     call BA_algth(N, K, sigma, mat_adj, M)
     
  else if (type .eq. 4) then
     
     ! Exponential random graph
     call EXP_algth(N, K, sigma, mat_adj, M)
     
  else if (type .eq. 5) then
     
     ! GBA random graph - power law with exponent = delta
     write(*,*) '    |->> Construindo Generalized Barabasi-Albert '
     call GBA_algth(N, K, sigma, mat_adj, M, delta)
     
  else if (type .eq. 6) then
     
     ! Watts Strogatz
     call WS_algth(N, K, sigma, mat_adj, M)
     
  else if (type .eq. 7) then
     
     ! Scale free Uncorrelated CM
     write(*,*) 'Rede UCM chamada!'
     call SFUCM_algth(N, K, sigma, mat_adj, M, delta)
     
     
  else if ( type .eq. -1 ) then
     
     ! Erdos-Renyi mean-field
     call ER_meanfield(N, K, sigma, mat_adj, M)
     
  else if ( type .eq. 8 ) then
     
     ! Square lattice
     call square_lattice(N, sigma, mat_adj, M)
     
  end if
  
  return
end subroutine connect_network_selection



!
! Erdos_Renyi
!
subroutine connect_network(N, K, sigma, mat_adj)

  integer, intent(in)    :: N
  real(8), intent(in)    :: sigma, K
  integer                :: contador, aux, NN, j, p, pp, M, progressbar
  real(8)                :: prob, auxp, pmax
  real(8), intent(inout) :: mat_adj(N, 0:5*int(K))
  logical                :: isntYetNeighbor

! Funcoes
  real(8)                :: ran1

!  CALL RANDOM_SEED (SIZE=pr)
!  seed = sizer + 29 * (/ (i - 1, i = 1, pr) /)
!  CALL RANDOM_SEED (PUT=seed)

  M = 5*int(K)

! Numero de possiveis vizinhos de um dado noh
  NN = N*int(K)/2

! Calculando pmax
  pmax = 2*sigma/K

! Inicia o progress bar!
  write(*,"(A)",advance='no') " >> ["

! Inicializa os elementos 0 da matriz de adjacencias
  do j = 1, N
     mat_adj(j, 0) = 0d0

     do p = 1, M
        mat_adj(j, p) = 0d0
     end do

  end do
  write(*,"(A)",advance="no") "="

! Inicializando o valor de j para o proximo loop
  j = 1

! Inicializando o contador e o progressbar
  contador = 0
  progressbar = N/1
  
! Algoritmo para gerar o Erdos-Renyie
  do while ( j .le. NN )

     pp = int( N*ran1() + 1 )
     p = int( N*ran1() + 1 )

!     write(*,*) pp, p

     if ( ( p .ne. pp ) .and. isntYetNeighbor(pp, p, N, 5*int(K), mat_adj) ) then
        
        ! Contabiliza um noh a mais para j e p
        mat_adj(pp, 0) = mat_adj(pp, 0) + 1d0
        mat_adj(p, 0) = mat_adj(p, 0) + 1d0
        
        ! Marca o peso da adjacencia usando a seguinte
        ! indexcao: X . Y (parte inteira . parte decimal)
        ! X = indice do vizinho
        ! Y = peso da ligacao
        
        ! Uniform distribution
        auxp = pmax*ran1()
        
        ! All weights exactly equal
        !auxp = sigma/K
        
        mat_adj( p , int(mat_adj(p,0)) ) = auxp + dfloat(pp)
        mat_adj( pp , int(mat_adj(pp,0)) ) = auxp + dfloat(p)
        
! Debuggin...
!           write(*,*) mod(p + auxp,1.0), auxp + p


        ! Iterando a variavel do loop.        
        j = j + 1

        ! Atualizando a barra de progresso.
        CALL progressbar_update(contador, progressbar)
        
     end if
     
  end do
  
  write(*,"(A)") "] "
  return
end subroutine connect_network





!
! Square Lattice implementation
!
! The basic idea is to separate in two parts the algorithm
! 1. Connect every node in the first row to its right neighbor.
! 2. For every other node, connect it to its right neighbor and
! is bottom neighor. Last node in each row must be treated as
! a special case.
! This way no double edges are introduced and everything seems to
! run pretty fast.
!
subroutine square_lattice(N, sigma, mat_adj, M)
  
  ! Model variables
  integer, intent(in)    :: N
  real(8), intent(inout) :: mat_adj(N, 0:M)  
  real(8), intent(in)    :: sigma
  
  ! Number of columns in adjacency matrix
  integer, intent(in)    :: M

  
  ! Important variables
  integer                :: K, L
  real(8)                :: prob, auxp, pmax
  integer                :: aux, LIN

  ! Index variables
  integer                :: j, p, pp
  integer                :: linha, from, to

  integer                :: contador, NN

  ! External functions
  real(8)                :: ran1
  logical                :: isntYetNeighbor
  integer                :: progressbar
  
  
! Side of the square
  L = dsqrt( dfloat(N) )
  
  
! Conectividade serah 4
  K = 4

! Last inner neuron index
  LIN = L - 1

! Calculando pmax
  pmax = 2.d0*sigma/K
  write(*,*) 'pmax', pmax

! Inicia o progress bar!
  write(*,"(A)",advance='no') " >> ["

! Inicializa os elementos 0 da matriz de adjacencias
  do j = 1, N
     mat_adj(j, 0) = 0d0

     do p = 1, M
        mat_adj(j, p) = 0d0
     end do

  end do
  
  
! Show that the adjacency matrix has being initialized.
  write(*,"(A)",advance="no") "="  
  
  
! Inicializando o contador e o progressbar
  contador = 0
  progressbar = N/1
  
  
  ! Let's start connecting the first row
  ! of neurons only to their right neighbors
  do j = 1, LIN
     
     from = j
     
     ! ===
     ! Right neighbor
     to = from + 1
     
     auxp = pmax*ran1()
     
     mat_adj(from,0) = mat_adj(from,0) + 1.d0
     mat_adj(to,0) = mat_adj(to,0) + 1.d0
     
     mat_adj( from , int(mat_adj(from,0)) ) = auxp + dfloat(to)
     mat_adj( to , int(mat_adj(to,0)) ) = auxp + dfloat(from)
     
  enddo
  
  
  
  ! Now the rest
  do j = 2, L
     
     linha = L*(j-1)
     
     do p = 1, LIN
        
        from = linha + p
        
        ! ===
        ! Bottom neighbour
        to = from - L
        auxp = pmax*ran1()
        
        mat_adj(from,0) = mat_adj(from,0) + 1.d0
        mat_adj(to,0) = mat_adj(to,0) + 1.d0
        
        mat_adj( from , int(mat_adj(from,0)) ) = auxp + dfloat(to)
        mat_adj( to , int(mat_adj(to,0)) ) = auxp + dfloat(from)
        
        !===
        ! Right neighbour
        to = from + 1
        auxp = pmax*ran1()
        
        mat_adj(from,0) = mat_adj(from,0) + 1.d0
        mat_adj(to,0) = mat_adj(to,0) + 1.d0
        
        mat_adj( from , int(mat_adj(from,0)) ) = auxp + dfloat(to)
        mat_adj( to , int(mat_adj(to,0)) ) = auxp + dfloat(from)
        
     enddo
     
     ! Connecting the last of a row to
     ! its upper neighbor
     from = from + 1
     to = from - L
     auxp = pmax*ran1()
     mat_adj(from,0) = mat_adj(from,0) + 1.d0
     mat_adj(to,0) = mat_adj(to,0) + 1.d0
     mat_adj( from , int(mat_adj(from,0)) ) = auxp + dfloat(to)
     mat_adj( to , int(mat_adj(to,0)) ) = auxp + dfloat(from)
     
  enddo
  
  
  write(*,"(A)") "] "
  return
end subroutine square_lattice





!
! Watts Strogatz
!
! Detalhes sobre a implementacao:
! Usualmente, descreve-se o algoritmo como partindo de uma rede regular de 
! z = K/2 vizinhos. No entanto, se isso for feito, a matriz de adjacencia pode
! ficar com buracos ao retirar ligacoes para realoca-las de forma nao local.
! Ou isso, ou nos teriamos que apagar a ligacao e recolocar as outras todas
! uma posicao para tras. Uma opcao bem melhor que isso eh simplesmente gerar
! apenas conexoes que vao se manter. Entao, seguirei a seguinte ideia:
! -- para todo ij satisfazendo |i-j| < z, verifico se vai ter ligacao
! -- se tiver, escrevo na matriz de adjacencia
! -- se nao, jah procuro quem serah o vizinho nao local
! Tem um problema com esta logica: o ideal eh que a geracao seguinte jah soubesse
! que ha ligacoes na rede regular. Entao....... PENSAR MAIS.
!
!
subroutine WS_algth(N, K, sigma, mat_adj, M)

  integer, intent(in)    :: N
  real(8), intent(in)    :: sigma, K
  integer                :: aux, NN, j, p, pp, M, z
  real(8)                :: prob, auxp, pmax
  real(8), intent(inout) :: mat_adj(N, 0:M)
  logical                :: isntYetNeighbor

  integer                :: progressbar, contador

! Funcoes
  real(8)                :: ran1


! Calculando pmax
  pmax = 2*sigma/K

! Variavel auxiliar
  NN = N - 1

! Inicia o progress bar!
  write(*,"(A)",advance='no') " >> ["

! Inicializa os elementos 0 da matriz de adjacencias
  do j = 1, N
     
     mat_adj(j, 0) = 0d0

     do p = 1, M
        mat_adj(j, p) = 0d0
     end do

  end do
  write(*,"(A)",advance="no") "="


! Numero de vizinhos aa direita
  z = int(K/2.d0)

  write(*,"(A)",advance="no") "="


! Inicializando o contador e o progressbar
  contador = 0
  progressbar = N/15

! Loop para ir gerando a rede
  do j = 1, N
     do jj = 1, NN
        do jind = 1, z
     
           if (ran1() .lt. beta) then

              ! Primeiramente vamos apagar esta conexão
              p = int( N*ran1() + 1 )
              
              if ( ( j .ne. p ) .and. isntYetNeighbor(j, p, N, M, mat_adj) ) then
                 
                 ! Contabiliza um noh a mais para j e p
                 mat_adj(j, 0) = mat_adj(j, 0) + 1d0
                 mat_adj(p, 0) = mat_adj(p, 0) + 1d0
                 
                 ! Settando os pesos
                 auxp = pmax*ran1()
                 mat_adj( j , int(mat_adj(j,0)) ) = auxp + dfloat(p)
                 mat_adj( p , int(mat_adj(p,0)) ) = auxp + dfloat(j)
                 
! Debuggin...
!           write(*,*) mod(p + auxp,1.0), auxp + p
                 
                 ! Atualizando a barra de progresso.
                 CALL progressbar_update(contador, progressbar)
                 
              end if
              
           else
              
              ! Contabiliza um noh a mais para j e p
              mat_adj(jj, 0) = mat_adj(jj, 0) + 1d0
              mat_adj(j, 0) = mat_adj(j, 0) + 1d0
              
              ! Settando os pesos
              auxp = pmax*ran1()
              mat_adj( j , int(mat_adj(j,0)) ) = auxp + dfloat(jj)
              mat_adj( jj , int(mat_adj(jj,0)) ) = auxp + dfloat(j)
              
           end if

        end do
     end do
  end do
  
  write(*,"(A)") "] "
  return
end subroutine WS_algth






!
!  Barabasi-Albert
!
!
subroutine BA_algth(N, K, sigma, mat_adj, Mad)

  integer, intent(in)    :: N, Mad
  real(8), intent(in)    :: K, sigma
  real(8), intent(inout) :: mat_adj(N, 0:Mad)

  integer                :: NN, progressbar
  integer                :: m, m0, m0aux, sorteado
  integer                :: SOMA_DEGREES, SOMA_DEGREES_aux, no_conn
  real(8)                :: prob, auxp, CORR, pjj, pmax

! Variaveis para gerar a distribuicao gaussiana truncada
  real(8)                :: tg_z, tg_sigma, tg_mu, tg_E0, tg_E1, tg_pi
  integer                :: tg_j


! Indices disponiveis
  integer                :: contador, aux, j, jj, jjj
  logical                :: continua

! Funcoes
  real(8)                :: ran1, inverf, dinvnorm
  logical                :: isntYetNeighbor



! Numero de possiveis vizinhos de um dado noh
  NN = N - 1

! Correcao para construcao dos pesos
!  CORR = 2*sigma/K - 1
  pmax = 2*sigma/K
  
! Inicia o progress bar!
  write(*,"(A)",advance='no') " >> ["
  
! Inicializa os elementos 0 da matriz de adjacencias
  do j = 1, N
     do jj = 0, Mad
        mat_adj(j, jj) = 0d0
     end do
  end do

  write(*,"(A)",advance="no") "="

! Inicializando o contador e o progressbar
  contador = 0
  progressbar = N/30

! Numero de vizinhos iniciais e numero de vizinhancas
! a cada passo temporal.
  m = int(K/2)
  m0 = int(K/2) + 5
  
! Condicao inicial do algor. BA: comecaremos com dois nohs,
! ambos conectados um ao outro.
  m0aux = m0 - 1
  do j = 1, m0
     mat_adj(j, 0) = m0aux
  enddo

! Inserindo os pesos dessas duas conexoes iniciais
  do j = 1, m0aux

     ! Inicializando o indice
     jjj = j

     do jj = j + 1, m0

        ! calculo do peso
        auxp = pmax*ran1()

        ! Associacao das ligacoes
        mat_adj(j, jjj) = jj + auxp
        mat_adj(jj, j) = j + auxp

        jjj = jjj + 1

     enddo
  enddo

!  do j = 1, m0
!     do jj = 1, m0aux
!        write(*,*) mat_adj(j,jj)
!     end do
!  end do
!  read(*,*)

! Soma dos graus
  SOMA_DEGREES = m0*m0aux

! Zerando variavel auxiliar: esta variavel serve
! para nao precisar ficar calculando a soma dos graus
! a cada iteracao
  SOMA_DEGREES_aux = 0


! Algoritmo para gerar o Erdos-Renyie
  do j = m0 + 1, N
     
     ! Auxilia na contagem de conexoes
     no_conn = 0
     
     do while ( no_conn .lt. m )
        
        ! Faz o sorteio do noh que deverah ser ligado
        sorteado = int( SOMA_DEGREES*ran1() + 1 )
        
        ! Iniciando as variaveis do loop seguinte
        continua = .TRUE.
        jj = 1
        pjj = 0.d0
 
        ! Busca o indice do sorteado
        do while ( continua )
           
           ! Calculando a probabilidade de o novo noh j se conectar
           ! com o noh jj
           pjj = mat_adj(jj, 0) + pjj

           if ( sorteado .le. pjj ) then
              continua = .FALSE.
           else if ( jj .lt. j ) then 
              jj = jj + 1
           else
              jj = j
              continua = .FALSE.
           end if
           
        enddo

        
        ! Removendo nohs nao elegiveis (duplas conexoes e mesmos nohs)
        if ( isntYetNeighbor(jj, j, N, Mad, mat_adj) .and. ( jj .ne. j ) ) then
           
           ! Contabiliza um noh a mais para j e p
           mat_adj(jj, 0) = mat_adj(jj, 0) + 1d0
           mat_adj(j, 0) = mat_adj(j, 0) + 1d0
           
           !           write(*,*) mat_adj(jj, 0), mat_adj(j+1,0)          
           
           ! Atualizando a normalizacao do preferential_probability
           SOMA_DEGREES = SOMA_DEGREES + 2
           
           ! Marca o peso da adjacencia usando a seguinte
           ! indexcao: X . Y (parte inteira . parte decimal)
           ! X = indice do vizinho
           ! Y = peso da ligacao
           
           ! Drawing weight from uniform distribution
           auxp = pmax*ran1()
           
           ! Drawing from a "pseudo-gaussian"
           !tg_mu = sigma/K
           !tg_sigma = 0.1d0
           !tg_E0 = erf( ( 0 - tg_mu )/sqrt(2.d0)/tg_sigma )
           !tg_E1 = erf( ( 1 - tg_mu )/sqrt(2.d0)/tg_sigma )
           !tg_z = 0.5d0*( 1.d0 + tg_E0 + ran1()*( tg_E1 - tg_E0  ) )
           !pi = 3.14159265d0
           !auxp = inverf( (tg_z) )
           !write(*,*) tg_z, auxp, erf(auxp), tg_z - erf(auxp)
           !
           !write(*,*) dinvnorm( erf(.1d0) )
           !read(*,*)
           
           ! Drawing from a Beta distribution
           
              
           
           ! Fixed degree weights
           !auxp = sigma/K
           
           mat_adj( j , int(mat_adj(j, 0)) ) = auxp + int(jj)
           mat_adj( jj , int(mat_adj(jj, 0)) ) = auxp + int(j)

!           if ( j .eq. 0 .or. jj .eq. 0) then 
!             write(*,*) "chicken basket"
!           end if
           no_conn = no_conn + 1
           
        end if
        
     enddo

     CALL progressbar_update(contador, progressbar)

  enddo

  write(*,"(A)") "] "
  return
end subroutine BA_algth





!
!  Exponential random graph
!
subroutine EXP_algth(N, K, sigma, mat_adj, Mad)

  integer, intent(in)    :: N, Mad
  real(8), intent(in)    :: K, sigma
  integer                :: NN, progressbar
  integer                :: m, m0, no_conn
  real(8)                :: prob, auxp, CORR, pjj, pmax
  real(8), intent(inout) :: mat_adj(N, 0:Mad)

! Indices disponiveis
  integer                :: contador, aux, j, jj, jjj

! Funcoes
  real(8)                :: ran1
  logical                :: isntYetNeighbor



! Numero de possiveis vizinhos de um dado noh
  NN = N - 1

! Correcao para construcao dos pesos
!  CORR = 2*sigma/K - 1
  pmax = 2*sigma/K

! Inicia o progress bar!
  write(*,"(A)",advance='no') " >> ["

! Inicializa os elementos 0 da matriz de adjacencias
  do j = 1, N

     mat_adj(j, 0) = 0d0

     do jj = 1, Mad
        mat_adj(j, jj) = 0d0
     end do

  end do

  write(*,"(A)",advance="no") "="

! Inicializando o contador e o progressbar
  contador = 0
  progressbar = N/15

! Numero de vizinhos iniciais e numero de vizinhancas
! a cada passo temporal.
  m = int(K/2)
  m0 = int(K/2)
  

! Algoritmo para gerar o grafo exponencial
  do j = m0 + 1, N

     ! Auxilia na contagem de conexoes
     no_conn = 0

     do while ( no_conn .lt. m )

        ! Sorteio uniforme de dois nos
        jj = int(j*ran1()) + 1
        jjj = int(j*ran1()) + 1


        ! Verificando se jah estao conectados
        if ( isntYetNeighbor(jj, jjj, N, M, mat_adj) .and. ( jj .ne. jjj ) ) then

           ! Contabiliza um noh a mais para j e p
           mat_adj(jj, 0) = mat_adj(jj, 0) + 1.d0
           mat_adj(jjj, 0) = mat_adj(jjj, 0) + 1.d0
                      
           ! Marca o peso da adjacencia usando a seguinte
           ! indexcao: X . Y (parte inteira . parte decimal)
           ! X = indice do vizinho
           ! Y = peso da ligacao
           
           auxp = pmax*ran1()
           mat_adj( jj , int(mat_adj(jj, 0)) ) = auxp + dfloat(jjj)
           mat_adj( jjj , int(mat_adj(jjj, 0)) ) = auxp + dfloat(jj)

           no_conn = no_conn + 1
           
        end if

     enddo

     CALL progressbar_update(contador, progressbar)

  end do
  
  write(*,"(A)") "] "
  return
end subroutine EXP_algth





!
!  Barabasi-Albert
!
!
subroutine GBA_algth(N, K, sigma, mat_adj, Mad, delta)
  
  integer, intent(in)    :: N, Mad
  real(8), intent(in)    :: K, sigma, delta
  integer                :: NN, progressbar
  integer                :: m, m0, m0aux, c
  integer                :: SOMA_DEGREES, SOMA_DEGREES_aux, no_conn
  real(8)                :: prob, auxp, CORR, pjj, pmax, Attract
  real(8), intent(inout) :: mat_adj(N, 0:Mad)
  
! Indices disponiveis
  integer                :: contador, aux, j, jj, jjj, sorteado
  logical                :: continua
  
! Funcoes
  real(8)                :: ran1
  logical                :: isntYetNeighbor
  
  
!  write(*,*)
!  write(*,*) "Greetings from Dorogovtsev-Mendes generator"
!  write(*,*)
  
  
! Numero de possiveis vizinhos de um dado noh
  NN = N - 1
  
! Correcao para construcao dos pesos
!  CORR = 2*sigma/ - 1
  pmax = 2*sigma/k
  
! Inicia o progress bar!
  write(*,"(A)",advance='no') " >> ["
  
! Inicializa os elementos 0 da matriz de adjacencias
  do j = 1, N
     
     mat_adj(j, 0) = 0d0
     
     do jj = 1, Mad
        mat_adj(j, jj) = 0d0
     end do
     
  end do
  
  write(*,"(A)",advance="no") "="
  
! Inicializando o contador e o progressbar
  contador = 0
  progressbar = N/15
  
! Numero de vizinhos iniciais e numero de vizinhancas
! a cada passo temporal.
  m = int(K/2.d0)
  m0 = int(K/2) + 5
  
  
! Verifying attractiveness
  Attract = -(delta + 3.d0)*dfloat(m)
  if ( Attract .lt. -m ) then
     write(*,*) "ERROR - Atrractiveness out of range!"
     stop
  end if
  
! Condicao inicial do algor. BA: comecaremos com dois nohs,
! ambos conectados um ao outro.
  m0aux = m0 - 1
  do j = 1, m0
     mat_adj(j, 0) = dfloat(m0aux)
  enddo


! Inserindo os pesos dessas duas conexoes iniciais
  do j = 1, m0aux

     ! Inicializando o indice
     jjj = j

     do jj = j + 1, m0

        ! calculo do peso
        auxp = pmax*ran1()

        ! Associacao das ligacoes
        mat_adj(j, jjj) = jj + auxp
        mat_adj(jj, j) = j + auxp

        jjj = jjj + 1

     enddo
  enddo

  
  
! Soma dos graus
  SOMA_DEGREES = m0*m0aux
  
! Zerando variavel auxiliar: esta variavel serve
! para nao precisar ficar calculando a soma dos graus
! a cada iteracao
  SOMA_DEGREES_aux = 0
  
  
! Algoritmo para gerar o Erdos-Renyie
  do j = m0, NN
     
     ! Auxilia na contagem de conexoes
     no_conn = 0
     
     do while ( no_conn .lt. m )
        
        ! Faz o sorteio do noh que deverah ser ligado
        sorteado = int( (SOMA_DEGREES + j*Attract)*ran1() + 1)
        
        ! Iniciando as variaveis do loop seguinte
        continua = .TRUE.
        jj = 1
        pjj = 0.d0
        
        ! Busca o indice do sorteado
        do while ( continua )
           
           ! Calculando a probabilidade de o novo noh j se conectar
           ! com o noh jj
           pjj = mat_adj(jj, 0) + Attract + pjj
           
           if ( sorteado .le. pjj ) then
              continua = .FALSE.
           else if ( jj .lt. j ) then 
              jj = jj + 1
           else
              jj = j
              continua = .FALSE.
           end if
           
        enddo
        
        
        ! Verificando se conecta
        if ( ( ii .ne. jjj ) .and. isntYetNeighbor(jj, jjj, N, Mad, mat_adj) ) then

           jjj = j + 1
           
           ! Contabiliza um noh a mais para j e p
           mat_adj(jj, 0) = mat_adj(jj, 0) + 1d0
           mat_adj(jjj, 0) = mat_adj(jjj, 0) + 1d0
           
           !           write(*,*) mat_adj(jj, 0), mat_adj(j+1,0)          
           
           ! Atualizando a normalizacao do preferential_probability
           SOMA_DEGREES = SOMA_DEGREES + 2
           
           ! Marca o peso da adjacencia usando a seguinte
           ! indexcao: X . Y (parte inteira . parte decimal)
           ! X = indice do vizinho
           ! Y = peso da ligacao
           
           auxp = pmax*ran1()
           mat_adj( jjj , int(mat_adj(jjj, 0)) ) = auxp + int(jj)
           mat_adj( jj , int(mat_adj(jj, 0)) ) = auxp + int(jjj)
           
           
           if ( mat_adj( jjj , int(mat_adj(jjj, 0)) ) .eq. 0.d0 ) then
              write(*,*) 'oi?????????????????????????'
           end if
           if ( mat_adj( jj , int(mat_adj(jj, 0)) ) .eq. 0.d0 ) then
              write(*,*) 'oi?????????????????????????'
           end if

           
           no_conn = no_conn + 1
           
        end if

     enddo

     CALL progressbar_update(contador, progressbar)

!     SOMA_DEGREES = SOMA_DEGREES + SOMA_DEGREES_aux

!     SOMA_DEGREES_aux = 0

  end do
  
  
  write(*,"(A)") "] "
  return
end subroutine GBA_algth






!
! Undirected Cayley Tree construction
! 
subroutine Cayley_tree(N, K, sigma, mat_adj)

  integer, intent(in)    :: N
  real(8), intent(in)    :: sigma, K
  integer                :: contador, aux, auxj, NN, j, p, pp, M, progressbar
  real(8)                :: prob, auxp, pmax
  real(8), intent(inout) :: mat_adj(N, 0:int(K)+1)
  logical                :: isntYetNeighbor

! Funcoes
  real(8)                :: ran1


! Numero de possiveis vizinhos de um dado noh
  NN = N*int(K)/2
  M = int(K) + 1

! Calculando pmax
  pmax = 2*sigma/K

! Inicia o progress bar!
  write(*,"(A)",advance='no') " >> ["

! Inicializa os elementos 0 da matriz de adjacencias
  do j = 1, N
     mat_adj(j, 0) = K + 1

     do p = 1, M
        mat_adj(j, p) = 0d0
     end do

  end do
  write(*,"(A)",advance="no") "="

! Inicializando o contador e o progressbar
  contador = 0
  progressbar = N/15


! Fazendo o primeiro noh na mao
  do aux = 1, int(K)
     auxp = pmax*ran1()

     mat_adj(1, aux) = aux + 1 + auxp
     mat_adj(aux, 1) = 1 + auxp
  enddo


! Inicializando o valor de j para o proximo loop
  j = K + 2

! Algoritmo para ir construindo as arvores
  do aux = 2, N

     ! Escolhendo os seus K vizinhos
     do auxj = 1, int(K)

        ! Indice em que pararam as atribuicoes
        j = j + 1
        
        ! Gerando o peso da vizinhanca
        auxp = pmax*ran1()
        
        ! Atribuindo a vizinhanda
        mat_adj(aux, j) = j + auxp
        mat_adj(j, 1) = aux + auxp

     enddo

        
! Debuggin...
!           write(*,*) mod(p + auxp,1.0), auxp + p

     CALL progressbar_update(contador, progressbar)
     
  end do
  
  write(*,"(A)") "] "
  return
end subroutine Cayley_Tree





!
! Scale-Free Uncorrelated Configuration Model
!
!
subroutine SFUCM_algth(N, K, sigma, mat_adj, Mad, delta)

  integer, intent(in)    :: N, Mad
  real(8), intent(in)    :: sigma, K
  integer                :: contador, aux, NN, j, p, pp, M, progressbar, listaKis(0:N), NK
  integer                :: r_Zipf, min_Zipf, min_Zp, conta_NK
  real(8)                :: prob, auxp, pmax, normaZipf, delta
  real(8), intent(inout) :: mat_adj(N, 0:Mad)
  logical                :: isntYetNeighbor
  real(8), allocatable   :: CDF(:)
  
! Indices utilizados para construcao das coisas
  integer                :: pc, tp, tpp
  
! Funcoes
  real(8)                :: ran1
  

  min_Zipf = int( K/2.d0 )
  min_Zp = min_Zipf + 1
  
! Numero de possiveis vizinhos de um dado noh
  NN = int( sqrt( dfloat(N) ) )
  
! Definindo a distribuicao cumulativa da distribuicao alvo
  allocate( CDF(min_Zipf:NN) )
  
! Condicoes iniciais
  normaZipf = dfloat(min_Zipf)**(-delta)
  CDF(min_Zipf) = normaZipf
  
  do j = min_Zp, NN
     auxp = j**(-delta)
     normaZipf = normaZipf + auxp
     CDF(j) = auxp + CDF(j-1)
  enddo
  

! Normalizando 
! Ao somar auxp, estavamos somando as "probabilidades sem estarem normalizadas
! Por isso, precisamos agora renormalizar - na verdade, estamos dividindo
! termo a termo da soma anterior pelo fator apropriado.
  do j = min_Zipf, NN
     CDF(j) = 1 - CDF(j)/normaZipf
  enddo

! ***
! Comment:
! CDF jah foi testada!
! ***


! Gerando a lista de stubs
  do j = 1, N
     listaKis(j) = 0
  enddo  
  
  listaKis(0) = 0
  do j = 1, N
     
     auxp = ran1()
     p = min_Zipf
     r_Zipf = 0.d0
     
     do while ( p .le. NN .and. r_Zipf .eq. 0)
        
        ! Buscando a inversa
        if ( auxp .gt. CDF(p) ) then
           r_Zipf = p
        end if
        
        ! Iterando p
        p = p + 1
        
     enddo
     
     listaKis(j) = r_Zipf
  enddo
  
! Compondo a lista de forma a facilitar a busca
  do j = 1, N
     listaKis(j) = listaKis(j) + listaKis(j-1)
  enddo



! ***
! Comment:
! Distribuicao das listasKis estah ficando uma lei de potencia
! como desejada!!
! ***


! Numero total de conexoes:
  NK = listaKis(N)


! Calculando pmax
  pmax = 2*sigma/K

! Inicia o progress bar!
  write(*,"(A)",advance='no') " >> ["

! Inicializa os elementos 0 da matriz de adjacencias
  do j = 1, N
     mat_adj(j, 0) = 0d0

     do p = 1, Mad
        mat_adj(j, p) = 0d0
     end do

  end do
  write(*,"(A)",advance="no") "="
  
! Inicializando o valor de j para o proximo loop
  j = 1

! Inicializando o contador e o progressbar
  contador = 0
  progressbar = N/1
  
  conta_NK = NK
  if ( MOD(conta_NK, 2) .ne. 0 ) then
     conta_NK = conta_NK - 1
  end if

  do while ( conta_NK .ne. 0 )
     
     ! Estes indices sao com respeito aos open stubs (pontas abertas) 
     ! saindo de cada um dos nohs
     pp = int( conta_NK*ran1() + 1 )
     p = int( conta_NK*ran1() + 1 )
     
     ! Precisamos agora identificar o indice deste noh na matriz de incidencias
     tp = 0
     tpp = 0
     pc = 1
     do while ( pc .le. N .and. ( tpp .eq. 0 .or. tp .eq. 0 ) )
        
        if ( pp .le. listaKis(pc) .and. tpp .eq. 0 ) then
           tpp = pc
        end if
        
        if ( p .le. listaKis(pc) .and. tp .eq. 0 ) then
           tp = pc
        end if
        
        ! Iterando p
        pc = pc + 1
     enddo   
     
     
     if ( ( tp .ne. tpp ) .and. isntYetNeighbor(tpp, tp, N, Mad, mat_adj) ) then
        
        ! Contabiliza um noh a mais para j e p
        mat_adj(tpp, 0) = mat_adj(tpp, 0) + 1d0
        mat_adj(tp, 0) = mat_adj(tp, 0) + 1d0
        
        ! Marca o peso da adjacencia usando a seguinte
        ! indexcao: X . Y (parte inteira . parte decimal)
        ! X = indice do vizinho
        ! Y = peso da ligacao
        
        auxp = pmax*ran1()
        
        if ( ( mat_adj(tp, 0) .ge. Mad ) .or. ( mat_adj(tpp, 0) .ge. Mad ) ) then
           write(*,*) 'Explodiu!'
           return
        end if

        mat_adj( tp , int(mat_adj(tp,0)) ) = auxp + dfloat(tpp)
        mat_adj( tpp , int(mat_adj(tpp,0)) ) = auxp + dfloat(tp)
        
        ! Iterando a variavel do loop.        
        j = j + 1

        ! Atualizando a barra de progresso.
        CALL progressbar_update(contador, progressbar)
        
        ! Apagando essas duas arestas omo soltas
        do i = tp, N
           listaKis(i) = listaKis(i) - 1
        enddo
        do i = tpp, N
           listaKis(i) = listaKis(i) - 1
        enddo        
        conta_NK = conta_NK - 2
        
     end if
     
  end do
  
  write(*,"(A)") "] "
  return
end subroutine SFUCM_algth






subroutine ER_meanfield(N, K, sigma, mat_adj, M)

  integer, intent(in)    :: N, M
  real(8), intent(in)    :: sigma, K
  integer                :: contador, aux, NN, j, p, pp, progressbar
  real(8)                :: prob, aij
  real(8), intent(inout) :: mat_adj(N, 0:M)
  logical                :: isntYetNeighbor

! Funcoes
  real(8)                :: ran1

!  CALL RANDOM_SEED (SIZE=pr)
!  seed = sizer + 29 * (/ (i - 1, i = 1, pr) /)
!  CALL RANDOM_SEED (PUT=seed)

! Numero de possiveis vizinhos de um dado noh
  NN = N*K/2

! Calculando os pesos
  aij = sigma/K

! Inicia o progress bar!
  write(*,"(A)",advance='no') " >> ["

! Inicializa os elementos 0 da matriz de adjacencias
  do j = 1, N
     mat_adj(j, 0) = 0d0

     do p = 1, M
        mat_adj(j, p) = 0d0
     end do

  end do
  write(*,"(A)",advance="no") "="

! Inicializando o valor de j para o proximo loop
  j = 1

! Inicializando o contador e o progressbar
  contador = 0
  progressbar = N/15

! Algoritmo para gerar o Erdos-Renyie
  do while ( j .le. NN )

     pp = int( N*ran1() + 1 )
     p = int( N*ran1() + 1 )

!     write(*,*) pp, p

     if ( ( p .ne. pp ) .and. isntYetNeighbor(pp, p, N, 5*int(K), mat_adj) ) then

        ! Contabiliza um noh a mais para j e p
        mat_adj(pp, 0) = mat_adj(pp, 0) + 1d0
        mat_adj(p, 0) = mat_adj(p, 0) + 1d0
        
        ! Marca o peso da adjacencia usando a seguinte
        ! indexcao: X . Y (parte inteira . parte decimal)
        ! X = indice do vizinho 
        ! Y = peso da ligacao

        mat_adj( p , int(mat_adj(p,0)) ) = aij + dfloat(pp)
        mat_adj( pp , int(mat_adj(pp,0)) ) = aij + dfloat(p)
        
! Debuggin...
!           write(*,*) mod(p + auxp,1.0), auxp + p
        
        j = j+1
        CALL progressbar_update(contador, progressbar)
        
     end if
     
  end do
  
  write(*,"(A)") "] "
  return
end subroutine ER_meanfield





logical function isntYetNeighbor(j, p, N, M, mat_adj) result(isit)

  integer    :: j, p, N
  real(8)    :: mat_adj(N, 0:M)
  integer    :: l, auxl

! Numero de vizinhos a verificar
  auxl = int( mat_adj(j,0) )

  do l = 1, auxl

     if ( int(mat_adj(j,l)) .eq. p ) then

        isit = .FALSE.
        return

     end if

  end do
  
  isit = .TRUE.

  return
end function isntYetNeighbor



subroutine connect_network_m2(N, K, mat_adj, pmax)

  integer, intent(in)    :: N, K
  real(8), intent(in)    :: pmax
  integer                :: contador, aux, NN, j, p, M, progressbar
  real(8)                :: prob, auxp
  real(8), intent(inout) :: mat_adj(N, 0:5*K)

  real(8), external      :: ran1

  M = 5*K

! Numero de possiveis vizinhos de um dado noh
  NN = N - 1

! Probabilidade de dados dois nos, haver um vertice
! que os liga tal que a media seja K
  prob = dfloat(K)/dfloat(N - 1)


! Inicia o progress bar!
  write(*,"(A)",advance='no') " >> ["

! Inicializa os elementos 0 da matriz de adjacencias
  do j = 1, N
     mat_adj(j, 0) = 0

     do p = 1, M
        mat_adj(j, p) = 0d0
     end do

  end do
  write(*,"(A)",advance="no") "="

  contador = 0
  progressbar = N/30

! Algoritmo para gerar o Erdos-Renyie
  do j = 1, N
     do p = 1, NN

        if ( ran1() < prob ) then

           ! Contabiliza um noh a mais para j e p
           mat_adj(j, 0) = mat_adj(j, 0) + 1
           mat_adj(p, 0) = mat_adj(p, 0) + 1

           ! Marca o peso da adjacencia usando a seguinte
           ! indexcao: X . Y (parte inteira . parte decimal)
           ! X = indice do vizinho
           ! Y = peso da ligacao

           auxp = pmax*ran1()
           mat_adj( j , int(mat_adj(j,0)) ) = auxp + p
           mat_adj( p , int(mat_adj(p,0)) ) = auxp + j

! Debuggin...
!           write(*,*) mod(p + auxp,1.0), auxp + p

        end if

     end do

     contador = contador + 1
     if ( contador .eq. progressbar ) then
        contador = 0
        write(*,"(A)",advance='no') "="
     end if

  end do
  
  write(*,"(A)") "] "
  return
end subroutine connect_network_m2





subroutine selectPairs(N, vec, Nnodes)
  
  integer, intent(in)    :: N, Nnodes
  integer, intent(inout) :: vec(Nnodes)
  
  integer                :: j, p, pp
  
! Funcoes
  real(8)                :: ran1
  
  do j = 1, Nnodes
     vec(j) = 0
  enddo
  
  j = 0
  do while ( j .le. N )

     pp = int( N*ran1() + 1 )
     p = int( pp*ran1() + 1 )
     
     if ( ( p .ne. pp )  ) then
        vec(pp) = p
        j = j + 1
     end if

  enddo
  
  return
end subroutine selectPairs
