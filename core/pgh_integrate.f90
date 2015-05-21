!
! mean_activity
!
! Funcao que retorna a atividade media apos um periodo dado
! de iteracoes
!
real(8) function mean_activity_wm(N, Madj, netstate, mat_adj, m, h, Nsteps) result(media)

  integer, intent(in)    :: N, Nsteps, m, Madj
  real(8), intent(in)    :: h, mat_adj(N, 0:Madj)
  integer, intent(inout) :: netstate(N)
  
  integer                :: netstate_new(N)
  real(8)                :: SynWeight

  
! Characteristic time for the probability decay
  real(8)                :: tau
  real(8)                :: instSynWeight
  
! Last time at which each neuron received inputs
  integer                :: lasttime(N)
  
! Membrane potential rising at the last time received input
  real(8)                :: ltPotential(N)
  
  
! Indices
  integer                :: j, jj, NN, p, aux, idx
  integer                :: contador, progressbar, chave
  
! Funcoes externas
  real(8)                :: getActivity, ran1
  logical                :: isDead


  
! Setting a "good" value for tau
! Should be placed as input argument
  tau = 1.d0/20.d0
  
  
! Como o algoritmo eh sincrono, entao eh necessario
! guardar os estados atualizados em um novo vetor.
  netstate_new = netstate

! Inicializando a media
  media = 0d0
  
  
! Initializing lasttime and ltPotential arrays
  do j = 1, N
     lasttime(j) = 0
     ltPotential(j) = 0.d0
  enddo
  
  
  
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
              
              ! Checks if idx is "still" excitable
              if ( netstate( idx ) + netstate_new( idx ) .eq. 0 ) then
                 
                 ! This is the stantaneous synaptic weight, renormalized as to give an
                 ! impression of integration
                 SynWeight = mod( mat_adj(p,jj), 1d0)
                 instSynWeight = min( 1.d0 , SynWeight + 2*ltPotential(idx)*exp( (lasttime(idx) - j)*tau ) )
                 
                 ! Checks if p's action potential will excite idx
                 if ( ran1() .lt. instSynWeight ) then
                    
                    netstate_new( idx ) = 1
                    
                    ! Since a spike is triggered in idx, ltPotnetial
                    ! goes to zero.
                    ltPotential(idx) = 0.d0
                    
                 else
                    
                    ltPotential(idx) = ltPotential(idx)*exp( (lasttime(idx) - j)*tau ) + SynWeight
                    
                 endif
                 
                 ! Updating the lasttime arrival
                 lasttime(idx) = j
                 
              endif
              
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
     
  end do

  media = media/float(Nsteps)

  return
end function mean_activity_wm
