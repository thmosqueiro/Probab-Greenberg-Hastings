!
! Retorna a atividade por neuronio
!
real(8) function getActivity(N, V) result(soma)

  integer, intent(in)    :: N, V(N)
  real(8)                :: delta

  integer                :: j

  soma = 0d0

  do j = 1, N
     soma = soma + delta(1, V(j))
  end do

  soma = soma/dfloat(N)
     
  return
end function getActivity



!
! Retorna a atividade por neuronio
!
real(8) function getActivity_subsampling(N, M, V, EI) result(soma)

  integer, intent(in)    :: N, V(N), M, EI(M)
  real(8)                :: delta

  integer                :: j

  soma = 0d0

  do j = 1, M
     soma = soma + delta(1, V( EI(j) )  )
  end do

  soma = soma/dfloat(M)
     
  return
end function getActivity_subsampling



!
! Retorna numero de neuronios atualmente excitados
!
integer function getNofExcitedNeurons(N, V) result(soma)

  integer, intent(in)    :: N, V(N)
  real(8)                :: delta

  soma = 0d0

  do j = 1, N
     soma = soma + int(delta(1, V(j)))
  end do
     
  return
end function getNofExcitedNeurons



!
! Verifica se a rede está morta.
!
! Mais rapido que getActivity em caso medio, compensa se nao for necessario
! o valor da atividade (c.c., usar getActivity).
!
logical function isDead(N, V) result(isit)

  integer, intent(in)    :: N, V(N)
  integer                :: j

  isit = .true.

  j = 1

  do while ( isit .and. ( j .le. N ) )

     if ( V(j) .eq. 1 ) then
        isit = .false.
     end if

     j = j + 1
     
  end do
     
  return
end function isDead


!
!
!
logical function isDead_loose(N, V, pl) result(isit)

  ! Input variables
  integer, intent(in)    :: N, V(N)
  real(8), intent(in)    :: pl

  ! Internal variables
  integer                :: j

  ! External function
  real(8)                :: ran1

  isit = .true.

  do j = 1, N

     if ( V(j) .eq. 1 ) then
        if ( ran1() .ge. pl ) then
           isit = .false.
           return
        end if
     end if
     
  end do
     
  return
end function isDead_loose





!
! Verifica se a parte da rede sendo provada está morta.
!
! Mais rapido que getActivity em caso medio, compensa se nao for necessario
! o valor da atividade (c.c., usar getActivity).
! *  V( EI(j) ) significa a leitura do j-esimon eletrodo.
!
logical function isDead_subsampling(N, V, nume, EI) result(isit)

  integer, intent(in)    :: N, V(N), nume, EI(nume)
  integer                :: j

  isit = .true.

  j = 1

  do while ( isit .and. ( j .le. nume ) )

     if ( V( EI(j) ) .eq. 1 ) then
        isit = .false.
     end if

     j = j + 1
     
  end do
     
  return
end function isDead_subsampling




!
! Calcula a conectividade media da rede usando a matriz de adjacencia.
!
!
real(8) function getMeanConnectivity(N, M, mat_adj) result(mc)

  integer, intent(in)    :: N, M
  real(8), intent(in)    :: mat_adj(N,0:M)
  integer                :: j

  mc = 0d0

  do j = 1, N
     mc = mc + mat_adj(j,0)
  end do

  mc = mc/dfloat(N)
     
  return
end function getMeanConnectivity



!
! Retorna o indice do neuronio mais conectado.
!
integer function get_mostConnected(N, Madj, mat_adj) result(k)

  integer, intent(in)   :: N, Madj
  real(8), intent(in)   ::  mat_adj(N, 0:Madj)
  integer               :: j, aux


  k = 1
  aux = 0
  do j = 1, N
     if ( mat_adj(j,0) .gt. aux ) then
        k = j
        aux = mat_adj(j,0)
     end if
  end do
  
  return
end function get_mostConnected


!
! Retorna o numero de conexoes.
!
integer function getMeanNofConnections(N, M, mat_adj) result(mc)

  integer, intent(in)    :: N, M
  real(8), intent(in)    :: mat_adj(N,0:M)
  integer                :: j

  mc = 0d0

  do j = 1, N
     mc = mc + int( mat_adj(j,0) )
  end do

  mc = mc/2
     
  return
end function getMeanNofConnections


!
! Calcula o average branching ratio
!
real(8) function getMeanBR(N, M, mat_adj) result(mc)

  integer, intent(in)    :: N, M
  real(8), intent(in)    :: mat_adj(N,0:M)
  integer                :: j, p, auxp

  mc = 0d0

  do j = 1, N

     auxp = int( mat_adj(j,0) )

     do p = 1, auxp
        mc = mc + mod(mat_adj(j,p), 1d0)
     end do

  end do

  mc = mc/dfloat(N)
     
  return
end function getMeanBR




!
! Calcula o average branching ratio
!
real(8) function getLargestEigenvalue(N, M, mat_adj) result(ev)

  integer, intent(in)    :: N, M
  real(8), intent(in)    :: mat_adj(N,0:M)
  integer                :: j, p, auxp
  real(8)                :: mc, mc2, mcaux

! Zerando as variaveis que vou usar pros momentos (1o e 2o).
  mc = 0d0
  mc2 = 0d0

! Iteracoes
  do j = 1, N

     auxp = int( mat_adj(j,0) )

     mcaux = 0.d0

     do p = 1, auxp
        mcaux = mcaux + mod(mat_adj(j,p), 1d0)
     end do
     
     mc = mcaux + mc
     mc2 = mcaux**2 + mc2
     
  end do
  
! Finalizando o calculo do autovalor mais alto
  ev = mc2/mc
  
  return
end function getLargestEigenvalue




!
! Calcula o average branching ratio
!
real(8) function getLargestEigenvalue_old(N, M, mat_adj, pmax) result(ev)

  integer, intent(in)    :: N, M
  real(8), intent(in)    :: mat_adj(N,0:M), pmax
  integer                :: j, p, auxp
  real(8)                :: mc, mc2, mcaux, dmax

! Zerando as variaveis que vou usar pros momentos (1o e 2o).
  mc = 0d0
  mc2 = 0d0

  dmax = 0.d0

! Iteracoes
  do j = 1, N

     auxp = int( mat_adj(j,0) )

     mcaux = 0.d0
     do p = 1, auxp

        mcaux = mcaux + mod(mat_adj(j,p), 1d0)

     end do

     mc = mcaux + mc
     mc2 = mcaux**2 + mc2

!     mc = mat_adj(j,0) + mc
!     mc2 = mat_adj(j,0)**2 + mc2


  end do

! Finalizando o calculo do autovalor mais alto
  ev = mc2/mc
  
  return
end function getLargestEigenvalue_old



!
! Calcula o branching ratio do j-esimo neuronio
!
real(8) function getBranchingRatio(N, M, mat_adj, j) result(mc)

  integer, intent(in)    :: N, M, j
  real(8), intent(in)    :: mat_adj(N,0:M)
  integer                :: p, auxp

  mc = 0d0

  auxp = int( mat_adj(j,0) )

  do p = 1, auxp
     mc = mc + mod( mat_adj(j,p), 1d0 )
  end do
     
  return
end function getBranchingRatio


!
! getActiveNeuronsNoRepeat
!
!
! This function returns the number of activie neurons that
! have never been active since simulation started.
! Matrix norep in updated during the execution of this function,
! then it is set as inout dummy variable.
!
integer function getActiveNeuronsNoRepeat(N, V, norep) result(soma)

    integer, intent(in)    :: N, V(N)
    integer, intent(inout) :: norep(N)
    real(8)                :: delta
    integer                :: t

    soma = 0

    do j = 1, N
        t = int( delta(1, V(j)) )
        if ( ( t .eq. 1 ) .and. ( norep(j) .eq. 0 ) ) then
            soma = soma + t
            norep(j) = 1
        end if
    end do

    return
end function getActiveNeuronsNoRepeat



!
! getActiveNeuronsNoRepeat
!
!
! This function returns the number of activie neurons that
! have never been active since simulation started.
! Matrix norep in updated during the execution of this function,
! then it is set as inout dummy variable.
!
subroutine getActiveNeuronsNoRepeat_loose(N, V, norep, pl, size, isDead)
  
  integer, intent(in)    :: N, V(N)
  integer, intent(inout) :: norep(N), size
  logical, intent(inout) :: isDead
  real(8), intent(in)    :: pl
  real(8)                :: delta
  integer                :: t
  
  ! External function
  real(8)                :: ran1
  
  isDead = .true.
  
  do j = 1, N

     ! Retrieve the neuronal state at this moment
     t = V(j)

     ! Now verify if it is in the excited state and, concomitantly,
     ! if its signal is not lost.
     if ( ( t .eq. 1 ) .and. ( ran1() .gt. pl ) ) then

        ! This is enough to define isDead
        isDead = .false.

        ! To increase the size, we need to be sure this neuron
        ! was never counted in the past.
        if ( ( norep(j) .eq. 0 ) ) then
           size = size + 1
           norep(j) = 1
        end if

     end if
  end do
  
  return
end subroutine getActiveNeuronsNoRepeat_loose




!
! getActiveNeurons
!
! Similar to getActivity, but returning a integer
! number that accounts for the number of active
! neurons.
!
integer function getActiveNeurons(N, V) result(soma)

    integer, intent(in)    :: N, V(N)
    real(8)                :: delta

    soma = 0

    do j = 1, N
            soma = soma + int( delta(1, V(j)) )
    end do
    
    return
end function getActiveNeurons







! ================
!
! Phase code
!
! ===============


!
! 
!
real(8) function syncOverlap(N, X, phi) result(soma)

  integer, intent(in)    :: N       ! numero de neuronios
  integer, intent(in)    :: X(N)    ! estado dos neuronios
  integer, intent(in)    :: phi(N)  ! vetor de fases
  
  ! Variaveis auxiliares
  real(8)                :: I, constante
  
  ! Variavel indice
  integer                :: j, k, jaux
  
  ! External function
  real(8)                :: delta
  
  
  ! Zerando a soma
  soma = 0d0
  
  ! Normalizacao
  constante = 1.d0/N/N
  
  
  ! Somando
  do j = 1, N
     
     I = 1.d0
     
     jaux = j-1
     do k = 1, jaux
        I = I + delta( 1, X(k) )*2.d0*cos( dfloat( phi(j) - phi(k) ) )
     enddo
     
     soma = soma + I*delta( 1, X(j) )
     
  end do

  soma = soma*constante
     
  return
end function syncOverlap











!
! Inverse of Error Function
!
real(8) function inverf(p) result(z)
  
  real(8), intent(in) :: p
  real*8 p_low,p_high
  real*8 a1,a2,a3,a4,a5,a6
  real*8 b1,b2,b3,b4,b5
  real*8 c1,c2,c3,c4,c5,c6
  real*8 d1,d2,d3,d4
  real*8 q,r
  
  a1=-39.6968302866538
  a2=220.946098424521
  a3=-275.928510446969
  a4=138.357751867269
  a5=-30.6647980661472
  a6=2.50662827745924
  b1=-54.4760987982241
  b2=161.585836858041
  b3=-155.698979859887
  b4=66.8013118877197
  b5=-13.2806815528857
  c1=-0.00778489400243029
  c2=-0.322396458041136
  c3=-2.40075827716184
  c4=-2.54973253934373
  c5=4.37466414146497
  c6=2.93816398269878
  d1=0.00778469570904146
  d2=0.32246712907004
  d3=2.445134137143
  d4=3.75440866190742
  
  p_low = 0.02425
  p_high = 1.d0 - p_low


  if ( p .le. p_low ) then
     
     q = dsqrt( - 2.d0*dlog(p) )
     
     z = ( ( ( ( (c1*q+c2)*q + c3 )*q + c4)*q + c5)*q + c6)/( ( ( ( d1*q + d2 )*q + d3 )*q + d4 )*q + 1 )
     
  else if ( ( p .ge. p_low ) .and. ( p .le. p_high ) ) then
     
     q = p - 0.5d0
     r = q**2
     z = ( ( ( ( ( a1*r + a2 )*r + a3 )*r + a4 )*r + a5 )*r + a6 )*q/( ( ( ( ( b1*r + b2 )*r + b3 )*r + b4 )*r + b5 )*r + 1 )
        
  else if ( p .gt. p_high ) then
        
     q = dsqrt( - 2.d0*dlog(1 - p) )
     z = - ( ( ( ( (c1*q + c2 )*q + c3 )*q + c4 )*q + c5 )*q + c6 )/( ( ( ( d1*q + d2 )*q + d3 )*q + d4 )*q + 1 )
     
  end if
  
  return
end function inverf



real*8 function dinvnorm(p)
      real*8 p,p_low,p_high
      real*8 a1,a2,a3,a4,a5,a6
      real*8 b1,b2,b3,b4,b5
      real*8 c1,c2,c3,c4,c5,c6
      real*8 d1,d2,d3,d4
      real*8 z,q,r
      a1=-39.6968302866538
      a2=220.946098424521
      a3=-275.928510446969
      a4=138.357751867269
      a5=-30.6647980661472
      a6=2.50662827745924
      b1=-54.4760987982241
      b2=161.585836858041
      b3=-155.698979859887
      b4=66.8013118877197
      b5=-13.2806815528857
      c1=-0.00778489400243029
      c2=-0.322396458041136
      c3=-2.40075827716184
      c4=-2.54973253934373
      c5=4.37466414146497
      c6=2.93816398269878
      d1=0.00778469570904146
      d2=0.32246712907004
      d3=2.445134137143
      d4=3.75440866190742
      p_low=0.02425
      p_high=1-p_low
      if(p.lt.p_low) goto 201
      if(p.ge.p_low) goto 301
201   q=dsqrt(-2*dlog(p))
      z=(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/ &
           ((((d1*q+d2)*q+d3)*q+d4)*q+1)
      goto 204
301   if((p.ge.p_low).and.(p.le.p_high)) goto 202
      if(p.gt.p_high) goto 302
202   q=p-0.5
      r=q*q
      z=(((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q/&
           (((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1)
      goto 204
302   if((p.gt.p_high).and.(p.lt.1)) goto 203
203   q=dsqrt(-2*dlog(1-p))
      z=-(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/&
           ((((d1*q+d2)*q+d3)*q+d4)*q+1)
204   dinvnorm=z
      return
    end function dinvnorm
