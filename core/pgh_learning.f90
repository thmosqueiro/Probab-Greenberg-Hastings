

! 
! subroutine PHASE_CODING_SL
! 
! This routine receives all neuron phases and connections and 
! updates (through function dJ_stdp) all connections weighs.
! 
subroutine phase_coding_sl(N, M, sigma, matadj, phases)

  integer, intent(in)    :: N, M
  real(8), intent(inout) :: matadj(N, 0:M), sigma
  integer, intent(in)    :: phases(N)
  
  integer                :: tmax, to, dphi
  real(8)                :: weight, mbr, delta, dj
  
! Loop variables
  integer                :: ncon, j, p, l
  integer                :: contador

! Funcoes
  real(8)                :: ran1, dJ_stdp, STDP, getMeanBR
  
  
! Time for learning process
  tmax = int( 1d3 )
  
  
  write(*,*) '>> Rede aprendendo...'
  
  do j = 1, N
     
     ncon = int( matadj(j, 0) )
     write(*,*) 'Rodando para ', j
     
     do p = 1, ncon
        to = int( matadj(j, p) )
        weight = matadj(j, p) - dfloat(to)
        
        !matadj(j, p) = matadj(j, p) + dJ_stdp(weight, j, to, tmax, phases(j), phases(to))
        dj = dJ_stdp(weight, j, to, tmax, phases(j), phases(to))
        
        ! Verificando se pode passar de um!
        if ( dabs(dj) + sigma/3.60d0 .gt. 0.9 ) write(*,*) 'cacete de aguia!!'
        
        ! Atualizando o peso sinaptico
        matadj(j, p) = dfloat(to) + sigma/3.60d0 + dj
                     
     enddo
     
  enddo
  
  
  write(*,*) '>> Atualizando mbr...'
  
  mbr = getMeanBR(N, M, matadj)
  delta = sigma - mbr
  write(*,*) 'Modified mbr: ', mbr, delta
  do j = 1, N
     
     ncon = int( matadj(j, 0) )
     
     do p = 1, ncon
        matadj(j, p) = matadj(j, p) + delta
     enddo
     
  enddo  
  
  mbr = getMeanBR(N, M, matadj)
  write(*,*) 'Final mbr: ', mbr
  
  
  return
end subroutine phase_coding_sl



! 
! function dJ_stdp
! 
! This functions implements the phase-coding associative memory
! learning stage, in which the synaptic weight difference is 
! evaluated.
! 
real(8) function dJ_stdp(Jjp, j, p, tmax, phi1, phi2) result(dJ)
  
  integer, intent(in)    :: j, P, tmax, phi1, phi2
  real(8), intent(in)    :: Jjp
  
  ! Period
  integer                :: T
  
  ! Unity gauge
  real(8)                :: gauge
  
  ! Loop variables for convolution
  integer                :: k,l
  
  ! External functions
  real(8)                :: STDP, fdeterm
  
  
! Initial value of dJ
  dJ = 0.d0
  
! Setting the unity gauge
  gauge = 25.d0
  
! Period
  T = int( 10d1 )
  
  do k = 1, tmax
     do l = 1, tmax
        dJ = dJ + &
             fdeterm( k + phi1, T )*   &
             STDP( dfloat( k - l )  )* &
             fdeterm( l + phi2 , T )
     enddo
  enddo
  
  dJ = dJ/dfloat(tmax)*gauge
  
  return
end function dJ_stdp



! 
! function FDETERM
! 
! During the learning stage of the phase-code associative memory,
! this is the assumed activity of each neuron, given its specific
! phase.
! 
real(8) function fdeterm(phi, T) result(f)
  
  integer, intent(in) :: T, phi
  
  if ( mod(phi, T) .eq. 0 ) then
     f = 1.d0
     return
  else
     f = 0.d0
     return
  endif
  
  
  return
end function fdeterm


!
! function STDP
!
! This function was tested and graph can be found in Textos.
!
real(8) function STDP(dt)
  
  real(8), intent(in) :: dt
  
  !Parameters
  real(8)             :: ap, ad, Tp, Td, gamma, eta
  
  ! Auxiliary variables
  real(8)             :: C
  
! Normalization condition
  C = 1./80.
  
! Values taken from Scarpetta et al EPL 2009 
  Tp = 10.2d0
  Td = 28.6d0
  eta = 4.d0
  gamma = 42.d0
  
  ap = gamma*( 1.0/Tp + eta/Td )**(-1.)
  ad = gamma*( eta/Tp + 1.0/Td )**(-1.)
  
! Calculating the value
  if ( dt >= 0 ) then
     STPD = C*( ap*exp(-dt/Tp) - ad*exp(-eta*dt/Td) )
     return
  else
     STDP = C*( ap*exp(eta*dt/Tp) - ad*exp(dt/Td) )
     return
  end if
  
end function STDP
