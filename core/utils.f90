!
! Gerador de logs
!
!  ----------------------> Sob teste
!
subroutine logGenerate(file_name, Nneurons, m, sigma, sigmae, K, Ke, h, Tsimul, nruns, date, logfile)
  
  character(len=100), intent(in) :: file_name, logfile, date
  integer, intent(in)            :: Nneurons, m, Tsimul, nruns, K
  real(8), intent(in)            :: h, sigma, sigmae, Ke
  
  open (unit=1, file=logfile, status='old', action='write', access='append')
  
  write(1,*)
  write(1,*) '--'
  write(1,*) file_name
  write(1,101) date
  write(1,*)
  write(1,102) Nneurons, m
  write(1,103) Tsimul, nruns
  write(1,104) sigma, sigmae
  write(1,105) K, Ke
  write(1,106) h
  write(1,*)
  write(1,*) '--'
  write(1,*)
  write(1,*)

  close(1)

101 format(" - Data: ", a)
102 format(" - N. neuronios: ", I6, t40, " - N. est. refrat.: ", I2)
103 format(" - N. max. de time steps: ", I6, t40, " - Num de repeticoes: ", I8)
104 format(" - Avg branching ratio (presc): ", F5.2, t40, " --> Medido: ", F4.2)
105 format(" - Media de conexoes (presc): ", I3, t40, " --> Medido: ", F6.3)
106 format(" - Eta: ", F6.5)
  
  return
end subroutine logGenerate




!
! Implementacao da funcao delta
!
real(8) function delta(x0, x) result(d)

  integer, intent(in) :: x, x0

  if (x .eq. x0) then
     d = 1d0
  else
     d = 0d0
  end if

  return
end function delta


subroutine progressbar_update(contador, progressbar)

  integer, intent(inout) :: contador
  integer, intent(in)    :: progressbar

  contador = contador + 1
  
  if ( contador .eq. progressbar ) then
     contador = 0
     write(*,"(A)",advance='no') "="
  end if

  return
end subroutine progressbar_update



! Funacao RAN1()
!
! Usando a funcao rotina mais nova do Fortran 95 (ou mais)
! para gerar numeros aleatorios
!
real(8) function ran1() result(x)
  CALL random_number(x)
  return
end function ran1

