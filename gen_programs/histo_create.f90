!
! Creating histograms
!
! thmosqueiro @ 06/2011
!
! EXEMPLO: histo_create linear input 10 output
!
!

program histogram

  implicit real*8(a-h,o-z)

  character(len=200)     :: bugger, input_file, output_file, metodo
  integer                :: Ntrials, nbins
  real(8)                :: Smax, Smin

  write(*,*) 
  write(*,*) "---"
  write(*,*) 
  write(*,*) "Creating histograms"
  write(*,*) 
  write(*,*) 


  ! Lendo o metodo do histograma
  call getarg(1, bugger)
  read(bugger,*) metodo

  ! Lendo o arquivo de dados
  call getarg(2, bugger)
  read(bugger,*) input_file

  ! Lendo o numero de bins
  call getarg(3, bugger)
  read(bugger,*) nbins

  ! Lendo o arquivo de dados
  call getarg(4, bugger)
  read(bugger,*) output_file


  call find_max_min(Ntrials, Smin, Smax, input_file)

  if ( metodo .eq. 'linear' ) then 
     CALL linear_histogram(Ntrials, nbins, Smax, input_file, output_file)
  else if ( metodo .eq. 'cumulative' ) then
     CALL cumulative_histogram(Ntrials, nbins, Smax, input_file, output_file)
  else if ( metodo .eq. 'log-cumulative' ) then
     CALL log_cumulative_histogram(Ntrials, nbins, Smax, input_file, output_file)
  else
     write(*,*) 'Option not recognized...'
  end if



  write(*,*) ""
  write(*,*) 'The end, my friend'

  stop
end program histogram



subroutine find_max_min(Ntrials, Smin, Smax, input_file)

  character(len=200)     :: input_file
  integer, intent(inout) :: Ntrials
  real(8), intent(inout) :: Smax, Smin
  real(8)                :: dt
  integer                :: counter, Reason
  
  ! Abrindo o arquivo de dados
  open(unit=27, file=input_file, action="read")

  counter = 0

  read(27, *, IOSTAT=Reason) dt

  if (Reason .lt. 0) then
     write(*,*) "ARQUIVO VAZIO!!"
     stop
  end if

  Smax = dt
  Smin = dt

  counter = counter + 1


  do
     read(27, *, IOSTAT=Reason) dt

     if ( Reason .lt. 0 ) then

        Ntrials = counter

        close(27)

        return
     else
        counter = counter + 1
        
        if (dt .gt. Smax) then
           Smax = dt
        end if
        if (dt .lt. Smin) then
           Smin = dt
        end if
        
     end if

  end do
  
  return
end subroutine find_max_min


!
! Subrotina linear_histogram
!
! Constroi a densidade de probabilidades a partir de dados passados 
! em um arquivo externo. O resultado jah sai normalizado.
!
subroutine linear_histogram(Ntrials, nbins, Smax, input_file, output_file)

  implicit real*8(a-h,o-z)

  character(len=200), intent(in) :: input_file, output_file
  integer, intent(in)            :: Ntrials, nbins
  real(8), intent(in)            :: Smax
  integer                        :: aux
  real(8)                        :: dbins, soma, cons_norm, dt
  real(8), allocatable           :: hist(:)

  write(*,*) 
  write(*,*) 'Criando o histograma normal'
  write(*,*) 

  ! Tamanhho do bin calculado da forma mais obvia
  dbins = Smax/dfloat(nbins)

  ! Alocando o vetor com o histograma
  allocate(hist(1:nbins))

  ! Setando o vetor com o histograma
  do j = 1, nbins
     hist(j) = 0.d0
  enddo


  ! Normalizacao
  cons_norm = 1.d0/dfloat(Ntrials)/dbins
  

  ! Abrindo o arquivo de dados
  open(unit=27, file=input_file, action="read")

! Neste laco do fazemos a contagem de cada possivel resultado
! para gerar o histograma, normalizado para Ntrials
  do j = 1, Ntrials

     ! Lendo do arquivo de dados
     read(27,*) dt

! Serve para debug
!     write(*,*) j

     ! Posicao em que dt deve se encaixar
     aux = int(dt/dbins) + 1

     ! Fazendo a contagem de forma a jah obter algo normalizado
     hist(aux) = hist(aux) + 1.d0

  enddo

  ! Fechando arquivo
  close(27)


  ! Abrindo o arquivo de saída
  open(unit=20, file=output_file, action="write")

  ! Escrevendo os resultados no arquivo
  do j = 1, nbins
     write(20,*) (j)*dbins - dbins/2.d0, hist(j)*cons_norm
  enddo

  ! Fechando o arquivo
  close(20)
  

! Compensa verificar a normalizacao
  soma = 0.d0
  do j = 1, nbins
     soma = hist(j)*dbins + soma
  enddo
  write(*,"('Normalizacao: ', F6.4)") soma*cons_norm

end subroutine linear_histogram





!
! Subrotina cumulative_histogram
!
! Constroi a densidade de probabilidades a partir de dados passados 
! em um arquivo externo. O resultado jah sai normalizado.
!
subroutine cumulative_histogram(Ntrials, nbins, Smax, input_file, output_file)

  implicit real*8(a-h,o-z)

  character(len=200), intent(in) :: input_file, output_file
  integer, intent(in)            :: Ntrials, nbins
  real(8), intent(in)            :: Smax
  integer                        :: aux
  real(8)                        :: dbins, soma, cons_norm, dt
  real(8), allocatable           :: hist(:)

  write(*,*) 
  write(*,*) 'Criando prob. cumulativa'
  write(*,*) 

  ! Tamanhho do bin calculado da forma mais obvia
  dbins = Smax/dfloat(nbins)

  ! Alocando o vetor com o histograma
  allocate(hist(1:nbins))

  ! Setando o vetor com o histograma
  do j = 1, nbins
     hist(j) = 0.d0
  enddo

  
  ! Abrindo o arquivo de dados
  open(unit=27, file=input_file, action="read")

! Neste laco do fazemos a contagem de cada possivel resultado
! para gerar o histograma, normalizado para Ntrials
  do j = 1, Ntrials

     ! Lendo do arquivo de dados
     read(27,*) dt

     ! Posicao em que dt deve se encaixar
     aux = int(dt/dbins) + 1


!     write(*,*) dt, aux, (aux-1)*dbins, dbins
!     read(*,*)

     do k = 1, aux
        hist(k) = hist(k) + 1.d0
     enddo

  enddo

  ! Fechando arquivo
  close(27)


! Contabilizando a normalizacao
  soma = 0.d0
  do j = 1, nbins
     soma = hist(j)*dbins + soma
  enddo
  cons_norm = 1.d0/dfloat(Ntrials)


  ! Abrindo o arquivo de saída
  open(unit=20, file=output_file, action="write")

  ! Escrevendo os resultados no arquivo
  do j = 1, nbins
     write(20,*) (j-1)*dbins, hist(j)*cons_norm
  enddo

  ! Fechando o arquivo
  close(20)
  

end subroutine cumulative_histogram





!
! Subrotina cumulative_histogram
!
! Constroi a densidade de probabilidades a partir de dados passados 
! em um arquivo externo. O resultado jah sai normalizado.
!
subroutine log_cumulative_histogram(Ntrials, nbins, Smax, input_file, output_file)

  implicit real*8(a-h,o-z)

  character(len=200), intent(in) :: input_file, output_file
  integer, intent(in)            :: Ntrials, nbins
  real(8), intent(in)            :: Smax
  integer                        :: aux
  real(8)                        :: dbins, soma, cons_norm, dt
  real(8), allocatable           :: hist(:)

  write(*,*) 
  write(*,*) 'Criando prob. cumulativa'
  write(*,*) 

  ! Tamanhho do bin calculado da forma mais obvia
  dbins = dlog(Smax)/dfloat(nbins)

  ! Alocando o vetor com o histograma
  allocate(hist(1:nbins))

  ! Setando o vetor com o histograma
  do j = 1, nbins
     hist(j) = 0.d0
  enddo

  
  ! Abrindo o arquivo de dados
  open(unit=27, file=input_file, action="read")

! Neste laco do fazemos a contagem de cada possivel resultado
! para gerar o histograma, normalizado para Ntrials
  do j = 1, Ntrials

     ! Lendo do arquivo de dados
     read(27,*) dt

     ! Posicao em que dt deve se encaixar
     aux = int(log(dt)/dbins) + 1

     do k = 1, aux
        hist(k) = hist(k) + 1.d0
     enddo

  enddo

  ! Fechando arquivo
  close(27)

  ! Acertando a normalizacao
  cons_norm = 1.d0/dfloat(Ntrials)

  ! Abrindo o arquivo de saída
  open(unit=20, file=output_file, action="write")

  ! Escrevendo os resultados no arquivo
  do j = 1, nbins
     write(20,*) dexp(dbins*j), hist(j)*cons_norm
  enddo

  ! Fechando o arquivo
  close(20)
  

end subroutine log_cumulative_histogram
