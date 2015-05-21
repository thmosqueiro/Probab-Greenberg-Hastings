program eval_entropy
!
! Purpose: To evaluate the entropy of a given probability
! distribution function. The binning should be linearly
! graduated.
!
! Record of revisions:
!  Date         Who                Description
! ==========   ===============    =====================================
! 9-12-2011     ThMosqueiro        Original code.
!
!

  implicit none

  character(len=100)     :: bugger, inputfile
  real(8)                :: Entropy, dx, logbase

! Volatile variables
  real(8)                :: x, prob


! Purely auxiliary variables
  real(8)                :: x0, x1, prob0, prob1
  integer                :: Reason
  


! Reading input file name into inputfile variable.
  call getarg(1, bugger)
  read(bugger,*) inputfile

! This makes the correct base for entropy
  logbase = 1.d0/dlog(2.d0)

! Openning as read-only the inputfile
  open(unit=27, file=inputfile, action="read", status="old")

! Initial condition for variables
  Entropy = 0.d0

  read(27,*) x0, prob0
  read(27,*) x1, prob1

  dx = dabs(x0 - x1)
  
  if ( prob0 .ne. 0.d0 ) then
     Entropy = - prob0*dx*dlog(prob0*dx)*logbase + Entropy
  end if
  if ( prob1 .ne. 0.d0 ) then
     Entropy = - prob1*dx*dlog(prob1*dx)*logbase + Entropy
  end if
  

! Now read every line of the file and evaluate its contribution to
! the entropy.
  do
     read(27, *, IOSTAT=Reason) x, prob

     if ( Reason .lt. 0 ) then
        close(27)
        exit 
     else

        if ( prob .ne. 0.d0 ) then
           Entropy = - prob*dx*dlog(prob*dx)*logbase + Entropy
        end if
        
     end if
     
  end do


  close(27)

  write(*,*) Entropy

  stop
end program eval_entropy
