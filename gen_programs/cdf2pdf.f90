program cdf2pdf
!
! Purpose: To evaluate the empirical probability density function
! starting from a cumulative distribution function. 
!
! Notes: 
!  - Forward differentiation is the best choice?
!
! Record of revisions:
!  Date         Who                Description
! ==========   ===============    =====================================
! 9-12-2011     ThMosqueiro        Original code.
!
!

  implicit none

  character(len=100)     :: bugger, inputfile, outputfile

! Volatile variables
  real(8)                :: x, xold, cdf, cdfold


! Purely auxiliary variables
  integer                :: Reason
  


! Reading input file name into inputfile variable.
  call getarg(1, bugger)
  read(bugger,*) inputfile

! Openning as read-only the inputfile
  open(unit=27, file=inputfile, action="read", status="old")

! Reading output file name into outputfile variable.
  call getarg(2, bugger)
  read(bugger,*) outputfile

! Openning as read-only the inputfile
  open(unit=37, file=outputfile, action="write", status="new")  

  read(27, *) xold, cdfold

! Now read every line of the file and evaluate its contribution to
! the entropy.
  do
     read(27, *, IOSTAT=Reason) x, cdf

     if ( Reason .lt. 0 ) then
        close(27)
        close(37)
        exit 
     else
        
        write(37,*) xold, -(cdf - cdfold)/(x - xold)
        
        xold = x
        cdfold = cdf

     end if
     
  end do

  stop
end program cdf2pdf
