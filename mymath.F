c Simple mathematic functions
c
ccccccccccccccccccccccccccccccccccccc
c Factrization
c
      integer function myfact(n)
      implicit none
      integer,intent(in) :: n
      integer :: i
      myfact=1
      do i=1,n
        myfact=myfact*i
      enddo
      end function myfact
ccccccccccccccccccccccccccccccccccccc
c Permutation
c
      integer function myperm(n,m)
      implicit none
      integer,intent(in) :: n,m
      integer :: i
      myperm=1
      do i=0,m-1
        myperm=myperm*(n-i)
      enddo
      end function myperm
ccccccccccccccccccccccccccccccccccccc      
c Combination
c
      integer function mycomb(n,m)
      implicit none
      integer,intent(in) :: n,m
      integer :: i,nmm
c  function
      integer :: myperm,myfact

      nmm=m
      if (n-m.lt.m) nmm=n-m

      mycomb=myperm(n,nmm)/myfact(nmm)
      end function mycomb
cccccccccccccccccccccccccccccccccccc      
c Calculate m-th derivative of x^n
c (d/dx^m) x^n = n*(n-1)*...*(n-m+1)*x^(n-m)
c
      double precision function mydif(x,n,m)
      implicit none
      double precision,intent(in) :: x
      integer,intent(in) :: n,m
c function
      integer :: myperm

      if (m.gt.n) then
        mydif=0
        return
      endif

      mydif=dble(myperm(n,m))*x**(n-m)
      end function mydif

