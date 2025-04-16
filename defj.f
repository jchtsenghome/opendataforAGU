      subroutine defj(rj,n)
c
c define the J matrix in MDS
c
      integer n
      real rj(n,n)
c
      call zilch(rj,n*n)
c
      do i=1,n
        rj(i,i)=1.0
      enddo
      do i=1,n
      do j=1,n
        rj(i,j)=rj(i,j)-1./float(n)
      enddo
      enddo
c
      return
      end
