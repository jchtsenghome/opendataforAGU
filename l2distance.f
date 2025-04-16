      subroutine l2distance(y,md,newmd,nt,l2dist,nl2kind)
c
c Purpose: calculate the L2 Distance based on different time 
c 
c
      integer md,newmd,nt
      real y(md,nt),l2dist(nt,nt)
c
c the kind of L2 distance
c    nl2kind=0 use square L2 distance
c    nl2kind=1 use sqrt L2 distance
c
      integer nl2kind
c=====================================================================
c
      call zilch(l2dist,nt*nt)
c
      do i=1,nt
      do j=i+1,nt
        do 1 m=1,newmd
          if (y(m,i).eq.-999.) goto 1
          l2dist(i,j)=l2dist(i,j)+(y(m,i)-y(m,j))**2
          l2dist(j,i)=l2dist(i,j)
    1   continue
      enddo
      enddo
c
      if (nl2kind.eq.1) then
      do i=1,nt
      do j=1,nt
      l2dist(i,j)=sqrt(l2dist(i,j))
      enddo
      enddo
      endif
c
      return
      end
c
c
      subroutine zilchbck (x,m)
c
c  subroutine to zero out arrays
c
c  ***input***
c
c  x: input array to zilch
c  m: number of elements to zilch
c
c  ***output***
c
c  x: zeroed array
c
c ****************************************************************
c
      real x(m)
      do 1 i=1,m
      x(i)= 0.0
    1 continue
      return
      end
c
c
