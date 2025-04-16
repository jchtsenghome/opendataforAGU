      program calipv5disomap
c
c Purpose: Calculate ISOMAP 
c          based on IPV pentad data
c
c Call: L2DISTANCE
c Call: RS, TRED1, TQLRAT, TRED2, TQL2, PYTHAG
c       from EISPACK softwares
c
      include 'omp_lib.h'
c      
c      parameter (nx=1536,my=nx/2,nt=124,nvar=65,lev=16)
c      parameter (nx=180,my=89,nt=482,nvar=65,lev=16)
c whole year
c      parameter (nx=144,my=73,nt=3053,nvar=17,lev=17)  ! 198001-202212
c JJA
      parameter (nx=144,my=73,nt=774,nvar=17,lev=17)  ! 198001-202212
c      parameter (nv=12)
      parameter (md=nx*my,npca=nt)
c
      real*4 gh(nx,my)
      real y(md,nt),yt(nt,md),cove(nt,nt)
      real ytemp(nt),ave,var
c      real avemon(md,71) ! for whole year pentad
      real avemon(md,18) ! for whole year pentad
c 
c ISOMAP variables
c
      real*4 l2dist(nt,nt),l2check(nt,nt)
      real rndtmp(nt),rndist(nt,nt)
      real rj(nt,nt)
      integer indx(nt),nv,niter
c
c end isomap set
c
      real reof(nt,md),treof(md,nt)
c      real reof(10,md)
      real pca(npca)
      integer matz,ierr
      real weigen(nt),zeigen(nt,nt),fv1(nt),fv2(nt)
      real ze10(10,nt)
      real zeall(nt,nt),zeall2(nt,nt)
      character keydate*12
      real rlon1,rlon2,rlat1,rlat2
c
c the kind of L2 distance
c    nl2kind=0 use square L2 distance
c    nl2kind=1 use sqrt L2 distance
c
      integer nl2kind
c
      integer levpncep(lev)       
      data levpncep/1000, 925, 850, 700, 600, 500, 400, 300, 250, 200,
     >              150, 100, 70, 50, 30, 20, 10/
      data nv/12/
      data niter/10/
c
c==========================================================================
c
      character*50 eofpcatail
      data eofpcatail/'_EOF_sst40years.dat'/
      character*120 svtfsano_ana2_file,eofano_ana2_file,
     >              eofratio_ana2_file,eofpca_ana2_file,
     >              eofev_ana2_file
      character*120 sstavemonfile
      namelist /calipveoflst/ keydate,svtfsano_ana2_file,
     >                     eofano_ana2_file,eofratio_ana2_file,
     >                     eofpca_ana2_file,eofev_ana2_file,
     >                     rlon1,rlon2,rlat1,rlat2,
c the number of nearest neighbor
     >                     nl2kind,nv,niter
c
c===============================================================================
c
c
      open(66,file='calipv5disomap.nmlst',form='formatted')
      read(66,calipveoflst,end=110)
  110 continue
      print calipveoflst
c
      nxmy4=nx*my*4
      open(11,file=trim(svtfsano_ana2_file),access='direct',
     &     form='unformatted',recl=nxmy4,status='old')
c      open(12,file=trim(eofano_ana2_file),access='direct',
c     &     form='unformatted',recl=nxmy4,status='unknown')
      open(13,file=trim(eofratio_ana2_file),
     &     form='formatted',status='unknown')
      nt4=npca*4
      open(14,file=trim(eofpca_ana2_file),access='direct',
     &     form='unformatted',recl=nt4,status='unknown')
      open(15,file=trim(eofev_ana2_file),
     &     form='formatted',status='unknown')
c
      if (nl2kind.eq.0) then
c      open(33,file='l2dist_sst_square.dat',
      open(33,file='l2dist_NHipv5d_JJA_square_JJA_1980_2022.dat',
     &     form='formatted',status='unknown')
      else if (nl2kind.eq.1) then
c      open(33,file='l2dist_sst_sqrt.dat',
c      open(33,file='l2dist_EAsst_sqrt.dat',
c     &     form='formatted',status='unknown')
c      open(33,file='l2dist_NHipv5d_cosf_sqrt_JJA_1980_2022.dat',
      open(33,file='l2dist_NHipv5d_sqrt_EA_JJA_1980_2022.txt.dat',
c     &     form='unformatted',status='unknown')
     &     form='formatted',status='unknown')
      else
      print*,'You forgot to define nl2kind in namelist file!'
      endif
c
      open(34,file='rndist.dat',
     &     form='formatted',status='unknown')
      open(35,file='cove.dat',
     &     form='formatted',status='unknown')
c
c Pick up the selected domain
c
      rlondiff=360./float(nx)
      rlatdiff=180./float(my-1)
      npx1=nint(rlon1/rlondiff)+1
      npx2=nint(rlon2/rlondiff)+1
      mpy1=nint((rlat1+90.)/rlatdiff)+1
      mpy2=nint((rlat2+90.)/rlatdiff)+1
      print *,'npx1,npx2= ',npx1,npx2
      print *,'mpy1,mpy2= ',mpy1,mpy2
      newpx=(npx2-npx1+1)
      newpy=(mpy2-mpy1+1)
      print *,'newpx,newpy= ',newpx,newpy
      newmd=newpx*newpy
      print *,'new spatial dimension, newmd = ',newmd
c
      pi=4.*atan(1.0)
      sscale=1.e6
c
      call zilch(y,md*nt)
      call zilch(yt,nt*md)
c
c      nvarh=1
      lev500=8  ! for ipv 340K
c
      do 20 n=1,nt
c
c read NCEP r4 monthly data format, the 500 hPa layer
      ngh=n
c      print *,'n,ngh= ',n,ngh
      call reads(11,ngh,nx,my,gh)
c
c      read(11,rec=n) gh
c
      m=0
      do j=mpy1,mpy2
      do i=npx1,npx2
        m=m+1
        cosf=cos(((j-1)*rlatdiff-90.)*pi/180.)
        y(m,n)=gh(i,j)*cosf
c        y(m,n)=gh(i,j)
cc        yt(n,m)=gh(i,j)
      enddo
      enddo
      if (n.le.3) then
      print *,'m= ',m
      print *,'gh(npx1+9,1)= ',gh(npx1+9,1)
      print *,'n,y(10,n)= ',n,y(10,n)
      print *,'n,y(2831,n)= ',n,y(2831,n)
      endif
c
   20 continue
c
      call avemonthly(y,md,newmd,nt,avemon,18)
c
      sstavemonfile=
     & '/geps1/chtseng/ncepdata/19802022/ipv5d_340k_clim_EA_JJA_1980_202
     &2.dat'
      call outavemon(sstavemonfile,avemon,18,md,newmd,gh,newpx,newpy)
c
      do m=1,newmd
      do n=1,nt
      yt(n,m)=y(m,n)
      enddo
      enddo 
c
c calculate L2 Distance
c
c      call l2distance(y,md,newmd,nt,l2dist)
      call l2distance(y,md,newmd,nt,l2dist,nl2kind)
      do i=1,nt
      do j=1,nt
c      print *,'l2dist=',j,l2dist(1,j)
      l2check(i,j)=l2dist(i,j)
      enddo
      enddo
      call findmaxinl2(l2dist,nt,rmaxl2)
      print *,'The Maximal Distance is: ',rmaxl2
      call Infsetup(rndist,nt*nt,rmaxl2)
c
c calculate the shortest path by Dijkstra's alogrithm
c
      print *,'begin to count Dijkstra!'
c
      do inode=1,nt
      initnode=inode
c      call dijkstra_distance ( initnode, nt, l2check, rndtmp )
c      if (inode.eq.1) then
c      do iv=1,nt
c      write(*,'(a,i4,2x,a,f20.6)') 'rndtmp(',iv,')=',rndtmp(iv)
c      enddo
c      endif
c
c search all nodes (verteces)
c
c      do iv=1,nt
c        rndist(inode,iv)=rndtmp(iv)
c      enddo
c
c      nv=12  ! now set to 12 months defined in parameter
c
      do j=1,nt
      rndtmp(j)=l2check(inode,j)
      enddo
c
      call indexx(nt,rndtmp,indx)
c      print *,'arrange knn inode= ',inode
      do iv=1,nv+1
        rndist(inode,indx(iv))=rndtmp(indx(iv))
c        l2check(inode,indx(iv))=1000.*rmaxl2*nt
c        l2check(inode,indx(iv))=rmaxl2
      enddo
c
c      if (inode.eq.1) then
c      do iv=1,nt
c        if (rndist(inode,iv).ne.0.0) then
c        print *,'rndist(inode,iv)= ',inode,iv,rndist(inode,iv)
c        endif
c      enddo
c      print *,''
c      do iv=1,nt
c      print *,'indx(iv)= ',iv,indx(iv)
c      enddo
c      endif
c
      enddo ! end loop of inode
c
      do j=1,nt
      do i=1,nt
      rndist(i,j)=min(rndist(i,j),rndist(j,i))
c      l2check(i,j)=rndist(i,j)
      enddo
      enddo
c
c calculate the shortest path again by Floyd's Algorithm
c
      call r8mat_floyd ( nt, rndist )
c
c end the second shortest path run
c
      do i=1,nt
      do j=1,nt
      rndist(i,j)=rndist(i,j)*rndist(i,j)
      enddo
      enddo
c
      print *,' '
      print *,'check the shortest distance if it is symmetric!'
      print *,'rndist(1,100)= ',rndist(1,100) 
      print *,'rndist(100,1)= ',rndist(100,1) 
      print *,'rndist(301,100)= ',rndist(301,100) 
      print *,'rndist(100,301)= ',rndist(100,301) 
      print *,' '
c      
      do i=1,nt
      rndist(i,i)=0.0
      enddo
c
c      write(33) l2dist
c
      do i=1,nt
      write(33,333) (l2dist(i,j),j=1,nt)
  333 format(1x,10(f12.7,2x))
      write(34,344) (rndist(i,j),j=1,nt)
  344 format(1x,10(f20.6,2x))
      enddo
c      stop
c
c traditional PCA
c
c      print *,'y(5,101)= ',y(5,101)
c      print *,'yt(101,5)= ',yt(5,101)
c      call mtxmlp (yt,y,cove,nt,newmd,md,nt)
c      print *,'cove(5,101)= ',cove(5,101)
c      print *,'cove(101,5)= ',cove(5,101)
c
c define J matrix in MDS
c
      call defj(rj,nt)
      call mtxmlp(rndist,rj,cove,nt,nt,nt,nt)
      call mtxmlp(rj,cove,rndist,nt,nt,nt,nt)
      do i=1,nt
      do j=1,nt
      cove(i,j)=-0.5*rndist(i,j)
      enddo
      enddo
      do j=1,nt
      write(35,*) (cove(i,j),i=1,nt)
      enddo
      call findmaxinl2(cove,nt,rmaxl2)
      print *,'the maximux value in cove is: ',rmaxl2
c
      matz=1 ! need eigenvectors
      call rs(nt,nt,cove,weigen,matz,zeigen,fv1,fv2,IERR)
      print *,'IERR= ',ierr
      sumw=0.0
      do n=1,nt
c        if (n<=10) print *,'n,EigenValues =',n,weigen(n)
        if (n>=nt-9) print *,'n,EigenValues =',n,weigen(n)
c        if (weigen(n).gt.0.0) then
c        sumw=sumw+sqrt(weigen(n))
c        endif 
c        if (weigen(n).ge.0.0) then
        sumw=sumw+weigen(n)
c        endif
      enddo
      print *,'begin to calculate EOF!'
      print *,'sumw= ',sumw
c      do n=nt,nt-9,-1
      do n=nt,1,-1
      do i=1,nt
      nn=nt-n+1
c      ze10(nn,i)=zeigen(i,n)
c old safety calculation
c      zeall(nn,i)=zeigen(i,n)  ! old safety calculation
c
      zeall2(i,nn)=zeigen(i,n)
      enddo
      enddo
c
c      print *,'!!!!!old calculation!!!!'
c      safety calculation
c      especially using call mtxmlp2 (zeall,nt,nt,yt,nt,newmd,reof)
c      now using  
c      call mtxmlp2 (zeall,nt,nt,nt,nt,yt,nt,nt,md,newmd,reof)
c
c      print *,'!!!!!new calculation!!!!'
      call mtxmlp2 (y,md,newmd,nt,nt,zeall2,nt,nt,nt,nt,treof)
c
c      call mtxmlp2 (ze10,10,nt,yt,nt,newmd,reof)
c      call mtxmlp2 (zeigen,nt,nt,yt,nt,newmd,reof)
c
      nrec=0
c      do n=nt,1,-1     
c      do n=1,10     
      do n=1,nt
      nrec=nrec+1
c      if (weigen(nt-n+1).gt.0.0) then 
c      rratio=sqrt(weigen(nt-n+1))/sumw
c      else
c      rratio=0.
c      endif
      rratio=weigen(nt-n+1)/sumw
      write(13,133) rratio 
  133 format(1x,f12.8)
      write(15,155) weigen(nt-n+1)
  155 format(1x,f20.8)
c
c      if (n<=npca) then  ! npca=30
      do np=1,npca ! npca=30
c      pca(np)=zeall(np,n)
c      pca(np)=zeall2(n,np)
      if (weigen(nt-np+1).le.0.) then
      pca(np)=0.
      else
      pca(np)=zeall2(n,np)*sqrt(weigen(nt-np+1))
      endif
      enddo
c      endif
      call outpca(14,pca,npca,n)
c
      enddo ! end of loop nt 
c
c old calculation
c      call outflds(eofano_ana2_file,reof,nt,newmd,gh,newpx,newpy)
c new calculation
      print *,'treof(1,2)= ',treof(1,2)
      call outflds2(eofano_ana2_file,treof,nt,md,newmd,gh,newpx,newpy)
c
      close(11)
      close(13)
      close(14)
      stop
      end
c
c
      subroutine reads(nunit,nrec,nx,my,ary4)
c      real*8 ary8(nx,my)
      real*4 tary4(nx,my)
      real*4 ary4(nx,my)
      integer nunit,nrec
c      ary4=ary8
      read(nunit,rec=nrec)tary4
      do j=1,my
      do i=1,nx
      ary4(i,j)=tary4(i,j) 
      enddo
      enddo
      return
      end
c
c
      subroutine mtxmlp (a,b,c,lm,ln,mlarge,nlarge)
c      subroutine mtxmlp (a,b,c,lm,ln)
c
c  matrix multiply routine
c
c  ***input***
c
c  a:  first array in matrix product
c  b:  second array in matrix product
c  lm: order of matrices
c  ln: order of matrices
c
c  mlarge,nlarge: total dimension declaration in main program
c          other than working dimension
c          for preventing Fortran get the wrong memory address
c
c
c ***output***
c
c  c: matrix product
c
c **************************************************************
c
      integer lm,ln,mlarge,nlarge
c      real a(lm,ln),b(ln,lm),c(lm,lm)
      dimension a(lm,mlarge),b(mlarge,lm),c(lm,lm)
      call zilch(c,lm*lm)
c
      print *,'a(5,101)= ',a(5,101)
      print *,'b(101,5)= ',b(101,5)
c
      do 3 k=1,ln
      do 2 j=1,lm
      do 1 i=1,lm
c      if (a(i,k).eq.-999.) goto 1
c      if (b(k,j).eq.-999.) goto 2
      c(i,j)= c(i,j)+a(i,k)*b(k,j)
    1 continue
    2 continue
    3 continue
      return 
      end
c
c
      subroutine mtxmlp2 (a,malarge,ma,nalarge,na,
     &                    b,mblarge,mb,nblarge,nb,
     &                    c)
c
c  matrix multiply routine
c
c  ***input***
c
c  a:  first array in matrix product
c  ma: row order of matrix a
c  na: column order of matrix a
c  b:  second array in matrix product
c  mb: row order of matrix b
c  nb: column order of matrix b
c
c  NOTE: na must be equal to mb
c
c  malarge,nalarge,mblarge,nblarge
c        : total dimension declaration in main program
c          other than working dimension
c          for preventing Fortran get the wrong memory address
c
c
c ***output***
c
c  c: matrix product
c
c **************************************************************
c
      integer ma,na,mb,nb
      integer malarge,nalarge,mblarge,nblarge
c      dimension a(ma,na),b(mb,nb),c(ma,nb)
      dimension a(malarge,nalarge),b(mblarge,nblarge),
     &          c(malarge,nblarge)
c
c      call zilch(c,ma*nb)
      call zilch(c,malarge*nblarge)
c
c check!
c
      if (na.ne.mb) then
        print *,'Cannot execute matrix mulitiplication!'
        print *,'Check the dimension of Matrix!'
        stop
      endif
c
      do 3 k=1,na
      do 2 j=1,nb
      do 1 i=1,ma
c      if (a(i,k).eq.-999.) goto 1
c      if (b(k,j).eq.-999.) goto 2
      if (a(i,k).ne.-999.and.b(k,j).ne.-999.) then
        c(i,j)= c(i,j)+a(i,k)*b(k,j)
      endif
    1 continue
    2 continue
    3 continue
      return 
      end
c
c
      subroutine zilch (x,m)
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
      subroutine Infsetup (x,m,rmaxl2)
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
      real x(m),rmaxl2
      do 1 i=1,m
      x(i)= 1000.*rmaxl2*sqrt(1.*m)
c      x(i)= 1.*rmaxl2
    1 continue
      return
      end
c
c
      SUBROUTINE avevar(data,n,ave,var)
      INTEGER n
      REAL ave,var,data(n)
      INTEGER j
      REAL s,ep
      ave=0.0
      do 11 j=1,n
        ave=ave+data(j)
11    continue
      ave=ave/n
      var=0.0
      ep=0.0
      do 12 j=1,n
        s=data(j)-ave
        ep=ep+s
        var=var+s*s
12    continue
      var=(var-ep**2/n)/(n-1)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software $!6)$3D#21)8.
c
      subroutine avemonthly(ywrk,md,newmd,nt,avemon,mon)
c
c calculate monthly average and return anomaly of y
c
c mon=12 for monthly
c mon=18 for JJA pentad
c mon=71 for whole year pentad
c
      integer md,nt,mon
      real ywrk(md,nt),avemon(md,mon)
      integer nrecmon(mon)
c
      print *,'in subroutine avemonthly'
      print *,'md= ',md
      print *,'newmd= ',newmd
      print *,'ywrk(10,1)= ',ywrk(10,1)
      print *,'ywrk(10,2)= ',ywrk(10,2)
      do i=1,mon
      do j=1,newmd
      avemon(j,i)=0.0
      enddo
      nrecmon(i)=0
      enddo
c
      do n=1,nt
      imon=mod(n,mon)
      if(imon.eq.0) then
        imon=mon
        nrecmon(mon)=nrecmon(mon)+1
      else
        nrecmon(imon)=nrecmon(imon)+1
      endif
c
      do 1 m=1,newmd
      if (ywrk(m,n).eq.-9.96921e+36) goto 1
      avemon(m,imon)=avemon(m,imon)+ywrk(m,n)
    1 continue
c      print *,'The 10 point, year, ywrk(10,n)= ',n,ywrk(10,n)
      enddo
      print *,'The sum of 10 point, Jun, avemon(10,1)= ',avemon(10,1)
      print *,'The sum of 30 point, Jun, avemon(30,1)= ',avemon(30,1)
c
      do i=1,mon
c      print *,'nrecmon= ',nrecmon(i)
      do j=1,newmd
c      print *,'i,j,avemon= ',i,j,avemon(j,i)
      avemon(j,i)=avemon(j,i)/nrecmon(i)
      enddo
      enddo
      print *,'The ave of 10 point, Jun, avemon(10,1)= ',avemon(10,1)
      print *,'The ave of 30 point, Jun, avemon(30,1)= ',avemon(30,1)
c
      do n=1,nt
      imon=mod(n,mon)
      if(imon.eq.0) then
        imon=mon
      endif
      do 2 m=1,newmd
      if (ywrk(m,n).eq.-9.96921e+36) goto 2
      ywrk(m,n)=ywrk(m,n)-avemon(m,imon)
    2 continue
      enddo
      print *,'finish avemonthly'
c
      return
      end
c
c
      subroutine outflds(eofano_ana2_file,reof,nt,md,gh,nx,my)
      integer nt,md,nx,my
      real reof(nt,md),gh(nx,my)
      character*120 eofano_ana2_file
c
      nxmy4=nx*my*4
      open(12,file=trim(eofano_ana2_file),access='direct',
     &     form='unformatted',recl=nxmy4,status='unknown')
c
      do n=1,nt
      nrec=n
      m=0
      do j=1,my
      do i=1,nx 
        m=m+1
        gh(i,j)=reof(n,m)
      enddo
      enddo
c      print *,'write out record nrec= ',nrec
c      print *,'gh(100,50)= ',gh(100,50)
      write(12,rec=nrec) gh
      enddo ! end loop of n, time series
c
      close(12)
      return 
      end
c
c
      subroutine outflds2(eofano_ana2_file,reof,nt,md,newmd,gh,nx,my)
      integer nt,md,nx,my
      real reof(md,nt),gh(nx,my)
      real xx(newmd,nt)
      character*120 eofano_ana2_file
c
      nxmy4=nx*my*4
      open(12,file=trim(eofano_ana2_file),access='direct',
     &     form='unformatted',recl=nxmy4,status='unknown')
c
      do n=1,nt
      do m=1,newmd
      xx(m,n)=reof(m,n)
      enddo
      enddo
      print *,'reof(1,2)=',reof(1,2)
      print *,'xx(1,2)=',xx(1,2)
c
      do n=1,nt
      nrec=n
      m=0
      do j=1,my
      do i=1,nx 
        m=m+1
        gh(i,j)=xx(m,n)
      enddo
      enddo
c      print *,'write out record nrec= ',nrec
c      print *,'gh(100,50)= ',gh(100,50)
      write(12,rec=nrec) gh
      enddo ! end loop of n, time series
c
      close(12)
      return 
      end
c
c
      subroutine outpca(ifid,pca,nt,ncheck)
      integer ncheck,nt,ifid
      real pca(nt)
c
      write(ifid,rec=ncheck) pca
c
      return
      end
c
c
      subroutine outavemon(eofano_ana2_file,reof,nt,md,newmd,gh,nx,my)
c output monthly average quantity
      integer nt,md,nx,my
      real reof(md,nt),gh(nx,my)
      real xx(newmd,nt)
      character*120 eofano_ana2_file
c
      nxmy4=nx*my*4
      open(12,file=trim(eofano_ana2_file),access='direct',
     &     form='unformatted',recl=nxmy4,status='unknown')
c
      do n=1,nt
      do m=1,newmd
      xx(m,n)=reof(m,n)
      enddo
      enddo
      print *,'reof(1,2)=',reof(1,2)
      print *,'xx(1,2)=',xx(1,2)
c
      do n=1,nt
      nrec=n
      m=0
      do j=1,my
      do i=1,nx
        m=m+1
        gh(i,j)=xx(m,n)
      enddo
      enddo
c      print *,'write out record nrec= ',nrec
c      print *,'gh(100,50)= ',gh(100,50)
      write(12,rec=nrec) gh
      enddo ! end loop of n, time series
c
      close(12)
      return
      end
c
c
      subroutine findmaxinl2(l2dist,nt,rmaxl2)
      integer nt
      real l2dist(nt,nt)
      real rmaxl2
c
      rmaxl2=-9999.
      do i=1,nt
      do j=1,nt
      if (l2dist(i,j).gt.rmaxl2) then
        rmaxl2=l2dist(i,j)
      endif
      enddo
      enddo
c
      return
      end
