      program getcontlengthv5
c
c Purpose: trace the contour line 
c          based on conrec marching triangular algorithm
c
      parameter (nx=144,my=73,lev=17,nt=365)
      parameter (nc=10,nx1=nx+1)
c
      integer yyear,dyend(12),totdate(12),npentad,nlon
      real cwbplev(lev)
      real wlon,pvu
      real distchk
c
      real*4 xx(nx,my),axx(nx,my)
      real*4 xx2(nx,my)
      real sxx(nx,my)
      real diffpv1(nx,my),diffpv2(nx,my),diffgradpv(nx,my)
      real pv2count(nx,my)
      real totzmeandist(nx),totdist(nx)
      real*4 sumcount(nx,my)
      real*4 pv2lat(nx,20)
      integer njwcount(nx)
      integer indx(nx1*my*2),mcon(nx1*my*2),kc(nx1*my*2)
c
      real*4 d(nx1,my),z(nc),x(nx1),y(my)
      real xcord(nx1*my*2),ycord(nx1*my*2)
      real xcord1(nx1*my*2),ycord1(nx1*my*2)
      real chk(nx1*my*2),chk1(nx1*my*2)
      real glon(nx),glat(my)
      real sumd(nx)
      integer numc
c
      character infilename*180,oufilename1*180,oufilename2*180,
     >          oufilename3*180,contourfilename*180,
     >          contintwofilename*180
c
      data dyend /31,28,31,30,31,30,31,31,30,31,30,31/
      data cwbplev/1000., 925., 850., 700., 600., 500., 400., 300.,
     >              250., 200., 150., 100., 70., 50., 30., 20., 10./
c
c      data thelev/270., 280., 290., 300., 310., 320., 330., 340.,
c     >            350., 360., 370., 380., 390., 400., 410., 420., 430./
c
      namelist /getcontlengthlst/yyear,infilename,
     &                           oufilename1,
     &                           wlon,npentad,pvu,
     &                           diffpvu,gradpvu
c
c======================================================================
c
      open(66,file='getcontlengthv5.nmlst',form='formatted')
      read(66,getcontlengthlst,end=110)
  110 continue
c
      nxmy4=nx*my*4      
      open(12,file=trim(infilename),access='direct',
     &     form='unformatted',recl=nxmy4,status='old')
      open(13,file=trim(oufilename1),access='direct',
     &     form='unformatted',recl=nxmy4,status='unknown')
c      open(14,file=trim(oufilename2),access='direct',
c     &     form='unformatted',recl=nxmy4,status='unknown')
c      open(15,file=trim(oufilename3),access='direct',
c     &     form='unformatted',recl=nxmy4,status='unknown')
c
c Pick up the selected domain
c
      rlondiff=360./float(nx)
      rlatdiff=180./float(my-1)
      nlon=wlon/rlondiff
      do i=1,nx
      glon(i)=real(i-1)*360./real(nx)-0.
      enddo
      do j=1,my
      glat(j)=real(j-1)*180./real(my-1)-90.
      enddo

c
      sscale=1.e6
c
      call ifleap(yyear,dyend(2))
      if (dyend(2)==29) then
        nday=366
      else
        nday=365
      endif
c
      totdate(1)=dyend(1)
      do n=2,12
        totdate(n)=totdate(n-1)+dyend(n)
      enddo
c
      print *,'total days in this year: ',totdate(12)
c
c 
c get lat=45 average distance
c
          call distgreatcir(0.,90.,0.,45.,distchk)
c
      do 1100 n=1,nday
c      do 1100 n=305,305
c
      do k=1,lev
c      do k=6,6
      call zilch(sxx,nx*my) ! every k should be
      call zilch(totzmeandist,nx*1)
      call zilch(totdist,nx*1)
      call zilch(pv2count,nx*my)
      call zilch(sumcount,nx*my)
      call zilch99(pv2lat,nx*20)
c
c count the distance
c
        irec=(n-1)*lev+k
        read(12,rec=irec) xx
c
        do j=1,my
        do i=1,nx
        d(i,j)=sscale*xx(i,j)
        enddo
        d(nx+1,j)=sscale*xx(1,j)
        enddo
c
        zmin =  1.0e30
        zmax = -1.0e30
        do i=1,nx
        do j=1,my
          zmin = min(zmin,d(i,j))
          zmax = max(zmax,d(i,j))
        enddo
        enddo
        zmin=1.
        zmax=10.
c
c     Set coordinates in Y array suitable for
c     automatic plotting on the graphics screen
c
        pymin=-90.
        pymax=90.
        do j=1,my
c          y(j) = pymax - j * (pymax - pymin) / float(my-1)
          y(j) = (j-1)* (pymax - pymin) / float(my-1)+pymin
        enddo
c
c     Set coordinates in X array suitable for
c     automatic plotting on the graphics screen
c
        pxmin=0.
        pxmax=360.
        do i=1,nx1
          x(i) = (i-1) * (pxmax - pxmin) / float(nx) + pxmin
        enddo
c
c     Set a full contingent of contour levels
c
        do i=1,nc
c          z(i) = i * (zmax - zmin) / (nc + 1)
          z(i) = float(i)
        enddo
        
c
        call conrecjohn(d,1,nx1,1,my,x,y,nc,z,xcord,ycord,chk,numc)
c
c        print *,'the number of contour point is: ',numc
        numc=0
        ii1=0
        do ii=1,nx1*my*2
          if (xcord(ii).ne.-999999) then
c          write(22,*) xcord(ii),ycord(ii),chk(ii)
c           print *,xcord(ii),ycord(ii),chk(ii)
          numc=numc+1
          ii1=ii1+1
          xcord1(ii1)=xcord(ii)
          ycord1(ii1)=ycord(ii)
          chk1(ii1)=chk(ii)
          endif
        enddo
c        print *,'the number of contour point is: ',numc
c
c
c        rewind (22)
c        nnn=0
c        do ii=1,nx*my*2
c        read(22,*,end=9999) xcord(ii),ycord(ii)
c        nnn=nnn+1
c        enddo
c 9999   continue
c        print *,'nnn= ',nnn
c
c
c        print *,'numc= ',numc
c
        call zilch(sumd,nx*1)
c
        do jdiv=3,nx-2  ! loop of jdiv
        glonmin=glon(jdiv)-5.0
        glonmax=glon(jdiv)+5.0
c
        do ii=1,numc-1,2
        if (chk1(ii).eq.0) then
        xini=xcord1(ii)
        yini=ycord1(ii)
        endif
        if (chk1(ii+1).eq.1) then
        xend=xcord1(ii+1)
        yend=ycord1(ii+1)
        endif
c
        if (xini.ge.glonmin.and.xini.le.glonmax) then
          call distgreatcir(xini,yini,xend,yend,dist)
          sumd(jdiv)=sumd(jdiv)+dist
        endif
c
        enddo
        enddo ! end loop of 3, nx-2
c
        do jdiv=1,2
        glonmin=glon(jdiv)-5.0
        glonmin=360.+glonmin
        glonmax=glon(jdiv)+5.0
        do ii=1,numc-1,2
        if (chk1(ii).eq.0) then
        xini=xcord1(ii)
        yini=ycord1(ii)
        endif
        if (chk1(ii+1).eq.1) then
        xend=xcord1(ii+1)
        yend=ycord1(ii+1)
        endif
c
        if ( (xini.ge.glonmin).or.
     &       (xini.le.glonmax) ) then
          call distgreatcir(xini,yini,xend,yend,dist)
          sumd(jdiv)=sumd(jdiv)+dist
        endif
c
        enddo
        enddo ! end loop of 1,2
c
        do jdiv=nx-1,nx
        glonmin=glon(jdiv)-5.0
        glonmax=glon(jdiv)+5.0
        glonmax=glonmax-360.
        do ii=1,numc-1,2
        if (chk1(ii).eq.0) then
        xini=xcord1(ii)
        yini=ycord1(ii)
        endif
        if (chk1(ii+1).eq.1) then
        xend=xcord1(ii+1)
        yend=ycord1(ii+1)
        endif
c
        if ( (xini.le.glonmax).or.
     &       (xini.ge.glonmin) ) then
          call distgreatcir(xini,yini,xend,yend,dist)
          sumd(jdiv)=sumd(jdiv)+dist
        endif
c
        enddo
        enddo ! end loop of nx-1, nx
c
c        do jdiv=1,nx
c        print *,'sumd= ',sumd(jdiv)
c        enddo ! end loop of jdiv
c
        do j=1,my
        do i=1,nx
        sumcount(i,j)=sumd(i)
        enddo
        enddo
c
        ireco=(n-1)*lev+k
        write(13,rec=ireco) sumcount
c
      enddo ! end loop of lev
c
 1100 continue
c
      close(12)
      close(13)
c
      stop
      end
c
c
      subroutine ifleap(yy,dyend)
      integer yy,dyend
      dyend=28
      if     ( mod(yy,400) == 0 )then
        dyend=29
      else if( mod(yy,100) == 0 )then
        dyend=28
      else if( mod(yy,4  ) == 0 )then
        dyend=29
      endif
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
      subroutine zilint (x,m)
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
      integer x(m)
      do 1 i=1,m
      x(i)= 0
    1 continue
      return
      end
c
c
      subroutine zilch99 (x,m)
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
      x(i)= -999.0
    1 continue
      return
      end
c
c
      subroutine rlonmean(sxx,nx,my,nlon,xx)
c 
c take longitudinal running mean
c nx,my: data dimensions
c nlon: the number of longitudinal grids
c sxx: input data
c xx: output data
c
      integer nx,my,nlon,igl
      real sxx(nx,my)
      real*4 xx(nx,my)
c
      do j=1,my
c
      do i=1,nlon/2
        slon=0.0
        do il=i-nlon/2,i+nlon/2
        igl=il
        if (il.le.0) igl=nx-il+1
        slon=slon+sxx(igl,j)
        enddo
c
        xx(i,j)=slon/(nlon+1)
      enddo
c
      do i=nlon/2+1,nx-nlon/2
c
        slon=0.0
        do il=i-nlon/2,i+nlon/2
        slon=slon+sxx(il,j)
        enddo
c
        xx(i,j)=slon/(nlon+1)
      enddo
c
      do i=nx-nlon/2+1,nx
        slon=0.0
        do il=i-nlon/2,i+nlon/2
        igl=il
        if (il.gt.nx) igl=il-nx
        slon=slon+sxx(igl,j)
        enddo
c
        xx(i,j)=slon/(nlon+1)
      enddo
c
      enddo  ! end loop of my
c
      return
      end
c
c
      subroutine distgreatcir(rlon1,rlat1,rlon2,rlat2,dist)
c
c Purpose: calculate great circle distance
c
      real rlon1,rlat1,rlon2,rlat2
      real pi,rad,torad,a,c,dist
c
      pi=4.0*atan(1.0)
      rad=6370. ! earth radius in km
      torad=pi/180.
      dlon=(rlon2-rlon1)*torad
      dlat=(rlat2-rlat1)*torad
      a=sin(dlat/2)*sin(dlat/2)+cos(rlat1*torad)*cos(rlat2*torad)
     >                         *sin(dlon/2)*sin(dlon/2)
      c=2*atan2(sqrt(a), sqrt(1-a))
      dist=rad*c
c      print *,'The distance between (lon1,lat1) and (lon2,lat2)= ',dist
c
      return
      end
c
c
      subroutine countlength(pv2count,nx,my,pv2length)
      integer nx
      real pv2count(nx,my),diffv(my)
c
      do j=1,my
        diffv(j)=999.
      enddo 
c
      do j=1,my
        diffv(j)=pv2count(1,j)-pv2count(nx,j)
      enddo
c
      do i=1,nx
c
      enddo
c
      return
      end
c
c
c======================================================================
c
c     CONREC is a contouring subroutine for rectangularily spaced data.
c
c     It emits calls to a line drawing subroutine supplied by the user
c     which draws a contour map corresponding to real*4 data on a randomly
c     spaced rectangular grid. The coordinates emitted are in the same
c     units given in the x() and y() arrays.
c
c     Any number of contour levels may be specified but they must be 
c     in order of increasing value.
c
c     subroutine conrec(d,ilb,iub,jlb,jub,x,y,nc,z)
c     real*4 d(ilb:iub,jlb:jub)  ! matrix of data to contour
c     integer ilb,iub,jlb,jub    ! index bounds of data matrix
c     real*4 x(ilb:iub)          ! data matrix column coordinates
c     real*4 y(jlb,jub)          ! data matrix row coordinates
c     integer nc                 ! number of contour levels
c     real*4 z(1:nc)             ! contour levels in increasing order
c
      subroutine conrecjohn(d,ilb,iub,jlb,jub,x,y,nc,z,xcord,ycord,
     &                      chk,numc)
      real*4 d(ilb:iub,jlb:jub)
      integer ilb,iub,jlb,jub
      real*4 x(ilb:iub)
      real*4 y(jlb:jub)
      integer nc,numc
      real*4 z(1:nc)
c
      real undef
      real xcord(iub*jub*2),ycord(iub*jub*2)
      real chk(iub*jub*2)
c
c     Local declarations
c
      real*4 h(0:4)
      integer sh(0:4)
      real*4 xh(0:4),yh(0:4)
      integer im(1:4),jm(1:4)
      integer case
      integer castab(-1:1,-1:1,-1:1)
      integer p1,p2
c
c     Data
c
      data im/0,1,1,0/
      data jm/0,0,1,1/
      data castab/0,0,9,  0,1,5,  7,4,8,
     1            0,3,6,  2,3,2,  6,3,0,
     2            8,4,7,  5,1,0,  9,0,0/
      data undef /-999./
c
c     Use statement functions for the line intersections
c
      xsect(p1,p2) = (h(p2)*xh(p1)-h(p1)*xh(p2))/(h(p2)-h(p1))
      ysect(p1,p2) = (h(p2)*yh(p1)-h(p1)*yh(p2))/(h(p2)-h(p1))
c
      do ii=1,iub*jub*2
        xcord(ii)=-999999.
        ycord(ii)=-999999.
      enddo
c
c     Scan the arrays, left to right, up down within rows
c
      ireg=1
c
c20    do 100 j=jub-1,jlb,-1
c         do 90 i=ilb,iub-1
20    do 100 i=ilb,iub-1
         do 90 j=jub-1,jlb,-1
            dmin = min(d(i,j),d(i,j+1),d(i+1,j),d(i+1,j+1))
            dmax = max(d(i,j),d(i,j+1),d(i+1,j),d(i+1,j+1))
            if (dmax.ge.z(1) .and. dmin.le.z(nc)) then
c               do 80 k=1,nc
               do 80 k=2,2
                  if (z(k).ge.dmin .and. z(k).le.dmax) then
                     do 22 m=4,0,-1
                        if (m.gt.0) then
                           h(m)=d(i+im(m),j+jm(m))-z(k)
                           xh(m)=x(i+im(m))
                           yh(m)=y(j+jm(m))
                        else
                           h(0)=0.25*(h(1)+h(2)+h(3)+h(4))
                           xh(0)=0.5*(x(i)+x(i+1))
                           yh(0)=0.5*(y(j)+y(j+1))
                        endif
                        if (h(m).gt.0.0) then
                           sh(m)=+1
                        else if (h(m).lt.0.0) then
                           sh(m)=-1
                        else
                           sh(m)=0
                        endif
22                   continue
c
c     Note: at this stage the relative heights of the corners and the
c     centre are in the h array, and the corresponding coordinates are
c     in the xh and yh arrays. The centre of the box is indexed by 0
c     and the 4 corners by 1 to 4 as shown below.
c     Each triangle is then indexed by the parameter m, and the 3
c     vertices of each triangle are indexed by parameters m1,m2,and m3.
c     It is assumed that the centre of the box is always vertex 2 though
c     this isimportant only when all 3 vertices lie exactly on the same
c     contour level, in which case only the side of the box is drawn.
c
c
c           vertex 4 +-------------------+ vertex 3
c                    | \               / |
c                    |   \    m-3    /   |
c                    |     \       /     |
c                    |       \   /       |
c                    |  m=2    X   m=2   |       the centre is vertex 0
c                    |       /   \       |
c                    |     /       \     |
c                    |   /    m=1    \   |
c                    | /               \ |
c           vertex 1 +-------------------+ vertex 2
c
c
c
c                    Scan each triangle in the box
c
c
                     do 60 m=1,4
c
                        ireg=ireg+1
c
                        m1=m
                        m2=0
                        if (m.ne.4) then
                           m3=m+1
                        else
                           m3=1
                        endif
                        case = castab(sh(m1),sh(m2),sh(m3))
                        if (case.ne.0) then
                           goto (31,32,33,34,35,36,37,38,39),case
c
c     Case 1 - Line between vertices 1 and 2
c
31                            x1=xh(m1)
                              y1=yh(m1)
                              x2=xh(m2)
                              y2=yh(m2)
                              goto 40
c
c     Case 2 - Line between vertices 2 and 3
c
32                            x1=xh(m2)
                              y1=yh(m2)
                              x2=xh(m3)
                              y2=yh(m3)
                              goto 40
c
c     Case 3 - Line between vertices 3 and 1
c
33                            x1=xh(m3)
                              y1=yh(m3)
                              x2=xh(m1)
                              y2=yh(m1)
                              goto 40
c
c     Case 4 - Line between vertex 1 and side 2-3
c
34                            x1=xh(m1)
                              y1=yh(m1)
                              x2=xsect(m2,m3)
                              y2=ysect(m2,m3)
                              goto 40
c
c     Case 5 - Line between vertex 2 and side 3-1
c
35                            x1=xh(m2)
                              y1=yh(m2)
                              x2=xsect(m3,m1)
                              y2=ysect(m3,m1)
                              goto 40
c
c     Case 6 - Line between vertex 3 and side 1-2
c
36                            x1=xh(m3)
                              y1=yh(m3)
                              x2=xsect(m1,m2)
                              y2=ysect(m1,m2)
                              goto 40
c
c     Case 7 - Line between sides 1-2 and 2-3
c
37                            x1=xsect(m1,m2)
                              y1=ysect(m1,m2)
                              x2=xsect(m2,m3)
                              y2=ysect(m2,m3)
                              goto 40
c
c     Case 8 - Line between sides 2-3 and 3-1
c
38                            x1=xsect(m2,m3)
                              y1=ysect(m2,m3)
                              x2=xsect(m3,m1)
                              y2=ysect(m3,m1)
                              goto 40
c
c     Case 9 - Line between sides 3-1 and 1-2
c
39                            x1=xsect(m3,m1)
                              y1=ysect(m3,m1)
                              x2=xsect(m1,m2)
                              y2=ysect(m1,m2)
                              goto 40
c
c40                         call vecout(x1,y1,x2,y2,z(k))
c
   40                       continue
c                            print *,'x1= ',x1,'  y1= ',y1
c                            print *,'x2= ',x2,'  y2= ',y2
c                            print *,x1,y1
c                            print *,x2,y2
                            xcord(ireg)=x1
                            ycord(ireg)=y1
                            chk(ireg)=0
                            ireg=ireg+1
                            xcord(ireg)=x2
                            ycord(ireg)=y2
                            chk(ireg)=1
                        endif
60                   continue
                  endif
80             continue
            endif
c                     print *,undef,undef
c                            xcord(ireg)=undef
c                            ycord(ireg)=undef
c                            ireg=ireg+1
c                            xcord(ireg)=undef
c                            ycord(ireg)=undef
90       continue
100   continue
c
      numc=0
      do ii=1,iub*jub*2
        if (xcord(ii).ne.-999999.) numc=numc+1 
      enddo
c      print *,'in conrec, numc= ',numc
c
      return
      end
c
c======================================================================
c
c     This is a sample vector output routine. For a local environment
c     either replace the VECOUT call in the main line, or better
c     replace the contents of this subroutine between the *'s shown.
c
c     There is often the requirement to distinguish each contour
c     line with a different colour or a different line style. This
c     can be done in many ways using the contour values z for a
c     particular line segment.
c
      subroutine vecout(x1,y1,x2,y2,z)
      implicit none
      real*4 x1,y1,x2,y2,z
c
c***** Replace from here
c
c     The following should be ignored since it is specific to
c     the version of FORTRAN running on the Macintosh microcomputer
c     on which this particular example was written.
c
      INTEGER LINETO
      PARAMETER (LINETO=89109000)
      INTEGER MOVETO
      PARAMETER (MOVETO=89309000)
c      call toolbx(MOVETO,nint(x1),nint(y1))
c      call toolbx(LINETO,nint(x2),nint(y2))
c
c***** To here
c
      return
      end
c
c
      subroutine getonegrid(xx,nx,my,xxv
     >                     ,glon,glat,stalon,stalat)
c
c get Radiosound's longitude and latitude
c
      parameter (nxt=128,myt=64)      ! NCEP and new CWB
c      parameter (nxt=64,myt=32)      ! CWB
      real*4 xx(nx,my)
      real xxv,stalon,stalat
      real glon(nx),glat(my)
      real indx(nxt),indy(myt)
      real x1a,x1b,y1a,y1b,t,u
      integer nx,my,lev
c
      if (stalon.lt.0) then
        stalon=360.+stalon
      endif
c      print *,'station lon,lat= ',stalon,stalat
c
c Interpolation
c
      imin=1
      diffmin=9999.
      do i=1,nx
        diff=abs(stalon-glon(i))
        if (diff.lt.diffmin) then
          diffmin=diff
          imin=i
        endif
      enddo
      if (stalon.lt.glon(imin)) then
        imin1=imin-1
        imin2=imin
      else if (stalon.ge.glon(imin)) then
        imin1=imin
        imin2=imin+1
        if (imin2.gt.nx) imin2=imin2-nx
      endif
c
      jmin=1
      diffmin=9999.
      do j=1,my
        diff=abs(stalat-glat(j))
        if (diff.lt.diffmin) then
          diffmin=diff
          jmin=j
        endif
      enddo
      if (stalat.lt.glat(jmin)) then
        jmin1=jmin-1
        jmin2=jmin
      else if (stalat.ge.glat(jmin)) then
        jmin1=jmin
        jmin2=jmin+1
      endif
c      print *,'imin, jmin = ',imin,jmin
c      print *,'imin1, jmin1 = ',imin1,jmin1
c      print *,'imin1, jmin2= ',imin1,jmin2
c      print *,'imin2, jmin2= ',imin2,jmin2
c      print *,'imin2, jmin1= ',imin2,jmin1
c-------------------------------------------------------------
c      xxv=(xx(imin1,jmin1)+xx(imin1,jmin2)+
c     >        xx(imin2,jmin1)+xx(imin2,jmin2))/4.
      x1a=(imin1-1)*2.5
      x1b=(imin2-1)*2.5
      y1a=(jmin1-1)*2.5-90.
      y1b=(jmin2-1)*2.5-90.
      t=(stalon-x1a)/(x1b-x1a)
      u=(stalat-y1a)/(y1b-y1a)
      xxv=(1.0-t)*(1.0-u)*xx(imin1,jmin1)+t*(1.0-u)*xx(imin2,jmin1)
     >   +t*u*xx(imin2,jmin2)+(1.0-t)*u*xx(imin1,jmin2)
      return
      end
c
c
      SUBROUTINE indexx(n,arr,indx)
      INTEGER n,indx(n),M,NSTACK
      REAL arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      REAL a
      do 11 j=1,n
        indx(j)=j
11    continue
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do 12 i=j-1,1,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=0
2         indx(i+1)=indxt
13      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(l)))then
          itemp=indx(l+1)
          indx(l+1)=indx(l)
          indx(l)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(jstack.gt.NSTACK) then
          print *, 'NSTACK too small in indexx'
          stop
        endif
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END
c
c
      SUBROUTINE indexxZ(n,arr,indx)
      INTEGER n,indx(n),M,NSTACK
      INTEGER arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      INTEGER a
      do 11 j=1,n
        indx(j)=j
11    continue
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do 12 i=j-1,1,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=0
2         indx(i+1)=indxt
13      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(l)))then
          itemp=indx(l+1)
          indx(l+1)=indx(l)
          indx(l)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(jstack.gt.NSTACK) then
          print *, 'NSTACK too small in indexx'
          stop
        endif
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END
C  (C) Copr. 1986-92 Numerical Recipes Software $!6)$3D#21)8.
