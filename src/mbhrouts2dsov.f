cc Copyright (C) 2017: Travis Askham
cc Contact: askhamwhat@gmail.com
cc 
cc This program is free software; you can redistribute it and/or modify 
cc it under the terms of the GNU General Public License as published by 
cc the Free Software Foundation; either version 2 of the License, or 
cc (at your option) any later version.  This program is distributed in 
cc the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
cc even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
cc PARTICULAR PURPOSE.  See the GNU General Public License for more 
cc details. You should have received a copy of the GNU General Public 
cc License along with this program; 
cc if not, see <http://www.gnu.org/licenses/>.

      subroutine mbh2dsov_circvals(u,up,source,ifcharge,charge,
     1     ifdipole,dipstr,dipvec,ifquad,quadstr,quadvec,
     2     ifoct,octstr,octvec,ns,center,rad,beta,nterms,zcirc)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     INPUT:
c
c     source(2,ns)    : coordinates of sources
c     ifcharge        : flag, IFCHARGE = 1 means there are charges
c     charge(ns)      : charge strengths
c     ifdipole        : flag, IFDIPOLE = 1 means there are dipoles
c     dipstr(ns)      : dipole strengths
c     dipvec(2,ns)    : dipole directions
c     ifquad          : flag, IFQUAD = 1 means there are quadrupoles
c     quadstr(ns)     : quadrupole strengths
c     quadvec(3,ns)   : quadrupole directions
c     ifoct           : flag, IFOCT = 1 means there are octopoles
c     octstr(ns)      : octopole strengths
c     octvec(4,ns)    : octopole directions
c     ns              : number of sources
c     center(2)       : expansion center
c     rad             : radius of circle where field values are matched
c     zcirc(2,nterms) : precomputed evenly-spaced points on unit circle
c     beta            : the modified biharmonic parameter
c     nterms          : number of points on circle
c
c     OUTPUT:
c
c     u               : values of field induced by the charges 
c                       on circle centered at CENTER with radius RAD
c     up              : values of normal derivative of field
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global 
      integer ns, nterms
      integer ifcharge, ifdipole, ifquad, ifoct
      real *8 source(2,*), charge(*), center(2), rad, beta
      real *8 u(*), up(*), dipstr(*), quadstr(*), octstr(*)
      real *8 dipvec(2,*), quadvec(3,*), octvec(4,*)
      real *8 zcirc(2,*)
c     local
      integer i, j
      integer ifpot, ifgrad, ifhess, ifder3, ifder4, ifder5
      real *8 pot, grad(2), hess(3), der3(4), der4(5), der5(6)
      real *8 zx(2), zy(2)

      ifpot = 1
      ifgrad = 1
      ifhess = 1
      ifder3 = 1
      ifder4 = 1
      ifder5 = 0

c     note: it's important to recenter so that the points truly
c     lie on a circle

      do i = 1,nterms
         zx(1) = rad*zcirc(1,i)
         zx(2) = rad*zcirc(2,i)
         u(i) = 0.0d0
         up(i) = 0.0d0
         do j = 1,ns
            zy(1) = source(1,j)-center(1)
            zy(2) = source(2,j)-center(2)
            call modbhgreen_all(beta,zx,zy,ifpot,pot,
     1           ifgrad,grad,ifhess,hess,ifder3,der3,ifder4,der4,
     2           ifder5,der5)
            if (ifcharge .eq. 1) then
               u(i) = u(i) + charge(j)*pot
               up(i) = up(i) + charge(j)*(zcirc(1,i)*grad(1) 
     1              + zcirc(2,i)*grad(2))
            endif
            if (ifdipole .eq. 1) then
               u(i) = u(i) - dipstr(j)*(dipvec(1,j)*grad(1) 
     1              + dipvec(2,j)*grad(2))
               up(i) = up(i) - dipstr(j)*(dipvec(1,j)*
     1              (zcirc(1,i)*hess(1)+zcirc(2,i)*hess(2))+dipvec(2,j)*
     2              (zcirc(1,i)*hess(2)+zcirc(2,i)*hess(3)))
            endif
            if (ifquad .eq. 1) then
               u(i) = u(i) + quadstr(j)*(quadvec(1,j)*hess(1) + 
     1              quadvec(2,j)*hess(2) + quadvec(3,j)*hess(3))
               up(i) = up(i) + quadstr(j)*(quadvec(1,j)*
     1              (zcirc(1,i)*der3(1)+zcirc(2,i)*der3(2)) + 
     2              quadvec(2,j)*(zcirc(1,i)*der3(2) + 
     3              zcirc(2,i)*der3(3)) + quadvec(3,j)*
     4              (zcirc(1,i)*der3(3) + zcirc(2,i)*der3(4)))
            endif
            if (ifoct .eq. 1) then
               u(i) = u(i) - octstr(j)*(octvec(1,j)*der3(1) + 
     1              octvec(2,j)*der3(2) + octvec(3,j)*der3(3) +
     2              octvec(4,j)*der3(4))
               up(i) = up(i) - octstr(j)*(octvec(1,j)*
     1              (zcirc(1,i)*der4(1)+zcirc(2,i)*der4(2)) + 
     2              octvec(2,j)*(zcirc(1,i)*der4(2) + 
     3              zcirc(2,i)*der4(3)) + octvec(3,j)*
     4              (zcirc(1,i)*der4(3) + zcirc(2,i)*der4(4))
     5              + octvec(4,j)*(zcirc(1,i)*der4(4) + 
     6              zcirc(2,i)*der4(5)))
            endif
         enddo
      enddo

      return
      end
      
      subroutine mbh2dsov_dtocmp(beta,rad,rscale,u,up,npts,mbhmp,
     1     ymp,wsave,work)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     This is a wrapper routine. The real work is carried out by 
c     MBH2DSOV_DTOC1. 
c
c     INPUT
c
c     beta - REAL *8, the parameter of the modified biharmonic equation
c     rad - REAL *8, the radius of the circle on which u and up are 
c           sampled
c     rscale - REAL *8, the scaling parameter for the basis functions
c     u - REAL *8 array of length NPTS, the values of the field at 
c         equispaced points on a circle of radius rad
c     up - REAL *8 array of length NPTS, the values of the derivative
c          of the field in the radial direction at the same points as u
c     npts - INTEGER, as described above
c     wsave - REAL *8 array of length NPTS, the precomputed array 
c             returned by DFFTI from the FFT library FFTPACK
c     work - REAL *8 work array. recommended length of work array:
c            14*npts + 500
c
c     OUTPUT
c
c     mbhmp, ymp - REAL *8 arrays of length npts which 
c                        give the coefficients for the basis 
c                        functions which approximate this field 
c                        outside the disk
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    
      implicit none
      complex *16 mbhmp(0:(npts+1)/2),ymp(0:(npts+1)/2)
      real *8 beta, rad, u(*), up(*), work(*)
      real *8 wsave(*), rscale
      integer npts
c     local variables
      integer idiffs, iders, ihvec, ihder, ikvec, ikder, icu, icup
      integer ldiffs, lhvec, lkvec, lcu, l, lused

      l = (npts+1)/2

      ldiffs = l+1
      lhvec = 2*(l+6)
      lkvec = l+6
      lcu = npts

      idiffs = 1
      iders = idiffs + ldiffs
      ihvec = iders + ldiffs
      ihder = ihvec + lhvec
      ikvec = ihder + lhvec
      ikder = ikvec + lkvec
      icu = ikder + lkvec
      icup = icu + lcu
      lused = icup + lcu

      call mbh2dsov_dtocmp1(beta,rad,rscale,u,up,npts,
     1     mbhmp,ymp,wsave,work(idiffs),work(iders),
     2     work(ihvec),work(ihder),work(ikvec),work(ikder),
     3     work(icu),work(icup))


      return
      end
      
      subroutine mbh2dsov_dtocmp1(beta,rad,rscale,u,up,npts,
     1     mbhmp,ymp,wsave,diffs,ders,hvec,hder,
     2     kvec,kder,cu,cup)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     see MBH2DSOV_DTOC for description
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      complex *16 hvec(0:1), hder(0:1), mbhmp(0:(npts+1)/2)
      complex *16 ymp(0:(npts+1)/2)
      real *8 kvec(0:1), kder(0:1), wsave(*), rscale
      real *8 beta, rad, u(*), up(*)
      real *8 diffs(0:1), ders(0:1), cu(*), cup(*)
      integer npts
c     local variables
      integer i, l, ifder, ifders, j, ierge
      real *8 pih, norm
      real *8 a, b, c, d, alphan, betan, alphad, betad
      real *8 cu2j, cup2j, cu2jp1, cup2jp1, gamma
      real *8 amat(2,2), x(2), y(2)
      complex *16 z, eye
      data eye /(0.0d0,1.0d0) /

      do i = 1,npts
         cu(i) = u(i)
         cup(i) = up(i)
      enddo

      call dfftf(npts,cu,wsave)
      call dfftf(npts,cup,wsave)

      norm = 2.0d0/npts

      cu(1) = cu(1)/npts
      cup(1) = cup(1)/npts

      do i = 2,npts
         cu(i) = cu(i)*norm
         cup(i) = cup(i)*norm
      enddo

      l = (npts+1)/2

      z = eye*rad*beta
      ifder = 1
      call h2dall(l+5,z,rscale,hvec,ifder,hder)

      pih = 2.0d0*atan(1.0d0)

      do j = 0,l+1,4
         kvec(j) = -dimag(hvec(j))*pih
         kvec(j+1) = -dreal(hvec(j+1))*pih
         kvec(j+2) = dimag(hvec(j+2))*pih
         kvec(j+3) = dreal(hvec(j+3))*pih
         kder(j) = -dreal(hder(j))*pih*beta
         kder(j+1) = dimag(hder(j+1))*pih*beta
         kder(j+2) = dreal(hder(j+2))*pih*beta
         kder(j+3) = -dimag(hder(j+3))*pih*beta
      enddo

      ifders = 1
      call diffslogbk_fast(rad,beta,rscale,diffs,ifders,ders,kvec,l)

      amat(1,1) = diffs(0)
      amat(2,1) = ders(0)
      amat(1,2) = kvec(0)
      amat(2,2) = kder(0)
      
      y(1) = cu(1)
      y(2) = cup(1)
      call mbh2dsov_ge22cp(ierge,amat,x,y)  

      mbhmp(0) = x(1)
      ymp(0) = x(2)

      do j = 1,l-1

         amat(1,1) = diffs(j)
         amat(2,1) = ders(j)
         amat(1,2) = kvec(j)
         amat(2,2) = kder(j)

         y(1) = cu(2*j)
         y(2) = cup(2*j)
         call mbh2dsov_ge22cp(ierge,amat,x,y)  

         mbhmp(j) = x(1)
         ymp(j) = x(2)

         y(1) = cu(2*j+1)
         y(2) = cup(2*j+1)
         call mbh2dsov_ge22cp(ierge,amat,x,y)  

         mbhmp(j) = mbhmp(j) - eye*x(1)
         ymp(j) = ymp(j) - eye*x(2)

      enddo
      

      return
      end

      subroutine mbh2dsov_dtocta(beta,rad,rscale,u,up,npts,
     1     mbhloc,lloc,wsave,work)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     This is a wrapper routine. The real work is carried out by 
c     MBH2DSOV_DTOCTA1. 
c
c     INPUT
c
c     beta - REAL *8, the parameter of the modified biharmonic equation
c     rad - REAL *8, the radius of the circle on which u and up are 
c           sampled
c     rscale - REAL *8, the scaling parameter for the basis functions
c     u - REAL *8 array of length NPTS, the values of the field at 
c         equispaced points on a circle of radius rad
c     up - REAL *8 array of length NPTS, the values of the derivative
c          of the field in the radial direction at the same points as u
c     npts - INTEGER, as described above
c     wsave - REAL *8 array of length NPTS, the precomputed array 
c             returned by DFFTI from the FFT library FFTPACK
c     work - REAL *8 work array. recommended length of work array:
c            14*npts + 500
c
c     OUTPUT
c
c     mbhloc          : coefficients for difference-type functions
c     lloc            : coefficients for z^k functions
c    
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    
      implicit none
      complex *16 mbhloc(0:(npts+1)/2), lloc(0:(npts+1)/2)
      real *8 beta, rad, u(*), up(*), work(*)
      real *8 wsave(*), rscale
      integer npts
c     local variables
      integer idiffs, iders, ifjs, ifjder, iivec, iider, icu, icup
      integer iiscale
      integer ldiffs, lwfjs, lfjder, livec, l, lused, lcu, liscale

      l = (npts+1)/2

      ldiffs = l+1
      lwfjs = 2*(l + 5 + 4*l + 100)
      liscale = l + 5 + 4*l + 100
      lfjder = 2*(l+6)
      livec = l+6
      lcu = npts

      idiffs = 1
      iders = idiffs + ldiffs
      ifjs = iders + ldiffs
      ifjder = ifjs + lwfjs
      iivec = ifjder + lfjder
      iider = iivec + livec
      icu = iider + livec
      icup = icu + lcu
      iiscale = icup + lcu
      lused = iiscale + liscale

      call mbh2dsov_dtocta1(beta,rad,rscale,u,up,npts,
     1     mbhloc,lloc,wsave,work(idiffs),work(iders),
     2     lwfjs,work(ifjs),work(ifjder),work(iiscale),
     3     work(iivec),work(iider),
     4     work(icu),work(icup))

      return
      end

      subroutine mbh2dsov_dtocta1(beta,rad,rscale,u,up,npts,
     1     mbhloc,lloc,wsave,diffs,ders,lwfjs,fjs,fjder,
     2     iscale,pow,dpow,cu,cup)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     see MBH2DSOV_DTOCTA for description
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      complex *16 fjs(0:1), fjder(0:1)
      complex *16 mbhloc(0:(npts+1)/2), lloc(0:(npts+1)/2)
      real *8 pow(0:1), dpow(0:1), wsave(*), rscale
      real *8 beta, rad, u(*), up(*)
      real *8 diffs(0:1), ders(0:1), cu(*), cup(*)
      integer iscale(*)
      integer npts, lwfjs
c     local variables
      integer i, l, ifder, ifders, j, ntop, ier, ierge
      real *8 pih, norm
      real *8 a, b, c, d, alphan, betan, alphad, betad
      real *8 rtemp
      real *8 cu2j, cup2j, cu2jp1, cup2jp1, gamma
      real *8 amat(2,2), x(2), y(2)
      complex *16 z, eye
      data eye /(0.0d0,1.0d0) /

      do i = 1,npts
         cu(i) = u(i)
         cup(i) = up(i)
      enddo

c     fourier transform u and up values

      call dfftf(npts,cu,wsave)
      call dfftf(npts,cup,wsave)

      norm = 2.0d0/npts

      cu(1) = cu(1)/npts
      cup(1) = cup(1)/npts

      do i = 2,npts
         cu(i) = cu(i)*norm
         cup(i) = cup(i)*norm
      enddo

c     obtain values and derivatives of basis functions
c     for radius of circle

      l = (npts+1)/2

      ifders = 1
      call diffszkik_fast(rad,beta,rscale,diffs,ifders,ders,pow,l)

      call mbh2dsov_rks(pow,dpow,rad,beta,rscale,l)


c     compute zero-th mode coefficients

      amat(1,1) = diffs(0)
      amat(2,1) = ders(0)
      amat(1,2) = pow(0)
      amat(2,2) = dpow(0)

      y(1) = cu(1)
      y(2) = cup(1)
      call mbh2dsov_ge22cp(ierge,amat,x,y)  

      mbhloc(0) = x(1)
      lloc(0) = x(2)

      do j = 1,l-1

c     compute coefficients for higher modes
c     these are separated into cos(j theta) and 
c     sin(j theta) parts by DFFTF

         amat(1,1) = diffs(j)
         amat(2,1) = ders(j)
         amat(1,2) = pow(j)
         amat(2,2) = dpow(j)

         y(1) = cu(2*j)
         y(2) = cup(2*j)

         call mbh2d_ge22cp(ierge,amat,x,y)  

         mbhloc(j) = x(1)
         lloc(j) = x(2)

         y(1) = cu(2*j+1)
         y(2) = cup(2*j+1)

         call mbh2d_ge22cp(ierge,amat,x,y)  

         mbhloc(j) = mbhloc(j) - eye*x(1)
         lloc(j) = lloc(j) - eye*x(2)

      enddo

      return
      end
      
      subroutine mbh2dsov_dtocmp_naive(beta,rad,rscale,rscalelap,u,up,
     1     npts,lmp,ymp,wsave,work)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     This is a wrapper routine. The real work is carried out by 
c     MBH2DSOV_DTOC1. 
c
c     INPUT
c
c     beta - REAL *8, the parameter of the modified biharmonic equation
c     rad - REAL *8, the radius of the circle on which u and up are 
c           sampled
c     rscale - REAL *8, the scaling parameter for the basis functions
c     u - REAL *8 array of length NPTS, the values of the field at 
c         equispaced points on a circle of radius rad
c     up - REAL *8 array of length NPTS, the values of the derivative
c          of the field in the radial direction at the same points as u
c     npts - INTEGER, as described above
c     wsave - REAL *8 array of length NPTS, the precomputed array 
c             returned by DFFTI from the FFT library FFTPACK
c     work - REAL *8 work array. recommended length of work array:
c            14*npts + 500
c
c     OUTPUT
c
c     mbhmp, ymp - REAL *8 arrays of length npts which 
c                        give the coefficients for the basis 
c                        functions which approximate this field 
c                        outside the disk
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    
      implicit none
      complex *16 lmp(0:(npts+1)/2),ymp(0:(npts+1)/2)
      real *8 beta, rad, u(*), up(*), work(*)
      real *8 wsave(*), rscale, rscalelap
      integer npts
c     local variables
      integer idiffs, iders, ihvec, ihder, ikvec, ikder, icu, icup
      integer ldiffs, lhvec, lkvec, lcu, l, lused

      l = (npts+1)/2

      ldiffs = l+1
      lhvec = 2*(l+6)
      lkvec = l+6
      lcu = npts

      idiffs = 1
      iders = idiffs + ldiffs
      ihvec = iders + ldiffs
      ihder = ihvec + lhvec
      ikvec = ihder + lhvec
      ikder = ikvec + lkvec
      icu = ikder + lkvec
      icup = icu + lcu
      lused = icup + lcu

      call mbh2dsov_dtocmp1_naive(beta,rad,rscale,rscalelap,u,up,npts,
     1     lmp,ymp,wsave,work(idiffs),work(iders),
     2     work(ihvec),work(ihder),work(ikvec),work(ikder),
     3     work(icu),work(icup))


      return
      end
      
      subroutine mbh2dsov_dtocmp1_naive(beta,rad,rscale,rscalelap,u,up,
     1     npts,lmp,ymp,wsave,pow,dpow,hvec,hder,
     2     kvec,kder,cu,cup)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     see MBH2DSOV_DTOC for description
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      complex *16 hvec(0:1), hder(0:1), ymp(0:(npts+1)/2)
      complex *16 lmp(0:(npts+1)/2)
      real *8 kvec(0:1), kder(0:1), wsave(*), rscale, rscalelap
      real *8 beta, rad, u(*), up(*)
      real *8 pow(0:1), dpow(0:1), cu(*), cup(*)
      integer npts
c     local variables
      integer i, l, ifder, ifders, j, ierge
      real *8 pih, norm, pi2
      real *8 a, b, c, d, alphan, betan, alphad, betad
      real *8 cu2j, cup2j, cu2jp1, cup2jp1, gamma
      real *8 amat(2,2), x(2), y(2)
      complex *16 z, eye
      data eye /(0.0d0,1.0d0) /

      do i = 1,npts
         cu(i) = u(i)
         cup(i) = up(i)
      enddo

      call dfftf(npts,cu,wsave)
      call dfftf(npts,cup,wsave)

      norm = 2.0d0/npts

      cu(1) = cu(1)/npts
      cup(1) = cup(1)/npts

      do i = 2,npts
         cu(i) = cu(i)*norm
         cup(i) = cup(i)*norm
      enddo

      l = (npts+1)/2

      z = eye*rad*beta
      ifder = 1
      call h2dall(l+5,z,rscale,hvec,ifder,hder)

      pih = 2.0d0*atan(1.0d0)

      do j = 0,l+1,4
         kvec(j) = -dimag(hvec(j))*0.25d0
         kvec(j+1) = -dreal(hvec(j+1))*0.25d0
         kvec(j+2) = dimag(hvec(j+2))*0.25d0
         kvec(j+3) = dreal(hvec(j+3))*0.25d0
         kder(j) = -dreal(hder(j))*beta*0.25d0
         kder(j+1) = dimag(hder(j+1))*beta*0.25d0
         kder(j+2) = dreal(hder(j+2))*beta*0.25d0
         kder(j+3) = -dimag(hder(j+3))*beta*0.25d0
      enddo

      call mbh2dsov_rkinv(pow,dpow,rad,rscalelap,l)

      pi2 = 8.0d0*datan(1.0d0)
      
      do j = 0,l
         pow(j) = pow(j)/pi2
         dpow(j) = dpow(j)/pi2
      enddo
      
      amat(1,1) = pow(0)
      amat(2,1) = dpow(0)
      amat(1,2) = kvec(0)
      amat(2,2) = kder(0)
      
      y(1) = cu(1)
      y(2) = cup(1)
      call mbh2dsov_ge22cp(ierge,amat,x,y)  

      lmp(0) = x(1)
      ymp(0) = x(2)

      do j = 1,l-1

         amat(1,1) = pow(j)
         amat(2,1) = dpow(j)
         amat(1,2) = kvec(j)
         amat(2,2) = kder(j)

         y(1) = cu(2*j)
         y(2) = cup(2*j)
         call mbh2dsov_ge22cp(ierge,amat,x,y)  

         lmp(j) = -x(1)
         ymp(j) = x(2)

         y(1) = cu(2*j+1)
         y(2) = cup(2*j+1)
         call mbh2dsov_ge22cp(ierge,amat,x,y)  

         lmp(j) = lmp(j) + eye*x(1)
         ymp(j) = ymp(j) - eye*x(2)

      enddo
      

      return
      end

      subroutine mbh2dsov_dtocta_naive(beta,rad,rscale,rscalelap,u,up,
     1     npts,iloc,lloc,wsave,work)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     This is a wrapper routine. The real work is carried out by 
c     MBH2DSOV_DTOCTA1. 
c
c     INPUT
c
c     beta - REAL *8, the parameter of the modified biharmonic equation
c     rad - REAL *8, the radius of the circle on which u and up are 
c           sampled
c     rscale - REAL *8, the scaling parameter for the basis functions
c     u - REAL *8 array of length NPTS, the values of the field at 
c         equispaced points on a circle of radius rad
c     up - REAL *8 array of length NPTS, the values of the derivative
c          of the field in the radial direction at the same points as u
c     npts - INTEGER, as described above
c     wsave - REAL *8 array of length NPTS, the precomputed array 
c             returned by DFFTI from the FFT library FFTPACK
c     work - REAL *8 work array. recommended length of work array:
c            14*npts + 500
c
c     OUTPUT
c
c     iloc          : coefficients for difference-type functions
c     lloc            : coefficients for z^k functions
c    
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    
      implicit none
      complex *16 iloc(0:(npts+1)/2), lloc(0:(npts+1)/2)
      real *8 beta, rad, u(*), up(*), work(*)
      real *8 wsave(*), rscale, rscalelap
      integer npts
c     local variables
      integer idiffs, iders, ifjs, ifjder, iivec, iider, icu, icup
      integer iiscale
      integer ldiffs, lwfjs, lfjder, livec, l, lused, lcu, liscale

      l = (npts+1)/2

      ldiffs = l+1
      lwfjs = 2*(l + 5 + 4*l + 100)
      liscale = l + 5 + 4*l + 100
      lfjder = 2*(l+6)
      livec = l+6
      lcu = npts

      idiffs = 1
      iders = idiffs + ldiffs
      ifjs = iders + ldiffs
      ifjder = ifjs + lwfjs
      iivec = ifjder + lfjder
      iider = iivec + livec
      icu = iider + livec
      icup = icu + lcu
      iiscale = icup + lcu
      lused = iiscale + liscale

      call mbh2dsov_dtocta1_naive(beta,rad,rscale,rscalelap,u,up,npts,
     1     iloc,lloc,wsave,work(idiffs),work(iders),
     2     lwfjs,work(ifjs),work(ifjder),work(iiscale),
     3     work(iivec),work(iider),
     4     work(icu),work(icup))

      return
      end

      subroutine mbh2dsov_dtocta1_naive(beta,rad,rscale,rscalelap,u,
     1     up,npts,iloc,lloc,wsave,pow,dpow,lwfjs,fjs,fjder,
     2     iscale,ival,ider,cu,cup)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     see MBH2DSOV_DTOCTA for description
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      complex *16 fjs(0:1), fjder(0:1)
      complex *16 iloc(0:(npts+1)/2), lloc(0:(npts+1)/2)
      real *8 ival(0:1), ider(0:1), wsave(*), rscale, rscalelap
      real *8 beta, rad, u(*), up(*)
      real *8 pow(0:1), dpow(0:1), cu(*), cup(*)
      integer iscale(*)
      integer npts, lwfjs
c     local variables
      integer i, l, ifder, ifders, j, ntop, ier, ierge, n
      real *8 pih, norm, pi2
      real *8 a, b, c, d, alphan, betan, alphad, betad
      real *8 rtemp
      real *8 cu2j, cup2j, cu2jp1, cup2jp1, gamma
      real *8 amat(2,2), x(2), y(2)
      complex *16 z, eye, zk
      data eye /(0.0d0,1.0d0) /

      do i = 1,npts
         cu(i) = u(i)
         cup(i) = up(i)
      enddo

c     fourier transform u and up values

      call dfftf(npts,cu,wsave)
      call dfftf(npts,cup,wsave)

      norm = 2.0d0/npts

      cu(1) = cu(1)/npts
      cup(1) = cup(1)/npts

      do i = 2,npts
         cu(i) = cu(i)*norm
         cup(i) = cup(i)*norm
      enddo

c     obtain values and derivatives of basis functions
c     for radius of circle

      l = (npts+1)/2

      zk = eye*beta
      z = zk*rad
      ifder=1
      call jfuns2d(ier,l+1,z,rscale,fjs,ifder,fjder,
     1           lwfjs,iscale,ntop)

      write(*,*) 'ier = ', ier
      
c     convert Bessel J_n to modified Bessel I_n
      do n = 0,l,4
         ival(n) = dreal(fjs(n))*0.25d0
         ival(n+1) = dimag(fjs(n+1))*0.25d0
         ival(n+2) = -dreal(fjs(n+2))*0.25d0
         ival(n+3) = -dimag(fjs(n+3))*0.25d0
         ider(n) = -dimag(fjder(n))*beta*0.25d0
         ider(n+1) = dreal(fjder(n+1))*beta*0.25d0
         ider(n+2) = dimag(fjder(n+2))*beta*0.25d0
         ider(n+3) = -dreal(fjder(n+3))*beta*0.25d0  
      enddo
      

      call mbh2dsov_rk(pow,dpow,rad,rscalelap,l)
      
      pi2 = 8.0d0*datan(1.0d0)
      
      do i = 0,l
         pow(i) = pow(i)/pi2
         dpow(i) = dpow(i)/pi2
      enddo


c     compute zero-th mode coefficients

      amat(1,1) = ival(0)
      amat(2,1) = ider(0)
      amat(1,2) = pow(0)
      amat(2,2) = dpow(0)

      y(1) = cu(1)
      y(2) = cup(1)
      call mbh2dsov_ge22cp(ierge,amat,x,y)  

      iloc(0) = x(1)
      lloc(0) = -x(2)

      do j = 1,l-1

c     compute coefficients for higher modes
c     these are separated into cos(j theta) and 
c     sin(j theta) parts by DFFTF

         amat(1,1) = ival(j)
         amat(2,1) = ider(j)
         amat(1,2) = pow(j)
         amat(2,2) = dpow(j)

         y(1) = cu(2*j)
         y(2) = cup(2*j)

         call mbh2d_ge22cp(ierge,amat,x,y)  

         iloc(j) = x(1)
         lloc(j) = -x(2)

         y(1) = cu(2*j+1)
         y(2) = cup(2*j+1)

         call mbh2d_ge22cp(ierge,amat,x,y)  

         iloc(j) = iloc(j) - eye*x(1)
         lloc(j) = lloc(j) - eye*x(2)

      enddo

      return
      end

      subroutine mbh2dsov_matsmp(beta,rad,rscale,npts,amats,
     1     work)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     This is a wrapper routine. The real work is carried out by 
c     MBH2DSOV_MATS1. 
c
c     INPUT
c
c     beta - REAL *8, the parameter of the modified biharmonic equation
c     rad - REAL *8, the radius of the circle on which u and up are 
c           sampled
c     rscale - REAL *8, the scaling parameter for the basis functions
c     u - REAL *8 array of length NPTS, the values of the field at 
c         equispaced points on a circle of radius rad
c     up - REAL *8 array of length NPTS, the values of the derivative
c          of the field in the radial direction at the same points as u
c     npts - INTEGER, as described above
c     wsave - REAL *8 array of length NPTS, the precomputed array 
c             returned by DFFTI from the FFT library FFTPACK
c     work - REAL *8 work array. recommended length of work array:
c            14*npts + 500
c
c     OUTPUT
c
c     mbhmp, ymp - REAL *8 arrays of length npts which 
c                        give the coefficients for the basis 
c                        functions which approximate this field 
c                        outside the disk
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    
      implicit none
      real *8 amats(2,2,0:(npts+1)/2)
      real *8 beta, rad, work(*)
      real *8 rscale
      integer npts
c     local variables
      integer idiffs, iders, ihvec, ihder, ikvec, ikder, icu, icup
      integer ldiffs, lhvec, lkvec, lcu, l, lused

      l = (npts+1)/2

      ldiffs = l+1
      lhvec = 2*(l+6)
      lkvec = l+6
      lcu = npts

      idiffs = 1
      iders = idiffs + ldiffs
      ihvec = iders + ldiffs
      ihder = ihvec + lhvec
      ikvec = ihder + lhvec
      ikder = ikvec + lkvec
      icu = ikder + lkvec
      icup = icu + lcu
      lused = icup + lcu

      call mbh2dsov_matsmp1(beta,rad,rscale,npts,
     1     amats,work(idiffs),work(iders),
     2     work(ihvec),work(ihder),work(ikvec),work(ikder))


      return
      end
      
      subroutine mbh2dsov_matsmp1(beta,rad,rscale,npts,
     1     amats,diffs,ders,hvec,hder,kvec,kder)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     see MBH2DSOV_MATS for description
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      complex *16 hvec(0:1), hder(0:1)
      real *8 amats(2,2,0:(npts+1)/2)
      real *8 kvec(0:1), kder(0:1), rscale
      real *8 beta, rad
      real *8 diffs(0:1), ders(0:1)
      integer npts
c     local variables
      integer i, l, ifder, ifders, j, ierge
      real *8 pih, norm
      real *8 a, b, c, d, alphan, betan, alphad, betad
      real *8 cu2j, cup2j, cu2jp1, cup2jp1, gamma
      real *8 amat(2,2), x(2), y(2)
      complex *16 z, eye
      data eye /(0.0d0,1.0d0) /

      norm = 2.0d0/npts

      l = (npts+1)/2

      z = eye*rad*beta
      ifder = 1
      call h2dall(l+5,z,rscale,hvec,ifder,hder)

      pih = 2.0d0*atan(1.0d0)

      do j = 0,l+1,4
         kvec(j) = -dimag(hvec(j))*pih
         kvec(j+1) = -dreal(hvec(j+1))*pih
         kvec(j+2) = dimag(hvec(j+2))*pih
         kvec(j+3) = dreal(hvec(j+3))*pih
         kder(j) = -dreal(hder(j))*pih*beta
         kder(j+1) = dimag(hder(j+1))*pih*beta
         kder(j+2) = dreal(hder(j+2))*pih*beta
         kder(j+3) = -dimag(hder(j+3))*pih*beta
      enddo

      ifders = 1
      call diffslogbk_fast(rad,beta,rscale,diffs,ifders,ders,kvec,l)

      amats(1,1,0) = diffs(0)
      amats(2,1,0) = ders(0)
      amats(1,2,0) = kvec(0)
      amats(2,2,0) = kder(0)
      
      do j = 1,l-1

         amats(1,1,j) = diffs(j)
         amats(2,1,j) = ders(j)
         amats(1,2,j) = kvec(j)
         amats(2,2,j) = kder(j)

      enddo
      

      return
      end

      subroutine mbh2dsov_matsta(beta,rad,rscale,npts,
     1     amats,work)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     This is a wrapper routine. The real work is carried out by 
c     MBH2DSOV_MATSTA1. 
c
c     INPUT
c
c     beta - REAL *8, the parameter of the modified biharmonic equation
c     rad - REAL *8, the radius of the circle on which u and up are 
c           sampled
c     rscale - REAL *8, the scaling parameter for the basis functions
c     u - REAL *8 array of length NPTS, the values of the field at 
c         equispaced points on a circle of radius rad
c     up - REAL *8 array of length NPTS, the values of the derivative
c          of the field in the radial direction at the same points as u
c     npts - INTEGER, as described above
c     wsave - REAL *8 array of length NPTS, the precomputed array 
c             returned by DFFTI from the FFT library FFTPACK
c     work - REAL *8 work array. recommended length of work array:
c            14*npts + 500
c
c     OUTPUT
c
c     mbhloc          : coefficients for difference-type functions
c     lloc            : coefficients for z^k functions
c    
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    
      implicit none
      real *8 amats(2,2,(npts+1)/2)
      real *8 beta, rad, work(*)
      real *8 rscale
      integer npts
c     local variables
      integer idiffs, iders, ifjs, ifjder, iivec, iider, icu, icup
      integer iiscale
      integer ldiffs, lwfjs, lfjder, livec, l, lused, lcu, liscale

      l = (npts+1)/2

      ldiffs = l+1
      lwfjs = 2*(l + 5 + 4*l + 100)
      liscale = l + 5 + 4*l + 100
      lfjder = 2*(l+6)
      livec = l+6
      lcu = npts

      idiffs = 1
      iders = idiffs + ldiffs
      ifjs = iders + ldiffs
      ifjder = ifjs + lwfjs
      iivec = ifjder + lfjder
      iider = iivec + livec
      icu = iider + livec
      icup = icu + lcu
      iiscale = icup + lcu
      lused = iiscale + liscale

      call mbh2dsov_matsta1(beta,rad,rscale,npts,
     1     amats,work(idiffs),work(iders),
     2     lwfjs,work(ifjs),work(ifjder),work(iiscale),
     3     work(iivec),work(iider))

      return
      end

      subroutine mbh2dsov_matsta1(beta,rad,rscale,npts,
     1     amats,diffs,ders,lwfjs,fjs,fjder,iscale,pow,dpow)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     see MBH2DSOV_MATSTA for description
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      complex *16 fjs(0:1), fjder(0:1)
      real *8 amats(2,2,0:(npts+1)/2)
      real *8 pow(0:1), dpow(0:1), rscale
      real *8 beta, rad
      real *8 diffs(0:1), ders(0:1)
      integer iscale(*)
      integer npts, lwfjs
c     local variables
      integer i, l, ifder, ifders, j, ntop, ier, ierge
      real *8 pih, norm
      real *8 a, b, c, d, alphan, betan, alphad, betad
      real *8 rtemp
      real *8 cu2j, cup2j, cu2jp1, cup2jp1, gamma
      real *8 amat(2,2), x(2), y(2)
      complex *16 z, eye
      data eye /(0.0d0,1.0d0) /

c     obtain values and derivatives of basis functions
c     for radius of circle

      l = (npts+1)/2

      ifders = 1
      call diffszkik_fast(rad,beta,rscale,diffs,ifders,ders,pow,l)

      call mbh2dsov_rks(pow,dpow,rad,beta,rscale,l)


c     compute zero-th mode coefficients

      amats(1,1,0) = diffs(0)
      amats(2,1,0) = ders(0)
      amats(1,2,0) = pow(0)
      amats(2,2,0) = dpow(0)

      do j = 1,l-1

c     compute coefficients for higher modes
c     these are separated into cos(j theta) and 
c     sin(j theta) parts by DFFTF

         amats(1,1,j) = diffs(j)
         amats(2,1,j) = ders(j)
         amats(1,2,j) = pow(j)
         amats(2,2,j) = dpow(j)

      enddo

      return
      end
      
      subroutine mbh2dsov_matsmp_naive(beta,rad,rscale,rscalelap,
     1     npts,amats,work)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     This is a wrapper routine. The real work is carried out by 
c     MBH2DSOV_MATS1. 
c
c     INPUT
c
c     beta - REAL *8, the parameter of the modified biharmonic equation
c     rad - REAL *8, the radius of the circle on which u and up are 
c           sampled
c     rscale - REAL *8, the scaling parameter for the basis functions
c     u - REAL *8 array of length NPTS, the values of the field at 
c         equispaced points on a circle of radius rad
c     up - REAL *8 array of length NPTS, the values of the derivative
c          of the field in the radial direction at the same points as u
c     npts - INTEGER, as described above
c     wsave - REAL *8 array of length NPTS, the precomputed array 
c             returned by DFFTI from the FFT library FFTPACK
c     work - REAL *8 work array. recommended length of work array:
c            14*npts + 500
c
c     OUTPUT
c
c     mbhmp, ymp - REAL *8 arrays of length npts which 
c                        give the coefficients for the basis 
c                        functions which approximate this field 
c                        outside the disk
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    
      implicit none
      real *8 amats(2,2,0:(npts+1)/2)
      real *8 beta, rad, work(*)
      real *8 rscale, rscalelap
      integer npts
c     local variables
      integer idiffs, iders, ihvec, ihder, ikvec, ikder, icu, icup
      integer ldiffs, lhvec, lkvec, lcu, l, lused

      l = (npts+1)/2

      ldiffs = l+1
      lhvec = 2*(l+6)
      lkvec = l+6
      lcu = npts

      idiffs = 1
      iders = idiffs + ldiffs
      ihvec = iders + ldiffs
      ihder = ihvec + lhvec
      ikvec = ihder + lhvec
      ikder = ikvec + lkvec
      icu = ikder + lkvec
      icup = icu + lcu
      lused = icup + lcu

      call mbh2dsov_matsmp1_naive(beta,rad,rscale,rscalelap,npts,
     1     amats,work(idiffs),work(iders),
     2     work(ihvec),work(ihder),work(ikvec),work(ikder))


      return
      end
      
      subroutine mbh2dsov_matsmp1_naive(beta,rad,rscale,rscalelap,
     1     npts,amats,pow,dpow,hvec,hder,kvec,kder)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     see MBH2DSOV_MATS for description
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      complex *16 hvec(0:1), hder(0:1)
      real *8 amats(2,2,0:(npts+1)/2)
      real *8 kvec(0:1), kder(0:1), rscale, rscalelap
      real *8 beta, rad
      real *8 pow(0:1), dpow(0:1)
      integer npts
c     local variables
      integer i, l, ifder, ifders, j, ierge
      real *8 pih, norm, pi2
      real *8 a, b, c, d, alphan, betan, alphad, betad
      real *8 cu2j, cup2j, cu2jp1, cup2jp1, gamma
      real *8 amat(2,2), x(2), y(2)
      complex *16 z, eye
      data eye /(0.0d0,1.0d0) /

      l = (npts+1)/2

      z = eye*rad*beta
      ifder = 1
      call h2dall(l+5,z,rscale,hvec,ifder,hder)

      pih = 2.0d0*atan(1.0d0)

      do j = 0,l+1,4
         kvec(j) = -dimag(hvec(j))*0.25d0
         kvec(j+1) = -dreal(hvec(j+1))*0.25d0
         kvec(j+2) = dimag(hvec(j+2))*0.25d0
         kvec(j+3) = dreal(hvec(j+3))*0.25d0
         kder(j) = -dreal(hder(j))*beta*0.25d0
         kder(j+1) = dimag(hder(j+1))*beta*0.25d0
         kder(j+2) = dreal(hder(j+2))*beta*0.25d0
         kder(j+3) = -dimag(hder(j+3))*beta*0.25d0
      enddo

      call mbh2dsov_rkinv(pow,dpow,rad,rscalelap,l)

      pi2 = 8.0d0*datan(1.0d0)
      
      do j = 0,l
         pow(j) = pow(j)/pi2
         dpow(j) = dpow(j)/pi2
      enddo
      
      amats(1,1,0) = pow(0)
      amats(2,1,0) = dpow(0)
      amats(1,2,0) = kvec(0)
      amats(2,2,0) = kder(0)

      do j = 1,l-1

         amats(1,1,j) = pow(j)
         amats(2,1,j) = dpow(j)
         amats(1,2,j) = kvec(j)
         amats(2,2,j) = kder(j)

      enddo
      

      return
      end

      subroutine mbh2dsov_matsta_naive(beta,rad,rscale,rscalelap,
     1     npts,amats,work)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     This is a wrapper routine. The real work is carried out by 
c     MBH2DSOV_MATSTA1. 
c
c     INPUT
c
c     beta - REAL *8, the parameter of the modified biharmonic equation
c     rad - REAL *8, the radius of the circle on which u and up are 
c           sampled
c     rscale - REAL *8, the scaling parameter for the basis functions
c     u - REAL *8 array of length NPTS, the values of the field at 
c         equispaced points on a circle of radius rad
c     up - REAL *8 array of length NPTS, the values of the derivative
c          of the field in the radial direction at the same points as u
c     npts - INTEGER, as described above
c     wsave - REAL *8 array of length NPTS, the precomputed array 
c             returned by DFFTI from the FFT library FFTPACK
c     work - REAL *8 work array. recommended length of work array:
c            14*npts + 500
c
c     OUTPUT
c
c     iloc          : coefficients for difference-type functions
c     lloc            : coefficients for z^k functions
c    
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    
      implicit none
      real *8 amats(2,2,0:(npts+1)/2)
      real *8 beta, rad, work(*)
      real *8 rscale, rscalelap
      integer npts
c     local variables
      integer idiffs, iders, ifjs, ifjder, iivec, iider, icu, icup
      integer iiscale
      integer ldiffs, lwfjs, lfjder, livec, l, lused, lcu, liscale

      l = (npts+1)/2

      ldiffs = l+1
      lwfjs = 2*(l + 5 + 4*l + 100)
      liscale = l + 5 + 4*l + 100
      lfjder = 2*(l+6)
      livec = l+6
      lcu = npts

      idiffs = 1
      iders = idiffs + ldiffs
      ifjs = iders + ldiffs
      ifjder = ifjs + lwfjs
      iivec = ifjder + lfjder
      iider = iivec + livec
      icu = iider + livec
      icup = icu + lcu
      iiscale = icup + lcu
      lused = iiscale + liscale

      call mbh2dsov_matsta1_naive(beta,rad,rscale,rscalelap,npts,
     1     amats,work(idiffs),work(iders),
     2     lwfjs,work(ifjs),work(ifjder),work(iiscale),
     3     work(iivec),work(iider))

      return
      end

      subroutine mbh2dsov_matsta1_naive(beta,rad,rscale,rscalelap,
     1     npts,amats,pow,dpow,lwfjs,fjs,fjder,iscale,ival,ider)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     see MBH2DSOV_MATSTA for description
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      complex *16 fjs(0:1), fjder(0:1)
      real *8 amats(2,2,0:(npts+1)/2)
      real *8 ival(0:1), ider(0:1), rscale, rscalelap
      real *8 beta, rad
      real *8 pow(0:1), dpow(0:1)
      integer iscale(*)
      integer npts, lwfjs
c     local variables
      integer i, l, ifder, ifders, j, ntop, ier, ierge, n
      real *8 pih, norm, pi2
      real *8 a, b, c, d, alphan, betan, alphad, betad
      real *8 rtemp
      real *8 cu2j, cup2j, cu2jp1, cup2jp1, gamma
      real *8 amat(2,2), x(2), y(2)
      complex *16 z, eye, zk
      data eye /(0.0d0,1.0d0) /

c     obtain values and derivatives of basis functions
c     for radius of circle

      l = (npts+1)/2

      zk = eye*beta
      z = zk*rad
      ifder=1
      call jfuns2d(ier,l+1,z,rscale,fjs,ifder,fjder,
     1           lwfjs,iscale,ntop)

c     convert Bessel J_n to modified Bessel I_n
      do n = 0,l,4
         ival(n) = dreal(fjs(n))*0.25d0
         ival(n+1) = dimag(fjs(n+1))*0.25d0
         ival(n+2) = -dreal(fjs(n+2))*0.25d0
         ival(n+3) = -dimag(fjs(n+3))*0.25d0
         ider(n) = -dimag(fjder(n))*beta*0.25d0
         ider(n+1) = dreal(fjder(n+1))*beta*0.25d0
         ider(n+2) = dimag(fjder(n+2))*beta*0.25d0
         ider(n+3) = -dreal(fjder(n+3))*beta*0.25d0  
      enddo
      

      call mbh2dsov_rk(pow,dpow,rad,rscalelap,l)
      
      pi2 = 8.0d0*datan(1.0d0)
      
      do i = 0,l
         pow(i) = pow(i)/pi2
         dpow(i) = dpow(i)/pi2
      enddo


c     compute zero-th mode coefficients

      amats(1,1,0) = ival(0)
      amats(2,1,0) = ider(0)
      amats(1,2,0) = pow(0)
      amats(2,2,0) = dpow(0)

      do j = 1,l-1

c     compute coefficients for higher modes
c     these are separated into cos(j theta) and 
c     sin(j theta) parts by DFFTF

         amats(1,1,j) = ival(j)
         amats(2,1,j) = ider(j)
         amats(1,2,j) = pow(j)
         amats(2,2,j) = dpow(j)

      enddo

      return
      end

      subroutine mbh2dsov_circpts(z,npts)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     INPUT:
c
c     npts       : the number of points on circle
c
c     OUTPUT:
c
c     z(2,npts)  : equally spaced points on unit circle
c                  z(1,i) = cos((i-1)*2*pi/npts)
c                  z(2,i) = sin((i-1)*2*pi/npts)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global 
      integer npts
      real *8 z(2,*)
c     local
      real *8 pi2
      integer i

      pi2 = 8.0d0*datan(1.0d0)
      
      do i = 1,npts
         z(1,i) = dcos((i-1)*pi2/npts)
         z(2,i) = dsin((i-1)*pi2/npts)
      enddo

      return
      end


      subroutine mbh2dsov_ge22cp(ier,amat,x,y)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     gaussian elimination with complete pivoting for 2x2 matrices
c     
c     solves:   amat*x = y
c
c     if the matrix is singular, the least squares solution is returned
c
c     INPUT:
c
c     a(2,2)        : system matrix a
c     y(2)          : right-hand side
c
c     OUTPUT:
c
c     x(2)          : solution
c     ier           : flag, IER = 0 means success
c                     IER = 1 means a is very small in norm
c                     IER = 2 means a is nearly not invertible
c                     IER = 3 means a is not invertible
c                      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global
      real *8 amat(2,2), x(2), y(2)
      integer ier
c     local
      integer imax, jmax, i, j
      real *8 a, b, c, d, e, f, g, h, gamma, dtemp, a2, b2, coa
      real *8 epssmall
      parameter (epssmall = 1.0d-200)
      
      imax = 1
      jmax = 1
      a = amat(1,1)
      if (dabs(amat(2,1)) .gt. dabs(a)) then
         imax = 2
         jmax = 1
         a = amat(2,1)
      endif
      if (dabs(amat(1,2)) .gt. dabs(a)) then
         imax = 1
         jmax = 2
         a = amat(1,2)
      endif
      if (dabs(amat(2,2)) .gt. dabs(a)) then
         imax = 2
         jmax = 2
         a = amat(2,2)
      endif

c     special case, norm of amat is tiny

      if (dabs(a) .eq. 0.0d0) then
         ier = 3
         x(1) = 0.0d0
         x(2) = 0.0d0
         return
      endif

      if (dabs(a) .lt. epssmall) ier = 1

      i = mod(imax,2)+1
      j = mod(jmax,2)+1

c     grab other entries in pivoted matrix (a,b;c,d)

      b = amat(imax,j)
      c = amat(i,jmax)
      d = amat(i,j)

c     grab pivoted right hand side

      g = y(imax)
      h = y(i)

c     perform elimination

      dtemp = d-b*c/a

c     special case, nearly not invertible or not invertible
      

c     not invertible, return least squares solution
      if (dabs(dtemp) .eq. 0.0d0) then
         ier = 3
         a2 = a**2
         b2 = b**2
         coa = c/a
         gamma = (g+h*coa)/(a+c*coa)
         e = gamma*a2/(a2+b2)
         f = gamma*b/a
         x(jmax) = e
         x(j) = f
         return
      endif

      if (dabs(dtemp) .lt. epssmall) ier = 2

      h = h-g*c/a

      f = h/dtemp
      e = (g-b*f)/a

c     copy solution to correct entries

      x(jmax) = e
      x(j) = f

      return
      end


      subroutine mbh2dsov_rks(pow,dpow,r,beta,rscale,nterms)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     INPUT:
c
c     r               : radius
c     rscale          : scaling factor
c     nterms          : number of terms in the expansion
c     
c     OUTPUT:
c
c     pow(0:nterms)   : pow(i) = (r/rscale)^i
c     dpow(0:nterms)  : dpow(i) = i*(r/rscale)^(i-1)/rscale
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global
      real *8 pow(0:nterms), dpow(0:nterms), r, rscale, beta
      integer nterms
c     local
      real *8 dtemp1, dtemp2
      integer i

      dtemp2 = r*beta/(rscale)
      dtemp1 = 1.0d0

      pow(0) = 1.0d0
      dpow(0) = 0.0d0

      do i = 1,nterms
         dpow(i) = i*beta*dtemp1/(rscale)
         dtemp1 = dtemp1*dtemp2
         if (dtemp1 .gt. 1.0d250) dtemp1 = 0.0d0         
         pow(i) = dtemp1
      enddo

      return
      end

      subroutine mbh2dsov_rk(pow,dpow,r,rscale,nterms)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     INPUT:
c
c     r               : radius
c     rscale          : scaling factor
c     nterms          : number of terms in the expansion
c     
c     OUTPUT:
c
c     pow(0:nterms)   : pow(i) = (r/rscale)^i
c     dpow(0:nterms)  : dpow(i) = i*(r/rscale)^(i-1)/rscale
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global
      real *8 pow(0:nterms), dpow(0:nterms), r, rscale
      integer nterms
c     local
      real *8 dtemp1, dtemp2
      integer i

      dtemp2 = r/(rscale)
      dtemp1 = 1.0d0

      pow(0) = 1.0d0
      dpow(0) = 0.0d0

      do i = 1,nterms
         dpow(i) = i*dtemp1/(rscale)
         dtemp1 = dtemp1*dtemp2
         if (dtemp1 .gt. 1.0d250) dtemp1 = 0.0d0         
         pow(i) = dtemp1
      enddo

      return
      end

      subroutine mbh2dsov_rkinv(pow,dpow,r,rscale,nterms)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     INPUT:
c
c     r               : radius
c     rscale          : scaling factor
c     nterms          : number of terms in the expansion
c     
c     OUTPUT:
c
c     pow(0:nterms)   : pow(i) = (r/rscale)^-i
c     dpow(0:nterms)  : dpow(i) = -i*(r/rscale)^(-i-1)/rscale
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global
      real *8 pow(0:nterms), dpow(0:nterms), r, rscale
      integer nterms
c     local
      real *8 dtemp1, dtemp2
      integer i

      dtemp2 = rscale/(r)
      dtemp1 = dtemp2

      pow(0) = dlog(r)
      dpow(0) = 1/r

      do i = 1,nterms
         pow(i) = dtemp1
         dtemp1 = dtemp1*dtemp2
         if (dtemp1 .gt. 1.0d250) dtemp1 = 0.0d0         
         dpow(i) = -i*dtemp1/(rscale)
      enddo

      return
      end
