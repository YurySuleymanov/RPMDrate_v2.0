
      subroutine dvdir(j,rinv)
c      subroutine dvdir(j)

c  evaluate the energy and gradients

      include 'traj.inc'

      dimension dvt(maxf,nbondm),ei(maxf),dei(maxf,nintm),vt1(maxf),
     .         rinv(nbondm),sumdvt(nmax,nbondm),vt(maxf),
     .         z(maxf,nintm),dtdr(maxf,nbondm),tsb(nbondm)


      dimension zed2(maxf), zedp(maxf), zed3(maxf), vtone(maxf)

c this version of dvdi uses a new two part primitive weight with a
c confidence "radius" rads(k) or sumb for each data point.


c  zero those who need it

      en(j)=0.d0
      totsum=0.d0


c      do i=1,nbond
c         dvr(j,i)=0.d0
c         sumdvt(j,i)=0.d0
c         dvt(maxf,i)=0.0d0
c         rinv(i)=1.d0/r(j,i)
c      enddo

      do k=1,nforc
         vt(k)=0.d0
         vt1(k)=0.d0
      enddo

      ip2=ipow*ipow

c  calc the raw wts and totsum

      do i=1,nbond
	  dvr(j,i)=0.d0
          sumdvt(j,i)=0.d0
	  dvt(maxf,i)=0.0d0
      do k=1,nforc
         i1=mda(k)
         i2=ip(i,nda(k))
         temp= rinv(i)-rda(i1,i2)
         vt1(k)= temp*temp/sumb(i1,i2) + vt1(k)
c         vtone(k)=temp*temp + vtone(k)
      enddo
      enddo

      do k=1,nforc
         zedp(k)=vt1(k)**ipow
         zed2(k)=vt1(k)**lowp
c         zed3(k)=vtone(k)**ipow
         vt(k)=1.d0/(zed2(k)+zedp(k))
         totsum=totsum+vt(k)
      enddo

c      do k=1,nforc
c         totsum=totsum+vt(k)
c      enddo

c  calculate the relative weights
  
c      do k=1,nforc
c         wt(k)=vt(k)/totsum
c      enddo

c  calc the derivatives of the raw weights v

      do k=1,nforc
         wt(k)=vt(k)/totsum
         temp=2.d0*wt(k)*vt(k)*(lowp*zed2(k)+ipow*zedp(k))/vt1(k)
         i1=mda(k)
         do i=1,nbond
         i2=ip(i,nda(k))
            dvt(k,i)=temp*(rinv(i)-rda(i1,i2))*rinv(i)*rinv(i)
     .               /sumb(i1,i2)
         enddo
      enddo

c dvt is actually the derivative of the primitive weight, divided by
c the sum of the primitive weights


      do i=1,nbond
      do k=1,nforc
         sumdvt(j,i)=sumdvt(j,i)+dvt(k,i)
      enddo
      enddo

c  calc the derivatives of the relative weights w

      do i=1,nbond
      do k=1,nforc
         dvt(k,i)= -wt(k)*sumdvt(j,i) + dvt(k,i)
      enddo
      enddo

c  check that the sum over weights vanishes

c      write(6,*)
c      write(6,*) '  weights : '
c      test=0.d0
c      do k=1,nforc
c         test=wt(k)+test
c         write(6,*) k,wt(k)
c      enddo
c      write(6,*) ' sum of weights = ',test 

c      do i=1,nbond
c         sumwt=0.d0
c         do k=1,nforc
c            sumwt=sumwt+dvt(k,i)
c         enddo
c         write(6,*)i,sumwt
c      enddo

c  construct the normal local coords for each neighbour
c  recalling that intern calculates bond lengths NOT inverses

      do k=1,nforc
      k1=mda(k)
      do i=1,nint
         z(k,i)=0.d0
         do l=1,nbond
             l1=ip(l,nda(k))
             z(k,i) = ut(k1,i,l1)*rinv(l) + z(k,i)
         enddo
             z(k,i)=z(k,i)-zda(k1,i)
      enddo
      enddo

c  degub 

c     do k=1,nforc
c        write(6,*)
c        write(6,*) ' bond lengths ',k,imol(k)
c        write(6,1971) (r(imol(k),i),i=1,nbond)
c        write(6,*) ' normal local coords z(k,nint)',k
c        write(6,1971) (z(k,i),i=1,nint)
c        write(6,*) ' data local coords zda',k
c        write(6,1971) (zda(mda(k),i),i=1,nint)
c        write(6,*) ' data bondlengths rda',k
c        write(6,1971) (rda(mda(k),ip(i,nda(k))),i=1,nbond)
c     enddo
1971  format(2x,6g11.4)

c  calculate the potential energy: scalar term 

c      do k=1,nforc
c         ei(k)=v0(mda(k))
c      enddo

c  linear term in energy plus diagonal part of quadratic term
c  first term in derivative of taylor poly wrt local coords

      do k=1,nforc
         ei(k)=v0(mda(k))
         k1=mda(k)
      do i=1,nint
c	 ei(k) = z(k,i)*v1(k1,i) 
c     .   + 0.5d0*v2(k1,i)*z(k,i)*z(k,i) + ei(k)

         zv2=z(k,i)*v2(k1,i)

	 ei(k) = z(k,i)*(v1(k1,i)+0.5d0*zv2) + ei(k)
	 
         dei(k,i) = zv2 + v1(k1,i)

      enddo
      enddo

c  the O(n^2) loop doesn't apply because we diagonalised d2v!!!

c  transform derivatives of taylor poly wrt normal local coords into 
c  derivs wrt bondlengths NOT inverses

!ocl vi(k)
      do i=1,nbond
         t=-rinv(i)*rinv(i)
         do k=1,nforc
            dtdr(k,i)=0.d0
            k1=mda(k)
            i1=ip(i,nda(k))
            do l=1,nint

               dtdr(k,i)=dei(k,l)*ut(k1,l,i1) + dtdr(k,i)

            enddo
            dtdr(k,i)=dtdr(k,i)*t
         enddo
      enddo

c  sum over neighbours

      do k=1,nforc
c         write(46,*) ei(k),v0(mda(k))
c         write(46,*) v0(mda(k))
c         ei(k)=v0(mda(k))
         en(j)=ei(k)*wt(k)+en(j)
c         write(46,*) k,ei(k),wt(k)
c        write(46,*) (1.0d0/rda(mda(k),ip(i,nda(k))),i=1,nbond)
      enddo

      do i=1,nbond
      do k=1,nforc
         dvr(j,i)=ei(k)*dvt(k,i)+wt(k)*dtdr(k,i)+dvr(j,i)
      enddo
      enddo

c variance in the energy

      rms=0.d0

      do k=1,nforc
         rms=wt(k)*(ei(k)-en(j))*(ei(k)-en(j)) + rms
      enddo

      rms=dsqrt(rms)

c        write(6,*) rms

c  temporarily (10/2/98) evaluate the rms using only the one-part
c  weight function

c     totsum1=0.d0
c     do k=1,nforc
c     totsum1=1.d0/zed3(k)+totsum1
c     enddo

c     do k=1,nforc
c      rms=rms+(1.d0/zed3(k)/totsum1)*(ei(k)-en(j))**2
c     enddo

c     rms=dsqrt(rms)

c variance in the gradients

c     tdv=0.d0
c     do i=1,nbond
c        tvg(i)=0.d0
c        tdv=dvr(j,i)**2 + tdv
c     enddo

c      grms=0.d0
c      do k=1,nforc
c      do i=1,nbond
c         t=dtdr(k,i)-dvr(j,i)
c         t=t*t
c        tvg(i)=wt(k)*t + tvg(i)
c         grms=wt(k)*t + grms
c      enddo
c      enddo

c      grms=dsqrt(grms)
     

c  look out for loony geometries

      if (nlowstop.eq.0) goto 1975

      ect=vmin-0.0005d0

      if (en(j).lt.ect) then
         if (nloop.gt.1) then
            do n=1,natom
               write(8,*)(c(j,n,i),i=1,3)
            enddo
            write(8,*)en(j),rms,totsum
      sumr=0.0d0
      do n=1,nbond
      tsb(n)=0.d0
      enddo
      do k=1,nforc
      sumr=sumr+wt(k)*rads(mda(k))
      do n=1,nbond
      tsb(n)=tsb(n)+wt(k)*sumb(mda(k),ip(n,nda(k)))
      enddo
      enddo
      write(17,*)sumr,(tsb(n),n=1,nbond)

         endif
         nfin(j)=3
         nstop=2
      endif

1975  continue


      return
      end

