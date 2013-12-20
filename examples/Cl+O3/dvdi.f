
      subroutine dvdi(j,rinv)

c  evaluate the energy and gradients

      include 'traj.inc'

      dimension dvt(maxf,nbondm),ei(maxf),dei(maxf,nintm),vt1(maxf),
     .        rinv(nbondm),sumdvt(nmax,nbondm),vt(maxf),
     .        z(maxf,nintm),dtdr(maxf,nbondm)

c  zero those who need it

      en(j)=0.d0
      totsum=0.d0

      do i=1,nbond
         dvr(j,i)=0.d0
         sumdvt(j,i)=0.d0
         dvt(maxf,i)=0.0d0
      enddo

      do k=1,nforc
         vt(k)=0.d0
         vt1(k)=0.d0
      enddo
      vt(maxf)=0.0d0

c  we will need to have some quantities for the gradients
c  of the weights
c  Note an approx here, due to neglect of small weights

      ip2=ipow*2

c  calc the raw wts and totsum

      do i=1,nbond
      do k=1,nforc
c         temp= 1.d0/r(j,i)-1.d0/rda(mda(k),ip(i,nda(k)))
         temp= rinv(i)-1.0d0/rda(mda(k),ip(i,nda(k)))
c        temp= r(j,i)-rda(mda(k),ip(i,nda(k)))
         vt1(k)= temp*temp + vt1(k)
      enddo
      enddo

      do k=1,nforc
         vt(k)=1.d0/vt1(k)**ipow
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
         temp=ip2*wt(k)/vt1(k)
         do i=1,nbond
            dvt(k,i)=temp*(rinv(i)-1.0/rda(mda(k),ip(i,nda(k))))*rinv(i)**2
         enddo
      enddo

c      do k=1,nforc
c        temp=-ip2*wt(k)/vt1(k)
c        do i=1,nbond
c           dvt(k,i)=temp*(r(j,i)-1.0/rda(mda(k),ip(i,nda(k))))
c        enddo
c      enddo

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

c      do i=1,nbond
c         rinv(i)=1.d0/r(j,i)
c      enddo

c  !ocl vi(k)
      do k=1,nforc
      do i=1,nint
         z(k,i)=0.d0
         do l=1,nbond
            z(k,i) = ut(mda(k),i,ip(l,nda(k)))*rinv(l) + z(k,i)
         enddo
c         z(k,i)=tanh(z(k,i)-zda(mda(k),i))
          z(k,i)=(z(k,i)-zda(mda(k),i))
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

      do k=1,nforc
         ei(k)=v0(mda(k))
      enddo

c  linear term in energy plus diagonal part of quadratic term
c  first term in derivative of taylor poly wrt local coords

      do i=1,nint
      do k=1,nforc

	 ei(k) = z(k,i)*v1(mda(k),i) 
     .   + 0.5d0*z(k,i)**2*v2(mda(k),i) + ei(k)

	 dei(k,i) = z(k,i)*v2(mda(k),i) + v1(mda(k),i)

      enddo
      enddo

c  the O(n^2) loop doesn't apply because we diagonalised d2v!!!

c  transform derivatives of taylor poly wrt normal local coords into 
c  derivs wrt bondlengths NOT inverses

c  !ocl vi(k)
      do i=1,nbond
         t=-rinv(i)**2
         do k=1,nforc
            dtdr(k,i)=0.d0
            do l=1,nint

c                sechsq=1.d0
c this was commented originally
c               sechsq=1.d0 - z(k,l)**2

c            dtdr(k,i)=sechsq*dei(k,l)*ut(mda(k),l,ip(i,nda(k))) + dtdr(k,i)

             dtdr(k,i)=dei(k,l)*ut(mda(k),l,ip(i,nda(k))) + dtdr(k,i)

            enddo
            dtdr(k,i)=dtdr(k,i)*t
         enddo
      enddo

c  scale by relative weights, and their derivatives

      do k=1,nforc
     
c        write(6,*) k,ei(k),wt(k)
         en(j)=ei(k)*wt(k)+en(j)
      enddo

      do i=1,nbond
      do k=1,nforc
         dvr(j,i)=ei(k)*dvt(k,i)+wt(k)*dtdr(k,i)+dvr(j,i)
      enddo
      enddo

c  smoothing attempt: <\delta T^2>

      rms=0.d0
      grms=0.d0
      do k=1,nforc
         rms=wt(k)*(ei(k)-en(j))**2 + rms
         do i=1,nbond
         grms=grms+wt(k)*(dtdr(k,i)-dvr(j,i))**2
         enddo
      enddo
      rms=sqrt(rms)
      totrms=totrms+rms

c  sum over interactions for each traj

c     do k=1,nforc
c        en(j)=en(j)+ei(k)
c     enddo

c     do i=1,nbond
c     do k=1,nforc
c        dvr(j,i)=dvr(j,i)+dtdr(k,i)
c     enddo
c     enddo

c  look out for loony geometries

      if (nlowstop.eq.0) goto 1975

      ect=vmin-0.0005d0

      if (en(j).lt.ect) then
         if (nloop.gt.1) then
            do n=1,natom
               write(8,*)(c(j,n,i),i=1,3)
            enddo
            write(8,*)en(j),rms,totsum
         endif
         nfin(j)=3
         nstop=2
      endif

1975  continue

      return
      end

