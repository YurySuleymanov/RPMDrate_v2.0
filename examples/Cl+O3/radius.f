
      subroutine radius

      include 'traj.inc'

c  this subroutine looks at the data set and finds the "distance"
c  of the closest data point to each data point. The radius of 
c  "confidence" for each data point is set to ntimes times this 
c  closest distance

      dimension dists(ndatam*ngroupm),
     .          smdis(ndatam*ngroupm),z(nintm)

      dimension vt(ndatam,ngroupm),vt1(ndatam,ngroupm),mdat(maxf),ndat(maxf)
     .       ,wdat(maxf),dwt(maxf,nbondm),sumdvt(nbondm),dei(nintm)

      dimension dtdr(nbondm),dvkr(nbondm),grad(nbondm)

c  nneigh is now the number of neighbours in the POT file
c  who are to lie within rads (this is how rads is determined).

c  check that nneigh is not too big

      call lisym(6,'-')
      write(6,*)
      write(6,*) ' nneigh = ',nneigh
      write(6,*)
!      write(6,*)' k rads(k)'
      call lisym(11,'-')
      write(11,*) ' ***  subroutine radius *** '
      write(11,*)
      write(11,*) ' nneigh = ',nneigh
      write(11,*)
      if (nneigh.gt.(ndatam*ngroup)) then
         write(11,*)'  nneigh has been set too large for the data set'
         stop
      endif

      do k=1,ndata

         dmin=1.d+06
         totsum=0.d0
         ic=0

         do m=1,ndata
            if (m.eq.k) go to 1
            do n=1,ngroup

               sum=0.d0

               do i=1,nbond
                  t=rda(k,i)-rda(m,ip(i,n))
                  t2=t*t
                  sum=t2+sum
               enddo

               ic=ic+1
               dists(ic)=sum

c              totsum=totsum+1.d0/sum**lowp

               if (sum.lt.dmin) then
                  dmin=sum
               endif

            enddo
1           continue
         enddo

c guard against identical data points

         if (dmin.lt.1.d-10) dmin=1.d-4

         dmin2=2.d0*dmin
50       continue

c find all the distances less than 2*dmin

         ndist=0

         do i=1,ic
            if (dists(i).le.dmin2) then
               ndist=ndist+1
               smdis(ndist)=dists(i)
            endif
         enddo

         if (ndist.lt.(nneigh+1)) then
            dmin2=dmin2*2.d0
c           write(11,500) k,ndist
            go to 50
         endif

500   format(2x,'increase dmin2 for k = ',i5,2x,i10)

c  we now have a set of at least nneigh small distances
c  we sort them

         call sort(ndist,smdis)

         rads(k)=smdis(nneigh)

c  temp patch 12/8/97

c         rads(k)=rads(k)/10.d0

c end of patch 12/8/97

c       write(6,3) k, rads(k)
3        format(2x,i5,g14.5)

      enddo

      avrads=0.d0
      do k=1,ndata
         avrads=avrads+rads(k)
      enddo
      avrads=avrads/float(ndata)

c temp patch 3/9/97 to put rads very small

c      do k=1,ndata
c       rads(k)=avrads/1000.d0
c      enddo

c fixrad ?   
c     av=0.04d0

c     do k=1,ndata
c        if (rads(k).lt.0.025d0) rads(k)=0.025d0
c        rads(k)=av
c     enddo

      write(11,2) avrads
      write(11,*)
2     format('  <rads(k)> = ',g12.5)
c     write(6,*)
c     write(6,2) avrads
c     stop
      return

c ---------------- now we interpolate -----------------------------

      call lisym(6,'-')
      write(6,*)
      write(6,*) ' nneigh = ',nneigh
      write(6,*)

c  build neighbour list for each data point
c  first calculate raw weights and totsum

      enmx=0.d0
      gdmx=0.d0
      averr=0.d0
      avgrad=0.d0
      avneigh=0.d0

      do k=104,ndata

         totsum=0.d0
         do m=1,ndata
            if (m.eq.k) goto 27
            do n=1,ngroup
         
               sum=0.d0

               do i=1,nbond
                  t=rda(k,i)-rda(m,ip(i,n))
                  sum=t*t + sum
               enddo
               sum=sum/rads(m)

               vt1(m,n)=sum

               vt(m,n)=1.d0/(sum**lowp + sum**ipow)

               totsum = totsum + vt(m,n)

            enddo
27          continue
         enddo

c    calculate neighbour list 

         nforc=0
         do m=1,ndata
            if (m.eq.k) goto 29
            do n=1,ngroup
 
               wt1=vt(m,n)/totsum

               if (wt1.gt.wtol) then

                  nforc=nforc+1
                  mdat(nforc)=m 
                  ndat(nforc)=n

               endif

            enddo
29          continue
         enddo

         avneigh=avneigh + nforc

c  recalculate totsum for neighbours

         totsum=0.d0
         do j=1,nforc
            totsum=totsum + vt(mdat(j),ndat(j))
         enddo

c  calculate relative weights 

         do j=1,nforc
            wdat(j) = vt(mdat(j),ndat(j))/totsum
         enddo

c  derivatives of relative weight wrt bond lengths
          
         do i=1,nbond
            sumdvt(i)=0.d0
         enddo
 
         do j=1,nforc
            tv=vt1(mdat(j),ndat(j))
            z1=tv**ipow
            z2=tv**lowp
            vtt=1.d0/(z1 + z2)
            t=2.d0*wdat(j)*vtt*(lowp*z2+ipow*z1)/tv/rads(mdat(j))
            do i=1,nbond
               dwt(j,i) = t*(rda(k,i)-rda(mdat(j),ip(i,ndat(j))))*rda(k,i)**2
               sumdvt(i) = sumdvt(i) + dwt(j,i)
            enddo
         enddo

         do j=1,nforc
         do i=1,nbond
            dwt(j,i)= -wdat(j)*sumdvt(i) + dwt(j,i)
         enddo
         enddo

c    now we have a neighbour list, calculate interpolated potential

         shep=0.d0

         do i=1,nbond
            dvkr(i)=0.d0
         enddo

         do j=1,nforc

c    local coords
               
            do i=1,nint
               z(i)=0.d0
               do l=1,nbond
c                 t = rda(k,l) - rda(mdat(j),ip(l,ndat(j)))
                  z(i) = ut(mdat(j),i,ip(l,ndat(j)))*rda(k,l) + z(i)
               enddo
               z(i) = z(i) - zda(mdat(j),i) 
            enddo

c    taylor polys and 1st derivatives thereof

            ei=v0(mdat(j))

            do i=1,nint
               ei = z(i)*v1(mdat(j),i) + 0.5d0*v2(mdat(j),i)*z(i)**2 + ei
               dei(i) = v1(mdat(j),i) + v2(mdat(j),i)*z(i)
            enddo

c    shephard interpolation of energy
   
            shep = wdat(j)*ei + shep

c    transformation of d(ei)/dz to d(ei)/dr

            do i=1,nbond
               dtdr(i)=0.d0
               do l=1,nint
                  dtdr(i) = dei(l)*ut(mdat(j),l,ip(i,ndat(j))) + dtdr(i)
               enddo 
               dtdr(i)=-dtdr(i)*rda(k,i)**2
            enddo

c    derivative of interpolated potential at rda(k)

            do i=1,nbond
               dvkr(i) = ei*dwt(j,i) + wdat(j)*dtdr(i) + dvkr(i)
            enddo

c  close loop over neighbours

         enddo

c  energy error at rda(k) - in kJ/mol

         err = 100.d0*dabs(v0(k)-shep) 
         averr = averr + err

         if (err.gt.enmx) enmx=err

c  transform dv/dz(k) -> dv/dr|_k

         do i=1,nbond
            grad(i)=0.d0
            do j=1,nint
               grad(i) = v1(k,j)*ut(k,j,i) + grad(i)
            enddo
            grad(i) = -grad(i)*rda(k,i)**2
         enddo

c  gradient error at rda(k)

         gg=0.d0
         t1=0.d0
         t2=0.d0
         do i=1,nbond
            t1=grad(i)**2 + t1
            t2=dvkr(i)**2 + t2
            gg=(grad(i)-dvkr(i))**2 + gg
         enddo
         gg=100.d0*2.d0*dsqrt(gg)/(dsqrt(t1) + dsqrt(t2))

         avgrad = avgrad + gg

         if (gg.gt.gdmx) gdmx=gg
            
c close outer loop over data points

c        write(6,99) k,err,gg
99       format(2x,i5,4x,g12.5,4x,g12.5)

      enddo

      denom=1.d0/float(ndata-1)
      averr=averr*denom
      avgrad=avgrad*denom
      avneigh=avneigh*denom

      write(11,*) ' nneigh : avneigh :   <enerr> : max enerr :',
     .   ' <graderr> : max graderr '
      write(11,*) '        :         :    kJ/mol :  kJ/mol   :',
     .   '    %      :     %    ' 
      write(11,*)
      write(11,23) nneigh,avneigh,averr,enmx,avgrad,gdmx
23    format('qt ',i3,2x,f7.2,4x,4(1x,g12.5))

      write(11,*)
      return
      end
