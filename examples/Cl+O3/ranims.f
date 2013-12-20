      real*8 function usran(ir)
c
c   this subroutine generates random values between 0.0 and 1.0 using
c   an integer seed
c   it is based on the imsl routine ggubs.
c
c   real*8 version
c
      implicit real*8 (a-h,o-z)
      parameter(da=16807.d0,db=2147483647.d0,dc=2147483648.d0)
      ir=abs(mod(da*ir,db)+0.5d0)
      usran=dfloat(ir)/dc
      return
      end
