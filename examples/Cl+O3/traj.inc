
      parameter(natomm=4)
      parameter(nbondm=natomm*(natomm-1)/2)
      parameter(nmax=25000)

      parameter (ndatam=2100)
      parameter (ngroupm=120)
      parameter (maxfo=49600,maxf=48000)
      parameter (na3=3*natomm,nintm=na3-6)

c      implicit real*8(a-h,o-z)
      implicit double precision (a-h,o-z)

      character*20 msunit,enunit,lnunit,tmunit

      character*80 title,comlin

      common /txt/ title,comlin

      common /sclfaa/ msunit,enunit,lnunit,tmunit

      common /sclfac/ amsc,esc,rsc,tsc,tups,scalef,
     .  eukjpm,ccmps,ccmptu,hjs,hkjpms,heutu

      common /ctl/ delta,pi,todeg,toang
      common/sig/ sigma,sigma2

   
      common /nvar/ nloop,nstep,nprint,nmol,ir,np,nq,nchannel(nmax),nstop

      common/shepdata/ v0(ndatam),v1(ndatam,nintm),v2(ndatam,nintm),
     .  ut(ndatam,nintm,nbondm),rda(ndatam,nbondm),zda(ndatam,nintm)
     .  ,rads(ndatam),wtol,outer,totsum,avrads,potmax,sumb(ndatam,nbondm)
     .  ,wt(maxf)

      common/shepuhf/ v0u(ndatam),v1u(ndatam,nintm),v2u(ndatam,nintm),
     .  utu(ndatam,nintm,nbondm),zdau(ndatam,nintm)

      common /cart/ amas(natomm),c(nmax,natomm,3),cold(nmax,natomm,3),
     .  v(nmax,natomm,3),vold(nmax,natomm,3),vm2(nmax,natomm,3),
     .  mb(nbondm),nb(nbondm) 

      common /enmol/ en(nmax),enprev(nmax)
      common /enmolu/ enu(nmax),enprevu(nmax)

      common /dvi/ r(nmax,nbondm),dr(nmax,nbondm,natomm,3)
     .  ,dv(nmax,natomm,3),dvr(nmax,nbondm)
     
      common /dviu/ dvru(nmax,nbondm)

      common /oup/ echeck(nmax),ek1(nmax),ek2(nmax),econ,tvg(nbondm),
     .  errck(nmax),nfin(nmax),ntest(nmax),nwrite

      common/ints/ ndata,nforc,mda(maxf),nda(maxf),lowp,nneigh,ipart,
     . ipow,neighci,neighco,maxnb,nlowstop,nforco,mdao(maxfo),ndao(maxfo),
     .  neighc

      common/prout/sample,rms,grms,totsm

      common/perm/ip(nbondm,ngroupm),nperm(ngroupm,natomm)

      common/addpot/ vmin,nsel,nadd

      common /global/ preact

      common /fragin/ nfraga,nfragb,nbfraga,nbfragb,
     .  ifraga(natomm),ifragb(natomm),ibfraga(nbondm),ibfragb(nbondm)

      common/frag/ ifrag(natomm),kf(natomm),map(natomm,natomm),
     .               ibfrag(natomm,nbondm),kb(natomm),nfrag,
     .             mfrag(nmax),ichansave(nmax,natomm)

      common/fragment/ rmin(natomm,natomm), bipmax, distab


        common/actual/ natom, nbond, n3, nint, ngroup, nset

        character*2 lab(natomm)
        common/elabels/ lab

