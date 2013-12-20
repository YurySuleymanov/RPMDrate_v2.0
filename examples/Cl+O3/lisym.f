      subroutine lisym(iwr,csym)
c
c   18/02/89: edited nsym=80 -> 78
c
      character*1 csym
      nsym=78
      write(iwr,50)(csym(1:1),i=1,nsym)
   50 format(1x,78a1)
      return
      end
