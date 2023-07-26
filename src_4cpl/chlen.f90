      subroutine chlen (char,n,lench)
!
      character*1 char(n)
!
!  search for first blank in character string
!
      do 10 k=1,n
      if (char(k).eq.' ') then
        lench=k-1
        return
      endif
   10 continue
!
!  no blank found, string must be full
!
      lench= n
!
      return
      end
