c --- external product 
      subroutine cross_x(a,b,c)
      implicit real*8 (a-h,p-z), integer(o)
      dimension a(3),b(3),c(3)
      c(1)=a(2)*b(3)-a(3)*b(2)
      c(2)=a(3)*b(1)-a(1)*b(3)
      c(3)=a(1)*b(2)-a(2)*b(1)
      return
      end
c
      subroutine dinv33x(plat,qlat)
C- This is a replacement of dinv33 of Ferdi's GW  => dinv33(plat,1,qlat,det) --------------
Cr THIS IS the SAME as the one of dinv33 in extens.f in ferdi/lmto/extens.f
      implicit none
      double precision plat(3,3),qlat(3,3),det
      call cross_x(plat(1,2),plat(1,3), qlat     )
      call cross_x(plat(1,3),plat     , qlat(1,2))
      call cross_x(plat     ,plat(1,2), qlat(1,3))
      det  = sum( plat(1:3,1)*qlat(1:3,1) )
      qlat = qlat/det
      end

      subroutine dinv33y(plat,qlat,det)
C- This is a replacement of dinv33 of Ferdi's GW  => dinv33(plat,1,qlat,det) --------------
Cr THIS IS the SAME as the one of dinv33 in extens.f in ferdi/lmto/extens.f
      implicit none
      double precision plat(3,3),qlat(3,3),det
      call cross_x(plat(1,2),plat(1,3), qlat     )
      call cross_x(plat(1,3),plat     , qlat(1,2))
      call cross_x(plat     ,plat(1,2), qlat(1,3))
      det  = sum( plat(1:3,1)*qlat(1:3,1) )
      qlat = qlat/det
      end
