c
c 
      subroutine efld_cd(dz,ls,vnorm,flux1,flux2,elparnw,flux0) !for cqlpmod.eq."enabled"
      implicit integer (i-n), real*8 (a-h,o-z)

c.......................................................................
c     Calculates the electric field (V/cm) and flux0 (1/(cm**2*sec)), 
c     to maintain constant flux (current) 
c     with no external electric field, such as in a RF driven torus.
c     (The internal electric field is generated by charge separation.)
c     Uses fluxes from the h and g functions.
c.......................................................................

      include 'param.h'

c..... input: ls,dz,flux0,flux1,flux2

      dimension dz(1:ls),flux1(0:ls+1),flux2(0:ls+1)
         !YuP[2019-05-30] corrected sub.efld_cd
         !                so that dz argument starts with index 1.
    
c.....input/output: 
      dimension elparnw(0:ls+1)

c.......................................................................
c     
c.......................................................................
      

      sum1=0.d0
      sum2=0.d0
      do kk=1,ls
         sum1=sum1 + flux1(kk)*dz(kk)/flux2(kk)
         sum2=sum2 + dz(kk)/flux2(kk)
      enddo

      flux0=vnorm*sum1/sum2

      do kk=1,ls
         delparnw=(flux0-flux1(kk)*vnorm)/(flux2(kk)*vnorm)
         elparnw(kk)=elparnw(kk) +300.d0*delparnw
      enddo

      return
      end  
c
c
      real*8 function fluxpar(kopt,x,coss,cynt2,cint2,f,iy_l,jx)
      !YuP Note a usual call: fluxpar(1,x,coss(1:iy_(l_),l_),
      !                 cynt2(1:iy_(l_),l_),cint2,temp1,iy_(l_),jx)
      !Renamed iy to iy_l, for convenience (input/local)
      implicit integer (i-n), real*8 (a-h,o-z)

c..................................................................
c     Compute the parallel flux ---
c     parallel flux (particles/sec per cm^2)=fluxpar*vnorm,
c     where vnorm is the velocity normalization in cm/sec, 
c     i.e., fluxpar has the units of 1/cm^3.
c.................................................................
c     kopt=1 the flux from the total function f
c     kopt=2 the flux from half of f  at v_par >0
c     kopt=3 the flux from half of f  at v_par <0
c..................................................................
      integer iy_l
      real*8 x(jx),coss(iy_l),cynt2(iy_l),cint2(jx),
     +                 f(0:iy_l+1,0:jx+1)

      fluxpar=0.d0

      if (kopt.eq.1) then
        do 20 i=1,iy_l
          dum=cynt2(i)*coss(i)         
          do 10 j=1,jx          
            fluxpar=fluxpar+f(i,j)*cint2(j)*x(j)*dum         
 10        continue         
 20      continue
      endif

      if (kopt.eq.2) then
        ihy=iy_l/2
        do i=1,ihy
          dum=cynt2(i)*coss(i)         
          do j=1,jx          
            fluxpar=fluxpar+f(i,j)*cint2(j)*x(j)*dum         
           enddo         
         enddo
      endif

      if (kopt.eq.3) then
        ihy=iy_l/2
        do i=ihy+1,iy_l
          dum=cynt2(i)*coss(i)         
          do j=1,jx          
            fluxpar=fluxpar+f(i,j)*cint2(j)*x(j)*dum         
           enddo         
         enddo
      endif

      return
      end
