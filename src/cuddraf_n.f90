       subroutine cuddraf_n                                 &      
     &   ( nxj, klon,     klev,    lddraf                      &
     &   , ptenh,    pqenh,    puen,     pven             &             
     &   , pgeo,     pgeoh,    paph,     prfl             &             
     &   , ptd,      pqd,      pud,      pvd,      pmfu   &           
     &   , pmfd,     pmfds,    pmfdq,    pdmfdp,   pmfdde_rate &
     &   )    
                                                               
!          this routine calculates cumulus downdraft descent              
                                                                 
!          m.tiedtke         e.c.m.w.f.    12/86 modif. 12/89             
                                                                   
!          purpose.                                                       
!          --------                                                       
!          to produce the vertical profiles for cumulus downdrafts        
!          (i.e. t,q,u and v and fluxes)                                  
                                                                    
!          interface                                                      
!          ---------                                                      
                                                                    
!          this routine is called from *cumastr*.                         
!          input is t,q,p,phi,u,v at half levels.                         
!          it returns fluxes of s,q and evaporation rate                  
!          and u,v at levels where downdraft occurs                       
                                                                     
!          method.                                                       
!          --------                                                       
!          calculate moist descent for entraining/detraining plume by     
!          a) moving air dry-adiabatically to next level below and        
!          b) correcting for evaporation to obtain saturated state.       
                                                                 
!     parameter     description                                   units   
!     ---------     -----------                                   -----   
!     input parameters (integer):                                         
                                                                  
!    *klon*         number of grid points per packet                      
!    *klev*         number of levels                                      
                                                                    
!    input parameters (logical):                                          
                                                                  
!    *lddraf*       .true. if downdrafts exist                            
                                                                 
!    input parameters (real):                                             
                                                                 
!    *ptenh*        env. temperature (t+1) on half levels          k      
!    *pqenh*        env. spec. humidity (t+1) on half levels     kg/kg    
!    *puen*         provisional environment u-velocity (t+1)      m/s     
!    *pven*         provisional environment v-velocity (t+1)      m/s     
!    *pgeo*         geopotential                                  m2/s2   
!    *pgeoh*        geopotential on half levels                  m2/s2    
!    *paph*         provisional pressure on half levels           pa      
!    *pmfu*         massflux updrafts                           kg/(m2*s) 
                                                                 
!    updated parameters (real):                                           
                                                           
!    *prfl*         precipitation rate                           kg/(m2*s)
                                                                
!    output parameters (real):                                            
                                                       
!    *ptd*          temperature in downdrafts                      k      
!    *pqd*          spec. humidity in downdrafts                 kg/kg    
!    *pud*          u-velocity in downdrafts                      m/s     
!    *pvd*          v-velocity in downdrafts                      m/s     
!    *pmfd*         massflux in downdrafts                       kg/(m2*s)!    *pmfds*        flux of dry static energy in downdrafts       j/(m2*s)
!    *pmfdq*        flux of spec. humidity in downdrafts         kg/(m2*s)
!    *pdmfdp*       flux difference of precip. in downdrafts     kg/(m2*s)
                                                               
!          externals                                                      
!          ---------                                                      
!          *cuadjtq* for adjusting t and q due to evaporation in          
!          saturated descent                                              
!----------------------------------------------------------------------
  USE mo_constants,    ONLY: rd,      &! gas constant for dry air
                             cpd,     &! specific heat at constant pressure
                             rcpd,    &! rcpd=1./cpd
                             vtmpc1,  &! vtmpc1=rv/rd-1
                             g,       &! gravity acceleration
                             zrg        ! 1.0/g
  USE mo_cumulus_flux, ONLY: entrdd,  &! entrainment rate for cumulus downdrafts
                             cmfcmin, &! minimum massflux value (for safety)
                             lmfdudv   ! true if cumulus friction is switched on
  
!---------------------------------------------------------------------- 
      implicit none
      
      integer  klev,klon,nxj
      real     ptenh(klon,klev),       pqenh(klon,klev),   &   
     &         puen(klon,klev),        pven(klon,klev),    & 
     &         pgeoh(klon,klev+1),     paph(klon,klev+1),  &  
     &         pgeo(klon,klev),        pmfu(klon,klev)             
                                                                   
      real     ptd(klon,klev),         pqd(klon,klev),     &    
     &         pud(klon,klev),         pvd(klon,klev),     &     
     &         pmfd(klon,klev),        pmfds(klon,klev),   &    
     &         pmfdq(klon,klev),       pdmfdp(klon,klev),  &    
     &         prfl(klon)     
      real     pmfdde_rate(klon,klev)   
      logical  lddraf(klon)   
                                                         
      real     zdmfen(klon),           zdmfde(klon),       &   
     &         zcond(klon),            zoentr(klon),       & 
     &         zbuoy(klon)                              
      real     zph(klon)                         
      logical  llo2(klon)                                   
      logical  llo1
! local variables
      integer  jl,jk
      integer  is,ik,icall,ike, itopde(klon)
      real     zentr,zdz,zzentr,zseen,zqeen,zsdde,zqdde,zdmfdp
      real     zmfdsk,zmfdqk,zbuo,zrain,zbuoyz,zmfduk,zmfdvk
                                                                
!----------------------------------------------------------------------   
!     1.           calculate moist descent for cumulus downdraft by       
!                     (a) calculating entrainment/detrainment rates,      
!                         including organized entrainment dependent on    
!                         negative buoyancy and assuming                  
!                         linear decrease of massflux in pbl              
!                     (b) doing moist descent - evaporative cooling       
!                         and moistening is calculated in *cuadjtq*       
!                     (c) checking for negative buoyancy and              
!                         specifying final t,q,u,v and downward fluxes    
!                    -------------------------------------------------    
       do jl = 1, nxj         
        zoentr(jl)=0.
        zbuoy(jl)=0.
        zdmfen(jl)=0.
        zdmfde(jl)=0.
!xb110>
        itopde(jl)=klev
        zcond(jl) = 0.
!xb110<
      enddo

      do jk=klev,1,-1              
       do jl = 1, nxj         
         pmfdde_rate(jl,jk) = 0.
         if((paph(jl,klev+1)-paph(jl,jk)).lt. 60.e2) itopde(jl)=jk
       end do
      end do
     
      do jk=3,klev
      is=0
       do jl = 1, nxj         
      zph(jl)=paph(jl,jk)
      llo2(jl)=lddraf(jl).and.pmfd(jl,jk-1).lt.0.
      if(llo2(jl)) then
         is=is+1
      endif
      end do

      if(is.eq.0) cycle
       do jl = 1, nxj         
      if(llo2(jl)) then
         zentr = entrdd*pmfd(jl,jk-1)*(pgeoh(jl,jk-1)-pgeoh(jl,jk))*zrg
         zdmfen(jl)=zentr
         zdmfde(jl)=zentr
      end if
      end do

       do jl = 1, nxj         
         if(llo2(jl)) then
         if(jk.gt.itopde(jl)) then
            zdmfen(jl)=0.
            zdmfde(jl)=pmfd(jl,itopde(jl))*      &
            (paph(jl,jk)-paph(jl,jk-1))/     &
            (paph(jl,klev+1)-paph(jl,itopde(jl)))
         end if
         end if
        end do

         do jl = 1,nxj
          if(llo2(jl)) then
          if(jk.le.itopde(jl)) then
            zdz=-(pgeoh(jl,jk-1)-pgeoh(jl,jk))*zrg
            zzentr=zoentr(jl)*zdz*pmfd(jl,jk-1)
            zdmfen(jl)=zdmfen(jl)+zzentr
            zdmfen(jl)=max(zdmfen(jl),0.3*pmfd(jl,jk-1))
            zdmfen(jl)=max(zdmfen(jl),-0.75*pmfu(jl,jk)-   &
   &         (pmfd(jl,jk-1)-zdmfde(jl)))
            zdmfen(jl)=min(zdmfen(jl),0.)
          endif
         endif
        enddo

       do jl = 1, nxj         
         if(llo2(jl)) then
            pmfd(jl,jk)=pmfd(jl,jk-1)+zdmfen(jl)-zdmfde(jl)
            zseen=(cpd*ptenh(jl,jk-1)+pgeoh(jl,jk-1))*zdmfen(jl)
            zqeen=pqenh(jl,jk-1)*zdmfen(jl)
            zsdde=(cpd*ptd(jl,jk-1)+pgeoh(jl,jk-1))*zdmfde(jl)
            zqdde=pqd(jl,jk-1)*zdmfde(jl)
            zmfdsk=pmfds(jl,jk-1)+zseen-zsdde
            zmfdqk=pmfdq(jl,jk-1)+zqeen-zqdde
            pqd(jl,jk)=zmfdqk*(1./min(-cmfcmin,pmfd(jl,jk)))
            ptd(jl,jk)=(zmfdsk*(1./min(-cmfcmin,pmfd(jl,jk)))- &
                       pgeoh(jl,jk))*rcpd
            ptd(jl,jk)=min(400.,ptd(jl,jk))
            ptd(jl,jk)=max(100.,ptd(jl,jk))
            zcond(jl)=pqd(jl,jk)
         end if
      end do

      ik=jk
      icall=2
      call cuadjtq_n(nxj, klon,klev,ik,zph,ptd,pqd,llo2,icall)
       do jl = 1, nxj         
         if(llo2(jl)) then
            zcond(jl)=zcond(jl)-pqd(jl,jk)
            zbuo=ptd(jl,jk)*(1.+vtmpc1*pqd(jl,jk))- &
                 ptenh(jl,jk)*(1.+vtmpc1*pqenh(jl,jk))
            if(prfl(jl).gt.0..and.pmfu(jl,jk).gt.0.) then
              zrain=prfl(jl)/pmfu(jl,jk)
              zbuo=zbuo-ptd(jl,jk)*zrain
            endif
            if(zbuo.ge.0..or.prfl(jl).le.(pmfd(jl,jk)*zcond(jl))) then
               pmfd(jl,jk)=0.
               zbuo=0.
            end if
            pmfds(jl,jk)=(cpd*ptd(jl,jk)+pgeoh(jl,jk))*pmfd(jl,jk)
            pmfdq(jl,jk)=pqd(jl,jk)*pmfd(jl,jk)
            zdmfdp=-pmfd(jl,jk)*zcond(jl)
            pdmfdp(jl,jk-1)=zdmfdp
            prfl(jl)=prfl(jl)+zdmfdp
! compute organized entrainment for use at next level
            zbuoyz=zbuo/ptenh(jl,jk)
            zbuoyz=min(zbuoyz,0.0)
            zdz=-(pgeo(jl,jk-1)-pgeo(jl,jk))
            zbuoy(jl)=zbuoy(jl)+zbuoyz*zdz
            zoentr(jl)=g*zbuoyz*0.5/(1.+zbuoy(jl))
            pmfdde_rate(jl,jk) = -zdmfde(jl)
          end if                                              
        end do
                                                      
      end do                                          
      return            
      end subroutine cuddraf_n
