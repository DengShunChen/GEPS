      subroutine cudlfs_n   &                              
     &    (nxj,      klon,     klev,  &                       
     &     kcbot,    kctop,    lndj,   ldcum,   &              
     &     ptenh,    pqenh,    puen,     pven,     &          
     &     pten,     pqsen,    pgeo,                &            
     &     pgeoh,    paph,     ptu,      pqu,      plu,&        
     &     puu,      pvu,      pmfub,    prfl,          &        
     &     ptd,      pqd,      pud,      pvd,            &        
     &     pmfd,     pmfds,    pmfdq,    pdmfdp,          &     
     &     kdtop,    lddraf)                                      
                                                                   
!          this routine calculates level of free sinking for                     
!          cumulus downdrafts and specifies t,q,u and v values                   
                                                                   
!          m.tiedtke         e.c.m.w.f.    12/86 modif. 12/89                    
                                                                     
!          purpose.                                                              
!          --------                                                              
!          to produce lfs-values for cumulus downdrafts                          
!          for massflux cumulus parameterization                                 
                                                                  
!          interface                                                             
!          ---------                                                             
!          this routine is called from *cumastr*.                                
!          input are environmental values of t,q,u,v,p,phi                       
!          and updraft values t,q,u and v and also                               
!          cloud base massflux and cu-precipitation rate.                        
!          it returns t,q,u and v values and massflux at lfs.                    
!          method.                                                               
!          check for negative buoyancy of air of equal parts of                  
!          moist environmental air and cloud air.                                
!     parameter     description                                   units          
!     ---------     -----------                                   -----          
!     input parameters (integer):                                                
!    *klon*         number of grid points per packet                             
!    *klev*         number of levels                                             
!    *kcbot*        cloud base level                                             
!    *kctop*        cloud top level                                              
!    input parameters (logical):                                                 
!    *lndj*       land sea mask (1 for land)                              
!    *ldcum*        flag: .true. for convective points                           
!    input parameters (real):                                                    
!    *ptenh*        env. temperature (t+1) on half levels          k             
!    *pqenh*        env. spec. humidity (t+1) on half levels     kg/kg           
!    *puen*         provisional environment u-velocity (t+1)      m/s            
!    *pven*         provisional environment v-velocity (t+1)      m/s            
!    *pten*         provisional environment temperature (t+1)       k            
!    *pqsen*        environment spec. saturation humidity (t+1)   kg/kg          
!    *pgeo*         geopotential                                  m2/s2          
!    *pgeoh*        geopotential on half levels                  m2/s2           
!    *paph*         provisional pressure on half levels           pa             
!    *ptu*          temperature in updrafts                        k             
!    *pqu*          spec. humidity in updrafts                   kg/kg           
!    *plu*          liquid water content in updrafts             kg/kg           
!    *puu*          u-velocity in updrafts                        m/s            
!    *pvu*          v-velocity in updrafts                        m/s            
!    *pmfub*        massflux in updrafts at cloud base           kg/(m2*s)       
!    updated parameters (real):                                                  
!    *prfl*         precipitation rate                           kg/(m2*s)       
!    output parameters (real):                                                   
!    *ptd*          temperature in downdrafts                      k             
!    *pqd*          spec. humidity in downdrafts                 kg/kg           
!    *pud*          u-velocity in downdrafts                      m/s            
!    *pvd*          v-velocity in downdrafts                      m/s            
!    *pmfd*         massflux in downdrafts                       kg/(m2*s)       
!    *pmfds*        flux of dry static energy in downdrafts       j/(m2*s)       
!    *pmfdq*        flux of spec. humidity in downdrafts         kg/(m2*s)       
!    *pdmfdp*       flux difference of precip. in downdrafts     kg/(m2*s)       
!    output parameters (integer):                                                
!    *kdtop*        top level of downdrafts                                      
!    output parameters (logical):                                                
!    *lddraf*       .true. if downdrafts exist                                   
!          externals                                                             
!          ---------                                                             
!          *cuadjtq* for calculating wet bulb t and q at lfs                     
!----------------------------------------------------------------------
  USE shr_kind_mod, only: r8 => shr_kind_r8
  USE mo_constants,    ONLY: vtmpc1,  &! vtmpc1=rv/rd-1
                             cpd       ! specific heat at constant pressure
  USE mo_cumulus_flux, ONLY: lmfdd,   &! true if cumulus downdr. is switched on
                             cmfdeps, &! fractional massflux for downdr. at lfs
                             lmfdudv   ! true if cumulus friction is switched on
         
!---------------------------------------------------------------------- 
      implicit none                                              
                   
      integer  klev,klon,nxj                                            
      real     ptenh(klon,klev),       pqenh(klon,klev), &     
     &         puen(klon,klev),        pven(klon,klev),  &            
     &         pten(klon,klev),        pqsen(klon,klev),  &            
     &         pgeo(klon,klev),                       &               
     &         pgeoh(klon,klev+1),     paph(klon,klev+1),&             
     &         ptu(klon,klev),         pqu(klon,klev),   &             
     &         puu(klon,klev),         pvu(klon,klev),   &           
     &         plu(klon,klev),                          &            
     &         pmfub(klon),            prfl(klon)                      
                                                                      
      real     ptd(klon,klev),         pqd(klon,klev),   &            
     &         pud(klon,klev),         pvd(klon,klev),    &     
     &         pmfd(klon,klev),        pmfds(klon,klev),  &        
     &         pmfdq(klon,klev),       pdmfdp(klon,klev)             
      integer  kcbot(klon),            kctop(klon),       &         
     &         kdtop(klon),            ikhsmin(klon)                
      logical  ldcum(klon),                              &
     &         lddraf(klon)   
      integer  lndj(klon)                                            
                                                                      
      real     ztenwb(klon,klev),      zqenwb(klon,klev), &          
     &         zcond(klon),            zph(klon),         &           
     &         zhsmin(klon)                                           
      logical  llo2(klon)                                             
! local variables
      integer  jl,jk
      integer  is,ik,icall,ike
      real     zhsk,zttest,zqtest,zbuo,zmftop,foelhm
!----------------------------------------------------------------------          
!     1.           set default values for downdrafts                             
!                  ---------------------------------                            
!xb110>
      zcond = 0.
!xb110< 
       do jl = 1, nxj         
        lddraf(jl)=.false.                                              
        kdtop(jl)=klev+1      
        ikhsmin(jl)=klev+1
        zhsmin(jl)=1.e8                                                 
      enddo                                                            
!----------------------------------------------------------------------          
                                                                        
!     2.           determine level of free sinking:                              
!                  downdrafts shall start at model level of minimum              
!                  of saturation moist static energy or below                    
!                  respectively                                                  
                                                                       
!                  for every point and proceed as follows:                       
                                                                       
!                    (1) determine level of minimum of hs                        
!                    (2) determine wet bulb environmental t and q                
!                    (3) do mixing with cumulus cloud air                        
!                    (4) check for negative buoyancy                             
!                    (5) if buoyancy>0 repeat (2) to (4) for next                
!                        level below                                             
                                                                       
!                  the assumption is that air of downdrafts is mixture           
!                  of 50% cloud air + 50% environmental air at wet bulb          
!                  temperature (i.e. which became saturated due to               
!                  evaporation of rain and cloud water)                          
!                  ----------------------------------------------------          
      do jk=3,klev-2
       do jl = 1, nxj         
           zhsk=cpd*pten(jl,jk)+pgeo(jl,jk) +  &
     &         foelhm(pten(jl,jk))*pqsen(jl,jk)
           if(zhsk .lt. zhsmin(jl)) then
              zhsmin(jl) = zhsk
              ikhsmin(jl)= jk
           end if
         end do
      end do


      ike=klev-3                                                       
      do jk=3,ike                                                      
                                                                       
!     2.1          calculate wet-bulb temperature and moisture                   
!                  for environmental air in *cuadjtq*                            
!                  -------------------------------------------                   
        is=0                                                          
       do jl = 1, nxj         
          ztenwb(jl,jk)=ptenh(jl,jk)                                  
          zqenwb(jl,jk)=pqenh(jl,jk)                                  
          zph(jl)=paph(jl,jk)                                           
          llo2(jl)=ldcum(jl).and.prfl(jl).gt.0..and..not.lddraf(jl).and. &      
     &     (jk.lt.kcbot(jl).and.jk.gt.kctop(jl)).and. jk.ge.ikhsmin(jl) 
          if(llo2(jl))then                                              
           is=is+1                                                     
          endif                                                         
        enddo                                                           
        if(is.eq.0) cycle                                               
                                                                       
        ik=jk                                                           
        icall=2                                                         
        call cuadjtq_n                                           &       
     &   ( nxj, klon, klev, ik, zph, ztenwb, zqenwb, llo2, icall)            
                                                                       
!     2.2          do mixing of cumulus and environmental air                    
!                  and check for negative buoyancy.                              
!                  then set values for downdraft at lfs.                         
!                  ----------------------------------------                      
       do jl = 1, nxj         
          if(llo2(jl)) then                                             
            zttest=0.5*(ptu(jl,jk)+ztenwb(jl,jk))                      
            zqtest=0.5*(pqu(jl,jk)+zqenwb(jl,jk))                      
            zbuo=zttest*(1.+vtmpc1  *zqtest)-                    &     
     &       ptenh(jl,jk)*(1.+vtmpc1  *pqenh(jl,jk))                 
            zcond(jl)=pqenh(jl,jk)-zqenwb(jl,jk)                      
            zmftop=-cmfdeps*pmfub(jl)                               
            if(zbuo.lt.0..and.prfl(jl).gt.10.*zmftop*zcond(jl)) then  
              kdtop(jl)=jk                                         
              lddraf(jl)=.true.                                       
              ptd(jl,jk)=zttest                                      
              pqd(jl,jk)=zqtest                                     
              pmfd(jl,jk)=zmftop                                    
              pmfds(jl,jk)=pmfd(jl,jk)*(cpd*ptd(jl,jk)+pgeoh(jl,jk))   
              pmfdq(jl,jk)=pmfd(jl,jk)*pqd(jl,jk)                      
              pdmfdp(jl,jk-1)=-0.5*pmfd(jl,jk)*zcond(jl)             
              prfl(jl)=prfl(jl)+pdmfdp(jl,jk-1)                      
            endif                                                     
          endif                                                     
        enddo                                                        
                                                                     
      enddo                                                          
                                                                     
      return             
      end subroutine cudlfs_n 
