       subroutine cuflx_n         & 
     &  (  nxj,     klon,     klev,     ztmst &  
     &  ,  pten,     pqen,     pqsen,    ptenh,    pqenh   &    
     &  ,  paph,     pap,      pgeoh,    lndj,   ldcum   &   
     &  ,  kcbot,    kctop,    kdtop,    ktopm2            &     
     &  ,  ktype,    lddraf                                &   
     &  ,  pmfu,     pmfd,     pmfus,    pmfds             &    
     &  ,  pmfuq,    pmfdq,    pmful,    plude             &   
     &  ,  pdmfup,   pdmfdp,   pdpmel,   plglac            &    
     &  ,  prain,    pmfdde_rate, pmflxr, pmflxs )
 
!          m.tiedtke         e.c.m.w.f.     7/86 modif. 12/89                  
!          purpose                                                             
!          -------                                                             
!          this routine does the final calculation of convective               
!          fluxes in the cloud layer and in the subcloud layer                 
!          interface                                                           
!          ---------                                                           
!          this routine is called from *cumastr*.                              
!     parameter     description                                   units        
!     ---------     -----------                                   -----        
!     input parameters (integer):                                              
!    *klon*         number of grid points per packet                           
!    *klev*         number of levels                                           
!    *kcbot*        cloud base level                                           
!    *kctop*        cloud top level                                            
!    *kdtop*        top level of downdrafts                                    
!    input parameters (logical):                                               
!    *lndj*       land sea mask (1 for land)                            
!    *ldcum*        flag: .true. for convective points                         
!    input parameters (real):                                                  
!    *ptsphy*       time step for the physics                       s          
!    *pten*         provisional environment temperature (t+1)       k          
!    *pqen*         provisional environment spec. humidity (t+1)  kg/kg        
!    *pqsen*        environment spec. saturation humidity (t+1)   kg/kg        
!    *ptenh*        env. temperature (t+1) on half levels           k          
!    *pqenh*        env. spec. humidity (t+1) on half levels      kg/kg        
!    *paph*         provisional pressure on half levels            pa          
!    *pap*          provisional pressure on full levels            pa          
!    *pgeoh*        geopotential on half levels                   m2/s2        
!    updated parameters (integer):                                             
!    *ktype*        set to zero if ldcum=.false.                               
!    updated parameters (logical):                                             
!    *lddraf*       set to .false. if ldcum=.false. or kdtop<kctop             
!    updated parameters (real):                                                
!    *pmfu*         massflux in updrafts                          kg/(m2*s)    
!    *pmfd*         massflux in downdrafts                        kg/(m2*s)    
!    *pmfus*        flux of dry static energy in updrafts          j/(m2*s)    
!    *pmfds*        flux of dry static energy in downdrafts        j/(m2*s)    
!    *pmfuq*        flux of spec. humidity in updrafts            kg/(m2*s)    
!    *pmfdq*        flux of spec. humidity in downdrafts          kg/(m2*s)    
!    *pmful*        flux of liquid water in updrafts              kg/(m2*s)    
!    *plude*        detrained liquid water                        kg/(m3*s)    
!    *pdmfup*       flux difference of precip. in updrafts        kg/(m2*s)    
!    *pdmfdp*       flux difference of precip. in downdrafts      kg/(m2*s)    
!    output parameters (real):                                                 
!    *pdpmel*       change in precip.-fluxes due to melting       kg/(m2*s)    
!    *plglac*       flux of frozen cloud water in updrafts        kg/(m2*s)    
!    *pmflxr*       convective rain flux                          kg/(m2*s)    
!    *pmflxs*       convective snow flux                          kg/(m2*s)    
!    *prain*        total precip. produced in conv. updrafts      kg/(m2*s)    
!                   (no evaporation in downdrafts)                             
!          externals                                                           
!          ---------                                                           
!          none                                                                
!----------------------------------------------------------------------
  USE shr_kind_mod, only: r8 => shr_kind_r8
  USE mo_constants,     ONLY: cpd,   &! specific heat at constant pressure
                              alf,   &! latent heat for vaporisation
                              g,     &! gravity acceleration
                              tmelt   ! temperature of fusion of ice
        
!----------------------------------------------------------------------
      implicit none                                                     

      integer  klev,klon,ktopm2,nxj 
      real     foelhm,foeewm,foealfa
      real     pten(klon,klev),        ptenh(klon,klev),          &
     &         pqen(klon,klev),        pqsen(klon,klev),          & 
     &         pqenh(klon,klev),       pap(klon,klev),            & 
     &         paph(klon,klev+1),      pgeoh(klon,klev+1),        & 
     &         plglac(klon,klev)                                  
                                                                    
      real     pmfu(klon,klev),        pmfd(klon,klev),           & 
     &         pmfus(klon,klev),       pmfds(klon,klev),          & 
     &         pmfuq(klon,klev),       pmfdq(klon,klev),          &  
     &         pdmfup(klon,klev),      pdmfdp(klon,klev),         & 
     &         pdpmel(klon,klev),      prain(klon),               & 
     &         pmful(klon,klev),       plude(klon,klev),          &
     &         pmflxr(klon,klev+1),    pmflxs(klon,klev+1)  
      real     pmfdde_rate(klon,klev)               
      integer  kcbot(klon),            kctop(klon),               &  
     &         kdtop(klon),            ktype(klon)        
      logical  lddraf(klon),                                &
     &         ldcum(klon)  
      integer  lndj(klon)
! local variables
      integer  jl,jk
      integer  is,ik,icall,ike,ikb
      real     ztmst,ztaumel,zcons1a,zcons1,zcons2,zcucov,zcpecons
      real     zalfaw,zrfl,zdrfl1,zrnew,zrmin,zrfln,zdrfl,zdenom
      real     zpdr,zpds,zzp,zfac,zsnmlt
      real     rhevap(klon)
      integer  idbas(klon)
      logical  llddraf
!--------------------------------------------------------------------          
!*             specify constants  

      ztaumel=18000.
      zcons1a=cpd/(alf*g*ztaumel)
      zcons2=3./(g*ztmst)
      zcucov=0.05
      zcpecons=5.44e-4/g

!*    1.0          determine final convective fluxes                           
!                  ---------------------------------                           
       do jl = 1, nxj         
        prain(jl)=0.
        rhevap(jl)=0.          !xb110
        if(.not.ldcum(jl).or.kdtop(jl).lt.kctop(jl)) lddraf(jl)=.false.
        if(.not.ldcum(jl)) ktype(jl)=0
        idbas(jl) = klev
        if(lndj(jl) .eq. 1) then
          rhevap(jl) = 0.7
        else
          rhevap(jl) = 0.9
        end if
      enddo

!xb110>
      pmflxr = 0.
      pmflxs = 0.
      pdpmel = 0.
!xb110<

      ktopm2= 2
      do jk=ktopm2,klev
        ikb = min(jk+1,klev)
       do jl = 1, nxj        
!xb110>
!          pmflxr(jl,jk) = 0.
!          pmflxs(jl,jk) = 0.
!          pdpmel(jl,jk) = 0.
!xb110<
          if(ldcum(jl).and.jk.ge.kctop(jl)) then
            pmfus(jl,jk)=pmfus(jl,jk)-pmfu(jl,jk)*           &
     &       (cpd*ptenh(jl,jk)+pgeoh(jl,jk))
            pmfuq(jl,jk)=pmfuq(jl,jk)-pmfu(jl,jk)*pqenh(jl,jk)
            plglac(jl,jk)=pmfu(jl,jk)*plglac(jl,jk)
            llddraf = lddraf(jl) .and. jk >= kdtop(jl)
            if ( llddraf .and.jk.ge.kdtop(jl)) then
              pmfds(jl,jk) = pmfds(jl,jk)-pmfd(jl,jk) * &
                           (cpd*ptenh(jl,jk)+pgeoh(jl,jk))
              pmfdq(jl,jk) = pmfdq(jl,jk)-pmfd(jl,jk)*pqenh(jl,jk)
            else
              pmfd(jl,jk) = 0.
              pmfds(jl,jk) = 0.
              pmfdq(jl,jk) = 0.
              pdmfdp(jl,jk-1) = 0.
            end if
            if ( llddraf .and. pmfd(jl,jk) < 0. .and. &
               abs(pmfd(jl,ikb)) < 1.e-20 ) then
               idbas(jl) = jk
            end if
          else
            pmfu(jl,jk)=0.
            pmfd(jl,jk)=0.
            pmfus(jl,jk)=0.
            pmfds(jl,jk)=0.
            pmfuq(jl,jk)=0.
            pmfdq(jl,jk)=0.
            pmful(jl,jk)=0.
            plglac(jl,jk)=0.
            pdmfup(jl,jk-1)=0.
            pdmfdp(jl,jk-1)=0.
            plude(jl,jk-1)=0.
          endif
        enddo
      enddo

       do jl = 1, nxj         
        pmflxr(jl,klev+1) = 0.
        pmflxs(jl,klev+1) = 0.
      end do
       do jl = 1, nxj         
        if(ldcum(jl)) then
          ikb=kcbot(jl)
          ik=ikb+1
          zzp=((paph(jl,klev+1)-paph(jl,ik))/                  &
     &     (paph(jl,klev+1)-paph(jl,ikb)))
          if(ktype(jl).eq.3) then
            zzp=zzp**2
          endif
          pmfu(jl,ik)=pmfu(jl,ikb)*zzp
          pmfus(jl,ik)=(pmfus(jl,ikb)-                         &
     &         foelhm(ptenh(jl,ikb))*pmful(jl,ikb))*zzp
          pmfuq(jl,ik)=(pmfuq(jl,ikb)+pmful(jl,ikb))*zzp
          pmful(jl,ik)=0.
        endif
      enddo

      do jk=ktopm2,klev
       do jl = 1, nxj         
          if(ldcum(jl).and.jk.gt.kcbot(jl)+1) then
            ikb=kcbot(jl)+1
            zzp=((paph(jl,klev+1)-paph(jl,jk))/                &
     &       (paph(jl,klev+1)-paph(jl,ikb)))
            if(ktype(jl).eq.3) then
              zzp=zzp**2
            endif
            pmfu(jl,jk)=pmfu(jl,ikb)*zzp
            pmfus(jl,jk)=pmfus(jl,ikb)*zzp
            pmfuq(jl,jk)=pmfuq(jl,ikb)*zzp
            pmful(jl,jk)=0.
          endif
          ik = idbas(jl)
          llddraf = lddraf(jl) .and. jk > ik .and. ik < klev
          if ( llddraf .and. ik == kcbot(jl)+1 ) then
           zzp = ((paph(jl,klev+1)-paph(jl,jk))/(paph(jl,klev+1)-paph(jl,ik)))
           if ( ktype(jl) == 3 ) zzp = zzp*zzp
            pmfd(jl,jk) = pmfd(jl,ik)*zzp
            pmfds(jl,jk) = pmfds(jl,ik)*zzp
            pmfdq(jl,jk) = pmfdq(jl,ik)*zzp
            pmfdde_rate(jl,jk) = -(pmfd(jl,jk-1)-pmfd(jl,jk))
           else if ( llddraf .and. ik /= kcbot(jl)+1 .and. jk == ik+1 ) then
            pmfdde_rate(jl,jk) = -(pmfd(jl,jk-1)-pmfd(jl,jk))
           end if
        enddo
      enddo
!*    2.            calculate rain/snow fall rates                             
!*                  calculate melting of snow                                  
!*                  calculate evaporation of precip                            
!                   -------------------------------                            

        do jk=ktopm2,klev
       do jl = 1, nxj         
          if(ldcum(jl) .and. jk >=kctop(jl)-1 ) then
            prain(jl)=prain(jl)+pdmfup(jl,jk)
            if(pmflxs(jl,jk).gt.0..and.pten(jl,jk).gt.tmelt) then
              zcons1=zcons1a*(1.+0.5*(pten(jl,jk)-tmelt))
              zfac=zcons1*(paph(jl,jk+1)-paph(jl,jk))
              zsnmlt=min(pmflxs(jl,jk),zfac*(pten(jl,jk)-tmelt))
              pdpmel(jl,jk)=zsnmlt
              pqsen(jl,jk)=foeewm(pten(jl,jk)-zsnmlt/zfac)/pap(jl,jk)
            endif
            zalfaw=foealfa(pten(jl,jk))
          !
          ! No liquid precipitation above melting level
          !
            if ( pten(jl,jk) < tmelt .and. zalfaw > 0. ) then
              plglac(jl,jk) = plglac(jl,jk)+zalfaw*(pdmfup(jl,jk)+pdmfdp(jl,jk))
              zalfaw = 0.
            end if
            pmflxr(jl,jk+1)=pmflxr(jl,jk)+zalfaw*               &
     &       (pdmfup(jl,jk)+pdmfdp(jl,jk))+pdpmel(jl,jk)
            pmflxs(jl,jk+1)=pmflxs(jl,jk)+(1.-zalfaw)*          &
     &       (pdmfup(jl,jk)+pdmfdp(jl,jk))-pdpmel(jl,jk)
            if(pmflxr(jl,jk+1)+pmflxs(jl,jk+1).lt.0.0) then
              pdmfdp(jl,jk)=-(pmflxr(jl,jk)+pmflxs(jl,jk)+pdmfup(jl,jk))
              pmflxr(jl,jk+1)=0.0
              pmflxs(jl,jk+1)=0.0
              pdpmel(jl,jk)  =0.0
            else if ( pmflxr(jl,jk+1) < 0. ) then
              pmflxs(jl,jk+1) = pmflxs(jl,jk+1)+pmflxr(jl,jk+1)
              pmflxr(jl,jk+1) = 0.
            else if ( pmflxs(jl,jk+1) < 0. ) then
              pmflxr(jl,jk+1) = pmflxr(jl,jk+1)+pmflxs(jl,jk+1)
              pmflxs(jl,jk+1) = 0.
            end if
          endif
        enddo
      enddo
      do jk=ktopm2,klev
       do jl = 1, nxj         
          if(ldcum(jl).and.jk.ge.kcbot(jl)) then
            zrfl=pmflxr(jl,jk)+pmflxs(jl,jk)
            if(zrfl.gt.1.e-20) then
              zdrfl1=zcpecons*max(0.,pqsen(jl,jk)-pqen(jl,jk))*zcucov* &
     &         (sqrt(paph(jl,jk)/paph(jl,klev+1))/5.09e-3*             &
     &         zrfl/zcucov)**0.5777*                                   &
     &         (paph(jl,jk+1)-paph(jl,jk))
              zrnew=zrfl-zdrfl1
              zrmin=zrfl-zcucov*max(0.,rhevap(jl)*pqsen(jl,jk) &
     &              -pqen(jl,jk)) *zcons2*(paph(jl,jk+1)-paph(jl,jk))
              zrnew=max(zrnew,zrmin)
              zrfln=max(zrnew,0.)
              zdrfl=min(0.,zrfln-zrfl)
              zdenom=1./max(1.e-20,pmflxr(jl,jk)+pmflxs(jl,jk))
              zalfaw=foealfa(pten(jl,jk))
              if ( pten(jl,jk) < tmelt ) zalfaw = 0.
              zpdr=zalfaw*pdmfdp(jl,jk)
              zpds=(1.-zalfaw)*pdmfdp(jl,jk)
              pmflxr(jl,jk+1)=pmflxr(jl,jk)+zpdr         &
     &         +pdpmel(jl,jk)+zdrfl*pmflxr(jl,jk)*zdenom
              pmflxs(jl,jk+1)=pmflxs(jl,jk)+zpds         &
     &         -pdpmel(jl,jk)+zdrfl*pmflxs(jl,jk)*zdenom
              pdmfup(jl,jk)=pdmfup(jl,jk)+zdrfl
              if ( pmflxr(jl,jk+1)+pmflxs(jl,jk+1) < 0. ) then
                pdmfup(jl,jk) = pdmfup(jl,jk)-(pmflxr(jl,jk+1)+pmflxs(jl,jk+1))
                pmflxr(jl,jk+1) = 0.
                pmflxs(jl,jk+1) = 0.
                pdpmel(jl,jk)   = 0.
              else if ( pmflxr(jl,jk+1) < 0. ) then
                pmflxs(jl,jk+1) = pmflxs(jl,jk+1)+pmflxr(jl,jk+1)
                pmflxr(jl,jk+1) = 0.
              else if ( pmflxs(jl,jk+1) < 0. ) then
                pmflxr(jl,jk+1) = pmflxr(jl,jk+1)+pmflxs(jl,jk+1)
                pmflxs(jl,jk+1) = 0.
              end if
            else
              pmflxr(jl,jk+1)=0.0
              pmflxs(jl,jk+1)=0.0
              pdmfdp(jl,jk)=0.0
              pdpmel(jl,jk)=0.0
            endif
          endif
        enddo
      enddo
                                      
      return
      end subroutine cuflx_n 
