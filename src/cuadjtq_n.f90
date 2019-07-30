
      subroutine cuadjtq_n                 &                           
     &    (nxj, klon, klev, kk, psp, pt, pq, ldflag,  kcall)                
!          m.tiedtke         e.c.m.w.f.     12/89     
!          purpose.                                                                 
!          --------                                                                 
!          to produce t,q and l values for cloud ascent                             
                                                                       
!          interface                                                                
!          ---------                                                                
!          this routine is called from subroutines:                                 
!              *cond*     (t and q at condensation level)                           
!              *cubase*   (t and q at condensation level)                           
!              *cuasc*    (t and q at cloud levels)                                 
!              *cuini*    (environmental t and qs values at half levels)            
!          input are unadjusted t and q values,                                     
!          it returns adjusted values of t and q                                    
                                                                        
!     parameter     description                                   units             
!     ---------     -----------                                   -----             
!     input parameters (integer):                                                   
                                                                       
!    *klon*         number of grid points per packet                                
!    *klev*         number of levels                                                
!    *kk*           level                                                           
!    *kcall*        defines calculation as                                          
!                      kcall=0  env. t and qs in*cuini*                             
!                      kcall=1  condensation in updrafts  (e.g. cubase, cuasc)      
!                      kcall=2  evaporation in downdrafts (e.g. cudlfs,cuddraf)     
!     input parameters (real):                                                      
!    *psp*          pressure                                        pa              
!     updated parameters (real):                                                    
!    *pt*           temperature                                     k               
!    *pq*           specific humidity                             kg/kg             
!          externals                                                                
!          ---------                                                                
!          for condensation calculations.                                           
!          the tables are initialised in *suphec*.                                  
!----------------------------------------------------------------------             
  USE mo_constants,      ONLY: vtmpc1,   &
                               r5alvcp,  &
                               r5alscp,  &
                               c4les,    &
                               c4ies,    &
                               c3ies,    &
                               tmelt,    &
                               c3les,    &
                               c2es   
  USE mo_convect_tables, ONLY: tlucua,   & ! table a
                               tlucub,   & ! table b
                               tlucuc,   & ! table c
                               jptlucu1, jptlucu2, &
                               lookuperror, lookupoverflow

      implicit none
 
      integer  klev,klon   
                 
      real     pt(klon,klev),          pq(klon,klev),  &               
     &         psp(klon)                                               
      logical  ldflag(klon)                                           
! local variables
      integer  jl,jk,nxj           
      integer  isum,kcall,kk
      real     zqmax,zqsat,zcor,zqp,zcond,zcond1,zl,zi,zf
      real     foealfa,foeldcpm,foeewm,foedem     
!------------- ---------------------------------------------------------             
!     1.           define constants                                                 
!                  ----------------                                                 
      zqmax=0.5  
!     2.           calculate condensation and adjust t and q accordingly            
!                  -----------------------------------------------------            
      if ( kcall == 1 ) then
       do jl = 1, nxj           
        if ( ldflag(jl) ) then
          zqp = 1./psp(jl)
          zl = 1./(pt(jl,kk)-c4les)
          zi = 1./(pt(jl,kk)-c4ies)
          zqsat = c2es*(foealfa(pt(jl,kk))*exp(c3les*(pt(jl,kk)-tmelt)*zl) + &
                (1.-foealfa(pt(jl,kk)))*exp(c3ies*(pt(jl,kk)-tmelt)*zi))
          zqsat = zqsat*zqp
          zqsat = min(0.5,zqsat)
          zcor = 1. - vtmpc1*zqsat
          zf = foealfa(pt(jl,kk))*r5alvcp*zl**2 + &
               (1.-foealfa(pt(jl,kk)))*r5alscp*zi**2
          zcond = (pq(jl,kk)*zcor**2-zqsat*zcor)/(zcor**2+zqsat*zf)
          if ( zcond > 0. ) then
            pt(jl,kk) = pt(jl,kk) + foeldcpm(pt(jl,kk))*zcond
            pq(jl,kk) = pq(jl,kk) - zcond
            zl = 1./(pt(jl,kk)-c4les)
            zi = 1./(pt(jl,kk)-c4ies)
            zqsat = c2es*(foealfa(pt(jl,kk)) * &
              exp(c3les*(pt(jl,kk)-tmelt)*zl)+(1.-foealfa(pt(jl,kk))) * &
              exp(c3ies*(pt(jl,kk)-tmelt)*zi))
            zqsat = zqsat*zqp
            zqsat = min(0.5,zqsat)
            zcor = 1. - vtmpc1*zqsat
            zf = foealfa(pt(jl,kk))*r5alvcp*zl**2 + &
                 (1.-foealfa(pt(jl,kk)))*r5alscp*zi**2
            zcond1 = (pq(jl,kk)*zcor**2-zqsat*zcor)/(zcor**2+zqsat*zf)
            if ( abs(zcond) < 1.e-20 ) zcond1 = 0.
            pt(jl,kk) = pt(jl,kk) + foeldcpm(pt(jl,kk))*zcond1
            pq(jl,kk) = pq(jl,kk) - zcond1
          end if
        end if
      end do
      elseif ( kcall == 2 ) then
      do jl = 1, nxj           
        if ( ldflag(jl) ) then
          zqp = 1./psp(jl)
          zqsat = foeewm(pt(jl,kk))*zqp
          zqsat = min(0.5,zqsat)
          zcor = 1./(1.-vtmpc1*zqsat)
          zqsat = zqsat*zcor
          zcond = (pq(jl,kk)-zqsat)/(1.+zqsat*zcor*foedem(pt(jl,kk)))
          zcond = min(zcond,0.)
          pt(jl,kk) = pt(jl,kk) + foeldcpm(pt(jl,kk))*zcond
          pq(jl,kk) = pq(jl,kk) - zcond
          zqsat = foeewm(pt(jl,kk))*zqp
          zqsat = min(0.5,zqsat)
          zcor = 1./(1.-vtmpc1*zqsat)
          zqsat = zqsat*zcor
          zcond1 = (pq(jl,kk)-zqsat)/(1.+zqsat*zcor*foedem(pt(jl,kk)))
          if ( abs(zcond) < 1.e-20 ) zcond1 = min(zcond1,0.)
          pt(jl,kk) = pt(jl,kk) + foeldcpm(pt(jl,kk))*zcond1
          pq(jl,kk) = pq(jl,kk) - zcond1
        end if
      end do
      else if ( kcall == 0 ) then
       do jl = 1, nxj              
        zqp = 1./psp(jl)
        zqsat = foeewm(pt(jl,kk))*zqp
        zqsat = min(0.5,zqsat)
        zcor = 1./(1.-vtmpc1*zqsat)
        zqsat = zqsat*zcor
        zcond1 = (pq(jl,kk)-zqsat)/(1.+zqsat*zcor*foedem(pt(jl,kk)))
        pt(jl,kk) = pt(jl,kk) + foeldcpm(pt(jl,kk))*zcond1
        pq(jl,kk) = pq(jl,kk) - zcond1
        zqsat = foeewm(pt(jl,kk))*zqp
        zqsat = min(0.5,zqsat)
        zcor = 1./(1.-vtmpc1*zqsat)
        zqsat = zqsat*zcor
        zcond1 = (pq(jl,kk)-zqsat)/(1.+zqsat*zcor*foedem(pt(jl,kk)))
        pt(jl,kk) = pt(jl,kk) + foeldcpm(pt(jl,kk))*zcond1
        pq(jl,kk) = pq(jl,kk) - zcond1
      end do
      end if

      return
      end subroutine cuadjtq_n
