  subroutine check_global_values(part,k,vname)
    use param
    use mpe
    use rank
    use index
    use const, only: RTYPE
    use mod_stochastic_physics, only : avevar_sppt2d
    implicit none

    real(kind=RTYPE) :: glob(nx,my)
    real :: aves,vars,stds

    real(kind=RTYPE), intent(in) :: part(nxp,my_max)
    integer, intent(in) :: k
    character(len=*), intent(in) :: vname
    character(len=120) :: fml

    call unify_reduceintp(nx,my,my_max,part,glob)
    if (myrank.eq.0) then
      write(fml,'(a2,i2.2,a16)')'(a',len(vname),',i3,3(1x,e15.6))'
      call avevar_sppt2d(glob,nx,my,aves,vars,stds)
      write(6,trim(fml))trim(vname),k,aves,vars,stds
    endif

  end subroutine check_global_values

