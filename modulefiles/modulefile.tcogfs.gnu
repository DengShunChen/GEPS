#%Module####################################################
## TCoGFS thin module for GEPS_COMPILER=gnu
## Expects GEPS_LIB geps/1.0 (or equivalent) already loaded.
##__________________________________________________________

module-whatis "TCoGFS build env aliases for GNU/OpenMPI (via GEPS_LIB)"

setenv GEPS_COMPILER gnu

if {[info exists ::env(MPIFC)] && $::env(MPIFC) ne ""} {
  setenv COMP_MP $::env(MPIFC)
} else {
  setenv COMP_MP mpifort
}
if {[info exists ::env(MPICC)] && $::env(MPICC) ne ""} {
  setenv C_COMP_MP $::env(MPICC)
} else {
  setenv C_COMP_MP mpicc
}
if {[info exists ::env(MPICXX)] && $::env(MPICXX) ne ""} {
  setenv CXX_COMP_MP $::env(MPICXX)
} else {
  setenv CXX_COMP_MP mpicxx
}
