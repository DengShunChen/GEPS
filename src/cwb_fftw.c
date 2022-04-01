#include <stdio.h>
#include <stdlib.h>
#include "fftw3.h"
/*
  call dfftw_init_threads( )
  call dfftw_plan_with_nthreads()
  call dfftw_make_planner_thread_safe( )

  the above calls have linking trouble when calling in fortran,
  so I make this c code for fortran to call.
*/

/* ------------------------------------------------------- */
void cwb_fftw_init_threads_( )
{
  fftw_init_threads( );
  return ;
}

/* ------------------------------------------------------- */
void cwb_fftw_plan_with_nthreads_( )
{
#include <omp.h>
  fftw_plan_with_nthreads(omp_get_max_threads());
/*
  fftw_plan_with_nthreads( );
*/
  return ;
}

/* ------------------------------------------------------- */
void cwb_fftw_make_planner_thread_safe_( )
{
  fftw_make_planner_thread_safe( );
  return ;
}

/* ------------------------------------------------------- */
void cwb_fftw_import_wisdom_from_filename_(char *cwbfile)
{
   int is;
   is=fftw_import_wisdom_from_filename(cwbfile);
   if(is == 1){
     printf("import wisdom file OK !\n");
   }
   return ;
}

/* ------------------------------------------------------- */
void cwb_fftw_export_wisdom_to_filename_(char *cwbfile)
{
   int is;
   is=fftw_export_wisdom_to_filename(cwbfile);
   if(is == 1){
     printf("export wisdom file OK !\n");
   }
   return ;
}

/* ------------------------------------------------------- */
void cwb_fftw_cleanup_threads_( )
{
  fftw_cleanup_threads( );
  return ;
}
