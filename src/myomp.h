#ifdef _OPENMP
  #include <omp.h>
#else
  // for machines with compilers void of openmp support
  #define omp_set_dynamic(a)      0;
  #define omp_set_num_threads(a)  1;
#endif
