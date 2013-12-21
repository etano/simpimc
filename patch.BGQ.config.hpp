--- config.hpp	2013-11-03 17:30:39.145454000 -0800
+++ ../../../src/armadillo-3.920.2/include/armadillo_bits/config.hpp	2013-11-03 17:32:10.609758000 -0800
@@ -22,15 +22,15 @@
 //// Without BLAS, matrix multiplication will still work, but might be slower.
 #endif
 
-/* #undef ARMA_USE_WRAPPER */
+#undef ARMA_USE_WRAPPER
 //// Comment out the above line if you're getting linking errors when compiling your programs,
 //// or if you prefer to directly link with LAPACK and/or BLAS.
 //// You will then need to link your programs directly with -llapack -lblas instead of -larmadillo
 
-// #define ARMA_BLAS_CAPITALS
+// #define ARMA_BLAS_CAPITALS
 //// Uncomment the above line if your BLAS and LAPACK libraries have capitalised function names (eg. ACML on 64-bit Windows)
 
-#define ARMA_BLAS_UNDERSCORE
+//#define ARMA_BLAS_UNDERSCORE
 //// Uncomment the above line if your BLAS and LAPACK libraries have function names with a trailing underscore.
 //// Conversely, comment it out if the function names don't have a trailing underscore.
 
