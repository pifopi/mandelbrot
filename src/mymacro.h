/* ----------------- */
/* --- mymacro.h --- */
/* ----------------- */


#ifndef __MY_MACRO_H__
#define __MY_MACRO_H__

// _WIN32 = we're in windows
#ifdef _WIN32
#include <intrin.h>
#else
#include <x86intrin.h>
#endif
//change for MSVC
//#include <x86intrin.h>

// activer ou desactiver le define ci-dessous pour passer du mode de mise au point
// au mode benchmark

#define ENABLE_BENCHMARK

// -------------------------------------------
// -- ne rien ecrire en dessous de cette ligne
// -------------------------------------------

#ifdef ENABLE_BENCHMARK

#pragma message("ENABLE_BENCHMARK is ON")

//#define CHRONO(X,t)  tmin = 1e38; for(run=0; run<nrun; run++) { t0 = (double) _rdtsc(); for(iter=0; iter<niter; iter++) { X; } t1 = (double) _rdtsc(); dt=t1-t0; if(dt<tmin) tmin = dt; } t = tmin / (double) niter

#define CHRONO(X,t)                       \
    tmin = 1e38;                          \
    for(run=0; run<nrun; run++) {         \
        t0 = (double) __rdtsc();           \
        for(iter=0; iter<niter; iter++) { \
            X;                            \
        }                                 \
        t1 = (double) __rdtsc();           \
        dt=t1-t0; if(dt<tmin) tmin = dt;  \
    }                                     \
    t = tmin / (double) niter


#define BENCH(X) X
#define DEBUG(X)
#define VERBOSE(X)

#else

#pragma message("ENABLE_BENCHMARK is OFF")

#define CHRONO(X,t)  X
#define BENCH(X)
#define DEBUG(X) X
//#define VERBOSE(X) X
#define VERBOSE(X)
#endif

#endif // __MY_MACRO_H__
