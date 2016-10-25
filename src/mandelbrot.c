/* -------------------- */
/* --- mandelbrot.c --- */
/* -------------------- */

#include <math.h>

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <ctype.h> /* isdigit */
#include <string.h> /* memcpy */

#ifdef OPENMP
#include <omp.h>
#endif

#include "nrdef.h"
#include "nrutil.h"

#include "vnrdef.h"
#include "vnrutil.h"

#include "mutil.h"
#include "simd_macro.h"
#include "mymacro.h"

//#include "ia32intrin.h" // si compilateur Intel

#define OPENMP
// #define OMP_PARAM 
// #define OMP_PARAM schedule (static)
#define OMP_PARAM schedule (dynamic)
// #define OMP_PARAM schedule (dynamic, 30)
// #define OMP_PARAM schedule (guided)
// #define OMP_PARAM schedule (runtime)
// #define OMP_PARAM for schedule (auto)

// --------------------------------------------------
int mandelbrot_scalar (float a, float b, int max_iter)
// --------------------------------------------------
{
	int iter = 0;
	float x = 0;
	float y = 0;

	while (iter != max_iter)
	{
		float x_temp = x * x;
		float y_temp = y * y;
		float mult_temp = x * y;

		x = x_temp - y_temp + a;
		y = 2 * mult_temp + b;

		float z = x * x + y * y;

		if (z > 4)
		{
			break;
		}

		iter++;
	}
	return iter;
}
// ------------------------------
void test_mandelbrot_scalar (void)
// ------------------------------
{
	int iter, max_iter;
	float a, b;

	DEBUG (puts ("------------------------------"));
	DEBUG (puts ("--- test_mandelbrot_scalar ---"));
	DEBUG (puts ("------------------------------"));
	// tests unitaire pour valider
	max_iter = 15;
	DEBUG (printf ("max_iter\t%d\n", max_iter));
	a = -0.8; b = +0.3; iter = mandelbrot_scalar (a, b, max_iter); DEBUG (printf ("(%4.2f %4.2f) -> %3d\n", a, b, iter));
	a = -0.7; b = +0.2; iter = mandelbrot_scalar (a, b, max_iter); DEBUG (printf ("(%4.2f %4.2f) -> %3d\n", a, b, iter));
	a = -0.8; b = -0.3; iter = mandelbrot_scalar (a, b, max_iter); DEBUG (printf ("(%4.2f %4.2f) -> %3d\n", a, b, iter));
	a = -0.7; b = -0.2; iter = mandelbrot_scalar (a, b, max_iter); DEBUG (printf ("(%4.2f %4.2f) -> %3d\n", a, b, iter));
	DEBUG (puts (""));

	max_iter = 20;
	DEBUG (printf ("max_iter\t%d\n", max_iter));
	a = -1.00; b = 0.40; iter = mandelbrot_scalar (a, b, max_iter); DEBUG (printf ("(%4.2f %4.2f) -> %3d\n", a, b, iter));
	a = -0.90; b = 0.30; iter = mandelbrot_scalar (a, b, max_iter); DEBUG (printf ("(%4.2f %4.2f) -> %3d\n", a, b, iter));
	a = -0.80; b = 0.30; iter = mandelbrot_scalar (a, b, max_iter); DEBUG (printf ("(%4.2f %4.2f) -> %3d\n", a, b, iter));
	a = -0.70; b = 0.10; iter = mandelbrot_scalar (a, b, max_iter); DEBUG (printf ("(%4.2f %4.2f) -> %3d\n", a, b, iter));
	DEBUG (puts (""));
}
// --------------------------------------------------------------
vuint32 mandelbrot_SIMD_F32 (vfloat32 a, vfloat32 b, int max_iter)
// --------------------------------------------------------------
{
	// version avec test de sortie en float
	vfloat32 iter = _mm_set_ps1 (0);
	vfloat32 x = _mm_set_ps1 (0);// On initialise x0
	vfloat32 y = _mm_set_ps1 (0);// On initialise y0

	int i = 0;
	vfloat32 seuil = _mm_set_ps1 (4);// On initialise la condition d'arrêt
	vfloat32 un = _mm_set_ps1 (1);// On initialise l'incrémentation qui sera masquée par la suite
	vfloat32 incr = _mm_set_ps1 (1);

	while (i != max_iter)
	{
		vfloat32 x_temp = _mm_mul_ps (x, x);// x_n²
		vfloat32 y_temp = _mm_mul_ps (y, y);// y_n²
		vfloat32 mult_temp = _mm_mul_ps (x, y);// x_n * y_n

		x = _mm_add_ps (_mm_sub_ps (x_temp, y_temp), a);// x_n² - y_n² + a
		y = _mm_add_ps (_mm_add_ps (mult_temp, mult_temp), b);// x_n * y_n + x_n * y_n + b

		vfloat32 z = _mm_add_ps (_mm_mul_ps (x, x), _mm_mul_ps (y, y));// x_n+1² + y_n+1²

		vfloat32 masque = _mm_cmple_ps (z, seuil);// On regarde quelles sont les valeurs encore inférieures au seuil
		if (_mm_movemask_ps (masque) == 0)// Si il n'y a aucune valeur inférieure au seuil, nous quittons la boucle
		{
			break;
		}

		iter = _mm_add_ps (iter, incr);// On incrémente les valeurs d'iter nécessaires
		i++;

		incr = _mm_and_ps (masque, un);// On crée le vecteur d'incrémentation pour le prochain tour de boucle
	}

	// Si on n'a pas atteint max_iter, nous sommes allés un cran trop loin, il faut donc décrementer d'un les valeurs trop augmentées
	vfloat32 masque = _mm_cmplt_ps (iter, _mm_set_ps1 (max_iter));
	incr = _mm_and_ps (masque, _mm_set_ps1 (1));
	iter = _mm_sub_ps (iter, incr);
	return _mm_cvttps_epi32 (iter);
}
// --------------------------------------------------------------
vuint32 mandelbrot_SIMD_I32 (vfloat32 a, vfloat32 b, int max_iter)
// --------------------------------------------------------------
{
	// version avec test de sortie en int
	vuint32 iter = _mm_set1_epi32 (0);
	vfloat32 x = _mm_set_ps1 (0);// On initialise x0
	vfloat32 y = _mm_set_ps1 (0);// On initialise y0

	int i = 0;
	vfloat32 seuil = _mm_set_ps1 (4);// On initialise la condition d'arrêt
	vuint32 un = _mm_set1_epi32 (1);// On initialise l'incrémentation qui sera masquée par la suite
	vuint32 incr = _mm_set1_epi32 (1);

	while (i != max_iter)
	{
		vfloat32 x_temp = _mm_mul_ps (x, x);// x_n²
		vfloat32 y_temp = _mm_mul_ps (y, y);// y_n²
		vfloat32 mult_temp = _mm_mul_ps (x, y);// x_n * y_n

		x = _mm_add_ps (_mm_sub_ps (x_temp, y_temp), a);// x_n² - y_n² + a
		y = _mm_add_ps (_mm_add_ps (mult_temp, mult_temp), b);// x_n * y_n + x_n * y_n + b

		vfloat32 z = _mm_add_ps (_mm_mul_ps (x, x), _mm_mul_ps (y, y));// x_n+1² + y_n+1²

		vfloat32 masque = _mm_cmple_ps (z, seuil);// On regarde quelles sont les valeurs encore inférieures au seuil
		if (_mm_movemask_ps (masque) == 0)// Si il n'y a aucune valeur inférieure au seuil, nous quittons la boucle
		{
			break;
		}

		iter = _mm_add_epi32 (iter, incr);// On incrémente les valeurs d'iter nécessaires
		i++;

		//incr = _mm_castps_si128(masque) & un;// On crée le vecteur d'incrémentation pour le prochain tour de boucle  BEWARE MSVC
		incr = _mm_and_si128 (_mm_castps_si128 (masque), un);
	}

	// Si on n'a pas atteint max_iter, nous sommes allés un cran trop loin, il faut donc décrementer d'un les valeurs trop augmentées
	vuint32 masque = _mm_cmplt_epi32 (iter, _mm_set1_epi32 (max_iter));
	incr = _mm_and_si128 (masque, _mm_set1_epi32 (1));
	iter = _mm_sub_epi32 (iter, incr);
	return iter;
}
// ----------------------------
void test_mandelbrot_SIMD (void)
// ----------------------------
{
	int max_iter = 20;
	vuint32 iter;
	vfloat32 a, b;

	DEBUG (puts ("----------------------------"));
	DEBUG (puts ("--- test_mandelbrot_SIMD ---"));
	DEBUG (puts ("----------------------------"));

	DEBUG (puts ("mandelbrot_SIMD_F32"));
	a = _mm_setr_ps (-0.8, -0.7, -0.8, -0.7);
	b = _mm_setr_ps (+0.3, +0.2, -0.3, -0.2);
	max_iter = 15;
	iter = mandelbrot_SIMD_F32 (a, b, max_iter);
	DEBUG (printf ("max_iter\t%d\n", max_iter));
	DEBUG (display_vfloat32 (a, "%10.2f ", "a\t")); DEBUG (puts (""));
	DEBUG (display_vfloat32 (b, "%10.2f ", "b\t")); DEBUG (puts (""));
	DEBUG (display_vuint32 (iter, "%10d ", "iter")); DEBUG (puts (""));
	DEBUG (puts (""));

	a = _mm_setr_ps (-1.00, -0.90, -0.80, -0.70);
	b = _mm_setr_ps (+0.40, +0.30, +0.30, +0.10);
	max_iter = 20;
	iter = mandelbrot_SIMD_F32 (a, b, max_iter);
	DEBUG (printf ("max_iter\t%d\n", max_iter));
	DEBUG (display_vfloat32 (a, "%10.2f ", "a\t")); DEBUG (puts (""));
	DEBUG (display_vfloat32 (b, "%10.2f ", "b\t")); DEBUG (puts (""));
	DEBUG (display_vuint32 (iter, "%10d ", "iter")); DEBUG (puts (""));
	DEBUG (puts (""));

	DEBUG (puts ("mandelbrot_SIMD_I32"));
	a = _mm_setr_ps (-0.8, -0.7, -0.8, -0.7);
	b = _mm_setr_ps (+0.3, +0.2, -0.3, -0.2);
	max_iter = 15;
	iter = mandelbrot_SIMD_I32 (a, b, max_iter);
	DEBUG (printf ("max_iter\t%d\n", max_iter));
	DEBUG (display_vfloat32 (a, "%10.2f ", "a\t")); DEBUG (puts (""));
	DEBUG (display_vfloat32 (b, "%10.2f ", "b\t")); DEBUG (puts (""));
	DEBUG (display_vuint32 (iter, "%10d ", "iter")); DEBUG (puts (""));
	DEBUG (puts (""));

	a = _mm_setr_ps (-1.00, -0.90, -0.80, -0.70);
	b = _mm_setr_ps (+0.40, +0.30, +0.30, +0.10);
	max_iter = 20;
	iter = mandelbrot_SIMD_I32 (a, b, max_iter);
	DEBUG (printf ("max_iter\t%d\n", max_iter));
	DEBUG (display_vfloat32 (a, "%10.2f ", "a\t")); DEBUG (puts (""));
	DEBUG (display_vfloat32 (b, "%10.2f ", "b\t")); DEBUG (puts (""));
	DEBUG (display_vuint32 (iter, "%10d ", "iter")); DEBUG (puts (""));
	DEBUG (puts (""));
}

// --------------------------------------------------------------------------------------------------------
void calc_mandelbrot_scalar (uint32 **M, int h, int w, float a0, float a1, float b0, float b1, int max_iter)
// --------------------------------------------------------------------------------------------------------
{
	// intervale de valeurs: [a0:a1]x[b0:b1]

	// la seule chose a modifier dans cette fonction est la ligne de pragma OpenMP
	float da = (a1 - a0) / w;
	float db = (b1 - b0) / h;

	int i;
	int j;

#ifdef OPENMP
#pragma omp parallel for private (i, j) OMP_PARAM
#endif

	for (i = 0; i<h; i++) {
		for (j = 0; j<w; j++) {

			// conversion (i,j) -> (x,y)
			float32 a = a0 + i * da;
			float32 b = b0 + j * db;

			uint32 iter = mandelbrot_scalar (a, b, max_iter);

			M[i][j] = iter;
		}
	}
}
// -----------------------------------------------------------------------------------------------------------
void calc_mandelbrot_SIMD_F32 (vuint32 **M, int h, int w, float a0, float a1, float b0, float b1, int max_iter)
// -----------------------------------------------------------------------------------------------------------
{
	float da = (a1 - a0) / w;
	float db = (b1 - b0) / h;

	int i;
	int j;

#ifdef OPENMP
#pragma omp parallel for private (i, j) OMP_PARAM
#endif
	for (i = 0; i<h; i++) {
		for (j = 0; j<w / 4; j++) {

			// conversion (i,j) -> (x,y)
			float sa = a0 + i * da;
			float sb = b0 + j * db * 4;

			vfloat32 b = _mm_setr_ps (sb, sb + db, sb + 2 * db, sb + 3 * db);
			vfloat32 a = _mm_set1_ps (sa);

			vuint32 iter = mandelbrot_SIMD_F32 (a, b, max_iter);
			M[i][j] = iter;
		}
	}
}
// -----------------------------------------------------------------------------------------------------------
void calc_mandelbrot_SIMD_I32 (vuint32 **M, int h, int w, float a0, float a1, float b0, float b1, int max_iter)
// -----------------------------------------------------------------------------------------------------------
{
	float da = (a1 - a0) / w;
	float db = (b1 - b0) / h;

	int i;
	int j;

#ifdef OPENMP
#pragma omp parallel for private (i, j) OMP_PARAM
#endif

	for (i = 0; i<h; i++) {
		for (j = 0; j<w / 4; j++) {

			// conversion (i,j) -> (x,y)
			float sa = a0 + i * da;
			float sb = b0 + j * db * 4;

			vfloat32 b = _mm_setr_ps (sb, sb + db, sb + 2 * db, sb + 3 * db);
			vfloat32 a = _mm_set1_ps (sa);

			vuint32 iter = mandelbrot_SIMD_I32 (a, b, max_iter);
			M[i][j] = iter;
		}
	}
}
// -----------------------------------------------------------------------------------
convert_ui32matrix_ui8matrix (uint32 **m32, int i0, int i1, int j0, int j1, uint8 **m8)
// -----------------------------------------------------------------------------------
{
	int i, j;
	for (i = i0; i <= i1; i++) {
		for (j = j0; j <= j1; j++) {
			m8[i][j] = (uint8)m32[i][j];
		}
	}
}
// ----------------------------------------------
void bench_mandelbrot_scalar (int n, int max_iter)
// ----------------------------------------------
{
	// ne rien modifier dans cette fonction

	int h, w;
	int i0, i1, j0, j1;
	float a0, a1, b0, b1;
	uint32 **M32;
	uint8 **M8;


	// chronometrie
	int iter, niter = 4;
	int run, nrun = 5;
	double t0, t1, dt, tmin, t;
	double cycles;

	//puts("-----------------------------");
	//puts("-- bench_mandelbrot_scalar --");
	//puts("-----------------------------");

	h = w = n;
	M32 = ui32matrix (0, h - 1, 0, w - 1);
	M8 = ui8matrix (0, h - 1, 0, w - 1);

	a0 = -1.5; a1 = +0.5;
	b0 = -1.0; b1 = +1.0;

	CHRONO (calc_mandelbrot_scalar (M32, h, w, a0, a1, b0, b1, max_iter), cycles);  printf ("scalar:   %10.2f\n", cycles / (n*n));

	DEBUG (convert_ui32matrix_ui8matrix (M32, 0, h - 1, 0, w - 1, M8));
	DEBUG (SavePGM_ui8matrix (M8, 0, h - 1, 0, w - 1, "M_scalar.pgm"));

	free_ui32matrix (M32, 0, h - 1, 0, w - 1);
	free_ui8matrix (M8, 0, h - 1, 0, w - 1);
}
// --------------------------------------------
void bench_mandelbrot_SIMD (int n, int max_iter)
// --------------------------------------------
{
	// ne rien modifier dans cette fonction

	int h, w;
	float a0, a1, b0, b1;
	vuint32 **M32;
	uint32 **wM32;
	uint8 **M8;

	// chronometrie
	int iter, niter = 4;
	int run, nrun = 5;
	double t0, t1, dt, tmin, t;
	double cycles;

	//puts("---------------------------");
	//puts("-- bench_mandelbrot_SIMD --");
	//puts("---------------------------");
	h = w = n;

	M32 = vui32matrix (0, h - 1, 0, w / 4 - 1);
	M8 = ui8matrix (0, h - 1, 0, w - 1);
	wM32 = (uint32**)M32;

	// ne pas changer
	a0 = -1.5; a1 = +0.5;
	b0 = -1.0; b1 = +1.0;

	CHRONO (calc_mandelbrot_SIMD_F32 (M32, h, w, a0, a1, b0, b1, max_iter), cycles);
	printf ("SIMD F32: %10.2f\n", cycles / (n*n));

	CHRONO (calc_mandelbrot_SIMD_I32 (M32, h, w, a0, a1, b0, b1, max_iter), cycles); // facultatif
	printf ("SIMD I32: %10.2f\n\n", cycles / (n*n));

	DEBUG (convert_ui32matrix_ui8matrix (wM32, 0, h - 1, 0, w - 1, M8));
	DEBUG (SavePGM_ui8matrix (M8, 0, h - 1, 0, w - 1, "M_v.pgm"));

	free_vui32matrix (M32, 0, h - 1, 0, w / 4 - 1);
	free_ui8matrix (M8, 0, h - 1, 0, w - 1);
}

// =========================================
int main_mandelbrot (int argc, char * argv[])
// =========================================
{
	int n, max_iter; // pour avoir les meme param en scalar et SIMD ...

	test_mandelbrot_scalar ();
	test_mandelbrot_SIMD ();

	// n = 512; max_iter = 256;
	n = 1024; max_iter = 256;
	// n = 2048; max_iter = 256;
	// n = 4096; max_iter = 256;
	// n = 8192; max_iter = 256;

	printf ("n = %4d max_iter = %d\n", n, max_iter);
	bench_mandelbrot_scalar (n, max_iter);
	bench_mandelbrot_SIMD (n, max_iter);

	return 0;
}
