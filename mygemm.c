#include "mygemm.h"

/**
 * 
 * Implement all functions here in this file.
 * Do NOT change input parameters and return type.
 * 
 **/

void dgemm0(const double* A, const double* B, double* C, const int n)
{
printf( "n size %d\n",n); 
printf( "dgemm0");  
for (int i = 0; i < n; i++)
	for (int j = 0; j < n; j++)
		for (int k = 0; k < n; k++)
		C[i*n + j] += A[i*n + k] * B[k*n + j]; 

}

void dgemm1(const double *A, const double *B, double *C, const int n) 
{
printf( "n size %d\n",n); 
printf( "Time of dgemm1 ...");  
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++) 
		{
			register double r = C[i*n + j];
			for (int k = 0; k < n; k++)
				r += A[i*n + k] * B[k*n + j];
				C[i*n + j] = r;
		}
	}
}

void dgemm2(const double *A, const double *B, double *C, const int n) 
{
	printf( "n size %d\n",n);  
	printf( "Time of dgemm2 ...");  
	for (int i = 0; i < n; i += 2)
    {
        for (int j = 0; j < n; j += 2)
        {
            register double c00 = C[i*n + j];
            register double c01 = C[i*n + j+1];
            register double c10 = C[(i+1)*n + j];
            register double c11 = C[(i+1)*n + j+1];

            for (int k = 0; k < n; k += 2)
            {
                register double a00 = A[i*n + k];
                register double a01 = A[i*n + k+1];
                register double a10 = A[(i+1)*n + k];
                register double a11 = A[(i+1)*n + k+1];
                register double b00 = B[k*n + j];
                register double b01 = B[k*n + j+1];
                register double b10 = B[(k+1)*n + j];
                register double b11 = B[(k+1)*n + j+1];

                c00 += a00 * b00 + a01 * b10;
                c01 += a00 * b01 + a01 * b11;
                c10 += a10 * b00 + a11 * b10;
                c11 += a10 * b01 + a11 * b11;
            }
			
            C[i*n + j] = c00;
            C[i*n + j+1] = c01;
            C[(i+1)*n + j] = c10;
            C[(i+1)*n + j+1] = c11;
        }
    }
}

void dgemm3(const double *A, const double *B, double *C, const int n) 
{
printf( "n size %d\n",n);  
printf( "Time of dgemm3 ...");  
for (int i = 0; i < n; i += 3)
    {
        for (int j = 0; j < n; j += 3)
        {
            register double c00 = C[i*n + j];
            register double c01 = C[i*n + j+1];
            register double c02 = C[i*n + j+2];
            register double c10 = C[(i+1)*n + j];
            register double c11 = C[(i+1)*n + j+1];
            register double c12 = C[(i+1)*n + j+2];
            register double c20 = C[(i+2)*n + j];
            register double c21 = C[(i+2)*n + j+1];
            register double c22 = C[(i+2)*n + j+2];

            for (int k = 0; k < n; k += 3)
            {
                for (int l = 0; l < 3; l++)
                {
                    register double a0 = A[i*n + k+l];
                    register double a1 = A[(i+1)*n + k+l];
                    register double a2 = A[(i+2)*n + k+l];
                    register double b0 = A[(k+l)*n + j];
                    register double b1 = B[(k+l)*n + j+1];
                    register double b2 = B[(k+l)*n + j+2];

                    c00 += a0 * b0;
                    c01 += a0 * b1;
                    c02 += a0 * b2;
                    c10 += a1 * b0;
                    c11 += a1 * b1;
                    c12 += a1 * b2;
                    c20 += a2 * b0;
                    c21 += a2 * b1;
                    c22 += a2 * b2;
                }
            }
            C[i*n + j] = c00;
            C[i*n + j+1] = c01;
            C[i*n + j+2] = c02;
            C[(i+1)*n + j] = c10;
            C[(i+1)*n + j+1] = c11;
            C[(i+1)*n + j+2] = c12;
            C[(i+2)*n + j] = c20;
            C[(i+2)*n + j+1] = c21;
            C[(i+2)*n + j+2] = c22;
        }
    }


}

void ijk(const double *A, const double *B, double *C, const int n) 
{
printf( "n size %d\n",n);  
printf( "Time of ijk...");  
	 for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            register double r = C[i * n + j];
            for (int k = 0; k < n; k++) {
                r += A[i * n + k] * B[k * n + j];
            }
            C[i * n + j] = r;
        }
    }


}

void bijk(const double *A, const double *B, double *C, const int n, const int b) 
{
printf( "n size %d\t",n); 
printf( "b size %d\n",b); 
printf( "Time of blockijk...");  
for (int i = 0; i < n; i += b) {
        for (int j = 0; j < n; j += b) {
            for (int k = 0; k < n; k += b) {
                for (int i1 = i; i1 < i + b && i1 < n; i1++) {
                    for (int j1 = j; j1 < j + b && j1 < n; j1++) {
                        register double r = C[i1 * n + j1];
                        for (int k1 = k; k1 < k + b && k1 < n; k1++) {
                            r += A[i1 * n + k1] * B[k1 * n + j1];
                        }
                        C[i1 * n + j1] = r;
                    }
                }
            }
        }
    }

}

void jik(const double *A, const double *B, double *C, const int n) 
{
printf( "n size %d\n",n);  
printf( "Time of jik...");  
	for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            register double r = C[i * n + j];
            for (int k = 0; k < n; k++) {
                r += A[i * n + k] * B[k * n + j];
            }
            C[i * n + j] = r;
        }
    }

}

void bjik(const double *A, const double *B, double *C, const int n, const int b) 
{
printf( "n size %d\t",n); 
printf( "b size %d\n",b); 
printf( "Time of blockjik...");  
for (int j = 0; j < n; j += b) {
        for (int i = 0; i < n; i += b) {
            for (int k = 0; k < n; k += b) {
                for (int j1 = j; j1 < j + b && j1 < n; j1++) {
                    for (int i1 = i; i1 < i + b && i1 < n; i1++) {
                        register double r = C[i1 * n + j1];
                        for (int k1 = k; k1 < k + b && k1 < n; k1++) {
                            r += A[i1 * n + k1] * B[k1 * n + j1];
                        }
                        C[i1 * n + j1] = r;
                    }
                }
            }
        }
    }

}

void kij(const double *A, const double *B, double *C, const int n) 
{
printf( "n size %d\n",n); 
printf( "Time of kij...");  
for (int k = 0; k < n; k++) {
        for (int i = 0; i < n; i++) {
            register double r = A[i * n + k];
            for (int j = 0; j < n; j++) {
                C[i * n + j] += r * B[k * n + j];
            }
        }
    }

}

void bkij(const double *A, const double *B, double *C, const int n, const int b) 
{
printf( "n size %d\t",n); 
printf( "b size %d\n",b); 
printf( "Time of blockkij...");  
for (int k = 0; k < n; k += b) {
        for (int i = 0; i < n; i += b) {
            for (int j = 0; j < n; j += b) {
                for (int k1 = k; k1 < k + b && k1 < n; k1++) {
                    for (int i1 = i; i1 < i + b && i1 < n; i1++) {
                        register double r = A[i1 * n + k1];
                        for (int j1 = j; j1 < j + b && j1 < n; j1++) {
                            C[i1 * n + j1] += r * B[k1 * n + j1];
                        }
                    }
                }
            }
        }
    }


}


void ikj(const double *A, const double *B, double *C, const int n) 
{
printf( "n size %d\n",n); 
printf( "Time of ikj...");  
	 for (int i = 0; i < n; i++) {
        for (int k = 0; k < n; k++) {
            register double r = A[i * n + k];
            for (int j = 0; j < n; j++) {
                C[i * n + j] += r * B[k * n + j];
            }
        }
    }
}

void bikj(const double *A, const double *B, double *C, const int n, const int b) 
{
printf( "n size %d\t",n); 
printf( "b size %d\n",b); 
printf( "Time of blockikj...");  
for (int i = 0; i < n; i += b) {
        for (int k = 0; k < n; k += b) {
            for (int j = 0; j < n; j += b) {
                for (int i1 = i; i1 < i + b && i1 < n; i1++) {
                    for (int k1 = k; k1 < k + b && k1 < n; k1++) {
                        register double r = A[i1 * n + k1];
                        for (int j1 = j; j1 < j + b && j1 < n; j1++) {
                            C[i1 * n + j1] += r * B[k1 * n + j1];
                        }
                    }
                }
            }
        }
    }

}

void jki(const double *A, const double *B, double *C, const int n) 
{
printf( "n size %d\n",n); 
printf( "Time of jki...");  
	for (int j = 0; j < n; j++) {
        for (int k = 0; k < n; k++) {
            register double r = B[k * n + j];
            for (int i = 0; i < n; i++) {
                C[i * n + j] += r * A[i * n + k];
            }
        }
    }

}

void bjki(const double *A, const double *B, double *C, const int n, const int b) 
{
printf( "n size %d\t",n); 
printf( "b size %d\n",b); 
printf( "Time of blockjki...");  
for (int j = 0; j < n; j += b) {
        for (int k = 0; k < n; k += b) {
            for (int i = 0; i < n; i += b) {
                for (int j1 = j; j1 < j + b && j1 < n; j1++) {
                    for (int k1 = k; k1 < k + b && k1 < n; k1++) {
                        register double r = B[k1 * n + j1];
                        for (int i1 = i; i1 < i + b && i1 < n; i1++) {
                            C[i1 * n + j1] += r * A[i1 * n + k1];
                        }
                    }
                }
            }
        }
    }

}

void kji(const double *A, const double *B, double *C, const int n) 
{
printf( "n size %d\n",n); 
printf( "Time of kji...");  
	for (int k = 0; k < n; k++) {
        for (int j = 0; j < n; j++) {
            register double r = B[k * n + j];
            for (int i = 0; i < n; i++) {
                C[i * n + j] += r * A[i * n + k];
            }
        }
    }

}

void bkji(const double *A, const double *B, double *C, const int n, const int b) 
{
printf( "n size %d\t",n); 
printf( "b size %d\n",b); 
printf( "Time of blockkji...");
for (int k = 0; k < n; k += b) {
        for (int j = 0; j < n; j += b) {
            for (int i = 0; i < n; i += b) {
                for (int k1 = k; k1 < k + b && k1 < n; k1++) {
                    for (int j1 = j; j1 < j + b && j1 < n; j1++) {
                        register double r = B[k1 * n + j1];
                        for (int i1 = i; i1 < i + b && i1 < n; i1++) {
                            C[i1 * n + j1] += r * A[i1 * n + k1];
                        }
                    }
                }
            }
        }
    }

}

void optimal(const double* A, const double* B, double *C, const int n, const int b)
{
printf( "n size %d\t",n); 
printf( "b size %d\n",b); 
printf( "Time of optimal");
  for (int i = 0; i < n; i += b) {
        for (int j = 0; j < n; j += b) {
            for (int k = 0; k < n; k += b) {
                int i1 = i, j1 = j, k1 = k;
                int ni = i + b > n ? n : i + b;
                int nj = j + b > n ? n : j + b;
                int nk = k + b > n ? n : k + b;

                for (i1 = i; i1 < ni; i1 += 3) {
                    for (j1 = j; j1 < nj; j1 += 3) {
                        int h = i1 * n + j1;
                        int hh = h + n;
                        int hhh = hh + n;
                        register double c00 = C[h];
                        register double c01 = C[h + 1];
                        register double c02 = C[h + 2];
                        register double c10 = C[hh];
                        register double c11 = C[hh + 1];
                        register double c12 = C[hh + 2];
                        register double c20 = C[hhh];
                        register double c21 = C[hhh + 1];
                        register double c22 = C[hhh + 2];

                        for (k1 = k; k1 < nk; k1 += 3) {
                            for (int l = 0; l < 3; l++) {
                                int ha = i1 * n + k1 + l;
                                int hha = ha + n;
                                int hhha = hha + n;
                                int hb = k1 * n + j1 + l * n;
                                register double a0 = A[ha];
                                register double a1 = A[hha];
                                register double a2 = A[hhha];
                                register double b0 = B[hb];
                                register double b1 = B[hb + 1];
                                register double b2 = B[hb + 2];

                                c00 += a0 * b0;
                                c01 += a0 * b1;
                                c02 += a0 * b2;
                                c10 += a1 * b0;
                                c11 += a1 * b1;
                                c12 += a1 * b2;
                                c20 += a2 * b0;
                                c21 += a2 * b1;
                                c22 += a2 * b2;
                            }
                        }
                        C[h] = c00;
                        C[h + 1] = c01;
                        C[h + 2] = c02;
                        C[hh] = c10;
                        C[hh + 1] = c11;
                        C[hh + 2] = c12;
                        C[hhh] = c20;
                        C[hhh + 1] = c21;
                        C[hhh + 2] = c22;

                    }
                }
            }
        }
    }

}