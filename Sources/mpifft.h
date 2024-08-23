#ifndef __MPIFFT_H__
#define __MPIFFT_H__

/****************** Sequential functions ********************************/
void ufft(double *x, int n, int sign, double *w);
void ufft_init(int n, double *w);
void twiddle(double *x, int n, int sign, double *w);
void twiddle_init(int n, double alpha, int *rho, double  *w);
void permute(double *x, int n, int *sigma);
void bitrev_init(int n, int *rho);

/****************** Parallel functions ********************************/
int k1_init(int n, int p);
void mpiredistr(double *x, int n, int p, int s, int c0, int c1,
                char rev, int *rho_p);
void mpifft(double *x, int n, int p, int s, int sign, double *w0, double *w,
            double *tw, int *rho_np, int *rho_p);
void mpifft_init(int n, int p, int s, double *w0, double *w, double *tw,
                 int *rho_np, int *rho_p);
void MPI_FFT(double *x, unsigned int n,unsigned int s, unsigned int p,int dir);


#endif 

