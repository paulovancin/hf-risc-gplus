#ifndef __MATRIX_H
#define __MATRIX_H
#define NMAX 25

#include "setType.h"
#include "libaux.h"

class Matrix {
  private:
		int row;
		int column;
		Complex str[NMAX][NMAX];
	public:
		void set_values(int r, int c, Complex s[]);
		void print_matrix();
		int get_row();
		int get_column();
		Complex get_value(int rw, int col);
    void transposed(Matrix *t);
    void sum(Matrix *m1, Matrix *m2);
    void subtraction(Matrix *m1, Matrix *m2);
    void multiplication(Matrix *m1, Matrix *m2);
    void multE(Matrix *m, Complex cte);
    void divE(Matrix *m, Complex cte);
    void sumE(Matrix *m, Complex cte);
    Complex determinant(int order);
    void cofactor(Matrix *out, int rw, int col, int order);
    Complex auxDet(int order);
    void inverse(Matrix *inv);
    void setEye(int order);
    void zeros(int rw, int col);
    void ones(int rw, int col);
    void crossProduct3(Complex x, Complex y, Complex z);
    void blkdiag2(Matrix *m1, Matrix *m2);
    void blkdiag3(Matrix *m1, Matrix *m2, Matrix *m3);
    void blkdiag4(Matrix *m1, Matrix *m2, Matrix *m3, Matrix *m4);
    void customMat(int n_mats, int conf[], ...);
    void copy(Matrix *o);
    int rank();
    int checkSymmetry();
    void luDecomposition(Matrix *L, Matrix *U);
    void qr(Matrix *L, Matrix *U);
    void diag(Matrix *out);
    void get_part(int ri, int rf, int ci, int cf, Matrix *in);
    void householder(Matrix *H, Matrix *Q);
    void eig(Matrix *in);
    void triu(Matrix *in);
    void tril(Matrix *in);
    Complex trace();
    void eig22(Matrix *in);
    Complex twoNorm();
    void riccati(Matrix *A, Matrix *B, Matrix *R, Matrix *Q);

};


#endif
