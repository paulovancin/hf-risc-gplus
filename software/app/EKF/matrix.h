#include "setType.h"
#define NMAX 30


struct Matrix {
  int line;
  int column;
  typ_var str[NMAX][NMAX];
};

struct ekf_out {
  struct Matrix states;
  struct Matrix MatP;
};


struct Matrix set_values(int l, int c, typ_var s[]);
void print_matrix(struct Matrix m);
int get_line(struct Matrix m);
int get_column(struct Matrix m);
typ_var get_value(struct Matrix m, int lin, int col);
struct Matrix transposed(struct Matrix m);
struct Matrix sum(struct Matrix m1, struct Matrix m2);
struct Matrix subtraction(struct Matrix m1, struct Matrix m2);
struct Matrix multiplication(struct Matrix m1, struct Matrix m2);
struct Matrix multE(struct Matrix m, typ_var cte);
struct Matrix divE(struct Matrix m, typ_var cte);
struct Matrix setEye(int order);
struct Matrix zeros(int lin, int col);
struct Matrix ones(int lin, int col);
struct Matrix crossProduct3(typ_var x, typ_var y, typ_var z);
struct Matrix blkdiag2(struct Matrix m1, struct Matrix m2);
struct Matrix customMat(int n_mats, int conf[], ...);
struct Matrix copy(struct Matrix m);
struct Matrix Invert(struct Matrix In);
