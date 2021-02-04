#include "setType.h"
#include "Matrix.h"
#include "quaternion.h"


int main(int argc, char *argv[]){

  unsigned int mem_allocmat1[NEW_SIZE];
  unsigned int mem_allocmat2[NEW_SIZE];
  unsigned int mem_allocmat3[NEW_SIZE];
  unsigned int mem_allocmat4[NEW_SIZE];
  unsigned int mem_allocmat5[NEW_SIZE];
  unsigned int mem_allocinv[NEW_SIZE];
  unsigned int mem_allocQ[NEW_SIZE];
  unsigned int mem_allocR[NEW_SIZE];
  unsigned int mem_allocC[NEW_SIZE];
  unsigned int mem_allocQ2[NEW_SIZE];
  unsigned int mem_allocH[NEW_SIZE];
  unsigned int mem_allocEIGEN[NEW_SIZE];
  unsigned int mem_allocEIGEN2[NEW_SIZE];
  unsigned int mem_allocEIGEN3[NEW_SIZE];
  unsigned int mem_allocTu[NEW_SIZE];
  unsigned int mem_allocTl[NEW_SIZE];
  unsigned int mem_allocU[NEW_SIZE];
  unsigned int mem_allocS[NEW_SIZE];
  unsigned int mem_allocA[NEW_SIZE];
  unsigned int mem_allocB[NEW_SIZE];
  unsigned int mem_allocQa[NEW_SIZE];
  unsigned int mem_allocRa[NEW_SIZE];
  unsigned int mem_allocX[NEW_SIZE];


  Matrix* mat1 = new (mem_allocmat1) Matrix;
  Matrix* mat2 = new (mem_allocmat2) Matrix;
  Matrix* mat3 = new (mem_allocmat3) Matrix;
  Matrix* mat4 = new (mem_allocmat4) Matrix;
  Matrix* mat5 = new (mem_allocmat5) Matrix;
  Matrix* inv = new (mem_allocinv) Matrix;
  Matrix* Q = new (mem_allocQ) Matrix;
  Matrix* R = new (mem_allocR) Matrix;
  Matrix* C = new (mem_allocC) Matrix;
  Matrix* Q2 = new (mem_allocQ2) Matrix;
  Matrix* H = new (mem_allocH) Matrix;
  Matrix* EIGEN = new (mem_allocEIGEN) Matrix;
  Matrix* EIGEN2 = new (mem_allocEIGEN2) Matrix;
  Matrix* EIGEN3 = new (mem_allocEIGEN3) Matrix;
  Matrix* Tu = new (mem_allocTu) Matrix;
  Matrix* Tl = new (mem_allocTl) Matrix;
  Matrix* U = new (mem_allocU) Matrix;
  Matrix* S = new (mem_allocS) Matrix;
  Matrix* A = new (mem_allocA) Matrix;
  Matrix* B = new (mem_allocB) Matrix;
  Matrix* Qa = new (mem_allocQa) Matrix;
  Matrix* Ra = new (mem_allocRa) Matrix;
  Matrix* X = new (mem_allocX) Matrix;


  Complex val1[16] = {val(1),val(2),val(3),val(4),val(5),val(6),val(7),val(8),val(9),val(10),val(11),val(12),val(13),val(14),val(15),val(16)};
  Complex val2[16] = {val(2),val(3),val(4),val(5),val(6),val(7),val(8),val(9),val(10),val(11),val(12),val(13),val(14),val(15),val(16),val(17)};
  Complex val3[16] = {val(19),val(-12),val(-14),val(8),val(17),val(-10),val(-14),val(8),val(12),val(-9),val(-9),val(7),val(13),val(-10),val(-12),val(10)};
  Complex val4[16] = {val(-1.2075),val(1.0347),val(-0.7873),val(-0.8095),val(0.7172),val(0.7269),val(0.8884),val(-2.9443),val(1.6302),val(-0.3034),val(-1.1471),val(1.4384),val(0.4889),val(0.2939),val(-1.0689),val(0.3252)};
  Complex val5[9] = {val(3),val(2),val(1),val(4),val(2),val(1),val(4),val(4),val(0)};
  Complex valA[36] = {val(0.987),val(1),val(1.0094),val(1),val(0),val(0),val(1),val(0.9884),val(1),val(1.0076),val(0),val(0),val(1),val(1),val(0.9907),val(1),
  val(0),val(0),val(1),val(1),val(1),val(0.9925),val(0),val(0),val(-1),val(0),val(0),val(0),val(0),val(0),val(0),val(-1),val(0),val(0),val(0),val(0)};
  Complex valB[12] = {val(-0.0547),val(-0.0452),val(-0.0429),val(-0.0571),val(-0.0644),val(-0.0354),val(-0.0309),val(-0.0691),val(0),val(0),val(0),val(0)};
  Complex valR[4] = {val(0.01),val(0),val(0),val(0.01)};
  Complex valQ[36] = {val(20),val(0),val(0),val(0),val(0),val(0),val(0),val(20),val(0),val(0),val(0),val(0),val(0),val(0),val(20),val(0),
  val(0),val(0),val(0),val(0),val(0),val(20),val(0),val(0),val(0),val(0),val(0),val(0),val(20),val(0),val(0),val(0),val(0),val(0),val(0),val(20)};
  mat1->set_values(4, 4, val1);
  mat2->set_values(4, 4, val2);
  mat3->set_values(4, 4, val3);
  mat4->set_values(4, 4, val4);
  mat5->set_values(3, 3, val5);
  A->set_values(6, 6, valA);
  B->set_values(6, 2, valB);
  Ra->set_values(2, 2, valR);
  Qa->set_values(6, 6, valQ);
  Complex det, tr, nor;
  Complex  x, y, z;


//--------------------------FUNCTIONS TESTS (MATRIX)----------------------------

  printf("Matrix 1 is:\n");
  mat1->print_matrix();
  printf("Matrix 2 is:\n");
  mat2->print_matrix();
  printf("Matrix 3 is:\n");
  mat3->print_matrix();
  printf("Matrix 5 is:\n");
  mat5->print_matrix();
  det = mat3->determinant(4);
  printf("Determinant of Matrix 3 is: ");
  show(det);
  // printf("\n");
  // x = val(2.5);
  // y = val(2);
  // z = square_root(x);
  // show(z);
  // printf("\n");
  mat3->inverse(inv);
  printf("Inverse of Matrix 3 is:\n");
  inv->print_matrix();
  mat3->qr(Q,R);
  printf("Q of Matrix 3 is:\n");
  Q->print_matrix();
  printf("R of Matrix 3 is:\n");
  R->print_matrix();
  // printf("A Chunk of the Matrix 3 is:\n");
  // C->get_part(0, 3, 3, 3, mat3);
  // C->print_matrix();
  // mat3->householder(H,Q2);
  // printf("Householder H of Matrix 3 is:\n");
  // H->print_matrix();
  // printf("Householder Q of Matrix 3 is:\n");
  // Q2->print_matrix();
  // printf("Eigenvalues of Matrix 3 are:\n");
  // EIGEN->eig(mat3);
  // EIGEN->print_matrix();
  // printf("Matrix 1 is:\n");
  // mat1->print_matrix();
  // printf("Zeros Below the Diagonal of Matrix 1:\n");
  // Tu->triu(mat1);
  // Tu->print_matrix();
  // printf("Zeros Above the Diagonal of Matrix 1:\n");
  // Tl->tril(mat1);
  // Tl->print_matrix();
  // tr = mat1->trace();
  // printf("Trace of Matrix 1 is: ");
  // show(tr);
  // printf("\n");
  // printf("Eigenvalues of Matrix 1 are:\n");
  // EIGEN2->eig(mat1);
  // EIGEN2->print_matrix();
  // printf("Matrix 4 is:\n");
  // mat4->print_matrix();
  // printf("Eigenvalues of Matrix 4 are:\n");
  // EIGEN3->eig(mat4);
  // EIGEN3->print_matrix();
  // nor = mat1->twoNorm();
  // printf("2-Norm of Matrix 1 is: ");
  // show(nor);
  // printf("\n");
  // printf("\n");
  // printf("\n");
  // printf("----------RICCATI EQUATION SOLVER SECTION------------");
  // printf("\n");
  // printf("Matrix A is:\n");
  // A->print_matrix();
  // printf("Matrix B is:\n");
  // B->print_matrix();
  // printf("Matrix R is:\n");
  // Ra->print_matrix();
  // printf("Matrix Q is:\n");
  // Qa->print_matrix();
  // printf("Matrix X is the Riccati Solution:\n");
  // X->riccati(A,B,Ra,Qa);
  // X->print_matrix();

}
