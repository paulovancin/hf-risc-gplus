#include "setType.h"
#include "matrix.h"



int main(int argc, char* argv[]) {


  struct Matrix M1, M2, M3, M4, M5, M6, M7;
  typ_var Multiplicator = val(5);

        typ_var val1[16] = {val(1),val(2),val(3),val(4),val(5),val(6),val(7),val(8),val(9),val(10),val(11),val(12),val(13),val(14),val(15),val(16)};
        typ_var val2[16] = {val(2),val(3),val(4),val(5),val(6),val(7),val(8),val(9),val(10),val(11),val(12),val(13),val(14),val(15),val(16),val(17)};
        typ_var val3[16] = {val(4),val(2),val(2),val(2),val(0),val(1),val(-3),val(3),val(0),val(-1),val(3),val(3),val(0),val(3),val(1),val(1)};

        M1 = set_values(4, 4, val1);
        M2 = set_values(4, 4, val2);
        M3 = set_values(4, 4, val3);

//----------------------------FUNCTIONS TESTS-----------------------------------


        printf("Matrix 1 is:\n");
        print_matrix(M1);
        printf("Matrix 2 is:\n");
        print_matrix(M2);
        printf("Matrix 3 is:\n");
        print_matrix(M3);
        printf("Transposed of Matrix 1 is:\n");
        M4 = transposed(M1);
        print_matrix(M4);
        printf("Sum of Matrix 1 and 2 is:\n");
        M5 = sum(M1, M2);
        print_matrix(M5);
        printf("Subtraction of Matrix 1 and 2 is:\n");
        M5 = subtraction(M1, M2);
        print_matrix(M5);
        printf("Multiplication of Matrix 1 and 2 is:\n");
        M5 = multiplication(M1, M2);
        print_matrix(M5);
        printf("Multiplication of Matrix 1 by ");
        SHOW(Multiplicator);
        printf(":\n");
        M5 = multE(M1, Multiplicator);
        print_matrix(M5);
        printf("Division of Matrix 1 by ");
        SHOW(Multiplicator);
        printf(":\n");
        M5 = divE(M1, Multiplicator);
        print_matrix(M5);
        printf("Inverse of Matrix 3 is:\n");
        M5 = Invert(M3);
        print_matrix(M5);
        printf("Matrix of zeros of size %d by %d:\n", 4, 4);
        M5 = zeros(4, 4);
        print_matrix(M5);
        printf("Identity Matrix of size %d by %d:\n", 4, 4);
        M5 = setEye(4);
        print_matrix(M5);
        printf("Matrix of ones of size %d by %d:\n", 3, 3);
        M6 = ones(3, 3);
        print_matrix(M6);
        printf("Diagonal Block Concatenation of Two Matrices:\n");
        M7 = blkdiag2(M1, M2);
        print_matrix(M7);
        int conf[4] = {2,2,8,8};
        M7 =  customMat(4, conf, M1, M2, M3, M5);
        printf("Custom Matrix Made of Four Matrices:\n");
        print_matrix(M7);




}
