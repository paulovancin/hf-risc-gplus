#include "matrix.h"

//----------------------------->  SET MATRIX  <---------------------------------

struct Matrix set_values(int l, int c, typ_var s[])
{
	int k = 0;
  struct Matrix out;
  out.line = l;
  out.column = c;

	for(int i = 0; i < out.line; i++){
		for(int j = 0; j < out.column; j++){
			out.str[i][j] = s[k];
			k++;
		}
	}

  return out;
}


//--------------------------->  PRINT MATRIX  <---------------------------------

void print_matrix(struct Matrix m)
{
    	for(int i = 0; i < m.line; i++){
				for(int j = 0; j < m.column; j++){
					SHOW(m.str[i][j]);
				}
				printf("\n");
			}
			printf("\n");
}

//----------------------------->  GET LINE  <-----------------------------------

int get_line(struct Matrix m)
{
	return m.line;
}

//---------------------------->  GET COLUMN  <----------------------------------

int get_column(struct Matrix m)
{
	return m.column;
}

//------------------------->  GET A SINGLE VALUE  <-----------------------------

typ_var get_value(struct Matrix m, int lin, int col)
{
    typ_var VAL;
    VAL = m.str[lin][col];
		return VAL;
}


//-------------------------> TRANSPOSED MATRIX <--------------------------------

struct Matrix transposed(struct Matrix m)
{

  struct Matrix out;
  out.line = m.column;
  out.column = m.line;

	for (int i = 0; i < out.line; i++){
		for (int j = 0; j < out.column; j++){
			out.str[i][j] = m.str[j][i] ;
		}
	}

  return out;

}

//-------------------------> MATRIX SUM <---------------------------------------

struct Matrix sum(struct Matrix m1, struct Matrix m2)
{

  struct Matrix out;
  out.line = m1.line;
  out.column = m2.column;

	for (int i = 0; i < m1.line; i++){
		for (int j = 0; j < m1.column; j++){
			out.str[i][j] = m1.str[i][j] + m2.str[i][j];
		}
	}

  return out;

}


 //-------------------------> MATRIX SUBTRACTION <------------------------------

struct Matrix subtraction(struct Matrix m1, struct Matrix m2)
{

  struct Matrix out;
  out.line = m1.line;
  out.column = m1.column;

	for (int i = 0; i < out.line; i++){
		for (int j = 0; j < out.column; j++){
			out.str[i][j] = m1.str[i][j] - m2.str[i][j];
		}
	}

  return out;

}

//-----------------------> MATRIX MULTIPLICATION <------------------------------

struct Matrix multiplication(struct Matrix m1, struct Matrix m2)
{

  struct Matrix out;
  out.line = m1.line;
  out.column = m2.column;
  typ_var auxm = val(0.0);

	if (m1.column == m2.line){
		for (int i = 0; i < m1.line; i++){
			for (int j = 0; j < m2.column; j++){
				for (int x = 0; x < m2.line; x++){
					auxm = auxm + mul(m1.str[i][x], m2.str[x][j]);
				}
				out.str[i][j] = auxm;
				auxm = 0;
			}
		}

    return out;
	}

		else
		{
			printf("Impossible to Multiplicate Those Matrices!!!\n");
		}

}

//----------------> MULTIPLICATION NUMBER = CTE * MATRIX  <---------------------

struct Matrix multE(struct Matrix m, typ_var cte)
{
  struct Matrix out;
  out.line = m.line;
  out.column = m.column;

	for (int i = 0; i < out.line; i++){
		for (int j = 0; j < out.column; j++){
			out.str[i][j] = mul(m.str[i][j], cte);
		}
	}

  return out;

}


//--------------------> DIVISION = MATRIX / CTE  <------------------------------

struct Matrix divE(struct Matrix m, typ_var cte)
{
  struct Matrix out;
  out.line = m.line;
  out.column = m.column;

	for (int i = 0; i < out.line; i++){
		for (int j = 0; j < out.column; j++){
			out.str[i][j] = div(m.str[i][j], cte);
		}
	}

  return out;

}


//---------------------------> SET IDENTITY MATRIX <----------------------------

struct Matrix setEye(int order)
{

  struct Matrix out;
  out.line = order;
  out.column = order;

	for (int i = 0; i < order; i++){
		for (int j = 0; j < order; j++){
			if (j == i){
				out.str[i][j] = val(1);
			}
			else{
				out.str[i][j] = val(0);
			}
		}
	}

    return out;

}

//---------------------------> MATRIX OF ZEROS <--------------------------------

struct Matrix zeros(int lin, int col)
{
  struct Matrix out;
  out.line = lin;
  out.column = col;

	for(int i = 0; i < lin; i++){
		for(int j = 0; j < col; j++){
			out.str[i][j] = val(0);
		}
	}

  return out;
}


//---------------------------> MATRIX OF ONES <--------------------------------

struct Matrix ones(int lin, int col)
{
  struct Matrix out;
  out.line = lin;
  out.column = col;

	for(int i = 0; i < out.line; i++){
		for(int j = 0; j < out.column; j++){
			out.str[i][j] = val(1);
		}
	}

  return out;
}


//------------------> MATRIX - CROSS PRODUCT 3X3 <------------------------------

struct Matrix crossProduct3(typ_var x, typ_var y, typ_var z)
{
  struct Matrix out;
  out.line = 3;
  out.column = 3;
  out.str[0][0] = val(0);
  out.str[0][1] = mul(val(-1),z);
  out.str[0][2] = y;
  out.str[1][0] = z;
  out.str[1][1] = val(0);
  out.str[1][2] = mul(val(-1),x);
  out.str[2][0] = mul(val(-1),y);
  out.str[2][1] = x;
  out.str[2][2] = val(0);

  return out;
}

//-------------------> BLOCK DIAGONAL CONCATENATION <---------------------------


//                               |A 0 0 .. 0|
// blkdiagN(A,B,C,...) produces  |0 B 0 .. 0|, CAUTION: Matrices must be square!
//                               |0 0 C .. 0|

//-------------> BLOCK DIAGONAL CONCATENATION: INPUT 2 MATRICES <---------------

struct Matrix blkdiag2(struct Matrix m1, struct Matrix m2)
{
  struct Matrix aux;
  struct Matrix out;
  out.line = m1.line + m2.line;
  out.column = m1.column + m2.column;

  aux = zeros(out.line, out.column);

  for(int i = 0; i < m1.line; i++){
    for(int j = 0; j < m1.column; j++){
      aux.str[i][j] = m1.str[i][j];
    }
  }

  for(int i = 0; i < m2.line; i++){
    for(int j = 0; j < m2.column; j++){
      aux.str[i+m1.line][j+m1.column] = m2.str[i][j];
    }
  }

  for(int i = 0; i < out.line; i++){
    for(int j = 0; j < out.column; j++){
      out.str[i][j] = aux.str[i][j];
    }
  }

 return out;

}

//-------------------> CUSTOM MATRIX OF MATRICES <------------------------------
//
//                    customMat(n_mats, vet_conf, A,B,C,D..)
//                where, vet_conf = [confX, confY, line, column]
//                produces:
//
//    |A B .. |
//    |C D .. |, confY and confX indicate how the matrix will be build. In case
//    |   ..  |  of confX = 2 and confY = 2, the matrices A, B, C and D will be
//               placed in a 2x2 grid. The parameters line and column will
//               indicate the total size of the matrix. Be aware of how the
//               matrices will be placed due to the differces of sizes, for now
//               will not be developed a restriction for placement.


struct Matrix customMat(int n_mats, int conf[], ...)
{
  struct Matrix aux;
  struct Matrix out;

  int auxX = 0;
  int auxY = 0;

  out.line = conf[2];
  out.column = conf[3];


  va_list mats;
  va_start(mats, n_mats);


    for (int x = 0; x < conf[0]; x++){
      for (int y = 0; y < conf[1]; y++){
        aux = va_arg(mats, struct Matrix);
          for (int i = 0; i < aux.line; i++){
            for (int j = 0; j < aux.column; j++){
              out.str[i + auxX][j + auxY] = aux.str[i][j];
            }
          }

          auxY = auxY + aux.column;
          if (auxY >= out.column){
          auxY = 0;
          }
        }

        auxX = auxX + aux.line;
        if (auxX >= out.line){
        auxX = 0;
        }
      }



  va_end(mats);
  return out;

}

//-------------------> COPY OF A MATRIX <------------------------------

struct Matrix copy(struct Matrix m)
{
  struct Matrix out;
  out.line = m.line;
  out.column = m.column;

	for (int i = 0; i < out.line; i++){
		for (int j = 0; j < out.column; j++){
			out.str[i][j] = m.str[i][j] ;
		}
	}
  return out;
}

//----------------> INVERSE FUNCTION - ONLY 4X4 MATRICES <----------------------

struct Matrix Invert(struct Matrix In)
{
	  typ_var inv[16], m[16], invOut[16];
    typ_var det;
    int i;
		int count = 0;
    struct Matrix out;
    out.line = In.line;
    out.column = In.column;


		for(int x = 0; x < 4; x++){
			for(int y = 0; y < 4; y++){
				m[count] = In.str[x][y];
				count++;
			}
		}

   inv[0] = mul(mul(m[5],m[10]),m[15]) -
            mul(mul(m[5],m[11]),m[14]) -
            mul(mul(m[9],m[6]),m[15]) +
            mul(mul(m[9],m[7]),m[14]) +
            mul(mul(m[13],m[6]),m[11]) -
            mul(mul(m[13],m[7]),m[10]);


   inv[4] = mul(mul(mul(m[4],val(-1)), m[10]),m[15]) +
            mul(mul(m[4],m[11]),m[14]) +
            mul(mul(m[8],m[6]),m[15]) -
            mul(mul(m[8],m[7]),m[14]) -
            mul(mul(m[12],m[6]),m[11]) +
            mul(mul(m[12],m[7]),m[10]);

   inv[8] = mul(mul(m[4],m[9]),m[15]) -
            mul(mul(m[4],m[11]),m[13]) -
            mul(mul(m[8],m[5]),m[15]) +
            mul(mul(m[8],m[7]),m[13]) +
            mul(mul(m[12],m[5]),m[11]) -
            mul(mul(m[12],m[7]),m[9]);

   inv[12] = mul(mul(mul(m[4],val(-1)), m[9]),m[14]) +
             mul(mul(m[4],m[10]),m[13]) +
             mul(mul(m[8],m[5]),m[14]) -
             mul(mul(m[8],m[6]),m[13]) -
             mul(mul(m[12],m[5]),m[11]) +
             mul(mul(m[12],m[6]),m[9]);

   inv[1] = mul(mul(mul(m[1],val(-1)), m[10]),m[15]) +
            mul(mul(m[1],m[11]),m[14]) +
            mul(mul(m[9],m[2]),m[15]) -
            mul(mul(m[9],m[3]),m[14]) -
            mul(mul(m[13],m[2]),m[11]) +
            mul(mul(m[13],m[3]),m[10]);

   inv[5] = mul(mul(m[0],m[10]),m[15]) -
            mul(mul(m[0],m[11]),m[14]) -
            mul(mul(m[8],m[2]),m[15]) +
            mul(mul(m[8],m[3]),m[14]) +
            mul(mul(m[12],m[2]),m[11]) -
            mul(mul(m[12],m[3]),m[10]);

   inv[9] = mul(mul(mul(m[0],val(-1)), m[9]),m[15]) +
            mul(mul(m[0],m[11]),m[13]) +
            mul(mul(m[8],m[1]),m[15]) -
            mul(mul(m[8],m[3]),m[13]) -
            mul(mul(m[12],m[1]),m[11]) +
            mul(mul(m[12],m[3]),m[9]);

   inv[13] = mul(mul(m[0],m[9]),m[14]) -
             mul(mul(m[0],m[10]),m[13]) -
             mul(mul(m[8],m[1]),m[14]) +
             mul(mul(m[8],m[2]),m[13]) +
             mul(mul(m[12],m[1]),m[10]) -
             mul(mul(m[12],m[2]),m[9]);

   inv[2] = mul(mul(m[1],m[6]),m[15]) -
            mul(mul(m[1],m[7]),m[14]) -
            mul(mul(m[5],m[2]),m[15]) +
            mul(mul(m[5],m[3]),m[14]) +
            mul(mul(m[13],m[2]),m[7]) -
            mul(mul(m[13],m[3]),m[6]);

   inv[6] = mul(mul(mul(m[0],val(-1)), m[6]),m[15]) +
            mul(mul(m[0],m[7]),m[14]) +
            mul(mul(m[4],m[2]),m[15]) -
            mul(mul(m[4],m[3]),m[14]) -
            mul(mul(m[12],m[2]),m[7]) +
            mul(mul(m[12],m[3]),m[6]);

  inv[10] = mul(mul(m[0],m[5]),m[15]) -
            mul(mul(m[0],m[7]),m[13]) -
            mul(mul(m[4],m[1]),m[15]) +
            mul(mul(m[4],m[3]),m[13]) +
            mul(mul(m[12],m[1]),m[7]) -
            mul(mul(m[12],m[3]),m[5]);

  inv[14] = mul(mul(mul(m[0],val(-1)), m[5]),m[14]) +
            mul(mul(m[0],m[6]),m[13]) +
            mul(mul(m[4],m[1]),m[14]) -
            mul(mul(m[4],m[2]),m[13]) -
            mul(mul(m[12],m[1]),m[6]) +
            mul(mul(m[12],m[2]),m[5]);

  inv[3] = mul(mul(mul(m[1],val(-1)), m[6]),m[11]) +
           mul(mul(m[1],m[7]),m[10]) +
           mul(mul(m[5],m[2]),m[11]) -
           mul(mul(m[5],m[3]),m[10]) -
           mul(mul(m[9],m[2]),m[7]) +
           mul(mul(m[9],m[3]),m[6]);

  inv[7] = mul(mul(m[0],m[6]),m[11]) -
           mul(mul(m[0],m[7]),m[10]) -
           mul(mul(m[4],m[2]),m[11]) +
           mul(mul(m[4],m[3]),m[10]) +
           mul(mul(m[8],m[2]),m[7]) -
           mul(mul(m[8],m[3]),m[6]);

  inv[11] = mul(mul(mul(m[0],val(-1)), m[5]),m[11]) +
            mul(mul(m[0],m[7]),m[9]) +
            mul(mul(m[4],m[1]),m[11]) -
            mul(mul(m[4],m[3]),m[9]) -
            mul(mul(m[8],m[1]),m[7]) +
            mul(mul(m[8],m[3]),m[5]);

  inv[15] = mul(mul(m[0],m[5]),m[10]) -
            mul(mul(m[0],m[6]),m[9]) -
            mul(mul(m[4],m[1]),m[10]) +
            mul(mul(m[4],m[2]),m[9]) +
            mul(mul(m[8],m[1]),m[6]) -
            mul(mul(m[8],m[2]),m[5]);


    det = mul(m[0], inv[0]) + mul(m[1], inv[4]) + mul(m[2], inv[8]) + mul(m[3], inv[12]);
    det = div(val(1), det);



    for (i = 0; i < 16; i++){
        invOut[i] = mul(inv[i], det);
    }

    int k = 0;

    for(int i = 0; i < out.line; i++){
      for(int j = 0; j < out.column; j++){
        out.str[i][j] = invOut[k];
        k++;
      }
    }

      return out;
}
