#include "Matrix.h"



//----------------------------->  SET MATRIX  <---------------------------------

void Matrix::set_values(int r, int c, Complex s[])
{
	row = r;
	column = c;
	int k = 0;

	for(int i = 0; i < row; i++){
		for(int j = 0; j < column; j++){
			str[i][j] = s[k];
			k++;
		}
	}
}


//--------------------------->  PRINT MATRIX  <---------------------------------

void Matrix::print_matrix()
{
    	for(int i = 0; i < row; i++){
				for(int j = 0; j < column; j++){
					show(str[i][j]);
					printf("     ");
				}
				printf("\n");
			}
			printf("\n");
}

//----------------------------->  GET ROW  <-----------------------------------

int Matrix::get_row()
{
	return row;
}

//---------------------------->  GET COLUMN  <----------------------------------

int Matrix::get_column()
{
	return column;
}

//------------------------->  GET A SINGLE VALUE  <-----------------------------

Complex Matrix::get_value(int rw, int col)
{
    Complex VAL;
    VAL = str[rw][col];
		return VAL;
}


//-------------------------> TRANSPOSED MATRIX <--------------------------------

void Matrix::transposed(Matrix *t)
{
  column = t->row;
  row = t->column;

	for (int i = 0; i < row; i++){
		for (int j = 0; j < column; j++){
			str[i][j] = t->str[j][i];
		}
	}

}

//-------------------------> MATRIX SUM <---------------------------------------

void Matrix::sum(Matrix *m1, Matrix *m2)
{
	row = m1->row;
	column = m1->column;

	for (int i = 0; i < m1->row; i++){
		for (int j = 0; j< m1->column; j++){
			str[i][j] = m1->str[i][j] + m2->str[i][j];
		}
	}
}


 //-------------------------> MATRIX SUBTRACTION <------------------------------

void Matrix::subtraction(Matrix *m1, Matrix *m2)
{
row = m1->row;
column = m1->column;

	for (int i = 0; i < m1->row; i++){
		for (int j = 0; j< m1->column; j++){
			str[i][j] = m1->str[i][j] - m2->str[i][j];
		}
	}

}

//-----------------------> MATRIX MULTIPLICATION <------------------------------

void Matrix::multiplication(Matrix *m1, Matrix *m2)
{

		Complex auxm = 0.0;
		row = m1->row;
		column = m2->column;


		if (m1->column == m2->row){
			for (int i = 0; i < m1->row; i++){
				for (int j = 0; j < m2->column; j++){
					for (int x = 0; x < m2->row; x++){
						auxm = auxm + mul(m1->str[i][x], m2->str[x][j]);
					}
					str[i][j] = auxm;
					auxm = 0;
				}
			}
		}

		else
		{
			printf("Impossible to Multiplicate Those Matrices!!!\n");
		}

}

//----------------> MULTIPLICATION NUMBER = CTE * MATRIX  <---------------------

void Matrix::multE(Matrix *m, Complex cte)
{
  row = 	m->row;
  column = m->column;

	for (int i = 0; i < row; i++){
		for (int j = 0; j < column; j++){
			str[i][j] = mul(m->str[i][j], cte);
		}
	}

}


//--------------------> DIVISION = MATRIX / CTE  <------------------------------

void Matrix::divE(Matrix *m, Complex cte)
{
  row = 	m->row;
  column = m->column;

	for (int i = 0; i < row; i++){
		for (int j = 0; j < column; j++){
			str[i][j] = division(m->str[i][j], cte);
		}
	}

}

//------------------> SUM NUMBER = CTE * MATRIX  <-----------------------------

void Matrix::sumE(Matrix *m, Complex cte)
{
  row = 	m->row;
  column = m->column;

	for (int i = 0; i < row; i++){
		for (int j = 0; j < column; j++){
			str[i][j] = m->str[i][j] + cte;
		}
	}

}


//-------------------------> MATRIX DETERMINANT  <------------------------------

Complex Matrix::determinant(int order)
{

  unsigned int mem_allocCoFactor[NEW_SIZE];
  Matrix* CoFactor = new (mem_allocCoFactor) Matrix;

 Complex D = val(0), D_aux;
 Complex sign = val(1);

	if (order == 1){
		  return (str[0][0]);
  }

  if (order > 1 && order <= 3){
    	D = auxDet(order);
  }

  if (order >= 4){
		for (int f = 0; f < order; f++){
			cofactor(CoFactor, 0, f, order);
			D_aux = mul(mul(sign, str[0][f]), CoFactor->determinant((order - 1)));
			D = D + D_aux;
			sign = mul(sign,val(-1));
		}
  }

 return D;

}

//--------------------------> MATRIX CO-FACTOR  <-------------------------------

void Matrix::cofactor(Matrix *out, int rw, int col, int order)
{

	int i = 0, j = 0;
	Complex temp[order][order];

	for (int g = 0; g < order; g++){
		for(int h = 0; h < order; h++){
			if(g != rw && h != col){
				temp[i][j++] = str[g][h];
					if(j == order - 1){
						j = 0;
						i++;
					}
			}
		}
	}

  out->row = order;
	out->column = order;

	for(int h = 0; h < order; h++){
		for(int g = 0; g < order; g++){
			out->str[h][g] = temp[h][g];
		}
	}

}

//------------------> AUXILIARY CALCULATION FOR DETERMINANT <-------------------

Complex Matrix::auxDet(int order)
{
	Complex D = 0, D_aux1, D_aux2, D_aux3, D_aux4, D_aux5, D_aux6;
	Complex D_aux7, D_aux8, D_aux9, D_aux10, D_aux11, D_aux12;

	switch (order) {

		case 1:
		D = str[0][0];
		return D;
		break;

		case 2:
		D = mul(str[0][0], str[1][1]) - mul(str[0][1], str[1][0]);
		return D;
		break;

		case 3:
		D_aux1 = (mul(mul(str[0][0],str[1][1]), str[2][2])) - (mul(mul(str[0][0],str[1][2]), str[2][1])) + (mul(mul(str[0][1],str[1][2]), str[2][0]));
		D_aux2 = mul(val(-1),(mul(mul(str[0][1],str[1][0]), str[2][2]))) + (mul(mul(str[0][2],str[1][0]), str[2][1])) - (mul(mul(str[0][2],str[1][1]), str[2][0]));
		D = D_aux1 + D_aux2;
		return D;
		break;

	}

}

//---------------------------> MATRIX INVERSE  <--------------------------------

void Matrix::inverse(Matrix *inv)
{
  unsigned int mem_allocb[NEW_SIZE];
  Matrix* b = new (mem_allocb) Matrix;

  unsigned int mem_allocfac[NEW_SIZE];
  Matrix* fac = new (mem_allocfac) Matrix;

 	int p,q,m,n,i,j;
	int order;
	Complex d;

	order = row;

	Complex Aux[order][order];
  inv->zeros(order, order);
	b->zeros(order, order);
	fac->zeros(order, order);

 	for (q = 0; q < order; q++){
   	for (p = 0; p < order; p++){
     	m = 0;
    	n = 0;
     	for (i = 0; i < order; i++){
       	for (j = 0; j < order; j++){
          if (i != q && j != p){
            b->str[n][m] = str[j][i];
            	if (n < (order-2)){
             		n++;
							}
           		else{
               n=0;
               m++;
              }
          }
        }
      }

			if (p == 0 && q == 0){
				fac->str[p][q] = b->determinant(order-1);
			}
			else {

				if ((p + q) % 2 == 0){
					fac->str[p][q] = mul(val(1.0), b->determinant(order-1));
				}
      	else {
					fac->str[p][q] = mul(val(-1.0), b->determinant(order-1));
				}
			}
     }
  	}

		for (int g = 0; g < order; g++){
			for (int h = 0; h < order; h++){
				Aux[h][g] = fac->str[g][h];
			}
		}

		d = determinant(order);


		for (int g = 0; g < order; g++){
			for (int h = 0; h < order; h++){
				inv->str[h][g] = division(Aux[h][g], d);
			}
		}

}



//---------------------------> SET IDENTITY MATRIX <----------------------------

void Matrix::setEye(int order)
{

	row = order;
	column = order;

	for (int i = 0; i < order; i++){
		for (int j = 0; j < order; j++){
			if (j == i){
				str[i][j] = val(1);
			}
			else{
				str[i][j] = val(0);
			}
		}
	}

}

//---------------------------> MATRIX OF ZEROS <--------------------------------

void Matrix::zeros(int rw, int col)
{
	row = rw;
	column = col;

	for(int i = 0; i < rw; i++){
		for(int j = 0; j < col; j++){
			str[i][j] = val(0);
		}
	}
}


//---------------------------> MATRIX OF ONES <--------------------------------

void Matrix::ones(int rw, int col)
{
	row = rw;
	column = col;

	for(int i = 0; i < rw; i++){
		for(int j = 0; j < col; j++){
			str[i][j] = val(1);
		}
	}
}


//------------------> MATRIX - CROSS PRODUCT 3X3 <------------------------------

void Matrix::crossProduct3(Complex x, Complex y, Complex z)
{
  row = 3;
  column = 3;
  str[0][0] = val(0);
  str[0][1] = mul(val(-1),z);
  str[0][2] = y;
  str[1][0] = z;
  str[1][1] = val(0);
  str[1][2] = mul(val(-1),x);
  str[2][0] = mul(val(-1),y);
  str[2][1] = x;
  str[2][2] = val(0);
}

//-------------------> BLOCK DIAGONAL CONCATENATION <---------------------------


//                               |A 0 0 .. 0|
// blkdiagN(A,B,C,...) produces  |0 B 0 .. 0|, CAUTION: Matrices must be square!
//                               |0 0 C .. 0|

//-------------> BLOCK DIAGONAL CONCATENATION: INPUT 2 MATRICES <---------------

void Matrix::blkdiag2(Matrix *m1, Matrix *m2)
{

  unsigned int mem_allocaux[NEW_SIZE];
  Matrix* aux = new (mem_allocaux) Matrix;

  row = m1->row + m2->row;
  column = m1->column + m2->column;

  aux->zeros(row, column);

  for(int i = 0; i < m1->row; i++){
    for(int j = 0; j < m1->column; j++){
      aux->str[i][j] = m1->str[i][j];
    }
  }

  for(int i = 0; i < m2->row; i++){
    for(int j = 0; j < m2->column; j++){
      aux->str[i+m1->row][j+m1->column] = m2->str[i][j];
    }
  }

  for(int i = 0; i < row; i++){
    for(int j = 0; j < column; j++){
      str[i][j] = aux->str[i][j];
    }
  }

}

//-------------> BLOCK DIAGONAL CONCATENATION: INPUT 3 MATRICES <---------------

void Matrix::blkdiag3(Matrix *m1, Matrix *m2, Matrix *m3)
{

  unsigned int mem_allocaux[NEW_SIZE];
  Matrix* aux = new (mem_allocaux) Matrix;

  row = m1->row + m2->row + m3->row;
  column = m1->column + m2->column + m3->column;

  aux->zeros(row, column);

  for(int i = 0; i < m1->row; i++){
    for(int j = 0; j < m1->column; j++){
      aux->str[i][j] = m1->str[i][j];
    }
  }

  for(int i = 0; i < m2->row; i++){
    for(int j = 0; j < m2->column; j++){
      aux->str[i+m1->row][j+m1->column] = m2->str[i][j];
    }
  }

  for(int i = 0; i < m3->row; i++){
    for(int j = 0; j < m3->column; j++){
      aux->str[i+m1->row+m2->row][j+m1->column+m2->column] = m3->str[i][j];
    }
  }

  for(int i = 0; i < row; i++){
    for(int j = 0; j < column; j++){
      str[i][j] = aux->str[i][j];
    }
  }

}

//-------------> BLOCK DIAGONAL CONCATENATION: INPUT 4 MATRICES <---------------

void Matrix::blkdiag4(Matrix *m1, Matrix *m2, Matrix *m3, Matrix *m4)
{

  unsigned int mem_allocaux[NEW_SIZE];
  Matrix* aux = new (mem_allocaux) Matrix;

  row = m1->row + m2->row + m3->row + m4->row;
  column = m1->column + m2->column + m3->column + m4->column;

  aux->zeros(row, column);

  for(int i = 0; i < m1->row; i++){
    for(int j = 0; j < m1->column; j++){
      aux->str[i][j] = m1->str[i][j];
    }
  }

  for(int i = 0; i < m2->row; i++){
    for(int j = 0; j < m2->column; j++){
      aux->str[i+m1->row][j+m1->column] = m2->str[i][j];
    }
  }

  for(int i = 0; i < m3->row; i++){
    for(int j = 0; j < m3->column; j++){
      aux->str[i+m1->row+m2->row][j+m1->column+m2->column] = m3->str[i][j];
    }
  }

  for(int i = 0; i < m4->row; i++){
    for(int j = 0; j < m4->column; j++){
      aux->str[i+m1->row+m2->row+m3->row][j+m1->column+m2->column+m3->column] = m4->str[i][j];
    }
  }

  for(int i = 0; i < row; i++){
    for(int j = 0; j < column; j++){
      str[i][j] = aux->str[i][j];
    }
  }

}


//-------------------> CUSTOM MATRIX OF MATRICES <------------------------------
//
//                    customMat(n_mats, vet_conf, A,B,C,D..)
//                where, vet_conf = [confX, confY, row, column]
//                produces:
//
//    |A B .. |
//    |C D .. |, confY and confX indicate how the matrix will be build. In case
//    |   ..  |  of confX = 2 and confY = 2, the matrices A, B, C and D will be
//               placed in a 2x2 grid. The parameters row and column will
//               indicate the total size of the matrix. Be aware of how the
//               matrices will be placed due to the differces of sizes, for now
//               will not be developed a restriction for placement.


void Matrix::customMat(int n_mats, int conf[], ...)
{

  unsigned int mem_allocaux[NEW_SIZE];
  Matrix* aux = new (mem_allocaux) Matrix;

  int auxX = 0;
  int auxY = 0;

  row = conf[2];
  column = conf[3];


  va_list mats;
  va_start(mats, n_mats);


    for (int x = 0; x < conf[0]; x++){
      for (int y = 0; y < conf[1]; y++){
        aux = va_arg(mats, Matrix*);
          for (int i = 0; i < aux->row; i++){
            for (int j = 0; j < aux->column; j++){
              str[i + auxX][j + auxY] = aux->str[i][j];
            }
          }

          auxY = auxY + aux->column;
          if (auxY >= column){
          auxY = 0;
          }
        }

        auxX = auxX + aux->row;
        if (auxX >= row){
        auxX = 0;
        }
      }



  va_end(mats);

}

//-------------------> COPY OF A MATRIX <------------------------------

void Matrix::copy(Matrix *o)
{
  row = o->row;
  column = o->column;

	for (int i = 0; i < row; i++){
		for (int j = 0; j < column; j++){
			str[i][j] = o->str[i][j] ;
		}
	}

}


//-------------------------> RANK OF A MATRIX <---------------------------------

int Matrix::rank()
{
  int rank = column;
  Complex mult, temp;

  for (int row = 0; row < rank; row++){
    if (str[row][row].re() != val(0)){
      for (int col = 0; col < row; col++){
        if (col != row){
          mult = division(str[col][row], str[row][row]);
          for (int i = 0; i < rank; i++){
            str[col][i] -= mul(mult, str[row][i]);
          }
        }
      }
    }
    else{
      bool reduce = true;
      for (int i = row + 1; i < row; i++){
        if (str[i][row].re() != val(0)){
          for (int j = 0; j < rank; j++){
            temp = str[row][j];
            str[row][j] = str[i][j];
            str[i][j] = temp;
          }
          reduce = false;
          break;
        }
      }
      if (reduce){
        rank--;
        for (int i = 0; i < row; i++){
          str[i][row] = str[i][rank];
        }
      }
      row--;
    }
  }

  return rank;
}

//---------------------> CHECK SYMMETRY OF A MATRIX <---------------------------

int Matrix::checkSymmetry()
{
  int sym;
  int sum = 0, check;

  unsigned int mem_alloccp[NEW_SIZE];
  Matrix* cp = new (mem_alloccp) Matrix;

  unsigned int mem_alloctr[NEW_SIZE];
  Matrix* tr = new (mem_alloctr) Matrix;

  cp->zeros(row, column);
  tr->zeros(row, column);

  if (row != column){
    sym = 0;
  }

  else{
    for(int i = 0; i < row; i++){
      for(int j = 0; j < column; j++){
        cp->str[i][j] = str[i][j];
      }
    }

    check = row * column;
    tr->transposed(cp);



    for(int i = 0; i < row; i++){
      for(int j = 0; j < column; j++){
        if (cp->str[i][j] == tr->str[i][j]){
          sum++;
        }
      }
    }

    if (check == sum){
      sym = 1;
    }

    else{
      sym = 0;
    }
  }

  if (sym == 0){
    printf("---->Input Matrix is NOT Symmetric!!!<----\n");
  }

  else{
    printf("---->Input Matrix IS Symmetric!!!<----\n");
  }

  return sym;
}

//--------> DECOMPOSE A MATRIX INTO LOWER AND UPPER TRINGULAR MATRICES <--------

void Matrix::luDecomposition(Matrix *L, Matrix *U)
{
  L->zeros(row, column);
  U->zeros(row, column);

   int i = 0, j = 0, k = 0;
   for (i = 0; i < row; i++) {
      for (j = 0; j < row; j++) {
         if (j < i)
         L->str[j][i] = val(0);
         else {
            L->str[j][i] = str[j][i];
            for (k = 0; k < i; k++) {
               L->str[j][i] = L->str[j][i] - mul(L->str[j][k], U->str[k][i]);
            }
         }
      }
      for (j = 0; j < row; j++) {
         if (j < i)
         U->str[i][j] = val(0);
         else if (j == i)
         U->str[i][j] = val(1);
         else {
            U->str[i][j] = division(str[i][j], L->str[i][i]);
            for (k = 0; k < i; k++) {
               U->str[i][j] = U->str[i][j] - division(mul(L->str[i][k], U->str[k][j]), L->str[i][i]);
            }
         }
      }
   }
}



//-------------> QR DECOMPOSITION OF A SQUARE MATRIX <--------------------------


void Matrix::qr(Matrix *Q, Matrix *R){

unsigned int mem_allocA[NEW_SIZE];
Matrix* A = new (mem_allocA) Matrix;

unsigned int mem_allocz[NEW_SIZE];
Matrix* z = new (mem_allocz) Matrix;

unsigned int mem_allocv[NEW_SIZE];
Matrix* v = new (mem_allocv) Matrix;

unsigned int mem_allocvt[NEW_SIZE];
Matrix* vt = new (mem_allocvt) Matrix;

unsigned int mem_allocP[NEW_SIZE];
Matrix* P = new (mem_allocP) Matrix;

unsigned int mem_allocP1[NEW_SIZE];
Matrix* P1 = new (mem_allocP1) Matrix;

unsigned int mem_allocP2[NEW_SIZE];
Matrix* P2 = new (mem_allocP2) Matrix;

unsigned int mem_allocP3[NEW_SIZE];
Matrix* P3 = new (mem_allocP3) Matrix;

unsigned int mem_allocP4[NEW_SIZE];
Matrix* P4 = new (mem_allocP4) Matrix;

unsigned int mem_allocI[NEW_SIZE];
Matrix* I = new (mem_allocI) Matrix;

unsigned int mem_allocRc[NEW_SIZE];
Matrix* Rc = new (mem_allocRc) Matrix;

unsigned int mem_allocRaux[NEW_SIZE];
Matrix* Raux = new (mem_allocRaux) Matrix;

unsigned int mem_allocQc[NEW_SIZE];
Matrix* Qc = new (mem_allocQc) Matrix;

unsigned int mem_allocQaux[NEW_SIZE];
Matrix* Qaux = new (mem_allocQaux) Matrix;

unsigned int mem_allocQt[NEW_SIZE];
Matrix* Qt = new (mem_allocQt) Matrix;

unsigned int mem_allocRt[NEW_SIZE];
Matrix* Rt = new (mem_allocRt) Matrix;


Complex v1, normz, P2inv;

A->row = row;
A->column = column;

  for (int h = 0; h < row; h++){
    for (int g = 0; g < column; g++){
      A->str[h][g] = str[h][g] ;
    }
  }

  Q->setEye(A->get_row());
  R->copy(A);

  for (int i = 0; i < A->get_column(); i++){
    z->get_part(i, (A->get_row()-1), i, i, R);
    Complex vectorz[z->get_row()];

    for (int j = 0; j < z->get_row(); j++){
      vectorz[j] = z->str[j][0];
    }

    normz = norm_vector(z->get_row(),vectorz);
    v1 = mul(mul(val(-1),sign(z->str[0][0])),normz) - z->str[0][0];

    v->row = z->get_row();
    v->column = z->get_column();

    for (int k = 0; k < z->get_row(); k++){
      if (k == 0){
        v->str[k][0] = v1;
      }
      else{
        v->str[k][0] = mul(val(-1),z->str[k][0]);
      }
    }


    vt->transposed(v);
    P1->multiplication(v,vt);
    P2->multiplication(vt,v);
    P2inv = P2->get_value(0,0);
    P3->divE(P1,P2inv);
    P4->multE(P3, val(2));
    I->setEye(A->get_row()-i+1);
    P->subtraction(I,P4);


    Rc->get_part(i,A->get_row(),0,A->get_column(),R);
    Raux->multiplication(P,Rc);
    Qc->get_part(i,A->get_row(),0,A->get_column(),Q);
    Qaux->multiplication(P,Qc);

    for(int g = i; g < A->get_row(); g++){
      for(int h = 0; h < A->get_column(); h++){
    //   R->str[g][h] = Raux->str[g-i][h]; // --> In eig, this is fucking up!¹
      }
    }
	//
    for(int g = i; g < A->get_row(); g++){
      for(int h = 0; h < A->get_column(); h++){
    //   Q->str[g][h] = Qaux->str[g-i][h]; //¹ --> Same thing!
      }
    }
  }
  Qt->transposed(Q);
  Q->copy(Qt);
  Rt->triu(R);
  R->copy(Rt);
}


//------------------------> GET MATRIX DIAGONAL <-------------------------------

void Matrix::diag(Matrix *out){

  out->zeros(row, 1);
  out->row = row;

  for(int i = 0; i < row; i++){
    out->str[i][0] = str[i][i];
  }

}


//------------------------> GET PART OF A MATRIX <------------------------------

void Matrix::get_part(int ri, int rf, int ci, int cf, Matrix *in){

  row = rf - ri + 1;
  column = cf - ci + 1;

  for (int i = 0; i < row; i++){
    for (int j = 0; j < column; j++){
      str[i][j] = in->str[i+ri][j+ci] ;
    }
  }

}

//----------------------> HOUSEHOLDER DECOMPOSITION <---------------------------


void Matrix::householder(Matrix *H, Matrix *Q){

  unsigned int mem_allocx[NEW_SIZE];
  Matrix* x = new (mem_allocx) Matrix;

  unsigned int mem_allocA[NEW_SIZE];
  Matrix* A = new (mem_allocA) Matrix;

  unsigned int mem_allocA2[NEW_SIZE];
  Matrix* A2 = new (mem_allocA2) Matrix;

  unsigned int mem_allocA3[NEW_SIZE];
  Matrix* A3 = new (mem_allocA3) Matrix;

  unsigned int mem_allocek[NEW_SIZE];
  Matrix* ek = new (mem_allocek) Matrix;

  unsigned int mem_allocu[NEW_SIZE];
  Matrix* u = new (mem_allocu) Matrix;

  unsigned int mem_allocut[NEW_SIZE];
  Matrix* ut = new (mem_allocut) Matrix;

  unsigned int mem_allocuut[NEW_SIZE];
  Matrix* uut = new (mem_allocuut) Matrix;

  unsigned int mem_alloctwo_uut[NEW_SIZE];
  Matrix* two_uut = new (mem_alloctwo_uut) Matrix;

  unsigned int mem_allocu_aux1[NEW_SIZE];
  Matrix* u_aux1 = new (mem_allocu_aux1) Matrix;

  unsigned int mem_allocu_aux2[NEW_SIZE];
  Matrix* u_aux2 = new (mem_allocu_aux2) Matrix;

  unsigned int mem_allocI[NEW_SIZE];
  Matrix* I = new (mem_allocI) Matrix;

  unsigned int mem_allocIi[NEW_SIZE];
  Matrix* Ii = new (mem_allocIi) Matrix;

  unsigned int mem_allocP[NEW_SIZE];
  Matrix* P = new (mem_allocP) Matrix;

  unsigned int mem_allocQt[NEW_SIZE];
  Matrix* Qt = new (mem_allocQt) Matrix;

  unsigned int mem_allocQtt[NEW_SIZE];
  Matrix* Qtt = new (mem_allocQtt) Matrix;

  unsigned int mem_allocQt2[NEW_SIZE];
  Matrix* Qt2 = new (mem_allocQt2) Matrix;

Complex sig1, nx, nu;

Q->setEye(column);

A->row = row;
A->column = column;

  for (int h = 0; h < row; h++){
    for (int g = 0; g < column; g++){
      A->str[h][g] = str[h][g] ;
    }
  }

  for(int i = 0; i < (column-2); i++){
    x->get_part(i+1,column-1,i,i,A);
    int s = x->get_row();
    ek->zeros(s,1);
    ek->str[0][0] = val(1);

    Complex vector_x[s];
    int vectorx_size = 0;
    for (int a = 0; a < x->row; a++){
      for (int b = 0; b < x->column; b++){
        vector_x[vectorx_size] = x->str[a][b];
        vectorx_size += 1;
      }
    }

    sig1 = sign(x->str[0][0]);
    nx = norm_vector(s, vector_x);



    u_aux1->multE(ek, mul(sig1,nx));
    u_aux2->sum(x,u_aux1);


    int t = x->get_row();
    Complex vector_u[t];
    int vectoru_size = 0;
    for (int a = 0; a < u_aux2->row; a++){
      for (int b = 0; b < u_aux2->column; b++){
        vector_u[vectoru_size] = u_aux2->str[a][b];
        vectoru_size += 1;
      }
    }

    nu = norm_vector(t, vector_u);
    u->divE(u_aux2, nu);


    I->setEye(s);
    ut->transposed(u);
    uut->multiplication(u,ut);
    two_uut->multE(uut, val(2));
    P->subtraction(I, two_uut);


    Ii->setEye(i+1);
    Qt->blkdiag2(Ii,P);
    Qt2->multiplication(Qt,Q);
    Qtt->transposed(Qt);

    A2->multiplication(Qtt,A);
    A3->multiplication(A2,Qt);
    A->copy(A3);
    Q->copy(Qt2);

  }

  H->copy(A);

}


//--------------------------> GET EIGENVALUES <---------------------------------

void Matrix::eig(Matrix *in){

  int n = in->row;

  unsigned int mem_allocH[NEW_SIZE];
  Matrix* H = new (mem_allocH) Matrix;

  unsigned int mem_allocq[NEW_SIZE];
  Matrix* q = new (mem_allocq) Matrix;

  unsigned int mem_allocr[NEW_SIZE];
  Matrix* r = new (mem_allocr) Matrix;

  unsigned int mem_alloceaux_in[NEW_SIZE];
  Matrix* eaux_in = new (mem_alloceaux_in) Matrix;

  unsigned int mem_alloceaux_out[NEW_SIZE];
  Matrix* eaux_out = new (mem_alloceaux_out) Matrix;

  H->copy(in);

  for (int i = 0; i < 100; i++){
    H->qr(q,r);
    H->multiplication(r,q);

  }

  int counter = 0;
  row = n;
  column = 1;



  for (int j = 0; j < n; j++){
    if (counter < n){
      if (H->str[counter+1][counter].re() <= val(0.000015)){
        str[counter][0] = H->str[counter][counter];
        counter = counter + 1;
      }
      else{
        eaux_in->row = 2;
        eaux_in->column = 2;
        eaux_in->str[0][0] = H->str[counter][counter];
        eaux_in->str[0][1] = H->str[counter][counter+1];
        eaux_in->str[1][0] = H->str[counter+1][counter];
        eaux_in->str[1][1] = H->str[counter+1][counter+1];
        eaux_out->eig22(eaux_in);
        str[counter][0] = eaux_out->str[0][0];
        str[counter+1][0] = eaux_out->str[1][0];
        counter = counter + 2;

      }
    }
  }



}

//-------------> ZEROS ALL THE ELEMENTS BELOW THE DIAGONAL <-------------------------

void Matrix::triu(Matrix *in){

  row = in->row;
  column = in->column;

  for (int i = 0; i < row; i++){
    for (int j = 0; j < column; j++){
      if (i > j){
        str[i][j] = val(0);
      }
      else {
        str[i][j] = in->str[i][j];
      }
    }
  }

}


//-------------> ZEROS ALL THE ELEMENTS ABOVE THE DIAGONAL <-------------------------

void Matrix::tril(Matrix *in){

  row = in->row;
  column = in->column;

  for (int i = 0; i < row; i++){
    for (int j = 0; j < column; j++){
      if (i < j){
        str[i][j] = val(0);
      }
      else {
        str[i][j] = in->str[i][j];
      }
    }
  }


}



//---------------------------> TRACE OF A MATRIX <--------------------------------

Complex Matrix::trace(){

  Complex t;
  int n;

  n = row;
  t = val(0);

  for (int i = 0; i < n; i++){
    t = t + str[i][i];
  }

  return t;


}


//---------------------------> EIGENVALUES OF A 2X2 MATRIX (AUXILIARY) <--------------------------------

void Matrix::eig22(Matrix *in){

row = 2;
column = 1;

Complex t = in->trace();
Complex t2 = power(t,val(2));
Complex det = in->determinant(2);
Complex y = t2-mul(4,det);
Complex x1 = division((t+square_root(y)),val(2));
Complex x2 = division((t-square_root(y)),val(2));

str[0][0] = x1;
str[1][0] = x2;

}


//-----------------------------> 2-NORM OF A MATRIX  <-------------------------------------------------

Complex Matrix::twoNorm(){

  Complex sum, x, out;

  for(int i = 0; i < row; i++){
    for(int j = 0; j < column; j++){
      x = Complex(fabs(str[i][j]),val(0));
      sum += power(x,val(2));
    }
  }

  out = square_root(sum);
  return out;

}


//-----------------> DISCRETE-TIME ALGEBRAIC RICCATI EQUATION SOLVER <------------------------------
//            Finds X as X = A'*X*A - (B'*X*A)' * (R+B'*X*B)^-1 * B'*X*A + Q

void Matrix::riccati(Matrix *A, Matrix *B, Matrix *R, Matrix *Q)
{

  unsigned int mem_allocHk[NEW_SIZE];
  Matrix* Hk = new (mem_allocHk) Matrix;

  unsigned int mem_allocGk[NEW_SIZE];
  Matrix* Gk = new (mem_allocGk) Matrix;

  unsigned int mem_allocGk1[NEW_SIZE];
  Matrix* Gk1 = new (mem_allocGk1) Matrix;

  unsigned int mem_allocBt[NEW_SIZE];
  Matrix* Bt = new (mem_allocBt) Matrix;

  unsigned int mem_allocAt[NEW_SIZE];
  Matrix* At = new (mem_allocAt) Matrix;

  unsigned int mem_allocRi[NEW_SIZE];
  Matrix* Ri = new (mem_allocRi) Matrix;

  unsigned int mem_allocAk[NEW_SIZE];
  Matrix* Ak = new (mem_allocAk) Matrix;

  unsigned int mem_allocI[NEW_SIZE];
  Matrix* I = new (mem_allocI) Matrix;

  unsigned int mem_allocGH[NEW_SIZE];
  Matrix* GH = new (mem_allocGH) Matrix;

  unsigned int mem_allocGHI[NEW_SIZE];
  Matrix* GHI = new (mem_allocGHI) Matrix;

  unsigned int mem_allocIGHI[NEW_SIZE];
  Matrix* IGHI = new (mem_allocIGHI) Matrix;

  unsigned int mem_allocAIGHI[NEW_SIZE];
  Matrix* AIGHI = new (mem_allocAIGHI) Matrix;

  unsigned int mem_allocGAT[NEW_SIZE];
  Matrix* GAT = new (mem_allocGAT) Matrix;

  unsigned int mem_allocAIGHIAT[NEW_SIZE];
  Matrix* AIGHIAT = new (mem_allocAIGHIAT) Matrix;

  unsigned int mem_allocHk_one[NEW_SIZE];
  Matrix* Hk_one = new (mem_allocHk_one) Matrix;

  unsigned int mem_allocHIGHI[NEW_SIZE];
  Matrix* HIGHI = new (mem_allocHIGHI) Matrix;

  unsigned int mem_allocHIGHIA[NEW_SIZE];
  Matrix* HIGHIA = new (mem_allocHIGHIA) Matrix;

  unsigned int mem_allocAHIGHIA[NEW_SIZE];
  Matrix* AHIGHIA = new (mem_allocAHIGHIA) Matrix;

  unsigned int mem_allocGk_one[NEW_SIZE];
  Matrix* Gk_one = new (mem_allocGk_one) Matrix;

  unsigned int mem_allocAk_one[NEW_SIZE];
  Matrix* Ak_one = new (mem_allocAk_one) Matrix;

  unsigned int mem_allocHn[NEW_SIZE];
  Matrix* Hn = new (mem_allocHn) Matrix;

  int n = A->row;
  Complex norms = Complex(val(1), val(0));
  Complex n1, n2;
  typ_var epslon = val(0.0001);

  R->inverse(Ri);
  Bt->transposed(B);
  Gk1->multiplication(B,Ri);

  Gk->multiplication(Gk1,Bt);
  Ak->copy(A);
  Hk->copy(Q);
  I->setEye(n);

  while (norms.re() >= epslon){
    GH->multiplication(Gk,Hk);
    GHI->sum(I,GH);
    GHI->inverse(IGHI);
    AIGHI->multiplication(Ak,IGHI);
    Ak_one->multiplication(AIGHI,Ak);
    At->transposed(Ak);
    GAT->multiplication(Gk,At);
    AIGHIAT->multiplication(AIGHI,GAT);
    Gk_one->sum(Gk,AIGHIAT);
    HIGHI->multiplication(Hk,IGHI);
    HIGHIA->multiplication(HIGHI,Ak);
    AIGHIAT->multiplication(At,HIGHIA);
    Hk_one->sum(Hk,AIGHIAT);

    Hn->subtraction(Hk_one,Hk);
    n1 = Hn->twoNorm();
    n2 = Hk_one->twoNorm();
    norms = division(n1,n2);

    Hk->copy(Hk_one);
    Gk->copy(Gk_one);
    Ak->copy(Ak_one);

  }


  row = Hk->row;
  column = Hk->column;

  for(int i = 0; i < row; i++){
    for(int j = 0; j < column; j++){
      str[i][j] = Hk->str[i][j];
    }
  }


}
