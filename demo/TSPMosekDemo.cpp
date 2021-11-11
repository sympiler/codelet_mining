//
// Created by kazem on 11/10/21.
//
#include <iostream>
#include <list>
#include <vector>

#include <Input.h>
#include <def.h>
#include <sparse_io.h>
#include <sparse_utilities.h>
#include "fusion.h"
#include "MosekTSP.h"


using namespace mosek::fusion;
using namespace monty;

int main(int argc, char *argv[]){
 auto config = DDT::parseInput(argc, argv);
 auto A = sym_lib::read_mtx(config.matrixPath);
 sym_lib::CSC *A_full = NULLPNTR;
 sym_lib::CSR *B = NULLPNTR, *L_csr = NULLPNTR;
 if (A->stype < 0) {
  A_full = sym_lib::make_full(A);
  B = sym_lib::csc_to_csr(A_full);
 } else {
  A_full = A;
  B = sym_lib::csc_to_csr(A);
 }
 //auto A_i = new_array_ptr<int, 1>({0, 1, 2, 3, 1, 0, 2, 0});
 //auto A_j = new_array_ptr<int, 1>({1, 2, 3, 0, 0, 2, 1, 3});
 auto A_j = new_array_ptr<int, 1>(A_full->nnz);
 auto A_i = new_array_ptr<int, 1>(A_full->nnz);
 auto A_x = new_array_ptr<double, 1>(A_full->nnz);

 for (int i = 0; i < A->n; ++i) {
  for (int j = A_full->p[i]; j < A_full->p[i+1]; ++j) {
   (*A_i)(j) = i;
   (*A_j)(j) = A_full->i[j];
   (*A_x)(j) = A->x[j];
  }
 }
 //auto A_x = new_array_ptr<double, 1>({1., 1., 1., 1., 0.1, 0.1, 0.1, 0.1});
 //int n = 4;
 int n = A_full->n;
 MOSEK::tsp(n, Matrix::sparse(n, n, A_i, A_j, 1.),
            Matrix::sparse(n, n, A_i,A_j,A_x),
            true, false);
 MOSEK::tsp(n, Matrix::sparse(n, n, A_i, A_j, 1.),
            Matrix::sparse(n, n, A_i,A_j,A_x),
            true, true);

 delete A_full;
 delete B;
 return 0;
}