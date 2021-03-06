/*
 * =====================================================================================
 *
 *       Filename:  Executor.cpp
 *
 *    Description:  Executes patterns found in codes from a differentiated matrix 
 *
 *        Version:  1.0
 *        Created:  2021-07-13 09:25:02 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Zachary Cetinic, 
 *   Organization:  University of Toronto
 *
 * =====================================================================================
 */

#include "DDT.h"
#include "Executor.h"
#include "Inspector.h"
#include "SpMVGenericCode.h"
#include "SpMMGenericCode.h"
#include "GenericCodeletsSpTRSV.h"
#include "SparseMatrixIO.h"

#include <chrono>
#include <vector>

namespace DDT {
void executeSpTRSCodelets(const std::vector<DDT::Codelet*>* cl, const DDT::Config& c) {
    // Read matrix
    auto m = readSparseMatrix<CSR>(c.matrixPath);

    // Setup memory
    auto x = new double[m.c]();
    for (int i = 0; i < m.c; i++) {
        x[i] = 1;
    }

    // Execute SpMV
    DDT::sptrsv_generic(m.r, m.Lp, m.Li, m.Lx, x, cl[0], c);

    // Clean up memory
    delete[] x;
}

void executeSPMVCodelets(const std::vector<DDT::Codelet*>* cl, const DDT::Config& c) {
  // Read matrix
  auto m = readSparseMatrix<CSR>(c.matrixPath);

  // Setup memory
  auto x = new double[m.c]();
  for (int i = 0; i < m.c; i++) {
      x[i] = i;
  }
  auto y = new double[m.r]();

  // Execute SpMV
  DDT::spmv_generic(m.r, m.Lp, m.Li, m.Lx, x, y, cl, c);

  // Clean up memory
  delete[] x;
  delete[] y;
}

void executeSPMMCodelets(const std::vector<DDT::Codelet*>* cl, const
DDT::Config& cfg, const int r, const int* Lp, const int *Li, const double*Ax,
const double* Bx, double* Cx, int bRows, int bCols) {
    // Execute SpMV
    DDT::spmm_generic(r, Lp, Li, Ax, Bx, Cx, bRows, bCols, cl, cfg);
}

void executeSPSPMMCodelets(const std::vector<DDT::Codelet*>* cl, const
DDT::Config& cfg, const int r, const int* Lp, const int *Li, const double*Lx,
const double* x, double* y) {
    // Execute SpMV
//    spspmm_generic(r, Lp, Li, Lx, x, y, cl, cfg);
}


void executeSPMVCodelets(const std::vector<DDT::Codelet*>* cl, const
DDT::Config& cfg, const int r, const int* Lp, const int *Li, const double*Lx,
const double* x, double* y) {
  // Execute SpMV
  spmv_generic(r, Lp, Li, Lx, x, y, cl, cfg);
 }

 void executeSPTRSVCodelets(const std::vector<DDT::Codelet*>* cl, const
 DDT::Config& cfg, const int r, const int* Lp, const int *Li, const double*Lx,
                             double* x, double* y) {
     sptrsv_generic(r, Lp, Li, Lx, x, cl[0], cfg);
 }

 void executeParallelSPTRSVCodelets(const DDT::GlobalObject& d, const DDT::Config& cfg, int r, const int* Lp, const int* Li, const double* Lx, double* x) {
     for (int i = 0; i < d.sm->_final_level_no; ++i) {
#pragma omp parallel for num_threads(cfg.nThread)
         for (int j = 0; j < d.sm->_wp_bounds[i]; ++j) {
             auto& cc = d.sm->_cl[i][j];
             sptrsv_generic(r, Lp, Li, Lx, x, cc, cfg);
         }
     }
 }


/**
 * @brief Executes codelets found in a matrix performing a computation
 *
 * @param cl List of codelets to perform computation on
 * @param c  Configuration object for setting up executor
 */
void executeCodelets(const std::vector<DDT::Codelet*>* cl, const DDT::Config& cfg) {
  switch (cfg.op) {
    case DDT::OP_SPMV:
      executeSPMVCodelets(cl, cfg);
      break;
    case DDT::OP_SPTRS:
      executeSpTRSCodelets(cl, cfg);
    default:
      break;
  }
}
 void executeCodelets(const std::vector<DDT::Codelet*>* cl, const DDT::Config&
 cfg, Args& args) {
  switch (cfg.op) {
   case DDT::OP_SPMV:
    executeSPMVCodelets(cl, cfg,args.r, args.Lp, args.Li, args.Lx, args.x, args
    .y);
    break;
   case DDT::OP_SPTRS:
       executeSPTRSVCodelets(cl, cfg, args.r, args.Lp, args.Li, args.Lx, args.x, args.y);
   default:
    break;
   case DDT::OP_SPMM:
       executeSPMMCodelets(cl, cfg, args.r, args.Lp, args.Li, args.Ax, args.Bx, args
       .Cx, args.bRows, args.bCols);
       break;
  }
 }

    void executeParallelCodelets(const DDT::GlobalObject& d, const DDT::Config&
    cfg, Args& args) {
        switch (cfg.op) {
            case DDT::OP_SPMV:
//                executeSPMVCodelets(cl, cfg,args.r, args.Lp, args.Li, args.Lx, args.x, args
//                        .y);
                break;
            case DDT::OP_SPTRS:
                executeParallelSPTRSVCodelets(d, cfg, args.r, args.Lp, args.Li, args.Lx, args.x);
            default:
                break;
        }
    }



}
