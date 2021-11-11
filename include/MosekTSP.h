//
// Created by kazem on 11/10/21.
//

#ifndef DDT_MOSEKTSP_H
#define DDT_MOSEKTSP_H

using namespace mosek::fusion;
using namespace monty;
namespace MOSEK{

 void tsp(int n, Matrix::t A, Matrix::t C, bool remove_1_hop_loops, bool
 remove_2_hop_loops);
}
#endif //DDT_MOSEKTSP_H
