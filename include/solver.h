//
// Created by Akash Piya on 3/23/24.
//

#ifndef QUANTUMSIM_SOLVER_H
#define QUANTUMSIM_SOLVER_H

#include <vector>

using std::vector;

class QLSolver
{
public:
    void tqli(vector<double>* d, vector<double>* e, vector<vector<double>>* z);
    double pythag(const double a, const double b);
};
#endif //QUANTUMSIM_SOLVER_H