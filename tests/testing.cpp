//
// Created by Akash Piya on 3/23/24.
//

#include <catch2/catch_test_macros.hpp>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include "solver.h"

using std::vector;

void print_output(vector<double> d, vector<vector<double>> z) {
    //Eigenvalues stored in d1 and vectors in z1
    std::cout << "Eigenvalues:" <<std::endl;
    for (double val : d)
        std::cout << val << " ";
    std::cout << std::endl;

    //Printing Eigenvectors
    std::cout << "Eigenvectors:" << std::endl;
    for (vector<double> a : z) {
        for (double val : a)
            std::cout << val << " ";
        std::cout << std::endl;
    }
}

bool check_equal1D(vector<double>& a, vector<double>& b) {
    const double EPS = 0.00001;
    if (a.size() != b.size())
        return false;
    vector<double> v1(a.begin(), a.end());
    vector<double> v2(b.begin(), b.end());
    std::sort(v1.begin(), v1.end());
    std::sort(v2.begin(), v2.end());

    int i;
    for (int i=0; i<a.size(); i++) {
        if (!(std::abs(v1[i] - v2[i])) < EPS) {
            return false;
        }
    }
    return true;
}

bool check_equal2d(vector<vector<double>>& a, vector<vector<double>>& b) {
}

TEST_CASE("QLSolver Test 0") {
    // 2 0
    // 0 1 ---> Should have eigenvalues 1, 2 with corresponding eigenvectors
    vector<vector<double>> z1 = QLSolver::create_identity(2); // QLSolver requires
    // vector<vector<double>> z1(2, vector<double>(2, 0.0));
    // the diagonal entries, and the off-diagonal elements and an input identity matrix
    vector<double> d1 = {2, 1};
    vector<double> e1 = {0, 0};

    QLSolver::tqli(d1, e1, z1);

    print_output(d1, z1);
}

TEST_CASE("QLSolver 3x3 Test 1") {
    // 2 0
    // 0 1 ---> Should have eigenvalues 1, 2 with corresponding eigenvectors
    vector<vector<double>> z1 = QLSolver::create_identity(2); // QLSolver requires
    // vector<vector<double>> z1(2, vector<double>(2, 0.0));
    // the diagonal entries, and the off-diagonal elements and an input identity matrix
    vector<double> d1 = {2, 1};
    vector<double> e1 = {0, 0};

    QLSolver::tqli(d1, e1, z1);

    print_output(d1, z1);
}