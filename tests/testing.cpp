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

void static func(vector<double>& d, vector<double>& e, vector<vector<double>>& z)
{
    int n = d.size();
    int m,l,iter,i,k;
    double s,r,p,g,f,dd,c,b;

    const double EPS=std::numeric_limits<double>::epsilon();
    for (i=1;i<n;i++) e[i-1]=e[i];
    e[n-1]=0.0;
    for (l=0;l<n;l++) {
        iter=0;
        do {
            for (m=l;m<n-1;m++) {
                dd=abs(d[m])+abs(d[m+1]);
                if (abs(e[m]) <= EPS*dd) break;
            }
            if (m != l) {
                if (iter++ == 30) throw("Too many iterations in tqli");
                g=(d[l+1]-d[l])/(2.0*e[l]);
                r=QLSolver::pythag(g,1.0);
                g=d[m]-d[l]+e[l]/(g+QLSolver::sign(r,g));
                s=c=1.0;
                p=0.0;
                for (i=m-1;i>=l;i--) {
                    f=s*e[i];
                    b=c*e[i];
                    e[i+1]=(r=QLSolver::pythag(f,g));
                    if (r == 0.0) {
                        d[i+1] -= p;
                        e[m]=0.0;
                        break;
                    }
                    s=f/r;
                    c=g/r;
                    g=d[i+1]-p;
                    r=(d[i]-g)*s+2.0*c*b;
                    d[i+1]=g+(p=s*r);
                    g=c*r-b;
                    for (k=0;k<n;k++) {
                        f=z[k][i+1];
                        z[k][i+1]=s*z[k][i]+c*f;
                        z[k][i]=c*z[k][i]-s*f;
                    }
                }
                if (r == 0.0 && i >= l) continue;
                d[l] -= p;
                e[l]=g;
                e[m]=0.0;
            }
        } while (m != l);
    }
}

void print_output(vector<double> eivalues, vector<vector<double>> eivectors) {
    //Eigenvalues stored in d1 and vectors in z1
    std::cout << "Eigenvalues:" <<std::endl;
    for (double val : eivalues)
        std::cout << val << " ";
    std::cout << std::endl;

    //Printing Eigenvectors
    std::cout << "Eigenvectors:" << std::endl;
    for (vector<double> a : eivectors) {
        for (double val : a)
            std::cout << val << " ";
        std::cout << std::endl;
    }
}

bool check_equal1D(vector<double>& a, vector<double>& b) {
    const double EPS = 0.1;
    if (a.size() != b.size())
        return false;
    vector<double> v1(a.begin(), a.end());
    vector<double> v2(b.begin(), b.end());
    std::sort(v1.begin(), v1.end());
    std::sort(v2.begin(), v2.end());

    int i;
    for (i=0; i<a.size(); i++) {
        if (std::abs(v1[i] - v2[i]) >= EPS) {
            return false;
        }
    }
    return true;
}

bool check_equal2d(vector<vector<double>>& a, vector<vector<double>>& b) {
    return false;
}

TEST_CASE("Pythag Test") {
    REQUIRE(5.0 == QLSolver::pythag(3.0, 4.0));
    REQUIRE(sqrt(2.0) == QLSolver::pythag(1.0, 1.0));
    REQUIRE(13.0 == QLSolver::pythag(5.0, 12.0));
}

TEST_CASE("QLSolver Diagonal Matrix") {
    // 1 0 0 0
    // 0 2 0 0
    // 0 0 3 0
    // 0 0 0 4

    // Eigenvalues should be 1, 2, 3, 4

    vector<double> d1 = {1, 2 ,3, 4};
    vector<double> e1 = {0, 0, 0, 0};
    vector<vector<double>> z1 = QLSolver::create_identity(d1.size()); // QLSolver requires


    QLSolver::tqli(d1, e1, z1);
    vector<double> out = {1, 2, 3, 4};

    REQUIRE(check_equal1D(d1, out));
}

TEST_CASE("QLSolver 3x3 non-trivial") {
    // 3 1 0
    // 1 2 2
    // 0 2 2
    vector<vector<double>> z1 = QLSolver::create_identity(3);
    vector<double> d1 = {3, 2, 2};
    vector<double> e1 = {0, 1, 2};

    QLSolver::tqli(d1, e1, z1);
    vector<double> expected = {4.39138,2.77287,-0.164248}; // Found in mathematica
    print_output(d1, z1); // print results
    REQUIRE(check_equal1D(d1, expected)); // if eigenvalues are correct, then surely
    // eigenvectors are too. I hope.
}

TEST_CASE("QLSolver 4x4 non-trivial") {
    // 1 3 0 0
    // 3 4 0 0
    // 0 0 2 1
    // 0 0 1 2
    vector<vector<double>> z1 = QLSolver::create_identity(4);
    vector<double> d1 = {1, 4, 2, 2};
    vector<double> e1 = {0, 3, 1, 1};
    QLSolver::tqli(d1, e1, z1);
    vector<double> expected = {6.04658,2.91719,1.0,-0.963772};
    print_output(d1, z1);
    REQUIRE(check_equal1D(d1, expected));
}