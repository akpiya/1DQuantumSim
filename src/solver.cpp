//
// Created by Akash Piya on 3/23/24.
//

#include "solver.h"
#include <vector>
#include <cmath>
using std::vector;

double QLSolver::pythag(const double a, const double b) {
    double absa=abs(a), absb=abs(b);
    if (absa > absb) {
        return absa * sqrt(1.0 + (absb/absa) * (absb/absa));
    } else {
        if (absa == 0.0) {
            return 0.0;
        } else {
            return absb * sqrt(1.0 + (absa/absb) * (absa/absb));
        }
    }
}

void QLSolver::tqli(vector<double>* d_ptr, vector<double>* e_ptr, vector<vector<double>>* z_ptr) {
    vector<double>& d = *d_ptr;
    vector<double>& e = *e_ptr;
    vector<double>& z = *z_ptr;

    int n = d.size();
    int m,l,iter,i,k;
    double s,r,p,g,f,dd,c,b;

    const double EPS = std::numeric_limits<double>::epsilon();

    //Assume d is of n length and e is of n length but only the first n-1 values matter
    for (l=0; l<n; l++) {
        iter=0;
        do {
            for (m=l;m<=n-1;m++) {
                dd = abs(d[m]) + abs(d[m+1]);
                if (abs(e[m]) <= EPS * dd)
                    break;
            }
            if (m!=l) {
                if (iter++ == 30) throw("Too many iterations");
                g = (d[l+1] - d[l]) / (2.0 * e[l]);
                r = pythag(g, 1.0);
                g = d[m] - d[l] + e[l]/(g+SIGN(r,g));

                s=c=1.0;
                p=0.0;
                for (i=m-1;i>=1;i--) {
                    f=s*e[i];
                    b=c*e[i];
                    e[i+1] = (r=pythag(f,g));
                    if (r == 0.0) {
                        d[i+1] -= p;
                        e[m] = 0.0;
                        break;
                    }
                    s = f/r;
                    c = g/r;
                    g = d[i+1] - p;
                    r = (d[i]-g)*s + 2.0*c*b;
                    d[i+1] = g+(p=s*r);
                    g = c*r - b;
                    for (k=0; k < n; k++) {
                        f = z[k][i+1];
                        z[k][i+1] = s*z[k][i] + c * f;
                        z[k][i] = c * z[k][i] - s * f;
                    }
                }
            }
            if (r == 0.0 && i >= 1) continue;
            d[l] -= p;
            e[l] = g;
            e[m] = 0.0;
        } while (m != l);
    }
}