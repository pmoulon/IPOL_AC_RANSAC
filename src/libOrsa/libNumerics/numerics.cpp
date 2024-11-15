//SPDX-License-Identifier: LGPL-3.0-or-later
/**
 * @file numerics.cpp
 * @brief Linear algebra basics
 * @author Pascal Monasse, Pierre Moulon
 */

#include "numerics.h"
#include <cmath>
#include <vector>
#include <limits>
#include <algorithm>
#include <cassert>

namespace libNumerics {

inline flnum ABS(flnum x)
{ return (x >= 0)? x: -x; }

/// Resolution by LU decomposition with pivot.
bool solveLU(const matrix<flnum>& A, const vector<flnum>& B, vector<flnum>& X)
{
    X = B;
    return solveLU(A, X);
}

/// Replace X by A^{-1}X, by LU solver.
bool solveLU(matrix<flnum> A, vector<flnum>& X)
{
    assert(A.nrow() == A.ncol());
    int	n = A.nrow();
    vector<flnum> rowscale(n); // Implicit scaling of each row
    std::vector<int> permut(n,0); // Permutation of rows

    // Get the implicit scaling information of each row
    for(int i=0; i< n; i++) {
        flnum max = 0.0;
        for(int j=0; j< n; j++) {
            flnum tmp = ABS(A(i,j));
            if (tmp> max)
                max = tmp;
        }
        if(max == 0.0)
            return false;
        rowscale(i) = 1.0/max;
    }
    // Perform the decomposition
    for(int k=0; k < n; k++) {
        // Search for largest pivot element
        flnum max = rowscale(k)*ABS(A(k,k));
        int imax = k;
        for(int i=k+1; i < n; i++) {
            flnum tmp = rowscale(i)*ABS(A(i,k));
            if(tmp > max) {
                max = tmp;
                imax = i;
            }
        }
        if(max == 0.0)
            return false;

        // Interchange rows if needed
        if(k != imax) {
            A.swapRows(k, imax);
            rowscale(imax) = rowscale(k); // Scale of row k no longer needed
        }
        permut[k] = imax; // permut(k) was not initialized before
        flnum Akk = 1/A(k,k);
        for(int i=k+1; i < n; i++) {
            flnum tmp = A(i,k) *= Akk; // Divide by pivot
            for (int j=k+1;j < n; j++) // Reduce the row
                A(i,j) -= tmp*A(k,j);
        }
    }
    // Forward substitution
    for (int k=0; k < n; k++) {
        flnum sum = X(permut[k]);
        X(permut[k]) = X(k);
        for(int j = 0; j < k; j++)
            sum -= A(k,j)*X(j);
        X(k) = sum;
    }
    // Backward substitution
    for(int k=n-1; k >= 0; k--) {
        flnum sum = X(k);
        for(int j=k+1; j < n; j++)
            sum -= A(k,j)*X(j);
        X(k) = sum/A(k,k);
    }
    return true;
}

/// Decompose A into U diag(D) V^T with U(m,m) and V(n,n) orthonormal matrices.
SVD::SVD(const matrix<flnum>& A, SVD::Mode mode)
: U(), V(), D(std::min(A.nrow(),A.ncol())), iter(0)
{
    int p = std::min(A.nrow(),A.ncol());
    if((mode & noU) == 0)
        U = matrix<flnum>(A.nrow(), (mode & compact)? p: A.nrow());
    if((mode & noV) == 0)
        V = matrix<flnum>(A.ncol(), (mode & compact)? p: A.ncol());
    D.fill(0);
    iter = A.SVD(U, D, V);
    sort();
}

/// Recompose from SVD. This should be the initial matrix, except if U or V
/// was not computed, in which case a 0x0 matrix is returned..
matrix<flnum> SVD::compose() const
{
    return U * D.diag(U.ncol(), V.ncol()) * V.t();
}

/// Return ith greathest singular value. It is safer than picking directly in D
/// as it checks against out of bounds.
flnum SVD::sv(int i) const
{
    assert(0<=i && (i<U.nrow() || i<V.ncol()));
    if(i<D.nrow())
        return D(i);
    return 0;
}

/// Functor to sort elements of a vector by decreasing value.
class SVDElement {
public:
    SVDElement(const vector<flnum>& D, int i)
    : m_val(D(i)), m_i(i) {}
    bool operator<(const SVDElement& e) const
    { return (m_val>e.m_val); }

    flnum m_val;
    int m_i;
};

/// Sort SVD by decreasing order of singular value.
void SVD::sort()
{
    std::vector<SVDElement> vec;
    for(int i=0; i < D.nrow(); i++)
        vec.push_back( SVDElement(D, i) );
    std::sort(vec.begin(), vec.end());
    // Apply permutation
    for(int i=D.nrow()-1; i >=0; i--)
        if(vec[i].m_i != i) { // Find cycle of i
            vector<flnum> colU; if(U.nElements()) colU = U.col(i);
            vector<flnum> colV; if(V.nElements()) colV = V.col(i);
            const flnum w = D(i);
            int j = i;
            while(vec[j].m_i != i) {
                if(U.nElements())
                    U.paste(0,j, U.col(vec[j].m_i));
                if(V.nElements())
                    V.paste(0,j, V.col(vec[j].m_i));
                D(j) = D(vec[j].m_i);
                std::swap(j,vec[j].m_i);
            }
            vec[j].m_i = j;
            if(U.nElements())
                U.paste(0,j, colU);
            if(V.nElements())
                V.paste(0,j, colV);
            D(j) = w;
        }
}

/// Solve the linear system Ax = 0 with ||x|| = 1.0 via SVD.
/// Check if there is a single singular value such that SV/maxSV<ratioExtremes.
/// However if the two lowest SV satisfy this, accept if minSV/minSV2<ratio2Min.
/// This indicates whether the solution is unique.
bool SVD::Nullspace(const Mat& A, Vec *nullspace,
                    double ratioExtremes, double ratio2Min)
{
    if(A.nrow()+1<A.ncol()) // rank<n-1
        return false;
    SVD svd(A, SVD::noU);
    int last=A.lastCol();
    assert(last>=1);
    flnum minSV1=std::max(svd.sv(last),
			  std::numeric_limits<flnum>::epsilon()*svd.sv(0));
    flnum minSV2=std::max(svd.sv(last-1),
			  std::numeric_limits<flnum>::epsilon()*svd.sv(0));    
    if(minSV1>=ratioExtremes*svd.sv(0) || // No low SV
       (minSV2<ratioExtremes*svd.sv(0) && // 2+ low SV...
        minSV1>=ratio2Min*minSV2))        // ... and rather close
        return false;
    *nullspace = svd.V.col(svd.V.lastCol());
    return true;
}

/// Inverse of norm-2 condition value (ratio of extreme singular values)
flnum SVD::InvCond(const Mat& A) {
    assert(A.nrow()==A.ncol());
    SVD svd(A, SVD::noU | SVD::noV);
    return svd.D(svd.D.lastRow())/svd.D(0);
}

/// Make square matrix of rank<=2.
void SVD::EnforceRank2_3x3(const Mat& A, Mat *ARank)
{
    assert(A.nrow()==3 && A.ncol()==3);
    libNumerics::SVD svd(A);
    svd.D(svd.D.lastRow()) = 0;
    *ARank = svd.compose();
}

/// Save the two last nullspace vectors as 3x3 matrices.
void SVD::Nullspace2_Remap33(const Mat &A, Mat& f1, Mat& f2) {
    assert(A.ncol()==9);
    assert(f1.nrow()==3 && f1.ncol()==3);
    assert(f2.nrow()==3 && f2.ncol()==3);      
    SVD svd(A, SVD::noU);
    f1.read( svd.V.col(svd.V.lastCol()-1) );
    f2.read( svd.V.col(svd.V.lastCol()) );
}

} // namespace libNumerics
