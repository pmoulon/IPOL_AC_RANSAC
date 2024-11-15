//SPDX-License-Identifier: LGPL-3.0-or-later
/**
* @file svd.cpp
* @brief Clean implementation of SVD algorithm.
* @author Pascal Monasse
* 
* Copyright (c) 2024 Pascal Monasse
*/

#ifdef MATRIX_H // Do nothing if not included from matrix.h

#include <cmath>
#include <cassert>
#include <limits>

// --- Part 1: Bidiagonalization by Householder reflections ---

// Compute the Householder vector of x.
// The Householder matrix is H=I-2/|v|^2 v v^t with v = x +/- |x| e1.
// Then Hx = -/+ e1. The sign is chosen so that no cancellation occurs in v:
// We can assume + if we note |x| the actual signed value  sign(x1) |x|.
// Noting alpha = 2/|v|^2, we get alpha = 1 / (|x|^2 + |x|x1).
// However, we store the 1-normalized vector w=v/v1. Since v1=x1+|x|,
// beta = v1^2 alpha = (x1 + |x|)/|x| = 1+x1/|x| and H=I - beta w w^t.
// Return -|x| and x is overwritten: x[0]=beta and w = [1 x[1:]].
template <typename T>
T householder(T* x, int n, int stride) {
    T s=0;
    for(int i=0; i<n; i++) {
        T* p = x+i*stride;
        s += *p * *p;
    }
    s = sqrt(s);
    if(x[0]<0) s = -s;
    T beta = 0;
    if(s!=0) {
        beta = 1 + x[0]/s;
        T t = 1/(x[0]+s);
        for(int i=1; i<n; i++)
            x[i*stride] *= t;
    }
    x[0] = beta;
    return -s;
}

// Left-multiply matrix a by Householder matrix of vector v. First component of
// v is implicitly 1 and v[0] is actually beta.
template <typename T>
void left_mult_householder(const T* v, int dv,
                           T* a, int m, int n, int dm, int dn) {
    if(n<=0) return;
    T beta = *v;
    if(beta == 0)
        return;
    for(int j=0; j<n; j++) {
        T* p = a+j*dn;
        T s = p[0]; // Compute dot product of v and column j of a
        for(int i=1; i<m; i++)
            s += v[i*dv] * p[i*dm];
        s *= beta;
        p[0] -= s; // First component of v is implicitly 1
        for(int i=1; i<m; i++)
            p[i*dm] -= s*v[i*dv];
    }
}

// Multiply a on the right by Householder matrix.
template <typename T>
void right_mult_householder(const T* v, int dv,
                            T* a, int m, int n, int dm, int dn) {
    left_mult_householder(v, dv, a, n, m, dn, dm);
}

// Make matrix A bidiagonal by alternating Householder matrix multiplications
// on the left and right. At output, a is overwritten with Householder vectors,
// d (length n) stores the diagonal and e (length n-1) the upper-diagonal
// coefficients.
template <typename T>
void bidiagonalize(T* a, int m, int n, T* d, T* e) {
    for(int j=0; j<n; j++) {
        T* p = a+j*(n+1); 
        d[j] = householder(p, m-j, n);
        left_mult_householder(p, n, p+1, m-j, n-j-1, n, 1);
        if(j+1<n) {
            ++p;
            e[j] = householder(p, n-j-1, 1);
            right_mult_householder(p, 1, p+n, m-j-1, n-j-1, n, 1);
        }
    }
}

// Fill matrix U composed of the left multiplications encoded in columns of A.
// U is of size m x nu with nu=n (compact SVD) or nu=m (full SVD).
template <typename T>
void fill_householder_u(const T* a, int m, int n, T* u, int nu) {
    T* q=u;
    for(int i=0,l=m*nu; i<l; ++i) *q++ = 0.; // u <- 0
    for(int i=0; i<nu; i++) u[i*(nu+1)] = 1.; // diag(u) <- 1
    for(int i=n-1; i>=0; i--) {
        const T* p = a+i*(n+1);
        q = u+i*(nu+1);
        left_mult_householder(p, n, q, m-i, nu-i, nu, 1);
    }
}

// Fill matrix V composed of the right multiplications encoded in rows of A.
template <typename T>
void fill_householder_v(const T* a, int n, T* v) {
    T* q=v;
    for(int i=0,l=n*n; i<l; ++i) *q++ = 0.; // u <- 0
    for(int i=0; i<n; i++) v[i*(n+1)] = 1.; // diag(u) <- 1
    for(int i=0; i+1<n; i++) {
        const T* p = a+i*(n+1)+1;
        q = v+n+i+1; // First row and column of V always remain (1 0...0)
        right_mult_householder(p, 1, q, n-1, n-i-1, n, 1);
    }
}

// --- Part 2: QR algorithm with Givens rotations ---

// Givens rotation G=(c,-s; s,c) bringing u=(a,b)^t to v=(|u|,0)^t: v=Gu.
template <typename T>
static void givens(T a, T b, T& c, T& s) {
    if(b==0) { // Handle in particular a=b=0
        c=1; s=0;
        return;
    }
    T norm = hypot(a,b);
    c = a/norm;
    s = -b/norm;
}

// Put at 0 small values in the upper diagonal.
template <typename T>
static void chop_small_elements(const T* d, T* e, int n) {
    const T eps = std::numeric_limits<T>::epsilon();
    T d_i = d[0];
    for(int i=0; i+1<n; i++) {
        double d_ip1 = d[i+1];
        if(std::abs(e[i]) < eps*(std::abs(d_i)+std::abs(d_ip1)))
            e[i] = 0;
        d_i = d_ip1;
    }
}

// Put 0 in row k0 of bidiagonal matrix when diagonal element at k0 is 0:
// ( . . 0 0 )    ( . . 0 0 )
// ( 0 0 . 0 ) -> ( 0 0 0 0 )
// ( 0 0 . . )    ( 0 0 . . )
// ( 0 0 0 . )    ( 0 0 0 . )
// Done by applying Givens rotations on the left (row combinations).
template <typename T>
static void chase_out_intermediate_zero(T* d, T* e, int n, T* U, int m, int nu,
                                        int k0) {
    assert(d[k0]==0);
    T x=e[k0], y=d[k0+1], c, s;
    e[k0]=0; // Ensured by first Givens rotation    
    for(int k=k0; k+1<n; k++) {
        givens(y, -x, c, s);
        // U <- U G
        if(U)
            for(int i=0,j0=k0,j=k+1; i<m; i++, j0+=nu, j+=nu) {
                T w   = c*U[j0]-s*U[j];
                U[j]  = s*U[j0]+c*U[j];
                U[j0] = w;
            }
        // B <- G^t B
        d[k+1] = s*x+c*y;
        if(k+2<n) {
            T z=e[k+1];
            e[k+1] = c*z;
            x = -s*z; y=d[k+2];
        }
    }
}

// Put 0 in last column of bidiagonal matrix when bottom-right element is 0:
// ( . . 0 )    ( . . 0 )
// ( 0 . . ) -> ( 0 . 0 )
// ( 0 0 0 )    ( 0 0 0 )
// Done by applying Givens rotations on the right (column combinations).
template <typename T>
static void chase_out_trailing_zero(T* d, T* e, int n, T* V, int nv) {
    assert(d[n-1]==0);
    T x=d[n-2], y=e[n-2];
    e[n-2]=0; // Ensured by first Givens rotation
    for(int k=n-2; k>=0; k--) {
        T c,s;
        givens(x,y,c,s);
        // V <- V G
        if(V)
            for(int i=0,j0=k,j=n-1; i<nv; i++, j0+=nv, j+=nv) {
                T w  = c*V[j0]-s*V[j];
                V[j] = s*V[j0]+c*V[j];
                V[j0] = w;
            }
        // B <- B G
        d[k] = c*x-s*y;
        if(k>0) {
            T z = e[k-1];
            e[k-1] = c*z;
            x = d[k-1]; y = s*z;
        }
    }
}

// Given 3x2 matrix A, compute eigenvalue of A^t A closest to its last element:
//     (fa  0)
// A = (da fb) --> A^t A = (da2+fa2 da*fb  ) = (ta  tab)
//     (0  db)             (da*fb   db2+fb2)   (tab tb )
template <typename T>
static T trailing_eigenvalue(const T* d, const T* e, int n) {
    T da=d[n-2], db=d[n-1], fa=(n>2)?e[n-3]:0, fb=e[n-2]; // Coeff of A
    T da2=da*da, db2=db*db;
    T fa2=fa*fa, fb2=fb*fb;
    T ta = da2+fa2; // Coeff of A^t A
    T tb = db2+fb2;
    T tab = da*fb;
    T dt = (ta-tb)/2;

    T S = ta+tb; // Trace
    T P = da2*db2 + fa2*db2 + fa2*fb2; // Determinant
    T D = hypot(dt,tab); // Half square root of discriminant of X2-SX+P=0
    T r1 = S/2 + D;

    T mu=r1; // Larger root
    if (dt>=0) /* tb < ta, choose smaller root */
        mu = (r1>0)?  P/r1: 0;
    return mu;
}

// The bidiagonal matrix has non-zero values in the upper diagonal.
// Apply one QR step and get again a bidiagonal matrix: after the Givens
// rotation of first two columns, we get a bulge coefficient in subdiagonal,
// which moves above the upper diagonal though rotation of first two rows.
// Iterations amount to chasing the bulge.
// d (diagonal) has size n, e (upper-diagonal) has size n-1.
// U is m x nu
// V is nv x nv
template <typename T>
static void qr_step(T* d, T* e, int n, T* U, int m, int nu, T* V, int nv) {
    // If a 0 on the diagonal, chase out the upper-diagonal coeff on its right
    for(int i=0; i+1<n; i++)
        if(d[i]==0) {
            chase_out_intermediate_zero(d, e, n, U, m, nu, i);
            return;
        }
    // If a 0 at last position, put a 0 at the upper-diagonal coeff above
    if(d[n-1]==0) {
        chase_out_trailing_zero(d, e, n, V, nv);
        return;
    }
    T d0=d[0], e0=e[0], d1=d[1];
    T mu = trailing_eigenvalue(d, e, n);
    T y=d0*d0-mu, z=d0*e0;
    T bk=0, ap=d0, bp=e0, aq=d1;

    for(int k=0; k+1<n; k++) {
        // At k>0, we have bidiagonal matrix with bulge z:
        // ( ak bk  z )
        // (  0 ap bp ). (i) combine cols 2 and 3 to put z at 0
        // (  0  0 aq ). (ii) combine rows 2 and 3 to move bulge right of bp
        // Notice ak is not used.
        T c,s;
        givens(y, z, c, s); // (i) col combinations
        // V <- V G
        if(V)
            for(int i=0,j=k; i<nv; i++, j+=nv) {
                T w    = c*V[j]-s*V[j+1];
                V[j+1] = s*V[j]+c*V[j+1];
                V[j] = w;
            }
        // B <- B G
        if(k>0)
            e[k-1] = c*bk-s*z;
        y  = c*ap-s*bp;
        z = -s*aq;
        bk = s*ap+c*bp;
        ap = c*aq;
        bp = (k+2<n)? e[k+1]: 0;
        givens(y, z, c, s); // (ii) row combinations
        // U <- U G
        if(U)
            for(int i=0,j=k; i<m; i++, j+=nu) {
                T w    = c*U[j]-s*U[j+1];
                U[j+1] = s*U[j]+c*U[j+1];
                U[j] = w;
            }
        // B <- G^t B
        d[k] = c*y-s*z;
        T b=bk;
        bk = c*b-s*ap;
        ap = s*b+c*ap;
        y = bk;
        z = -s*bp;
        bp = c*bp;
        aq = (k+2<n)? d[k+2]: 0;
    }
    e[n-2]=bk;
    d[n-1]=ap;
}

// Apply implicit QR algo of a real bidiagonal matrix and update the
// orthogonal transformation matrices U and V.
// - d diagonal: vector of length n
// - e upper-diagonal: vector of length n-1
// - U matrix of size mxnu (can be null)
// - V matrix of size nxn (can be null)
template <typename T>
static int qr_bidiag(T* d, T* e, int n, T* U, int m, int nu, T* V) {
    chop_small_elements(d, e, n);
    int iter=0;
    int maxIter=100*n;
    for(int b=n-1; b>0 && iter<maxIter; iter++) {
        if(e[b-1] == 0) {
            --b;
            --iter; // Do not count as iteration
            continue;
        }
        assert(b>0);
        int a = b-1;
        while(a>0 && e[a-1]!=0)
            --a;
        int n2 = b-a+1;
        qr_step(d+a, e+a, n2, U?U+a:0, m, nu, V?V+a:0, n);
        chop_small_elements(d+a, e+a, n2);
    }
    return iter;
}

// Compute the singular value decomposition A=U*diag(D,nu,n)*V'.
// A is mxn (m>=n) -> U mxnu, D n and V nxn. nu=m (full SVD) or nu=n (compact).
// u and/or v can be null, if the reconstruction is not desired.
// Return the number of qrstep iterations.
template <typename T>
int svd(T* a, int m, int n, T* d, T* u, int nu, T* v) {
    assert(m>=n);

    T* e = new T[n-1];
    bidiagonalize(a, m, n, d, e);
    if(u)
        fill_householder_u(a, m, n, u, nu);
    if(v)
        fill_householder_v(a, n, v);
    int iter = qr_bidiag(d,e,n,u,m,nu,v);
    delete [] e;

    for(int i=0; i<n; ++i) { // Make positive values on diagonal
        if(d[i]<0.) {
            d[i] = -d[i];
            if(v)
                for(int j=0; j<n; ++j) {
                    T* p = v+n*j+i;
                    *p = -*p;
                }
        }
    }
    return iter;
}

#endif // MATRIX_H
