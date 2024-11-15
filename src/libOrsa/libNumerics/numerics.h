//SPDX-License-Identifier: LGPL-3.0-or-later
/**
 * @file numerics.h
 * @brief Linear algebra: system solving by LU decomposition and SVD
 * @author Pascal Monasse
 * 
 * Copyright (c) 2010-2012, 2024 Pascal Monasse
 */

#ifndef NUMERICS_H
#define NUMERICS_H

#include "matrix.h"

namespace libNumerics {
    typedef double flnum;

    /// Solve system AX = B.
    bool solveLU(const matrix<flnum>& A, const vector<flnum>& B,
                 vector<flnum>& X);
    bool solveLU(matrix<flnum> A, vector<flnum>& B);

    /// Singular Value Decomposition: U diag(D) V, U in O(m), V in O(n), D>=0.
    class SVD {
    public:
        typedef int Mode;
        static const Mode full=0x0, noU=0x1, noV=0x2, compact=0x4;
        SVD(const matrix<flnum>& A, Mode mode=full);
        matrix<flnum> compose() const;
        flnum sv(int i) const;

        matrix<flnum> U, V;
        vector<flnum> D;
        int iter; ///< Number of qrsteps iterations to compute SVD.

        typedef matrix<flnum> Mat;
        typedef vector<flnum> Vec;
        // Static functions related to SVD
        static bool Nullspace(const Mat& A, Vec* nullspace,
                              double ratioExtremes=1e-2, double ratio2Min=.5);
        static flnum InvCond(const Mat& A);
        static void EnforceRank2_3x3(const Mat& A, Mat* ARank);
        static void Nullspace2_Remap33(const Mat& A, Mat& f1, Mat& f2);
    private:
        void sort();
    };

} // namespace libNumerics

#endif
