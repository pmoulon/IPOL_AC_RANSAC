/**
 * @file cubicRoots.h
 * @brief Find roots of a cubic polynomial
 * @author Pascal Monasse
 *
 * Copyright (c) 2015 Pascal Monasse
 * All rights reserved.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef CUBICROOTS_H_
#define CUBICROOTS_H_

#define _USE_MATH_DEFINES //For Windows (M_PI)
#include <cmath>
#include <limits>

namespace orsa {

/// Find roots of the cubic polynomial
/// \f[ x^3 + a*x^2 + b*x + c = 0. \f]
/// The number of roots is returned.
template<typename Real>
int SolveCubicPolynomial(Real a, Real b, Real c, Real x[3]) {
    const Real eps = std::numeric_limits<Real>::epsilon();
    a /= 3;
    Real p = (b-3*a*a)/3;
    Real q = (2*a*a*a - a*b + c)/2;
    Real d = q*q+p*p*p;
    Real tolq = std::max(std::abs(2*a*a*a),std::max(std::abs(a*b),std::abs(c)));
    Real tolp = std::max(std::abs(b),std::abs(3*a*a));
    int n = (d>eps*std::max(p*p*tolp,std::abs(q)*tolq)? 1: 3);
    if(n==1) { // Single root: Cardano's formula
        d = std::pow(std::abs(q)+std::sqrt(d), 1/(Real)3);
        x[0] = d - p/d;
        if(q>0)
            x[0] = -x[0];
    } else { // Three roots: Viete's formula
        if(3*p>=-eps*tolp) { // p=0 and d<=0 implies q=0: triple root 0
            n = 1;
            x[0] = 0;
        } else {
            p = std::sqrt(-p);
            q /= p*p*p;
            d = Real((q<=-1)? M_PI: (q>=1)? 0: std::acos(q));
            for(int i=0; i<3; i++)
                x[i] = Real(-2*p*std::cos((d+2*M_PI*i)/3));
        }
    }
    for(int i=0; i<n; i++)
        x[i] -= a;
    return n;
}

/// Find real roots of cubic polynomial equation.
/// \f[ coeffs[3] x^3 + coeffs[2] x^2 + coeffs[1] x + coeffs[0] = 0. \f]
/// The coefficients are in ascending powers, i.e. coeffs[i]*x^i.
/// Return the number of solutions (1 or 3). A triple root is reported only once
/// but a double root is reported twice (so 3 solutions).
template<typename Real>
int CubicRoots(const Real coeffs[4], Real roots[3]) {
    if(coeffs[0] == 0.0)
        return 0;
    Real a = coeffs[2] / coeffs[3];
    Real b = coeffs[1] / coeffs[3];
    Real c = coeffs[0] / coeffs[3];
    return SolveCubicPolynomial(a, b, c, roots);
}

}  // namespace orsa

#endif  // LIBS_NUMERIC_POLY_H_
