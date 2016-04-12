#define _USE_MATH_DEFINES //For Windows (M_PI)
#include <algorithm>
#include "CppUnitLite/TestHarness.h"
#include "cubicRoots.h"
#include <cerrno>
#include <cstring>

static const int ITER=100000; // Number of tests
static const float RANGE=50.0f; // Amplitude of random roots
static const int PERCENT_FAIL=1; // Max failure rate

// Random number in [-RANGE/2,RANGE/2]
float genRand() {
    return (std::rand()/(float)RAND_MAX -.5f) * RANGE;
}

// Evaluate polynomial at x
float eval(float coeffs[4], float x) {
    x = ((coeffs[3]*x+coeffs[2])*x
         +coeffs[1])*x + coeffs[0];
    return x;
}

// Error tolerance near a single root
float bound(float coeffs[4], float x) {
    float delta = std::numeric_limits<float>::epsilon()*RANGE/2*
        (std::abs(coeffs[0]) * 3*RANGE/2*RANGE/2 +
         std::abs(x) * 6*RANGE/2 +
         std::abs(x*x) * 3);
    float deriv = coeffs[1] + 2.0f*coeffs[2]*x + 3.0f*coeffs[3]*x*x;
    return delta/std::abs(deriv);
}

// Test for polynomials with a single real root
TEST(CubicRoots, SingleRoot) {
    errno = 0;
    int fails=0;
    for(int i=0; i<ITER; i++) {
        float x = genRand();
        float a = genRand(); // Real part of complex root
        // Ensure that a*a+b*b > a*a with fp precision, otherwise a is also root
        float minb = 0.01f*std::abs(a);
        float b = std::rand()/(float)RAND_MAX * (RANGE-minb) + minb;

        // Roots: x, a+ib, a-ib
        float coeffs[4] = {
            -x*(a*a+b*b),
            x*(a+a)+(a*a+b*b),
            -(x+a+a),
            1.0f};
        float roots[3];
        int n = orsa::CubicRoots(coeffs, roots);
        CHECK_EQUAL(1,n);
        float err = std::abs(x-roots[0]);
        std::cout << i << ' ' << x << ' ' << roots[0] << ' ' << err;
        if(err >= bound(coeffs,x)) {
            ++fails;
            std::cout << " (fail)";
        }
        std::cout << std::endl;
    }
    CHECK(fails*100 <= ITER*PERCENT_FAIL);
    CHECK(errno==0);
    if(errno)
        std::cout << strerror(errno) << std::endl;
}

// Test for polynomials with three real roots
TEST(CubicRoots, ThreeRoots) {
    errno = 0;
    int fails=0;
    for(int i=0; i<ITER; i++) {
        float x[3] = { genRand(), genRand(), genRand() };
        float coeffs[4] = {
            -x[0]*x[1]*x[2],
            x[0]*x[1]+x[0]*x[2]+x[1]*x[2],
            -(x[0]+x[1]+x[2]),
            1.0f};

        float roots[3];
        int n = orsa::CubicRoots(coeffs, roots);
        std::sort(roots, roots+n);
        std::sort(x, x+3);
        for(int j=0; j<n; j++) {
            float err = std::abs(x[j]-roots[j]);
            std::cout << i << ' ' << x[j] << ' ' << roots[j] << ' ' << err;
            if(err >= bound(coeffs,x[j])) {
                ++fails;
                std::cout << " (fail)";
            }
            std::cout << std::endl;
        }
        if(n!=3) {
            std::cout << "Failure: 1 root found instead of 3" << std::endl;
            ++fails;
        }
    }
    CHECK(fails*100 <= ITER*PERCENT_FAIL);
    CHECK(errno==0);
    if(errno)
        std::cout << strerror(errno) << std::endl;
}

// Test for polynomials with a triple root. Finding more than 1 root is accepted
// but additional roots must evaluate to a small value.
TEST(CubicRoots, TripleRoot) {
    errno = 0;
    int fails=0;
    for(int i=0; i<ITER; i++) {
        float x = genRand();
        float coeffs[4] = {
            0, //-x*x*x,
            3*x*x,
            -3*x,
            1.0f};
        coeffs[0] = -eval(coeffs,x); // Adjusting poly to yield 0 at triple root
        float roots[3];
        int n = orsa::CubicRoots(coeffs, roots);
        for(int j=0; j<n; j++) {
            float err = std::abs(x-roots[j]);
            std::cout << i << ' ' << x << ' ' << roots[j] << ' ' << err;
            if(eval(coeffs,roots[j]) >= 1E-2f) {
                ++fails;
                std::cout << " (fail)";
            }
            std::cout << std::endl;
        }
    }
    CHECK(fails*100 <= ITER*PERCENT_FAIL);
    CHECK(errno==0);
    if(errno)
        std::cout << strerror(errno) << std::endl;
}

// Test for polynomials with one double real root and another root.
// This is tricky, because it is the limit case between 1 and 3 solutions, and
// fp rounding could tip the balance one side or the other.
// For this reason, the test is less strict: we require that one returned root
// is close to the ground truth single root, and other possibly found roots must
// evaluate to a small value.
TEST(CubicRoots, DoubleRoot) {
    errno = 0;
    int fails=0;
    for(int i=0; i<ITER; i++) {
        float x[3] = { genRand(), genRand() };
        x[2] = x[1];
        float coeffs[4] = {
            0, //-x[0]*x[1]*x[2],
            x[0]*x[1]+x[0]*x[2]+x[1]*x[2],
            -(x[0]+x[1]+x[2]),
            1.0f};
        coeffs[0] = -eval(coeffs,x[1]); // Adjusting to yield 0 at double root

        float roots[3];
        int n = orsa::CubicRoots(coeffs, roots);
        // x[0] must be found
        int jmin = 0;
        for(int j=1; j<n; j++)
            if(std::abs(x[0]-roots[j]) < std::abs(x[0]-roots[jmin]))
                jmin = j;
        float err = std::abs(x[0]-roots[jmin]);
        std::cout << i << ' ' << x[0] << ' ' << roots[jmin] << ' ' << err;
        if(err >= bound(coeffs,x[0])) {
            ++fails;
            std::cout << " (fail)";
        }
        std::cout << std::endl;
        // Other found roots must be close to double root x[1]=x[2]
        for(int k=1; k<n; k++) {
            float r = roots[(jmin+k)%n];
            float err = std::abs(x[1]-r);
            std::cout << i << ' ' << x[1] << ' ' << r << ' ' << err;
            if(eval(coeffs,r) >= 1E-2f) {
                ++fails;
                std::cout << " (fail)";
            }
            std::cout << std::endl;
        }
    }
    CHECK(fails*100 <= ITER*PERCENT_FAIL);
    CHECK(errno==0);
    if(errno)
        std::cout << strerror(errno) << std::endl;
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
