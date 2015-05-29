#include <algorithm>
#include "CppUnitLite/TestHarness.h"
#include "libOrsa/cubicRoots.h"
#include "testing/testing.h"

static const int ITER=100000; // Number of tests
static const float RANGE=50.0f; // Amplitude of random roots
static const int PERCENT_FAIL=1; // Max failure rate

float genRand() {
    return (std::rand()/(float)RAND_MAX -.5f) * RANGE;
}

float eval(float coeffs[4], float x) {
    x = ((coeffs[3]*x+coeffs[2])*x
         +coeffs[1])*x + coeffs[0];
    return x;
}

float bound(float coeffs[4], float x) {
    float delta = 1E-3f *
        (std::abs(coeffs[0]) +
         std::abs(coeffs[1]*x) +
         std::abs(coeffs[2]*x*x));
    float deriv = coeffs[1] + 2.0f*coeffs[2]*x + 3.0f*coeffs[3]*x*x;
    return delta/std::abs(deriv);
}

// Test for polynomials with a single real root
TEST(CubicRoots, SingleRoot) {
    int fails=0;
    for(int i=0; i<ITER; i++) {
        float x = genRand();
        float a = genRand(); // Real part of complex root
        // Ensure that a*a+b*b > a*a with fp precision, otherwise a is also root
        float minb = 0.01f*std::abs(a);
        float b = std::rand()/(float)RAND_MAX * (RANGE-minb) + minb; // imaginary part

        // Roots: x, a+ib, a-ib
        float coeffs[4] = {
            -x*(a*a+b*b),
            x*(a+a)+(a*a+b*b),
            -(x+a+a),
            1.0f};
        float roots[3];
        int n = orsa::CubicRoots(coeffs, roots);
        CHECK_EQUAL(1,n);
        std::cout << i << ' ' << x << ' ' << roots[0];
        if(std::abs(x-roots[0]) >= bound(coeffs,x)) {
            ++fails;
            std::cout << " (fail)";
        }
        std::cout << std::endl;
    }
    EXPECT_TRUE(fails*100 <= ITER*PERCENT_FAIL);
}

// Test for polynomials with three real roots
TEST(CubicRoots, ThreeRoots) {
    int fails=0;
    for(int i=0; i<ITER; i++) {
        float x[3] = { genRand(), genRand(), genRand() };
        float coeffs[4] = {
            -x[0]*x[1]*x[2],
            x[0]*x[1]+x[0]*x[2]+x[1]*x[2],
            -(x[0]+x[1]+x[2]),
            1.0f};

        // Make sure that fp rounding does not perturbate the three roots
        float d = coeffs[2]*coeffs[2]-3*coeffs[1];
        bool ok = (d>0);
        if(ok) {
            float min = (-coeffs[2]-std::sqrt(d))/3,
                  max = (-coeffs[2]+std::sqrt(d))/3;
            ok = (eval(coeffs,min)> RANGE*1.0e-5) &&
                 (eval(coeffs,max)<-RANGE*1.0e-5);
        }
        if(! ok) { // Regenerate in case of bad rounding
            --i;
            continue;
        }

        float roots[3];
        int n = orsa::CubicRoots(coeffs, roots);
        std::sort(roots, roots+n);
        std::sort(x, x+3);
        for(int j=0; j<n; j++) {
            std::cout << i << ' ' << x[j] << ' ' << roots[j];
            if(std::abs(x[j]-roots[j]) >= bound(coeffs,x[j])) {
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
    EXPECT_TRUE(fails*100 <= ITER*PERCENT_FAIL);
}

// Test for polynomials with a triple root
TEST(CubicRoots, TripleRoot) {
    int fails=0;
    for(int i=0; i<ITER; i++) {
        float x = genRand();
        float coeffs[4] = {
            -x*x*x,
            3*x*x,
            -3*x,
            1.0f};
        float roots[3];
        int n = orsa::CubicRoots(coeffs, roots);
        CHECK_EQUAL(1,n);
        std::cout << i << ' ' << x << ' ' << roots[0];
        if(std::abs(x-roots[0]) >= bound(coeffs,x)) {
            ++fails;
            std::cout << " (fail)";
        }
        std::cout << std::endl;
    }
    EXPECT_TRUE(fails*100 <= ITER*PERCENT_FAIL);
}

// Test for polynomials with one double real root and another root.
// This is tricky, because it is the limit case between 1 and 3 solutions, and
// fp rounding could tip the balance one side or the other.
// For this reason, the test is less strict: we require that one returned root
// is close to the ground truth single root.
TEST(CubicRoots, DoubleRoot) {
    int fails=0;
    for(int i=0; i<ITER; i++) {
        float x[3] = { genRand(), genRand() };
        x[2] = x[1];
        float coeffs[4] = {
            -x[0]*x[1]*x[2],
            x[0]*x[1]+x[0]*x[2]+x[1]*x[2],
            -(x[0]+x[1]+x[2]),
            1.0f};

        // Make sure that fp rounding does not make a triple root
        float d = coeffs[2]*coeffs[2]-3*coeffs[1];
        bool ok = (d>0);
        if(ok) {
            float min = (-coeffs[2]-std::sqrt(d))/3;
            float max = (-coeffs[2]+std::sqrt(d))/3;
            float bound = std::max(RANGE,std::abs(coeffs[0]))*1e-5f;
            ok = (std::abs(eval(coeffs,min))> bound) ||
                 (std::abs(eval(coeffs,max))> bound);
        }
        if(! ok) { // Regenerate in case of bad rounding
            --i;
            continue;
        }     

        float roots[3];
        int n = orsa::CubicRoots(coeffs, roots);
        // x[0] must be found
        int jmin = 0;
        for(int j=1; j<n; j++)
            if(std::abs(x[0]-roots[j]) < std::abs(x[0]-roots[jmin]))
                jmin = j;
        std::cout << i << ' ' << x[0] << ' ' << roots[jmin];
        if(std::abs(x[0]-roots[jmin]) >= bound(coeffs,x[0])) {
            ++fails;
            std::cout << " (fail)";
        }
        std::cout << std::endl;
    }
    EXPECT_TRUE(fails*100 <= ITER*PERCENT_FAIL);
}

/* ************************************************************************* */
int main() { TestResult tr; return TestRegistry::runAllTests(tr);}
/* ************************************************************************* */
