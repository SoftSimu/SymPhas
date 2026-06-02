// Autogen stencil reduction & correctness regression test.
// See testautogenstencil.h.

#include "testautogenstencil.h"

#include "symphas.h"

#include <cmath>
#include <cstdio>
#include <type_traits>
#include <typeinfo>

namespace {

static int g_failures = 0;
static int g_subtests = 0;

inline void check_close(double got, double want, double tol, const char* tag) {
    ++g_subtests;
    const double diff = std::fabs(got - want);
    const double scale = std::max({std::fabs(got), std::fabs(want), 1.0});
    if (!std::isfinite(got) || diff / scale > tol) {
        std::fprintf(stderr,
                     "  [FAIL] %-60s  got %.6g  want %.6g  rel %.3e\n",
                     tag, got, want, diff / scale);
        ++g_failures;
    }
}

template <typename T>
inline bool eval_returns_symbol() {
    using ret_t = decltype(std::declval<T const&>().eval(0));
    return std::is_same_v<ret_t, expr::symbols::Symbol>;
}

// SFINAE probe: does Gen accept (T*, Stride, double) and produce double?
template <typename Gen, typename Stride, typename = void>
struct gen_invokable : std::false_type {};

template <typename Gen, typename Stride>
struct gen_invokable<
    Gen, Stride,
    std::void_t<decltype(std::declval<Gen const&>()(
        std::declval<double const*>(),
        std::declval<Stride const&>(),
        std::declval<double>()))>>
    : std::is_same<decltype(std::declval<Gen const&>()(
                       std::declval<double const*>(),
                       std::declval<Stride const&>(),
                       std::declval<double>())),
                   double> {};

// --------------------------------------------------------------------------
// 1D derivative stencils.  For (O, OA=2, D=1):
//   stencil_apply_type<O, 2, 1> { v[-R..+R], stride, divh } -> derivative
// We probe by constructing a 1D buffer holding p(x) = sum_k a_k x^k and
// applying the stencil at the center.
// --------------------------------------------------------------------------

// Evaluate polynomial p(x) = c0 + c1 x + c2 x^2 + c3 x^3 + c4 x^4 at x.
inline double poly1d(double x, const double (&c)[5]) {
    double y = 0.0;
    double xp = 1.0;
    for (int k = 0; k < 5; ++k) { y += c[k] * xp; xp *= x; }
    return y;
}

// Analytic derivative of degree-4 polynomial: d^O p / dx^O at x.
inline double poly1d_dn(int O, double x, const double (&c)[5]) {
    double y = 0.0;
    for (int k = O; k < 5; ++k) {
        double fact = 1.0;
        for (int j = 0; j < O; ++j) fact *= (k - j);
        y += fact * c[k] * std::pow(x, k - O);
    }
    return y;
}

template <size_t O>
void test_1d_axial_derivative(double h, const double (&c)[5], const char* tag) {
    // Build the autogen stencil dictionary for d^O / dx^O, OA=2, D=1.
    auto dict = expr::get_central_space_stencil<O, 2, 1>();

    // GeneratedStencilApply wraps the dictionary so we can call it on a buffer.
    using gen_t = symphas::internal::GeneratedStencilApply<decltype(dict)>;

    // We need enough samples to support O+2 stencil reach on each side.
    // For OA=2 the half-width R = ceil(O/2)+1 worst-case 3, so allocate 9.
    constexpr int N = 9;
    double buf[N] = {0};
    // Place the center sample at index 4.
    const int c0 = N / 2;
    for (int i = 0; i < N; ++i) {
        const double x = (i - c0) * h;
        buf[i] = poly1d(x, c);
    }

    const double divh = 1.0 / h;
    char buftag[160];
    std::snprintf(buftag, sizeof(buftag), "%s d^%zu/dx^%zu (1D, OA=2)", tag, O, O);
    if constexpr (gen_invokable<gen_t, len_type>::value) {
        const double got = gen_t{}(buf + c0, /*stride=*/len_type{1}, divh);
        const double want = poly1d_dn(O, 0.0, c);
        check_close(got, want, 1e-9, buftag);
    } else {
        ++g_subtests; ++g_failures;
        std::fprintf(stderr, "  [FAIL] %-60s  Symbol leak (dictionary not reduced)\n", buftag);
    }
}

// --------------------------------------------------------------------------
// 2D axial derivative.  For (O, OA=2, D=2) the dictionary covers
// S2_symbol<i,j> and is applied with stride[2].
// --------------------------------------------------------------------------

inline double poly2d(double x, double y,
                     const double (&c)[5][5]) {
    double yval = 0.0;
    for (int kx = 0; kx < 5; ++kx) {
        double xp = std::pow(x, kx);
        for (int ky = 0; ky < 5; ++ky) {
            yval += c[kx][ky] * xp * std::pow(y, ky);
        }
    }
    return yval;
}

inline double poly2d_dxn_dyn(int Ox, int Oy, double x, double y,
                              const double (&c)[5][5]) {
    double yval = 0.0;
    for (int kx = Ox; kx < 5; ++kx) {
        double fx = 1.0;
        for (int j = 0; j < Ox; ++j) fx *= (kx - j);
        for (int ky = Oy; ky < 5; ++ky) {
            double fy = 1.0;
            for (int j = 0; j < Oy; ++j) fy *= (ky - j);
            yval += fx * fy * c[kx][ky] *
                    std::pow(x, kx - Ox) * std::pow(y, ky - Oy);
        }
    }
    return yval;
}

template <size_t O>
void test_2d_axial_derivative(double h, const double (&c)[5][5],
                              const char* tag) {
    auto dict = expr::get_central_space_stencil<O, 2, 2>();
    using gen_t = symphas::internal::GeneratedStencilApply<decltype(dict)>;

    constexpr int N = 9;
    double buf[N * N] = {0};
    const int c0 = N / 2;
    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < N; ++i) {
            const double x = (i - c0) * h;
            const double y = (j - c0) * h;
            buf[j * N + i] = poly2d(x, y, c);
        }
    }
    const double divh = 1.0 / h;
    const int center = c0 * N + c0;
    const len_type stride[2] = {1, N};

    char buftag[160];
    std::snprintf(buftag, sizeof(buftag), "%s d^%zu (2D, OA=2)", tag, O);
    if constexpr (!gen_invokable<gen_t, len_type[2]>::value) {
        ++g_subtests; ++g_failures;
        std::fprintf(stderr, "  [FAIL] %-60s  Symbol leak (dictionary not reduced)\n", buftag);
        return;
    } else {
    const double got = gen_t{}(buf + center, stride, divh);
    // get_central_space_stencil<O,2,2> is "axial-X" because of the
    // builder's chosen d_op: for O==1, d/dx; for even O it's pow<O/2>(d2x+d2y)
    // which is the laplacian-based operator -> not strictly d^O/dx^O.
    // We special-case below.
    double want;
    if (O == 1) {
        want = poly2d_dxn_dyn(1, 0, 0, 0, c);
    } else if (O % 2 == 0) {
        // pow<O/2>(d2x + d2y) of polynomial at origin
        // Apply (d2x + d2y) O/2 times to poly.
        // For 4th-order poly we can compute by formula for small O.
        // Use a numerical reference: build it analytically by recursion.
        double work[5][5];
        for (int i = 0; i < 5; ++i) for (int j = 0; j < 5; ++j) work[i][j] = c[i][j];
        for (size_t step = 0; step < O / 2; ++step) {
            double next[5][5] = {{0}};
            for (int kx = 0; kx < 5; ++kx) for (int ky = 0; ky < 5; ++ky) {
                if (kx >= 2) next[kx - 2][ky] += kx * (kx - 1) * work[kx][ky];
                if (ky >= 2) next[kx][ky - 2] += ky * (ky - 1) * work[kx][ky];
            }
            for (int i = 0; i < 5; ++i) for (int j = 0; j < 5; ++j) work[i][j] = next[i][j];
        }
        want = work[0][0];
    } else {
        // odd O > 1: d_op = d/dx * pow<O/2>(d2x + d2y).  For the
        // polynomial bases we use (max degree 4), the maximum relevant O
        // is 3.  For O=3, this is d/dx * (d2x + d2y) of the poly.
        double work[5][5];
        for (int i = 0; i < 5; ++i) for (int j = 0; j < 5; ++j) work[i][j] = c[i][j];
        for (size_t step = 0; step < O / 2; ++step) {
            double next[5][5] = {{0}};
            for (int kx = 0; kx < 5; ++kx) for (int ky = 0; ky < 5; ++ky) {
                if (kx >= 2) next[kx - 2][ky] += kx * (kx - 1) * work[kx][ky];
                if (ky >= 2) next[kx][ky - 2] += ky * (ky - 1) * work[kx][ky];
            }
            for (int i = 0; i < 5; ++i) for (int j = 0; j < 5; ++j) work[i][j] = next[i][j];
        }
        // one final d/dx
        double next[5][5] = {{0}};
        for (int kx = 1; kx < 5; ++kx) for (int ky = 0; ky < 5; ++ky)
            next[kx - 1][ky] += kx * work[kx][ky];
        want = next[0][0];
    }

    char buftag[160];
    std::snprintf(buftag, sizeof(buftag), "%s d^%zu (2D, OA=2)", tag, O);
    check_close(got, want, 1e-9, buftag);
    }
}

// --------------------------------------------------------------------------
// 2D mixed derivative.  For (O1, O2) with OA=2, D=2.
// --------------------------------------------------------------------------

template <size_t O1, size_t O2>
void test_2d_mixed_derivative(double h, const double (&c)[5][5],
                              const char* tag) {
    auto dict = expr::get_central_space_mixed_stencil<2>(
        std::index_sequence<O1, O2>{});
    using gen_t = symphas::internal::GeneratedStencilApply<decltype(dict)>;

    constexpr int N = 9;
    double buf[N * N] = {0};
    const int c0 = N / 2;
    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < N; ++i) {
            const double x = (i - c0) * h;
            const double y = (j - c0) * h;
            buf[j * N + i] = poly2d(x, y, c);
        }
    }
    const double divh = 1.0 / h;
    const int center = c0 * N + c0;
    const len_type stride[2] = {1, N};

    char buftag[160];
    std::snprintf(buftag, sizeof(buftag),
                  "%s d^%zu/dx^%zu d^%zu/dy^%zu (2D mixed)",
                  tag, O1, O1, O2, O2);
    if constexpr (!gen_invokable<gen_t, len_type[2]>::value) {
        ++g_subtests; ++g_failures;
        std::fprintf(stderr, "  [FAIL] %-60s  Symbol leak (dictionary not reduced)\n", buftag);
    } else {
        const double got = gen_t{}(buf + center, stride, divh);
        const double want = poly2d_dxn_dyn(O1, O2, 0, 0, c);
        check_close(got, want, 1e-9, buftag);
    }
}

}  // namespace

int testautogenstencil() {
    g_failures = 0;
    g_subtests = 0;
    std::fprintf(stderr, "--- testautogenstencil: autogen reduction regression ---\n");

    const double h = 0.5;

    // For each order O, use a polynomial of exact degree O so the
    // OA=2 central-difference stencil reproduces the derivative with
    // zero truncation error.

    {
        // d/dx: use p(x) = x  -> p'(0) = 1
        const double c[5] = {0.0, 1.0, 0.0, 0.0, 0.0};
        test_1d_axial_derivative<1>(h, c, "p(x)=x");
    }
    {
        // d²/dx²: use p(x) = x²/2  -> p''(0) = 1
        const double c[5] = {0.0, 0.0, 0.5, 0.0, 0.0};
        test_1d_axial_derivative<2>(h, c, "p(x)=x^2/2");
    }
    {
        // d³/dx³: use p(x) = x³/6  -> p'''(0) = 1.  OA=2 d³ stencil
        // has truncation O(h²) so leave one extra degree in the poly.
        const double c[5] = {0.0, 0.0, 0.0, 1.0 / 6.0, 0.0};
        test_1d_axial_derivative<3>(h, c, "p(x)=x^3/6");
    }
    {
        // d⁴/dx⁴: use p(x) = x⁴/24  -> p''''(0) = 1
        const double c[5] = {0.0, 0.0, 0.0, 0.0, 1.0 / 24.0};
        test_1d_axial_derivative<4>(h, c, "p(x)=x^4/24");
    }

    // 2D axial derivatives.  Each polynomial chosen so the OA=2 stencil
    // is exact (degree matches operator order).
    {
        double c2d[5][5] = {{0}};
        // p(x,y) = x  -> d/dx p = 1
        c2d[1][0] = 1.0;
        test_2d_axial_derivative<1>(h, c2d, "p=x");
    }
    {
        double c2d[5][5] = {{0}};
        // For O=2 even, builder uses pow<1>(d2x+d2y).
        // p = x²/2  -> lap p = 1
        c2d[2][0] = 0.5;
        test_2d_axial_derivative<2>(h, c2d, "p=x^2/2");
    }
    {
        double c2d[5][5] = {{0}};
        // O=3: builder uses d/dx * pow<1>(d2x+d2y) so
        // p = x³/6 -> (d/dx)(d²p/dx² + d²p/dy²) = (d/dx)(x) = 1
        c2d[3][0] = 1.0 / 6.0;
        test_2d_axial_derivative<3>(h, c2d, "p=x^3/6");
    }
    {
        double c2d[5][5] = {{0}};
        // O=4 even: builder uses pow<2>(d2x+d2y) = bilaplacian.
        // For p = x⁴/24: lap p = x² /2, lap(lap p) = 1.
        c2d[4][0] = 1.0 / 24.0;
        test_2d_axial_derivative<4>(h, c2d, "p=x^4/24");
    }

    // 2D mixed derivatives.
    {
        double c[5][5] = {{0}};
        c[1][1] = 1.0;          // p = xy -> d²/(dx dy) = 1
        test_2d_mixed_derivative<1, 1>(h, c, "p=xy");
    }
    {
        double c[5][5] = {{0}};
        c[2][1] = 0.5;          // p = x²y/2 -> d²/dx² d/dy = 1
        test_2d_mixed_derivative<2, 1>(h, c, "p=x^2y/2");
    }
    {
        double c[5][5] = {{0}};
        c[1][2] = 0.5;          // p = xy²/2 -> d/dx d²/dy² = 1
        test_2d_mixed_derivative<1, 2>(h, c, "p=xy^2/2");
    }
    {
        double c[5][5] = {{0}};
        c[2][2] = 0.25;         // p = x²y²/4 -> d²/dx² d²/dy² = 1
        test_2d_mixed_derivative<2, 2>(h, c, "p=x^2y^2/4");
    }
    {
        double c[5][5] = {{0}};
        c[3][0] = 1.0 / 6.0;
        test_2d_mixed_derivative<3, 0>(h, c, "p=x^3/6");
    }
    {
        double c[5][5] = {{0}};
        c[0][3] = 1.0 / 6.0;
        test_2d_mixed_derivative<0, 3>(h, c, "p=y^3/6");
    }
    {
        double c[5][5] = {{0}};
        c[3][1] = 1.0 / 6.0;
        test_2d_mixed_derivative<3, 1>(h, c, "p=x^3y/6");
    }
    {
        double c[5][5] = {{0}};
        c[1][3] = 1.0 / 6.0;
        test_2d_mixed_derivative<1, 3>(h, c, "p=xy^3/6");
    }

    std::fprintf(stderr,
                 "--- testautogenstencil: %d failures out of %d subtests ---\n",
                 g_failures, g_subtests);
    return g_failures;
}
