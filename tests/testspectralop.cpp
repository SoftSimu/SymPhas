#include "testspectralop.h"

#include "expressiontypeincludes.h"
#include "modelmacros.h"
#include "modelarray.h"
#include "solverinclude.h"
#include "stencilincludes.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <unistd.h>

// solverinclude.h already pulls in solversp.h / spectrallib.h when
// USING_FFTW; including them again here would double-register the solver.

namespace {

#ifdef USING_FFTW

template <typename E>
std::string capture_printe(E const& e, const char* label) {
  fflush(stdout);
  char buf[8192] = {0};
  int saved = dup(fileno(stdout));
  FILE* tmp = tmpfile();
  dup2(fileno(tmp), fileno(stdout));
  expr::printe(e, label);
  fflush(stdout);
  fseek(tmp, 0, SEEK_SET);
  size_t n = fread(buf, 1, sizeof(buf) - 1, tmp);
  buf[n] = '\0';
  fclose(tmp);
  dup2(saved, fileno(stdout));
  close(saved);
  return std::string(buf);
}

// Compute (i, j) wave-vector components the same way k_field_entry does.
// Used to predict the closed-form L(k) we expect for Cahn-Hilliard.
struct kpoint {
  int n;          // linear grid index
  double kx, ky;  // physical wave numbers
};

kpoint kindex(int i, int j, const len_type dims[2], const double h[2]) {
  len_type L = dims[0], M = dims[1];
  double dk_i = 2 * symphas::PI / (h[0] * L);
  double dk_j = 2 * symphas::PI / (h[1] * M);
  double kx = (i < L / 2) ? i * dk_i : (i - L) * dk_i;
  double ky = (j < M / 2) ? j * dk_j : (j - M) * dk_j;
  return {i + j * (int)L, kx, ky};
}

#endif  // USING_FFTW

}  // namespace

void testspectralop() {
  fprintf(stdout, "\n========== SP linear-operator regression ==========\n");

#ifndef USING_FFTW
  fprintf(stdout, "  (USING_FFTW not defined -- skipped)\n");
  return;
#else

  constexpr size_t D = 2;
  constexpr size_t Z = 0;
  constexpr len_type N = 32;
  len_type dims[D] = {N, N};
  double h[D] = {1.0, 1.0};
  double dt = 0.005;

  // Construct symbolic Model B RHS exactly the way the model macros would.
  using Sp = SolverFT<Stencil2d2h<9, 6, 13>>;
  Sp dummy_solver{dims, 1.0};

  BoundaryGrid<scalar_t, D> psi_grid(dims);
  auto psi = expr::make_term<Z>(symphas::ref<BoundaryGrid<scalar_t, D>>(psi_grid));
  auto lap_op = expr::make_operator_derivative<2>(dummy_solver);
  auto bilap_op = expr::make_operator_derivative<4>(dummy_solver);

  // Model B: dpsi/dt = -bilap(psi) - lap(psi) + lap(psi^3)
  //                  = -V^4 psi - V^2 psi - V^2 (-psi^3)
  auto rhs = expr::apply_operators(
      -bilap_op(psi) - lap_op(psi) + lap_op(psi * psi * psi));

  std::string s_rhs = capture_printe(rhs, "given equation");
  fprintf(stdout, "[STEP 0] symbolic RHS\n  %s", s_rhs.c_str());

  // Split linear / nonlinear in the spectral pipeline. Same calls as
  // SolverSP::form_expr_one (examples/solvers/solversp.h:88-94).
  auto&& [linear, nonlinear] = expr::split::by_linear(rhs);
  auto&& [lin_Z, lin_nonZ] = expr::split::separate_var<Z>(linear);
  auto&& [l_op, non_op] = solver_sp::get_l_op<Z>(lin_Z, h);

  std::string s_lop = capture_printe(l_op, "l_op (linear symbol)");
  fprintf(stdout, "[STEP 1] symbolic l_op (Fourier symbol of linear part)\n");
  fprintf(stdout, "  %s", s_lop.c_str());

  auto A_expr = solver_sp::form_A_op<D>(l_op, dt, dims);
  auto B_expr = solver_sp::form_B_op<D>(l_op, dt, dims);

  std::string s_A = capture_printe(A_expr, "A_expr");
  std::string s_B = capture_printe(B_expr, "B_expr");
  fprintf(stdout, "[STEP 2] A = exp(dt * l_op)\n  %s", s_A.c_str());
  fprintf(stdout, "[STEP 3] B = (A - 1) / l_op\n  %s", s_B.c_str());

  // --- Numeric probe ----------------------------------------------------
  // Evaluate l_op AT the k-grid indices and compare to the closed-form
  // Cahn-Hilliard symbol L(k) = -|k|^4 + |k|^2 (Fourier symbol of
  // -bilap - lap, where bilap -> +|k|^4 and lap -> -|k|^2).
  fprintf(stdout, "[STEP 4] numeric l_op(n) vs closed-form L(k)\n");
  fprintf(stdout,
          "  %-6s %-10s %-10s %-15s %-15s %-15s\n",
          "(i,j)", "kx", "ky", "l_op.eval(n)", "expected L(k)", "delta");

  int fail = 0;
  int probes[][2] = {{0, 0}, {1, 0}, {0, 1}, {2, 2}, {4, 0},
                      {8, 8}, {N / 2, 0}, {N / 2, N / 2}, {1, N / 2}};
  for (auto& p : probes) {
    int i = p[0], j = p[1];
    auto k = kindex(i, j, dims, h);
    // k_field_entry substitutes a different epsilon per derivative order
    // when a component is exactly zero, so compute k^2 and k^4
    // separately, mirroring that substitution. For O=2 the zero-axis
    // floor is pow(EPS, 1/2); for O=4 it is pow(EPS, 1/4).
    double kx2 = (k.kx == 0) ? std::pow(symphas::EPS, 1.0 / 2.0) : k.kx;
    double ky2 = (k.ky == 0) ? std::pow(symphas::EPS, 1.0 / 2.0) : k.ky;
    double kx4 = (k.kx == 0) ? std::pow(symphas::EPS, 1.0 / 4.0) : k.kx;
    double ky4 = (k.ky == 0) ? std::pow(symphas::EPS, 1.0 / 4.0) : k.ky;
    double k2_o2 = kx2 * kx2 + ky2 * ky2;
    double k2_o4 = kx4 * kx4 + ky4 * ky4;
    double expected = -k2_o4 * k2_o4 + k2_o2;  // -|k|^4 + |k|^2
    auto val = l_op.eval(k.n);
    double got = static_cast<double>(val);
    double delta = got - expected;
    fprintf(stdout, "  (%2d,%2d) %-10.4f %-10.4f %-15.6f %-15.6f %-15.3e%s\n",
            i, j, k.kx, k.ky, got, expected, delta,
            (std::fabs(delta) > 1e-9 * (1.0 + std::fabs(expected))) ? "  FAIL"
                                                                   : "");
    if (std::fabs(delta) > 1e-9 * (1.0 + std::fabs(expected))) fail++;
  }

  if (fail == 0) {
    fprintf(stdout,
            "\nSP-regression: l_op numeric values MATCH closed-form L(k).\n"
            "  Conclusion: printed symbol is misleading but algebra is correct.\n");
  } else {
    fprintf(stderr,
            "\nSP-regression: %d numeric mismatches -- real algebra bug.\n",
            fail);
    std::exit(1);
  }

  // --- Full scheme construction --------------------------------------------
  // Build the exact same nonlinear_scheme that solversp.h does, dump its
  // printed form and numerically check the coefficient applied to F{psi^3}.
  // For Model B with linear L(k) = -|k|^4 + |k|^2 and nonlinear term
  // lap(-psi^3), the ETD1 nonlinear coefficient should be
  //   C(k) = (A(k) - 1) / L(k) * (-|k|^2)
  // multiplied onto F{psi^3} (note the leading + because lap(-psi^3) carries
  // its own -1, so the F{psi^3} coefficient is +(A-1)/L * (-|k|^2) * -1
  //                                        = (A-1)/L * |k|^2).
  fprintf(stdout, "\n[STEP 5] full nonlinear-scheme construction\n");

  // Replicate the linear_in_nonZ + non_op + nonlinear sum that
  // form_expr_one feeds into construct_nonlinear.
  auto nonlin_input = lin_nonZ + non_op + nonlinear;
  std::string s_nlin = capture_printe(nonlin_input, "nonlin_input");
  fprintf(stdout, "  %s", s_nlin.c_str());

  // construct_nonlinear needs a "systems" tuple. Mimic with a 1-tuple
  // referencing our psi_grid (the only system).
  auto systems = std::make_tuple(symphas::ref<BoundaryGrid<scalar_t, D>>(psi_grid));
  auto nlin_scheme = solver_sp::construct_nonlinear<Z, D>(
      systems, B_expr, nonlin_input, h, dims);
  std::string s_nls = capture_printe(nlin_scheme, "nlin_scheme");
  fprintf(stdout, "  %s", s_nls.c_str());

  // Numeric evaluation of nlin_scheme is complex (it carries FFT'd data
  // and grid-cached terms via NamedData). We rely on the symbolic dump
  // above to inspect what coefficient is applied to F{psi^3}.
  fprintf(stdout, "  (nlin_scheme numeric probe skipped -- see symbolic form)\n");
  (void)probes;

  // --- Trace each substep of construct_nonlinear<OpDerivative> --------
  // This mirrors spectrallib.h:738-749 directly so we can see where the
  // |k|^2 collapses to literal 1.
  fprintf(stdout, "\n[STEP 6] manual trace of construct_nonlinear<lap(-psi^3)>\n");
  // For Model B, `nonlinear` is `lap(-psi^3)`. Extract it from the rhs.
  auto lap_neg_p3 = lap_op(-(psi * psi * psi));
  std::string s_lnp3 = capture_printe(lap_neg_p3, "lap(-psi^3)");
  fprintf(stdout, "  %s", s_lnp3.c_str());

  auto [op_inner, en_inner] = expr::split::separate_operator(lap_neg_p3);
  std::string s_op = capture_printe(op_inner, "op (lap operator)");
  std::string s_en = capture_printe(en_inner, "en (-psi^3)");
  fprintf(stdout, "  %s", s_op.c_str());
  fprintf(stdout, "  %s", s_en.c_str());

  auto to_ft_lap = expr::transform::to_ft<D>(op_inner, h, dims);
  std::string s_tft = capture_printe(to_ft_lap, "to_ft(lap)");
  fprintf(stdout, "  %s", s_tft.c_str());

  auto B_x_lap = B_expr * to_ft_lap;
  std::string s_Bxl = capture_printe(B_x_lap, "B_expr * to_ft(lap)");
  fprintf(stdout, "  %s", s_Bxl.c_str());

  // Also test the multiplication on l_op directly. If l_op * WaveVector<2>
  // collapses incorrectly, that's a symbolic-algebra bug in the
  // WaveVector*WaveVector operator overload.
  auto l_x_lap = l_op * to_ft_lap;
  std::string s_lxl = capture_printe(l_x_lap, "l_op * to_ft(lap)");
  fprintf(stdout, "  %s", s_lxl.c_str());

  // And test plain |k|^2 * |k|^2 (should be |k|^4 / 1 / something sane):
  auto k2 = expr::make_term(k_grid_type<2, D>(dims, h));
  auto k4 = expr::make_term(k_grid_type<4, D>(dims, h));
  std::string s_k2 = capture_printe(k2, "k2 = WaveVector<2>");
  std::string s_k4 = capture_printe(k4, "k4 = WaveVector<4>");
  fprintf(stdout, "  %s  %s", s_k2.c_str(), s_k4.c_str());

  auto k2_x_k2 = k2 * k2;
  std::string s_k2k2 = capture_printe(k2_x_k2, "k2 * k2");
  fprintf(stdout, "  %s", s_k2k2.c_str());

  // The killer test: (-k4 - k2) * k2 -- this is L * k2.
  // If it produces (-k6 - k4) → correct.
  // If it produces (-k4 - 1) → BUG.
  auto neg_l = -k4 - k2;  // mimic the symbolic L
  std::string s_negl = capture_printe(neg_l, "-k4 - k2");
  fprintf(stdout, "  %s", s_negl.c_str());

  auto L_times_k2 = neg_l * k2;
  std::string s_Ltk2 = capture_printe(L_times_k2, "(-k4 - k2) * k2");
  fprintf(stdout, "  %s", s_Ltk2.c_str());

  // And the reciprocal direction: k2 / (-k4 - k2).
  auto k2_div_L = k2 / neg_l;
  std::string s_div = capture_printe(k2_div_L, "k2 / (-k4 - k2)");
  fprintf(stdout, "  %s", s_div.c_str());

  // Numeric check that the printed simplification matches the actual
  // mathematical value of k^2 / (k^2 - k^4) at several k-points.
  // Closed-form value: k2/(-k4+k2) where k2 = -|k|^2 (the WaveVector
  // runtime value) and k4 = +|k|^4. So expression is
  //   -|k|^2 / (-|k|^4 + |k|^2)  = -|k|^2 / (|k|^2(1 - |k|^2))
  //                              = -1 / (1 - |k|^2)
  // (well-defined for |k| != 1).
  fprintf(stdout,
          "  numeric check k2 / (-k4 - k2)  vs  closed-form -1/(1-|k|^2):\n");
  int div_fail = 0;
  // Avoid probes where either axis is zero — `k_field_entry` substitutes
  // `pow(EPS, 1/O)` for zero components and the substitution differs
  // between O=2 and O=4, breaking a closed-form scalar comparison.
  int div_probes[][2] = {{2, 2}, {3, 5}, {8, 8}, {16, 16}};
  for (auto& p : div_probes) {
    int i = p[0], j = p[1];
    auto k = kindex(i, j, dims, h);
    double k2v = k.kx * k.kx + k.ky * k.ky;
    if (std::fabs(1.0 - k2v) < 1e-9) continue;  // skip singular
    double expected = -1.0 / (1.0 - k2v);
    double got_real = complex_t(k2_div_L.eval(k.n)).real();
    double rel = std::fabs(got_real - expected) /
                 (1.0 + std::fabs(expected));
    bool ok = rel < 1e-9;
    fprintf(stdout, "    (%2d,%2d) |k|^2=%-8.4f got=%-13.6g expected=%-13.6g  %s\n",
            i, j, k2v, got_real, expected, ok ? "OK" : "FAIL");
    if (!ok) div_fail++;
  }
  if (div_fail > 0) {
    fprintf(stderr,
            "\nDIVISION FAIL: k2/(-k4-k2) does not equal -1/(1-|k|^2) "
            "at %d points -- factor extraction is buggy.\n",
            div_fail);
    std::exit(1);
  } else {
    fprintf(stdout, "  DIVISION OK: factor extraction is sound.\n");
  }

  // Additional factor-extraction probes covering common simplification
  // patterns. Each must round-trip cleanly to the analytic value.
  fprintf(stdout, "\n[STEP 7] additional factor-extraction probes\n");

  // (k4) / (k2) should be k2 (i.e. -|k|^2). Single-term division.
  auto k4_div_k2 = k4 / k2;
  fprintf(stdout, "  k4 / k2  ->  ");
  fflush(stdout);
  expr::printe(k4_div_k2, "");
  // Numeric: k4 = |k|^4, k2 = -|k|^2, k4/k2 = -|k|^2 = k2 = -k.kx^2-k.ky^2.
  int factor_fail = 0;
  for (auto& p : div_probes) {
    auto k = kindex(p[0], p[1], dims, h);
    double k2v = k.kx * k.kx + k.ky * k.ky;
    double expected = -k2v;  // |k|^4 / -|k|^2 = -|k|^2
    double got = complex_t(k4_div_k2.eval(k.n)).real();
    double rel = std::fabs(got - expected) / (1.0 + std::fabs(expected));
    if (rel >= 1e-9) {
      fprintf(stdout, "    (%2d,%2d) k4/k2 got=%g expected=%g FAIL\n", p[0],
              p[1], got, expected);
      factor_fail++;
    }
  }

  // (k4 + k2) / k2 should be k2 + 1 (uniform factor extraction case).
  auto sum_div = (k4 + k2) / k2;
  fprintf(stdout, "  (k4 + k2) / k2  ->  ");
  fflush(stdout);
  expr::printe(sum_div, "");
  for (auto& p : div_probes) {
    auto k = kindex(p[0], p[1], dims, h);
    double k2v = k.kx * k.kx + k.ky * k.ky;
    // (k4 + k2)/k2 = (|k|^4 + -|k|^2) / -|k|^2 = (|k|^4)/(-|k|^2) + 1
    //              = -|k|^2 + 1
    double expected = -k2v + 1.0;
    double got = complex_t(sum_div.eval(k.n)).real();
    double rel = std::fabs(got - expected) / (1.0 + std::fabs(expected));
    if (rel >= 1e-9) {
      fprintf(stdout,
              "    (%2d,%2d) (k4+k2)/k2 got=%g expected=%g FAIL\n", p[0],
              p[1], got, expected);
      factor_fail++;
    }
  }

  // (k4 + 1) / k2 -- heterogeneous, factor extraction must NOT succeed
  // (k2 doesn't divide the constant 1). Result must equal the raw value
  // (|k|^4 + 1)/-|k|^2.
  auto het_div = (k4 + expr::make_literal(1.0)) / k2;
  fprintf(stdout, "  (k4 + 1) / k2  ->  ");
  fflush(stdout);
  expr::printe(het_div, "");
  for (auto& p : div_probes) {
    auto k = kindex(p[0], p[1], dims, h);
    double k2v = k.kx * k.kx + k.ky * k.ky;
    double k4v = k2v * k2v;
    double expected = (k4v + 1.0) / (-k2v);
    double got = complex_t(het_div.eval(k.n)).real();
    double rel = std::fabs(got - expected) / (1.0 + std::fabs(expected));
    if (rel >= 1e-9) {
      fprintf(stdout,
              "    (%2d,%2d) (k4+1)/k2 got=%g expected=%g FAIL\n", p[0],
              p[1], got, expected);
      factor_fail++;
    }
  }

  if (factor_fail > 0) {
    fprintf(stderr,
            "\nFACTOR PROBE FAIL: %d factor-extraction probes incorrect.\n",
            factor_fail);
    std::exit(1);
  } else {
    fprintf(stdout, "  ALL FACTOR PROBES OK\n");
  }

  // --- STEP 8: PFC-shape distribution probe -----------------------------
  // The PFC dynamic equation is constructed as
  //   ∇² * (bulk_polynomial + coupled_linear)
  // where coupled_linear = α + β(ν² + ∇²)² n.
  // If apply_operators distributes the outer ∇² across the OpAdd, the
  // resulting tree contains ∇²(∇⁴ n) which the FTCS scheme evaluates via
  // apply<6> -- a chain of three Laplacians that is *anti-diffusive at
  // Nyquist*, blowing up PFC at h=1. The "good" form keeps the outer ∇²
  // around the whole sum so the solver can choose a stable apply<2>
  // chain.
  //
  // Build a minimal stand-in:  ∇²(psi + psi^3 + lap(psi)) and inspect
  // whether apply_operators distributes it.
  fprintf(stdout, "\n[STEP 8] PFC-shape: apply_operators distribution\n");
  auto pfc_like = lap_op(psi + psi * psi * psi + lap_op(psi));
  std::string s_pfc_pre = capture_printe(pfc_like, "PFC-like pre-apply");
  fprintf(stdout, "  pre :  %s", s_pfc_pre.c_str());
  auto pfc_post = expr::apply_operators(pfc_like);
  std::string s_pfc_post = capture_printe(pfc_post, "PFC-like post-apply");
  fprintf(stdout, "  post:  %s", s_pfc_post.c_str());

  // Heuristic: count how many "V^2(" occurrences appear post-apply.
  // 1 = single outer ∇² preserved (good).
  // >1 = outer ∇² distributed (bad — produces ∇²(∇²·) chains).
  size_t pos = 0, n_lap = 0;
  while ((pos = s_pfc_post.find("V^2(", pos)) != std::string::npos) {
    ++n_lap;
    ++pos;
  }
  fprintf(stdout, "  V^2( count post-apply = %zu\n", n_lap);
  if (n_lap >= 2) {
    fprintf(stderr,
            "  PFC distribution: apply_operators DID distribute ∇² over OpAdd.\n"
            "  This produces ∇²(∇²·) chains that go through apply<4> -> chained Laplacians\n"
            "  rather than the isotropic apply<4> override. Stable for ∇² and ∇⁴ but ∇²(∇⁴·) -> apply<6>\n"
            "  fires the anti-diffusive triple-Laplacian path at Nyquist.\n");
  } else {
    fprintf(stdout, "  PFC distribution: outer ∇² preserved (good).\n");
  }

  // --- STEP 9: numeric evaluation of each PFC term ----------------------
  // Build expressions matching what the PFC model prints, evaluate each
  // term ON a known plane-wave field, and compare to the analytic
  // continuum derivative.
  //
  // We fill psi_grid with sin(kx_p * x) (cell-centered: x = i*h) for
  // probe wavenumbers kx_p. Then for each "term":
  //   lap(psi)     -> continuum -kx^2 * sin(kx*x)
  //   bilap(psi)   -> continuum +kx^4 * sin(kx*x)
  //   hexlap(psi)  -> continuum -kx^6 * sin(kx*x)
  //   lap(lap(psi)) (chained)
  //   lap(lap(lap(psi))) (chained, what apply_operators may leave PFC as)
  //
  // We then check value at a single point against analytic expectation,
  // and also compute the L2 norm of (numerical - analytic) / |analytic|.
  fprintf(stdout, "\n[STEP 9] numeric PFC term evaluation on plane wave\n");

  // Choose a few probe wavenumbers, including low-k, mid-k, and Nyquist.
  struct kprobe { double k; const char* tag; };
  kprobe probes_k[] = {
      {2.0 * symphas::PI / (N * h[0]) * 4,  "low-k (4 cycles)"},
      {2.0 * symphas::PI / (N * h[0]) * 16, "mid-k (16 cycles)"},
      {2.0 * symphas::PI / (N * h[0]) * 32, "high-k (32 cycles = N/4)"},
      {symphas::PI / h[0],                  "Nyquist (k = pi/h)"},
  };

  // Build the expressions. Distribute_operator suppression: build each
  // applied directly without an outer OpAdd to keep tree small.
  auto term_lap     = lap_op(psi);
  auto term_bilap   = bilap_op(psi);
  // Single hexalap via two-stage: use lap of bilap, which apply_operators
  // should fold into derivative order 6.
  auto term_hexlap_via_bb = expr::apply_operators(lap_op(bilap_op(psi)));
  std::string s_h = capture_printe(term_hexlap_via_bb, "lap(bilap(psi)) post-apply");
  fprintf(stdout, "  %s", s_h.c_str());
  // Chained lap-of-lap, lap-of-lap-of-lap WITHOUT apply_operators
  // (mimics what would happen if symbolic algebra never folded the
  // chain into a single OpDerivative<order=4 or 6>).
  // But OpOperatorDerivative * OpOperatorDerivative is currently composed
  // into OpOperatorDerivative<sum_order>, so we can't easily build the
  // un-folded chain. Instead, evaluate the printed-PFC term 4 expression:
  //   lap(lap(1+lap)(psi)) = lap(lap(psi + lap(psi)))
  auto term4_raw = lap_op(lap_op(psi + lap_op(psi)));
  std::string s_t4_pre  = capture_printe(term4_raw, "term4 pre-apply");
  fprintf(stdout, "  %s", s_t4_pre.c_str());
  auto term4_applied = expr::apply_operators(term4_raw);
  std::string s_t4_post = capture_printe(term4_applied, "term4 post-apply");
  fprintf(stdout, "  %s", s_t4_post.c_str());

  // Evaluate each on a plane wave. BoundaryGrid stores a (N+2*BD)² buffer
  // with BOUNDARY_DEPTH halo cells on every side. We fill the entire raw
  // buffer with sin(k * x_raw) so the stencil-evaluator sees a self-
  // consistent plane wave including the halo.
  // BoundaryGrid here is constructed with dims = {N, N} and has no halo
  // (verified: psi_grid.len = N*N). Stride for the stencil is dims[0] = N.
  // Sample point must stay ≥ 3 cells from any edge so the 7-wide hexlap
  // stencil doesn't read out of bounds.
  fprintf(stdout, "  BOUNDARY_DEPTH=%d, N=%d, psi_grid.len=%d\n",
          (int)BOUNDARY_DEPTH, (int)N, (int)psi_grid.len);
  auto fill_plane_wave = [&](double k) {
    for (iter_type n = 0; n < N * N; ++n) {
      int i = n % N;
      double x = i * h[0];
      psi_grid[n] = std::cos(k * x);
    }
  };
  // Fill the 2D checkerboard mode  (-1)^(i+j) = cos(pi*i)*cos(pi*j).
  auto fill_2d_nyquist = [&]() {
    for (iter_type n = 0; n < N * N; ++n) {
      int i = n % N;
      int j = n / N;
      psi_grid[n] = ((i + j) & 1) ? -1.0 : 1.0;
    }
  };
  auto interior_idx = [&](int i, int j) -> iter_type {
    return i + j * N;
  };

  for (auto& kp : probes_k) {
    double k = kp.k;
    fill_plane_wave(k);
    int i0 = N / 2 + 3;  // off-axis to avoid zeros of cos
    int j0 = N / 2 + 1;
    double x = i0 * h[0];
    double psi_val = std::cos(k * x);
    // Sanity: check the grid value at the sample point matches cos(k*x)
    iter_type p = interior_idx(i0, j0);
    double grid_at_p = static_cast<double>(psi_grid[p]);
    if (std::fabs(grid_at_p - psi_val) > 1e-10) {
      fprintf(stdout,
              "    [indexing mismatch] grid[p]=%+.6g  cos(k*x)=%+.6g\n",
              grid_at_p, psi_val);
    }
    // Continuum: lap[cos(kx)] = -k^2 cos(kx), bilap = +k^4 cos, hexlap = -k^6 cos
    double cont_lap   = -k*k * psi_val;
    double cont_bilap =  k*k*k*k * psi_val;
    double cont_hex   = -k*k*k*k*k*k * psi_val;

    double n_lap   = static_cast<double>(term_lap.eval(p));
    double n_bilap = static_cast<double>(term_bilap.eval(p));
    double n_hex   = static_cast<double>(term_hexlap_via_bb.eval(p));
    double n_t4    = static_cast<double>(term4_applied.eval(p));

    fprintf(stdout,
            "  %-26s k=%.4f   sin=%+.4f\n",
            kp.tag, k, psi_val);
    fprintf(stdout,
            "    lap     got=%+13.6g cont=%+13.6g  rel=%+.2e\n",
            n_lap, cont_lap,
            (cont_lap == 0) ? 0 : (n_lap - cont_lap) / std::fabs(cont_lap));
    fprintf(stdout,
            "    bilap   got=%+13.6g cont=%+13.6g  rel=%+.2e\n",
            n_bilap, cont_bilap,
            (cont_bilap == 0) ? 0 : (n_bilap - cont_bilap) / std::fabs(cont_bilap));
    fprintf(stdout,
            "    hexlap  got=%+13.6g cont=%+13.6g  rel=%+.2e\n",
            n_hex, cont_hex,
            (cont_hex == 0) ? 0 : (n_hex - cont_hex) / std::fabs(cont_hex));
    // Term 4 analytic: lap(lap(psi + lap(psi))) on plane wave
    //   = lap(lap(sin(kx) - k^2 sin(kx)))
    //   = lap(lap(sin(kx)(1 - k^2)))
    //   = (1 - k^2) * lap(lap(sin(kx)))
    //   = (1 - k^2) * k^4 * sin(kx)
    double cont_t4 = (1.0 - k*k) * k*k*k*k * psi_val;
    fprintf(stdout,
            "    term4   got=%+13.6g cont=%+13.6g  rel=%+.2e\n",
            n_t4, cont_t4,
            (cont_t4 == 0) ? 0 : (n_t4 - cont_t4) / std::fabs(cont_t4));
  }

  // 2D corner-Nyquist probe: checkerboard (-1)^(i+j). Eigenvalues:
  //   continuum: lap=-(pi^2+pi^2)=-2pi^2, bilap=+4pi^4, hexlap=-8pi^6
  //   So on field=+1 at (i+j) even: cont_lap=-2pi^2, cont_bilap=+4pi^4, cont_hex=-8pi^6.
  fprintf(stdout, "  2D corner-Nyquist checkerboard (kx=ky=pi):\n");
  fill_2d_nyquist();
  int ic = N / 2 + 2;  // i+j even => field=+1
  int jc = N / 2 + 2;
  iter_type pc = interior_idx(ic, jc);
  double fc = static_cast<double>(psi_grid[pc]);
  double pi = symphas::PI;
  double clap   = -2*pi*pi * fc;
  double cbilap =  4*pi*pi*pi*pi * fc;
  double chex   = -8*pi*pi*pi*pi*pi*pi * fc;
  double nlap   = static_cast<double>(term_lap.eval(pc));
  double nbilap = static_cast<double>(term_bilap.eval(pc));
  double nhex   = static_cast<double>(term_hexlap_via_bb.eval(pc));
  double nt4    = static_cast<double>(term4_applied.eval(pc));
  fprintf(stdout,
          "    field=%+.1f   lap  got=%+9.4g cont=%+9.4g    bilap got=%+9.4g cont=%+9.4g\n"
          "                   hex  got=%+9.4g cont=%+9.4g    term4 got=%+9.4g cont=%+9.4g\n",
          fc, nlap, clap, nbilap, cbilap, nhex, chex, nt4,
          (1.0 - 1.0) * 4*pi*pi*pi*pi * fc + chex);  // term4 cont = (1-k^2)k^4·fc + (-k^6·fc) ... but k=(pi,pi) so |k|^2=2pi^2; term4=lap(lap(psi+lap(psi))) on 2D nyquist = lap(lap((1-2pi^2)·fc)) = (1-2pi^2)·(2pi^2)^2 ·fc
  // Better: compute analytically:
  //   psi = fc
  //   lap psi   = -(|k|^2) psi
  //   psi+lap psi = (1-(|k|^2)) psi
  //   lap(psi+lap psi) = -(|k|^2)(1-(|k|^2)) psi
  //   lap(lap(...)) = (|k|^2)^2 (1-(|k|^2)) psi
  double k2_2d = 2*pi*pi;
  double t4_cont_2d = k2_2d*k2_2d * (1 - k2_2d) * fc;
  fprintf(stdout,
          "    [analytic] term4 cont (2D corner) = %+9.4g\n", t4_cont_2d);

  // --- STEP 10: chained-Laplacian dispatch test ---------------------
  // Apply the 9-pt isotropic Laplacian THREE times sequentially via
  // explicit grids: g1 = lap(psi), g2 = lap(g1), g3 = lap(g2). Compare
  // g3 at the 2D-corner Nyquist mode against:
  //   - the unified apply<6> result (already shown above as `hex`)
  //   - the continuum -|k|^6 value
  // The hypothesis we are testing: does the working binary's chained-
  // Laplacian path produce a DIFFERENT eigenvalue at corner Nyquist than
  // the unified 29-pt apply<6> override?
  fprintf(stdout, "\n[STEP 10] chained apply<2>(apply<2>(apply<2>(.))) on checkerboard\n");
  fill_2d_nyquist();

  // Build two temp grids same dims; will use term_lap on g_in to fill g_out.
  BoundaryGrid<scalar_t, D> g1(dims), g2(dims), g3(dims);
  auto psi_ref = expr::make_term<0>(symphas::ref<BoundaryGrid<scalar_t, D>>(psi_grid));
  auto g1_ref  = expr::make_term<1>(symphas::ref<BoundaryGrid<scalar_t, D>>(g1));
  auto g2_ref  = expr::make_term<2>(symphas::ref<BoundaryGrid<scalar_t, D>>(g2));
  auto lap_psi_expr = expr::make_operator_derivative<2>(dummy_solver)(psi_ref);
  auto lap_g1_expr  = expr::make_operator_derivative<2>(dummy_solver)(g1_ref);
  auto lap_g2_expr  = expr::make_operator_derivative<2>(dummy_solver)(g2_ref);

  // Materialize g1 = lap(psi) at every interior point [3..N-3]
  for (int j = 3; j < N - 3; ++j)
    for (int i = 3; i < N - 3; ++i)
      g1[i + j * N] = static_cast<double>(lap_psi_expr.eval(i + j * N));
  // g2 = lap(g1)
  for (int j = 3; j < N - 3; ++j)
    for (int i = 3; i < N - 3; ++i)
      g2[i + j * N] = static_cast<double>(lap_g1_expr.eval(i + j * N));
  // g3 = lap(g2)
  for (int j = 3; j < N - 3; ++j)
    for (int i = 3; i < N - 3; ++i)
      g3[i + j * N] = static_cast<double>(lap_g2_expr.eval(i + j * N));

  // Sample well inside (avoid the 3-cell unfilled ring on every side
  // since we only updated [3..N-3] and the lap stencil reads ±1).
  int ic2 = N / 2;
  int jc2 = N / 2;
  double fc2 = static_cast<double>(psi_grid[ic2 + jc2 * N]);
  double g1v = static_cast<double>(g1[ic2 + jc2 * N]);
  double g2v = static_cast<double>(g2[ic2 + jc2 * N]);
  double g3v = static_cast<double>(g3[ic2 + jc2 * N]);
  double hex_unified = static_cast<double>(term_hexlap_via_bb.eval(ic2 + jc2 * N));
  fprintf(stdout,
          "  field=%+.1f  lap1=%+9.4g  lap2=%+9.4g  lap3=%+9.4g     hex_unified(V^6)=%+9.4g\n",
          fc2, g1v, g2v, g3v, hex_unified);
  fprintf(stdout,
          "  ratio lap3/field = %+9.4g   (this is the chained-Lap eigenvalue at corner Nyquist)\n",
          fc2 == 0 ? 0 : g3v / fc2);
  fprintf(stdout,
          "  ratio hex_unified/field = %+9.4g\n",
          fc2 == 0 ? 0 : hex_unified / fc2);
  fprintf(stdout, "  continuum -|k|^6 at corner = %+9.4g\n", -8.0 * std::pow(symphas::PI, 6));

  // --- Negative control -------------------------------------------------
  // Build a deliberately WRONG RHS where lap(psi) has the opposite sign
  // (-bilap(psi) + lap(psi) instead of -bilap(psi) - lap(psi)) and verify
  // that l_op.eval() now DIFFERS from the correct Cahn-Hilliard symbol.
  // If the test still "passes" against the original reference, then it
  // would be insensitive to a real algebra bug.
  fprintf(stdout, "\n[NEG CTRL] same probe with WRONG RHS (-bilap + lap)\n");
  auto rhs_wrong = expr::apply_operators(
      -bilap_op(psi) + lap_op(psi) + lap_op(psi * psi * psi));
  auto&& [lin_w, nl_w] = expr::split::by_linear(rhs_wrong);
  auto&& [lin_w_Z, lin_w_nonZ] = expr::split::separate_var<Z>(lin_w);
  auto&& [l_op_w, non_op_w] = solver_sp::get_l_op<Z>(lin_w_Z, h);
  std::string s_lop_w = capture_printe(l_op_w, "wrong l_op");
  fprintf(stdout, "  %s", s_lop_w.c_str());

  int neg_disagree = 0;
  for (auto& p : probes) {
    int i = p[0], j = p[1];
    auto k = kindex(i, j, dims, h);
    double kx2 = (k.kx == 0) ? std::pow(symphas::EPS, 1.0 / 2.0) : k.kx;
    double ky2 = (k.ky == 0) ? std::pow(symphas::EPS, 1.0 / 2.0) : k.ky;
    double kx4 = (k.kx == 0) ? std::pow(symphas::EPS, 1.0 / 4.0) : k.kx;
    double ky4 = (k.ky == 0) ? std::pow(symphas::EPS, 1.0 / 4.0) : k.ky;
    double k2_o2 = kx2 * kx2 + ky2 * ky2;
    double k2_o4 = kx4 * kx4 + ky4 * ky4;
    double correct_L = -k2_o4 * k2_o4 + k2_o2;  // what MB *should* give
    double got = static_cast<double>(l_op_w.eval(k.n));
    double delta = got - correct_L;
    if (std::fabs(delta) > 1e-9 * (1.0 + std::fabs(correct_L))) neg_disagree++;
    fprintf(stdout,
            "  (%2d,%2d) got=%-13.6f correct_MB=%-13.6f delta=%-12.3e%s\n", i,
            j, got, correct_L, delta, (std::fabs(delta) > 1e-9) ? "  diff" : "");
  }
  if (neg_disagree > 0) {
    fprintf(stdout,
            "\nNEG CTRL OK: %d/%lu probes disagree with correct symbol --\n"
            "  the test IS sensitive to a real algebra bug.\n",
            neg_disagree, sizeof(probes) / sizeof(probes[0]));
  } else {
    fprintf(stderr,
            "\nNEG CTRL FAIL: wrong RHS gave same numeric l_op as correct --\n"
            "  the test would not catch a real algebra bug.\n");
    std::exit(1);
  }

#endif  // USING_FFTW
}
