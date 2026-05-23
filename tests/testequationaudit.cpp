// Audit printed-equation form for suspicious tokens.  See testequationaudit.h.

#include "testequationaudit.h"

#include "symphas.h"

#define MODEL_SET_1
#define PoissonSolver(E) expr::poisson_solver(E)

#include "modeldefinitions.h"
#include "modelinclude.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <regex>
#include <string>
#include <vector>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

namespace {

struct StdoutCapture {
  int saved_stdout = -1;
  int tmp_fd = -1;
  std::string path = "/tmp/sym_audit_capture.XXXXXX";
  void start() {
    std::fflush(stdout);
    saved_stdout = ::dup(fileno(stdout));
    std::vector<char> buf(path.begin(), path.end());
    buf.push_back('\0');
    tmp_fd = ::mkstemp(buf.data());
    path.assign(buf.data());
    ::dup2(tmp_fd, fileno(stdout));
  }
  std::string stop() {
    std::fflush(stdout);
    ::dup2(saved_stdout, fileno(stdout));
    ::close(saved_stdout);
    ::close(tmp_fd);
    FILE* f = std::fopen(path.c_str(), "r");
    std::string out;
    if (f) {
      char buf[4096];
      size_t n;
      while ((n = std::fread(buf, 1, sizeof(buf), f)) > 0)
        out.append(buf, n);
      std::fclose(f);
    }
    ::unlink(path.c_str());
    return out;
  }
};

std::vector<std::string> extract_equations(std::string const& captured) {
  std::vector<std::string> eqs;
  std::regex re("given equation\\s*->\\s*([^\\n]+)");
  auto begin = std::sregex_iterator(captured.begin(), captured.end(), re);
  auto end = std::sregex_iterator();
  for (auto it = begin; it != end; ++it) {
    std::string line = (*it)[1].str();
    while (!line.empty() &&
           std::isspace(static_cast<unsigned char>(line.back())))
      line.pop_back();
    eqs.push_back(line);
  }
  return eqs;
}

static int g_pass = 0;
static int g_fail = 0;

// Build model, capture stdout, return printed equations.
template <typename ModelT>
std::vector<std::string> instantiate_and_capture(size_t num_fields,
                                                 double const* coeffs,
                                                 size_t n_coeffs) {
  symphas::problem_parameters_type pp{static_cast<int>(num_fields)};
  symphas::interval_element_type intvl;
  intvl.set_count(0.0, 32.0, 32);
  std::vector<symphas::interval_data_type> vdata(
      num_fields, symphas::interval_data_type(2, intvl));
  std::vector<symphas::b_data_type> bdata(
      num_fields, symphas::b_data_type(2, BoundaryType::PERIODIC));
  std::vector<symphas::init_data_type> tdata(
      num_fields, symphas::init_data_type(Inside::UNIFORM, {-0.1, 0.1}));
  pp.set_boundary_data(bdata.data(), num_fields);
  pp.set_initial_data(tdata.data(), num_fields);
  pp.set_interval_data(vdata.data(), num_fields);
  pp.set_time_step(0.001);

  StdoutCapture cap;
  cap.start();
  try {
    if (coeffs == nullptr) {
      ModelT m{pp};
      (void)m;
    } else {
      ModelT m{coeffs, n_coeffs, pp};
      (void)m;
    }
  } catch (...) {
  }
  return extract_equations(cap.stop());
}

// Report a finding.  `ok` is true if the audit passed.
void report(const char* model_tag, const char* check, bool ok,
            std::string const& evidence = "") {
  if (ok) {
    ++g_pass;
    std::fprintf(stderr, "  [OK]   %s :: %s\n", model_tag, check);
  } else {
    ++g_fail;
    std::fprintf(stderr, "  [FAIL] %s :: %s\n", model_tag, check);
    if (!evidence.empty())
      std::fprintf(stderr, "         %s\n", evidence.c_str());
  }
}

// Dump every captured equation to stderr so a human can audit.
void dump(const char* tag, std::vector<std::string> const& eqs) {
  std::fprintf(stderr, "\n[%s]  (%zu equation%s)\n", tag, eqs.size(),
               eqs.size() == 1 ? "" : "s");
  for (size_t i = 0; i < eqs.size(); ++i)
    std::fprintf(stderr, "  eq[%zu] %s\n", i, eqs[i].c_str());
}

bool contains(std::string const& hay, const char* needle) {
  return hay.find(needle) != std::string::npos;
}

// Concatenate all equations into one string for easier substring scanning.
std::string joined(std::vector<std::string> const& eqs) {
  std::string s;
  for (auto& e : eqs) {
    s += e;
    s += "\n";
  }
  return s;
}

}  // namespace

int testequationaudit() {
  g_pass = 0;
  g_fail = 0;
  std::fprintf(stderr,
               "--- testequationaudit: printed-equation suspicious-token "
               "audit ---\n");

  using Sp = SolverFT<Stencil2d2h<>>;

  // ============================================================
  // ModelA -- nothing exotic.  Reference baseline.
  // ============================================================
  {
    auto eqs = instantiate_and_capture<model_MA_t<2, Sp>>(1, nullptr, 0);
    dump("model_MA", eqs);
    auto j = joined(eqs);
    report("MA", "no trig of spatial coord (no sin(x)/cos(x))",
           !contains(j, "sin(x)") && !contains(j, "cos(x)") &&
               !contains(j, "sin(y)") && !contains(j, "cos(y)"),
           j);
    report("MA", "contains psi", contains(j, "psi"));
    report("MA", "contains laplacian V^2", contains(j, "V^2(psi)"));
  }

  // ============================================================
  // ModelB -- conserved.
  // ============================================================
  {
    auto eqs = instantiate_and_capture<model_MB_t<2, Sp>>(1, nullptr, 0);
    dump("model_MB", eqs);
    auto j = joined(eqs);
    report("MB", "no trig of spatial coord",
           !contains(j, "sin(x)") && !contains(j, "cos(x)") &&
               !contains(j, "sin(y)") && !contains(j, "cos(y)"),
           j);
    // Conserved => bilaplacian (V^4) must appear.
    report("MB", "contains V^4 (bilaplacian)", contains(j, "V^4"));
  }

  // ============================================================
  // ModelC -- coupled scalar.
  // ============================================================
  {
    auto eqs = instantiate_and_capture<model_MC_t<2, Sp>>(2, nullptr, 0);
    dump("model_MC", eqs);
    auto j = joined(eqs);
    report("MC", "no trig of spatial coord",
           !contains(j, "sin(x)") && !contains(j, "cos(x)") &&
               !contains(j, "sin(y)") && !contains(j, "cos(y)"),
           j);
    report("MC", "two equations", eqs.size() == 2);
  }

  // ============================================================
  // ModelH -- scalar + vector.
  // ============================================================
  {
    auto eqs = instantiate_and_capture<model_MH_t<2, Sp>>(2, nullptr, 0);
    dump("model_MH", eqs);
    auto j = joined(eqs);
    report("MH", "no trig of spatial coord",
           !contains(j, "sin(x)") && !contains(j, "cos(x)") &&
               !contains(j, "sin(y)") && !contains(j, "cos(y)"),
           j);
    report("MH", "two equations", eqs.size() == 2);
  }

  // ============================================================
  // MagneticPFC2013 -- Faghihi et al. 2013 paper coefficients.
  // The model definition uses `e(x)` and `e(y)` to represent unit
  // basis vectors.  But `e(...)` is defined as
  // `make_unit_vector<Dm>(...)` which interprets its argument as an
  // *angle in radians* -- and the symbol `x` in the EVOLUTION block
  // is the spatial coordinate (from `auto [x, y, z] = make_coords`).
  // So `e(x)` builds the position-dependent vector field
  // [cos(x); sin(x)] instead of the constant basis vector [1; 0].
  // This audit catches that.
  // ============================================================
  {
    double cs[] = {0.02,  0.98, 0.5, 1.0 / 3.0, 0.0,
                   1.0,   0.01, 0.04, 1.0,       0.001};
    auto eqs =
        instantiate_and_capture<model_MagneticPFC2013_t<2, Sp>>(2, cs, 10);
    dump("model_MagneticPFC2013", eqs);
    auto j = joined(eqs);

    // Structural sanity.
    report("MagneticPFC2013", "two equations", eqs.size() == 2);
    report("MagneticPFC2013", "contains rho", contains(j, "rho"));
    report("MagneticPFC2013", "contains V^4 (PFC stiffness chain)",
           contains(j, "V^4"));
    report("MagneticPFC2013", "contains psi^2 (coupling)",
           contains(j, "psi^2") || contains(j, "psi*psi"));

    // The suspicious-token check.  The model definition contains no
    // call to sin() or cos(); any sin/cos of the spatial coordinate
    // in the printed equation indicates the `e(x)`/`e(y)` macro is
    // producing a coordinate-dependent vector instead of a constant
    // basis vector.
    bool no_trig = !contains(j, "sin(x)") && !contains(j, "cos(x)") &&
                   !contains(j, "sin(y)") && !contains(j, "cos(y)");
    report("MagneticPFC2013",
           "no trig of spatial coord (e(x)/e(y) -> [1;0]/[0;1])", no_trig,
           j);

    // 2D curl(m) z-component is dm_y/dx - dm_x/dy.  The Poisson chain
    // prints as `d/d{x,y}(curl(m))`; if the curl was implemented as
    // dm_x/dx - dm_y/dy (a div-like fragment), the printed inner
    // expression is `drho_x/dx - drho_y/dy`.  Detect that explicitly.
    bool curl_ok = !contains(j, "drho_x/dx - drho_y/dy");
    report("MagneticPFC2013",
           "2D curl renders as dm_y/dx - dm_x/dy (not dm_x/dx - dm_y/dy)",
           curl_ok, j);

    // Cross-equation leak check.  The density equation dop(1) in the
    // model has NO PoissonSolver invocation -- the magnetic Poisson
    // chain is only in dop(2).  If `V^*` (the print glyph for
    // PoissonSolver) appears in eq[0], the algebra has leaked m-eqn
    // terms into the rho-eqn.  Setting c(10)=0 zeroes the runtime
    // contribution, but the structural leak is a real bug.  This is
    // CURRENTLY KNOWN-FAILING -- left in the suite as a regression
    // tracker until the algebra rule is found and fixed.
    if (!eqs.empty()) {
      bool no_leak = eqs[0].find("V^*") == std::string::npos;
      report("MagneticPFC2013",
             "no PoissonSolver leak into rho-eqn (KNOWN-FAILING)", no_leak,
             eqs[0]);
    }

    // Conservation: the rho-equation (dop(1)) is V^2(...) of something,
    // so the printed equation should NOT have a bare V^0 term outside a
    // V^2 wrapper.  Look for the outer laplacian.
    if (!eqs.empty()) {
      report("MagneticPFC2013", "rho-eqn outer V^2 (Model B conserved)",
             contains(eqs[0], "V^2("), eqs[0]);
    }
  }

  // ============================================================
  // FMPFCLinearField -- has provisional vars (var(1), var(2), var(3))
  // and *also* uses `e(x)` / `e(y)` for the unit basis vectors.
  // Same `e(...)` audit as MagneticPFC2013.
  // ============================================================
  {
    double cs[] = {0.02, 0.98, 0.5, 1.0 / 3.0, 0.0,  1.0,
                   0.01, 0.04, 1.0, 0.001,     0.1,  10.0};
    auto eqs =
        instantiate_and_capture<model_FMPFCLinearField_t<2, Sp>>(2, cs, 12);
    dump("model_FMPFCLinearField", eqs);
    auto j = joined(eqs);

    bool no_trig = !contains(j, "sin(x)") && !contains(j, "cos(x)") &&
                   !contains(j, "sin(y)") && !contains(j, "cos(y)");
    report("FMPFCLinearField",
           "no trig of spatial coord (e(x)/e(y) -> [1;0]/[0;1])", no_trig,
           j);
  }

  std::fprintf(stderr,
               "\n--- testequationaudit: %d check(s) passed, %d failed "
               "---\n",
               g_pass, g_fail);
  return g_fail;
}
