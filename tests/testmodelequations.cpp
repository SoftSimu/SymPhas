// Regression test for printed equation form of each registered model.
// See testmodelequations.h for design.

#include "testmodelequations.h"

#include "symphas.h"

// Pull in the standard model set (MA, MB, MC, MH).
#define MODEL_SET_1
#define PoissonSolver(E) expr::poisson_solver(E)

#include "modeldefinitions.h"
#include "modelinclude.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <regex>
#include <unistd.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>

namespace {

// Lightweight stdout capture: dup, redirect to a temp file, then
// read back into a string when stop_capture() is called.  We can't
// use freopen because SYMPHAS_LOG is a macro that expands to stdout,
// which is referenced through the file descriptor table.
struct StdoutCapture {
    int saved_stdout = -1;
    int tmp_fd = -1;
    std::string path = "/tmp/sym_equation_capture.XXXXXX";

    void start() {
        std::fflush(stdout);
        saved_stdout = ::dup(fileno(stdout));
        // mkstemp wants a writable buffer.
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

// Extract every "given equation -> ..." line into a vector.
std::vector<std::string> extract_equations(std::string const& captured) {
    std::vector<std::string> eqs;
    std::regex re("given equation\\s*->\\s*([^\\n]+)");
    auto begin = std::sregex_iterator(captured.begin(), captured.end(), re);
    auto end = std::sregex_iterator();
    for (auto it = begin; it != end; ++it) {
        std::string line = (*it)[1].str();
        // strip trailing whitespace
        while (!line.empty() && std::isspace(static_cast<unsigned char>(line.back())))
            line.pop_back();
        eqs.push_back(line);
    }
    return eqs;
}

// Normalise small differences (extra spaces, trailing zeros) so the
// comparison is structural.
std::string normalize(std::string const& s) {
    std::string out;
    out.reserve(s.size());
    bool prev_space = false;
    for (char c : s) {
        if (std::isspace(static_cast<unsigned char>(c))) {
            if (!prev_space && !out.empty()) out.push_back(' ');
            prev_space = true;
        } else {
            out.push_back(c);
            prev_space = false;
        }
    }
    while (!out.empty() && out.back() == ' ') out.pop_back();
    return out;
}

static int g_pass = 0;
static int g_fail = 0;

template <typename ModelT>
void check_model(const char* tag,
                 std::vector<const char*> const& expected,
                 size_t num_fields = 1,
                 double const* coeffs = nullptr, size_t n_coeffs = 0) {
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
    pp.set_time_step(0.01);

    StdoutCapture cap;
    cap.start();
    try {
        if (coeffs == nullptr) {
            ModelT model{pp};
            (void)model;
        } else {
            ModelT model{coeffs, n_coeffs, pp};
            (void)model;
        }
    } catch (...) {
    }
    std::string captured = cap.stop();
    auto eqs = extract_equations(captured);

    std::fprintf(stderr, "\n[%s]\n", tag);
    bool model_ok = true;
    for (size_t i = 0; i < std::max(eqs.size(), expected.size()); ++i) {
        std::string got = i < eqs.size() ? normalize(eqs[i]) : std::string{"<missing>"};
        std::string want = i < expected.size() ? normalize(expected[i]) : std::string{"<unexpected>"};
        bool ok = (got == want);
        if (!ok) model_ok = false;
        std::fprintf(stderr, "  eq[%zu] %s\n", i, ok ? "OK" : "FAIL");
        if (!ok) {
            std::fprintf(stderr, "    got:  %s\n", got.c_str());
            std::fprintf(stderr, "    want: %s\n", want.c_str());
        } else {
            std::fprintf(stderr, "    %s\n", got.c_str());
        }
    }
    if (model_ok) ++g_pass; else ++g_fail;
}

}  // namespace

int testmodelequations() {
    g_pass = 0;
    g_fail = 0;
    std::fprintf(stderr, "--- testmodelequations: printed equation regression ---\n");

    using Sp = SolverFT<Stencil2d2h<>>;

    // Discover phase: emit every model's equations and skip assertions.
    // Useful when intentionally changing the print form.
    if (std::getenv("SYMPHAS_TEST_EQ_DISCOVER")) {
        check_model<model_MA_t<2, Sp>>("model_MA discover", {}, 1);
        check_model<model_MB_t<2, Sp>>("model_MB discover", {}, 1);
        check_model<model_MC_t<2, Sp>>("model_MC discover", {}, 2);
        check_model<model_MH_t<2, Sp>>("model_MH discover", {}, 2);
        double cs[] = {0.02, 0.98, 0.5, 1.0/3.0, 0.0, 1.0, 0.01, 0.04, 1.0, 0.001};
        check_model<model_MagneticPFC2013_t<2, Sp>>(
            "model_MagneticPFC2013 discover", {}, 2, cs, 10);
        std::fprintf(stderr,
                     "\n--- discover-only run, no assertions made ---\n");
        return 0;
    }

    // Canonical regression contracts.  Locked in from the discover-phase
    // output once verified by inspection.  Defaults: c(N) = 1.

    // ModelA -- dpsi = lap(psi) + (c1 - 4 c2 psi^2) psi
    check_model<model_MA_t<2, Sp>>(
        "model_MA",
        {"V^2(psi) + psi - psi^3"},
        1);

    // ModelB -- conserved: dpsi = -bilap(psi) - lap((c1 - c2 psi^2) psi)
    check_model<model_MB_t<2, Sp>>(
        "model_MB",
        {"-V^4(psi) - V^2(psi) - V^2(-psi^3)"},
        1);

    // ModelC -- coupled scalar (MA + MB)
    check_model<model_MC_t<2, Sp>>(
        "model_MC",
        {
          "-V^4(psi) - V^2(-psi^3) - V^2(m^2) - V^2(psi)",
          "2.psi*m + V^2(m) + m - m^3",
        },
        2);

    // ModelH -- scalar + vector, Hohenberg-Halperin form
    check_model<model_MH_t<2, Sp>>(
        "model_MH",
        {
          "-V^2(m) - V^2(V^2(m)) - V^2(-m^3) - d/dy(m*j_y) - d/dx(m*j_x)",
          "(-m)*([1;0]d/dx(V^2(m)) + [0;1]d/dy(V^2(m)) + [1;0]d/dx(-m^3) + [0;1]d/dy(-m^3) + [1;0]dm/dx + [0;1]dm/dy) + [0;1]V^2(j_y) + [1;0]V^2(j_x)",
        },
        2);

    // MagneticPFC2013 -- the full model with PoissonSolver chain.
    // Coefficients are the paper values from magnetic_pfc_ref.
    {
        StdoutCapture cap; cap.start();
        symphas::problem_parameters_type pp{2};
        symphas::interval_element_type iv; iv.set_count(0.0, 32.0, 32);
        symphas::interval_data_type vd(2, iv);
        symphas::interval_data_type vds[2] = {vd, vd};
        symphas::b_data_type bd(2, BoundaryType::PERIODIC);
        symphas::b_data_type bds[2] = {bd, bd};
        symphas::init_data_type td(Inside::UNIFORM, {-0.1, 0.1});
        symphas::init_data_type tds[2] = {td, td};
        pp.set_boundary_data(bds, 2); pp.set_initial_data(tds, 2);
        pp.set_interval_data(vds, 2); pp.set_time_step(0.001);
        double cs[] = {0.02, 0.98, 0.5, 1.0/3.0, 0.0, 1.0, 0.01, 0.04, 1.0, 0.001};
        try { model_MagneticPFC2013_t<2, Sp> m{cs, 10, pp}; (void)m; }
        catch (...) {}
        auto eqs = extract_equations(cap.stop());
        std::fprintf(stderr, "\n[model_MagneticPFC2013]\n");
        for (size_t i = 0; i < eqs.size(); ++i)
            std::fprintf(stderr, "  eq[%zu]  %s\n", i, eqs[i].c_str());
        bool ok = (eqs.size() == 2);
        auto contains = [&](size_t i, const char* needle) {
            if (i >= eqs.size()) return false;
            return eqs[i].find(needle) != std::string::npos;
        };
        // dop(1): must contain laplacian-of-bilaplacian, lap-of-psi^3,
        // psi*rho coupling terms.
        ok &= contains(0, "V^2(");
        ok &= contains(0, "V^4(psi)");
        ok &= contains(0, "psi^3");
        ok &= contains(0, "rho");
        // dop(2): vector RHS with laplacian of magnetization components,
        // psi^2 coupling, and the PoissonSolver chain (renders as V^*).
        ok &= contains(1, "V^2(rho_");
        ok &= contains(1, "psi^2");
        ok &= contains(1, "drho_x/dx");
        ok &= contains(1, "drho_y/dy");
        if (ok) ++g_pass;
        else { ++g_fail; std::fprintf(stderr, "  FAIL: missing canonical pieces\n"); }
    }

    std::fprintf(stderr,
                 "\n--- testmodelequations: %d pass, %d fail ---\n",
                 g_pass, g_fail);
    return g_fail;
}
