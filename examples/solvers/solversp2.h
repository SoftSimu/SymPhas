#include "solver.h"
#ifdef USING_FFTW
#include <memory>
#include "gridfunctions.h"
#include "spectrallib.h"
#include "stencilincludes.h"
#ifdef USING_MPI
#include "spsmpi.h"
#endif
namespace solver_sp2_detail {
  template <typename T> struct is_derivative : std::false_type {};
  template <typename Dd, typename V, typename E, typename Sp>
  struct is_derivative<OpDerivative<Dd, V, E, Sp>> : std::true_type {};
  template <typename A1, typename A2, typename E>
  struct is_derivative<OpChain<A1, A2, E>> : std::true_type {};
  template <typename A1, typename A2, typename E>
  struct is_derivative<OpCombination<A1, A2, E>> : std::true_type {};

  // Check if an expression contains any derivative operators (for OpAdd sums).
  template <typename T> struct contains_derivative : std::false_type {};
  template <typename Dd, typename V, typename E, typename Sp>
  struct contains_derivative<OpDerivative<Dd, V, E, Sp>> : std::true_type {};
  template <typename A1, typename A2, typename E>
  struct contains_derivative<OpChain<A1, A2, E>> : std::true_type {};
  template <typename A1, typename A2, typename E>
  struct contains_derivative<OpCombination<A1, A2, E>> : std::true_type {};
  template <typename... Es>
  struct contains_derivative<OpAdd<Es...>> : std::bool_constant<(contains_derivative<Es>::value || ...)> {};

  // Check if ALL terms in an OpAdd are derivatives (for factoring out common operator).
  template <typename T> struct all_terms_derivative : is_derivative<T> {};
  template <typename... Es>
  struct all_terms_derivative<OpAdd<Es...>> : std::bool_constant<(is_derivative<Es>::value && ...)> {};

  // Factor out the common derivative operator from an OpAdd where every term
  // is a derivative. Uses separate_operator on each summand and collects the
  // inner expressions with expr::add_all. The operator is taken from the first term.
  template <typename... Es, size_t... Is>
  auto factor_derivative_sum_impl(OpAdd<Es...> const& e, std::index_sequence<Is...>) {
    auto ops = std::make_tuple(expr::split::separate_operator(expr::get<Is>(e))...);
    return std::make_pair(
        std::get<0>(std::get<0>(ops)),
        expr::add_all(std::get<1>(std::get<Is>(ops))...));
  }
  template <typename... Es>
  auto factor_derivative_sum(OpAdd<Es...> const& e) {
    return factor_derivative_sum_impl(e, std::make_index_sequence<sizeof...(Es)>{});
  }
}
template <typename NL_type>
struct SpectralDataSP2 {
  std::shared_ptr<complex_t[]> A_data;
  std::shared_ptr<complex_t[]> B_data;
  std::shared_ptr<scalar_t[]> nl_scratch;
  len_type klen;
  len_type rlen;
  NL_type nl_expr;
  // Cross-field spectral contributions: one kernel per coupled field.
  // Each entry stores B(k)*D(k)*L_cross_j(k) and a pointer to field j's frame_t.
  static constexpr size_t MAX_CROSS = 8;
  std::shared_ptr<complex_t[]> cross_kernels[MAX_CROSS];
  complex_t* cross_frame_ts[MAX_CROSS];
  size_t num_cross;
#ifdef USING_CUDA
  // Device mirrors of the spectral operators for the GPU spectral path.
  // Lazily allocated/filled on the first CUDA equation() call (the A/B
  // arrays are constant in time). Stored as void* so this header needs no
  // cuFFT types in non-CUDA translation units; reinterpret_cast at use.
  mutable void* A_dev = nullptr;
  mutable void* B_dev = nullptr;
  mutable bool dev_ready = false;
#endif
  SpectralDataSP2(complex_t* A, complex_t* B, len_type klen, len_type rlen, NL_type const& nl)
      : A_data{A, std::default_delete<complex_t[]>()},
        B_data{B, std::default_delete<complex_t[]>()},
        nl_scratch{new scalar_t[rlen], std::default_delete<scalar_t[]>()},
        klen{klen}, rlen{rlen}, nl_expr{nl},
        cross_kernels{}, cross_frame_ts{}, num_cross{0} {}
  void update() { expr::prune::update(nl_expr); }
};
START_NEW_SOLVER(SolverSP2)
double h[3];
SolverSP2(const double* h = nullptr, double dt = 1, size_t dim = 0)
    : parent_type(dt), h{0} { std::copy(h, h + dim, this->h); }
template <typename S> void step(S&&) const {}

template <size_t Z, typename S, typename NL_type>
void equation(std::pair<Variable<Z, symphas::ref<S>>, SpectralDataSP2<NL_type>>& r) const {
  SYMPHAS_MPI_PROFILE_SCOPE("sp2_equation");
  auto& [sys, data] = r;
  auto& s = sys.get();
  {
    SYMPHAS_MPI_PROFILE_SCOPE("sp2_data_update");
    data.update();
  }
  len_type rlen = grid::length(s);
#if defined(USING_FFTW) && defined(USING_MPI) && defined(USING_FFTW_MPI)
  if constexpr (is_spectral_mpi_system_v<S>) {
    constexpr size_t Dm = grid::dimension_of<S>::value;
    // MPI path: save/restore only local slab, no Allgatherv.
    len_type lstart = s.local_real_start();
    len_type llen = s.local_real_len();
    // Save local slab of ψ.
    {
      SYMPHAS_MPI_PROFILE_SCOPE("sp2_save_psi");
      std::copy(s.values + lstart, s.values + lstart + llen, data.nl_scratch.get());
    }
    // Evaluate NL only on this rank's local slab via region_interval
    // overload of expr::result. The non-local rows of s.values are stale
    // (neighbors' data from a prior step), but they are not read here
    // because the expression iterator steps only over [lstart, lstart+llen).
    // Only local entries of ψ are valid; scatter_real_to_work reads
    // only the local slab of the result. The slab is the trailing axis
    // (y in 2D, z in 3D), so only that axis is restricted to the rank range.
    {
      SYMPHAS_MPI_PROFILE_SCOPE("sp2_nl_eval");
      len_type global_dims[Dm];
      len_type intervals[Dm][2];
      for (size_t d = 0; d < Dm; ++d) {
        global_dims[d] = s.dims[d];
        intervals[d][0] = 0;
        intervals[d][1] = s.dims[d];
      }
      intervals[Dm - 1][0] = static_cast<len_type>(s.local_0_start);
      intervals[Dm - 1][1] =
          static_cast<len_type>(s.local_0_start + s.local_n0);
      grid::region_interval<Dm> local_region{global_dims, intervals};
      expr::result(data.nl_expr, s.values, local_region);
    }
    // Scatter local NL into padded workspace and forward FFT.
    {
      SYMPHAS_MPI_PROFILE_SCOPE("sp2_scatter_to_work");
      s.scatter_real_to_work();
    }
    // Restore local slab of ψ.
    {
      SYMPHAS_MPI_PROFILE_SCOPE("sp2_restore_psi");
      std::copy(data.nl_scratch.get(), data.nl_scratch.get() + llen, s.values + lstart);
    }
    {
      SYMPHAS_MPI_PROFILE_SCOPE("sp2_fwd_fft");
      symphas::dft::fftw_execute(s.p_to_t);
    }
    {
      SYMPHAS_MPI_PROFILE_SCOPE("sp2_spectral_mul");
      // A/B are global arrays; offset to this rank's local k-space slab.
      // The k-space slab stride per slab index is (Nx/2+1) in 2D and
      // Ny*(Nx/2+1) in 3D (the trailing real axis is half-spectrum).
      len_type tlen = s.transformed_len;
      len_type kx_len = s.dims[0] / 2 + 1;
      len_type kplane = (Dm == 2) ? kx_len : s.dims[1] * kx_len;
      len_type offset = static_cast<len_type>(s.local_0_start) * kplane;
      auto* A = data.A_data.get() + offset;
      auto* B = data.B_data.get() + offset;
      for (len_type i = 0; i < tlen; ++i)
        s.dframe[i] = A[i] * s.frame_t[i] + B[i] * s.dframe[i];
    }
  } else
#endif
#if defined(USING_CUDA)
  if constexpr (is_spectral_cuda_system_v<S>) {
    // GPU spectral path. All cuFFT calls and kernel launches are encapsulated
    // in SolverSystemSpectralCUDA methods (defined in solversystem.cuh, only
    // under nvcc) so this shared header carries no device syntax.
    s.upload_operators(data.A_dev, data.B_dev, data.dev_ready,
                       data.A_data.get(), data.B_data.get());
    // Save field, evaluate NL into the device field, transform + spectral
    // step, then restore the field for the next iteration's NL evaluation.
    // The NL expression must be evaluated on the GPU: select the kernel eval
    // handler (KernelEvalHandler for CUDA storage) exactly as the FD CUDA
    // solver does, so the symbolic evaluation runs in device kernels rather
    // than the host loop (which would dereference device memory on the host).
    s.save_field();
    {
      using nl_t = std::decay_t<decltype(data.nl_expr)>;
      expr::eval_handler_type<nl_t> handler;
      // Pass the device GridCUDA object (not a raw pointer): the CUDA
      // evaluate_expression_trait constructor binds GridCUDA<T,D>& and writes
      // its .values on device. grid::length(s) selects the contiguous-length
      // overload (full-grid evaluation).
      handler.result(data.nl_expr, s.as_grid(), grid::length(s));
    }
    s.solve_spectral_step(data.A_dev, data.B_dev);
    s.restore_field();
  } else
#endif
  {
    // Serial path.
    std::copy(s.values, s.values + rlen, data.nl_scratch.get());
    expr::result(data.nl_expr, s.values, rlen);
    constexpr size_t D = grid::dimension_of<S>::value;
    symphas::dft::arrange_fftw_stip<D>(s.values, reinterpret_cast<scalar_t*>(s.dframe), s.dims);
    std::copy(data.nl_scratch.get(), data.nl_scratch.get() + rlen, s.values);
    symphas::dft::fftw_execute(s.p_to_t);
    len_type tlen = s.transformed_len;
    auto* A = data.A_data.get();
    auto* B = data.B_data.get();
    for (len_type i = 0; i < tlen; ++i)
      s.dframe[i] = A[i] * s.frame_t[i] + B[i] * s.dframe[i];
    // Add cross-field spectral contributions.
    for (size_t c = 0; c < data.num_cross; ++c) {
      auto* K = data.cross_kernels[c].get();
      auto* Ft = data.cross_frame_ts[c];
      for (len_type i = 0; i < tlen; ++i)
        s.dframe[i] += K[i] * Ft[i];
    }
  }
}
template <size_t En, typename... Ss, size_t Z, typename S, typename E,
          typename T_src = typename grid::value_type_of<S>::type,
          size_t D = grid::dimension_of<S>::value>
decltype(auto) form_expr_one(std::tuple<Ss...> const& systems,
    std::pair<Variable<Z, symphas::ref<S>>, E> const& e) const {
  auto&& [sys, equation] = e;
  auto&& [linear, nonlinear] = expr::split::by_linear(expr::apply_operators(equation));
  auto&& [linear_in_Z, linear_in_nonZ] = expr::split::separate_var<Z>(linear);
  auto&& [l_op, non_op] = solver_sp::get_l_op<Z>(linear_in_Z, h);
  auto&& A_expression = solver_sp::form_A_op<D>(l_op, dt, sys.get().dims);
  auto&& B_expression = solver_sp::form_B_op<D>(l_op, dt, sys.get().dims);
  len_type tlen = sys.get().transformed_len;
  len_type ab_len = tlen;
#if defined(USING_MPI) && defined(USING_FFTW_MPI)
  if constexpr (is_spectral_mpi_system_v<S>) {
    // A/B(k) are built for the FULL global k-space (each rank slices its slab
    // at evaluation time). Global k-space size: 2D = Ny*(Nx/2+1),
    // 3D = Nz*Ny*(Nx/2+1); the trailing real axis Nx is half-spectrum.
    len_type kx_len = sys.get().dims[0] / 2 + 1;
    if constexpr (D == 2) {
      ab_len = static_cast<len_type>(sys.get().dims[1] * kx_len);
    } else {
      ab_len = static_cast<len_type>(
          sys.get().dims[2] * sys.get().dims[1] * kx_len);
    }
  }
#endif
  complex_t* A_data = new complex_t[ab_len];
  complex_t* B_data = new complex_t[ab_len];
  {
    auto A_grid = expr::transform::to_grid(solver_sp::sthc_apply_on_scalar<T_src>(A_expression));
    auto B_grid = expr::transform::to_grid(solver_sp::sthc_apply_on_scalar<T_src>(B_expression));
    if constexpr (std::is_arithmetic_v<decltype(A_grid)>) {
      for (len_type i = 0; i < ab_len; ++i) { A_data[i] = complex_t(A_grid, 0); B_data[i] = complex_t(B_grid, 0); }
    } else {
      std::copy(A_grid.values, A_grid.values + ab_len, A_data);
      std::copy(B_grid.values, B_grid.values + ab_len, B_data);
    }
  }
  expr::printe(A_expression, "SP2 A(k)");

  // ---------- Nonlinear handling ----------
  // The NL has three components with potentially different derivative structures:
  //   nonlinear:      may have a common derivative from conserved dynamics (e.g., ∇²)
  //   linear_in_nonZ: cross-field linear terms, may have different derivative orders
  //   non_op:         typically OpVoid, leftover self-linear terms
  // Handle nonlinear and cross-field SEPARATELY to preserve correct operators.

  // Step 1: Factor any common derivative from the self-field nonlinear part.
  // For conserved PFC: nonlinear = ∇²(n² + n³ + ...) → factor ∇², fold into B(k).
  // For non-conserved: nonlinear = n² + n³ + ... → no derivative to factor.
  using nonlinear_type = std::decay_t<decltype(nonlinear)>;
  auto factor_nl = [&]() {
    if constexpr (solver_sp2_detail::is_derivative<nonlinear_type>::value) {
      return expr::split::separate_operator(nonlinear);
    } else if constexpr (solver_sp2_detail::all_terms_derivative<nonlinear_type>::value) {
      return solver_sp2_detail::factor_derivative_sum(nonlinear);
    } else {
      return std::make_pair(OpIdentity{}, nonlinear);
    }
  };
  auto [nl_deriv_op, nl_pointwise] = factor_nl();

  // Apply the extracted derivative to B(k) spectrally.
  if constexpr (!std::is_same_v<std::decay_t<decltype(nl_deriv_op)>, OpIdentity>) {
    auto deriv_k = expr::transform::to_ft<D>(nl_deriv_op, h, sys.get().dims);
    auto deriv_grid = expr::transform::to_grid(solver_sp::sthc_apply_on_scalar<T_src>(deriv_k));
    if constexpr (!std::is_arithmetic_v<decltype(deriv_grid)>) {
      for (len_type i = 0; i < ab_len; ++i) B_data[i] *= deriv_grid.values[i];
    }
    expr::printe(nl_deriv_op, "SP2 factored derivative (applied to B(k))");
  }

  // Build the real-space NL expression: non_op + factored nonlinear inner.
  auto nl_real = non_op + nl_pointwise;

  // Step 2: Handle cross-field linear terms (linear_in_nonZ).
  // For multi-field systems, linear_in_nonZ contains derivative operators on
  // other fields (e.g., (1+∇²)²n_other). Extract the spectral operators
  // for each other field and precompute cross kernels B(k)*D(k)*L_cross(k).
  using cross_type = std::decay_t<decltype(linear_in_nonZ)>;
  constexpr bool has_cross = !std::is_same_v<cross_type, OpVoid> && sizeof...(Ss) > 1;

  if constexpr (has_cross) {
    // Use get_l_op on linear_in_nonZ with each other variable index to extract
    // the k-space operator. For 2-field systems, all cross terms reference one
    // other field. For N>2, separate_var first isolates each field's terms.
    auto data = SpectralDataSP2<decltype(nl_real)>(A_data, B_data, tlen,
        grid::length(sys.get()), nl_real);

    // For each other field J != Z: extract spectral operator, build cross kernel.
    // Use separate_var<J> to isolate field J's terms, then get_l_op<J>.
    auto add_cross_for = [&]<size_t J>() {
      auto&& [terms_J, terms_rest] = expr::split::separate_var<J>(linear_in_nonZ);
      using tj = std::decay_t<decltype(terms_J)>;
      if constexpr (!std::is_same_v<tj, OpVoid>) {
        auto&& [lop_J, non_J] = solver_sp::get_l_op<J>(terms_J, h);
        using lj = std::decay_t<decltype(lop_J)>;
        if constexpr (!std::is_same_v<lj, OpVoid>) {
          // Build cross kernel: B(k) * D(k) * lop_J(k)
          auto ck_expr = [&]() {
            if constexpr (!std::is_same_v<std::decay_t<decltype(nl_deriv_op)>, OpIdentity>) {
              auto dk = expr::transform::to_ft<D>(nl_deriv_op, h, sys.get().dims);
              return B_expression * dk * lop_J;
            } else {
              return B_expression * lop_J;
            }
          }();
          auto cg = expr::transform::to_grid(solver_sp::sthc_apply_on_scalar<T_src>(ck_expr));
          complex_t* ck = new complex_t[ab_len];
          if constexpr (std::is_arithmetic_v<decltype(cg)>) {
            for (len_type i = 0; i < ab_len; ++i) ck[i] = complex_t(cg, 0);
          } else {
            std::copy(cg.values, cg.values + ab_len, ck);
          }
          size_t idx = data.num_cross++;
          data.cross_kernels[idx] = std::shared_ptr<complex_t[]>(ck, std::default_delete<complex_t[]>());
          // frame_t is complex_t* for host systems and cufftDoubleComplex* for
          // the GPU spectral system; both are layout-compatible double2. The
          // GPU path targets single-field models (num_cross stays 0), so this
          // cross-field branch is only exercised by the host solvers, but the
          // assignment must still compile when instantiated for either type.
          data.cross_frame_ts[idx] =
              reinterpret_cast<complex_t*>(std::get<J>(systems).frame_t);
          expr::printe(lop_J, "SP2 cross-field (spectral)");
        }
      }
    };

    // Expand for each other field. For 2-field models we only instantiate one.
    // For N-field models, the fold expands at compile time for each J != Z.
    [&]<size_t... Is>(std::index_sequence<Is...>) {
      auto maybe_add = [&]<size_t J>(std::integral_constant<size_t, J>) {
        if constexpr (J != Z) {
          add_cross_for.template operator()<J>();
        }
      };
      (maybe_add(std::integral_constant<size_t, Is>{}), ...);
    }(std::make_index_sequence<sizeof...(Ss)>{});

    expr::printe(nl_real, "SP2 nonlinear (real-space)");
    return std::make_pair(sys, data);
  } else {
    // No cross-field terms: pure pointwise NL.
    expr::printe(nl_real, "SP2 nonlinear (real-space)");
    auto data = SpectralDataSP2<decltype(nl_real)>(A_data, B_data, tlen,
        grid::length(sys.get()), nl_real);
    return std::make_pair(sys, data);
  }
}
static auto make_solver(symphas::problem_parameters_type const& parameters) {
  size_t dim = parameters.get_dimension();
  double* h = new double[dim];
  for (iter_type i = 0; i < dim; ++i)
    h[i] = parameters.get_interval_data()[0].at(symphas::index_to_axis(i)).width();
  auto s = SolverSP2{h, parameters.get_time_step(), dim};
  delete[] h;
  return s;
}
END_SOLVER

#if defined(USING_CUDA)
// GPU spectral build: the SP2 solver uses the device-resident spectral system
// (cuFFT transforms + device spectral algebra). Manually written specialization
// (as the ASSOCIATE_SOLVER_SYSTEM_TYPE macro would generate) forwarding the
// simulation dimension D.
namespace symphas::internal {
template <>
struct solver_system_type_match<solver_id_type_SolverSP2, 0> {
  template <typename Ty, size_t D>
  using type = SolverSystemSpectralCUDA<D>;
};
constexpr solver_count_index<solver_id_type_SolverSP2, 1>
    solver_counter(solver_count_index<solver_id_type_SolverSP2, 1>);
}
#elif defined(USING_MPI) && defined(USING_FFTW_MPI)
// SolverSystemSpectralMPI<D> is selected (instead of the serial spectral
// system) when MPI is active. We manually write the specialization that the
// ASSOCIATE_SOLVER_SYSTEM_TYPE macro would generate, forwarding the simulation
// dimension D so the distributed-FFT system is 2D or 3D as appropriate.
namespace symphas::internal {
template <>
struct solver_system_type_match<solver_id_type_SolverSP2, 0> {
  template <typename Ty, size_t D>
  using type = SolverSystemSpectralMPI<D>;
};
constexpr solver_count_index<solver_id_type_SolverSP2, 1>
    solver_counter(solver_count_index<solver_id_type_SolverSP2, 1>);
}
#else
ASSOCIATE_SOLVER_SYSTEM_TYPE(SolverSP2, SolverSystemSpectral)
#endif

ASSOCIATE_PROVISIONAL_SYSTEM_TYPE(SolverSP2, ProvisionalSystemSpectral)
SYMPHAS_SOLVER_ALL_SUPPORTED(SolverSP2)
#endif
