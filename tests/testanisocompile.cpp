// Force compile-time instantiation of every AnisotropicFMPFC test-model
// variant.  At runtime this just touches a static counter so the linker
// keeps the model classes; the value is whether the file *compiles*.

#include "testanisocompile.h"

#include <cstddef>

namespace {

// Pull each ANISO_T* model into a small static so the model class is
// fully instantiated even if no one constructs it.  sizeof on the model
// type is enough to force the templates to expand.
template <typename T>
constexpr std::size_t force_instantiate_one() {
    return sizeof(T);
}

}  // namespace

void testanisocompile() {
    using Sp = SolverFT<Stencil2d2h<>>;
    std::size_t total = 0;
    total += force_instantiate_one<model_AnisoT1_t<2, Sp>>();
    total += force_instantiate_one<model_AnisoT2_t<2, Sp>>();
    total += force_instantiate_one<model_AnisoT3_t<2, Sp>>();
    total += force_instantiate_one<model_AnisoT4_t<2, Sp>>();
    total += force_instantiate_one<model_AnisoT5_t<2, Sp>>();
    total += force_instantiate_one<model_AnisoT6_t<2, Sp>>();
    total += force_instantiate_one<model_AnisoT3a_t<2, Sp>>();
    total += force_instantiate_one<model_AnisoT3b_t<2, Sp>>();
    total += force_instantiate_one<model_AnisoT3c_t<2, Sp>>();
    total += force_instantiate_one<model_AnisoT3d_t<2, Sp>>();
    total += force_instantiate_one<model_AnisoT3e_t<2, Sp>>();
    total += force_instantiate_one<model_AnisoT7_t<2, Sp>>();
    total += force_instantiate_one<model_AnisoT8_t<2, Sp>>();
    total += force_instantiate_one<model_AnisoL6_t<2, Sp>>();
    int n_models = 14;
#ifdef SYMPHAS_TEST_BUG4_HEAVY
    total += force_instantiate_one<model_AnisoL7_t<2, Sp>>();
    total += force_instantiate_one<model_AnisoL8_t<2, Sp>>();
    n_models += 2;
#endif
#ifdef SYMPHAS_TEST_ANISO_FMPFC_ORIGINAL
    total += force_instantiate_one<model_AnisotropicFMPFC_t<2, Sp>>();
    n_models += 1;
#endif
    std::printf("--- testanisocompile: %d model variants instantiated, "
                "total sizeof = %zu bytes ---\n",
                n_models, total);
}
