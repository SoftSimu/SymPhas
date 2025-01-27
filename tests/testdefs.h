#pragma once

#include "modelmacros.h"

// a test case to make sure smoothing is done correctly
MODEL_TYPE(SMOOTHING, SCALAR) NO_PROVISIONAL
EVOLUTION(
	dpsi = smoothing(psi) - psi)
LINK_WITH_NAME(SMOOTHING, SMOOTHING)

// ensure that convolution maintains a correct indexing pattern for
// a complex valued order parameter
MODEL_TYPE(COMPLEX_CONVOLUTION, COMPLEX) NO_PROVISIONAL
EVOLUTION(
	dpsi = smoothing(psi) - psi)
LINK_WITH_NAME(COMPLEX_CONVOLUTION, COMPLEX_CONVOLUTION)

/* based on model B
 * the state is in the laplacian because an expression is passed
 * through tests (on the main computer), I found this has equivalent
 * performance to the regular modelb
 */
MODEL_TYPE(STATE_SPEED, SCALAR) NO_PROVISIONAL
EVOLUTION(
	dpsi = -bilap(psi) + lap(-(A - lit(4.) * B * psi * psi) * psi))
LINK_WITH_NAME(STATE_SPEED, STATE_SPEED)


// Model B modified, performing laplacian in provisional equation so
// that it doesn't use the bilaplacian
MODEL_TYPE(BB, SCALAR) PROVISIONAL(SCALAR)
PROVISIONAL_DEF(
var(1) = lap(psi) + (A - lit(4.) * B * psi * psi) * psi)
EVOLUTION(
	dpsi = -lap(var(1))
)
LINK_WITH_NAME(BB, MODELBB)

// Model H as the way written in original manuscript
MODEL_TYPE(H2, SCALAR, VECTOR) PROVISIONAL(SCALAR)
PROVISIONAL_DEF(
	var(1) = -(A - lit(4.) * B * psi * psi) * psi)
EVOLUTION(
	dpsi = -bilap(psi) + lap(var(1)) - E * grad(psi) * rho,
	drho = lap(rho) + E * grad(psi) * (lap(psi) - var(1))
)
LINK_WITH_NAME(H2, MODELH2)

