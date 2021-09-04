
/* ***************************************************************************
 * This file is part of the SymPhas library, a framework for implementing
 * solvers for phase-field problems with compile-time symbolic algebra.
 *
 * Copyright (c) 2018-2021 by Steven A. Silber and Mikko Karttunen
 *
 * SymPhas is free software, which can be redistributed or modified under
 * the terms of the GNU Lesser General Public License (LGPL) as published
 * by the Free Software Foundation; LGPL version 3, or later versions at
 * your choice.
 *
 * SymPhas is distributed with the faith that it will be helpful and
 * practical but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser
 * General Public License for more details.
 *
 * ***************************************************************************
 *
 * MODULE:  expr
 * PURPOSE: Forward declaration of all expression objects.
 *
 * ***************************************************************************
 */

#pragma once

#include "spslib.h"

/*!
 * \defgroup Op Symbolic Algebra Objects
 * @{
 */

struct OpIdentity;
struct OpNegIdentity;
struct OpVoid;
template<typename T>
struct OpLiteral;


template<typename E1, typename E2>
struct OpBinaryAdd;
template<typename E1, typename E2>
struct OpBinarySub;
template<typename E1, typename E2>
struct OpBinaryMul;
template<typename E1, typename E2>
struct OpBinaryDiv;



template<typename E>
struct OpExpression;
template<typename E>
struct OpOperator;


template<size_t Z, typename G = OpVoid>
struct Variable;
template<typename G>
struct NamedData;
template<typename T, typename G>
struct OpLVariable;
template<typename T, typename... Gs>
struct OpNLVariable;

template<typename Dd, typename V, typename E, typename T>
struct OpFuncDerivative;
template<size_t O, typename V, typename T>
struct OpOperatorDerivative;
template<Axis ax, size_t O, typename V, typename Sp>
struct OpOperatorDirectionalDerivative;

template<typename V, typename E1, typename E2>
struct OpFuncConvolution;
template<size_t D>
struct GaussianSmoothing;

template<typename A1, typename A2, typename E>
struct OpCombination;
template<typename A1, typename A2>
struct OpOperatorCombination;
template<typename A1, typename A2, typename E>
struct OpChain;
template<typename A1, typename A2>
struct OpOperatorChain;
template<typename E>
struct OpOperator;
template<typename V, typename E>
struct OpExponential;
template<typename G, typename V, typename E>
struct OpMap;

template<typename V, typename E, typename F, typename Arg, typename... Args>
struct OpFunc;
template<auto f, typename V, typename E>
struct OpFuncApply;


template<size_t D, Axis ax>
struct GridAxis;


//! @}
