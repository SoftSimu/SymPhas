
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
 * MODULE:  lib
 * PURPOSE: Defines Fourier transform algorithms that can be used in solvers.
 *
 * ***************************************************************************
 */

#pragma once


#include "dft.h"


struct fftw_plan_s;
typedef double fftw_complex[2];
typedef fftw_plan_s* fftw_plan;

#ifdef EXECUTION_HEADER_AVAILABLE
#include <execution>
#endif

namespace symphas::dft
{

	/* Subset of FFTW functions used in SymPhas. FFTW is not exposed to
	 * the whole program, but a limited number of functions are specified
	 * and the rest are adapted to the specific workings of SymPhas.
	 */

	//! Executes the FFTW plan.
	void fftw_execute(fftw_plan p);
	//! Destroys the FFTW plan.
	void fftw_destroy_plan(fftw_plan p);
	//! Allocates memory aligned complex array of the prescribed size.
	fftw_complex* fftw_alloc_complex(size_t);
	//! Allocates memory aligned real-valued array of the prescribed size.
	scalar_t* fftw_alloc_real(size_t);
	//! Frees the memory associated with an aligned FFTW array.
	void fftw_free(fftw_complex*&);



	//! Creates a new FFTW plan from the given types and dimension.
	/*!
	 * Using the given dimension, construct a new FFTW plan based on
	 * transforming a source array of type `T_src` to a transform target
	 * of type `T_target`. The template parameters are used to specialize the
	 * plan type which is selected.
	 * 
	 * \tparam D The dimension of the system that is transformed.
	 * \tparam T_src The source system value type.
	 * \tparam T_target The type of the Fourier transform output.
	 */
	template<size_t D, typename T_src, typename T_target>
	struct new_fftw_plan;


	/* specialized plans, in particular implementing real to complex and vice versa
	 * note that these plans (when only complex values are involved) are purely forward, hence
	 * if a "backward" fourier transform needs to be performed, the output grid needs to be passed
	 * as the second argument
	 */

	//! Specialization based on symphas::dft::new_fftw_plan.
	template<>
	struct new_fftw_plan<1, complex_t, complex_t>
	{
		//! Creates a new complex to complex FFTW plan.
		/*!
		 * Creates a plan based on transforming \p src to \p target.
		 * The direction of the complex to complex plan is considered to be
		 * forward by default, unless \p backward is provided equal to true.
		 * 
		 * \param src The source array.
		 * \param target The target array.
		 * \param dims The dimension of the source array.
		 * \param estimate Whether or not to estimate the plan. If the plan
		 * is not estimated, the source and target arrays will be modified.
		 * \param backward Whether or not the plan refers to taking a forward
		 * Fourier transform, or backward.
		 */
		fftw_plan operator()(fftw_complex* src, fftw_complex* target, const len_type* dims, bool estimate = false, bool backward = false);
		//! Forwards arbitrary arguments to create a specific FFTW plan.
		/*!
		 * Creates a plan based on transforming \p src to \p target.
		 *
		 * \param src The source array.
		 * \param target The target array.
		 * \param dims The dimension of the source array.
		 * \param estimate Whether or not to estimate the plan. If the plan
		 * is not estimated, the source and target arrays will be modified.
		 * \param backward Whether or not the plan refers to taking a forward
		 * Fourier transform, or backward.
		 */
		template<typename T1, typename T2>
		fftw_plan operator()(T1* src, T2* target, const len_type* dims, bool estimate = false, bool backward = false)
		{
			return new_fftw_plan<1, complex_t, complex_t>{}(reinterpret_cast<fftw_complex*>(src), reinterpret_cast<fftw_complex*>(target), dims, estimate, backward);
		}
	};

	//! Specialization based on symphas::dft::new_fftw_plan.
	template<>
	struct new_fftw_plan<2, complex_t, complex_t>
	{
		//! Creates a new complex to complex FFTW plan.
		/*!
		 * Creates a plan based on transforming \p src to \p target.
		 * The direction of the complex to complex plan is considered to be
		 * forward by default, unless \p backward is provided equal to true.
		 *
		 * \param src The source array.
		 * \param target The target array.
		 * \param dims The dimension of the source array.
		 * \param estimate Whether or not to estimate the plan. If the plan
		 * is not estimated, the source and target arrays will be modified.
		 * \param backward Whether or not the plan refers to taking a forward
		 * Fourier transform, or backward.
		 */
		fftw_plan operator()(fftw_complex* src, fftw_complex* target, const len_type* dims, bool estimate = false, bool backward = false);
		//! Forwards arbitrary arguments to create a specific FFTW plan.
		/*!
		 * Creates a plan based on transforming \p src to \p target.
		 *
		 * \param src The source array.
		 * \param target The target array.
		 * \param dims The dimension of the source array.
		 * \param estimate Whether or not to estimate the plan. If the plan
		 * is not estimated, the source and target arrays will be modified.
		 * \param backward Whether or not the plan refers to taking a forward
		 * Fourier transform, or backward.
		 */
		template<typename T1, typename T2>
		fftw_plan operator()(T1* src, T2* target, const len_type* dims, bool estimate = false, bool backward = false)
		{
			return new_fftw_plan<2, complex_t, complex_t>{}(reinterpret_cast<fftw_complex*>(src), reinterpret_cast<fftw_complex*>(target), dims, estimate, backward);
		}
	};

	//! Specialization based on symphas::dft::new_fftw_plan.
	template<>
	struct new_fftw_plan<3, complex_t, complex_t>
	{
		//! Creates a new complex to complex FFTW plan.
		/*!
		 * Creates a plan based on transforming \p src to \p target.
		 * The direction of the complex to complex plan is considered to be
		 * forward by default, unless \p backward is provided equal to true.
		 *
		 * \param src The source array.
		 * \param target The target array.
		 * \param dims The dimension of the source array.
		 * \param estimate Whether or not to estimate the plan. If the plan
		 * is not estimated, the source and target arrays will be modified.
		 * \param backward Whether or not the plan refers to taking a forward
		 * Fourier transform, or backward.
		 */
		fftw_plan operator()(fftw_complex* src, fftw_complex* target, const len_type* dims, bool estimate = false, bool backward = false);
		//! Forwards arbitrary arguments to create a specific FFTW plan.
		/*!
		 * Creates a plan based on transforming \p src to \p target.
		 *
		 * \param src The source array.
		 * \param target The target array.
		 * \param dims The dimension of the source array.
		 * \param estimate Whether or not to estimate the plan. If the plan
		 * is not estimated, the source and target arrays will be modified.
		 * \param backward Whether or not the plan refers to taking a forward
		 * Fourier transform, or backward.
		 */
		template<typename T1, typename T2>
		fftw_plan operator()(T1* src, T2* target, const len_type* dims, bool estimate = false, bool backward = false)
		{
			return new_fftw_plan<3, complex_t, complex_t>{}(reinterpret_cast<fftw_complex*>(src), reinterpret_cast<fftw_complex*>(target), dims, estimate, backward);
		}
	};

	//! Specialization based on symphas::dft::new_fftw_plan.
	template<>
	struct new_fftw_plan<1, scalar_t, complex_t>
	{
		//! Creates a new real to complex FFTW plan.
		/*!
		 * Creates a plan based on transforming \p src to \p target.
		 * The direction of the real to complex plan is always forward.
		 *
		 * \param src The source array.
		 * \param target The target array.
		 * \param dims The dimension of the source array.
		 * \param estimate Whether or not to estimate the plan. If the plan
		 * is not estimated, the source and target arrays will be modified.
		 */
		fftw_plan operator()(scalar_t* src, fftw_complex* target, const len_type* dims, bool estimate = false);
		//! Forwards arbitrary arguments to create a specific FFTW plan.
		/*!
		 * Creates a plan based on transforming \p src to \p target.
		 *
		 * \param src The source array.
		 * \param target The target array.
		 * \param dims The dimension of the source array.
		 * \param estimate Whether or not to estimate the plan. If the plan
		 * is not estimated, the source and target arrays will be modified.
		 */
		template<typename T1, typename T2>
		fftw_plan operator()(T1* src, T2* target, const len_type* dims, bool estimate = false)
		{
			return new_fftw_plan<1, scalar_t, complex_t>{}(reinterpret_cast<scalar_t*>(src), reinterpret_cast<fftw_complex*>(target), dims, estimate);
		}
	};

	//! Specialization based on symphas::dft::new_fftw_plan.
	template<>
	struct new_fftw_plan<2, scalar_t, complex_t>
	{
		//! Creates a new real to complex FFTW plan.
		/*!
		 * Creates a plan based on transforming \p src to \p target.
		 * The direction of the real to complex plan is always forward.
		 *
		 * \param src The source array.
		 * \param target The target array.
		 * \param dims The dimension of the source array.
		 * \param estimate Whether or not to estimate the plan. If the plan
		 * is not estimated, the source and target arrays will be modified.
		 */
		fftw_plan operator()(scalar_t* src, fftw_complex* target, const len_type* dims, bool estimate = false);
		//! Forwards arbitrary arguments to create a specific FFTW plan.
		/*!
		 * Creates a plan based on transforming \p src to \p target.
		 *
		 * \param src The source array.
		 * \param target The target array.
		 * \param dims The dimension of the source array.
		 * \param estimate Whether or not to estimate the plan. If the plan
		 * is not estimated, the source and target arrays will be modified.
		 */
		template<typename T1, typename T2>
		fftw_plan operator()(T1* src, T2* target, const len_type* dims, bool estimate = false)
		{
			return new_fftw_plan<2, scalar_t, complex_t>{}(reinterpret_cast<scalar_t*>(src), reinterpret_cast<fftw_complex*>(target), dims, estimate);
		}
	};

	//! Specialization based on symphas::dft::new_fftw_plan.
	template<>
	struct new_fftw_plan<3, scalar_t, complex_t>
	{
		//! Creates a new real to complex FFTW plan.
		/*!
		 * Creates a plan based on transforming \p src to \p target.
		 * The direction of the real to complex plan is always forward.
		 *
		 * \param src The source array.
		 * \param target The target array.
		 * \param dims The dimension of the source array.
		 * \param estimate Whether or not to estimate the plan. If the plan
		 * is not estimated, the source and target arrays will be modified.
		 */
		fftw_plan operator()(scalar_t* src, fftw_complex* target, const len_type* dims, bool estimate = false);
		//! Forwards arbitrary arguments to create a specific FFTW plan.
		/*!
		 * Creates a plan based on transforming \p src to \p target.
		 *
		 * \param src The source array.
		 * \param target The target array.
		 * \param dims The dimension of the source array.
		 * \param estimate Whether or not to estimate the plan. If the plan
		 * is not estimated, the source and target arrays will be modified.
		 */
		template<typename T1, typename T2>
		fftw_plan operator()(T1* src, T2* target, const len_type* dims, bool estimate = false)
		{
			return new_fftw_plan<3, scalar_t, complex_t>{}(reinterpret_cast<scalar_t*>(src), reinterpret_cast<fftw_complex*>(target), dims, estimate);
		}
	};

	//! Specialization based on symphas::dft::new_fftw_plan.
	template<>
	struct new_fftw_plan<1, complex_t, scalar_t>
	{
		//! Creates a new complex to real FFTW plan.
		/*!
		 * Creates a plan based on transforming \p src to \p target.
		 * The direction of the complex to real plan is always backward.
		 *
		 * \param src The source array.
		 * \param target The target array.
		 * \param dims The dimension of the source array.
		 * \param estimate Whether or not to estimate the plan. If the plan
		 * is not estimated, the source and target arrays will be modified.
		 */
		fftw_plan operator()(fftw_complex* src, scalar_t* target, const len_type* dims, bool estimate = false);
		//! Forwards arbitrary arguments to create a specific FFTW plan.
		/*!
		 * Creates a plan based on transforming \p src to \p target.
		 *
		 * \param src The source array.
		 * \param target The target array.
		 * \param dims The dimension of the source array.
		 * \param estimate Whether or not to estimate the plan. If the plan
		 * is not estimated, the source and target arrays will be modified.
		 */
		template<typename T1, typename T2>
		fftw_plan operator()(T1* src, T2* target, const len_type* dims, bool estimate = false)
		{
			return new_fftw_plan<1, complex_t, scalar_t>{}(reinterpret_cast<fftw_complex*>(src), reinterpret_cast<scalar_t*>(target), dims, estimate);
		}
	};

	//! Specialization based on symphas::dft::new_fftw_plan.
	template<>
	struct new_fftw_plan<2, complex_t, scalar_t>
	{
		//! Creates a new complex to real FFTW plan.
		/*!
		 * Creates a plan based on transforming \p src to \p target.
		 * The direction of the complex to real plan is always backward.
		 *
		 * \param src The source array.
		 * \param target The target array.
		 * \param dims The dimension of the source array.
		 * \param estimate Whether or not to estimate the plan. If the plan
		 * is not estimated, the source and target arrays will be modified.
		 */
		fftw_plan operator()(fftw_complex* src, scalar_t* target, const len_type* dims, bool estimate = false);
		//! Forwards arbitrary arguments to create a specific FFTW plan.
		/*!
		 * Creates a plan based on transforming \p src to \p target.
		 *
		 * \param src The source array.
		 * \param target The target array.
		 * \param dims The dimension of the source array.
		 * \param estimate Whether or not to estimate the plan. If the plan
		 * is not estimated, the source and target arrays will be modified.
		 */
		template<typename T1, typename T2>
		fftw_plan operator()(T1* src, T2* target, const len_type* dims, bool estimate = false)
		{
			return new_fftw_plan<2, complex_t, scalar_t>{}(reinterpret_cast<fftw_complex*>(src), reinterpret_cast<scalar_t*>(target), dims, estimate);
		}
	};

	//! Specialization based on symphas::dft::new_fftw_plan.
	template<>
	struct new_fftw_plan<3, complex_t, scalar_t>
	{
		//! Creates a new complex to real FFTW plan.
		/*!
		 * Creates a plan based on transforming \p src to \p target.
		 * The direction of the complex to real plan is always backward.
		 *
		 * \param src The source array.
		 * \param target The target array.
		 * \param dims The dimension of the source array.
		 * \param estimate Whether or not to estimate the plan. If the plan
		 * is not estimated, the source and target arrays will be modified.
		 */
		fftw_plan operator()(fftw_complex* src, scalar_t* target, const len_type* dims, bool estimate = false);
		//! Forwards arbitrary arguments to create a specific FFTW plan.
		/*!
		 * Creates a plan based on transforming \p src to \p target.
		 *
		 * \param src The source array.
		 * \param target The target array.
		 * \param dims The dimension of the source array.
		 * \param estimate Whether or not to estimate the plan. If the plan
		 * is not estimated, the source and target arrays will be modified.
		 */
		template<typename T1, typename T2>
		fftw_plan operator()(T1* src, T2* target, const len_type* dims, bool estimate = false)
		{
			return new_fftw_plan<3, complex_t, scalar_t>{}(reinterpret_cast<fftw_complex*>(src), reinterpret_cast<scalar_t*>(target), dims, estimate);
		}
	};


	// **************************************************************************************

	//! Arrange FFTW data from sequential to in place format.
	/*!
	 * From source, copy the data to the target array. This arranges the data
	 * in target so that the FFTW transform for an in place r2c transform
	 * in target can produce the correct results.
	 * 
	 * This algorithm supports the case when src and target are the same array.
	 * 
	 * See arrange_fftw_ipts, which performs the opposite transposition.
	 *
	 * \param src The array which is arranged sequentially, and will be copied
	 * from into target.
	 * \param target The data here is copied from src according to an 
	 * arrangement that allows FFTW to perform an in place transform.
	 * \param dims The underlying dimensions of the complex data.
	 */
	template<size_t D>
	void arrange_fftw_stip(const scalar_t* src, scalar_t* target, const len_type* dims) = delete;

	//! Specialization based on symphas::dft::arrange_fftw_stip.
	template<>
	inline void arrange_fftw_stip<1>(const scalar_t* src, scalar_t* target, const len_type* dims)
	{
		iter_type ii = grid::length<1>(dims) - 1; //!< The index in the source array.
		iter_type c = 2 - dims[0] % 2; //!< The jump to the next row as according to the FFTW arrangement.
		iter_type n = (dims[0] + c) - 1; //!< Index in the transformed array.

		n -= c;

		for (iter_type i = 0; i < dims[0]; ++i)
		{
			target[n--] = src[ii--];
		}
	}

	//! Specialization based on symphas::dft::arrange_fftw_stip.
	template<>
	inline void arrange_fftw_stip<2>(const scalar_t* src, scalar_t* target, const len_type* dims)
	{
		iter_type ii = grid::length<2>(dims) - 1;
		iter_type c = 2 - dims[0] % 2;
		iter_type n = (dims[0] + c) * dims[1] - 1;

		for (iter_type j = 0; j < dims[1]; ++j)
		{
			n -= c;
			for (iter_type i = 0; i < dims[0]; ++i)
			{
				target[n--] = src[ii--];
			}
		}

	}

	//! Specialization based on symphas::dft::arrange_fftw_stip.
	template<>
	inline void arrange_fftw_stip<3>(const scalar_t* src, scalar_t* target, const len_type* dims)
	{
		iter_type ii = grid::length<3>(dims) - 1;
		iter_type c = 2 - dims[0] % 2;
		iter_type n = (dims[0] + c) * dims[1] * dims[2] - 1;

		for (iter_type k = 0; k < dims[2]; ++k)
		{
			for (iter_type j = 0; j < dims[1]; ++j)
			{
				n -= c;
				for (iter_type i = 0; i < dims[0]; ++i)
				{
					target[n--] = src[ii--];
				}
			}
		}
	}

	//! Arrange FFTW data from in place format to sequential format.
	/*!
	 * From source, copy the data to the target array. This arranges the data
	 * in target so that it is sequential, as it assumes that the source
	 * array is in the format given by an in place FFTW transform.
	 *
	 * This algorithm supports the case when src and target are the same array.
	 * 
	 * See arrange_fftw_stip, which performs the opposite transposition.
	 *
	 * \param src The array which is arranged according to the FFTW format for
	 * in place Fourier transforms.
	 * \param target The data here is copied from src into a sequential
	 * arrangement suitable for regular array use cases.
	 * \param dims The underlying dimensions of the complex data.
	 */
	template<size_t D>
	void arrange_fftw_ipts(scalar_t* src, scalar_t* target, const len_type* dims) = delete;

	//! Specialization based on symphas::dft::arrange_fftw_ipts.
	template<>
	inline void arrange_fftw_ipts<1>(scalar_t* src, scalar_t* target, const len_type* dims)
	{
		iter_type ii = 0;
		iter_type n = 0;

		for (iter_type i = 0; i < dims[0]; ++i)
		{
			target[ii++] = src[n++];
		}
	}

	//! Specialization based on symphas::dft::arrange_fftw_ipts.
	template<>
	inline void arrange_fftw_ipts<2>(scalar_t* src, scalar_t* target, const len_type* dims)
	{
		iter_type ii = 0;
		iter_type n = 0;
		iter_type c = 2 - dims[0] % 2;

		for (iter_type j = 0; j < dims[1]; ++j)
		{
			for (iter_type i = 0; i < dims[0]; ++i)
			{
				target[ii++] = src[n++];
			}
			n += c;
		}
	}

	//! Specialization based on symphas::dft::arrange_fftw_ipts.
	template<>
	inline void arrange_fftw_ipts<3>(scalar_t* src, scalar_t* target, const len_type* dims)
	{
		iter_type ii = 0;
		iter_type n = 0;
		iter_type c = 2 - dims[0] % 2;

		for (iter_type k = 0; k < dims[2]; ++k)
		{
			for (iter_type j = 0; j < dims[1]; ++j)
			{
				for (iter_type i = 0; i < dims[0]; ++i)
				{
					target[ii++] = src[n++];
				}
				n += c;
			}
		}
	}



	//! Arranges FFTW half complex output into a (full) sequential array.
	/*!
	 * When an FFTW r2c plan has been executed, the complex array storing the
	 * result will only store the non-redundant half. In some cases, the whole
	 * complex array should be represented, and this algorithm rearranges the
	 * values of the array to fit from the result into the desired usable shape.
	 * 
	 * The source and target array may be the same
	 *
	 * \param src The source array which is in half complex output.
	 * \param target The array that is rearranged so that it is expanded into
	 * the full dimensions.0
	 * \param dims The underlying dimensions of the data.
	 */
	template<size_t D>
	void arrange_fftw_hcts(complex_t* src, complex_t* target, const len_type* dims) = delete;

	//! Specialization based on symphas::dft::arrange_fftw_hcts.
	template<>
	inline void arrange_fftw_hcts<1>(complex_t* src, complex_t* target, const len_type* dims)
	{
		iter_type dn = dims[0] / 2 + 1;
#		pragma omp parallel for
		for (iter_type i = dn; i < dims[0]; ++i)
		{
			target[i] = std::conj(src[dn - i + dn - 2]);
		}
	}

	//! Specialization based on symphas::dft::arrange_fftw_hcts.
	template<>
	inline void arrange_fftw_hcts<2>(complex_t* src, complex_t* target, const len_type* dims)
	{
		iter_type dn = dims[0] / 2 + 1;
		iter_type ii = grid::length<2>(dims) - dims[0] + dn;
		iter_type n = dn * dims[1];

		/* Iterate over the hermitian half.
		 */
		for (iter_type j = 0; j < dims[1]; ++j)
		{
			for (iter_type i = 0; i < dn; ++i)
			{
				target[--ii] = std::conj(src[--n]);
			}
			ii -= dims[0] - dn;
		}

		n = dims[0] * (dims[1] - 1) + dn - 2;
		ii = dims[0] + dn;
		for (iter_type j = 1; j < dims[1]; ++j)
		{
			for (iter_type i = 0; i < dims[0] - dn; ++i)
			{
				target[ii++] = std::conj(target[n--]);
			}
			n -= dn;
			ii += dn;
		}
		ii = dn;
		for (iter_type i = 0; i < dims[0] - dn; ++i)
		{
			target[ii++] = std::conj(target[n--]);
		}
	}

	//! Specialization based on symphas::dft::arrange_fftw_hcts.
	template<>
	inline void arrange_fftw_hcts<3>(complex_t* src, complex_t* target, const len_type* dims)
	{
		iter_type dn = dims[0] / 2 + 1;
		iter_type ii = grid::length<3>(dims) - dims[0] + dn;
		iter_type n = dn * dims[1] * dims[2];

		for (iter_type k = 0; k < dims[2]; ++k)
		{
			for (iter_type j = 0; j < dims[1]; ++j)
			{
				for (iter_type i = 0; i < dn; ++i)
				{
					target[--ii] = std::conj(src[--n]);
				}
				ii -= dims[0] - dn;
			}
		}

		iter_type layer_offset = 0,
			layer_length = dims[0] * dims[1],
			last_layer_offset = layer_length * (dims[2] - 1);
		for (iter_type k = 1; k < dims[2]; ++k)
		{
			n = (last_layer_offset - layer_offset) + dims[0] * (dims[1] - 1) + dn - 2;
			ii = layer_length + layer_offset + dims[0] + dn;
			for (iter_type j = 1; j < dims[1]; ++j)
			{
				for (iter_type i = 0; i < dims[0] - dn; ++i)
				{
					target[ii++] = std::conj(target[n--]);
				}
				n -= dn;
				ii += dn;
			}
			ii = layer_length + layer_offset + dn;
			for (iter_type i = 0; i < dims[0] - dn; ++i)
			{
				target[ii++] = std::conj(target[n--]);
			}
			layer_offset += layer_length;
		}
		n = dims[0] * (dims[1] - 1) + dn - 2;
		ii = dims[0] + dn;
		for (iter_type j = 1; j < dims[1]; ++j)
		{
			for (iter_type i = 0; i < dims[0] - dn; ++i)
			{
				target[ii++] = std::conj(target[n--]);
			}
			n -= dn;
			ii += dn;
		}
		ii = dn;
		for (iter_type i = 0; i < dims[0] - dn; ++i)
		{
			target[ii++] = std::conj(target[n--]);
		}
	}


	//! Arranges a (full) sequential array into FFTW half complex output.
	/*!
	 * When an FFTW r2c plan has been executed, the complex array storing the
	 * result will only store the non-redundant half. In some cases, the whole
	 * complex array should be represented, and this algorithm rearranges the
	 * values of the array to fit from the result into the desired usable shape.
	 *
	 * The source and target array may be the same
	 *
	 * \param src The source array which is the source complex output.
	 * \param target The array that is rearranged so that it is expanded into
	 * the full dimensions.
	 * \param dims The underlying dimensions of the data.
	 */
	template<size_t D>
	void arrange_fftw_sthc(complex_t* src, complex_t* target, const len_type* dims) = delete;
	
	//! Specialization based on symphas::dft::arrange_fftw_sthc.
	template<>
	inline void arrange_fftw_sthc<1>(complex_t* src, complex_t* target, const len_type* dims) {}

	//! Specialization based on symphas::dft::arrange_fftw_sthc.
	template<>
	inline void arrange_fftw_sthc<2>(complex_t* src, complex_t* target, const len_type* dims)
	{
		iter_type dn = dims[0] / 2 + 1;
		iter_type ii = dn;

		for (iter_type j = 1; j < dims[1]; ++j)
		{
			iter_type n = dims[0] * j;
			for (iter_type i = 0; i < dn; ++i)
			{
				target[ii++] = src[n++];
			}
		}
	}

	//! Specialization based on symphas::dft::arrange_fftw_sthc.
	template<>
	inline void arrange_fftw_sthc<3>(complex_t* src, complex_t* target, const len_type* dims)
	{
		iter_type dn = dims[0] / 2 + 1;
		iter_type ii = dn;

		for (iter_type j = 1; j < dims[1]; ++j)
		{
			iter_type n = dims[0] * j;
			for (iter_type i = 0; i < dn; ++i)
			{
				target[ii++] = src[n++];
			}
		}
		for (iter_type k = 1; k < dims[2]; ++k)
		{
			for (iter_type j = 0; j < dims[1]; ++j)
			{
				iter_type n = dims[0] * j + dims[0] * dims[1] * k;
				for (iter_type i = 0; i < dn; ++i)
				{
					target[ii++] = src[n++];
				}
			}
		}
	}




	// **************************************************************************************


	namespace
	{


		//! Executes a function with consideration of an FFTW arrangement.
		/*!
		 * When using FFTW transformations between complex- and real-valued
		 * arrays, the data of the real-valued grid corresponds to only the
		 * non-duplicated half of the complex (Fourier space) array. Thus, an
		 * iteration must iterate to skip over the duplicated values when
		 * iterating over the complex array.
		 * 
		 * The template parameter specializes for either complex or real types,
		 * and when the real type is provided, the iteration will skip the
		 * duplicated part, and otherwise there will be no skipping.
		 * 
		 * The function given to the operator member is expected to have two
		 * indices, the first is the index in the complex array which will be
		 * skipped on duplicated rows when the given type is ::scalar_t.
		 * The second index is the index is always sequential.
		 * 
		 * If the iteration is for a complex type, it will be entirely
		 * sequential (both indices will always match).
		 * 
		 * \tparam T The type of the grid values.
		 * \tparam D The dimension of the grid.
		 */
		template<typename T, size_t D>
		struct iterate_complex_dup
		{
			template<typename F>
			void operator()(F f, const len_type* dim)
			{
				len_type len = grid::length<D>(dim);
#				pragma omp parallel for
				for (iter_type i = 0; i < len; ++i)
				{
					f(i, i);
				}
			}

			template<typename F>
			void operator()(T* values, F f, const len_type* dim)
			{
				len_type len = grid::length<D>(dim);
#				pragma omp parallel for
				for (iter_type i = 0; i < len; ++i)
				{
					values[i] = f(i);
				}
			}
		};

		//! Specialization of iterate_complex_dup.
		template<>
		struct iterate_complex_dup<scalar_t, 2>
		{
			template<typename F>
			void operator()(F&& f, const len_type* dim)
			{
				len_type Lh = dim[0] / 2 + 1;

#				pragma omp parallel for
				for (iter_type j = 0; j < dim[1]; ++j)
				{
					iter_type skip_index = j * dim[0];	// Iterated across the whole dimension.
					iter_type seq_index = j * Lh;		// Iterated sequentially.
					for (iter_type i = 0; i < Lh; ++i)
					{
						f(skip_index + i, seq_index + i);
					}
				}
			}

			template<typename F>
			void operator()(scalar_t* values, F&& f, const len_type* dim)
			{
				const len_type Lh = dim[0] / 2 + 1;

#				pragma omp parallel for shared(f)
				for (iter_type j = 0; j < dim[1]; ++j)
				{
					iter_type skip_index = j * dim[0];
					iter_type seq_index = j * Lh;
					for (iter_type i = 0; i < Lh; ++i)
					{
						values[seq_index + i] = f(skip_index + i);
					}
				}
			}
		};

		//! Specialization of iterate_complex_dup.
		template<>
		struct iterate_complex_dup<scalar_t, 3>
		{
			template<typename F>
			void operator()(F&& f, const len_type* dim)
			{
				len_type Lh = dim[0] / 2 + 1;

#				pragma omp parallel for
				for (iter_type kj = 0; kj < dim[1] * dim[2]; ++kj)
				{
					iter_type skip_index = kj * dim[0];
					iter_type seq_index = kj * Lh;
					for (iter_type i = 0; i < Lh; ++i)
					{
						f(skip_index + i, seq_index + i);
					}
				}
			}

			template<typename F>
			void operator()(scalar_t* values, F&& f, const len_type* dim)
			{
				len_type Lh = dim[0] / 2 + 1;

#				pragma omp parallel for
				for (iter_type kj = 0; kj < dim[1] * dim[2]; ++kj)
				{
					iter_type skip_index = kj * dim[0];
					iter_type seq_index = kj * Lh;
					for (iter_type i = 0; i < Lh; ++i)
					{
						values[seq_index + i] = f(skip_index + i);
					}
				}
			}
		};
	}

	//! Computes the number of array elements for the given type.
	/*!
	 * Using the FFTW algorithm means that the number of elements in the 
	 * complex array may be shortened since the transform will be duplicated
	 * by way of Hermitian conjugate.
	 * 
	 * When the given type is complex, the true length of the data array
	 * represented by the given dimensions is returned. When the given type
	 * is real, then only the length of the non-duplicated half is returned.
	 * 
	 * \param dims The dimensions of the grid which is transformed.
	 * 
	 * \tparam T The type which the length is adjusted to. A real type will
	 * return a length half of the total dimensions.
	 * \tparam D the dimension of the corresponding dimension.
	 */
	template<typename T, size_t D>
	len_type len(len_type const* dims) = delete;


	//! Specialization of symphas::dft::len().
	template<>
	inline len_type len<scalar_t, 1>(len_type const* dims)
	{
		return (dims[0] / 2 + 1);
	}

	//! Specialization of symphas::dft::len().
	template<>
	inline len_type len<scalar_t, 2>(len_type const* dims)
	{
		return (dims[0] / 2 + 1) * dims[1];
	}

	//! Specialization of symphas::dft::len().
	template<>
	inline len_type len<scalar_t, 3>(len_type const* dims)
	{
		return (dims[0] / 2 + 1) * dims[1] * dims[2];
	}

	//! Specialization of symphas::dft::len().
	template<>
	inline len_type len<complex_t, 1>(len_type const* dims)
	{
		return grid::length<1>(dims);
	}

	//! Specialization of symphas::dft::len().
	template<>
	inline len_type len<complex_t, 2>(len_type const* dims)
	{
		return grid::length<2>(dims);
	}

	//! Specialization of symphas::dft::len().
	template<>
	inline len_type len<complex_t, 3>(len_type const* dims)
	{
		return grid::length<3>(dims);
	}




	//! Applies the binary function over matching indices.
	/*!
	 * See symphas::dft::iterate_complex_dup. Iterates through by skipping over
	 * duplicated data in each row, which is chosen by the source type.
	 */
	template<typename T, size_t D, typename F>
	void iterate_rc(F&& f, len_type const* dims)
	{
		iterate_complex_dup<T, D>{}(std::forward<F>(f), dims);
	}

	//! Applies the binary function over matching indices.
	/*!
	 * See symphas::dft::iterate_complex_dup. Iterates through by skipping over
	 * duplicated data in each row, which is chosen by the source type.
	 */
	template<typename T, size_t D, typename F>
	void iterate_rc(T* values, F&& f, len_type const* dims)
	{
		iterate_complex_dup<T, D>{}(values, std::forward<F>(f), dims);
	}

	//! Executes the given plan to compute the Fourier transform.
	/*!
	 * An overload that takes an instance of the FFTW plan type and executes it
	 * to compute the Fourier transform.
	 * 
	 * \param p An instance of the FFTW plan type.
	 */
	inline void dft(fftw_plan p)
	{
		symphas::dft::fftw_execute(p);
	}
	
	//! Compute the Fouier transform of a `D`-dimensional system.
	/*!
	 * Compute the Fourier transform of the given real-valued system and store
	 * the result in the given complex array. The given data is a sequential
	 * array representing a `D`-dimensional grid with dimensions \p dims. The
	 * input and output both have the same dimensions.
	 * 
	 * The FFTW routine is used by creating an estimated plan, and then
	 * executed. Since the input is real-valued, the result then has to be
	 * rearranged to represent the whole system that will include duplicated
	 * values. The plan will always do a forward type transformation. The plan 
	 * is then destroyed.
	 * 
	 * \param[in] in The input to the Fourier transform, as real-valued data.
	 * \param[out] out The output array of the Fourier transform, where the
	 * result is written to. This has dimensions equal to the input.
	 * \param dims The dimensions of the input array of length `D`.
	 * 
	 * \tparam D The dimension of the input system.
	 */
	template<size_t D>
	void dft(const scalar_t* in, complex_t* out, const len_type* dims, bool)
	{
		fftw_plan p = new_fftw_plan<D, scalar_t, complex_t>{}(const_cast<scalar_t*>(in), out, dims, true);
		symphas::dft::fftw_execute(p);
		arrange_fftw_hcts<D>(out, out, dims);
		symphas::dft::fftw_destroy_plan(p);
	}

	//! Compute the Fouier transform of a `D`-dimensional system.
	/*!
	 * Compute the Fourier transform of the given complex-valued system and 
	 * store the result in the given complex array. The given data is a 
	 * sequential array representing a `D`-dimensional grid with dimensions 
	 * \p dims. The input and output both have the same dimensions.
	 *
	 * The FFTW routine is used by creating an estimated plan, and then
	 * executed. The plan is then destroyed.
	 *
	 * \param[in] in The input to the Fourier transform.
	 * \param[out] out The output array of the Fourier transform, where the
	 * result is written to. This has dimensions equal to the input.
	 * \param dims The dimensions of the input array of length `D`.
	 *
	 * \tparam D The dimension of the input system.
	 */
	template<size_t D>
	void dft(const complex_t* in, complex_t* out, const len_type* dims, bool backward = false)
	{
		fftw_plan p = new_fftw_plan<D, complex_t, complex_t>{}(const_cast<complex_t*>(in), out, dims, true, backward);
		symphas::dft::fftw_execute(p);
		symphas::dft::fftw_destroy_plan(p);
	}

	//! Compute the Fouier transform of a `D`-dimensional system.
	/*!
	 * Compute the Fourier transform of the given complex-valued system and 
	 * store the result in the given real-valued array. The given data is a 
	 * sequential array representing a `D`-dimensional grid with dimensions 
	 * \p dims. The input and output both have the same dimensions.
	 *
	 * The FFTW routine is used by creating an estimated plan, and then
	 * executed. Since the output is real-valued, the input is first rearranged
	 * because this type of transformation assumes that the input array has
	 * values which repeat in the two halves. Once the plan is executed, the
	 * input is restored to the original arrangement. Therefore, this plan
	 * always assumes a backward type transformation. That is, it is the
	 * opposite of symphas::dft::dft(const scalar_t*, complex_t*, const len_type*).
	 * The plan is subsequently destroyed.
	 * 
	 * **NOTE**: If the input array does not have repeated values between each
	 * half, then this routine will destroy the original data. 
	 
	 * \param[in] in The input to the Fourier transform, as complex-valued data.
	 * \param[out] out The output array of the Fourier transform, where the
	 * result is written to. This has dimensions equal to the input.
	 * \param dims The dimensions of the input array of length `D`.
	 *
	 * \tparam D The dimension of the input system.
	 */
	template<size_t D>
	void dft(const complex_t* in, scalar_t* out, const len_type* dims, bool)
	{
		complex_t* hf = const_cast<complex_t*>(in);
		arrange_fftw_sthc<D>(hf, hf, dims);
		fftw_plan p = new_fftw_plan<D, complex_t, scalar_t>{}(hf, out, dims, true);
		symphas::dft::fftw_execute(p);
		symphas::dft::fftw_destroy_plan(p);
		arrange_fftw_hcts<D>(hf, hf, dims);
	}

	//! An arbitrary Fourier transform chosen based on the argument types.
	/*!
	 * The type of transform is chosen based on the argument types, which
	 * must be compatible with the ::complex_t and ::scalar_t types.
	 * 
	 * \param[in] in The data which is transformed.
	 * \param[out] out The result of the Fourier transform of \p in.
	 * \param dims The dimensions of the input array.
	 * 
	 * \tparam D The dimension of the system.
	 */
	template<size_t D, typename T, typename C>
	void dft(T&& in, C&& out, std::initializer_list<len_type> dims, bool backward = false)
	{
		dft<D>(std::forward<T>(in), std::forward<C>(out), dims.begin(), backward);
	}

	//! An arbitrary Fourier transform chosen based on the argument types.
	/*!
	 * The type of transform is chosen based on the argument types, which
	 * must be compatible with the ::complex_t and ::scalar_t types. This
	 * transform takes the input system to be a 1-dimensional system.
	 *
	 * \param[in] in The data which is transformed.
	 * \param[out] out The result of the Fourier transform of \p in.
	 * \param L The length of the input array.
	 */
	template<typename T, typename C>
	void dft(T&& in, C&& out, len_type L, bool backward = false)
	{
		dft<1>(std::forward<T>(in), std::forward<C>(out), &L, backward);
	}

	//! An arbitrary Fourier transform chosen based on the argument types.
	/*!
	 * The type of transform is chosen based on the argument types, which
	 * must be compatible with the ::complex_t and ::scalar_t types. This
	 * transform takes the input system to be a 2-dimensional system.
	 *
	 * \param[in] in The data which is transformed.
	 * \param[out] out The result of the Fourier transform of \p in.
	 * \param L The length of the \f$x\f$-axis.
	 * \param M The length of the \f$y\f$-axis.
	 */
	template<typename T, typename C>
	void dft(T&& in, C&& out, len_type L, len_type M, bool backward = false)
	{
		dft<2>(std::forward<T>(in), std::forward<C>(out), { L, M }, backward);
	}

	//! An arbitrary Fourier transform chosen based on the argument types.
	/*!
	 * The type of transform is chosen based on the argument types, which
	 * must be compatible with the ::complex_t and ::scalar_t types. This
	 * transform takes the input system to be a 3-dimensional system.
	 *
	 * \param[in] in The data which is transformed.
	 * \param[out] out The result of the Fourier transform of \p in.
	 * \param L The length of the \f$x\f$-axis.
	 * \param M The length of the \f$y\f$-axis.
	 * \param N The length of the \f$z\f$-axis.
	 */
	template<typename T, typename C>
	void dft(T&& in, C&& out, len_type L, len_type M, len_type N, bool backward = false)
	{
		dft<3>(std::forward<T>(in), std::forward<C>(out), { L, M, N }, backward);
	}


}









