
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
 * PURPOSE: Defines almost all the rules present in the symbolic algebra
 * functionality.
 * 
 * ***************************************************************************
 */

#pragma once

#include <cassert>

#include "expressionlogic.h"

namespace expr
{
	//! Specialization based on SymbolID.
	template<typename A>
	struct SymbolID<symphas::ref<A>>
	{
		static decltype(auto) get(symphas::ref<A> const& a)
		{
			return SymbolID<A>::get(a);
		}
	};

	//! Specialization based on SymbolID.
	template<size_t Z, typename G>
	struct SymbolID<Variable<Z, G>>
	{
		static decltype(auto) get(Variable<Z, G> const& a)
		{
			return SymbolID<G>::get(a);
		}
	};



	//! Specialization based on BaseData.
	template<typename A>
	struct BaseData<symphas::ref<A>>
	{
		static decltype(auto) get(symphas::ref<A> const& data, iter_type n)
		{
			return BaseData<A>::get(data, n);
		}
		static decltype(auto) get(symphas::ref<A> const& data)
		{
			return BaseData<A>::get(data);
		}
		static decltype(auto) get(symphas::ref<A>& data, iter_type n)
		{
			return BaseData<A>::get(data, n);
		}
		static decltype(auto) get(symphas::ref<A>& data)
		{
			return BaseData<A>::get(data);
		}
	};

	//! Specialization based on BaseData.
	template<size_t Z, typename T>
	struct BaseData<Variable<Z, T>>
	{
		static decltype(auto) get(Variable<Z, T> const& data, iter_type n)
		{
			return BaseData<T>::get(data, n);
		}
		static decltype(auto) get(Variable<Z, T> const& data)
		{
			return BaseData<T>::get(data);
		}
		static decltype(auto) get(Variable<Z, T>& data, iter_type n)
		{
			return BaseData<T>::get(data, n);
		}
		static decltype(auto) get(Variable<Z, T>& data)
		{
			return BaseData<T>::get(data);
		}
	};
}


//! Specialization based on expr::base_data_type.
template<typename A>
struct expr::base_data_type<symphas::ref<A>>
{
	using type = typename expr::base_data_type<A>::type;
};

//! Specialization based on expr::base_data_type.
template<size_t Z, typename T>
struct expr::base_data_type<Variable<Z, T>>
{
	using type = typename expr::base_data_type<T>::type;
};

//! Specialization based on expr::base_data_type.
template<typename A>
struct expr::base_data_type<const A>
{
	using type = typename expr::base_data_type<A>::type;
};

//! Specialization based on expr::base_data_type.
template<typename A>
struct expr::base_data_type<A&>
{
	using type = typename expr::base_data_type<A>::type;
};



// ******************************************************************************************


namespace expr
{
	//! Indicates whether a term can be collected into like terms.
	/*!
	 * Predicate indicating whether a data type can combine with another data type.
	 * Combining means that two data are considered "equivalent" and can be put
	 * together into like terms.
	 *
	 * Grid combination requires a level of indirection in order
	 * to allow variables of all types to be combined. This means that some kind
	 * of subtype needs to be compared, as the encapsulating types can vary.
	 * By default, variables are allowed to be combined into like terms.
	 *
	 * \tparam E The type which is checked if it satisfies the predicate.
	 */
	template<typename E, typename Enable = void>
	constexpr bool is_combinable = grid_can_combine<typename original_data_type<E>::type>::value;

	//! Specialization based on expr::is_combinable.
	template<size_t Z, typename G>
	constexpr bool is_combinable<Variable<Z, G>> = true;


	//! Indicates whether a term can be multiplied.
	/*!
	 * Predicate indicating whether a data type can combine with another data 
	 * type through multiplication.
	 * Combining means that two data are considered "equivalent" and can be put
	 * together into a nonlinear variable.
	 *
	 * It may be desirable to ensure that some terms cannot be multiplied
	 * into nonlinear variables, in which case, only a OpBinaryMul will be
	 * used to combine them.
	 *
	 * \tparam E The type which is checked if it satisfies the predicate.
	 */
	template<typename E>
	constexpr bool is_multiplies = expr::grid_can_multiply<typename expr::base_data_type<E>::type>::value;

}

namespace symphas::internal
{

	//! Recursively checks for the correct type to remove.
	/*!
	 * The specified template type `G0` is compared to the first type that appears
	 * in the second list. If it matches, then a new list is returned as the
	 * concatenation of the first and second without the matched type. If the type
	 * does not match, then the first type is moved into the first list and the
	 * check is performed recursively.
	 *
	 * \tparam G0 The first type that is removed.
	 */
	template<typename G0, typename Gn1, typename... Gcs, typename... Gns,
		typename std::enable_if_t<!std::is_same<G0, Gn1>::value && (sizeof...(Gns) == 0), int> = 0>
		constexpr auto subtract_common_type(std::tuple<Gcs*...>, std::tuple<Gn1*, Gns*...>)
	{
		return std::tuple<>{};
	}

	template<typename G0, typename Gn1, typename... Gcs, typename... Gns, 
		typename std::enable_if_t<std::is_same<G0, Gn1>::value, int> = 0>
	constexpr auto subtract_common_type(std::tuple<Gcs*...>, std::tuple<Gn1*, Gns*...>)
	{
		return std::tuple<Gcs*..., Gns*...>{};
	}

	template<typename G0, typename Gn1, typename... Gcs, typename... Gns, 
		typename std::enable_if_t<!std::is_same<G0, Gn1>::value && (sizeof...(Gns) > 0), int> = 0>
	constexpr auto subtract_common_type(std::tuple<Gcs*...>, std::tuple<Gn1*, Gns*...>)
	{
		return subtract_common_type<G0>(std::tuple_cat(std::tuple<Gcs*...>{}, std::tuple<Gn1*>{}), std::tuple<Gns*...>{});
	}


	//! Removes the specified type from the type list.
	/*!
	 * The given template parameter type is recursively searched for in the given
	 * list and the list is returned with the type removed. If the type is not
	 * found, the same list is returned.
	 *
	 * \param list The list of types from which the given template parameter type
	 * is taken away from.
	 *
	 * \tparam G0 The type that is removed from the list.
	 */
	template<typename G0, typename Gn, typename... Gcs, typename std::enable_if_t<std::is_same<G0, Gn>::value, int> = 0>
	constexpr auto subtract_common_type([[maybe_unused]] std::tuple<Gn*, Gcs*...> list)
	{
		return std::tuple<Gcs*...>{};
	}

	template<typename G0, typename Gn, typename... Gcs, typename std::enable_if_t<!std::is_same<G0, Gn>::value, int> = 0>
	constexpr auto subtract_common_type([[maybe_unused]] std::tuple<Gn*, Gcs*...> list)
	{
		return subtract_common_type<G0>(std::tuple<Gn*>{}, std::tuple<Gcs*...>{});
	}

}


namespace expr
{
	//! Compares whether two multi-variables can be combined as like terms.
	/*!
	 * Implements a compile time algorithm that compares the types of two lists.
	 * Like types are eliminated by calling the function subtract_common_type
	 * and checks whether they all satisfy grid_can_combine.
	 *
	 * If there is no common type, then subtract_common_type returns an empty
	 * tuple and the result is that the two multi-variables cannot be combined.
	 *
	 * G0 is compared to all the types that are template parameters from the
	 * nested struct, and the common type is eliminated from the list using
	 * recursion. The result of one iteration is a new nl_can_combine, which
	 * repeats the procedure.
	 */
	template<typename... Gs>
	struct nl_can_combine;
}

namespace symphas::internal
{


	template<typename... Gs, typename g0, typename... gns>
	constexpr static bool next_type_check(std::tuple<Gs*...>, std::tuple<g0*, gns*...>)
	{
		static_assert(sizeof...(gns) + 1 == sizeof...(Gs));
		return expr::nl_can_combine<Gs...>::template with<g0, gns...>::value;
	}

	template<typename... Gs>
	constexpr static bool next_type_check(std::tuple<Gs*...>, std::tuple<>)
	{
		return false;
	}

	template<typename G0, typename... Gs, typename g0, typename g1, typename... gs, 
		typename std::enable_if_t<(1 + sizeof...(Gs) == 2 + sizeof...(gs)), int> = 0>
	constexpr bool nlv_combines_pred(std::tuple<G0*, Gs*...>, std::tuple<g0*, g1*, gs*...>)
	{
		return expr::is_combinable<G0>
			&& next_type_check(std::tuple<Gs*...>{}, symphas::internal::subtract_common_type<G0>(std::tuple<g0*, g1*, gs*...>{}));
	}

	template<typename... Gs, typename... gs,
		typename std::enable_if_t<(sizeof...(Gs) != sizeof...(gs)), int> = 0>
	constexpr bool nlv_combines_pred(std::tuple<Gs*...>, std::tuple<gs*...>)
	{
		return false;
	}

	template<typename... Gs>
	struct nlv_pred_outer
	{
		template<typename... gs>
		struct nlv_pred_inner
		{
			static const bool value = nlv_combines_pred(
				std::tuple<Gs*...>{},
				std::tuple<gs*...>{});
		};
	};

}

namespace expr
{

	//! Specialization based on expr::nl_can_combine.
	template<typename G0>
	struct nl_can_combine<G0>
	{
		template<typename g>
		struct with
		{
			static const bool value = expr::is_combinable<G0> && std::is_same<G0, g>::value;
		};
	};

	//! Specialization based on expr::nl_can_combine.
	template<typename G0, typename... Gs>
	struct nl_can_combine<G0, Gs...>
	{
		template<typename g0, typename g1, typename... gs>
		struct with
		{
			static const bool value = symphas::internal::nlv_pred_outer<G0, Gs...>::template nlv_pred_inner<g0, g1, gs...>::value;
		};

	};


	//! Compares whether two multi-variables can be multiplied.
	/*!
	 * Compares the types of two lists to check if two multi-variables can be
	 * multiplied together and combined into a single multi-variable.
	 *
	 * Like types are eliminated by calling the function subtract_common_type.
	 * If there is no common type, then subtract_common_type returns an empty
	 * tuple, and the result is that the two multi-variables cannot be
	 * multiplied together.
	 *
	 * The first type is compared to all the types that are template parameters to the 
	 * nested `with' type, and the common type is eliminated from the list.
	 * Recursion is performed, creating a new expr::nl_can_multiply with the 
	 * remaining types.
	 * 
	 * \param Gs... The list of types that are checked against another list, 
	 * to see if is possible to come the lists as a new multivariable.
	 */
	template<typename... Gs>
	struct nl_can_multiply
	{
		//! Specifies the list of variables to compare with the first.
		/*!
		 * Specifies the list of variables to compare with the first.
		 * 
		 * \tparam gs... Given list to compare with.
		 */
		template<typename... gs>
		struct with;
	};
}



/* an object can be multiplied into others if it satisfies one of two conditions
 * 1) all the objects in the test object support multiplication AND the subject supports multiplication
 * 2) the subject and the test objects are the same
 */

 //! Specialization based on expr::nl_can_multiply.
template<typename... Gs>
template<typename g>
struct expr::nl_can_multiply<Gs...>::with<g>
{
	static const bool value =
		((expr::is_multiplies<Gs> && ... && expr::is_multiplies<g>)) ||
		((std::is_same<typename expr::base_data_type<Gs>::type, typename expr::base_data_type<g>::type>::value && ...));
};

//! Specialization based on expr::nl_can_multiply.
template<typename... Gs>
template<typename g0, typename g1, typename... gs>
struct expr::nl_can_multiply<Gs...>::with<g0, g1, gs...>
{
protected:

	constexpr static bool get_value()
	{
		return expr::nl_can_multiply<Gs...>::template with<g0>::value 
			&& expr::nl_can_multiply<Gs...>::template with<g1, gs...>::value;
	}

public:

	static const bool value = get_value();

};


// ******************************************************************************************



namespace expr
{

	//! Apply an inverse to a scalar value.
	inline auto inverse(scalar_t e)
	{
		return symphas::lib::get_identity<scalar_t>() / e;
	}

	//! Apply an inverse to an integer value.
	inline auto inverse(int e)
	{
		return symphas::lib::get_identity<scalar_t>() / e;
	}

	//! Apply an inverse to a complex value.
	inline auto inverse(complex_t const& e)
	{
		return symphas::lib::get_identity<complex_t>() / e;
	}

	//! Apply an inverse to an expression.
	template<typename E>
	auto inverse(OpExpression<E> const& e)
	{
		return expr::make_div(OpIdentity{}, (*static_cast<const E*>(&e)));
	}

	//! Apply an inverse to an expression literal.
	template<typename T>
	auto inverse(OpLiteral<T> const& e)
	{
		return expr::make_literal(inverse(e.value));
	}

}





// ******************************************************************************************


template<typename A1, typename A2, typename B1, typename B2>
auto operator*(OpBinarySub<A1, A2> const& a, OpBinarySub<B1, B2> const& b);
template<typename A1, typename A2, typename B1, typename B2>
auto operator*(OpBinaryAdd<A1, A2> const& a, OpBinarySub<B1, B2> const& b);
template<typename A1, typename A2, typename B1, typename B2>
auto operator*(OpBinarySub<A1, A2> const& a, OpBinaryAdd<B1, B2> const& b);
template<typename A1, typename A2, typename B1, typename B2>
auto operator*(OpBinaryAdd<A1, A2> const& a, OpBinaryAdd<B1, B2> const& b);
template<typename A1, typename B1, typename B2>
auto operator*(OpExpression<A1> const& a, OpBinaryAdd<B1, B2> const& b);
template<typename B2, typename A1, typename A2>
auto operator*(OpBinaryAdd<A1, A2> const& a, OpExpression<B2> const& b);
template<typename A1, typename B1, typename B2>
auto operator*(OpExpression<A1> const& a, OpBinarySub<B1, B2> const& b);
template<typename B2, typename A1, typename A2>
auto operator*(OpBinarySub<A1, A2> const& a, OpExpression<B2> const& b);
template<typename T, typename B1, typename B2>
auto operator*(OpLiteral<T> const& a, OpBinaryMul<B1, B2> const& b);
template<typename T, typename B1, typename B2>
auto operator*(OpLiteral<T> const& a, OpBinaryDiv<B1, B2> const& b);


template<typename A1, typename A2, typename B1, typename B2>
auto operator+(OpBinaryAdd<A1, A2> const& a, OpBinaryAdd<B1, B2> const& b);
template<typename A1, typename A2, typename B1, typename B2>
auto operator+(OpBinaryAdd<A1, A2> const& a, OpBinarySub<B1, B2> const& b);
template<typename A1, typename A2, typename B1, typename B2>
auto operator+(OpBinarySub<A1, A2> const& a, OpBinaryAdd<B1, B2> const& b);
template<typename A1, typename A2, typename B1, typename B2>
auto operator+(OpBinarySub<A1, A2> const& a, OpBinarySub<B1, B2> const& b);
template<typename A1, typename A2, typename B1, typename B2>
auto operator-(OpBinaryAdd<A1, A2> const& a, OpBinaryAdd<B1, B2> const& b);
template<typename A1, typename A2, typename B1, typename B2>
auto operator-(OpBinaryAdd<A1, A2> const& a, OpBinarySub<B1, B2> const& b);
template<typename A1, typename A2, typename B1, typename B2>
auto operator-(OpBinarySub<A1, A2> const& a, OpBinaryAdd<B1, B2> const& b);
template<typename A1, typename A2, typename B1, typename B2>
auto operator-(OpBinarySub<A1, A2> const& a, OpBinarySub<B1, B2> const& b);

template<typename A1, typename A2, typename E>
auto operator+(OpExpression<E> const& a, OpBinaryAdd<A1, A2> const& b);
template<typename A1, typename A2, typename E>
auto operator+(OpBinaryAdd<A1, A2> const& a, OpExpression<E> const& b);
template<typename A1, typename A2, typename E>
auto operator+(OpExpression<E> const& a, OpBinarySub<A1, A2> const& b);
template<typename A1, typename A2, typename E>
auto operator+(OpBinarySub<A1, A2> const& a, OpExpression<E> const& b);
template<typename A1, typename A2, typename E>
auto operator-(OpExpression<E> const& a, OpBinaryAdd<A1, A2> const& b);
template<typename A1, typename A2, typename E>
auto operator-(OpBinaryAdd<A1, A2> const& a, OpExpression<E> const& b);
template<typename A1, typename A2, typename E>
auto operator-(OpExpression<E> const& a, OpBinarySub<A1, A2> const& b);
template<typename A1, typename A2, typename E>
auto operator-(OpBinarySub<A1, A2> const& a, OpExpression<E> const& b);




template<typename A1, typename A2, typename T, typename G>
auto operator*(OpBinaryMul<A1, A2> const& a, OpLVariable<T, G> const& b);
template<typename T, typename G, typename B1, typename B2>
auto operator*(OpLVariable<T, G> const& a, OpBinaryMul<B1, B2> const& b);
template<typename A1, typename A2, typename T, typename... Gs>
auto operator*(OpBinaryMul<A1, A2> const& a, OpNLVariable<T, Gs...> const& b);
template<typename T, typename... Gs, typename B1, typename B2>
auto operator*(OpNLVariable<T, Gs...> const& a, OpBinaryMul<B1, B2> const& b);
template<typename A1, typename A2, typename B1, typename B2>
auto operator*(OpBinaryMul<A1, A2> const& a, OpBinaryMul<B1, B2> const& b);


template<typename A1, typename A2, typename E>
auto operator*(OpBinaryDiv<A1, A2> const& a, OpExpression<E> const& b);
template<typename E, typename B1, typename B2>
auto operator*(OpExpression<E> const& a, OpBinaryDiv<B1, B2> const& b);
template<typename A1, typename A2, typename B1, typename B2>
auto operator*(OpBinaryDiv<A1, A2> const& a, OpBinaryDiv<B1, B2> const& b);
template<typename A1, typename A2, typename B1, typename B2>
auto operator*(OpBinaryDiv<A1, A2> const& a, OpBinaryAdd<B1, B2> const& b);
template<typename A1, typename A2, typename B1, typename B2>
auto operator*(OpBinaryDiv<A1, A2> const& a, OpBinarySub<B1, B2> const& b);
template<typename A1, typename A2, typename B1, typename B2>
auto operator*(OpBinaryAdd<A1, A2> const& a, OpBinaryDiv<B1, B2> const& b);
template<typename A1, typename A2, typename B1, typename B2>
auto operator*(OpBinarySub<A1, A2> const& a, OpBinaryDiv<B1, B2> const& b);





/* Handling recursion at the point of an inverse.
 */

template<typename A2, typename B1, typename B2>
auto operator*(OpBinaryDiv<OpNegIdentity, A2> const& a, OpBinaryDiv<B1, B2> const& b);
template<typename A2, typename B1, typename B2>
auto operator*(OpBinaryDiv<OpIdentity, A2> const& a, OpBinaryDiv<B1, B2> const& b);
template<typename T, typename A2, typename B1, typename B2>
auto operator*(OpBinaryDiv<OpLiteral<T>, A2> const& a, OpBinaryDiv<B1, B2> const& b);
template<typename A1, typename A2, typename B2>
auto operator*(OpBinaryDiv<A1, A2> const& a, OpBinaryDiv<OpNegIdentity, B2> const& b);
template<typename A1, typename A2, typename B2>
auto operator*(OpBinaryDiv<A1, A2> const& a, OpBinaryDiv<OpIdentity, B2> const& b);
template<typename A1, typename A2, typename T, typename B2>
auto operator*(OpBinaryDiv<A1, A2> const& a, OpBinaryDiv<OpLiteral<T>, B2> const& b);

template<typename V1, typename V2, typename Dd, typename E, typename Sp>
auto operator/(OpFuncDerivative<Dd, V1, E, Sp> const& a, OpFuncDerivative<Dd, V2, E, Sp> const& b);










/*
 *
 * Multiplcation between basic values, particularly identities.
 *
 ******************************************************************************/

 //! Ensure a value constant is always multiplied on the left by anything else.
template<typename E, typename T>
auto operator*(OpExpression<E> const& a, OpLiteral<T> const& b)
{
	return b * (*static_cast<const E*>(&a));
}


//! Multiplication between integer and complex types.
inline auto operator*(int const a, complex_t const& b)
{
	return static_cast<double>(a) * b;
}

//! Multiplication between complex and integer types.
inline auto operator*(complex_t const& a, int const b)
{
	return a * static_cast<double>(b);
}

//! Multiplication between the identity expression and a scalar.
inline auto operator*(OpIdentity const, scalar_t const& b)
{
	return b;
}

//! Multiplication between the identity expression and a complex value.
inline auto operator*(OpIdentity const, complex_t const& b)
{
	return b;
}

//! Multiplication between the identity expression and an integer value.
inline auto operator*(OpIdentity const, int b)
{
	return b;
}

//! Multiplication between the identity expression and a scalar.
inline auto operator*(scalar_t const& a, OpIdentity const)
{
	return a;
}

//! Multiplication between the identity expression and a complex value.
inline auto operator*(complex_t const& a, OpIdentity const)
{
	return a;
}

//! Multiplication between the identity expression and an integer value.
inline auto operator*(int a, OpIdentity const)
{
	return a;
}

//! Multiplication between the identity expression and anything.
template<typename T>
auto operator*(OpIdentity const, OpLiteral<T> const& b)
{
	return b;
}

//! Multiplication between anything and the identity expression.
template<typename T>
auto operator*(OpLiteral<T> const& a, OpIdentity)
{
	return a;
}

//! Multiplication between the identity expression and anything.
template<typename E>
decltype(auto) operator*(OpIdentity, E const& b)
{
	return b;
}

//! Multiplication between anything and the identity expression.
template<typename E>
decltype(auto) operator*(E const& a, OpIdentity)
{
	return OpIdentity{} * a;
}

//! Multiplication between the negative identity expression and anything.
template<typename E>
decltype(auto) operator*(OpNegIdentity const, E const& b)
{
	return -(OpIdentity{} * b);
}

//! Multiplication between the negative identity expression and anything.
template<typename E>
decltype(auto) operator*(E const& a, OpNegIdentity const)
{
	return (OpNegIdentity{} * a);
}

//! Multiplication between two identity expressions.
inline auto operator*(OpIdentity const, OpIdentity const)
{
	return OpIdentity{};
}

//! Multiplication between two negative identity expressions.
inline auto operator*(OpNegIdentity const, OpNegIdentity const)
{
	return OpIdentity{};
}

//! Multiplication between the negative identity and identity expressions.
inline auto operator*(OpNegIdentity const, OpIdentity const)
{
	return OpNegIdentity{};
}

//! Multiplication between the identity and negative identity expressions.
inline auto operator*(OpIdentity const, OpNegIdentity const)
{
	return OpNegIdentity{};
}

//! Multiplication between anything and the 0 identity.
template<typename E>
decltype(auto) operator*(E const&, OpVoid const)
{
	return OpVoid{};
}

//! Multiplication between the 0 identity and anything.
template<typename E>
decltype(auto) operator*(OpVoid const, E const&)
{
	return OpVoid{};
}

//! Multiplication between the identity expression and 0 identity.
inline auto operator*(OpIdentity const, OpVoid const)
{
	return OpVoid{};
}

//! Multiplication between the 0 identity and identity expression.
inline auto operator*(OpVoid const, OpIdentity const)
{
	return OpVoid{};
}

//! Multiplication between the negative identity and 0 identity.
inline auto operator*(OpNegIdentity const, OpVoid const)
{
	return OpVoid{};
}

//! Multiplication between the 0 identity and negative identity.
inline auto operator*(OpVoid const, OpNegIdentity const)
{
	return OpVoid{};
}

//! Multiplication between two 0 identities.
inline auto operator*(OpVoid const, OpVoid const)
{
	return OpVoid{};
}





/*
 *
 * Division between basic values, particularly identities.
 *
 ******************************************************************************/

 //! Division between anything and the identity expression.
template<typename E>
decltype(auto) operator/(E const& a, OpIdentity const)
{
	return a;
}

//! Division between the identity expression and anything.
template<typename E>
decltype(auto) operator/(OpIdentity const, E const& b)
{
	return expr::inverse(b);
}


//! Division between the and anything and the negative identity.
template<typename E>
decltype(auto) operator/(E const& a, OpNegIdentity const)
{
	return -a;
}

//! Division between the negative identity expression and anything.
template<typename E>
decltype(auto) operator/(OpNegIdentity const, E const& b)
{
	return expr::inverse(-b);
}

//! Division between two identities.
inline auto operator/(OpIdentity const, OpIdentity const)
{
	return OpIdentity{};
}

//! Division between two negative identities.
inline auto operator/(OpNegIdentity const, OpNegIdentity const)
{
	return OpIdentity{};
}

//! Division between the identity and the negative identity.
inline auto operator/(OpIdentity const, OpNegIdentity const)
{
	return OpNegIdentity{};
}

//! Division between the negative identity and the identity.
inline auto operator/(OpNegIdentity const, OpIdentity const)
{
	return OpNegIdentity{};
}

//! Division between anything and the 0 identity.
template<typename E>
auto operator/(E&&, OpVoid const) = delete;

//! Division between 0 and 0.
inline auto operator/(OpVoid const, OpVoid const) = delete;


//! Division between 0 identity and anything.
template<typename E>
decltype(auto) operator/(OpVoid const, E const&)
{
	return OpVoid{};
}

//! Division between 0 identity and identity expression.
inline auto operator/(OpVoid const, OpIdentity const)
{
	return OpVoid{};
}

//! Division between 0 identity and negative identity expression.
inline auto operator/(OpVoid const, OpNegIdentity const)
{
	return OpVoid{};
}

//! Division between a value constant and anything.
template<typename T, typename E>
decltype(auto) operator/(OpLiteral<T> const& a, E const& b)
{
	return expr::make_div(a, b);
}

//! Division between two value constants.
template<typename T, typename S>
auto operator/(OpLiteral<T> const& a, OpLiteral<S> const& b)
{
	return expr::make_literal(a.value / b.value);
}

//! Division between anything and a value constant.
template<typename T, typename E>
decltype(auto) operator/(E const& a, OpLiteral<T> const& b)
{
	return expr::inverse(b) * a;
}






/*
 *
 * Addition between basic values, particularly identities.
 *
 ******************************************************************************/


 //! Addition between integer and complex types.
inline auto operator+(int const a, complex_t const& b)
{
	return static_cast<double>(a) + b;
}

//! Addition between complex and integer types.
inline auto operator+(complex_t const& a, int const b)
{
	return a + static_cast<double>(b);
}



//! Addition between anything and the 0 identity.
template<typename E>
decltype(auto) operator+(E const& a, OpVoid const)
{
	return a;
}

//! Addition between the 0 identity and anything.
template<typename E>
decltype(auto) operator+(OpVoid const, E const& b)
{
	return b;
}

//! Addition between two 0 identities.
inline auto operator+(OpVoid const, OpVoid const)
{
	return OpVoid{};
}

//! Addition between two identities.
inline auto operator+(OpIdentity const, OpIdentity const)
{
	return expr::make_literal(2 * OpIdentity{}.eval());
}

//! Addition between two negative identities.
inline auto operator+(OpNegIdentity const, OpNegIdentity const)
{
	return expr::make_literal(2 * OpNegIdentity{}.eval());
}

//! Addition between the identity expression and negative identity.
inline OpVoid operator+(OpIdentity const, OpNegIdentity const)
{
	return OpVoid{};
}

//! Addition between the negative identity expression and identity.
inline OpVoid operator+(OpNegIdentity const, OpIdentity const)
{
	return OpVoid{};
}




//! Addition between a primitive value and the identity expression.
inline auto operator+(scalar_t const a, OpIdentity const)
{
	return a + OpIdentity{}.eval();
}

//! Addition between the identity expression and a primitive value.
inline auto operator+(OpIdentity const, scalar_t const b)
{
	return OpIdentity{}.eval() + b;
}

//! Addition between a primitive value and the identity expression.
inline auto operator+(int const a, OpIdentity const)
{
	return a + OpIdentity{}.eval();
}

//! Addition between the identity expression and a primitive value.
inline auto operator+(OpIdentity const, int const b)
{
	return OpIdentity{}.eval() + b;
}

//! Addition between a primitive complex type and the identity expression.
inline auto operator+(complex_t const& a, OpIdentity const)
{
	return a + static_cast<complex_t>(OpIdentity{}.eval());
}

//! Addition between the identity expression and a primitive complex type.
inline auto operator+(OpIdentity const, complex_t const b)
{
	return static_cast<complex_t>(OpIdentity{}.eval()) + b;
}




//! Addition between a primitive value and the negative identity expression.
inline auto operator+(scalar_t const a, OpNegIdentity const)
{
	return a + OpNegIdentity{}.eval();
}

//! Addition between the negative identity expression and a primitive value.
inline auto operator+(OpNegIdentity const, scalar_t const b)
{
	return OpNegIdentity{}.eval() + b;
}

//! Addition between a primitive value and the negative identity expression.
inline auto operator+(int const a, OpNegIdentity const)
{
	return a + OpNegIdentity{}.eval();
}

//! Addition between the negative identity expression and a primitive value.
inline auto operator+(OpNegIdentity const, int const b)
{
	return OpNegIdentity{}.eval() + b;
}

//! Addition between a primitive complex type and the negative identity.
inline auto operator+(complex_t const& a, OpNegIdentity const)
{
	return a + static_cast<complex_t>(OpNegIdentity{}.eval());
}

//! Addition between the negative identity and a primitive complex type.
inline auto operator+(OpNegIdentity const, complex_t const b)
{
	return static_cast<complex_t>(OpNegIdentity{}.eval()) + b;
}



/*
 *
 * Subtraction between basic values, particularly identities.
 *
 ******************************************************************************/

 //! Subtraction between anything and the 0 identity.
inline auto operator-(int const a, complex_t const& b)
{
	return static_cast<double>(a) - b;
}

//! Subtraction between complex and integer types.
inline auto operator-(complex_t const& a, int const b)
{
	return a - static_cast<double>(b);
}




//! Subtraction between anything and the 0 identity.
template<typename E>
decltype(auto) operator-(E const& a, OpVoid const)
{
	return a;
}

//! Subtraction between the 0 identity and anything.
template<typename E>
decltype(auto) operator-(OpVoid const, E const& b)
{
	return -b;
}

//! Subtraction between two 0 identities.
inline auto operator-(OpVoid const, OpVoid const)
{
	return OpVoid{};
}

//! Subtraction between two identities.
inline auto operator-(OpIdentity const, OpIdentity const)
{
	return OpVoid{};
}

//! Subtraction between two negative identities.
inline auto operator-(OpNegIdentity const, OpNegIdentity const)
{
	return OpVoid{};
}





//! Subtraction between a primitive value and the identity expression.
inline auto operator-(scalar_t const a, OpIdentity const)
{
	return a - OpIdentity{}.eval();
}

//! Subtraction between the identity expression and a primitive value.
inline auto operator-(OpIdentity const, scalar_t const b)
{
	return OpIdentity{}.eval() - b;
}

//! Subtraction between a primitive value and the identity expression.
inline auto operator-(int const a, OpIdentity const)
{
	return a - OpIdentity{}.eval();
}

//! Subtraction between the identity expression and a primitive value.
inline auto operator-(OpIdentity const, int const b)
{
	return OpIdentity{}.eval() - b;
}

//! Subtraction between a primitive complex type and the identity expression.
inline auto operator-(complex_t const a, OpIdentity const)
{
	return a - static_cast<complex_t>(OpIdentity{}.eval());
}

//! Subtraction between the identity expression and a primitive complex type.
inline auto operator-(OpIdentity const, complex_t const b)
{
	return static_cast<complex_t>(OpIdentity{}.eval()) - b;
}





//! Subtraction between a primitive value and the negative identity expression.
inline auto operator-(scalar_t const a, OpNegIdentity const)
{
	return a - OpNegIdentity{}.eval();
}

//! Subtraction between the negative identity expression and a primitive value.
inline auto operator-(OpNegIdentity const, scalar_t const b)
{
	return OpNegIdentity{}.eval() - b;
}

//! Subtraction between a primitive value and the negative identity expression.
inline auto operator-(int const a, OpNegIdentity const)
{
	return a - OpNegIdentity{}.eval();
}

//! Subtraction between the negative identity expression and a primitive value.
inline auto operator-(OpNegIdentity const, int const b)
{
	return OpNegIdentity{}.eval() - b;
}

//! Subtraction between a primitive complex type and the negative identity.
inline auto operator-(complex_t const& a, OpNegIdentity const)
{
	return a - static_cast<complex_t>(OpNegIdentity{}.eval());
}

//! Subtraction between the negative identity and a primitive complex type.
inline auto operator-(OpNegIdentity const, complex_t const& b)
{
	return static_cast<complex_t>(OpNegIdentity{}.eval()) - b;
}



template<typename E>
inline auto operator-(OpExpression<E> const& e, OpNegIdentity const)
{
	return e + OpIdentity{};
}

template<typename E>
inline auto operator+(OpExpression<E> const& e, OpNegIdentity const)
{
	return e - OpIdentity{};
}




// ******************************************************************************************





/*
 *
 * Collection of like terms.
 *
 ******************************************************************************/

namespace expr
{
	//! Returns true if the base type matches between two types.
	/*!
	 * The base type of a type is evaluated by expr::base_data_type, and 
	 * this will return true if the base type of two types is the same.
	 */
	template<typename G1, typename G2>
	constexpr bool is_same_base = std::is_same<
		typename expr::base_data_type<G1>::type,
		typename expr::base_data_type<G2>::type>::value;
}

//! Addition of two variables with data that can be combined.
template<typename A, typename B, typename G, typename std::enable_if<expr::is_combinable<G>, int>::type = 0>
auto operator+(OpLVariable<A, G> const& a, OpLVariable<B, G> const& b)
{
	return OpLVariable(a.value + b.value, a.data);
}

//! Subtraction of two variables with data that can be combined.
template<typename A, typename B, typename G, typename std::enable_if<expr::is_combinable<G>, int>::type = 0>
auto operator-(OpLVariable<A, G> const& a, OpLVariable<B, G> const& b)
{
	return OpLVariable(a.value - b.value, a.data);
}

//! Addition of two multi variables with data that can be combined.
template<typename A, typename B, typename... G1s, typename... G2s,
	typename std::enable_if<expr::nl_can_combine<G1s...>::template with<G2s...>::value, int>::type = 0>
	auto operator+(OpNLVariable<A, G1s...> const& a, OpNLVariable<B, G2s...> const& b)
{
	return OpNLVariable(a.value + b.value, a.datas);
}

//! Subtraction of two multi variables with data that can be combined.
template<typename A, typename B, typename... G1s, typename... G2s,
	typename std::enable_if<expr::nl_can_combine<G1s...>::template with<G2s...>::value, int>::type = 0>
	auto operator-(OpNLVariable<A, G1s...> const& a, OpNLVariable<B, G2s...> const& b)
{
	return OpNLVariable(a.value - b.value, a.datas);
}



//! Subtraction of two variables with data that cancel.
template<typename G, typename std::enable_if<expr::is_combinable<G>, int>::type = 0>
auto operator-(OpLVariable<OpIdentity, G> const&, OpLVariable<OpIdentity, G> const&)
{
	return OpVoid{};
}

//! Subtraction of two variables with data that cancel.
template<typename G, typename std::enable_if<expr::is_combinable<G>, int>::type = 0>
auto operator-(OpLVariable<OpNegIdentity, G> const&, OpLVariable<OpNegIdentity, G> const&)
{
	return OpVoid{};
}

//! Addition of two variables with data that cancel.
template<typename G, typename std::enable_if<expr::is_combinable<G>, int>::type = 0>
auto operator+(OpLVariable<OpNegIdentity, G> const&, OpLVariable<OpIdentity, G> const&)
{
	return OpVoid{};
}

//! Addition of two variables with data that cancel.
template<typename G, typename std::enable_if<expr::is_combinable<G>, int>::type = 0>
auto operator+(OpLVariable<OpIdentity, G> const&, OpLVariable<OpNegIdentity, G> const&)
{
	return OpVoid{};
}

//! Subtraction of two multi variables with data that cancel.
template<typename... G1s, typename... G2s,
	typename std::enable_if<expr::nl_can_combine<G1s...>::template with<G2s...>::value, int>::type = 0>
	auto operator-(OpNLVariable<OpIdentity, G1s...> const&, OpNLVariable<OpIdentity, G2s...> const&)
{
	return OpVoid{};
}

//! Subtraction of two multi variables with data that cancel.
template<typename... G1s, typename... G2s,
	typename std::enable_if<expr::nl_can_combine<G1s...>::template with<G2s...>::value, int>::type = 0>
	auto operator-(OpNLVariable<OpNegIdentity, G1s...> const&, OpNLVariable<OpNegIdentity, G2s...> const&)
{
	return OpVoid{};
}

//! Addition of two multi variables with data that cancel.
template<typename... G1s, typename... G2s,
	typename std::enable_if<expr::nl_can_combine<G1s...>::template with<G2s...>::value, int>::type = 0>
	auto operator+(OpNLVariable<OpNegIdentity, G1s...> const&, OpNLVariable<OpIdentity, G2s...> const&)
{
	return OpVoid{};
}

//! Addition of two multi variables with data that cancel.
template<typename... G1s, typename... G2s,
	typename std::enable_if<expr::nl_can_combine<G1s...>::template with<G2s...>::value, int>::type = 0>
	auto operator+(OpNLVariable<OpIdentity, G1s...> const&, OpNLVariable<OpNegIdentity, G2s...> const&)
{
	return OpVoid{};
}




/*
 *
 * Overloads to simplify expressions for division.
 *
 ******************************************************************************/


 //! Applies an identity to the division between two variables.
template<typename T1, typename T2, typename G1, typename G2,
	typename std::enable_if<expr::d_idy<G1, G2>::enable, int>::type = 0>
	auto operator/(OpLVariable<T1, G1> const& a, OpLVariable<T2, G2> const& b)
{
	return expr::d_idy<G1, G2>::template apply(a, b);
}

//! Division between two variables that can be reduced.
/*!
 * A fundamental rule for division, needs to be implemented for each new
 * expression term individually, such that the same type will reduce to a
 * literal formed by the division of their leading coefficients combinations
 * of terms (addition and subtraction) do not satisfy this rule.
 */
template<typename T1, typename T2, typename G,
	typename std::enable_if<(expr::is_combinable<G> && !expr::d_idy<G, G>::enable), int>::type = 0>
	auto operator/(OpLVariable<T1, G> const& a, OpLVariable<T2, G> const& b)
{
	return expr::make_literal(a.value * (OpIdentity{} / b.value));
}

//! Perform a usual division when two variables have no identity.
template<typename T1, typename T2, typename G1, typename G2,
	typename std::enable_if<(!expr::is_same_base<G1, G2> && !expr::d_idy<G1, G2>::enable), int>::type = 0>
	auto operator/(OpLVariable<T1, G1> const& a, OpLVariable<T2, G2> const& b)
{
	return expr::make_div(a, b);
}





/*
 *
 * Distributivity rules.
 *
 ******************************************************************************/

 //! Distributing subtraction expression into a subtraction expression.
template<typename A1, typename A2, typename B1, typename B2>
auto operator*(OpBinarySub<A1, A2> const& a, OpBinarySub<B1, B2> const& b)
{
	auto v1 = (*static_cast<const A1*>(&a.a)) * (*static_cast<const B1*>(&b.a));
	auto v2 = (*static_cast<const A2*>(&a.b)) * (*static_cast<const B1*>(&b.a));
	auto v3 = (*static_cast<const A1*>(&a.a)) * (*static_cast<const B2*>(&b.b));
	auto v4 = (*static_cast<const A2*>(&a.b)) * (*static_cast<const B2*>(&b.b));
	return v1 - v2 - v3 + v4;
}

//! Distributing addition expression into a subtraction expression.
template<typename A1, typename A2, typename B1, typename B2>
auto operator*(OpBinaryAdd<A1, A2> const& a, OpBinarySub<B1, B2> const& b)
{
	auto v1 = (*static_cast<const A1*>(&a.a)) * (*static_cast<const B1*>(&b.a));
	auto v2 = (*static_cast<const A2*>(&a.b)) * (*static_cast<const B1*>(&b.a));
	auto v3 = (*static_cast<const A1*>(&a.a)) * (*static_cast<const B2*>(&b.b));
	auto v4 = (*static_cast<const A2*>(&a.b)) * (*static_cast<const B2*>(&b.b));
	return v1 + v2 - v3 - v4;
}

//! Distributing subtraction expression into a addition expression.
template<typename A1, typename A2, typename B1, typename B2>
auto operator*(OpBinarySub<A1, A2> const& a, OpBinaryAdd<B1, B2> const& b)
{
	auto v1 = (*static_cast<const A1*>(&a.a)) * (*static_cast<const B1*>(&b.a));
	auto v2 = (*static_cast<const A2*>(&a.b)) * (*static_cast<const B1*>(&b.a));
	auto v3 = (*static_cast<const A1*>(&a.a)) * (*static_cast<const B2*>(&b.b));
	auto v4 = (*static_cast<const A2*>(&a.b)) * (*static_cast<const B2*>(&b.b));
	return v1 - v2 + v3 - v4;
}

//! Distributing addition expression into a addition expression.
template<typename A1, typename A2, typename B1, typename B2>
auto operator*(OpBinaryAdd<A1, A2> const& a, OpBinaryAdd<B1, B2> const& b)
{
	auto v1 = (*static_cast<const A1*>(&a.a)) * (*static_cast<const B1*>(&b.a));
	auto v2 = (*static_cast<const A2*>(&a.b)) * (*static_cast<const B1*>(&b.a));
	auto v3 = (*static_cast<const A1*>(&a.a)) * (*static_cast<const B2*>(&b.b));
	auto v4 = (*static_cast<const A2*>(&a.b)) * (*static_cast<const B2*>(&b.b));
	return v1 + v2 + v3 + v4;
}

//! Distributing an LHS expression between operands in addition expression.
template<typename A1, typename B1, typename B2>
auto operator*(OpExpression<A1> const& a, OpBinaryAdd<B1, B2> const& b)
{
	auto v1 = (*static_cast<const A1*>(&a)) * (*static_cast<const B1*>(&b.a));
	auto v2 = (*static_cast<const A1*>(&a)) * (*static_cast<const B2*>(&b.b));
	return v1 + v2;
}

//! Distributing an RHS expression between operands in addition expression.
template<typename B2, typename A1, typename A2>
auto operator*(OpBinaryAdd<A1, A2> const& a, OpExpression<B2> const& b)
{
	auto v1 = (*static_cast<const A1*>(&a.a)) * (*static_cast<const B2*>(&b));
	auto v2 = (*static_cast<const A2*>(&a.b)) * (*static_cast<const B2*>(&b));
	return v1 + v2;
}

//! Distributing an LHS expression between operands in subtraction expression.
template<typename A1, typename B1, typename B2>
auto operator*(OpExpression<A1> const& a, OpBinarySub<B1, B2> const& b)
{
	auto v1 = (*static_cast<const A1*>(&a)) * (*static_cast<const B1*>(&b.a));
	auto v2 = (*static_cast<const A1*>(&a)) * (*static_cast<const B2*>(&b.b));
	return v1 - v2;
}

//! Distributing an RHS expression between operands in subtraction expression.
template<typename B2, typename A1, typename A2>
auto operator*(OpBinarySub<A1, A2> const& a, OpExpression<B2> const& b)
{
	auto v1 = (*static_cast<const A1*>(&a.a)) * (*static_cast<const B2*>(&b));
	auto v2 = (*static_cast<const A2*>(&a.b)) * (*static_cast<const B2*>(&b));
	return v1 - v2;
}


//! Distributing a literal into the first operand in multiplication expression.
template<typename T, typename B1, typename B2>
auto operator*(OpLiteral<T> const& a, OpBinaryMul<B1, B2> const& b)
{
	return (a * b.a) * b.b;
}

//! Distributing a literal into the first operand in division expression.
template<typename T, typename B1, typename B2>
auto operator*(OpLiteral<T> const& a, OpBinaryDiv<B1, B2> const& b)
{
	return (a * b.a) / b.b;
}




namespace symphas::internal
{
	using expr::has_nmi_coeff;

	/*
	 *
	 * Overloads to simplify expressions for addition and subtraction.
	 *
	 ******************************************************************************/

	 //! End of the recursion for addition.
	 /*!
	  * Terminate the addition recursion, which groups like terms and return the
	  * final expression. Completes the addition by adding the *second* term of the
	  * addition with the second parameter.
	  */
	template<typename E1, typename E2, typename std::enable_if_t<!has_nmi_coeff<E1>, int> = 0>
	auto terminate_add_2(OpExpression<E1> const& a, OpExpression<E2> const& b)
	{
		return expr::make_add(*static_cast<E1 const*>(&a), *static_cast<E2 const*>(&b));
	}

	template<typename E>
	auto terminate_add_2(OpExpression<E> const& a, OpVoid const)
	{
		return *static_cast<E const*>(&a);
	}

	template<typename E>
	auto terminate_add_2(OpVoid const, OpExpression<E> const& b)
	{
		return *static_cast<E const*>(&b);
	}

	//! End of the recursion for subtraction.
	/*!
	 * See terminate_add_2(OpExpression const&, OpExpression const& b).
	 */
	template<typename E1, typename E2, typename std::enable_if_t<!has_nmi_coeff<E1>, int> = 0>
	auto terminate_sub_2(OpExpression<E1> const& a, OpExpression<E2> const& b)
	{
		return expr::make_sub(*static_cast<E1 const*>(&a), *static_cast<E2 const*>(&b));
	}

	template<typename E>
	auto terminate_sub_2(OpExpression<E> const& a, OpVoid const)
	{
		return *static_cast<E const*>(&a);
	}

	template<typename E>
	auto terminate_sub_2(OpVoid const, OpExpression<E> const& b)
	{
		return -*static_cast<E const*>(&b);
	}


	//! Overloaded end of the recursion for subtraction, removing double minus.
	/*
	 * See terminate_add_2(OpExpression const&, OpExpression const& b).
	 */
	template<typename E, typename E1, typename E2, typename std::enable_if_t<has_nmi_coeff<E1>, int> = 0>
	auto terminate_sub_2(OpExpression<E> const& a, OpBinaryAdd<E1, E2> const& b)
	{
		return expr::make_add(*static_cast<E const*>(&a), -b.a - b.b);
	}
	
	template<typename E, typename E1, typename E2, typename std::enable_if_t<has_nmi_coeff<E1>, int> = 0>
	auto terminate_sub_2(OpExpression<E> const& a, OpBinarySub<E1, E2> const& b)
	{
		return expr::make_add(*static_cast<E const*>(&a), -b.a + b.b);
	}

	template<typename E1, typename E2, typename std::enable_if_t<has_nmi_coeff<E1>, int> = 0>
	auto terminate_sub_2(OpExpression<E1> const& a, OpExpression<E2> const& b)
	{
		return expr::make_add(*static_cast<E1 const*>(&a), -*static_cast<E2 const*>(&b));
	}


	//! Overloaded end of the recursion for addition, removing leading minus.
	/*
	 * See terminate_add_2(OpExpression const&, OpExpression const& b).
	 */
	template<typename E, typename E1, typename E2, typename std::enable_if_t<has_nmi_coeff<E1>, int> = 0>
	auto terminate_add_2(OpExpression<E> const& a, OpBinaryAdd<E1, E2> const& b)
	{
		return expr::make_sub(*static_cast<E const*>(&a), -b.a - b.b);
	}

	template<typename E, typename E1, typename E2, typename std::enable_if_t<has_nmi_coeff<E1>, int> = 0>
	auto terminate_add_2(OpExpression<E> const& a, OpBinarySub<E1, E2> const& b)
	{
		return expr::make_sub(*static_cast<E const*>(&a), -b.a + b.b);
	}

	template<typename E1, typename E2, typename std::enable_if_t<has_nmi_coeff<E1>, int> = 0>
	auto terminate_add_2(OpExpression<E1> const& a, OpExpression<E2> const& b)
	{
		return expr::make_sub(*static_cast<E1 const*>(&a), -*static_cast<E2 const*>(&b));
	}





	//! Last step before end of recursion for addition.
	/*!
	 * Used in reorganizing the addition between two expressions before ending the
	 * recursion. See terminate_add_2(OpExpression const&, OpExpression const& b).
	 */
	template<typename E1, typename E2>
	auto terminate_add(OpExpression<E1> const& a, OpExpression<E2> const& b)
	{
		return terminate_add_2(*static_cast<E1 const*>(&a), *static_cast<E2 const*>(&b));
	}

	template<typename A1, typename A2, typename E>
	auto terminate_add(OpBinaryAdd<A1, A2> const& a, OpExpression<E> const& b)
	{
		return terminate_add_2(a.a, a.b + *static_cast<E const*>(&b));
	}

	template<typename A1, typename A2, typename E>
	auto terminate_add(OpBinarySub<A1, A2> const& a, OpExpression<E> const& b)
	{
		return terminate_sub_2(a.a, a.b - *static_cast<E const*>(&b));
	}


	//! Last step before end of recursion for subtraction.
	/*!
	 * Used in reorganizing the addition between two expressions before ending the
	 * recursion. See terminate_sub_2(OpExpression const&, OpExpression const& b).
	 */
	template<typename E1, typename E2>
	auto terminate_sub(OpExpression<E1> const& a, OpExpression<E2> const& b)
	{
		return terminate_sub_2(*static_cast<E1 const*>(&a), *static_cast<E2 const*>(&b));
	}

	template<typename A1, typename A2, typename E>
	auto terminate_sub(OpBinaryAdd<A1, A2> const& a, OpExpression<E> const& b)
	{
		return terminate_add_2(a.a, a.b - *static_cast<E const*>(&b));
	}

	template<typename A1, typename A2, typename E>
	auto terminate_sub(OpBinarySub<A1, A2> const& a, OpExpression<E> const& b)
	{
		return terminate_sub_2(a.a, a.b + *static_cast<E const*>(&b));
	}




	//! Terminates the multiplication recursion of expressions.
	/*!
	 * These multiplication rules take into account multiplication between variables
	 * and existing multiplication objects, so that variables which don't multiply
	 * with each other can be multiplied into a new variable which allows putting
	 * them into a multi variable.
	 *
	 * Consequently this assumes all multiplication is associative.
	 *
	 * Another assumption is that multiplication between objects that don't
	 * multiply with each other is commutative. i.e.:
	 * consider objects \f$A\f$, \f$B\f$ and \f$C\f$. \f$A\f$ and \f$B\f$ can be
	 * multiplied with each other to produce a multi variable variable, but \f$C\f$
	 * can only be multiplied with itself to produce a multi variable. Then
	 * \f$A \cdot B = AB\f$, but \f$A \cdot C = A \cdot C\f$ and
	 * \f$B \cdot C = B \cdot C\f$ however, the last two expressions are assumed to
	 * commute, e.g. \f$A \cdot C = C \cdot A\f$.
	 */
	template<typename E1, typename E2>
	auto terminate_mul(OpExpression<E1> const& a, OpExpression<E2> const& b)
	{
		return expr::reevaluate(expr::make_mul(*static_cast<E1 const*>(&a), *static_cast<E2 const*>(&b)));
	}

	template<typename A1, typename A2, typename E>
	auto terminate_mul(OpBinaryMul<A1, A2> const& a, OpExpression<E> const& b)
	{
		return expr::reevaluate(expr::make_mul(a.a * (*static_cast<E const*>(&b)), a.b));
	}

}


/*
 *
 *
 *
 * Recursive behaviour of addition; adding many terms together requires
 * simplifying if possible and collection of like terms if possible. This is
 * performed by adding together newly added terms with all other ones in the
 * sum, recursively having all terms be added or subtracted with each other.
 *
 ******************************************************************************/

template<typename A1, typename A2, typename E>
auto operator+(OpExpression<E> const& a, OpBinaryAdd<A1, A2> const& b)
{
	return symphas::internal::terminate_add(b.a + *static_cast<E const*>(&a), b.b);
}

template<typename A1, typename A2, typename E>
auto operator+(OpBinaryAdd<A1, A2> const& a, OpExpression<E> const& b)
{
	return symphas::internal::terminate_add(a.a + *static_cast<E const*>(&b), a.b);
}

template<typename A1, typename A2, typename E>
auto operator+(OpExpression<E> const& a, OpBinarySub<A1, A2> const& b)
{
	return symphas::internal::terminate_sub(b.a + *static_cast<E const*>(&a), b.b);
}

template<typename A1, typename A2, typename E>
auto operator+(OpBinarySub<A1, A2> const& a, OpExpression<E> const& b)
{
	return symphas::internal::terminate_sub(a.a + *static_cast<E const*>(&b), a.b);
}

template<typename A1, typename A2, typename E>
auto operator-(OpExpression<E> const& a, OpBinaryAdd<A1, A2> const& b)
{
	return symphas::internal::terminate_sub(-b.a + *static_cast<E const*>(&a), b.b);
}

template<typename A1, typename A2, typename E>
auto operator-(OpBinaryAdd<A1, A2> const& a, OpExpression<E> const& b)
{
	return symphas::internal::terminate_add(a.a - *static_cast<E const*>(&b), a.b);
}

template<typename A1, typename A2, typename E>
auto operator-(OpExpression<E> const& a, OpBinarySub<A1, A2> const& b)
{
	return symphas::internal::terminate_add(-b.a + *static_cast<E const*>(&a), b.b);
}

template<typename A1, typename A2, typename E>
auto operator-(OpBinarySub<A1, A2> const& a, OpExpression<E> const& b)
{
	return symphas::internal::terminate_sub(a.a - *static_cast<E const*>(&b), a.b);
}


template<typename A1, typename A2, typename B1, typename B2>
auto operator+(OpBinaryAdd<A1, A2> const& a, OpBinaryAdd<B1, B2> const& b)
{
	return a.a + (a.b + b);
}

template<typename A1, typename A2, typename B1, typename B2>
auto operator+(OpBinaryAdd<A1, A2> const& a, OpBinarySub<B1, B2> const& b)
{
	return a.a + (a.b + b);
}

template<typename A1, typename A2, typename B1, typename B2>
auto operator+(OpBinarySub<A1, A2> const& a, OpBinaryAdd<B1, B2> const& b)
{
	return a.a - (a.b - b);
}

template<typename A1, typename A2, typename B1, typename B2>
auto operator+(OpBinarySub<A1, A2> const& a, OpBinarySub<B1, B2> const& b)
{
	return a.a - (a.b - b);
}

template<typename A1, typename A2, typename B1, typename B2>
auto operator-(OpBinaryAdd<A1, A2> const& a, OpBinaryAdd<B1, B2> const& b)
{
	return a.a + (a.b - b);
}

template<typename A1, typename A2, typename B1, typename B2>
auto operator-(OpBinaryAdd<A1, A2> const& a, OpBinarySub<B1, B2> const& b)
{
	return a.a + (a.b - b);
}

template<typename A1, typename A2, typename B1, typename B2>
auto operator-(OpBinarySub<A1, A2> const& a, OpBinaryAdd<B1, B2> const& b)
{
	return a.a - (a.b + b);
}

template<typename A1, typename A2, typename B1, typename B2>
auto operator-(OpBinarySub<A1, A2> const& a, OpBinarySub<B1, B2> const& b)
{
	return a.a - (a.b + b);
}


/*
 *
 * Overloads to simplify expressions for multiplication.
 *
 ******************************************************************************/


template<typename A1, typename A2, typename T, typename G>
auto operator*(OpBinaryMul<A1, A2> const& a, OpLVariable<T, G> const& b)
{
	return a.a * (a.b * b);
}

template<typename T, typename G, typename B1, typename B2>
auto operator*(OpLVariable<T, G> const& a, OpBinaryMul<B1, B2> const& b)
{
	return symphas::internal::terminate_mul(a * b.a, b.b);
}

template<typename A1, typename A2, typename T, typename... Gs>
auto operator*(OpBinaryMul<A1, A2> const& a, OpNLVariable<T, Gs...> const& b)
{
	return a.a * (a.b * b);
}

template<typename T, typename... Gs, typename B1, typename B2>
auto operator*(OpNLVariable<T, Gs...> const& a, OpBinaryMul<B1, B2> const& b)
{
	return symphas::internal::terminate_mul(a * b.a, b.b);
}

template<typename A1, typename A2, typename B1, typename B2>
auto operator*(OpBinaryMul<A1, A2> const& a, OpBinaryMul<B1, B2> const& b)
{
	return a.a * (a.b * b);
}






namespace symphas::internal
{

	//! Terminates the recursion of the division operation.
	/*!
	 * Endpoint of the recursion for applying the division rules, which multiplies
	 * two reduced terms, with a key emphasis on reduced terms. there are several
	 * rules for handling various reduced forms. There is also the capability to
	 * automatically reduce terms. Makes the assumption that a division is always
	 * the first term given to the function.
	 */
	template<typename A2, typename E>
	auto terminate_div(OpBinaryDiv<OpIdentity, A2> const& a, OpExpression<E> const& b)
	{
		return expr::make_div((*static_cast<const E*>(&b)), (*static_cast<const A2*>(&a.b)));
	}

	template<typename A2, typename E>
	auto terminate_div(OpBinaryDiv<OpNegIdentity, A2> const& a, OpExpression<E> const& b)
	{
		return expr::make_div(-(*static_cast<const E*>(&b)), (*static_cast<const A2*>(&a.b)));
	}

	template<typename T, typename A2, typename E>
	auto terminate_div(OpBinaryDiv<OpLiteral<T>, A2> const& a, OpExpression<E> const& b)
	{
		return expr::make_div(a.a * (*static_cast<const E*>(&b)), (*static_cast<const A2*>(&a.b)));
	}

	template<typename A1, typename A2, typename E>
	auto terminate_div(OpBinaryDiv<A1, A2> const& a, OpExpression<E> const& b)
	{
		return expr::make_mul(a, (*static_cast<const E*>(&b)));
	}

	template<typename A1, typename A2, typename B1, typename B2>
	auto terminate_div(OpBinaryDiv<A1, A2> const& a, OpBinaryDiv<B1, B2> const& b)
	{
		return expr::make_div(
			(*static_cast<const A1*>(&a.a)) * (*static_cast<const B1*>(&b.a)),
			(*static_cast<const A2*>(&a.b)) * (*static_cast<const B2*>(&b.b)));
	}


	template<typename A1, typename A2, typename T>
	auto terminate_div(OpBinaryDiv<A1, A2> const& a, OpLiteral<T> const& b)
	{
		return b * a;
	}

	template<typename A1, typename A2>
	auto terminate_div(OpBinaryDiv<A1, A2> const& a, OpIdentity const)
	{
		return a;
	}

	template<typename A1, typename A2>
	auto terminate_div(OpBinaryDiv<A1, A2> const& a, OpNegIdentity const)
	{
		return -a;
	}


	template<typename E, typename T>
	auto terminate_div(OpLiteral<T> const& a, OpExpression<E> const& b)
	{
		return (*static_cast<const E*>(&b)) * a;
	}

	template<typename E>
	auto terminate_div(OpIdentity const, OpExpression<E> const& b)
	{
		return (*static_cast<const E*>(&b));
	}

	template<typename E>
	auto terminate_div(OpNegIdentity const, OpExpression<E> const& b)
	{
		return -(*static_cast<const E*>(&b));
	}

	template<typename T, typename G, typename E>
	auto terminate_div(OpLVariable<T, G> const& a, OpExpression<E> const& b)
	{
		return expr::make_mul(a, (*static_cast<const E*>(&b)));
	}

	template<typename T, typename G, typename S>
	auto terminate_div(OpLVariable<T, G> const& a, OpLiteral<S> const& b)
	{
		return b * a;
	}

	template<typename T, typename G>
	auto terminate_div(OpLVariable<T, G> const& a, OpIdentity const)
	{
		return a;
	}

	template<typename T, typename G>
	auto terminate_div(OpLVariable<T, G> const& a, OpNegIdentity const)
	{
		return -a;
	}



	template<typename A1, typename A2, typename E>
	auto terminate_div(OpBinaryMul<A1, A2> const& a, OpExpression<E> const& b)
	{
		return expr::make_mul(a, (*static_cast<const E*>(&b)));
	}

	template<typename A1, typename A2, typename T>
	auto terminate_div(OpBinaryMul<A1, A2> const& a, OpLiteral<T> const& b)
	{
		return b * a;
	}

	template<typename A1, typename A2>
	auto terminate_div(OpBinaryMul<A1, A2> const& a, OpIdentity const)
	{
		return a;
	}

	template<typename A1, typename A2>
	auto terminate_div(OpBinaryMul<A1, A2> const& a, OpNegIdentity const)
	{
		return -a;
	}

	/* handles a special case; puts the 2nd parameter into the denominator of the
	 * first
	 */

	template<typename E1, typename E2>
	auto to_denom(OpExpression<E1> const& a, OpExpression<E2> const& b)
	{
		return expr::make_div(*static_cast<const E1*>(&a), (*static_cast<const E2*>(&b)));
	}

	template<typename A1, typename A2, typename E>
	auto to_denom(OpBinaryDiv<A1, A2> const& a, OpExpression<E> const& b)
	{
		return expr::make_div(*static_cast<const A1*>(&a.a), (*static_cast<const A2*>(&a.b) * (*static_cast<const E*>(&b))));
	}

	template<typename E>
	auto to_denom(OpNegIdentity const, OpExpression<E> const& b)
	{
		return inverse(-*static_cast<E const*>(&b));
	}

	template<typename E>
	auto to_denom(OpIdentity const, OpExpression<E> const& b)
	{
		return inverse(*static_cast<E const*>(&b));
	}

	template<typename T, typename E>
	auto to_denom(OpLiteral<T> const& a, OpExpression<E> const& b)
	{
		return a * inverse(*static_cast<E const*>(&b));
	}
}



/*
 * Recursively handle division between multiplications and divisions.
 */

template<typename A1, typename A2, typename E,
	typename std::enable_if<expr::factor_list_all<OpBinaryMul<A1, A2>, E>::value == 0, int>::type = 0>
	auto operator/(OpBinaryMul<A1, A2> const& a, OpExpression<E> const& b)
{
	return ((*static_cast<const A1*>(&a.a)) / (*static_cast<const E*>(&b)))
		* (*static_cast<const A2*>(&a.b));
}

template<typename B1, typename B2, typename E,
	typename std::enable_if<expr::factor_list_all<E, OpBinaryMul<B1, B2>>::value == 0, int>::type = 0>
	auto operator/(OpExpression<E> const& a, OpBinaryMul<B1, B2> const& b)
{
	return ((*static_cast<const E*>(&a)) / (*static_cast<const B1*>(&b.a)))
		* expr::inverse(*static_cast<const B2*>(&b.b));
}

template<typename A1, typename A2, typename B1, typename B2,
	typename std::enable_if<expr::factor_list_all<OpBinaryMul<A1, A2>, OpBinaryMul<B1, B2>>::value == 0, int>::type = 0>
	auto operator/(OpBinaryMul<A1, A2> const& a, OpBinaryMul<B1, B2> const& b)
{
	return ((*static_cast<const A1*>(&a.a)) / (*static_cast<const B1*>(&b.a)))
		* (*static_cast<const A2*>(&a.b) / *static_cast<const B2*>(&b.b));
}

template<typename A1, typename A2, typename E,
	typename std::enable_if<expr::factor_list_all<OpBinaryDiv<A1, A2>, E>::value == 0, int>::type = 0>
	auto operator/(OpBinaryDiv<A1, A2> const& a, OpExpression<E> const& b)
{
	return symphas::internal::to_denom(
		(*static_cast<const A1*>(&a.a)) / (*static_cast<const E*>(&b)),
		*static_cast<const A2*>(&a.b));
}

template<typename B1, typename B2, typename E,
	typename std::enable_if<expr::factor_list_all<E, OpBinaryDiv<B1, B2>>::value == 0, int>::type = 0>
	auto operator/(OpExpression<E> const& a, OpBinaryDiv<B1, B2> const& b)
{
	return (*static_cast<const E*>(&a))
		* ((*static_cast<const B2*>(&b.b)) / (*static_cast<const B1*>(&b.a)));
}

template<typename A1, typename A2, typename B1, typename B2,
	typename std::enable_if<expr::factor_list_all<OpBinaryDiv<A1, A2>, OpBinaryDiv<B1, B2>>::value == 0, int>::type = 0>
	auto operator/(OpBinaryDiv<A1, A2> const& a, OpBinaryDiv<B1, B2> const& b)
{
	return ((*static_cast<const A1*>(&a.a)) / (*static_cast<const B1*>(&b.a)))
		* ((*static_cast<const B2*>(&b.b)) / (*static_cast<const A2*>(&a.b)));
}



/*
 *
 *
 * Multiplication rules in the context of division, used for canceling like 
 * terms. All the divisions that are parameters to these functions are assumed 
 * to be in reduced form. 
 *
 ******************************************************************************/

template<typename A1, typename A2, typename E>
auto operator*(OpBinaryDiv<A1, A2> const& a, OpExpression<E> const& b)
{
	return symphas::internal::terminate_div(
		(*static_cast<const E*>(&b)) / (*static_cast<const A2*>(&a.b)), 
		*static_cast<const A1*>(&a.a));
}

template<typename E, typename B1, typename B2>
auto operator*(OpExpression<E> const& a, OpBinaryDiv<B1, B2> const& b)
{
	return b * (*static_cast<const E*>(&a));
}

template<typename A1, typename A2, typename B1, typename B2>
auto operator*(OpBinaryDiv<A1, A2> const& a, OpBinaryDiv<B1, B2> const& b)
{
	return symphas::internal::terminate_div(
		(*static_cast<const A1*>(&a.a)) / (*static_cast<const B2*>(&b.b)), 
		(*static_cast<const B1*>(&b.a)) / (*static_cast<const A2*>(&a.b)));
}



template<typename A1, typename A2, typename B1, typename B2>
auto operator*(OpBinaryDiv<A1, A2> const& a, OpBinaryAdd<B1, B2> const& b)
{
	return (*static_cast<const A1*>(&a.a)) 
		* ((*static_cast<const B1*>(&b.a)) / (*static_cast<const A2*>(&a.b)) 
		+ (*static_cast<const B2*>(&b.b)) / (*static_cast<const A2*>(&a.b)));
}

template<typename A1, typename A2, typename B1, typename B2>
auto operator*(OpBinaryDiv<A1, A2> const& a, OpBinarySub<B1, B2> const& b)
{
	return (*static_cast<const A1*>(&a.a))
		* ((*static_cast<const B1*>(&b.a)) / (*static_cast<const A2*>(&a.b)) 
		- (*static_cast<const B2*>(&b.b)) / (*static_cast<const A2*>(&a.b)));
}

template<typename A1, typename A2, typename B1, typename B2>
auto operator*(OpBinaryAdd<A1, A2> const& a, OpBinaryDiv<B1, B2> const& b)
{
	return b * a;
}

template<typename A1, typename A2, typename B1, typename B2>
auto operator*(OpBinarySub<A1, A2> const& a, OpBinaryDiv<B1, B2> const& b)
{
	return b * a;
}

/* Handling recursion at the point of an inverse.
 */

template<typename A2, typename B1, typename B2>
auto operator*(OpBinaryDiv<OpNegIdentity, A2> const& a, OpBinaryDiv<B1, B2> const& b)
{
	return b * a;
}

template<typename A2, typename B1, typename B2>
auto operator*(OpBinaryDiv<OpIdentity, A2> const& a, OpBinaryDiv<B1, B2> const& b)
{
	return b * a;
}

template<typename T, typename A2, typename B1, typename B2>
auto operator*(OpBinaryDiv<OpLiteral<T>, A2> const& a, OpBinaryDiv<B1, B2> const& b)
{
	return b * a;
}

template<typename A1, typename A2, typename B2>
auto operator*(OpBinaryDiv<A1, A2> const& a, OpBinaryDiv<OpNegIdentity, B2> const& b)
{
	return symphas::internal::terminate_div(
		(*static_cast<const A1*>(&a.a)) / (*static_cast<const B2*>(&b.b)), 
		expr::inverse(*static_cast<const A2*>(&a.b)));
}

template<typename A1, typename A2, typename B2>
auto operator*(OpBinaryDiv<A1, A2> const& a, OpBinaryDiv<OpIdentity, B2> const& b)
{
	return symphas::internal::terminate_div(
		(*static_cast<const A1*>(&a.a)) / (*static_cast<const B2*>(&b.b)), 
		expr::inverse(*static_cast<const A2*>(&a.b)));
}

template<typename A1, typename A2, typename T, typename B2>
auto operator*(OpBinaryDiv<A1, A2> const& a, OpBinaryDiv<OpLiteral<T>, B2> const& b)
{
	return b.a * symphas::internal::terminate_div(
		(*static_cast<const A1*>(&a.a)) / (*static_cast<const B2*>(&b.b)), 
		expr::inverse(*static_cast<const A2*>(&a.b)));
}




/* Division rules for some other objects but will certainly not be used in 
 * practice.
 */

template<typename V1, typename V2, typename Dd, typename E, typename Sp>
auto operator/(OpFuncDerivative<Dd, V1, E, Sp> const& a, OpFuncDerivative<Dd, V2, E, Sp> const& b)
{
	return expr::make_literal(a.value * (OpIdentity{} / b.value));
}



// ******************************************************************************************


template<typename E1, typename E2>
auto OpBinaryAdd<E1, E2>::operator-() const
{
	return -a - b;
}

template<typename E1, typename E2>
auto OpBinarySub<E1, E2>::operator-() const
{
	return -a + b;
}

template<typename E1, typename E2>
auto OpBinaryMul<E1, E2>::operator-() const
{
	return -a * b;
}

template<typename E1, typename E2>
auto OpBinaryDiv<E1, E2>::operator-() const
{
	return -a / b;
}



template<typename E1>
auto operator+(OpExpression<E1> const& a, double b)
{
	return *static_cast<const E1*>(&a) + expr::make_literal(b);
}

template<typename E2>
auto operator+(double a, OpExpression<E2> const& b)
{
	return expr::make_literal(a) + *static_cast<const E2*>(&b);
}


template<typename E1>
auto operator-(OpExpression<E1> const& a, double b)
{
	return *static_cast<const E1*>(&a) - expr::make_literal(b);
}

template<typename E2>
auto operator-(double a, OpExpression<E2> const& b)
{
	return expr::make_literal(a) - *static_cast<const E2*>(&b);
}


template<typename E1>
auto operator*(OpExpression<E1> const& a, double b)
{
	return *static_cast<const E1*>(&a) * expr::make_literal(b);
}

template<typename E2>
auto operator*(double a, OpExpression<E2> const& b)
{
	return expr::make_literal(a) * *static_cast<const E2*>(&b);
}

template<typename E1>
auto operator/(OpExpression<E1> const& a, double b)
{
	return *static_cast<const E1*>(&a) / expr::make_literal(b);
}

template<typename E2>
auto operator/(double a, OpExpression<E2> const& b)
{
	return expr::make_literal(a) / *static_cast<const E2*>(&b);
}



template<typename E1>
auto operator+(OpExpression<E1> const& a, int b)
{
	return *static_cast<const E1*>(&a) + expr::make_literal(b);
}

template<typename E2>
auto operator+(int a, OpExpression<E2> const& b)
{
	return expr::make_literal(a) + *static_cast<const E2*>(&b);
}


template<typename E1>
auto operator-(OpExpression<E1> const& a, int b)
{
	return *static_cast<const E1*>(&a) - expr::make_literal(b);
}

template<typename E2>
auto operator-(int a, OpExpression<E2> const& b)
{
	return expr::make_literal(a) - *static_cast<const E2*>(&b);
}


template<typename E1>
auto operator*(OpExpression<E1> const& a, int b)
{
	return *static_cast<const E1*>(&a) * expr::make_literal(b);
}

template<typename E2>
auto operator*(int a, OpExpression<E2> const& b)
{
	return expr::make_literal(a) * *static_cast<const E2*>(&b);
}

template<typename E1>
auto operator/(OpExpression<E1> const& a, int b)
{
	return *static_cast<const E1*>(&a) / expr::make_literal(b);
}

template<typename E2>
auto operator/(int a, OpExpression<E2> const& b)
{
	return expr::make_literal(a) / *static_cast<const E2*>(&b);
}





