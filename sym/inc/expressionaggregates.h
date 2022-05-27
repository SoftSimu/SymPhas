
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
 * PURPOSE: Implements the Variable, NamedData and binary operators.
 *
 * ***************************************************************************
 */

#pragma once

#include <tuple>
#include <functional>
#include <vector>

#include "expressionrules.h"
#include "expressionsprint.h"

//! \cond
#ifdef EXPR_EXPORTS
#define DLLEXPR DLLEXPORT
#else
#define DLLEXPR DLLIMPORT
#endif
//! \endcond


//! Associates an identifier with a piece of data.
/*!
 * Associates a template constant number with a given
 * object; primarily used for parsing functions. Used in particular when 
 * parsing or pruning equations. Acts as a data identifier.
 */
template<size_t Z, typename G>
struct Variable : G
{
	using G::G;

	Variable(G const& data) : G(data) {}
	Variable(G&& data) noexcept : G(std::forward<G>(data)) {}
};


//! Gives a name to data, to be used in expressions.
/*!
 * By inheriting from an object that acts like data for an OpLVariable, this
 * construct will add a name that can be interpreted by the get_op_name()
 * function. When used in conjunction with Variable, Variable will take this as
 * the template parameter, not the other way around.
 * 
 * \tparam G The object from which the NamedData inherits. Since the NamedData
 * gives a piece of data a name, the template refers to the data type, not
 * any expression type or Variable.
 */
template<typename G>
struct NamedData : G
{
#ifdef PRINTABLE_EQUATIONS

	//! Constructs the named instance using the given data.
	/*!
	 * Constructs the named instance using the given data, which is passed
	 * to the constructor of the underlying class.
	 * 
	 * \param data The data that is named.
	 * \param name The name of the data.
	 */
	NamedData(G data, std::string name) : G(data), name{ name } {}


	//! Constructs the named instance represented by the given expression.
	/*!
	 * Constructs the named instance with a name that is derived from the
	 * given instance. That is, it is assumed that the given data is
	 * equal to the given expression, and the name derived from the print
	 * output of that expression.
	 *
	 * \param data The data that is named.
	 * \param expr From where the name is derived.
	 */
	template<typename E>
	NamedData(G data, OpExpression<E> const& e) : G(data), name{ "" }
	{
		char* str = new char[e.print_length() + 1];
		e.print(str);
		name = std::string(str);
		delete[] str;
	}

	std::string name;	//!< Name given to the data.

#else

	NamedData(G data, ...) : G(data) {}

#endif


};

//! Gives a name to a data array, to be used in expressions.
/*!
 * Provides a name to a data array, and is a specialization of ::NamedData<G>.
 * The specialization is necessary because a pointer type cannot be inherited,
 * so the data object must mimic a pointer instead, and store the pointer.
 * To this end, deallocating the original pointer array would also make this
 * object invalid.
 *
 * \tparam G The object from which the NamedData inherits. Since the NamedData
 * gives a piece of data a name, the template refers to the data type, not
 * any expression type or Variable.
 */
template<typename G>
struct NamedData<G*>
{

#ifdef PRINTABLE_EQUATIONS

	//! Constructs the named instance using the given data.
	/*!
	 * Constructs the named instance using the given data, which is passed
	 * to the constructor of the underlying class.
	 *
	 * \param data The data that is named.
	 * \param name The name of the data.
	 */
	NamedData(G* data, std::string name) : data{ data }, name{ name } {}

	//! Constructs the named instance represented by the given expression.
	/*!
	 * Constructs the named instance with a name that is derived from the
	 * given instance. That is, it is assumed that the given data is
	 * equal to the given expression, and the name derived from the print
	 * output of that expression.
	 *
	 * \param data The data that is named.
	 * \param expr From where the name is derived.
	 */
	template<typename E>
	NamedData(G* data, OpExpression<E> const& e) : data{ data }, name{ "" }
	{
		char* str = new char[e.print_length() + 1];
		e.print(str);
		name = std::string(str);
		delete[] str;
	}

	std::string name;	//!< Name given to the data.

#else

	NamedData(G* data, ...) : data{ data } {}

#endif


	operator G* ()
	{
		return data;
	}

	operator G* const () const
	{
		return data;
	}

	G* data;			//!< Pointer to the data.
};


//! Gives a name to a double value, to be used in expressions.
/*!
 * Provides a name to a data array, and is a specialization of ::NamedData<G>.
 * The specialization is necessary because a POD cannot be inherited,
 * so the data object must mimic the ::Scalar_t instead instead, and store the pointer.
 * To this end, deallocating the original pointer array would also make this
 * object invalid.
 */
template<>
struct NamedData<scalar_t>
{
#ifdef PRINTABLE_EQUATIONS

	//! Constructs the named instance using the given data.
	/*!
	 * Constructs the named instance using the given data, which is passed
	 * to the constructor of the underlying class.
	 *
	 * \param data The data that is named.
	 * \param name The name of the data.
	 */
	NamedData(scalar_t data, std::string name) : data{ data }, name{ name } {}

	//! Constructs the named instance represented by the given expression.
	/*!
	 * Constructs the named instance with a name that is derived from the
	 * given instance. That is, it is assumed that the given data is
	 * equal to the given expression, and the name derived from the print
	 * output of that expression.
	 *
	 * \param data The data that is named.
	 * \param expr From where the name is derived.
	 */
	template<typename E>
	NamedData(scalar_t data, OpExpression<E> const& e) : data{ data }, name{ "" }
	{
		char* str = new char[e.print_length() + 1];
		e.print(str);
		name = std::string(str);
		delete[] str;
	}

	std::string name;	//!< Name given to the data.

#else

	NamedData(scalar_t data, ...) : data{ data } {}

#endif

	operator scalar_t()
	{
		return data;
	}

	operator const scalar_t() const
	{
		return data;
	}

	scalar_t data;		//!< Data value.
};

template<typename G>
NamedData(G&, std::string)->NamedData<symphas::ref<G>>;
template<typename G>
NamedData(G*, std::string)->NamedData<G*>;


namespace expr
{
	//! The given parameter is created as a Variable.
	/*!
	 * The given parameter is created and returned as a Variable.
	 * 
	 * \param g The object which is created as a Variable.
	 */
	template<size_t Z, typename G>
	auto as_variable(G&& g)
	{
		return Variable<Z, G>(std::forward<G>(g));
	}

	//! The given parameter is associated to a Variable.
	/*!
	 * The given parameter is passed by reference, so its contents are assumed
	 * to be associated outside of the scope of a Variable. In that case, the
	 * given parameter is associated by a reference wrapper to a new Variable
	 * instance which is returned as the result.
	 *
	 * \param g The object which is associated to a Variable.
	 */
	template<size_t Z, typename G>
	auto as_variable(G& g)
	{
		return Variable<Z, symphas::ref<G>>(g);
	}

#ifdef PRINTABLE_EQUATIONS

	DLLEXPR extern int NAME_PTR_POS;					//!< Current position in selecting name for arbitrary data pointers.
	DLLEXPR extern std::vector<const void*> NAME_PTRS;	//!< List of all data pointers with names associated with them.
	DLLEXPR extern std::vector<char*> MORE_NAMES;		//!< List of overflow names for data pointers.


	//! Gets the string name associated with the data.
	template<typename A>
	const char* get_op_name(symphas::ref<A> const& a);

	//! Gets the string name associated with the data.
	template<typename T>
	const char* get_op_name(T* ptr);


	//! Gets the string name associated with the data.
	template<>
	inline const char* get_op_name(char* a)
	{
		return a;
	}

	//! Gets the string name associated with the data.
	template<>
	inline const char* get_op_name(const char* a)
	{
		return a;
	}


	//! Gets the string name associated with the data.
	template<typename A>
	const char* get_op_name(A const& a)
	{
		return get_op_name(expr::SymbolID<A>::get(a));
	}

	//! Gets the string name associated with the data.
	template<typename A>
	const char* get_op_name(symphas::ref<A> const& a)
	{
		return get_op_name(expr::SymbolID<A>::get(a));
	}

	//! Gets the string name associated with the data.
	template<size_t N>
	const char* get_op_name(Variable<N> const& a)
	{
		static size_t NN = 0;
		static char** names;
		const char prefix[] = "var";
		if (N >= NN)
		{
			char** new_names = new char* [N + 1];
			for (iter_type i = 0; i < NN; ++i)
			{
				new_names[i] = new char[std::strlen(names[i]) + 1];
				std::strcpy(new_names[i], names[i]);
			}
			for (size_t i = NN; i <= N; ++i)
			{
				new_names[i] = new char[STR_ARR_LEN(prefix) + symphas::lib::num_digits(N)];
				sprintf(new_names[i], "%s%zd", prefix, N);
			}
			delete[] names;
			names = new_names;
			return names[N];
		}
		else
		{
			return names[N];
		}
	}

	//! Gets the string name associated with the data.
	template<typename T>
	const char* get_op_name(T* ptr)
	{
		if (!ptr)
		{
			return "?";
		}
		else
		{
			const void* ptr_cmp = static_cast<const void*>(ptr);
			constexpr size_t MAX_NAME_COUNT = sizeof(VARIABLE_NAMES) / sizeof(*VARIABLE_NAMES);

			for (iter_type i = 0; i < NAME_PTR_POS; ++i)
			{
				if (NAME_PTRS[i] == ptr_cmp)
				{
					if (i < MAX_NAME_COUNT)
					{
						return VARIABLE_NAMES[i];
					}
					else
					{
						return MORE_NAMES[i - MAX_NAME_COUNT];
					}
				}
			}
			NAME_PTRS.push_back(ptr_cmp);

			if (NAME_PTR_POS < MAX_NAME_COUNT)
			{
				return VARIABLE_NAMES[NAME_PTR_POS++];
			}
			else
			{
				char* name = new char[BUFFER_LENGTH_R4];
				sprintf(name, VARIABLE_NAME_EXTRA_FMT, NAME_PTR_POS);
				MORE_NAMES.push_back(name);
				return MORE_NAMES[NAME_PTR_POS++ - MAX_NAME_COUNT];
			}
		}
	}

#endif
}


// **************************************************************************************

namespace symphas::internal
{

	//! Helper in order to construct an OpLVariable.
	/*!
	 * Allows the type to be delegated to the correct function, so that
	 * the correct OpLVariable object is generated.
	 * In particular, it is important that if a reference is made, that
	 * an additional reference is not nested.
	 * 
	 * The usual template will assume that an rvalue is given, and 
	 * therefore the data can be moved into the OpLVariable.
	 * 
	 * \tparam A The type being given to create an OpLVariable.
	 */
	template<typename A>
	struct make_oplvariable
	{
		template<typename S>
		static auto get(S value, A data)
		{
			return OpLVariable(value, data);
		}

		static auto get(A data)
		{
			return OpLVariable(data);
		}

		template<size_t Z, typename S>
		static auto get(S value, A data)
		{
			return OpLVariable(value, Variable<Z, A>(data));
		}

		template<size_t Z>
		static auto get(A data)
		{
			return OpLVariable(Variable<Z, A>(data));
		}
	};


	//! Helper in order to construct an OpLVariable.
	/*!
	 * Specialization for an lvalue; this means that a reference should be 
	 * created.
	 */
	template<typename A>
	struct make_oplvariable<A&>
	{
		template<typename S>
		static auto get(S value, A& data)
		{
			return OpLVariable(value, symphas::ref<A>(data));
		}

		static auto get(A& data)
		{
			return OpLVariable(symphas::ref<A>(data));
		}

		template<size_t Z, typename S>
		static auto get(S value, A& data)
		{
			return OpLVariable(value, Variable<Z, symphas::ref<A>>(data));
		}

		template<size_t Z>
		static auto get(A& data)
		{
			return OpLVariable(Variable<Z, symphas::ref<A>>(data));
		}
	};

	//! Helper in order to construct an OpLVariable.
	/*!
	 * Specialization for a pointer; this means the pointer should be copied to
	 * the variable.
	 */
	template<typename A>
	struct make_oplvariable<A*>
	{
		template<typename S>
		static auto get(S value, A* data)
		{
			return OpLVariable(value, data);
		}

		static auto get(A* data)
		{
			return OpLVariable(data);
		}

		template<size_t Z, typename S>
		static auto get(S value, A* data)
		{
			return get(value, data);
		}

		template<size_t Z>
		static auto get(A* data)
		{
			return get(data);
		}
	};

	//! Helper in order to construct an OpLVariable.
	/*!
	 * The data is passed as a reference, rather than copied.
	 */
	template<size_t Z, typename A>
	struct make_oplvariable<Variable<Z, A>>
	{
		template<typename S>
		static auto get(S value, Variable<Z, A>& data)
		{
			return make_oplvariable<A>::template get<Z>(value, data);
		}

		static auto get(Variable<Z, A>& data)
		{
			return make_oplvariable<A>::template get<Z>(data);
		}

		template<size_t Y, typename S>
		static auto get(S value, Variable<Z, A>& data)
		{
			return make_oplvariable<A>::template get<Y>(value, data);
		}

		template<size_t Y>
		static auto get(Variable<Z, A>& data)
		{
			return make_oplvariable<A>::template get<Y>(data);
		}
	};

	//! Helper in order to construct an OpLVariable.
	/*!
	 * The data is passed as a reference, rather than copied.
	 */
	template<typename A>
	struct make_oplvariable<symphas::ref<A>>
	{
		template<typename S>
		static auto get(S value, symphas::ref<A> data)
		{
			return make_oplvariable<A>::template get(value, data.get());
		}

		static auto get(symphas::ref<A> data)
		{
			return make_oplvariable<A>::get(data.get());
		}

		template<size_t Z, typename S>
		static auto get(S value, symphas::ref<A> data)
		{
			return make_oplvariable<A>::template get<Z>(value, data.get());
		}

		template<size_t Z>
		static auto get(symphas::ref<A> data)
		{
			return make_oplvariable<A>::template get<Z>(data.get());
		}
	};


	//! Helper in order to construct an OpLVariable.
	/*!
	 * When these types are given by lvalue, then the original functions need to
	 * be applied, not the functions corresponding to the simple lvalue
	 * specialization.
	 */
	template<typename A>
	struct make_oplvariable<A*&> : make_oplvariable<A*> {};

	template<typename A>
	struct make_oplvariable<symphas::ref<A>&> : make_oplvariable<symphas::ref<A>> {};

	template<size_t Z, typename A>
	struct make_oplvariable<Variable<Z, A>&> : make_oplvariable<Variable<Z, A>> {};
	template<typename A>
	struct make_oplvariable<NamedData<A>&> : make_oplvariable<NamedData<A>> {};
}


namespace expr
{

	//! Delegate to create an OpLVariable with the correct template. 
	/*!
	 * Delegate to create an OpLVariable with the correct template. For
	 * more information on specific implementation, see symphas::internal::make_oplvariable.
	 *
	 * \param data The data from which to create an OpLVariable.
	 *
	 * \tparam A The type of the data to make an OpLVariable.
	 */
	template<typename A>
	auto make_op(A&& data)
	{
		return symphas::internal::make_oplvariable<A>::get(std::forward<A>(data));
	}

	//! Delegate to create an OpLVariable with the correct template. 
	/*!
	 * Delegate to create an OpLVariable with the correct template. For
	 * more information on specific implementation, see symphas::internal::make_oplvariable.
	 * This overload will use the given index to create an OpLVariable of a
	 * Variable.
	 *
	 * \param data The data from which to create an OpLVariable.
	 *
	 * \tparam Z The index of the Variable to create.
	 * \tparam A The type of the data to make an OpLVariable.
	 */
	template<size_t Z, typename A>
	auto make_op(A&& data)
	{
		return symphas::internal::make_oplvariable<A>::template get<Z>(std::forward<A>(data));
	}

	//! Delegate to create an OpLVariable with the correct template. 
	/*!
	 * Delegate to create an OpLVariable with the correct template. For
	 * more information on specific implementation, see symphas::internal::make_oplvariable.
	 *
	 * \param value The coefficient of the OpLVariable.
	 * \param data The data from which to create an OpLVariable.
	 *
	 * \param A The type of the data to make an OpLVariable.
	 */
	template<typename S, typename A>
	auto make_op(S&& value, A&& data)
	{
		return symphas::internal::make_oplvariable<A>::template get(std::forward<S>(value), std::forward<A>(data));
	}

	//! Delegate to create an OpLVariable with the correct template. 
	/*!
	 * Delegate to create an OpLVariable with the correct template. For
	 * more information on specific implementation, see symphas::internal::make_oplvariable.
	 * This overload will use the given index to create an OpLVariable of a
	 * Variable.
	 *
	 * \param value The coefficient of the OpLVariable.
	 * \param data The data from which to create an OpLVariable.
	 *
	 * \tparam Z The index of the Variable to create.
	 * \tparam A The type of the data to make an OpLVariable.
	 */
	template<size_t Z, typename S, typename A>
	auto make_op(S&& value, A&& data)
	{
		return symphas::internal::make_oplvariable<A>::template get<Z>(std::forward<S>(value), std::forward<A>(data));
	}


}

// *************************************************************************************

//! An expression for the left hand side of an equation.
/*!
 * In the context of the program implementation, the left hand side of an   
 * equation can only be a single data. There are no other expressions 
 * on the left hand side of an 
 * equals sign and no coefficients either, therefore the LHS object simply 
 * inherits from the data that it represents and implements the assignment
 * functions so that it can be combined with an expression.
 * 
 * \tparam G The type of the left hand side of the equation.
 */
template<typename G>
struct OpLHS : G
{
	using G::G;

	//! Construct the left hand side expression.
	/*!
	 * Construct an expression for the left hand side of an equals sign.
	 * 
	 * \param data The data that is stored on the left hand side.
	 */
	OpLHS(G const& data) : G(data) {}
	OpLHS(G&& data) noexcept : G(data) {}

	//! Combines this data with an expression.
	template<typename E>
	auto operator=(OpExpression<E>&& expr) const
	{
		return std::pair<G, E>(*this, *static_cast<E const*>(&expr));
	}

	//! Combines this data with an expression.
	template<typename E>
	auto operator=(OpExpression<E> const& expr) const
	{
		return std::pair<G, E>(*this, *static_cast<E const*>(&expr));
	}
};




// *************************************************************************************


//! A expression term that manages data, representing a variable.
/*!
 * Manages data to be used in an expression, representing a variable because
 * the underlying data can change.
 * The correct data type has
 * to be provided, requiring this object to be created with expr::make_op().
 * 
 * \tparam T The coefficient type of this term.
 * \tparam G The type of the data managed by this term.
 */
template<typename T, typename G>
struct OpLVariable : OpExpression<OpLVariable<T, G>>
{
	//! Generates a new variable term.
	/*!
	 * Generates a new variable term.
	 * 
	 * \param value The coefficient value.
	 * \param data The data managed by the variable.
	 */
	OpLVariable(T value, G data) : value{ value }, data{ data } {}

	inline auto eval(iter_type n) const
	{
		return value * expr::BaseData<G>::get(data, n);
	}

	template<size_t N>
	auto pow() const;

	auto operator-() const
	{
		auto v = -value;
		return OpLVariable<decltype(v), G>(-value, data);
	}

#ifdef PRINTABLE_EQUATIONS

	size_t print(FILE* out) const
	{
		return expr::print_with_coeff(out, expr::get_op_name(data), value);
	}

	size_t print(char* out) const
	{
		return expr::print_with_coeff(out, expr::get_op_name(data), value);
	}

	size_t print_length() const
	{
		return expr::coeff_print_length(value) + std::strlen(expr::get_op_name(data));
	}

#endif

	T value;			//!< Coefficient of the variable term.
	G data;				//!< The data which this variable represents.
};

//! A expression term that manages data, representing a variable.
/*!
 * Specialization based on having a coefficient equal to the multiplicative
 * identity.
 * 
 * \tparam G The type of the data managed by this term.
 */
template<typename G>
struct OpLVariable<OpIdentity, G> : OpExpression<OpLVariable<OpIdentity, G>>
{
	//! Generates a new variable term.
	/*!
	 * Generates a new variable term.
	 *
	 * \param data The data managed by the variable.
	 */
	OpLVariable(OpIdentity, G data) : data{ data } {}
	OpLVariable(G data) : data{ data } {}

	inline auto eval(iter_type n) const
	{
		return expr::BaseData<G>::get(data, n);
	}

	template<size_t N>
	auto pow() const;

	auto operator-() const
	{
		return OpLVariable<OpNegIdentity, G>(OpNegIdentity{}, data);
	}

#ifdef PRINTABLE_EQUATIONS

	size_t print(FILE* out) const
	{
		return expr::print_with_coeff(out, expr::get_op_name(data), OpIdentity{});
	}

	size_t print(char* out) const
	{
		return expr::print_with_coeff(out, expr::get_op_name(data), OpIdentity{});
	}

	size_t print_length() const
	{
		return expr::coeff_print_length(value) + std::strlen(expr::get_op_name(data));
	}

#endif

	/* when the opvariable is the special type that it has an identity as its constant, we
	 * can use it as an LHS
	 */

	template<typename E>
	auto operator=(OpExpression<E>&& expr) const
	{
		return (OpLHS(data) = expr);
	}

	template<typename E>
	auto operator=(OpExpression<E> const& expr) const
	{
		return (OpLHS(data) = expr);
	}

	OpIdentity value;	//!< Member is added to standardize member access across all OpLVariable specializations.
	G data;				//!< The data which this variable represents.
};

//! A expression term that manages data, equivalent to a variable.
/*!
 * Specialization based on the coefficient equal to -1.
 */
template<typename G>
struct OpLVariable<OpNegIdentity, G> : OpExpression<OpLVariable<OpNegIdentity, G>>
{
	//! Generates a new variable term.
	/*!
	 * Generates a new variable term.
	 *
	 * \param data The data managed by the variable.
	 */
	OpLVariable(OpNegIdentity, G data) : data{ data } {}
	OpLVariable(G data) : data{ data } {}

	inline auto eval(iter_type n) const
	{
		return -expr::BaseData<G>::get(data, n);
	}

	template<size_t N>
	auto pow() const;

	auto operator-() const
	{
		return OpLVariable<OpIdentity, G>(OpIdentity{}, data);
	}

#ifdef PRINTABLE_EQUATIONS

	size_t print(FILE* out) const
	{
		return expr::print_with_coeff(out, expr::get_op_name(data), OpNegIdentity{});
	}

	size_t print(char* out) const
	{
		return expr::print_with_coeff(out, expr::get_op_name(data), OpNegIdentity{});
	}

	size_t print_length() const
	{
		return expr::coeff_print_length(value) + std::strlen(expr::get_op_name(data));
	}

#endif

	/* when the opvariable is the special type that it has an identity as its constant, we
	 * can use it as an LHS
	 */

	template<typename E>
	auto operator=(OpExpression<E>&& expr) const
	{
		return (OpLHS(data) = expr);
	}

	template<typename E>
	auto operator=(OpExpression<E> const& expr) const
	{
		return (OpLHS(data) = expr);
	}

	OpNegIdentity value;	//!< Member is added to standardize member access across all OpLVariable specializations.
	G data;					//!< The data which this variable represents.
};


template<typename G>
OpLVariable(G)->OpLVariable<OpIdentity, G>;


//! Alias for variable type with identity coefficient.
template<typename G>
using OpILVariable = OpLVariable<OpIdentity, G>;

template<typename G2, typename T1, typename T2>
auto operator*(OpLiteral<T1> const& a, OpLVariable<T2, G2> const& b)
{
	return OpLVariable<mul_result_t<T1, T2>, G2>(a.value * b.value, b.data);
}


namespace symphas::internal
{

	template<size_t I>
	struct data_to_mul;

	template<>
	struct data_to_mul<0>
	{
		template<typename T, typename... Gs>
		auto operator()(T value, std::tuple<Gs...> const& datas)
		{
			return OpLVariable(value, std::get<0>(datas));
		}
	};

	template<size_t I>
	struct data_to_mul
	{
		template<typename T, typename... Gs>
		auto operator()(T value, std::tuple<Gs...> const& datas)
		{
			return expr::make_mul(data_to_mul<I - 1>(value, datas), OpLVariable(std::get<I>(datas)));
		}
	};


#ifdef PRINTABLE_EQUATIONS

	template<typename type, typename count, typename... Gs>
	struct print_op;

	template<typename... Gs>
	struct print_op<std::tuple<>, std::index_sequence<>, Gs...>
	{
		template<typename parent_type, typename os_type, size_t... Is>
		size_t operator()(parent_type* p, os_type* out, std::index_sequence<Is...>)
		{
			return 0;
		}
	};

	template<typename uniq0, typename... uniqs, size_t pow0, size_t... pows, typename... Gs>
	struct print_op<std::tuple<uniq0, uniqs...>, std::index_sequence<pow0, pows...>, Gs...>
	{
		template<typename parent_type, typename os_type>
		size_t operator()(parent_type* p, os_type* out, std::index_sequence<>)
		{
			return 0;
		}

		template<typename parent_type, typename os_type, size_t I0, size_t... Is>
		size_t operator()(parent_type* p, os_type* out, std::index_sequence<I0, Is...>)
		{
			using G0 = typename symphas::lib::type_of_index<I0, Gs...>::type;
			constexpr bool last_print = (sizeof...(uniqs) == 0 || sizeof...(Is) == 0);

			if constexpr (std::is_same<G0, uniq0>::value && expr::is_combinable<G0>)
			{
				auto [n, offset] = p->template print_op_apply<I0, pow0>(out);
				if constexpr (!last_print)
				{
					auto [n0, offset0] = p->print_sep_apply(out + offset);
					n += n0;
					offset += offset0;
				}
				return print_op<std::tuple<uniqs...>, std::index_sequence<pows...>, Gs...>{}(p, out + offset, std::index_sequence<Is...>{}) + n;
			}
			else if constexpr (std::is_same<G0, uniq0>::value)
			{
				auto [n, offset] = p->template print_op_apply<I0, 1>(out);
				if constexpr (!last_print)
				{
					auto [n0, offset0] = p->print_sep_apply(out + offset);
					n += n0;
					offset += offset0;
				}
				return print_op<std::tuple<uniq0, uniqs...>, std::index_sequence<pow0, pows...>, Gs...>{}(p, out + offset, std::index_sequence<Is...>{}) + n;
			}
			else if constexpr (expr::is_combinable<G0> && !((false || ... || std::is_same<G0, uniqs>::value)))
			{
				return print_op<std::tuple<uniq0, uniqs...>, std::index_sequence<pow0, pows...>, Gs...>{}(p, out, std::index_sequence<Is...>{});
			}
			else
			{
				return print_op<std::tuple<uniqs..., uniq0>, std::index_sequence<pows..., pow0>, Gs...>{}(p, out, std::index_sequence<I0, Is...>{});
			}
		}
	};

#endif

}



// *************************************************************************************

//! A expression term equivalent to the product of two or more OpLVariable.
/*!
 * Represents the product of two or more variables. Also known as a
 * multivariable.
 * Used as a standard way
 * to evaluate the product of multiple terms. Also used in the symbolic
 * algebra rules because it is preferred over 
 * 
 * \tparam T The coefficient of the nonlinear variable.
 * \tparam Gs... The data types that the variable represents.
 */
template<typename T, typename... Gs>
struct OpNLVariable : OpExpression<OpNLVariable<T, Gs...>>
{
	//! Create a multivariable from the given coefficient and data.
	/*!
	 * The multivariable will be created with the given coefficient value
	 * and manage the given list of data, and will represent the product
	 * of variables that represent that data.
	 * 
	 * \param value The coefficient.
	 * \param datas... The list of data that is put together as a product.
	 */
	OpNLVariable(T value, Gs... datas) : datas{ datas... }, value{ value } {}

	//! Create a multivariable from the given coefficient and data.
	/*!
	 * The multivariable will be created with the given coefficient value
	 * and manage the given list of data, and will represent the product
	 * of variables that represent that data.
	 *
	 * \param value The coefficient.
	 * \param datas The list of data that is put together as a product.
	 */
	OpNLVariable(T value, std::tuple<Gs...> const& datas) : datas{ datas }, value{ value } {}

	//! Create a multivariable from the given coefficient and data.
	/*!
	 * The multivariable will be created without a coefficient
	 * and manage the given list of data, and will represent the product
	 * of variables that represent that data.
	 *
	 * \param datas... The list of data that is put together as a product.
	 */
	OpNLVariable(Gs... datas) : datas{ datas... }, value{ OpIdentity{} } {}

	//! Create a multivariable from the given coefficient and data.
	/*!
	 * The multivariable will be created without a coefficient
	 * and manage the given list of data, and will represent the product
	 * of variables that represent that data.
	 *
	 * \param datas The list of data that is put together as a product.
	 */
	OpNLVariable(std::tuple<Gs...> const& datas) : datas{ datas }, value{ OpIdentity{} } {}

	inline auto eval(iter_type n) const
	{
		return _eval(n, std::make_index_sequence<sizeof...(Gs)>{});
	}
	auto operator-() const
	{
		auto v = -value;
		return OpNLVariable<decltype(v), Gs...>(-value, datas);
	}

#ifdef PRINTABLE_EQUATIONS

	size_t print(FILE* out) const
	{
		return _print(out, std::make_index_sequence<sizeof...(Gs)>{});
	}

	size_t print(char* out) const
	{
		return _print(out, std::make_index_sequence<sizeof...(Gs)>{});
	}

	size_t print_length() const
	{
		return _print_length(std::make_index_sequence<sizeof...(Gs)>{});
	}

#endif

	auto as_mul() const
	{
		return _as_mul();
	}

	std::tuple<Gs...> datas;
	T value;

protected:

	template<size_t I = 0>
	auto _as_mul() const
	{
		symphas::internal::data_to_mul<sizeof...(Gs) - 1>(value, datas);
	}

	template<size_t... Is>
	inline auto _eval(iter_type n, std::index_sequence<Is...>) const
	{
		return ((value * ... * expr::BaseData<Gs>::get(std::get<Is>(datas), n)));
	}

#ifdef PRINTABLE_EQUATIONS

	template<typename type, typename count, typename... G0s>
	friend struct symphas::internal::print_op;

	template<size_t I, size_t pow>
	std::pair<size_t, size_t> print_op_apply(FILE* out) const
	{
		if constexpr (pow == 1)
		{
			return { fprintf(out, "%s", expr::get_op_name(std::get<I>(datas))), 0 };
		}
		else
		{
			return { fprintf(out, "%s" SYEX_POW_SEP_A "%zd" SYEX_POW_SEP_B, expr::get_op_name(std::get<I>(datas)), pow), 0 };
		}
	}

	template<size_t I, size_t pow>
	std::pair<size_t, size_t> print_op_apply(char* out) const
	{
		if constexpr (pow == 1)
		{
			size_t n = sprintf(out, "%s", expr::get_op_name(std::get<I>(datas)));
			return { n, n };
		}
		else
		{
			size_t n = sprintf(out, "%s" SYEX_POW_SEP_A "%zd" SYEX_POW_SEP_B, expr::get_op_name(std::get<I>(datas)), pow);
			return { n, n };
		}
	}

	template<size_t I, size_t pow>
	std::pair<size_t, size_t> print_op_apply(int* out) const
	{
		if constexpr (I == 1)
		{
			return { std::strlen(expr::get_op_name(std::get<I>(datas))), 0 };
		}
		else
		{
			return { std::strlen(expr::get_op_name(std::get<I>(datas)))
				+ symphas::lib::num_digits<pow>()
				+ STR_ARR_LEN(SYEX_POW_SEP_A SYEX_POW_SEP_B) - 1,
				0 };
		}
	}


	inline std::pair<size_t, size_t> print_sep_apply(char* out) const
	{
		size_t n = sprintf(out, SYEX_NLV_SEP_OP);
		return { n, n };
	}

	inline std::pair<size_t, size_t> print_sep_apply(FILE* out) const
	{
		return { fprintf(out, SYEX_NLV_SEP_OP), 0 };
	}

	inline std::pair<size_t, size_t> print_sep_apply(int* out) const
	{
		return { STR_ARR_LEN(SYEX_NLV_SEP_OP) - 1, 0 };
	}



	using cc = symphas::lib::cc_like_types<Gs...>;
	using pr = symphas::internal::print_op<typename cc::type, typename cc::count, Gs...>;

	template<size_t... Is>
	size_t _print(FILE* out, std::index_sequence<Is...>) const
	{
		size_t n = 0;
		n += expr::print_with_coeff(out, value);
		n += pr{}(this, out, std::make_index_sequence<sizeof...(Gs)>{});
		return n;
	}

	template<size_t... Is>
	size_t _print(char* out, std::index_sequence<Is...>) const
	{
		size_t n = 0;
		n += expr::print_with_coeff(out, value);
		n += pr{}(this, out + n, std::make_index_sequence<sizeof...(Gs)>{});
		return n;
	}

	template<size_t... Is>
	size_t _print_length(std::index_sequence<Is...>) const
	{
		return expr::coeff_print_length(value)
			+ pr{}(this, (int*)(0), std::make_index_sequence<sizeof...(Gs)>{});
	}

#endif

};

template<typename... Gs>
OpNLVariable(Gs...)->OpNLVariable<OpIdentity, Gs...>;
template<typename... Gs>
OpNLVariable(std::tuple<Gs...>)->OpNLVariable<OpIdentity, Gs...>;


template<typename T1, typename T2, typename... E>
auto operator*(OpLiteral<T1> const& a, OpNLVariable<T2, E...> const& b)
{
	return OpNLVariable(a.value * b.value, b.datas);
}



template<typename T1, typename T2, typename... G1s, typename... G2s, 
	typename std::enable_if<expr::nl_can_multiply<G1s...>::template with<G2s...>::value, int>::type = 0>
auto operator*(OpNLVariable<T1, G1s...> const& a, OpNLVariable<T2, G2s...> const& b)
{
	return symphas::internal::NL_rules(
		OpNLVariable(
			a.value * b.value,
			symphas::internal::NL_sort(std::tuple_cat(a.datas, b.datas))
		));
}


template<typename T1, typename T2, typename G2, typename... Gs, 
	typename std::enable_if<expr::nl_can_multiply<Gs...>::template with<G2>::value, int>::type = 0>
auto operator*(OpNLVariable<T1, Gs...> const& a, OpLVariable<T2, G2> const& b)
{
	return symphas::internal::NL_rules(
		OpNLVariable(
			a.value * b.value, 
			symphas::internal::NL_sort(std::tuple_cat(a.datas, std::tuple<G2>(b.data)))
		));
}

template<typename T1, typename G1, typename T2, typename... Gs, 
	typename std::enable_if<expr::nl_can_multiply<G1>::template with<Gs...>::value, int>::type = 0>
auto operator*(OpLVariable<T1, G1> const& a, OpNLVariable<T2, Gs...> const& b)
{
	return symphas::internal::NL_rules(
		OpNLVariable(
			b.value * a.value, 
			symphas::internal::NL_sort(std::tuple_cat(std::tuple<G1>(a.data), b.datas))
		));
}



/*
 * self multiplication operator for linear variable transforming to nonlinear variable
 */


template<typename T1, typename G1, typename T2, typename G2, 
	typename std::enable_if<expr::m_idy<G1, G2>::enable, int>::type = 0>
auto operator*(OpLVariable<T1, G1> const& a, OpLVariable<T2, G2> const& b)
{
	return expr::m_idy<G1, G2>::template apply(a, b);
}

template<typename T1, typename G1, typename T2, typename G2, 
	typename std::enable_if<!expr::m_idy<G1, G2>::enable && expr::nl_can_multiply<G1>::template with<G2>::value, int>::type = 0>
auto operator*(OpLVariable<T1, G1> const& a, OpLVariable<T2, G2> const& b)
{
	return symphas::internal::NL_rules(
		OpNLVariable(
			a.value * b.value, 
			symphas::internal::NL_sort(symphas::lib::make_tuple(a.data, b.data))
		));
}



template<typename T, typename G>
template<size_t N>
auto OpLVariable<T, G>::pow() const
{
	return expr::pow<N>(*this);
}

template<typename G>
template<size_t N>
auto OpLVariable<OpIdentity, G>::pow() const
{
	return expr::pow<N>(*this);
}

template<typename G>
template<size_t N>
auto OpLVariable<OpNegIdentity, G>::pow() const
{
	return expr::pow<N>(*this);
}






