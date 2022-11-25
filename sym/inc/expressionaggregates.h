
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

#include <vector>

#include "expressions.h"
#include "expressionlogic.h"
#include "expressionsprint.h"

//! \cond
#ifdef EXPR_EXPORTS
#define DLLEXPR DLLEXPORT
#else
#define DLLEXPR DLLIMPORT
#endif
//! \endcond

namespace symphas::internal
{
	template<typename G>
	struct make_data
	{
		G operator()()
		{
			return { 0 };
		}
	};

	template<typename G>
	struct make_data<symphas::ref<G>>
	{
		symphas::ref<G> operator()()
		{
			G data;
			return std::ref(data);
		}
	};
}

//! Gives a name to data, to be used in expressions.
/*!
 * Adds a name that can be interpreted by the get_op_name()
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

	NamedData() : G(symphas::internal::make_data<G>{}() ), name{ "" } {}

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
		char* str = new char[(*static_cast<E const*>(&e)).print_length() + 1];
		(*static_cast<E const*>(&e)).print(str);
		name = std::string(str);
		delete[] str;
	}

	std::string name;	//!< Name given to the data.

#else

	NamedData() : G(data) {}

	template<typename T>
	NamedData(G data, T&&) : G(data) {}

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
	G* data;			//!< Pointer to the data.

#ifdef PRINTABLE_EQUATIONS

	NamedData() : data{}, name{ "" } {}

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

	NamedData() : data{} {}

	template<typename T>
	NamedData(G* data, T&&) : data{ data } {}

#endif


	operator G* ()
	{
		return data;
	}

	operator G* const () const
	{
		return data;
	}

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
	scalar_t data;		//!< Data value.

#ifdef PRINTABLE_EQUATIONS

	NamedData() : data{}, name{ "" } {}

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

	NamedData() : data{} {}

	template<typename T>
	NamedData(scalar_t data, T&&) : data{ data } {}

#endif

	operator scalar_t()
	{
		return data;
	}

	operator const scalar_t() const
	{
		return data;
	}

};

#ifdef PRINTABLE_EQUATIONS

template<typename G>
NamedData(G&, std::string)->NamedData<symphas::ref<G>>;
template<typename G>
NamedData(G&&, std::string)->NamedData<G>;
template<typename G>
NamedData(G*, std::string)->NamedData<G*>;


#else

template<typename G, typename T>
NamedData(G&, T&&)->NamedData<symphas::ref<G>>;
template<typename G, typename T>
NamedData(G&&, T&&)->NamedData<G>;
template<typename G, typename T>
NamedData(G*, T&&)->NamedData<G*>;

#endif



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

	Variable() : G() {}
	Variable(G const& data) : G(data) {}
	Variable(G&& data) noexcept : G(std::move(data)) {}
};


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

	//! The given parameter is created as a Variable.
	/*!
	 * The given parameter is created and returned as a Variable.
	 *
	 * \param g The object which is created as a Variable.
	 */
	template<Axis ax, typename G>
	auto as_component(G&& g)
	{
		return VectorComponent<ax, G>(std::forward<G>(g));
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
	template<Axis ax, typename G>
	auto as_component(G& g)
	{
		return VectorComponent<ax, symphas::ref<G>>(g);
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
	template<Axis ax, size_t Z, typename G>
	auto as_component(Variable<Z, G> g)
	{
		return Variable<Z, VectorComponent<ax, G>>(*static_cast<G*>(&g));
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

	//! Specialization based on SymbolID.
	template<Axis ax, typename G>
	const char* get_op_name(VectorComponent<ax, G> const& a)
	{
		static std::map<std::string, char *> map;
		const char* name = get_op_name(SymbolID<G>::get(*static_cast<G const*>(&a)));
		
		auto with_component = map.find(name);
		if (with_component == map.end())
		{
			map[name] = new char[std::strlen(name) + 3];

			std::strcpy(map[name], name);
			map[name][std::strlen(name)] = '_';
			map[name][std::strlen(name) + 1] = (ax == Axis::X) ? 'x' : (ax == Axis::Y) ? 'y' : (ax == Axis::Z) ? 'z' : '?';
			map[name][std::strlen(name) + 2] = '\0';
		}

		return map[name];
	};

	//! Specialization based on SymbolID.
	template<size_t Z, Axis ax, typename G>
	const char* get_op_name(Variable<Z, VectorComponent<ax, G>> const& a)
	{
		return get_op_name(*static_cast<VectorComponent<ax, G> const*>(&a));
	};

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



	template<typename T>
	auto get_fourier_name(T const& t)
	{
		char* name = expr::get_op_name(std::forward<T>(t));
		return std::string(SYEX_FT_OF_OP_FMT_A) + std::string(name) + std::string(SYEX_FT_OF_OP_FMT_B);
	}

#endif
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


template<typename G, expr::exp_key_t X>
struct Term : G
{
	using G::G;

	Term() : G() {}
	Term(G const& data) : G(data) {}
	Term(G&& data) noexcept : G(std::forward<G>(data)) {}

	inline auto eval(iter_type n) const
	{
		if constexpr (expr::is_symbol<G>)
		{
			return G{};
		}
		else if constexpr (expr::_Xk<X> == 1)
		{
			return expr::BaseData<G>::get(data(), n);
		}
		else
		{
			using expr::pow;
			using symphas::math::pow;

			if constexpr (expr::_Xk_t<X>::D == 1)
			{
				auto result = pow<expr::_Xk_t<X>::N>(expr::BaseData<G>::get(data(), n));
				if constexpr (expr::_Xk_t<X>::sign)
				{
					return 1.0 / result;
				}
				else
				{
					return result;
				}
			}
			else
			{
				return pow(expr::BaseData<G>::get(data(), n), expr::_Xk<X>);
			}
		}
	}

#ifdef PRINTABLE_EQUATIONS

	size_t print(FILE* out) const
	{
		if constexpr (expr::_Xk<X> == 1)
		{
			return fprintf(out, "%s", expr::get_op_name(data()));
		}
		else if constexpr (expr::_Xk_t<X>::D == 1)
		{
			return fprintf(out, "%s" SYEX_POW_SEP_A "%u" SYEX_POW_SEP_B, expr::get_op_name(data()), expr::_Xk_t<X>::N);
		}
		else
		{
			return fprintf(out, "%s" SYEX_POW_SEP_A "%u/%u" SYEX_POW_SEP_B, expr::get_op_name(data()), expr::_Xk_t<X>::N, expr::_Xk_t<X>::D);
		}
	}

	size_t print(char* out) const
	{
		if constexpr (expr::_Xk<X> == 1)
		{
			return sprintf(out, "%s", expr::get_op_name(*static_cast<G const*>(this)));
		}
		else if constexpr (expr::_Xk_t<X>::D == 1)
		{
			return sprintf(out, "%s" SYEX_POW_SEP_A "%u" SYEX_POW_SEP_B, expr::get_op_name(data()), expr::_Xk_t<X>::N);
		}
		else
		{
			return sprintf(out, "%s" SYEX_POW_SEP_A "%u/%u" SYEX_POW_SEP_B, expr::get_op_name(data()), expr::_Xk_t<X>::N, expr::_Xk_t<X>::D);
		}
	}

	size_t print_length() const
	{
		// + std::strlen(expr::get_op_name(data))
		size_t n = 0;
		if constexpr (expr::_Xk_t<X>::N > 1)
		{
			n += STR_ARR_LEN(SYEX_POW_SEP_A SYEX_POW_SEP_B)
				+ symphas::lib::num_digits<expr::_Xk_t<X>::N>();
			
			if constexpr (expr::_Xk_t<X>::D > 1)
			{
				n += symphas::lib::num_digits<expr::_Xk_t<X>::D>()
					+ 1;
			}
		}
		n += std::strlen(expr::get_op_name(data()));
		return n;
	}

#endif

	//! Return a new term with the opposite sign of the exponent.
	auto operator~() const;

	//! Apply the exponent, resulting in a new term.
	template<size_t N0>
	auto pow() const;

	//! Take the root of the given order, resulting in a new term.
	template<size_t N0>
	auto root() const;

	auto& data()
	{
		return *static_cast<G*>(this);
	}

	const auto& data() const
	{
		return *static_cast<G const*>(this);
	}
};

template<typename T, expr::exp_key_t X>
struct Term<T*, X> : Term<symphas::ref<T>, X>
{
	using parent_type = Term<symphas::ref<T>, X>;
	using parent_type::parent_type;

	Term(T* data) : parent_type(symphas::ref<T>(*data)) {}
};

namespace symphas::internal
{
	
	template<size_t N, typename E>
	struct Nth_type_of_terms;

	template<size_t N, typename T0, typename... Ts>
	struct Nth_type_of_terms<N, OpTerms<T0, Ts...>>
	{
		using type = typename Nth_type_of_terms<N - 1, OpTerms<Ts...>>::type;
	};

	template<size_t N>
	struct Nth_type_of_terms<N, OpTerms<>>
	{
		using type = OpIdentity;
	};

	template<typename T0, typename... Ts>
	struct Nth_type_of_terms<0, OpTerms<T0, Ts...>>
	{
		using type = OpTerms<T0, Ts...>;
	};


	template<typename estream, typename E>
	size_t print_one_term(estream* out, E const& e)
	{
		return expr::print_with_coeff(out, "", e);
	}

	template<typename estream, typename E>
	size_t print_one_term(estream* out, OpExpression<E> const& e)
	{
		return (*static_cast<E const*>(&e)).print(out);
	}

	template<typename estream, typename G, expr::exp_key_t X>
	size_t print_one_term(estream* out, Term<G, X> const& term)
	{
		return term.print(out);
	}

	template<typename E>
	size_t print_length_one_term(OpExpression<E> const& e)
	{
		return expr::coeff_print_length(*static_cast<E const*>(&e));
	}

	template<typename G, expr::exp_key_t X>
	size_t print_length_one_term(Term<G, X> const& term)
	{
		return term.print_length();
	}
}


template<typename G>
Term(G)->Term<G, expr::Xk<1>>;
template<typename G>
Term(G*)->Term<symphas::ref<G*>, expr::Xk<1>>;

namespace expr
{

	template<size_t N, typename V, typename... Gs, exp_key_t... Xs>
	const auto& get(OpTerms<V, Term<Gs, Xs>...> const& e)
	{
		using Nth_type = typename symphas::internal::Nth_type_of_terms<N, OpTerms<V, Term<Gs, Xs>...>>::type;
		return static_cast<Nth_type const*>(&e)->term;
	}

	template<size_t N, typename V, typename... Gs, exp_key_t... Xs>
	auto& get(OpTerms<V, Term<Gs, Xs>...>& e)
	{
		using Nth_type = typename symphas::internal::Nth_type_of_terms<N, OpTerms<V, Term<Gs, Xs>...>>::type;
		return static_cast<Nth_type*>(&e)->term;
	}

	template<typename V, typename... Gs, exp_key_t... Xs>
	const auto& terms_after_first(OpTerms<V, Term<Gs, Xs>...> const& e)
	{
		return *static_cast<OpTerms<Term<Gs, Xs>...> const*>(&e);
	}

	template<typename V, typename... Gs, exp_key_t... Xs>
	auto& terms_after_first(OpTerms<V, Term<Gs, Xs>...>& e)
	{
		return *static_cast<OpTerms<Term<Gs, Xs>...>*>(&e);
	}

	template<typename V, typename G, exp_key_t X>
	const auto& data(OpTerms<V, Term<G, X>> const& e)
	{
		return static_cast<OpTerms<Term<G, X>> const*>(&e)->term.data();
	}

	template<typename V, typename G, exp_key_t X>
	auto& data(OpTerms<V, Term<G, X>>& e)
	{
		return static_cast<OpTerms<Term<G, X>>*>(&e)->term.data();
	}

	template<typename E>
	constexpr bool is_term = false;

	template<typename G, exp_key_t X>
	constexpr bool is_term<Term<G, X>> = true;
}


template<typename G, expr::exp_key_t X1, expr::exp_key_t X2, typename std::enable_if_t<!expr::has_identity<G, G>, int> = 0>
auto operator*(Term<G, X1> const& a, Term<G, X2> const& b)
{
	return Term<G, expr::XXk_t<X1, X2>::value>(*static_cast<G const*>(&a));
}


template<typename G, expr::exp_key_t X>
auto Term<G, X>::operator~() const
{
	constexpr expr::exp_key_t X2 = X ^ symphas::internal::xsm;
	return Term<G, X2>(*static_cast<G const*>(this));
}

template<typename G, expr::exp_key_t X>
template<size_t N0>
auto Term<G, X>::pow() const
{
	if constexpr (N0 == 0)
	{
		return (*this) * ~(*this);
	}
	else if constexpr (N0 == 1)
	{
		return *this;
	}
	else if constexpr (N0 == 2)
	{
		return (*this) * (*this);
	}
	else
	{
		constexpr size_t N1 = N0 / 2;
		if constexpr (N0 % 2 == 1)
		{
			return pow<N1>() * pow<N1>() * (*this);
		}
		else
		{
			return pow<N1>() * pow<N1>();
		}
	}
}

template<typename G, expr::exp_key_t X>
template<size_t N0>
auto Term<G, X>::root() const
{
	if constexpr (N0 == 0)
	{
		return (*this) * ~(*this);
	}
	else if constexpr (N0 == 1)
	{
		return *this;
	}
	else if constexpr (N0 == 2)
	{
		return (*this) * (*this);
	}
	else
	{
		constexpr size_t N1 = N0 / 2;
		if constexpr (N0 % 2 == 1)
		{
			return pow<N1>() * pow<N1>() * (*this);
		}
		else
		{
			return pow<N1>() * pow<N1>();
		}
	}
}

template<>
struct OpTerms<> : OpExpression<OpTerms<>>
{
#ifdef PRINTABLE_EQUATIONS

	auto eval(iter_type n) const
	{
		return OpIdentity{};
	}

	size_t print(FILE* out) const
	{
		return 0;
	}

	size_t print(char* out) const
	{
		return 0;
	}

	size_t print_length() const
	{
		return 0;
	}

#endif
};

template<typename V, typename... Gs, expr::exp_key_t... Xs>
struct OpExpression<OpTerms<V, Term<Gs, Xs>...>> : OpTerms<Term<Gs, Xs>...>
{
	using parent_type = OpTerms<Term<Gs, Xs>...>;
	using parent_type::parent_type;
	using E = OpTerms<V, Term<Gs, Xs>...>;

	OpExpression(OpTerms<Term<Gs, Xs>...> const& terms) : parent_type(terms) {}

	auto operator()(iter_type n) const
	{
		return cast().eval(n);
	}

	template<typename E0>
	auto operator()(OpExpression<E0> const& e) const
	{
		return cast() * (*static_cast<E0 const*>(&e));
	}

	auto& cast() const
	{
		return *static_cast<E const*>(this);
	}

	symphas::internal::expression_iterator<E> begin() const
	{
		return symphas::internal::expression_iterator<E>(cast());
	}

	symphas::internal::expression_iterator<E> end(len_type len) const
	{
		return symphas::internal::expression_iterator<E>(cast(), len);
	}
};


//! A expression term equivalent to the product of two or more OpTerm.
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
template<typename V, typename... Gs, expr::exp_key_t... Xs>
struct OpTerms<V, Term<Gs, Xs>...> : OpExpression<OpTerms<V, Term<Gs, Xs>...>>
{
	using parent_type = OpExpression<OpTerms<V, Term<Gs, Xs>...>>;

	OpTerms() : parent_type(), term{ V{} } {}

	OpTerms(V const& term0, Term<Gs, Xs> const&... terms) 
		: parent_type(terms...), term{ term0 } {}
	
	OpTerms(V const& term0, OpTerms<Term<Gs, Xs>...> const& terms) 
		: parent_type(terms), term{ term0 } {}
	
	OpTerms(OpTerms<V, Term<Gs, Xs>...> const& terms) 
		: parent_type(expr::terms_after_first(terms)), term{ terms.term } {}

	template<typename data_type = V, typename nested_type = OpTerms<Term<Gs, Xs>...>,
		typename std::enable_if_t<(
			std::is_default_constructible<data_type>::value &&
			std::is_default_constructible<nested_type>::value), int> = 0>
	OpTerms() : parent_type(), term{} {}

	template<typename E = V, typename std::enable_if_t<expr::is_term<E>, int> = 0>
	auto _eval(iter_type n) const
	{
		return term.eval(n) * parent_type::eval(n);
	}

	template<typename E = V, typename std::enable_if_t<!expr::is_term<E>, int> = 0>
	auto _eval(iter_type n) const
	{
		return expr::eval(term) * parent_type::eval(n);
	}

	inline auto eval(iter_type n) const
	{
		return _eval(n);
	}

	template<typename E, size_t N = sizeof...(Gs), typename std::enable_if_t<(N == 1), int> = 0>
	void operator=(OpExpression<E> const& e)
	{
		expr::result(
			*static_cast<E const*>(&e),
			static_cast<OpTerms<Term<Gs, Xs>...>*>(this)->term);
	}
	
	auto operator-() const;

#ifdef PRINTABLE_EQUATIONS

	size_t print(FILE* out) const
	{
		size_t n = symphas::internal::print_one_term(out, term);
		return n + parent_type::print(out);
	}

	size_t print(char* out) const
	{
		size_t n = symphas::internal::print_one_term(out, term);
		return n + parent_type::print(out + n);
	}

	size_t print_length() const
	{
		size_t n = 0;
		if constexpr (expr::is_term<V>)
		{
			n += symphas::internal::print_length_one_term(term);
		}
		else
		{
			n += symphas::internal::print_length_one_term(expr::make_literal(term));
		}
		return n + parent_type::print_length();
	}

#endif

	V term;
};



template<typename V, typename... Gs, expr::exp_key_t... Xs>
OpTerms(V, Term<Gs, Xs>...)->OpTerms<V, Term<Gs, Xs>...>;

template<typename V, typename... Gs, expr::exp_key_t... Xs>
OpTerms(OpLiteral<V>, Term<Gs, Xs>...)->OpTerms<V, Term<Gs, Xs>...>;

template<typename V, typename... Gs, expr::exp_key_t... Xs>
OpTerms(V, OpTerms<Term<Gs, Xs>...>)->OpTerms<V, Term<Gs, Xs>...>;

template<typename V, typename... Gs, expr::exp_key_t... Xs>
OpTerms(OpLiteral<V>, OpTerms<Term<Gs, Xs>...>)->OpTerms<V, Term<Gs, Xs>...>;

template<typename V, typename... Gs, expr::exp_key_t... Xs, typename G0, expr::exp_key_t X0>
OpTerms(OpTerms<V, Term<Gs, Xs>...>, Term<G0, X0>)->OpTerms<V, Term<G0, X0>, Term<Gs, Xs>...>;



template<typename V, typename... Gs, expr::exp_key_t... Xs>
auto OpTerms<V, Term<Gs, Xs>...>::operator-() const
{
	return ::OpTerms(-term, *static_cast<OpTerms<Term<Gs, Xs>...> const*>(this));
}



namespace symphas::internal
{
	template<expr::exp_key_t X, typename G>
	auto to_term_element(G const& data)
	{
		return Term<G, X>(data);
	}

	template<typename G>
	auto to_term_element(G const& data)
	{
		return to_term_element<expr::Xk<1>>(data);
	}

	template<size_t Z, typename G>
	auto to_variable_element(G const& data)
	{
		return Variable<Z, G>(data);
	}


	template<Axis ax, typename G>
	auto to_vector_component(G const& data)
	{
		return VectorComponent<ax, G>(*static_cast<G const*>(&data));
	}

	template<Axis ax, size_t Z, typename G>
	auto to_vector_component(Variable<Z, G> const& data)
	{
		return expr::as_variable<Z>(to_vector_component<ax>(*static_cast<G const*>(&data)));;
	}

	template<Axis ax, typename G, expr::exp_key_t X>
	auto to_vector_component(Term<G, X> const& term)
	{
		return to_term_element<X>(to_vector_component<ax>(*static_cast<G const*>(&term)));
	}

	template<Axis ax, typename... Gs, expr::exp_key_t... Xs, size_t... Ns, size_t M0, size_t... Ms>
	auto terms_make_component(OpTerms<Term<Gs, Xs>...> const& b, std::index_sequence<Ns...>, std::index_sequence<M0, Ms...>)
	{
		return OpTerms(OpIdentity{}, expr::get<Ns>(b)..., to_vector_component<ax>(expr::get<M0>(b)), expr::get<Ms>(b)...);
	}

	template<Axis ax, typename... Gs, expr::exp_key_t... Xs, size_t... Is>
	auto terms_make_component(OpTerms<Term<Gs, Xs>...> const& b, std::index_sequence<Is...>)
	{
		using mask_t = symphas::lib::seq_join_t<
			std::index_sequence<>,
			std::conditional_t<
				(expr::eval_type<Term<Gs, Xs>>::rank == 0),
				std::index_sequence<>,
				std::index_sequence<Is>>
			...>;
		constexpr size_t M0 = symphas::lib::seq_index_value<mask_t::size() - 1, mask_t>::value;

		return terms_make_component<ax>(b, 
			std::make_index_sequence<M0>{}, 
			symphas::lib::seq_add_t<
				symphas::lib::seq_repeating_value_t<sizeof...(Gs) - M0, size_t, M0>, 
				std::make_index_sequence<sizeof...(Gs) - M0>
			>{});
	}
}

template<typename coeff_t, typename V, typename... Gs, expr::exp_key_t... Xs, 
	typename std::enable_if_t<(expr::is_coeff<coeff_t> && !expr::is_tensor<V>), int> = 0>
auto operator*(coeff_t const& value, OpTerms<V, Term<Gs, Xs>...> const& b)
{
	return OpTerms(value * expr::coeff(b), *static_cast<OpTerms<Term<Gs, Xs>...> const*>(&b));
}

template<typename coeff_t, typename tensor_t, typename... Gs, expr::exp_key_t... Xs,
	typename std::enable_if_t<(expr::is_coeff<coeff_t> && expr::is_tensor<tensor_t>), int> = 0>
auto operator*(coeff_t const& value, OpTerms<tensor_t, Term<Gs, Xs>...> const& b)
{
	return (value * expr::coeff(b)) * OpTerms(OpIdentity{}, *static_cast<OpTerms<Term<Gs, Xs>...> const*>(&b));
}


template<typename T, typename V, typename... Gs, expr::exp_key_t... Xs,
	size_t N0, size_t N1, size_t N,
	size_t D = expr::eval_type<OpTerms<Term<Gs, Xs>...>>::rank,
	typename std::enable_if_t<(!expr::is_tensor<V> && expr::eval_type<OpTerms<Term<Gs, Xs>...>>::rank > 0), int> = 0>
auto operator*(OpTensor<T, N0, N1, N, D> const& tensor, OpTerms<V, Term<Gs, Xs>...> const& b)
{
	constexpr Axis ax = (N1 == 0) ? Axis::X : (N1 == 1) ? Axis::Y : Axis::Z;
	if constexpr (N == 1)
	{
		auto coeff = symphas::internal::tensor_cast::cast(tensor) * expr::coeff(b);
		return coeff * symphas::internal::terms_make_component<ax>(
			*static_cast<OpTerms<Term<Gs, Xs>...> const*>(&b), std::make_index_sequence<sizeof...(Gs)>{});
	}
	else
	{
		auto coeff = expr::make_tensor<N0, N>(symphas::internal::tensor_cast::cast(tensor)) * expr::coeff(b);
		return coeff * symphas::internal::terms_make_component<ax>(
			*static_cast<OpTerms<Term<Gs, Xs>...> const*>(&b), std::make_index_sequence<sizeof...(Gs)>{});
	}

}


template<typename T, typename V, typename... Gs, expr::exp_key_t... Xs,
	size_t N0, size_t N1, size_t N, size_t M,
	size_t D = expr::eval_type<OpTerms<Term<Gs, Xs>...>>::dimension,
	typename = std::enable_if_t<(D > 1 && M != D), int>>
auto operator*(OpTensor<T, N0, N1, N, M> const& tensor, OpTerms<V, Term<Gs, Xs>...> const& b) = delete;




template<typename E, typename G0, expr::exp_key_t X0, typename... Gs, expr::exp_key_t... Xs>
auto operator*(OpExpression<E> const& a, OpTerms<Term<G0, X0>, Term<Gs, Xs>...> const& b)
{
	return (*static_cast<E const*>(&a)) * OpTerms(OpIdentity{}, b);
}

template<typename E, typename G0, expr::exp_key_t X0, typename... Gs, expr::exp_key_t... Xs>
auto operator*(OpTerms<Term<G0, X0>, Term<Gs, Xs>...> const& a, OpExpression<E> const& b)
{
	return OpTerms(OpIdentity{}, a) * (*static_cast<E const*>(&b));
}


namespace symphas::internal
{

	template<typename V, typename... Gs, expr::exp_key_t... Xs, size_t... Is>
	auto invert_terms(OpTerms<V, Term<Gs, Xs>...> const& b, std::index_sequence<Is...>)
	{
		return OpTerms(~expr::get<Is>(b)...);
	}

	inline auto combine_terms(std::index_sequence<>, OpTerms<> const& terms_num, OpTerms<> const& terms_den)
	{
		return OpIdentity{};
	}

	template<typename... Gs, expr::exp_key_t... Xs, typename G, expr::exp_key_t X, size_t... P0s, size_t... P1s>
	auto place_at(OpTerms<Term<Gs, Xs>...> const& terms, Term<G, X> const& term0, std::index_sequence<P0s...>, std::index_sequence<P1s...>)
	{
		return OpTerms(expr::get<P0s>(terms)..., term0, expr::get<P1s + sizeof...(P0s)>(terms)...);
	}

	template<size_t P, typename... Gs, expr::exp_key_t... Xs, typename G, expr::exp_key_t X, size_t... Is>
	auto place_at(OpTerms<Term<Gs, Xs>...> const& terms, Term<G, X> const& term0)
	{
		return place_at(terms, term0, std::make_index_sequence<P>{}, std::make_index_sequence<sizeof...(Gs) - P>{});
	}

	template<typename G02, expr::exp_key_t X02, typename... G2s, expr::exp_key_t... X2s>
	auto combine_terms(std::index_sequence<>, OpTerms<> const& terms_num, OpTerms<Term<G02, X02>, Term<G2s, X2s>...> const& terms_den)
	{
		return expr::make_div(OpIdentity{}, OpTerms(OpIdentity{}, terms_den));
	}

	template<typename G01, expr::exp_key_t X01, typename... G1s, expr::exp_key_t... X1s>
	auto combine_terms(std::index_sequence<>, OpTerms<Term<G01, X01>, Term<G1s, X1s>...> const& terms_num, OpTerms<> const& terms_den)
	{
		return OpTerms(OpIdentity{}, terms_num);
	}

	template<typename G01, expr::exp_key_t X01, typename... G1s, expr::exp_key_t... X1s, typename G02, expr::exp_key_t X02, typename... G2s, expr::exp_key_t... X2s>
	auto combine_terms(std::index_sequence<>, OpTerms<Term<G01, X01>, Term<G1s, X1s>...> const& terms_num, OpTerms<Term<G02, X02>, Term<G2s, X2s>...> const& terms_den)
	{
		return OpBinaryDiv(OpTerms(OpIdentity{}, terms_num), OpTerms(OpIdentity{}, terms_den));
	}

	template<size_t P0, size_t... Ps, typename... G1s, expr::exp_key_t... X1s, typename... G2s, expr::exp_key_t... X2s, typename G2, typename... Ts>
	auto combine_terms(std::index_sequence<P0, Ps...>, OpTerms<Term<G1s, X1s>...> const& terms_num, OpTerms<Term<G2s, X2s>...> const& terms_den, Term<G2, 0> const& term0, Ts&&... rest)
	{
		return combine_terms(std::index_sequence<Ps...>{}, terms_num, terms_den, std::forward<Ts>(rest)...);
	}

	template<size_t P0, size_t... Ps, typename... G1s, expr::exp_key_t... X1s, typename... G2s, expr::exp_key_t... X2s, typename G2, expr::exp_key_t X2, typename... Ts>
	auto combine_terms(std::index_sequence<P0, Ps...>, OpTerms<Term<G1s, X1s>...> const& terms_num, OpTerms<Term<G2s, X2s>...> const& terms_den, Term<G2, X2> const& term0, Ts&&... rest)
	{
		if constexpr (expr::_Xk_t<X2>::sign)
		{
			return combine_terms(std::index_sequence<Ps...>{}, terms_num, OpTerms(~term0, terms_den), std::forward<Ts>(rest)...);
		}
		else
		{
			return combine_terms(std::index_sequence<Ps...>{}, place_at<P0>(terms_num, term0), terms_den, std::forward<Ts>(rest)...);
		}
	}

	template<typename... G1s, expr::exp_key_t... X1s, typename... G2s, expr::exp_key_t... X2s, size_t... Is, size_t... Js, size_t... Ms, size_t... Ns>
	auto combine_terms(OpTerms<Term<G1s, X1s>...> const& a, OpTerms<Term<G2s, X2s>...> const& b, 
		std::index_sequence<Is...>, std::index_sequence<Js...>, 
		std::index_sequence<Ms...>, std::index_sequence<Ns...>)
	{
		if constexpr ((expr::_Xk_t<X2s>::sign && ...))
		{
			return combine_terms(std::index_sequence<Ms...>{}, OpTerms(expr::get<Is>(a)...), OpTerms(~expr::get<Js>(b)...), (expr::get<Ms>(a) * expr::get<Ns>(b))...);
		}
		else
		{
			return combine_terms(std::index_sequence<Ms...>{}, OpTerms(expr::get<Is>(a)..., expr::get<Js>(b)...), OpTerms<>{}, (expr::get<Ms>(a) * expr::get<Ns>(b))...);
		}
	}

	template<typename G, typename T1>
	constexpr bool commutes_through = false;

	template<typename G, typename... Gs>
	constexpr bool commutes_through<G, symphas::lib::types_list<Gs...>> = (expr::is_commutable<G, Gs> && ...);

	template<typename G0, typename G1, typename T>
	constexpr bool commutes_to = false;

	//! Indicates whether G0 can reach a G1 in Gs list through commuting.
	/*!
	 * Returns true if the `G0` type can reach the `G1` type by commuting through
	 * the `Gs` list. Two conditions must effectively be met:
	 * 1. `G1` must exist in `Gs`; if not then this will evaluate to false;
	 * 2. `G0` must commute with everything in `Gs`, up to the occurence of a `G1`. 
	 * If both conditions are met, then this statement will evaluate true.
	 */
	template<typename G0, typename G1, typename... Gs>
	constexpr bool commutes_to<G0, G1, types_list<Gs...>> = 
		(index_of_type<G1, Gs...> >= 0) 
		? commutes_through<G0, types_before_index<fixed_max<0, index_of_type<G1, Gs...>>, Gs...>>
		: false;

	template<typename G1, typename T1, typename T2>
	constexpr bool satisfies_combination = false;

	template<typename G1, typename... Gs, typename... G2s>
	constexpr bool satisfies_combination<G1, types_list<Gs...>, types_list<G2s...>> = (
		(expr::is_combinable<G1> && commutes_to<G1, G1, types_list<Gs..., G2s...>>)	// will combine with another G1 in G2s list
		|| ((expr::has_identity<G1, G2s> && commutes_to<G1, G2s, types_list<Gs..., G2s...>>) && ...)); // will combine with a G2 if there is an identity


	template<typename... G1s, expr::exp_key_t... X1s, typename... G2s, expr::exp_key_t... X2s, size_t... Is, size_t... Js>
	auto combine_terms(OpTerms<Term<G1s, X1s>...> const& a, OpTerms<Term<G2s, X2s>...> const& b, std::index_sequence<Is...>, std::index_sequence<Js...>)
	{
		using seq_a_mask = std::integer_sequence<bool, (
			satisfies_combination<G1s, 
				types_after_at_index<Is + 1, G1s...>, 
				types_list<G2s...>>
			)...>;
		using seq_b_mask = std::integer_sequence<bool, (
			satisfies_combination<G2s,
				reverse_types_list<types_before_index<Js, G2s...>>,
				reverse_types_list<types_list<G1s...>>>
			)...>;

		using seq_a_pick = symphas::lib::seq_join_t<std::index_sequence<>, std::conditional_t<symphas::lib::seq_index_value<Is, seq_a_mask>::value, std::index_sequence<Is>, std::index_sequence<>>...>;
		using seq_a_pick_ = symphas::lib::seq_join_t<std::index_sequence<>, std::conditional_t<!symphas::lib::seq_index_value<Is, seq_a_mask>::value, std::index_sequence<Is>, std::index_sequence<>>...>;

		using seq_b_pick = symphas::lib::seq_join_t<std::index_sequence<>, std::conditional_t<symphas::lib::seq_index_value<Js, seq_b_mask>::value, std::index_sequence<Js>, std::index_sequence<>>...>;
		using seq_b_pick_ = symphas::lib::seq_join_t<std::index_sequence<>, std::conditional_t<!symphas::lib::seq_index_value<Js, seq_b_mask>::value, std::index_sequence<Js>, std::index_sequence<>>...>;

		return combine_terms(a, b, seq_a_pick_{}, seq_b_pick_{}, seq_a_pick{}, seq_b_pick{});
	}


}


//! Addition of two multi variables with data that can be combined.
template<typename V1, typename... G1s, expr::exp_key_t... X1s, typename V2, typename... G2s, expr::exp_key_t... X2s>
auto operator*(OpTerms<V1, Term<G1s, X1s>...> const& a, OpTerms<V2, Term<G2s, X2s>...> const& b)
{
	return expr::coeff(a) * expr::coeff(b) *
		symphas::internal::combine_terms(
			expr::terms_after_first(a), 
			expr::terms_after_first(b), 
			std::make_index_sequence<sizeof...(G1s)>{}, 
			std::make_index_sequence<sizeof...(G2s)>{});
}


//! Addition of two multi variables with data that can be combined.
template<typename V1, typename... G1s, expr::exp_key_t... X1s, typename V2, typename... G2s, expr::exp_key_t... X2s>
auto operator/(OpTerms<V1, Term<G1s, X1s>...> const& a, OpTerms<V2, Term<G2s, X2s>...> const& b)
{
	return expr::coeff(a) * (OpIdentity{} / expr::coeff(b)) *
		symphas::internal::combine_terms(
			expr::terms_after_first(a), 
			symphas::internal::invert_terms(expr::terms_after_first(b), std::make_index_sequence<sizeof...(G2s)>{}),
			std::make_index_sequence<sizeof...(G1s)>{},
			std::make_index_sequence<sizeof...(G2s)>{});
}


namespace symphas::internal
{

	using expr::Xk;

	template<typename A>
	struct construct_term
	{
		template<typename S>
		static auto get(S value, A data)
		{
			if constexpr (expr::is_coeff<A>)
			{
				return value * data;
			}
			else
			{
				return OpTerms(value, Term<A, Xk<1>>(data));
			}
		}

		static auto get(A data)
		{
			if constexpr (expr::is_coeff<A>)
			{
				return data;
			}
			else
			{
				return OpTerms(OpIdentity{}, Term<A, Xk<1>>(data));
			}
		}

		template<size_t Z, typename S>
		static auto get(S value, A data)
		{
			if constexpr (expr::is_coeff<A>)
			{
				return value * data;
			}
			else
			{
				return OpTerms(value, Term<Variable<Z, A>, Xk<1>>(data));
			}
		}

		template<size_t Z>
		static auto get(A data)
		{
			if constexpr (expr::is_coeff<A>)
			{
				return data;
			}
			else
			{
				return OpTerms(OpIdentity{}, Term<Variable<Z, A>, Xk<1>>(data));
			}
		}
	};

	template<>
	struct construct_term<scalar_t>
	{
		template<typename S>
		static auto get(S value, scalar_t data)
		{
			return OpTerms(value, Term<OpLiteral<scalar_t>, Xk<1>>(data));
		}

		static auto get(scalar_t data)
		{
			return OpTerms(OpIdentity{}, Term<OpLiteral<scalar_t>, Xk<1>>(data));
		}

		template<size_t Z, typename S>
		static auto get(S value, scalar_t data)
		{
			return OpTerms(value, Term<Variable<Z, OpLiteral<scalar_t>>, Xk<1>>(data));
		}

		template<size_t Z>
		static auto get(scalar_t data)
		{
			return OpTerms(OpIdentity{}, Term<Variable<Z, OpLiteral<scalar_t>>, Xk<1>>(data));
		}
	};

	template<>
	struct construct_term<complex_t>
	{
		template<typename S>
		static auto get(S value, complex_t data)
		{
			return OpTerms(value, Term<OpLiteral<complex_t>, Xk<1>>(data));
		}

		static auto get(complex_t data)
		{
			return OpTerms(OpIdentity{}, Term<OpLiteral<complex_t>, Xk<1>>(data));
		}

		template<size_t Z, typename S>
		static auto get(S value, complex_t data)
		{
			return OpTerms(value, Term<Variable<Z, OpLiteral<complex_t>>, Xk<1>>(data));
		}

		template<size_t Z>
		static auto get(complex_t data)
		{
			return OpTerms(OpIdentity{}, Term<Variable<Z, OpLiteral<complex_t>>, Xk<1>>(data));
		}
	};

	template<typename T, size_t D>
	struct construct_term<any_vector_t<T, D>>
	{
		template<typename S>
		static auto get(S value, any_vector_t<T, D> data)
		{
			return OpTerms(value, Term<OpLiteral<any_vector_t<T, D>>, Xk<1>>(data));
		}

		static auto get(any_vector_t<T, D> data)
		{
			return OpTerms(OpIdentity{}, Term<OpLiteral<any_vector_t<T, D>>, Xk<1>>(data));
		}

		template<size_t Z, typename S>
		static auto get(S value, any_vector_t<T, D> data)
		{
			return OpTerms(value, Term<Variable<Z, OpLiteral<any_vector_t<T, D>>>, Xk<1>>(data));
		}

		template<size_t Z>
		static auto get(any_vector_t<T, D> data)
		{
			return OpTerms(OpIdentity{}, Term<Variable<Z, OpLiteral<any_vector_t<T, D>>>, Xk<1>>(data));
		}
	};


	//! Helper in order to construct an OpTerms.
	/*!
	 * Specialization for an lvalue; this means that a reference should be
	 * created.
	 */
	template<typename A>
	struct construct_term<A&>
	{
		template<typename S>
		static auto get(S value, A& data)
		{
			return construct_term<symphas::ref<A>>::get(value, data);
		}

		static auto get(A& data)
		{
			return construct_term<symphas::ref<A>>::get(data);
		}

		template<size_t Z, typename S>
		static auto get(S value, A& data)
		{
			return construct_term<symphas::ref<A>>::template get<Z>(value, data);
		}

		template<size_t Z>
		static auto get(A& data)
		{
			return construct_term<symphas::ref<A>>::template get<Z>(data);
		}
	};

	//! Helper in order to construct an OpTerms.
	/*!
	 * Specialization for a pointer; this means the pointer should be copied to
	 * the variable.
	 */
	template<typename A>
	struct construct_term<A*>
	{
		template<typename S>
		static auto get(S value, A* data)
		{
			return OpTerms(value, Term<A*, Xk<1>>(data));
		}

		static auto get(A* data)
		{
			return OpTerms(OpIdentity{}, Term<A*, Xk<1>>(data));
		}

		template<size_t Z, typename S>
		static auto get(S value, A* data)
		{
			return OpTerms(value, Term<Variable<Z, A*>, Xk<1>>(data));
		}

		template<size_t Z>
		static auto get(A* data)
		{
			return OpTerms(OpIdentity{}, Term<Variable<Z, A*>, Xk<1>>(data));
		}
	};

	//! Helper in order to construct an OpTerms.
	/*!
	 * The data is passed as a reference, rather than copied.
	 */
	template<size_t Z, typename A>
	struct construct_term<Variable<Z, A>>
	{
		template<typename S>
		static auto get(S value, Variable<Z, A>& data)
		{
			return construct_term<A>::template get<Z>(value, data);
		}

		static auto get(Variable<Z, A>& data)
		{
			return construct_term<A>::template get<Z>(data);
		}

		template<size_t Y, typename S>
		static auto get(S value, Variable<Z, A>& data)
		{
			return construct_term<A>::template get<Y>(value, data);
		}

		template<size_t Y>
		static auto get(Variable<Z, A>& data)
		{
			return construct_term<A>::template get<Y>(data);
		}
	};

	//! Helper in order to construct an OpTerms.
	/*!
	 * The data is passed as a reference, rather than copied.
	 */
	template<typename A>
	struct construct_term<symphas::ref<A>>
	{
		template<typename S>
		static auto get(S value, symphas::ref<A> data)
		{
			return OpTerms(value, Term<symphas::ref<A>, Xk<1>>(data));
		}

		static auto get(symphas::ref<A> data)
		{
			return OpTerms(OpIdentity{}, Term<symphas::ref<A>, Xk<1>>(data));
		}

		template<size_t Z, typename S>
		static auto get(S value, symphas::ref<A> data)
		{
			return OpTerms(value, Term<Variable<Z, symphas::ref<A>>, Xk<1>>(data));
		}

		template<size_t Z>
		static auto get(symphas::ref<A> data)
		{
			return OpTerms(OpIdentity{}, Term<Variable<Z, symphas::ref<A>>, Xk<1>>(data));
		}
	};


	//! Helper in order to construct an OpTerms.
	/*!
	 * When these types are given by lvalue, then the original functions need to
	 * be applied, not the functions corresponding to the simple lvalue
	 * specialization.
	 */
	template<typename A>
	struct construct_term<A*&> : construct_term<A*> {};
	template<typename A>
	struct construct_term<symphas::ref<A>&> : construct_term<symphas::ref<A>> {};
	template<size_t Z, typename A>
	struct construct_term<Variable<Z, A>&> : construct_term<Variable<Z, A>> {};
	template<typename A>
	struct construct_term<NamedData<A>&> : construct_term<NamedData<A>> {};

	template<>
	struct construct_term<NamedData<OpIdentity>> : construct_term<OpIdentity> {};
	template<>
	struct construct_term<NamedData<OpNegIdentity>> : construct_term<OpNegIdentity> {};
	template<size_t N, size_t D>
	struct construct_term<NamedData<OpFractionLiteral<N, D>>> : construct_term<OpFractionLiteral<N, D>> {};
	template<size_t N, size_t D>
	struct construct_term<NamedData<OpNegFractionLiteral<N, D>>> : construct_term<OpNegFractionLiteral<N, D>> {};
	template<>
	struct construct_term<NamedData<OpVoid>> : construct_term<OpVoid> {};
}


namespace expr
{

	//! Delegate to create an OpTerm with the correct template. 
	/*!
	 * Delegate to create an OpTerm with the correct template. For
	 * more information on specific implementation, see symphas::internal::make_termlvariable.
	 *
	 * \param data The data from which to create an OpTerm.
	 *
	 * \tparam A The type of the data to make an OpTerm.
	 */
	template<typename A>
	auto make_term(A&& data)
	{
		return symphas::internal::construct_term<A>::get(std::forward<A>(data));
	}

	//! Delegate to create an OpTerm with the correct template. 
	/*!
	 * Delegate to create an OpTerm with the correct template. For
	 * more information on specific implementation, see symphas::internal::make_termlvariable.
	 * This overload will use the given index to create an OpTerm of a
	 * Variable.
	 *
	 * \param data The data from which to create an OpTerm.
	 *
	 * \tparam Z The index of the Variable to create.
	 * \tparam A The type of the data to make an OpTerm.
	 */
	template<size_t Z, typename A>
	auto make_term(A&& data)
	{
		return symphas::internal::construct_term<A>::template get<Z>(std::forward<A>(data));
	}

	//! Delegate to create an OpTerm with the correct template. 
	/*!
	 * Delegate to create an OpTerm with the correct template. For
	 * more information on specific implementation, see symphas::internal::make_termlvariable.
	 *
	 * \param value The coefficient of the OpTerm.
	 * \param data The data from which to create an OpTerm.
	 *
	 * \param A The type of the data to make an OpTerm.
	 */
	template<typename V, typename A>
	auto make_term(V&& value, A&& data)
	{
		return symphas::internal::construct_term<A>::template get(expr::make_literal(std::forward<V>(value)), std::forward<A>(data));
	}

	//! Delegate to create an OpTerm with the correct template. 
	/*!
	 * Delegate to create an OpTerm with the correct template. For
	 * more information on specific implementation, see symphas::internal::make_termlvariable.
	 * This overload will use the given index to create an OpTerm of a
	 * Variable.
	 *
	 * \param value The coefficient of the OpTerm.
	 * \param data The data from which to create an OpTerm.
	 *
	 * \tparam Z The index of the Variable to create.
	 * \tparam A The type of the data to make an OpTerm.
	 */
	template<size_t Z, typename V, typename A>
	auto make_term(V&& value, A&& data)
	{
		return symphas::internal::construct_term<A>::template get<Z>(expr::make_literal(std::forward<V>(value)), std::forward<A>(data));
	}

	//! Creates a symbolic variable with no data. 
	/*!
	 * For symbolic algebra manipulation, create a variable with no data, simply a symbolic type,
	 * in order to represent a symbol.
	 *
	 * \tparam Z The index of the Variable to create.
	 */
	template<size_t Z>
	auto make_term()
	{
		return OpTerms(OpIdentity{}, Term<Variable<Z>, Xk<1>>{});
	}

}


