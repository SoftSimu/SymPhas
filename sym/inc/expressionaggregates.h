
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


namespace symphas::internal
{

	template<typename G>
	struct make_data
	{
		template<typename GG = G, std::enable_if_t<std::is_default_constructible_v<GG>, int> = 0>
		G operator()()
		{
			return {};
		}

		template<typename GG = G, std::enable_if_t<!std::is_default_constructible_v<GG>, int> = 0>
		G operator()()
		{
			return { 0 };
		}
	};

	template<typename G, size_t D>
	struct make_data<GridData<G, D>>
	{
		GridData<G, D> operator()()
		{
			return {};
		}
	};

	template<typename G>
	struct make_data<DynamicVariable<G>>
	{
		DynamicVariable<G> operator()()
		{
			return {};
		}
	};

	template<typename G>
	struct make_data<symphas::ref<G>>
	{
		symphas::ref<G> operator()()
		{
			auto data = make_data<G>{}();
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
	using G::G;

#ifdef PRINTABLE_EQUATIONS


	//! Constructs the named instance using the given data.
	/*!
	 * Constructs the named instance using the given data, which is passed
	 * to the constructor of the underlying class.
	 *
	 * \param data The data that is named.
	 * \param name The name of the data.
	 */
	NamedData(G const& data, const char* name) : 
		G(data), name{ (name != nullptr && std::strlen(name) > 0) ? new char[std::strlen(name) + 1] : nullptr }
	{
		if (this->name != nullptr)
		{
			std::strcpy(this->name, name);
		}
	}

	NamedData(G const& data, std::string name) : NamedData(data, name.c_str()) {}

	NamedData(G data = symphas::internal::make_data<G>{}()) : NamedData(data, nullptr) {}

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
	NamedData(G const& data, OpExpression<E> const& e) : G(data), name{}
	{
		name = new char[(*static_cast<E const*>(&e)).print_length() + 1];
		(*static_cast<E const*>(&e)).print(name);
		//name = std::string(str);
		//delete[] str;
	}

	NamedData(NamedData<G> const& other) :
		NamedData(*static_cast<G const*>(&other), other.name) 
	{}

	NamedData(NamedData<G>&& other) : NamedData()
	{
		swap(*this, other);
	}

	NamedData<G> operator=(NamedData<G> other)
	{
		swap(*this, other);
		return *this;
	}

	friend void swap(NamedData<G>& first, NamedData<G>& second)
	{
		using std::swap;
		swap(static_cast<G&>(first), static_cast<G&>(second));
		swap(first.name, second.name);
	}



	~NamedData()
	{
		delete[] name;
	}

	char* name;

#else


	template<typename T>
	NamedData(G data, T&&) : G(data) {}
	NamedData(G data = symphas::internal::make_data<G>{}()) : NamedData(data, nullptr) {}

#endif

	NamedData(len_type) : NamedData() {}


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


	//! Constructs the named instance using the given data.
	/*!
	 * Constructs the named instance using the given data, which is passed
	 * to the constructor of the underlying class.
	 *
	 * \param data The data that is named.
	 * \param name The name of the data.
	 */
	NamedData(G* data, const char* name) :
		data{ data }, name{(name != nullptr && std::strlen(name) > 0) ? new char[std::strlen(name) + 1] : nullptr}
	{
		if (this->name != nullptr)
		{
			std::strcpy(this->name, name);
		}
	}

	NamedData(G* data, std::string name) : NamedData(data, name.c_str()) {}

	NamedData(G* data = nullptr) : NamedData(data, nullptr) {}

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
	NamedData(G* data, OpExpression<E> const& e) : data{ data }, name{}
	{
		name = new char[(*static_cast<E const*>(&e)).print_length() + 1];
		(*static_cast<E const*>(&e)).print(name);
	}

	NamedData(NamedData<G*> const& other) :
		NamedData(other.data, other.name)
	{}

	NamedData(NamedData<G*>&& other) : NamedData()
	{
		swap(*this, other);
	}

	NamedData<G*> operator=(NamedData<G*> other)
	{
		swap(*this, other);
		return *this;
	}

	friend void swap(NamedData<G*>& first, NamedData<G*>& second)
	{
		using std::swap;
		swap(first.data, second.data);
		swap(first.name, second.name);
	}

	~NamedData()
	{
		delete[] name;
	}

	char* name;


#else

	NamedData() : data{} {}

	template<typename T = char>
	NamedData(G* data, T&& = '\0') : data{data} {}

#endif


	operator G* ()
	{
		return data;
	}

	operator G* const () const
	{
		return data;
	}

protected:

	NamedData(len_type) : NamedData() {}

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

	NamedData(scalar_t = 0) : data{}, name{ "" } {}

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
	NamedData(scalar_t data, OpOperator<E> const& e) : data{ data }, name{ "" }
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
NamedData(G&, std::string) -> NamedData<symphas::ref<G>>;
template<typename G>
NamedData(G&&, std::string) -> NamedData<G>;
template<typename G>
NamedData(G*, std::string) -> NamedData<G*>;

template<typename G>
NamedData(G&, const char*) -> NamedData<symphas::ref<G>>;
template<typename G>
NamedData(G&&, const char*) -> NamedData<G>;
template<typename G>
NamedData(G*, const char*) -> NamedData<G*>;

template<typename G, typename E>
NamedData(G&, OpExpression<E>) -> NamedData<symphas::ref<G>>;
template<typename G, typename E>
NamedData(G&&, OpExpression<E>) -> NamedData<G>;
template<typename G, typename E>
NamedData(G*, OpExpression<E>) -> NamedData<G*>;

template<typename G, typename E>
NamedData(G&, OpOperator<E>) -> NamedData<symphas::ref<G>>;
template<typename G, typename E>
NamedData(G&&, OpOperator<E>) -> NamedData<G>;
template<typename G, typename E>
NamedData(G*, OpOperator<E>) -> NamedData<G*>;

#else

template<typename G, typename T>
NamedData(G&, T&&) -> NamedData<symphas::ref<G>>;
template<typename G, typename T>
NamedData(G&&, T&&) -> NamedData<G>;
template<typename G, typename T>
NamedData(G*, T&&) -> NamedData<G*>;

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

	constexpr Variable() : G() {}
	Variable(G const& data) : G(data) {}
	Variable(G&& data) noexcept : G(std::move(data)) {}
};

template<size_t Z, typename G>
struct Variable<Z, G*>
{
	constexpr Variable() : data{} {}
	Variable(G* data) : data(data) {}

	operator const G* () const 
	{
		return data;
	}

	operator G* ()
	{
		return data;
	}

	G* data;
};

//! Associates an identifier with a piece of data.
/*!
 * Associates a template constant number with a given
 * object; primarily used for parsing functions. Used in particular when
 * parsing or pruning equations. Acts as a data identifier.
 */
template<size_t Z, typename G>
struct Variable<Z, symphas::ref<G>> : symphas::ref<G>
{
	using parent_type = symphas::ref<G>;
	using parent_type::parent_type;

	Variable(symphas::ref<G> const& data) : parent_type(data) {}
	Variable(symphas::ref<G>&& data) noexcept : parent_type(std::move(data)) {}
	Variable(G&& data) noexcept : parent_type(data) {}
	constexpr Variable() : Variable(((G*)0)[0]) {}
};

template<typename G>
struct DynamicVariable : DynamicIndex
{
	constexpr DynamicVariable() : DynamicIndex(), data{ nullptr } {}
	DynamicVariable(DynamicIndex const& index, G* data) : DynamicIndex(index), data{ data } {}

	const auto& get() const
	{
		return data[*DynamicIndex::data];
	}

	auto& get()
	{
		return data[*DynamicIndex::data];
	}

	G* data;
};

template<typename G>
struct DynamicVariable<NamedData<G*>> : DynamicIndex
{
	constexpr DynamicVariable() : DynamicIndex(), data{ nullptr } {}
	DynamicVariable(DynamicIndex const& index, NamedData<G*> const& data) : DynamicIndex(index), data{ data } {}

	const auto& get() const
	{
		return data.data[*DynamicIndex::data];
	}

	auto& get()
	{
		return data.data[*DynamicIndex::data];
	}

	NamedData<G*> data;
};

template<size_t Z, typename G>
struct DynamicVariable<Variable<Z, G*>> : DynamicIndex
{
	constexpr DynamicVariable() : DynamicIndex(), data{ nullptr } {}
	DynamicVariable(DynamicIndex const& index, Variable<Z, G*> const& data) : DynamicIndex(index), data{ data } {}

	const auto& get() const
	{
		return data.data[*DynamicIndex::data];
	}

	auto& get()
	{
		return data.data[*DynamicIndex::data];
	}

	Variable<Z, G*> data;
};

template<typename G>
DynamicVariable(DynamicIndex, NamedData<G*>) -> DynamicVariable<NamedData<G*>>;
template<size_t Z, typename G>
DynamicVariable(DynamicIndex, Variable<Z, G*>) -> DynamicVariable<Variable<Z, G*>>;


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
	template<typename G>
	auto as_variable(G&& g)
	{
		return DynamicVariable<G>(std::forward<G>(g));
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
	template<typename G>
	auto as_variable(G& g)
	{
		return DynamicVariable<symphas::ref<G>>(g);
	}

	////! The given parameter is created as a Variable.
	///*!
	// * The given parameter is created and returned as a Variable.
	// *
	// * \param g The object which is created as a Variable.
	// */
	//template<Axis ax, typename G>
	//auto as_component(G& g)
	//{
	//	return VectorComponent<ax, G>(g);
	//}

	//! The given parameter is created as a Variable.
	/*!
	 * The given parameter is created and returned as a Variable.
	 *
	 * \param g The object which is created as a Variable.
	 */
	template<Axis ax, typename G>
	auto as_component(G&& g)
	{
		if constexpr (expr::is_expression<G>)
		{
			return to_axis<ax>(std::forward<G>(g));;
		}
		else
		{
			return VectorComponent<ax, G>(g);
		}
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
	auto as_component(DynamicVariable<G> g)
	{
		return DynamicVariable<VectorComponent<ax, G>>(*static_cast<G*>(&g));
	}

}


#ifdef PRINTABLE_EQUATIONS
template<typename E>
auto symphas::internal::set_var_string(OpExpression<E> const& var)
{
	char* buffer = new char[static_cast<E const*>(&var)->print_length() + 1];
	static_cast<E const*>(&var)->print(buffer);
	return buffer;
}

template<typename E>
auto symphas::internal::set_var_string(OpOperator<E> const& var)
{
	char* buffer = new char[static_cast<E const*>(&var)->print_length() + 1];
	static_cast<E const*>(&var)->print(buffer);
	return buffer;
}

template<typename G>
auto symphas::internal::set_var_string(G const& var)
{
	const char* name = expr::get_op_name(var);
	char* buffer = new char[std::strlen(name) + 1];
	std::strcpy(buffer, name);
	return buffer;
}
#endif


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

	//! Combines this data with an expression.
	template<typename E>
	auto operator=(OpOperator<E> const& expr) const
	{
		return std::pair<G, E>(*this, *static_cast<E const*>(&expr));
	}
};


// *************************************************************************************


template<typename G, expr::exp_key_t X>
struct Term : G
{
	using G::G;

	constexpr Term() : G() {}
	Term(G const& data) : G(data) {}
	Term(G&& data) noexcept : G(std::move(data)) {}

	inline auto eval(iter_type n) const
	{
		if constexpr (expr::is_symbol<G>)
		{
			return expr::BaseData<G>::get(G{});
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
				auto result = pow<size_t(expr::_Xk_t<X>::N)>(expr::BaseData<G>::get(data(), n));
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
			return sprintf(out, "%s" SYEX_POW_SEP_A "%u" SYEX_POW_DIV_SEP "%u" SYEX_POW_SEP_B, expr::get_op_name(data()), expr::_Xk_t<X>::N, expr::_Xk_t<X>::D);
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

	Term<G, X> const& operator*(OpIdentity) const
	{
		return *this;
	}

	Term<G, X>& operator*(OpIdentity)
	{
		return *this;
	}

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
struct Term<T*, X> : Term<symphas::pointer_type<T>, X>
{
	using parent_type = Term<symphas::pointer_type<T>, X>;
	using parent_type::parent_type;

	Term(T* data) : parent_type(symphas::pointer_type<T>(data)) {}

	Term<T*, X> const& operator*(OpIdentity) const
	{
		return *this;
	}

	Term<T*, X>& operator*(OpIdentity)
	{
		return *this;
	}
};

template<typename... Es>
struct OpTermsList;

template<>
struct OpTermsList<>
{
	inline auto _eval(iter_type) const
	{
		return OpIdentity{};
	}

#ifdef PRINTABLE_EQUATIONS

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
		using type = OpTermsList<>;
	};

	template<typename T0, typename... Ts>
	struct Nth_type_of_terms<0, OpTerms<T0, Ts...>>
	{
		using type = OpTermsList<T0, Ts...>;
	};


	template<typename cast_type>
	struct cast_term;

	template<typename... T0s>
	struct cast_term<OpTerms<T0s...>>
	{
		using cast_type = OpTerms<T0s...>;

		template<typename... Ts>
		static cast_type cast(OpTerms<Ts...> const& terms)
		{
			return *static_cast<OpTermsList<T0s...> const*>(&terms);
		}

		template<typename... Ts>
		static cast_type cast(OpTerms<Ts...>& terms)
		{
			return *static_cast<OpTermsList<T0s...>*>(&terms);
		}
	};

	template<typename... T0s>
	struct cast_term<OpTermsList<T0s...>>
	{
		using cast_type = OpTermsList<T0s...>;

		template<typename... Ts>
		static cast_type const& cast(OpTerms<Ts...> const& terms)
		{
			return *static_cast<OpTermsList<T0s...> const*>(&terms);
		}

		template<typename... Ts>
		static cast_type& cast(OpTerms<Ts...>& terms)
		{
			return *static_cast<OpTermsList<T0s...>*>(&terms);
		}

		template<typename... Ts>
		static cast_type const& cast(OpTermsList<Ts...> const& terms)
		{
			return *static_cast<OpTermsList<T0s...> const*>(&terms);
		}

		template<typename... Ts>
		static cast_type& cast(OpTermsList<Ts...>& terms)
		{
			return *static_cast<OpTermsList<T0s...>*>(&terms);
		}
	};

	template<typename... Gs, expr::exp_key_t... Xs>
	struct cast_term<symphas::lib::types_list<Term<Gs, Xs>...>>
	{
		template<typename V, typename... G0s, expr::exp_key_t... X0s>
		static OpTerms<Term<Gs, Xs>...> cast(OpTerms<V, Term<G0s, X0s>...> const& terms)
		{
			return *static_cast<OpTermsList<Term<Gs, Xs>...> const*>(&terms);
		}

		template<typename V, typename... G0s, expr::exp_key_t... X0s>
		static OpTerms<Term<Gs, Xs>...> cast(OpTerms<V, Term<G0s, X0s>...>& terms)
		{
			return *static_cast<OpTermsList<Term<Gs, Xs>...>*>(&terms);
		}

		template<typename V, typename... G0s, expr::exp_key_t... X0s>
		static OpTermsList<Term<Gs, Xs>...> const& cast(OpTermsList<V, Term<G0s, X0s>...> const& terms)
		{
			return *static_cast<OpTermsList<Term<Gs, Xs>...> const*>(&terms);
		}

		template<typename V, typename... G0s, expr::exp_key_t... X0s>
		static OpTermsList<Term<Gs, Xs>...>& cast(OpTermsList<V, Term<G0s, X0s>...>& terms)
		{
			return *static_cast<OpTermsList<Term<Gs, Xs>...>*>(&terms);
		}
	};

#ifdef PRINTABLE_EQUATIONS

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

	template<typename E>
	size_t print_length_one_term(OpOperator<E> const& e)
	{
		return expr::coeff_print_length(*static_cast<E const*>(&e));
	}

	template<typename G, expr::exp_key_t X>
	size_t print_length_one_term(Term<G, X> const& term)
	{
		return term.print_length();
	}

#endif



	inline auto to_term(OpTermsList<>)
	{
		return OpIdentity{};
	}

	template<typename T>
	auto to_term(T const& e)
	{
		return OpTerms(e);
	}

	template<typename T>
	auto to_term(T& e)
	{
		return OpTerms(e);
	}
}


template<typename G>
Term(G) -> Term<G, expr::Xk<1>>;
template<typename G>
Term(G*) -> Term<symphas::ref<G*>, expr::Xk<1>>;


namespace expr
{
	template<typename E>
	constexpr bool is_term = false;

	template<typename G, exp_key_t X>
	constexpr bool is_term<Term<G, X>> = true;

	template<size_t N, typename V, typename... Gs, exp_key_t... Xs>
	const auto& get(OpTerms<V, Term<Gs, Xs>...> const& e)
	{
		using Nth_type = typename symphas::internal::Nth_type_of_terms<N, OpTerms<V, Term<Gs, Xs>...>>::type;
		return symphas::internal::cast_term<Nth_type>::cast(e).term;
	}

	template<size_t N, typename V, typename... Gs, exp_key_t... Xs>
	auto& get(OpTerms<V, Term<Gs, Xs>...>& e)
	{
		using Nth_type = typename symphas::internal::Nth_type_of_terms<N, OpTerms<V, Term<Gs, Xs>...>>::type;
		return symphas::internal::cast_term<Nth_type>::cast(e).term;
	}

	template<size_t N, typename V, typename... Gs, exp_key_t... Xs>
	const auto& get(OpTermsList<V, Term<Gs, Xs>...> const& e)
	{
		using Nth_type = typename symphas::internal::Nth_type_of_terms<N, OpTerms<V, Term<Gs, Xs>...>>::type;
		return symphas::internal::cast_term<Nth_type>::cast(e).term;
	}

	template<size_t N, typename V, typename... Gs, exp_key_t... Xs>
	auto& get(OpTermsList<V, Term<Gs, Xs>...>& e)
	{
		using Nth_type = typename symphas::internal::Nth_type_of_terms<N, OpTerms<V, Term<Gs, Xs>...>>::type;
		return symphas::internal::cast_term<Nth_type>::cast(e).term;
	}

	//template<typename G0, exp_key_t X0, typename G1, exp_key_t X1, typename... Gs, exp_key_t... Xs>
	//const auto& terms_after_first(OpTerms<Term<G0, X0>, Term<G1, X1>, Term<Gs, Xs>...> const& e)
	//{
	//    using cast_to_t = OpTerms<Term<G0, X0>, Term<G1, X1>, Term<Gs, Xs>...>;
	//	//return symphas::internal::cast_term<Nth_type>::cast(e).term;
	//    return e;
	//    //*static_cast<OpTerms<Term<G0, X0>, Term<G1, X1>, Term<Gs, Xs>...> const*>(&e);
	//}

	//template<typename G0, exp_key_t X0, typename G1, exp_key_t X1, typename... Gs, exp_key_t... Xs>
	//auto& terms_after_first(OpTerms<Term<G0, X0>, Term<G1, X1>, Term<Gs, Xs>...>& e)
	//{
	//    return e;
	//	//return *static_cast<OpTerms<Term<G0, X0>, Term<G1, X1>, Term<Gs, Xs>...>*>(&e);
	//}

	template<typename G0, exp_key_t X0>
	const auto& terms_after_first(OpTerms<Term<G0, X0>> const& e)
	{
		return e;
	}

	template<typename G0, exp_key_t X0>
	auto& terms_after_first(OpTerms<Term<G0, X0>>& e)
	{
		return e;
	}

	template<typename V, typename G0, exp_key_t X0, typename... Gs, exp_key_t... Xs>
	decltype(auto) terms_after_first(OpTerms<V, Term<G0, X0>, Term<Gs, Xs>...> const& e)
	{
		if constexpr (is_term<V>)
		{
			return e;
		}
		else
		{
			using cast_to_t = OpTerms<Term<G0, X0>, Term<Gs, Xs>...>;
			return symphas::internal::cast_term<cast_to_t>::cast(e);
		}
	}

	template<typename V, typename G0, exp_key_t X0, typename... Gs, exp_key_t... Xs>
	decltype(auto) terms_after_first(OpTerms<V, Term<G0, X0>, Term<Gs, Xs>...>& e)
	{
		if constexpr (is_term<V>)
		{
			return e;
		}
		else
		{
			using cast_to_t = OpTerms<Term<G0, X0>, Term<Gs, Xs>...>;
			return symphas::internal::cast_term<cast_to_t>::cast(e);
		}
	}

	template<typename V>
	decltype(auto) terms_after_first(OpTerms<V> const& e)
	{
		using cast_to_t = OpTerms<>;
		return symphas::internal::cast_term<cast_to_t>::cast(e);
	}

	template<typename V>
	decltype(auto) terms_after_first(OpTerms<V>& e)
	{
		using cast_to_t = OpTerms<>;
		return symphas::internal::cast_term<cast_to_t>::cast(e);
	}

	template<size_t N, typename V, typename... Gs, exp_key_t... Xs>
	decltype(auto) terms_after_n(OpTermsList<V, Term<Gs, Xs>...> const& e)
	{
		using ts = symphas::lib::types_after_at_index<N + 1, V, Term<Gs, Xs>...>;
		if constexpr (N + 1 == sizeof...(Gs))
		{
			return symphas::internal::cast_term<ts>::cast(e);
		}
		else if constexpr (N + 1 < sizeof...(Gs) + 1)
		{
			return symphas::internal::cast_term<ts>::cast(e);
		}
		else
		{
			return OpTermsList<>{};
		}
	}

	template<size_t N, typename V, typename... Gs, exp_key_t... Xs>
	decltype(auto) terms_after_n(OpTermsList<V, Term<Gs, Xs>...>& e)
	{
		using ts = symphas::lib::types_after_at_index<N + 1, V, Term<Gs, Xs>...>;
		if constexpr (N + 1 == sizeof...(Gs))
		{
			return symphas::internal::cast_term<ts>::cast(e);
		}
		else if constexpr (N + 1 < sizeof...(Gs) + 1)
		{
			return symphas::internal::cast_term<ts>::cast(e);
		}
		else
		{
			return OpTermsList<>{};
		}
	}

	template<typename... Ts>
	decltype(auto) terms_after_first(OpTermsList<Ts...> const& e)
	{
		if constexpr (sizeof...(Ts) > 1)
		{
			return terms_after_n<0>(e);
		}
		else
		{
			return OpIdentity{};
		}
	}

	template<typename... Ts>
	decltype(auto) terms_after_first(OpTermsList<Ts...>& e)
	{
		if constexpr (sizeof...(Ts) > 1)
		{
			return terms_after_n<0>(e);
		}
		else
		{
			return OpIdentity{};
		}
	}


	template<size_t N, typename V, typename... Gs, exp_key_t... Xs>
	decltype(auto) terms_after_n(OpTerms<V, Term<Gs, Xs>...> const& e)
	{
		return symphas::internal::to_term(terms_after_n<N>(*static_cast<OpTermsList<V, Term<Gs, Xs>...> const*>(&e)));
	}

	template<size_t N, typename V, typename... Gs, exp_key_t... Xs>
	decltype(auto) terms_after_n(OpTerms<V, Term<Gs, Xs>...>& e)
	{
		return symphas::internal::to_term(terms_after_n<N>(*static_cast<OpTermsList<V, Term<Gs, Xs>...> const*>(&e)));
	}

	template<typename V, typename G, exp_key_t X>
	const auto& data(OpTerms<V, Term<G, X>> const& e)
	{
		if constexpr (is_term<V>)
		{
			return expr::get<0>(e).data();
		}
		else
		{
			return expr::get<1>(e).data();
		}
	}

	template<typename V, typename G, exp_key_t X>
	auto& data(OpTerms<V, Term<G, X>>& e)
	{
		if constexpr (is_term<V>)
		{
			return expr::get<0>(e).data();
		}
		else
		{
			return expr::get<1>(e).data();
		}
	}

}

namespace std
{
	template<typename V, typename G>
	struct tuple_size<OpTerm<V, G>> : integral_constant<size_t, 2> {};

	template<size_t I, typename V, typename G>
	struct tuple_element<I, OpTerm<V, G>> : tuple_element<I, tuple<V, G>> {};
}



template<size_t I, typename V, typename G>
std::tuple_element_t<I, OpTerm<V, G>>& get(OpTerm<V, G>& e)
{
	return expr::get<I>(e);
}

template<size_t I, typename V, typename G>
const std::tuple_element_t<I, OpTerm<V, G>>& get(OpTerm<V, G> const& e)
{
	return expr::get<I>(e);
}



ALLOW_COMBINATION((Axis ax, typename G), (VectorComponent<ax, G>))

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

template<typename V, typename... Gs, expr::exp_key_t... Xs>
struct OpTermsList<V, Term<Gs, Xs>...> : OpTermsList<Term<Gs, Xs>...>
{
	using parent_type = OpTermsList<Term<Gs, Xs>...>;

	OpTermsList() : parent_type(), term{ } {}
	OpTermsList(V const& e, Term<Gs, Xs> const&... terms) : parent_type(terms...), term{ e } {}
	OpTermsList(V const& e, OpTerms<Term<Gs, Xs>...> const& rest) :
		parent_type(*static_cast<OpTermsList<Term<Gs, Xs>...> const*>(&rest)), term{ e } {}
	OpTermsList(OpTerms<V, Term<Gs, Xs>...> const& list) :
		parent_type(*static_cast<OpTermsList<Term<Gs, Xs>...> const*>(&list)), term{ list.term } {}
	//OpTermsList(OpTermsList<V, Term<Gs, Xs>...> const& list) : 
	//    parent_type(*static_cast<OpTermsList<Term<Gs, Xs>...> const*>(&list)), term{ list.term } {}

	template<typename G0, expr::exp_key_t X0>
	auto _eval(Term<G0, X0> const& term0, iter_type n) const
	{
		return term0.eval(n) * parent_type::_eval(n);
	}

	template<typename V0>
	auto _eval(V0 const& value, iter_type n) const
	{
		return expr::eval(value) * parent_type::_eval(n);
	}

	inline auto _eval(iter_type n) const
	{
		return _eval(term, n);
	}


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
struct OpTerms<V, Term<Gs, Xs>...> : OpExpression<OpTerms<V, Term<Gs, Xs>...>>, OpTermsList<V, Term<Gs, Xs>...>
{
	using parent_type = OpTermsList<V, Term<Gs, Xs>...>;
	using parent_type::term;

#ifdef PRINTABLE_EQUATIONS
	using parent_type::print;
	using parent_type::print_length;
#endif

	constexpr OpTerms() : parent_type() {}

	OpTerms(V const& term0, Term<Gs, Xs> const&... terms)
		: parent_type(term0, terms...) {}
	OpTerms(V const& term0, OpTerms<Term<Gs, Xs>...> const& terms)
		: parent_type(term0, terms) {}
	//OpTerms(OpTerms<V, Term<Gs, Xs>...> const& terms) 
	//	: parent_type(terms) {}
	OpTerms(OpTermsList<V, Term<Gs, Xs>...> const& terms)
		: parent_type(terms) {}

	template<typename E, size_t N = sizeof...(Gs), typename std::enable_if_t<(N == 1), int> = 0>
	void operator=(OpExpression<E> const& e)
	{
		expr::result(
			*static_cast<E const*>(&e),
			static_cast<OpTermsList<Term<Gs, Xs>...> const*>(this)->term);
	}

	auto operator-() const;

	auto eval(iter_type n = 0) const
	{
		return parent_type::_eval(n);
	}

};

template<>
struct OpTerms<> : OpTermsList<> {};

template<typename V, typename... Gs, expr::exp_key_t... Xs>
OpTerms(V, Term<Gs, Xs>...) -> OpTerms<V, Term<Gs, Xs>...>;

template<typename V, typename... Gs, expr::exp_key_t... Xs>
OpTerms(V, OpTerms<Term<Gs, Xs>...>) -> OpTerms<V, Term<Gs, Xs>...>;

template<typename V, typename... Gs, expr::exp_key_t... Xs, typename G0, expr::exp_key_t X0>
OpTerms(OpTerms<V, Term<Gs, Xs>...>, Term<G0, X0>) -> OpTerms<V, Term<G0, X0>, Term<Gs, Xs>...>;

template<typename... Ts>
OpTerms(OpTermsList<Ts...>) -> OpTerms<Ts...>;

OpTerms() -> OpTerms<>;


template<typename V, typename... Gs, expr::exp_key_t... Xs>
auto OpTerms<V, Term<Gs, Xs>...>::operator-() const
{
	return ::OpTerms(-term, expr::terms_after_n<0>(*this));
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

	template<Axis ax, typename G>
	auto to_vector_component(DynamicVariable<G> const& data)
	{
		return expr::as_variable(to_vector_component<ax>(*static_cast<G const*>(&data)));;
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
	size_t R = expr::eval_type<OpTerms<Term<Gs, Xs>...>>::rank,
	typename std::enable_if_t<(expr::is_coeff<coeff_t> && !expr::is_tensor<V>
		&& ((R > 1) ? R != expr::eval_type<coeff_t>::template rank_<1> : true)), int> = 0>
	auto operator*(coeff_t const& value, OpTerms<V, Term<Gs, Xs>...> const& b)
{
	return OpTerms(value * expr::coeff(b), expr::terms_after_first(b));
}

template<typename coeff_t, typename tensor_t, typename... Gs, expr::exp_key_t... Xs,
	typename std::enable_if_t<(expr::is_coeff<coeff_t>&& expr::is_tensor<tensor_t>), int> = 0>
auto operator*(coeff_t const& value, OpTerms<tensor_t, Term<Gs, Xs>...> const& b)
{
	return (value * expr::coeff(b)) * OpTerms(OpIdentity{}, expr::terms_after_first(b));
}

template<typename T, typename V, typename... Gs, expr::exp_key_t... Xs,
	size_t N0, size_t N1, size_t N,
	size_t D = expr::eval_type<OpTerms<Term<Gs, Xs>...>>::rank,
	typename std::enable_if_t<(!expr::is_tensor<V>&& expr::eval_type<OpTerms<Term<Gs, Xs>...>>::rank > 0), int> = 0>
	auto operator*(OpTensor<T, N0, N1, N, D> const& tensor, OpTerms<V, Term<Gs, Xs>...> const& b)
{
	constexpr Axis ax = (N1 == 0) ? Axis::X : (N1 == 1) ? Axis::Y : (N1 == 2) ? Axis::Z : Axis::NONE;
	if constexpr (N == 1)
	{
		auto coeff = symphas::internal::tensor_cast::cast(tensor) * expr::coeff(b);
		return coeff * symphas::internal::terms_make_component<ax>(
			expr::terms_after_first(b), std::make_index_sequence<sizeof...(Gs)>{});
	}
	else
	{
		auto coeff = expr::make_tensor<N0, N>(symphas::internal::tensor_cast::cast(tensor))* expr::coeff(b);
		return coeff * symphas::internal::terms_make_component<ax>(
			expr::terms_after_first(b), std::make_index_sequence<sizeof...(Gs)>{});
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

template<typename E, typename G0, expr::exp_key_t X0, typename... Gs, expr::exp_key_t... Xs>
auto operator*(OpTerms<Term<G0, X0>, Term<Gs, Xs>...> const& a, OpOperator<E> const& b)
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

	template<size_t P, typename... Gs, typename G, expr::exp_key_t X, expr::exp_key_t... Xs>
	auto place_at(OpTerms<Term<Gs, Xs>...> const& terms, Term<G, X> const& term0)
	{
		return place_at(terms, term0, std::make_index_sequence<P>{}, std::make_index_sequence<sizeof...(Gs) - P>{});
	}

    template<typename G0, expr::exp_key_t X0, typename... Gs, expr::exp_key_t... Xs>
    auto combine_terms(Term<G0, X0> const& term0, Term<Gs, Xs> const& ...terms)
    {
        return OpTerms(term0, terms...);
    }
    
    inline auto combine_terms()
    {
        return OpTerms<>();
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

	template<size_t P0, size_t... Ps, typename... G1s, expr::exp_key_t... X1s, typename... G2s, expr::exp_key_t... X2s, typename G2, expr::exp_key_t X2, typename... Ts>
	auto combine_terms(std::index_sequence<P0, Ps...>, OpTerms<Term<G1s, X1s>...> const& terms_num, OpTerms<Term<G2s, X2s>...> const& terms_den, Term<G2, X2> const& term0, Ts&&... rest)
	{
		if constexpr (X2 == 0)
		{
			return combine_terms(std::index_sequence<Ps...>{}, terms_num, terms_den, std::forward<Ts>(rest)...);
		}
		else
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
	}

	template<typename... G1s, expr::exp_key_t... X1s, typename... G2s, expr::exp_key_t... X2s, size_t... Is, size_t... Js, size_t... Ms, size_t... Ns>
	auto combine_terms(OpTerms<Term<G1s, X1s>...> const& a, OpTerms<Term<G2s, X2s>...> const& b,
		std::index_sequence<Is...>, std::index_sequence<Js...>,
		std::index_sequence<Ms...>, std::index_sequence<Ns...>)
	{
		if constexpr ((expr::_Xk_t<X2s>::sign && ...))
		{
			return combine_terms(std::index_sequence<Ms...>{}, combine_terms(expr::get<Is>(a)...), combine_terms(~expr::get<Js>(b)...), (expr::get<Ms>(a) * expr::get<Ns>(b))...);
		}
		else
		{
			return combine_terms(std::index_sequence<Ms...>{}, combine_terms(expr::get<Is>(a)..., expr::get<Js>(b)...), OpTerms<>{}, (expr::get<Ms>(a) * expr::get<Ns>(b))...);
		}
	}

	template<typename G, typename T1>
	constexpr bool commutes_through = false;

	template<typename G, typename... Gs>
	constexpr bool commutes_through<G, symphas::lib::types_list<Gs...>> = (true && ... && expr::is_commutable<G, Gs>);


	template<typename G0, typename G1, typename T>
	constexpr int commutes_to = -1;

	//! Indicates whether G0 can reach a G1 in Gs list through commuting.
	/*!
	 * Returns true if the `G0` type can reach the `G1` type by commuting through
	 * the `Gs` list. Two conditions must effectively be met:
	 * 1. `G1` must exist in `Gs`; if not then this will evaluate to false;
	 * 2. `G0` must commute with everything in `Gs`, up to the occurence of a `G1`.
	 * If both conditions are met, then this statement will evaluate true.
	 */
	template<typename G0, typename G1, typename... Gs>
	constexpr int commutes_to<G0, G1, types_list<Gs...>> =
		(commutes_through<G0, types_before_index<fixed_max<0, index_of_type<G1, Gs...>>, Gs...>>)
		? index_of_type<G1, Gs...>
		: -1;

	template<typename G1, typename T1, typename T2>
	constexpr int satisfies_identity = -1;

	template<typename G1, typename... Gs, typename G2, typename... G2s>
	constexpr int satisfies_identity<G1, types_list<Gs...>, types_list<G2, G2s...>> = 
		(expr::has_identity<G1, G2>) 
		? commutes_to<G1, G2, types_list<Gs..., G2s...>> - int(sizeof...(Gs))
		: satisfies_identity<G1, types_list<Gs...>, types_list<G2s...>>; // will combine with a G2 if there is an identity

	template<typename G1, typename... Gs>
	constexpr int satisfies_identity<G1, types_list<Gs...>, types_list<>> = -1; // will combine with a G2 if there is an identity

	template<typename G1, typename T1, typename T2>
	constexpr int satisfies_combination = -1;

	template<typename G1, typename... Gs, typename... G2s>
	constexpr int satisfies_combination<G1, types_list<Gs...>, types_list<G2s...>> =
		(expr::is_combinable<G1>)			// checks if G1 combines with another G1 in G2s list
		? commutes_to<G1, G1, types_list<Gs..., G2s...>> - int(sizeof...(Gs))
		: satisfies_identity<G1, types_list<Gs...>, types_list<G2s...>>;


	template<typename... G1s, expr::exp_key_t... X1s, typename... G2s, expr::exp_key_t... X2s, size_t... Is, size_t... Js>
	auto combine_terms(OpTerms<Term<G1s, X1s>...> const& a, OpTerms<Term<G2s, X2s>...> const& b, std::index_sequence<Is...>, std::index_sequence<Js...>)
	{
		using seq_a_match_b = std::integer_sequence<int, (
			(satisfies_combination<G1s,
				types_after_at_index<Is + 1, G1s...>,
				types_list<G2s...>>)
			)...>;

		using seq_a_pick = symphas::lib::seq_join_t<std::index_sequence<>,
			std::conditional_t<
				(symphas::lib::seq_index_value<Is, seq_a_match_b>::value >= 0),
				std::index_sequence<Is>,
				std::index_sequence<>>
			...>;
		using seq_b_pick = symphas::lib::seq_join_t<std::index_sequence<>,
			std::conditional_t<
				(symphas::lib::seq_index_value<Is, seq_a_match_b>::value >= 0),
				std::index_sequence<(size_t)symphas::lib::seq_index_value<Is, seq_a_match_b>::value>,
				std::index_sequence<>>
			...>;

		using seq_a_pick_ = symphas::lib::filter_seq_t<std::index_sequence<Is...>, seq_a_pick>;
		using seq_b_pick_ = symphas::lib::filter_seq_t<std::index_sequence<Js...>, seq_b_pick>;

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

		template<typename S>
		static auto get(S value, Variable<Z, A> const& data)
		{
			return construct_term<A>::template get<Z>(value, data);
		}

		static auto get(Variable<Z, A> const& data)
		{
			return construct_term<A>::template get<Z>(data);
		}

		template<size_t Y, typename S>
		static auto get(S value, Variable<Z, A> const& data)
		{
			return construct_term<A>::template get<Y>(value, data);
		}

		template<size_t Y>
		static auto get(Variable<Z, A> const& data)
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

	template<typename G, expr::exp_key_t X>
	struct construct_term<Term<G, X>>
	{
		template<typename S>
		static auto construct(S value, Term<G, X> data)
		{
			return OpTerms<S, Term<G, X>>(value, data);
		}

		template<typename S>
		static auto get(S value, Term<G, X> data)
		{
			return construct(value, data);
		}

		static auto get(Term<G, X> data)
		{
			return construct(OpIdentity{}, data);
		}

		template<size_t Z, typename S>
		static auto get(S value, Term<G, X> data)
		{
			return construct(value, expr::as_variable<Z>(data));
		}

		template<size_t Z>
		static auto get(Term<G, X> data)
		{
			return construct(OpIdentity{}, expr::as_variable<Z>(data));
		}
	};


	template<typename... Vs>
	struct construct_term<OpTerms<Vs...>>
	{
		template<typename S, typename... Gs, expr::exp_key_t... Xs>
		static auto construct(S value, S, OpTerms<Term<Gs, Xs>...> data)
		{
			return OpTerms<S, Term<Gs, Xs>...>(value, data);
		}

		template<typename G0, expr::exp_key_t X0, typename... Gs, expr::exp_key_t... Xs>
		static auto construct(OpIdentity, Term<G0, X0> term, OpTerms<Term<Gs, Xs>...> data)
		{
            auto terms = OpTermsList<Term<G0, X0>, Term<Gs, Xs>...>(term, data);
			return OpTerms<OpIdentity, Term<G0, X0>, Term<Gs, Xs>...>(OpIdentity{}, terms);
		}
        
		template<typename G0, expr::exp_key_t X0>
		static auto construct(OpIdentity, Term<G0, X0> term, OpIdentity)
		{
            auto terms = OpTermsList<Term<G0, X0>>(term);
			return OpTerms<OpIdentity, Term<G0, X0>>(OpIdentity{}, terms);
		}

		template<typename S>
		static auto get(S value, OpTerms<Vs...> data)
		{
			return value * construct(expr::coeff(data), expr::get<0>(data), expr::terms_after_n<0>(data));
		}

		static auto get(OpTerms<Vs...> data)
		{
			return construct(expr::coeff(data), expr::get<0>(data), expr::terms_after_n<0>(data));
		}

		template<size_t Z, typename S>
		static auto get(S value, OpTerms<Vs...> data)
		{
			return value * construct(expr::coeff(data), expr::get<0>(data), expr::terms_after_n<0>(data));
		}

		template<size_t Z>
		static auto get(OpTerms<Vs...> data)
		{
			return construct(expr::coeff(data), expr::get<0>(data), expr::terms_after_n<0>(data));
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
	struct construct_term<const A> : construct_term<A> {};
	template<typename A>
	struct construct_term<symphas::ref<A>&> : construct_term<symphas::ref<A>> {};
	template<size_t Z, typename A>
	struct construct_term<Variable<Z, A>&> : construct_term<Variable<Z, A>> {};
	template<typename A>
	struct construct_term<DynamicVariable<A>&> : construct_term<DynamicVariable<A>> {};
	template<typename A>
	struct construct_term<NamedData<A>&> : construct_term<NamedData<A>> {};
	template<typename... Ts>
	struct construct_term<OpTerms<Ts...>&> : construct_term<OpTerms<Ts...>> {};
	template<>
	struct construct_term<expr::symbols::Symbol&> : construct_term<expr::symbols::Symbol> {};

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
	template<typename... Ts>
	struct construct_term<SymbolicCase<Ts...>&> : construct_term<SymbolicCase<Ts...>> {};
	template<typename T, typename I>
	struct construct_term<OpCoeff<T, I>&> : construct_term<OpCoeff<T, I>> {};
	template<typename G, expr::exp_key_t X>
	struct construct_term<Term<G, X>&> : construct_term<Term<G, X>> {};
	template<typename G, expr::exp_key_t X>
	struct construct_term<Term<G, X> const&> : construct_term<Term<G, X>> {};
	template<typename... Ts>
	struct construct_term<OpTermsList<Ts...>> : construct_term<OpTerms<Ts...>> {};
	template<typename... Ts>
	struct construct_term<OpTermsList<Ts...>&> : construct_term<OpTermsList<Ts...>> {};
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

	
	template<typename V, typename A>
	auto make_term_dynamic(DynamicIndex const& index, V&& value, A&& data)
	{
		return expr::make_term(std::forward<V>(value), DynamicVariable(index, std::forward<A>(data)));
	}

	template<typename A>
	auto make_term_dynamic(DynamicIndex const& index, A&& data)
	{
		return expr::make_term_dynamic(index, OpIdentity{}, std::forward<A>(data));
	}

}


namespace symphas::internal
{

#ifdef PRINTABLE_EQUATIONS

	template<typename T>
	void update_dynamic_op_name(T* data, const char* const* names, len_type len)
	{
		for (iter_type i = 0; i < len; ++i)
		{
			char*& ptr = expr::get_op_name_reference(data);
			if (std::strcmp(names[i], ptr) != 0)
			{
				char* cpy = new char[std::strlen(names[i]) + 1];
				std::strcpy(cpy, names[i]);
				std::swap(cpy, ptr);
				delete[] cpy;
			}
		}
	}

	template<typename T>
	char** format_dynamic_op_name(T* data, const char* format, DynamicIndex index)
	{
		char** names = new char* [index.end() - index.start() + 1];
		for (iter_type i = index.start(); i <= index.end(); ++i)
		{
			names[i] = new char[std::strlen(format) - 1 + symphas::lib::num_digits(i)];
			sprintf(names[i], format, i);
		}
		update_dynamic_op_name(data, names, index.end() - index.start() + 1);
		return names;
	}

#endif

}


namespace expr
{

#ifdef PRINTABLE_EQUATIONS

	//! Gets the string name associated with the data.
	template<typename G>
	const char* get_op_name(DynamicVariable<NamedData<G*>> const& a)
	{
		char* s = std::strstr(a.data.name, "%d");
		if (s != NULL)
		{
			static len_type start = a.start();
			static len_type end = a.end();
			static char** str = symphas::internal::format_dynamic_op_name(a.data.data, a.data.name, a);

			if (a.index() > end || a.index() < start)
			{
				for (iter_type i = start; i <= end; ++i)
				{
					delete[] str[i];
				}
				delete[] str;

				start = std::min(start, a.index());
				end = std::max(end, a.index());
				str = symphas::internal::format_dynamic_op_name(a.data.data, a.data.name, DynamicIndex(0, start, end));
			}

			return str[a.index()];
		}
		else
		{
			return a.data.name;
		}
	}

#endif

}




