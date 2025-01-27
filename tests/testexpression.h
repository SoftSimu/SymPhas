#pragma once

#include "conf.h"
#include "solver.h"
#include "stencilincludes.h"

#include "expressiontypeincludes.h"

void testexpressioneval();
void testexpressionspeed();
void testexpressionmodel();
void testexpressionmodelspeed();
void testexpressionoperator();
void testexpressionoperatorspeed();
void testexpressionstate();
void testexpressioncollect();
void testexpressionconvolution();
void testexpressionvariable();
void testexpressionswap();
void testexpressionsort();
void testexpressiondivision();
void testexpressionexponential();
void testexpressionevaluation();
void testexpressionkgrid();
void testderivativefactor();
void testexpressiondistribution();
void test12varexpression();
void testccliketypes();
void testcoefftype();


inline void testexpression()
{
	printf("\n- test testexpressioneval ----------------------------------\n");
	testexpressioneval();
	printf("\n- test testexpressionmodel ---------------------------------\n");
	/*testexpressionmodel();
	printf("\n- test testexpressionoperator ------------------------------\n");*/
	testexpressionoperator();
	printf("\n- test testexpressionstate ---------------------------------\n");
	testexpressionstate();
	printf("\n- test testexpressioncollect -------------------------------\n");
	testexpressioncollect();
	printf("\n- test testexpressionconvolution ---------------------------\n");
	testexpressionconvolution();
	printf("\n- test testexpressionvariable ------------------------------\n");
	testexpressionvariable();
	printf("\n- test testexpressionswap ----------------------------------\n");
	testexpressionswap();
	printf("\n- test testexpressionsort ----------------------------------\n");
	testexpressionsort();
	printf("\n- test testexpressiondivision ------------------------------\n");
	testexpressiondivision();
	printf("\n- test testexpressionexponential ---------------------------\n");
	testexpressionexponential();
	printf("\n- test testexpressionevaluation ----------------------------\n");
	testexpressionevaluation();
	printf("\n- test testderivativefactor --------------------------------\n");
	testderivativefactor();
	printf("\n- test testexpressiondistribution --------------------------\n");
	testexpressiondistribution();
	printf("\n- test test12varexpression ---------------------------------\n");
	test12varexpression();

#ifdef SPEED_TEST
	printf("\n- test expressionspeed -------------------------------------\n");
	testexpressionspeed();
	printf("\n- test testexpressionmodelspeed ----------------------------\n");
	testexpressionmodelspeed();
	printf("\n- test testexpressionoperatorspeed -------------------------\n");
	testexpressionoperatorspeed();
#endif
}

/* test the speed of generating a random number (evaluating a complicated expression) and
 * taking the derivative (central space approximation)
 */

template<typename T>
double neighbour_average(T* v, int n)
{
	return v[-n] + v[n] - 2 * v[0];
}

template<typename T>
struct MyEquation
{
	T constant;
	int n;
	T eval(Grid<T, 1>& v, int i)
	{
		return rand() / RAND_MAX + neighbour_average(v.values + i, n) + constant;
	}
};

// ************************************************************************************

/* testing the speed of wrapping a function in an object that then calls it
 * on the stored objects
 *
 * compared with the MyEquation object for checking how long different methods take
 */

template<typename T, typename R, typename A>
struct MyWrapper
{
	MyWrapper(T& t, R r, A a) : a{ a }, t{ t }, r{ r } {}
	auto eval(int i)
	{
		return (t.*r)(a, i);
	}

protected:
	A a;
	T& t;
	R r;
};




// ************************************************************************************



// objects used to test a more complicated expression

struct Number;

template<size_t Z>
struct NumberSpecialized;

struct Eta;

template<size_t Z>
struct EtaSpecialized;



// realiasing to distinguish individual variables

using Number_uk = NumberSpecialized<0>;
using Number_vmk = NumberSpecialized<1>;
using Number_umk = NumberSpecialized<2>;
using Number_vk = NumberSpecialized<3>;
using Number_udk = NumberSpecialized<4>;
using Number_vdmk = NumberSpecialized<5>;
using Number_udmk = NumberSpecialized<6>;
using Number_vdk = NumberSpecialized<7>;

using Number_vk_0 = NumberSpecialized<8>;
using Number_umk_0 = NumberSpecialized<9>;
using Number_udmk_0 = NumberSpecialized<10>;
using Number_vdk_0 = NumberSpecialized<11>;

using Eta_dk = EtaSpecialized<0>;
using Eta_mk = EtaSpecialized<1>;
using Eta_dmk = EtaSpecialized<2>;
using Eta_k = EtaSpecialized<3>;


// ********************************************************************

// variables used in test12varexpression function

extern complex_t* u_data;
extern complex_t* v_data;
extern double* thetak;

extern size_t k_len;
extern size_t t_len;


extern Number_uk uk;
extern Number_vmk vmk;
extern Number_umk umk;
extern Number_vk vk;
extern Number_udk udk;
extern Number_vdmk vdmk;
extern Number_udmk udmk;
extern Number_vdk vdk;
extern Eta_dk etadk;
extern Eta_mk etamk;
extern Eta_dmk etadmk;
extern Eta_k etak;

extern Number_vk_0 vk0;
extern Number_umk_0 umk0;
extern Number_udmk_0 udmk0;
extern Number_vdk_0 vdk0;


// ********************************************************************


struct Number
{
	Number() : values{ new complex_t[k_len] } {}
	complex_t* values;

	const complex_t& operator[](iter_type n) const
	{
		return values[n];
	}

	complex_t& operator[](iter_type n)
	{
		return values[n];
	}

	~Number()
	{
		delete[] values;
	}
};

template<size_t Z>
struct NumberSpecialized : Number
{
	NumberSpecialized() : Number()
	{
		fill_data(*this);
	}

	void update_data()
	{
		fill_data(*this);
	}
};

struct Eta
{
	Eta() : x{ 0 } {}
	int x;
};

template<size_t Z>
struct EtaSpecialized : Eta
{
	using Eta::Eta;
};

// add these types to the expression logic and define a couple rules

DEFINE_BASE_TYPE((size_t Z), (NumberSpecialized<Z>), (Number))
DEFINE_BASE_TYPE((size_t Z), (EtaSpecialized<Z>), (Eta))

ADD_EXPR_TYPE((), (Number), values)
ADD_EXPR_TYPE((size_t Z), (NumberSpecialized<Z>), values)
ADD_EXPR_TYPE_POINT((), (Eta), x)
ADD_EXPR_TYPE_POINT((size_t Z), (EtaSpecialized<Z>), x)

ALLOW_COMBINATION((size_t Z), (NumberSpecialized<Z>))
ALLOW_COMBINATION((size_t Z), (EtaSpecialized<Z>))

RESTRICT_COMMUTATIVITY((), (Eta))


// ********************************************************************


template<size_t Z>
void fill_data(NumberSpecialized<Z>&);

template<>
inline void fill_data(Number_uk& data)
{
	for (iter_type i = 0; i < k_len; ++i)
	{
		data.values[i] = u_data[i];
	}
}

template<>
inline void fill_data(Number_umk& data)
{
	for (iter_type i = 0; i < k_len; ++i)
	{
		data.values[i] = u_data[k_len - 1 - i];
	}
}

template<>
inline void fill_data(Number_vk& data)
{
	for (iter_type i = 0; i < k_len; ++i)
	{
		data.values[i] = v_data[i];
	}
}

template<>
inline void fill_data(Number_vmk& data)
{
	for (iter_type i = 0; i < k_len; ++i)
	{
		data.values[i] = v_data[k_len - 1 - i];
	}
}

template<>
inline void fill_data(Number_udk& data)
{
	for (iter_type i = 0; i < k_len; ++i)
	{
		data.values[i] = std::conj(u_data[i]);
	}
}

template<>
inline void fill_data(Number_udmk& data)
{
	for (iter_type i = 0; i < k_len; ++i)
	{
		data.values[i] = std::conj(u_data[k_len - 1 - i]);
	}
}

template<>
inline void fill_data(Number_vdk& data)
{
	for (iter_type i = 0; i < k_len; ++i)
	{
		data.values[i] = std::conj(v_data[i]);
	}
}

template<>
inline void fill_data(Number_vdmk& data)
{
	for (iter_type i = 0; i < k_len; ++i)
	{
		data.values[i] = std::conj(v_data[k_len - 1 - i]);
	}
}

template<>
inline void fill_data(Number_vk_0& data)
{
	for (iter_type i = 0; i < k_len; ++i)
	{
		data.values[i] = v_data[i];
	}
}

template<>
inline void fill_data(Number_umk_0& data)
{
	for (iter_type i = 0; i < k_len; ++i)
	{
		data.values[i] = u_data[k_len - 1 - i];
	}
}

template<>
inline void fill_data(Number_udmk_0& data)
{
	for (iter_type i = 0; i < k_len; ++i)
	{
		data.values[i] = std::conj(u_data[k_len - 1 - i]);
	}
}

template<>
inline void fill_data(Number_vdk_0& data)
{
	for (iter_type i = 0; i < k_len; ++i)
	{
		data.values[i] = std::conj(v_data[i]);
	}
}



// ********************************************************************

// custom identities

template<typename A, typename B>
struct eta_identity
{
	template<typename N1, typename N2>
	static auto apply(N1, N2)
	{
		return OpVoid{};
	}
};


template<>
struct eta_identity<Eta_dk, Eta_dmk>
{
	template<typename N1, typename N2>
	static auto apply(N1, N2)
	{
		return expr::make_term(vk0) * expr::make_term(umk0);
	}
};

template<>
struct eta_identity<Eta_dmk, Eta_dk>
{
	template<typename N1, typename N2>
	static auto apply(N1 a, N2 b)
	{
		return -eta_identity<Eta_dk, Eta_dmk>::template apply(a, b);
	}
};


template<>
struct eta_identity<Eta_dk, Eta_k>
{
	template<typename N1, typename N2>
	static auto apply(N1, N2)
	{
		return expr::make_term(vk0) * expr::make_term(vk0);
	}
};

template<>
struct eta_identity<Eta_k, Eta_dk>
{
	template<typename N1, typename N2>
	static auto apply(N1 a, N2 b)
	{
		return OpIdentity{} -eta_identity<Eta_dk, Eta_k>::template apply(a, b);
	}
};



template<>
struct eta_identity<Eta_mk, Eta_dmk>
{
	template<typename N1, typename N2>
	static auto apply(N1, N2)
	{
		return expr::make_term(umk0) * expr::make_term(umk0);
	}
};

template<>
struct eta_identity<Eta_dmk, Eta_mk>
{
	template<typename N1, typename N2>
	static auto apply(N1 a, N2 b)
	{
		return OpIdentity{} - eta_identity<Eta_mk, Eta_dmk>::template apply(a, b);
	}
};


template<>
struct eta_identity<Eta_mk, Eta_k>
{
	template<typename N1, typename N2>
	static auto apply(N1, N2)
	{
		return expr::make_term(udmk0) * expr::make_term(vdk0);
	}
};

template<>
struct eta_identity<Eta_k, Eta_mk>
{
	template<typename N1, typename N2>
	static auto apply(N1 a, N2 b)
	{
		return -eta_identity<Eta_mk, Eta_k>::template apply(a, b);
	}
};

// ********************************************************************


template<typename A, typename B>
using etai = eta_identity<typename expr::original_data_type<A>::type, typename expr::original_data_type<B>::type>;

template<typename A, typename B>
inline auto eta_identity_apply(A const& a, B const& b)
{
	return etai<A, B>::template apply(a, b);
}

//
//
//template<>
//struct NL_applied_rule<Eta, Eta, Eta, Eta>
//{
//	static const bool enable = true;
//	template<typename... Gs>
//	static auto apply(std::tuple<Gs...> const& datas)
//	{
//		auto [a, b, c, d] = datas;
//
//		auto ab = eta_identity_apply(a, b);
//		auto cd = eta_identity_apply(c, d);
//		auto ac = eta_identity_apply(a, c);
//		auto bd = eta_identity_apply(b, d);
//		auto ad = eta_identity_apply(a, d);
//		auto bc = eta_identity_apply(b, c);
//
//		return ab * cd - ac * bd + ad * bc;
//	}
//};


/* custom pruning algorithm to apply the dual identity to a formed expression
 */


template<typename E>
auto eta_identity_prune(OpExpression<E> const& e)
{
	return *static_cast<E const*>(&e);
}

/*
template<typename T, typename A, typename B, typename std::enable_if<(std::is_same<typename expr::base_data_type<A>::type, Eta>::value && std::is_same<typename expr::base_data_type<B>::type, Eta>::value), int>::type = 0>
auto eta_identity_prune(OpTerms<T, A, B> const& e)
{
	return eta_identity_apply(std::get<0>(e.datas), std::get<1>(e.datas));
}*/



template<typename A1, typename A2>
auto eta_identity_prune(OpAdd<A1, A2> const& e)
{
	return eta_identity_prune(e.a) + eta_identity_prune(e.b);
}

template<typename A1, typename A2>
auto eta_identity_prune(OpBinaryMul<A1, A2> const& e)
{
	return eta_identity_prune(e.a) * eta_identity_prune(e.b);
}

template<typename A1, typename A2>
auto eta_identity_prune(OpBinaryDiv<A1, A2> const& e)
{
	return eta_identity_prune(e.a) / eta_identity_prune(e.b);
}




