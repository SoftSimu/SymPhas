//#include "symphas.h"
#include <iostream>
#include <utility>
#include <tuple>

#pragma once

template<typename... Ts>
struct variadic_type {};

template<typename T>
struct variadic_entry {};

template<int I0, int I1>
struct variadic_index {};

template<typename T0, typename T1>
struct select_type {};

template<typename T0, typename T10, typename T11>
struct select_type<T0, variadic_type<variadic_entry<T10>, T11>> 
{
    select_type(T0 = {}, variadic_type<variadic_entry<T10>, T11> = {}) {}
};

template<typename T0, typename T10, int I10, int I11>
struct select_type<T0, variadic_type<variadic_entry<T10>, variadic_index<I10, I11>>> 
{
    select_type(T0 = {}, variadic_type<variadic_entry<T10>, variadic_index<I10, I11>> = {}) {}
};

template<typename T0, typename T10, typename T11>
select_type(T0, variadic_type<variadic_entry<T10>, T11>) -> select_type<T0, variadic_type<variadic_entry<T10>, T11>>;

template<typename T0, typename T10, int I10, int I11>
select_type(T0, variadic_type<variadic_entry<T10>, variadic_index<I10, I11>>) -> select_type<T0, variadic_type<variadic_entry<T10>, variadic_index<I10, I11>>>;
    
template<int... I0s, typename... Ts, typename T0, typename T10, typename T11>
auto var_(std::integer_sequence<int, I0s...>, std::tuple<Ts...>, select_type<T0, variadic_type<variadic_entry<T10>, T11>> const&)
{
    return 0;
}

template<int... I0s, typename... Ts, typename T0, typename T10, int I10, int I11>
auto var_(std::integer_sequence<int, I0s...>, std::tuple<Ts...>, select_type<T0, variadic_type<variadic_entry<T10>, variadic_index<I10, I11>>> const&)
{
    return 1;
}

int nothing()
{
    
    auto seq = std::integer_sequence<int, -1, 1, 0>{};
    auto tupl = std::make_tuple(0, "asdf", '6');
    auto v0 = select_type<int, variadic_type<variadic_entry<int>, variadic_entry<void>>>{};
    auto v1 = select_type<int, variadic_type<variadic_entry<int>, variadic_index<0, 0>>>{};
    auto f0 = var_(seq, tupl, v0);
    auto f1 = var_(seq, tupl, v1);
    
    printf("%d %d\n", f0, f1);
    //symphas::b_data_type bdata;
    //symphas::interval_element_type interval;
    //interval.set_interval_count(0, 1, 64);
    //printf("%lf \n", interval.width());
}


