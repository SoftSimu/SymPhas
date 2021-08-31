
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
 */


#include "symphasthread.h"


void ThreadPool::main_wait()
{
	std::unique_lock<std::mutex> lck(mm);
	cm.wait(lck, [&] { return pred(); });
}

void ThreadPool::init_false()
{
	for (iter_type i = 0; i < THREADS; ++i)
	{
		status[i] = false;
	}
}

void ThreadPool::notify()
{
	std::lock_guard<std::mutex> lck(mt);
	++pr;
	ct.notify_all();
}

void ThreadPool::sync(size_t i, size_t p)
{
	{
		std::lock_guard<std::mutex> lck(mm);
		status[i] = true;
		cm.notify_one();
	}

	std::unique_lock<std::mutex> lck(mt);
	ct.wait(lck, [&] { return pr == p; });
}

void ThreadPool::main_thr_sync()
{
	main_wait();
	init_false();
	notify();
}

void ThreadPool::init_thr()
{
	/*!
		* The pr variable is used to count thread function completion
		* (in this method they are update, equation, step) by incrementing
		* pr once all threads have have finished it, acting as a sync for
		* combined threads to continue (threads wait for increment).
		*/
	pr = 0;
	init_false();

#ifdef DEBUG
	//Time t("initializing %d threads", THREADS);
#endif
}

void ThreadPool::close_thr()
{
#ifdef DEBUG
	//Time t("joining %d threads", THREADS);
#endif

	for (iter_type i = 0; i < THREADS; ++i)
	{
		thr[i].join();
	}
}

