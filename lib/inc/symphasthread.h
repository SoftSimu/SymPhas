#pragma once

#include <mutex>
#include <condition_variable>
#include <algorithm>
#include <thread>

#include "timer.h"


//! Manages the ability to execute parallel jobs.
/*!
 * Functions can be started on multiple threads, and all of them joined
 * together from the main (calling) thread in order to continue execution.
 */
struct ThreadPool
{

protected:
	std::thread thr[THREADS];
	std::mutex mm, mt;
	std::condition_variable cm, ct;

	bool status[THREADS];
	size_t pr;


	inline bool pred() const
	{
		for (iter_type i = 0; i < THREADS; ++i)
		{
			if (!status[i]) return false;
		}
		return true;
	}

	/*!
	 * Wait until all the threads have finished by checking
	 * that all status return true.
	 */
	void main_wait();

	/*!
	 * Set all threads to false, to reinitialize their jobs.
	 */
	void init_false();

	/*!
	 * Notify worker threads so they can resume their work after sync()
	 * call.
	 */
	void notify();

public:
	ThreadPool(ThreadPool const&) = delete;
	ThreadPool operator=(ThreadPool &) = delete;
	ThreadPool() : status{ false }, pr { 0 } {}

	/*! 
	 * In the case when a series of threads must be synchronized, the 
	 * main_thr_sync method will be called in the main thread and main_thr_sync  
	 * will make the main thread BLOCK until the workers finish.
	 * 
	 * Once the thread finishes its work, this function must be called with an 
	 * incremented time step, initially starting from 0, and will BLOCK until 
	 * all threads are done once this function is called by all threads, the 
	 * main thread will finally  UNBLOCK, and each worker will also UNBLOCK.
	 * 
	 * \param i Thread index.
	 * \param p Step number.
	 */
	void sync(size_t i, size_t p);

	/*!
	 * Blocking call from the main thread, that will wait for all the workers
	 * to finish before starting them again.
	 * the template allows passing of a function that will be called once the 
	 * threads have all finished, but before they are restarted.
	 */
	template<typename F>
	void main_thr_sync(F fn)
	{
		main_wait();
		fn();

		init_false();
		notify();

	}

	void main_thr_sync();
	/*!
	 * Initializes the threading instance.
	 */
	void init_thr();

	/*!
	 * Starts the thread. The threads are immediately initialized.
	 */
	template<typename F, typename  ...Args>
	void start_thr(F fn, Args&& ...args)
	{
		init_thr();

		/*!
		 * Once each status is false, it is safe to begin threading.
		 */
		for (iter_type i = 0; i < THREADS; ++i)
		{
			thr[i] = std::thread(fn, std::forward<Args>(args)..., i);
		}
	}


	void close_thr();


	//! Return in the index of the thread.
	/*!
	 * For the index within the total length, compute the thread index
	 * to which it is assigned.
	 * 
	 * \param i The index in the list.
	 * \param len The total length of the list.
	 */
	static inline iter_type thr_start_i(iter_type i, len_type len)
	{
		iter_type
			count = len / THREADS,
			rem = len - (count * THREADS);

		return count * i + std::min(i, rem);
	}
	static inline iter_type thr_start_i(size_t i, size_t len)
	{
		return thr_start_i(static_cast<iter_type>(i), static_cast<len_type>(len));
	}
	static inline iter_type thr_start_i(size_t i, len_type len)
	{
		return thr_start_i(static_cast<iter_type>(i), len);
	}


	//! Return the compile time constant value of the thread index.
	/*!
	 * For the index within the total length, compute the thread index
	 * to which it is assigned.
	 *
	 * \tparam I The index in the list.
	 * \tparam N The total length of the list.
	 */
	template<len_type N, iter_type I>
	static iter_type constexpr thr_start_i()
	{
		iter_type constexpr
			count = N / THREADS,
			rem = N % THREADS;

		return count * I + fixed_min<I, rem>;
	}


};
