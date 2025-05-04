


/* convenience definitions to standardize the names of the algorithms as they appear in the objects
 */

#define ALG_NAME_VECTOR(NAME) Vector ## NAME
#define ALG_NAME_SCALAR(NAME) Scalar ## NAME
#define ALG_NAME_POINT(NAME) Point ## NAME
#define ALG_NAME_DYNAMIC(NAME) Dynamic ## NAME


 /* these definitions are used to prototype or otherwise declare the function under which the algorithm
  * is actually implemented
  */

#define _ALGORITHM_1D_DECLARATION(NAME, TYPE) \
template<typename T> template<typename model_type> auto __ ## NAME::TYPE ## NAME<axis_1d_type*, T*>::get_data( \
	[[maybe_unused]] const axis_1d_type* data_x, [[maybe_unused]] T* data_y, [[maybe_unused]] len_type len, \
	[[maybe_unused]] len_type L, [[maybe_unused]] model_type const& model)
#define _ALGORITHM_2D_DECLARATION(NAME, TYPE) \
template<typename T> template<typename model_type> auto __ ## NAME::TYPE ## NAME<axis_2d_type*, T*>::get_data( \
	[[maybe_unused]] const axis_2d_type* data_x, [[maybe_unused]] T* data_y, [[maybe_unused]] len_type len, \
	[[maybe_unused]] len_type L, [[maybe_unused]] len_type M, [[maybe_unused]] model_type const& model)
#define _ALGORITHM_3D_DECLARATION(NAME, TYPE) \
template<typename T> template<typename model_type> auto __ ## NAME::TYPE ## NAME<axis_3d_type*, T*>::get_data( \
	[[maybe_unused]] const axis_3d_type* data_x, [[maybe_unused]] T* data_y,  \
	[[maybe_unused]] len_type len, [[maybe_unused]] len_type L, [[maybe_unused]] len_type M, [[maybe_unused]] len_type N, \
	[[maybe_unused]] model_type const& model)

//! Function header for defining the algorithm for a 1D vector process.
/*!
 * Provides the following parameters:
 * > **data_x, data_y, len, L, model**
 * 
 * This process should return data of the type vector_data<Y, 1>.
 * 
 * \param NAME The name of the process, which must match all other \p NAME arguments
 * passed to the other macros associated to the definition of this data process.
 */
#define ALGORITHM_VECTOR_1D_DECLARATION(NAME) _ALGORITHM_1D_DECLARATION(NAME, Vector)
//! Function header for defining the algorithm for a 2D vector process.
/*!
 * Provides the following parameters:
 * **data_x, data_y, len, L, M, model**
 * 
 * This process should return data of the type vector_data<Y, 2>.
 * 
 * \param NAME The name of the process, which must match all other \p NAME arguments
 * passed to the other macros associated to the definition of this data process.
 */
#define ALGORITHM_VECTOR_2D_DECLARATION(NAME) _ALGORITHM_2D_DECLARATION(NAME, Vector)
//! Function header for defining the algorithm for a 3D vector process.
/*!
 * Provides the following parameters:
 * > **data_x, data_y, len, L, M, N, model**
 * 
 * This process should return data of the type vector_data<Y, 3>.
 * 
 * \param NAME The name of the process, which must match all other \p NAME arguments
 * passed to the other macros associated to the definition of this data process.
 */
#define ALGORITHM_VECTOR_3D_DECLARATION(NAME) _ALGORITHM_3D_DECLARATION(NAME, Vector)


//! Function header for defining the algorithm for a 1D scalar process.
/*!
 * Provides the following parameters:
 * > **data_x, data_y, len, L, model**
 * 
 * This process should return data of the type scalar_data<Y>.
 * 
 * \param NAME The name of the process, which must match all other \p NAME arguments
 * passed to the other macros associated to the definition of this data process.
 */
#define ALGORITHM_SCALAR_1D_DECLARATION(NAME) _ALGORITHM_1D_DECLARATION(NAME, Scalar)
//! Function header for defining the algorithm for a 2D scalar process.
/*!
 * Provides the following parameters:
 * > **data_x, data_y, len, L, M, model**
 * 
 * This process should return data of the type scalar_data<Y>.
 * 
 * \param NAME The name of the process, which must match all other \p NAME arguments
 * passed to the other macros associated to the definition of this data process.
 */
#define ALGORITHM_SCALAR_2D_DECLARATION(NAME) _ALGORITHM_2D_DECLARATION(NAME, Scalar)
//! Function header for defining the algorithm for a 3D scalar process.
/*!
 * Provides the following parameters:
 * > **data_x, data_y, len, L, M, N, model**
 * 
 * This process should return data of the type scalar_data<Y>.
 * 
 * \param NAME The name of the process, which must match all other \p NAME arguments
 * passed to the other macros associated to the definition of this data process.
 */
#define ALGORITHM_SCALAR_3D_DECLARATION(NAME) _ALGORITHM_3D_DECLARATION(NAME, Scalar)


//! Function header for defining the algorithm for a 1D point process.
/*!
 * Provides the following parameters:
 * > **data_x, data_y, len, L, model**
 * 
 * This process should return data of the type point_data<Y>.
 * 
 * \param NAME The name of the process, which must match all other \p NAME arguments
 * passed to the other macros associated to the definition of this data process.
 */
#define ALGORITHM_POINT_1D_DECLARATION(NAME) _ALGORITHM_1D_DECLARATION(NAME, Point)
//! Function header for defining the algorithm for a 2D point process.
/*!
 * Provides the following parameters:
 * > **data_x, data_y, len, L, M, model**
 * 
 * This process should return data of the type point_data<Y>.
 * 
 * \param NAME The name of the process, which must match all other \p NAME arguments
 * passed to the other macros associated to the definition of this data process.
 */
#define ALGORITHM_POINT_2D_DECLARATION(NAME) _ALGORITHM_2D_DECLARATION(NAME, Point)
//! Function header for defining the algorithm for a 3D point process.
/*!
 * Provides the following parameters:
 * > **data_x, data_y, len, L, M, N, model**
 * 
 * This process should return data of the type point_data<Y>.
 * 
 * \param NAME The name of the process, which must match all other \p NAME arguments
 * passed to the other macros.
 */
#define ALGORITHM_POINT_3D_DECLARATION(NAME) _ALGORITHM_3D_DECLARATION(NAME, Point)

 //! Function header for defining the algorithm of a dynamic process.
 /*!
  * Provides the following parameters:
  * > **data_x, data_y, len**
  * 
  * The template parameter `Y` refers to the type of the output of the other data processes
  * which feed into this process. For example, if scalar-type data processes were defined,
  * their outputs would be aggregated and then fed into the dynamic algorithm as a pointer array
  * of those Fields (`F**` type parameter), and `Y` would be a pointer type. It is derived
  * from the result of symphas::field_y_t.
  *
  * This process can return any field type.
  *
  * \param NAME The name of the process, which must match all other \p NAME arguments
  * passed to the other macros.
  */
#define ALGORITHM_DYNAMIC_DECLARATION(NAME) \
template<typename F> template<typename model_type, size_t D> auto __ ## NAME :: ALG_NAME_DYNAMIC(NAME)<iter_type*, F**>::get_data( \
	const iter_type* data_x, F** data_y, len_type len, [[maybe_unused]] model_type const& model)



  // ***************************************************************************************************************************************


  /* this defines the general object and prototypes all the functions for each of 1,2 and 3 dimensions
   * in the case that one of them is not defined, then the compiler will issue an error, thus the algorithm
   * has to be defined for all the dimensions
   */

//! Defines a data process returning a vector type.
/*!
 * Introduces a new data process which generates vector-valued data from its input. This
 * macro defines the class used to implement and call the data process as part of the 
 * data process control flow.
 * 
 * \param NAME The name of the process, which must match all other \p NAME arguments
 * passed to the other macros associated to the definition of this data process.
 */
#define DEFINE_ALGORITHM_VECTOR(NAME) \
namespace __ ## NAME { \
template<typename X, typename Y> struct ALG_NAME_VECTOR(NAME) : AlgorithmVectorType<ALG_NAME_VECTOR(NAME), X, Y> {}; \
template<typename T> struct ALG_NAME_VECTOR(NAME)<axis_1d_type*, T*> : AlgorithmVectorType<ALG_NAME_VECTOR(NAME), axis_1d_type*, T*> \
{ \
	using Y = std::remove_const_t<T>; \
	using X = axis_1d_type; \
	template<typename M> \
	static auto get_data(const axis_1d_type* data_x, T* data_y, len_type len, len_type l, M const& model); \
}; \
template<typename T> struct ALG_NAME_VECTOR(NAME)<axis_2d_type*, T*> : AlgorithmVectorType<ALG_NAME_VECTOR(NAME), axis_2d_type*, T*> \
{ \
	using Y = std::remove_const_t<T>; \
	using X = axis_2d_type; \
	template<typename M> \
	static auto get_data(const axis_2d_type* data_x, T* data_y, len_type len, len_type l, len_type m, M const& model); \
}; \
template<typename T> struct ALG_NAME_VECTOR(NAME)<axis_3d_type*, T*> : AlgorithmVectorType<ALG_NAME_VECTOR(NAME), axis_3d_type*, T*> \
{ \
	using Y = std::remove_const_t<T>; \
	using X = axis_3d_type; \
	template<typename M> \
	static auto get_data(const axis_3d_type* data_x, T* data_y, len_type len, len_type l, len_type m, len_type n, M const& model); \
}; \
template<typename X, typename Y> using vector_algorithm = ALG_NAME_VECTOR(NAME)<X, Y>; } \
template<typename X, typename Y> using ALG_NAME_VECTOR(NAME) = __ ## NAME::ALG_NAME_VECTOR(NAME)<X, Y>;


//! Defines a data process returning a scalar type.
/*!
 * Introduces a new data process which generates scalar-valued data from its input. This
 * macro defines the class used to implement and call the data process as part of the 
 * data process control flow.
 * 
 * \param NAME The name of the process, which must match all other \p NAME arguments
 * passed to the other macros associated to the definition of this data process.
 */
#define DEFINE_ALGORITHM_SCALAR(NAME) \
namespace __ ## NAME { \
template<typename X, typename Y> struct ALG_NAME_SCALAR(NAME) : AlgorithmScalarType<ALG_NAME_SCALAR(NAME), X, Y> {}; \
template<typename T> struct ALG_NAME_SCALAR(NAME)<axis_1d_type*, T*> : AlgorithmScalarType<ALG_NAME_SCALAR(NAME), axis_1d_type*, T*> \
{ \
	using Y = std::remove_const_t<T>; \
	template<typename M> \
	static auto get_data(const axis_1d_type* data_x, T* data_y, len_type len, len_type l, M const& model); \
}; \
template<typename T> struct ALG_NAME_SCALAR(NAME)<axis_2d_type*, T*> : AlgorithmScalarType<ALG_NAME_SCALAR(NAME), axis_2d_type*, T*> \
{ \
	using Y = std::remove_const_t<T>; \
	template<typename M> \
	static auto get_data(const axis_2d_type* data_x, T* data_y, len_type len, len_type l, len_type m, M const& model); \
}; \
template<typename T> struct ALG_NAME_SCALAR(NAME)<axis_3d_type*, T*> : AlgorithmScalarType<ALG_NAME_SCALAR(NAME), axis_3d_type*, T*> \
{ \
	using Y = std::remove_const_t<T>; \
	template<typename M> \
	static auto get_data(const axis_3d_type* data_x, T* data_y, len_type len, len_type l, len_type m, len_type n, M const& model); \
}; \
template<typename X, typename Y> using scalar_algorithm = ALG_NAME_SCALAR(NAME)<X, Y>; } \
template<typename X, typename Y> using ALG_NAME_SCALAR(NAME) = __ ## NAME::ALG_NAME_SCALAR(NAME)<X, Y>;


//! Defines a data process returning a point type.
/*!
 * Introduces a new data process which generates point-valued data from its input. This
 * macro defines the class used to implement and call the data process as part of the 
 * data process control flow.
 * 
 * \param NAME The name of the process, which must match all other \p NAME arguments
 * passed to the other macros associated to the definition of this data process.
 */
#define DEFINE_ALGORITHM_POINT(NAME) \
namespace __ ## NAME { \
template<typename X, typename Y> struct ALG_NAME_POINT(NAME) : AlgorithmPointType<ALG_NAME_POINT(NAME), X, Y> {}; \
template<typename T> struct ALG_NAME_POINT(NAME)<axis_1d_type*, T*> : AlgorithmPointType<ALG_NAME_POINT(NAME), axis_1d_type*, T*> \
{ \
	using Y = std::remove_const_t<T>; \
	template<typename M> \
	static auto get_data(const axis_1d_type* data_x, T* data_y, len_type len, len_type l, M const& model); \
}; \
template<typename T> struct ALG_NAME_POINT(NAME)<axis_2d_type*, T*> : AlgorithmPointType<ALG_NAME_POINT(NAME), axis_2d_type*, T*> \
{ \
	using Y = std::remove_const_t<T>; \
	template<typename M> \
	static auto get_data(const axis_2d_type* data_x, T* data_y, len_type len, len_type l, len_type m, M const& model); \
}; \
template<typename T> struct ALG_NAME_POINT(NAME)<axis_3d_type*, T*> : AlgorithmPointType<ALG_NAME_POINT(NAME), axis_3d_type*, T*> \
{ \
	using Y = std::remove_const_t<T>; \
	template<typename M> \
	static auto get_data(const axis_3d_type* data_x, T* data_y, len_type len, len_type l, len_type m, len_type n, M const& model); \
}; \
template<typename X, typename Y> using point_algorithm = ALG_NAME_POINT(NAME)<X, Y>; } \
template<typename X, typename Y> using ALG_NAME_POINT(NAME) = __ ## NAME::ALG_NAME_POINT(NAME)<X, Y>;


//! Defines a data process that takes place across multiple calls.
/*!
 * Introduces a new data process which generates data of any value type from its input. This
 * macro defines the class used to implement and call the data process as part of the 
 * data process control flow. Dynamic means that this data process will aggregate data so that
 * multiple indices can be processed.
 * 
 * \param NAME The name of the process, which must match all other \p NAME arguments
 * passed to the other macros associated to the definition of this data process.
 */
#define DEFINE_ALGORITHM_DYNAMIC(NAME) \
namespace __ ## NAME { \
template<typename X, typename Y> struct ALG_NAME_DYNAMIC(NAME) : AlgorithmDynamicType<ALG_NAME_DYNAMIC(NAME), X, Y> {}; \
template<typename F> struct ALG_NAME_DYNAMIC(NAME)<iter_type*, F**> : AlgorithmDynamicType<ALG_NAME_DYNAMIC(NAME), iter_type*, F**> \
{ \
	using X = symphas::field_x_t<F>; \
	using Y = symphas::field_y_t<F>; \
	template<typename M, size_t D = model_dimension<M>::value> \
	static auto get_data(const iter_type* data_x, F** data_y, len_type len, M const& model); \
}; \
template<typename X, typename Y> using dynamic_algorithm = ALG_NAME_DYNAMIC(NAME)<X, Y>; } \
template<typename X, typename Y> using ALG_NAME_DYNAMIC(NAME) = __ ## NAME::ALG_NAME_DYNAMIC(NAME)<X, Y>;


  // ***************************************************************************************************************************************



/* the definition which brings together all of the data definitions under one convenient object
 * this way the postprocessing utility knows which algorithms a specific postprocessing part will use
 */

//! Define a new data collection group of the given name.
/*!
 * Creates a new data collection group with the given name, which is subsequently accompanied by
 * #DEFINE_ALGORITHM_VECTOR, #DEFINE_ALGORITHM_SCALAR, #DEFINE_ALGORITHM_POINT and/or 
 * #DEFINE_ALGORITHM_DYNAMIC, with the same name provided to each. The implementation of
 * a data process is then complete using the algorithm declaration macros, e.g. 
 * #ALGORITHM_VECTOR_1D_DECLARATION.
 */
#define DEFINE_DATA(NAME, DEFAULT_PROC, ...) \
namespace __ ## NAME { \
inline std::vector<DataParams> params{}; \
struct DataApplied : Data<__VA_ARGS__> \
{ \
	using Data<__VA_ARGS__>::Data; \
	static bool add_default(const char* name, ProcessType t) { symphas::internal::default_procs_map[name] = t; return true; } \
	static bool add_process(const char* name) { symphas::internal::data_map[name] = &params; return true; } \
}; \
} \
using Data ## NAME = __ ## NAME::DataApplied; \
inline auto __DATA_DEFAULT_PROC_ ## NAME = ProcessType:: DEFAULT_PROC; \
inline auto __ADD_DEFAULT_ ## NAME = Data ## NAME::add_default(#NAME, __DATA_DEFAULT_PROC_ ## NAME); \
inline auto __ADD_PROCESS_ ## NAME = Data ## NAME::add_process(#NAME);

#define DEFINE_DATA_ALIAS_NAME(NAME, ALIAS) \
inline auto __ADD_DEFAULT_ ## ALIAS ## _ ## NAME = Data ## NAME::add_default(#ALIAS, __DATA_DEFAULT_PROC_ ## NAME); \
inline auto __ADD_PROCESS_ ## ALIAS ## _ ## NAME = Data ## NAME::add_process(#ALIAS);




/* these are the VA_ARGS passed to the above definition in order to hook the data in with
 * the specific algorithms that were implemented for the postprocessing function
 */

#define ALG_POINT point_algorithm
#define ALG_VECTOR vector_algorithm
#define ALG_SCALAR scalar_algorithm
#define ALG_DYNAMIC dynamic_algorithm


#define _RUN_ALGORITHM_1D_ARGS(NAME, TYPE) TYPE ## NAME <axis_1d_type*, T*> :: get_data
#define _RUN_ALGORITHM_2D_ARGS(NAME, TYPE) TYPE ## NAME <axis_2d_type*, T*> :: get_data
#define _RUN_ALGORITHM_3D_ARGS(NAME, TYPE) TYPE ## NAME <axis_3d_type*, T*> :: get_data

#define ALGORITHM_ARGS_1D data_x, data_y, len, L, model
#define ALGORITHM_ARGS_2D data_x, data_y, len, L, M, model
#define ALGORITHM_ARGS_3D data_x, data_y, len, L, M, N, model

#define ALGORITHM_ARGS_DATA_Y_1D(DATA_Y) data_x, DATA_Y, len, L, model
#define ALGORITHM_ARGS_DATA_Y_2D(DATA_Y) data_x, DATA_Y, len, L, M, model
#define ALGORITHM_ARGS_DATA_Y_3D(DATA_Y) data_x, DATA_Y, len, L, M, N, model

#define ALGORITHM_ARGS_DATA_XY_1D(DATA_X, DATA_Y) DATA_X, DATA_Y, len, L, model
#define ALGORITHM_ARGS_DATA_XY_2D(DATA_X, DATA_Y) DATA_X, DATA_Y, len, L, M, model
#define ALGORITHM_ARGS_DATA_XY_3D(DATA_X, DATA_Y) DATA_X, DATA_Y, len, L, M, N, model

#define RUN_ALGORITHM_VECTOR_1D(NAME) _RUN_ALGORITHM_1D_ARGS(NAME, Vector)
#define RUN_ALGORITHM_VECTOR_2D(NAME) _RUN_ALGORITHM_2D_ARGS(NAME, Vector)
#define RUN_ALGORITHM_VECTOR_3D(NAME) _RUN_ALGORITHM_3D_ARGS(NAME, Vector)
#define RUN_ALGORITHM_SCALAR_1D(NAME) _RUN_ALGORITHM_1D_ARGS(NAME, Scalar)
#define RUN_ALGORITHM_SCALAR_2D(NAME) _RUN_ALGORITHM_2D_ARGS(NAME, Scalar)
#define RUN_ALGORITHM_SCALAR_3D(NAME) _RUN_ALGORITHM_3D_ARGS(NAME, Scalar)
#define RUN_ALGORITHM_POINT_1D(NAME) _RUN_ALGORITHM_1D_ARGS(NAME, Point)
#define RUN_ALGORITHM_POINT_2D(NAME) _RUN_ALGORITHM_2D_ARGS(NAME, Point)
#define RUN_ALGORITHM_POINT_3D(NAME) _RUN_ALGORITHM_3D_ARGS(NAME, Point)



// ***************************************************************************************************************************************


/* the post processing needs to be joined under one roof; for the data processing, the data saving and the
 * processing type
 *
 * there is a definition for both the static and dynamic types of processes
 *
 * moreover it will automatically be able to connect with the data variable which is automatically populated by the
 * dataconfig functionality on program startup
 */

#define DEFINE_COLLECTION_GENERIC(TYPE, NAME) \
template<typename X, typename Y> \
struct Collect ## NAME : Process ## TYPE<X, Y, Save ## NAME, Data ## NAME> \
{ \
	Collect ## NAME() : Process ## TYPE<X, Y, Save ## NAME, Data ## NAME>(*symphas::internal::get_data_map()[#NAME]) {} \
}; \



#define DEFINE_COLLECTION_STATIC(NAME) DEFINE_COLLECTION_GENERIC(Static, NAME)
#define DEFINE_COLLECTION_DYNAMIC(NAME) DEFINE_COLLECTION_GENERIC(Dynamic, NAME)


