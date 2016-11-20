/*
 * Copyright (c) 2016 Felix Lelchuk
 *
 * This software is provided 'as-is', without any express or implied
 * warranty. In no event will the authors be held liable for any damages
 * arising from the use of this software.
 *
 * Permission is granted to anyone to use this software for any purpose,
 * including commercial applications, and to alter it and redistribute it
 * freely, subject to the following restrictions:
 *
 * 1. The origin of this software must not be misrepresented; you must not
 *    claim that you wrote the original software. If you use this software
 *    in a product, an acknowledgment in the product documentation would be
 *    appreciated but is not required.
 * 2. Altered source versions must be plainly marked as such, and must not be
 *    misrepresented as being the original software.
 * 3. This notice may not be removed or altered from any source distribution.
 **/

/**
 * @file matio_cpp.h
 * @brief C++11 interface to MAT I/O 
 * @author Felix Lelchuk
 * @version 20/11/2016
 */

#ifndef MATIO_CPP_H
#define MATIO_CPP_H

#include <matio.h>

#include <cassert>
#include <cinttypes>

#include <string>
#include <array>
#include <algorithm>
#include <numeric>
#include <memory>

/** @brief Namespace containing the C++ 11 interface to MAT I/O  */
namespace MatioCPP
{
   using SizeType = size_t;

   /**
    * @brief Base class for n-dimensional array memory ordering.
    * @tparam N Number of dimensions
    *
    * This class pre-calculates strides for the addressing of a multi-dimensional
    * array's memory. The resulting strides a stores inside the class.
    */
   template <size_t N>
   class MemOrder
   {
   public:
      static const size_t Dims = N; /**< Holds the dimension template parameter N. */

      /**
       * @brief Calculates the strides for the given dimensions.
       * @param dims std::array of the sizes in each dimension.
       */
      virtual void calcStrides(const std::array<SizeType, N>& dims) = 0;

      /**
       * @brief Returns the calculated strides.
       * @return const_iterator to the head of the stride list.
       * @pre calcStrides() was called.
       */
      inline typename std::array<SizeType, N>::const_iterator enumStride() const
      {
         return _strides.cbegin();
      }

   protected:
      std::array<SizeType, N> _strides;   /**< Stide list */
   };

   /**
    * @brief Class implementing row-major memory ordering.
    * @tparam N Number of dimensions
    *
    * This class pre-calculates strides for row-major addressing and stores them.
    * Row-major ordering is the default memory layout for C/C++ an many other languages.
    */
   template <size_t N>
   class RowMajor : public MemOrder<N>
   {
   public:
      /**
       * @brief Calculates the strides for the given dimensions in a row-major fashion.
       * @param dims std::array of the sizes in each dimension.
       */
      virtual void calcStrides(const std::array<SizeType, N>& dims)
      {
         _strides.back() = 1;

         if (N > 1)
            std::partial_sum(dims.cbegin() + 1, dims.cend(), _strides.rbegin() + 1, std::multiplies<SizeType>());
      }
   };

   /**
    * @brief Class implementing column-major memory ordering.
    * @tparam N Number of dimensions
    *
    * This class pre-calculates strides for column-major addressing and stores them.
    * Color-major ordering is the default memory layout for MATLAB and Fortran.
    */
   template <size_t N>
   class ColumnMajor : public MemOrder<N>
   {
   public:
      /**
       * @brief Calculates the strides for the given dimensions in a column-major fashion.
       * @param dims std::array of the sizes in each dimension.
       */
      virtual void calcStrides(const std::array<SizeType, N>& dims)
      {
         if (N > 1)
            std::partial_sum(dims.cbegin(), dims.cend() - 1, this->_strides.begin() + 1, std::multiplies<SizeType>());

         _strides.front() = 1;
      }
   };

   /**
    * @brief Index into a multi-dimensional array
    * @tparam N Number of dimensions
    */
   template <size_t N>
   class Index;

   template <>
   class Index<1>
   {
   public:
      static const size_t Dims = 1;

      inline Index(int last) : _pos(last)
      {
      }

      template <typename StrideList>
      inline int toLinear(StrideList strideList) const
      {
         return _pos * (*strideList);
      }

      template <typename DimList>
      inline void checkBounds(DimList dims) const
      {
         assert(_pos >= 0);
         assert(_pos < *dims);
      }

   private:
      int _pos;
   };

   template <size_t N>
   class Index
   {
   public:
      static const size_t Dims = N;

      template <typename... Rest>
      Index(int first, Rest... rest) :
         _inner(rest...),
         _pos(first)
      {
      }

      template <typename StrideList>
      inline int toLinear(StrideList strideList) const
      {
         return _pos * (*strideList) + _inner.toLinear(++strideList);
      }

      template <typename DimList>
      inline void checkBounds(DimList dims) const
      {
         assert(_pos >= 0);
         assert(_pos < *dims);

         _inner.checkBounds(++dims);
      }

   private:
      Index<N - 1> _inner;
      int _pos;
   };

   /**
    * @brief Order-independent base class of MultiArray.
    * @tparam T Element type.
    * @tparam N Number of dimensions.
    */
   template <typename T, size_t N>
   class ArrayBase
   {
   public:
      /**
       * @brief Returns a pointer the contiguous memory backing the array.
       * @return Pointer to the first element.
       */
      inline const T* data() const
      {
         return _data;
      }

      /**
       * @brief Returns a pointer the contiguous memory backing the array.
       * @return Pointer to the first element.
       */
      inline T* data()
      {
         return _data;
      }

      /**
       * @brief Returns the total number of elements, in all dimensions.
       * @return Overall number of elements.
       */
      inline SizeType size() const
      {
         return _size;
      }

      /**
       * @brief Returns the number of elements in a certain dimension.
       * @param n Zero-based number of the dimension.
       * @return Number of elements in the dimension of interest.
       */
      inline SizeType dim(SizeType n) const
      {
         return _dims.at(n);
      }

      /**
       * @brief Returns a N-element array containing the number of elements in each dimension.
       * @return std::array with the number of elements for each dimension.
       */
      inline const std::array<SizeType, N>& dims() const
      {
         return _dims;
      }

      /**
       * @brief Sanity-check if the linear index is in range.
       * @param index Index to check.
       * 
       * This method uses C's assert() to report errors.
       */
      void checkBounds(int index) const
      {
         assert(index >= 0);
         assert(index < size());
      }

      /**
       * @brief Return an immutable reference to an array element by its linear index.
       * @param Linear index of the element.
       * @return Constant element reference.
       */
      inline const T& linear(int index) const
      {
         //checkBounds(index);
         return data()[index];
      }

      /**
       * @brief Return a mutable reference to an array element by its linear index.
       * @param Linear index of the element.
       * @return Mutable element reference.
       */
      inline T& linear(int index)
      {
         //checkBounds(index);
         return data()[index];
      }

      /**
       * @brief Return an immutable reference to the first array element.
       * @param Linear index of the element.
       * @return Constant element reference.
       */
      inline const T& value() const
      {
         return linear(0);
      }

      /**
       * @brief Return a mutable reference to the first array element.
       * @param Linear index of the element.
       * @return Constant element reference.
       */
      inline T& value()
      {
         return linear(0);
      }

   protected:
      /**
       * @internal
       * @brief Constructs an object for the given dimensions.
       * @param dims Number of elements in each of the N dimensions.
       *
       * Internal to the implementation. Does not allocate memory.
       * This is left to the MultiArray.
       */
      ArrayBase(const std::array<SizeType, N>& dims) :
         _data(nullptr),
         _dims(dims),
         _ownData(true)
      {
         _size = std::accumulate(dims.cbegin(), dims.cend(), SizeType(1), std::multiplies<SizeType>());
      }

      /**
       * @brief Destroys the array, freeing owned memory through setData().
       */
      ~ArrayBase()
      {
         setData(nullptr);
      }

   protected:
      /**
       * @brief Assign a new backing memory buffer to this array.
       * @param data New memory to use.
       * @postcondition The old memory is freed (if previously owned).
       * @postcondition The newly assigned data is not owned.
       *
       * Frees the old memory if owned. Assigns new memory which is NOT
       * owned by this object.
       */
      inline void setData(T *data)
      {
         if (data == _data)
            return;

         if (_data && _ownData)
         {
            delete[] _data;
         }

         _ownData = false;
         _data = data;
      }

      /**
       * @brief Creates a new memory region for backing and claims ownership thereof.
       * @param numElements Overall number of elements.
       * @postcondition setData() was called, see its postconditions!
       * @postcondition This object holds ownership of the newly allocated memory.
       */
      inline void createData(SizeType numElements)
      {
         T* t = new T[numElements]();
         setData(t);
         _ownData = true;
      }

   protected:
      T *_data;                        /**< Pointer to the backing memory. */
      SizeType _size;                  /**< Overall number of elements in the array (all dimension). */
      std::array<SizeType, N> _dims;   /**< Numbers of elements in each of the N dimensions. */
      bool _ownData;                   /**< Free backing memory when no longer used? */
   };

   /**
    * @brief Multi-dimensional array with custom memory ordering.
    * @tparam T Element type.
    * @tparam N Number of dimensions.
    * @tparam Order The Order class to use.
    */
   template <typename T, size_t N, template <size_t> class Order>
   class MultiArray : public ArrayBase<T, N>
   {
   public:
      /**
       * @brief Constructs a multi-dimensional array with the given number of elements in each dimension.
       * @param first Number of elements in the first dimension.
       * @param rest Number of elements in the following dimensions.
       *
       * This constructor allows to supply the element numbers as arguments rather than using std::array.
       */
      template <typename... Rest>
      MultiArray(SizeType first, Rest... rest)
         : ArrayBase<T, N>({ first, rest... })
      {
         init();
      }

      /**
       * @brief Creates a copy of the supplies MultiArray.
       * @param other Original array.
       *
       * The newly created copy has the same dimensions as the original and contains an independent
       * copy of the original array data.
       */
      MultiArray(const MultiArray<T, N, Order>& other)
         : ArrayBase<T, N>(other._dims)
      {
         init();
         std::copy(other._data, other._data + other.size(), _data);
      }

      /**
       * @brief Constructs a multi-dimensional array with the given number of elements in each dimension.
       * @param dims std::array of the element numbers per dimension.
       */
      MultiArray(const std::array<SizeType, N>& dims) : ArrayBase<T, N>(dims)
      {
         init();
      }

      /**
       * @brief Constructs a MultiArray with the given number of elements in each dimension reusing existing memory.
       * @param dims std::array of the element numbers per dimension.
       * @param data Pointer to the exisiting memory.
       * @postcondition The memory pointed to by the data argument may must persist until destruction of this object.
       * @postcondition The memory is NOT owned by this object and must be freed elsewhere.
       *
       * Create a MultiArray based on existing memory. E.g. allows severals MultiArrays to share the same memory.
       */
      MultiArray(const std::array<SizeType, N>& dims, T* data) : ArrayBase<T, N>(dims)
      {
         init(data);
      }

      /**
       * @brief Assigns the other MultiArray object to this object.
       * @param other Original array.
       * @return Immutable reference to this object.
       * @postcondition This object's size is adjusted to the argument's.
       * @postcondition This object holds an independent copy of the argument's data.
       */
      const MultiArray<T, N, Order>& operator=(const MultiArray<T, N, Order>& other)
      {
         if (&other == this)
            return;

         init();
         std::copy(other._data, other._data + other.size(), _data);
         return *this;
      }

      /**
       * @brief Returns a constant reference to the element at the given Index.
       * @param index N-dimensional Index of the element.
       * @return Constant reference to the element.
       * @precondition The index argument must be within the bounds of the array.
       */
      inline const T& at(const Index<N> &index) const
      {
         //index.checkBounds(_dims.cbegin());
         return linear(index.toLinear(_order.enumStride()));
      }

      /**
       * @brief Returns a mutable reference to the element at the given Index.
       * @param index N-dimensional Index of the element.
       * @return Mutable reference to the element.
       * @precondition The index argument must be within the bounds of the array.
       */
      inline T& at(const Index<N> &index)
      {
         //index.checkBounds(_dims.cbegin());
         return linear(index.toLinear(_order.enumStride()));
      }

      /**
       * @brief Returns a constant reference to the element at the given Index.
       * @param first Index into the first dimension.
       * @param rest Index into the following dimensions.
       * @return Constant reference to the element.
       * @precondition The index must be within the bounds of the array.
       *
       * This version of the at() method allows to pass the
       * sub-indices as parameters, rather than using an Index object.
       */
      template <typename... Rest>
      inline const T& at(int first, Rest... rest) const
      {
         return at(Index<N>(first, rest...));
      }

      /**
       * @brief Returns a mutable reference to the element at the given Index.
       * @param first Index into the first dimension.
       * @param rest Index into the following dimensions.
       * @return Mutable reference to the element.
       * @precondition The index must be within the bounds of the array.
       *
       * This version of the at() method allows to pass the
       * sub-indices as parameters, rather than using an Index object.
       */
      template <typename... Rest>
      inline T& at(int first, Rest... rest)
      {
         return at(Index<N>(first, rest...));
      }

   private:

      /**
       * @internal
       * @brief Allocates memory for the array and calculates strides.
       * @postcondition Postconditions of createData() apply.
       * @postcondition The array received its very own backing memory buffer.
       * @postcondition The strides were calculated according to the Order class.
       */
      inline void init()
      {
         SizeType sz = size();
         createData(sz);
         _order.calcStrides(_dims);
      }

      /**
       * @internal
       * @brief Assigns existing memory for the array and calculates strides.
       * @postcondition Postconditions of setData() apply.
       * @postcondition The array is configured to use (but not own) the given memory.
       * @postcondition The strides were calculated according to the Order class.
       */
      inline void init(T *data)
      {
         setData(data);
         _order.calcStrides(_dims);
      }

   private:
      Order<N> _order;  /**< Instance of the Order class. Contains stides for addressing. */
   };

   /** @internal */
   template <typename T>
   struct MatioCompat
   {
      static const matio_classes cls = MAT_C_EMPTY;
      static const matio_types typ = MAT_T_UNKNOWN;
   };

   template <>
   struct MatioCompat<uint8_t>
   {
      static const matio_classes cls = MAT_C_UINT8;
      static const matio_types typ = MAT_T_UINT8;
   };

   template <>
   struct MatioCompat<int8_t>
   {
      static const matio_classes cls = MAT_C_INT8;
      static const matio_types typ = MAT_T_INT8;
   };

   template <>
   struct MatioCompat<uint16_t>
   {
      static const matio_classes cls = MAT_C_UINT16;
      static const matio_types typ = MAT_T_UINT16;
   };

   template <>
   struct MatioCompat<int16_t>
   {
      static const matio_classes cls = MAT_C_INT16;
      static const matio_types typ = MAT_T_INT16;
   };

   template <>
   struct MatioCompat<uint32_t>
   {
      static const matio_classes cls = MAT_C_UINT32;
      static const matio_types typ = MAT_T_UINT32;
   };
   
   template <>
   struct MatioCompat<int32_t>
   {
      static const matio_classes cls = MAT_C_INT32;
      static const matio_types typ = MAT_T_INT32;
   };

   template <>
   struct MatioCompat<uint64_t>
   {
      static const matio_classes cls = MAT_C_UINT64;
      static const matio_types typ = MAT_T_UINT64;
   };

   template <>
   struct MatioCompat<int64_t>
   {
      static const matio_classes cls = MAT_C_INT64;
      static const matio_types typ = MAT_T_INT64;
   };

   template <>
   struct MatioCompat<float>
   {
      static const matio_classes cls = MAT_C_SINGLE;
      static const matio_types typ = MAT_T_SINGLE;
   };

   template <>
   struct MatioCompat<double>
   {
      static const matio_classes cls = MAT_C_DOUBLE;
      static const matio_types typ = MAT_T_DOUBLE;
   };

   /**
    * @brief Copy-on-write smart pointer.
    *
    * The held object can be pointed to by several ImplicitlySharedDataPtr.
    * However, each write access (non-const dereferencing) of the smart pointer
    * will cause the object to be duplicated, if more than one
    * ImplicitlySharedDataPtr holds the object.
    */
   template <class T>
   class ImplicitlySharedDataPtr
   {
   public:
      typedef std::shared_ptr<T> RefPtr;  /**< Implementation of the shared_ptr. */

   private:
      RefPtr _ref;                        /**< Contained shared_ptr. */

      /**
       * @brief Method to duplicate the held object, if required.
       * @postcondition The newly held object may now be a copy of the previously held.
       */
      void detach()
      {
         T* tmp = _ref.get();
         if (tmp != nullptr && !_ref.unique())
         {
            _ref = RefPtr(new T(*tmp));
         }
      }

   public:
      /**
       * @brief Constructs the ImplicitlySharedDataPtr to hold the object pointed to by the argument.
       * @param t Object to hold and take ownership of.
       * @postcondition Takes ownership of the given object.
       */
      ImplicitlySharedDataPtr(T* t) :
         _ref(t)
      {
      }

      /**
       * @brief Constructs the ImplicitlySharedDataPtr to hold the object pointed to by the argument.
       * @param t Object to hold.
       */
      ImplicitlySharedDataPtr(const RefPtr& refptr) :
         _ref(refptr)
      {
      }

      /**
       * @brief Check if any object is being held.
       * @return false if the object is nullptr, true otherwise.
       */
      inline operator bool() const
      {
         return (bool)_ref;
      }
      
      /**
       * @brief Dereference immutably.
       * @return Constant reference to the held object.
       * @precondition Must hold an object, not nullptr.
       * @see operator bool()
       */
      inline const T& operator*() const
      {
         return *_ref;
      }
      
      /**
       * @brief Dereference mutably.
       * @return Mutable reference to the held object.
       * @precondition Must hold an object, not nullptr.
       * @see operator bool()
       */
      inline T& operator*()
      {
         detach();
         return *_ref;
      }

      /**
       * @brief Immutable pointer access.
       * @return Constant pointer to the held object.
       * @see operator bool()
       */
      inline const T* operator->() const
      {
         return _ref.operator->();
      }

      /**
       * @brief Mutable pointer access.
       * @return Mutable pointer to the held object.
       * @see operator bool()
       */
      inline T* operator->()
      {
         detach();
         return _ref.operator->();
      }
   };

   /** @brief Base class for MATLAB variables. */
   template <typename T, size_t N, matio_classes MC = MatioCompat<T>::cls, matio_types MT = MatioCompat<T>::typ>
   class MatVarBase
   {
   public:
      using ArrayType = MultiArray<T, N, ColumnMajor>;   /**< Type of the underlying MultiArray template class. */

   protected:
      /** @internal */
      struct PrivateData
      {
         PrivateData() :
            var(nullptr)
         {
         }

         PrivateData(matvar_t *variable) :
            var(variable)
         {
         }

         virtual ~PrivateData()
         {
            if (var) Mat_VarFree(var);
         }

         matvar_t *var;
      };

   public:
      /**
       * @brief Constructs an invalid object.
       * @postcondition isValid() returns false on this object.
       */
      MatVarBase() :
         _private(nullptr)
      {
      }

      /**
       * @brief Check if the object is valid.
       * @return true if this is valid, false otherwise.
       */
      inline bool isValid() const
      {
         return (bool)_private;
      }

      /**
       * @brief Returns the variable's name.
       * @return std::string containing the variable name.
       * @precondition This object is valid.
       * @see isValid()
       */
      inline std::string name() const
      {
         if (_private->var->name)
            return _private->var->name;
         else
            return std::string();
      }

      /**
       * @brief Returns the variable's name.
       * @param name std::string containing the variable name.
       * @precondition This object is valid.
       * @postcondition The string holding the old name was freed.
       * @postcondition The name was changed.
       * @see isValid()
       */
      inline void setName(const std::string& name)
      {
         char* &n = _private->var->name;

         if (n)
         {
            free(n);
         }

         n = malloc(name.length() + 1);
         size_t copied = name.copy(n, name.length());
         n[copied] = '\0';
      }

      /**
       * Stores this variable to a matio file.
       * @param file The matio file to write to.
       * @param compress matio compression flag. Defaults to no compression.
       * @precondition This object is valid.
       * @postcondition The variable was written to the given matio file.
       * @see isValid()
       */
      inline void write(mat_t *file, matio_compression compress = MAT_COMPRESSION_NONE) const
      {
         Mat_VarWrite(file, _private->var, compress);
      }

   protected:
      /** @internal */
      MatVarBase(const ImplicitlySharedDataPtr<PrivateData> &priv) :
         _private(priv)
      {
      }

      static matvar_t* readBasic(mat_t *file, const std::string& name)
      {
         matvar_t *var = Mat_VarReadInfo(file, name.c_str());
         if (!var) goto leave;

         // Check if the variable we just read from the file is compatible with this array.
         if (var->rank != N) goto del;
         if (var->class_type != MC) goto del;
         if (var->data_type != MT) goto del;

         return var;

      del:
         Mat_VarFree(var);
      leave:
         return nullptr;
      }

      static bool readData(mat_t *file, matvar_t *var)
      {
         return Mat_VarReadDataAll(file, var) == 0;
      }

   protected:
      ImplicitlySharedDataPtr<PrivateData> _private;  /**< Smart pointer to private data structure. */
   };

   /**
    * @brief Multidimensional non-complex MATLAB variable.
    * @tparam T Element type.
    * @tparam N Number of dimensions.
    * @tparam MC matio_classes enum value for given T.
    * @tparam MT matio_types enum value for given T.
    */
   template <typename T, size_t N, matio_classes MC = MatioCompat<T>::cls, matio_types MT = MatioCompat<T>::typ>
   class MatVar : public MatVarBase<T, N, MC, MT>
   {
   protected:
      using Base = MatVarBase<T, N, MC, MT>; /**< Parameterized type of the MatVarBase base class. */

      /** @internal */
      struct PrivateData : Base::PrivateData
      {
         PrivateData(const std::string& name, const std::array<SizeType, N> &dims) :
            arr(new ArrayType(dims))
         {
            var = Mat_VarCreate(name.c_str(), MC, MT, N, const_cast<SizeType *>(dims.data()), arr->data(), MAT_F_DONT_COPY_DATA);
            assert(var != nullptr);
         }

         // NOTE: The PrivateData takes ownership of the variable.
         PrivateData(const std::array<SizeType, N> &dims, matvar_t *variable) :
            Base::PrivateData(variable),
            arr(new ArrayType(dims, static_cast<T*>(variable->data)))
         {
            // Here, the variable reserved the memory.
            // Have the Array use it (NOT own it).
         }

         PrivateData(const PrivateData& other)
         {
            if (other.arr)
            {
               // Duplicate the array and the variable.
               // The variable will access (NOT own) the Array's data.
               arr.swap(std::unique_ptr<ArrayType>(new ArrayType(*other.arr)));

               var = Mat_VarDuplicate(other.var, 0);  // 0 = don't copy data
               assert(var != nullptr);
               var->data = arr->data();               // replace the variable's data pointer but...
               var->mem_conserve = 1;                 // ...tell the variable that memory is NOT owned
            }
         }

         std::unique_ptr<ArrayType> arr;
      };

   public:
      /**
       * @brief Constructs an invalid MatVar.
       */
      MatVar()
      {
      }

      /**
       * @brief Constructs a MatVar with given name and dimensional extents.
       * @param name Name to assign to the variable.
       * @param dims Number of elements in each of the N dimensions.
       */
      MatVar(const std::string &name, std::array<SizeType, N> dims) :
         Base(new PrivateData(name, dims))
      {
      }

      /**
       * @brief Constructs a MatVar with given name and dimensional extents.
       * @param name Name to assign to the variable.
       * @param dim1 Number of elements in the first dimension.
       * @param rest Number of elements in the other N-1 dimensions.
       */
      template <typename... Rest>
      MatVar(const std::string &name, SizeType dim1, Rest... rest) :
         MatVar(name, { dim1, rest... })
      {
      }

      /**
       * @brief Reads a MatVar with given name from the given matio file.
       * @param file The matio file to read from.
       * @param name Name of the variable to read.
       * @postcondition The MatVar is invalid if the variable is not found in the matio file.
       * @postcondition The MatVar is invalid if the variable's rank is not N.
       * @postcondition The MatVar is invalid if the variable's matio_classes class is not MC.
       * @postcondition The MatVar is invalid if the variable's matio_types type is not MT.
       * @postcondition The MatVar is invalid if the variable is complex.
       * @postcondition The MatVar's dimensions are according to the variable stored in the matio file.
       * @postcondition The MatVar contains the data from the matio file.
       */
      MatVar(mat_t *file, const std::string &name)
      {
         std::array<SizeType, N> dims;

         // Read basic variable info from the file.
         matvar_t *var = Base::readBasic(file, name);
         if (!var) return;

         // Check additional constraints
         if (var->isComplex) goto del;

         // Extract dimensions from the variable.
         std::copy(var->dims, var->dims + N, dims.begin());

         if (!Base::readData(file, var)) goto del;
         _private = new PrivateData(dims, var);
         return;
         
      del:
         Mat_VarFree(var); // This will leave _private uninitialized and thus the object invalid.
      }

      /**
        * @brief Returns a constant reference to the underlying MultiArray.
        * @precondition This object is valid.
        * @see isValid()
        */
      inline const ArrayType& array() const
      {
         return *static_cast<const PrivateData&>(*_private).arr;
      }

      /**
       * @brief Returns a mutable reference to the underlying MultiArray.
       * @precondition This object is valid.
       * @see isValid()
       */
      inline ArrayType& array()
      {
         return *static_cast<PrivateData&>(*_private).arr;
      }

      /**
       * @brief Returns the number of elements in the given dimension.
       * @param n Zero-based dimension number.
       * @return Number of elements in the n'th dimension.
       * @precondition This object is valid.
       * @precondition n is in bounds (0 <= n < N)
       * @see isValid()
       */
      inline SizeType dim(SizeType n) const
      {
         return array().dim(n);
      }

      /**
       * @brief Returns true if any dimension has zero elements. False otherwise.
       * @param 
       * @precondition This object is valid.
       * @see isValid()
       */
      inline bool isEmpty() const
      {
         const std::array<SizeType, N> &dims = array().dims();
         return std::find(dims.cbegin(), dims.cend(), SizeType(0)) != dims.cend();
      }

      /**
       * @brief Constant pointer access. Returns a constant pointer to the underlying MultiArray.
       * @precondition This object is valid.
       * @see isValid()
       */
      inline const ArrayType* operator ->() const
      {
         return static_cast<const PrivateData&>(*_private).arr.get();
      }

      /**
       * @brief Mutable pointer access. Returns a mutable pointer to the underlying MultiArray.
       * @precondition This object is valid.
       * @see isValid()
       */
      inline ArrayType* operator ->()
      {
         return static_cast<PrivateData&>(*_private).arr.get();
      }
   };

   /**
    * @brief Non-complex MATLAB row vector.
    * @tparam T Element type.
    * @tparam MC matio_classes enum value for given T.
    * @tparam MT matio_types enum value for given T.
    */
   template <typename T, matio_classes MC = MatioCompat<T>::cls, matio_types MT = MatioCompat<T>::typ>
   class MatVec : public MatVar<T, 2, MC, MT>
   {
   public:
      using Parent = MatVar<T, 2, MC, MT>; /**< Parameterized type of the MatVar parent class. */

      /**
       * @brief Constructs an invalid MatVec.
       */
      MatVec()
      {
      }

      /**
       * @brief Constructs a MatVec with given name and length.
       * @param name Name to assign to the variable.
       * @param length Number of elements in the first dimension.
       */
      MatVec(const std::string &name, SizeType length) :
         Parent(name, length, 1)
      {
      }

      /**
       * @brief Reads a MatVec with given name from the given matio file.
       * @param file The matio file to read from.
       * @param name Name of the variable to read.
       * @postcondition The MatVec is invalid if the variable is not found in the matio file.
       * @postcondition The MatVec is invalid if the variable's rank is not 2.
       * @postcondition The MatVec is invalid if the variable's matio_classes class is not MC.
       * @postcondition The MatVec is invalid if the variable's matio_types type is not MT.
       * @postcondition The MatVec is invalid if the variable is complex.
       * @postcondition The MatVec's dimensions are according to the variable stored in the matio file.
       * @postcondition The MatVec contains the data from the matio file.
       */
      MatVec(mat_t *file, const std::string &name)
      {
         std::array<SizeType, 2> dims;

         // Read basic variable info from the file.
         matvar_t *var = Parent::readBasic(file, name);
         if (!var) return;

         // Check additional constraints
         if (var->isComplex) goto del;
         if (var->dims[1] != 1) goto del;

         // Extract dimensions from the variable.
         std::copy(var->dims, var->dims + 2, dims.begin());

         if (!Parent::readData(file, var)) goto del;
         _private = new PrivateData(dims, var);
         return;

      del:
         Mat_VarFree(var);
      }
   };

   /**
    * @brief Non-complex MATLAB scalar.
    * @tparam T Element type.
    * @tparam MC matio_classes enum value for given T.
    * @tparam MT matio_types enum value for given T.
    */
   template <typename T, matio_classes MC = MatioCompat<T>::cls, matio_types MT = MatioCompat<T>::typ>
   class MatScalar : public MatVar<T, 1, MC, MT>
   {
   public:
      using Parent = MatVar<T, 2, MC, MT>;

      /**
       * @brief Constructs an invalid MatScalar.
       */
      MatScalar()
      {
      }

      /**
       * @brief Constructs a MatScalar with given name.
       * @param name Name to assign to the variable.
       */
      MatScalar(const std::string &name) :
         Parent(name, 1, 1)
      {
      }

      /**
       * @brief Reads a MatScalar with given name from the given matio file.
       * @param file The matio file to read from.
       * @param name Name of the variable to read.
       * @postcondition The MatScalar is invalid if the variable is not found in the matio file.
       * @postcondition The MatScalar is invalid if the variable's rank is not 2.
       * @postcondition The MatScalar is invalid if the variable's matio_classes class is not MC.
       * @postcondition The MatScalar is invalid if the variable's matio_types type is not MT.
       * @postcondition The MatScalar is invalid if the variable is complex.
       * @postcondition The MatScalar's dimensions are according to the variable stored in the matio file.
       * @postcondition The MatScalar contains the data from the matio file.
       */
      MatScalar(mat_t *file, const std::string &name)
      {
         std::array<SizeType, 2> dims;

         // Read basic variable info from the file.
         matvar_t *var = Parent::readBasic(file, name);
         if (!var) return;

         // Check additional constraints
         if (var->isComplex) goto del;
         if (var->dims[0] != 1) goto del;
         if (var->dims[1] != 1) goto del;

         // Extract dimensions from the variable.
         std::copy(var->dims, var->dims + 2, dims.begin());

         if (!Parent::readData(file, var)) goto del;
         _private = new PrivateData(dims, var);
         return;

      del:
         Mat_VarFree(var);
      }
   };

   /**
    * @brief Multidimensional complex MATLAB variable.
    * @tparam T Element type.
    * @tparam N Number of dimensions.
    * @tparam MC matio_classes enum value for given T.
    * @tparam MT matio_types enum value for given T.
    */
   template <typename T, size_t N, matio_classes MC = MatioCompat<T>::cls, matio_types MT = MatioCompat<T>::typ>
   class MatVarComplex : public MatVarBase<T, N, MC, MT>
   {
   protected:
      using Base = MatVarBase<T, N, MC, MT>; /**< Parameterized type of the MatVarBase base class. */

      /** @internal */
      struct PrivateData : Base::PrivateData
      {
         PrivateData(const std::string& name, const std::array<SizeType, N> &dims) :
            real(new ArrayType(dims)),
            imag(new ArrayType(dims))
         {
            cplx.Re = real->data();
            cplx.Im = imag->data();
            var = Mat_VarCreate(name.c_str(), MC, MT, N, const_cast<SizeType *>(dims.data()), &cplx, MAT_F_COMPLEX | MAT_F_DONT_COPY_DATA);
            assert(var != nullptr);
         }

         // NOTE: The PrivateData takes ownership of the variable.
         PrivateData(const std::array<SizeType, N> &dims, matvar_t *variable) :
            Base::PrivateData(variable),
            real(new ArrayType(dims, castRe(variable->data))),
            imag(new ArrayType(dims, castIm(variable->data)))
         {
            // Here, the variable reserved the memory.
            // Have the Array use it (NOT own it).
         }

         PrivateData(const PrivateData& other)
         {
            if (other.real && other.imag)
            {
               // Duplicate the array and the variable.
               // The variable will access (NOT own) the Array's data.
               real.swap(std::unique_ptr<ArrayType>(new ArrayType(*other.real)));
               imag.swap(std::unique_ptr<ArrayType>(new ArrayType(*other.imag)));

               var = Mat_VarDuplicate(other.var, 0);  // 0 = don't copy data
               assert(var != nullptr);
               
               cplx.Re = real->data();                // Point the complex struct members...
               cplx.Im = imag->data();                // ... to the both arrays

               var->data = &cplx;                     // replace the variable's data pointer but...
               var->mem_conserve = 1;                 // ...tell the variable that memory is NOT owned
            }
         }

         std::unique_ptr<ArrayType> real;
         std::unique_ptr<ArrayType> imag;
         mat_complex_split_t cplx;

      private:
         inline static T* castRe(void *data)
         {
            return static_cast<T*>(static_cast<mat_complex_split_t*>(data)->Re);
         }

         inline static T* castIm(void *data)
         {
            return static_cast<T*>(static_cast<mat_complex_split_t*>(data)->Im);
         }
      };

   public:

      /**
       * @brief Constructs an invalid MatVarComplex.
       */
      MatVarComplex()
      {
         // invalid = no private data
      }

      /**
       * @brief Constructs a MatVarComplex with given name and dimensional extents.
       * @param name Name to assign to the variable.
       * @param dims Number of elements in each of the N dimensions.
       */
      MatVarComplex(const std::string &name, std::array<SizeType, N> dims) :
         Base(new PrivateData(name, dims))
      {
      }

      /**
       * @brief Constructs a MatVarComplex with given name and dimensional extents.
       * @param name Name to assign to the variable.
       * @param dim1 Number of elements in the first dimension.
       * @param rest Number of elements in the other N-1 dimensions.
       */
      template <typename... Rest>
      MatVarComplex(const std::string &name, SizeType dim1, Rest... rest) :
         MatVarComplex(name, { dim1, rest... })
      {
      }

      /**
       * @brief Reads a MatVarComplex with given name from the given matio file.
       * @param file The matio file to read from.
       * @param name Name of the variable to read.
       * @postcondition The MatVarComplex is invalid if the variable is not found in the matio file.
       * @postcondition The MatVarComplex is invalid if the variable's rank is not N.
       * @postcondition The MatVarComplex is invalid if the variable's matio_classes class is not MC.
       * @postcondition The MatVarComplex is invalid if the variable's matio_types type is not MT.
       * @postcondition The MatVarComplex is invalid if the variable is not complex.
       * @postcondition The MatVarComplex's dimensions are according to the variable stored in the matio file.
       * @postcondition The MatVarComplex contains the data from the matio file.
       */
      MatVarComplex(mat_t *file, const std::string &name)
      {
         std::array<SizeType, N> dims;

         // Read basic variable info from the file.
         matvar_t *var = Base::readBasic(file, name);
         if (!var) return;

         // Check additional constraints
         if (!var->isComplex) goto del;

         // Extract dimensions from the variable.
         std::copy(var->dims, var->dims + N, dims.begin());

         if (!Base::readData(file, var)) goto del;
         _private = new PrivateData(dims, var);
         return;

      del:
         Mat_VarFree(var); // This will leave _private uninitialized and thus the object invalid.
      }

      /**
       * @brief Returns a constant reference to the underlying MultiArray of real values.
       * @precondition This object is valid.
       * @see isValid()
       */
      inline const ArrayType& real() const
      {
         return *static_cast<const PrivateData&>(*_private).real;
      }

      /**
       * @brief Returns a mutable reference to the underlying MultiArray of real values.
       * @precondition This object is valid.
       * @see isValid()
       */
      inline ArrayType& real()
      {
         return *static_cast<PrivateData&>(*_private).real;
      }

      /**
       * @brief Returns a constant reference to the underlying MultiArray of imaginary values.
       * @precondition This object is valid.
       * @see isValid()
       */
      inline const ArrayType& imag() const
      {
         return *static_cast<const PrivateData&>(*_private).imag;
      }

      /**
       * @brief Returns a mutable reference to the underlying MultiArray of imaginary values.
       * @precondition This object is valid.
       * @see isValid()
       */
      inline ArrayType& imag()
      {
         return *static_cast<PrivateData&>(*_private).imag;
      }

      /**
       * @brief Returns the number of elements in the given dimension.
       * @param n Zero-based dimension number.
       * @return Number of elements in the n'th dimension.
       * @precondition This object is valid.
       * @precondition n is in bounds (0 <= n < N)
       * @see isValid()
       */
      inline SizeType dim(SizeType n) const
      {
         return real().dim(n);
      }

      /**
      * @brief Returns true if any dimension has zero elements. False otherwise.
      * @param
      * @precondition This object is valid.
      * @see isValid()
      */
      inline bool isEmpty() const
      {
         const std::array<SizeType, N> &dims = array().dims();
         return std::find(dims.cbegin(), dims.cend(), SizeType(0)) != dims.cend();
      }
   };

   /**
    * @brief Complex MATLAB row vector.
    * @tparam T Element type.
    * @tparam MC matio_classes enum value for given T.
    * @tparam MT matio_types enum value for given T.
    */
   template <typename T, matio_classes MC = MatioCompat<T>::cls, matio_types MT = MatioCompat<T>::typ>
   class MatVecComplex : public MatVarComplex<T, 1, MC, MT>
   {
   public:
      using Parent = MatVarComplex<T, 1, MC, MT>;

      /**
       * @brief Constructs an invalid MatVecComplex.
       */
      MatVecComplex()
      {
      }

      /**
       * @brief Constructs a MatVecComplex with given name and length.
       * @param name Name to assign to the variable.
       * @param length Number of elements in the first dimension.
       */
      MatVecComplex(const std::string &name, SizeType length) :
         Parent(name, length)
      {
      }

      /**
       * @brief Reads a MatVecComplex with given name from the given matio file.
       * @param file The matio file to read from.
       * @param name Name of the variable to read.
       * @postcondition The MatVecComplex is invalid if the variable is not found in the matio file.
       * @postcondition The MatVecComplex is invalid if the variable's rank is not 2.
       * @postcondition The MatVecComplex is invalid if the variable's matio_classes class is not MC.
       * @postcondition The MatVecComplex is invalid if the variable's matio_types type is not MT.
       * @postcondition The MatVecComplex is invalid if the variable is not complex.
       * @postcondition The MatVecComplex's dimensions are according to the variable stored in the matio file.
       * @postcondition The MatVecComplex contains the data from the matio file.
       */
      MatVecComplex(mat_t *file, const std::string &name)
      {
         std::array<SizeType, 2> dims;

         // Read basic variable info from the file.
         matvar_t *var = Parent::readBasic(file, name);
         if (!var) return;

         // Check additional constraints
         if (!var->isComplex) goto del;
         if (var->dims[1] != 1) goto del;

         // Extract dimensions from the variable.
         std::copy(var->dims, var->dims + 2, dims.begin());

         if (!Parent::readData(file, var)) goto del;
         _private = new PrivateData(dims, var);
         return;

      del:
         Mat_VarFree(var);
      }
   };

   /**
    * @brief Complex MATLAB scalar.
    * @tparam T Element type.
    * @tparam MC matio_classes enum value for given T.
    * @tparam MT matio_types enum value for given T.
    */
   template <typename T, matio_classes MC = MatioCompat<T>::cls, matio_types MT = MatioCompat<T>::typ>
   class MatScalarComplex : public MatVarComplex<T, 1, MC, MT>
   {
   public:
      using Parent = MatVar<T, 2, MC, MT>;

      /**
       * @brief Constructs an invalid MatScalarComplex.
       */
      MatScalarComplex()
      {
      }

      /**
       * @brief Constructs a MatScalarComplex with given name.
       * @param name Name to assign to the variable.
       */
      MatScalarComplex(const std::string &name) :
         Parent(name, 1, 1)
      {
      }

      /**
       * @brief Reads a MatScalarComplex with given name from the given matio file.
       * @param file The matio file to read from.
       * @param name Name of the variable to read.
       * @postcondition The MatScalarComplex is invalid if the variable is not found in the matio file.
       * @postcondition The MatScalarComplex is invalid if the variable's rank is not 2.
       * @postcondition The MatScalarComplex is invalid if the variable's matio_classes class is not MC.
       * @postcondition The MatScalarComplex is invalid if the variable's matio_types type is not MT.
       * @postcondition The MatScalarComplex is invalid if the variable is not complex.
       * @postcondition The MatScalarComplex's dimensions are according to the variable stored in the matio file.
       * @postcondition The MatScalarComplex contains the data from the matio file.
       */
      MatScalarComplex(mat_t *file, const std::string &name)
      {
         std::array<SizeType, 2> dims;

         // Read basic variable info from the file.
         matvar_t *var = Parent::readBasic(file, name);
         if (!var) return;

         // Check additional constraints
         if (!var->isComplex) goto del;
         if (var->dims[0] != 1) goto del;
         if (var->dims[1] != 1) goto del;

         // Extract dimensions from the variable.
         std::copy(var->dims, var->dims + 2, dims.begin());

         if (!Parent::readData(file, var)) goto del;
         _private = new PrivateData(dims, var);
         return;

      del:
         Mat_VarFree(var);
      }
   };

   /**
    * @brief Construct object and assign to reference.
    * @tparam T Some constructible type.
    * @param var Reference to the value of type T to assign to.
    * @param args Arguments to pass to the constructor of T.
    * @postcondition A new object of type T was constructed with the given arguments.
    * @postcondition The referenced var is assigned to newly created object.
    *
    * Compact notation due to parameter deduction.
    */
   template <typename T, typename... Args>
   inline void Construct(T& var, Args&&... args)
   {
      var = T(args...);
   }

   /**
    * @brief Assigns a reference its default-constructed type.
    * @tparam T Some default-constructible type.
    * @param var Reference to the value of type T to assign to.
    * @postcondition A new object of type T was default-constructed.
    * @postcondition The referenced var is assigned to newly created object.
    *
    * Compact notation due to parameter deduction.
    */
   template <typename T>
   inline void SetToDefault(T& var)
   {
      Construct(var);
   }
}

#endif
