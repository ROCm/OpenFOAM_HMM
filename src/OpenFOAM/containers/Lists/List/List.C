/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2022 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "List.H"
#include "ListLoopM.H"
#include "FixedList.H"
#include "PtrList.H"
#include "SLList.H"
#include "contiguous.H"
#include <utility>

#include <stdlib.h>  //LG1 AMD

#include <type_traits>



#ifdef USE_ROCTX
#include <roctracer/roctx.h>
#endif

#ifdef USE_MEM_POOL
void * provide_umpire_pool(size_t N);
void free_umpire_pool( void * data);
bool is_umpire_pool_ptr(void *ptr);

//#define USE_MEM_POOL
#endif


#ifndef OMP_UNIFIED_MEMORY_REQUIRED
#pragma omp requires unified_shared_memory
#define OMP_UNIFIED_MEMORY_REQUIRED
#endif


//forward declaration of Foam::Vector<>

namespace Foam 
{
template<typename T> class Vector;
template<typename T> class Tensor;
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class T>
void Foam::List<T>::doResize(const label len)
{
    if (len == this->size_)
    {
        return;
    }

    if (len > 0)
    {
        // With sign-check to avoid spurious -Walloc-size-larger-than
        
	
        #ifdef USE_ROCTX
	if (len  > 10000){
	  char roctx_name[128];
	  sprintf(roctx_name,"resizing_%zu",sizeof(T)*len);
          roctxRangePush(roctx_name);
	}
        #endif
        
        T *nv;

	   //fprintf(stderr,"doResize: calling mem allocator len = %d\n",len);

	   #ifdef USE_MEM_POOL

	   if (len > 10000 && is_contiguous<T>::value ){
	       void * tmp_ptr = provide_umpire_pool(sizeof(T)*len);
               nv = new (tmp_ptr) T[len]; //use placement new	    
	   }
	   else {
		size_t alignement = 16;
                size_t bytes_needed = sizeof(T)*len;
                if (bytes_needed > 2*100){ //LG1 AMD
                   alignement = 256;       //LG1 AMD
                }
                nv = new (std::align_val_t( alignement)) T[len];
	   }
           #else
                size_t alignement = 16;
                size_t bytes_needed = sizeof(T)*len;
                if (bytes_needed > 2*100){ //LG1 AMD
                   alignement = 256;       //LG1 AMD
                }
                nv = new (std::align_val_t( alignement)) T[len];

           #endif
	   //if (nv == NULL) fprintf(stderr,"nv is NULL\n");


        #ifdef USE_ROCTX
	if (len  > 10000){
           roctxRangePop();
        }
        #endif


        const label overlap = Foam::min(this->size_, len);

        if (overlap)
        {

            #ifdef USEMEMCPY
            if (is_contiguous<T>::value)
            {
                std::memcpy
                (
                    static_cast<void*>(nv), this->v_, overlap*sizeof(T)
                );
            }
            else
            #endif
            {
                List_ACCESS(T, *this, vp);
                for (label i = 0; i < overlap; ++i)
                {
                    nv[i] = std::move(vp[i]);
                }
            }
        }

        clear();
        this->size_ = len;
        this->v_ = nv;
    }
    else
    {
        // Or only #ifdef FULLDEBUG
        if (len < 0)
        {
            FatalErrorInFunction
                << "bad size " << len
                << abort(FatalError);
        }
        // #endif

        clear();
    }
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class T>
Foam::List<T>::List(const label len)
:
    UList<T>(nullptr, len)
{
    if (len < 0)
    {
        FatalErrorInFunction
            << "bad size " << len
            << abort(FatalError);
    }

    doAlloc();
}


template<class T>
Foam::List<T>::List(const label len, const T& val)
:
    UList<T>(nullptr, len)
{
    if (len < 0)
    {
        FatalErrorInFunction
            << "bad size " << len
            << abort(FatalError);
    }

    if (len)
    {

	#ifdef USE_ROCTX
        roctxRangePush("List::List_ref");
        #endif

    	doAlloc();

        //PRINTS
          // if  ( !(std::is_same<T,scalar>() || std::is_same<T,int>() || std::is_same<T,unsigned int>() || std::is_same_v<T,Foam::Vector<scalar>> || std::is_same_v<T,Foam::Tensor<scalar>> ) && len>10000 ) fprintf(stderr,"List:not scalar/vector/tensor line=%d\n",__LINE__);


	if constexpr ( std::is_same<T,scalar>() ) {
           scalar * __restrict__ vp_ptr = (*this).begin();
           #pragma omp target teams distribute parallel for if(target:len>20000)
           for (label i=0; i < len; ++i)
           {
              vp_ptr[i] = val;
           }
	}
	else if constexpr ( std::is_same<T,int>() ) {
           int * __restrict__ vp_ptr = (*this).begin();
           #pragma omp target teams distribute parallel for if(target:len>20000)
           for (label i=0; i < len; ++i)
           {
              vp_ptr[i] = val;
           }
        }
        else if constexpr ( std::is_same<T,unsigned int>() ) {
           unsigned int * __restrict__ vp_ptr = (*this).begin();
           #pragma omp target teams distribute parallel for if(target:len>20000)
           for (label i=0; i < len; ++i)
           {
              vp_ptr[i] = val;
           }
        }
	else if (std::is_same_v<T,Foam::Vector<scalar>> || std::is_same_v<T,Foam::Tensor<scalar>> ) {
           T * __restrict__ vp_ptr = (*this).begin();
           #pragma omp target teams distribute parallel for if(target:len>20000)
           for (label i=0; i < len; ++i)
           {
              vp_ptr[i] = val;
           }
         }
	else{
//	test_type<T>;

          List_ACCESS(T, (*this), vp);

/*
lnInclude/List.C:254:26: error: cannot pass non-trivial object of type 'Foam::SphericalTensor<double>' to variadic function; expected type from format string was 'unsigned int' [-Wnon-pod-varargs]
lnInclude/List.C:254:26: error: cannot pass non-trivial object of type 'Foam::Vector<double>' to variadic function; expected type from format string was 'unsigned int' [-Wnon-pod-varargs]
lnInclude/List.C:254:26: error: cannot pass non-trivial object of type 'Foam::Vector<double>' to variadic function; expected type from format string was 'unsigned int' [-Wnon-pod-varargs]
lnInclude/List.C:254:26: error: cannot pass non-trivial object of type 'Foam::Tensor<double>' to variadic function; expected type from format string was 'unsigned int' [-Wnon-pod-varargs]
lnInclude/List.C:254:26: error: cannot pass non-trivial object of type 'Foam::Tensor<double>' to variadic function; expected type from format string was 'unsigned int' [-Wnon-pod-varargs]
*/

          for (label i=0; i < len; ++i)
          {
            vp[i] = val;
          }
	}

	#ifdef USE_ROCTX
        roctxRangePop();
        #endif
    }
}


template<class T>
Foam::List<T>::List(const label len, const Foam::zero)
:
    UList<T>(nullptr, len)
{
    if (len < 0)
    {
        FatalErrorInFunction
            << "bad size " << len
            << abort(FatalError);
    }

    if (len)
    {
        doAlloc();


        List_ACCESS(T, (*this), vp);
        #pragma omp target teams distribute parallel for if(target:len>20000)
        for (label i=0; i < len; ++i)
        {
            vp[i] = Zero;
        }
    }
}


template<class T>
Foam::List<T>::List(const Foam::one, const T& val)
:
    UList<T>(new T[1], 1)
{
    this->v_[0] = val;
}


template<class T>
Foam::List<T>::List(const Foam::one, T&& val)
:
    UList<T>(new T[1], 1)
{
    this->v_[0] = std::move(val);
}


template<class T>
Foam::List<T>::List(const Foam::one, const Foam::zero)
:
    UList<T>(new T[1], 1)
{
    this->v_[0] = Zero;
}


template<class T>
Foam::List<T>::List(const UList<T>& a)
:
    UList<T>(nullptr, a.size_)
{
    const label len = this->size_;

    if (len)
    {
        doAlloc();

        #ifdef USEMEMCPY
        if (is_contiguous<T>::value)
        {
            std::memcpy
            (
                static_cast<void*>(this->v_), a.v_, this->size_bytes()
            );
        }
        else
        #endif
        {
            List_ACCESS(T, (*this), vp);
            List_CONST_ACCESS(T, a, ap);
            for (label i = 0; i < len; ++i)
            {
                vp[i] = ap[i];
            }
        }
    }
}


template<class T>
Foam::List<T>::List(const List<T>& a)
:
    UList<T>(nullptr, a.size_)
{
    const label len = this->size_;

    if (len)
    {
        #ifdef USE_ROCTX
        roctxRangePush("List::List");
        #endif

        doAlloc();

        #ifdef USEMEMCPY
        if (is_contiguous<T>::value)
        {
            std::memcpy
            (
                static_cast<void*>(this->v_), a.v_, this->size_bytes()
            );
        }
        else
        #endif
        {
            //PRINTS
           //if  ( !(std::is_same<T,scalar>() || std::is_same<T,int>() || std::is_same<T,unsigned int>() || std::is_same<T,Foam::Vector<scalar>>() ) && len>10000 ) fprintf(stderr,"List:not scalar/Vector line=%d\n",__LINE__);

            if constexpr ( std::is_same<T,scalar>() ) {
              scalar * __restrict__ vp_ptr = (*this).begin();
              const scalar * __restrict__ ap_ptr = a.begin();
              #pragma omp target teams distribute parallel for if(target:len > 10000) 
              for (label i = 0; i < len; ++i)
              {
                vp_ptr[i] = ap_ptr[i];
              }
	    }
            else if constexpr ( std::is_same<T,int>() ) {
              int * __restrict__ vp_ptr = (*this).begin();
              const int * __restrict__ ap_ptr = a.begin();
              #pragma omp target teams distribute parallel for if(target:len > 10000) 
              for (label i = 0; i < len; ++i)
              {
                vp_ptr[i] = ap_ptr[i];
              }
            }
            else if constexpr ( std::is_same<T,unsigned int>() ) {
              unsigned int * __restrict__ vp_ptr = (*this).begin();
              const unsigned int * __restrict__ ap_ptr = a.begin();
              #pragma omp target teams distribute parallel for if(target:len > 10000) 
              for (label i = 0; i < len; ++i)
              {
                vp_ptr[i] = ap_ptr[i];
              }
            }
	    else if constexpr ( std::is_same<T,Foam::Vector<scalar>>() ) {
              T * __restrict__ vp_ptr = (*this).begin();
              const T * __restrict__ ap_ptr = a.begin();
              #pragma omp target teams distribute parallel for if(target:len > 10000)
              for (label i = 0; i < len; ++i)
              {
                vp_ptr[i] = ap_ptr[i];
              }
            }


	    else{
              List_ACCESS(T, (*this), vp);
              List_CONST_ACCESS(T, a, ap);
              for (label i = 0; i < len; ++i)
              {
                vp[i] = ap[i];
              }
	    }
        }
        #ifdef USE_ROCTX
        roctxRangePop();
        #endif
    }
}


template<class T>
Foam::List<T>::List(List<T>& a, bool reuse)
:
    UList<T>(nullptr, a.size_)
{
    if (reuse)
    {
        // Steal content
        this->v_ = a.v_;
        a.v_ = nullptr;
        a.size_ = 0;
        return;
    }

    const label len = this->size_;

    if (len)
    {
        doAlloc();

        #ifdef USEMEMCPY
        if (is_contiguous<T>::value)
        {
            std::memcpy
            (
                static_cast<void*>(this->v_), a.v_, this->size_bytes()
            );
        }
        else
        #endif
        {
	    if constexpr ( std::is_same<T,scalar>() || std::is_same<T,int>() || std::is_same<T,unsigned int>() || std::is_same<T,Foam::Vector<scalar>>() ) {
              T * __restrict__ vp_ptr = (*this).begin();
              const T * __restrict__ ap_ptr = a.begin();
              #pragma omp target teams distribute parallel for if(target:len > 10000)
              for (label i = 0; i < len; ++i)
              {
                vp_ptr[i] = ap_ptr[i];
              }
            }
	    else{
              List_ACCESS(T, (*this), vp);
              List_CONST_ACCESS(T, a, ap);
              for (label i = 0; i < len; ++i)
              {
                  vp[i] = ap[i];
              }
	    }
        }
    }
}


template<class T>
Foam::List<T>::List(const UList<T>& list, const labelUList& indices)
:
    UList<T>(nullptr, indices.size())
{
    const label len = indices.size();

    if (len)
    {
        doAlloc();
        //AMD LG used in snappy Mesh a lot

        List_ACCESS(T, (*this), vp);

        for (label i=0; i < len; ++i)
        {
            vp[i] = list[indices[i]];
        }
    }
}


template<class T>
template<unsigned N>
Foam::List<T>::List
(
    const UList<T>& list,
    const FixedList<label,N>& indices
)
:
    UList<T>(nullptr, label(N))
{
    const label len = label(N);

    doAlloc();

    fprintf(stderr,"line=%d copy\n",__LINE__);


    List_ACCESS(T, (*this), vp);

    for (label i=0; i < len; ++i)
    {
        vp[i] = list[indices[i]];
    }
}


template<class T>
template<unsigned N>
Foam::List<T>::List(const FixedList<T, N>& list)
:
    UList<T>(nullptr, label(N))
{
    doAlloc();
    copyList(list);
}


template<class T>
Foam::List<T>::List(const PtrList<T>& list)
:
    UList<T>(nullptr, list.size())
{
    doAlloc();
    copyList(list);
}


template<class T>
Foam::List<T>::List(const SLList<T>& list)
:
    List<T>(list.begin(), list.end(), list.size())
{}


template<class T>
template<class Addr>
Foam::List<T>::List(const IndirectListBase<T, Addr>& list)
:
    UList<T>(nullptr, list.size())
{
    doAlloc();
    copyList(list);
}


template<class T>
Foam::List<T>::List(std::initializer_list<T> list)
:
    List<T>(list.begin(), list.end(), list.size())
{}


template<class T>
Foam::List<T>::List(List<T>&& list)
:
    UList<T>()
{
    // Can use transfer or swap to manage content
    transfer(list);
}


template<class T>
template<int SizeMin>
Foam::List<T>::List(DynamicList<T, SizeMin>&& list)
:
    UList<T>()
{
    transfer(list);
}


template<class T>
Foam::List<T>::List(SLList<T>&& list)
:
    UList<T>()
{
    operator=(std::move(list));
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

template<class T>
Foam::List<T>::~List()
{
    if (this->v_)
    {

        #ifdef USE_ROCTX
        if (this->size_  > 10000){
          char roctx_name[128];
          sprintf(roctx_name,"deleting_%zu",sizeof(T)*this->size_);
          roctxRangePush(roctx_name);
	}
        #endif

	#ifdef USE_MEM_POOL
        bool umpr_ptr =  is_umpire_pool_ptr( reinterpret_cast<void*>(this->v_) );
        if (umpr_ptr == true  ){
           free_umpire_pool( reinterpret_cast<void*>(this->v_) ); //AMD
        }
        else    
        #endif
           delete[] this->v_;


	#ifdef USE_ROCTX
	if (this->size_  > 10000){
           roctxRangePop();
	}
        #endif

    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
void Foam::List<T>::resize(const label len, const T& val)
{
    label idx = this->size_;
    this->doResize(len);

    List_ACCESS(T, *this, vp);
    while (idx < len)
    {
        vp[idx] = val;
        ++idx;
    }
}


template<class T>
void Foam::List<T>::transfer(List<T>& list)
{
    if (this == &list)
    {
        return;  // Self-assignment is a no-op
    }

    // Clear and swap - could also check for self assignment
    //fprintf(stderr,"calling clear()\n");
    clear();
    //fprintf(stderr,"after clear()\n");

    this->size_ = list.size_;
    this->v_ = list.v_;

    list.size_ = 0;
    list.v_ = nullptr;
}


template<class T>
template<int SizeMin>
void Foam::List<T>::transfer(DynamicList<T, SizeMin>& list)
{
    // Shrink the allocated space to the number of elements used
    list.shrink();
    transfer(static_cast<List<T>&>(list));

    // Ensure DynamicList has proper capacity=0 too
    list.clearStorage();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
void Foam::List<T>::operator=(const UList<T>& a)
{
    if (this == &a)
    {
        return;  // Self-assignment is a no-op
    }

    reAlloc(a.size_);

    const label len = this->size_;

    if (len)
    {
        //fprintf(stderr,"line=%d copy\n",__LINE__);
        #ifdef USE_ROCTX
        roctxRangePush("List_equal_A");
        #endif

        #ifdef USEMEMCPY
        if (is_contiguous<T>::value)
        {
            std::memcpy
            (
                static_cast<void*>(this->v_), a.v_, this->size_bytes()
            );
        }
        else
        #endif
        {
           //PRINTS
//           if  ( !(std::is_same<T,scalar>() || std::is_same<T,int>() || std::is_same<T,unsigned int>() || std::is_same<T,Foam::Vector<scalar>>() || std::is_same<T,Foam::Tensor<scalar>>() ) && len>10000 ) 
//		   fprintf(stderr,"List:not scalar or vector<scalar>, its %s  line=%d\n", typeid(T).name(),  __LINE__);


           if constexpr ( std::is_same<T,scalar>() ) {
              scalar * __restrict__ vp_ptr = (*this).begin();
              const scalar * __restrict__ ap_ptr = a.begin();
              #pragma omp target teams distribute parallel for if(target:len > 10000) 
              for (label i = 0; i < len; ++i)
              {
                vp_ptr[i] = ap_ptr[i];
              }
	   } 
	   else if constexpr ( std::is_same<T,int>() ) {
              int * __restrict__ vp_ptr = (*this).begin();
              const int * __restrict__ ap_ptr = a.begin();
              #pragma omp target teams distribute parallel for if(target:len > 10000)
              for (label i = 0; i < len; ++i)
              {
                vp_ptr[i] = ap_ptr[i];
              }
           }
	   else if constexpr ( std::is_same<T,unsigned int>() ) {
              unsigned int * __restrict__ vp_ptr = (*this).begin();
              const unsigned int * __restrict__ ap_ptr = a.begin();
              #pragma omp target teams distribute parallel for if(target:len > 10000)
              for (label i = 0; i < len; ++i)
              {
                vp_ptr[i] = ap_ptr[i];
              }
           }
	   else if constexpr ( std::is_same<T,Foam::Vector<scalar>>() || std::is_same<T,Foam::Tensor<scalar>>() ) {
              T * __restrict__ vp_ptr = (*this).begin();
              const T * __restrict__ ap_ptr = a.begin();
              #pragma omp target teams distribute parallel for if(target:len > 10000)
              for (label i = 0; i < len; ++i)
              {
                vp_ptr[i] = ap_ptr[i];
              }
           }
	   else{
              List_ACCESS(T, (*this), vp);
              List_CONST_ACCESS(T, a, ap);
              for (label i = 0; i < len; ++i)
              {
                vp[i] = ap[i];
              }
	   }
        }
	#ifdef USE_ROCTX
        roctxRangePop();
        #endif
    }
}


template<class T>
void Foam::List<T>::operator=(const List<T>& list)
{
    if (this == &list)
    {
        return;  // Self-assignment is a no-op
    }

    operator=(static_cast<const UList<T>&>(list));
}


template<class T>
void Foam::List<T>::operator=(const SLList<T>& list)
{
    const label len = list.size();

    reAlloc(len);

    if (len)
    {
        T* iter = this->begin();

        for (const T& val : list)
        {
            *iter = val;
            ++iter;
        }
    }
}


template<class T>
template<unsigned N>
void Foam::List<T>::operator=(const FixedList<T, N>& list)
{
    reAlloc(static_cast<label>(N));

    T* iter = this->begin();

    for (const T& val : list)
    {
        *iter = val;
        ++iter;
    }
}


template<class T>
template<class Addr>
void Foam::List<T>::operator=(const IndirectListBase<T, Addr>& list)
{
    const label len = list.size();

    reAlloc(len);

    if (len)
    {
        List_ACCESS(T, (*this), vp);

        for (label i=0; i < len; ++i)
        {
            vp[i] = list[i];
        }
    }
}


template<class T>
void Foam::List<T>::operator=(std::initializer_list<T> list)
{
    const label len = list.size();

    reAlloc(len);

    if (len)
    {
        T* iter = this->begin();

        for (const T& val : list)
        {
            *iter = val;
            ++iter;
        }
    }
}


template<class T>
void Foam::List<T>::operator=(List<T>&& list)
{
    if (this == &list)
    {
        return;  // Self-assignment is a no-op
    }

    transfer(list);
}


template<class T>
template<int SizeMin>
void Foam::List<T>::operator=(DynamicList<T, SizeMin>&& list)
{
    transfer(list);
}


template<class T>
void Foam::List<T>::operator=(SLList<T>&& list)
{
    label len = list.size();

    reAlloc(len);

    for (T* iter = this->begin(); len--; ++iter)
    {
        *iter = std::move(list.removeHead());
    }

    list.clear();
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class T>
Foam::labelList Foam::sortedOrder
(
    const UList<T>& list
)
{
    labelList order;
    Foam::sortedOrder(list, order, typename UList<T>::less(list));
    return order;
}


template<class T>
void Foam::sortedOrder
(
    const UList<T>& list,
    labelList& order
)
{
    Foam::sortedOrder(list, order, typename UList<T>::less(list));
}


template<class T, class ListComparePredicate>
void Foam::sortedOrder
(
    const UList<T>& list,
    labelList& order,
    const ListComparePredicate& comp
)
{
    // List lengths must be identical. Old content is overwritten
    order.resize_nocopy(list.size());

    // Same as std::iota and ListOps::identity
    label value = 0;
    for (label& item : order)
    {
        item = value;
        ++value;
    }

    Foam::stableSort(order, comp);
}


// ************************************************************************* //
