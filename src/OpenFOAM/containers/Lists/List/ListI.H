/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2021 OpenCFD Ltd.
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


#ifdef USE_ROCTX
#include <roctracer/roctx.h>
#endif

#ifdef USE_MEM_POOL
void * provide_umpire_pool(size_t N);
void free_umpire_pool( void * data);
bool is_umpire_pool_ptr(void *ptr);

//#define USE_MEM_POOL
#endif

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class T>
inline void Foam::List<T>::doAlloc()
{
    if (this->size_ > 0)
    {
        // With sign-check to avoid spurious -Walloc-size-larger-than
        #if 1
          size_t alignement = 16;
          size_t bytes_needed = sizeof(T)*this->size_;
          if (bytes_needed > 2*100){ //LG1 AMD
           alignement = 256;       //LG1 AMD        
          }
	 
//	   fprintf(stderr,"doAlloc: calling mem allocator size_ = %d\n",this->size_);

           #ifdef USE_MEM_POOL
           if (this->size_ >= 5000 && is_contiguous<T>::value){
	       void * tmp_ptr = provide_umpire_pool(sizeof(T)*this->size_);
               this->v_ = new (tmp_ptr) T[this->size_]; //use placement new 
	   }
           else
               this->v_ = new (std::align_val_t( alignement)) T[this->size_];
           #else
               this->v_ = new (std::align_val_t( alignement)) T[this->size_];
           #endif


        #else
          this->v_ = new T[this->size_];
        #endif

    }
}


template<class T>
inline void Foam::List<T>::reAlloc(const label len)
{
    if (this->size_ != len)
    {
        clear();
        this->size_ = len;
        doAlloc();
    }
}


template<class T>
template<class List2>
inline void Foam::List<T>::copyList(const List2& list)
{
    const label len = this->size_;


    for (label i=0; i<len; ++i)
    {
        this->operator[](i) = list[i];
    }
}


template<class T>
template<class InputIterator>
inline Foam::List<T>::List
(
    InputIterator begIter,
    InputIterator endIter,
    const label len
)
:
    UList<T>(nullptr, len)
{
    if (this->size_)
    {
        doAlloc();

        #ifdef USE_ROCTX
        if (len  > 10000){
         char roctx_name[128];
         sprintf(roctx_name,"copy_iterator_%zu",sizeof(T)*len);
         roctxRangePush(roctx_name);
        }
        #endif

        InputIterator iter = begIter;
        for (label i = 0; i < len; ++i)
        {
            this->operator[](i) = *iter;
            ++iter;
        }

        #ifdef USE_ROCTX
        if (len  > 10000){
           roctxRangePop();
        }
        #endif


    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
inline constexpr Foam::List<T>::List() noexcept
{}


template<class T>
inline Foam::autoPtr<Foam::List<T>> Foam::List<T>::clone() const
{
    return autoPtr<List<T>>::New(*this);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline const Foam::List<T>& Foam::List<T>::null()
{
    return NullObjectRef<List<T>>();
}


template<class T>
inline void Foam::List<T>::clear()
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


        if (umpr_ptr == true) //AMD
           free_umpire_pool( reinterpret_cast<void*> (this->v_) );
        else    
           delete[] this->v_;

        #else
	    delete[] this->v_;
        #endif


        this->v_ = nullptr;

        #ifdef USE_ROCTX
        if (this->size_  > 10000){
           roctxRangePop(); 
	}
        #endif

    }
    this->size_ = 0;

}


namespace Foam
{
    // Template specialization for bool. Fills new entries with false
    template<>
    inline void List<bool>::resize(const label newLen)
    {
        this->resize(newLen, false);
    }
}


template<class T>
inline void Foam::List<T>::resize(const label len)
{
    this->doResize(len);
}


template<class T>
inline void Foam::List<T>::resize_nocopy(const label len)
{
    this->reAlloc(len);
}


template<class T>
inline T& Foam::List<T>::newElmt(const label i)
{
    label n = this->size();

    if (i >= n)
    {
        if (!n) n = 1;  // Avoid dead-lock when starting from zero-sized

        do
        {
            n *= 2;
        }
        while (i >= n);

        resize(n);
    }

    return UList<T>::operator[](i);
}


template<class T>
inline void Foam::List<T>::append(const T& val)
{
    const label idx = this->size();
    resize(idx + 1);

    this->operator[](idx) = val;  // copy element
}


template<class T>
inline void Foam::List<T>::append(T&& val)
{
    const label idx = this->size();
    resize(idx + 1);

    this->operator[](idx) = std::move(val);  // move assign element
}


template<class T>
inline void Foam::List<T>::append(const UList<T>& list)
{
    if (this == &list)
    {
        FatalErrorInFunction
            << "Attempted appending to self" << abort(FatalError);
    }

    label idx = this->size();
    const label n = list.size();

    resize(idx + n);

    #ifdef USE_ROCTX
    if (n  > 10000){
         char roctx_name[128];
         sprintf(roctx_name,"append1_%zu",sizeof(T)*n);
         roctxRangePush(roctx_name);
    }
    #endif

    for (label i=0; i<n; ++i)
    {
        this->operator[](idx++) = list[i];  // copy element
    }

    #ifdef USE_ROCTX
    if (n  > 10000){
           roctxRangePop();
    }
    #endif

}


template<class T>
template<class Addr>
inline void Foam::List<T>::append(const IndirectListBase<T, Addr>& list)
{
    label idx = this->size();
    const label n = list.size();

    resize(idx + n);

    #ifdef USE_ROCTX
    if (n  > 10000){
         char roctx_name[128];
         sprintf(roctx_name,"append2_%zu",sizeof(T)*n);
         roctxRangePush(roctx_name);
    }
    #endif

    for (label i=0; i<n; ++i)
    {
        this->operator[](idx++) = list[i];  // copy element
    }

    #ifdef USE_ROCTX
    if (n  > 10000){
           roctxRangePop();
    }
    #endif

}


template<class T>
inline Foam::label Foam::List<T>::appendUniq(const T& val)
{
    if (this->found(val))
    {
        return 0;
    }
    else
    {
        this->append(val);
        return 1;  // Increased list length by one
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
inline void Foam::List<T>::operator=(const T& val)
{
    UList<T>::operator=(val);
}


template<class T>
inline void Foam::List<T>::operator=(const Foam::zero)
{
    UList<T>::operator=(Zero);
}


// ************************************************************************* //
