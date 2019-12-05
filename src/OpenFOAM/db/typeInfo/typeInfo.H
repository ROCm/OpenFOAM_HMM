/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2019 OpenCFD Ltd.
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

Typedef
    Foam::typeInfo

Description
    Basic run-time type information using word as the type's name.
    Used to enhance the standard RTTI to cover I/O.

    The user can get the type's type name using the type info access function
    \code
        type()
    \endcode

    The reference type cast template function:
    \code
        refCast<T>(r)
    \endcode

    wraps dynamic_cast to handle the bad_cast exception and generate a
    FatalError.

    The isA function:
    \code
        isA<T>(r)
    \endcode

    returns true if r is of type T or derived from type T.

\*---------------------------------------------------------------------------*/

#ifndef typeInfo_H
#define typeInfo_H

#include "error.H"
#include "className.H"
#include <typeinfo>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Declarations (for use in header files)

//- Declare a ClassNameNoDebug() with extra virtual type info
#define TypeNameNoDebug(TypeNameString)                                        \
    ClassNameNoDebug(TypeNameString);                                          \
    virtual const word& type() const { return typeName; }

//- Declare a ClassName() with extra virtual type info
#define TypeName(TypeNameString)                                               \
    ClassName(TypeNameString);                                                 \
    virtual const word& type() const { return typeName; }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

//- Reference type cast template function,
//- wraps dynamic_cast to handle bad_cast exception and generate a FatalError.
template<class To, class From>
inline To& dynamicCast(From& r)
{
    try
    {
        return dynamic_cast<To&>(r);
    }
    catch (const std::bad_cast&)
    {
        FatalErrorInFunction
            << "Attempt to cast type " << typeid(r).name()
            << " to type " << typeid(To).name()
            << abort(FatalError);

        return dynamic_cast<To&>(r);
    }
}


//- Reference type cast template function,
//- wraps dynamic_cast to handle bad_cast exception and generate a FatalError.
template<class To, class From>
inline To& dynamicCast(From& r, const dictionary& d)
{
    try
    {
        return dynamic_cast<To&>(r);
    }
    catch (const std::bad_cast&)
    {
        FatalIOErrorInFunction(d)
            << "Attempt to cast type " << typeid(r).name()
            << " to type " << typeid(To).name()
            << abort(FatalIOError);

        return dynamic_cast<To&>(r);
    }
}


//- Reference type cast template function.
//  As per dynamicCast, but handles type names for the error messages
//  via the virtual type() method.
template<class To, class From>
inline To& refCast(From& r)
{
    try
    {
        return dynamic_cast<To&>(r);
    }
    catch (const std::bad_cast&)
    {
        FatalErrorInFunction
            << "Attempt to cast type " << r.type()
            << " to type " << To::typeName
            << abort(FatalError);

        return dynamic_cast<To&>(r);
    }
}


//- Reference type cast template function.
//  As per dynamicCast, but handles type names for the error messages
//  via the virtual type() method.
//  Can use index to convey the context.
template<class To, class From>
inline To& refCast(From& r, const label index)
{
    try
    {
        return dynamic_cast<To&>(r);
    }
    catch (const std::bad_cast&)
    {
        FatalErrorInFunction
            << "Attempt to cast type " << r.type()
            << " to type " << To::typeName
            << " at index " << index
            << abort(FatalError);

        return dynamic_cast<To&>(r);
    }
}


//- Reference type cast template function.
//  As per dynamicCast, but handles type names for the error messages
//  via the virtual type() method.
template<class To, class From>
inline To& refCast(From& r, const dictionary& d)
{
    try
    {
        return dynamic_cast<To&>(r);
    }
    catch (const std::bad_cast&)
    {
        FatalIOErrorInFunction(d)
            << "Attempt to cast type " << r.type()
            << " to type " << To::typeName
            << abort(FatalIOError);

        return dynamic_cast<To&>(r);
    }
}


//- Check if dynamic_cast to TargetType is possible
template<class TargetType, class Type>
inline const TargetType* isA(const Type& t)
{
    const Type* p = &t;
    return dynamic_cast<const TargetType*>(p);
}


//- Check if dynamic_cast to TargetType is possible, as a functor
template<class TargetType>
struct isAOp
{
    template<class Type>
    inline bool operator()(const Type& t) const
    {
        return isA<TargetType,Type>(t);
    }
};


//- Check is typeid is identical to the TargetType
template<class TargetType, class Type>
inline bool isType(const Type& t)
{
    return typeid(TargetType) == typeid(t);
}


//- Check is typeid is identical to the TargetType, as a functor
template<class TargetType>
struct isTypeOp
{
    template<class Type>
    inline bool operator()(const Type& t) const
    {
        return isType<TargetType,Type>(t);
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
