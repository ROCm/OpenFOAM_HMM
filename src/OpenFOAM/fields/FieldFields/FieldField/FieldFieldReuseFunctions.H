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

\*---------------------------------------------------------------------------*/

#ifndef FieldFieldReuseFunctions_H
#define FieldFieldReuseFunctions_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<template<class> class Field, class TypeR, class Type1>
struct reuseTmpFieldField
{
    static tmp<FieldField<Field, TypeR>> New
    (
        const tmp<FieldField<Field, Type1>>& tf1
    )
    {
        return FieldField<Field, TypeR>::NewCalculatedType(tf1());
    }
};


template<template<class> class Field, class TypeR>
struct reuseTmpFieldField<Field, TypeR, TypeR>
{
    static tmp<FieldField<Field, TypeR>> New
    (
        const tmp<FieldField<Field, TypeR>>& tf1,
        const bool initCopy = false
    )
    {
        if (tf1.isTmp())
        {
            return tf1;
        }

        auto rtf = FieldField<Field, TypeR>::NewCalculatedType(tf1());

        if (initCopy)
        {
            rtf.ref() = tf1();
        }

        return rtf;
    }
};


//- Global function forwards to reuseTmpFieldField::New
template<template<class> class Field, class TypeR>
tmp<FieldField<Field, TypeR>> New
(
    const tmp<FieldField<Field, TypeR>>& tf1,
    const bool initCopy = false
)
{
    return reuseTmpFieldField<Field, TypeR, TypeR>::New(tf1, initCopy);
}


template
<
    template<class> class Field,
    class TypeR,
    class Type1,
    class Type12,
    class Type2
>
struct reuseTmpTmpFieldField
{
    static tmp<FieldField<Field, TypeR>> New
    (
        const tmp<FieldField<Field, Type1>>& tf1,
        const tmp<FieldField<Field, Type2>>& tf2
    )
    {
        return FieldField<Field, TypeR>::NewCalculatedType(tf1());
    }
};


template<template<class> class Field, class TypeR, class Type1, class Type12>
struct reuseTmpTmpFieldField<Field, TypeR, Type1, Type12, TypeR>
{
    static tmp<FieldField<Field, TypeR>> New
    (
        const tmp<FieldField<Field, Type1>>& tf1,
        const tmp<FieldField<Field, TypeR>>& tf2
    )
    {
        if (tf2.isTmp())
        {
            return tf2;
        }

        return FieldField<Field, TypeR>::NewCalculatedType(tf1());
    }
};


template<template<class> class Field, class TypeR, class Type2>
struct reuseTmpTmpFieldField<Field, TypeR, TypeR, TypeR, Type2>
{
    static tmp<FieldField<Field, TypeR>> New
    (
        const tmp<FieldField<Field, TypeR>>& tf1,
        const tmp<FieldField<Field, Type2>>& tf2
    )
    {
        if (tf1.isTmp())
        {
            return tf1;
        }

        return FieldField<Field, TypeR>::NewCalculatedType(tf1());
    }
};


template<template<class> class Field, class TypeR>
struct reuseTmpTmpFieldField<Field, TypeR, TypeR, TypeR, TypeR>
{
    static tmp<FieldField<Field, TypeR>> New
    (
        const tmp<FieldField<Field, TypeR>>& tf1,
        const tmp<FieldField<Field, TypeR>>& tf2
    )
    {
        if (tf1.isTmp())
        {
            return tf1;
        }
        else if (tf2.isTmp())
        {
            return tf2;
        }

        return FieldField<Field, TypeR>::NewCalculatedType(tf1());
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
