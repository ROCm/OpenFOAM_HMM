/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "symmTransformField.H"
#include "FieldM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class Type>
void transform
(
    Field<Type>& rtf,
    const symmTensorField& trf,
    const Field<Type>& tf
)
{
    if (trf.size() == 1)
    {
        return transform(rtf, trf[0], tf);
    }
    else
    {
        TFOR_ALL_F_OP_FUNC_F_F
        (
            Type, rtf, =, transform, symmTensor, trf, Type, tf
        )
    }
}


template<class Type>
tmp<Field<Type>> transform
(
    const symmTensorField& trf,
    const Field<Type>& tf
)
{
    auto tresult = tmp<Field<Type>>::New(tf.size());
    transform(tresult.ref(), trf, tf);
    return tresult;
}


template<class Type>
tmp<Field<Type>> transform
(
    const symmTensorField& trf,
    const tmp<Field<Type>>& ttf
)
{
    tmp<Field<Type>> tresult = New(ttf);
    transform(tresult.ref(), trf, ttf());
    ttf.clear();
    return tresult;
}


template<class Type>
tmp<Field<Type>> transform
(
    const tmp<symmTensorField>& ttrf,
    const Field<Type>& tf
)
{
    auto tresult = tmp<Field<Type>>::New(tf.size());
    transform(tresult.ref(), ttrf(), tf);
    ttrf.clear();
    return tresult;
}


template<class Type>
tmp<Field<Type>> transform
(
    const tmp<symmTensorField>& ttrf,
    const tmp<Field<Type>>& ttf
)
{
    tmp<Field<Type>> tresult = New(ttf);
    transform(tresult.ref(), ttrf(), ttf());
    ttf.clear();
    ttrf.clear();
    return tresult;
}


template<class Type>
void transform
(
    Field<Type>& rtf,
    const symmTensor& t,
    const Field<Type>& tf
)
{
    TFOR_ALL_F_OP_FUNC_S_F(Type, rtf, =, transform, tensor, t, Type, tf)
}


template<class Type>
tmp<Field<Type>> transform
(
    const symmTensor& t,
    const Field<Type>& tf
)
{
    auto tresult = tmp<Field<Type>>::New(tf.size());
    transform(tresult.ref(), t, tf);
    return tresult;
}


template<class Type>
tmp<Field<Type>> transform
(
    const symmTensor& t,
    const tmp<Field<Type>>& ttf
)
{
    tmp<Field<Type>> tresult = New(ttf);
    transform(tresult.ref(), t, ttf());
    ttf.clear();
    return tresult;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
