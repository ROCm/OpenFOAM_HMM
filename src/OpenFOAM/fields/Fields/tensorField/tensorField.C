/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "tensorField.H"
#include "transformField.H"

#define TEMPLATE
#include "FieldFunctionsM.C"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

UNARY_FUNCTION(scalar, tensor, tr)
UNARY_FUNCTION(sphericalTensor, tensor, sph)
UNARY_FUNCTION(symmTensor, tensor, symm)
UNARY_FUNCTION(symmTensor, tensor, twoSymm)
UNARY_FUNCTION(tensor, tensor, skew)
UNARY_FUNCTION(tensor, tensor, dev)
UNARY_FUNCTION(tensor, tensor, dev2)
UNARY_FUNCTION(scalar, tensor, det)
UNARY_FUNCTION(tensor, tensor, cof)

void inv(Field<tensor>& tf, const UList<tensor>& tf1)
{
    if (tf.empty())
    {
        return;
    }

    scalar scale = magSqr(tf1[0]);
    Vector<bool> removeCmpts
    (
        magSqr(tf1[0].xx())/scale < SMALL,
        magSqr(tf1[0].yy())/scale < SMALL,
        magSqr(tf1[0].zz())/scale < SMALL
    );

    if (removeCmpts.x() || removeCmpts.y() || removeCmpts.z())
    {
        tensorField tf1Plus(tf1);

        if (removeCmpts.x())
        {
            tf1Plus += tensor(1,0,0,0,0,0,0,0,0);
        }

        if (removeCmpts.y())
        {
            tf1Plus += tensor(0,0,0,0,1,0,0,0,0);
        }

        if (removeCmpts.z())
        {
            tf1Plus += tensor(0,0,0,0,0,0,0,0,1);
        }

        TFOR_ALL_F_OP_FUNC_F(tensor, tf, =, inv, tensor, tf1Plus)

        if (removeCmpts.x())
        {
            tf -= tensor(1,0,0,0,0,0,0,0,0);
        }

        if (removeCmpts.y())
        {
            tf -= tensor(0,0,0,0,1,0,0,0,0);
        }

        if (removeCmpts.z())
        {
            tf -= tensor(0,0,0,0,0,0,0,0,1);
        }
    }
    else
    {
        TFOR_ALL_F_OP_FUNC_F(tensor, tf, =, inv, tensor, tf1)
    }
}

tmp<tensorField> inv(const UList<tensor>& tf)
{
    auto tres = tmp<tensorField>::New(tf.size());
    inv(tres.ref(), tf);
    return tres;
}

tmp<tensorField> inv(const tmp<tensorField>& tf)
{
    auto tres = New(tf);
    inv(tres.ref(), tf());
    tf.clear();
    return tres;
}

UNARY_FUNCTION(vector, symmTensor, eigenValues)
UNARY_FUNCTION(tensor, symmTensor, eigenVectors)


template<>
tmp<Field<tensor>> transformFieldMask<tensor>
(
    const symmTensorField& stf
)
{
    auto tres = tmp<tensorField>::New(stf.size());
    auto& res = tres.ref();
    TFOR_ALL_F_OP_F(tensor, res, =, symmTensor, stf)
    return tres;
}

template<>
tmp<Field<tensor>> transformFieldMask<tensor>
(
    const tmp<symmTensorField>& tstf
)
{
    tmp<Field<tensor>> ret = transformFieldMask<tensor>(tstf());
    tstf.clear();
    return ret;
}


// * * * * * * * * * * * * * * * global operators  * * * * * * * * * * * * * //

UNARY_OPERATOR(vector, tensor, *, hdual)
UNARY_OPERATOR(tensor, vector, *, hdual)

BINARY_OPERATOR(vector, vector, tensor, /, divide)
BINARY_TYPE_OPERATOR(vector, vector, tensor, /, divide)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "undefFieldFunctionsM.H"

// ************************************************************************* //
