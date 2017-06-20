/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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

#include "axesRotation.H"
#include "dictionary.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(axesRotation, 0);
    addToRunTimeSelectionTable
    (
        coordinateRotation,
        axesRotation,
        dictionary
    );
    addToRunTimeSelectionTable
    (
        coordinateRotation,
        axesRotation,
        objectRegistry
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::axesRotation::axesRotation()
:
    R_(sphericalTensor::I),
    Rtr_(sphericalTensor::I)
{}


Foam::axesRotation::axesRotation(const axesRotation& r)
:
    R_(r.R_),
    Rtr_(r.Rtr_)
{}


Foam::axesRotation::axesRotation(const tensor& R)
:
    R_(R),
    Rtr_(R_.T())
{}


Foam::axesRotation::axesRotation
(
    const vector& axis,
    const vector& dir,
    const axisOrder& order
)
:
    R_(sphericalTensor::I),
    Rtr_(sphericalTensor::I)
{
    setTransform(axis, dir, order);
}


Foam::axesRotation::axesRotation
(
    const vector& axis
)
:
    R_(sphericalTensor::I),
    Rtr_(sphericalTensor::I)
{
    direction maxCmpt = 0, dirCmpt = 1;

    scalar maxVal = mag(axis[maxCmpt]);
    bool negative = (axis[maxCmpt] < 0);

    for (direction cmpt = 1; cmpt < vector::nComponents; ++cmpt)
    {
        const scalar val = mag(axis[cmpt]);

        if (maxVal < val)
        {
            maxVal  = val;
            maxCmpt = cmpt;
            dirCmpt = maxCmpt+1;
            negative = (axis[cmpt] < 0);

            if (dirCmpt >= vector::nComponents)
            {
                dirCmpt = 0;
            }
        }
    }

    vector dir = Zero;
    dir.component(dirCmpt) = (negative ? -1 : 1);

    setTransform(axis, dir, E3_E1);
}


Foam::axesRotation::axesRotation
(
    const dictionary& dict
)
:
    R_(sphericalTensor::I),
    Rtr_(sphericalTensor::I)
{
    operator=(dict);
}


Foam::axesRotation::axesRotation
(
    const dictionary& dict,
    const objectRegistry&
)
:
    axesRotation(dict)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::axesRotation::setTransform
(
    const vector& axis1,
    const vector& axis2,
    const axisOrder& order
)
{
    const vector a = axis1/mag(axis1);
    vector b = axis2;

    b = b - (b & a)*a;

    if (mag(b) < SMALL)
    {
        FatalErrorInFunction
            << "axis1, axis2 appear to be co-linear: "
            << axis1 << ", " << axis2 << endl
            << abort(FatalError);
    }

    b = b/mag(b);
    const vector c = a^b;

    // Global->local transformation
    switch (order)
    {
        case E1_E2:
        {
            Rtr_ = tensor(a, b, c);
            break;
        }
        case E2_E3:
        {
            Rtr_ = tensor(c, a, b);
            break;
        }
        case E3_E1:
        {
            Rtr_ = tensor(b, c, a);
            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unhandled axes specification" << endl
                << abort(FatalError);

            break;
        }
    }

    // Local->global transformation
    R_ = Rtr_.T();
}


const Foam::tensorField& Foam::axesRotation::Tr() const
{
    NotImplemented;
    return NullObjectRef<tensorField>();
}


Foam::tmp<Foam::vectorField> Foam::axesRotation::transform
(
    const vectorField& st
) const
{
    return (R_ & st);
}


Foam::vector Foam::axesRotation::transform(const vector& st) const
{
    return (R_ & st);
}


Foam::tmp<Foam::vectorField> Foam::axesRotation::invTransform
(
    const vectorField& st
) const
{
    return (Rtr_ & st);
}


Foam::vector Foam::axesRotation::invTransform(const vector& st) const
{
    return (Rtr_ & st);
}


Foam::tmp<Foam::tensorField> Foam::axesRotation::transformTensor
(
    const tensorField& st
) const
{
    NotImplemented;
    return tmp<tensorField>(nullptr);
}


Foam::tensor Foam::axesRotation::transformTensor
(
    const tensor& st
) const
{
    return (R_ & st & Rtr_);
}


Foam::tmp<Foam::tensorField> Foam::axesRotation::transformTensor
(
    const tensorField& st,
    const labelList& cellMap
) const
{
    NotImplemented;
    return tmp<tensorField>(nullptr);
}


Foam::tmp<Foam::symmTensorField> Foam::axesRotation::transformVector
(
    const vectorField& st
) const
{
    tmp<symmTensorField> tfld(new symmTensorField(st.size()));
    symmTensorField& fld = tfld.ref();

    forAll(fld, i)
    {
        fld[i] = transformPrincipal(R_, st[i]);
    }
    return tfld;
}


Foam::symmTensor Foam::axesRotation::transformVector
(
    const vector& st
) const
{
    return transformPrincipal(R_, st);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::axesRotation::operator=(const dictionary& dict)
{
    vector axis1, axis2;

    if (dict.readIfPresent("e1", axis1) && dict.readIfPresent("e2", axis2))
    {
        setTransform(axis1, axis2, E1_E2);
    }
    else if (dict.readIfPresent("e2", axis1) && dict.readIfPresent("e3", axis2))
    {
        setTransform(axis1, axis2, E2_E3);
    }
    else if (dict.readIfPresent("e3", axis1) && dict.readIfPresent("e1", axis2))
    {
        setTransform(axis1, axis2, E3_E1);
    }
    else if (dict.found("axis") || dict.found("direction"))
    {
        // Both "axis" and "direction" are required
        // If one is missing the appropriate error message will be generated
        dict.lookup("axis") >> axis1;
        dict.lookup("direction") >> axis2;

        setTransform(axis1, axis2, E3_E1);
    }
    else
    {
        FatalErrorInFunction
            << "not entry of the type (e1, e2) or (e2, e3) or (e3, e1) "
            << "found "
            << exit(FatalError);
    }
}


// ************************************************************************* //
