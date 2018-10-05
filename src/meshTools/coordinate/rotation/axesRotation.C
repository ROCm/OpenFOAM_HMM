/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017-2018 OpenCFD Ltd.
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


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::tensor Foam::axesRotation::rotation
(
    const vector& axis1,
    const vector& axis2,
    axisOrder order
)
{
    const scalar magAxis1(mag(axis1));
    scalar magAxis2(mag(axis2));

    if (magAxis1 < ROOTVSMALL)
    {
        FatalErrorInFunction
            << "Dominant coordinate axis cannot have zero length"
            << nl << endl
            << abort(FatalError);
    }

    const vector ax1(axis1 / magAxis1);  // normalise
    vector ax2(axis2);

    if (magAxis2 < ROOTVSMALL)
    {
        // axis2 == Zero : Use best-guess for the second axis.
        ax2 = findOrthogonal(axis1);
    }

    // Remove colinear component
    ax2 -= ((ax1 & ax2) * ax1);

    magAxis2 = mag(ax2);

    if (magAxis2 < SMALL)
    {
        WarningInFunction
            << "axis1, axis2 appear to be co-linear: "
            << axis1 << ", " << axis2 << "  Revert to guessing axis2"
            << nl << endl;

        ax2 = findOrthogonal(axis1);

        // Remove colinear component
        ax2 -= ((ax1 & ax2) * ax1);

        magAxis2 = mag(ax2);

        if (magAxis2 < SMALL)
        {
            FatalErrorInFunction
                << "Could not find an appropriate second axis"
                << nl << endl
                << abort(FatalError);
        }
    }

    ax2 /= magAxis2;  // normalise


    // The local axes are columns of the rotation matrix

    tensor rotTensor;

    switch (order)
    {
        case E1_E2:
        {
            rotTensor.col<0>(ax1);
            rotTensor.col<1>(ax2);
            rotTensor.col<2>(ax1^ax2);
            break;
        }
        case E2_E3:
        {
            rotTensor.col<0>(ax1^ax2);
            rotTensor.col<1>(ax1);
            rotTensor.col<2>(ax2);
            break;
        }
        case E3_E1:
        case E3_E1_COMPAT:
        {
            rotTensor.col<0>(ax2);
            rotTensor.col<1>(ax1^ax2);
            rotTensor.col<2>(ax1);
            break;
        }
    }

    return rotTensor;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::axesRotation::read(const dictionary& dict)
{
    vector axis1, axis2;
    axisOrder order = E3_E1;

    if
    (
        dict.readIfPresent("e1", axis1)
     && dict.readIfPresent("e2", axis2)
    )
    {
        order = E1_E2;
    }
    else if
    (
        dict.readIfPresent("e2", axis1)
     && dict.readIfPresent("e3", axis2)
    )
    {
        order = E2_E3;
    }
    else if
    (
        dict.readIfPresent("e3", axis1)
     && dict.readIfPresent("e1", axis2)
    )
    {
        order = E3_E1;
    }
    else if
    (
        dict.readIfPresent("axis", axis1)
     && dict.readIfPresent("direction", axis2)
    )
    {
        order = E3_E1_COMPAT;
    }
    else
    {
        FatalErrorInFunction
            << "No entries of the type (e1, e2) or (e2, e3) or (e3, e1) found"
            << exit(FatalError);
    }

    R_ = rotation(axis1, axis2, order);
    Rtr_ = R_.T();
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
    R_(rotation(axis, dir, order)),
    Rtr_(R_.T())
{}


Foam::axesRotation::axesRotation
(
    const vector& axis
)
:
    R_(rotation(axis, findOrthogonal(axis), E3_E1)),
    Rtr_(R_.T())
{}


Foam::axesRotation::axesRotation
(
    const dictionary& dict
)
:
    R_(sphericalTensor::I),
    Rtr_(sphericalTensor::I)
{
    read(dict);
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
    return nullptr;
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
    return nullptr;
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
    read(dict);
}


// ************************************************************************* //
