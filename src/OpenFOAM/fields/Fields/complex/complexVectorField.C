/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "complexVectorField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineCompoundTypeName(List<complexVector>, complexVectorList);
    addCompoundToRunTimeSelectionTable(List<complexVector>, complexVectorList);
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

Foam::complexVectorField Foam::ComplexField
(
    const UList<vector>& re,
    const UList<vector>& im
)
{
    complexVectorField cvf(re.size());

    for (direction cmpt=0; cmpt<vector::nComponents; ++cmpt)
    {
        forAll(cvf, i)
        {
            cvf[i].component(cmpt).Re() = re[i].component(cmpt);
            cvf[i].component(cmpt).Im() = im[i].component(cmpt);
        }
    }

    return cvf;
}


Foam::complexVectorField Foam::ReComplexField(const UList<vector>& re)
{
    complexVectorField cvf(re.size());

    for (direction cmpt=0; cmpt<vector::nComponents; ++cmpt)
    {
        forAll(cvf, i)
        {
            cvf[i].component(cmpt).Re() = re[i].component(cmpt);
            cvf[i].component(cmpt).Im() = 0.0;
        }
    }

    return cvf;
}


Foam::complexVectorField Foam::ImComplexField(const UList<vector>& im)
{
    complexVectorField cvf(im.size());

    for (direction cmpt=0; cmpt<vector::nComponents; ++cmpt)
    {
        forAll(cvf, i)
        {
            cvf[i].component(cmpt).Re() = 0.0;
            cvf[i].component(cmpt).Im() = im[i].component(cmpt);
        }
    }

    return cvf;
}


Foam::vectorField Foam::ReImSum(const UList<complexVector>& cvf)
{
    vectorField vf(cvf.size());

    for (direction cmpt=0; cmpt<vector::nComponents; ++cmpt)
    {
        forAll(cvf, i)
        {
            vf[i].component(cmpt) =
                cvf[i].component(cmpt).Re() + cvf[i].component(cmpt).Im();
        }
    }

    return vf;
}


Foam::vectorField Foam::Re(const UList<complexVector>& cvf)
{
    vectorField vf(cvf.size());

    for (direction cmpt=0; cmpt<vector::nComponents; ++cmpt)
    {
        forAll(cvf, i)
        {
            vf[i].component(cmpt) = cvf[i].component(cmpt).Re();
        }
    }

    return vf;
}


Foam::vectorField Foam::Im(const UList<complexVector>& cvf)
{
    vectorField vf(cvf.size());

    for (direction cmpt=0; cmpt<vector::nComponents; ++cmpt)
    {
        forAll(cvf, i)
        {
            vf[i].component(cmpt) = cvf[i].component(cmpt).Im();
        }
    }

    return vf;
}


Foam::complexVectorField Foam::operator^
(
    const UList<vector>& vf,
    const UList<complexVector>& cvf
)
{
    return ComplexField(vf^Re(cvf), vf^Im(cvf));
}


// ************************************************************************* //
