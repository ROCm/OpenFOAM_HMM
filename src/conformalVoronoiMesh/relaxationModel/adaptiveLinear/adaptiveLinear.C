/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "adaptiveLinear.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(adaptiveLinear, 0);
addToRunTimeSelectionTable(relaxationModel, adaptiveLinear, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

adaptiveLinear::adaptiveLinear
(
    const dictionary& relaxationDict,
    const conformalVoronoiMesh& cvMesh
)
:
    relaxationModel(typeName, relaxationDict, cvMesh),
    relaxationStart_(readScalar(coeffDict().lookup("relaxationStart"))),
    relaxationEnd_(readScalar(coeffDict().lookup("relaxationEnd"))),
    lastTimeValue_(cvMesh_.time().timeOutputValue()),
    relaxation_(relaxationStart_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar adaptiveLinear::relaxation()
{
    if (cvMesh_.time().timeOutputValue() > lastTimeValue_)
    {
        relaxation_ -=
        (relaxation_ - relaxationEnd_)
       /(
            (
                cvMesh_.time().endTime().value()
              - cvMesh_.time().timeOutputValue()
            )
           /(cvMesh_.time().timeOutputValue() - lastTimeValue_)
          + 1
        );

        lastTimeValue_ = cvMesh_.time().timeOutputValue();
    }

    return relaxation_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
