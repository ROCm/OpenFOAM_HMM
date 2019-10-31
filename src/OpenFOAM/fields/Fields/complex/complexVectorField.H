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

Typedef
    Foam::complexVectorField

Description
    Specialisation of Field\<T\> for complexVector.

SourceFiles
    complexVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef complexVectorField_H
#define complexVectorField_H

#include "complex.H"
#include "complexVector.H"
#include "vectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef Field<complexVector> complexVectorField;

//- Zip up two lists of values into a list of complex
complexVectorField ComplexField
(
    const UList<vector>& re,
    const UList<vector>& im
);


//- Create complex field from a list of real (using imag == 0)
complexVectorField ReComplexField(const UList<vector>& re);

//- Create complex field from a list of imag (using real == 0)
complexVectorField ImComplexField(const UList<vector>& im);

//- Extract real component
vectorField Re(const UList<complexVector>& cvf);

//- Extract imag component
vectorField Im(const UList<complexVector>& cvf);

//- Sum real and imag components
vectorField ReImSum(const UList<complexVector>& cvf);

complexVectorField operator^
(
    const UList<vector>& vf,
    const UList<complexVector>& cvf
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
