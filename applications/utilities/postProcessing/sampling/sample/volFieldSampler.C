/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

Description

\*---------------------------------------------------------------------------*/

#include "volFieldSampler.H"
#include "volPointInterpolation.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
template <class Type>
Foam::volFieldSampler<Type>::volFieldSampler
(
    const volPointInterpolation& pInterp,
    const dictionary& interpolationSchemes,
    const GeometricField<Type, fvPatchField, volMesh>& field,
    const PtrList<sampleSet>& samplers
)
:
    List<Field<Type> >(samplers.size()),
    name_(field.name())
{
    autoPtr<interpolation<Type> > interpolator
    (
        interpolation<Type>::New(interpolationSchemes, pInterp, field)
    );

    forAll(samplers, setI)
    {
        Field<Type>& values = this->operator[](setI);

        const sampleSet& samples = samplers[setI];

        values.setSize(samples.size());
        forAll(samples, sampleI)
        {
            const point& samplePt = samples[sampleI];
            label cellI = samples.cells()[sampleI];
            label faceI = samples.faces()[sampleI];

            values[sampleI] =
                interpolator().interpolate
                (
                    samplePt,
                    cellI,
                    faceI
                );
        }
    }
}


// Construct from components
template <class Type>
Foam::volFieldSampler<Type>::volFieldSampler
(
    const List<Field<Type> >& values,
    const word& name
)
:
    List<Field<Type> >(values),
    name_(name)
{}


// ************************************************************************* //
