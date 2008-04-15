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

Class
    Foam::volFieldSampler

Description
    Samples a volField given a list of sampleSets (i.e. points).

SourceFiles
    volFieldSampler.C

\*---------------------------------------------------------------------------*/

#ifndef volFieldSampler_H
#define volFieldSampler_H

#include "sampleSet.H"
#include "List.H"
#include "PtrList.H"
#include "volFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class volPointInterpolation;
class dictionary;

/*---------------------------------------------------------------------------*\
                           Class volFieldSampler Declaration
\*---------------------------------------------------------------------------*/

template <class Type>
class volFieldSampler
:
    public List<Field<Type> >
{
    //- Name of this collection of values
    const word name_;

public:

    // Constructors

        //- Construct from sampleSet, interpolator, Field
        volFieldSampler
        (
            const volPointInterpolation&,
            const dictionary& interpolationSchemes,
            const GeometricField<Type, fvPatchField, volMesh>& field,
            const PtrList<sampleSet>&
        );

        //- Construct from components
        volFieldSampler
        (
            const List<Field<Type> >& values,
            const word& name
        );

    // Member functions

        const word& name() const
        {
            return name_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "volFieldSampler.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
