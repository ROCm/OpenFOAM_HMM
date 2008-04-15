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
    Foam::fieldsCache

Description
        
SourceFiles
    fieldsCache.C

\*---------------------------------------------------------------------------*/

#ifndef fieldsCache_H
#define fieldsCache_H

#include "volFields.H"
#include "pointFields.H"
#include "HashPtrTable.H"
#include "interpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class fieldsCache Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class fieldsCache
:
    public HashPtrTable<GeometricField<Type, fvPatchField, volMesh> >
{

private:

    //- Cached pointFields
    mutable HashPtrTable<GeometricField<Type, pointPatchField, pointMesh> >
        pointFields_;

    //- Cached interpolators
    mutable HashPtrTable<interpolation<Type> > interpolators_;


public:

    // Constructors

        //- Construct null
        fieldsCache();


    // Member Functions

        //- Get pointField of named field
        const GeometricField<Type, pointPatchField, pointMesh>& pointField
        (
            const word&,
            const volPointInterpolation& pInterp
        ) const;

        //- Get interpolator for named field
        const interpolation<Type>& interpolator
        (
            const word&,
            const volPointInterpolation& pInterp,
            const dictionary& interpolationSchemes
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "fieldsCache.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
