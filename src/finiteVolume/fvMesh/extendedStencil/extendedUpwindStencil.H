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
    Foam::extendedUpwindStencil

Description

SourceFiles
    extendedUpwindStencil.C
    extendedUpwindStencilTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef extendedUpwindStencil_H
#define extendedUpwindStencil_H

#include "extendedStencil.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class faceStencil;

/*---------------------------------------------------------------------------*\
                           Class extendedUpwindStencil Declaration
\*---------------------------------------------------------------------------*/

class extendedUpwindStencil
:
    public extendedStencil
{
    // Private data

        //- Swap map for getting neigbouring data
        autoPtr<mapDistribute> ownMapPtr_;
        autoPtr<mapDistribute> neiMapPtr_;

        //- Per face the stencil.
        labelListList ownStencil_;
        labelListList neiStencil_;



    // Private Member Functions

        //- Find most 'opposite' faces of cell
        void selectOppositeFaces
        (
            const boolList& nonEmptyFace,
            const scalar minOpposedness,
            const label faceI,
            const label cellI,
            DynamicList<label>& oppositeFaces
        ) const;

        //- Transport (centred) face stencil to 'opposite' face.
        void transportStencil
        (
            const boolList& nonEmptyFace,
            const labelListList& faceStencil,
            const scalar minOpposedness,
            const label faceI,
            const label cellI,
            const bool stencilHasNeighbour,

            DynamicList<label>& oppositeFaces,
            labelHashSet& faceStencilSet,
            labelList& transportedStencil
        ) const;

        //- Transport (centred) face stencil to 'opposite' faces.
        void transportStencils
        (
            const labelListList& faceStencil,
            const scalar minOpposedness,
            labelListList& ownStencil,
            labelListList& neiStencil
        );


        //- Disallow default bitwise copy construct
        extendedUpwindStencil(const extendedUpwindStencil&);

        //- Disallow default bitwise assignment
        void operator=(const extendedUpwindStencil&);


public:

    // Constructors

        //- Construct from mesh and uncompacted face stencil
        extendedUpwindStencil
        (
            const faceStencil&,
            const scalar minOpposedness
        );


    // Member Functions

        //- Return reference to the parallel distribution map
        const mapDistribute& ownMap() const
        {
            return ownMapPtr_();
        }

        //- Return reference to the parallel distribution map
        const mapDistribute& neiMap() const
        {
            return neiMapPtr_();
        }

        //- Return reference to the stencil
        const labelListList& ownStencil() const
        {
            return ownStencil_;
        }

        //- Return reference to the stencil
        const labelListList& neiStencil() const
        {
            return neiStencil_;
        }

        //- Sum vol field contributions to create face values
        template<class Type>
        tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > weightedSum
        (
            const surfaceScalarField& phi,
            const GeometricField<Type, fvPatchField, volMesh>& fld,
            const List<List<scalar> >& ownWeights,
            const List<List<scalar> >& neiWeights
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "extendedUpwindStencilTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
