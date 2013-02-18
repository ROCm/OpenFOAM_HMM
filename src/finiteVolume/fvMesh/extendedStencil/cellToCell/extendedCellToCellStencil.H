/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
     \\/     M anipulation  |
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

Class
    Foam::extendedCellToCellStencil

Description
    Baseclass for cell-to-cell stencils

SourceFiles
    extendedCellToCellStencil.C
    extendedCellToCellStencilTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef extendedCellToCellStencil_H
#define extendedCellToCellStencil_H

#include "mapDistribute.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class globalIndex;

/*---------------------------------------------------------------------------*\
                  Class extendedCellToCellStencil Declaration
\*---------------------------------------------------------------------------*/

class extendedCellToCellStencil
{
protected:

    // Protected data

        const polyMesh& mesh_;

public:

    // Constructors

        //- Construct from mesh
        explicit extendedCellToCellStencil(const polyMesh&);


    // Member Functions

//        //- Use map to get the data into stencil order
//        template<class T>
//        static void collectData
//        (
//            const mapDistribute& map,
//            const labelListList& stencil,
//            const GeometricField<T, fvsPatchField, surfaceMesh>& fld,
//            List<List<T> >& stencilFld
//        );
//
        //- Sum surface field contributions to create cell values
        template<class Type>
        static tmp<GeometricField<Type, fvPatchField, volMesh> > weightedSum
        (
            const mapDistribute& map,
            const labelListList& stencil,
            const GeometricField<Type, fvPatchField, volMesh>& fld,
            const List<List<scalar> >& stencilWeights
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "extendedCellToCellStencilTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
