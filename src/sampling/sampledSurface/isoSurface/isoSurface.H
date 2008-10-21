/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
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
    Foam::isoSurface

Description
    A surface formed by the iso value.
    After "Polygonising A Scalar Field Using Tetrahedrons", Paul Bourke

SourceFiles
    isoSurface.C

\*---------------------------------------------------------------------------*/

#ifndef isoSurface_H
#define isoSurface_H

#include "triSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;

/*---------------------------------------------------------------------------*\
                       Class isoSurface Declaration
\*---------------------------------------------------------------------------*/

class isoSurface
:
    public triSurface
{
    // Private data

        //- Reference to mesh
        const polyMesh& mesh_;

        //- Reference to cell values
        const scalarField& cellValues_;

        //- Reference to point values
        const scalarField& pointValues_;

        //- Isosurface value
        const scalar iso_;


        //- For every triangle the original cell in mesh
        labelList meshCells_;

        //- For every unmerged triangle point the point in the triSurface
        labelList triPointMergeMap_;


    // Private Member Functions

        template<class T>
        static T vertexInterp
        (
            const scalar iso,
            const T& p0,
            const T& p1,
            const scalar s0,
            const scalar s1
        );

        template<class T>
        static void vertexInterp
        (
            const scalar iso,
            const scalar s0,
            const scalar s1,
            const scalar s2,
            const scalar s3,

            const T& p0,
            const T& p1,
            const T& p2,
            const T& p3,

            DynamicList<T>& points
        );

public:

    //- Runtime type information
    TypeName("isoSurface");


    // Constructors

        //- Construct from dictionary
        isoSurface
        (
            const polyMesh& mesh,
            const scalarField& cellValues,
            const scalarField& pointValues,
            const scalar iso
        );


    // Member Functions

        //- For every face original cell in mesh
        const labelList& meshCells() const
        {
            return meshCells_;
        }

        //- For every unmerged triangle point the point in the triSurface
        const labelList triPointMergeMap() const
        {
            return triPointMergeMap_;
        }


        //- sample field on faces
        template <class Type>
        tmp<Field<Type> > sample
        (
            const Field<Type>& sampleCellValues
        ) const;

        //- interpolate field to points
        template <class Type>
        tmp<Field<Type> >
        interpolate
        (
            const Field<Type>& sampleCellValues,
            const Field<Type>& samplePointValues
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "isoSurfaceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
