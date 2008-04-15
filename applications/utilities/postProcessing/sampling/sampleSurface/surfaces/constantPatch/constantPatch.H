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
    Foam::constantPatch

Description

SourceFiles
    constantPatch.C

\*---------------------------------------------------------------------------*/

#ifndef constantPatch_H
#define constantPatch_H

#include "surface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//class meshSearch;
class dictionary;
class plane;
class volPointInterpolation;
template<class Type> class fieldsCache;

/*---------------------------------------------------------------------------*\
                           Class constantPatch Declaration
\*---------------------------------------------------------------------------*/

class constantPatch
:
    public surface
{
    // Private data

        //- Name of patch
        const word patchName_;

        //- Index of patch in boundaryMesh
        const label patchIndex_;

        //- Make triangles or keep faces
        const bool triangulate_;

        //- Zero size or copy of patch.localPoints()
        pointField points_;

        //- Generated triangles
        faceList faces_;

        //- Original patch face labels
        labelList patchFaceLabels_;

    // Private Member Functions

        //- Triangulate patch
        void makeTriangles();

        void copyFaces();

        //- Do all to construct geometry.
        void createGeometry();

        //- Interpolate field onto faces
        template <class Type>
        tmp<Field<Type> > interpolate
        (
            const word&,
            const fieldsCache<Type>&
        ) const;


public:

    //- Runtime type information
    TypeName("constantPatch");


    // Constructors

        //- Construct from components
        constantPatch
        (
            const polyMesh& mesh,
            meshSearch& searchEngine,
            const word& name,
            const word& patchName,
            const bool triangulate = true
        );

        //- Construct from dictionary
        constantPatch
        (
            const polyMesh& mesh,
            meshSearch& searchEngine,
            const dictionary& dict          
        );


    // Destructor

        virtual ~constantPatch();


    // Member Functions

        const word patchName() const
        {
            return patchName_;
        }

        label patchIndex() const
        {
            return patchIndex_;
        }

        const labelList& patchFaceLabels() const
        {
            return patchFaceLabels_;
        }

        //- Points of surface
        virtual const pointField& points() const
        {
            return points_;
        }

        //- Faces of surface
        virtual const faceList& faces() const
        {
            return faces_;
        }

        //- Correct for mesh movement and/or field changes
        virtual void correct
        (
            const bool meshChanged,
            const volPointInterpolation&,
            const dictionary& interpolationSchemes,
            const fieldsCache<scalar>&
        );

        //- interpolate field to surface
        virtual tmp<scalarField> interpolate
        (
            const word&,
            const fieldsCache<scalar>&,
            const volPointInterpolation&,
            const dictionary& interpolationSchemes
        ) const;

        //- interpolate field to surface
        virtual tmp<vectorField> interpolate
        (
            const word&,
            const fieldsCache<vector>&,
            const volPointInterpolation&,
            const dictionary& interpolationSchemes
        ) const;

        //- interpolate field to surface
        virtual tmp<sphericalTensorField> interpolate
        (
            const word&,
            const fieldsCache<sphericalTensor>&,
            const volPointInterpolation&,
            const dictionary& interpolationSchemes
        ) const;

        //- interpolate field to surface
        virtual tmp<symmTensorField> interpolate
        (
            const word&,
            const fieldsCache<symmTensor>&,
            const volPointInterpolation&,
            const dictionary& interpolationSchemes
        ) const;

        //- interpolate field to surface
        virtual tmp<tensorField> interpolate
        (
            const word&,
            const fieldsCache<tensor>&,
            const volPointInterpolation&,
            const dictionary& interpolationSchemes
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
