/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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
    Foam::cellShapeControl

Description

SourceFiles
    cellShapeControlI.H
    cellShapeControl.C

\*---------------------------------------------------------------------------*/

#ifndef cellShapeControl_H
#define cellShapeControl_H

#include "dictionary.H"
#include "autoPtr.H"
#include "tensor.H"
#include "point.H"
#include "primitiveFieldsFwd.H"
#include "pointFieldFwd.H"
#include "triadField.H"
#include "Time.H"
#include "searchableSurfaces.H"
#include "conformationSurfaces.H"
#include "cellAspectRatioControl.H"
#include "cellSizeAndAlignmentControls.H"
#include "cellShapeControlMesh.H"
#include "backgroundMeshDecomposition.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                     Class cellShapeControl Declaration
\*---------------------------------------------------------------------------*/

class cellShapeControl
:
    public dictionary
{
    // Private data

        const Time& runTime_;

        const searchableSurfaces& allGeometry_;

        const conformationSurfaces& geometryToConformTo_;

        const scalar defaultCellSize_;

        cellShapeControlMesh shapeControlMesh_;

        cellAspectRatioControl aspectRatio_;

        cellSizeAndAlignmentControls sizeAndAlignment_;


    // Private Member Functions

        template <class Triangulation, class Type>
        tmp<Field<Type> > filterFarPoints
        (
            const Triangulation& mesh,
            const Field<Type>& field
        );

        template <class Triangulation>
        autoPtr<mapDistribute> buildMap
        (
            const Triangulation& mesh,
            labelListList& pointPoints
        );

        template <class Triangulation>
        autoPtr<mapDistribute> buildReferredMap
        (
            const Triangulation& mesh,
            labelList& indices
        );

        template <class Triangulation>
        tmp<triadField> buildAlignmentField(const Triangulation& mesh);

        template <class Triangulation>
        tmp<pointField> buildPointField(const Triangulation& mesh);

        //- Disallow default bitwise copy construct
        cellShapeControl(const cellShapeControl&);

        //- Disallow default bitwise assignment
        void operator=(const cellShapeControl&);


public:

    //- Runtime type information
    ClassName("cellShapeControl");


    // Constructors

        //- Construct from dictionary and references to conformalVoronoiMesh and
        //  searchableSurfaces
        cellShapeControl
        (
            const Time& runTime,
            const dictionary& cellShapeControlDict,
            const searchableSurfaces& allGeometry,
            const conformationSurfaces& geometryToConformTo
        );


    //- Destructor
    ~cellShapeControl();


    // Member Functions

        // Access

            inline const scalar& defaultCellSize() const;

            inline cellShapeControlMesh& shapeControlMesh();

            inline const cellShapeControlMesh& shapeControlMesh() const;

            inline const cellAspectRatioControl& aspectRatio() const;

            inline const cellSizeAndAlignmentControls& sizeAndAlignment() const;


        // Query

            //- Return the cell size at the given location
            scalar cellSize(const point& pt) const;

            scalarField cellSize(const pointField& pts) const;

            //- Return the cell alignment at the given location
            tensor cellAlignment(const point& pt) const;

            void cellSizeAndAlignment
            (
                const point& pt,
                scalar& size,
                tensor& alignment
            ) const;


        // Edit

            void initialMeshPopulation
            (
                const autoPtr<backgroundMeshDecomposition>& decomposition
            );

            label refineMesh
            (
                const autoPtr<backgroundMeshDecomposition>& decomposition
            );

            void smoothMesh();

            //- Add a control point with a specified size and alignment
//            virtual void addControlPoint
//            (
//                const point& pt,
//                const scalar& size,
//                const tensor& alignment
//            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cellShapeControlI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
