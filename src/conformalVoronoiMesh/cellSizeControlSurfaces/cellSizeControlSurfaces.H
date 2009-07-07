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

Class
    Foam::cellSizeControlSurfaces

Description

SourceFiles
    cellSizeControlSurfacesI.H
    cellSizeControlSurfaces.C

\*---------------------------------------------------------------------------*/

#ifndef cellSizeControlSurfaces_H
#define cellSizeControlSurfaces_H

#include "searchableSurfaces.H"
#include "searchableSurfacesQueries.H"
#include "triSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class conformalVoronoiMesh;

class cellSizeFunction;


/*---------------------------------------------------------------------------*\
                   Class cellSizeControlSurfaces Declaration
\*---------------------------------------------------------------------------*/

class cellSizeControlSurfaces
{
    // Private data

        //- Reference to the conformalVoronoiMesh holding this object
        const conformalVoronoiMesh& cvMesh_;

        //- Reference to the searchableSurfaces object holding all geometry data
        const searchableSurfaces& allGeometry_;

        //- Indices of surfaces in allGeometry that are to be conformed to
        labelList surfaces_;

        //- A list of all of the cellSizeFunction objects
        PtrList<cellSizeFunction> cellSizeFunctions_;

        //- Regions where no cell size function is specified use defaultCellSize
        scalar defaultCellSize_;

        //- Assigning a priority to all requests for cell sizes, the highest
        //  overrules.  Regions requesting with the same priority take the
        //  average
        label defaultPriority_;


    // Private Member Functions

        //- Helper function to decide whether to terminate searching
        //  cellSizeFunctions by determining if the next cellSizeFunction has
        //  the same priority, includes looking off the end of the
        //  cellSizeFunctions_ list at the defaultPriority_ value
        bool samePriorityNext(label i) const;

        //- Disallow default bitwise copy construct
        cellSizeControlSurfaces(const cellSizeControlSurfaces&);

        //- Disallow default bitwise assignment
        void operator=(const cellSizeControlSurfaces&);


public:

    // Constructors

        //- Construct from dictionary and references to conformalVoronoiMesh and
        //  searchableSurfaces
        cellSizeControlSurfaces
        (
            const conformalVoronoiMesh& cvMesh,
            const searchableSurfaces& allGeometry,
            const dictionary& motionControlDict
        );

    //- Destructor
    ~cellSizeControlSurfaces();


    // Member Functions

        // Access

            //- Return reference to the searchableSurfaces object containing all
            //  of the geometry
            inline const searchableSurfaces& geometry() const;

            //- Return the surface indices
            inline const labelList& surfaces() const;

        // Query

            //- Return the cell size at the given location
            scalar cellSize(const point& pt) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cellSizeControlSurfacesI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
