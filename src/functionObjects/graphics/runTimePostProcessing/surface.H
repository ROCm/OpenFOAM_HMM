/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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
    Foam::functionObjects::runTimePostPro::surface

Description

SourceFiles
    surface.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_runTimePostPro_surface_H
#define functionObjects_runTimePostPro_surface_H

#include "geometryBase.H"
#include "Enum.H"
#include "runTimeSelectionTables.H"

#include "vtkSmartPointer.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class vtkActor;
class vtkRenderer;
class vtkPolyData;

namespace Foam
{
namespace functionObjects
{
namespace runTimePostPro
{

/*---------------------------------------------------------------------------*\
                          Class surface Declaration
\*---------------------------------------------------------------------------*/

class surface
:
    public geometryBase
{
public:

    // Public enumerations

        enum representationType
        {
            rtNone,
            rtWireframe,
            rtSurface,
            rtSurfaceWithEdges,
            rtGlyph
        };

        static const Enum<representationType> representationTypeNames;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        surface(const surface&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const surface&) = delete;


protected:

    // Protected data

        //- Representation type
        representationType representation_;

        //- Activate feature edges
        bool featureEdges_;

        //- Surface colour
        autoPtr<Function1<vector>> surfaceColour_;

        //- Edge colour
        autoPtr<Function1<vector>> edgeColour_;

        //- Surface actor
        vtkSmartPointer<vtkActor> surfaceActor_;

        //- Edge actor
        vtkSmartPointer<vtkActor> edgeActor_;

        //- Max glyph length for representation type rtGlyph
        scalar maxGlyphLength_;


    // Protected Member Functions

        //- Set the representation
        void setRepresentation(vtkActor* actor) const;

        //- Add feature edges to scene
        void addFeatureEdges
        (
            vtkRenderer* renderer,
            vtkPolyData* data
        ) const;


public:

    //- Run-time type information
    TypeName("surface");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            surface,
            dictionary,
            (
                const runTimePostProcessing& parent,
                const dictionary& dict,
                const HashPtrTable<Function1<vector>, word>& colours
            ),
            (parent, dict, colours)
        );


    // Constructors

        //- Construct from dictionary
        surface
        (
            const runTimePostProcessing& parent,
            const dictionary& dict,
            const HashPtrTable<Function1<vector>, word>& colours
        );


    // Selectors

        //- Return a reference to the selected RAS model
        static autoPtr<surface> New
        (
            const runTimePostProcessing& parent,
            const dictionary& dict,
            const HashPtrTable<Function1<vector>, word>& colours,
            const word& surfaceName
        );


    //- Destructor
    virtual ~surface();


    // Member Functions

        //- Update the actors
        virtual void updateActors(const scalar position);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace runTimePostPro
} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
