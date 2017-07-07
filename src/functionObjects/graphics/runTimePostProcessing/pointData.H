/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2016 OpenCFD Ltd.
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
    Foam::functionObjects::runTimePostPro::pointData

Description

SourceFiles
    pointData.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_runTimePostPro_pointData_H
#define functionObjects_runTimePostPro_pointData_H

#include "geometryBase.H"
#include "Enum.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class vtkActor;
class vtkPolyData;
class vtkPolyDataMapper;

namespace Foam
{
namespace functionObjects
{
namespace runTimePostPro
{

/*---------------------------------------------------------------------------*\
                          Class pointData Declaration
\*---------------------------------------------------------------------------*/

class pointData
:
    public geometryBase
{
public:

    // Public enumerations

        enum representationType
        {
            rtSphere,               //!< Sphere
            rtVector                //!< Vector
        };

        static const Enum<representationType> representationTypeNames;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        pointData(const pointData&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const pointData&) = delete;


protected:

    // Protected data

        //- Representation type
        representationType representation_;

        //- Max glyph length
        scalar maxGlyphLength_;

        //- Point colour
        autoPtr<Function1<vector>> pointColour_;


    // Protected Member Functions

        //- Add the point data to the renderer
        void addPoints
        (
            const label framei,
            vtkActor* actor,
            vtkPolyDataMapper* mapper,
            vtkPolyData* data
        ) const;


public:

    //- Run-time type information
    TypeName("pointData");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            pointData,
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
        pointData
        (
            const runTimePostProcessing& parent,
            const dictionary& dict,
            const HashPtrTable<Function1<vector>, word>& colours
        );


    // Selectors

        //- Return a reference to the selected RAS model
        static autoPtr<pointData> New
        (
            const runTimePostProcessing& parent,
            const dictionary& dict,
            const HashPtrTable<Function1<vector>, word>& colours,
            const word& pointDataName
        );


    //- Destructor
    virtual ~pointData();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace runTimePostPro
} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
