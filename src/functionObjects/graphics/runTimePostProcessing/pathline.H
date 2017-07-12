/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd
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
    Foam::pathline

Description

SourceFiles
    pathline.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_runTimePostPro_pathline_H
#define functionObjects_runTimePostPro_pathline_H

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
                          Class pathline Declaration
\*---------------------------------------------------------------------------*/

class pathline
:
    public geometryBase
{
public:

    // Public enumerations

        enum representationType
        {
            rtNone,
            rtLine,
            rtTube,
            rtVector
        };

        static const Enum<representationType> representationTypeNames;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        pathline(const pathline&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const pathline&) = delete;


protected:

    // Protected data

        //- Representation type
        representationType representation_;

        //- Radius for the tube filter
        scalar tubeRadius_;

        //- Line colour
        autoPtr<Function1<vector>> lineColour_;


    // Protected Member Functions

        //- Add the pathlines to the renderer
        void addLines
        (
            const label framei,
            vtkActor* actor,
            vtkPolyData* data
        ) const;


public:

    //- Run-time type information
    TypeName("pathline");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            pathline,
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
        pathline
        (
            const runTimePostProcessing& parent,
            const dictionary& dict,
            const HashPtrTable<Function1<vector>, word>& colours
        );


    // Selectors

        //- Return a reference to the selected RAS model
        static autoPtr<pathline> New
        (
            const runTimePostProcessing& parent,
            const dictionary& dict,
            const HashPtrTable<Function1<vector>, word>& colours,
            const word& pathlineName
        );


    //- Destructor
    virtual ~pathline();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace runTimePostPro
} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
