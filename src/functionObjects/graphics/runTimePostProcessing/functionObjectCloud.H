/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2018 OpenCFD Ltd.
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
    Foam::functionObjects::runTimePostPro::functionObjectCloud

Description
    Visualisation of cloud data from function object output

SourceFiles
    functionObjectCloud.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_runTimePostPro_functionObjectCloud_H
#define functionObjects_runTimePostPro_functionObjectCloud_H

#include "pointData.H"
#include "functionObjectBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace runTimePostPro
{

/*---------------------------------------------------------------------------*\
                     Class functionObjectCloud Declaration
\*---------------------------------------------------------------------------*/

class functionObjectCloud
:
    public pointData,
    public functionObjectBase
{
protected:

    // Protected Data

        //- Name of functionObjectCloud
        word cloudName_;

        //- The input filename used
        fileName inputFileName_;

        //- Name of field to colour by
        word colourFieldName_;

        //- Actor
        vtkSmartPointer<vtkActor> actor_;


    // Protected Member Functions

        //- No copy construct
        functionObjectCloud(const functionObjectCloud&) = delete;

        //- No copy assignment
        void operator=(const functionObjectCloud&) = delete;


public:

    //- Run-time type information
    TypeName("functionObjectCloud");


    // Constructors

        //- Construct from dictionary
        functionObjectCloud
        (
            const runTimePostProcessing& parent,
            const dictionary& dict,
            const HashPtrTable<Function1<vector>>& colours
        );


    //- Destructor
    virtual ~functionObjectCloud();


    // Member Functions

        //- Add tube(s) to scene
        virtual void addGeometryToScene
        (
            const scalar position,
            vtkRenderer* renderer
        );

        //- Update actors
        virtual void updateActors(const scalar position);

        //- Clear files used to create the object(s)
        virtual bool clear();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace runTimePostPro
} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
