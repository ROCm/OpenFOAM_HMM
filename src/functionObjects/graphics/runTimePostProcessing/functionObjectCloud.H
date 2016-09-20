/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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
    Foam::functionObjectCloud

Description

SourceFiles
    functionObjectCloud.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjectCloud_H
#define functionObjectCloud_H

#include "pointData.H"
#include "fieldVisualisationBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class functionObjectCloud Declaration
\*---------------------------------------------------------------------------*/

class functionObjectCloud
:
    public pointData,
    public fieldVisualisationBase
{
private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        functionObjectCloud(const functionObjectCloud&);

        //- Disallow default bitwise assignment
        void operator=(const functionObjectCloud&);


protected:

    // Protected data

        //- Name of cloud
        word cloudName_;

        //- Name of cloud function object result to render
        word functionObject_;

        //- Name of field to colour by
        word colourFieldName_;

        //- Actor
        vtkSmartPointer<vtkActor> actor_;


public:

    //- Run-time type information
    TypeName("functionObjectCloud");


    // Constructors

        //- Construct from dictionary
        functionObjectCloud
        (
            const runTimePostProcessing& parent,
            const dictionary& dict,
            const HashPtrTable<Function1<vector>, word>& colours
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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
