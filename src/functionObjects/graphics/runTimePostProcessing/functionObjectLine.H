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
    Foam::functionObjects::runTimePostPro::functionObjectLine

Description
    Visualisation of line data from function object output

SourceFiles
    functionObjectLine.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_runTimePostPro_functionObjectLine_H
#define functionObjects_runTimePostPro_functionObjectLine_H

#include "pathline.H"
#include "functionObjectBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace runTimePostPro
{

/*---------------------------------------------------------------------------*\
                     Class functionObjectLine Declaration
\*---------------------------------------------------------------------------*/

class functionObjectLine
:
    public pathline,
    public functionObjectBase
{
private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        functionObjectLine(const functionObjectLine&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const functionObjectLine&) = delete;


protected:

    // Protected data

        //- Actor
        vtkSmartPointer<vtkActor> actor_;


public:

    //- Run-time type information
    TypeName("functionObjectLine");


    // Constructors

        //- Construct from dictionary
        functionObjectLine
        (
            const runTimePostProcessing& parent,
            const dictionary& dict,
            const HashPtrTable<Function1<vector>, word>& colours
        );


    //- Destructor
    virtual ~functionObjectLine();


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
