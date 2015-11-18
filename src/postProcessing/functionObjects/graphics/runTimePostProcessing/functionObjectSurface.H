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
    Foam::functionObjectSurface

Description

SourceFiles
    functionObjectSurface.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjectSurface_H
#define functionObjectSurface_H

#include "geometrySurface.H"
#include "fieldVisualisationBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class functionObjectSurface Declaration
\*---------------------------------------------------------------------------*/

class functionObjectSurface
:
    public geometrySurface,
    public fieldVisualisationBase
{


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        functionObjectSurface(const functionObjectSurface&);

        //- Disallow default bitwise assignment
        void operator=(const functionObjectSurface&);


protected:

    // Protected data

        //- Name of function object result to render
        word functionObject_;


public:

    //- Run-time type information
    TypeName("functionObject");


    // Constructors

        //- Construct from dictionary
        functionObjectSurface
        (
            const runTimePostProcessing& parent,
            const dictionary& dict,
            const HashPtrTable<DataEntry<vector>, word>& colours
        );


    //- Destructor
    virtual ~functionObjectSurface();


    // Member Functions

        //- Add surface(s) to scene
        virtual void addGeometryToScene
        (
            const scalar position,
            vtkRenderer* renderer
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
