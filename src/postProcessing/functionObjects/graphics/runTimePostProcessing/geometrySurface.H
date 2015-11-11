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
    Foam::geometrySurface

Description

SourceFiles
    geometrySurface.C

\*---------------------------------------------------------------------------*/

#ifndef geometrySurface_H
#define geometrySurface_H

#include "surface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class vtkPolyData;

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class geometrySurface Declaration
\*---------------------------------------------------------------------------*/

class geometrySurface
:
    public surface
{
private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        geometrySurface(const geometrySurface&);

        //- Disallow default bitwise assignment
        void operator=(const geometrySurface&);


protected:

    // Protected data

        //- File names
        List<fileName> fileNames_;


    // Protected Member Functions

        //- Add surface (file) to scene
        void addGeometryToScene
        (
            const scalar position,
            vtkRenderer* renderer,
            const fileName& fName
        ) const;


public:

    //- Run-time type information
    TypeName("geometry");


    // Constructors

        //- Construct from dictionary
        geometrySurface
        (
            const runTimePostProcessing& parent,
            const dictionary& dict,
            const HashPtrTable<DataEntry<vector>, word>& colours
        );

        //- Construct from components
        geometrySurface
        (
            const runTimePostProcessing& parent,
            const dictionary& dict,
            const HashPtrTable<DataEntry<vector>, word>& colours,
            const List<fileName>& fileNames
        );


    //- Destructor
    virtual ~geometrySurface();


    // Member Functions

        //- Add surface(s) to scene
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
