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
    Foam::text

Description

SourceFiles
    text.C

\*---------------------------------------------------------------------------*/

#ifndef text_H
#define text_H

#include "geometryBase.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class vtkRenderer;

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class text Declaration
\*---------------------------------------------------------------------------*/

class text
:
    public geometryBase
{
private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        text(const text&);

        //- Disallow default bitwise assignment
        void operator=(const text&);


protected:

    // Protected data

        //- Text
        string string_;

        //- Position
        Tuple2<scalar, scalar> position_;

        //- Size
        scalar size_;

        //- Colour
        autoPtr<DataEntry<vector> > colour_;

        //- Bold flag
        bool bold_;


public:

    // Constructors

        //- Construct from dictionary
        text
        (
            const runTimePostProcessing& parent,
            const dictionary& dict,
            const HashPtrTable<DataEntry<vector>, word>& colours
        );


    //- Destructor
    virtual ~text();


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
