/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2018 OpenCFD Ltd.
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
    Foam::cylindricalCS

Description
    Cylindrical coordinate system, always in radians.

SourceFiles
    cylindricalCS.C

\*---------------------------------------------------------------------------*/

#ifndef cylindricalCS_H
#define cylindricalCS_H

#include "coordinateSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class cylindricalCS Declaration
\*---------------------------------------------------------------------------*/

class cylindricalCS
:
    public coordinateSystem
{
protected:

    // Protected Member Functions


        //- Convert from local coordinate system to the global Cartesian system
        //  with optional translation for the origin
        virtual vector localToGlobal(const vector&, bool translate) const;

        //- Convert from local coordinate system to the global Cartesian system
        //  with optional translation for the origin
        virtual tmp<vectorField> localToGlobal
        (
            const vectorField&,
            bool translate
        ) const;

        //- Convert from global Cartesian system to the local coordinate system
        //  with optional translation for the origin
        virtual vector globalToLocal(const vector&, bool translate) const;

        //- Convert from global Cartesian system to the local coordinate system
        //  with optional translation for the origin
        virtual tmp<vectorField> globalToLocal
        (
            const vectorField&,
            bool translate
        ) const;


public:


    // Constructors

        //- Construct null
        cylindricalCS();

        //- Construct copy
        cylindricalCS
        (
            const coordinateSystem&
        );

        //- Construct copy with a different name
        cylindricalCS
        (
            const word& name,
            const coordinateSystem&
        );

        //- Construct from origin and rotation
        cylindricalCS
        (
            const word& name,
            const point& origin,
            const coordinateRotation&
        );

        //- Construct from origin and 2 axes
        cylindricalCS
        (
            const word& name,
            const point& origin,
            const vector& axis,
            const vector& dirn
        );

        //- Construct from dictionary and name
        cylindricalCS(const word&, const dictionary&);

        //- Construct from dictionary and objectRegistry
        cylindricalCS(const objectRegistry&, const dictionary&);


    //- Destructor
    virtual ~cylindricalCS() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
