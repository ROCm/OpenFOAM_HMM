/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017-2018 OpenCFD Ltd.
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
    Foam::coordinateRotations::euler

Description
    A coordinateRotation defined in the z-x-z (intrinsic) Euler convention.

    The 3 rotations are defined in the Euler intrinsic convention
    (around Z, around X' and around Z'').
    The order of the parameter arguments matches this rotation order.

    For reference and illustration, see
    https://en.wikipedia.org/wiki/Euler_angles

    \verbatim
    coordinateRotation
    {
        type    euler;
        angles  (0 0 180);
    }
    \endverbatim

    \heading Dictionary entries
    \table
        Property    | Description                           | Required | Default
        type        | Type name: euler (or EulerRotation)   | yes   |
        angles      | The z-x-z rotation angles             | yes   |
        degrees     | Angles are in degrees                 | no    | true
    \endtable

SourceFiles
    EulerCoordinateRotation.C

\*---------------------------------------------------------------------------*/

#ifndef coordinateRotations_euler_H
#define coordinateRotations_euler_H

#include "coordinateRotation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace coordinateRotations
{

/*---------------------------------------------------------------------------*\
                 Class coordinateRotations::euler Declaration
\*---------------------------------------------------------------------------*/

class euler
:
    public coordinateRotation
{
    // Private Data

        //- The rotation angles
        vector angles_;

        //- Angles measured in degrees
        bool degrees_;


public:

    //- Runtime type information
    TypeNameNoDebug("euler");


    // Constructors

        //- Construct null - an identity transform
        euler();

        //- Copy construct
        euler(const euler& crot);

        //- Construct from Euler rotation angles (z-x-z)
        euler(const vector& phiThetaPsi, bool degrees);

        //- Construct from Euler rotation angles (z-x-z)
        euler(scalar phi, scalar theta, scalar psi, bool degrees);

        //- Construct from dictionary
        explicit euler(const dictionary& dict);

        //- Return clone
        autoPtr<coordinateRotation> clone() const
        {
            return
                autoPtr<coordinateRotation>::NewFrom
                <coordinateRotations::euler>(*this);
        }


    //- Destructor
    virtual ~euler() = default;


    // Static Member Functions

        //- The rotation tensor calculated for the specified Euler angles
        //- interpreted as phi/theta/psi (z-x-z order)
        static tensor rotation(const vector& angles, bool degrees);


    // Member Functions

        //- Reset specification
        virtual void clear();

        //- The rotation tensor calculated for the specified Euler angles.
        virtual tensor R() const;

        //- Write information
        virtual void write(Ostream& os) const;

        //- Write dictionary entry
        virtual void writeEntry(const word& keyword, Ostream& os) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace coordinateRotations


//- Compatibility typedef 1806
typedef coordinateRotations::euler EulerCoordinateRotation;

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
