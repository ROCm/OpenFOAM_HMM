/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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
    Foam::interpolateSolidThermo

Description
    The thermophysical properties of a interpolateSolidThermo

SourceFiles
    interpolateSolidThermo.C

\*---------------------------------------------------------------------------*/

#ifndef interpolateSolidThermo_H
#define interpolateSolidThermo_H

#include "basicSolidThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class basicSolidThermo;

Ostream& operator<<
(
     Ostream&,
     const basicSolidThermo&
);


/*---------------------------------------------------------------------------*\
                            Class interpolateSolidThermo Declaration
\*---------------------------------------------------------------------------*/

class interpolateSolidThermo
:
    public basicSolidThermo
{
    // Private data

        //- Temperature points for which there are values
        const Field<scalar> TValues_;

        //- Density at given temperatures
        const Field<scalar> rhoValues_;

        const Field<scalar> cpValues_;

        const Field<symmTensor> KValues_;

        const Field<scalar> HfValues_;

        const Field<scalar> emissivityValues_;

public:

    //- Runtime type information
    TypeName("interpolateSolidThermo");


    // Constructors

        //- Construct from mesh
        interpolateSolidThermo(const fvMesh& mesh, const dictionary& dict);

    // Destructor

        virtual ~interpolateSolidThermo();


    // Member Functions

        //- Update properties
        virtual void correct();


        // I-O

            //- Write the interpolateSolidThermo properties
            virtual void write(Ostream& os) const;


        // Ostream Operator

        friend Ostream& operator<<
        (
            Ostream& os,
            const interpolateSolidThermo& s
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
