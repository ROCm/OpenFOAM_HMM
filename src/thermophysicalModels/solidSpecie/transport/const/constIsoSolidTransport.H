/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::constIsoSolidTransport

Description
    Constant properties Transport package.
    Templated into a given thermodynamics package (needed for thermal
    conductivity).

SourceFiles
    constIsoSolidTransportI.H
    constIsoSolidTransport.C

\*---------------------------------------------------------------------------*/

#ifndef constIsoSolidTransport_H
#define constIsoSolidTransport_H

#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class thermo> class constIsoSolidTransport;

template<class thermo>
inline constIsoSolidTransport<thermo> operator*
(
    const scalar,
    const constIsoSolidTransport<thermo>&
);

template<class thermo>
Ostream& operator<<
(
    Ostream&,
    const constIsoSolidTransport<thermo>&
);


/*---------------------------------------------------------------------------*\
                   Class constIsoSolidTransport Declaration
\*---------------------------------------------------------------------------*/

template<class thermo>
class constIsoSolidTransport
:
    public thermo
{
    // Private data

        //- Constant isotropic thermal conductivity
        scalar kappa_;


    // Private Member Functions

        //- Construct from components
        inline constIsoSolidTransport(const thermo& t, const scalar kappa);


public:

    // Constructors

        //- Construct as named copy
        inline constIsoSolidTransport
        (
            const word&,
            const constIsoSolidTransport&
        );

        //- Construct from Istream
        constIsoSolidTransport(const dictionary& dict);


    // Member functions

        //- Isotropic thermal conductivity [W/mK]
        inline scalar kappa(const scalar T) const;

        //- Un-isotropic thermal conductivity [W/mK]
        inline vector Kappa(const scalar T) const;

        //- Write to Ostream
        void write(Ostream& os) const;


    // Member operators

        inline constIsoSolidTransport& operator=
        (
            const constIsoSolidTransport&
        );
        inline void operator+=(const constIsoSolidTransport&);
        inline void operator-=(const constIsoSolidTransport&);


    // Friend operators

        friend constIsoSolidTransport operator* <thermo>
        (
            const scalar,
            const constIsoSolidTransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <thermo>
        (
            Ostream&,
            const constIsoSolidTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "constIsoSolidTransportI.H"

#ifdef NoRepository
#   include "constIsoSolidTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
