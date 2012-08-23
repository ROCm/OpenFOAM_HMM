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
    Foam::constAnIsoSolidTransport

Description
    Constant properties Transport package.
    Templated into a given thermodynamics package (needed for thermal
    conductivity).

SourceFiles
    constAnIsoSolidTransportI.H
    constAnIsoSolidTransport.C

\*---------------------------------------------------------------------------*/

#ifndef constAnIsoSolidTransport_H
#define constAnIsoSolidTransport_H
#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class thermo> class constAnIsoSolidTransport;

template<class thermo>
inline constAnIsoSolidTransport<thermo> operator*
(
    const scalar,
    const constAnIsoSolidTransport<thermo>&
);

template<class thermo>
Ostream& operator<<
(
    Ostream&,
    const constAnIsoSolidTransport<thermo>&
);


/*---------------------------------------------------------------------------*\
                     Class constAnIsoSolidTransport Declaration
\*---------------------------------------------------------------------------*/

template<class thermo>
class constAnIsoSolidTransport
:
    public thermo
{
    // Private data

        //- Constant anisotropic thermal conductivity.
        vector kappa_;


    // Private Member Functions

        //- Construct from components
        inline constAnIsoSolidTransport
        (
            const thermo& t,
            const vector kappa
        );


public:

    // Constructors

        //- Construct as named copy
        inline constAnIsoSolidTransport
        (
            const word&,
            const constAnIsoSolidTransport&
        );

        //- Construct from Istream
        //constAnIsoSolidTransport(Istream&);
        constAnIsoSolidTransport(const dictionary&);


    // Member functions

        //- Isotropic thermal conductivity [W/mK]
        inline scalar kappa(const scalar T) const;

        //- Un-isotropic thermal conductivity [W/mK]
        inline vector Kappa(const scalar T) const;


    // Member operators

        inline constAnIsoSolidTransport& operator=
        (
            const constAnIsoSolidTransport&
        );
        inline void operator+=(const constAnIsoSolidTransport&);
        inline void operator-=(const constAnIsoSolidTransport&);

         // Friend operators

        friend constAnIsoSolidTransport operator* <thermo>
        (
            const scalar,
            const constAnIsoSolidTransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <thermo>
        (
            Ostream&,
            const constAnIsoSolidTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "constAnIsoSolidTransportI.H"

#ifdef NoRepository
#   include "constAnIsoSolidTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
