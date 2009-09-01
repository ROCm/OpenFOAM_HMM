/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::RanzMarshallBirdCorrection

Description
    The Ranz-Marshall correlation for heat transfer with the Bird correction
    to account for the local shielding effect due to emitted species.

\*---------------------------------------------------------------------------*/

#ifndef RanzMarshallBirdCorrection_H
#define RanzMarshallBirdCorrection_H

#include "RanzMarshall.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                 Class RanzMarshallBirdCorrection Declaration
\*---------------------------------------------------------------------------*/

template <class CloudType>
class RanzMarshallBirdCorrection
:
    public RanzMarshall<CloudType>
{

public:

    //- Runtime type information
    TypeName("RanzMarshallBirdCorrection");


    // Constructors

        //- Construct from dictionary
        RanzMarshallBirdCorrection
        (
            const dictionary& dict,
            CloudType& cloud
        );


    //- Destructor
    virtual ~RanzMarshallBirdCorrection();


    // Member Functions

        //- Return heat transfer coefficient
        virtual scalar htc
        (
            const scalar dp,
            const scalar Re,
            const scalar Pr,
            const scalar kappa,
            const scalar NCpW
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "RanzMarshallBirdCorrection.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
