/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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
    coalParcel

Description


SourceFiles
    coalParcel.C
    coalParcelIO.C

\*---------------------------------------------------------------------------*/

#ifndef coalParcel_H
#define coalParcel_H

#include "ReactingMultiphaseParcel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class coalParcel Declaration
\*---------------------------------------------------------------------------*/

class coalParcel
:
    public ReactingMultiphaseParcel<coalParcel>
{
public:

    // Run-time type information
    TypeName("coalParcel");

    // Constructors

        //- Construct from owner, position, and cloud owner
        //  Other properties initialised as null
        coalParcel
        (
             ReactingMultiphaseCloud<coalParcel>& owner,
             const vector& position,
             const label cellI
        );

        //- Construct from components
        coalParcel
        (
             ReactingMultiphaseCloud<coalParcel>& owner,
             const vector& position,
             const label cellI,
             const label typeId,
             const scalar nParticle0,
             const scalar d0,
             const vector& U0,
             const scalarField& YMixture0,
             const scalarField& YGas0,
             const scalarField& YLiquid0,
             const scalarField& YSolid0,
             const constantProperties& constProps
        );

        //- Construct from Istream
        coalParcel
        (
            const Cloud<coalParcel>& c,
            Istream& is,
            bool readFields = true
        );

        //- Construct and return a clone
        autoPtr<coalParcel> clone() const
        {
            return autoPtr<coalParcel>(new coalParcel(*this));
        }


    //- Destructor
    virtual ~coalParcel();
};


template<>
inline bool contiguous<coalParcel>()
{
    return false; // Now have scalar lists/fields (mass fractions)
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
