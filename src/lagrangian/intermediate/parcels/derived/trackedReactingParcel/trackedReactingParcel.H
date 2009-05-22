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
    Foam::trackedReactingParcel

Description


SourceFiles
    trackedReactingParcel.C
    trackedReactingParcelIO.C

\*---------------------------------------------------------------------------*/

#ifndef trackedReactingParcel_H
#define trackedReactingParcel_H

#include "TrackedReactingParcel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class trackedReactingParcel Declaration
\*---------------------------------------------------------------------------*/

class trackedReactingParcel
:
    public TrackedReactingParcel<trackedReactingParcel>
{

public:

    //- Run-time type information
    TypeName("trackedReactingParcel");

    // Constructors

        //- Construct from components
        trackedReactingParcel
        (
            ReactingCloud<trackedReactingParcel>& owner,
            const vector& position,
            const label cellI,
            const label typeId,
            const scalar nParticle0,
            const scalar d0,
            const vector& U0,
            const scalarField& Y0,
            const constantProperties& constProps
        );

        //- Construct from Istream
        trackedReactingParcel
        (
            const Cloud<trackedReactingParcel>& c,
            Istream& is,
            bool readFields = true
        );

        //- Construct as a copy
        trackedReactingParcel(const trackedReactingParcel& p);

        //- Construct and return a clone
        autoPtr<trackedReactingParcel> clone() const
        {
            return
                autoPtr<trackedReactingParcel>
                (
                    new trackedReactingParcel(*this)
                );
        }


    //- Destructor
    virtual ~trackedReactingParcel();
};


template<>
inline bool contiguous<trackedReactingParcel>()
{
    return false; // Now have scalar lists/fields (mass fractions)
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
