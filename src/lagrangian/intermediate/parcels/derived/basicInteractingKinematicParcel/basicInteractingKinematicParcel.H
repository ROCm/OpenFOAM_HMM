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
    Foam::basicInteractingKinematicParcel

Description


SourceFiles
    basicInteractingKinematicParcel.C

\*---------------------------------------------------------------------------*/

#ifndef basicInteractingKinematicParcel_H
#define basicInteractingKinematicParcel_H

#include "InteractingKinematicParcel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class basicInteractingKinematicParcel Declaration
\*---------------------------------------------------------------------------*/

class basicInteractingKinematicParcel
:
    public InteractingKinematicParcel<basicInteractingKinematicParcel>
{

public:

    //- Run-time type information
    TypeName("basicInteractingKinematicParcel");


    // Constructors

        //- Construct from owner, position, and cloud owner
        //  Other properties initialised as null
        basicInteractingKinematicParcel
        (
            InteractingKinematicCloud<basicInteractingKinematicParcel>& owner,
            const vector& position,
            const label cellI
        );

        //- Construct from components
        basicInteractingKinematicParcel
        (
            InteractingKinematicCloud<basicInteractingKinematicParcel>& owner,
            const vector& position,
            const label cellI,
            const label typeId,
            const scalar nParticle0,
            const scalar d0,
            const vector& U0,
            const vector& f0,
            const vector& pi0,
            const vector& tau0,
            const constantProperties& constProps
        );

        //- Construct from Istream
        basicInteractingKinematicParcel
        (
            const Cloud<basicInteractingKinematicParcel>& c,
            Istream& is,
            bool readFields = true
        );

        //- Construct as a copy
        basicInteractingKinematicParcel
        (
            const basicInteractingKinematicParcel& p
        );

        //- Construct and return a clone
        autoPtr<basicInteractingKinematicParcel> clone() const
        {
            return
                autoPtr<basicInteractingKinematicParcel>
                (
                    new basicInteractingKinematicParcel(*this)
                );
        }


    //- Destructor
    virtual ~basicInteractingKinematicParcel();
};


template<>
inline bool contiguous<basicInteractingKinematicParcel>()
{
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
