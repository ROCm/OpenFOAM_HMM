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
    Foam::trackedReactingCloud

Description
    Reacting cloud templated on the reacting parcel

SourceFiles
    trackedReactingCloud.C

\*---------------------------------------------------------------------------*/

#ifndef trackedReactingCloud_H
#define trackedReactingCloud_H

#include "ReactingCloud.H"
#include "trackedReactingParcel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class trackedReactingCloud Declaration
\*---------------------------------------------------------------------------*/

class trackedReactingCloud
:
    public ReactingCloud<trackedReactingParcel>
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        trackedReactingCloud(const trackedReactingCloud&);

        //- Disallow default bitwise assignment
        void operator=(const trackedReactingCloud&);


public:

    //- Runtime type information
    TypeName("trackedReactingCloud");


    // Constructors

        //- Construct given carrier gas fields
        trackedReactingCloud
        (
            const word& cloudName,
            const volScalarField& rho,
            const volVectorField& U,
            const dimensionedVector& g,
            hCombustionThermo& thermo,
            PtrList<specieReactingProperties>& gases
        );


    //- Destructor
    ~trackedReactingCloud();


    // Member Functions

        //- Write fields
        void writeFields() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
