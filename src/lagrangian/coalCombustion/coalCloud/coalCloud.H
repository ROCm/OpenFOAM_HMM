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
    coalCloud

Description

SourceFiles
    coalCloud.C

\*---------------------------------------------------------------------------*/

#ifndef coalCloud_H
#define coalCloud_H

#include "ReactingMultiphaseCloud.H"
#include "coalParcel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class coalCloud Declaration
\*---------------------------------------------------------------------------*/

class coalCloud
:
    public ReactingMultiphaseCloud<coalParcel>
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        coalCloud(const coalCloud&);

        //- Disallow default bitwise assignment
        void operator=(const coalCloud&);


public:

    //-Runtime type information
    TypeName("coalCloud");


    // Constructors

        //- Construct given carrier gas fields
        coalCloud
        (
            const word& cloudType,
            const volScalarField& rho,
            const volVectorField& U,
            const dimensionedVector& g,
            hCombustionThermo& thermo,
            PtrList<specieReactingProperties>& gases
        );


    //- Destructor
    ~coalCloud();


    // Member Functions

        //- Write fields
        void writeFields() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
