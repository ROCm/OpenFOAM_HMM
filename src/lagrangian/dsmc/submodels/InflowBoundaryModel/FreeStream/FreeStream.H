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
    Foam::FreeStream

Description
    Inserting new particles across the faces of a specified patch for a free
    stream.  Uniform values of temperature, velocity and number densities

\*---------------------------------------------------------------------------*/

#ifndef FreeStream_H
#define FreeStream_H

#include "InflowBoundaryModel.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                        Class FreeStream Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class FreeStream
:
    public InflowBoundaryModel<CloudType>
{

    // Private data

        //- Index of patch to introduce particles across
        label patchIndex_;

        //- Temperature of the free stream
        scalar temperature_;

        //- Velocity of the free stream
        vector velocity_;

        //- The molecule types to be introduced
        List<label> moleculeTypeIds_;

        //- The number density of the species in the inflow
        Field<scalar> numberDensities_;

        //- A List of Fields, one Field for every species to be introduced, each
        // field entry corresponding to a face on the patch to be injected
        // across.
        List<Field<scalar> > particleFluxAccumulators_;

public:

    //- Runtime type information
    TypeName("FreeStream");


    // Constructors

        //- Construct from dictionary
        FreeStream
        (
            const dictionary& dict,
            CloudType& cloud
        );


    // Destructor
    virtual ~FreeStream();


    // Member Functions

        //- Introduce particles
        virtual void inflow();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "FreeStream.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
