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
    Foam::ManualInjection

Description
    Manual injection

    For manual injection,
    - User specifies
      -# Total mass to inject
      -# Parcel positions in file \<positionsFile\>
      -# Initial parcel velocity
    - Parcel diameters obtained by PDF model
    - All parcels introduced at the start of the calculation

    NOTE - not suitable for 2-D slab/wedge simulations unless the positions
           file describes 2-D data!

SourceFiles
    ManualInjection.C

\*---------------------------------------------------------------------------*/

#ifndef manualInjection_H
#define manualInjection_H

#include "InjectionModel.H"
#include "pdf.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class ManualInjection Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ManualInjection
:
    public InjectionModel<CloudType>
{

    // Private data

        //- Coefficients dictionary
        dictionary coeffDict_;

        //- Time to introduce parcels
        const scalar injectionTime_;

        //- Name of file containing positions data
        const word positionsFile_;

        //- Field of parcel positions
        vectorIOField positions_;

        //- Field of parcel diameters
        scalarField diameters_;

        //- Initial parcel velocity
        const vector U0_;

        //- Total volume of parcels to introduce [m^3]
        scalar volumeTotal_;

        //- Parcel size PDF model
        const autoPtr<pdf> parcelPDF_;

        //- Number of particles represented by each parcel
        scalar nParticlesPerParcel_;


public:

    //- Runtime type information
    TypeName("ManualInjection");


    // Constructors

        //- Construct from dictionary
        ManualInjection
        (
            const dictionary& dict,
            CloudType& owner
        );


    // Destructor

        ~ManualInjection();


    // Member Functions

        bool active() const;

        scalar timeStart() const;

        scalar timeEnd() const;

        vector position
        (
            const label iParcel,
            const scalar time,
            const polyMeshInfo& meshInfo,
            Random&
        ) const;

        label nParcelsToInject
        (
            const label,
            const scalar time0,
            const scalar time1
        ) const;

        scalar volume
        (
            const scalar,
            const scalar,
            const polyMeshInfo&
        ) const;

        scalar volumeFraction
        (
            const scalar time0,
            const scalar time1
        ) const;

        scalar d0
        (
            const label iParcel,
            const scalar
        ) const;

        vector velocity
        (
            const label,
            const scalar
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ManualInjection.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
