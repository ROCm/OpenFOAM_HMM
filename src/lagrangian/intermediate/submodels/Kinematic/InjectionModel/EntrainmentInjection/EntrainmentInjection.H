/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
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
    Foam::EntrainmentInjection

Description
    Injection at specified positions, with the condition that for injection
    to be allowed

        rhol <= c_*0.5*rhoc*|Uc|^2

    where:
      - rhol is the effective density (of all clouds)
      - c_ is a modelling constant [s2/m2]
      - rhoc is the carrier phase density
      - Uc is the carrier phase velocity

SourceFiles
    EntrainmentInjection.C

\*---------------------------------------------------------------------------*/

#ifndef EntrainmentInjection_H
#define EntrainmentInjection_H

#include "InjectionModel.H"
#include "pdf.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class EntrainmentInjection Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class EntrainmentInjection
:
    public InjectionModel<CloudType>
{
    // Private data

        // Model parameters

            //- Modelling constant
            const scalar c_;

            //- Carrier phase density field
            const volScalarField& rhoc_;

            //- Lagrangian phase density field
            const scalarField& rhol_;

            //- Carrier phase velocity field
            const vectorField& Uc_;


        // Injector properties

            //- Name of file containing positions data
            const word positionsFile_;

            //- Field of injector (x,y,z) positions
            vectorIOField positions_;

            //- Field of cell labels corresoponding to injector positions
            labelField injectorCells_;

            //- nParcels per injector
            const label nParcelsPerInjector_;

            //- Field of number of parcels injected for each injector
            labelField nParcelsInjected_;


        // Parcel properties

            //- Initial parcel velocity
            const vector U0_;

            //- Field of parcel diameters
            scalarField diameters_;

            //- Parcel size PDF model
            const autoPtr<pdf> parcelPDF_;


protected:

    // Protected member functions

        //- Number of parcels to introduce over the time step
        label parcelsToInject
        (
            const scalar time0,
            const scalar time1
        ) const;

        //- Volume of parcels to introduce over the time step
        scalar volumeToInject
        (
            const scalar time0,
            const scalar time1
        ) const;


public:

    //- Runtime type information
    TypeName("EntrainmentInjection");


    // Constructors

        //- Construct from dictionary
        EntrainmentInjection
        (
            const dictionary& dict,
            CloudType& owner
        );


    // Destructor
    virtual ~EntrainmentInjection();


    // Member Functions

        //- Flag to indicate whether model activates injection model
        bool active() const;

        //- Return the end-of-injection time
        scalar timeEnd() const;


        // Injection geometry

            //- Set the injection position and owner cell
            void setPositionAndCell
            (
                const label parcelI,
                const scalar time,
                vector& position,
                label& cellOwner
            );

            //- Return the velocity of the parcel to introduce at a time
            vector velocity
            (
                const label parcelI,
                const scalar time
            );

            //- Return the diameter of the parcel to introduce at a time
            scalar d0
            (
                const label parcelI,
                const scalar time
            ) const;

            //- Return flag to identify whether or not injection in cellI is
            //  permitted
            bool validInjection(const label parcelI);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "EntrainmentInjection.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
