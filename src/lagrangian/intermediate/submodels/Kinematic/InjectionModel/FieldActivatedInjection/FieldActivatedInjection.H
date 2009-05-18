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
    Foam::FieldActivatedInjection

Description
    Injection at specified positions, with the condition that for injection
    to be allowed

        referenceField[cellI] >= thresholdField[cellI]

    where:
      - referenceField is the field used to supply the look-up values
      - thresholdField supplies the values beyond which the injection is
        permitted

SourceFiles
    FieldActivatedInjection.C

\*---------------------------------------------------------------------------*/

#ifndef FieldActivatedInjection_H
#define FieldActivatedInjection_H

#include "InjectionModel.H"
#include "pdf.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class FieldActivatedInjection Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class FieldActivatedInjection
:
    public InjectionModel<CloudType>
{
    // Private data

        // Model parameters

            //- Reference field
            const volScalarField& referenceField_;

            //- Threshold field
            const volScalarField& thresholdField_;


        // Injector properties

            //- Name of file containing positions data
            const word positionsFile_;

            //- Field of injector (x,y,z) positions
            vectorIOField positions_;

            //- Field of cell labels corresoponding to injector positions
            labelField injectorCells_;

            //- Number of parcels per injector
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
    TypeName("FieldActivatedInjection");


    // Constructors

        //- Construct from dictionary
        FieldActivatedInjection
        (
            const dictionary& dict,
            CloudType& owner
        );


    // Destructor
    virtual ~FieldActivatedInjection();


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
#   include "FieldActivatedInjection.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
