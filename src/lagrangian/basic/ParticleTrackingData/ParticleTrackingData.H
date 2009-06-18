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
    Foam::ParticleTrackingData

Description
    Class to provide additional properties to allow construction of
    particle tracks

SourceFiles
    ParticleTrackingData.C
    ParticleTrackingDataIO.C

\*---------------------------------------------------------------------------*/

#ifndef ParticleTrackingData_H
#define ParticleTrackingData_H

#include "Cloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes and friend functions
template<class ParcelType>
class ParticleTrackingData;


template<class ParcelType>
Ostream& operator<<
(
    Ostream&,
    const ParticleTrackingData<ParcelType>&
);


/*---------------------------------------------------------------------------*\
                   Class ParticleTrackingData Declaration
\*---------------------------------------------------------------------------*/

template<class ParcelType>
class ParticleTrackingData
{
    // Private data

        //- Reference to the cloud
        const Cloud<ParcelType>& cloud_;

        //- Originating processor id
        label origProc_;

        //- Local particle id
        label id_;

        //- Cumulative particle count used for particle id
        static label PARTICLE_COUNT;


    // Private member functions

        //- Write properties - particle count
        static void writeProperties(const Cloud<ParcelType>& cloud);

        //- Read properties - particle count
        static void readProperties(const Cloud<ParcelType>& cloud);


public:

    // Constructors

        //- Construct from cloud
        ParticleTrackingData(const Cloud<ParcelType>& cloud);

        //- Construct copy
        ParticleTrackingData(const ParticleTrackingData& ptd);

        //- Construct from Istream and mesh
        ParticleTrackingData
        (
            const Cloud<ParcelType>& cloud,
            Istream& is,
            bool readFields
        );


    //- Destructor
    ~ParticleTrackingData();


    // Member functions

        // Access

            //- Return const access to the cloud
            inline const Cloud<ParcelType>& cloud() const;

            //- Return const access to the originating processor id
            inline label origProc() const;

            //- Return const access to the local particle id
            inline label id() const;


    // I-O

        //- Read fields
        static void readFields(Cloud<ParcelType>& c);

        //- Write fields
        static void writeFields(const Cloud<ParcelType>& c);


    // Ostream Operator

        friend Ostream& operator<< <ParcelType>
        (
            Ostream&,
            const ParticleTrackingData<ParcelType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ParticleTrackingDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ParticleTrackingData.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
