/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2011 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::SprayParcel

Description
    Derived from ReacingParcel, with added functionality
    for atomization and breakup

\*---------------------------------------------------------------------------*/

#ifndef SprayParcel_H
#define SprayParcel_H

#include "ReactingParcel.H"
#include "SprayCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class ParcelType>
class SprayParcel;

template<class ParcelType>
Ostream& operator<<
(
    Ostream&,
    const SprayParcel<ParcelType>&
);

/*---------------------------------------------------------------------------*\
                          Class SprayParcel Declaration
\*---------------------------------------------------------------------------*/

template<class ParcelType>
class SprayParcel
:
    public ReactingParcel<ParcelType>
{
protected:

    // protected properties

        // spray parcel properties

            //- initial droplet diameter
            scalar d0_;

            //- injection position
            vector position0_;

            //- part of liquid core ( > 0.5-liquid, < 0.5-droplet )
            scalar liquidCore_;

            //- index for KH Breakup
            scalar KHindex_;

            //- spherical deviation
            scalar y_;

            //- rate of change of spherical deviation
            scalar yDot_;

            //- characteristic time (used in atomization and/or breakup model)
            scalar tc_;

            //- stripped parcel mass due to breakup
            scalar ms_;

            //- injected from injector (needed for ex. for calculating distance from injector)
            scalar injector_;

            //- momentum relaxation time (needed for calculating parcel acc.)
            scalar tMom_;

            //- passive scalar (extra variable to be defined by user)
            scalar user_;

public:

    //- Class to hold reacting particle constant properties
    class constantProperties
    :
        public ReactingParcel<ParcelType>::constantProperties
    {
    public:

        //- Constructor
        constantProperties(const dictionary& parentDict);
   };


    //- Class used to pass reacting tracking data to the trackToFace function
    class trackData
    :
        public ReactingParcel<ParcelType>::trackData
    {

        // Private data

            //- Reference to the cloud containing this particle
            SprayCloud<ParcelType>& cloud_;

            //- Particle constant properties
            const constantProperties& constProps_;

    public:

        typedef SprayCloud<ParcelType> cloudType;


        // Constructors

            //- Construct from components
            inline trackData
            (
                SprayCloud<ParcelType>& cloud,
                const constantProperties& constProps,
                const interpolation<scalar>& rhoInterp,
                const interpolation<vector>& UInterp,
                const interpolation<scalar>& muInterp,
                const interpolation<scalar>& TInterp,
                const interpolation<scalar>& CpInterp,
                const interpolation<scalar>& pInterp,
                const vector& g
            );


        // Member functions

            //- Return access to the owner cloud
            inline SprayCloud<ParcelType>& cloud();

            //- Return const access to the constant properties
            inline const constantProperties& constProps() const;

    };


    // Static data members

        //- String representation of properties
        static string propHeader;

        //- Runtime type information
        TypeName("SprayParcel");


    friend class Cloud<ParcelType>;


    // Constructors

        //- Construct from owner, position, and cloud owner
        //  Other properties initialised as null
        inline SprayParcel
        (
            SprayCloud<ParcelType>& owner,
            const vector& position,
            const label cellI
        );

        //- Construct from components
        inline SprayParcel
        (
            SprayCloud<ParcelType>& owner,
            const vector& position,
            const label cellI,
            const label typeId,
            const scalar nParticle0,
            const scalar d0,
            const vector& U0,
            const scalarField& Y0,
            const scalar& dInit,
            const scalar& liquidCore,
            const scalar& KHindex,
            const scalar& y,
            const scalar& yDot,
            const scalar& tc,
            const scalar& ms,
            const scalar& injector,
            const scalar& tMom,
            const scalar& user,
            const constantProperties& constProps
        );

        //- Construct from Istream
        SprayParcel
        (
            const Cloud<ParcelType>& c,
            Istream& is,
            bool readFields = true
        );

        //- Construct as a copy
        SprayParcel(const SprayParcel& p);

        //- Construct and return a clone
        autoPtr<SprayParcel> clone() const
        {
            return autoPtr<SprayParcel>(new SprayParcel(*this));
        }

    // Member Functions

        // Access

            //- Return const access to initial droplet diameter
            inline const scalar& d0() const;

            //- Return const access to initial droplet position
            inline const vector& position0() const;

            //- Return const access to liquid core
            inline const scalar& liquidCore() const;

            //- Return const access to Kelvin-Helmholtz breakup index
            inline const scalar& KHindex() const;

            //- Return const access to spherical deviation
            inline const scalar& y() const;

            //- Return const access to rate of change of spherical deviation
            inline const scalar& yDot() const;

            //- Return const access to atomization characteristic time
            inline const scalar& tc() const;

            //- Return const access to stripped parcel mass
            inline const scalar& ms() const;

            //- Return const access to injector id
            inline const scalar& injector() const;

            //- Return const access to momentum relaxation time
            inline const scalar& tMom() const;

            //- Return const access to passive user scalar
            inline const scalar& user() const;

        // Edit

            //- Return access to initial droplet diameter
            inline scalar& d0();

            //- Return access to initial droplet position
            inline vector& position0();

            //- Return access to liquid core
            inline scalar& liquidCore();

            //- Return access to Kelvin-Helmholtz breakup index
            inline scalar& KHindex();

            //- Return access to spherical deviation
            inline scalar& y();

            //- Return access to rate of change of spherical deviation
            inline scalar& yDot();

            //- Return access to atomization characteristic time
            inline scalar& tc();

            //- Return access to stripped parcel mass
            inline scalar& ms();

            //- Return access to injector id
            inline scalar& injector();

            //- Return access to momentum relaxation time
            inline scalar& tMom();

            //- Return access to passive user scalar
            inline scalar& user();


        // Main calculation loop

            //- Set cell values
            template<class TrackData>
            void setCellValues
            (
                TrackData& td,
                const scalar dt,
                const label cellI
            );

            //- Correct parcel properties according to atomization model
            template<class TrackData>
            void calcAtomization
            (
                TrackData& td,
                const scalar dt,
                const label cellI
            );

            //- Correct parcel properties according to breakup model
            template<class TrackData>
            void calcBreakup
            (
                TrackData& td,
                const scalar dt,
                const label cellI
            );

            //- Correct cell values using latest transfer information
            template<class TrackData>
            void cellValueSourceCorrection
            (
                TrackData& td,
                const scalar dt,
                const label cellI
            );

            //- Correct surface values due to emitted species
            template<class TrackData>
            void correctSurfaceValues
            (
                TrackData& td,
                const label cellI,
                const scalar T,
                const scalarField& Cs,
                scalar& rhos,
                scalar& mus,
                scalar& Pr,
                scalar& kappa
            );

            //- Update parcel properties over the time interval
            template<class TrackData>
            void calc
            (
                TrackData& td,
                const scalar dt,
                const label cellI
            );

            //- Calculate the chi-factor for flash-boiling for the atomization model
            template<class TrackData>
            scalar chi
            (
                TrackData& td,
                const scalarField& X
            ) const;


            //- Solve the TAB equation
            template<class TrackData>
            void solveTABEq
            (
                TrackData& td,
                const scalar& dt
            );

        // I-O

            //- Read
            static void readFields(Cloud<ParcelType>& c);

            //- Write
            static void writeFields(const Cloud<ParcelType>& c);


    // Ostream Operator

        friend Ostream& operator<< <ParcelType>
        (
            Ostream&,
            const SprayParcel<ParcelType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "SprayParcelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SprayParcel.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
