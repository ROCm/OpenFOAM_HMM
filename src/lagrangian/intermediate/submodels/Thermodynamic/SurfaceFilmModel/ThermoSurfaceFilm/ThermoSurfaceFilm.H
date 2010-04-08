/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
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
    Foam::ThermoSurfaceFilm

Description
    Thermo parcel surface film model.

SourceFiles
    ThermoSurfaceFilm.C
    ThermoSurfaceFilmI.H

\*---------------------------------------------------------------------------*/

#ifndef ThermoSurfaceFilm_H
#define ThermoSurfaceFilm_H

#include "SurfaceFilmModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class ThermoSurfaceFilm Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ThermoSurfaceFilm
:
    public SurfaceFilmModel<CloudType>
{
protected:

    // Protected data

        //- Convenience typedef to the cloud's parcel type
        typedef typename CloudType::parcelType parcelType;


       // Cached injector fields per film patch

            //- Film temperature / patch face
            scalarList TFilmPatch_;

            //- Film specific heat capacity / patch face
            scalarList cpFilmPatch_;


    // Protected functions

        //- Cache the film fields in preparation for injection
        virtual void cacheFilmFields
        (
            const label filmPatchI,
            const mapDistribute& distMap,
            const surfaceFilmModels::surfaceFilmModel& filmModel
        );

        //- Set the individual parcel properties
        virtual void setParcelProperties
        (
            parcelType& p,
            const label filmFaceI
        ) const;


public:

    //- Runtime type information
    TypeName("ThermoSurfaceFilm");


    // Constructors

        //- Construct from components
        ThermoSurfaceFilm
        (
            const dictionary& dict,
            CloudType& owner,
            const dimensionedVector& g
        );


    //- Destructor
    virtual ~ThermoSurfaceFilm();


    // Member Functions

        // Evaluation

            //- Flag to indicate whether model activates surface film model
            virtual bool active() const;

            //- Transfer parcel from cloud to surface film
            //  Returns true if parcel is to be transferred
            virtual bool transferParcel
            (
                const parcelType& p,
                const label patchI
            );


        // I-O

            //- Write surface film info to stream
            virtual void info(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ThermoSurfaceFilm.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
