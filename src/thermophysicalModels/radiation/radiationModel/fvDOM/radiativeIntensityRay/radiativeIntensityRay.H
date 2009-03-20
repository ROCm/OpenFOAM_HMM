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
    Foam::radiation::radiativeIntensityRay

Description

SourceFiles
    radiativeIntensityRay.C

\*---------------------------------------------------------------------------*/

#ifndef radiativeIntensityRay_H
#define radiativeIntensityRay_H

#include "absorptionEmissionModel.H"
#include "blackBodyEmission.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{


/*---------------------------------------------------------------------------*\
                           Class radiativeIntensityRay Declaration
\*---------------------------------------------------------------------------*/
class fvDOM;


class radiativeIntensityRay
{

    // Private data

        //- absorptionEmissionModel
        const absorptionEmissionModel& absEmmModel_;

        //- Temperature
        const fvMesh& mesh_;

        // black Body
        const blackBodyEmission&  blackBody_;

        //- Total Radiative Intensity on i direction  / [W/m2]
        volScalarField I_;

        // - Total  radiative heat flux on boundary on i direction
        mutable volScalarField Qri_;

        //- Direction i
        vector Si_;

        //- theta angle of direction i
        scalar thetai_;

        //- phi angle of direction i
        scalar phii_;

        //- solid angle
        scalar omegai_;

        //- Number of bands on i direction
        label Nlambdaj_;

        //- average vector inside the solid angle
        vector Di_;

        //- List of Pointers to Radiative Intensity wave-length
        PtrList<volScalarField> Ilambdaj_;

    // Private member functions

        //- Disallow default bitwise copy construct
        radiativeIntensityRay(const radiativeIntensityRay&);

        //- Disallow default bitwise assignment
        void operator=(const radiativeIntensityRay&);

        //- Integrate Intensity on this direction
//        void IntegrateRadiativeIntensity(const label i);
public:

    // Static data members
        static label rayId;

    // Constructors

        //- Null constructor
        radiativeIntensityRay
        (
            scalar& phii,
            scalar& thetai,
            scalar& deltaPhi,
            scalar& deltaTheta,
            label& lambdaj,
            const fvMesh& mesh_,
            const absorptionEmissionModel& absEmmModel_,
            const blackBodyEmission& blackBody
        );


    // Destructor

        ~radiativeIntensityRay();


    // Member functions

        // Edit

            // Update radiative intensity on i direction
            scalar correct(fvDOM*);

            // init the ray on i direction
            void init
            (
                scalar phii, scalar thetai, scalar deltaPhi,scalar
                deltaTheta, scalar lambda
            );

            // add radiative heat flux on walls from the boundary patch
            void add(const scalarField&, label) const;

            // add Radiative intensities from all the bands
            void addIntensity();

        // Access

            //- Return Intensity on i direction
            inline const volScalarField& I() const
            {
                return I_;
            }

            //- Return heat flux on boundary on i direction
            inline const volScalarField& Qri() const
            {
                return Qri_;
            }

            inline const vector Si() const
            {
                return Si_;
            }

            inline const vector Di() const
            {
                return Di_;
            }

            scalar lambdaj() const
            {
                return  Nlambdaj_;
            }

            scalar phii() const
            {
                return  phii_;
            }

            scalar thetai() const
            {
                return  thetai_;
            }

            scalar omegai() const
            {
                return  omegai_;
            }


            const volScalarField& Ilambdaj(label i) const
            {
                return  Ilambdaj_[i];
            }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
