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

// Forward declaration of classes
class fvDOM;

/*---------------------------------------------------------------------------*\
                    Class radiativeIntensityRay Declaration
\*---------------------------------------------------------------------------*/

class radiativeIntensityRay
{
    // Private data

        //- Absorption/emission model
        const absorptionEmissionModel& absEmmModel_;

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Black body
        const blackBodyEmission&  blackBody_;

        //- Total radiative intensity in i direction  / [W/m2]
        volScalarField I_;

        //- Total radiative heat flux on boundary in i direction
        mutable volScalarField Qri_;

        //- Direction i
        vector Si_;

        //- Theta angle of direction i
        scalar thetai_;

        //- Phi angle of direction i
        scalar phii_;

        //- Solid angle of direction i
        scalar omegai_;

        //- Number of bands on i direction
        label nLambdaj_;

        //- Average vector inside the solid angle
        vector Di_;

        //- List of pointers to radiative intensity wave-length
        PtrList<volScalarField> Ilambdaj_;


    // Private member functions

        //- Disallow default bitwise copy construct
        radiativeIntensityRay(const radiativeIntensityRay&);

        //- Disallow default bitwise assignment
        void operator=(const radiativeIntensityRay&);


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

            //- Update radiative intensity on i direction
            scalar correct(fvDOM*);

            //- Initialise the ray in i direction
            void init
            (
                const scalar phii,
                const scalar thetai,
                const scalar deltaPhi,
                const scalar deltaTheta,
                const scalar lambda
            );

            //- Add radiative heat flux on walls from the boundary patch
            void add(const scalarField&, const label) const;

            //- Add Radiative intensities from all the bands
            void addIntensity();


        // Access

            //- Return intensity in i direction
            inline const volScalarField& I() const
            {
                return I_;
            }

            //- Return heat flux on boundary in i direction
            inline const volScalarField& Qri() const
            {
                return Qri_;
            }

            //- Return direction i
            inline const vector Si() const
            {
                return Si_;
            }

            //- Return the average vector inside the solid angle
            inline const vector Di() const
            {
                return Di_;
            }

            //- Return the number of bands on i direction
            scalar nLambdaj() const
            {
                return  nlambdaj_;
            }

            //- Return the phi angle of direction i
            scalar phii() const
            {
                return  phii_;
            }

            //- Return the theta angle of direction i
            scalar thetai() const
            {
                return  thetai_;
            }

            //- Return the solid angle of direction i
            scalar omegai() const
            {
                return  omegai_;
            }

            //- Return the list of pointers to radiative intensity wave-length
            const volScalarField& Ilambdaj(const label i) const
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
