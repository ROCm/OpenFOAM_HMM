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
    Foam::radiation::blackBodyEmission

Description
    Class Black Body Emission Intensity.

SourceFiles
    blackBodyEmission.C


\*---------------------------------------------------------------------------*/

#ifndef blackModyEmission_H
#define blackModyEmission_H

#include "volFields.H"
#include "dimensionedScalar.H"
#include "mathematicalConstants.H"
#include "radiationConstants.H"
#include "interpolationLookUpTable.H"
#include "Vector2D.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{
/*---------------------------------------------------------------------------*\
                           Class blackBodyEmission Declaration
\*---------------------------------------------------------------------------*/
class blackBodyEmission
{
    // Private data

        mutable interpolationLookUpTable<scalar> blackBodyEmissiveTable_;

        scalar flambdaT(const scalar lambdaT) const;

        const dimensionedScalar C1_;

        const dimensionedScalar C2_;

        // Ptr List for enregy black body emission
        PtrList<volScalarField> bj_;

        // Reference to the temperature field
        const volScalarField& T_;


public:

    // Constructors

        //- Construct from components
        blackBodyEmission
        (
            const fileName& fn,
            const word& instance,
            label lambdaj,
            const volScalarField& T
        );


    // Destructor
    ~blackBodyEmission();


    // Member functions

        // Access

            //- Black body spectrum
            inline const volScalarField& bj(label i) const
            {
                return bj_[i];
            }

            //- Spectral emission for the black body at T and lambda
            inline dimensionedScalar EblambdaT
            (
                const dimensionedScalar T,
                const scalar lambda
            ) const
            {
                return (C1_/(pow5(lambda)*(exp(C2_/(lambda*T)) - 1)));
            }

            //- Integral energy at T from lambda1 to lambda2
            tmp<Foam::volScalarField> EbDeltaLambdaT
            (
                const volScalarField& T,
                const Vector2D<scalar>& band
            ) const;


    // Edit

        // Update black body emission
        void correct(label j, const Vector2D<scalar>& band);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
} // End namespace radiation

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
