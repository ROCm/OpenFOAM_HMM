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
    Foam::liquidEvaporation

Description
    Liquid evaporation model

\*---------------------------------------------------------------------------*/

#ifndef liquidEvaporation_H
#define liquidEvaporation_H

#include "PhaseChangeModel.H"
#include "liquidMixture.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                       Class liquidEvaporation Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class liquidEvaporation
:
    public PhaseChangeModel<CloudType>
{
protected:

    // Protected data

        //- Global liquid properties data
        autoPtr<liquidMixture> liquids_;

        //- Vaporisation temperature [K]
        //  Local droplet temperature must exceed Tvap before evaporation is
        //  allowed
        scalar Tvap_;

        //- Component diffusion coefficients - one per specie
        //  Note: need to be in same order as defined by liquids in thermo
        //        props
        scalarList Dm_;

        //- Component saturation pressure vs temperature
        List<List<Tuple2<scalar, scalar> > > pSat_;


    // Protected member functions

        //- Sherwood number as a function of Reynolds and Schmidt numbers
        scalar Sh(const scalar Re, const scalar Sc) const;

        //- Return the saturation pressure for species i at temperature, T
        scalar pSat(const label i, const scalar T) const;


public:

    //- Runtime type information
    TypeName("liquidEvaporation");


    // Constructors

        //- Construct from dictionary
        liquidEvaporation
        (
            const dictionary& dict,
            CloudType& cloud
        );


    //- Destructor
    virtual ~liquidEvaporation();


    // Member Functions

        //- Flag to indicate whether model activates phase change model
        bool active() const;

        //- Update model
        scalar calculate
        (
            const scalar T,
            const scalar d,
            const scalarField& Xc,
            const scalarField& dMassMT,
            const vector Ur,
            const scalar Tc,
            const scalar pc,
            const scalar nuc,
            const scalar dt
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "liquidEvaporation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
