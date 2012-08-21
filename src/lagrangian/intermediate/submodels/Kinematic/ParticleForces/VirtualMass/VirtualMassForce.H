/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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
    Foam::VirtualMassForce

Description
    Calculates particle virtual mass force

SourceFiles
    VirtualMassForceI.H
    VirtualMassForce.C

\*---------------------------------------------------------------------------*/

#ifndef VirtualMassForce_H
#define VirtualMassForce_H

#include "ParticleForce.H"
#include "volFields.H"
#include "interpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class VirtualMassForce Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class VirtualMassForce
:
    public ParticleForce<CloudType>
{
    // Private data

        //- Name of velocity field
        const word UName_;

        //- Virtual mass coefficient - typically 0.5
        scalar Cvm_;

        //- Rate of change of carrier phase velocity
        volVectorField* DUcDtPtr_;

        //- Rate of change of carrier phase velocity interpolator
        autoPtr<interpolation<vector> > DUcDtInterpPtr_;


public:

    //- Runtime type information
    TypeName("virtualMass");


    // Constructors

        //- Construct from mesh
        VirtualMassForce
        (
            CloudType& owner,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Construct copy
        VirtualMassForce(const VirtualMassForce& pgf);

        //- Construct and return a clone
        virtual autoPtr<ParticleForce<CloudType> > clone() const
        {
            return autoPtr<ParticleForce<CloudType> >
            (
                new VirtualMassForce<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~VirtualMassForce();


    // Member Functions

        // Access

            //- Return the rate of change of carrier phase velocity
            inline const volVectorField& DUcDt() const;

            //- Return the rate of change of carrier phase velocity interpolator
            inline const interpolation<vector>& DUcDtInterp() const;


        // Evaluation

            //- Cache fields
            virtual void cacheFields(const bool store);

            //- Calculate the non-coupled force
            virtual forceSuSp calcCoupled
            (
                const typename CloudType::parcelType& p,
                const scalar dt,
                const scalar mass,
                const scalar Re,
                const scalar muc
            ) const;

            //- Return the added mass
            virtual scalar massAdd
            (
                const typename CloudType::parcelType& p,
                const scalar mass
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "VirtualMassForceI.H"

#ifdef NoRepository
    #include "VirtualMassForce.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
