/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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
    Foam::LESModels::kOmegaSSTDES

Group
    grpDESTurbulence

Description
    k-omega-SST DES turbulence model for incompressible and compressible flows

    Reference:
    \verbatim
        Strelets, M. (2001)
        Detached Eddy Simulation of Massively Separated Flows,
        39th AIAA Aerospace Sciences Meeting and Exhibit, Reno, NV
    \endverbatim

SourceFiles
    kOmegaSSTDES.C

\*---------------------------------------------------------------------------*/

#ifndef kOmegaSSTDES_H
#define kOmegaSSTDES_H

#include "DESModel.H"
#include "kOmegaSSTBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                        class kOmegaSSTDES Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class kOmegaSSTDES
:
    public kOmegaSSTBase<DESModel<BasicTurbulenceModel> >
{
    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        kOmegaSSTDES(const kOmegaSSTDES&);
        kOmegaSSTDES& operator=(const kOmegaSSTDES&);


protected:

    // Protected data

        // Model coefficients

            dimensionedScalar kappa_;
            dimensionedScalar CDESkom_;
            dimensionedScalar CDESkeps_;


    // Protected Member Functions

        //- Blending for CDES parameter
        virtual tmp<volScalarField> CDES(const volScalarField& F1) const
        {
            return this->blend(F1, CDESkom_, CDESkeps_);
        }

        virtual void correctNut(const volScalarField& S2);
        virtual void correctNut();

        //- Length scale
        virtual tmp<volScalarField> dTilda
        (
            const volScalarField& magGradU,
            const volScalarField& CDES
        ) const;


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("kOmegaSSTDES");


    // Constructors

        //- Construct from components
        kOmegaSSTDES
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~kOmegaSSTDES()
    {}


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();

        //- Return the LES field indicator
        virtual tmp<volScalarField> LESRegion() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
#   include "kOmegaSSTDES.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //
