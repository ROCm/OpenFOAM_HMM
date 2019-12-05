/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::LESModels::kOmegaSSTDDES

Group
    grpDESTurbulence

Description
    k-omega-SST DDES turbulence model for incompressible and compressible flows

    Reference:
    \verbatim
        Gritskevich, M.S., Garbaruk, A.V., Schuetze, J., Menter, F.R. (2011)
        Development of DDES and IDDES Formulations for the k-omega
        Shear Stress Transport Model, Flow, Turbulence and Combustion,
        pp. 1-19
    \endverbatim

SourceFiles
    kOmegaSSTDDES.C

\*---------------------------------------------------------------------------*/

#ifndef kOmegaSSTDDES_H
#define kOmegaSSTDDES_H

#include "kOmegaSSTDES.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                        Class kOmegaSSTDDES Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class kOmegaSSTDDES
:
    public kOmegaSSTDES<BasicTurbulenceModel>
{
    // Private Member Functions

        tmp<volScalarField> fd(const volScalarField& magGradU) const;

        tmp<volScalarField> rd(const volScalarField& magGradU) const;

        //- No copy construct
        kOmegaSSTDDES(const kOmegaSSTDDES&) = delete;

        //- No copy assignment
        void operator=(const kOmegaSSTDDES&) = delete;


protected:

    // Protected data

        // Model coefficients

            dimensionedScalar Cd1_;
            dimensionedScalar Cd2_;


    // Protected Member Functions

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
    TypeName("kOmegaSSTDDES");


    // Constructors

        //- Construct from components
        kOmegaSSTDDES
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
    virtual ~kOmegaSSTDDES() = default;


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "kOmegaSSTDDES.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
