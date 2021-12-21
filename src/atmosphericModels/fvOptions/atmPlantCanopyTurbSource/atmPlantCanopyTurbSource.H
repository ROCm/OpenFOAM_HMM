/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 ENERCON GmbH
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::fv::atmPlantCanopyTurbSource

Group
    grpFvOptionsSources

Description
    Applies sources on either \c epsilon or \c omega to incorporate effects
    of plant canopy for atmospheric boundary layer modelling.

    Corrections applied to either of the below, if exist:
    \verbatim
      epsilon   | Turbulent kinetic energy dissipation rate [m2/s3]
      omega     | Specific dissipation rate                 [1/s]
    \endverbatim

    Required fields:
    \verbatim
      epsilon/omega   | Dissipation rate OR Spec. dissipation rate [m2/s3]/[1/s]
      plantCd         | Plant canopy drag coefficient              [-]
      leafAreaDensity | Leaf area density                          [1/m]
    \endverbatim

    References:
    \verbatim
        Influence of forest (tag:SP):
            Sogachev, A., & Panferov, O. (2006).
            Modification of two-equation models to account for plant drag.
            Boundary-Layer Meteorology, 121(2), 229-266.
            DOI:10.1007/s10546-006-9073-5
    \endverbatim

Usage
    Example by using \c constant/fvOptions:
    \verbatim
    atmPlantCanopyTurbSource1
    {
        // Mandatory entries (unmodifiable)
        type             atmPlantCanopyTurbSource;

        atmPlantCanopyTurbSourceCoeffs
        {
            // Mandatory (inherited) entries (unmodifiable)
            selectionMode    all;

            // Optional entries (unmodifiable)
            rho          rho;
        }

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                         | Type   | Req'd  | Dflt
      type      | Type name: atmPlantCanopyTurbSource | word   | yes    | -
      rho       | Name of density field               | word   | no     | rho
    \endtable

    The inherited entries are elaborated in:
      - \link fvOption.H \endlink
      - \link cellSetOption.H \endlink

SourceFiles
    atmPlantCanopyTurbSource.C
    atmPlantCanopyTurbSourceTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef fv_atmPlantCanopyTurbSource_H
#define fv_atmPlantCanopyTurbSource_H

#include "cellSetOption.H"
#include "turbulentTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
               Class atmPlantCanopyTurbSource Declaration
\*---------------------------------------------------------------------------*/

class atmPlantCanopyTurbSource
:
    public fv::cellSetOption
{
    // Private Data

        //- Internal flag to determine the working field is epsilon or omega
        Switch isEpsilon_;

        //- Name of density field
        const word rhoName_;

        //- Required turbulence model coefficients (copied from turb model)
        dimensionedScalar Cmu_;
        dimensionedScalar C1_;
        dimensionedScalar C2_;


        // Fields

            //- Plant canopy drag coefficient field [-]
            volScalarField plantCd_;

            //- Leaf area density field [1/m]
            volScalarField leafAreaDensity_;


     // Private Member Functions

        //- Return the modifier for plant canopy effects
        tmp<volScalarField::Internal> calcPlantCanopyTerm
        (
            const volVectorField::Internal& U
        ) const;

        //- Apply atmPlantCanopyTurbSource to epsilon
        template<class AlphaFieldType, class RhoFieldType>
        void atmPlantCanopyTurbSourceEpsilon
        (
            const AlphaFieldType& alpha,
            const RhoFieldType& rho,
            fvMatrix<scalar>& eqn,
            const label fieldi
        ) const;

        //- Apply atmPlantCanopyTurbSource to omega
        template<class AlphaFieldType, class RhoFieldType>
        void atmPlantCanopyTurbSourceOmega
        (
            const AlphaFieldType& alpha,
            const RhoFieldType& rho,
            fvMatrix<scalar>& eqn,
            const label fieldi
        ) const;


public:

    //- Runtime type information
    TypeName("atmPlantCanopyTurbSource");


    // Constructors

        //- Construct from explicit source name and mesh
        atmPlantCanopyTurbSource
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        atmPlantCanopyTurbSource(const atmPlantCanopyTurbSource&) = delete;

        //- No copy assignment
        void operator=(const atmPlantCanopyTurbSource&) = delete;


    // Member Functions

        //- Add explicit contribution to epsilon or omega equation
        //- for incompressible flow computations
         virtual void addSup
        (
            fvMatrix<scalar>& eqn,
            const label fieldi
        );

        //- Add explicit contribution to epsilon or omega equation
        //- for compressible flow computations
        virtual void addSup
        (
            const volScalarField& rho,
            fvMatrix<scalar>& eqn,
            const label fieldi
        );

        //- Add explicit contribution to epsilon or omega equation
        //- for multiphase flow computations
        virtual void addSup
        (
            const volScalarField& alpha,
            const volScalarField& rho,
            fvMatrix<scalar>& eqn,
            const label fieldi
        );

        //- Read source dictionary (effectively no-op)
        virtual bool read(const dictionary& dict)
        {
            return true;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "atmPlantCanopyTurbSourceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
