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
    Foam::fv::atmPlantCanopyUSource

Group
    grpFvOptionsSources

Description
    Applies sources on velocity, i.e. \c U, to incorporate effects
    of plant canopy for atmospheric boundary layer modelling.

    Corrections applied to:
    \verbatim
      U                | Velocity                               [m/s]
    \endverbatim

    Required fields:
    \verbatim
      U                | Velocity                               [m/s]
      plantCd          | Plant canopy drag coefficient          [-]
      leafAreaDensity  | Leaf area density                      [1/m]
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
    atmPlantCanopyUSource1
    {
        // Mandatory entries (unmodifiable)
        type             atmPlantCanopyUSource;

        atmPlantCanopyUSourceCoeffs
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
      Property  | Description                      | Type   | Req'd  | Dflt
      type      | Type name: atmPlantCanopyUSource | word   | yes    | -
      rho       | Name of density field            | word   | no     | rho
    \endtable

    The inherited entries are elaborated in:
      - \link fvOption.H \endlink
      - \link cellSetOption.H \endlink

SourceFiles
    atmPlantCanopyUSource.C
    atmPlantCanopyUSourceTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef fv_atmPlantCanopyUSource_H
#define fv_atmPlantCanopyUSource_H

#include "cellSetOption.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
               Class atmPlantCanopyUSource Declaration
\*---------------------------------------------------------------------------*/

class atmPlantCanopyUSource
:
    public fv::cellSetOption
{
    // Private Data

        //- Name of density field
        const word rhoName_;


        // Fields

            //- Plant canopy drag coefficient field [-]
            volScalarField plantCd_;

            //- Leaf area density field [1/m]
            volScalarField leafAreaDensity_;


public:

    //- Runtime type information
    TypeName("atmPlantCanopyUSource");


    // Constructors

        //- Construct from explicit source name and mesh
        atmPlantCanopyUSource
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        atmPlantCanopyUSource(const atmPlantCanopyUSource&) = delete;

        //- No copy assignment
        void operator=(const atmPlantCanopyUSource&) = delete;


    // Member Functions

        //- Add explicit contribution to momentum equation
        virtual void addSup
        (
            fvMatrix<vector>& eqn,
            const label fieldi
        );

        //- Add explicit contribution to compressible momentum equation
        virtual void addSup
        (
            const volScalarField& rho,
            fvMatrix<vector>& eqn,
            const label fieldi
        );

        //- Add explicit contribution to phase momentum equation
        virtual void addSup
        (
            const volScalarField& alpha,
            const volScalarField& rho,
            fvMatrix<vector>& eqn,
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

#endif

// ************************************************************************* //
