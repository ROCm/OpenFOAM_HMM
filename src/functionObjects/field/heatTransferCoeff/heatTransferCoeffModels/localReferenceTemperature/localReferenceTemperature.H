/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2020 OpenCFD Ltd.
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
    Foam::heatTransferCoeffModels::localReferenceTemperature

Description
    Heat transfer coefficient calculation that employs
    the patch internal field as the reference temperature.

    The heat transfer coefficient is specified by:

    \f[
        h = \frac{q}{T_c - T_w}
    \f]

    where
    \vartable
      h       | Convective heat transfer coefficient of the flow  [W/(m2 K)]
      q       | Heat flux                                         [W/m2]
      T_{ref} | Reference temperature of patch internal field     [K]
      T_w     | Patch temperature                                 [K]
    \endvartable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    heatTransferCoeff1
    {
        // Mandatory and other optional entries
        ...
        htcModel        localReferenceTemperature;
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                            | Type | Reqd | Dflt
      type      | Model name: localReferenceTemperature  | word | yes  | -
    \endtable

See also
  - Foam::heatTransferCoeffModel
  - Foam::heatTransferCoeffModels::fixedReferenceTemperature
  - Foam::heatTransferCoeffModels::ReynoldsAnalogy

SourceFiles
    localReferenceTemperature.C

\*---------------------------------------------------------------------------*/

#ifndef heatTransferCoeffModels_localReferenceTemperature_H
#define heatTransferCoeffModels_localReferenceTemperature_H

#include "heatTransferCoeffModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace heatTransferCoeffModels
{

/*---------------------------------------------------------------------------*\
                   Class localReferenceTemperature Declaration
\*---------------------------------------------------------------------------*/

class localReferenceTemperature
:
    public heatTransferCoeffModel
{
protected:

    // Protected Member Functions

        //- Set the heat transfer coefficient
        virtual void htc
        (
            volScalarField& htc,
            const FieldField<Field, scalar>& q
        );

        //- No copy construct
        localReferenceTemperature(const localReferenceTemperature&) = delete;

        //- No copy assignment
        void operator=(const localReferenceTemperature&) = delete;


public:

    //- Runtime type information
    TypeName("localReferenceTemperature");


    // Constructors

        //- Construct from components
        localReferenceTemperature
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const word& TName
        );


    //- Destructor
    virtual ~localReferenceTemperature() = default;


    // Member Functions

        //- Read from dictionary
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace heatTransferCoeffModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
