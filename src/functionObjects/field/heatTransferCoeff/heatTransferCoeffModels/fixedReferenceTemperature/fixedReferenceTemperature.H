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
    Foam::heatTransferCoeffModels::fixedReferenceTemperature

Description
    Heat transfer coefficient calculation
    that employs a fixed reference temperature.

    The heat transfer coefficient is specified by:

    \f[
        h = \frac{q}{T_{ref} - T_w}
    \f]

    where
    \vartable
      h       | Convective heat transfer coefficient of the flow  [W/(m2 K)]
      q       | Heat flux                                         [W/m2]
      T_{ref} | Reference temperature of surrounding fluid        [K]
      T_w     | Patch temperature                                 [K]
    \endvartable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    heatTransferCoeff1
    {
        // Mandatory and other optional entries
        ...
        htcModel        fixedReferenceTemperature;

        // Conditional mandatory entries (runtime modifiable)
        Tref            0;
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                            | Type | Reqd | Dflt
      type      | Model name: fixedReferenceTemperature  | word | yes  | -
      Tref      | Reference temperature of surrounding fluid | scalar | yes | -
    \endtable

See also
  - Foam::heatTransferCoeffModel
  - Foam::heatTransferCoeffModels::localReferenceTemperature
  - Foam::heatTransferCoeffModels::ReynoldsAnalogy

SourceFiles
    fixedReferenceTemperature.C

\*---------------------------------------------------------------------------*/

#ifndef heatTransferCoeffModels_fixedReferenceTemperature_H
#define heatTransferCoeffModels_fixedReferenceTemperature_H

#include "heatTransferCoeffModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace heatTransferCoeffModels
{

/*---------------------------------------------------------------------------*\
                   Class fixedReferenceTemperature Declaration
\*---------------------------------------------------------------------------*/

class fixedReferenceTemperature
:
    public heatTransferCoeffModel
{
protected:

    // Protected Data

        //- Reference temperature
        scalar TRef_;


    // Protected Member Functions

        //- Set the heat transfer coefficient
        virtual void htc
        (
            volScalarField& htc,
            const FieldField<Field, scalar>& q
        );


public:

    //- Runtime type information
    TypeName("fixedReferenceTemperature");


    // Constructors

        //- Construct from components
        fixedReferenceTemperature
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const word& TName
        );

        //- No copy construct
        fixedReferenceTemperature(const fixedReferenceTemperature&) = delete;

        //- No copy assignment
        void operator=(const fixedReferenceTemperature&) = delete;


    //- Destructor
    virtual ~fixedReferenceTemperature() = default;


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
