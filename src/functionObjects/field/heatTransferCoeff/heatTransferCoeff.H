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
    Foam::functionObjects::heatTransferCoeff

Group
    grpFieldFunctionObjects

Description
    Computes the heat transfer coefficient as a \c volScalarField
    for a set of patches.

    The field is stored on the mesh database so that it can be retrieved and
    used for other applications.

    Operands:
    \table
      Operand           | Type           | Location
      input             | -              | -
      output file       | -              | -
      output field      | volScalarField | $FOAM_CASE/\<time\>/\<outField\>
    \endtable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    heatTransferCoeff1
    {
        // Mandatory entries (unmodifiable)
        type        heatTransferCoeff;
        libs        (fieldFunctionObjects);

        field       T;
        patches     ("walls.*");

        htcModel    ReynoldsAnalogy;
        UInf        (20 0 0);
        Cp          CpInf;
        CpInf       1000;
        rho         rhoInf;
        rhoInf      1.2;
    }
    \endverbatim

    Example usage for mode \c ReynoldsAnalogy for incompressible case:
    \verbatim
    heatTransferCoeff1
    {
        // Mandatory entries (unmodifiable)
        type        heatTransferCoeff;
        libs        (fieldFunctionObjects);

        field       T;
        patches     ("walls.*");

        htcModel    ReynoldsAnalogy;
        UInf        (20 0 0);
        Cp          CpInf;
        CpInf       1000;
        rho         rhoInf;
        rhoInf      1.2;
    }
    \endverbatim

    Example usage for mode \c ReynoldsAnalogy for compressible case:
    \verbatim
    htc
    {
        type        heatTransferCoeff;
        libs        (fieldFunctionObjects);

        field       T;
        patches     ("walls.*");

        htcModel    ReynoldsAnalogy;
        UInf        (20 0 0);
    }
    \endverbatim

    Example usage for mode \c localReferenceTemperature for compressible case:
    \verbatim
    htc
    {
        type        heatTransferCoeff;
        libs        (fieldFunctionObjects);

        field       T;
        patches     ("walls.*");
        htcModel    localReferenceTemperature;
    }
    \endverbatim

    Example usage for mode \c fixedReferenceTemperature for compressible case:
    \verbatim
    htc
    {
        type        heatTransferCoeff;
        libs        (fieldFunctionObjects);

        field       T;
        patches     ("walls.*");
        htcModel    fixedReferenceTemperature;
        TRef        300;
    }
    \endverbatim

    Options for the \c htcModel entry:
    \verbatim
        ReynoldsAnalogy           | Reynold's analogy
        localReferenceTemperature | Local reference temperature
        fixedReferenceTemperature | Specified reference temperature
    \endverbatim

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link fieldExpression.H \endlink

    Usage by the \c postProcess utility is not available.

See also
    - Foam::functionObject
    - Foam::functionObjects::fieldExpression
    - Foam::heatTransferCoeffModels::fixedReferenceTemperature
    - Foam::heatTransferCoeffModels::localReferenceTemperature
    - ExtendedCodeGuide::functionObjects::field::heatTransferCoeff

SourceFiles
    heatTransferCoeff.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_heatTransferCoeff_H
#define functionObjects_heatTransferCoeff_H

#include "fieldExpression.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class heatTransferCoeffModel;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                      Class heatTransferCoeff Declaration
\*---------------------------------------------------------------------------*/

class heatTransferCoeff
:
    public fieldExpression
{
    // Private Data

        //- Heat transfer coefficient model
        autoPtr<heatTransferCoeffModel> htcModelPtr_;


protected:

        //- Calculate the heat transfer coefficient field and return true
        //- if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("heatTransferCoeff");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        heatTransferCoeff
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        heatTransferCoeff(const heatTransferCoeff&) = delete;

        //- No copy assignment
        void operator=(const heatTransferCoeff&) = delete;


    //- Destructor
    virtual ~heatTransferCoeff() = default;


    // Member Functions

        //- Read the heatTransferCoeff data
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
