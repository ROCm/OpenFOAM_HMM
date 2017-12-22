/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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
    Foam::functionObjects::heatTransferCoeff

Group
    grpFieldFunctionObjects

Description
    This function object calculates and writes the heat transfer coefficient
    as a volScalarField for a set of patches.

    The field is stored on the mesh database so that it can be retrieved and
    used for other applications.  Heat transfer coefficient, htc [W/m2/K]
    can be evaluated using one of the following modes:
    - ReynoldsAnalogy: Reynold's analogy
    - localReferenceTemperature: local reference temperature
    - fixedReferenceTemperature: specified reference temperature

Usage
    Example usage for mode 'ReynoldsAnalogy' for incompressible case
    \verbatim
    htc
    {
        type        heatTransferCoeff;
        libs        ("libfieldFunctionObjects.so");

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

    Example usage for mode 'ReynoldsAnalogy' for compressible case
    \verbatim
    htc
    {
        type        heatTransferCoeff;
        libs        ("libfieldFunctionObjects.so");

        field       T;
        patches     ("walls.*");

        htcModel    ReynoldsAnalogy;
        UInf        (20 0 0);
    }
    \endverbatim

    Example usage for mode 'localReferenceTemperature' for compressible case
    \verbatim
    htc
    {
        type        heatTransferCoeff;
        libs        ("libfieldFunctionObjects.so");

        field       T;
        patches     ("walls.*");
        htcModel    local;
    }
    \endverbatim

    Example usage for mode 'fixedReferenceTemperature' for compressible case
    \verbatim
    htc
    {
        type        heatTransferCoeff;
        libs        ("libfieldFunctionObjects.so");

        field       T;
        patches     ("walls.*");
        htcModel    local;
        TRef        300;
    }
    \endverbatim

See also
    Foam::functionObjects::fieldExpression
    Foam::heatTransferCoeffModels::fixedReferenceTemperature
    Foam::heatTransferCoeffModels::localReferenceTemperature

SourceFiles
    heatTransferCoeff.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_heatTransferCoeff_H
#define functionObjects_heatTransferCoeff_H

#include "fieldExpression.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

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

private:

    // Private data

        //- Heat transfer coefficient model
        autoPtr<heatTransferCoeffModel> htcModelPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        heatTransferCoeff(const heatTransferCoeff&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const heatTransferCoeff&) = delete;


protected:

        //- Calculate the heat transfer coefficient field and return true
        // if successful
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


    //- Destructor
    virtual ~heatTransferCoeff();


    // Member Functions

        //- Read the heatTransferCoeff data
        virtual bool read(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
