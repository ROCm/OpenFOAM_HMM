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
    Foam::heatTransferCoeffModels::fixedReferenceTemperature

Description
    Heat transfer coefficient calculation that employs a fixed reference
    temperature

    The heat transfer coefficient is specified by:
    
    \f[
        h = \frac{q}{T_{ref} - T_w}
    \f]

Usage
    Example of function object specification:
    \verbatim
        type            heatTransferCoeff;
        libs            ("libfieldFunctionObjects.so");
        ...
    	htcModel	    fixedReferenceTemperature;
	    TRef		    300;
        ...
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description             | Required    | Default value
        type         | type name: heatTransferCoeff | yes    |
        htcModel     | selected htc model      | yes         |
        TRef         | reference temperature   | yes |
    \endtable

SourceFiles
    fixedReferenceTemperature.C

SeeAlso
    Foam::heatTransferCoeffModel

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
    // Private Member Functions

        //- Disallow copy construct
        fixedReferenceTemperature(const fixedReferenceTemperature&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const fixedReferenceTemperature&) = delete;


protected:

    // Protected data

        //- Reference tempetaure
        scalar TRef_;


    // Protected Member Functions

        //- Set the heat transfer coefficient
        virtual void htc(volScalarField& htc);


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


    //- Destructor
    virtual ~fixedReferenceTemperature()
    {}


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
