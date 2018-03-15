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
    Foam::heatTransferCoeffModels::localReferenceTemperature

Description
    Heat transfer coefficient calculation that employs the patch internal
    field as the reference temperature

    The heat transfer coefficient is specified by:
    
    \f[
        h = \frac{q}{T_c - T_w}
    \f]

Usage
    Example of function object specification:
    \verbatim
        type            heatTransferCoeff;
        libs            ("libfieldFunctionObjects.so");
        ...
    	htcModel	    localReferenceTemperature;
        ...
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description             | Required    | Default value
        type         | type name: heatTransferCoeff | yes    |
        htcModel     | selected htc model      | yes         |
    \endtable

SourceFiles
    localReferenceTemperature.C

SeeAlso
    Foam::heatTransferCoeffModel

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
    // Private Member Functions

        //- Disallow copy construct
        localReferenceTemperature(const localReferenceTemperature&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const localReferenceTemperature&) = delete;


protected:

    // Protected Member Functions

        //- Set the heat transfer coefficient
        virtual void htc(volScalarField& htc);


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
    virtual ~localReferenceTemperature()
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
