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
    Foam::heatTransferCoeffModels::ReynoldsAnalogy

Description
    Heat transfer coefficient calculation based on Reynolds Analogy

    The heat transfer coefficient is derived from the skin friction
    coefficient:

    \f[
        C_f = \frac{\tau_w}{0.5 \rho_\infty |U|^2}
    \f]

    as:

    \f[
        h = 0.5 \rho_\infty \C_{p,\infty} |U_{\infty}| C_f
    \f]

Usage
    Example of function object specification:
    \verbatim
        type            heatTransferCoeff;
        libs            ("libfieldFunctionObjects.so");
        ...
    	htcModel	    ReynoldsAnalogy;
        UInf            (20 0 0);
        Cp              CpInf;
        CpInf           1005;
        ...
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description             | Required    | Default value
        type         | type name: heatTransferCoeff | yes    |
        htcModel     | selected htc model      | yes         |
        UInf         | reference velocity      | yes |
        Cp           | specific heat capacity field name | no |
        rho          | density field name      | no |
    \endtable

    Note:
    - to use a reference \c Cp, set \c Cp to \c CpInf
    - to use a reference \c rho, set \c rho to \c rhoInf

SourceFiles
    ReynoldsAnalogy.C

SeeAlso
    Foam::heatTransferCoeffModel

\*---------------------------------------------------------------------------*/

#ifndef heatTransferCoeffModels_ReynoldsAnalogy_H
#define heatTransferCoeffModels_ReynoldsAnalogy_H

#include "heatTransferCoeffModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace heatTransferCoeffModels
{

/*---------------------------------------------------------------------------*\
                   Class ReynoldsAnalogy Declaration
\*---------------------------------------------------------------------------*/

class ReynoldsAnalogy
:
    public heatTransferCoeffModel
{
    // Private Member Functions

        //- Disallow copy construct
        ReynoldsAnalogy(const ReynoldsAnalogy&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const ReynoldsAnalogy&) = delete;


protected:

    // Protected data

        //- Name of velocity field
        word UName_;

        //- Reference velocity
        vector URef_;

        //- Name of density field
        word rhoName_;

        //- Reference density
        scalar rhoRef_;

        //- Name of specific heat capacity field
        word CpName_;

        //- Reference specific heat capacity
        scalar CpRef_;


    // Protected Member Functions

        virtual tmp<Field<scalar>> rho(const label patchi) const;

        virtual tmp<Field<scalar>> Cp(const label patchi) const;

        virtual tmp<volSymmTensorField> devReff() const;

        tmp<FieldField<Field, scalar>> Cf() const;

        //- Set the heat transfer coefficient
        virtual void htc(volScalarField& htc);


public:

    //- Runtime type information
    TypeName("ReynoldsAnalogy");


    // Constructors

        //- Construct from components
        ReynoldsAnalogy
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const word& TName
        );


    //- Destructor
    virtual ~ReynoldsAnalogy()
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
