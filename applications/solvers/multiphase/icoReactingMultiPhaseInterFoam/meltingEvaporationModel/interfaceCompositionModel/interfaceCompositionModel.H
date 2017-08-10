/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenFOAM Foundation
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
    Foam::interfaceCompositionModel

Description
    Generic base class for interface composition models. These models describe
    the composition in phase 1 of the supplied pair at the interface with phase
    2.

SourceFiles
    interfaceCompositionModel.C

\*---------------------------------------------------------------------------*/

#ifndef interfaceCompositionModel_H
#define interfaceCompositionModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "dictionary.H"
#include "hashedWordList.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

class phaseModel;
class phasePair;

/*---------------------------------------------------------------------------*\
                        Class interfaceCompositionModel Declaration
\*---------------------------------------------------------------------------*/

class interfaceCompositionModel
{
public:

    // Public type

        //- Enumeration for variable based mass transfer models
        enum modelVariable
        {
            T,  /* temperature based */
            P,  /* pressure based */
            Y   /* mass fraction based */
        };

        static const NamedEnum<modelVariable, 3> modelVariableNames;

        //- Enumeration for model variables
        modelVariable modelVariable_;

        //- Semi-Implicit or Explicit model
        bool semiImplicit_;


protected:

    // Protected data

        //- Phase pair
        const phasePair& pair_;

        //- Names of the transferring species
        const hashedWordList speciesNames_;

        //- Reference to mesh
        const fvMesh& mesh_;


public:

    //- Runtime type information
    TypeName("interfaceCompositionModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            interfaceCompositionModel,
            dictionary,
            (
                const dictionary& dict,
                const phasePair& pair
            ),
            (dict, pair)
        );


    // Constructors

        //- Construct from a dictionary and a phase pair
        interfaceCompositionModel
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~interfaceCompositionModel();


    // Selectors

        static autoPtr<interfaceCompositionModel> New
        (
            const dictionary& dict,
            const phasePair& pair
        );


    // Member Functions


        //- Return the transferring species names
        const hashedWordList& species() const;

        //- Return pair
        const phasePair& pair() const;

        //- Returns whether the species is transported by the model and
        //  provides the name of the diffused species
        bool transports(word& speciesName) const;

        //- Interface mass fraction
        virtual tmp<volScalarField> Yf
        (
            const word& speciesName,
            const volScalarField& Tf
        ) const = 0;

        //- Mass fraction difference between the interface and the field
        virtual tmp<volScalarField> dY
        (
            const word& speciesName,
            const volScalarField& Tf
        ) const = 0;

        //- Mass diffusivity
        virtual tmp<volScalarField> D
        (
            const word& speciesName
        ) const = 0;

        //- Latent heat (delta Hc)
        virtual tmp<volScalarField> L
        (
            const word& speciesName,
            const volScalarField& Tf
        ) const = 0;

        //- Heat of formation (Hc) from source phase
        virtual tmp<volScalarField> HcSource
        (
             const word& speciesName
        ) const = 0;

        //- Explicit mass transfer coefficient
        virtual tmp<volScalarField> Kexp
        (
            label modelVariable,
            const volScalarField& field
        ) const = 0;

        //- Implicit mass transfer coefficient of the shape
        //  Kimp*(variable - refValue)
        virtual tmp<volScalarField> Kimp
        (
            label modelVariable,
            const volScalarField& field
        ) const = 0;

        //- Explicit mass transfer for energy
        virtual tmp<volScalarField> KexpEnergy
        (
            label modelVariable,
            const volScalarField& field
        ) const;

        //- Reference value
        virtual const dimensionedScalar& Tactivate() const = 0;

        //- Returns the variable on which the model is based
        const word variable() const;

        //- Returns type of model. false is fully Explicit
        // true is treated semiImplicit as :
        // Sk-i = Kimp()(variable - Tactivate) for positive k to i
        // Si-k = Kimp()(Tactivate - variable) for positive i to k
        bool semiImplicit() const;

        //- Derivate sign of the source term:
        // d(Sk-i)/d(variable) > 0 or d(Si-k)/d(variable) < 0
        virtual label dSdVariable() = 0;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
