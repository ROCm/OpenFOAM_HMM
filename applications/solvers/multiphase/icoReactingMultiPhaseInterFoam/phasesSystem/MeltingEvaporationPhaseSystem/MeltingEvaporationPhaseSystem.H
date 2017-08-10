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
    Foam::MeltingEvaporationPhaseSystem

Description
    Holds mass transfer models

SourceFiles
    MeltingEvaporationPhaseSystem.C

\*---------------------------------------------------------------------------*/

#ifndef MeltingEvaporationPhaseSystem_H
#define MeltingEvaporationPhaseSystem_H

#include "phaseSystem.H"
#include "HashPtrTable.H"
#include "interfaceCompositionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class MeltingEvaporationPhaseSystem Declaration
\*---------------------------------------------------------------------------*/

template<class BasePhaseSystem>
class MeltingEvaporationPhaseSystem
:
    public BasePhaseSystem
{
protected:

    // Protected typedefs

        typedef
            HashTable
            <
                autoPtr<interfaceCompositionModel>,
                phasePairKey,
                phasePairKey::hash
            >
            meltingEvaporationModelTable;


       typedef
            HashPtrTable
            <
                volScalarField,
                Pair<word>,
                FixedList<word, 2>::Hash<>
            >
            dmdtYiTable;


       typedef
            HashTable
            <
                dmdtYiTable,
                phasePairKey,
                phasePairKey::hash
            >
            dmdtPhaseVarYiTable;


       typedef
            HashPtrTable
            <
                volScalarField,
                phasePairKey,
                phasePairKey::hash
            >
            dmdtTable;


    // Protected data

        //- Overall inter-phase mass transfer rates
        dmdtTable dmdt_;

        //- Mass transfer models
        meltingEvaporationModelTable meltingEvaporationModels_;

        //- Variable based inter-phase mass transfer rates
        dmdtPhaseVarYiTable dmdtYi_;


public:

    // Constructors

        //- Construct from fvMesh
        MeltingEvaporationPhaseSystem(const fvMesh&);


    //- Destructor
    virtual ~MeltingEvaporationPhaseSystem();


    // Member Functions

        // Total mass transfer per variable

            //- Return total interfacial mass flow rate
            virtual tmp<volScalarField> dmdt(const phasePairKey& key) const;


            //- Return species interfacial source mass rate per species
            //  driven by pressure or temperature
            virtual tmp<volScalarField> dmdtYi
            (
                const word& phaseSpeciesName
            ) const;


            //- Return total interfacial mass driven by modelVariable
            tmp<volScalarField> dmdtYi
            (
                const phasePairKey& key,
                word modelVar,
                word specieName
            ) const;


            //- Return non-const reference total interfacial mass driven
            //  by modelVariable
            volScalarField& dmdtYi
            (
                const phasePairKey& key,
                word modelVar,
                word specieName
            );


        // Explicit and Implicit mass transfer matrices

            //- Return the heat transfer matrix
            virtual  tmp<fvScalarMatrix> heatTransfer
            (
                const volScalarField& T
            );

            //- Return semi-implicit mass transfer matrix for
            //  multicomponent phases based on mass fractions
            virtual  autoPtr<phaseSystem::massTransferTable> massTransfer
            (
                const volScalarField& T
            );


        //- Return phase inverse density
        virtual tmp<volScalarField> coeffs(const word& key) const;

        //- Correct the phase properties other than the thermo and turbulence
        virtual void correct();

        //- Read base phaseProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "MeltingEvaporationPhaseSystem.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
