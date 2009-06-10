/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::chemistryModel

Description
    Foam::chemistryModel

SourceFiles
    chemistryModelI.H
    chemistryModel.C
    newChemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef chemistryModel_H
#define chemistryModel_H

#include "IOdictionary.H"
#include "Switch.H"
#include "scalarField.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "hCombustionThermo.H"
//#include "basicMultiComponentMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                       class chemistryModel Declaration
\*---------------------------------------------------------------------------*/

class chemistryModel
:
    public IOdictionary
{
    // Private Member Functions

        //- Construct as copy (not implemented)
        chemistryModel(const chemistryModel&);

        //- Disallow default bitwise assignment
        void operator=(const chemistryModel&);


protected:

    // Protected data

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Thermo package
        autoPtr<hCombustionThermo> thermo_;

        //- Chemistry activation switch
        Switch chemistry_;

        //- Latest estimation of integration step
        scalarField deltaTChem_;


    // Protected member functions

        //- Return non-const access to the latest estimation of integration
        //  step, e.g. for multi-chemistry model
        scalarField& deltaTChem();


public:

    //- Runtime type information
    TypeName("chemistryModel");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        chemistryModel,
        fvMesh,
        (
            const fvMesh& mesh
        ),
        (mesh)
    );


    // Constructors

        //- Construct from mesh
        chemistryModel(const fvMesh& mesh);


    //- Selector
    static autoPtr<chemistryModel> New(const fvMesh& mesh);

    //- Destructor
    virtual ~chemistryModel();


    // Member Functions

        //- Return const access to the mesh database
        inline const fvMesh& mesh() const;

        //- Return access to the thermo package
        inline hCombustionThermo& thermo();

        //- Return const access to the thermo package
        inline const hCombustionThermo& thermo() const;

        //- Chemistry activation switch
        inline Switch chemistry() const;

        //- Return the latest estimation of integration step
        inline const scalarField& deltaTChem() const;


        // Functions to be derived in derived classes

            // Fields

                //- Return const access to chemical source terms
                virtual tmp<volScalarField> RR(const label i) const = 0;


            // Chemistry solution

                //- Solve the reaction system for the given start time and
                //  timestep and return the characteristic time
                virtual scalar solve(const scalar t0, const scalar deltaT) = 0;

                //- Return the chemical time scale
                virtual tmp<volScalarField> tc() const = 0;

                //- Return the heat release
                virtual tmp<volScalarField> dQ() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "chemistryModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
