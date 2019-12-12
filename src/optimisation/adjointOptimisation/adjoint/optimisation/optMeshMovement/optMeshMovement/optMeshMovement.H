/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019 PCOpt/NTUA
    Copyright (C) 2013-2019 FOSS GP
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::optMeshMovement

Description
    Abstract base class for translating an update of the design variables
    into mesh movement

SourceFiles
    optMeshMovement.C

\*---------------------------------------------------------------------------*/

#ifndef optMeshMovement_H
#define optMeshMovement_H

#include "displacementMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class optMeshMovement Declaration
\*---------------------------------------------------------------------------*/

class optMeshMovement
{
private:

    // Private Member Functions

        //- Max allowed boundary displacement for the first optimisation cycle
        autoPtr<scalar> maxAllowedDisplacement_;


    // Private Member Functions

        //- No copy construct
        optMeshMovement(const optMeshMovement&) = delete;

        //- No copy assignment
        void operator=(const optMeshMovement&) = delete;


protected:

    // Protected data

        fvMesh& mesh_;
        const dictionary& dict_;

        //- Correction of design variables
        scalarField correction_;

        //- IDs of patches to be moved
        labelList patchIDs_;

        //- Fall back points in case line-search is used
        vectorField pointsInit_;

        //- Mesh movement engine and interface for applying mesh movement
        //- boundary conditions
        autoPtr<displacementMethod> displMethodPtr_;

        //- Whether to write the mesh quality metrics to files each time the
        //- mesh is updated
        bool writeMeshQualityMetrics_;


    // Protected Member Functions

        //- Get maxAllowedDisplacement, is set
        scalar getMaxAllowedDisplacement() const;


public:

    //- Runtime type information
    TypeName("optMeshMovement");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            optMeshMovement,
            dictionary,
            (
                fvMesh& mesh,
                const dictionary& dict,
                const labelList& patchIDs
            ),
            (
                mesh,
                dict,
                patchIDs
            )
        );


    // Constructors

        //- Construct from components
        optMeshMovement
        (
            fvMesh& mesh,
            const dictionary& dict,
            const labelList& patchIDs
        );


    // Selectors

        static autoPtr<optMeshMovement> New
        (
            fvMesh& mesh,
            const dictionary& dict,
            const labelList& patchIDs
        );


    //- Destructor
    virtual ~optMeshMovement() = default;


    // Member Functions

       //- Set design variable correction
       void setCorrection(const scalarField& correction);

       //- Calculates mesh movemnt based on the correction of the design
       //- variables
       virtual void moveMesh();

       //- Return displacementMethod
       autoPtr<displacementMethod>& returnDisplacementMethod();

       //- Return patchIDs
       const labelList& getPatchIDs();

       //- Write mesh quality metrics
       void writeMeshQualityMetrics();

       //- Store design variables and mesh, to act as the starting point of
       //- line search
       virtual void storeDesignVariables();

       //- Reset to starting point of line search
       virtual void resetDesignVariables();

       //- Compute eta value based on max displacement
       virtual scalar computeEta(const scalarField& correction) = 0;

       //- Whether maxAllowedDisplacement has been set
       bool maxAllowedDisplacementSet() const;

       //- Return active design variables.
       //  Implemented only for certain parametetisations
       virtual labelList getActiveDesignVariables() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
