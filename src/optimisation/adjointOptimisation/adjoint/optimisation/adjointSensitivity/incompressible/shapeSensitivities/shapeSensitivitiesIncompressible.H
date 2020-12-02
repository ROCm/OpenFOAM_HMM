/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 PCOpt/NTUA
    Copyright (C) 2020 FOSS GP
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
    Foam::incompressible::shapeSensitivitiesBase

Description
    Base class supporting shape sensitivity derivatives for
    incompressible flows

SourceFiles
    shapeSensitivitiesBase.C

\*---------------------------------------------------------------------------*/

#ifndef shapeSensitivitiesIncompressible_H
#define shapeSensitivitiesIncompressible_H

#include "adjointSensitivityIncompressible.H"
#include "shapeSensitivitiesBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

/*---------------------------------------------------------------------------*\
                      Class shapeSensitivities Declaration
\*---------------------------------------------------------------------------*/

class shapeSensitivities
:
    public adjointSensitivity,
    public shapeSensitivitiesBase
{
protected:

    // Protected data

        //- Fields related to direct sensitivities
        autoPtr<boundaryVectorField> dSfdbMult_;
        autoPtr<boundaryVectorField> dnfdbMult_;
        autoPtr<boundaryVectorField> dxdbDirectMult_;
        autoPtr<boundaryVectorField> bcDxDbMult_;


    // Protected Member Fuctions

        //- Accumulate direct sensitivities
        virtual void accumulateDirectSensitivityIntegrand(const scalar dt);

        //- Accumulate sensitivities enamating from the boundary conditions
        virtual void accumulateBCSensitivityIntegrand(const scalar dt);

        //- Compute multiplier of dv_i/db
        tmp<boundaryVectorField> dvdbMult() const;


private:

    // Private Member Functions

        //- No copy construct
        shapeSensitivities(const shapeSensitivities&) = delete;

        //- No copy assignment
        void operator=(const shapeSensitivities&) = delete;


public:

    //- Runtime type information
    TypeName("shapeSensitivities");


    // Constructors

        //- Construct from components
        shapeSensitivities
        (
            const fvMesh& mesh,
            const dictionary& dict,
            incompressibleVars& primalVars,
            incompressibleAdjointVars& adjointVars,
            objectiveManager& objectiveManager
        );


    //- Destructor
    virtual ~shapeSensitivities() = default;


    // Member Functions

        //- Accumulate sensitivity integrands
        virtual void accumulateIntegrand(const scalar dt) = 0;

        //- Assemble sensitivities
        virtual void assembleSensitivities() = 0;

        //- Zero sensitivity fields and their constituents
        virtual void clearSensitivities();

        //- Write sensitivity fields.
        virtual void write(const word& baseName = word::null);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
