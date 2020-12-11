/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2020 PCOpt/NTUA
    Copyright (C) 2013-2020 FOSS GP
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
    Foam::incompressible::sensitivityBezierFI

Description
    Calculation of adjoint based sensitivities for Bezier control points
    using the FI appoach

SourceFiles
    sensitivityBezierFI.C

\*---------------------------------------------------------------------------*/

#ifndef sensitivityBezierFIIncompressible_H
#define sensitivityBezierFIIncompressible_H

#include "primitiveFieldsFwd.H"
#include "volFieldsFwd.H"
#include "pointFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "volPointInterpolation.H"
#include "FIBaseIncompressible.H"
#include "PrimitivePatchInterpolation.H"
#include "deltaBoundary.H"
#include "Bezier.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace incompressible
{

/*---------------------------------------------------------------------------*\
                     Class sensitivityBezierFI Declaration
\*---------------------------------------------------------------------------*/

class sensitivityBezierFI
:
    public FIBase
{
protected:

    // Protected data
        Bezier Bezier_;

        //- Flow related term
        scalarField flowSens_;

        //- Term depending on delta(n dS)/delta b
        scalarField dSdbSens_;

        //- Term depending on delta(n)/delta b
        scalarField dndbSens_;

        //- Term depending on delta(x)/delta b for objectives that directly
        //- depend on x
        scalarField dxdbDirectSens_;

        //- Term depending on delta(V)/delta b
        scalarField dVdbSens_;

        //- Term depending on distance differentiation
        scalarField distanceSens_;

        //- Term depending on fvOptions
        scalarField optionsSens_;

        //- Term depending on the differenation of boundary conditions
        scalarField bcSens_;

        fileName derivativesFolder_;

        label meshMovementIters_;
        scalar meshMovementResidualLimit_;
        volVectorField dxdb_;

        void read();

        tmp<volVectorField> solveMeshMovementEqn
        (
            const label iCP,
            const label idir
        );


private:

    // Private Member Functions

        //- No copy construct
        sensitivityBezierFI(const sensitivityBezierFI&) = delete;

        //- No copy assignment
        void operator=(const sensitivityBezierFI&) = delete;


public:

    //- Runtime type information
    TypeName("BezierFI");


    // Constructors

        //- Construct from components
        sensitivityBezierFI
        (
            const fvMesh& mesh,
            const dictionary& dict,
            incompressibleVars& primalVars,
            incompressibleAdjointVars& adjointVars,
            objectiveManager& objectiveManager
        );


    //- Destructor
    virtual ~sensitivityBezierFI() = default;


    // Member Functions

        //- Assemble sensitivities
        virtual void assembleSensitivities();

        //- Zero sensitivity fields and their constituents
        virtual void clearSensitivities();

        //- Write sensitivities to file
        virtual void write(const word& baseName = word::null);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
