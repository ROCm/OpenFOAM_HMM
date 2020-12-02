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
    Foam::incompressible::sensitivityBezier

Description
    Calculation of adjoint based sensitivities for Bezier control points

SourceFiles
    sensitivityBezier.C

\*---------------------------------------------------------------------------*/

#ifndef sensitivityBezierIncompressible_H
#define sensitivityBezierIncompressible_H

#include "primitiveFieldsFwd.H"
#include "volFieldsFwd.H"
#include "pointFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "volPointInterpolation.H"
#include "SIBaseIncompressible.H"
#include "PrimitivePatchInterpolation.H"
#include "PrimitivePatch.H"
#include "deltaBoundary.H"
#include "Bezier.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace incompressible
{

/*---------------------------------------------------------------------------*\
                      Class sensitivityBezier Declaration
\*---------------------------------------------------------------------------*/

class sensitivityBezier
:
    public SIBase
{
protected:

    // Protected data

        Bezier Bezier_;

        vectorField sens_;
        vectorField flowSens_;
        vectorField dSdbSens_;
        vectorField dndbSens_;
        vectorField dxdbDirectSens_;
        vectorField bcSens_;

        fileName derivativesFolder_;


private:

    // Private Member Functions

        //- No copy construct
        sensitivityBezier(const sensitivityBezier&) = delete;

        //- No copy assignment
        void operator=(const sensitivityBezier&) = delete;


public:

    //- Runtime type information
    TypeName("Bezier");


    // Constructors

        //- Construct from components
        sensitivityBezier
        (
            const fvMesh& mesh,
            const dictionary& dict,
            incompressibleVars& primalVars,
            incompressibleAdjointVars& adjointVars,
            objectiveManager& objectiveManager
        );


    //- Destructor
    virtual ~sensitivityBezier() = default;


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
