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
    Foam::incompressible::sensitivityVolBSplinesFI

Description
    Calculation of adjoint based sensitivities at vol B-Splines control points
    using the FI approach.

SourceFiles
    sensitivityVolBSplinesFI.C

\*---------------------------------------------------------------------------*/

#ifndef sensitivityVolBSplinesFIIncompressible_H
#define sensitivityVolBSplinesFIIncompressible_H

#include "FIBaseIncompressible.H"
#include "volBSplinesBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace incompressible
{

/*---------------------------------------------------------------------------*\
                  Class sensitivityVolBSplinesFI Declaration
\*---------------------------------------------------------------------------*/

class sensitivityVolBSplinesFI
:
    public FIBase
{
protected:

    // Protected data

        //- Reference to underlaying volumetric B-Splines morpher
        volBSplinesBase& volBSplinesBase_;

        //- Flow related term
        vectorField flowSens_;

        //- Term depending on delta(n dS)/delta b
        vectorField dSdbSens_;

        //- Term depending on delta(n)/delta b
        vectorField dndbSens_;

        //- Term depending on delta(x)/delta b for objectives that directly
        //- depend on x
        vectorField dxdbDirectSens_;

        //- Term depending on delta(V)/delta b
        vectorField dVdbSens_;

        //- Term depending on distance differentiation
        vectorField distanceSens_;

        //- Term depending on fvOptions
        vectorField optionsSens_;

        //- Term depending on the differentiation of boundary conditions
        vectorField bcSens_;

        fileName derivativesFolder_;


private:

    // Private Member Functions

        //- No copy construct
        sensitivityVolBSplinesFI(const sensitivityVolBSplinesFI&) = delete;

        //- No copy assignment
        void operator=(const sensitivityVolBSplinesFI&) = delete;


public:

    //- Runtime type information
    TypeName("volumetricBSplinesFI");


    // Constructors

        //- Construct from components
        sensitivityVolBSplinesFI
        (
            const fvMesh& mesh,
            const dictionary& dict,
            incompressibleVars& primalVars,
            incompressibleAdjointVars& adjointVars,
            objectiveManager& objectiveManager
        );


    //- Destructor
    virtual ~sensitivityVolBSplinesFI() = default;


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
