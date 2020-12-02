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
    Foam::incompressible::SIBase

Description
    Base class for Surface Integral-based sensitivity derivatives

SourceFiles
    SIBase.C

\*---------------------------------------------------------------------------*/

#ifndef SIBaseIncompressible_H
#define SIBaseIncompressible_H

#include "shapeSensitivitiesIncompressible.H"
#include "sensitivitySurfaceIncompressible.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace incompressible
{

/*---------------------------------------------------------------------------*\
                            Class SIBase Declaration
\*---------------------------------------------------------------------------*/

class SIBase
:
    public shapeSensitivities
{
protected:

    // Protected data

        //- Surface sensitivities
        sensitivitySurface surfaceSensitivity_;

        //- Whether to include direct sensitivities or not.
        //  Used to avoid double contributions from both here and the
        //  sensitivitySurface object which might have already accounted for
        //  them
        bool includeObjective_;

        //- Write sensitivity map upon write
        bool writeSensitivityMap_;


    // Protected Member Functions

        //- Read options from dict
        void read();


private:

    // Private Member Functions

        //- No copy construct
        SIBase(const SIBase&) = delete;

        //- No copy assignment
        void operator=(const SIBase&) = delete;


public:

    //- Runtime type information
    TypeName("volumetricBSplinesFI");


    // Constructors

        //- Construct from components
        SIBase
        (
            const fvMesh& mesh,
            const dictionary& dict,
            incompressibleVars& primalVars,
            incompressibleAdjointVars& adjointVars,
            objectiveManager& objectiveManager
        );


    //- Destructor
    virtual ~SIBase() = default;


    // Member Functions

        //- Read dict if changed
        virtual bool readDict(const dictionary& dict);

        //- Accumulate sensitivity integrands
        virtual void accumulateIntegrand(const scalar dt);

        //- Assemble sensitivities
        virtual void assembleSensitivities() = 0;

        //- Zero sensitivity fields and their constituents
        virtual void clearSensitivities();

        //- Return reference to underlaying surface sensitivities
        const sensitivitySurface& getSurfaceSensitivities() const;

        //- Write sensitivity map
        virtual void write(const word& baseName = word::null);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
