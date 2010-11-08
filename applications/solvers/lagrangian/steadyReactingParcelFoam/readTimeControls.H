/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2010 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

// Maximum flow Courant number
scalar maxCo(readScalar(runTime.controlDict().lookup("maxCo")));

// Maximum time scale
scalar maxDeltaT = readScalar(runTime.controlDict().lookup("maxDeltaT"));

// Smoothing parameter (0-1) when smoothing iterations > 0
scalar alphaTauSmooth
(
    runTime.controlDict().lookupOrDefault("alphaTauSmooth", 0.1)
);

// Maximum change in cell density per iteration (relative to previous value)
scalar alphaTauRho
(
    runTime.controlDict().lookupOrDefault("alphaTauRho", 0.05)
);

// Maximum change in cell velocity per iteration (relative to previous value)
scalar alphaTauU
(
    runTime.controlDict().lookupOrDefault("alphaTauU", 0.05)
);

// Maximum change in cell temperature per iteration (relative to previous value)
scalar alphaTauTemp
(
    runTime.controlDict().lookupOrDefault("alphaTauTemp", 0.05)
);

// Max specie mass fraction that can be consumed/gained per chemistry
// integration step
scalar alphaTauSpecie
(
    runTime.controlDict().lookupOrDefault("alphaTauSpecie", 0.05)
);

// Maximum unboundedness allowed (fraction of 1)
scalar specieMaxUnbound
(
    runTime.controlDict().lookupOrDefault("specieMaxUnbound", 0.01)
);


// ************************************************************************* //
