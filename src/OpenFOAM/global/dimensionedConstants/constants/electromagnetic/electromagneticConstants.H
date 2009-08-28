/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2009 OpenCFD Ltd.
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

Namespace
    Foam::constant::em

Description
    Electromagnetic constants

\*---------------------------------------------------------------------------*/

#ifndef electromagneticConstants_H
#define electromagneticConstants_H

#include "dimensionedScalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace constant
{
namespace electromagnetic
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    //- Group name for electromagnetic constants
    extern const char* group;

    //- Magnetic constant/permeability of free space: default SI units: [H/m]
    extern const dimensionedScalar mu0;

    //- Electric constant: default SI units: [F/m]
    extern const dimensionedScalar epsilon0;

    //- Characteristic impedance of a vacuum: default SI units: [ohm]
    extern const dimensionedScalar Z0;

    //- Coulomb constant: default SI units: [N.m2/C2]
    extern const dimensionedScalar kappa;

    //- Conductance quantum: default SI units: [S]
    extern const dimensionedScalar G0;

    //- Josephson constant: default SI units: [Hz/V]
    extern const dimensionedScalar KJ;

    //- Magnetic flux quantum: default SI units: [Wb]
    extern const dimensionedScalar phi0;

    //- von Klitzing constant: default SI units: [ohm]
    extern const dimensionedScalar RK;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace electromagnetic
} // end namespace constant
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //


