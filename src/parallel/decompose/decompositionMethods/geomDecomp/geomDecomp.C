/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "geomDecomp.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::geomDecomp::readCoeffs()
{
    coeffsDict_.readIfPresent("delta", delta_);

    coeffsDict_.lookup("n") >> n_;

    // Verify that the input makes sense
    if (nDomains_ != n_.x()*n_.y()*n_.z())
    {
        FatalErrorInFunction
            << "Wrong number of domain divisions in geomDecomp:" << nl
            << "Number of domains    : " << nDomains_ << nl
            << "Wanted decomposition : " << n_
            << exit(FatalError);
    }

    const scalar d = 1 - 0.5*delta_*delta_;
    const scalar d2 = sqr(d);

    const scalar a = delta_;
    const scalar a2 = sqr(a);

    rotDelta_ = tensor
    (
        d2,         -a*d,         a,
        a*d - a2*d,  a*a2 + d2,  -2*a*d,
        a*d2 + a2,   a*d - a2*d,  d2 - a2
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::geomDecomp::geomDecomp
(
    const word& derivedType,
    const dictionary& decompDict,
    int select
)
:
    decompositionMethod(decompDict),
    coeffsDict_(findCoeffsDict(derivedType + "Coeffs", select)),
    n_(1,1,1),
    delta_(0.001),
    rotDelta_(I)
{
    readCoeffs();
}


Foam::geomDecomp::geomDecomp
(
    const word& derivedType,
    const dictionary& decompDict,
    const word& regionName,
    int select
)
:
    decompositionMethod(decompDict, regionName),
    coeffsDict_(findCoeffsDict(derivedType + "Coeffs", select)),
    n_(1,1,1),
    delta_(0.001),
    rotDelta_(I)
{
    readCoeffs();
}


// ************************************************************************* //
