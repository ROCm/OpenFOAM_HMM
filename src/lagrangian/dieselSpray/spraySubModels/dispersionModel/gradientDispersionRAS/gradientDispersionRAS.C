/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "error.H"

#include "gradientDispersionRAS.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(gradientDispersionRAS, 0);

    addToRunTimeSelectionTable
    (
        dispersionModel,
        gradientDispersionRAS,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::gradientDispersionRAS::gradientDispersionRAS
(
    const dictionary& dict,
    spray& sm
)
:
    dispersionRASModel(dict, sm)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::gradientDispersionRAS::~gradientDispersionRAS()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::gradientDispersionRAS::disperseParcels() const
{

    const scalar cps = 0.16432;

    scalar dt = spray_.runTime().deltaTValue();
    const volScalarField& k = turbulence().k();
    const volVectorField gradk(fvc::grad(k));
    const volScalarField& epsilon = turbulence().epsilon();
    const volVectorField& U = spray_.U();

    forAllIter(spray, spray_, elmnt)
    {
        const label cellI = elmnt().cell();
        scalar UrelMag = mag(elmnt().U() - U[cellI] - elmnt().Uturb());

        scalar Tturb = min
        (
            k[cellI]/epsilon[cellI],
            cps*pow(k[cellI], 1.5)/epsilon[cellI]/(UrelMag + SMALL)
        );
        // parcel is perturbed by the turbulence
        if (dt < Tturb)
        {
            elmnt().tTurb() += dt;

            if (elmnt().tTurb() > Tturb)
            {
                elmnt().tTurb() = 0.0;

                scalar sigma = sqrt(2.0*k[cellI]/3.0);
                vector dir = -gradk[cellI]/(mag(gradk[cellI]) + SMALL);

                // numerical recipes... Ch. 7. Random Numbers...
                scalar x1 = 0.0;
                scalar x2 = 0.0;
                scalar rsq = 10.0;
                while ((rsq > 1.0) || (rsq == 0.0))
                {
                    x1 = 2.0*spray_.rndGen().sample01<scalar>() - 1.0;
                    x2 = 2.0*spray_.rndGen().sample01<scalar>() - 1.0;
                    rsq = x1*x1 + x2*x2;
                }

                scalar fac = sqrt(-2.0*log(rsq)/rsq);

                // in 2D calculations the -grad(k) is always
                // away from the axis of symmetry
                // This creates a 'hole' in the spray and to
                // prevent this we let x1 be both negative/positive
                if (spray_.twoD())
                {
                    fac *= x1;
                }
                else
                {
                    fac *= mag(x1);
                }

                elmnt().Uturb() = sigma*fac*dir;
            }
        }
        else
        {
            elmnt().tTurb() = GREAT;
            elmnt().Uturb() = vector::zero;
        }
    }
}


// ************************************************************************* //
