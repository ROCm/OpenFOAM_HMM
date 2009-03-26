/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "blackBodyEmission.H"
#include "dimensionedConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::scalar
Foam::radiation::blackBodyEmission::
emissivePowerTable[Foam::radiation::blackBodyEmission::nTableEntries][2] =
{
    {1000.0, 0.00},
    {1100.0, 0.00},
    {1200.0, 0.00},
    {1300.0, 0.00},
    {1400.0, 0.01},
    {1500.0, 0.01},
    {1600.0, 0.02},
    {1700.0, 0.03},
    {1800.0, 0.04},
    {1900.0, 0.05},
    {2000.0, 0.07},
    {2100.0, 0.08},
    {2200.0, 0.10},
    {2300.0, 0.12},
    {2400.0, 0.14},
    {2500.0, 0.16},
    {2600.0, 0.18},
    {2700.0, 0.21},
    {2800.0, 0.23},
    {2900.0, 0.25},
    {3000.0, 0.27},
    {3100.0, 0.30},
    {3200.0, 0.32},
    {3300.0, 0.34},
    {3400.0, 0.36},
    {3500.0, 0.38},
    {3600.0, 0.40},
    {3700.0, 0.42},
    {3800.0, 0.44},
    {3900.0, 0.46},
    {4000.0, 0.48},
    {4100.0, 0.50},
    {4200.0, 0.52},
    {4300.0, 0.53},
    {4400.0, 0.55},
    {4500.0, 0.56},
    {4600.0, 0.58},
    {4700.0, 0.59},
    {4800.0, 0.61},
    {4900.0, 0.62},
    {5000.0, 0.63},
    {5100.0, 0.65},
    {5200.0, 0.66},
    {5300.0, 0.67},
    {5400.0, 0.68},
    {5500.0, 0.69},
    {5600.0, 0.70},
    {5700.0, 0.71},
    {5800.0, 0.72},
    {5900.0, 0.73},
    {6000.0, 0.74},
    {6100.0, 0.75},
    {6200.0, 0.75},
    {6300.0, 0.76},
    {6400.0, 0.77},
    {6500.0, 0.78},
    {6600.0, 0.78},
    {6700.0, 0.79},
    {6800.0, 0.80},
    {6900.0, 0.80},
    {7000.0, 0.81},
    {7100.0, 0.81},
    {7200.0, 0.82},
    {7300.0, 0.82},
    {7400.0, 0.83},
    {7500.0, 0.83},
    {7600.0, 0.84},
    {7700.0, 0.84},
    {7800.0, 0.85},
    {7900.0, 0.85},
    {8000.0, 0.86},
    {8100.0, 0.86},
    {8200.0, 0.86},
    {8300.0, 0.87},
    {8400.0, 0.87},
    {8500.0, 0.87},
    {8600.0, 0.88},
    {8700.0, 0.88},
    {8800.0, 0.88},
    {8900.0, 0.89},
    {9000.0, 0.89},
    {9100.0, 0.89},
    {9200.0, 0.90},
    {9300.0, 0.90},
    {9400.0, 0.90},
    {9500.0, 0.90},
    {9600.0, 0.91},
    {9700.0, 0.91},
    {9800.0, 0.91},
    {9900.0, 0.91},
    {10000.0, 0.92}
};


// * * * * * * * * * * * * * Private member funstions  * * * * * * * * * * * //

Foam::List<Foam::Tuple2<Foam::scalar, Foam::scalar> >
Foam::radiation::blackBodyEmission::tableToList() const
{
    List<Tuple2<scalar, scalar> > newList(nTableEntries);

    forAll(newList, i)
    {
        newList[i].first() = emissivePowerTable[i][0];
        newList[i].second() = emissivePowerTable[i][1];
    }

    return newList;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::blackBodyEmission::blackBodyEmission
(
    const label nLambda,
    const volScalarField& T
)
:
    table_
    (
        tableToList(),
        interpolationTable<scalar>::CLAMP,
        "blackBodyEmissivePower"
    ),
    C1_("C1", dimensionSet(1, 4, 3, 0, 0, 0, 0), 3.7419e-16),
    C2_("C2", dimensionSet(0, 1, 0, 1, 0, 0, 0), 14.388e-6),
    bLambda_(nLambda),
    T_(T)
{
    forAll(bLambda_, lambdaI)
    {
        bLambda_.set
        (
            lambdaI,
            new volScalarField
            (
                IOobject
                (
                    "bLambda_" + Foam::name(lambdaI) ,
                    T.mesh().time().timeName(),
                    T.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                radiation::sigmaSB*pow4(T)
            )
        );

    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::blackBodyEmission::~blackBodyEmission()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::radiation::blackBodyEmission::fLambdaT
(
    const scalar lambdaT
) const
{
    return  table_(lambdaT*1.0e6);
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::blackBodyEmission::EbDeltaLambdaT
(
    const volScalarField& T,
    const Vector2D<scalar>& band
) const
{
    tmp<volScalarField> Eb
    (
        new volScalarField
        (
            IOobject
            (
                "Eb",
                T.mesh().time().timeName(),
                T.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            radiation::sigmaSB*pow4(T)
        )
    );


    if (band == Vector2D<scalar>::one)
    {
        return Eb;
    }
    else
    {
        forAll(T, i)
        {
            scalar T1 = fLambdaT(band[1]*T[i]);
            scalar T2 = fLambdaT(band[0]*T[i]);
            dimensionedScalar fLambdaDelta
            (
                "fLambdaDelta",
                dimless,
                T1 - T2
            );
            Eb()[i] = Eb()[i]*fLambdaDelta.value();
        }
        return Eb;
    }
}


void Foam::radiation::blackBodyEmission::correct
(
    const label lambdaI,
    const Vector2D<scalar>& band
)
{
    bLambda_[lambdaI] = EbDeltaLambdaT(T_, band);
}


// ************************************************************************* //
