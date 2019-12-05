/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

Description
    Test minMax

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "Time.H"
#include "BitOps.H"
#include "HashOps.H"
#include "ListOps.H"
#include "scalarField.H"
#include "MinMax.H"
#include "dimensionedScalar.H"
#include "dimensionedMinMax.H"

using namespace Foam;


template<class T>
Ostream& printInfo(const MinMax<T>& range)
{
    Info<< range << " valid=" << range.valid();

    return Info;
}


dimensionedScalarMinMax rhoLimit(const dictionary& dict)
{
    Info<< "From " << dict;

    dimensionedScalarMinMax range =
        makeDimensionedMinMax<scalar>
        (
            "rhoLimit", dimDensity, scalarMinMax{Zero, GREAT}, dict,
            "rhoMin", "rhoMax"
        );

    Info<< "=> " << range << nl;

    return range;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noCheckProcessorDirectories();

    #include "setRootCase.H"

    Info<< "Test min/max " << nl;


    {
        scalarMinMax range1(10, 20);
        scalarMinMax range2(40, 50);
        Info<< range1 << " + " << range2 << " = " << (range1 + range2) <<nl;
    }

    {
        Info<< "Dimensioned range : "
            << dimensioned<scalarMinMax>("velrange", dimVelocity, {1, 20})
            << nl;

        dimensioned<scalarMinMax> range1("a", dimVelocity, {10, 20});
        dimensioned<scalarMinMax> range2("b", dimVelocity, {40, 50});

        Info<< "Dimensioned range : " << (range1 + range2) << endl;
    }

    {
        Info<< nl << "makeDimensionedMinMax:" << nl << nl;

        Info
            << makeDimensionedMinMax<scalar>("rhoa", dimDensity, 1, 20)
            << nl;

        {
            dimensionedScalar minval("min", dimDensity, 0.3);
            dimensionedScalar maxval("max", dimDensity, 0.5);

            Info
                << makeDimensionedMinMax<scalar>(minval, maxval)
                << nl;

            Info
                << makeDimensionedMinMax<scalar>("rhob", minval, maxval)
                << nl;

        }


        {
            dictionary dict1, dict2, dict3, dict4;

            dict1.add("rhoMin", dimensionedScalar("", dimDensity, 0.1));
            dict2.add("rhoMax", dimensionedScalar("", dimDensity, 20));

            dict3.add("rhoMin", dimensionedScalar("", dimDensity, 0.3));
            dict3.add("rhoMax", dimensionedScalar("", dimDensity, 30));

            dict4.add
            (
                "rhoLimit",
                dimensionedScalarMinMax("", dimDensity, scalarMinMax(0.4, 40))
            );

            rhoLimit(dict1);
            rhoLimit(dict2);
            rhoLimit(dict3);
            rhoLimit(dict4);
        }
    }


    return 0;
}


// ************************************************************************* //
