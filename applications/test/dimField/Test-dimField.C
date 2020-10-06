/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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

Application
    Test-dimField

Description
    Simple tests for DimensionedField

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "GeometricFields.H"

// #undef TEST_UINT8_FIELD

#ifdef TEST_UINT8_FIELD
namespace Foam
{
    // Something like an internal state field. Probably only dimensionless
    typedef DimensionedField<uint8_t, volMesh> dimUint8Field;

    defineTemplate2TypeNameAndDebug(dimUint8Field, 0);

} // End namespace Foam
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    {
        Info<< "Tensor field\n" << endl;
        DimensionedField<tensor, volMesh> tensorfld
        (
            IOobject
            (
                "tensor",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<tensor>(dimless, tensor(1,2,3,4,5,6,7,8,9))
        );

        Info().beginBlock("transformed")
            << tensorfld.T() << nl;
        Info().endBlock();
    }

    #ifdef TEST_UINT8_FIELD
    {
        Info<< "uint8 field\n" << endl;
        DimensionedField<uint8_t, volMesh> statefld
        (
            IOobject
            (
                "state",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<uint8_t>(dimless, uint8_t{100})
        );

        Info().beginBlock("stateField")
            << statefld << nl;
        Info().endBlock();
    }
    #endif


    Info<< "\nEnd\n" << nl;

    return 0;
}


// ************************************************************************* //
