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

#include "variableHeatTransfer.H"
#include "IObasicSourceList.H"
#include "turbulenceModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(variableHeatTransfer, 0);
    addToRunTimeSelectionTable
    (
        basicSource,
        variableHeatTransfer,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::variableHeatTransfer::variableHeatTransfer
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    interRegionHeatTransferModel(name, modelType, dict, mesh),
    a_(0),
    b_(0),
    c_(0),
    ds_(0),
    Pr_(0),
    area_()
{
    if (master_)
    {
        a_ = readScalar(dict_.lookup("a"));
        b_ = readScalar(dict_.lookup("b"));
        c_ = readScalar(dict_.lookup("c"));
        ds_ = readScalar(dict_.lookup("ds"));
        Pr_ = readScalar(dict_.lookup("Pr"));
        area_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "area",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::variableHeatTransfer::~variableHeatTransfer()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


const Foam::tmp<Foam::volScalarField> Foam::variableHeatTransfer::
calculateHtc()
{

    const fvMesh& secondaryMesh =
        mesh_.time().lookupObject<fvMesh>(mapRegionName());

    const compressible::turbulenceModel& turb =
        secondaryMesh.lookupObject<compressible::turbulenceModel>
        (
            "turbulenceModel"
        );

    const basicThermo& secondaryThermo =
        secondaryMesh.lookupObject<basicThermo>
        (
            "thermophysicalProperties"
        );

    const volVectorField& U =
        secondaryMesh.lookupObject<volVectorField>("U");

    const volScalarField Re
    (
        mag(U)*ds_*secondaryThermo.rho()/turb.mut()
    );

    const volScalarField Nu(a_*pow(Re, b_)*pow(Pr_, c_));

    const volScalarField K(turb.alphaEff()*secondaryThermo.Cp());

    scalarField htcMapped(htc_.internalField().size(), 0.0);

    secondaryToPrimaryInterpPtr_->interpolateInternalField
    (
        htcMapped,
        Nu*K/ds_,
        meshToMesh::MAP,
        eqOp<scalar>()
    );

    htc_.internalField() = htcMapped*area_/mesh_.V();

    return htc_;
}


void Foam::variableHeatTransfer::writeData(Ostream& os) const
{
    os  << indent << token::BEGIN_BLOCK << incrIndent << nl;
    interRegionHeatTransferModel::writeData(os);

    os.writeKeyword("a") << a_ << token::END_STATEMENT << nl;
    os.writeKeyword("b") << b_ << token::END_STATEMENT << nl;
    os.writeKeyword("c") << c_ << token::END_STATEMENT << nl;
    os.writeKeyword("ds") << ds_ << token::END_STATEMENT << nl;
    os.writeKeyword("Pr") << Pr_ << token::END_STATEMENT << nl;

    os << indent << "variableHeatTransfer";

    dict_.write(os);

    os << decrIndent << indent << token::END_BLOCK << endl;
}


bool Foam::variableHeatTransfer::read(const dictionary& dict)
{
    if (basicSource::read(dict))
    {

        const dictionary& sourceDict = dict.subDict(name());
        const dictionary& subDictCoeffs =
            sourceDict.subDict(typeName + "Coeffs");
        subDictCoeffs.readIfPresent("a", a_);
        subDictCoeffs.readIfPresent("b", b_);
        subDictCoeffs.readIfPresent("c", c_);
        subDictCoeffs.readIfPresent("ds", ds_);
        subDictCoeffs.readIfPresent("Pr", Pr_);

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //