/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

\*----------------------------------------------------------------------------*/

#include "porousZone.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "oneField.H"

// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //

// adjust negative resistance values to be multiplier of max value
void Foam::porousZone::adjustNegativeResistance(vector& resist)
{
    scalar maxCmpt = max(0, cmptMax(resist));

    for (label cmpt=0; cmpt < vector::nComponents; ++cmpt)
    {
        if (resist[cmpt] < 0)
        {
            resist[cmpt] *= -maxCmpt;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porousZone::porousZone
(
    const fvMesh& mesh,
    const word& name,
    const dictionary& dict
)
:
    mesh_(mesh),
    name_(name),
    dict_(dict),
    cellZoneID_(mesh_.cellZones().findZoneID(name)),
    coordSys_(dict),
    porosity_(1),
    C0_(0),
    C1_(0),
    D_("D", dimensionSet(0, -2, 0, 0, 0), tensor::zero),
    F_("F", dimensionSet(0, -1, 0, 0, 0), tensor::zero)
{
    if (cellZoneID_ == -1 && !Pstream::parRun())
    {
        FatalErrorIn
        (
            "Foam::porousZone::porousZone"
            "(const fvMesh&, const Istream&)"
        )   << "cannot find porous cellZone " << name_
            << exit(FatalError);
    }

    // local-to-global transformation tensor
    const tensor& E = coordSys_.R();

    // porosity
    if (dict_.found("porosity"))
    {
        dict_.lookup("porosity") >> porosity_;

        if (porosity_ <= 0.0 || porosity_ > 1.0)
        {
            FatalIOErrorIn
            (
                "Foam::porousZone::porousZone(const fvMesh&, const Istream&)",
                dict_
            )
                << "out-of-range porosity value " << porosity_
                << exit(FatalIOError);
        }
    }

    // powerLaw coefficients
    if (dict_.found("powerLaw"))
    {
        const dictionary& subDict = dict_.subDict("powerLaw");
        if (subDict.found("C0"))
        {
            subDict.lookup("C0") >> C0_;
        }
        if (subDict.found("C1"))
        {
            subDict.lookup("C1") >> C1_;
        }
    }

    // Darcy-Forchheimer coefficients
    if (dict_.found("Darcy"))
    {
        const dictionary& subDict = dict_.subDict("Darcy");

        dimensionedVector d("d", D_.dimensions(), vector::zero);
        dimensionedVector f("f", F_.dimensions(), vector::zero);

        if (subDict.found("d"))
        {
            // d = dimensionedVector("d", subDict.lookup("d"));
            subDict.lookup("d") >> d;
            adjustNegativeResistance(d.value());

        }
        if (subDict.found("f"))
        {
            // f = dimensionedVector("f", subDict.lookup("f"));
            subDict.lookup("f") >> f;
            adjustNegativeResistance(f.value());
        }

        if (D_.dimensions() != d.dimensions())
        {
            FatalIOErrorIn
            (
                "Foam::porousZone::porousZone(const fvMesh&, const Istream&)",
                dict_
            )   << "incorrect dimensions for d: " << d.dimensions()
                << " should be " << D_.dimensions()
                << exit(FatalIOError);
        }

        if (F_.dimensions() != f.dimensions())
        {
            FatalIOErrorIn
            (
                "Foam::porousZone::porousZone(const fvMesh&, const Istream&)",
                dict_
            )   << "incorrect dimensions for f: " << f.dimensions()
                << " should be " << F_.dimensions()
                << exit(FatalIOError);
        }

        D_.value().xx() = d.value().x();
        D_.value().yy() = d.value().y();
        D_.value().zz() = d.value().z();
        D_.value() = (E & D_ & E.T()).value();

        // leading 0.5 is from 1/2 * rho
        F_.value().xx() = 0.5*f.value().x();
        F_.value().yy() = 0.5*f.value().y();
        F_.value().zz() = 0.5*f.value().z();
        F_.value() = (E & F_ & E.T()).value();
    }

    // provide some feedback for the user
    // writeDict(Info, false);

    // it is an error not to define anything
    if
    (
        C0_ <= VSMALL
     && magSqr(D_.value()) <= VSMALL
     && magSqr(F_.value()) <= VSMALL
    )
    {
        FatalIOErrorIn
        (
            "Foam::porousZone::porousZone(const fvMesh&, const Istream&)",
            dict_
        )   << "neither powerLaw (C0/C1) "
               "nor Darcy-Forchheimer law (d/f) specified"
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::porousZone::addResistance(fvVectorMatrix& UEqn) const
{
    if (cellZoneID_ == -1)
    {
        return;
    }

    bool compressible = false;
    if (UEqn.dimensions() == dimensionSet(1, 1, -2, 0, 0))
    {
        compressible = true;
    }

    const labelList& cells = mesh_.cellZones()[cellZoneID_];
    const scalarField& V = mesh_.V();
    scalarField& Udiag = UEqn.diag();
    vectorField& Usource = UEqn.source();
    const vectorField& U = UEqn.psi();

    if (C0_ > VSMALL)
    {
        if (compressible)
        {
            addPowerLawResistance
            (
                Udiag,
                cells,
                V,
                mesh_.lookupObject<volScalarField>("rho"),
                U
            );
        }
        else
        {
            addPowerLawResistance
            (
                Udiag,
                cells,
                V,
                oneField(),
                U
            );
        }
    }

    const tensor& D = D_.value();
    const tensor& F = F_.value();

    if (magSqr(D) > VSMALL || magSqr(F) > VSMALL)
    {
        if (compressible)
        {
            addViscousInertialResistance
            (
                Udiag,
                Usource,
                cells,
                V,
                mesh_.lookupObject<volScalarField>("rho"),
                mesh_.lookupObject<volScalarField>("mu"),
                U
            );
        }
        else
        {
            addViscousInertialResistance
            (
                Udiag,
                Usource,
                cells,
                V,
                oneField(),
                mesh_.lookupObject<volScalarField>("nu"),
                U
            );
        }
    }
}


void Foam::porousZone::addResistance
(
    const fvVectorMatrix& UEqn,
    volTensorField& AU,
    bool correctAUprocBC
) const
{
    if (cellZoneID_ == -1)
    {
        return;
    }

    bool compressible = false;
    if (UEqn.dimensions() == dimensionSet(1, 1, -2, 0, 0))
    {
        compressible = true;
    }

    const labelList& cells = mesh_.cellZones()[cellZoneID_];
    const vectorField& U = UEqn.psi();

    if (C0_ > VSMALL)
    {
        if (compressible)
        {
            addPowerLawResistance
            (
                AU,
                cells,
                mesh_.lookupObject<volScalarField>("rho"),
                U
            );
        }
        else
        {
            addPowerLawResistance
            (
                AU,
                cells,
                oneField(),
                U
            );
        }
    }

    const tensor& D = D_.value();
    const tensor& F = F_.value();

    if (magSqr(D) > VSMALL || magSqr(F) > VSMALL)
    {
        if (compressible)
        {
            addViscousInertialResistance
            (
                AU,
                cells,
                mesh_.lookupObject<volScalarField>("rho"),
                mesh_.lookupObject<volScalarField>("mu"),
                U
            );
        }
        else
        {
            addViscousInertialResistance
            (
                AU,
                cells,
                oneField(),
                mesh_.lookupObject<volScalarField>("nu"),
                U
            );
        }
    }

    if (correctAUprocBC)
    {
        // Correct the boundary conditions of the tensorial diagonal to ensure
        // processor boundaries are correctly handled when AU^-1 is interpolated
        // for the pressure equation.
        AU.correctBoundaryConditions();
    }
}


void Foam::porousZone::writeDict(Ostream& os, bool subDict) const
{
    if (subDict)
    {
        os  << indent << token::BEGIN_BLOCK << incrIndent << nl;
        os.writeKeyword("name") << zoneName() << token::END_STATEMENT << nl;
    }
    else
    {
        os  << indent << zoneName() << nl
            << indent << token::BEGIN_BLOCK << incrIndent << nl;
    }

    coordSys_.writeDict(os, true);

    if (dict_.found("porosity"))
    {
        os.writeKeyword("porosity")
            << porosity()
            << token::END_STATEMENT << nl;
    }

    if (dict_.found("powerLaw"))
    {
        const dictionary& subDict = dict_.subDict("powerLaw");

        os << indent << "powerLaw";
        subDict.write(os);
    }

    if (dict_.found("Darcy"))
    {
        const dictionary& subDict = dict_.subDict("Darcy");

        os << indent << "Darcy";
        subDict.write(os);
    }

    os << decrIndent << indent << token::END_BLOCK << endl;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const porousZone& pZone)
{
    pZone.writeDict(os);
    return os;
}

// ************************************************************************* //
