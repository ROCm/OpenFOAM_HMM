/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
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

#include "molecule.H"
#include "IOstreams.H"
#include "moleculeCloud.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::molecule::molecule
(
    const Cloud<molecule>& cloud,
    Istream& is,
    bool readFields
)
:
    Particle<molecule>(cloud, is),
    R_(tensor::zero),
    v_(vector::zero),
    a_(vector::zero),
    omega_(vector::zero),
    alpha_(vector::zero),
    siteForces_(List<vector>(0,vector::zero)),
    sitePositions_(List<vector>(0,vector::zero)),
    specialPosition_(vector::zero),
    potentialEnergy_(0.0),
    rf_(tensor::zero),
    special_(0),
    id_(0)
{
    Info<< "Set sizes of siteForces_ and sitePositions_ "
        << "from molCloud reference if possible" << endl;

    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            is >> R_;
            is >> v_;
            is >> a_;
            is >> omega_;
            is >> alpha_;
            is >> siteForces_;
            is >> sitePositions_;
            is >> specialPosition_;
            potentialEnergy_ = readScalar(is);
            is >> rf_;
            special_ = readLabel(is);
            id_ = readLabel(is);
        }
        else
        {
            is.read
            (
                reinterpret_cast<char*>(&R_),
                sizeof(R_)
                + sizeof(v_)
                + sizeof(a_)
                + sizeof(omega_)
                + sizeof(alpha_)
                + sizeof(siteForces_)
                + sizeof(sitePositions_)
                + sizeof(specialPosition_)
                + sizeof(rf_)
                + sizeof(special_)
                + sizeof(id_)
            );
        }
    }

    // Check state of Istream
    is.check
    (
        "Foam::molecule::molecule"
        "(const Cloud<molecule>& cloud, Foam::Istream&), bool"
    );
}


void Foam::molecule::readFields(moleculeCloud& mC)
{
    if (!mC.size())
    {
        return;
    }

    IOField<tensor> R(mC.fieldIOobject("R", IOobject::MUST_READ));
    mC.checkFieldIOobject(mC, R);

    IOField<vector> v(mC.fieldIOobject("v", IOobject::MUST_READ));
    mC.checkFieldIOobject(mC, v);

    IOField<vector> a(mC.fieldIOobject("a", IOobject::MUST_READ));
    mC.checkFieldIOobject(mC, a);

    IOField<vector> omega(mC.fieldIOobject("omega", IOobject::MUST_READ));
    mC.checkFieldIOobject(mC, omega);

    IOField<vector> alpha(mC.fieldIOobject("alpha", IOobject::MUST_READ));
    mC.checkFieldIOobject(mC, alpha);

    IOField<vector> specialPosition
    (
        mC.fieldIOobject("specialPosition", IOobject::MUST_READ)
    );
    mC.checkFieldIOobject(mC, specialPosition);

    IOField<label> special(mC.fieldIOobject("special", IOobject::MUST_READ));
    mC.checkFieldIOobject(mC, special);

    IOField<label> id(mC.fieldIOobject("id", IOobject::MUST_READ));
    mC.checkFieldIOobject(mC, id);

    label i = 0;
    forAllIter(moleculeCloud, mC, iter)
    {
        molecule& mol = iter();

        mol.R_ = R[i];
        mol.v_ = v[i];
        mol.a_ = a[i];
        mol.omega_ = omega[i];
        mol.alpha_ = alpha[i];
        mol.specialPosition_ = specialPosition[i];
        mol.special_ = special[i];
        mol.id_ = id[i];
        i++;
    }
}


void Foam::molecule::writeFields(const moleculeCloud& mC)
{
    Particle<molecule>::writeFields(mC);

    label np = mC.size();

    IOField<tensor> R(mC.fieldIOobject("R", IOobject::NO_READ), np);
    IOField<vector> v(mC.fieldIOobject("v", IOobject::NO_READ), np);
    IOField<vector> a(mC.fieldIOobject("a", IOobject::NO_READ), np);
    IOField<vector> omega(mC.fieldIOobject("omega", IOobject::NO_READ), np);
    IOField<vector> alpha(mC.fieldIOobject("alpha", IOobject::NO_READ), np);
    IOField<vector> specialPosition
    (
        mC.fieldIOobject("specialPosition", IOobject::NO_READ),
        np
    );
    IOField<label> special(mC.fieldIOobject("special", IOobject::NO_READ), np);
    IOField<label> id(mC.fieldIOobject("id", IOobject::NO_READ), np);

    label i = 0;
    forAllConstIter(moleculeCloud, mC, iter)
    {
        const molecule& mol = iter();

        R[i] = mol.R_;
        v[i] = mol.v_;
        a[i] = mol.a_;
        omega[i] = mol.omega_;
        alpha[i] = mol.alpha_;
        specialPosition[i] = mol.specialPosition_;
        special[i] = mol.special_;
        id[i] = mol.id_;
        i++;
    }

    R.write();
    v.write();
    a.write();
    omega.write();
    alpha.write();
    specialPosition.write();
    special.write();
    id.write();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const molecule& mol)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << token::SPACE << static_cast<const Particle<molecule>&>(mol)
            << token::SPACE << mol.face()
            << token::SPACE << mol.stepFraction()
            << token::SPACE << mol.R_
            << token::SPACE << mol.v_
            << token::SPACE << mol.a_
            << token::SPACE << mol.omega_
            << token::SPACE << mol.alpha_
            << token::SPACE << mol.siteForces_
            << token::SPACE << mol.sitePositions_
            << token::SPACE << mol.specialPosition_
            << token::SPACE << mol.potentialEnergy_
            << token::SPACE << mol.rf_
            << token::SPACE << mol.special_
            << token::SPACE << mol.id_;
    }
    else
    {
        os  << static_cast<const Particle<molecule>&>(mol);
        os.write
        (
            reinterpret_cast<const char*>(&mol.R_),
            sizeof(mol.R_)
            + sizeof(mol.v_)
            + sizeof(mol.a_)
            + sizeof(mol.omega_)
            + sizeof(mol.alpha_)
            + sizeof(mol.siteForces_)
            + sizeof(mol.sitePositions_)
            + sizeof(mol.specialPosition_)
            + sizeof(mol.potentialEnergy_)
            + sizeof(mol.rf_)
            + sizeof(mol.special_)
            + sizeof(mol.id_)
        );
    }

    // Check state of Ostream
    os.check
    (
        "Foam::Ostream& Foam::operator<<"
        "(Foam::Ostream&, const Foam::molecule&)"
    );

    return os;
}


// ************************************************************************* //
