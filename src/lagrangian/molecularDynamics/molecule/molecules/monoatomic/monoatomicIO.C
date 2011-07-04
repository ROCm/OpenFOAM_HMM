/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2011 OpenCFD Ltd.
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

#include "monoatomic.H"
#include "IOstreams.H"
#include "Cloud.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::monoatomic::monoatomic
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields
)
:
    particle(mesh, is, readFields),
    v_(vector::zero),
    a_(vector::zero),
    specialPosition_(vector::zero),
    potentialEnergy_(0.0),
    rf_(tensor::zero),
    special_(0),
    id_(0),
    siteForces_(),
    sitePositions_()
{
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            is  >> v_;
            is  >> a_;
            is  >> siteForces_;
            potentialEnergy_ = readScalar(is);
            is  >> rf_;
            special_ = readLabel(is);
            id_ = readLabel(is);
            is  >> sitePositions_;
            is  >> specialPosition_;
        }
        else
        {
            is.read
            (
                reinterpret_cast<char*>(&v_),
                sizeof(v_)
              + sizeof(a_)
              + sizeof(specialPosition_)
              + sizeof(potentialEnergy_)
              + sizeof(rf_)
              + sizeof(special_)
              + sizeof(id_)
            );

            is  >> siteForces_ >> sitePositions_;
        }
    }

    // Check state of Istream
    is.check
    (
        "Foam::monoatomic::monoatomic"
        "("
            "const polyMesh& mesh,"
            "Istream& is,"
            "bool readFields"
        ")"
    );
}


void Foam::monoatomic::readFields(Cloud<monoatomic>& mC)
{
    if (!mC.size())
    {
        return;
    }

    particle::readFields(mC);

    IOField<vector> v(mC.fieldIOobject("v", IOobject::MUST_READ));
    mC.checkFieldIOobject(mC, v);

    IOField<vector> a(mC.fieldIOobject("a", IOobject::MUST_READ));
    mC.checkFieldIOobject(mC, a);

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

    forAllIter(typename Cloud<monoatomic>, mC, iter)
    {
        monoatomic& mol = iter();

        mol.v_ = v[i];
        mol.a_ = a[i];
        mol.specialPosition_ = specialPosition[i];
        mol.special_ = special[i];
        mol.id_ = id[i];
        i++;
    }
}


void Foam::monoatomic::writeFields(const Cloud<monoatomic>& mC)
{
    particle::writeFields(mC);

    label np = mC.size();

    IOField<vector> v(mC.fieldIOobject("v", IOobject::NO_READ), np);
    IOField<vector> a(mC.fieldIOobject("a", IOobject::NO_READ), np);
    IOField<vector> specialPosition
    (
        mC.fieldIOobject("specialPosition", IOobject::NO_READ),
        np
    );
    IOField<label> special(mC.fieldIOobject("special", IOobject::NO_READ), np);
    IOField<label> id(mC.fieldIOobject("id", IOobject::NO_READ), np);

    label i = 0;
    forAllConstIter(typename Cloud<monoatomic>, mC, iter)
    {
        const monoatomic& mol = iter();

        v[i] = mol.v_;
        a[i] = mol.a_;
        specialPosition[i] = mol.specialPosition_;
        special[i] = mol.special_;
        id[i] = mol.id_;
        i++;
    }

    v.write();
    a.write();
    specialPosition.write();
    special.write();
    id.write();
}


void Foam::monoatomic::info(monoatomic::trackingData& td)
{
    vector totalLinearMomentum(vector::zero);
    scalar maxVelocityMag = 0.0;
    scalar totalMass = 0.0;
    scalar totalLinearKE = 0.0;
    scalar totalPE = 0.0;
    scalar totalrDotf = 0.0;
    label nMols = td.cloud().size();

    forAllConstIter(typename Cloud<monoatomic>, td.cloud(), mol)
    {
        const label molId = mol().id();
        scalar molMass(td.cloud().constProps(molId).mass());
        totalMass += molMass;
    }

    forAllConstIter(typename Cloud<monoatomic>, td.cloud(), mol)
    {
        const label molId = mol().id();
        const monoatomic::constantProperties cP
        (
            td.cloud().constProps(molId)
        );
        scalar molMass(cP.mass());
        const vector& molV(mol().v());
        totalLinearMomentum += molV * molMass;
        if (mag(molV) > maxVelocityMag)
        {
            maxVelocityMag = mag(molV);
        }
        totalLinearKE += 0.5*molMass*magSqr(molV);
        totalPE += mol().potentialEnergy();
        totalrDotf += tr(mol().rf());
    }

    scalar meshVolume = sum(td.cloud().mesh().cellVolumes());

    if (Pstream::parRun())
    {
        reduce(totalLinearMomentum, sumOp<vector>());
        reduce(maxVelocityMag, maxOp<scalar>());
        reduce(totalMass, sumOp<scalar>());
        reduce(totalLinearKE, sumOp<scalar>());
        reduce(totalPE, sumOp<scalar>());
        reduce(totalrDotf, sumOp<scalar>());
        reduce(nMols, sumOp<label>());
        reduce(meshVolume, sumOp<scalar>());
    }

    if (nMols)
    {
        Info<< nl << "Number of molecules in " << td.cloud().name() << " = "
            << nMols << nl
            << "    Overall number density = "
            << nMols/meshVolume << nl
            << "    Overall mass density = "
            << totalMass/meshVolume << nl
            << "    Average linear momentum per molecule = "
            << totalLinearMomentum/nMols << ' '
            << mag(totalLinearMomentum)/nMols << nl
            << "    maximum |velocity| = "
            << maxVelocityMag << nl
            << "    Average linear KE per molecule = "
            << totalLinearKE/nMols << nl
            << "    Average angular KE per molecule = "
            << totalPE/nMols << nl
            << "    Average TE per molecule = "
            <<
            (
                  totalLinearKE
                + totalPE
            )
            /nMols
            << nl << endl;
    }
    else
    {
        Info<< nl << "No molecules in " << td.cloud().name() << endl;
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const monoatomic& mol)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << token::SPACE << static_cast<const particle&>(mol)
            << token::SPACE << mol.face()
            << token::SPACE << mol.stepFraction()
            << token::SPACE << mol.v_
            << token::SPACE << mol.a_
            << token::SPACE << mol.specialPosition_
            << token::SPACE << mol.potentialEnergy_
            << token::SPACE << mol.rf_
            << token::SPACE << mol.special_
            << token::SPACE << mol.id_
            << token::SPACE << mol.siteForces_
            << token::SPACE << mol.sitePositions_;
    }
    else
    {
        os  << static_cast<const particle&>(mol);
        os.write
        (
            reinterpret_cast<const char*>(&mol.v_),
            sizeof(mol.v_)
          + sizeof(mol.a_)
          + sizeof(mol.specialPosition_)
          + sizeof(mol.potentialEnergy_)
          + sizeof(mol.rf_)
          + sizeof(mol.special_)
          + sizeof(mol.id_)
        );
        os  << mol.siteForces_ << mol.sitePositions_;
    }

    // Check state of Ostream
    os.check
    (
        "Foam::Ostream& Foam::operator<<"
        "(Foam::Ostream&, const Foam::monoatomic&)"
    );

    return os;
}


// ************************************************************************* //
