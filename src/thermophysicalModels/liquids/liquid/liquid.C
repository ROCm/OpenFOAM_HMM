/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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

#include "liquid.H"
#include "HashTable.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(liquid, 0);
    defineRunTimeSelectionTable(liquid,);
    defineRunTimeSelectionTable(liquid, Istream);
    defineRunTimeSelectionTable(liquid, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::liquid::liquid
(
    scalar W,
    scalar Tc,
    scalar Pc,
    scalar Vc,
    scalar Zc,
    scalar Tt,
    scalar Pt,
    scalar Tb,
    scalar dipm,
    scalar omega,
    scalar delta
)
:
    W_(W),
    Tc_(Tc),
    Pc_(Pc),
    Vc_(Vc),
    Zc_(Zc),
    Tt_(Tt),
    Pt_(Pt),
    Tb_(Tb),
    dipm_(dipm),
    omega_(omega),
    delta_(delta)
{}


Foam::liquid::liquid(Istream& is)
:
    W_(readScalar(is)),
    Tc_(readScalar(is)),
    Pc_(readScalar(is)),
    Vc_(readScalar(is)),
    Zc_(readScalar(is)),
    Tt_(readScalar(is)),
    Pt_(readScalar(is)),
    Tb_(readScalar(is)),
    dipm_(readScalar(is)),
    omega_(readScalar(is)),
    delta_(readScalar(is))
{}


Foam::liquid::liquid(const dictionary& dict)
:
    W_(readScalar(dict.lookup("W"))),
    Tc_(readScalar(dict.lookup("Tc"))),
    Pc_(readScalar(dict.lookup("Pc"))),
    Vc_(readScalar(dict.lookup("Vc"))),
    Zc_(readScalar(dict.lookup("Zc"))),
    Tt_(readScalar(dict.lookup("Tt"))),
    Pt_(readScalar(dict.lookup("Pt"))),
    Tb_(readScalar(dict.lookup("Tb"))),
    dipm_(readScalar(dict.lookup("dipm"))),
    omega_(readScalar(dict.lookup("omega"))),
    delta_(readScalar(dict.lookup("delta")))
{}


Foam::liquid::liquid(const liquid& liq)
:
    W_(liq.W_),
    Tc_(liq.Tc_),
    Pc_(liq.Pc_),
    Vc_(liq.Vc_),
    Zc_(liq.Zc_),
    Tt_(liq.Tt_),
    Pt_(liq.Pt_),
    Tb_(liq.Tb_),
    dipm_(liq.dipm_),
    omega_(liq.omega_),
    delta_(liq.delta_)
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::liquid> Foam::liquid::New(Istream& is)
{
    if (debug)
    {
        Info<< "liquid::New(Istream&) : " << "constructing liquid" << endl;
    }

    const word liquidType(is);
    const word coeffs(is);

    if (coeffs == "defaultCoeffs")
    {
        ConstructorTable::iterator cstrIter =
            ConstructorTablePtr_->find(liquidType);

        if (cstrIter == ConstructorTablePtr_->end())
        {
            FatalErrorIn("liquid::New(Istream&)")
                << "Unknown liquid type "
                << liquidType << nl << nl
                << "Valid liquid types are:" << nl
                << ConstructorTablePtr_->sortedToc()
                << abort(FatalError);
        }

        return autoPtr<liquid>(cstrIter()());
    }
    else if (coeffs == "coeffs")
    {
        IstreamConstructorTable::iterator cstrIter =
            IstreamConstructorTablePtr_->find(liquidType);

        if (cstrIter == IstreamConstructorTablePtr_->end())
        {
            FatalErrorIn("liquid::New(Istream&)")
                << "Unknown liquid type "
                << liquidType << nl << nl
                << "Valid liquid types are:" << nl
                << IstreamConstructorTablePtr_->sortedToc()
                << abort(FatalError);
        }

        return autoPtr<liquid>(cstrIter()(is));
    }
    else
    {
        FatalErrorIn("liquid::New(Istream&)")
            << "liquid type " << liquidType
            << ", option " << coeffs << " given"
            << ", should be coeffs or defaultCoeffs"
            << abort(FatalError);

        return autoPtr<liquid>(NULL);
    }
}


Foam::autoPtr<Foam::liquid> Foam::liquid::New(const dictionary& dict)
{
    if (debug)
    {
        Info<< "liquid::New(const dictionary&) : " << "constructing liquid"
            << endl;
    }

    const word& liquidTypeName = dict.dictName();

    const Switch defaultCoeffs(dict.lookup("defaultCoeffs"));

    if (defaultCoeffs)
    {
        ConstructorTable::iterator cstrIter =
            ConstructorTablePtr_->find(liquidTypeName);

        if (cstrIter == ConstructorTablePtr_->end())
        {
            FatalErrorIn("liquid::New(const dictionary&, const word&)")
                << "Unknown liquid type "
                << liquidTypeName << nl << nl
                << "Valid liquid types are:" << nl
                << ConstructorTablePtr_->sortedToc()
                << abort(FatalError);
        }

        return autoPtr<liquid>(cstrIter()());
    }
    else
    {
        dictionaryConstructorTable::iterator cstrIter =
            dictionaryConstructorTablePtr_->find(liquidTypeName);

        if (cstrIter == dictionaryConstructorTablePtr_->end())
        {
            FatalErrorIn("liquid::New(const dictionary&, const word&)")
                << "Unknown liquid type "
                << liquidTypeName << nl << nl
                << "Valid liquid types are:" << nl
                << dictionaryConstructorTablePtr_->sortedToc()
                << abort(FatalError);
        }

        return autoPtr<liquid>
        (
            cstrIter()(dict.subDict(liquidTypeName + "Coeffs"))
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::liquid::rho(scalar p, scalar T) const
{
    notImplemented("Foam::scalar Foam::liquid::rho(scalar, scalar) const");
    return 0.0;
}


Foam::scalar Foam::liquid::pv(scalar p, scalar T) const
{
    notImplemented("Foam::scalar Foam::liquid::pv(scalar, scalar) const");
    return 0.0;
}


Foam::scalar Foam::liquid::hl(scalar p, scalar T) const
{
    notImplemented("Foam::scalar Foam::liquid::hl(scalar, scalar) const");
    return 0.0;
}


Foam::scalar Foam::liquid::Cp(scalar p, scalar T) const
{
    notImplemented("Foam::scalar Foam::liquid::Cp(scalar, scalar) const");
    return 0.0;
}


Foam::scalar Foam::liquid::h(scalar p, scalar T) const
{
    notImplemented("Foam::scalar Foam::liquid::h(scalar, scalar) const");
    return 0.0;
}


Foam::scalar Foam::liquid::Cpg(scalar p, scalar T) const
{
    notImplemented("Foam::scalar Foam::liquid::Cpg(scalar, scalar) const");
    return 0.0;
}


Foam::scalar Foam::liquid::mu(scalar p, scalar T) const
{
    notImplemented("Foam::scalar Foam::liquid::mu(scalar, scalar) const");
    return 0.0;
}


Foam::scalar Foam::liquid::mug(scalar p, scalar T) const
{
    notImplemented("Foam::scalar Foam::liquid::mug(scalar, scalar) const");
    return 0.0;
}


Foam::scalar Foam::liquid::K(scalar p, scalar T) const
{
    notImplemented("Foam::scalar Foam::liquid::K(scalar, scalar) const");
    return 0.0;
}


Foam::scalar Foam::liquid::Kg(scalar p, scalar T) const
{
    notImplemented("Foam::scalar Foam::liquid::Kg(scalar, scalar) const");
    return 0.0;
}


Foam::scalar Foam::liquid::sigma(scalar p, scalar T) const
{
    notImplemented("Foam::scalar Foam::liquid::sigms(scalar, scalar) const");
    return 0.0;
}


Foam::scalar Foam::liquid::D(scalar p, scalar T) const
{
    notImplemented("Foam::scalar Foam::liquid::D(scalar, scalar) const");
    return 0.0;
}


Foam::scalar Foam::liquid::D(scalar p, scalar T, scalar Wb) const
{
    notImplemented("Foam::scalar Foam::liquid::D(scalar, scalar) const");
    return 0.0;
}


void Foam::liquid::writeData(Ostream& os) const
{

    os  << W_ << token::SPACE
        << Tc_ << token::SPACE
        << Pc_ << token::SPACE
        << Vc_ << token::SPACE
        << Zc_ << token::SPACE
        << Tt_ << token::SPACE
        << Pt_ << token::SPACE
        << Tb_ << token::SPACE
        << dipm_ << token::SPACE
        << omega_<< token::SPACE
        << delta_;
}


// ************************************************************************* //
