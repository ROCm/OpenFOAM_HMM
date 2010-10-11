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
    dipm_(readScalar(dict.lookup("dipm"))),
    omega_(readScalar(dict.lookup("omega"))),
    delta_(readScalar(dict.lookup("delta")))
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


Foam::autoPtr<Foam::liquid> Foam::liquid::New
(
    const dictionary& dict,
    const word& name
)
{
    if (debug)
    {
        Info<< "liquid::New(const dictionary&) : " << "constructing liquid"
            << endl;
    }

    const Switch defaultCoeffs(dict.lookup("defaultCoeffs"));

    if (defaultCoeffs)
    {
        ConstructorTable::iterator cstrIter = ConstructorTablePtr_->find(name);

        if (cstrIter == ConstructorTablePtr_->end())
        {
            FatalErrorIn("liquid::New(const dictionary&, const word&)")
                << "Unknown liquid type "
                << name << nl << nl
                << "Valid liquid types are:" << nl
                << ConstructorTablePtr_->sortedToc()
                << abort(FatalError);
        }

        return autoPtr<liquid>(cstrIter()());
    }
    else
    {
        dictionaryConstructorTable::iterator cstrIter =
            dictionaryConstructorTablePtr_->find(name);

        if (cstrIter == dictionaryConstructorTablePtr_->end())
        {
            FatalErrorIn("liquid::New(const dictionary&, const word&)")
                << "Unknown liquid type "
                << name << nl << nl
                << "Valid liquid types are:" << nl
                << dictionaryConstructorTablePtr_->sortedToc()
                << abort(FatalError);
        }

        return autoPtr<liquid>(cstrIter()(dict.subDict(name + "Coeffs")));
    }
}


// ************************************************************************* //
