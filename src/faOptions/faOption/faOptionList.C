/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 OpenCFD Ltd.
------------------------------------------------------------------------------
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

#include "faOptionList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fa
{
    defineTypeNameAndDebug(optionList, 0);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

const Foam::dictionary& Foam::fa::optionList::optionsDict
(
    const dictionary& dict
) const
{
    if (dict.found("options"))
    {
        return dict.subDict("options");
    }
    else
    {
        return dict;
    }
}


bool Foam::fa::optionList::readOptions(const dictionary& dict)
{
    checkTimeIndex_ = mesh_.time().timeIndex() + 2;

    bool allOk = true;
    forAll(*this, i)
    {
        option& bs = this->operator[](i);
        bool ok = bs.read(dict.subDict(bs.name()));
        allOk = (allOk && ok);
    }
    return allOk;
}


void Foam::fa::optionList::checkApplied() const
{
    if (mesh_.time().timeIndex() == checkTimeIndex_)
    {
        forAll(*this, i)
        {
            const option& bs = this->operator[](i);
            bs.checkApplied();
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fa::optionList::optionList
(
    const fvPatch& p,
    const dictionary& dict
)
:
    PtrList<option>(),
    mesh_(p.boundaryMesh().mesh()),
    patch_(p),
    checkTimeIndex_(mesh_.time().startTimeIndex() + 2)
{
    reset(optionsDict(dict));
}


Foam::fa::optionList::optionList(const fvPatch& p)
:
    PtrList<option>(),
    mesh_(p.boundaryMesh().mesh()),
    patch_(p),
    checkTimeIndex_(mesh_.time().startTimeIndex() + 2)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fa::optionList::reset(const dictionary& dict)
{
    // Count number of active faOptions
    label count = 0;
    for (const entry& dEntry : dict)
    {
        if (dEntry.isDict())
        {
            ++count;
        }
    }

    this->resize(count);

    count = 0;
    for (const entry& dEntry : dict)
    {
        if (dEntry.isDict())
        {
            const word& name = dEntry.keyword();
            const dictionary& sourceDict = dEntry.dict();

            this->set
            (
                count++,
                option::New(name, sourceDict, patch_)
            );
        }
    }
}


bool Foam::fa::optionList::read(const dictionary& dict)
{
    return readOptions(optionsDict(dict));
}


bool Foam::fa::optionList::writeData(Ostream& os) const
{
    // Write list contents
    forAll(*this, i)
    {
        os  << nl;
        this->operator[](i).writeHeader(os);
        this->operator[](i).writeData(os);
        this->operator[](i).writeFooter(os);
    }

    // Check state of IOstream
    return os.good();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const fa::optionList& options)
{
    options.writeData(os);
    return os;
}


// ************************************************************************* //
