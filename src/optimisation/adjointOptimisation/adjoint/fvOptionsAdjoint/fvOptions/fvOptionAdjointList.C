/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2007-2019 PCOpt/NTUA
                            | Copyright (C) 2013-2019 FOSS GP
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

#include "fvOptionAdjointList.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(optionAdjointList, 0);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

const Foam::dictionary& Foam::fv::optionAdjointList::optionAdjointsDict
(
    const dictionary& dict
) const
{
    if (dict.found("optionAdjoints"))
    {
        return dict.subDict("optionAdjoints");
    }
    else
    {
        return dict;
    }
}


bool Foam::fv::optionAdjointList::readOptionAdjoints(const dictionary& dict)
{
    checkTimeIndex_ = mesh_.time().timeIndex() + 2;

    bool allOk = true;
    forAll(*this, i)
    {
        optionAdjoint& bs = this->operator[](i);
        bool ok = bs.read(dict.subDict(bs.name()));
        allOk = (allOk && ok);
    }
    return allOk;
}


void Foam::fv::optionAdjointList::checkApplied() const
{
    if (mesh_.time().timeIndex() == checkTimeIndex_)
    {
        forAll(*this, i)
        {
            const optionAdjoint& bs = this->operator[](i);
            bs.checkApplied();
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::optionAdjointList::optionAdjointList
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    PtrList<optionAdjoint>(),
    mesh_(mesh),
    checkTimeIndex_(mesh_.time().startTimeIndex() + 2)
{
    reset(optionAdjointsDict(dict));
}


Foam::fv::optionAdjointList::optionAdjointList(const fvMesh& mesh)
:
    PtrList<optionAdjoint>(),
    mesh_(mesh),
    checkTimeIndex_(mesh_.time().startTimeIndex() + 2)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::optionAdjointList::reset(const dictionary& dict)
{
    label count = 0;
    forAllConstIter(dictionary, dict, iter)
    {
        // safety:
        if (iter().isDict())
        {
            count++;
        }
    }

    this->setSize(count);
    label i = 0;
    forAllConstIter(dictionary, dict, iter)
    {
        if (iter().isDict())
        {
            const word& name = iter().keyword();
            const dictionary& sourceDict = iter().dict();

            this->set
            (
                i++,
                optionAdjoint::New(name, sourceDict, mesh_)
            );
        }
    }
}


bool Foam::fv::optionAdjointList::read(const dictionary& dict)
{
    return readOptionAdjoints(optionAdjointsDict(dict));
}


bool Foam::fv::optionAdjointList::writeData(Ostream& os) const
{
    // Write list contents
    forAll(*this, i)
    {
        os  << nl;
        this->operator[](i).writeData(os);
    }

    // Check state of IOstream
    return os.good();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

namespace Foam
{
    Ostream& operator<<
    (
        Ostream& os,
        const fv::optionAdjointList& optionAdjoints
    )
    {
        optionAdjoints.writeData(os);
        return os;
    }
}


// ************************************************************************* //
