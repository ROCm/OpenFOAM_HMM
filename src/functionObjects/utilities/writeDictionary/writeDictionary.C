/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2017 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

#include "writeDictionary.H"
#include "Time.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(writeDictionary, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        writeDictionary,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::writeDictionary::writeHeader()
{
    if (firstChange_)
    {
        Info<< type() << " " << name() << " write:" << nl << endl;

        IOobject::writeDivider(Info);
        Info<< endl;
        firstChange_ = false;
    }
}


void Foam::functionObjects::writeDictionary::checkDictionary
(
    const dictionary& dict,
    const label dicti
)
{
    if (dict.digest() != digests_[dicti])
    {
        writeHeader();

        digests_[dicti] = dict.digest();

        Info<< dict.dictName() << dict << nl;

        IOobject::writeDivider(Info);
        Info<< endl;
    }
}


bool Foam::functionObjects::writeDictionary::tryDirectory
(
    const word& location,
    const label dicti
)
{
    IOobject dictIO
    (
        dictNames_[dicti],
        location,
        obr_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    if (dictIO.typeHeaderOk<IOdictionary>(true))
    {
        checkDictionary(IOdictionary(dictIO), dicti);

        return true;
    }

    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::writeDictionary::writeDictionary
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    regionFunctionObject(name, runTime, dict),
    dictNames_(),
    digests_(),
    firstChange_(true)
{
    read(dict);
    execute();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::writeDictionary::read(const dictionary& dict)
{
    regionFunctionObject::read(dict);

    wordList dictNames(dict.get<wordList>("dictNames"));
    wordHashSet uniqueNames(dictNames);
    dictNames_ = uniqueNames.toc();

    digests_.setSize(dictNames_.size(), SHA1Digest());

    Info<< type() << " " << name() << ": monitoring dictionaries:" << nl;
    if (dictNames_.size())
    {
        for (const word & dictName : dictNames_)
        {
            Info<< "    " << dictName << endl;
        }
    }
    else
    {
        Info<< "    none" << nl;
    }
    Info<< endl;

    return true;
}


bool Foam::functionObjects::writeDictionary::execute()
{
    return true;
}


bool Foam::functionObjects::writeDictionary::write()
{
    firstChange_ = true;

    forAll(dictNames_, dicti)
    {
        const IOdictionary* dictptr =
            obr_.cfindObject<IOdictionary>(dictNames_[dicti]);

        if (dictptr)
        {
            checkDictionary(*dictptr, dicti);
        }
        else
        {
            bool processed = tryDirectory(obr_.time().timeName(), dicti);

            if (!processed)
            {
                processed = tryDirectory(obr_.time().constant(), dicti);
            }

            if (!processed)
            {
                processed = tryDirectory(obr_.time().system(), dicti);
            }

            if (!processed)
            {
                writeHeader();

                Info<< "    Unable to locate dictionary " << dictNames_[dicti]
                    << nl << endl;

                IOobject::writeDivider(Info);
                Info<< endl;
            }
        }
    }

    return true;
}


// ************************************************************************* //
