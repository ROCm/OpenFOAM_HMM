/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 OpenCFD Ltd.
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

#include "ensightWrite.H"
#include "Time.H"
#include "polyMesh.H"
#include "wordRes.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(ensightWrite, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        ensightWrite,
        dictionary
    );
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

int Foam::functionObjects::ensightWrite::process(const word& fieldName)
{
    int state = 0;

    writeVolField<scalar>(fieldName, state);
    writeVolField<vector>(fieldName, state);
    writeVolField<sphericalTensor>(fieldName, state);
    writeVolField<symmTensor>(fieldName, state);
    writeVolField<tensor>(fieldName, state);

    return state;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::ensightWrite::ensightWrite
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeOpts_
    (
        dict.found("format")
      ? IOstream::formatEnum(dict.lookup("format"))
      : runTime.writeFormat()
    ),
    caseOpts_(writeOpts_.format()),
    selectFields_(),
    dirName_("ensightWrite"),
    consecutive_(false)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::ensightWrite::~ensightWrite()
{
    if (ensCase_.valid())
    {
        // finalize case
        ensCase().write();
        ensCase_.clear();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::ensightWrite::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    //
    // writer options
    //
    writeOpts_.noPatches
    (
        dict.lookupOrDefault<Switch>("noPatches", false)
    );

    if (dict.found("patches"))
    {
        wordReList lst(dict.lookup("patches"));
        wordRes::inplaceUniq(lst);

        writeOpts_.patchSelection(lst);
    }

    if (dict.found("faceZones"))
    {
        wordReList lst(dict.lookup("faceZones"));
        wordRes::inplaceUniq(lst);

        writeOpts_.faceZoneSelection(lst);
    }


    //
    // case options
    //
    caseOpts_.width(dict.lookupOrDefault<label>("width", 8));

    // remove existing output directory
    caseOpts_.overwrite(dict.lookupOrDefault<Switch>("overwrite", false));


    //
    // other options
    //
    dict.readIfPresent("directory", dirName_);
    consecutive_ = dict.lookupOrDefault<Switch>("consecutive", false);


    //
    // output fields
    //
    dict.lookup("fields") >> selectFields_;
    wordRes::inplaceUniq(selectFields_);

    return true;
}


bool Foam::functionObjects::ensightWrite::execute()
{
    return true;
}


bool Foam::functionObjects::ensightWrite::write()
{
    const Time& t = obr_.time();

    if (!ensCase_.valid())
    {
        // Define sub-directory name to use for EnSight data.
        // The path to the ensight directory is at case level only
        // - For parallel cases, data only written from master

        fileName ensightDir = dirName_;
        if (!ensightDir.isAbsolute())
        {
            ensightDir = t.rootPath()/t.globalCaseName()/ensightDir;
        }

        ensCase_.reset
        (
            new ensightCase
            (
                ensightDir,
                t.globalCaseName(),
                caseOpts_
            )
        );
    }

    if (!ensMesh_.valid())
    {
        ensMesh_.reset(new ensightMesh(mesh_, writeOpts_));

        if (ensMesh_().needsUpdate())
        {
            ensMesh_().correct();
        }

        // assume static geometry - need to fix later
        autoPtr<ensightGeoFile> os = ensCase_().newGeometry(false);
        ensMesh_().write(os);
    }
    else if (ensMesh_().needsUpdate())
    {
        // appears to have moved
        ensMesh_().correct();

        autoPtr<ensightGeoFile> os = ensCase_().newGeometry(true);
        ensMesh_().write(os);
    }

    Log << type() << " " << name() << " write: (";

    if (consecutive_)
    {
        ensCase().nextTime(t.value());
    }
    else
    {
        ensCase().setTime(t.value(), t.timeIndex());
    }

    wordHashSet candidates(subsetStrings(selectFields_, mesh_.names()));
    DynamicList<word> missing(selectFields_.size());
    DynamicList<word> ignored(selectFields_.size());

    // check exact matches first
    forAll(selectFields_, i)
    {
        const wordRe& select = selectFields_[i];
        if (!select.isPattern())
        {
            const word& fieldName = static_cast<const word&>(select);

            if (!candidates.erase(fieldName))
            {
                missing.append(fieldName);
            }
            else if (process(fieldName) < 1)
            {
                ignored.append(fieldName);
            }
        }
    }

    forAllConstIter(wordHashSet, candidates, iter)
    {
        process(iter.key());
    }

    Log << " )" << endl;

    if (missing.size())
    {
        WarningInFunction
            << "Missing field " << missing << endl;
    }
    if (ignored.size())
    {
        WarningInFunction
            << "Unprocessed field " << ignored << endl;
    }

    return true;
}


bool Foam::functionObjects::ensightWrite::end()
{
    if (ensCase_.valid())
    {
        // finalize case
        ensCase().write();
        ensCase_.clear();
    }

    return true;
}


void Foam::functionObjects::ensightWrite::updateMesh(const mapPolyMesh& mpm)
{
    // fvMeshFunctionObject::updateMesh(mpm);

    if (ensMesh_.valid())
    {
        ensMesh_().expire();
    }
}


void Foam::functionObjects::ensightWrite::movePoints(const polyMesh& mpm)
{
    // fvMeshFunctionObject::updateMesh(mpm);

    if (ensMesh_.valid())
    {
        ensMesh_().expire();
    }
}


// ************************************************************************* //
