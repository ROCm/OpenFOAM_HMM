/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2018 OpenCFD Ltd.
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
        IOstreamOption::formatNames.lookupOrFailsafe
        (
            "format",
            dict,
            runTime.writeFormat()
        )
    ),
    caseOpts_(writeOpts_.format()),
    selectFields_(),
    dirName_("ensightWrite"),
    consecutive_(false)
{
    if (postProcess)
    {
        // Disable for post-process mode.
        // Emit as FatalError for the try/catch in the caller.
        FatalError
            << type() << " disabled in post-process mode"
            << exit(FatalError);
    }

    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::ensightWrite::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    //
    // writer options
    //
    writeOpts_.noPatches(dict.lookupOrDefault("noPatches", false));

    if (dict.found("patches"))
    {
        wordRes list(dict.lookup("patches"));
        list.uniq();

        writeOpts_.patchSelection(list);
    }

    if (dict.found("faceZones"))
    {
        wordRes list(dict.lookup("faceZones"));
        list.uniq();

        writeOpts_.faceZoneSelection(list);
    }


    //
    // case options
    //

    caseOpts_.nodeValues(dict.lookupOrDefault("nodeValues", false));

    caseOpts_.width(dict.lookupOrDefault<label>("width", 8));

    // remove existing output directory
    caseOpts_.overwrite(dict.lookupOrDefault("overwrite", false));

    //
    // other options
    //
    dict.readIfPresent("directory", dirName_);
    consecutive_ = dict.lookupOrDefault("consecutive", false);


    //
    // output fields
    //
    dict.lookup("fields") >> selectFields_;
    selectFields_.uniq();

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

    if (consecutive_)
    {
        ensCase().nextTime(t.value());
    }
    else
    {
        ensCase().setTime(t.value(), t.timeIndex());
    }

    bool writeGeom = false;
    if (!ensMesh_.valid())
    {
        writeGeom = true;
        ensMesh_.reset(new ensightMesh(mesh_, writeOpts_));
    }
    if (ensMesh_().needsUpdate())
    {
        writeGeom = true;
        ensMesh_().correct();
    }

    if (writeGeom)
    {
        // Treat all geometry as moving, since we do not know a priori
        // if the simulation has mesh motion later on.
        autoPtr<ensightGeoFile> os = ensCase_().newGeometry(true);
        ensMesh_().write(os);
    }

    Log << type() << " " << name() << " write: (";

    wordHashSet candidates(subsetStrings(selectFields_, mesh_.names()));
    DynamicList<word> missing(selectFields_.size());
    DynamicList<word> ignored(selectFields_.size());

    // check exact matches first
    for (const wordRe& select : selectFields_)
    {
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

    for (const word& cand : candidates)
    {
        process(cand);
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

    ensCase().write();  // Flush case information

    return true;
}


bool Foam::functionObjects::ensightWrite::end()
{
    return true;
}


void Foam::functionObjects::ensightWrite::updateMesh(const mapPolyMesh& mpm)
{
    // fvMeshFunctionObject::updateMesh(mpm);

    if (ensMesh_.valid())
    {
        ensMesh_->expire();
    }
}


void Foam::functionObjects::ensightWrite::movePoints(const polyMesh& mpm)
{
    // fvMeshFunctionObject::updateMesh(mpm);

    if (ensMesh_.valid())
    {
        ensMesh_->expire();
    }
}


// ************************************************************************* //
