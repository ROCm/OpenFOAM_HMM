/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
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

#include "particleDistribution.H"
#include "addToRunTimeSelectionTable.H"
#include "general.H"
#include "fvMesh.H"
#include "cloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(particleDistribution, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        particleDistribution,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::particleDistribution::particleDistribution
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(runTime, name, typeName, dict),
    cloudName_("unknown-cloudName"),
    fieldNames_(),
    tagFieldName_("unknown-tagFieldName"),
    distributionBinWidth_(0),
    rndGen_(1234, -1)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::particleDistribution::~particleDistribution()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::particleDistribution::read(const dictionary& dict)
{
    if (fvMeshFunctionObject::read(dict) && writeFile::read(dict))
    {
        dict.lookup("cloud") >> cloudName_;
        dict.lookup("fields") >> fieldNames_;
        dict.lookup("tagField") >> tagFieldName_;
        dict.lookup("distributionBinWidth") >> distributionBinWidth_;

        Info<< type() << " " << name() << " output:"
            << "    Processing cloud : " << cloudName_ << nl
            << endl;

        return true;
    }

    return false;
}


bool Foam::functionObjects::particleDistribution::execute()
{
    return true;
}


bool Foam::functionObjects::particleDistribution::write()
{
    Log << type() << " " << name() << " output:" << endl;

    if (!mesh_.foundObject<cloud>(cloudName_))
    {
        WarningInFunction
            << "Unable to find cloud " << cloudName_
            << " in the mesh database" << endl;

        return false;
    }

    const cloud& c = mesh_.lookupObject<cloud>(cloudName_);

    objectRegistry cloudObr
    (
        IOobject
        (
            name() & "CloudRegistry",
            mesh_.time().timeName(),
            cloud::prefix,
            mesh_.time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    );

    c.writeObjects(cloudObr);

    List<DynamicList<label>> tagAddr;
    if
    (
        tagFieldName_ != "none"
     && cloudObr.foundObject<IOField<scalar>>(tagFieldName_)
    )
    {
        // Tag field present - generate distribution per tag
        const IOField<label>& tag =
            cloudObr.lookupObject<IOField<label>>(tagFieldName_);
        const HashSet<label> tagMap(tag);
        const label tagMax = tagMap.size();

        List<DynamicList<label>> tagAddr(tagMax);
        forAll(tag, i)
        {
            label newTag = tagMap[tag[i]];
            tagAddr[newTag].append(i);
        }
    }

    bool ok = false;
    forAll(fieldNames_, i)
    {
        const word fName = fieldNames_[i];
        ok = ok || processField<scalar>(cloudObr, fName, tagAddr);
        ok = ok || processField<vector>(cloudObr, fName, tagAddr);
        ok = ok || processField<tensor>(cloudObr, fName, tagAddr);
        ok = ok || processField<sphericalTensor>(cloudObr, fName, tagAddr);
        ok = ok || processField<symmTensor>(cloudObr, fName, tagAddr);
        ok = ok || processField<tensor>(cloudObr, fName, tagAddr);

        if (log && !ok)
        {
            WarningInFunction
                << "Unable to find field " << fName << " in the "
                << cloudName_ << " cloud database" << endl;
        }
    }


    return true;
}


void Foam::functionObjects::particleDistribution::generateDistribution
(
    const word& fieldName,
    const scalarField& field,
    const label tag
)
{
    Ostream& os = file();

    if (tag != -1)
    {
        os  << tag << token::TAB;
    }

    distributionModels::general distribution
    (
        field,
        distributionBinWidth_,
        rndGen_
    );

    distribution.writeData(os);

    os  << endl;
}


// ************************************************************************* //
