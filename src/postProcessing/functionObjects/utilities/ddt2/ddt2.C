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

#include "ddt2.H"

#include "volFields.H"
#include "dictionary.H"
#include "FieldFunctions.H"
#include "fvcDdt.H"
#include "steadyStateDdtScheme.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ddt2, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::ddt2::checkFormatName(const word& str)
{
    if (str.find("@@") == string::npos)
    {
        WarningInFunction
            << "Bad result naming "
            << "(no '@@' token found), deactivating."
            << nl << endl;

        return false;
    }
    else if (str == "@@")
    {
        WarningInFunction
            << "Bad result naming "
            << "(only a '@@' token found), deactivating."
            << nl
            << endl;

        return false;
    }
    else
    {
        return true;
    }
}


void Foam::ddt2::uniqWords(wordReList& lst)
{
    boolList retain(lst.size());
    wordHashSet uniq;
    forAll(lst, i)
    {
        const wordRe& select = lst[i];

        retain[i] =
        (
            select.isPattern()
         || uniq.insert(static_cast<const word&>(select))
        );
    }

    inplaceSubset(retain, lst);
}


bool Foam::ddt2::accept(const word& fieldName) const
{
    // check input vs possible result names
    // to avoid circular calculations
    return !blacklist_.match(fieldName);
}


template<class FieldType>
int Foam::ddt2::apply
(
    const fvMesh& mesh,
    const word& inputName,
    int& state
)
{
    // state: return 0 (not-processed), -1 (skip), +1 ok

    // already done, or not available
    if (state || !mesh.foundObject<FieldType>(inputName))
    {
        return state;
    }

    const FieldType& input = mesh.lookupObject<FieldType>(inputName);

    word outputName(resultName_);
    outputName.replace("@@", inputName);

    results_.set(outputName);

    if (!mesh.foundObject<volScalarField>(outputName))
    {
        mesh.objectRegistry::store
        (
            new volScalarField
            (
                IOobject
                (
                    outputName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar
                (
                    "0",
                    (
                        mag_
                      ? mag(input.dimensions()/dimTime)
                      : magSqr(input.dimensions()/dimTime)
                    ),
                    Zero
                )
            )
        );
    }

    volScalarField& output = const_cast<volScalarField&>
    (
        mesh.lookupObject<volScalarField>(outputName)
    );

    if (mag_)
    {
        output = mag(fvc::ddt(input));
    }
    else
    {
        output = magSqr(fvc::ddt(input));
    }

    if (log_)
    {
        // could add additional statistics here
        Info<< type() << " " << name_
            << " field " << output
            << " average: " << gAverage(output) << endl;
    }

    state = +1;
    return state;
}


int Foam::ddt2::process(const fvMesh& mesh, const word& fieldName)
{
    if (!accept(fieldName))
    {
        return -1;
    }

    int state = 0;

    apply<volScalarField>(mesh, fieldName, state);
    apply<volVectorField>(mesh, fieldName, state);

    return state;
}


void Foam::ddt2::process(const fvMesh& mesh)
{
    results_.clear();

    wordHashSet candidates = subsetStrings(selectFields_, mesh.names());
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
            else if (process(mesh, fieldName) < 1)
            {
                ignored.append(fieldName);
            }
        }
    }

    forAllConstIter(wordHashSet, candidates, iter)
    {
        process(mesh, iter.key());
    }

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
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ddt2::ddt2
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    name_(name),
    obr_(obr),
    active_(true),
    selectFields_(),
    resultName_(word::null),
    blacklist_(),
    results_(),
    log_(true),
    mag_(dict.lookupOrDefault<Switch>("mag", false))
{
    // Check if the available mesh is an fvMesh and transient,
    // otherwise deactivate
    if (isA<fvMesh>(obr_))
    {
        const fvMesh& mesh = refCast<const fvMesh>(obr_);

        if (word(mesh.ddtScheme("default")) == "steadyState")
        {
            active_ = false;
            WarningInFunction
                << "Steady-state, deactivating." << nl
                << endl;
        }
    }
    else
    {
        active_ = false;
        WarningInFunction
            << "No fvMesh available, deactivating." << nl
            << endl;
    }

    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ddt2::~ddt2()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ddt2::read(const dictionary& dict)
{
    if (active_)
    {
        log_.readIfPresent("log", dict);

        dict.lookup("fields") >> selectFields_;
        uniqWords(selectFields_);

        resultName_ = dict.lookupOrDefault<word>
        (
            "result",
            ( mag_ ? "mag(ddt(@@))" : "magSqr(ddt(@@))" )
        );

        active_ = checkFormatName(resultName_);
        if (active_)
        {
            blacklist_.set
            (
                string::quotemeta<regExp>
                (
                    resultName_
                ).replace("@@", "(.+)")
            );
        }
        else
        {
            blacklist_.clear();
        }
    }
}


void Foam::ddt2::execute()
{
    results_.clear();
    if (active_)
    {
        process(refCast<const fvMesh>(obr_));
    }
}


void Foam::ddt2::write()
{
    if (active_)
    {
        // consistent output order
        const wordList outputList = results_.sortedToc();

        forAll(outputList, i)
        {
            const word& fieldName = outputList[i];

            if (obr_.foundObject<regIOobject>(fieldName))
            {
                const regIOobject& io =
                    obr_.lookupObject<regIOobject>(fieldName);

                io.write();
            }
        }
    }
}


// ************************************************************************* //
