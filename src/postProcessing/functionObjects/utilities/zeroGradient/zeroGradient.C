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

#include "zeroGradient.H"

#include "volFields.H"
#include "dictionary.H"
#include "symmetryFvPatchField.H"
#include "symmetryPlaneFvPatchField.H"
#include "zeroGradientFvPatchField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(zeroGradient, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::zeroGradient::checkFormatName(const word& str)
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


void Foam::zeroGradient::uniqWords(wordReList& lst)
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


template<class Type>
bool Foam::zeroGradient::accept
(
    const GeometricField<Type, fvPatchField, volMesh>& input
)
{
    const typename GeometricField<Type, fvPatchField, volMesh>
    ::GeometricBoundaryField& patches = input.boundaryField();

    forAll(patches, patchi)
    {
        const fvPatchField<Type>& p = patches[patchi];
        const polyPatch& pp = p.patch().patch();

        bool ignore =
        (
            isA<emptyPolyPatch>(pp)
         || isA<processorPolyPatch>(pp)

         || isA<zeroGradientFvPatchField<Type>>(p)
         || isA<symmetryFvPatchField<Type>>(p)
         || isA<symmetryPlaneFvPatchField<Type>>(p)
        );

        if (!ignore)
        {
            return true;
        }
    }

    return false;
}


template<class Type>
int Foam::zeroGradient::apply
(
    const fvMesh& mesh,
    const word& inputName,
    int& state
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> vfType;

    // state: return 0 (not-processed), -1 (skip), +1 ok

    // already done, or not available
    if (state || !mesh.foundObject<vfType>(inputName))
    {
        return state;
    }

    const vfType& input = mesh.lookupObject<vfType>(inputName);

    if (!returnReduce(accept(input), orOp<bool>()))
    {
        state = -1;
        return state;
    }

    word outputName(resultName_);
    outputName.replace("@@", inputName);

    // also save the field-type, just in case we want it later
    results_.set(outputName, vfType::typeName);

    if (!mesh.foundObject<vfType>(outputName))
    {
        mesh.objectRegistry::store
        (
            new vfType
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
                dimensioned<Type>("0", input.dimensions(), Zero),
                zeroGradientFvPatchField<Type>::typeName
            )
        );
    }

    vfType& output = const_cast<vfType&>(mesh.lookupObject<vfType>(outputName));
    output = input;
    output.correctBoundaryConditions();

    state = +1;
    return state;
}


int Foam::zeroGradient::process(const fvMesh& mesh, const word& fieldName)
{
    int state = 0;
    apply<scalar>(mesh, fieldName, state);
    apply<vector>(mesh, fieldName, state);
    apply<sphericalTensor>(mesh, fieldName, state);
    apply<symmTensor>(mesh, fieldName, state);
    apply<tensor>(mesh, fieldName, state);

    return state;
}


void Foam::zeroGradient::process(const fvMesh& mesh)
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

Foam::zeroGradient::zeroGradient
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
    resultName_(string::null),
    results_(),
    log_(false)
{
    // Check if the available mesh is an fvMesh, otherwise deactivate
    if (!isA<fvMesh>(obr_))
    {
        active_ = false;
        WarningInFunction
            << "No fvMesh available, deactivating." << nl
            << endl;
    }

    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::zeroGradient::~zeroGradient()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::zeroGradient::read(const dictionary& dict)
{
    if (active_)
    {
        log_ = dict.lookupOrDefault<Switch>("log", false);

        dict.lookup("fields") >> selectFields_;
        uniqWords(selectFields_);

        resultName_ = dict.lookupOrDefault<word>
        (
            "result",
            type() + "(@@)"
        );
        active_ = checkFormatName(resultName_);
    }
}


void Foam::zeroGradient::execute()
{
    results_.clear();
    if (active_)
    {
        process(refCast<const fvMesh>(obr_));
    }
}


void Foam::zeroGradient::write()
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

                if (log_)
                {
                    Info<< type() << " " << name_
                        << " output: writing field " << fieldName << endl;
                }

                io.write();
            }
        }
    }
}


// ************************************************************************* //
