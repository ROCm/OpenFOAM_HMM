/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2017 OpenCFD Ltd.
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

#include "functionObjectList.H"
#include "Time.H"
#include "mapPolyMesh.H"
#include "profiling.H"
#include "argList.H"
#include "timeControlFunctionObject.H"
//#include "IFstream.H"
#include "dictionaryEntry.H"
#include "stringOps.H"
#include "wordRes.H"
#include "Tuple2.H"
#include "etcFiles.H"
#include "IOdictionary.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

Foam::fileName Foam::functionObjectList::functionObjectDictPath
(
    "caseDicts/postProcessing"
);


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::functionObjectList::createStateDict() const
{
    // Cannot set the state dictionary on construction since Time has not
    // been fully initialised
    stateDictPtr_.reset
    (
        new IOdictionary
        (
            IOobject
            (
                "functionObjectProperties",
                time_.timeName(),
                "uniform"/word("functionObjects"),
                time_,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            )
        )
    );
}


Foam::functionObject* Foam::functionObjectList::remove
(
    const word& key,
    label& oldIndex
)
{
    functionObject* ptr = nullptr;

    // Find index of existing functionObject
    HashTable<label>::iterator fnd = indices_.find(key);

    if (fnd != indices_.end())
    {
        oldIndex = fnd();

        // Retrieve the pointer and remove it from the old list
        ptr = this->set(oldIndex, 0).ptr();
        indices_.erase(fnd);
    }
    else
    {
        oldIndex = -1;
    }

    return ptr;
}


void Foam::functionObjectList::listDir
(
    const fileName& dir,
    HashSet<word>& foMap
)
{
    // Search specified directory for functionObject configuration files
    {
        fileNameList foFiles(fileHandler().readDir(dir));
        forAll(foFiles, f)
        {
            if (foFiles[f].ext().empty())
            {
                foMap.insert(foFiles[f]);
            }
        }
    }

    // Recurse into sub-directories
    {
        fileNameList foDirs(fileHandler().readDir(dir, fileName::DIRECTORY));
        forAll(foDirs, fd)
        {
            listDir(dir/foDirs[fd], foMap);
        }
    }
}


void Foam::functionObjectList::list()
{
    HashSet<word> foMap;

    fileNameList etcDirs(findEtcDirs(functionObjectDictPath));

    forAll(etcDirs, ed)
    {
        listDir(etcDirs[ed], foMap);
    }

    Info<< nl
        << "Available configured functionObjects:"
        << foMap.sortedToc()
        << nl;
}


Foam::fileName Foam::functionObjectList::findDict(const word& funcName)
{
    // First check if there is a functionObject dictionary file in the
    // case system directory
    fileName dictFile = stringOps::expand("$FOAM_CASE")/"system"/funcName;

    if (isFile(dictFile))
    {
        return dictFile;
    }
    else
    {
        fileNameList etcDirs(findEtcDirs(functionObjectDictPath));

        forAll(etcDirs, i)
        {
            dictFile = search(funcName, etcDirs[i]);
            if (!dictFile.empty())
            {
                return dictFile;
            }
        }
    }

    return fileName::null;
}


bool Foam::functionObjectList::readFunctionObject
(
    const string& funcNameArgs,
    dictionary& functionsDict,
    HashSet<wordRe>& requiredFields,
    const word& region
)
{
    // Parse the optional functionObject arguments:
    //     'Q(U)' -> funcName = Q; args = (U); field = U
    //
    // Supports named arguments:
    //     'patchAverage(patch=inlet, p)' -> funcName = patchAverage;
    //         args = (patch=inlet, p); field = p

    word funcName(funcNameArgs);

    int argLevel = 0;
    wordReList args;

    List<Tuple2<word, string>> namedArgs;
    bool namedArg = false;
    word argName;

    word::size_type start = 0;
    word::size_type i = 0;

    for
    (
        word::const_iterator iter = funcNameArgs.begin();
        iter != funcNameArgs.end();
        ++iter
    )
    {
        char c = *iter;

        if (c == '(')
        {
            if (argLevel == 0)
            {
                funcName = funcNameArgs.substr(start, i - start);
                start = i+1;
            }
            ++argLevel;
        }
        else if (c == ',' || c == ')')
        {
            if (argLevel == 1)
            {
                if (namedArg)
                {
                    namedArgs.append
                    (
                        Tuple2<word, string>
                        (
                            argName,
                            funcNameArgs.substr(start, i - start)
                        )
                    );
                    namedArg = false;
                }
                else
                {
                    args.append
                    (
                        wordRe
                        (
                            word::validate
                            (
                                funcNameArgs.substr(start, i - start)
                            )
                        )
                    );
                }
                start = i+1;
            }

            if (c == ')')
            {
                if (argLevel == 1)
                {
                    break;
                }
                --argLevel;
            }
        }
        else if (c == '=')
        {
            argName = word::validate
            (
                funcNameArgs.substr(start, i - start)
            );

            start = i+1;
            namedArg = true;
        }

        ++i;
    }

    // Search for the functionObject dictionary
    fileName path = findDict(funcName);

    if (path == fileName::null)
    {
        WarningInFunction
            << "Cannot find functionObject file " << funcName << endl;
        return false;
    }

    // Read the functionObject dictionary
    //IFstream fileStream(path);
    autoPtr<ISstream> fileStreamPtr(fileHandler().NewIFstream(path));
    ISstream& fileStream = fileStreamPtr();

    dictionary funcsDict(fileStream);
    dictionary* funcDictPtr = &funcsDict;

    if (funcsDict.found(funcName) && funcsDict.isDict(funcName))
    {
        funcDictPtr = &funcsDict.subDict(funcName);
    }

    dictionary& funcDict = *funcDictPtr;

    // Insert the 'field' and/or 'fields' entry corresponding to the optional
    // arguments or read the 'field' or 'fields' entry and add the required
    // fields to requiredFields
    if (args.size() == 1)
    {
        funcDict.set("field", args[0]);
        funcDict.set("fields", args);
        requiredFields.insert(args[0]);
    }
    else if (args.size() > 1)
    {
        funcDict.set("fields", args);
        requiredFields.insert(args);
    }
    else if (funcDict.found("field"))
    {
        requiredFields.insert(wordRe(funcDict.lookup("field")));
    }
    else if (funcDict.found("fields"))
    {
        requiredFields.insert(wordReList(funcDict.lookup("fields")));
    }

    // Insert named arguments
    forAll(namedArgs, i)
    {
        IStringStream entryStream
        (
            namedArgs[i].first() + ' ' + namedArgs[i].second() + ';'
        );
        funcDict.set(entry::New(entryStream).ptr());
    }

    // Insert the region name if specified
    if (region != word::null)
    {
        funcDict.set("region", region);
    }

    // Merge this functionObject dictionary into functionsDict
    dictionary funcArgsDict;
    funcArgsDict.add(word::validate(funcNameArgs), funcDict);
    functionsDict.merge(funcArgsDict);

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjectList::functionObjectList
(
    const Time& runTime,
    const bool execution
)
:
    PtrList<functionObject>(),
    digests_(),
    indices_(),
    time_(runTime),
    parentDict_(runTime.controlDict()),
    stateDictPtr_(),
    execution_(execution),
    updated_(false)
{}


Foam::functionObjectList::functionObjectList
(
    const Time& runTime,
    const dictionary& parentDict,
    const bool execution
)
:
    PtrList<functionObject>(),
    digests_(),
    indices_(),
    time_(runTime),
    parentDict_(parentDict),
    stateDictPtr_(),
    execution_(execution),
    updated_(false)
{}


Foam::autoPtr<Foam::functionObjectList> Foam::functionObjectList::New
(
    const argList& args,
    const Time& runTime,
    dictionary& controlDict,
    HashSet<wordRe>& requiredFields
)
{
    autoPtr<functionObjectList> functionsPtr;

    controlDict.add
    (
        dictionaryEntry("functions", controlDict, dictionary::null)
    );

    dictionary& functionsDict = controlDict.subDict("functions");

    word region = word::null;

    // Set the region name if specified
    if (args.optionFound("region"))
    {
        region = args["region"];
    }

    if
    (
        args.optionFound("dict")
     || args.optionFound("func")
     || args.optionFound("funcs")
    )
    {
        if (args.optionFound("dict"))
        {
            controlDict.merge
            (
                IOdictionary
                (
                    IOobject
                    (
                        args["dict"],
                        runTime,
                        IOobject::MUST_READ_IF_MODIFIED
                    )
                )
            );
        }

        if (args.optionFound("func"))
        {
            readFunctionObject
            (
                args["func"],
                functionsDict,
                requiredFields,
                region
            );
        }

        if (args.optionFound("funcs"))
        {
            wordList funcs(args.optionLookup("funcs")());

            forAll(funcs, i)
            {
                readFunctionObject
                (
                    funcs[i],
                    functionsDict,
                    requiredFields,
                    region
                );
            }
        }

        functionsPtr.reset(new functionObjectList(runTime, controlDict));
    }
    else
    {
        functionsPtr.reset(new functionObjectList(runTime));
    }

    functionsPtr->start();

    return functionsPtr;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjectList::~functionObjectList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::functionObjectList::resetState()
{
    // Reset (re-read) the state dictionary
    stateDictPtr_.clear();
    createStateDict();
}


Foam::IOdictionary& Foam::functionObjectList::stateDict()
{
    if (!stateDictPtr_.valid())
    {
        createStateDict();
    }

    return stateDictPtr_();
}


const Foam::IOdictionary& Foam::functionObjectList::stateDict() const
{
    if (!stateDictPtr_.valid())
    {
        createStateDict();
    }

    return stateDictPtr_();
}


void Foam::functionObjectList::clear()
{
    PtrList<functionObject>::clear();
    digests_.clear();
    indices_.clear();
    updated_ = false;
}


Foam::label Foam::functionObjectList::findObjectID(const word& name) const
{
    forAll(*this, objectI)
    {
        if (operator[](objectI).name() == name)
        {
            return objectI;
        }
    }

    return -1;
}


void Foam::functionObjectList::on()
{
    execution_ = true;
}


void Foam::functionObjectList::off()
{
    // For safety, also force a read() when execution is turned back on
    updated_ = execution_ = false;
}


bool Foam::functionObjectList::status() const
{
    return execution_;
}


bool Foam::functionObjectList::start()
{
    return read();
}


bool Foam::functionObjectList::execute()
{
    bool ok = true;

    if (execution_)
    {
        if (!updated_)
        {
            read();
        }

        forAll(*this, objectI)
        {
            const word& objName = operator[](objectI).name();
            {
                addProfiling(fo, "functionObject::" + objName + "::execute");

                ok = operator[](objectI).execute() && ok;
            }

            {
                addProfiling(fo, "functionObject::" + objName + "::write");

                ok = operator[](objectI).write() && ok;
            }
        }
    }

    // Force writing of state dictionary after function object execution
    if (time_.writeTime())
    {
        label oldPrecision = IOstream::precision_;
        IOstream::precision_ = 16;

        stateDictPtr_->writeObject
        (
            IOstream::ASCII,
            IOstream::currentVersion,
            time_.writeCompression(),
            true
        );

        IOstream::precision_ = oldPrecision;
    }

    return ok;
}


bool Foam::functionObjectList::execute(const label subIndex)
{
    bool ok = execution_;

    if (ok)
    {
        forAll(*this, obji)
        {
            functionObject& funcObj = operator[](obji);

            ok = funcObj.execute(subIndex) && ok;
        }
    }

    return ok;
}


bool Foam::functionObjectList::execute
(
    const wordRes& functionNames,
    const label subIndex
)
{
    bool ok = execution_;

    if (ok && functionNames.size())
    {
        forAll(*this, obji)
        {
            functionObject& funcObj = operator[](obji);

            if (functionNames.match(funcObj.name()))
            {
                ok = funcObj.execute(subIndex) && ok;
            }
        }
    }

    return ok;
}


bool Foam::functionObjectList::end()
{
    bool ok = true;

    if (execution_)
    {
        if (!updated_)
        {
            read();
        }

        forAll(*this, objectI)
        {
            const word& objName = operator[](objectI).name();

            addProfiling(fo, "functionObject::" + objName + "::end");

            ok = operator[](objectI).end() && ok;
        }
    }

    return ok;
}


bool Foam::functionObjectList::adjustTimeStep()
{
    bool ok = true;

    if (execution_)
    {
        if (!updated_)
        {
            read();
        }

        forAll(*this, objectI)
        {
            const word& objName = operator[](objectI).name();

            addProfiling(fo, "functionObject::" + objName + "::adjustTimeStep");

            ok = operator[](objectI).adjustTimeStep() && ok;
        }
    }

    return ok;
}


bool Foam::functionObjectList::read()
{
    if (!stateDictPtr_.valid())
    {
        createStateDict();
    }

    bool ok = true;
    updated_ = execution_;

    // Avoid reading/initializing if execution is off
    if (!execution_)
    {
        return true;
    }

    // Update existing and add new functionObjects
    const entry* entryPtr = parentDict_.lookupEntryPtr
    (
        "functions",
        false,
        false
    );

    if (entryPtr)
    {
        PtrList<functionObject> newPtrs;
        List<SHA1Digest> newDigs;
        HashTable<label> newIndices;

        label nFunc = 0;

        addProfiling(fo,"functionObjects::read");

        if (!entryPtr->isDict())
        {
            FatalIOErrorInFunction(parentDict_)
                << "'functions' entry is not a dictionary"
                << exit(FatalIOError);
        }

        const dictionary& functionsDict = entryPtr->dict();

        const_cast<Time&>(time_).libs().open
        (
            functionsDict,
            "libs",
            functionObject::dictionaryConstructorTablePtr_
        );

        newPtrs.setSize(functionsDict.size());
        newDigs.setSize(functionsDict.size());

        forAllConstIter(dictionary, functionsDict, iter)
        {
            const word& key = iter().keyword();

            if (!iter().isDict())
            {
                if (key != "libs")
                {
                    IOWarningInFunction(parentDict_)
                        << "Entry " << key << " is not a dictionary" << endl;
                }

                continue;
            }

            const dictionary& dict = iter().dict();
            bool enabled = dict.lookupOrDefault("enabled", true);

            newDigs[nFunc] = dict.digest();

            label oldIndex;
            functionObject* objPtr = remove(key, oldIndex);

            if (objPtr)
            {
                if (enabled)
                {
                    // Dictionary changed for an existing functionObject
                    if (newDigs[nFunc] != digests_[oldIndex])
                    {
                        addProfiling
                        (
                            fo2,
                            "functionObject::" + objPtr->name() + "::read"
                        );

                        enabled = objPtr->read(dict);
                        ok = enabled && ok;
                    }
                }

                if (!enabled)
                {
                    // Delete the disabled/invalid(read) functionObject
                    delete objPtr;
                    objPtr = nullptr;
                    continue;
                }
            }
            else if (enabled)
            {
                autoPtr<functionObject> foPtr;

                // Throw FatalError, FatalIOError as exceptions
                const bool throwingError = FatalError.throwExceptions();
                const bool throwingIOerr = FatalIOError.throwExceptions();

                try
                {
                    // New functionObject
                    addProfiling
                    (
                        fo2,
                        "functionObject::" + key + "::new"
                    );
                    if (functionObjects::timeControl::entriesPresent(dict))
                    {
                        foPtr.reset
                        (
                            new functionObjects::timeControl(key, time_, dict)
                        );
                    }
                    else
                    {
                        foPtr = functionObject::New(key, time_, dict);
                    }
                }
                catch (Foam::IOerror& ioErr)
                {
                    Info<< ioErr << nl << endl;
                    ::exit(1);
                }
                catch (Foam::error& err)
                {
                    // Bit of trickery to get the original message
                    err.write(Warning, false);
                    InfoInFunction << nl << endl;
                }

                // Restore previous exception throwing state
                FatalError.throwExceptions(throwingError);
                FatalIOError.throwExceptions(throwingIOerr);

                // If one processor only has thrown an exception (so exited the
                // constructor) invalidate the whole functionObject
                if (returnReduce(foPtr.valid(), andOp<bool>()))
                {
                    objPtr = foPtr.ptr();
                }
                else
                {
                    ok = false;
                }
            }

            // Insert active functionObjects into the list
            if (objPtr)
            {
                newPtrs.set(nFunc, objPtr);
                newIndices.insert(key, nFunc);
                nFunc++;
            }
        }

        newPtrs.setSize(nFunc);
        newDigs.setSize(nFunc);

        // Updating the PtrList of functionObjects deletes any
        // existing unused functionObjects
        PtrList<functionObject>::transfer(newPtrs);
        digests_.transfer(newDigs);
        indices_.transfer(newIndices);
    }
    else
    {
        PtrList<functionObject>::clear();
        digests_.clear();
        indices_.clear();
    }

    return ok;
}


bool Foam::functionObjectList::filesModified() const
{
    bool ok = false;
    if (execution_)
    {
        forAll(*this, objectI)
        {
            bool changed = operator[](objectI).filesModified();
            ok = ok || changed;
        }
    }
    return ok;
}


void Foam::functionObjectList::updateMesh(const mapPolyMesh& mpm)
{
    if (execution_)
    {
        forAll(*this, objectI)
        {
            operator[](objectI).updateMesh(mpm);
        }
    }
}


void Foam::functionObjectList::movePoints(const polyMesh& mesh)
{
    if (execution_)
    {
        forAll(*this, objectI)
        {
            operator[](objectI).movePoints(mesh);
        }
    }
}


// ************************************************************************* //
