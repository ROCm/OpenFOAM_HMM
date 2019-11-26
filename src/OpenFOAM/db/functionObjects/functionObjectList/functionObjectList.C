/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2019 OpenCFD Ltd.
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
#include "dictionaryEntry.H"
#include "stringOps.H"
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


void Foam::functionObjectList::createOutputRegistry() const
{
    objectsRegistryPtr_.reset
    (
        new objectRegistry
        (
            IOobject
            (
                "functionObjectObjects",
                time_.timeName(),
                time_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            )
        )
    );
}


Foam::autoPtr<Foam::functionObject> Foam::functionObjectList::remove
(
    const word& key,
    label& oldIndex
)
{
    autoPtr<functionObject> oldptr;

    auto iter = indices_.find(key);  // Index of existing functionObject

    if (iter.found())
    {
        oldIndex = *iter;

        // Remove pointer from the old list
        oldptr = this->release(oldIndex);
        indices_.erase(iter);
    }
    else
    {
        oldIndex = -1;
    }

    return oldptr;
}


void Foam::functionObjectList::listDir
(
    const fileName& dir,
    wordHashSet& available
)
{
    // Search specified directory for functionObject configuration files
    for (const fileName& f : fileHandler().readDir(dir))
    {
        if (f.ext().empty())
        {
            available.insert(f);
        }
    }

    // Recurse into sub-directories
    for (const fileName& d : fileHandler().readDir(dir, fileName::DIRECTORY))
    {
        listDir(dir/d, available);
    }
}


void Foam::functionObjectList::list()
{
    wordHashSet available;

    for (const fileName& d : findEtcDirs(functionObjectDictPath))
    {
        listDir(d, available);
    }

    Info<< nl
        << "Available configured functionObjects:"
        << available.sortedToc()
        << nl;
}


Foam::fileName Foam::functionObjectList::findDict(const word& funcName)
{
    // First check for functionObject dictionary file in globalCase system/

    fileName dictFile = stringOps::expand("<system>")/funcName;

    if (isFile(dictFile))
    {
        return dictFile;
    }

    for (const fileName& d : findEtcDirs(functionObjectDictPath))
    {
        dictFile = search(funcName, d);
        if (!dictFile.empty())
        {
            return dictFile;
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
    bool hasNamedArg = false;
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
                if (hasNamedArg)
                {
                    namedArgs.append
                    (
                        Tuple2<word, string>
                        (
                            argName,
                            funcNameArgs.substr(start, i - start)
                        )
                    );
                    hasNamedArg = false;
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
            hasNamedArg = true;
        }

        ++i;
    }

    // Search for the functionObject dictionary
    fileName path = functionObjectList::findDict(funcName);

    if (path == fileName::null)
    {
        WarningInFunction
            << "Cannot find functionObject file " << funcName << endl;
        return false;
    }

    // Read the functionObject dictionary
    autoPtr<ISstream> fileStreamPtr(fileHandler().NewIFstream(path));
    ISstream& fileStream = fileStreamPtr();

    dictionary funcsDict(fileStream);
    dictionary* funcDictPtr = funcsDict.findDict(funcName);
    dictionary& funcDict = (funcDictPtr ? *funcDictPtr : funcsDict);


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
        requiredFields.insert(funcDict.get<wordRe>("field"));
    }
    else if (funcDict.found("fields"))
    {
        requiredFields.insert(funcDict.get<wordRes>("fields"));
    }

    // Insert named arguments
    for (const Tuple2<word, string>& namedArg : namedArgs)
    {
        IStringStream entryStream
        (
            namedArg.first() + ' ' + namedArg.second() + ';'
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
    objectsRegistryPtr_(),
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
    objectsRegistryPtr_(),
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
    // Merge any functions from the provided controlDict
    controlDict.add
    (
        dictionaryEntry("functions", controlDict, dictionary::null),
        true
    );

    dictionary& functionsDict = controlDict.subDict("functions");

    const word regionName = args.get<word>("region", "");

    bool modifiedControlDict = false;

    if (args.found("dict"))
    {
        modifiedControlDict = true;

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

    if (args.found("func"))
    {
        modifiedControlDict = true;

        readFunctionObject
        (
            args["func"],
            functionsDict,
            requiredFields,
            regionName
        );
    }

    if (args.found("funcs"))
    {
        modifiedControlDict = true;

        wordList funcNames = args.getList<word>("funcs");

        for (const word& funcName : funcNames)
        {
            readFunctionObject
            (
                funcName,
                functionsDict,
                requiredFields,
                regionName
            );
        }
    }


    autoPtr<functionObjectList> functionsPtr;

    if (modifiedControlDict)
    {
        functionsPtr.reset(new functionObjectList(runTime, controlDict));
    }
    else
    {
        functionsPtr.reset(new functionObjectList(runTime));
    }

    functionsPtr->start();

    return functionsPtr;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::functionObjectList::triggerIndex() const
{
    label triggeri = labelMin;
    stateDict().readIfPresent("triggerIndex", triggeri);

    return triggeri;
}


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

    return *stateDictPtr_;
}


const Foam::IOdictionary& Foam::functionObjectList::stateDict() const
{
    if (!stateDictPtr_.valid())
    {
        createStateDict();
    }

    return *stateDictPtr_;
}


Foam::objectRegistry& Foam::functionObjectList::storedObjects()
{
    if (!objectsRegistryPtr_.valid())
    {
        createOutputRegistry();
    }

    return *objectsRegistryPtr_;
}


const Foam::objectRegistry& Foam::functionObjectList::storedObjects() const
{
    if (!objectsRegistryPtr_.valid())
    {
        createOutputRegistry();
    }

    return *objectsRegistryPtr_;
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
    label id = 0;

    for (const functionObject& funcObj : functions())
    {
        if (funcObj.name() == name)
        {
            return id;
        }

        ++id;
    }

    return -1;
}


void Foam::functionObjectList::on()
{
    execution_ = true;
}


void Foam::functionObjectList::off()
{
    // For safety, also force a read() when execution is resumed
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

        for (functionObject& funcObj : functions())
        {
            const word& objName = funcObj.name();
            {
                addProfiling(fo, "functionObject::" + objName + "::execute");

                ok = funcObj.execute() && ok;
            }

            {
                addProfiling(fo, "functionObject::" + objName + "::write");

                ok = funcObj.write() && ok;
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
        for (functionObject& funcObj : functions())
        {
            ok = funcObj.execute(subIndex) && ok;
        }
    }

    return ok;
}


bool Foam::functionObjectList::execute
(
    const UList<wordRe>& functionNames,
    const label subIndex
)
{
    bool ok = execution_;

    if (ok && functionNames.size())
    {
        for (functionObject& funcObj : functions())
        {
            if (stringOps::match(functionNames, funcObj.name()))
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

        for (functionObject& funcObj : functions())
        {
            const word& objName = funcObj.name();

            addProfiling(fo, "functionObject::" + objName + "::end");

            ok = funcObj.end() && ok;
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

        for (functionObject& funcObj : functions())
        {
            const word& objName = funcObj.name();

            addProfiling(fo, "functionObject::" + objName + "::adjustTimeStep");

            ok = funcObj.adjustTimeStep() && ok;
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

    updated_ = execution_;

    // Avoid reading/initializing if execution is off
    if (!execution_)
    {
        return true;
    }

    // Update existing and add new functionObjects
    const entry* entryPtr =
        parentDict_.findEntry("functions", keyType::LITERAL);

    bool ok = true;

    if (!entryPtr)
    {
        // No functions
        PtrList<functionObject>::clear();
        digests_.clear();
        indices_.clear();
    }
    else if (!entryPtr->isDict())
    {
        // Bad entry type
        ok = false;
        FatalIOErrorInFunction(parentDict_)
            << "'functions' entry is not a dictionary"
            << exit(FatalIOError);
    }
    else
    {
        const dictionary& functionsDict = entryPtr->dict();

        PtrList<functionObject> newPtrs(functionsDict.size());
        List<SHA1Digest> newDigs(functionsDict.size());
        HashTable<label> newIndices;

        addProfiling(fo, "functionObjects::read");

        const_cast<Time&>(time_).libs().open
        (
            functionsDict,
            "libs",
            functionObject::dictionaryConstructorTablePtr_
        );

        label nFunc = 0;

        for (const entry& dEntry : functionsDict)
        {
            const word& key = dEntry.keyword();

            if (!dEntry.isDict())
            {
                if (key != "libs")
                {
                    IOWarningInFunction(parentDict_)
                        << "Entry " << key << " is not a dictionary" << endl;
                }

                continue;
            }

            const dictionary& dict = dEntry.dict();

            bool enabled = dict.lookupOrDefault("enabled", true);

            newDigs[nFunc] = dict.digest();

            label oldIndex = -1;
            autoPtr<functionObject> objPtr = remove(key, oldIndex);

            if (objPtr)
            {
                // Re-read if dictionary content changed for
                // existing functionObject
                if (enabled && newDigs[nFunc] != digests_[oldIndex])
                {
                    addProfiling
                    (
                        fo2,
                        "functionObject::" + objPtr->name() + "::read"
                    );

                    if (functionObjects::timeControl::entriesPresent(dict))
                    {
                        if (isA<functionObjects::timeControl>(objPtr()))
                        {
                            // Already a time control - normal read
                            enabled = objPtr->read(dict);
                        }
                        else
                        {
                            // Was not a time control - need to re-create
                            objPtr.reset
                            (
                                new functionObjects::timeControl
                                (
                                    key,
                                    time_,
                                    dict
                                )
                            );

                            enabled = true;
                        }
                    }
                    else
                    {
                        // Plain function object - normal read
                        enabled = objPtr->read(dict);
                    }

                    ok = enabled && ok;
                }

                if (!enabled)
                {
                    // Delete disabled or an invalid(read) functionObject
                    objPtr.clear();
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
                catch (const Foam::IOerror& ioErr)
                {
                    Info<< ioErr << nl << endl;
                    std::exit(1);
                }
                catch (const Foam::error& err)
                {
                    // Bit of trickery to get the original message
                    err.write(Warning, false);
                    InfoInFunction
                        << nl
                        << "--> while loading function object '" << key << "'"
                        << nl << endl;
                }

                // Restore previous exception throwing state
                FatalError.throwExceptions(throwingError);
                FatalIOError.throwExceptions(throwingIOerr);

                // Required functionObject to be valid on all processors
                if (returnReduce(foPtr.valid(), andOp<bool>()))
                {
                    objPtr.reset(foPtr.release());
                }
                else
                {
                    ok = false;
                }
            }

            // Insert active functionObject into the list
            if (objPtr)
            {
                newPtrs.set(nFunc, objPtr);
                newIndices.insert(key, nFunc);
                ++nFunc;
            }
        }

        newPtrs.resize(nFunc);
        newDigs.resize(nFunc);

        // Updating PtrList of functionObjects deletes any
        // existing unused functionObjects
        PtrList<functionObject>::transfer(newPtrs);
        digests_.transfer(newDigs);
        indices_.transfer(newIndices);
    }

    return ok;
}


bool Foam::functionObjectList::filesModified() const
{
    bool ok = false;
    if (execution_)
    {
        for (const functionObject& funcObj : functions())
        {
            bool changed = funcObj.filesModified();
            ok = ok || changed;
        }
    }
    return ok;
}


void Foam::functionObjectList::updateMesh(const mapPolyMesh& mpm)
{
    if (execution_)
    {
        for (functionObject& funcObj : functions())
        {
            funcObj.updateMesh(mpm);
        }
    }
}


void Foam::functionObjectList::movePoints(const polyMesh& mesh)
{
    if (execution_)
    {
        for (functionObject& funcObj : functions())
        {
            funcObj.movePoints(mesh);
        }
    }
}


// ************************************************************************* //
