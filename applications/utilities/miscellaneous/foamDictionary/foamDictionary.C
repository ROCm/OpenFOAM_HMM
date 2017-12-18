/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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

Application
    foamDictionary

Description
    Interrogates and manipulates dictionaries.

Usage
    \b foamDictionary [OPTION] dictionary

      - \par -entry \<name\>
        Selects an entry

      - \par -keywords
        Prints the keywords (of the selected entry or of the top level if
        no entry was selected

      - \par -add \<value\>
        Adds the entry (should not exist yet)

      - \par -set \<value\>
        Adds or replaces the entry

      - \par -remove
        Remove the selected entry

      - \par -diff \<dictionary\>
        Write differences with respect to the specified dictionary
        (or sub entry if -entry specified)

      - \par -diffEtc \<dictionary\>
        Write differences with respect to the specified dictionary
        (or sub entry if -entry specified)

      - \par -expand
        Read the specified dictionary file, expand the macros etc. and write
        the resulting dictionary to standard output.

      - \par -includes
        List the \c \#include and \c \#includeIfPresent files to standard output

      - \par -disableFunctionEntries
        Do not expand macros or directives (\#include etc)

    Example usage:
      - Change simulation to run for one timestep only:
        \verbatim
          foamDictionary system/controlDict -entry stopAt -set writeNow
        \endverbatim

      - Change solver:
        \verbatim
           foamDictionary system/fvSolution -entry solvers.p.solver -set PCG
        \endverbatim

      - Print bc type:
        \verbatim
           foamDictionary 0/U -entry boundaryField.movingWall.type
        \endverbatim

      - Change bc parameter:
        \verbatim
           foamDictionary 0/U -entry boundaryField.movingWall.value \
             -set "uniform (2 0 0)"
        \endverbatim

      - Change whole bc type:
        \verbatim
          foamDictionary 0/U -entry boundaryField.movingWall \
            -set "{type uniformFixedValue; uniformValue (2 0 0);}"
        \endverbatim

      - Write the differences with respect to a template dictionary:
        \verbatim
          foamDictionary 0/U -diffEtc templates/closedVolume/0/U
        \endverbatim

      - Write the differences in boundaryField with respect to a
        template dictionary:
        \verbatim
          foamDictionary 0/U -diffEtc templates/closedVolume/0/U \
            -entry boundaryField
        \endverbatim

      - Change patch type:
        \verbatim
          foamDictionary constant/polyMesh/boundary \
            -entry entry0.fixedWalls.type -set patch
        \endverbatim
        This uses special parsing of Lists which stores these in the
        dictionary with keyword 'entryDDD' where DDD is the position
        in the dictionary (after ignoring the FoamFile entry).

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "profiling.H"
#include "Time.H"
#include "Fstream.H"
#include "etcFiles.H"
#include "includeEntry.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Convert older ':' scope syntax to newer '.' scope syntax,
//  but leave anything with '/' delimiters untouched
bool upgradeScope(word& entryName)
{
    if
    (
        entryName.find('/') == string::npos
     && entryName.find(':') != string::npos
    )
    {
        const wordList names(fileName(entryName).components(':'));

        entryName.resize(0);

        for (const word& name : names)
        {
            if (entryName.size()) entryName.append(".");

            entryName.append(name);
        }

        return true;
    }

    // Nothing changed
    return false;
}


//- Split into dictionary name and the entry name
class dictAndKeyword
{
    word dict_;
    word key_;

public:
    dictAndKeyword(const word& scopedName)
    {
        string::size_type i = scopedName.rfind('/');
        if (i == string::npos)
        {
            i = scopedName.rfind('.');
        }

        if (i != string::npos)
        {
            dict_ = scopedName.substr(0, i);
            key_  = scopedName.substr(i+1);
        }
        else
        {
            key_  = scopedName;
        }
    }

    inline const word& dict() const
    {
        return dict_;
    }

    inline const word& key() const
    {
        return key_;
    }
};


const dictionary& lookupScopedDict
(
    const dictionary& dict,
    const word& subDictName
)
{
    if (subDictName.empty())
    {
        return dict;
    }

    const entry* eptr = dict.lookupScopedEntryPtr(subDictName, false, false);

    if (!eptr || !eptr->isDict())
    {
        FatalIOErrorInFunction(dict)
            << "keyword " << subDictName
            << " is undefined in dictionary "
            << dict.name() << " or is not a dictionary"
            << endl
            << "Valid keywords are " << dict.keys()
            << exit(FatalIOError);
    }

    return eptr->dict();
}


void removeDict(dictionary& dict, const dictionary& dictToRemove)
{
    for (const entry& refEntry : dictToRemove)
    {
        auto finder = dict.search(refEntry.keyword(), false, false);

        bool purge = false;

        if (finder.isDict())
        {
            if (refEntry.isDict())
            {
                removeDict(finder.dict(), refEntry.dict());

                // Purge if dictionary is empty
                purge = finder.dict().empty();
            }
        }
        else if (finder.found() && !refEntry.isDict())
        {
            // Purge if entries match
            purge = (finder.ref() == refEntry);
        }

        if (purge)
        {
            dict.remove(refEntry.keyword());
        }
    }
}


int main(int argc, char *argv[])
{
    argList::addNote("manipulates dictionaries");

    argList::noBanner();
    argList::noJobInfo();
    argList::addArgument("dictionary");
    argList::addBoolOption("keywords", "List keywords");
    argList::addOption("entry", "name", "report/select the named entry");
    argList::addBoolOption
    (
        "value",
        "Print entry value"
    );
    argList::addOption
    (
        "set",
        "value",
        "Set entry value or add new entry"
    );
    argList::addOption
    (
        "add",
        "value",
        "Add a new entry"
    );
    argList::addBoolOption
    (
        "remove",
        "Remove the entry."
    );
    argList::addOption
    (
        "diff",
        "dict",
        "Write differences with respect to the specified dictionary"
    );
    argList::addOption
    (
        "diffEtc",
        "dict",
        "As per -diff, but locate the file as per foamEtcFile"
    );
    argList::addBoolOption
    (
        "includes",
        "List the #include/#includeIfPresent files to standard output"
    );
    argList::addBoolOption
    (
        "expand",
        "Read the specified dictionary file, expand the macros etc. and write "
        "the resulting dictionary to standard output"
    );
    argList::addBoolOption
    (
        "disableFunctionEntries",
        "Disable expansion of dictionary directives - #include, #codeStream etc"
    );
    profiling::disable(); // Disable profiling (and its output)

    argList args(argc, argv);

    const bool listIncludes = args.optionFound("includes");

    if (listIncludes)
    {
        Foam::functionEntries::includeEntry::log = true;
    }

    const bool disableEntries = args.optionFound("disableFunctionEntries");
    if (disableEntries)
    {
        Info<< "Not expanding variables or dictionary directives" << endl;
        entry::disableFunctionEntries = true;
    }


    const fileName dictFileName(args[1]);

    autoPtr<IFstream> dictFile(new IFstream(dictFileName));
    if (!dictFile().good())
    {
        FatalErrorInFunction
            << "Cannot open file " << dictFileName
            << exit(FatalError, 1);
    }


    bool changed = false;

    // Read but preserve headers
    dictionary dict;
    dict.read(dictFile(), true);

    if (listIncludes)
    {
        return 0;
    }
    else if (args.optionFound("expand"))
    {
        IOobject::writeBanner(Info)
            <<"//\n// " << dictFileName << "\n//\n";
        dict.write(Info, false);
        IOobject::writeDivider(Info);

        return 0;
    }


    // Has "diff" or "diffEtc"
    bool optDiff = false;

    // Reference dictionary for -diff / -diffEtc
    dictionary diffDict;
    {
        fileName diffFileName;
        if (args.optionReadIfPresent("diff", diffFileName))
        {
            IFstream diffFile(diffFileName);
            if (!diffFile.good())
            {
                FatalErrorInFunction
                    << "Cannot open file " << diffFileName
                    << exit(FatalError, 1);
            }

            // Read but preserve headers
            diffDict.read(diffFile, true);
            optDiff = true;
        }
        else if (args.optionReadIfPresent("diffEtc", diffFileName))
        {
            fileName foundName = findEtcFile(diffFileName);
            if (foundName.empty())
            {
                FatalErrorInFunction
                    << "Cannot find etcFile " << diffFileName
                    << exit(FatalError, 1);
            }

            IFstream diffFile(foundName);
            if (!diffFile.good())
            {
                FatalErrorInFunction
                    << "Cannot open file " << foundName
                    << exit(FatalError, 1);
            }

            // Read but preserve headers
            diffDict.read(diffFile, true);
            optDiff = true;
        }
    }

    word scopedName;  // Actually fileName, since it can contain '/' scoping
    if (args.optionReadIfPresent("entry", scopedName))
    {
        upgradeScope(scopedName);

        string newValue;
        if
        (
            args.optionReadIfPresent("set", newValue)
         || args.optionReadIfPresent("add", newValue)
        )
        {
            const bool overwrite = args.optionFound("set");

            // Dictionary name and keyword
            const dictAndKeyword dAk(scopedName);

            // The context for the action
            const dictionary& d(lookupScopedDict(dict, dAk.dict()));

            // Create a new entry
            IStringStream str(string(dAk.key()) + ' ' + newValue + ';');
            entry* ePtr(entry::New(str).ptr());

            if (overwrite)
            {
                const_cast<dictionary&>(d).set(ePtr);
            }
            else
            {
                const_cast<dictionary&>(d).add(ePtr, false);
            }
            changed = true;

            // Print the changed entry
            const auto finder = dict.csearchScoped
            (
                scopedName,
                false,
                true  // Support wildcards
            );

            if (finder.found())
            {
                Info<< finder.ref();
            }
        }
        else if (args.optionFound("remove"))
        {
            // Dictionary name and keyword
            const dictAndKeyword dAk(scopedName);

            // The context for the action
            const dictionary& d(lookupScopedDict(dict, dAk.dict()));

            const_cast<dictionary&>(d).remove(dAk.key());
            changed = true;
        }
        else
        {
            // Optionally remove a second dictionary
            if (optDiff)
            {
                // Dictionary name and keyword
                const dictAndKeyword dAk(scopedName);

                const dictionary& d1(lookupScopedDict(dict, dAk.dict()));
                const dictionary& d2(lookupScopedDict(diffDict, dAk.dict()));

                const entry* e1Ptr = d1.lookupEntryPtr(dAk.key(), false, true);
                const entry* e2Ptr = d2.lookupEntryPtr(dAk.key(), false, true);

                if (e1Ptr && e2Ptr)
                {
                    if (*e1Ptr == *e2Ptr)
                    {
                        const_cast<dictionary&>(d1).remove(dAk.key());
                    }
                    else if (e1Ptr->isDict() && e2Ptr->isDict())
                    {
                        removeDict
                        (
                            const_cast<dictionary&>(e1Ptr->dict()),
                            e2Ptr->dict()
                        );
                    }
                }
            }

            const auto finder = dict.csearchScoped
            (
                scopedName,
                false,
                true  // Support wildcards
            );

            if (!finder.found())
            {
                FatalIOErrorInFunction(dictFile())
                    << "Cannot find entry " << scopedName
                    << exit(FatalIOError, 2);
            }
            else if (args.optionFound("keywords"))
            {
                for (const entry& e : finder.dict())
                {
                    Info<< e.keyword() << endl;
                }
            }
            else if (args.optionFound("value"))
            {
                if (finder.isDict())
                {
                    Info<< finder.dict();
                }
                else if (finder.ref().isStream())
                {
                    const tokenList& tokens = finder.ref().stream();
                    forAll(tokens, i)
                    {
                        Info<< tokens[i];
                        if (i < tokens.size() - 1)
                        {
                            Info<< token::SPACE;
                        }
                    }
                    Info<< endl;
                }
            }
            else
            {
                Info<< finder.ref();
            }
        }
    }
    else if (args.optionFound("keywords"))
    {
        for (const entry& e : dict)
        {
            Info<< e.keyword() << endl;
        }
    }
    else if (optDiff)
    {
        removeDict(dict, diffDict);
        dict.write(Info, false);
    }
    else
    {
        dict.write(Info, false);
    }

    if (changed)
    {
        dictFile.clear();
        OFstream os(dictFileName);
        IOobject::writeBanner(os);
        dict.write(os, false);
        IOobject::writeEndDivider(os);
    }

    return 0;
}


// ************************************************************************* //
