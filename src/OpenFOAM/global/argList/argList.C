/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2022 OpenCFD Ltd.
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

#include "argList.H"
#include "OSspecific.H"
#include "Switch.H"
#include "clock.H"
#include "dictionary.H"
#include "IOobject.H"
#include "JobInfo.H"
#include "labelList.H"
#include "IOobject.H"
#include "dynamicCode.H"
#include "simpleObjectRegistry.H"
#include "sigFpe.H"
#include "sigInt.H"
#include "sigQuit.H"
#include "sigSegv.H"
#include "foamVersion.H"
#include "stringOps.H"
#include "CStringList.H"
#include "stringListOps.H"
#include "fileOperation.H"
#include "fileOperationInitialise.H"

#include <cctype>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

bool Foam::argList::argsMandatory_ = true;
bool Foam::argList::checkProcessorDirectories_ = true;

Foam::SLList<Foam::string>    Foam::argList::validArgs;
Foam::HashSet<Foam::string>   Foam::argList::advancedOptions;
Foam::HashTable<Foam::string> Foam::argList::validOptions;
Foam::HashTable<Foam::string> Foam::argList::validParOptions;
Foam::HashTable<Foam::string, Foam::label, Foam::Hash<Foam::label>>
    Foam::argList::argUsage;
Foam::HashTable<Foam::string> Foam::argList::optionUsage;
Foam::HashTable<std::pair<Foam::word,int>> Foam::argList::validOptionsCompat;
Foam::HashTable<std::pair<bool,int>> Foam::argList::ignoreOptionsCompat;
Foam::SLList<Foam::string>    Foam::argList::notes;
Foam::word Foam::argList::postProcessOptionName("postProcess");

std::string::size_type Foam::argList::usageMin = 20;
std::string::size_type Foam::argList::usageMax = 80;

Foam::argList::initValidTables::initValidTables()
{
    argList::addOption
    (
        "case",
        "dir",
        "Specify case directory to use (instead of cwd)"
    );
    argList::addOption
    (
        "lib",
        "name",
        "Additional library or library list to load"
        " (can be used multiple times)",
        true  // advanced option
    );

    argList::addOption
    (
        "debug-switch",
        "name=val",
        "Specify the value of a registered debug switch."
        " Default is 1 if the value is omitted."
        " (Can be used multiple times)",
        true  // advanced option
    );

    argList::addOption
    (
        "info-switch",
        "name=val",
        "Specify the value of a registered info switch."
        " Default is 1 if the value is omitted."
        " (Can be used multiple times)",
        true  // advanced option
    );

    argList::addOption
    (
        "opt-switch",
        "name=val",
        "Specify the value of a registered optimisation switch."
        " Default is 1 if the value is omitted."
        " (Can be used multiple times)",
        true  // advanced option
    );

    argList::addBoolOption("parallel", "Run in parallel");
    validParOptions.set("parallel", "");
    argList::addOption
    (
        "roots",
        "(dir1 .. dirN)",
        "Subprocess root directories for distributed running",
        true  // advanced option
    );
    validParOptions.set
    (
        "roots",
        "(dir1 .. dirN)"
    );

    argList::addOption
    (
        "decomposeParDict",
        "file",
        "Use specified file for decomposePar dictionary"
    );
    argList::addOption
    (
        "hostRoots",
        "((host1 dir1) .. (hostN dirN))",
        "Per-subprocess root directories for distributed running."
        " The host specification can be a regex.",
        true  // advanced option
    );
    validParOptions.set
    (
        "hostRoots",
        "((host1 dir1) .. (hostN dirN))"
    );

    argList::addBoolOption
    (
        "noFunctionObjects",
        "Do not execute function objects",
        true  // advanced option
    );

    argList::addOption
    (
        "fileHandler",
        "handler",
        "Override the file handler type",
        true  // advanced option
    );

    argList::addOption
    (
        "world",
        "name",
        "Name of the local world for parallel communication",
        true  // advanced option
    );
    validParOptions.set
    (
        "world",
        "name"
    );


    // Some standard option aliases (with or without version warnings)
//     argList::addOptionCompat
//     (
//         "noFunctionObjects", {"no-function-objects", 0}
//     );

    Pstream::addValidParOptions(validParOptions);
}

Foam::argList::initValidTables dummyInitValidTables;


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Counted per machine name
// Does not include any sorting since we wish to know the ordering according to
// mpi rank.
//
// Always include the master too.
// This provides a better overview of the subscription
static void printHostsSubscription(const UList<string>& hostProcs)
{
    Info<< "Hosts  :\n(" << nl;

    std::string prev = Foam::hostName();
    int count = 1;

    for (const auto& str : hostProcs)
    {
        std::string curr(str.substr(0, str.rfind('.')));

        if (prev != curr)
        {
            if (count)
            {
                // Finish previous
                Info<< "    (" << prev.c_str() << ' ' << count << ')' << nl;
                count = 0;
            }

            prev = std::move(curr);
        }
        ++count;
    }

    if (count)
    {
        // Finished last one
        Info<< "    (" << prev.c_str() << ' ' << count << ')' << nl;
    }

    Info<< ')' << nl;
}

} // End namespace Foam


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

void Foam::argList::checkITstream(const ITstream& is, const label index)
{
    const label remaining = is.nRemainingTokens();

    if (remaining)
    {
        std::cerr
            << nl
            << "--> FOAM WARNING:" << nl
            << "Argument " << index << " has "
            << remaining << " excess tokens" << nl << nl;
    }
    else if (!is.size())
    {
        std::cerr
            << nl
            << "--> FOAM WARNING:" << nl
            << "Argument " << index << " had no tokens" << nl << nl;
    }
}


void Foam::argList::checkITstream(const ITstream& is, const word& optName)
{
    const label remaining = is.nRemainingTokens();

    if (remaining)
    {
        std::cerr
            << nl
            << "--> FOAM WARNING:" << nl
            << "Option -" << optName << " has "
            << remaining << " excess tokens" << nl << nl;
    }
    else if (!is.size())
    {
        std::cerr
            << nl
            << "--> FOAM WARNING:" << nl
            << "Option -" << optName << " had no tokens" << nl << nl;
    }
}


void Foam::argList::raiseBadInput(const word& optName) const
{
    // Can use FatalError
    // predicate checks are not used at the earliest stages
    FatalErrorIn(executable())
        << "Option -" << optName << " with invalid input" << nl
        << exit(FatalError);
}


void Foam::argList::addArgument
(
    const string& argName,
    const string& usage
)
{
    validArgs.append(argName);

    // The first program argument starts at 1 - obtain index after the append

    const label index = validArgs.size();

    if (usage.empty())
    {
        argUsage.erase(index);
    }
    else
    {
        argUsage.set(index, usage);
    }
}


void Foam::argList::addBoolOption
(
    const word& optName,
    const string& usage,
    bool advanced
)
{
    argList::addOption(optName, "", usage, advanced);
}


void Foam::argList::addOption
(
    const word& optName,
    const string& param,
    const string& usage,
    bool advanced
)
{
    validOptions.set(optName, param);
    if (!usage.empty())
    {
        optionUsage.set(optName, usage);
    }
    if (advanced)
    {
        advancedOptions.set(optName);
    }
}


void Foam::argList::setAdvanced(const word& optName, bool advanced)
{
    if (advanced && validOptions.found(optName))
    {
        advancedOptions.set(optName);
    }
    else
    {
        advancedOptions.erase(optName);
    }
}


void Foam::argList::addOptionCompat
(
    const word& optName,
    std::pair<const char*,int> compat
)
{
    validOptionsCompat.insert
    (
        compat.first,
        std::pair<word,int>(optName, compat.second)
    );
}


void Foam::argList::ignoreOptionCompat
(
    std::pair<const char*,int> compat,
    bool expectArg
)
{
    ignoreOptionsCompat.insert
    (
        compat.first,
        std::pair<bool,int>(expectArg, compat.second)
    );
}


void Foam::argList::addUsage
(
    const word& optName,
    const string& usage
)
{
    if (usage.empty())
    {
        optionUsage.erase(optName);
    }
    else
    {
        optionUsage.set(optName, usage);
    }
}


void Foam::argList::addNote(const string& note)
{
    if (!note.empty())
    {
        notes.append(note);
    }
}


void Foam::argList::removeOption(const word& optName)
{
    validOptions.erase(optName);
    optionUsage.erase(optName);
    advancedOptions.erase(optName);
}


void Foam::argList::noMandatoryArgs()
{
    argsMandatory_ = false;
}


bool Foam::argList::argsMandatory()
{
    return argsMandatory_;
}


void Foam::argList::noBanner()
{
    ::Foam::infoDetailLevel = 0;
}


bool Foam::argList::bannerEnabled()
{
    return (::Foam::infoDetailLevel > 0);
}


void Foam::argList::addDryRunOption
(
    const string& usage,
    bool advanced
)
{
    argList::addBoolOption("dry-run", usage, advanced);
}


void Foam::argList::addVerboseOption
(
    const string& usage,
    bool advanced
)
{
    argList::addBoolOption("verbose", usage, advanced);
}


void Foam::argList::noFunctionObjects(bool addWithOption)
{
    removeOption("noFunctionObjects");

    // Ignore this bool option without warning
    // - cannot tie to any particular version anyhow
    ignoreOptionCompat({"noFunctionObjects", 0}, false);

    if (addWithOption)
    {
        argList::addBoolOption
        (
            "withFunctionObjects",
            "Execute functionObjects",
            true  // advanced option
        );
    }
}


void Foam::argList::noJobInfo()
{
    JobInfo::disable();
}


void Foam::argList::noLibs()
{
    argList::addBoolOption
    (
        "no-libs",
        "Disable use of the controlDict libs entry",
        true  // advanced option
    );
}


void Foam::argList::noParallel()
{
    removeOption("parallel");
    removeOption("roots");
    removeOption("decomposeParDict");
    removeOption("hostRoots");
    removeOption("world");
    validParOptions.clear();
}


void Foam::argList::noCheckProcessorDirectories()
{
    checkProcessorDirectories_ = false;
}


bool Foam::argList::postProcess(int argc, char *argv[])
{
    for (int i=1; i<argc; ++i)
    {
        if (argv[i] == '-' + postProcessOptionName)
        {
            return true;
        }
    }

    return false;
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::word Foam::argList::envExecutable()
{
    return Foam::getEnv("FOAM_EXECUTABLE");
}


Foam::fileName Foam::argList::envGlobalPath()
{
    return Foam::getEnv("FOAM_CASE");
}


Foam::fileName Foam::argList::envRelativePath
(
    const fileName& input,
    const bool caseTag
)
{
    if (input.isAbsolute())
    {
        return input.relative(envGlobalPath(), caseTag);
    }

    return input;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::word Foam::argList::optionCompat(const word& optName)
{
    // NB: optName includes the leading '-' so that the return value
    // can be used directly

    if (!validOptionsCompat.empty())
    {
        const auto fnd = validOptionsCompat.cfind(optName.substr(1));

        if (fnd.found())
        {
            const auto& alt = fnd.val();

            // No error::master() guard - only called on master anyhow
            if (error::warnAboutAge(alt.second))
            {
                std::cerr
                    << "--> FOAM IOWarning :" << nl
                    << "    Found [v" << alt.second << "] '"
                    << optName << "' instead of '-"
                    << alt.first << "' option"
                    << nl
                    << std::endl;

                error::warnAboutAge("option", alt.second);
            }

            return "-" + alt.first;
        }
    }

    // Nothing found - pass through the original input
    return optName;
}


int Foam::argList::optionIgnore(const word& optName)
{
    // NB: optName is without the leading '-'

    if (!ignoreOptionsCompat.empty())
    {
        const auto fnd = ignoreOptionsCompat.cfind(optName);

        if (fnd.found())
        {
            const auto& alt = fnd.val();

            // Number to skip (including the option itself)
            // '-option ARG' or '-option'
            const int nskip = (alt.first ? 2 : 1);

            // No error::master() guard - only called on master anyhow
            if (error::warnAboutAge(alt.second))
            {
                std::cerr
                    << "--> FOAM IOWarning :" << nl
                    << "    Ignoring [v" << alt.second << "] '-"
                    << optName << (nskip > 1 ? " ARG" : "")
                    << "' option"
                    << nl
                    << std::endl;

                error::warnAboutAge("option", alt.second);
            }

            return nskip;
        }
    }

    return 0; // Do not skip
}


bool Foam::argList::regroupArgv(int& argc, char**& argv)
{
    int nArgs = 1;
    int ignore = 0;
    unsigned depth = 0;
    string group;  // For grouping ( ... ) arguments

    // Note: we rewrite directly into args_
    // and use a second pass to sort out args/options

    args_[0] = fileName(argv[0]);
    for (int argi = 1; argi < argc; ++argi)
    {
        if (strcmp(argv[argi], "(") == 0)
        {
            ++depth;
            group += '(';
        }
        else if (strcmp(argv[argi], ")") == 0)
        {
            if (depth)
            {
                --depth;
                group += ')';
                if (!depth)
                {
                    args_[nArgs++] = group;
                    group.clear();
                }
            }
            else
            {
                args_[nArgs++] = argv[argi];
            }
        }
        else if (depth)
        {
            // Quote each string element
            group += '"';
            group += argv[argi];
            group += '"';
        }
        else if (argv[argi][0] == '-')
        {
            // Appears to be an option
            const char *optName = &argv[argi][1];

            if (validOptions.found(optName))
            {
                // Known option name
                args_[nArgs++] = argv[argi];
            }
            else if ((ignore = optionIgnore(optName)) > 0)
            {
                // Option to be ignored (with/without an argument)
                if (ignore > 1)
                {
                    ++argi;
                }
            }
            else
            {
                // Try alias for the option name
                args_[nArgs++] = optionCompat(argv[argi]);
            }
        }
        else
        {
            args_[nArgs++] = argv[argi];
        }
    }

    if (group.size())
    {
        // Group(s) not closed, but flush anything still pending
        args_[nArgs++] = group;
    }

    args_.resize(nArgs);

    std::string::size_type len = (nArgs-1); // Spaces between args
    for (const auto& s : args_)
    {
        len += s.length();
    }

    // Length needed for regrouped command-line
    commandLine_.reserve(len);

    return nArgs < argc;
}


void Foam::argList::setCasePaths()
{
    fileName caseDir;

    const auto optIter = options_.cfind("case");  // [-case dir] specified?

    if (optIter.found())
    {
        caseDir = fileName::validate(optIter.val());  // includes 'clean'

        if (caseDir.empty() || caseDir == ".")
        {
            // Treat "", "." and "./" as if -case was not specified
            caseDir = cwd();
            options_.erase("case");
        }
        else
        {
            caseDir.expand();
            caseDir.toAbsolute();
        }
    }
    else
    {
        // Nothing specified, use the current dir
        caseDir = cwd();
    }

    // The caseDir is a cleaned, absolute path

    rootPath_   = caseDir.path();
    globalCase_ = caseDir.name();
    case_       = globalCase_;  // The (processor) local case name

    // OPENFOAM API
    setEnv("FOAM_API", std::to_string(foamVersion::api), true);

    // Global case (directory) and case-name as environment variables
    setEnv("FOAM_CASE", caseDir, true);
    setEnv("FOAM_CASENAME", globalCase_, true);

    // Executable name, unless already present in the environment
    setEnv("FOAM_EXECUTABLE", executable_, false);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::argList::argList
(
    int& argc,
    char**& argv,
    bool checkArgs,
    bool checkOpts,
    bool initialise
)
:
    args_(argc),
    options_(argc),
    libs_()
{
    // Check for -fileHandler, which requires an argument.
    word handlerType;
    for (int argi = argc-2; argi > 0; --argi)
    {
        if (argv[argi][0] == '-')
        {
            const char *optName = &argv[argi][1];

            if (strcmp(optName, "fileHandler") == 0)
            {
                handlerType = argv[argi+1];
                break;
            }
        }
    }
    if (handlerType.empty())
    {
        handlerType = Foam::getEnv("FOAM_FILEHANDLER");
        if (handlerType.empty())
        {
            handlerType = fileOperation::defaultFileHandler;
        }
    }

    // Detect any parallel options
    const bool needsThread = fileOperations::fileOperationInitialise::New
    (
        handlerType,
        argc,
        argv
    )().needsThreading();


    // Check if this run is a parallel run by searching for any parallel option
    // If found call runPar which might filter argv
    for (int argi = 1; argi < argc; ++argi)
    {
        if (argv[argi][0] == '-')
        {
            const char *optName = &argv[argi][1];

            if (validParOptions.found(optName))
            {
                runControl_.runPar(argc, argv, needsThread);
                break;
            }
        }
    }

    // Convert argv -> args_ and capture ( ... ) lists
    regroupArgv(argc, argv);
    commandLine_ += args_[0];

    // Set executable name immediately - useful when emitting errors.
    executable_ = fileName(args_[0]).name();

    // Count -dry-run and -verbose switches
    int numDryRun = 0, numVerbose = 0;

    // Check arguments and options, argv[0] was already handled
    int nArgs = 1;
    for (int argi = 1; argi < args_.size(); ++argi)
    {
        commandLine_ += ' ';
        commandLine_ += args_[argi];

        if (args_[argi][0] == '-')
        {
            const char *optName = &args_[argi][1];

            if (!*optName)
            {
                Warning
                    << "Ignoring lone '-' on the command-line" << endl;
                continue;
            }

            // Option known and expects an argument?
            // - use Switch for a tri-state
            //   True  : known option, expects a parameter
            //   False : known option, no parameter
            //   bad() : unknown option

            Switch wantArg(Switch::INVALID);
            auto optIter = validOptions.cfind(optName);
            if
            (
                optIter.found()
             || (optIter = validParOptions.cfind(optName)).found()
            )
            {
                wantArg = !optIter.val().empty();
            }


            if (wantArg)
            {
                // Known option and expects a parameter
                // - get it or emit a FatalError.

                ++argi;
                if (argi >= args_.size())
                {
                    foamVersion::printBuildInfo(Info.stdStream(), false);

                    Info<< nl
                        <<"Error: option '-" << optName
                        << "' requires an argument" << nl << nl
                        << "See '" << executable_ << " -help' for usage"
                        << nl << nl;

                    Pstream::exit(1);  // works for serial and parallel
                }

                commandLine_ += ' ';
                commandLine_ += args_[argi];

                //
                // Special handling of these options
                //

                if (strcmp(optName, "lib") == 0)
                {
                    // The '-lib' option:
                    // Append name(s) to libs for later opening
                    libs().append(this->getList<fileName>(argi));
                }
                else if (strcmp(optName, "debug-switch") == 0)
                {
                    // The '-debug-switch' option:
                    // change registered debug switch
                    DetailInfo << "debug-switch ";
                    debug::debugObjects()
                        .setNamedValue(args_[argi], 1, true);
                }
                else if (strcmp(optName, "info-switch") == 0)
                {
                    // The '-info-switch' option:
                    // change registered info switch
                    DetailInfo << "info-switch ";
                    debug::infoObjects()
                        .setNamedValue(args_[argi], 1, true);
                }
                else if (strcmp(optName, "opt-switch") == 0)
                {
                    // The '-opt-switch' option:
                    // change registered optimisation switch
                    DetailInfo << "opt-switch ";
                    debug::optimisationObjects()
                        .setNamedValue(args_[argi], 1, true);
                }
                else
                {
                    // Regular option (with a parameter):
                    // Duplicates handled by using the last -option specified
                    options_.set(optName, args_[argi]);
                }
            }
            else
            {
                // All other options (including unknown ones) are simply
                // registered as existing.

                options_.insert(optName, "");

                // Special increment handling for some known flags
                if (wantArg.good())
                {
                    if (strcmp(optName, "dry-run") == 0)
                    {
                        ++numDryRun;
                    }
                    else if (strcmp(optName, "verbose") == 0)
                    {
                        ++numVerbose;
                    }
                }
            }
        }
        else
        {
            if (nArgs != argi)
            {
                args_[nArgs] = args_[argi];
            }
            ++nArgs;
        }
    }

    // Commit number of -dry-run and -verbose flag occurrences
    runControl_.dryRun(numDryRun);
    runControl_.verbose(numVerbose);

    args_.resize(nArgs);

    parse(checkArgs, checkOpts, initialise);
}


Foam::argList::argList
(
    const argList& args,
    const HashTable<string>& options,
    bool checkArgs,
    bool checkOpts,
    bool initialise
)
:
    runControl_(args.runControl_),
    args_(args.args_),
    options_(options),
    libs_(),
    executable_(args.executable_),
    rootPath_(args.rootPath_),
    globalCase_(args.globalCase_),
    case_(args.case_),
    commandLine_(args.commandLine_)
{
    parse(checkArgs, checkOpts, initialise);
}


void Foam::argList::parse
(
    bool checkArgs,
    bool checkOpts,
    bool initialise
)
{
    // Help/documentation options:
    //   -doc         Display documentation in browser
    //   -doc-source  Display source code in browser
    //   -help        Display short help and exit
    //   -help-compat Display compatibility options
    //   -help-full   Display full help and exit
    {
        bool quickExit = false;

        // Display either application or source documentation, not both
        if (options_.found("doc"))
        {
            displayDoc(false);
            quickExit = true;
        }
        else if (options_.found("doc-source"))
        {
            displayDoc(true);
            quickExit = true;
        }

        // Display either short or full help, not both
        if (options_.found("help-full"))
        {
            printUsage(true);
            quickExit = true;
        }
        else if (options_.found("help-notes"))
        {
            printNotes();
            Info<< nl;
            quickExit = true;
        }
        else if (options_.found("help"))
        {
            printUsage(false);
            quickExit = true;
        }
        else if (options_.found("help-man"))
        {
            printMan();
            quickExit = true;
        }

        // Allow independent display of compatibility information
        if (options_.found("help-compat"))
        {
            printCompat();
            quickExit = true;
        }

        if (quickExit)
        {
            std::exit(0);
        }
    }

    // Print the collected error messages and exit if check fails
    if (!check(checkArgs, checkOpts))
    {
        foamVersion::printBuildInfo(Info.stdStream(), false);
        FatalError.write(Info, false);

        Pstream::exit(1);  // works for serial and parallel
    }

    if (initialise)
    {
        const string dateString = clock::date();
        const string timeString = clock::clockTime();

        // Print the banner once only for parallel runs
        if (Pstream::master() && bannerEnabled())
        {
            IOobject::writeBanner(Info, true)
                << "Build  : ";

            if (foamVersion::build.size())
            {
                Info<< foamVersion::build.c_str() << ' ';
            }

            Info<< "OPENFOAM=" << foamVersion::api;

            if (foamVersion::patched())
            {
                // Patch-level, when defined
                Info<< " patch=" << foamVersion::patch.c_str();
            }

            Info<< " version=" << foamVersion::version.c_str();

            Info<< nl
                << "Arch   : " << foamVersion::buildArch << nl
                << "Exec   : " << commandLine_.c_str() << nl
                << "Date   : " << dateString.c_str() << nl
                << "Time   : " << timeString.c_str() << nl
                << "Host   : " << Foam::hostName().c_str() << nl
                << "PID    : " << pid() << nl;
        }

        jobInfo.add("startDate", dateString);
        jobInfo.add("startTime", timeString);
        jobInfo.add("userName", userName());

        jobInfo.add("foamApi", foamVersion::api);
        jobInfo.add("foamVersion", word(foamVersion::version));

        // Add build information - only use the first word
        {
            std::string build(foamVersion::build);
            const auto space = build.find(' ');
            if (space != std::string::npos)
            {
                build.resize(space);
            }
            jobInfo.add("foamBuild", build);
        }

        jobInfo.add("code", executable_);
        jobInfo.add("argList", commandLine_);
        jobInfo.add("currentDir", cwd());
        jobInfo.add("PPID", ppid());
        jobInfo.add("PGID", pgid());

        // Load additional libraries (verbosity according to banner setting)
        libs().open(bannerEnabled());
    }


    // Set fileHandler. In increasing order of priority:
    // 1. default = uncollated
    // 2. env variable "FOAM_FILEHANDLER"
    // 3. etc/controlDict optimisationSwitches 'fileHandler'
    // 4. system/controlDict 'fileHandler' (not handled here; done in TimeIO.C)
    // 5. '-fileHandler' commmand-line option

    {
        word handlerType
        (
            options_.lookup("fileHandler", Foam::getEnv("FOAM_FILEHANDLER"))
        );

        if (handlerType.empty())
        {
            handlerType = fileOperation::defaultFileHandler;
        }

        Foam::fileHandler(fileOperation::New(handlerType, bannerEnabled()));
    }


    stringList hostMachine;
    stringList hostProcs;
    const int writeHostsSwitch = Foam::debug::infoSwitch("writeHosts", 1);

    // Collect machine/pid, and check that the build is identical
    if (runControl_.parRun())
    {
        if (Pstream::master())
        {
            hostMachine.resize(Pstream::nProcs()-1);
            hostProcs.resize(Pstream::nProcs()-1);
            string procBuild;
            label procPid;
            int proci = 0;
            for (const int subproci : Pstream::subProcs())
            {
                IPstream fromSubproc(Pstream::commsTypes::scheduled, subproci);

                fromSubproc >> procBuild >> hostMachine[proci] >> procPid;

                hostProcs[proci] = hostMachine[proci] + "." + name(procPid);
                ++proci;

                // Verify that all processors are running the same build
                if (procBuild != foamVersion::build)
                {
                    FatalErrorIn(executable())
                        << "Running build version " << foamVersion::build
                        << " but proc " << subproci << " is running "
                        << procBuild << nl
                        << exit(FatalError);
                }
            }
        }
        else
        {
            OPstream toMaster
            (
                Pstream::commsTypes::scheduled,
                Pstream::masterNo()
            );
            toMaster << foamVersion::build << Foam::hostName() << Foam::pid();
        }
    }


    // Case is a single processor run unless it is running parallel
    int nProcs = 1;

    // Roots if running distributed
    fileNameList roots;

    // If this actually is a parallel run
    if (runControl_.parRun())
    {
        // For the master
        if (Pstream::master())
        {
            // Establish rootPath_/globalCase_/case_ for master
            setCasePaths();

            // The system/decomposeParDict (or equivalent)
            fileName source;

            if (this->readIfPresent("decomposeParDict", source))
            {
                bool adjustOpt = false;

                if (isDir(source))
                {
                    source /= "decomposeParDict";
                    adjustOpt = true;
                }

                // Case-relative if not absolute and not "./" etc
                if (!source.isAbsolute() && !source.starts_with('.'))
                {
                    source = rootPath_/globalCase_/source;
                    adjustOpt = true;
                }

                // Could also check for absolute path, but shouldn't be needed
                if (adjustOpt)
                {
                    source.clean();  // Remove unneeded ".."
                    options_.set("decomposeParDict", source);
                }
            }

            // If running distributed (different roots for different procs)
            label dictNProcs = -1;
            if (this->readListIfPresent("roots", roots))
            {
                source = "-roots";
                runControl_.distributed(true);

                if (roots.empty())
                {
                    FatalErrorInFunction
                        << "The -roots option must contain values"
                        << exit(FatalError);
                }
                if (roots.size() > 1)
                {
                    dictNProcs = roots.size()+1;
                }
            }
            else if (options_.found("hostRoots"))
            {
                source = "-hostRoots";
                runControl_.distributed(true);

                ITstream is(this->lookup("hostRoots"));

                List<Tuple2<wordRe, fileName>> hostRoots(is);
                checkITstream(is, "hostRoots");

                if (hostRoots.empty())
                {
                    FatalErrorInFunction
                        << "The -hostRoots option must contain values"
                        << exit(FatalError);
                }

                // Match machine names to roots
                roots.resize(Pstream::nProcs()-1, fileName::null);
                for (const auto& hostRoot : hostRoots)
                {
                    labelList matched
                    (
                        findMatchingStrings(hostRoot.first(), hostMachine)
                    );
                    for (const label matchi : matched)
                    {
                        if (!roots[matchi].empty())
                        {
                            FatalErrorInFunction
                                << "Multiple matching roots for "
                                << hostMachine[matchi] << " in "
                                << hostRoots << nl
                                << exit(FatalError);
                        }

                        roots[matchi] = hostRoot.second();
                    }
                }

                // Check
                forAll(roots, hosti)
                {
                    if (roots[hosti].empty())
                    {
                        FatalErrorInFunction
                            << "No matching roots for "
                            << hostMachine[hosti] << " in "
                            << hostRoots << nl
                            << exit(FatalError);
                    }
                }

                if (roots.size() > 1)
                {
                    dictNProcs = roots.size()+1;
                }
            }
            else if (checkProcessorDirectories_ && Pstream::nProcs() > 1)
            {
                // Check values from decomposeParDict

                const bool useDefault = source.empty();
                if (useDefault)
                {
                    source = rootPath_/globalCase_/"system"/"decomposeParDict";
                }

                // Disable any parallel comms happening inside the fileHandler
                // since we are on master. This can happen e.g. inside
                // the masterUncollated/collated handler. Note that we
                // also have to protect the actual dictionary parsing since
                // it might trigger file access (e.g. #include, #codeStream)
                const bool oldParRun = Pstream::parRun(false);

                autoPtr<ISstream> dictStream
                (
                    fileHandler().NewIFstream(source)
                );

                if (dictStream && dictStream->good())
                {
                    dictionary decompDict(*dictStream);
                    bool nDomainsMandatory = false;

                    if (decompDict.getOrDefault("distributed", false))
                    {
                        nDomainsMandatory = true;
                        runControl_.distributed(true);
                        decompDict.readEntry("roots", roots);

                        if (roots.empty())
                        {
                            DetailInfo
                                << "WARNING: running distributed"
                                << " but did not specify roots!" << nl;
                        }
                    }

                    // Get numberOfSubdomains if it exists.
                    // - mandatory when running with distributed roots
                    decompDict.readEntry
                    (
                        "numberOfSubdomains",
                        dictNProcs,
                        keyType::LITERAL,
                        nDomainsMandatory
                    );
                }
                else
                {
                    if (useDefault)
                    {
                        // Optional if using default location
                        DetailInfo
                            << "WARNING: running without decomposeParDict "
                            << this->relativePath(source) << nl;
                    }
                    else
                    {
                        // Mandatory if specified as -decomposeParDict
                        FatalError
                            << "Cannot read decomposeParDict: "
                            << this->relativePath(source) << nl
                            << exit(FatalError);
                    }
                }

                Pstream::parRun(oldParRun);  // Restore parallel state

                if (Pstream::nProcs() == 1)
                {
                    Warning
                        << "Running parallel on single processor. This only"
                        << " makes sense for multi-world simulation" << endl;
                    dictNProcs = 1;
                }
            }

            // Convenience:
            // when a single root is specified, use it for all processes
            if (roots.size() == 1)
            {
                const fileName rootName(roots[0]);
                roots.resize(Pstream::nProcs()-1, rootName);

                // Adjust dictNProcs for command-line '-roots' option
                if (dictNProcs <= 0)
                {
                    dictNProcs = roots.size()+1;
                }
            }


            // Check number of processors.
            // nProcs     => number of actual procs
            // dictNProcs => number of procs specified in decompositionDict
            // nProcDirs  => number of processor directories
            //               (n/a when running distributed)
            //
            // - normal running : nProcs = dictNProcs = nProcDirs
            // - decomposition to more  processors : nProcs = dictNProcs
            // - decomposition to fewer processors : nProcs = nProcDirs
            if
            (
                checkProcessorDirectories_
             && Pstream::nProcs() > 1
             && dictNProcs > Pstream::nProcs()
            )
            {
                FatalError
                    << this->relativePath(source)
                    << " specifies " << dictNProcs
                    << " processors but job was started with "
                    << Pstream::nProcs() << " processors."
                    << exit(FatalError);
            }

            // Distributed data
            if (roots.size())
            {
                if (roots.size() != Pstream::nProcs()-1)
                {
                    FatalError
                        << "number of entries in roots "
                        << roots.size()
                        << " is not equal to the number of sub-processes "
                        << Pstream::nProcs()-1
                        << exit(FatalError);
                }

                for (fileName& dir : roots)
                {
                    dir.expand();
                }

                // Distribute the master's argument list (with new root)
                const bool hadCaseOpt = options_.found("case");
                for (const int subproci : Pstream::subProcs())
                {
                    options_.set("case", roots[subproci-1]/globalCase_);

                    OPstream toProc(Pstream::commsTypes::scheduled, subproci);

                    toProc
                        << args_ << options_
                        << runControl_.distributed()
                        << label(runControl_.dryRun())
                        << label(runControl_.verbose());
                }
                options_.erase("case");

                // Restore [-case dir]
                if (hadCaseOpt)
                {
                    options_.set("case", rootPath_/globalCase_);
                }
            }
            else
            {
                // Possibly going to fewer processors.
                // Check if all procDirs are there.
                if
                (
                    checkProcessorDirectories_
                 && Pstream::nProcs() > 1
                 && dictNProcs >= 1
                 && dictNProcs < Pstream::nProcs()
                )
                {
                    label nProcDirs = 0;
                    while
                    (
                        isDir
                        (
                            rootPath_/globalCase_
                          / ("processor" + Foam::name(++nProcDirs))
                        )
                    )
                    {}

                    if (nProcDirs != Pstream::nProcs())
                    {
                        FatalError
                            << "number of processor directories = "
                            << nProcDirs
                            << " is not equal to the number of processors = "
                            << Pstream::nProcs()
                            << exit(FatalError);
                    }
                }

                // Distribute the master's argument list (unaltered)
                for (const int subproci : Pstream::subProcs())
                {
                    OPstream toProc(Pstream::commsTypes::scheduled, subproci);

                    toProc
                        << args_ << options_
                        << runControl_.distributed()
                        << label(runControl_.dryRun())
                        << label(runControl_.verbose());
                }
            }
        }
        else
        {
            // Collect the master's argument list
            bool isDistributed;
            label numDryRun, numVerbose;

            IPstream fromMaster
            (
                Pstream::commsTypes::scheduled,
                Pstream::masterNo()
            );

            fromMaster
                >> args_ >> options_
                >> isDistributed
                >> numDryRun >> numVerbose;

            runControl_.distributed(isDistributed);
            runControl_.dryRun(numDryRun);
            runControl_.verbose(numVerbose);

            // Establish rootPath_/globalCase_/case_ for sub-process
            setCasePaths();
        }

        nProcs = Pstream::nProcs();
        if (Pstream::nProcs() > 1)
        {
            case_ = globalCase_/("processor" + Foam::name(Pstream::myProcNo()));
        }
        else
        {
            case_ = globalCase_;
        }
    }
    else
    {
        // Establish rootPath_/globalCase_/case_
        setCasePaths();
        case_ = globalCase_;   // Redundant, but extra safety?
    }

    // If needed, adjust fileHandler for distributed roots
    if (runControl_.distributed())
    {
        if (fileOperation::fileHandlerPtr_)
        {
            fileOperation::fileHandlerPtr_->distributed(true);
        }
    }

    // Keep/discard sub-process host/root information for reporting:
    if (Pstream::master() && runControl_.parRun())
    {
        if (!writeHostsSwitch)
        {
            // Clear here to ensures it doesn't show in the jobInfo
            hostProcs.clear();
        }
        if (!debug::infoSwitch("writeRoots", 1))
        {
            roots.clear();
        }
    }

    if (Pstream::master() && bannerEnabled())
    {
        Info<< "Case   : " << (rootPath_/globalCase_).c_str() << nl
            << "nProcs : " << nProcs << nl;

        if (runControl_.parRun())
        {
            if (hostProcs.size())
            {
                if (writeHostsSwitch == 1)
                {
                    // Compact output (see etc/controlDict)
                    printHostsSubscription(hostProcs);
                }
                else if (writeHostsSwitch)
                {
                    // Full output of "host.pid"
                    Info<< "Hosts  :\n(" << nl;

                    // Include master in the list
                    Info<< "    " << Foam::hostName().c_str() << '.'
                        << Foam::pid() << nl;

                    // Sub-processes
                    for (const auto& str : hostProcs)
                    {
                        Info<< "    " << str.c_str() << nl;
                    }
                    Info<< ')' << nl;
                }
            }
            if (roots.size())
            {
                Info<< "Roots  : " << roots << nl;
            }
            Info<< "Pstream initialized with:" << nl
                << "    floatTransfer      : " << Pstream::floatTransfer << nl
                << "    nProcsSimpleSum    : " << Pstream::nProcsSimpleSum << nl
                << "    commsType          : "
                << Pstream::commsTypeNames[Pstream::defaultCommsType] << nl
                << "    polling iterations : " << Pstream::nPollProcInterfaces
                << nl;
            if (UPstream::allWorlds().size() > 1)
            {
                Info<< "    worlds             : "
                    << flatOutput(UPstream::allWorlds()) << nl
                    << "    world              : " << UPstream::myWorld()
                    << nl;
            }
        }
    }

    if (initialise)
    {
        jobInfo.add("root", rootPath_);
        jobInfo.add("case", globalCase_);
        jobInfo.add("nProcs", nProcs);
        if (hostProcs.size())
        {
            jobInfo.add("hosts", hostProcs);
        }
        if (roots.size())
        {
            jobInfo.add("roots", roots);
        }
        jobInfo.write();

        // Switch on signal trapping. We have to wait until after Pstream::init
        // since this sets up its own ones.
        sigFpe::set(bannerEnabled());
        sigInt::set(bannerEnabled());
        sigQuit::set(bannerEnabled());
        sigSegv::set(bannerEnabled());

        if (Pstream::master() && bannerEnabled())
        {
            Info<< "fileModificationChecking : "
                << "Monitoring run-time modified files using "
                << IOobject::fileCheckTypesNames
                    [
                        IOobject::fileModificationChecking
                    ];
            if
            (
                IOobject::fileModificationChecking == IOobject::timeStamp
             || IOobject::fileModificationChecking == IOobject::timeStampMaster
            )
            {
                if (IOobject::maxFileModificationPolls == 1)
                {
                    Info<< " (fileModificationSkew "
                        << IOobject::fileModificationSkew
                        << ")";
                }
                else if (IOobject::maxFileModificationPolls > 1)
                {
                    Info<< " (fileModificationSkew "
                        << IOobject::fileModificationSkew
                        << ", maxFileModificationPolls "
                        << IOobject::maxFileModificationPolls
                        << ")";
                }
                else
                {
                    FatalErrorInFunction
                        << "Invalid setting for maxFileModificationPolls "
                        << IOobject::maxFileModificationPolls
                        << exit(FatalError);
                }
            }
            Info<< nl;

            Info<< "allowSystemOperations : ";
            if (dynamicCode::allowSystemOperations)
            {
                Info<< "Allowing";
            }
            else
            {
                Info<< "Disallowing";
            }
            Info<< " user-supplied system call operations" << nl
                << nl;
            IOobject::writeDivider(Info);
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::argList::~argList()
{
    jobInfo.stop();     // Normal job termination

    // Delete file handler to flush any remaining IO
    Foam::fileHandler(nullptr);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::argList::count(const UList<word>& optionNames) const
{
    label n = 0;
    for (const word& optName : optionNames)
    {
        if (options_.found(optName))
        {
            ++n;
        }
    }
    return n;
}


Foam::label Foam::argList::count
(
    std::initializer_list<word> optionNames
) const
{
    label n = 0;
    for (const word& optName : optionNames)
    {
        if (options_.found(optName))
        {
            ++n;
        }
    }
    return n;
}


bool Foam::argList::setOption(const word& optName, const string& param)
{
    // Some options are always protected
    if
    (
        optName == "case"
     || optName == "parallel"
     || optName == "roots"
    )
    {
        FatalErrorInFunction
            <<"Option: '" << optName << "' is protected" << nl
            << exit(FatalError);
        return false;
    }

    if (options_.found(optName) ? (options_[optName] != param) : true)
    {
        options_.set(optName, param);
        return true;
    }

    return false;
}


bool Foam::argList::unsetOption(const word& optName)
{
    // Some options are always protected
    if
    (
        optName == "case"
     || optName == "parallel"
     || optName == "roots"
     || optName == "hostRoots"
    )
    {
        FatalErrorInFunction
            <<"Option: '" << optName << "' is protected" << nl
            << exit(FatalError);
        return false;
    }

    // Remove the option, return true if state changed
    return options_.erase(optName);
}


void Foam::argList::displayDoc(bool source) const
{
    const dictionary& docDict = debug::controlDict().subDict("Documentation");
    fileNameList docDirs(docDict.get<fileNameList>("doxyDocDirs"));
    fileName docExt(docDict.get<fileName>("doxySourceFileExt"));

    // For source code: change xxx_8C.html to xxx_8C_source.html
    if (source)
    {
        docExt.replace(".", "_source.");
    }

    fileName url;

    for (const fileName& dir : docDirs)
    {
        // The http protocols are last in the list
        if (dir.starts_with("http:") || dir.starts_with("https:"))
        {
            url = dir/executable_ + docExt;
            break;
        }

        fileName docFile = stringOps::expand(dir/executable_ + docExt);

        if
        (
            docFile.starts_with("file://")
          ? isFile(docFile.substr(7))   // check part after "file://"
          : isFile(docFile)
        )
        {
            url = std::move(docFile);
            break;
        }
    }

    if (url.empty())
    {
        Info<< nl
            << "No documentation found for " << executable_
            << ", but you can use -help to display the usage\n" << endl;

        return;
    }

    string docBrowser = getEnv("FOAM_DOC_BROWSER");
    if (docBrowser.empty())
    {
        docDict.readEntry("docBrowser", docBrowser);
    }

    // Can use FOAM_DOC_BROWSER='application file://%f' if required
    if (docBrowser.find("%f") != std::string::npos)
    {
        docBrowser.replace("%f", url);
    }
    else
    {
        docBrowser += " " + url;
    }

    // Split on whitespace to use safer version of Foam::system()

    CStringList command(stringOps::splitSpace(docBrowser));

    Info
        << "OpenFOAM " << foamVersion::api << " documentation:" << nl
        << "    " << command << nl << endl;

    Foam::system(command, true);
}


bool Foam::argList::check(bool checkArgs, bool checkOpts) const
{
    bool ok = true;

    if (Pstream::master())
    {
        const label nargs = args_.size()-1;
        if (checkArgs && nargs != validArgs.size())
        {
            FatalError
                << "Expected " << validArgs.size()
                << " arguments but found " << nargs << endl;
            ok = false;
        }

        if (checkOpts)
        {
            forAllConstIters(options_, iter)
            {
                const word& optName = iter.key();
                if
                (
                    !validOptions.found(optName)
                 && !validParOptions.found(optName)
                )
                {
                    FatalError
                        << "Invalid option: -" << optName << endl;
                    ok = false;
                }
            }
        }

        if (!ok)
        {
            FatalError
                << nl
                << "See '" << executable_ << " -help' for usage"
                << nl << nl;
        }
    }

    return ok;
}


bool Foam::argList::checkRootCase() const
{
    if (!fileHandler().isDir(rootPath()))
    {
        FatalError
            << executable_
            << ": cannot open root directory " << rootPath()
            << endl;

        return false;
    }

    const fileName pathDir(fileHandler().filePath(path()));

    if (checkProcessorDirectories_ && pathDir.empty() && Pstream::master())
    {
        // Allow non-existent processor directories on sub-processes,
        // to be created later (e.g. redistributePar)
        FatalError
            << executable_
            << ": cannot open case directory " << path()
            << endl;

        return false;
    }

    return true;
}


// ************************************************************************* //
