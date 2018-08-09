/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2018 OpenCFD Ltd.
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
#include "clock.H"
#include "IFstream.H"
#include "dictionary.H"
#include "IOobject.H"
#include "JobInfo.H"
#include "labelList.H"
#include "regIOobject.H"
#include "dynamicCode.H"
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
#include "IOmanip.H"

#include <cctype>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

bool Foam::argList::argsMandatory_ = true;
bool Foam::argList::checkProcessorDirectories_ = true;
Foam::SLList<Foam::string>    Foam::argList::validArgs;
Foam::HashSet<Foam::string>   Foam::argList::advancedOptions;
Foam::HashTable<Foam::string> Foam::argList::validOptions;
Foam::HashTable<Foam::string> Foam::argList::validParOptions;
Foam::HashTable<Foam::string> Foam::argList::optionUsage;
Foam::HashTable<std::pair<Foam::word,int>> Foam::argList::validOptionsCompat;
Foam::HashTable<std::pair<bool,int>> Foam::argList::ignoreOptionsCompat;
Foam::SLList<Foam::string>    Foam::argList::notes;
Foam::string::size_type Foam::argList::usageMin = 20;
Foam::string::size_type Foam::argList::usageMax = 80;
Foam::word Foam::argList::postProcessOptionName("postProcess");

Foam::argList::initValidTables::initValidTables()
{
    argList::addOption
    (
        "case", "dir",
        "Specify case directory to use (instead of the cwd)"
    );
    argList::addBoolOption("parallel", "Run in parallel");
    validParOptions.set("parallel", "");
    argList::addOption
    (
        "roots", "(dir1 .. dirN)",
        "Slave root directories for distributed running",
        true  // advanced option
    );
    validParOptions.set("roots", "(dir1 .. dirN)");

    argList::addOption
    (
        "decomposeParDict", "file",
        "Use specified file for decomposePar dictionary"
    );
    argList::addOption
    (
        "hostRoots", "(((host1 dir1) .. (hostN dirN))",
        "slave root directories (per host) for distributed running. "
        "The host specification can use a regex.",
        true  // advanced option
    );
    validParOptions.set("hostRoots", "((host1 dir1) .. (hostN dirN))");

    argList::addBoolOption
    (
        "noFunctionObjects",
        "Do not execute function objects"
    );

    argList::addOption
    (
        "fileHandler", "handler",
        "Override the file handler type",
        true  // advanced option
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

namespace
{

// Should issue warning if there is +ve versioning (+ve version number)
// and the this version number is older than the current OpenFOAM version
// as conveyed by the OPENFOAM compiler define.

static inline constexpr bool shouldWarnVersion(const int version)
{
    return (version > 0 && version < OPENFOAM);
}

} // End anonymous namespace


namespace Foam
{

// Counted per machine name
// Does not include any sorting since we wish to know the ordering according to
// mpi rank.
//
// Always include the master too.
// This provides a better overview of the subscription
static void printHostsSubscription(const UList<string>& slaveProcs)
{
    Info<< "Hosts  :" << nl << "(" << nl;

    std::string prev = hostName();
    int count = 1;

    for (const auto& str : slaveProcs)
    {
        std::string curr(str.substr(0, str.rfind('.')));

        if (prev != curr)
        {
            if (count)
            {
                // Finish previous
                Info<<"    (" << prev.c_str() << " " << count << ")" << nl;
                count = 0;
            }

            prev = std::move(curr);
        }
        ++count;
    }

    if (count)
    {
        // Finished last one
        Info<<"    (" << prev.c_str() << " " << count << ")" << nl;
    }

    Info<< ")" << nl;
}


// Print information about version, build, arch
static void printBuildInfo(const bool full=true)
{
    Info<<"Using: OpenFOAM-" << Foam::FOAMversion
        << " (see www.OpenFOAM.com)" << nl
        << "Build: " << Foam::FOAMbuild << nl;

    if (full)
    {
        Info << "Arch:  " << Foam::FOAMbuildArch << nl;
    }
}

} // End namespace Foam


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

void Foam::argList::warnTrailing(const ITstream& is, const label index)
{
    const label nExcess = is.nRemainingTokens();

    if (nExcess)
    {
        std::cerr
            << nl
            << "--> FOAM WARNING:" << nl
            << "argument " << index << " has "
            << nExcess << " excess tokens" << nl << nl;
    }
}


void Foam::argList::warnTrailing(const ITstream& is, const word& optName)
{
    const label nExcess = is.nRemainingTokens();

    if (nExcess)
    {
        std::cerr
            << nl
            << "--> FOAM WARNING:" << nl
            << "option -" << optName << " has "
            << nExcess << " excess tokens" << nl << nl;
    }
}


void Foam::argList::addArgument(const string& argName)
{
    validArgs.append(argName);
}


void Foam::argList::addBoolOption
(
    const word& optName,
    const string& usage,
    const bool advanced
)
{
    addOption(optName, "", usage, advanced);
}


void Foam::argList::addOption
(
    const word& optName,
    const string& param,
    const string& usage,
    const bool advanced
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
    const bool expectArg
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


void Foam::argList::nonMandatoryArgs()
{
    argsMandatory_ = false;
}


void Foam::argList::noBanner()
{
    ::Foam::infoDetailLevel = 0;
}


bool Foam::argList::bannerEnabled()
{
    return (::Foam::infoDetailLevel > 0);
}


void Foam::argList::noFunctionObjects(bool addWithOption)
{
    removeOption("noFunctionObjects");

    // Ignore this bool option without warning
    // - cannot tie to any particular version anyhow
    ignoreOptionCompat({"noFunctionObjects", 0}, false);

    if (addWithOption)
    {
        addBoolOption
        (
            "withFunctionObjects",
            "Execute functionObjects"
        );
    }
}


void Foam::argList::noJobInfo()
{
    JobInfo::writeJobInfo = false;
}


void Foam::argList::noLibs()
{
    addBoolOption
    (
        "no-libs",
        "Disable use of the controlDict libs entry",
        true  // advanced
    );
}


void Foam::argList::noParallel()
{
    removeOption("parallel");
    removeOption("roots");
    removeOption("decomposeParDict");
    removeOption("hostRoots");
    validParOptions.clear();
}


void Foam::argList::noCheckProcessorDirectories()
{
    checkProcessorDirectories_ = false;
}


void Foam::argList::printOptionUsage
(
    const label location,
    const string& str
)
{
    const string::size_type textWidth = usageMax - usageMin;
    const string::size_type strLen = str.size();

    if (strLen)
    {
        // Minimum of 2 spaces between option and usage:
        if (string::size_type(location) + 2 <= usageMin)
        {
            for (string::size_type i = location; i < usageMin; ++i)
            {
                Info<<' ';
            }
        }
        else
        {
            // or start a new line
            Info<< nl;
            for (string::size_type i = 0; i < usageMin; ++i)
            {
                Info<<' ';
            }
        }

        // Text wrap
        string::size_type pos = 0;
        while (pos != string::npos && pos + textWidth < strLen)
        {
            // Potential end point and next point
            string::size_type curr = pos + textWidth - 1;
            string::size_type next = string::npos;

            if (isspace(str[curr]))
            {
                // We were lucky: ended on a space
                next = str.find_first_not_of(" \t\n", curr);
            }
            else if (isspace(str[curr+1]))
            {
                // The next one is a space - so we are okay
                ++curr;  // otherwise the length is wrong
                next = str.find_first_not_of(" \t\n", curr);
            }
            else
            {
                // Search for end of a previous word break
                string::size_type prev = str.find_last_of(" \t\n", curr);

                // Reposition to the end of previous word if possible
                if (prev != string::npos && prev > pos)
                {
                    curr = prev;
                }
            }

            if (next == string::npos)
            {
                next = curr + 1;
            }

            // Indent following lines (not the first one)
            if (pos)
            {
                for (string::size_type i = 0; i < usageMin; ++i)
                {
                    Info<<' ';
                }
            }

            Info<< str.substr(pos, (curr - pos)).c_str() << nl;
            pos = next;
        }

        // Output the remainder of the string
        if (pos != string::npos)
        {
            // Indent following lines (not the first one)
            if (pos)
            {
                for (string::size_type i = 0; i < usageMin; ++i)
                {
                    Info<<' ';
                }
            }

            Info<< str.substr(pos).c_str() << nl;
        }
    }
    else
    {
        Info<< nl;
    }
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
            const auto& iter = *fnd;

            if (shouldWarnVersion(iter.second))
            {
                std::cerr
                    << "--> FOAM IOWarning :" << nl
                    << "    Found [v" << iter.second << "] '"
                    << optName << "' instead of '-"
                    << iter.first << "' option"
                    << nl
                    << std::endl;

                error::warnAboutAge("option", iter.second);
            }

            return "-" + iter.first;
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
            const auto& iter = *fnd;

            // Number to skip (including the option itself)
            // '-option ARG' or '-option'
            const int nskip = (iter.first ? 2 : 1);

            if (shouldWarnVersion(iter.second))
            {
                std::cerr
                    << "--> FOAM IOWarning :" << nl
                    << "    Ignoring [v" << iter.second << "] '-"
                    << optName << (nskip > 1 ? " ARG" : "")
                    << "' option"
                    << nl
                    << std::endl;

                error::warnAboutAge("option", iter.second);
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
    forAll(args_, argi)
    {
        len += args_[argi].size();
    }

    // Length needed for regrouped command-line
    argListStr_.reserve(len);

    return nArgs < argc;
}


void Foam::argList::getRootCase()
{
    fileName casePath;

    // [-case dir] specified
    const auto optIter = options_.cfind("case");

    if (optIter.found())
    {
        casePath = optIter.object();
        casePath.clean();

        if (casePath.empty() || casePath == ".")
        {
            // Handle degenerate form and '-case .' like no -case specified
            casePath = cwd();
            options_.erase("case");
        }
        else if (!casePath.isAbsolute() && casePath.name() == "..")
        {
            // Avoid relative cases ending in '..' - makes for very ugly names
            casePath = cwd()/casePath;
            casePath.clean();
        }
    }
    else
    {
        // Nothing specified, use the current dir
        casePath = cwd();
    }

    rootPath_   = casePath.path();
    globalCase_ = casePath.name();
    case_       = globalCase_;

    // The name of the executable, unless already present in the environment
    setEnv("FOAM_EXECUTABLE", executable_, false);

    // Set the case and case-name as an environment variable
    if (rootPath_.isAbsolute())
    {
        // Absolute path - use as-is
        setEnv("FOAM_CASE", rootPath_/globalCase_, true);
        setEnv("FOAM_CASENAME", globalCase_, true);
    }
    else
    {
        // Qualify relative path
        casePath = cwd()/rootPath_/globalCase_;
        casePath.clean();

        setEnv("FOAM_CASE", casePath, true);
        setEnv("FOAM_CASENAME", casePath.name(), true);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::argList::argList
(
    int& argc,
    char**& argv,
    bool checkArgs,
    bool checkOpts,
    const bool initialise
)
:
    args_(argc),
    options_(argc),
    distributed_(false)
{
    // Check for -fileHandler, which requires an argument.
    word handlerType(getEnv("FOAM_FILEHANDLER"));
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
        handlerType = fileOperation::defaultFileHandler;
    }

    // Detect any parallel options
    bool needsThread = fileOperations::fileOperationInitialise::New
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
                parRunControl_.runPar(argc, argv, needsThread);
                break;
            }
        }
    }

    // Convert argv -> args_ and capture ( ... ) lists
    regroupArgv(argc, argv);
    argListStr_ += args_[0];

    // Set executable name immediately - useful when emitting errors.
    executable_ = fileName(args_[0]).name();

    // Check arguments and options, argv[0] was already handled
    int nArgs = 1;
    for (int argi = 1; argi < args_.size(); ++argi)
    {
        argListStr_ += ' ';
        argListStr_ += args_[argi];

        if (args_[argi][0] == '-')
        {
            const char *optName = &args_[argi][1];

            if (!*optName)
            {
                Warning
                    <<"Ignoring lone '-' on the command-line" << endl;
            }
            else if
            (
                validOptions.lookup(optName, "").size()
             || validParOptions.lookup(optName, "").size()
            )
            {
                // If the option is known to require an argument,
                // get it or emit a FatalError.

                ++argi;
                if (argi >= args_.size())
                {
                    printBuildInfo(false);

                    Info<<nl
                        <<"Error: option '-" << optName
                        << "' requires an argument" << nl << nl
                        << "See '" << executable_ << " -help' for usage"
                        << nl << nl;

                    Pstream::exit(1); // works for serial and parallel
                }

                argListStr_ += ' ';
                argListStr_ += args_[argi];
                // Handle duplicates by taking the last -option specified
                options_.set(optName, args_[argi]);
            }
            else
            {
                // All other options (including unknown ones) are simply
                // registered as existing.

                options_.insert(optName, "");
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
    parRunControl_(args.parRunControl_),
    args_(args.args_),
    options_(options),
    executable_(args.executable_),
    rootPath_(args.rootPath_),
    globalCase_(args.globalCase_),
    case_(args.case_),
    argListStr_(args.argListStr_)
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
        else if
        (
            options_.found("doc-source")
         || options_.found("srcDoc")  // Compat 1706
        )
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
        else if (options_.found("help"))
        {
            printUsage(false);
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
            ::exit(0);
        }
    }

    // Print the collected error messages and exit if check fails
    if (!check(checkArgs, checkOpts))
    {
        printBuildInfo(false);
        FatalError.write(Info, false);

        Pstream::exit(1); // works for serial and parallel
    }

    if (initialise)
    {
        const string dateString = clock::date();
        const string timeString = clock::clockTime();

        // Print the banner once only for parallel runs
        if (Pstream::master() && bannerEnabled())
        {
            IOobject::writeBanner(Info, true)
                << "Build  : " << Foam::FOAMbuild << nl
                << "Arch   : " << Foam::FOAMbuildArch << nl
                << "Exec   : " << argListStr_.c_str() << nl
                << "Date   : " << dateString.c_str() << nl
                << "Time   : " << timeString.c_str() << nl
                << "Host   : " << hostName() << nl
                << "PID    : " << pid() << endl;
        }

        jobInfo.add("startDate", dateString);
        jobInfo.add("startTime", timeString);
        jobInfo.add("userName", userName());
        jobInfo.add("foamVersion", word(Foam::FOAMversion));
        jobInfo.add("code", executable_);
        jobInfo.add("argList", argListStr_);
        jobInfo.add("currentDir", cwd());
        jobInfo.add("PPID", ppid());
        jobInfo.add("PGID", pgid());

        // Add build information - only use the first word
        {
            std::string build(Foam::FOAMbuild);
            std::string::size_type space = build.find(' ');
            if (space != std::string::npos)
            {
                build.resize(space);
            }
            jobInfo.add("foamBuild", build);
        }
    }


    // Set fileHandler. In increasing order of priority:
    // 1. default = uncollated
    // 2. environment var FOAM_FILEHANDLER
    // 3. etc/controlDict optimisationSwitches 'fileHandler'
    // 4. system/controlDict 'fileHandler' (not handled here; done in TimeIO.C)
    // 5. '-fileHandler' commmand-line option

    {
        word fileHandlerName =
            options_.lookup("fileHandler", getEnv("FOAM_FILEHANDLER"));

        if (fileHandlerName.empty())
        {
            fileHandlerName = fileOperation::defaultFileHandler;
        }

        auto handler = fileOperation::New(fileHandlerName, bannerEnabled());
        Foam::fileHandler(handler);
    }


    stringList slaveProcs;
    stringList slaveMachine;
    const int writeHostsSwitch = debug::infoSwitch("writeHosts", 1);

    // Collect slave machine/pid, and check that the build is identical
    if (parRunControl_.parRun())
    {
        if (Pstream::master())
        {
            slaveProcs.resize(Pstream::nProcs()-1);
            slaveMachine.resize(Pstream::nProcs()-1);
            label proci = 0;
            for
            (
                int slave = Pstream::firstSlave();
                slave <= Pstream::lastSlave();
                slave++
            )
            {
                IPstream fromSlave(Pstream::commsTypes::scheduled, slave);

                string slaveBuild;
                label slavePid;
                fromSlave >> slaveBuild >> slaveMachine[proci] >> slavePid;

                slaveProcs[proci] = slaveMachine[proci] + "." + name(slavePid);
                proci++;

                // Verify that all processors are running the same build
                if (slaveBuild != Foam::FOAMbuild)
                {
                    FatalErrorIn(executable())
                        << "Master is running version " << Foam::FOAMbuild
                        << "; slave " << proci << " is running version "
                        << slaveBuild
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
            toMaster << string(Foam::FOAMbuild) << hostName() << pid();
        }
    }


    // Case is a single processor run unless it is running parallel
    int nProcs = 1;

    // Roots if running distributed
    fileNameList roots;

    // If this actually is a parallel run
    if (parRunControl_.parRun())
    {
        // For the master
        if (Pstream::master())
        {
            // Establish rootPath_/globalCase_/case_ for master
            getRootCase();

            // Establish location of decomposeParDict, allow override with
            // the -decomposeParDict option.
            fileName source = rootPath_/globalCase_/"system"/"decomposeParDict";
            if (options_.found("decomposeParDict"))
            {
                bool adjustOpt = false;

                source = options_["decomposeParDict"];
                if (isDir(source))
                {
                    source = source/"decomposeParDict";
                    adjustOpt = true;
                }

                // Case-relative if not absolute and not "./" etc
                if (!source.isAbsolute() && !source.startsWith("."))
                {
                    source = rootPath_/globalCase_/source;
                    adjustOpt = true;
                }

                // Could also check for absolute path, but shouldn't be needed
                if (adjustOpt)
                {
                    source.clean();
                    options_.set("decomposeParDict", source);
                }
            }

            // If running distributed (different roots for different procs)
            label dictNProcs = -1;
            if (this->readListIfPresent("roots", roots))
            {
                distributed_ = true;
                source = "-roots";
                if (roots.size() != 1)
                {
                    dictNProcs = roots.size()+1;
                }
            }
            else if (options_.found("hostRoots"))
            {
                roots.resize(Pstream::nProcs()-1, fileName::null);

                source = "-hostRoots";
                ITstream is = this->lookup("hostRoots");
                List<Tuple2<wordRe, fileName>> hostRoots(is);
                warnTrailing(is, "hostRoots");

                for (const auto& hostRoot : hostRoots)
                {
                    const wordRe& re = hostRoot.first();
                    labelList matched(findStrings(re, slaveMachine));
                    for (const label slavei : matched)
                    {
                        if (!roots[slavei].empty())
                        {
                            FatalErrorInFunction
                                << "Slave " << slaveMachine[slavei]
                                << " has multiple matching roots in "
                                << hostRoots << exit(FatalError);
                        }

                        roots[slavei] = hostRoot.second();
                    }
                }

                // Check
                forAll(roots, slavei)
                {
                    if (roots[slavei].empty())
                    {
                        FatalErrorInFunction
                            << "Slave " << slaveMachine[slavei]
                            << " has no matching roots in "
                            << hostRoots << exit(FatalError);
                    }
                }

                if (roots.size() != 1)
                {
                    dictNProcs = roots.size()+1;
                }
            }
            else if (checkProcessorDirectories_)
            {
                // Use values from decomposeParDict, the location was already
                // established above.

                IFstream decompDictStream(source);

                if (!decompDictStream.good())
                {
                    FatalError
                        << "Cannot read decomposeParDict from "
                        << decompDictStream.name()
                        << exit(FatalError);
                }

                dictionary decompDict(decompDictStream);

                dictNProcs = decompDict.get<label>("numberOfSubdomains");

                if (decompDict.lookupOrDefault("distributed", false))
                {
                    distributed_ = true;
                    decompDict.readEntry("roots", roots);
                }
            }

            // Convenience:
            // when a single root is specified, use it for all processes
            if (roots.size() == 1)
            {
                const fileName rootName(roots[0]);
                roots.resize(Pstream::nProcs()-1, rootName);

                // adjust dictNProcs for command-line '-roots' option
                if (dictNProcs < 0)
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
            if (checkProcessorDirectories_ && dictNProcs > Pstream::nProcs())
            {
                FatalError
                    << source
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
                        << " is not equal to the number of slaves "
                        << Pstream::nProcs()-1
                        << exit(FatalError);
                }

                forAll(roots, i)
                {
                    roots[i].expand();
                }

                // Distribute the master's argument list (with new root)
                const bool hadCaseOpt = options_.found("case");
                for
                (
                    int slave = Pstream::firstSlave();
                    slave <= Pstream::lastSlave();
                    slave++
                )
                {
                    options_.set("case", roots[slave-1]/globalCase_);

                    OPstream toSlave(Pstream::commsTypes::scheduled, slave);
                    toSlave << args_ << options_ << roots.size();
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
                 && dictNProcs < Pstream::nProcs()
                )
                {
                    label nProcDirs = 0;
                    while
                    (
                        isDir
                        (
                            rootPath_/globalCase_/"processor"
                          + name(++nProcDirs)
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
                for
                (
                    int slave = Pstream::firstSlave();
                    slave <= Pstream::lastSlave();
                    slave++
                )
                {
                    OPstream toSlave(Pstream::commsTypes::scheduled, slave);
                    toSlave << args_ << options_ << roots.size();
                }
            }
        }
        else
        {
            // Collect the master's argument list
            IPstream fromMaster
            (
                Pstream::commsTypes::scheduled,
                Pstream::masterNo()
            );
            fromMaster >> args_ >> options_ >> distributed_;

            // Establish rootPath_/globalCase_/case_ for slave
            getRootCase();
        }

        nProcs = Pstream::nProcs();
        case_ = globalCase_/(word("processor") + name(Pstream::myProcNo()));
    }
    else
    {
        // Establish rootPath_/globalCase_/case_
        getRootCase();
        case_ = globalCase_;
    }

    // Keep or discard slave and root information for reporting:
    if (Pstream::master() && parRunControl_.parRun())
    {
        if (!writeHostsSwitch)
        {
            // Clear here to ensures it doesn't show in the jobInfo
            slaveProcs.clear();
        }
        if (!debug::infoSwitch("writeRoots", 1))
        {
            roots.clear();
        }
    }

    if (Pstream::master() && bannerEnabled())
    {
        Info<< "Case   : " << (rootPath_/globalCase_).c_str() << nl
            << "nProcs : " << nProcs << endl;

        if (parRunControl_.parRun())
        {
            if (slaveProcs.size())
            {
                if (writeHostsSwitch == 1)
                {
                    // Compact output (see etc/controlDict)
                    printHostsSubscription(slaveProcs);
                }
                else
                {
                    // Full output of "slave.pid"
                    Info<< "Slaves : " << slaveProcs << nl;
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
                << endl;
        }
    }

    if (initialise)
    {
        jobInfo.add("root", rootPath_);
        jobInfo.add("case", globalCase_);
        jobInfo.add("nProcs", nProcs);
        if (slaveProcs.size())
        {
            jobInfo.add("slaves", slaveProcs);
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
                << regIOobject::fileCheckTypesNames
                    [
                        regIOobject::fileModificationChecking
                    ];
            if
            (
                (
                    regIOobject::fileModificationChecking
                 == regIOobject::timeStamp
                )
             || (
                    regIOobject::fileModificationChecking
                 == regIOobject::timeStampMaster
                )
            )
            {
                Info<< " (fileModificationSkew "
                    << regIOobject::fileModificationSkew << ")";
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
                << endl;
            IOobject::writeDivider(Info);
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::argList::~argList()
{
    jobInfo.end();

    // Delete file handler to flush any remaining IO
    autoPtr<fileOperation> dummy(nullptr);
    fileHandler(dummy);
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


void Foam::argList::printNotes() const
{
    // Output notes directly - no automatic text wrapping
    if (!notes.empty())
    {
        Info<< nl;
        forAllConstIters(notes, iter)
        {
            Info<< iter().c_str() << nl;
        }
    }
}


void Foam::argList::printUsage(bool full) const
{
    Info<< "\nUsage: " << executable_ << " [OPTIONS]";

    if (validArgs.size())
    {
        Info<< ' ';

        if (!argsMandatory_)
        {
            Info<< '[';
        }

        label i = 0;
        forAllConstIters(validArgs, iter)
        {
            if (i++) Info<< ' ';
            Info<< '<' << iter().c_str() << '>';
        }

        if (!argsMandatory_)
        {
            Info<< ']';
        }
    }

    Info<< "\noptions:\n";

    for (const word& optName : validOptions.sortedToc())
    {
        // Suppress advanced options for regular -help.
        if (advancedOptions.found(optName) && !full)
        {
            continue;
        }

        Info<< "  -" << optName;
        label len = optName.size() + 3;  // Length includes leading '  -'

        const auto optIter = validOptions.cfind(optName);
        if (optIter().size())
        {
            // Length includes space between option/param and '<>'
            len += optIter().size() + 3;
            Info<< " <" << optIter().c_str() << '>';
        }

        const auto usageIter = optionUsage.cfind(optName);
        if (usageIter.found())
        {
            printOptionUsage(len, usageIter());
        }
        else
        {
            Info<< nl;
        }
    }

    // Place documentation/help options at the end

    Info<< "  -doc";
    printOptionUsage(6, "Display documentation in browser");

    Info<< "  -doc-source";
    printOptionUsage(13, "Display source code in browser" );

    Info<< "  -help";
    printOptionUsage(7, "Display short help and exit");

    if
    (
        argList::validOptionsCompat.size()
      + argList::ignoreOptionsCompat.size()
    )
    {
        Info<< "  -help-compat";
        printOptionUsage(14, "Display compatibility options and exit");
    }

    Info<< "  -help-full";
    printOptionUsage(12, "Display full help and exit");

    printNotes();

    Info<< nl;
    printBuildInfo();
    Info<< endl;
}


void Foam::argList::printCompat() const
{
    const label nopt
    (
        argList::validOptionsCompat.size()
      + argList::ignoreOptionsCompat.size()
    );

    Info<< nopt << " compatibility options for " << executable_ << nl;
    if (nopt)
    {
        Info<< setf(ios_base::left) << setw(34) << "old option"
            << setf(ios_base::left) << setw(32) << "new option"
            << "version" << nl
            << setf(ios_base::left) << setw(34) << "=========="
            << setf(ios_base::left) << setw(32) << "=========="
            << "=======" << nl;

        for (const word& k : argList::validOptionsCompat.sortedToc())
        {
            const auto& iter = *argList::validOptionsCompat.cfind(k);

            Info<< "  -"
                << setf(ios_base::left) << setw(32) << k
                << " -"
                << setf(ios_base::left) << setw(30) << iter.first << ' ';

            if (iter.second > 0)
            {
                Info<< (iter.second);
            }
            else if (iter.second < 0)
            {
                Info<< (-iter.second);
            }
            else
            {
                Info<< '-';
            }
            Info<< nl;
        }

        for (const word& k : argList::ignoreOptionsCompat.sortedToc())
        {
            const auto& iter = *argList::ignoreOptionsCompat.cfind(k);

            if (iter.first)
            {
                Info<< "  -"
                    << setf(ios_base::left) << setw(32)
                    << (k + " <arg>").c_str();
            }
            else
            {
                Info<< "  -"
                    << setf(ios_base::left) << setw(32) << k;
            }

            Info<< setf(ios_base::left) << setw(32) << " removed" << ' ';
            if (iter.second > 0)
            {
                Info<< (iter.second);
            }
            else if (iter.second < 0)
            {
                Info<< (-iter.second);
            }
            else
            {
                Info<< '-';
            }
            Info<< nl;
        }
    }
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
        // http protocols are last in the list
        if (dir.startsWith("http:") || dir.startsWith("https:"))
        {
            url = dir/executable_ + docExt;
            break;
        }

        fileName docFile = stringOps::expand(dir/executable_ + docExt);

        if
        (
            docFile.startsWith("file://")
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

    Info<<"OpenFOAM-" << Foam::FOAMversion << " documentation:" << nl
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
        // Allow slaves on non-existing processor directories, created later
        // (e.g. redistributePar)
        FatalError
            << executable_
            << ": cannot open case directory " << path()
            << endl;

        return false;
    }

    return true;
}


// ************************************************************************* //
