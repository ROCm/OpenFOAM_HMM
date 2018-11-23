/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
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

#include "argList.H"
#include "foamVersion.H"
#include "stringOps.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// manpage Footer
static inline void printManFooter()
{
    Info<< ".SH \"SEE ALSO\"" << nl
        << "Online documentation "
        << "https://www.openfoam.com/documentation/" << nl
        << ".SH COPYRIGHT" << nl
        << "Copyright 2018 OpenCFD Ltd." << nl;
}


// Regular option
static void printManOption(const word& optName)
{
    Info<< ".TP\n\\fB\\-" << optName << "\\fR";

    // Option has arg?
    const auto optIter = argList::validOptions.cfind(optName);
    if (optIter().size())
    {
        Info<< " \\fI" << optIter().c_str() << "\\fR";
    }
    Info<< nl;

    // Option has usage information?

    const auto usageIter = argList::optionUsage.cfind(optName);
    if (usageIter.found())
    {
        stringOps::writeWrapped(Info, *usageIter, argList::usageMax, 0, true);
    }
    else
    {
        Info<< nl;
    }
}


// Simple, hard-coded option
static inline void printManOption(const char* optName, const char* optUsage)
{
    Info<< ".TP\n\\fB\\-" << optName << "\\fR" << nl
        << optUsage << nl;
}


static void printOptionUsage
(
    std::string::size_type start,
    const string& str
)
{
    if (str.empty())
    {
        Info<< nl;
        return;
    }

    // Indent the first line. Min 2 spaces between option and usage
    if (start + 2 > argList::usageMin)
    {
        Info<< nl;
        start = 0;
    }
    while (start < argList::usageMin)
    {
        Info<<' ';
        ++start;
    }

    stringOps::writeWrapped
    (
        Info,
        str,
        (argList::usageMax - argList::usageMin),
        argList::usageMin
    );
}


// Simple, hard-coded option
static inline void printOption(const char* optName, const char* optUsage)
{
    Info<< "  -" << optName;
    printOptionUsage(3 + strlen(optName), optUsage);
}

} // End namespace Foam


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::argList::printMan() const
{
    // .TH "<APPLICATION>" 1 "OpenFOAM-<version>" "source" "category"

    Info
        << ".TH" << token::SPACE
        // All uppercase and quoted
        << stringOps::upper(executable_) << token::SPACE
        << "\"1\"" << token::SPACE
        << token::DQUOTE << "OpenFOAM-v" << OPENFOAM << token::DQUOTE
        << token::SPACE
        << token::DQUOTE << "www.openfoam.com" << token::DQUOTE
        << token::SPACE
        << token::DQUOTE << "OpenFOAM Commands Manual" << token::DQUOTE
        << nl;


    // .SH NAME
    // <application> \- part of OpenFOAM (The Open Source CFD Toolbox).
    Info
        << ".SH \"NAME\"" << nl
        << executable_
        << " \\- part of \\fBOpenFOAM\\fR (The Open Source CFD Toolbox)."
        << nl;


    // .SH SYNOPSIS
    // .B command [OPTIONS] ...

    Info
        << ".SH \"SYNOPSIS\"" << nl
        << "\\fB" << executable_ << "\\fR [\\fIOPTIONS\\fR]";

    if (validArgs.size())
    {
        Info<< ' ';

        if (!argsMandatory_)
        {
            Info<< '[';
        }

        label i = 0;
        for (const std::string& argName : validArgs)
        {
            if (i++) Info<< ' ';
            Info << "\\fI" << argName.c_str() << "\\fR";
        }

        if (!argsMandatory_)
        {
            Info<< ']';
        }
    }
    Info<< nl;


    // .SH DESCRIPTION
    {
        Info
            << ".SH \"DESCRIPTION\"" << nl;

        Info<< ".nf" << nl;

        if (notes.empty())
        {
            Info<< "No description available\n";
        }
        else
        {
            Info<< nl;
            for (const std::string& note : notes)
            {
                if (note.empty())
                {
                    Info<< nl;
                }
                else
                {
                    stringOps::writeWrapped(Info, note, usageMax, 0, true);
                }
            }
        }
        Info<< ".fi" << nl;
    }


    // .SH "OPTIONS"
    Info
        << ".SH \"OPTIONS\"" << nl;

    for (const word& optName : validOptions.sortedToc())
    {
        // Normal options
        if (!advancedOptions.found(optName))
        {
            printManOption(optName);
        }
    }

    // Standard documentation/help options

    printManOption("doc", "Display documentation in browser");
    printManOption("doc-source", "Display source code in browser");
    printManOption("help", "Display short help and exit");
    printManOption("help-full", "Display full help and exit");


    // .SH "ADVANCED OPTIONS"
    Info
        << ".SH \"ADVANCED OPTIONS\"" << nl;

    for (const word& optName : validOptions.sortedToc())
    {
        // Advanced options
        if (advancedOptions.found(optName))
        {
            printManOption(optName);
        }
    }

    // Compatibility information
    if
    (
        argList::validOptionsCompat.size()
      + argList::ignoreOptionsCompat.size()
    )
    {
        printManOption("help-compat", "Display compatibility options and exit");
    }

    printManOption("help-man", "Display full help (manpage format) and exit");

    // Footer
    printManFooter();
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
        for (const std::string& argName : validArgs)
        {
            if (i++) Info<< ' ';
            Info<< '<' << argName.c_str() << '>';
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

        // Length includes leading '  -'
        label len = optName.size() + 3;

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

    printOption("doc", "Display documentation in browser");
    printOption("doc-source", "Display source code in browser");
    printOption("help", "Display short help and exit");

    if
    (
        argList::validOptionsCompat.size()
      + argList::ignoreOptionsCompat.size()
    )
    {
        printOption("help-compat", "Display compatibility options and exit");
    }

    if (full)
    {
        printOption("help-man", "Display full help (manpage format) and exit");
    }

    printOption("help-full", "Display full help and exit");

    printNotes();

    Info<< nl;
    foamVersion::printBuildInfo(true);
    Info<< endl;
}


void Foam::argList::printNotes() const
{
    // Output notes with automatic text wrapping
    if (!notes.empty())
    {
        Info<< nl;

        for (const std::string& note : notes)
        {
            if (note.empty())
            {
                Info<< nl;
            }
            else
            {
                stringOps::writeWrapped(Info, note, usageMax);
            }
        }
    }
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


// ************************************************************************* //
