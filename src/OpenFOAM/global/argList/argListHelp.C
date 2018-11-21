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
#include "stringOps.H"

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


static void printManOption(const word& optName)
{
    Info<< ".TP" << nl << ".B \\-" << optName;

    // Option has arg?
    const auto optIter = argList::validOptions.cfind(optName);
    if (optIter().size())
    {
        Info<< " <" << optIter().c_str() << '>';
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
        << " \\- part of OpenFOAM (The Open Source CFD Toolbox)."
        << nl;


    // .SH SYNOPSIS
    // .B command [OPTIONS] ...

    Info
        << ".SH \"SYNOPSIS\"" << nl
        << ".B " << executable_ << " [OPTIONS]";

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
            Info << '<' << argName.c_str() << '>';
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
            Info<<"No description available\n";
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

    Info<< ".TP" << nl << ".B \\-" << "doc" << nl
        <<"Display documentation in browser" << nl;

    Info<< ".TP" << nl << ".B \\-" << "doc-source" << nl
        << "Display source code in browser" << nl;

    Info<< ".TP" << nl << ".B \\-" << "help" << nl
        << "Display short help and exit" << nl;

    Info<< ".TP" << nl << ".B \\-" << "help-full" << nl
        << "Display full help and exit" << nl;


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
        Info<< ".TP" << nl << ".B \\-" << "help-compat" << nl
            << "Display compatibility options and exit" << nl;
    }

    Info<< ".TP" << nl << ".B \\-" << "help-man" << nl
        << "Display full help (manpage format) and exit" << nl;

    // Footer
    printManFooter();
}


// ************************************************************************* //
