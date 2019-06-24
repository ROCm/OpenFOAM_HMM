/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
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

Description
    Test atomic copyFile as per timeActivatedFileUpdate

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "OSspecific.H"
#include "Switch.H"

using namespace Foam;

#ifdef _WIN32
#undef DebugInfo    // Windows name clash with OpenFOAM messageStream
#define WIN32_LEAN_AND_MEAN
#include <windows.h>

// Prefix '\\?\' for extended-length path
inline std::string longName(const std::string& file)
{
    std::string longName;
    longName.reserve(4 + file.length());

    longName.append("\\\\?\\");
    longName.append(file);

    return longName;
}

bool win_cp(const fileName& src, const fileName& dst)
{
    const std::string srcName(longName(src));
    const std::string dstName(longName(dst));

    return ::CopyFile(srcName.c_str(), dstName.c_str(), false);
}

bool win_mv(const fileName& src, const fileName& dst)
{
    const std::string srcName(longName(src));
    const std::string dstName(longName(dst));

    return ::MoveFile(srcName.c_str(), dstName.c_str());
}

#else

bool win_cp(const fileName& src, const fileName& dst)
{
    Info<< "No Windows cp, using Foam::cp" << nl;
    return Foam::cp(src, dst);
}


bool win_mv(const fileName& src, const fileName& dst)
{
    Info<< "No Windows mv, using Foam::mv" << nl;
    return Foam::mv(src, dst);
}

#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noBanner();
    argList::noParallel();
    argList::noFunctionObjects();
    argList::addNote("Test atomic copyFile methods");

    #ifdef _WIN32
    argList::addBoolOption("win1", "Foam cp, Windows mv");
    argList::addBoolOption("win2", "Windows cp, Foam mv");
    argList::addBoolOption("win3", "Windows cp, Windows mv");
    #endif

    argList::addArgument("srcFile");
    argList::addArgument("dstFile");

    #include "setRootCase.H"

    const fileName srcFile(fileName::validate(args[1]));
    const fileName dstFile(fileName::validate(args[2]));
    const fileName tmpFile(dstFile + Foam::name(pid()));

    Info<< "src   : " << srcFile << nl
        << "tmp   : " << tmpFile << nl
        << "dst   : " << dstFile << nl
        << nl;

    bool cpOk = false, mvOk = false;


    if (args.found("win1"))
    {
        const auto usage = argList::optionUsage.cfind("win1");
        if (usage.good())
        {
            Info<< "    " << (*usage).c_str() << nl;
        }

        cpOk = Foam::cp(srcFile, tmpFile);
        mvOk = win_mv(tmpFile, dstFile);
    }
    else if (args.found("win2"))
    {
        const auto usage = argList::optionUsage.cfind("win2");
        if (usage.good())
        {
            Info<< "    " << (*usage).c_str() << nl;
        }

        cpOk = win_cp(srcFile, tmpFile);
        mvOk = Foam::mv(tmpFile, dstFile);
    }
    else if (args.found("win3"))
    {
        const auto usage = argList::optionUsage.cfind("win3");
        if (usage.good())
        {
            Info<< "    " << (*usage).c_str() << nl;
        }

        cpOk = win_cp(srcFile, tmpFile);
        mvOk = win_mv(tmpFile, dstFile);
    }
    else
    {
        Info<< "    Foam cp, Foam mv" << nl;

        cpOk = Foam::cp(srcFile, tmpFile);
        mvOk = Foam::mv(tmpFile, dstFile);
    }

    Info<< nl
        << "cp: " << Switch(cpOk) << nl
        << "mv: " << Switch(mvOk) << nl;

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
