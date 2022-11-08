/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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
    Test OFstream. Primarily atomic operations

\*---------------------------------------------------------------------------*/

#include "Fstream.H"
#include "IOstreams.H"
#include "OSspecific.H"
#include "argList.H"
#include "ListOps.H"

using namespace Foam;

void listFiles(const fileName& dir)
{
    wordList files = ListOps::create<word>
    (
        readDir(dir, fileName::FILE),
        nameOp<fileName>()
    );

    Info
        << nl
        << "files:" << nl
        << files << nl
        << "ls" << nl
        << "============" << endl;
    Foam::system("ls -al " + dir);
    Info<< "============" << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::addBoolOption("gz", "Use compression");
    argList::addBoolOption("append", "Use append mode");
    argList::addBoolOption("atomic", "Use atomic");
    argList::addBoolOption("keep", "Do not remove test directory");
    argList::addOption("write", "file", "test writing to file");

    #include "setRootCase.H"

    const fileName baseDir("Test-OFstream-directory");

    Foam::mkDir(baseDir);

    InfoErr<< "mkdir: " << baseDir << endl;

    IOstreamOption streamOpt;

    if (args.found("gz"))
    {
        streamOpt.compression(IOstreamOption::COMPRESSED);
    }

    IOstreamOption::appendType append =
    (
        args.found("append")
      ? IOstreamOption::APPEND
      : IOstreamOption::NON_APPEND
    );
    IOstreamOption::atomicType atomic =
    (
        args.found("atomic")
      ? IOstreamOption::ATOMIC
      : IOstreamOption::NON_ATOMIC
    );

    {
        OFstream(baseDir/"dummy")() << "Some file content" << endl;

        Foam::ln("dummy", baseDir/"Test2.txt");
        Foam::ln("dummy", baseDir/"Test3.txt");
        Foam::ln("dummy", baseDir/"Test4.txt");
        Foam::ln("dummy", baseDir/"Test4.txt.gz");
        Foam::ln("dummy", baseDir/"Test5.txt");
        Foam::ln("dummy", baseDir/"Test5.txt.gz");
    }

    {
        OFstream os
        (
            atomic,
            baseDir/"Test1.txt",
            streamOpt,
            append
        );

        os << "=========================" << endl;

        InfoErr<< "open: " << os.name() << endl;
        InfoErr<< "... sleep" << endl;

        listFiles(baseDir);

        sleep(2);

        os << "+++++++++++++++++++++++++++++++++++" << endl;
    }

    {
        OFstream os
        (
            atomic,
            baseDir/"Test2.txt",
            streamOpt
            // NON_APPEND
        );

        os << "=========================" << endl;

        InfoErr<< "open: " << os.name() << endl;
        InfoErr<< "... sleep" << endl;

        listFiles(baseDir);

        sleep(2);

        os << "+++++++++++++++++++++++++++++++++++" << endl;
    }
    {
        OFstream os
        (
            atomic,
            baseDir/"Test3.txt",
            streamOpt,
            IOstreamOption::APPEND
        );

        os << "=========================" << endl;

        InfoErr<< "open: " << os.name() << endl;
        InfoErr<< "... sleep" << endl;

        listFiles(baseDir);

        sleep(2);

        os << "+++++++++++++++++++++++++++++++++++" << endl;
    }
    {
        OFstream os
        (
            baseDir/"Test4.txt",
            IOstreamOption::ASCII,
            IOstreamOption::COMPRESSED
        );

        os << "=========================" << endl;

        InfoErr<< "open: " << os.name() << endl;
        InfoErr<< "... sleep" << endl;

        listFiles(baseDir);

        sleep(2);

        os << "+++++++++++++++++++++++++++++++++++" << endl;
    }
    {
        OFstream os
        (
            IOstreamOption::ATOMIC,
            baseDir/"Test5.txt"
            // ASCII UNCOMPRESSED NON_APPEND
        );

        os << "=========================" << endl;

        InfoErr<< "open: " << os.name() << endl;
        InfoErr<< "... sleep" << endl;

        listFiles(baseDir);

        sleep(2);

        os << "+++++++++++++++++++++++++++++++++++" << endl;
    }

    Info<< nl << "done:" << endl;

    listFiles(baseDir);

    if (args.found("keep"))
    {
        InfoErr<< "keep: " << baseDir << endl;
    }
    else
    {
        InfoErr<< "rmdir: " << baseDir << endl;
        Foam::rmDir(baseDir);
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
