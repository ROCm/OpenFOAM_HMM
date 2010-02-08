/*---------------------------------------------------------------------------*\
 =========                   |
 \\      /   F ield          | OpenFOAM: The Open Source CFD Toolbox
  \\    /    O peration      |
   \\  /     A nd            | Copyright (C) 2010-2010 OpenCFD Ltd.
    \\/      M anipulation   |
------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    wmkdepend

Description
    A fast dependency list generator that emulates the behaviour and
    output of cpp -M. However, the output contains no duplications and
    is ~40% faster than cpp.

    The algorithm uses flex to scan for includes and searches the files
    found.  Each file is entered into a hash table so that files are scanned
    only once.  This is why this program is faster than cpp.

Usage
    wmkdep [ -Idirectory ... -Idirectory ] filename

\*---------------------------------------------------------------------------*/

#include <cstdio>
#include <stdlib.h>
#include <cstring>

#include "wmkdependParser.h"

// Note: since we use the Coco/R default error messages, we must use
// wide streams for stderr.

void printUsage(const char* message = NULL)
{
    if (message)
    {
        fwprintf(stderr, L"\nError: %s\n\n", message);
    }

    fwprintf
    (
        stderr,
        L"Usage: wmkdepend [ -Idirectory ... -Idirectory ] filename\n"
    );
}


int main(int argc, char* argv[])
{
    if (argc == 1)
    {
        printUsage("Error: input file not supplied");
        ::exit(1);
    }

    for (int i=1; i < argc; i++)
    {
        if (strncmp(argv[i], "-I", 2) == 0 && strlen(argv[i]) > 2)
        {
            std::string dirName(argv[i] + 2);

            // add trailing slash if required
            if (dirName.rfind('/') != dirName.size()-1)
            {
                dirName += '/';
            }

            wmake::Parser::includeDirs.push_back(dirName);
        }
    }

    std::string sourceFile(argv[argc-1]);

    fwprintf
    (
        stderr,
        L"Making dependency list for source file %s\n",
        sourceFile.c_str()
    );

    std::string::size_type basePos = sourceFile.rfind('/');
    if (basePos == std::string::npos)
    {
        basePos = 0;
    }
    else
    {
        basePos++;
    }

    std::string::size_type dotPos = sourceFile.rfind('.');
    if
    (
        dotPos == std::string::npos
     || dotPos == sourceFile.size()-1
     || dotPos <= basePos
    )
    {
        fwprintf
        (
            stderr,
            L"Cannot find extension in source file name %s\n",
            sourceFile.c_str()
        );
        ::exit(1);
    }

    std::string depFile = sourceFile.substr(0, dotPos);
    depFile += ".dep";

    const std::string sourceExt = sourceFile.substr(dotPos);
    if (sourceExt == ".java")
    {
        // import directories to ignore
        wmake::Parser::ignoreDir("java.*");
        wmake::Parser::ignoreDir("org.*");
        wmake::Parser::ignoreDir("com.*");
        wmake::Parser::ignoreDir("sunw.*");
        wmake::Parser::ignoreDir("sun.*");
        wmake::Parser::ignoreDir("launcher.*");

        std::cout
            << "$(CLASSES_DIR)/"
            << sourceFile.substr(basePos, dotPos - basePos) << ".class: "
            << depFile << "\n";
    }
    else
    {
        std::cout
            << "$(OBJECTS_DIR)/"
            << sourceFile.substr(basePos, dotPos - basePos) << ".o: "
            << depFile << "\n";
    }


    wmake::Parser::sourceFile = sourceFile;
    wmake::Parser::depFile = depFile;

    wmake::Parser::includeFile(sourceFile);

    return 0;
}



/*****************************************************************************/
