/*----------------------------------------------------------------------------*\
 =========                   |
 \\      /   F ield          | OpenFOAM: The Open Source CFD Toolbox
  \\    /    O peration      |
   \\  /     A nd            | Copyright (C) 2011 OpenFOAM Foundation
    \\/      M anipulation   | Copyright (C) 2017 OpenCFD Ltd.
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
    dirToString

Description
    Converts a directory path into a camelCase string.
    e.g. dir1/dir2/dir3 becomes dir1Dir2Dir3

Usage
    echo dirName | dirToString

    e.g.
        using sh
        baseDirName=$(echo $dir | $bin/dirToString -strip)

        using csh
        set baseDirName=`echo $dir | $bin/dirToString -strip`

\*----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

/* The executable name (for messages), without requiring access to argv[] */
#define EXENAME  "dirToString"

int main(int argc, char* argv[])
{
    int c;

    if (argc > 1)
    {
        if (!strncmp(argv[1], "-h", 2))
        {
            /* Option: -h, -help */

            fputs
            (
                "\nUsage: " EXENAME
                " [-strip]\n\n"
                "  -strip    ignore leading [./] characters.\n\n"
                "Converts a directory path into a camelCase string\n\n",
                stderr
            );
            return 0;
        }

        if (!strcmp(argv[1], "-s") || !strcmp(argv[1], "-strip"))
        {
            /* Option: -s, -strip */

            while ((c=getchar()) != EOF && (c == '.' || c == '/'))
            {
                /* nop */
            }

            if (c == EOF)
            {
                return 0;
            }

            putchar(c);
        }
    }


    int nextUpper = 0;
    while ((c = getchar()) != EOF)
    {
        if (c == '/')
        {
            nextUpper = 1;
        }
        else if (nextUpper)
        {
            putchar(toupper(c));
            nextUpper = 0;
        }
        else
        {
            putchar(c);
        }
    }

    return 0;
}


/*****************************************************************************/
