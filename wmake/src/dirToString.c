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
    converts a directory path into a string with appropriate capitalisation
    e.g. dir1/dir2 becomes dir1Dir2

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

int main(int argc, char* argv[])
{
    int c;
    int nextupper = 0;

    if (argc > 1)
    {
        if (!strncmp(argv[1], "-h", 2))  /* -h, -help */
        {
            fprintf
            (
                stderr,
                "\nUsage: %s [-strip]\n\n",
                "dirToString"
            );

            fprintf
            (
                stderr,
                "  -strip    ignore leading [./] characters.\n\n"
                "Transform dir1/dir2 to camel-case dir1Dir2\n\n"
            );

            return 0;
        }

        if (!strcmp(argv[1], "-s") || !strcmp(argv[1], "-strip"))  /* -s, -strip */
        {
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


    while ((c=getchar()) != EOF)
    {
        if (c == '/')
        {
            nextupper = 1;
        }
        else
        {
            if (nextupper)
            {
                putchar(toupper(c));
                nextupper = 0;
            }
            else
            {
                putchar(c);
            }
        }
    }

    return 0;
}


/*****************************************************************************/
