/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    testHashing

Description

\*---------------------------------------------------------------------------*/

#include "IOstreams.H"
#include "IOobject.H"
#include "IFstream.H"

#include "stringList.H"
#include "labelList.H"
#include "labelPair.H"

#include "Hashing.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    IFstream is("hashingTests");


    while (is.good())
    {
        const word listType(is);

        Info<< endl;
        IOobject::writeDivider(Info);
        Info<< listType << endl;

        if (listType == "stringList")
        {
            Info<< "contiguous = " << contiguous<string>() << endl << endl;

            stringList lst(is);

            forAll(lst, i)
            {
                unsigned hash1 = Hashing::jenkins(lst[i]);
                unsigned hash2 = string::hash()(lst[i]);

                Info<< hex << hash1
                    << " (prev " << hash2 << ")"
                    << ": " << lst[i] << endl;
            }

        }
        else if (listType == "labelList")
        {
            Info<<"contiguous = " << contiguous<label>() << " "
                << "sizeof(label) " << unsigned(sizeof(label)) << endl << endl;

            labelList lst(is);

            unsigned hash4 = 0;

            forAll(lst, i)
            {
                unsigned hash1 = Hashing::intHash(lst[i]);

                unsigned hash2 = Hashing::jenkins
                (
                    reinterpret_cast<const char*>(&lst[i]),
                    sizeof(label)
                );

                unsigned hash3 = Hashing::jenkins
                (
                    reinterpret_cast<const unsigned*>(&lst[i]),
                    1
                );

                // incremental
                hash4 = Hashing::jenkins
                (
                    reinterpret_cast<const char*>(&lst[i]),
                    sizeof(label),
                    hash4
                );

                Info<< hex << hash1
                    << " (alt: " << hash2 << ")"
                    << " (alt: " << hash3 << ")"
                    << " (incr: " << hash4 << ")"
                    << ": " << dec << lst[i] << endl;
            }

            if (contiguous<label>())
            {
                unsigned hash1 = Hashing::jenkins
                (
                    lst.cdata(),
                    lst.size() * sizeof(label)
                );

                Info<<"contiguous hashed value " << hex << hash1 << endl;
            }
        }
        else if (listType == "labelListList")
        {
            List< List<label> > lst(is);

            forAll(lst, i)
            {
                unsigned hash1 = Hashing::jenkins
                (
                    reinterpret_cast<const char*>(lst[i].cdata()),
                    lst[i].size() * sizeof(label)
                );

                Info<< hex << hash1
                    << ": " << dec << lst[i] << endl;
            }

        }
    }

    return 0;
}


// ************************************************************************* //
