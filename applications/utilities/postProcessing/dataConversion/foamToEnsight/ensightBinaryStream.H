/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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

Class
    Foam::ensightBinaryStream

Description

SourceFiles
    ensightBinaryStream.C

\*---------------------------------------------------------------------------*/

#ifndef ensightBinaryStream_H
#define ensightBinaryStream_H

#include "ensightStream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class ensightBinaryStream Declaration
\*---------------------------------------------------------------------------*/

class ensightBinaryStream
:
    public ensightStream
{
    // Private data

        //- Output file stream
        autoPtr<std::ofstream> str_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        ensightBinaryStream(const ensightBinaryStream&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const ensightBinaryStream&) = delete;


public:

    // Constructors

        //- Construct from components
        ensightBinaryStream(const fileName& f)
        :
            ensightStream(f),
            str_
            (
                new std::ofstream
                (
                    f.c_str(),
                    ios_base::out | ios_base::binary | ios_base::trunc
                )
            )
        {}


    //- Destructor
    virtual ~ensightBinaryStream()
    {}


    // Member Functions

        virtual void write(const char* val)
        {
            char buffer[80];
            strncpy(buffer, val, 80);
            str_().write(buffer, 80*sizeof(char));
        }

        virtual void write(const int val)
        {
            str_().write(reinterpret_cast<const char*>(&val), sizeof(int));
        }

        virtual void write(const scalarField& sf)
        {
            if (sf.size())
            {
                List<float> temp(sf.size());

                forAll(sf, i)
                {
                    temp[i] = float(sf[i]);
                }

                str_().write
                (
                    reinterpret_cast<const char*>(temp.begin()),
                    sf.size()*sizeof(float)
                );
            }
        }

        virtual void write(const List<int>& sf)
        {
            str_().write
            (
                reinterpret_cast<const char*>(sf.begin()),
                sf.size()*sizeof(int)
            );
        }

        virtual void writePartHeader(const label partI)
        {
            write("part");
            write(partI);
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
