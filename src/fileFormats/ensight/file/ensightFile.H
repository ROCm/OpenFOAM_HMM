/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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
    Foam::ensightFile

Description
    Ensight output with specialized write() for strings, integers and floats.
    Correctly handles binary write as well.

\*---------------------------------------------------------------------------*/

#ifndef ensightFile_H
#define ensightFile_H

#include "OFstream.H"
#include "IOstream.H"

#include "ensightFileName.H"
#include "ensightVarName.H"
#include "UList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class ensightFile Declaration
\*---------------------------------------------------------------------------*/

class ensightFile
:
    public OFstream
{
    // Private data

        //- Allow undef in results
        static bool allowUndef_;

        //- Value to represent undef in results
        static scalar undefValue_;

        //- The '*' mask appropriate for subDir
        static string mask_;

        //- The printf format for zero-padded subdirectory numbers
        static string dirFmt_;


    // Private Member Functions

        //- Initialize by setting the ASCII output formatting
        void initialize();

        //- Disallow default bitwise assignment
        void operator=(const ensightFile&) = delete;

        //- Disallow default copy constructor
        ensightFile(const ensightFile&) = delete;

public:

    // Static Member Functions

        //- Return a null ensightFile
        inline static const ensightFile& null()
        {
            return NullObjectRef<ensightFile>();
        }


    // Constructors

        //- Construct from pathname.
        //  The entire pathname is checked for valid ensight naming.
        ensightFile
        (
            const fileName& pathname,
            IOstream::streamFormat format=IOstream::BINARY
        );

        //- Construct from path and name.
        //  Only the name portion is checked for valid ensight naming.
        ensightFile
        (
            const fileName& path,
            const fileName& name,
            IOstream::streamFormat format=IOstream::BINARY
        );


    //- Destructor
    ~ensightFile();


    // Access

        //- Return setting for whether 'undef' values are allowed in results
        static bool allowUndef();

        //- The '*' mask appropriate for subDir
        static string mask();

        //- Consistent zero-padded numbers for subdirectories
        static string subDir(const label);

        //- Set width of subDir and mask. Default width is 8 digits.
        //  Max width is 31 digits.
        static void subDirWidth(const label);

        //- Return current width of subDir and mask.
        static label subDirWidth();


    // Edit

        static bool allowUndef(bool);

        //- Assign the value to represent undef in the results
        //  Returns the previous value
        //  NB: do not use values larger than floatScalarVGREAT
        static scalar undefValue(const scalar);


    // Output

        //- Inherit write from Ostream
        using Ostream::write;

        //- Binary write
        virtual Ostream& write(const char*, std::streamsize count);

        //- Write element keyword with trailing newline, optionally with undef
        virtual Ostream& writeKeyword(const keyType&);

        //- Write "C Binary" for binary files (eg, geometry/measured)
        Ostream& writeBinaryHeader();

        //- Write undef value
        Ostream& writeUndef();

        //- Write C-string as "%79s" or as binary (max 80 chars)
        Ostream& write(const char*);

        //- Write string as "%79s" or as binary (max 80 chars)
        Ostream& write(const string&);

        //- Write integer as "%10d" or as binary
        Ostream& write(const label);

        //- Write integer with specified width or as binary
        Ostream& write(const label, const label fieldWidth);

        //- Write float as "%12.5e" or as binary
        Ostream& write(const scalar);

        //- Add carriage return to ascii stream
        void newline();


    // Convenience Output Methods

        //- Begin a part (0-based index).
        void beginPart(const label index);

        //- Begin a "particle coordinates" block (measured data)
        void beginParticleCoordinates(const label nparticles);

        //- Write a list of floats as "%12.5e" or as binary
        //  With carriage return after each value (ascii stream)
        void writeList(const UList<scalar>& field);

        //- Write an indirect list of scalars as "%12.5e" or as binary
        //  With carriage return after each value (ascii stream)
        void writeList
        (
            const UList<scalar>& field,
            const labelUList& idList
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
