/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
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
#include "ensightFileName.H"
#include "ensightVarName.H"
#include "IndirectListBase.H"

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
    // Private Data

        //- Allow undef in results
        static bool allowUndef_;

        //- Value to represent undef in results (default: 1e+37, floatVGREAT)
        static scalar undefValue_;

        //- The '*' mask appropriate for subDir
        static string mask_;

        //- The printf format for zero-padded subdirectory numbers
        static string dirFmt_;


    // Private Member Functions

        //- Initialize by setting the ASCII output formatting
        void initialize();

        //- No copy construct
        ensightFile(const ensightFile&) = delete;

        //- No copy assignment
        void operator=(const ensightFile&) = delete;


public:

    // Static Data Members

        //- The keyword "coordinates"
        static const char* const coordinates;


    // Static Functions

        //- Return a null ensightFile
        inline static const ensightFile& null()
        {
            return NullObjectRef<ensightFile>();
        }


    // Constructors

        //- Construct from pathName.
        //  The entire pathName is checked for valid ensight naming.
        explicit ensightFile
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
    ~ensightFile() = default;


    // Member Functions

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

        //- Enable/disable use of \c undef keyword and value
        static bool allowUndef(bool enabled);

        //- Assign the value to represent undef in the results
        //  Returns the previous value
        //  NB: do not use values larger than floatScalarVGREAT
        static scalar undefValue(const scalar value);


    // Output

        //- Write element keyword with trailing newline,
        //- optionally with undef and the value for undefined
        virtual Ostream& writeKeyword(const keyType& key);

        //- Write "C Binary" string for binary files (eg, geometry/measured)
        Ostream& writeBinaryHeader();

        //- Write C-string as "%79s" or as binary (max 80 chars)
        Ostream& writeString(const char* str);

        //- Write string as "%79s" or as binary (max 80 chars)
        Ostream& writeString(const std::string& str);

        //- Write undef value
        Ostream& writeUndef();


        //- Writing token does not make sense
        virtual bool write(const token&)
        {
            NotImplemented;
            return false;
        }

        //- Writing single character does not make sense
        virtual Ostream& write(const char)
        {
            NotImplemented;
            return *this;
        }

        //- Binary write
        virtual Ostream& write(const char* buf, std::streamsize count);

        //- Write C-string, uses writeString()
        virtual Ostream& write(const char* str);

        //- Write word, uses writeString()
        virtual Ostream& write(const word& str);

        //- Write string, uses writeString()
        virtual Ostream& write(const string& str);

        //- Write integer as "%10d" or as binary
        virtual Ostream& write(const int32_t val);

        //- Write integer as "%10d" or as binary
        virtual Ostream& write(const int64_t val);

        //- Write integer with specified width or as binary
        Ostream& write(const label value, const label fieldWidth);

        //- Write floating-point as "%12.5e" or as binary
        virtual Ostream& write(const floatScalar val);

        //- Write floating-point as "%12.5e" or as binary
        virtual Ostream& write(const doubleScalar val);

        //- Add carriage return to ascii stream
        void newline();


    // Convenience Output Methods

        //- Begin a part (0-based index internally).
        void beginPart(const label index);

        //- Begin a "particle coordinates" block (measured data)
        void beginParticleCoordinates(const label nparticles);

        //- Write a list of integers
        //  With carriage return after each value (ascii stream)
        void writeLabels(const UList<label>& list);

        //- Write a list of integers
        //  With carriage return after each value (ascii stream)
        template<class Addr>
        void writeLabels(const IndirectListBase<label, Addr>& list);

        //- Write a list of integers as float values
        //  With carriage return after each value (ascii stream)
        void writeList(const UList<label>& field);

        //- Write a list of floats as "%12.5e" or as binary
        //  With carriage return after each value (ascii stream)
        void writeList(const UList<scalar>& field);

        //- Write an indirect list of scalars as "%12.5e" or as binary
        //  With carriage return after each value (ascii stream)
        template<class Addr>
        void writeList(const IndirectListBase<scalar, Addr>& field);


    // Other Methods

        //- Check for any NaN in the field
        static bool isUndef(const UList<scalar>& field);

        //- Check for any NaN in the field
        template<class Addr>
        static bool isUndef(const IndirectListBase<scalar, Addr>& field);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ensightFileTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
