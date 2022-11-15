/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2022 OpenCFD Ltd.
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
    Foam::ifstreamPointer

Description
    A wrapped \c std::ifstream with possible compression handling
    (igzstream) that behaves much like a \c std::unique_ptr.

Note
    No <tt>operator bool</tt> to avoid inheritance ambiguity with
    <tt>std::ios::operator bool</tt>.

SourceFiles
    fstreamPointers.C

Class
    Foam::ofstreamPointer

Description
    A wrapped \c std::ofstream with possible compression handling
    (ogzstream) that behaves much like a \c std::unique_ptr.

Note
    No <tt>operator bool</tt> to avoid inheritance ambiguity with
    <tt>std::ios::operator bool</tt>.

SourceFiles
    fstreamPointers.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_fstreamPointer_H
#define Foam_fstreamPointer_H

#include "IOstreamOption.H"
#include "fileName.H"
#include <fstream>
#include <memory>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class ifstreamPointer Declaration
\*---------------------------------------------------------------------------*/

class ifstreamPointer
{
    // Private Data

        //- The stream pointer (ifstream or igzstream)
        std::unique_ptr<std::istream> ptr_;


protected:

    // Protected Member Functions

        //- Special 'rewind' method for compressed stream
        void reopen_gz(const std::string& pathname);


public:

    // Generated Methods

        //- Default construct (empty)
        ifstreamPointer() noexcept = default;

        //- No copy construct
        ifstreamPointer(const ifstreamPointer&) = delete;

        //- Move construct
        ifstreamPointer(ifstreamPointer&&) = default;

        //- No copy assignment
        void operator=(const ifstreamPointer&) = delete;

        //- Move assignment
        ifstreamPointer& operator=(ifstreamPointer&&) = default;

        //- Destructor
        ~ifstreamPointer() = default;


    // Constructors

        //- Construct from pathname.
        //  Attempts to read the specified file.
        //  If that fails, try as a compressed file (.gz ending).
        explicit ifstreamPointer(const fileName& pathname);


    // Member Functions

        //- True if compiled with libz support
        static bool supports_gz();


    // Access

        //- The stream pointer (ifstream or igzstream)
        std::istream* get() noexcept { return ptr_.get(); }

        //- The stream pointer (ifstream or igzstream)
        const std::istream* get() const noexcept { return ptr_.get(); }

        //- Which compression type?
        IOstreamOption::compressionType whichCompression() const;


    // Edit

        //- Return managed pointer and release ownership
        std::istream* release() noexcept { return ptr_.release(); }

        //- Replace the managed pointer
        void reset(std::istream* ptr) noexcept { ptr_.reset(ptr); }


    // Operators

        //- Reference to the stream (no nullptr checking)
        std::istream& operator*() { return *ptr_; }

        //- Const-reference to the stream (no nullptr checking)
        const std::istream& operator*() const { return *ptr_; }

        //- Pointer dereference
        std::istream* operator->() noexcept { return ptr_.get(); }

        //- Pointer dereference
        const std::istream* operator->() const noexcept { return ptr_.get(); }
};


/*---------------------------------------------------------------------------*\
                      Class ofstreamPointer Declaration
\*---------------------------------------------------------------------------*/

class ofstreamPointer
{
    // Private Data

        //- The stream pointer (ofstream | ogzstream | ocountstream)
        std::unique_ptr<std::ostream> ptr_;

        //- Atomic file creation
        bool atomic_;


protected:

    // Protected Member Functions

        //- Reopen for compressed/non-compressed
        void reopen(const std::string& pathname);

        //- Close stream and rename file
        void close(const std::string& pathname);


public:

    // Generated Methods

        //- No copy construct
        ofstreamPointer(const ofstreamPointer&) = delete;

        //- Move construct
        ofstreamPointer(ofstreamPointer&&) = default;

        //- No copy assignment
        void operator=(const ofstreamPointer&) = delete;

        //- Move assignment
        ofstreamPointer& operator=(ofstreamPointer&&) = default;

        //- Destructor
        ~ofstreamPointer() = default;


    // Constructors

        //- Default construct (empty)
        ofstreamPointer() noexcept;

        //- Construct as null output stream using Foam::ocountstream
        explicit ofstreamPointer(std::nullptr_t);

        //- Construct from pathname, compression, append, file handling atomic
        //  \param pathname The file name to open for writing
        //  \param comp     UNCOMPRESSED | COMPRESSED
        //  \param append   Open in append mode
        //  \param atomic   Write into temporary file (not target file).
        //      This option should only be used with a stream wrapper
        //      (eg, OFstream) that handles the final renaming.
        explicit ofstreamPointer
        (
            const fileName& pathname,
            IOstreamOption::compressionType comp = IOstreamOption::UNCOMPRESSED,
            const bool append = false,
            const bool atomic = false
        );


    // Member Functions

        //- True if compiled with libz support
        static bool supports_gz();


    // Access

        //- The stream pointer (ofstream or ogzstream)
        std::ostream* get() noexcept { return ptr_.get(); }

        //- The stream pointer (ofstream or ogzstream)
        const std::ostream* get() const noexcept { return ptr_.get(); }

        //- Which compression type?
        IOstreamOption::compressionType whichCompression() const;


    // Edit

        //- Return managed pointer and release ownership
        std::ostream* release() noexcept { return ptr_.release(); }

        //- Replace the managed pointer
        void reset(std::ostream* ptr) noexcept { ptr_.reset(ptr); }


    // Operators

        //- Reference to the stream (no nullptr checking)
        std::ostream& operator*() { return *ptr_; }

        //- Const-reference to the stream (no nullptr checking)
        const std::ostream& operator*() const { return *ptr_; }

        //- Pointer dereference
        std::ostream* operator->() noexcept { return ptr_.get(); }

        //- Pointer dereference
        const std::ostream* operator->() const noexcept { return ptr_.get(); }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
