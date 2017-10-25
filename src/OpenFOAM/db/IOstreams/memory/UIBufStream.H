/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 OpenCFD Ltd.
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

Class
    Foam::UIBufStream

Description
    Similar to IStringStream but using an externally managed buffer for its
    input. This allows the input buffer to be filled (and refilled) from
    various sources.

    Note that this stream will normally be used as a "one-shot" reader.
    Caution must be exercised that the referenced buffer remains valid and
    without any intermediate resizing for the duration of the stream's use.

    An example of possible use:
    \code
        DynamicList<char> buffer(4096);   // allocate some large buffer

        nread = something.read(buffer.data(),1024); // fill with content
        buffer.setSize(nread);            // content size

        // construct dictionary, or something else
        UIBufStream is(buffer)
        dictionary dict1(is);

        // sometime later
        nread = something.read(buffer.data(),2048); // fill with content
        buffer.setSize(nread);              // content size

        // without intermediate variable
        dictionary dict2(UIBufStream(buffer)());
    \endcode

\*---------------------------------------------------------------------------*/

#ifndef UIBufStream_H
#define UIBufStream_H

#include "BufStreamAllocator.H"
#include "ISstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class UIBufStream Declaration
\*---------------------------------------------------------------------------*/

class UIBufStream
:
    public BufStreamAllocator<std::istream,false>,
    public ISstream
{
    typedef BufStreamAllocator<std::istream,false> allocator_type;

public:

    // Constructors

        //- Construct using specified buffer and number of bytes
        UIBufStream
        (
            const char* buffer,
            size_t nbytes,
            streamFormat format=ASCII,
            versionNumber version=currentVersion,
            const Foam::string& name="input"
        )
        :
            allocator_type(const_cast<char*>(buffer), nbytes),
            ISstream(stream_, name, format, version)
        {}


        //- Construct using data area from a List and number of bytes
        UIBufStream
        (
            const UList<char>& buffer,
            label size,
            streamFormat format=ASCII,
            versionNumber version=currentVersion,
            const Foam::string& name="input"
        )
        :
            UIBufStream(buffer.cdata(), size, format,version,name)
        {}


        //- Construct using data area from a List and its inherent storage size
        //  Uses addressed size, thus no special treatment for a DynamicList
        UIBufStream
        (
            const UList<char>& buffer,
            streamFormat format=ASCII,
            versionNumber version=currentVersion,
            const Foam::string& name="input"
        )
        :
            UIBufStream(buffer.cdata(), buffer.size(), format,version,name)
        {}


    //- Destructor
    ~UIBufStream()
    {}


    // Member functions

        //- Return the current get position in the buffer
        std::streampos pos() const
        {
            return allocator_type::tellg();
        }


        //- Rewind the stream, clearing any old errors
        virtual void rewind()
        {
            allocator_type::rewind();
            setGood();  // resynchronize with internal state
        }


    // Member operators

        //- A non-const reference to const Istream
        //  Needed for read-constructors where the stream argument is temporary
        Istream& operator()() const
        {
            return const_cast<Istream&>(static_cast<const Istream&>(*this));
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
