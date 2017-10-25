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
    Foam::UOBufStream

Description
    Similar to OStringStream but using an externally managed buffer for
    its output.

    This allows the output buffer to be reused and can make it easier when
    writing out data.  It is the user's responsibility to ensure proper
    synchronization in the sizes. Provided that the external buffer is large
    enough that overflow does not occur, the following usage pattern
    works.

    \code
        DynamicList<char> buffer(4096);     // allocate some large buffer

        {
            UOBufStream os(buffer);
            os << "content1" << " and more content";
            buffer.setSize(os.size());      // synchronize sizes
        }

        something.write(buffer, buffer.size());
    \endcode

    Although the UOBufStream is quite lightweight, there may be cases
    where it is preferable to reuse the stream as well.
    \code
        DynamicList<char> buffer(4096);     // allocate some large buffer

        UOBufStream os(buffer);
        os << "content1" << " and more content";
        buffer.setSize(os.size());          // synchronize sizes

        something.write(buffer, buffer.size());

        os.rewind();
        os << "content2";
        buffer.setSize(os.size());      // synchronize sizes

        something.write(buffer, buffer.size());

        // or simply using the output size directly (without sync)
        os.rewind();
        os << "content3";

        something.write(buffer, os.size());
    \endcode

\*---------------------------------------------------------------------------*/

#ifndef UOBufStream_H
#define UOBufStream_H

#include "BufStreamAllocator.H"
#include "OSstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration
template<class T, int SizeMin> class DynamicList;


/*---------------------------------------------------------------------------*\
                         Class UOBufStream Declaration
\*---------------------------------------------------------------------------*/

//- An OSstream attached to an unallocated external buffer
class UOBufStream
:
    public BufStreamAllocator<std::ostream,false>,
    public OSstream
{
    typedef BufStreamAllocator<std::ostream,false> allocator_type;

public:

    // Constructors

        //- Construct using specified buffer and number of bytes
        UOBufStream
        (
            char* buffer,
            size_t nbytes,
            streamFormat format=ASCII,
            versionNumber version=currentVersion
        )
        :
            allocator_type(buffer, nbytes),
            OSstream(stream_, "output", format,version)
        {}


        //- Construct using data area from a List and number of bytes
        UOBufStream
        (
            UList<char>& buffer,
            size_t size,
            streamFormat format=ASCII,
            versionNumber version=currentVersion
        )
        :
            UOBufStream(buffer.data(), size, format,version)
        {}


        //- Construct using data area from a List and its inherent storage size
        UOBufStream
        (
            UList<char>& buffer,
            streamFormat format=ASCII,
            versionNumber version=currentVersion
        )
        :
            UOBufStream(buffer.data(), buffer.size(), format,version)
        {}


        //- Construct using data area from a DynamicList and its capacity
        template<int SizeMin>
        UOBufStream
        (
            DynamicList<char,SizeMin>& buffer,
            streamFormat format=ASCII,
            versionNumber version=currentVersion
        )
        :
            UOBufStream(buffer.data(), buffer.capacity(), format,version)
        {}


    //- Destructor
    ~UOBufStream()
    {}


    // Member functions

        //- Return the current output position in the buffer
        std::streampos size() const
        {
            return allocator_type::tellp();
        }


        //- Rewind the stream, clearing any old errors
        void rewind()
        {
            allocator_type::rewind();
            setGood();  // resynchronize with internal state
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
