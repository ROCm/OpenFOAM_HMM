/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2018 OpenCFD Ltd.
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
    Foam::UIListStream

Description
    Similar to IStringStream but using an externally managed buffer for its
    input. This allows the input buffer to be filled (and refilled) from
    various sources.

    Note that this stream will normally be used as a "one-shot" reader.
    Caution must be exercised that the referenced buffer remains valid and
    without any intermediate resizing for the duration of the stream's use.

    An example of possible use:
    \code
        DynamicList<char> buffer(4096);     // allocate some large buffer

        nread = something.read(buffer.data(),1024); // fill with content
        buffer.resize(nread);               // content size

        // construct dictionary, or something else
        UIListStream is(buffer)
        dictionary dict1(is);

        // sometime later
        nread = something.read(buffer.data(),2048); // fill with content
        buffer.resize(nread);               // content size

        // without intermediate variable
        dictionary dict2(UIListStream(buffer)());
    \endcode

See Also
    Foam::IListStream
    Foam::OListStream
    Foam::UOListStream

\*---------------------------------------------------------------------------*/

#ifndef UIListStream_H
#define UIListStream_H

#include "FixedList.H"
#include "UList.H"
#include "ISstream.H"
#include "memoryStreamBuffer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class uiliststream Declaration
\*---------------------------------------------------------------------------*/

//- Similar to std::istringstream, but with an externally managed input buffer.
//  This allows the input buffer to be filled or refilled from various sources
//  without copying.
class uiliststream
:
    virtual public std::ios,
    protected memorybuf::in,
    public std::istream
{
public:

    //- Construct for character array and number of bytes
    uiliststream(const char* buffer, size_t nbytes)
    :
        memorybuf::in(const_cast<char*>(buffer), nbytes),
        std::istream(static_cast<memorybuf::in*>(this))
    {}

    //- Reset buffer pointers
    inline void reset(char *buffer, size_t nbytes)
    {
        resetg(buffer, nbytes);
    }

    //- Rewind the stream, clearing any old errors
    void rewind()
    {
        this->pubseekpos(0, std::ios_base::in);
        clear(); // for safety, clear any old errors
    }
};


namespace Detail
{

/*---------------------------------------------------------------------------*\
                Class Detail::UIListStreamAllocator Declaration
\*---------------------------------------------------------------------------*/

//- An stream/stream-buffer input allocator for a externally allocated list
class UIListStreamAllocator
{
protected:

    // Protected data

        typedef std::istream stream_type;

        //- The stream buffer
        memorybuf::in buf_;

        //- The stream
        stream_type stream_;


    // Constructors

        //- Construct for character array and number of bytes
        UIListStreamAllocator(char *buffer, size_t nbytes)
        :
            buf_(buffer, nbytes),
            stream_(&buf_)
        {}


    // Protected Member Functions

        //- Reset buffer pointers
        inline void reset(char *buffer, size_t nbytes)
        {
            buf_.resetg(buffer, nbytes);
        }

        void printBufInfo(Ostream& os) const
        {
            buf_.printBufInfo(os);
        }

public:

    // Public Member Functions

        //- Const UList access to the input characters (shallow copy).
        inline const UList<char> list() const
        {
            return buf_.list();
        }

        //- Non-const UList access to the input characters (shallow copy).
        inline UList<char> list()
        {
            return buf_.list();
        }

        //- The list size
        inline label size() const
        {
            return buf_.capacity();
        }

        //- Position of the get buffer
        std::streampos tellg() const
        {
            return buf_.tellg();
        }

        //- Move to buffer start, clear errors
        void rewind()
        {
            buf_.pubseekpos(0, std::ios_base::in);
            stream_.clear(); // for safety, clear any old errors
        }
};

} // End namespace Detail


/*---------------------------------------------------------------------------*\
                         Class UIListStream Declaration
\*---------------------------------------------------------------------------*/

class UIListStream
:
    public Detail::UIListStreamAllocator,
    public ISstream
{
    typedef Detail::UIListStreamAllocator allocator_type;

public:

    // Constructors

        //- Construct using specified buffer and number of bytes
        UIListStream
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


        //- Construct using data area from a FixedList
        template<unsigned N>
        UIListStream
        (
            const FixedList<char, N>& buffer,
            streamFormat format=ASCII,
            versionNumber version=currentVersion,
            const Foam::string& name="input"
        )
        :
            UIListStream(buffer.cdata(), N, format, version, name)
        {}

        //- Construct using data area from a List and number of bytes
        UIListStream
        (
            const UList<char>& buffer,
            label size,
            streamFormat format=ASCII,
            versionNumber version=currentVersion,
            const Foam::string& name="input"
        )
        :
            UIListStream(buffer.cdata(), size, format, version, name)
        {}


        //- Construct using data area from a List and its inherent storage size
        //  Uses addressed size, thus no special treatment for a DynamicList
        UIListStream
        (
            const UList<char>& buffer,
            streamFormat format=ASCII,
            versionNumber version=currentVersion,
            const Foam::string& name="input"
        )
        :
            UIListStream(buffer.cdata(), buffer.size(), format, version, name)
        {}


    // Member Functions

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

        //- Print description to Ostream
        virtual void print(Ostream& os) const;


    // Member Operators

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
