/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2021 OpenCFD Ltd.
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

InClass
    Foam::StringStream

Description
    Input/output from string buffers.

SourceFiles
    StringStream.C

\*---------------------------------------------------------------------------*/

#ifndef StringStream_H
#define StringStream_H

#include "ISstream.H"
#include "OSstream.H"
#include <sstream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace Detail
{

/*---------------------------------------------------------------------------*\
                Class Detail::StringStreamAllocator Declaration
\*---------------------------------------------------------------------------*/

//- Allocator for variants of a std stringstream
template<class StreamType>
class StringStreamAllocator
{
protected:

    // Protected Member Data

        //- The stream type
        typedef StreamType stream_type;

        //- The input/output stream.
        stream_type stream_;


    // Constructors

        //- Default construct
        StringStreamAllocator() = default;

        //- Copy construct from string
        StringStreamAllocator(const std::string& s)
        :
            stream_(s)
        {}


public:

    // Member Functions

        //- Get the string - as Foam::string rather than std::string
        Foam::string str() const
        {
            return Foam::string(stream_.str());
        }

        //- Set the string
        void str(const std::string& s)
        {
            stream_.str(s);
        }
};

} // End namespace Detail


/*---------------------------------------------------------------------------*\
                        Class IStringStream Declaration
\*---------------------------------------------------------------------------*/

//- Input from string buffer, using a ISstream. Always UNCOMPRESSED.
class IStringStream
:
    public Detail::StringStreamAllocator<std::istringstream>,
    public ISstream
{
    typedef Detail::StringStreamAllocator<std::istringstream> allocator_type;

public:

    // Constructors

        //- Default construct or with specified stream option
        explicit IStringStream
        (
            IOstreamOption streamOpt = IOstreamOption()
        )
        :
            allocator_type(),
            ISstream(stream_, "input", streamOpt.format(), streamOpt.version())
        {}

        //- Construct from std::string
        explicit IStringStream
        (
            const std::string& s,
            IOstreamOption streamOpt = IOstreamOption()
        )
        :
            allocator_type(s),
            ISstream(stream_, "input", streamOpt.format(), streamOpt.version())
        {}

        //- Construct from char*
        explicit IStringStream
        (
            const char* s,
            IOstreamOption streamOpt = IOstreamOption()
        )
        :
            allocator_type(s),
            ISstream(stream_, "input", streamOpt.format(), streamOpt.version())
        {}

        //- Copy construct, copies content and format
        IStringStream(const IStringStream& str)
        :
            allocator_type(str.str()),
            ISstream(stream_, str.name(), str.format(), str.version())
        {}


    // Member Functions

        //- Reset the input buffer and rewind the stream
        virtual void reset(const std::string& s)
        {
            this->str(s);
            this->rewind();
        }

        //- Print stream description to Ostream
        virtual void print(Ostream& os) const;


    // Member Operators

        //- Return a non-const reference to const Istream
        //  Needed for read-constructors where the stream argument is temporary.
        Istream& operator()() const
        {
            // Could also rewind
            return const_cast<IStringStream&>(*this);
        }


    // Additional constructors and methods (as per v2012 and earlier)
    #ifdef Foam_IOstream_extras

        //- Default construct
        explicit IStringStream
        (
            IOstreamOption::streamFormat fmt,
            IOstreamOption::versionNumber ver = IOstreamOption::currentVersion
        )
        :
            IStringStream(IOstreamOption(fmt, ver))
        {}

        //- Construct from std::string
        IStringStream
        (
            const std::string& s,
            IOstreamOption::streamFormat fmt,
            IOstreamOption::versionNumber ver = IOstreamOption::currentVersion
        )
        :
            IStringStream(s, IOstreamOption(fmt, ver))
        {}

        //- Construct from char*
        IStringStream
        (
            const char* s,
            IOstreamOption::streamFormat fmt,
            IOstreamOption::versionNumber ver = IOstreamOption::currentVersion
        )
        :
            IStringStream(s, IOstreamOption(fmt, ver))
        {}

    #endif /* Foam_IOstream_extras */
};


/*---------------------------------------------------------------------------*\
                        Class OStringStream Declaration
\*---------------------------------------------------------------------------*/

//- Output to string buffer, using a OSstream. Always UNCOMPRESSED.
class OStringStream
:
    public Detail::StringStreamAllocator<std::ostringstream>,
    public OSstream
{
    typedef Detail::StringStreamAllocator<std::ostringstream> allocator_type;

public:

    // Constructors

        //- Default construct or with specified stream option
        explicit OStringStream
        (
            IOstreamOption streamOpt = IOstreamOption()
        )
        :
            allocator_type(),
            OSstream(stream_, "output", streamOpt.format(), streamOpt.version())
        {}

        //- Copy construct, copies content and format
        OStringStream(const OStringStream& str)
        :
            allocator_type(str.str()),
            OSstream(stream_, str.name(), str.format(), str.version())
        {}


    // Member Functions

        //- Reset the output buffer and rewind the stream
        void reset()
        {
            this->str("");  // No other way to reset the end
            this->rewind();
        }

        //- Rewind the output stream
        virtual void rewind()
        {
            // pubseekpos() instead of seekp() for symmetry with other classes
            stream_.rdbuf()->pubseekpos(0, std::ios_base::out);
        }

        //- Print stream description to Ostream
        virtual void print(Ostream& os) const;


    // Older style, without stream option (including 2012 release)
    #ifdef Foam_IOstream_extras

        //- Default construct
        explicit OStringStream
        (
            IOstreamOption::streamFormat fmt,
            IOstreamOption::versionNumber ver = IOstreamOption::currentVersion
        )
        :
            OStringStream(IOstreamOption(fmt, ver))
        {}

    #endif /* Foam_IOstream_extras */
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
