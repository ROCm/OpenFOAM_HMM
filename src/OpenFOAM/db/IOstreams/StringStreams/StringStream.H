/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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

/*---------------------------------------------------------------------------*\
                    Class StringStreamAllocator Declaration
\*---------------------------------------------------------------------------*/

//- Allocator for variants of a std stringstream
template<class StreamType>
class StringStreamAllocator
{
protected:

    // Member Data

        //- The stream type
        typedef StreamType stream_type;

        //- The input/output stream.
        stream_type stream_;


    // Constructors

        //- Construct null
        StringStreamAllocator()
        :
            stream_()
        {}


        //- Copy construct from string
        StringStreamAllocator(const std::string& buffer)
        :
            stream_(buffer)
        {}


public:

    // Public Member Functions

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


/*---------------------------------------------------------------------------*\
                        Class IStringStream Declaration
\*---------------------------------------------------------------------------*/

//- Input from string buffer, using a ISstream
class IStringStream
:
    public StringStreamAllocator<std::istringstream>,
    public ISstream
{
    typedef StringStreamAllocator<std::istringstream> allocator_type;

public:

    // Constructors

        //- Construct from string
        IStringStream
        (
            const Foam::string& buffer,
            streamFormat format=ASCII,
            versionNumber version=currentVersion,
            const Foam::string& name="input"
        )
        :
            allocator_type(buffer),
            ISstream(stream_, name, format, version)
        {}


        //- Construct from char*
        IStringStream
        (
            const char* buffer,
            streamFormat format=ASCII,
            versionNumber version=currentVersion,
            const Foam::string& name="input"
        )
        :
            allocator_type(buffer),
            ISstream(stream_, name, format, version)
        {}


        //- Construct as copy of content
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

        //- Print description to Ostream
        virtual void print(Ostream& os) const;


    // Member operators

        //- Return a non-const reference to const Istream
        //  Needed for read-constructors where the stream argument is temporary.
        Istream& operator()() const
        {
            return const_cast<IStringStream&>(*this);
        }

};


/*---------------------------------------------------------------------------*\
                        Class OStringStream Declaration
\*---------------------------------------------------------------------------*/

//- Output to string buffer, using a OSstream
class OStringStream
:
    public StringStreamAllocator<std::ostringstream>,
    public OSstream
{
    typedef StringStreamAllocator<std::ostringstream> allocator_type;

public:

    // Constructors

        //- Construct and set stream status
        OStringStream
        (
            streamFormat format=ASCII,
            versionNumber version=currentVersion
        )
        :
            allocator_type(),
            OSstream(stream_, "output", format, version)
        {}


        //- Construct as copy of content
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

        //- Print description to Ostream
        virtual void print(Ostream& os) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
