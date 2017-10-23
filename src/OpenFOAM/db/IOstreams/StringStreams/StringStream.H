/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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

        //- The allocated stream pointer.
        StreamType* allocatedPtr_;


    // Constructors

        //- Construct null
        StringStreamAllocator()
        :
            allocatedPtr_(new StreamType())
        {}

        //- Construct from pointer, taking ownership
        StringStreamAllocator(StreamType* ptr)
        :
            allocatedPtr_(ptr)
        {}

        //- Construct from string
        StringStreamAllocator(const std::string& buffer)
        :
            allocatedPtr_(new StreamType(buffer))
        {}


    //- Destructor
    ~StringStreamAllocator()
    {
        deallocate();
    }


    // Protected Member Functions

        //- Delete the stream pointer
        void deallocate()
        {
            if (allocatedPtr_)
            {
                delete allocatedPtr_;
                allocatedPtr_ = nullptr;
            }
        }


public:

    // Public Member Functions

        //- Get the string
        Foam::string str() const
        {
            return allocatedPtr_->str();
        }

        //- Set the string
        void str(const std::string& s)
        {
            allocatedPtr_->str(s);
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
            StringStreamAllocator<std::istringstream>(buffer),
            ISstream(*allocatedPtr_, name, format, version)
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
            StringStreamAllocator<std::istringstream>(buffer),
            ISstream(*allocatedPtr_, name, format, version)
        {}


    //- Destructor
    ~IStringStream()
    {}


    // Member Functions

        //- Print description to Ostream
        void print(Ostream& os) const;

        //- Reset the input buffer and rewind the stream
        void reset(const std::string& s)
        {
            this->str(s);
            this->rewind();
        }


    // Member operators

        //- Return a non-const reference to const Istream
        //  Needed for read-constructors where the stream argument is temporary:
        //  e.g. thing thisThing(IFstream("thingFileName")());
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
public:

    // Constructors

        //- Construct and set stream status
        OStringStream
        (
            streamFormat format=ASCII,
            versionNumber version=currentVersion
        )
        :
            StringStreamAllocator<std::ostringstream>(),
            OSstream(*allocatedPtr_, "output", format, version)
        {}

        //- Construct as copy
        OStringStream(const OStringStream& oss)
        :
            StringStreamAllocator<std::ostringstream>(oss.str()),
            OSstream(*allocatedPtr_, oss.name(), oss.format(), oss.version())
        {}


    //- Destructor
    ~OStringStream()
    {}


    // Member Functions

        //- Reset the output buffer and rewind the stream
        void reset()
        {
            this->str("");  // No other way to reset the end
            this->rewind();
        }

        //- Rewind the output stream
        void rewind()
        {
            // pubseekpos() instead of seekp() for symmetry with other classes
            allocatedPtr_->rdbuf()->pubseekpos(0, std::ios_base::out);
        }

        //- Print description to Ostream
        void print(Ostream& os) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
