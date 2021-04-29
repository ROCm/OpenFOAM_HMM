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

Class
    Foam::IListStream

Description
    An input stream that reads from a List and manages the List storage.
    Similar to IStringStream but with a List for its storage instead of
    as string to allow reuse of List contents without copying.

See Also
    Foam::OListStream
    Foam::UIListStream
    Foam::UOListStream

\*---------------------------------------------------------------------------*/

#ifndef IListStream_H
#define IListStream_H

#include "List.H"
#include "UIListStream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Detail
{

/*---------------------------------------------------------------------------*\
                Class Detail::IListStreamAllocator Declaration
\*---------------------------------------------------------------------------*/

//- An stream/stream-buffer input allocator with List storage
class IListStreamAllocator
:
    private List<char>,
    public UIListStreamAllocator
{
protected:

    // Constructors

        //- Default construct
        IListStreamAllocator()
        :
            List<char>(),
            UIListStreamAllocator(List<char>::data(), List<char>::size())
        {}

        //- Move construct from List
        IListStreamAllocator(List<char>&& buffer)
        :
            List<char>(std::move(buffer)),
            UIListStreamAllocator(List<char>::data(), List<char>::size())
        {}

        //- Move construct from DynamicList
        template<int SizeMin>
        IListStreamAllocator(DynamicList<char,SizeMin>&& buffer)
        :
            List<char>(std::move(buffer)),
            UIListStreamAllocator(List<char>::data(), List<char>::size())
        {}


    // Protected Member Functions

        //- Convenience method to address the underlying List storage
        inline void reset_gbuffer()
        {
            UIListStreamAllocator::reset
            (
                List<char>::data(),
                List<char>::size()
            );
        }

public:

    // Member Functions

        //- The current get position in the buffer
        using UIListStreamAllocator::size;

        //- Clear storage
        inline void clearStorage()
        {
            List<char>::clear();
            reset_gbuffer();
        }

        //- Transfer contents to other List
        inline void swap(List<char>& list)
        {
            List<char>::swap(list);
            reset_gbuffer();
        }
};

} // End namespace Detail


/*---------------------------------------------------------------------------*\
                         Class IListStream Declaration
\*---------------------------------------------------------------------------*/

//- An ISstream with internal List storage. Always UNCOMPRESSED.
class IListStream
:
    public Detail::IListStreamAllocator,
    public ISstream
{
    typedef Detail::IListStreamAllocator allocator_type;

public:

    // Constructors

        //- Default construct with an empty list
        explicit IListStream
        (
            IOstreamOption streamOpt = IOstreamOption()
        )
        :
            allocator_type(),
            ISstream(stream_, "input", streamOpt.format(), streamOpt.version())
        {}

        //- Move construct from List
        explicit IListStream
        (
            ::Foam::List<char>&& buffer,  // Fully qualify (issue #1521)
            IOstreamOption streamOpt = IOstreamOption()
        )
        :
            allocator_type(std::move(buffer)),
            ISstream(stream_, "input", streamOpt.format(), streamOpt.version())
        {}


        //- Move construct from DynamicList
        template<int SizeMin>
        explicit IListStream
        (
            DynamicList<char,SizeMin>&& buffer,
            IOstreamOption streamOpt = IOstreamOption()
        )
        :
            allocator_type(std::move(buffer)),
            ISstream(stream_, "input", streamOpt.format(), streamOpt.version())
        {}


    // Member Functions

        //- The current get position in the buffer
        using allocator_type::size;


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


        //- Print stream description to Ostream
        virtual void print(Ostream& os) const;


    // Member Operators

        //- A non-const reference to const Istream
        //  Needed for read-constructors where the stream argument is temporary
        Istream& operator()() const
        {
            // Could also rewind
            return const_cast<IListStream&>(*this);
        }


    // Additional constructors and methods (as per v2012 and earlier)
    #ifdef Foam_IOstream_extras

        //- Construct with an empty list
        explicit IListStream
        (
            IOstreamOption::streamFormat fmt,
            IOstreamOption::versionNumber ver = IOstreamOption::currentVersion
        )
        :
            IListStream(IOstreamOption(fmt, ver))
        {}


        //- Move construct from List
        IListStream
        (
            ::Foam::List<char>&& buffer,  // Fully qualify (issue #1521)
            IOstreamOption::streamFormat fmt,
            IOstreamOption::versionNumber ver = IOstreamOption::currentVersion
        )
        :
            IListStream(std::move(buffer), IOstreamOption(fmt, ver))
        {}


        //- Move construct from DynamicList
        template<int SizeMin>
        explicit IListStream
        (
            DynamicList<char,SizeMin>&& buffer,
            IOstreamOption::streamFormat fmt,
            IOstreamOption::versionNumber ver = IOstreamOption::currentVersion
        )
        :
            IListStream(std::move(buffer), IOstreamOption(fmt, ver))
        {}

    #endif /* Foam_IOstream_extras */
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
