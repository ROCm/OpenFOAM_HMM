/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017-2018 OpenCFD Ltd.
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

        //- Construct empty
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

        //- Transfer contents to other list
        inline void swap(List<char>& list)
        {
            List<char>::swap(list);
            reset_gbuffer();
        }
};

} // End namespace Detail


/*---------------------------------------------------------------------------*\
                         Class IListStream Declaration
\*----------------------------------------------d-----------------------------*/

//- An ISstream with internal List storage
class IListStream
:
    public Detail::IListStreamAllocator,
    public ISstream
{
    typedef Detail::IListStreamAllocator allocator_type;

public:

    // Constructors

        //- Construct with an empty list
        IListStream
        (
            streamFormat format=ASCII,
            versionNumber version=currentVersion,
            const Foam::string& name="input"
        )
        :
            allocator_type(),
            ISstream(stream_, name, format, version)
        {}


        //- Move construct from List
        IListStream
        (
            List<char>&& buffer,
            streamFormat format=ASCII,
            versionNumber version=currentVersion,
            const Foam::string& name="input"
        )
        :
            allocator_type(std::move(buffer)),
            ISstream(stream_, name, format, version)
        {}


        //- Move construct from DynamicList
        template<int SizeMin>
        IListStream
        (
            DynamicList<char,SizeMin>&& buffer,
            streamFormat format=ASCII,
            versionNumber version=currentVersion,
            const Foam::string& name="input"
        )
        :
            allocator_type(std::move(buffer)),
            ISstream(stream_, name, format, version)
        {}


    // Member functions

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


        //- Print description to Ostream
        virtual void print(Ostream& os) const;


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
