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

/*---------------------------------------------------------------------------*\
                    Class IListStreamAllocator Declaration
\*---------------------------------------------------------------------------*/

//- An stream/stream-buffer input allocator for with List storage
class IListStreamAllocator
:
    private List<char>,
    public UIListStreamAllocator
{
protected:

    // Constructors

        //- Construct with an empty list
        IListStreamAllocator()
        :
            List<char>(),
            UIListStreamAllocator(list_storage())
        {}

        //- Move construct from an existing List
        IListStreamAllocator(List<char>&& buffer)
        :
            List<char>(std::move(buffer)),
            UIListStreamAllocator(list_storage())
        {}


        //- Move construct from an existing DynamicList
        template<int SizeMin>
        IListStreamAllocator(DynamicList<char,SizeMin>&& buffer)
        :
            List<char>(std::move(buffer)),
            UIListStreamAllocator(list_storage())
        {}

        //- Transfer (move) construct
        IListStreamAllocator(const Xfer<List<char>>& buffer)
        :
            List<char>(buffer),
            UIListStreamAllocator(list_storage())
        {}


    // Protected Member Functions

        //- Convenience method to address the underlying List storage
        inline List<char>& list_storage()
        {
            return static_cast<List<char>&>(*this);
        }

public:

    // Member Functions

        //- The current get position in the buffer
        using UIListStreamAllocator::size;

        //- Clear storage
        inline void clearStorage()
        {
            list_storage().clear();
            sync_gbuffer_to_list();
        }

        //- Transfer contents to the Xfer container as a plain List
        inline Xfer<List<char>> xfer()
        {
            Xfer<List<char>> moved = list_storage().xfer();

            // Update buffer pointers for remaining (zero) size list
            sync_gbuffer_to_list(); // or rewind

            return moved;
        }
};


/*---------------------------------------------------------------------------*\
                         Class IListStream Declaration
\*----------------------------------------------d-----------------------------*/

//- An ISstream with internal List storage
class IListStream
:
    public IListStreamAllocator,
    public ISstream
{
    typedef IListStreamAllocator allocator_type;

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


        //- Move construct from an existing List
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


        //- Move construct from an existing DynamicList
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


        //- Transfer (move) construct
        IListStream
        (
            const Xfer<List<char>>& buffer,
            streamFormat format=ASCII,
            versionNumber version=currentVersion,
            const Foam::string& name="input"
        )
        :
            allocator_type(buffer),
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
