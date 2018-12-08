/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2018 OpenCFD Ltd.
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
    Foam::topoSetSource

Description
    Base class of a source for a topoSet.

    Implementer must modify the given set (see applyToSet) according to
    its function and the setAction (one of add/delete/new)

SourceFiles
    topoSetSource.C

\*---------------------------------------------------------------------------*/

#ifndef topoSetSource_H
#define topoSetSource_H

#include "pointField.H"
#include "labelList.H"
#include "faceList.H"
#include "typeInfo.H"
#include "autoPtr.H"
#include "Enum.H"
#include "HashTable.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class polyMesh;
class topoSet;

/*---------------------------------------------------------------------------*\
                           Class topoSetSource Declaration
\*---------------------------------------------------------------------------*/

class topoSetSource
{
public:

    // Public data types

        //- Enumeration defining the types of sources
        enum sourceType
        {
            CELLSETSOURCE = 0x1,        //!< Cells
            FACESETSOURCE = 0x2,        //!< Faces
            POINTSETSOURCE = 0x4,       //!< Points

            CELLZONESOURCE = 0x11,      //!< Cells as zone
            FACEZONESOURCE = 0x12,      //!< Faces as zone
            POINTZONESOURCE = 0x14,     //!< Points as zone
        };

        //- Enumeration defining the valid actions
        enum setAction
        {
            ADD,                //!< Add elements to the set
            SUBTRACT,           //!< Subtract elements from the set
            SUBSET,             //!< Subset with elements in the set
            INVERT,             //!< Invert the elements in the set
            CLEAR,              //!< Clear the set, possibly creating it
            NEW,                //!< Create a new set and ADD elements to it
            REMOVE,             //!< Remove the set (from the file system)
            LIST,               //!< Print contents of the set
            DELETE = SUBTRACT,  //!< \deprecated(2018-10) Alias for SUBTRACT
        };

        //- The setActions text representations
        static const Enum<setAction> actionNames;


protected:

        //- A table of usage strings
        static HashTable<string>* usageTablePtr_;

        //- Class with constructor to add usage string to table
        class addToUsageTable
        {
        public:

            addToUsageTable(const word& name, const string& msg)
            {
                if (!usageTablePtr_)
                {
                    usageTablePtr_ = new HashTable<string>();
                }
                usageTablePtr_->insert(name, msg);
            }

            ~addToUsageTable()
            {
                if (usageTablePtr_)
                {
                    delete usageTablePtr_;
                    usageTablePtr_ = nullptr;
                }
            }
        };


    // Protected Data

        //- Reference to the mesh
        const polyMesh& mesh_;

        //- Verbosity (default: true)
        bool verbose_;


    // Protected Member Functions

        //- Detect and remove any values less than 0 or ge maxLabel.
        //  \return false if invalid elements were detected (and removed)
        static bool check(labelList& list, const label maxLabel);

        //- Add or delete id from set. Add when 'add' is true
        void addOrDelete(topoSet& set, const label id, const bool add) const;

        //- Add or delete labels from set. Add when 'add' is true
        void addOrDelete
        (
            topoSet& set,
            const labelUList& labels,
            const bool add
        ) const;


private:

        static const string illegalSource_;


    // Private Member Functions

        //- No copy construct
        topoSetSource(const topoSetSource&) = delete;

        //- No copy assignment
        void operator=(const topoSetSource&) = delete;


public:

    //- Runtime type information
    TypeName("topoSetSource");


    // Static Functions

        //- Deprecated(2018-07) convert string to action
        //  \deprecated(2018-07) - use actionNames[] directly
        static setAction toAction(const word& actionName)
        {
            return actionNames[actionName];
        }

        //- Check state of stream.
        static Istream& checkIs(Istream& is);


    // Declare run-time constructor selection table

        // For the dictionary constructor
        declareRunTimeSelectionTable
        (
            autoPtr,
            topoSetSource,
            word,
            (
                const polyMesh& mesh,
                const dictionary& dict
            ),
            (mesh, dict)
        );

        // For the Istream constructor
        declareRunTimeSelectionTable
        (
            autoPtr,
            topoSetSource,
            istream,
            (
                const polyMesh& mesh,
                Istream& is
            ),
            (mesh, is)
        );


        //- Class used for the read-construction of
        //  PtrLists of topoSetSource
        class iNew
        {
            const polyMesh& mesh_;

        public:

            iNew(const polyMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<topoSetSource> operator()(Istream& is) const
            {
                const word sourceTypeName(is);
                dictionary dict(is);
                return topoSetSource::New(sourceTypeName, mesh_, dict);
            }
        };


        static const string& usage(const word& name)
        {
            if (!usageTablePtr_)
            {
                usageTablePtr_ = new HashTable<string>();
            }

            return usageTablePtr_->lookup(name, illegalSource_);
        }


    // Constructors

        //- Construct from components
        explicit topoSetSource(const polyMesh& mesh);

        //- Clone (disallowed)
        autoPtr<topoSetSource> clone() const
        {
            NotImplemented;
            return nullptr;
        }


    // Selectors

        //- Return a reference to the selected topoSetSource
        static autoPtr<topoSetSource> New
        (
            const word& topoSetSourceType,
            const polyMesh& mesh,
            const dictionary& dict
        );

        //- Return a reference to the selected topoSetSource
        static autoPtr<topoSetSource> New
        (
            const word& topoSetSourceType,
            const polyMesh& mesh,
            Istream& is
        );


    //- Destructor
    virtual ~topoSetSource() = default;


    // Member Functions

        //- Reference to the mesh
        const polyMesh& mesh() const
        {
            return mesh_;
        }

        //- Return the current verbosity
        bool verbose() const
        {
            return verbose_;
        }

        //- Set the current verbosity
        void verbose(bool on)
        {
            verbose_ = on;
        }


    // Member Functions

        //- The type of set being used by the source
        virtual sourceType setType() const = 0;

        //- Apply specified action to the topoSet
        virtual void applyToSet
        (
            const setAction action,
            topoSet& set
        ) const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
