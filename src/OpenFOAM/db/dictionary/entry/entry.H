/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

Class
    Foam::entry

Description
    A keyword and a list of tokens is an 'entry'.

    An entry can be read, written and printed, and the types and values of
    its tokens analysed.  An entry is a high-level building block for data
    description.  It is a front-end for the token parser. A list of entries
    can be used as a set of keyword syntax elements, for example.

SourceFiles
    entry.C
    entryIO.C

\*---------------------------------------------------------------------------*/

#ifndef entry_H
#define entry_H

#include "keyType.H"
#include "IDLList.H"
#include "fileName.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class ITstream;
class dictionary;

// Forward declaration of friend functions and operators

class entry;
Ostream& operator<<(Ostream& os, const entry& e);

/*---------------------------------------------------------------------------*\
                           Class entry Declaration
\*---------------------------------------------------------------------------*/

class entry
:
    public IDLList<entry>::link
{
public:

    // Enums/Typedefs

        //- The input mode options
        enum class inputMode
        {
            MERGE,      //!< Merge sub-dictionaries when possible
            OVERWRITE,  //!< Keep last entry. Silently remove previous ones
            PROTECT,    //!< Keep initial entry. Silently ignore subsequent ones
            WARN,       //!< Keep initial entry. Warn about subsequent ones.
            ERROR,      //!< FatalError for duplicate entries
            GLOBAL      //!< Use global value from #globalInputMode variable
        };


private:

    // Private data

        //- Keyword of entry
        keyType keyword_;


    // Private Member Functions

        //- Get the next valid keyword.
        //  \return True if it is a valid keyType.
        static bool getKeyword
        (
            keyType& keyword,
            token& keyToken,
            Istream& is
        );

        //- Get the next valid keyword.
        //  Warn when an invalid token is encountered that is not EOF or
        //  end-of-block
        //  \return True if it is a valid keyType.
        static bool getKeyword(keyType& keyword, Istream& is);


public:

    //- Enable or disable use of function entries and variable expansions.
    static int disableFunctionEntries;

    //- The current global input-mode.
    static inputMode globalInputMode;


    // Constructors

        //- Construct from keyword
        entry(const keyType& keyword);

        //- Construct as copy
        entry(const entry& e);

        //- Construct on freestore as copy with reference to the
        //  dictionary the copy belongs to
        virtual autoPtr<entry> clone
        (
            const dictionary& parentDict
        ) const = 0;

        //- Construct on freestore as copy
        //  Note: the parent directory is set to dictionary::null
        virtual autoPtr<entry> clone() const;

        //- Construct from Istream and insert into dictionary
        static bool New
        (
            dictionary& parentDict,
            Istream& is,
            const inputMode inMode = inputMode::GLOBAL
        );

        //- Construct on freestore from Istream and return
        static autoPtr<entry> New(Istream& is);

        //- Reset the #globalInputMode to %merge
        static void resetInputMode();


    //- Destructor
    virtual ~entry()
    {}


    // Member functions

        //- Return keyword
        const keyType& keyword() const
        {
            return keyword_;
        }

        //- Return non-const access to keyword
        keyType& keyword()
        {
            return keyword_;
        }

        //- Return the dictionary name
        virtual const fileName& name() const = 0;

        //- Return the dictionary name
        virtual fileName& name() = 0;

        //- Return line number of first token in dictionary
        virtual label startLineNumber() const = 0;

        //- Return line number of last token in dictionary
        virtual label endLineNumber() const = 0;


        //- Return true if this entry is a stream
        virtual bool isStream() const
        {
            return false;
        }

        //- Return token stream, if entry is a primitive entry
        virtual ITstream& stream() const = 0;


        //- Return true if this entry is a dictionary
        virtual bool isDict() const
        {
            return this->dictPtr();
        }

        //- Return pointer to dictionary, if entry is a dictionary.
        //  Return nullptr if the entry is not a dictionary.
        virtual const dictionary* dictPtr() const
        {
            return nullptr;
        }

        //- Return non-const pointer to dictionary, if entry is a dictionary
        //  Return nullptr if the entry is not a dictionary.
        virtual dictionary* dictPtr()
        {
            return nullptr;
        }

        //- Return dictionary, if entry is a dictionary
        virtual const dictionary& dict() const = 0;

        //- Return non-const access to dictionary, if entry is a dictionary
        virtual dictionary& dict() = 0;


        //- Write
        virtual void write(Ostream& os) const = 0;


    // Member operators

        void operator=(const entry& e);

        bool operator==(const entry& e) const;
        bool operator!=(const entry& e) const;


    // Ostream operator

        friend Ostream& operator<<(Ostream& os, const entry& e);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
