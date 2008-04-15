/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::DictionaryWriter.H

Description

SourceFiles
    DictionaryWriter.C

\*---------------------------------------------------------------------------*/

#ifndef DictionaryWriter_H
#define DictionaryWriter_H

// Foam header files.
#include "word.H"
#include "wordList.H"
#include "string.H"
#include "stringList.H"
#include "fileName.H"
#include "OFstream.H"

// FoamX header files.
#include "FoamX.H"
#include "FoamXErrors.H"
#include "FoamXString.H"
#include "FoamXStringList.H"
#include "FoamXWordList.H"
#include "FoamXFileNameList.H"
#include "FoamXAnyList.H"

// Namespaces
#include "FoamXNameSpaces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace FoamX
{

/*---------------------------------------------------------------------------*\
                           Class DictionaryWriter Declaration
\*---------------------------------------------------------------------------*/

class DictionaryWriter
{
    // Private data

        OFstream* pFileStream_;
        fileName  fileName_;
        fileName  caseRoot_;
        fileName  caseName_;
        fileName  instance_;
        fileName  object_;
        int       entryNameSize_;

 
    // Private member functions

        OFstream& file()
        {
            return (*pFileStream_);
        }

        void writeChars(char ch, int count);


public:

    // Constructors

        DictionaryWriter(const fileName& fileName);

        DictionaryWriter
        (
            const fileName& caseRoot,
            const fileName& caseName,
            const fileName& instance,
            const fileName& object
        );


    // Destructor

        ~DictionaryWriter();


    // Member Functions

        // Access

            const fileName& pathName() const
            {
                return fileName_;
            }


        // Header methods.

            void writeHeader
            (
                const Foam::string& title,
                const word& className
            );

            void writeSectionHeader(const Foam::string& name);


        // Entry name formatting methods.

            void writeKeyword(const word& name);
            void writeKeywordOnly(const word& name);


        // Entry value formatting methods.

            void writeValue(const FoamXString& value)
            {
                value.write(file());
            }

            void writeValue(const Foam::token& value)
            {
                file() << value;
            }

            void writeValue(const Foam::string& value)
            {
                file() << value;
            }

            void writeValue(const word& value)
            {
                file() << value;
            }

            void writeValue(const scalar& value)
            {
                file() << value;
            }

            void writeValue(const label& value)
            {
                file() << value;
            }

            void writeValue(const bool& value)
            {
                file() << value;
            }

            void writeValue(const FoamXAny& value)
            {
                value.write(file());
            }

            void writeValue
            (
                const FoamXServer::DimensionSet& value
            );


        // Entry methods.

            void writeEntry(const word& name, const char* value);
            void writeEntry(const word& name, const FoamXString& value);
            void writeEntry(const word& name, const Foam::string& value);
            void writeEntry(const word& name, const word& value);
            void writeEntry(const word& name, const scalar& value);
            void writeEntry(const word& name, const label& value);
            void writeEntry(const word& name, const bool& value);
            void writeEntry(const word& name, const FoamXAny& value);
            void writeEntry
            (
                const word& name,
                const FoamXServer::DimensionSet& value
            );


        // FixedList methods

            void startFixedList();
            void endFixedList();


        // List methods.

            void startList(const label& numElements);
            void endList();

            void writeEntry
            (
                const word& name,
                const FoamXServer::StringList& list
            );

            void writeEntry
            (
                const word& name,
                const FoamXStringList& list
            );

            void writeEntry
            (
                const word& name,
                const FoamXWordList& list
            );

            void writeEntry
            (
                const word& name,
                const FoamXAnyList& list
            );


        // Sub dictionary methods.

            void startSubDict();
            void startSubDict(const word& name);
            void endDict();
            void endSubDict();

            // Miscellaneous methods.
            void writeString(const Foam::string& text);
            void writeLine(const Foam::string& text);
            void writeComment(const Foam::string& text);
            void writeBar();
            void writeEndBar();
            void writeEndl();
            void endEntry();


        // General formatting

            void indent();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FoamX

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
