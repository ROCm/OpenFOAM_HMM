/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2016 OpenCFD Ltd.
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
    Foam::functionObjects::writeFile

Description
    functionObject base class for writing single files

See also
    Foam::functionObject
    Foam::functionObjects::logFiles

SourceFiles
    writeFile.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_writeFile_H
#define functionObjects_writeFile_H

#include "objectRegistry.H"
#include "OFstream.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class writeFile Declaration
\*---------------------------------------------------------------------------*/

class writeFile
{

protected:

    // Protected data

        //- Reference to the region objectRegistry
        const objectRegistry& fileObr_;

        //- Prefix
        const word prefix_;

        //- Name of file
        word fileName_;

        //- File pointer
        autoPtr<OFstream> filePtr_;

        //- Write precision
        label writePrecision_;

        //- Flag to enable/disable writing to file
        bool writeToFile_;

        //- Flag to identify whether the header has been written
        bool writtenHeader_;


    // Protected Member Functions

        //- Initialise the output stream for writing
        void initStream(Ostream& os) const;

        //- Return the base directory for output
        fileName baseFileDir() const;

        //- Return the base directory for the current time value
        fileName baseTimeDir() const;

        //- Return an autoPtr to a new file
        virtual autoPtr<OFstream> createFile(const word& name) const;

        //- Reset internal file pointer to new file with new name
        virtual void resetFile(const word& name);

        //- Return the value width when writing to stream with optional offset
        Omanip<int> valueWidth(const label offset = 0) const;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        writeFile(const writeFile&);

        //- Disallow default bitwise assignment
        void operator=(const writeFile&);


public:

    //- Directory prefix
    static const word outputPrefix;

    //- Additional characters for writing
    static label addChars;


    // Constructors

        //- Construct from objectRegistry and prefix
        writeFile
        (
            const objectRegistry& obr,
            const word& prefix
        );

        //- Construct from objectRegistry and prefix, and read options
        //  from dictionary
        writeFile
        (
            const objectRegistry& obr,
            const word& prefix,
            const word& fileName,
            const dictionary& dict
        );


    //- Destructor
    virtual ~writeFile();


    // Member Functions

        //- Read
        virtual bool read(const dictionary& dict);

        //- Return access to the file (if only 1)
        virtual OFstream& file();

        //- Flag to allow writing to file
        virtual bool writeToFile() const;

        //- Return width of character stream output
        virtual label charWidth() const;

        //- Write a commented string to stream
        virtual void writeCommented(Ostream& os, const string& str) const;

        //- Write a tabbed string to stream
        virtual void writeTabbed(Ostream& os, const string& str) const;

        //- Write a commented header to stream
        virtual void writeHeader(Ostream& os, const string& str) const;

        //- Write the current time to stream
        virtual void writeTime(Ostream& os) const;

        //- Write a break marker to the stream
        virtual void writeBreak(Ostream& os) const;

        //- Write a (commented) header property and value pair
        template<class Type>
        void writeHeaderValue
        (
            Ostream& os,
            const string& property,
            const Type& value
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "writeFileTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
