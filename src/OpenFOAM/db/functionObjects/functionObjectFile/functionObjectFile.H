/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015 OpenCFD Ltd
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
    Foam::functionObjectFile

Description
    Base class for output file data handling

See Also
    Foam::functionObject
    Foam::OutputFilterFunctionObject

SourceFiles
    functionObjectFile.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjectFile_H
#define functionObjectFile_H

#include "objectRegistry.H"
#include "OFstream.H"
#include "PtrList.H"
#include "HashSet.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                     Class functionObjectFile Declaration
\*---------------------------------------------------------------------------*/

class functionObjectFile
{
private:

    // Private data

        //- Reference to the database
        const objectRegistry& obr_;

        //- Prefix
        const word prefix_;

        //- Name of file
        word fileName_;

        //- File pointer
        autoPtr<OFstream> filePtr_;

        //- Write precision
        label writePrecision_;


protected:

    // Protected Data

        //- Flag to enable/disable writing to file
        bool writeToFile_;


    // Protected Member Functions

        //- Initialise the output stream for writing
        virtual void initStream(Ostream& os) const;

        //- Return the base directory for output
        virtual fileName baseFileDir() const;

        //- Return the base directory for the current time value
        virtual fileName baseTimeDir() const;

        //- Return an autoPtr to a new file
        virtual autoPtr<OFstream> createFile(const word& name) const;

        //- Reset internal file pointer to new file with new name
        virtual void resetFile(const word& name);

        //- Return the value width when writing to stream with optional offset
        virtual Omanip<int> valueWidth(const label offset = 0) const;

        //- Disallow default bitwise copy construct
        functionObjectFile(const functionObjectFile&);

        //- Disallow default bitwise assignment
        void operator=(const functionObjectFile&);


public:

    //- Directory prefix
    static const word outputPrefix;

    //- Additional characters for writing
    static label addChars;

    // Constructors

        //- Construct null
        functionObjectFile(const objectRegistry& obr, const word& prefix);

        //- Construct from components and read options from dictionary
        functionObjectFile
        (
            const objectRegistry& obr,
            const word& prefix,
            const word& fileName,
            const dictionary& dict
        );


    //- Destructor
    virtual ~functionObjectFile();


    // Member Functions

        //- Read
        void read(const dictionary& dict);

        //- Return access to the file (if only 1)
        OFstream& file();

        //- Return true if can write to file
        bool writeToFile() const;

        //- Write a commented string to stream
        void writeCommented
        (
            Ostream& os,
            const string& str
        ) const;

        //- Write a tabbed string to stream
        void writeTabbed
        (
            Ostream& os,
            const string& str
        ) const;

        //- Write a commented header to stream
        void writeHeader
        (
            Ostream& os,
            const string& str
        ) const;

        //- Write the current time to stream
        void writeTime(Ostream& os) const;

        //- Write a (commented) header property and value pair
        template<class Type>
        void writeHeaderValue
        (
            Ostream& os,
            const string& property,
            const Type& value
        ) const;

        //- Return width of character stream output
        label charWidth() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "functionObjectFileTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
