/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2011 OpenCFD Ltd.
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
    Foam::codeStreamTools

Description
    Base for all things on-the-fly from dictionary

SourceFiles
    codeStreamTools.C

\*---------------------------------------------------------------------------*/

#ifndef codeStreamTools_H
#define codeStreamTools_H

#include "Tuple2.H"
#include "Pair.H"
#include "SHA1Digest.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class ISstream;
class OSstream;

/*---------------------------------------------------------------------------*\
                        Class codeStreamTools Declaration
\*---------------------------------------------------------------------------*/

class codeStreamTools
{
public:
    typedef Tuple2<fileName, string> fileAndContent;

    //- Helper class for managing file and variables
    class fileAndVars
    :
        public HashTable<string>
    {
        // Private data
        fileName file_;

    public:
        //- Construct null
        fileAndVars()
        {}

        //- Return the file name
        const fileName& file() const
        {
            return file_;
        }

        //- Return the file name
        fileName& file()
        {
            return file_;
        }
    };


private:
    // Private data

        //- Name for underlying set
        word name_;

        //- Files to copy
        List<codeStreamTools::fileAndVars> copyFiles_;

        //- Direct contents for files
        List<fileAndContent> filesContents_;

protected:

        void copyAndExpand
        (
            ISstream&,
            OSstream&,
            const HashTable<string>& mapping
        ) const;

public:

    // Static data members

        //- Name of the code template environment variable
        //  Used to located the codeTemplateName
        const static word codeTemplateEnvName;

        //- Name of the code template sub-directory
        //  Used when locating the codeTemplateName via Foam::findEtcFile
        const static fileName codeTemplateDirName;


        static int allowSystemOperations;

    // Constructors

        //- Construct null
        codeStreamTools();

        //- Construct from dictionary
        codeStreamTools(const word& name, const dictionary& dict);

        //- Copy from components
        codeStreamTools
        (
            const word& name,
            const List<fileAndVars>&,
            const List<fileAndContent>&
        );

        //- Construct copy
        codeStreamTools(const codeStreamTools&);


    // Member functions

        const word& name() const
        {
            return name_;
        }

        const List<fileAndVars>& copyFiles() const
        {
            return copyFiles_;
        }

        const List<Tuple2<fileName, string> >& filesContents() const
        {
            return filesContents_;
        }

        bool copyFilesContents(const fileName& dir) const;


        //- Find a code-template via the codeTemplateEnvName
        //  alternatively in the codeTemplateDirName via Foam::findEtcFile
        static fileName findTemplate(const word& templateName);

        //- List searched locations in a format suitable for display an error
        static string searchedLocations();


        static void* findLibrary(const fileName& libPath);

        static bool writeDigest(const fileName& dir, const SHA1Digest& sha1);
        static SHA1Digest readDigest(const fileName& dir);
        static bool upToDate(const fileName& dir, const SHA1Digest& sha1);

        bool read(const dictionary&);

        void writeDict(Ostream&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
