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
    Foam::dynamicCode

Description
    Tools for handling dynamic code compilation

SourceFiles
    dynamicCode.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicCode_H
#define dynamicCode_H

#include "Tuple2.H"
#include "SHA1Digest.H"
#include "HashTable.H"
#include "DynamicList.H"
#include "dlLibraryTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class dynamicCodeContext;
class ISstream;
class OSstream;

/*---------------------------------------------------------------------------*\
                        Class dynamicCode Declaration
\*---------------------------------------------------------------------------*/

class dynamicCode
{
public:
    typedef Tuple2<fileName, string> fileAndContent;

private:
    // Private data

        //- Name for code
        word codeName_;

        //- Name for code subdirectory
        mutable word codeDirName_;

        //- Variables to use during filtering
        HashTable<string> filterVars_;

        //- Direct contents for files
        DynamicList<fileAndContent> createFiles_;

        //- Files to copy and filter
        DynamicList<fileName> filterFiles_;


protected:

        void copyAndFilter
        (
            ISstream&,
            OSstream&,
            const HashTable<string>& mapping
        ) const;

        //- Resolve code-templates via the codeTemplateEnvName
        //  alternatively in the codeTemplateDirName via Foam::findEtcFile
        static List<fileName> resolveTemplates
        (
            const UList<fileName>& names
        );


        bool writeDigest(const fileName& dir, const SHA1Digest& sha1) const;
        SHA1Digest readDigest(const fileName& dir) const;

public:

    // Static data members

        //- Name of the code template environment variable
        //  Used to located the codeTemplateName
        static const word codeTemplateEnvName;

        //- Name of the code template sub-directory
        //  Used when locating the codeTemplateName via Foam::findEtcFile
        static const fileName codeTemplateDirName;

        static int allowSystemOperations;


    // Static Member functions

        //- Check security for creating dynamic code
        static void checkSecurity
        (
            const char* title,
            const dictionary& dict
        );


    // Constructors

        //- Construct for a specified code name
        dynamicCode(const word& codeName);

        //- Construct for a specified code name and directory name
        dynamicCode(const word& codeName, const word& dirName);


    // Member functions

        //- Return the code-name
        const word& codeName() const
        {
            return codeName_;
        }

        //- Return the code-dirname
        const word& codeDirName() const
        {
            return codeDirName_;
        }

        //- Clear variables and files
        void clear();


        //- Add a file to create with its contents. Will not be filtered
        void addCreateFile(const fileName& name, const string& contents);

        //- Add a file template name, which will be found and filtered
        void addFilterFile(const fileName& name);

        //- Define filter variables for code, codeInclude, SHA1sum
        void setFilterContext(const dynamicCodeContext&);

        //- Define a filter variable
        void setFilterVariable(const word& key, const string& value);


        //- Local path for specified code name
        //  Expanded from \$FOAM_CASE/dynamicCode/codeDirName
        fileName codePath() const;

        //- Local library path for specified code name
        //  Expanded from \$FOAM_CASE/dynamicCode/platforms/\$WM_OPTIONS/lib
        fileName libPath() const;

        //- The library target path for Make/files
        string libTarget() const;


        //- Verify if the copied code is up-to-date
        bool upToDate(const dynamicCodeContext& context) const;

        //- Verify if the copied code is up-to-date
        bool upToDate(const SHA1Digest& sha1) const;

        //- Copy/create files prior to compilation
        bool copyFilesContents() const;

        //- Compile a libso
        bool wmakeLibso() const;

//        //- Open the libPath() library
//        bool openLibrary() const;
//
//        //- Close the libPath() library
//        bool closeLibrary() const;
//
//        //- Find the handle of the libPath() library
//        void* findLibrary() const;

//        bool read(const dictionary&);
//        void writeDict(Ostream&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
