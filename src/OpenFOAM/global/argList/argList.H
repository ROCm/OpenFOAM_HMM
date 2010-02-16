/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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
    Foam::argList

Description
    Extract command arguments and options from the supplied
    @a argc and @a argv parameters.

    Sequences with "(" ... ")" are transformed into a stringList.
    For example,
    @verbatim
        program -listFiles \( *.txt \)
    @endverbatim
    would create a stringList:
    @verbatim
        ( "file1.txt" "file2.txt" ... "fileN.txt" )
    @endverbatim
    The backslash-escaping has been used to avoid shell expansions.

    @par Default command-line options
    @param -case \<dir\> \n
        select an case directory instead of the current working directory
    @param -parallel \n
        specify case as a parallel job
    @param -doc \n
        display the documentation in browser
    @param -srcDoc \n
        display the source documentation in browser
    @param -help \n
       print the usage

    The environment variable @b FOAM_CASE is set to the path of the
    global case (same for serial and parallel jobs).
    The environment variable @b FOAM_CASENAME is set to the name of the
    global case.

Note
    - Adjustment of the valid (mandatory) arguments
      by directly manipulating the argList::validArgs static member.
    - Adjustment of the valid options
      via the addOption/removeOption static methods or by directly
      manipulating the argList::validOptions static member.

SourceFiles
    argList.C

\*---------------------------------------------------------------------------*/

#ifndef argList_H
#define argList_H

#include "stringList.H"
#include "SubList.H"
#include "SLList.H"
#include "HashTable.H"
#include "word.H"
#include "fileName.H"
#include "parRun.H"
#include "IStringStream.H"

#include "sigFpe.H"
#include "sigInt.H"
#include "sigQuit.H"
#include "sigSegv.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class argList
{
    // Private data
        static bool bannerEnabled;

        stringList args_;
        HashTable<string> options_;

        word executable_;
        fileName rootPath_;
        fileName globalCase_;
        fileName case_;

        ParRunControl parRunControl_;

        // Signal handlers
        sigFpe sigFpe_;
        sigInt sigInt_;
        sigQuit sigQuit_;
        sigSegv sigSegv_;


    // Private Member Functions

        //- Helper function for printUsage
        static void printOptionUsage
        (
            const label location,
            const string& str
        );

        //- get rootPath_ / globalCase_ from one of the following forms
        //   * [-case dir]
        //   * cwd
        //
        // Also export FOAM_CASE and FOAM_CASENAME environment variables
        // so they can be used immediately (eg, in decomposeParDict)
        void getRootCase();

        //- Transcribe argv into internal args_
        //  return true if any "(" ... ")" sequences were captured
        bool regroupArgv(int& argc, char**& argv);


public:

    // Static data members

        //- A list of valid (mandatory) arguments
        static SLList<string> validArgs;

        //- A list of valid options
        static HashTable<string> validOptions;

        //- A list of valid parallel options
        static HashTable<string> validParOptions;

        //- Short usage information for validOptions
        static HashTable<string> optionUsage;

        //- Additional notes for usage
        static SLList<string> notes;

        //- Min offset for displaying usage (default: 20)
        static string::size_type usageMin;

        //- Max screen width for displaying usage (default: 80)
        static string::size_type usageMax;

        //! @cond ignoreDocumentation
        class initValidTables
        {
        public:

            initValidTables();
        };
        //! @endcond ignoreDocumentation


    // Constructors

        //- Construct from argc and argv
        //  checking the arguments and options as requested
        argList
        (
            int& argc,
            char**& argv,
            bool checkArgs=true,
            bool checkOpts=true
        );


        //- Destructor
        virtual ~argList();


    // Member functions

        // Access

            //- Name of executable
            inline const word& executable() const;

            //- Return root path
            inline const fileName& rootPath() const;

            //- Return case name (parallel run) or global case (serial run)
            inline const fileName& caseName() const;

            //- Return case name
            inline const fileName& globalCaseName() const;

            //- Return the path to the caseName
            inline fileName path() const;

            //- Return arguments
            inline const stringList& args() const;

            //- Return the number of arguments
            inline label size() const;

            //- Return the argument corresponding to index.
            //  Index 0 corresponds to the name of the executable.
            //  Index 1 corresponds to the first argument.
            inline const string& operator[](const label index) const;

            //- Read a value from the argument at index.
            //  Index 0 corresponds to the name of the executable.
            //  Index 1 corresponds to the first argument.
            template<class T>
            inline T argRead(const label index) const;

            //- Return options
            inline const Foam::HashTable<string>& options() const;

            //- Return the argument string associated with the named option
            inline const string& option(const word& opt) const;

            //- Return true if the named option is found
            inline bool optionFound(const word& opt) const;

            //- Return an IStringStream from the named option
            inline IStringStream optionLookup(const word& opt) const;

            //- Read a value from the named option
            template<class T>
            inline T optionRead(const word& opt) const;

            //- Read a value from the named option if present.
            //  Return true if the named option was found.
            template<class T>
            inline bool optionReadIfPresent(const word& opt, T&) const;

            //- Read a value from the named option if present.
            //  Return true if the named option was found, otherwise
            //  use the supplied default and return false.
            template<class T>
            inline bool optionReadIfPresent
            (
                const word& opt,
                T&,
                const T& deflt
            ) const;

            //- Read a value from the named option if present.
            //  Return true if the named option was found.
            template<class T>
            inline T optionLookupOrDefault
            (
                const word& opt,
                const T& deflt
            ) const;

            //- Read a List of values from the named option
            template<class T>
            List<T> optionReadList(const word& opt) const
            {
                return readList<T>(optionLookup(opt)());
            }


        // Edit

            //- Add to a bool option to validOptions with usage information
            static void addBoolOption
            (
                const word& opt,
                const string& usage = ""
            );

            //- Add to an option to validOptions with usage information
            //  An option with an empty param is a bool option
            static void addOption
            (
                const word& opt,
                const string& param = "",
                const string& usage = ""
            );

            //- Add option usage information to optionUsage
            static void addUsage
            (
                const word& opt,
                const string& usage
            );

            //- Add extra notes for the usage information
            //  This string is used "as-is" without additional formatting
            static void addNote(const string&);

            //- Remove option from validOptions and from optionUsage
            static void removeOption(const word& opt);

            //- Disable emitting the banner information
            static void noBanner();

            //- Remove the parallel options
            static void noParallel();


        // Print

            //- Print usage
            void printUsage() const;

            //- Display documentation in browser
            //  Optionally display the application source code
            void displayDoc(bool source=false) const;


        // Check

            //- Check argument list
            bool check(bool checkArgs=true, bool checkOpts=true) const;

            //- Check root path and case path
            bool checkRootCase() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "argListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
