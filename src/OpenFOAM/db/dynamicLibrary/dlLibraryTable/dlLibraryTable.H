/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::dlLibraryTable

Description
    A table of dynamically loaded libraries

SourceFiles
    dlLibraryTable.C

\*---------------------------------------------------------------------------*/

#ifndef dlLibraryTable_H
#define dlLibraryTable_H

#include "label.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class dlLibraryTable Declaration
\*---------------------------------------------------------------------------*/

class dlLibraryTable
{
    // Private data

        DynamicList<void*> libPtrs_;

        DynamicList<fileName> libNames_;


    // Private Member Functions

        //- No copy construct
        dlLibraryTable(const dlLibraryTable&) = delete;

        //- No copy assignment
        void operator=(const dlLibraryTable&) = delete;


public:

    // Declare name of the class and its debug switch
    ClassName("dlLibraryTable");

    // Constructors

        //- Construct null
        dlLibraryTable();

        //- Open all libraries listed in the 'libsEntry' entry in the
        //- given dictionary.
        dlLibraryTable(const dictionary& dict, const word& libsEntry);


    //- Destructor
    ~dlLibraryTable();


    // Member Functions

        //- Open the named library, optionally with warnings if problems occur
        bool open(const fileName& libName, const bool verbose = true);

        //- Close the named library, optionally with warnings if problems occur
        bool close(const fileName& libName, const bool verbose = true);

        //- Find the handle of the named library
        void* findLibrary(const fileName& libName);

        //- Open all libraries listed in the 'libsEntry' entry in the
        //- given dictionary.
        bool open(const dictionary& dict, const word& libsEntry);

        //- Open all libraries listed in the 'libsEntry' entry in the
        //- given dictionary and check the additions
        //- to the given constructor table
        template<class TablePtr>
        bool open
        (
            const dictionary& dict,
            const word& libsEntry,
            const TablePtr& tablePtr
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "dlLibraryTableTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
