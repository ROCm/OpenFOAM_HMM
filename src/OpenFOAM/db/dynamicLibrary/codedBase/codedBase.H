/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::codedBase

Description
    Base class for function objects and boundary conditions using dynamic code
    that provides methods for managing loading/unloading/updating
    of a dynamic library. For these purposes, it uses a dynamicCodeContext
    object to maintain information about the state.

    For simple coded objects, the default state management is sufficient.
    When there are more complicated code segements
    (eg, functionObjects::codedFunctionObject), the state management
    must also register these elements as well, starting with an initial
    setCodeContext() call and followed by append() to register each element.

SourceFiles
    codedBase.C

\*---------------------------------------------------------------------------*/

#ifndef codedBase_H
#define codedBase_H

#include "dictionary.H"
#include "dynamicCodeContext.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class Ostream;
class dynamicCode;
class dlLibraryTable;

/*---------------------------------------------------------------------------*\
                       Class codedBase Declaration
\*---------------------------------------------------------------------------*/

class codedBase
{
    // Private Data

        //- The code context
        dynamicCodeContext context_;

        //- Previously loaded library
        mutable fileName oldLibPath_;


    // Private Member Functions

        //- Global loader/unloader function type
        typedef void (*loaderFunctionType)(bool);

        //- Load specified library and execute globalFuncName(true)
        void* loadLibrary
        (
            const fileName& libPath,
            const string& globalFuncName,
            const dynamicCodeContext& context
        ) const;

        //- Execute globalFuncName(false) and unload specified library
        void unloadLibrary
        (
            const fileName& libPath,
            const string& globalFuncName,
            const dynamicCodeContext& context
        ) const;

        //- Create library based on the dynamicCodeContext
        void createLibrary
        (
            dynamicCode& dynCode,
            const dynamicCodeContext& context
        ) const;

        //- No copy construct
        codedBase(const codedBase&) = delete;

        //- No copy assignment
        void operator=(const codedBase&) = delete;


protected:

        //- Write code-dictionary contents
        static void writeCodeDict(Ostream& os, const dictionary& dict);


    // Protected Member Functions

        //- Set code context from a dictionary
        void setCodeContext(const dictionary& dict);

        //- Add content to SHA1 hashing
        void append(const std::string& str);


        //- Update library as required, using the given context
        void updateLibrary
        (
            const word& name,
            const dynamicCodeContext& context
        ) const;

        //- Update library as required, using the given code dictionary
        //- to use for the context
        void updateLibrary
        (
            const word& name,
            const dictionary& dict
        ) const;

        //- Update library as required, using the predefined context
        //- or use the codeDict() to generate one
        void updateLibrary(const word& name) const;

        //- Get the loaded dynamic libraries
        virtual dlLibraryTable& libs() const = 0;

        //- Adapt the context for the current object
        virtual void prepare
        (
            dynamicCode& dynCode,
            const dynamicCodeContext& context
        ) const = 0;

        // Return a description (type + name) for the output
        virtual string description() const = 0;

        // Clear any redirected objects
        virtual void clearRedirect() const = 0;

        // Get the dictionary to initialize the codeContext
        virtual const dictionary& codeDict() const = 0;


public:

    //- Runtime type information
    ClassName("codedBase");


    // Constructors

        //- Construct null
        codedBase() = default;


    //- Destructor
    virtual ~codedBase() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
