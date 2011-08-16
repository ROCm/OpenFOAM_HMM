/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::codedFunctionObject

Description
    functionObject using dynamic code compilation.

SourceFiles
    codedFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef codedFunctionObject_H
#define codedFunctionObject_H

#include "pointFieldFwd.H"
#include "functionObject.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class mapPolyMesh;
class dynamicCode;
class dynamicCodeContext;
class IOdictionary;

/*---------------------------------------------------------------------------*\
                       Class codedFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class codedFunctionObject
:
    public functionObject
{
protected:

    // Protected data

        //- Reference to the time database
        const Time& time_;

        //- Input dictionary
        dictionary dict_;

        word redirectType_;

        string codeRead_;
        string codeExecute_;
        string codeEnd_;

        //- Previously loaded library
        mutable fileName oldLibPath_;

        //- Underlying functionObject
        mutable autoPtr<functionObject> redirectFunctionObjectPtr_;

    // Private Member Functions

        //- Global loader/unloader function type
        typedef void (*loaderFunctionType)(bool);

        //- Load specified library and execute globalFuncName(true)
        void* loadLibrary
        (
            const fileName& libPath,
            const string& globalFuncName,
            const dictionary& contextDict
        ) const;

        //- Execute globalFuncName(false) and unload specified library
        void unloadLibrary
        (
            const fileName& libPath,
            const string& globalFuncName,
            const dictionary& contextDict
        ) const;


        //- Create library based on the dynamicCodeContext
        void createLibrary(dynamicCode&, const dynamicCodeContext&) const;

        //- Update library as required
        void updateLibrary() const;

        //- Read relevant dictionary entries
        void readDict();

        //- Disallow default bitwise copy construct
        codedFunctionObject(const codedFunctionObject&);

        //- Disallow default bitwise assignment
        void operator=(const codedFunctionObject&);


public:

    //- Runtime type information
    TypeName("coded");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        codedFunctionObject
        (
            const word& name,
            const Time& time,
            const dictionary& dict
        );


    //- Destructor
    virtual ~codedFunctionObject();


    // Member Functions

        //- Dynamically compiled functionObject
        functionObject& redirectFunctionObject() const;

        //- Called at the start of the time-loop
        virtual bool start();

        //- Called at each ++ or += of the time-loop. forceWrite overrides the
        //  outputControl behaviour.
        virtual bool execute(const bool forceWrite);

        //- Called when Time::run() determines that the time-loop exits.
        //  By default it simply calls execute().
        virtual bool end();

        //- Read and set the function object if its data have changed
        virtual bool read(const dictionary&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
