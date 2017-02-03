/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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
    Foam::functionObjectState

Description
    Base class for function objects, adding functionality to read/write state
    information (data required for smooth restart behaviour) and results
    to/from the state dictionary

    Note: cannot access the state dictionary until after construction of the
    function objects, since the owner container functionObjectList is owned
    by time, and time owns the state dictionary.  I.e. need to wait for time
    to be fully constructed.

See also
    Foam::functionObject

SourceFiles
    stateFunctionObject.C
    stateFunctionObjectTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_stateFunctionObject_H
#define functionObjects_stateFunctionObject_H

#include "functionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class IOdictionary;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                     Class stateFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class stateFunctionObject
:
    public functionObject
{

private:

    // Private member data

        //- Name of the results dictionary
        static const word resultsName_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        stateFunctionObject(const stateFunctionObject&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const stateFunctionObject&) = delete;


protected:

    // Protected Member Data

        //- Reference to the time database
        const Time& time_;


    // Protacted Member Functions

        //- Return a const reference to the state dictionary
        const IOdictionary& stateDict() const;

        //- Return non-const access to the state dictionary
        IOdictionary& stateDict();



public:

    // Constructors

        //- Construct from components
        stateFunctionObject(const word& name, const Time& runTime);


    //- Destructor
    ~stateFunctionObject();


    // Member Functions

        //- Return access to the property dictionary
        dictionary& propertyDict();


        // Properties

            //- Return true if the property exists
            bool foundProperty(const word& entryName) const;

            //- Set dictionary, return true if set
            bool getDict
            (
                const word& entryName,
                dictionary& dict
            ) const;

            //- Set dictionary from named object, return true if set
            bool getObjectDict
            (
                const word& objectName,
                const word& entryName,
                dictionary& dict
            ) const;

            //- Retrieve generic property
            template<class Type>
            Type getProperty
            (
                const word& entryName,
                const Type& defaultValue = Type(Zero)
            ) const;

            //- Set generic property, return true if set
            template<class Type>
            bool getProperty(const word& entryName, Type& value) const;

            //- Add generic property
            template<class Type>
            void setProperty(const word& entryName, const Type& value);

            //- Retrieve generic property from named object
            template<class Type>
            Type getObjectProperty
            (
                const word& objectName,
                const word& entryName,
                const Type& defaultValue = Type(Zero)
            ) const;

            //- Set generic property from named object, return true if set
            template<class Type>
            bool getObjectProperty
            (
                const word& objectName,
                const word& entryName,
                Type& value
            ) const;

            //- Add generic property from named object
            template<class Type>
            void setObjectProperty
            (
                const word& objectName,
                const word& entryName,
                const Type& value
            );


        // Results

            //- Add result
            template<class Type>
            void setResult
            (
                const word& entryName,
                const Type& value
            );

            //- Add result from named object
            template<class Type>
            void setObjectResult
            (
                const word& objectName,
                const word& entryName,
                const Type& value
            );

            //- Retrieve result
            template<class Type>
            Type getResult
            (
                const word& entryName,
                const Type& defaultValue = Type(Zero)
            ) const;

            //- Retrieve result from named object
            template<class Type>
            Type getObjectResult
            (
                const word& objectName,
                const word& entryName,
                const Type& defaultValue = Type(Zero)
            ) const;

            //- Set result from named object, return true if set
            template<class Type>
            bool getObjectResult
            (
                const word& objectName,
                const word& entryName,
                Type& value
            ) const;

            //- Retrieve the result type
            word resultType(const word& entryName) const;

            //- Return the type of result
            word objectResultType
            (
                const word& objectName,
                const word& entryName
            ) const;

            //- Retrieve the result entries
            List<word> objectResultEntries() const;

            //- Return result entries for named object
            List<word> objectResultEntries(const word& objectName) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "stateFunctionObjectTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
