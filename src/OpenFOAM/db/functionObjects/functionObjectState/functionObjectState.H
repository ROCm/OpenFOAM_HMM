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

See Also
    Foam::functionObject

SourceFiles
    functionObjectState.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjectState_H
#define functionObjectState_H

#include "objectRegistry.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class functionObjectState Declaration
\*---------------------------------------------------------------------------*/

class functionObjectState
{
private:

    // Private data

        //- Name of the results dictionary
        static const word resultsName_;

        //- Reference to the database
        const objectRegistry& obr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        functionObjectState(const functionObjectState&);

        //- Disallow default bitwise assignment
        void operator=(const functionObjectState&);


protected:

    // Protected data

        //- Name of model
        const word name_;

        //- Flag to indicate whether the object is active
        bool active_;


    // Protacted Member Functions

        //- Return a const reference to the state dictionary
        const IOdictionary& stateDict() const;

        //- Return non-const access to the state dictionary
        IOdictionary& stateDict();



public:

    // Constructors

        //- Construct from components
        functionObjectState(const objectRegistry& obr, const word& name);


    //- Destructor
    virtual ~functionObjectState();


    // Member Functions

        //- Return the name
        const word& name() const;

        //- Return the active flag
        bool active() const;

        //- Return access to the property dictionary
        dictionary& propertyDict();

        //- Set the active status by querying objectRegistry type
        //  returns new active status
        template<class Type>
        bool setActive();


        // Properties

            //- Return true if the property exists
            bool foundProperty(const word& entryName) const;

            //- Retrieve generic property
            template<class Type>
            Type getProperty
            (
                const word& entryName,
                const Type& defaultValue = pTraits<Type>::zero
            ) const;

            //- Retrieve generic property
            template<class Type>
            void getProperty(const word& entryName, Type& value) const;

            //- Add generic property
            template<class Type>
            void setProperty(const word& entryName, const Type& value);

            //- Retrieve generic property from named object
            template<class Type>
            Type getObjectProperty
            (
                const word& objectName,
                const word& entryName,
                const Type& defaultValue = pTraits<Type>::zero
            ) const;

            //- Retrieve generic property from named object
            template<class Type>
            void getObjectProperty
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
                const Type& defaultValue = pTraits<Type>::zero
            ) const;

            //- Retrieve result from named object
            template<class Type>
            Type getObjectResult
            (
                const word& objectName,
                const word& entryName,
                const Type& defaultValue = pTraits<Type>::zero
            ) const;

            //- Retrieve result from named object
            template<class Type>
            void getObjectResult
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

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "functionObjectStateTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
