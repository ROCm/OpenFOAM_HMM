/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2017 OpenCFD Ltd.
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
    Foam::functionObjects::regionFunctionObject

Description
    Specialization of Foam::functionObject for a region and providing a
    reference to the region Foam::objectRegistry.
    Also provides support for referencing a sub-region, which is typically
    needed when dealing with surfMesh and their fields.

See also
    Foam::functionObjects::stateFunctionObject

SourceFiles
    regionFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_regionFunctionObject_H
#define functionObjects_regionFunctionObject_H

#include "stateFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                    Class regionFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class regionFunctionObject
:
    public stateFunctionObject
{

protected:

    // Protected member data

        //- Reference to the region objectRegistry
        const objectRegistry& obr_;

        //- Optional reference to the sub-region objectRegistry.
        //  If a sub-region is not in effect, this reference is identical
        //  to the usual region objectRegistry.
        const objectRegistry& subObr_;


    // Protected member functions

        //- Selector for alternative sub-registry,
        //  when the keyword %subRegion is present in the dictionary
        static const objectRegistry& whichSubRegistry
        (
            const objectRegistry& obr,
            const dictionary& dict
        );

        //- The region or sub-region registry being used
        virtual const objectRegistry& obr() const;


        //- Find object (eg, a field) in the (sub) objectRegistry
        template<class ObjectType>
        bool foundObject(const word& fieldName) const;

        //- Lookup and return object (eg, a field) from the (sub) objectRegistry
        template<class ObjectType>
        const ObjectType& lookupObject(const word& fieldName) const;

        //- Lookup and return object (eg, a field) from the (sub) objectRegistry
        template<class ObjectType>
        ObjectType& lookupObjectRef(const word& fieldName) const;

        //- Lookup and return pointer to the object,
        //  otherwise nullptr if the object was not found,
        //  or had the incorrect type.
        template<class ObjectType>
        const ObjectType* lookupObjectPtr(const word& fieldName) const;

        //- Lookup and return non-const pointer to the object,
        //  otherwise nullptr if the object was not found,
        //  or had the incorrect type.
        template<class ObjectType>
        ObjectType* lookupObjectRefPtr(const word& fieldName) const;

        //- Store the field in the (sub) objectRegistry under the given name
        //  Note: sets the fieldName to tfield().name() if not already set
        template<class ObjectType>
        bool store
        (
            word& fieldName,
            const tmp<ObjectType>& tfield,
            bool cacheable = false
        );

        //- Write field if present in the (sub) objectRegistry
        bool writeObject(const word& fieldName);

        //- Clear field from the (sub) objectRegistry if present
        bool clearObject(const word& fieldName);


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        regionFunctionObject(const regionFunctionObject&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const regionFunctionObject&) = delete;


public:

    //- Runtime type information
    TypeName("regionFunctionObject");


    // Constructors

        //- Construct from Time and dictionary.
        //  The region objectRegistry is looked-up runTime with the name
        //  looked-up from the dictionary (defaults to polyMesh::defaultRegion)
        regionFunctionObject
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Construct from the region objectRegistry and dictionary
        regionFunctionObject
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict
        );


    //- Destructor
    virtual ~regionFunctionObject();


    // Member Functions

        //- Read optional controls
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "regionFunctionObjectTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
