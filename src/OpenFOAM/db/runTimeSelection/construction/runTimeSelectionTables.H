/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

Description
    Macros to ease declaration of run-time selection tables.

    declareRunTimeSelectionTable is used to create a run-time selection table
    for a base-class which holds constructor pointers on the table.

    declareRunTimeNewSelectionTable is used to create a run-time selection
    table for a derived-class which holds "New" pointers on the table.

\*---------------------------------------------------------------------------*/

#include "token.H"

#ifndef runTimeSelectionTables_H
#define runTimeSelectionTables_H

#include "autoPtr.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Common infrastructure

// Not used directly: declare run-time selection (variables and methods)
#define declareRunTimeSelectionTableBase(returnType,prefix,argList)            \
                                                                               \
    /* Function pointer type: construct from argList */                        \
    typedef returnType (*prefix##Ptr)argList;                                  \
                                                                               \
    /* Function table type: construct from argList */                          \
    typedef ::Foam::HashTable                                                  \
    <                                                                          \
        prefix##Ptr,                                                           \
        ::Foam::word,                                                          \
        ::Foam::string::hasher                                                 \
    > prefix##TableType;                                                       \
                                                                               \
    /* Table singleton (storage) */                                            \
    static prefix##TableType* prefix##TablePtr_;                               \
                                                                               \
    /* Table construct/destruct helper */                                      \
    static void prefix##TablePtr_construct(bool load)


// Not used directly: storage and helper methods for runtime tables
#define defineRunTimeSelectionTableBase(baseType,prefix,Tspecialize)           \
                                                                               \
    /* Define table singleton (storage) */                                     \
    Tspecialize prefix##TableType* prefix##TablePtr_(nullptr);                 \
                                                                               \
    /* Table construct/destruct helper */                                      \
    Tspecialize void prefix##TablePtr_construct(bool load)                     \
    {                                                                          \
        static bool constructed = false;                                       \
        if (load)                                                              \
        {                                                                      \
            if (!constructed)                                                  \
            {                                                                  \
                prefix##TablePtr_ = new prefix##TableType;                     \
                constructed = true;                                            \
            }                                                                  \
        }                                                                      \
        else if (prefix##TablePtr_)                                            \
        {                                                                      \
            delete prefix##TablePtr_;                                          \
            prefix##TablePtr_ = nullptr;                                       \
            constructed = false;                                               \
        }                                                                      \
    }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//
// Declaration Macros
//

//- Declare a run-time selection (variables and adder classes)
//  The ptrWrapper is typically a Foam::autoPtr or a Foam::tmp container.
#define declareRunTimeSelectionTable\
(ptrWrapper,baseType,argNames,argList,parList)                                 \
                                                                               \
    declareRunTimeSelectionTableBase(                                          \
        ptrWrapper<baseType>, argNames##Constructor, argList);                 \
                                                                               \
    /* Helper to add constructor from argList to table */                      \
    template<class baseType##Type>                                             \
    struct add##argNames##ConstructorToTable                                   \
    {                                                                          \
        static ptrWrapper<baseType> New argList                                \
        {                                                                      \
            return ptrWrapper<baseType>(new baseType##Type parList);           \
        }                                                                      \
                                                                               \
        explicit add##argNames##ConstructorToTable                             \
        (                                                                      \
            const ::Foam::word& k = baseType##Type::typeName                   \
        )                                                                      \
        {                                                                      \
            argNames##ConstructorTablePtr_construct(true);                     \
            if (!argNames##ConstructorTablePtr_->insert(k, New))               \
            {                                                                  \
                std::cerr                                                      \
                    << "Duplicate entry " << k << " in runtime table "         \
                    << #baseType << std::endl;                                 \
                ::Foam::error::safePrintStack(std::cerr);                      \
            }                                                                  \
        }                                                                      \
                                                                               \
        ~add##argNames##ConstructorToTable()                                   \
        {                                                                      \
            argNames##ConstructorTablePtr_construct(false);                    \
        }                                                                      \
                                                                               \
        add##argNames##ConstructorToTable                                      \
            (const add##argNames##ConstructorToTable&) = delete;               \
                                                                               \
        void operator=                                                         \
            (const add##argNames##ConstructorToTable&) = delete;               \
    };                                                                         \
                                                                               \
    /* Helper to add constructor from argList to table */                      \
    /* Remove only the entry (not the table) upon destruction */               \
    template<class baseType##Type>                                             \
    struct addRemovable##argNames##ConstructorToTable                          \
    {                                                                          \
        const ::Foam::word name;  /* Lookup name for later removal */          \
                                                                               \
        static ptrWrapper<baseType> New argList                                \
        {                                                                      \
            return ptrWrapper<baseType>(new baseType##Type parList);           \
        }                                                                      \
                                                                               \
        explicit addRemovable##argNames##ConstructorToTable                    \
        (                                                                      \
            const ::Foam::word& k = baseType##Type::typeName                   \
        )                                                                      \
        :                                                                      \
            name(k)                                                            \
        {                                                                      \
            argNames##ConstructorTablePtr_construct(true);                     \
            argNames##ConstructorTablePtr_->set(k, New);                       \
        }                                                                      \
                                                                               \
        ~addRemovable##argNames##ConstructorToTable()                          \
        {                                                                      \
            if (argNames##ConstructorTablePtr_)                                \
            {                                                                  \
                argNames##ConstructorTablePtr_->erase(name);                   \
            }                                                                  \
        }                                                                      \
                                                                               \
        addRemovable##argNames##ConstructorToTable                             \
            (const addRemovable##argNames##ConstructorToTable&) = delete;      \
                                                                               \
        void operator=                                                         \
            (const addRemovable##argNames##ConstructorToTable&) = delete;      \
    };


//- Declare a run-time selection for derived classes
#define declareRunTimeNewSelectionTable\
(ptrWrapper,baseType,argNames,argList,parList)                                 \
                                                                               \
    declareRunTimeSelectionTableBase(                                          \
        ptrWrapper<baseType>,argNames##Constructor,argList);                   \
                                                                               \
    /* Helper to add constructor from argList to table */                      \
    template<class baseType##Type>                                             \
    struct add##argNames##ConstructorToTable                                   \
    {                                                                          \
        static ptrWrapper<baseType> New##baseType argList                      \
        {                                                                      \
            return ptrWrapper<baseType>(baseType##Type::New parList.ptr());    \
        }                                                                      \
                                                                               \
        explicit add##argNames##ConstructorToTable                             \
        (                                                                      \
            const ::Foam::word& k = baseType##Type::typeName                   \
        )                                                                      \
        {                                                                      \
            argNames##ConstructorTablePtr_construct(true);                     \
            if (!argNames##ConstructorTablePtr_->insert(k, New##baseType))     \
            {                                                                  \
                std::cerr                                                      \
                    << "Duplicate entry " << k << " in runtime table "         \
                    << #baseType << std::endl;                                 \
                ::Foam::error::safePrintStack(std::cerr);                      \
            }                                                                  \
        }                                                                      \
                                                                               \
        ~add##argNames##ConstructorToTable()                                   \
        {                                                                      \
            argNames##ConstructorTablePtr_construct(false);                    \
        }                                                                      \
                                                                               \
        add##argNames##ConstructorToTable                                      \
            (const add##argNames##ConstructorToTable&) = delete;               \
                                                                               \
        void operator=                                                         \
            (const add##argNames##ConstructorToTable&) = delete;               \
    };                                                                         \
                                                                               \
    /* Helper to add constructor from argList to table */                      \
    template<class baseType##Type>                                             \
    struct addRemovable##argNames##ConstructorToTable                          \
    {                                                                          \
        const ::Foam::word name;  /* Retain name for later removal */          \
                                                                               \
        static ptrWrapper<baseType> New##baseType argList                      \
        {                                                                      \
            return ptrWrapper<baseType>(baseType##Type::New parList.ptr());    \
        }                                                                      \
                                                                               \
        explicit addRemovable##argNames##ConstructorToTable                    \
        (                                                                      \
            const ::Foam::word& k = baseType##Type::typeName                   \
        )                                                                      \
        :                                                                      \
            name(k)                                                            \
        {                                                                      \
            argNames##ConstructorTablePtr_construct(true);                     \
            argNames##ConstructorTablePtr_->set(k, New##baseType);             \
        }                                                                      \
                                                                               \
        ~addRemovable##argNames##ConstructorToTable()                          \
        {                                                                      \
            if (argNames##ConstructorTablePtr_)                                \
            {                                                                  \
                argNames##ConstructorTablePtr_->erase(name);                   \
            }                                                                  \
        }                                                                      \
                                                                               \
        addRemovable##argNames##ConstructorToTable                             \
            (const addRemovable##argNames##ConstructorToTable&) = delete;      \
                                                                               \
        void operator=                                                         \
            (const addRemovable##argNames##ConstructorToTable&) = delete;      \
    };


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//
// Definition Macros
//

//- Define run-time selection table
#define defineRunTimeSelectionTable(baseType,argNames)                         \
                                                                               \
    defineRunTimeSelectionTableBase(                                           \
        baseType,baseType::argNames##Constructor,)


//- Define run-time selection table for template classes
//  use when baseType doesn't need a template argument (eg, is a typedef)
#define defineTemplateRunTimeSelectionTable(baseType,argNames)                 \
                                                                               \
    defineRunTimeSelectionTableBase(                                           \
        baseType,baseType::argNames##Constructor,template<>)


//- Define run-time selection table for template classes
//  use when baseType requires the Targ template argument
#define defineTemplatedRunTimeSelectionTable(baseType,argNames,Targ)           \
                                                                               \
    defineRunTimeSelectionTableBase(                                           \
        baseType,baseType<Targ>::argNames##Constructor,template<>)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
