/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2022 OpenCFD Ltd.
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

Application
    Test-objectRegistry2

Description
    Print objectRegistry information, with some additional tests.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvCFD.H"
#include "fvMesh.H"
#include "volFields.H"
#include "timeSelector.H"
#include "ReadFields.H"
#include "IOstreams.H"
#include "PtrListOps.H"
#include "IOobjectList.H"
#include "objectRegistry.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
bool loadField(fvMesh& mesh, const word& fieldName)
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> SurfaceFieldType;

    if (mesh.objectRegistry::found(fieldName))
    {
        // Info<< fieldName << " already in database" << endl;
        return false;
    }

    IOobject fieldHeader
    (
        fieldName,
        mesh.time().timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (fieldHeader.typeHeaderOk<VolFieldType>(true, true, false))
    {
        // Store field on mesh database
        VolFieldType* ptr = new VolFieldType(fieldHeader, mesh);
        mesh.objectRegistry::store(ptr);
        return true;
    }
    else if (fieldHeader.typeHeaderOk<SurfaceFieldType>(true, true, false))
    {
        // Store field on mesh database
        SurfaceFieldType* ptr = new SurfaceFieldType(fieldHeader, mesh);
        mesh.objectRegistry::store(ptr);
        return true;
    }

    return false;
}


bool loadField(fvMesh& mesh, const word& fieldName)
{
    return
    (
        !mesh.objectRegistry::found(fieldName)
    &&
        (
            loadField<scalar>(mesh, fieldName)
         || loadField<vector>(mesh, fieldName)
         || loadField<sphericalTensor>(mesh, fieldName)
         || loadField<symmTensor>(mesh, fieldName)
         || loadField<tensor>(mesh, fieldName)
        )
    );
}


void loadFields(fvMesh& mesh, const IOobjectList& objects)
{
    for (const word& fieldName : objects.names())
    {
        loadField(mesh, fieldName);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void report(const UPtrList<const Type>& objects)
{
    Info<< Type::typeName << " name/type:" << nl
        << objects.size() << nl << '(' << nl;

    for (const Type& obj : objects)
    {
        Info<< "  " << obj.name() << " : " << obj.type() << nl;
    }

    Info<< ')' << nl << endl;
}


void printRegistry
(
    Foam::Ostream& os,
    const Foam::objectRegistry& obr,
    Foam::label indent = 4
);


void printRegistry
(
    Foam::Ostream& os,
    const Foam::objectRegistry& obr,
    Foam::label indent
)
{
    UPtrList<const regIOobject> objects(obr.sorted());
    wordList regNames(obr.sortedNames<objectRegistry>());

    std::string prefix;
    for (label i=indent; i; --i)
    {
        prefix += ' ';
    }

    os  << '#' << prefix.c_str() << obr.name()
        << " parent:" << obr.parent().name() << nl;

    os  << ' ' << prefix.c_str() << "objects: "
        << flatOutput(PtrListOps::names(objects)) << nl;
    os  << ' ' << prefix.c_str() << "registries: "
        << flatOutput(regNames) << nl;


    // Print without expanding sub-registries
    for (const regIOobject& obj : objects)
    {
        os  << (isA<objectRegistry>(obj) ? '-' : ' ')
            << prefix.c_str() << obj.name() << " => " << obj.type() << nl;
    }
    for (label i=indent; i; --i)
    {
        os  << '-'; // divider
    }
    os  << '\n';

    // Now descend into the sub-registries
    for (const word& name : regNames)
    {
        const objectRegistry& next = obr.lookupObject<objectRegistry>
        (
            name,
            false // non-recursive
        );

        os  << prefix.c_str()
            << "current:" << obr.name() << " next:"
            << next.name() << " next-parent:" << next.parent().name() << nl;

        os  << prefix.c_str() << name << " => " << obr[name]->type();

        if ("dictionary" == obr[name]->type())
        {
            os  << " (skip dictionary)" << nl;
        }
        else
        {
            os  << nl;
            printRegistry(os, next, indent + 4);
        }
    }
}


template<class Type>
void filterTest(const objectRegistry& obr, const wordRe& re)
{
    Info<< nl << "Filter on names:" << nl;

    Info<< "Filter = " << re << nl;

    const word& typeName = Type::typeName;

    Info<< "    <" << typeName <<">(" << re << ") : "
        << obr.count<Type>(re) << nl
        << "    (" << typeName << "::typeName, " << re << ") : "
        << obr.count(typeName, re) << nl;

    Info<< "    <" << typeName << ">(" << re << ") : "
        << flatOutput(obr.sortedNames<Type>(re)) << nl
        // << flatOutput(obr.names<Type>(re)) << nl
        << "    (" << typeName << "::typeName, " << re << ") : "
        << flatOutput(obr.sortedNames(typeName, re)) << nl
        //<< flatOutput(obr.names(typeName, re)) << nl
        ;


    wordRe reClass("vol.*Field", wordRe::REGEX);
    wordRe re2(re, wordRe::REGEX_ICASE);

    Info<< "General" << nl
        << "    <void>(" << re << ") : "
        << flatOutput(obr.sortedNames<void>(re)) << nl
        << "    (" << reClass << ", " << re2 <<" ignore-case) : "
        << flatOutput(obr.sortedNames(reClass, re2)) << nl
        ;

    Info<< nl;
}


void registryTests(const objectRegistry& obr)
{
    Info<< nl << "Registry: " << obr.name() << nl
        << " names: " << flatOutput(obr.sortedNames()) << nl;

    Info<< "count" << nl
        << "    <void>()    : " << obr.count<void>() << nl
        << "    <labelList>()   : " << obr.count<labelList>() << nl
        << "    <labelList>(strict) : " << obr.count<labelList>(true) << nl
        << "    <scalarList>()   : " << obr.count<scalarList>() << nl
        << "    <scalarList>(strict) : " << obr.count<scalarList>(true) << nl;
    Info<< "    <volScalarField>()    : "
        << obr.count<volScalarField>() << nl
        << "    (volScalarField::typeName) : "
        << obr.count(volScalarField::typeName) << nl;
    Info<< "    <volVectorField>()    : "
        << obr.count<volVectorField>() << nl
        << "    (volVectorField::typeName) : "
        << obr.count(volVectorField::typeName) << nl;

    Info<< nl << "Filter on names:" << nl;

    filterTest<volScalarField>(obr, wordRe("[p-z].*", wordRe::DETECT));

    Info<< nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:

int main(int argc, char *argv[])
{
    argList::noBanner();
    //argList::noParallel();
//    argList::addOption
//    (
//        "filter",
//        "wordRes",
//        "filter keys with names or regexs"
//    );

    // timeSelector::addOptions();
    timeSelector::addOptions(true, true);

    #include "setRootCase.H"

//    wordRes matcher;
//    if (args.readListIfPresent<wordRe>("filter", matcher))
//    {
//        Info<<"limit names: " << matcher << nl;
//    }

    #include "createTime.H"
    #include "createMesh.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Info<< "Time: " << runTime.timeName() << endl;

        // Read objects in time directory
        IOobjectList objects(mesh, runTime.timeName());

        // Read volFields
        loadFields(mesh, objects);

        printRegistry(Info, mesh);

        registryTests(mesh);

        report(mesh.sorted<const volScalarField>());
        report(mesh.csorted<volVectorField>());

        Info<< nl;
    }


    Info<<"\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
