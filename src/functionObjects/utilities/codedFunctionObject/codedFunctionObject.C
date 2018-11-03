/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "codedFunctionObject.H"
#include "volFields.H"
#include "dictionary.H"
#include "Time.H"
#include "SHA1Digest.H"
#include "dynamicCode.H"
#include "dynamicCodeContext.H"
#include "stringOps.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(codedFunctionObject, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        codedFunctionObject,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::codedFunctionObject::prepare
(
    dynamicCode& dynCode,
    const dynamicCodeContext& context
) const
{
    // Set additional rewrite rules
    dynCode.setFilterVariable("typeName", name_);
    dynCode.setFilterVariable("codeData", codeData_);
    dynCode.setFilterVariable("codeRead", codeRead_);
    dynCode.setFilterVariable("codeExecute", codeExecute_);
    dynCode.setFilterVariable("codeWrite", codeWrite_);
    dynCode.setFilterVariable("codeEnd", codeEnd_);

    // Compile filtered C template
    dynCode.addCompileFile("functionObjectTemplate.C");

    // Copy filtered H template
    dynCode.addCopyFile("functionObjectTemplate.H");

    // Debugging: make BC verbose
    // dynCode.setFilterVariable("verbose", "true");
    // Info<<"compile " << name_ << " sha1: "
    //     << context.sha1() << endl;

    // Define Make/options
    dynCode.setMakeOptions
    (
        "EXE_INC = -g \\\n"
        "-I$(LIB_SRC)/finiteVolume/lnInclude \\\n"
        "-I$(LIB_SRC)/meshTools/lnInclude \\\n"
      + context.options()
      + "\n\nLIB_LIBS = \\\n"
      + "    -lOpenFOAM \\\n"
      + "    -lfiniteVolume \\\n"
      + "    -lmeshTools \\\n"
      + context.libs()
    );
}


Foam::dlLibraryTable& Foam::functionObjects::codedFunctionObject::libs() const
{
    return const_cast<Time&>(time_).libs();
}


Foam::string Foam::functionObjects::codedFunctionObject::description() const
{
    return "functionObject " + name();
}


void Foam::functionObjects::codedFunctionObject::clearRedirect() const
{
    redirectFunctionObjectPtr_.clear();
}


const Foam::dictionary&
Foam::functionObjects::codedFunctionObject::codeDict() const
{
    return dict_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::codedFunctionObject::codedFunctionObject
(
    const word& name,
    const Time& time,
    const dictionary& dict
)
:
    functionObject(name),
    codedBase(),
    time_(time),
    dict_(dict)
{
    read(dict_);

    updateLibrary(name_);
    redirectFunctionObject();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::codedFunctionObject::~codedFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::functionObject&
Foam::functionObjects::codedFunctionObject::redirectFunctionObject() const
{
    if (!redirectFunctionObjectPtr_.valid())
    {
        dictionary constructDict(dict_);
        constructDict.set("type", name_);

        redirectFunctionObjectPtr_ = functionObject::New
        (
            name_,
            time_,
            constructDict
        );
    }
    return *redirectFunctionObjectPtr_;
}


bool Foam::functionObjects::codedFunctionObject::execute()
{
    updateLibrary(name_);
    return redirectFunctionObject().execute();
}


bool Foam::functionObjects::codedFunctionObject::write()
{
    updateLibrary(name_);
    return redirectFunctionObject().write();
}


bool Foam::functionObjects::codedFunctionObject::end()
{
    updateLibrary(name_);
    return redirectFunctionObject().end();
}


bool Foam::functionObjects::codedFunctionObject::read(const dictionary& dict)
{
    functionObject::read(dict);

    dict.readCompat<word>("name", {{"redirectType", 1706}}, name_);

    const entry* dataPtr = dict.findEntry("codeData", keyType::LITERAL);

    if (dataPtr)
    {
        dataPtr->readEntry(codeData_);
        stringOps::inplaceExpand(codeData_, dict);
        dynamicCodeContext::addLineDirective
        (
            codeData_,
            dataPtr->startLineNumber(),
            dict.name()
        );
    }

    const entry* readPtr = dict.findEntry("codeRead", keyType::LITERAL);

    if (readPtr)
    {
        readPtr->readEntry(codeRead_);
        stringOps::inplaceExpand(codeRead_, dict);
        dynamicCodeContext::addLineDirective
        (
            codeRead_,
            readPtr->startLineNumber(),
            dict.name()
        );
    }

    const entry* execPtr = dict.findEntry("codeExecute", keyType::LITERAL);

    if (execPtr)
    {
        execPtr->readEntry(codeExecute_);
        stringOps::inplaceExpand(codeExecute_, dict);
        dynamicCodeContext::addLineDirective
        (
            codeExecute_,
            execPtr->startLineNumber(),
            dict.name()
        );
    }

    const entry* writePtr = dict.findEntry("codeWrite", keyType::LITERAL);

    if (writePtr)
    {
        writePtr->readEntry(codeWrite_);
        stringOps::inplaceExpand(codeWrite_, dict);
        dynamicCodeContext::addLineDirective
        (
            codeWrite_,
            writePtr->startLineNumber(),
            dict.name()
        );
    }

    const entry* endPtr = dict.findEntry("codeEnd", keyType::LITERAL);

    if (endPtr)
    {
        endPtr->readEntry(codeEnd_);
        stringOps::inplaceExpand(codeEnd_, dict);
        dynamicCodeContext::addLineDirective
        (
            codeEnd_,
            endPtr->startLineNumber(),
            dict.name()
        );
    }

    if (!dataPtr && !readPtr && !execPtr && !writePtr && !endPtr)
    {
        IOWarningInFunction(dict)
            << "No critical \"code\" prefixed keywords were found."
            << " Please check the code documentation for more details."
            << nl << endl;
    }

    updateLibrary(name_);
    return redirectFunctionObject().read(dict);
}


// ************************************************************************* //
