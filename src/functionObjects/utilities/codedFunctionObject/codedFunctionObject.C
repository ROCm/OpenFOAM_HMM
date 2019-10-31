/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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
#include "dynamicCode.H"
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
    dynCode.addCompileFile(codeTemplateC);

    // Copy filtered H template
    dynCode.addCopyFile(codeTemplateH);

    // Debugging: make verbose
    // dynCode.setFilterVariable("verbose", "true");
    // DetailInfo
    //     <<"compile " << name_ << " sha1: "
    //     << context.sha1() << endl;

    // Define Make/options
    dynCode.setMakeOptions
    (
        "EXE_INC = -g \\\n"
        "-I$(LIB_SRC)/finiteVolume/lnInclude \\\n"
        "-I$(LIB_SRC)/meshTools/lnInclude \\\n"
      + context.options()
      + "\n\nLIB_LIBS = \\\n"
        "    -lOpenFOAM \\\n"
        "    -lfiniteVolume \\\n"
        "    -lmeshTools \\\n"
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
    const Time& runTime,
    const dictionary& dict
)
:
    timeFunctionObject(name, runTime),
    codedBase(),
    dict_(dict)
{
    read(dict_);

    updateLibrary(name_);
    redirectFunctionObject();
}


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
    timeFunctionObject::read(dict);

    codedBase::setCodeContext(dict);

    dict.readCompat<word>("name", {{"redirectType", 1706}}, name_);

    label nKeywords = 0;

    const entry* eptr;

    codeData_.clear();
    codedBase::append("<codeData>");
    if ((eptr = dict.findEntry("codeData", keyType::LITERAL)) != nullptr)
    {
        eptr->readEntry(codeData_);
        dynamicCodeContext::inplaceExpand(codeData_, dict);
        codedBase::append(codeData_);

        dynamicCodeContext::addLineDirective
        (
            codeData_,
            eptr->startLineNumber(),
            dict.name()
        );

        ++nKeywords;
    }

    codeRead_.clear();
    codedBase::append("<codeRead>");
    if ((eptr = dict.findEntry("codeRead", keyType::LITERAL)) != nullptr)
    {
        eptr->readEntry(codeRead_);
        dynamicCodeContext::inplaceExpand(codeRead_, dict);
        codedBase::append(codeRead_);

        dynamicCodeContext::addLineDirective
        (
            codeRead_,
            eptr->startLineNumber(),
            dict.name()
        );

        ++nKeywords;
    }

    codeExecute_.clear();
    codedBase::append("<codeExecute>");
    if ((eptr = dict.findEntry("codeExecute", keyType::LITERAL)) != nullptr)
    {
        eptr->readEntry(codeExecute_);
        dynamicCodeContext::inplaceExpand(codeExecute_, dict);
        codedBase::append(codeExecute_);

        dynamicCodeContext::addLineDirective
        (
            codeExecute_,
            eptr->startLineNumber(),
            dict.name()
        );

        ++nKeywords;
    }

    codeWrite_.clear();
    codedBase::append("<codeWrite>");
    if ((eptr = dict.findEntry("codeWrite", keyType::LITERAL)) != nullptr)
    {
        eptr->readEntry(codeWrite_);
        dynamicCodeContext::inplaceExpand(codeWrite_, dict);
        codedBase::append(codeWrite_);

        dynamicCodeContext::addLineDirective
        (
            codeWrite_,
            eptr->startLineNumber(),
            dict.name()
        );

        ++nKeywords;
    }

    codeEnd_.clear();
    codedBase::append("<codeEnd>");
    if ((eptr = dict.findEntry("codeEnd", keyType::LITERAL)) != nullptr)
    {
        eptr->readEntry(codeEnd_);
        dynamicCodeContext::inplaceExpand(codeEnd_, dict);
        codedBase::append(codeEnd_);

        dynamicCodeContext::addLineDirective
        (
            codeEnd_,
            eptr->startLineNumber(),
            dict.name()
        );

        ++nKeywords;
    }

    if (!nKeywords)
    {
        IOWarningInFunction(dict)
            << "No critical \"code\" prefixed keywords found." << nl
            << "Please check the code documentation for more details." << nl
            << endl;
    }

    updateLibrary(name_);
    return redirectFunctionObject().read(dict);
}


// ************************************************************************* //
