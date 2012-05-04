/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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

#include "CodedSource.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "dynamicCode.H"
#include "dynamicCodeContext.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::CodedSource<Type>::prepare
(
    dynamicCode& dynCode,
    const dynamicCodeContext& context
) const
{
    word sourceType(pTraits<Type>::typeName);

    // Set additional rewrite rules
    dynCode.setFilterVariable("typeName", redirectType_);
    dynCode.setFilterVariable("TemplateType", sourceType);
    dynCode.setFilterVariable("SourceType", sourceType + "Source");

    //dynCode.removeFilterVariable("code");
    dynCode.setFilterVariable("codeCorrect", codeCorrect_);
    dynCode.setFilterVariable("codeAddSup", codeAddSup_);
    dynCode.setFilterVariable("codeSetValue", codeSetValue_);

    // compile filtered C template
    dynCode.addCompileFile("codedBasicSourceTemplate.C");

    // copy filtered H template
    dynCode.addCopyFile("codedBasicSourceTemplate.H");

    // debugging: make BC verbose
    //         dynCode.setFilterVariable("verbose", "true");
    //         Info<<"compile " << redirectType_ << " sha1: "
    //             << context.sha1() << endl;

    // define Make/options
    dynCode.setMakeOptions
        (
            "EXE_INC = -g \\\n"
            "-I$(LIB_SRC)/fieldSources/lnInclude \\\n"
            "-I$(LIB_SRC)/finiteVolume/lnInclude \\\n"
            "-I$(LIB_SRC)/meshTools/lnInclude \\\n"
            "-I$(LIB_SRC)/sampling/lnInclude \\\n"
            + context.options()
            + "\n\nLIB_LIBS = \\\n"
            + "    -lmeshTools \\\n"
            + "    -lfieldSources \\\n"
            + "    -lsampling \\\n"
            + "    -lfiniteVolume \\\n"
            + context.libs()
        );
}


template<class Type>
Foam::dlLibraryTable& Foam::CodedSource<Type>::libs() const
{
    return const_cast<Time&>(mesh_.time()).libs();
}


template<class Type>
Foam::string Foam::CodedSource<Type>::description() const
{
    return "basicSource " + name_;
}


template<class Type>
void Foam::CodedSource<Type>::clearRedirect() const
{
    redirectBasicSourcePtr_.clear();
}


template<class Type>
const Foam::dictionary& Foam::CodedSource<Type>::codeDict() const
{
    return coeffs_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::CodedSource<Type>::CodedSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    basicSource(name, modelType, dict, mesh)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::basicSource& Foam::CodedSource<Type>::redirectBasicSource()
const
{
    if (!redirectBasicSourcePtr_.valid())
    {
        dictionary constructDict(dict_);
        constructDict.set("type", redirectType_);

        redirectBasicSourcePtr_ = basicSource::New
        (
            redirectType_,
            constructDict,
            mesh_
        );
    }
    return redirectBasicSourcePtr_();
}


template<class Type>
void Foam::CodedSource<Type>::correct
(
    GeometricField<Type, fvPatchField, volMesh>& fld
)
{
    if (debug)
    {
        Info<< "CodedSource<"<< pTraits<Type>::typeName
            << ">::correct for source " << name_ << endl;
    }

    updateLibrary(redirectType_);
    redirectBasicSource().correct(fld);
}


template<class Type>
void Foam::CodedSource<Type>::addSup
(
    fvMatrix<Type>& eqn,
    const label fieldI
)
{
    if (debug)
    {
        Info<< "CodedSource<"<< pTraits<Type>::typeName
            << ">::addSup for source " << name_ << endl;
    }

    updateLibrary(redirectType_);
    redirectBasicSource().addSup(eqn, fieldI);
}


template<class Type>
void Foam::CodedSource<Type>::setValue
(
    fvMatrix<Type>& eqn,
    const label fieldI
)
{
    if (debug)
    {
        Info<< "CodedSource<"<< pTraits<Type>::typeName
            << ">::setValue for source " << name_ << endl;
    }

    updateLibrary(redirectType_);
    redirectBasicSource().setValue(eqn, fieldI);
}


// ************************************************************************* //
