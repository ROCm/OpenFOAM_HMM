/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    Selects a cell set through a dictionary.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "Time.H"
#include "fvMesh.H"
#include "topoSetSource.H"
#include "cellSet.H"
#include "volFields.H"

using namespace Foam;

template<class Type>
bool setFieldType
(
    const word& fieldTypeDesc,
    const fvMesh& mesh,
    const labelList& selectedCells,
    Istream& fieldValueStream
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    if (fieldTypeDesc != fieldType::typeName + "Value")
    {
        return false;
    }

    word fieldName(fieldValueStream);

    IOobject fieldHeader
    (
        fieldName,
        mesh.time().timeName(),
        mesh,
        IOobject::MUST_READ
    );

    // Check field exists
    if (fieldHeader.headerOk())
    {
        Info<< "    Setting " << fieldHeader.headerClassName()
            << " " << fieldName << endl;

        fieldType field(fieldHeader, mesh);

        const Type& value = pTraits<Type>(fieldValueStream);

        if (selectedCells.size() == field.size())
        {
            field.internalField() = value;
        }
        else
        {
            forAll(selectedCells, celli)
            {
                field[selectedCells[celli]] = value;
            }
        }

        forAll(field.boundaryField(), patchi)
        {
            field.boundaryField()[patchi] =
                field.boundaryField()[patchi].patchInternalField();
        }

        field.write();
    }
    else
    {
        WarningIn
        (
            "void setFieldType"
            "(const fvMesh& mesh, const labelList& selectedCells,"
            "Istream& fieldValueStream)"
        ) << "Field " << fieldName << " not found" << endl;
    }

    return true;
}


class setField
{

public:

    setField()
    {}

    autoPtr<setField> clone() const
    {
        return autoPtr<setField>(new setField());
    }

    class iNew
    {
        const fvMesh& mesh_;
        const labelList& selectedCells_;

    public:

        iNew(const fvMesh& mesh, const labelList& selectedCells)
        :
            mesh_(mesh),
            selectedCells_(selectedCells)
        {}

        autoPtr<setField> operator()(Istream& fieldValues) const
        {
            word fieldType(fieldValues);

            if
            (
               !(
                    setFieldType<scalar>
                        (fieldType, mesh_, selectedCells_, fieldValues)
                 || setFieldType<vector>
                        (fieldType, mesh_, selectedCells_, fieldValues)
                 || setFieldType<sphericalTensor>
                        (fieldType, mesh_, selectedCells_, fieldValues)
                 || setFieldType<symmTensor>
                        (fieldType, mesh_, selectedCells_, fieldValues)
                 || setFieldType<tensor>
                        (fieldType, mesh_, selectedCells_, fieldValues)
                )
            )
            {
                WarningIn("setField::iNew::operator()(Istream& is)")
                    << "field type " << fieldType << " not currently supported"
                    << endl;
            }

            return autoPtr<setField>(new setField());
        }
    };
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();

#   include "setRootCase.H"
#   include "createTime.H"

    // Get times list
    instantList timeDirs = timeSelector::select0(runTime, args);

#   include "createMesh.H"

    Info<< "Reading setFieldsDict\n" << endl;

    IOdictionary setFieldsDict
    (
        IOobject
        (
            "setFieldsDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    if (setFieldsDict.found("defaultFieldValues"))
    {
        Info<< "Setting field default values" << endl;
        PtrList<setField> defaultFieldValues
        (
            setFieldsDict.lookup("defaultFieldValues"),
            setField::iNew(mesh, labelList(mesh.nCells()))
        );
        Info<< endl;
    }


    Info<< "Setting field region values" << endl;

    PtrList<entry> regions(setFieldsDict.lookup("regions"));

    forAll(regions, regionI)
    {
        const entry& region = regions[regionI];

        autoPtr<topoSetSource> cellSelector =
            topoSetSource::New(region.keyword(), mesh, region.dict());

        cellSet selectedCellSet
        (
            mesh,
            "cellSet",
            mesh.nCells()/10+1  // Reasonable size estimate.
        );

        cellSelector->applyToSet
        (
            topoSetSource::NEW,
            selectedCellSet
        );

        PtrList<setField> fieldValues
        (
            region.dict().lookup("fieldValues"),
            setField::iNew(mesh, selectedCellSet.toc())
        );
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
