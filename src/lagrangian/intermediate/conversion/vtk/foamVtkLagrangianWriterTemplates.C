/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2019 OpenCFD Ltd.
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

#include "foamVtkLagrangianWriter.H"
#include "IOobjectList.H"
#include "IOField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::vtk::lagrangianWriter::write(const IOField<Type>& field)
{
    if (isState(outputState::CELL_DATA))
    {
        ++nCellData_;
    }
    else if (isState(outputState::POINT_DATA))
    {
        ++nPointData_;
    }
    else
    {
        FatalErrorInFunction
            << "Bad writer state (" << stateNames[state_]
            << ") - should be (" << stateNames[outputState::CELL_DATA]
            << ") or (" << stateNames[outputState::POINT_DATA]
            << ") for field " << field.name() << nl << endl
            << exit(FatalError);
    }

    static_assert
    (
        (
            std::is_same<label, typename pTraits<Type>::cmptType>::value
         || std::is_floating_point<typename pTraits<Type>::cmptType>::value
        ),
        "Label and Floating-point vector space only"
    );


    // Other integral types (eg, bool etc) would need cast/convert to label.
    // Similarly for labelVector etc.

    // // Ensure consistent output width
    // for (const Type& val : field)
    // {
    //     for (int cmpt=0; cmpt < nCmpt; ++cmpt)
    //     {
    //         format().write(label(component(val, cmpt)));
    //     }
    // }


    const direction nCmpt(pTraits<Type>::nComponents);

    label nVals = field.size();

    if (parallel_)
    {
        reduce(nVals, sumOp<label>());
    }

    if (format_)
    {
        // Non-legacy

        if (std::is_same<label, typename pTraits<Type>::cmptType>::value)
        {
            const uint64_t payLoad = vtk::sizeofData<label, nCmpt>(nVals);

            format().beginDataArray<label, nCmpt>(field.name());
            format().writeSize(payLoad);
        }
        else
        {
            const uint64_t payLoad = vtk::sizeofData<float, nCmpt>(nVals);

            format().beginDataArray<float, nCmpt>(field.name());
            format().writeSize(payLoad);
        }
    }

    if (parallel_)
    {
        vtk::writeListParallel(format_.ref(), field);
    }
    else
    {
        vtk::writeList(format(), field);
    }


    if (format_)
    {
        format().flush();
        format().endDataArray();
    }
}


template<class Type>
Foam::label Foam::vtk::lagrangianWriter::writeFields
(
    const wordList& fieldNames,
    bool verbose
)
{
    const fileName localDir(cloudDir());

    label nFields = 0;

    for (const word& fieldName : fieldNames)
    {
        // Globally the field is expected to exist (MUST_READ), but can
        // be missing on a local processor.
        //
        // However, constructing IOField with MUST_READ and valid=false fails.
        // Workaround: READ_IF_PRESENT and verify the header globally

        IOobject io
        (
            fieldName,
            mesh_.time().timeName(),
            localDir,
            mesh_,
            IOobject::READ_IF_PRESENT
        );

        // Check global existence to avoid any errors
        const bool ok =
            returnReduce
            (
                io.typeHeaderOk<IOField<Type>>(false),
                orOp<bool>()
            );

        if (!ok)
        {
            continue;
        }

        if (verbose)
        {
            if (!nFields)
            {
                Info<< "        " << pTraits<Type>::typeName << ":";
            }

            Info<< " " << fieldName;
        }

        IOField<Type> field(io);
        this->write<Type>(field);

        ++nFields;
    }

    if (verbose && nFields)
    {
        Info << endl;
    }

    return nFields;
}


template<class Type>
Foam::label Foam::vtk::lagrangianWriter::writeFields
(
    const IOobjectList& objects,
    const bool verbose
)
{
    return writeFields<Type>
    (
        objects.allNames<IOField<Type>>(),  // These are in sorted order
        verbose
    );
}


// ************************************************************************* //
