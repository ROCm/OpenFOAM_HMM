/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2018 OpenCFD Ltd.
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
#include "IOField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::vtk::lagrangianWriter::writeIOField(const wordList& fieldNames)
{
    const int nCmpt(pTraits<Type>::nComponents);

    const bool useIntField =
        std::is_integral<typename pTraits<Type>::cmptType>();

    const fileName cloudDir(cloud::prefix/cloudName_);

    for (const word& fldName : fieldNames)
    {
        // Globally the field is expected to exist (MUST_READ), but can
        // be missing on a local processor.
        //
        // However, constructing IOField with MUST_READ and valid=false fails.
        // Workaround: READ_IF_PRESENT and verify the header globally

        IOobject fieldObject
        (
            fldName,
            mesh_.time().timeName(),
            cloudDir,
            mesh_,
            IOobject::READ_IF_PRESENT
        );

        // Check global existence - could make an error
        const bool fieldExists =
            returnReduce
            (
                fieldObject.typeHeaderOk<IOField<Type>>(false),
                orOp<bool>()
            );

        if (!fieldExists)
        {
            continue;
        }

        IOField<Type> fld(fieldObject);

        // NOTE: Could skip if there are no local parcels...

        if (useIntField)
        {
            const uint64_t payLoad(fld.size() * nCmpt * sizeof(label));

            if (legacy_)
            {
                legacy::intField(os(), fldName, nCmpt, fld.size());
            }
            else
            {
                format().openDataArray<label, nCmpt>(fldName)
                    .closeTag();
            }

            format().writeSize(payLoad);

            // Ensure consistent output width
            for (const Type& val : fld)
            {
                for (int cmpt=0; cmpt < nCmpt; ++cmpt)
                {
                    format().write(label(component(val, cmpt)));
                }
            }
        }
        else
        {
            const uint64_t payLoad(fld.size() * nCmpt * sizeof(float));

            if (legacy_)
            {
                legacy::floatField(os(), fldName, nCmpt, fld.size());
            }
            else
            {
                format().openDataArray<float, nCmpt>(fldName)
                    .closeTag();
            }

            format().writeSize(payLoad);
            vtk::writeList(format(), fld);
        }

        format().flush();

        if (!legacy_)
        {
            format().endDataArray();
        }
    }
}


// ************************************************************************* //
