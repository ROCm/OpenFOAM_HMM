/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "foamVtkLagrangianWriter.H"
#include "IOField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::vtk::lagrangianWriter::writeIOField
(
    const wordList& objectNames
)
{
    const int nCmpt(pTraits<Type>::nComponents);

    const bool useIntField =
        std::is_integral<typename pTraits<Type>::cmptType>();

    for (const word& fldName : objectNames)
    {
        IOobject header
        (
            fldName,
            mesh_.time().timeName(),
            cloud::prefix/cloudName_,
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false  // no register
        );

        IOField<Type> fld(header);

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
