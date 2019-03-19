/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2019 OpenCFD Ltd.
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

#include "foamVtkPatchWriter.H"
#include "foamVtkOutput.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::vtk::patchWriter::writeUniform
(
    const word& fieldName,
    const Type& val
)
{
    if (isState(outputState::CELL_DATA))
    {
        ++nCellData_;
        vtk::fileWriter::writeUniform<Type>(fieldName, val, numberOfCells_);
    }
    else if (isState(outputState::POINT_DATA))
    {
        ++nPointData_;
        vtk::fileWriter::writeUniform<Type>(fieldName, val, numberOfPoints_);
    }
    else
    {
        WarningInFunction
            << "Ignore bad writer state (" << stateNames[state_]
            << ") for field " << fieldName << nl << endl
            << exit(FatalError);
    }
}


template<class Type, template<class> class PatchField>
void Foam::vtk::patchWriter::write
(
    const GeometricField<Type, PatchField, pointMesh>& field
)
{
    if (isState(outputState::POINT_DATA))
    {
        ++nPointData_;
    }
    else
    {
        FatalErrorInFunction
            << "Bad writer state (" << stateNames[state_]
            << ") - should be (" << stateNames[outputState::POINT_DATA]
            << ") for field " << field.name() << nl << endl
            << exit(FatalError);
    }

    const direction nCmpt(pTraits<Type>::nComponents);

    label nPoints = nLocalPoints_;

    if (parallel_)
    {
        reduce(nPoints, sumOp<label>());
    }


    if (format_)
    {
        if (legacy())
        {
            legacy::floatField<nCmpt>(format(), field.name(), nPoints);
        }
        else
        {
            const uint64_t payLoad =
                vtk::sizeofData<float, nCmpt>(nPoints);

            format().beginDataArray<float, nCmpt>(field.name());
            format().writeSize(payLoad);
        }
    }


    if (parallel_ ? Pstream::master() : true)
    {
        for (const label patchId : patchIDs_)
        {
            const auto& pfld = field.boundaryField()[patchId];

            vtk::writeList(format(), pfld.patchInternalField()());
        }
    }


    if (parallel_)
    {
        // Patch Ids are identical across all processes
        const label nPatches = patchIDs_.size();

        if (Pstream::master())
        {
            Field<Type> recv;

            // Receive each patch field and write
            for
            (
                int slave=Pstream::firstSlave();
                slave<=Pstream::lastSlave();
                ++slave
            )
            {
                IPstream fromSlave(Pstream::commsTypes::blocking, slave);

                for (label i=0; i < nPatches; ++i)
                {
                    fromSlave >> recv;

                    vtk::writeList(format(), recv);
                }
            }
        }
        else
        {
            // Send each patch field to master
            OPstream toMaster
            (
                Pstream::commsTypes::blocking,
                Pstream::masterNo()
            );

            for (const label patchId : patchIDs_)
            {
                const auto& pfld = field.boundaryField()[patchId];

                toMaster << pfld.patchInternalField()();
            }
        }
    }


    if (format_)
    {
        format().flush();
        format().endDataArray();
    }
}


template<class Type, template<class> class PatchField>
void Foam::vtk::patchWriter::write
(
    const GeometricField<Type, PatchField, volMesh>& field
)
{
    if (isState(outputState::CELL_DATA))
    {
        ++nCellData_;
    }
    else
    {
        FatalErrorInFunction
            << "Bad writer state (" << stateNames[state_]
            << ") - should be (" << stateNames[outputState::CELL_DATA]
            << ") for field " << field.name() << nl << endl
            << exit(FatalError);
    }

    const direction nCmpt(pTraits<Type>::nComponents);

    label nFaces = nLocalFaces_;

    if (parallel_)
    {
        reduce(nFaces, sumOp<label>());
    }


    if (format_)
    {
        if (legacy())
        {
            legacy::floatField<nCmpt>(format(), field.name(), nFaces);
        }
        else
        {
            const uint64_t payLoad =
                vtk::sizeofData<float, nCmpt>(nFaces);

            format().beginDataArray<float, nCmpt>(field.name());
            format().writeSize(payLoad);
        }
    }


    if (parallel_ ? Pstream::master() : true)
    {
        for (const label patchId : patchIDs_)
        {
            const auto& pfld = field.boundaryField()[patchId];

            if (useNearCellValue_)
            {
                vtk::writeList(format(), pfld.patchInternalField()());
            }
            else
            {
                vtk::writeList(format(), pfld);
            }
        }
    }

    if (parallel_)
    {
        // Patch Ids are identical across all processes
        const label nPatches = patchIDs_.size();

        if (Pstream::master())
        {
            Field<Type> recv;

            // Receive each patch field and write
            for
            (
                int slave=Pstream::firstSlave();
                slave<=Pstream::lastSlave();
                ++slave
            )
            {
                IPstream fromSlave(Pstream::commsTypes::blocking, slave);

                for (label i=0; i < nPatches; ++i)
                {
                    fromSlave >> recv;

                    vtk::writeList(format(), recv);
                }
            }
        }
        else
        {
            // Send each patch field to master
            OPstream toMaster
            (
                Pstream::commsTypes::blocking,
                Pstream::masterNo()
            );

            for (const label patchId : patchIDs_)
            {
                const auto& pfld = field.boundaryField()[patchId];

                if (useNearCellValue_)
                {
                    toMaster << pfld.patchInternalField()();
                }
                else
                {
                    toMaster << static_cast<const Field<Type>&>(pfld);
                }
            }
        }
    }


    if (format_)
    {
        format().flush();
        format().endDataArray();
    }
}


template<class Type>
void Foam::vtk::patchWriter::write
(
    const GeometricField<Type, fvPatchField, volMesh>& field,
    const PrimitivePatchInterpolation<primitivePatch>& pInter
)
{
    if (isState(outputState::POINT_DATA))
    {
        ++nPointData_;
    }
    else
    {
        FatalErrorInFunction
            << "Bad writer state (" << stateNames[state_]
            << ") - should be (" << stateNames[outputState::POINT_DATA]
            << ") for field " << field.name() << nl << endl
            << exit(FatalError);
    }

    const direction nCmpt(pTraits<Type>::nComponents);

    label nPoints = nLocalPoints_;

    if (parallel_)
    {
        reduce(nPoints, sumOp<label>());
    }


    if (format_)
    {
        if (legacy())
        {
            legacy::floatField<nCmpt>(format(), field.name(), nPoints);
        }
        else
        {
            const uint64_t payLoad =
                vtk::sizeofData<float, nCmpt>(nPoints);

            format().beginDataArray<float, nCmpt>(field.name());
            format().writeSize(payLoad);
        }
    }


    if (parallel_ ? Pstream::master() : true)
    {
        for (const label patchId : patchIDs_)
        {
            const auto& pfld = field.boundaryField()[patchId];

            if (useNearCellValue_)
            {
                auto tfield =
                    pInter.faceToPointInterpolate
                    (
                        pfld.patchInternalField()()
                    );

                vtk::writeList(format(), tfield());
            }
            else
            {
                auto tfield = pInter.faceToPointInterpolate(pfld);

                vtk::writeList(format(), tfield());
            }
        }
    }


    if (parallel_)
    {
        // Patch Ids are identical across all processes
        const label nPatches = patchIDs_.size();

        if (Pstream::master())
        {
            Field<Type> recv;

            // Receive each patch field and write
            for
            (
                int slave=Pstream::firstSlave();
                slave<=Pstream::lastSlave();
                ++slave
            )
            {
                IPstream fromSlave(Pstream::commsTypes::blocking, slave);

                for (label i=0; i < nPatches; ++i)
                {
                    fromSlave >> recv;

                    vtk::writeList(format(), recv);
                }
            }
        }
        else
        {
            // Send each patch field to master
            OPstream toMaster
            (
                Pstream::commsTypes::blocking,
                Pstream::masterNo()
            );

            for (const label patchId : patchIDs_)
            {
                const auto& pfld = field.boundaryField()[patchId];

                if (useNearCellValue_)
                {
                    auto tfield =
                        pInter.faceToPointInterpolate
                        (
                            pfld.patchInternalField()()
                        );

                    toMaster << tfield();
                }
                else
                {
                    auto tfield = pInter.faceToPointInterpolate(pfld);

                    toMaster << tfield();
                }
            }
        }
    }


    if (format_)
    {
        format().flush();
        format().endDataArray();
    }
}


// ************************************************************************* //
