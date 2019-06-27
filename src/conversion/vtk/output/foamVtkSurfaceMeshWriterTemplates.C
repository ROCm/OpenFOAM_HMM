/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2018 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011-2016 OpenFOAM Foundation
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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::vtk::surfaceMeshWriter::getFaceField
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& sfld
) const
{
    const polyBoundaryMesh& patches = sfld.mesh().boundaryMesh();

    const labelList& faceAddr = this->patch().addressing();

    auto tfld = tmp<Field<Type>>::New(faceAddr.size());
    auto iter = tfld.ref().begin();

    for (const label facei : faceAddr)
    {
        const label patchi = patches.whichPatch(facei);

        if (patchi == -1)
        {
            *iter = sfld[facei];
        }
        else
        {
            const label localFacei = facei - patches[patchi].start();
            *iter = sfld.boundaryField()[patchi][localFacei];
        }

        ++iter;
    }

    return tfld;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::vtk::surfaceMeshWriter::write
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& field
)
{
    if (notState(outputState::CELL_DATA))
    {
        FatalErrorInFunction
            << "Bad writer state (" << stateNames[state_]
            << ") - should be (" << stateNames[outputState::CELL_DATA]
            << ") for field " << field.name() << nl << endl
            << exit(FatalError);
    }

    this->indirectPatchWriter::write(field.name(), getFaceField(field)());
}


template<class Type>
void Foam::vtk::surfaceMeshWriter::write
(
    const GeometricField<Type, faPatchField, areaMesh>& field
)
{
    if (notState(outputState::CELL_DATA))
    {
        FatalErrorInFunction
            << "Bad writer state (" << stateNames[state_]
            << ") - should be (" << stateNames[outputState::CELL_DATA]
            << ") for field " << field.name() << nl << endl
            << exit(FatalError);
    }

    this->indirectPatchWriter::write(field.name(), field.primitiveField());
}


// ************************************************************************* //
