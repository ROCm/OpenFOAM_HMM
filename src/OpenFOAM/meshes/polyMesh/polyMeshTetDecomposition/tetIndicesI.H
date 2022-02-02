/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2021-2022 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline constexpr Foam::tetIndices::tetIndices() noexcept
:
    celli_(-1),
    facei_(-1),
    tetPti_(-1)
{}


inline constexpr Foam::tetIndices::tetIndices
(
    label celli,
    label facei,
    label tetPointi
) noexcept
:
    celli_(celli),
    facei_(facei),
    tetPti_(tetPointi)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::triFace Foam::tetIndices::faceTriIs
(
    const polyMesh& mesh,
    const bool warn
) const
{
    const Foam::face& f = mesh.faces()[face()];

    label faceBasePtI = mesh.tetBasePtIs()[face()];

    if (faceBasePtI < 0)
    {
        faceBasePtI = 0;

        if (warn && nWarnings_ < maxNWarnings)
        {
            ++nWarnings_;
            WarningInFunction
                << "No base point for face " << face() << ", " << f
                << ", produces a valid tet decomposition." << endl;
            if (nWarnings_ == maxNWarnings)
            {
                Warning
                    << "Suppressing any further warnings." << endl;
            }
        }
    }

    label facePtI = (tetPti_ + faceBasePtI) % f.size();
    label faceOtherPtI = f.fcIndex(facePtI);

    if (mesh.faceOwner()[face()] != cell())
    {
        std::swap(facePtI, faceOtherPtI);
    }

    return triFace(f[faceBasePtI], f[facePtI], f[faceOtherPtI]);
}


inline Foam::triFace Foam::tetIndices::triIs
(
    const polyMesh& mesh,
    const bool warn
) const
{
    const Foam::face& f = mesh.faces()[face()];

    label faceBasePtI = mesh.tetBasePtIs()[face()];

    if (faceBasePtI < 0)
    {
        faceBasePtI = 0;

        if (warn && nWarnings_ < maxNWarnings)
        {
            ++nWarnings_;
            WarningInFunction
                << "No base point for face " << face() << ", " << f
                << ", produces a valid tet decomposition." << endl;
            if (nWarnings_ == maxNWarnings)
            {
                Warning
                    << "Suppressing any further warnings." << endl;
            }
        }
    }

    label facePtI = (tetPti_ + faceBasePtI) % f.size();
    label faceOtherPtI = f.fcIndex(facePtI);

    if (mesh.faceOwner()[face()] != cell())
    {
        std::swap(facePtI, faceOtherPtI);
    }

    return triFace(faceBasePtI, facePtI, faceOtherPtI);
}


inline Foam::tetPointRef Foam::tetIndices::tet(const polyMesh& mesh) const
{
    const pointField& meshPoints = mesh.points();
    const triFace tri = faceTriIs(mesh);

    return tetPointRef
    (
        mesh.cellCentres()[cell()],
        meshPoints[tri[0]],
        meshPoints[tri[1]],
        meshPoints[tri[2]]
    );
}


inline Foam::triPointRef Foam::tetIndices::faceTri(const polyMesh& mesh) const
{
    const pointField& meshPoints = mesh.points();
    const triFace tri = faceTriIs(mesh);

    return triPointRef
    (
        meshPoints[tri[0]],
        meshPoints[tri[1]],
        meshPoints[tri[2]]
    );
}


inline Foam::triPointRef Foam::tetIndices::oldFaceTri
(
    const polyMesh& mesh
) const
{
    const pointField& meshOldPoints = mesh.oldPoints();
    const triFace tri = faceTriIs(mesh);

    return triPointRef
    (
        meshOldPoints[tri[0]],
        meshOldPoints[tri[1]],
        meshOldPoints[tri[2]]
    );
}


inline Foam::label Foam::tetIndices::compare
(
    const tetIndices& a,
    const tetIndices& b
) noexcept
{
    label diff;
    return
    (
        ((diff = (a.cell() - b.cell())) != 0) ? diff
      : ((diff = (a.face() - b.face())) != 0) ? diff
      : (a.tetPt() - b.tetPt())
    );
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

inline bool Foam::operator==(const tetIndices& a, const tetIndices& b) noexcept
{
    // Possibly slightly faster version than compare
    return
    (
        a.cell() == b.cell()
     && a.face() == b.face()
     && a.tetPt() == b.tetPt()
    );
}


inline bool Foam::operator!=(const tetIndices& a, const tetIndices& b) noexcept
{
    return !(a == b);
}


// ************************************************************************* //
