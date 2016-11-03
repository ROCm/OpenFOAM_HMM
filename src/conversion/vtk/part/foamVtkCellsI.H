/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
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

#include "foamVtkCells.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline bool Foam::foamVtkCells::decomp::requested() const
{
    return requested_;
}


inline bool Foam::foamVtkCells::decomp::used() const
{
    return !superCells_.empty();
}


inline const Foam::labelList&
Foam::foamVtkCells::decomp::superCells() const
{
    return superCells_;
}


inline const Foam::labelList&
Foam::foamVtkCells::decomp::addPointCellLabels() const
{
    return addPointCellLabels_;
}


inline const Foam::labelList&
Foam::foamVtkCells::decomp::pointMap() const
{
    return pointMap_;
}


inline bool Foam::foamVtkCells::decomposeRequested() const
{
    return decompose_.requested();
}


inline const Foam::labelList&
Foam::foamVtkCells::vertLabels() const
{
    return vertLabels_;
}


inline const Foam::labelList&
Foam::foamVtkCells::vertOffsets() const
{
    return vertOffset_;
}


inline const Foam::List<uint8_t>&
Foam::foamVtkCells::cellTypes() const
{
    return cellTypes_;
}


inline const Foam::labelList&
Foam::foamVtkCells::faceLabels() const
{
    return faceLabels_;
}


inline const Foam::labelList&
Foam::foamVtkCells::faceOffsets() const
{
    return faceOffset_;
}


inline const Foam::labelList&
Foam::foamVtkCells::addPointCellLabels() const
{
    return decompose_.addPointCellLabels();
}


inline const Foam::labelList&
Foam::foamVtkCells::superCells() const
{
    return decompose_.superCells();
}


inline Foam::label
Foam::foamVtkCells::nFieldCells() const
{
    return cellTypes_.size();
}


inline Foam::foamVtkCells::const_iterator
Foam::foamVtkCells::begin() const
{
    return const_iterator(*this);
}


inline Foam::foamVtkCells::const_iterator
Foam::foamVtkCells::end() const
{
    return const_iterator(*this, true);
}


// * * * * * * * * * * * * * * * * Iterators * * * * * * * * * * * * * * * * //

inline Foam::foamVtkCells::const_iterator::const_iterator
(
    const foamVtkCells& cells,
    bool isEnd
)
:
    parent_(cells),
    index_(0),
    begVert_(0),
    begFace_(0),
    legacy_()
{
    if (isEnd)
    {
        index_ = parent_.vertOffsets().size();
    }
}


inline
Foam::foamVtkCells::const_iterator&
Foam::foamVtkCells::const_iterator::operator++()
{
    ++index_;
    legacy_.clear();

    return *this;
}


inline
const Foam::UList<Foam::label>&
Foam::foamVtkCells::const_iterator::legacyCell() const
{
    if
    (
        legacy_.valid()
     || index_ >= parent_.vertOffsets().size()
    )
    {
        return legacy_();
    }

    const label endVert = parent_.vertOffsets()[index_];

    const label endFace =
    (
        parent_.faceOffsets().size()
      ? parent_.faceOffsets()[index_]
      : -1
    );

    if (endFace > 0)
    {
        // poly with face stream

        legacy_.reset
        (
            new SubList<label>
            (
                parent_.faceLabels(),
                endFace - begFace_,
                begFace_
            )
        );

        begFace_ = endFace;
    }
    else
    {
        // primitive without face stream
        legacy_.reset
        (
            new SubList<label>
            (
                parent_.vertLabels(),
                endVert - begVert_,
                begVert_
            )
        );
    }

    begVert_ = endVert;

    return legacy_();
}


inline bool
Foam::foamVtkCells::const_iterator::operator!=
(
    const const_iterator& rhs
) const
{
    return (index_ != rhs.index_);
}


// ************************************************************************* //
