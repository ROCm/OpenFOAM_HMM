/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2019 OpenCFD Ltd.
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

#include "mapDistribute.H"
#include "globalIndexAndTransform.H"
#include "transformField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mapDistribute, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<>
void Foam::mapDistribute::transform::operator()
(
    const vectorTensorTransform&,
    const bool,
    List<label>&
) const {}

template<>
void Foam::mapDistribute::transform::operator()
(
    const coupledPolyPatch&,
    UList<label>&
) const {}

template<>
void Foam::mapDistribute::transform::operator()
(
    const coupledPolyPatch&,
    Map<label>&
) const {}

template<>
void Foam::mapDistribute::transform::operator()
(
    const coupledPolyPatch&,
    EdgeMap<label>&
) const {}


template<>
void Foam::mapDistribute::transform::operator()
(
    const vectorTensorTransform&,
    const bool,
    List<scalar>&
) const {}

template<>
void Foam::mapDistribute::transform::operator()
(
    const coupledPolyPatch&,
    UList<scalar>&
) const {}

template<>
void Foam::mapDistribute::transform::operator()
(
    const coupledPolyPatch&,
    Map<scalar>&
) const {}

template<>
void Foam::mapDistribute::transform::operator()
(
    const coupledPolyPatch&,
    EdgeMap<scalar>&
) const {}


template<>
void Foam::mapDistribute::transform::operator()
(
    const vectorTensorTransform&,
    const bool,
    List<bool>&
) const {}

template<>
void Foam::mapDistribute::transform::operator()
(
    const coupledPolyPatch&,
    UList<bool>&
) const {}

template<>
void Foam::mapDistribute::transform::operator()
(
    const coupledPolyPatch&,
    Map<bool>&
) const {}

template<>
void Foam::mapDistribute::transform::operator()
(
    const coupledPolyPatch&,
    EdgeMap<bool>&
) const {}


void Foam::mapDistribute::printLayout(Ostream& os) const
{
    mapDistributeBase::printLayout(os);

    forAll(transformElements_, trafoI)
    {
        if (transformElements_[trafoI].size() > 0)
        {
            os  << "transform " << trafoI << ':' << endl
                << "    start : " << transformStart_[trafoI] << endl
                << "    size  : " << transformElements_[trafoI].size() << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mapDistribute::mapDistribute(const label comm)
:
    mapDistributeBase(comm)
{}


Foam::mapDistribute::mapDistribute(const mapDistribute& map)
:
    mapDistributeBase(map),
    transformElements_(map.transformElements_),
    transformStart_(map.transformStart_)
{}


Foam::mapDistribute::mapDistribute(mapDistribute&& map)
:
    mapDistribute()
{
    transfer(map);
}


Foam::mapDistribute::mapDistribute
(
    const label constructSize,
    labelListList&& subMap,
    labelListList&& constructMap,
    const bool subHasFlip,
    const bool constructHasFlip,
    const label comm
)
:
    mapDistributeBase
    (
        constructSize,
        std::move(subMap),
        std::move(constructMap),
        subHasFlip,
        constructHasFlip,
        comm
    )
{}


Foam::mapDistribute::mapDistribute
(
    const label constructSize,
    labelListList&& subMap,
    labelListList&& constructMap,
    labelListList&& transformElements,
    labelList&& transformStart,
    const bool subHasFlip,
    const bool constructHasFlip,
    const label comm
)
:
    mapDistributeBase
    (
        constructSize,
        std::move(subMap),
        std::move(constructMap),
        subHasFlip,
        constructHasFlip,
        comm
    ),
    transformElements_(std::move(transformElements)),
    transformStart_(std::move(transformStart))
{}


Foam::mapDistribute::mapDistribute
(
    const labelUList& sendProcs,
    const labelUList& recvProcs,
    const label comm
)
:
    mapDistributeBase(sendProcs, recvProcs, comm)
{}


Foam::mapDistribute::mapDistribute
(
    const globalIndex& globalNumbering,
    labelList& elements,
    List<Map<label>>& compactMap,
    const int tag,
    const label comm
)
:
    mapDistributeBase
    (
        globalNumbering,
        elements,
        compactMap,
        tag,
        comm
    )
{}


Foam::mapDistribute::mapDistribute
(
    const globalIndex& globalNumbering,
    labelListList& cellCells,
    List<Map<label>>& compactMap,
    const int tag,
    const label comm
)
:
    mapDistributeBase
    (
        globalNumbering,
        cellCells,
        compactMap,
        tag,
        comm
    )
{}


Foam::mapDistribute::mapDistribute
(
    const globalIndex& globalNumbering,
    labelList& elements,
    const globalIndexAndTransform& globalTransforms,
    const labelPairList& transformedElements,
    labelList& transformedIndices,
    List<Map<label>>& compactMap,
    const int tag,
    const label comm
)
:
    mapDistributeBase(comm)
{
    const label myRank = Pstream::myProcNo(comm);

    // Construct per processor compact addressing of the global elements
    // needed. The ones from the local processor are not included since
    // these are always all needed.
    calcCompactAddressing
    (
        globalNumbering,
        elements,
        compactMap
    );

    // Add all (non-local) transformed elements needed.
    forAll(transformedElements, i)
    {
        labelPair elem = transformedElements[i];
        label proci = globalTransforms.processor(elem);
        if (proci != myRank)
        {
            label index = globalTransforms.index(elem);
            label nCompact = compactMap[proci].size();
            compactMap[proci].insert(index, nCompact);
        }
    }


    // Exchange what I need with processor that supplies it. Renumber elements
    // into compact numbering
    labelList compactStart;
    exchangeAddressing
    (
        tag,
        globalNumbering,
        elements,
        compactMap,
        compactStart
    );


    // Renumber the transformed elements
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Count per transformIndex
    label nTrafo = globalTransforms.transformPermutations().size();
    labelList nPerTransform(nTrafo, Zero);
    forAll(transformedElements, i)
    {
        labelPair elem = transformedElements[i];
        label trafoI = globalTransforms.transformIndex(elem);
        nPerTransform[trafoI]++;
    }
    // Offset per transformIndex
    transformStart_.setSize(nTrafo);
    transformElements_.setSize(nTrafo);
    forAll(transformStart_, trafoI)
    {
        transformStart_[trafoI] = constructSize_;
        constructSize_ += nPerTransform[trafoI];
        transformElements_[trafoI].setSize(nPerTransform[trafoI]);
    }

    // Sort transformed elements into their new slot.
    nPerTransform = 0;

    transformedIndices.setSize(transformedElements.size());
    forAll(transformedElements, i)
    {
        labelPair elem = transformedElements[i];
        label proci = globalTransforms.processor(elem);
        label index = globalTransforms.index(elem);
        label trafoI = globalTransforms.transformIndex(elem);

        // Get compact index for untransformed element
        label rawElemI =
        (
            proci == myRank
          ? index
          : compactMap[proci][index]
        );

        label& n = nPerTransform[trafoI];
        // index of element to transform
        transformElements_[trafoI][n] = rawElemI;
        // destination of transformed element
        transformedIndices[i] = transformStart_[trafoI]+n;
        n++;
    }

    if (debug)
    {
        printLayout(Pout);
    }
}


Foam::mapDistribute::mapDistribute
(
    const globalIndex& globalNumbering,
    labelListList& cellCells,
    const globalIndexAndTransform& globalTransforms,
    const List<labelPairList>& transformedElements,
    labelListList& transformedIndices,
    List<Map<label>>& compactMap,
    const int tag,
    const label comm
)
:
    mapDistributeBase(comm)
{
    const label myRank = Pstream::myProcNo(comm_);

    // Construct per processor compact addressing of the global elements
    // needed. The ones from the local processor are not included since
    // these are always all needed.
    calcCompactAddressing
    (
        globalNumbering,
        cellCells,
        compactMap
    );

    // Add all (non-local) transformed elements needed.
    forAll(transformedElements, celli)
    {
        const labelPairList& elems = transformedElements[celli];

        forAll(elems, i)
        {
            label proci = globalTransforms.processor(elems[i]);
            if (proci != myRank)
            {
                label index = globalTransforms.index(elems[i]);
                label nCompact = compactMap[proci].size();
                compactMap[proci].insert(index, nCompact);
            }
        }
    }


    // Exchange what I need with processor that supplies it. Renumber elements
    // into compact numbering
    labelList compactStart;
    exchangeAddressing
    (
        tag,
        globalNumbering,
        cellCells,
        compactMap,
        compactStart
    );


    // Renumber the transformed elements
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Count per transformIndex
    label nTrafo = globalTransforms.transformPermutations().size();
    labelList nPerTransform(nTrafo, Zero);
    forAll(transformedElements, celli)
    {
        const labelPairList& elems = transformedElements[celli];

        forAll(elems, i)
        {
            label trafoI = globalTransforms.transformIndex(elems[i]);
            nPerTransform[trafoI]++;
        }
    }
    // Offset per transformIndex
    transformStart_.setSize(nTrafo);
    transformElements_.setSize(nTrafo);
    forAll(transformStart_, trafoI)
    {
        transformStart_[trafoI] = constructSize_;
        constructSize_ += nPerTransform[trafoI];
        transformElements_[trafoI].setSize(nPerTransform[trafoI]);
    }

    // Sort transformed elements into their new slot.
    nPerTransform = 0;

    transformedIndices.setSize(transformedElements.size());
    forAll(transformedElements, celli)
    {
        const labelPairList& elems = transformedElements[celli];
        transformedIndices[celli].setSize(elems.size());

        forAll(elems, i)
        {
            label proci = globalTransforms.processor(elems[i]);
            label index = globalTransforms.index(elems[i]);
            label trafoI = globalTransforms.transformIndex(elems[i]);

            // Get compact index for untransformed element
            label rawElemI =
            (
                proci == myRank
              ? index
              : compactMap[proci][index]
            );

            label& n = nPerTransform[trafoI];
            // index of element to transform
            transformElements_[trafoI][n] = rawElemI;
            // destination of transformed element
            transformedIndices[celli][i] = transformStart_[trafoI]+n;
            n++;
        }
    }

    if (debug)
    {
        printLayout(Pout);
    }
}


Foam::mapDistribute::mapDistribute
(
    labelListList&& subMap,
    const bool subHasFlip,
    const bool constructHasFlip,
    const label comm
)
:
    mapDistributeBase(std::move(subMap), subHasFlip, constructHasFlip, comm)
{}


Foam::mapDistribute::mapDistribute(Istream& is)
{
    is  >> *this;
}


Foam::autoPtr<Foam::mapDistribute> Foam::mapDistribute::clone() const
{
    return autoPtr<mapDistribute>::New(*this);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::label Foam::mapDistribute::whichTransform(const label index) const
{
    return findLower(transformStart_, index+1);
}


void Foam::mapDistribute::transfer(mapDistribute& rhs)
{
    if (this == &rhs)
    {
        // Self-assignment is a no-op
    }

    mapDistributeBase::transfer(rhs);
    transformElements_.transfer(rhs.transformElements_);
    transformStart_.transfer(rhs.transformStart_);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::mapDistribute::operator=(const mapDistribute& rhs)
{
    if (this == &rhs)
    {
        return;  // Self-assignment is a no-op
    }

    mapDistributeBase::operator=(rhs);
    transformElements_ = rhs.transformElements_;
    transformStart_ = rhs.transformStart_;
}


void Foam::mapDistribute::operator=(mapDistribute&& rhs)
{
    if (this != &rhs)
    {
        // Avoid self-assignment
        transfer(rhs);
    }
}


// * * * * * * * * * * * * * * Istream Operator  * * * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, mapDistribute& map)
{
    is.fatalCheck(FUNCTION_NAME);

    is  >> static_cast<mapDistributeBase&>(map)
        >> map.transformElements_ >> map.transformStart_;

    return is;
}


// * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const mapDistribute& map)
{
    os  << static_cast<const mapDistributeBase&>(map) << token::NL
        << map.transformElements_ << token::NL
        << map.transformStart_;

    return os;
}


// ************************************************************************* //
