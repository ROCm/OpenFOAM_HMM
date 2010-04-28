/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2010 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::label Foam::globalIndexAndTransform::encodeTransformIndex
(
    const List<label>& permutationIndices
) const
{
    if (permutationIndices.size() != transforms_.size())
    {
        FatalErrorIn
        (
            "Foam::label encodeTransformIndex"
            "("
                "const List<label>& permutationIndices,"
            ") const"
        )
            << "permutationIndices " << permutationIndices
            << "are of a different size to the number of independent transforms"
            << abort(FatalError);
    }

    label transformIndex = 0;

    label w = 1;

    for (label b = 0; b < transforms_.size(); b++)
    {
        if (mag(permutationIndices[b]) > 1)
        {
            FatalErrorIn
            (
                "Foam::label encodeTransformIndex"
                "("
                "const List<label>& permutationIndices,"
                ") const"
            )
                << "permutationIndices " << permutationIndices
                << "are illegal, they must all be only -1, 0 or +1"
                << abort(FatalError);
        }

        transformIndex += (permutationIndices[b] + 1)*w;

        w *= 3;
    }

    return transformIndex;
}


Foam::labelPair Foam::globalIndexAndTransform::encode
(
    const label index,
    const label transformIndex
)
{
    return encode(Pstream::myProcNo(), index, transformIndex);
}


Foam::labelPair Foam::globalIndexAndTransform::encode
(
    const label procI,
    const label index,
    const label transformIndex
)
{
    if (transformIndex < 0 || transformIndex >= base_)
    {
        FatalErrorIn
        (
            "Foam::labelPair Foam::globalIndexAndTransform::encode"
            "("
                "const label procI, "
                "const label index, "
                "const label transformIndex"
            ")"
        )
            << "TransformIndex " << transformIndex
            << " is outside allowed range of 0 to "
            << base_ - 1
            << abort(FatalError);
    }

    if (procI > labelMax/base_)
    {
        FatalErrorIn
        (
            "Foam::labelPair Foam::globalIndexAndTransform::encode"
            "("
                "const label procI, "
                "const label index, "
                "const label transformIndex"
            ")"
        )
            << "Overflow : encoding processor " << procI << " in base " << base_
            << " exceeds capability of label (" << labelMax
            << "). Please recompile with larger datatype for label."
            << exit(FatalError);
    }

    return labelPair
    (
        index,
        transformIndex + procI*base_
    );
}


Foam::label  Foam::globalIndexAndTransform::index
(
    const labelPair& globalIAndTransform
)
{
    return globalIAndTransform.first();
}


Foam::label Foam::globalIndexAndTransform::processor
(
    const labelPair& globalIAndTransform
)
{
    return globalIAndTransform.second()/base_;
}


Foam::label Foam::globalIndexAndTransform::transformIndex
(
    const labelPair& globalIAndTransform
)
{
    return globalIAndTransform.second() % base_;
}


Foam::label Foam::globalIndexAndTransform::nIndependentTransforms() const
{
    return transforms_.size();
}


const Foam::List<Foam::vectorTensorTransform>&
Foam::globalIndexAndTransform::transforms() const
{
    return transforms_;
}


const Foam::List<Foam::vectorTensorTransform>&
Foam::globalIndexAndTransform::transformPermutations() const
{
    return transformPermutations_;
}


const Foam::vectorTensorTransform& Foam::globalIndexAndTransform::transform
(
    label transformIndex
) const
{
    return transformPermutations_[transformIndex];
}


// ************************************************************************* //
