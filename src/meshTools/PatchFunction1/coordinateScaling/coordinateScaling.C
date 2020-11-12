/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2020 OpenCFD Ltd.
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

template<class Type>
Foam::coordinateScaling<Type>::coordinateScaling()
:
    coordSys_(nullptr),
    scale_(),
    active_(false)
{}


template<class Type>
Foam::coordinateScaling<Type>::coordinateScaling
(
    const objectRegistry& obr,
    const dictionary& dict
)
:
    coordSys_
    (
        dict.found(coordinateSystem::typeName_())
      ? coordinateSystem::New(obr, dict)
      : nullptr
    ),
    scale_(3),
    active_(bool(coordSys_))
{
    for (direction dir = 0; dir < vector::nComponents; ++dir)
    {
        const word key("scale" + Foam::name(dir+1));

        if (dict.found(key))
        {
            scale_.set(dir, Function1<Type>::New(key, dict));
            active_ = true;
        }
    }
}


template<class Type>
Foam::coordinateScaling<Type>::coordinateScaling(const coordinateScaling& rhs)
:
    coordSys_(rhs.coordSys_.clone()),
    scale_(rhs.scale_),
    active_(rhs.active_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::coordinateScaling<Type>::transform
(
    const pointField& pos,
    const Field<Type>& p0
) const
{
    auto tfld = tmp<Field<Type>>::New(p0);
    auto& fld = tfld.ref();

    if (coordSys_)
    {
        const vectorField local(coordSys_->localPosition(pos));
        for (direction dir = 0; dir < vector::nComponents; ++dir)
        {
            if (scale_.set(dir))
            {
                fld = cmptMultiply
                (
                    fld,
                    scale_[dir].value(local.component(dir))
                );
            }
        }

        return coordSys_->transform(pos, fld);
    }
    else if (scale_.size())
    {
        for (direction dir = 0; dir < vector::nComponents; ++dir)
        {
            if (scale_.set(dir))
            {
                fld = cmptMultiply
                (
                    fld,
                    scale_[dir].value(pos.component(dir))
                );
            }
        }
    }

    return tfld;
}


template<class Type>
void Foam::coordinateScaling<Type>::writeEntry(Ostream& os) const
{
    if (coordSys_)
    {
        coordSys_->writeEntry(coordinateSystem::typeName_(), os);
    }
    forAll(scale_, dir)
    {
        if (scale_.set(dir))
        {
            scale_[dir].writeData(os);
        }
    }
}


// ************************************************************************* //
