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

#include "ReactionList.H"
#include "IFstream.H"
#include "SLPtrList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::ReactionList<ThermoType>::ReactionList
(
    const speciesTable& species,
    const HashPtrTable<ThermoType>& thermoDb
)
:
    PtrList<Reaction<ThermoType> >(),
    species_(species),
    thermoDb_(thermoDb),
    dict_(dictionary::null)
{}


template<class ThermoType>
Foam::ReactionList<ThermoType>::ReactionList
(
    const speciesTable& species,
    const HashPtrTable<ThermoType>& thermoDb,
    const dictionary& dict
)
:
    PtrList<Reaction<ThermoType> >(),
    species_(species),
    thermoDb_(thermoDb),
    dict_(dict)
{
    readReactionDict();
}


template<class ThermoType>
Foam::ReactionList<ThermoType>::ReactionList
(
    const speciesTable& species,
    const HashPtrTable<ThermoType>& thermoDb,
    const fileName& fName
)
:
    PtrList<Reaction<ThermoType> >
    (
        dictionary(IFstream(fName)()).lookup("reactions"),
        Reaction<ThermoType>::iNew(species, thermoDb)
    ),
    species_(species),
    thermoDb_(thermoDb),
    dict_(dictionary::null)
{}


template<class ThermoType>
Foam::ReactionList<ThermoType>::ReactionList(const ReactionList& reactions)
:
    PtrList<Reaction<ThermoType> >(reactions),
    species_(reactions.species_),
    thermoDb_(reactions.thermoDb_),
    dict_(reactions.dict_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::ReactionList<ThermoType>::~ReactionList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
bool Foam::ReactionList<ThermoType>::readReactionDict()
{
    const dictionary& reactions(dict_.subDict("reactions"));
    PtrList<Reaction<ThermoType> > newPtrs(reactions.size());

    label i = 0;
    forAllConstIter(dictionary, reactions, iter)
    {
        newPtrs.set
        (
            i++,
            Reaction<ThermoType>::New
            (
                species_,
                thermoDb_,
                reactions.subDict(iter().keyword())
            )
        );
    }

    PtrList<Reaction<ThermoType> >::transfer(newPtrs);

    return true;
}


template<class ThermoType>
const Foam::PtrList<Foam::Reaction<ThermoType> >&
Foam::ReactionList<ThermoType>::reactions() const
{
    return *this;
}


// ************************************************************************* //
