/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

InClass
    Foam::basicMixture

Description
    Macros for creating 'basic' mixtures for basic fluid thermo packages

\*---------------------------------------------------------------------------*/

#ifndef makeBasicMixture_H
#define makeBasicMixture_H

#include "basicMixture.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeBasicMixture(Mixture,Transport,Type,Thermo,EqnOfState,Specie)     \
                                                                              \
typedef                                                                       \
    Mixture<Transport<species::thermo<Thermo<EqnOfState<Specie> >, Type> > >  \
    Mixture##Transport##Type##Thermo##EqnOfState##Specie;                     \
                                                                              \
defineTemplateTypeNameAndDebugWithName                                        \
    (Mixture##Transport##Type##Thermo##EqnOfState##Specie,                    \
    #Mixture"<"#Transport"<"#Thermo"<"#EqnOfState"<"#Specie">>,"#Type">>", 0)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
