/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM, licensed under GNU General Public License
    <http://www.gnu.org/licenses/>.

Namespace
    Foam::cellModeller

Description
    Compatibility definitions of static cellModel lookups.

    Superseded (NOV-2017) by cellModel methods.

\*---------------------------------------------------------------------------*/
#ifndef cellModeller_H
#define cellModeller_H

#include "cellModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace cellModeller
{

    //- Deprecated(2017-11) equivalent to cellModel::ptr static method.
    //  \deprecated(2017-11) use cellModel::ptr instead
    inline const cellModel* lookup(const word& modelName)
    {
        return cellModel::ptr(modelName);
    }

    //- Deprecated(2017-11) equivalent to cellModel::ptr static method.
    //  \deprecated(2017-11) use cellModel::ptr instead
    inline const cellModel* lookup(const label modelIndex)
    {
        return cellModel::ptr(modelIndex);
    }

} // End namespace cellModeller
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
