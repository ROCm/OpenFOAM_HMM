/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
     \\/     M anipulation  |
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

    //- Equivalent to cellModel::ptr static method.
    //  \deprecated use cellModel::ptr instead (NOV-2017)
    inline const cellModel* lookup(const word& modelName)
    {
        return cellModel::ptr(modelName);
    }

    //- Equivalent to cellModel::ptr static method.
    //  \deprecated use cellModel::ptr instead (NOV-2017)
    inline const cellModel* lookup(const label modelIndex)
    {
        return cellModel::ptr(modelIndex);
    }

} // End namespace cellModeller
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
