/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

    As a special exception, you have permission to link this program with the
    CGAL library and distribute executables, as long as you follow the
    requirements of the GNU GPL in regard to all of the software in the
    executable aside from CGAL.

\*---------------------------------------------------------------------------*/

#include "indexedCell.H"

// * * * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * //

template<class Gt, class Cb>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const InfoProxy<CGAL::indexedCell<Gt, Cb> >& p
)
{
    const CGAL::indexedCell<Gt, Cb>& iv = p.t_;

    os  << "Cell : index:" << iv.index_ << " filterCount:" << iv.filterCount_
        << nl;
    os  << "    " << iv.vertex(0)->info();
    os  << "    " << iv.vertex(1)->info();
    os  << "    " << iv.vertex(2)->info();
    os  << "    " << iv.vertex(3)->info();

    return os;
}


// ************************************************************************* //
