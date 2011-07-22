/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2011 OpenCFD Ltd.
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

#include "indexedVertex.H"
//#include "conformalVoronoiMesh.H"
#include "point.H"

// * * * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * //

template<class Gt, class Vb>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const InfoProxy<CGAL::indexedVertex<Gt, Vb> >& p
)
{
    const CGAL::indexedVertex<Gt, Vb>& iv = p.t_;

    if (iv.type_ == CGAL::indexedVertex<Gt, Vb>::vtNearBoundary)
    {
        os  << "internal near boundary point" << nl;
    }
    else if (iv.type_ == CGAL::indexedVertex<Gt, Vb>::vtInternal)
    {
        os  << "internal point" << nl;
    }
    else if (iv.type_ == CGAL::indexedVertex<Gt, Vb>::vtFar)
    {
        os  << "far point" << nl;
    }
    else if (iv.type_ > CGAL::indexedVertex<Gt, Vb>::vtFar && iv.type_ < 0)
    {
        if (iv.index_ >= 0)
        {
            os  << "referred (master) point from processor " << iv.procIndex()
                << nl;
        }
        else
        {
            os  << "referred (slave) point from processor " << iv.procIndex()
                << nl;
        }
    }
    else if (iv.type_ >= 0)
    {
        if (iv.ppMaster())
        {
            os  << "master of point pair. paired up with point " << iv.index_
            << nl;
        }
        else if (iv.ppSlave())
        {
            os  << "slave of point pair. paired up with point " << iv.index_
            << nl;
        }
        else
        {
            FatalErrorIn
            (
                "operator<<"
                "(Ostream&, const InfoProxy<CGAL::indexedVertex<Gt, Vb> >&)"
            )   << "unhandled type " << iv.type_ << " index " << iv.index_
                << abort(FatalError);
        }
    }
    else
    {
        FatalErrorIn
        (
            "operator<<"
            "(Ostream&, const InfoProxy<CGAL::indexedVertex<Gt, Vb> >&)"
        )   << "unhandled type " << iv.type_ << " index " << iv.index_
            << abort(FatalError);
    }
    const Foam::point pos
    (
        CGAL::to_double(iv.point().x()),
        CGAL::to_double(iv.point().y()),
        CGAL::to_double(iv.point().z())
    );

    os  << "    type:" << iv.type_ << " index:" << iv.index_
        << " at:" << pos << nl;

    return os;
}


// ************************************************************************* //
