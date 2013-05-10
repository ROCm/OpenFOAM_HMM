/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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

#include "DelaunayMesh.H"
#include "labelPair.H"
#include "PrintTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Triangulation>
Foam::DelaunayMesh<Triangulation>::DelaunayMesh()
:
    Triangulation(),
    vertexCount_(0),
    cellCount_(0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Triangulation>
Foam::DelaunayMesh<Triangulation>::~DelaunayMesh()
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Triangulation>
void Foam::DelaunayMesh<Triangulation>::reset()
{
    Info<< "Clearing triangulation" << endl;

    this->clear();

    resetVertexCount();
    resetCellCount();
}


template<class Triangulation>
void Foam::DelaunayMesh<Triangulation>::insertPoints(const List<Vb>& vertices)
{
    rangeInsertWithInfo
    (
        vertices.begin(),
        vertices.end(),
        true
    );
}


template<class Triangulation>
bool Foam::DelaunayMesh<Triangulation>::Traits_for_spatial_sort::Less_x_3::
operator()
(
    const Point_3& p,
    const Point_3& q
) const
{
    return typename Gt::Less_x_3()(*(p.first), *(q.first));
}

template<class Triangulation>
bool Foam::DelaunayMesh<Triangulation>::Traits_for_spatial_sort::Less_y_3::
operator()
(
    const Point_3& p,
    const Point_3& q
) const
{
    return typename Gt::Less_y_3()(*(p.first), *(q.first));
}

template<class Triangulation>
bool Foam::DelaunayMesh<Triangulation>::Traits_for_spatial_sort::Less_z_3::
operator()
(
    const Point_3& p,
    const Point_3& q
) const
{
    return typename Gt::Less_z_3()(*(p.first), *(q.first));
}

template<class Triangulation>
typename Foam::DelaunayMesh<Triangulation>::Traits_for_spatial_sort::Less_x_3
Foam::DelaunayMesh<Triangulation>::Traits_for_spatial_sort::less_x_3_object()
const
{
    return Less_x_3();
}

template<class Triangulation>
typename Foam::DelaunayMesh<Triangulation>::Traits_for_spatial_sort::Less_y_3
Foam::DelaunayMesh<Triangulation>::Traits_for_spatial_sort::less_y_3_object()
const
{
    return Less_y_3();
}

template<class Triangulation>
typename Foam::DelaunayMesh<Triangulation>::Traits_for_spatial_sort::Less_z_3
Foam::DelaunayMesh<Triangulation>::Traits_for_spatial_sort::less_z_3_object()
const
{
    return Less_z_3();
}


template<class Triangulation>
template<class PointIterator>
void Foam::DelaunayMesh<Triangulation>::rangeInsertWithInfo
(
    PointIterator begin,
    PointIterator end,
    bool printErrors
)
{
    typedef DynamicList
    <
        std::pair
        <
            const typename Triangulation::Point*,
            label
        >
    > vectorPairPointIndex;

    vectorPairPointIndex points;

    label count = 0;
    for (PointIterator it = begin; it != end; ++it)
    {
        points.append
        (
            std::make_pair(&(it->point()), count++)
        );
    }

    std::random_shuffle(points.begin(), points.end());

    spatial_sort
    (
        points.begin(),
        points.end(),
        Traits_for_spatial_sort()
    );

    Vertex_handle hint;

    for
    (
        typename vectorPairPointIndex::const_iterator p = points.begin();
        p != points.end();
        ++p
    )
    {
        const size_t checkInsertion = Triangulation::number_of_vertices();

        hint = this->insert(*(p->first), hint);

        const Vb& vert = *(begin + p->second);

        if (checkInsertion != Triangulation::number_of_vertices() - 1)
        {
            if (printErrors)
            {
                Vertex_handle nearV =
                    Triangulation::nearest_vertex(*(p->first));

                Pout<< "Failed insertion : " << vert.info()
                    << "         nearest : " << nearV->info();
            }
        }
        else
        {
            hint->index() = getNewVertexIndex();
            hint->type() = vert.type();
            hint->procIndex() = vert.procIndex();
            hint->targetCellSize() = vert.targetCellSize();
            hint->alignment() = vert.alignment();
        }
    }
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "DelaunayMeshIO.C"

// ************************************************************************* //
