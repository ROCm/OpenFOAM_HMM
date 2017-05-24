/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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

#include "triSurfaceLoader.H"
#include "fileNameList.H"
#include "Time.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::triSurfaceLoader::triSurfaceLoader(const fileName& directory)
:
    directory_(directory),
    available_(),
    selected_()
{
    readDir();
}


Foam::triSurfaceLoader::triSurfaceLoader(const Time& runTime)
:
    directory_(runTime.constantPath()/"triSurface"),
    available_(),
    selected_()
{
    readDir();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::triSurfaceLoader::~triSurfaceLoader()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::triSurfaceLoader::readDir()
{
    fileNameList files = Foam::readDir(directory_, fileName::FILE);

    // Will be using triSurface
    //
    // - filter according to what is supported
    //
    // Transform from fileName to word and eliminate duplicates
    // (eg, files with/without .gz)
    wordHashSet names(2*files.size());

    forAll(files, filei)
    {
        const fileName& f = files[filei];
        if (triSurface::canRead(f))
        {
            names.insert(f.name());
        }
    }
    available_ = names.sortedToc();     // Also hashes the names

    return available_.size();
}


Foam::label Foam::triSurfaceLoader::selectAll()
{
    selected_ = available_;
    return selected_.size();
}


Foam::label Foam::triSurfaceLoader::select(const word& name)
{
    if (available_.found(name))
    {
        selected_ = wordList{name};  // hashedWordList::operator[] is hidden!
    }
    else
    {
        selected_.clear();
    }

    return selected_.size();
}


Foam::label Foam::triSurfaceLoader::select(const wordRe& mat)
{
    DynamicList<label> foundIds(available_.size());

    if (mat.isPattern())
    {
        foundIds = findMatchingStrings(mat, available_);
        sort(foundIds);
    }
    else
    {
        const word& plain = static_cast<const word&>(mat);
        if (available_.found(plain))
        {
            foundIds.append(available_[plain]);
        }
        else
        {
            FatalErrorInFunction
                << "Specified the surfaces " << mat << nl
                << "  - but could not find it"
                << exit(FatalError);
        }
    }

    selected_ = wordList(available_, foundIds);
    return selected_.size();
}


Foam::label Foam::triSurfaceLoader::select(const wordReList& matcher)
{
    // Need to be more careful when select.
    // - preserve same order as the input matcher itself
    // - avoid duplicate selections
    // - flag explicitly named files as missing.
    //   (Eg, foo.stl must exist, but "foo.*\.stl" is optional)

    // Track which files have already been selected
    DynamicList<label> foundIds(available_.size());
    labelHashSet hashedFound(2*available_.size());

    DynamicList<word>  missing(matcher.size());
    wordHashSet hashedMissing(2*matcher.size());

    // Exact matches must exist
    forAll(matcher, i)
    {
        const wordRe& mat = matcher[i];

        if (mat.isPattern())
        {
            labelList indices = findMatchingStrings(mat, available_);
            sort(indices);

            forAll(indices, j)
            {
                const label idx = indices[j];
                if (hashedFound.insert(idx))
                {
                    foundIds.append(idx);
                }
            }
        }
        else
        {
            const word& plain = static_cast<const word&>(mat);
            if (available_.found(plain))
            {
                const label idx = available_[plain];
                if (hashedFound.insert(idx))
                {
                    foundIds.append(idx);
                }
            }
            else if (hashedMissing.insert(plain))
            {
                missing.append(plain);
            }
        }
    }

    if (missing.size())
    {
        FatalErrorInFunction
            << "Specified the surfaces " << flatOutput(matcher) << nl
            << "  - but could not find " << flatOutput(missing)
            << exit(FatalError);
    }

    selected_ = wordList(available_, foundIds);
    return selected_.size();
}


Foam::autoPtr<Foam::triSurface> Foam::triSurfaceLoader::load() const
{
    if (selected_.empty())
    {
        return autoPtr<triSurface>();
    }
    else if (selected_.size() == 1)
    {
        return autoPtr<triSurface>(new triSurface(directory_/selected_[0]));
    }

    List<labelledTri> faces;
    pointField points;
    label regoff = 0; // region offset

    // collect all patches, preserving names.
    // This will be horrible for output, but is good if we are relying on
    // the names for defining baffles.
    DynamicList<geometricSurfacePatch> patches(16*selected_.size());

    forAll(selected_, surfi)
    {
        triSurface addsurf(directory_/selected_[surfi]);

        List<labelledTri> addfaces(addsurf.xferFaces());
        List<point> addpoints(addsurf.xferPoints());
        patches.append(addsurf.patches());

        if (surfi)
        {
            const label ptoff = points.size(); // point offset

            forAll(addfaces, facei)
            {
                labelledTri& f = addfaces[facei];
                forAll(f, fi)
                {
                    f[fi] += ptoff;
                }
                f.region() += regoff;
            }

            faces.append(addfaces);
            points.append(addpoints);
        }
        else
        {
            faces.transfer(addfaces);
            points.transfer(addpoints);
        }

        regoff += addsurf.patches().size();
    }

    return autoPtr<triSurface>(new triSurface(faces, patches, points, true));
}


// ************************************************************************* //
