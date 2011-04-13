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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::pimpleControl::start()
{
    corr_ = 0;

    return true;
}


inline bool Foam::pimpleControl::loop()
{
    if (criteriaSatisfied())
    {
        Info<< "PIMPLE loop converged in " << corr_ << " iterations" << endl;
        return false;
    }
    else
    {
        if (finalIter())
        {
            mesh_.data::add("finalIteration", true);
        }

        if (corr_ < nCorr_)
        {
            Info<< "PIMPLE iteration " << corr_ + 1 << endl;
            return true;
        }
        else
        {
            Info<< "PIMPLE loop not converged within " << nCorr_
                << " iterations" << endl;
            return false;
        }
    }
}


inline bool Foam::pimpleControl::finalIter() const
{
    return corr_ == nCorr_-1;
}


inline void Foam::pimpleControl::operator++(int)
{
    if (finalIter())
    {
        mesh_.data::remove("finalIteration");
    }

    corr_++;
}


// ************************************************************************* //
