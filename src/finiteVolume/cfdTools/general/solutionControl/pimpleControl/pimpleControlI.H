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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::label Foam::pimpleControl::nCorrPIMPLE() const
{
    return nCorrPIMPLE_;
}


inline Foam::label Foam::pimpleControl::nCorrPISO() const
{
    return nCorrPISO_;
}


inline Foam::label Foam::pimpleControl::nCorrNonOrtho() const
{
    return nCorrPISO_;
}


inline Foam::label Foam::pimpleControl::corrPIMPLE() const
{
    return corrPIMPLE_;
}


inline Foam::label Foam::pimpleControl::corrPISO() const
{
    return corrPISO_;
}


inline Foam::label Foam::pimpleControl::corrNonOrtho() const
{
    return corrNonOrtho_;
}


inline bool Foam::pimpleControl::loop()
{
    corrPIMPLE_++;

    if (debug)
    {
        Info<< algorithmName_ << " loop: corrPIMPLE = " << corrPIMPLE_ << endl;
    }

    if (corrPIMPLE_ == nCorrPIMPLE_ + 1)
    {
        if ((!residualControl_.empty()) && (nCorrPIMPLE_ != 1))
        {
            Info<< algorithmName_ << ": not converged within "
                << nCorrPIMPLE_ << " iterations" << endl;
        }

        corrPIMPLE_ = 0;
        mesh_.data::remove("finalIteration");
        return false;
    }

    bool completed = false;
    if (criteriaSatisfied())
    {
        Info<< algorithmName_ << ": converged in " << corrPIMPLE_
            << " iterations" << endl;
        completed = true;
    }
    else
    {
        if (finalIter())
        {
            mesh_.data::add("finalIteration", true);
        }

        if (corrPIMPLE_ <= nCorrPIMPLE_)
        {
            if (nCorrPIMPLE_ != 1)
            {
                Info<< algorithmName_ << ": iteration " << corrPIMPLE_ << endl;
                storePrevIterFields();
            }

            completed = false;
        }
    }

    return !completed;
}


inline bool Foam::pimpleControl::correct()
{
    corrPISO_++;

    if (debug)
    {
        Info<< algorithmName_ << " correct: corrPISO = " << corrPISO_ << endl;
    }

    if (corrPISO_ <= nCorrPISO_)
    {
        return true;
    }
    else
    {
        corrPISO_ = 0;
        return false;
    }
}


inline bool Foam::pimpleControl::correctNonOrthogonal()
{
    corrNonOrtho_++;

    if (debug)
    {
        Info<< algorithmName_ << " correctNonOrthogonal: corrNonOrtho = "
            << corrNonOrtho_ << endl;
    }

    if (corrNonOrtho_ <= nCorrNonOrtho_ + 1)
    {
        return true;
    }
    else
    {
        corrNonOrtho_ = 0;
        return false;
    }
}


inline bool Foam::pimpleControl::storeInitialResiduals() const
{
    // start from second PIMPLE iteration
    return (corrPIMPLE_ == 2) && (corrPISO_ == 0) && (corrNonOrtho_ == 0);
}


inline bool Foam::pimpleControl::finalIter() const
{
    return corrPIMPLE_ == nCorrPIMPLE_;
}


inline bool Foam::pimpleControl::finalNonOrthogonalIter() const
{
    return corrNonOrtho_ == nCorrNonOrtho_ + 1;
}


inline bool Foam::pimpleControl::finalInnerIter() const
{
     return
        corrPIMPLE_ == nCorrPIMPLE_
     && corrPISO_ == nCorrPISO_
     && corrNonOrtho_ == nCorrNonOrtho_ + 1;
}


inline bool Foam::pimpleControl::turbCorr() const
{
    return !turbOnFinalIterOnly_ || finalIter();
}


// ************************************************************************* //
