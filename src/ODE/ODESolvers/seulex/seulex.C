/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "seulex.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(seulex, 0);

    addToRunTimeSelectionTable(ODESolver, seulex, dictionary);

    const scalar
        seulex::stepFactor1_ = 0.6,
        seulex::stepFactor2_ = 0.93,
        seulex::stepFactor3_ = 0.1,
        seulex::stepFactor4_ = 4.0,
        seulex::stepFactor5_ = 0.5,
        seulex::kFactor1_ = 0.7,
        seulex::kFactor2_ = 0.9;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::seulex::seulex(const ODESystem& ode, const dictionary& dict)
:
    ODESolver(ode, dict),
    nSeq_(iMaxx_),
    cost_(iMaxx_),
    dfdx_(n_),
    factrl_(iMaxx_),
    table_(kMaxx_, n_),
    fSave_((iMaxx_ - 1)*(iMaxx_ + 1)/2 + 2, n_),
    dfdy_(n_),
    calcJac_(false),
    dense_(false),
    a_(n_),
    coeff_(iMaxx_, iMaxx_),
    pivotIndices_(n_, 0.0),
    dxOpt_(iMaxx_),
    work_(iMaxx_),
    y0_(n_),
    ySequence_(n_),
    scale_(n_),
    del_(n_),
    yTemp_(n_),
    dyTemp_(n_),
    delTemp_(n_)
{
    const scalar costfunc = 1.0, costjac = 5.0, costlu = 1.0, costsolve = 1.0;

    jacRedo_ = min(1.0e-4, min(relTol_));
    theta_ = 2.0*jacRedo_;
    nSeq_[0] = 2;
    nSeq_[1] = 3;

    for (int i=2; i<iMaxx_; i++)
    {
        nSeq_[i] = 2*nSeq_[i-2];
    }
    cost_[0] = costjac + costlu + nSeq_[0]*(costfunc + costsolve);

    for (int k=0; k<kMaxx_; k++)
    {
        cost_[k+1] = cost_[k] + (nSeq_[k+1]-1)*(costfunc + costsolve) + costlu;
    }

    //NOTE: the first element of relTol_ and absTol_ are used here.
    scalar logfact =- log10(relTol_[0] + absTol_[0])*0.6 + 0.5;

    kTarg_ = min(1, min(kMaxx_ - 1, label(logfact)));

    for (int k=0; k<iMaxx_; k++)
    {
        for (int l=0; l<k; l++)
        {
            scalar ratio = scalar(nSeq_[k])/nSeq_[l];
            coeff_[k][l] = 1.0/(ratio - 1.0);
        }
    }

    factrl_[0] = 1.0;
    for (int k=0; k<iMaxx_ - 1; k++)
    {
        factrl_[k+1] = (k+1)*factrl_[k];
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool Foam::seulex::dy
(
    const scalar x,
    scalarField& y,
    const scalar dxTot,
    const label k,
    scalarField& yend,
    label& ipt,
    scalarField& scale
) const
{
    label nstep = nSeq_[k];
    scalar dx = dxTot/nstep;

    for (label i=0; i<n_; i++)
    {
        for (label j=0; j<n_; j++)
        {
            a_[i][j] = -dfdy_[i][j];
        }

        a_[i][i] += 1.0/dx;
    }

    LUDecompose(a_, pivotIndices_);

    scalar xnew = x + dx;

    odes_.derivatives(xnew, y, del_);

    yTemp_ = y;

    LUBacksubstitute(a_, pivotIndices_, del_);

    if (dense_ && nstep == k + 1)
    {
        ipt++;
        for (label i = 0; i < n_; i++)
        {
            fSave_[ipt][i] = del_[i];
        }
    }
    for (label nn = 1; nn < nstep; nn++)
    {
        for (label i=0;i<n_;i++)
        {
            yTemp_[i] += del_[i];
        }
        xnew += dx;

        odes_.derivatives(xnew, yTemp_, yend);

        if (nn == 1 && k<=1)
        {
            scalar del1=0.0;
            for (label i = 0; i < n_; i++)
            {
                del1 += sqr(del_[i]/scale[i]);
            }
            del1 = sqrt(del1);

            odes_.derivatives(x + dx, yTemp_, dyTemp_);
            for (label i=0;i<n_;i++)
            {
                del_[i] = dyTemp_[i] - del_[i]/dx;
            }

            LUBacksubstitute(a_, pivotIndices_, del_);

            scalar del2 = 0.0;
            for (label i = 0; i <n_ ; i++)
            {
                del2 += sqr(del_[i]/scale[i]);
            }
            del2 = sqrt(del2);
            theta_ = del2 / min(1.0, del1 + SMALL);

            if (theta_ > 1.0)
            {
                return false;
            }
        }

        delTemp_ = yend;
        LUBacksubstitute(a_, pivotIndices_, delTemp_);
        del_ = delTemp_;

        if (dense_ && nn >= nstep-k-1)
        {
            ipt++;
            for (label i=0;i<n_;i++)
            {
                fSave_[ipt][i]=del_[i];
            }
        }
    }

    yend = yTemp_ + del_;

    return true;
}


void Foam::seulex::polyextr
(
    const label k,
    scalarRectangularMatrix& table,
    scalarField& last
) const
{
    int l = last.size();
    for (int j=k - 1; j>0; j--)
    {
        for (label i=0; i<l; i++)
        {
            table[j-1][i] =
                table[j][i] + coeff_[k][j]*(table[j][i] - table[j-1][i]);
        }
    }

    for (int i=0; i<l; i++)
    {
        last[i] = table[0][i] + coeff_[k][0]*(table[0][i] - last[i]);
    }
}


void Foam::seulex::solve
(
    scalar& x,
    scalarField& y,
    scalar& dxTry
) const
{
    static bool firstStep = true, lastStep = false;
    static bool reject = false, prevReject = false;
    static scalar errold;

    work_[0] = GREAT;
    scalar dx = dxTry;
    dxTry = -VGREAT;

    bool forward = dx > 0 ? true : false;

    y0_ = y;

    if (dx != dxTry && !firstStep)
    {
        lastStep = true;
    }

    if (reject)
    {
        prevReject = true;
        lastStep = false;
        theta_=2.0*jacRedo_;
    }

    forAll(scale_, i)
    {
        scale_[i] = absTol_[i] + relTol_[i]*mag(y[i]);
    }

    reject = false;
    bool firstk = true;
    scalar dxNew = mag(dx);

    if (theta_ > jacRedo_ && !calcJac_)
    {
        odes_.jacobian(x, y, dfdx_, dfdy_);
        calcJac_ = true;
    }


    int k;
    while (firstk || reject)
    {
        dx = forward ? dxNew : -dxNew;
        firstk = false;
        reject = false;

        if (mag(dx) <= mag(x)*SMALL)
        {
             WarningIn("seulex::step(const scalar")
                    << "step size underflow in step :"  << dx << endl;
        }
        label ipt = -1;

        for (k = 0; k <= kTarg_ + 1; k++)
        {
            bool success = dy(x, y0_, dx, k, ySequence_, ipt, scale_);

            if (!success)
            {
                reject = true;
                dxNew = mag(dx)*stepFactor5_;
                break;
            }

            if (k == 0)
            {
                 y = ySequence_;
            }
            else
            {
                forAll(ySequence_, i)
                {
                    table_[k-1][i] = ySequence_[i];
                }
            }
            if (k != 0)
            {
                polyextr(k, table_, y);
                scalar err = 0.0;
                forAll(scale_, i)
                {
                    scale_[i] = absTol_[i] + relTol_[i]*mag(y0_[i]);
                    err += sqr((y[i] - table_[0][i])/scale_[i]);
                }
                err = sqrt(err/n_);
                if (err > 1.0/SMALL || (k > 1 && err >= errold))
                {
                    reject = true;
                    dxNew = mag(dx)*stepFactor5_;
                    break;
                }
                errold = min(4.0*err, 1.0);
                scalar expo = 1.0/(k + 1);
                scalar facmin = pow(stepFactor3_, expo);
                scalar fac;
                if (err == 0.0)
                {
                    fac = 1.0/facmin;
                }
                else
                {
                    fac = stepFactor2_/pow(err/stepFactor1_, expo);
                    fac = max(facmin/stepFactor4_, min(1.0/facmin, fac));
                }
                dxOpt_[k] = mag(dx*fac);
                work_[k] = cost_[k]/dxOpt_[k];

                if ((firstStep || lastStep) && err <= 1.0)
                {
                    break;
                }
                if (k == kTarg_-1 && !prevReject && !firstStep && !lastStep)
                {
                    if (err <= 1.0)
                    {
                        break;
                    }
                    else if (err > nSeq_[kTarg_]*nSeq_[kTarg_ + 1]*4.0)
                    {
                        reject = true;
                        kTarg_ = k;
                        if (kTarg_>1 && work_[k-1] < kFactor1_*work_[k])
                        {
                            kTarg_--;
                        }
                        dxNew = dxOpt_[kTarg_];
                        break;
                    }
                }
                if (k == kTarg_)
                {
                    if (err <= 1.0)
                    {
                        break;
                    }
                    else if (err > nSeq_[k + 1]*2.0)
                    {
                        reject = true;
                        if (kTarg_>1 && work_[k-1] < kFactor1_*work_[k])
                        {
                            kTarg_--;
                        }
                        dxNew = dxOpt_[kTarg_];
                        break;
                    }
                }
                if (k == kTarg_+1)
                {
                    if (err > 1.0)
                    {
                        reject = true;
                        if (kTarg_ > 1 && work_[kTarg_-1] < kFactor1_*work_[kTarg_])
                        {
                            kTarg_--;
                        }
                        dxNew = dxOpt_[kTarg_];
                    }
                    break;
                }
            }
        }
        if (reject)
        {
            prevReject = true;
            if (!calcJac_)
            {
                theta_ = 2.0*jacRedo_;

                if (theta_ > jacRedo_ && !calcJac_)
                {
                    odes_.jacobian(x, y, dfdx_, dfdy_);
                    calcJac_ = true;
                }
            }
        }
    }

    calcJac_ = false;

    x += dx;
    firstStep = false;

    label kopt;
    if (k == 1)
    {
        kopt = 2;
    }
    else if (k <= kTarg_)
    {
        kopt=k;
        if (work_[k-1] < kFactor1_*work_[k])
        {
            kopt = k - 1;
        }
        else if (work_[k] < kFactor2_*work_[k - 1])
        {
            kopt = min(k + 1, kMaxx_ - 1);
        }
    }
    else
    {
        kopt = k - 1;
        if (k > 2 && work_[k-2] < kFactor1_*work_[k - 1])
        {
            kopt = k - 2;
        }
        if (work_[k] < kFactor2_*work_[kopt])
        {
            kopt = min(k, kMaxx_ - 1);
        }
    }

    if (prevReject)
    {
        kTarg_ = min(kopt, k);
        dxNew = min(mag(dx), dxOpt_[kTarg_]);
        prevReject = false;
    }
    else
    {
        if (kopt <= k)
        {
            dxNew = dxOpt_[kopt];
        }
        else
        {
            if (k < kTarg_ && work_[k] < kFactor2_*work_[k - 1])
            {
                dxNew = dxOpt_[k]*cost_[kopt + 1]/cost_[k];
            }
            else
            {
                dxNew = dxOpt_[k]*cost_[kopt]/cost_[k];
            }
        }
        kTarg_ = kopt;
    }

    if (forward)
    {
        dxTry = dxNew;
    }
    else
    {
        dxTry = -dxNew;
    }
}


// ************************************************************************* //
