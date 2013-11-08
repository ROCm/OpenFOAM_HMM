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
        seulex::EPS = VSMALL,
        seulex::STEPFAC1 = 0.6,
        seulex::STEPFAC2 = 0.93,
        seulex::STEPFAC3 = 0.1,
        seulex::STEPFAC4 = 4.0,
        seulex::STEPFAC5 = 0.5,
        seulex::KFAC1 = 0.7,
        seulex::KFAC2 = 0.9;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::seulex::seulex(const ODESystem& ode, const dictionary& dict)
:
    ODESolver(ode, dict),
    ode_(ode),
    nSeq_(IMAXX),
    cost_(IMAXX),
    dfdx_(n_),
    factrl_(IMAXX),
    table_(KMAXX,n_),
    fSave_((IMAXX-1)*(IMAXX+1)/2 + 2,n_),
    dfdy_(n_),
    calcJac_(false),
    dense_(false),
    a_(n_),
    coeff_(IMAXX,IMAXX),
    pivotIndices_(n_, 0.0),
    hOpt_(IMAXX),
    work_(IMAXX),
    ySaved_(n_),
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

    for (label i = 2; i < IMAXX; i++)
    {
        nSeq_[i] = 2*nSeq_[i-2];
    }
    cost_[0] = costjac + costlu + nSeq_[0]*(costfunc + costsolve);

    for (label k=0;k<KMAXX;k++)
    {
        cost_[k+1] = cost_[k] + (nSeq_[k+1]-1)*(costfunc + costsolve) + costlu;
    }

    hnext_ =- VGREAT;

    //NOTE: the first element of relTol_ and absTol_ are used here.
    scalar logfact =- log10(relTol_[0] + absTol_[0])*0.6 + 0.5;

    kTarg_ = min(1, min(KMAXX - 1, label(logfact)));

    for (label k = 0; k < IMAXX; k++)
    {
        for (label l = 0; l < k; l++)
        {
            scalar ratio = scalar(nSeq_[k])/nSeq_[l];
            coeff_[k][l] = 1.0/(ratio - 1.0);
        }
    }

    factrl_[0] = 1.0;
    for (label k = 0; k < IMAXX - 1; k++)
    {
        factrl_[k+1] = (k+1)*factrl_[k];
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


void Foam::seulex::solve
(
    const ODESystem& ode,
    scalar& x,
    scalarField& y,
    scalar& dxTry
) const
{
    step(dxTry, x, y);
    dxTry = hnext_;
}


void Foam::seulex::step(const scalar& htry, scalar& x,  scalarField& y) const
{
    static bool first_step = true, last_step = false;
    static bool forward, reject = false, prev_reject = false;
    static scalar errold;
    label i, k;
    scalar fac, h, hnew, err;
    bool firstk;

    work_[0] = 1.e30;
    h = htry;
    forward = h > 0 ? true : false;

    ySaved_ = y;

    if (h != hnext_ && !first_step)
    {
        last_step = true;
    }

    if (reject)
    {
        prev_reject = true;
        last_step = false;
        theta_=2.0*jacRedo_;
    }

    for (i = 0; i < n_; i++)
    {
        scale_[i] = absTol_[i] + relTol_[i]*mag(y[i]);
    }

    reject = false;
    firstk = true;
    hnew = mag(h);

    if (theta_ > jacRedo_ && !calcJac_)
    {
        ode_.jacobian(x, y, dfdx_, dfdy_);
        calcJac_ = true;
    }

    while (firstk || reject)
    {
        h = forward ? hnew : -hnew;
        firstk = false;
        reject = false;

        if (mag(h) <= mag(x)*EPS)
        {
             WarningIn("seulex::step(const scalar")
                    << "step size underflow in step :"  << h << endl;
        }
        label ipt=-1;

        for (k = 0; k <= kTarg_ + 1; k++)
        {
            bool success = dy(x, ySaved_, h, k, ySequence_, ipt, scale_);

            if (!success)
            {
                reject = true;
                hnew = mag(h)*STEPFAC5;
                break;
            }

            if (k == 0)
            {
                 y = ySequence_;
            }
            else
            {
                for (i = 0; i < n_; i++)
                    table_[k-1][i] = ySequence_[i];
            }
            if (k != 0)
            {
                polyextr(k, table_, y);
                err = 0.0;
                for (i = 0; i < n_; i++)
                {
                    scale_[i] = absTol_[i] + relTol_[i]*mag(ySaved_[i]);
                    err += sqr((y[i] - table_[0][i])/scale_[i]);
                }
                err = sqrt(err/n_);
                if (err > 1.0/EPS || (k > 1 && err >= errold))
                {
                    reject = true;
                    hnew = mag(h)*STEPFAC5;
                    break;
                }
                errold = min(4.0*err, 1.0);
                scalar expo = 1.0/(k + 1);
                scalar facmin = pow(STEPFAC3, expo);
                if (err == 0.0)
                {
                    fac = 1.0/facmin;
                }
                else
                {
                    fac = STEPFAC2/pow(err/STEPFAC1, expo);
                    fac = max(facmin/STEPFAC4, min(1.0/facmin, fac));
                }
                hOpt_[k] = mag(h*fac);
                work_[k] = cost_[k]/hOpt_[k];

                if ((first_step || last_step) && err <= 1.0)
                {
                    break;
                }
                if (k == kTarg_-1 && !prev_reject && !first_step && !last_step)
                {
                    if (err <= 1.0)
                    {
                        break;
                    }
                    else if (err > nSeq_[kTarg_]*nSeq_[kTarg_ + 1]*4.0)
                    {
                        reject = true;
                        kTarg_ = k;
                        if (kTarg_>1 && work_[k-1] < KFAC1*work_[k])
                        {
                            kTarg_--;
                        }
                        hnew = hOpt_[kTarg_];
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
                        if (kTarg_>1 && work_[k-1] < KFAC1*work_[k])
                        {
                            kTarg_--;
                        }
                        hnew = hOpt_[kTarg_];
                        break;
                    }
                }
                if (k == kTarg_+1)
                {
                    if (err > 1.0)
                    {
                        reject = true;
                        if (kTarg_ > 1 && work_[kTarg_-1] < KFAC1*work_[kTarg_])
                        {
                            kTarg_--;
                        }
                        hnew = hOpt_[kTarg_];
                    }
                    break;
                }
            }
        }
        if (reject)
        {
            prev_reject = true;
            if (!calcJac_)
            {
                theta_ = 2.0*jacRedo_;

                if (theta_ > jacRedo_ && !calcJac_)
                {
                    ode_.jacobian(x, y, dfdx_, dfdy_);
                    calcJac_ = true;
                }
            }
        }
    }

    calcJac_ = false;

    x += h;
    hdid_ = h;
    first_step = false;

    label kopt;
    if (k == 1)
    {
        kopt = 2;
    }
    else if (k <= kTarg_)
    {
        kopt=k;
        if (work_[k-1] < KFAC1*work_[k])
        {
            kopt = k - 1;
        }
        else if (work_[k] < KFAC2*work_[k - 1])
        {
            kopt = min(k + 1, KMAXX - 1);
        }
    }
    else
    {
        kopt = k - 1;
        if (k > 2 && work_[k-2] < KFAC1*work_[k - 1])
        {
            kopt = k - 2;
        }
        if (work_[k] < KFAC2*work_[kopt])
        {
            kopt = min(k, KMAXX - 1);
        }
    }

    if (prev_reject)
    {
        kTarg_ = min(kopt, k);
        hnew = min(mag(h), hOpt_[kTarg_]);
        prev_reject = false;
    }
    else
    {
        if (kopt <= k)
        {
            hnew = hOpt_[kopt];
        }
        else
        {
            if (k < kTarg_ && work_[k] < KFAC2*work_[k - 1])
            {
                hnew = hOpt_[k]*cost_[kopt + 1]/cost_[k];
            }
            else
            {
                hnew = hOpt_[k]*cost_[kopt]/cost_[k];
            }
        }
        kTarg_ = kopt;
    }

    if (forward)
    {
        hnext_ = hnew;
    }
    else
    {
        hnext_ =- hnew;
    }
}


bool Foam::seulex::dy
(
    const scalar& x,
    scalarField& y,
    const scalar htot,
    const label k,
    scalarField& yend,
    label& ipt,
    scalarField& scale
) const
{
    label nstep = nSeq_[k];
    scalar h = htot/nstep;
    for (label i = 0; i < n_; i++)
    {
        for (label j = 0; j < n_; j++) a_[i][j] = -dfdy_[i][j];
        a_[i][i] += 1.0/h;
    }
    LUDecompose(a_, pivotIndices_);
    scalar xnew = x + h;

    ode_.derivatives(xnew, y, del_);

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
        xnew += h;

        ode_.derivatives(xnew, yTemp_, yend);

        if (nn == 1 && k<=1)
        {
            scalar del1=0.0;
            for (label i = 0; i < n_; i++)
            {
                del1 += sqr(del_[i]/scale[i]);
            }
            del1 = sqrt(del1);

            ode_.derivatives(x+h, yTemp_, dyTemp_);
            for (label i=0;i<n_;i++)
            {
                del_[i] = dyTemp_[i] - del_[i]/h;
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
    label l=last.size();
    for (label j = k - 1; j > 0; j--)
    {
        for (label i=0; i<l; i++)
        {
            table[j-1][i] =
                table[j][i] + coeff_[k][j]*(table[j][i] - table[j-1][i]);
        }
    }
    for (label i = 0; i < l; i++)
    {
        last[i] = table[0][i] + coeff_[k][0]*(table[0][i] - last[i]);
    }
}

// ************************************************************************* //
