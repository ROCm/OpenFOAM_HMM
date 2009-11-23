/*---------------------------------------------------------------------------*\
 =========                   |
 \\      /   F ield          | OpenFOAM: The Open Source CFD Toolbox
  \\    /    O peration      |
   \\  /     A nd            | Copyright (C) 2008-2009 OpenCFD Ltd.
    \\/      M anipulation   |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    DistributionTest

Description
    Test the Distribution class

    Plot normal distribution test in gnuplot using:

    @verbatim
    normalDistribution(mean, sigma, x) = \
        sqrt(1.0/(2.0*pi*sigma**2))*exp(-(x - mean)**2.0/(2.0*sigma**2))

    plot normalDistribution(8.5, 2.5, x)
    @endverbatim

\*---------------------------------------------------------------------------*/

#include "vector.H"
#include "labelVector.H"
#include "tensor.H"
#include "Distribution.H"
#include "Random.H"
#include "dimensionedTypes.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

int main(int argc, char *argv[])
{
    Random R(918273);

    {
        // scalar
        label randomDistributionTestSize = 50000000;

        Distribution<scalar> dS(scalar(5e-2));

        Info<< nl << "Distribution<scalar>" << nl
            << "Sampling "
            << randomDistributionTestSize
            << " times from GaussNormal distribution."
            << endl;

        for (label i = 0; i < randomDistributionTestSize; i++)
        {
            dS.add(2.5*R.GaussNormal() + 8.5);
        }

        Info<< "Mean " << dS.mean() << nl
            << "Median " << dS.median()
            << endl;

        dS.write("Distribution_scalar_test", dS.normalised());
    }

    {
        // vector
        Distribution<vector> dV(vector(0.1, 0.05, 0.15));

        label randomDistributionTestSize = 1000000;

        Info<< nl << "Distribution<vector>" << nl
            << "Sampling "
            << randomDistributionTestSize
            << " times from uniform and GaussNormal distribution."
            << endl;

        for (label i = 0; i < randomDistributionTestSize; i++)
        {
            dV.add(R.vector01());

            // Adding separate GaussNormal components with component
            // weights

            dV.add
            (
                vector
                (
                    R.GaussNormal()*3.0 + 1.5,
                    R.GaussNormal()*0.25 + 4.0,
                    R.GaussNormal()*3.0 - 1.5
                ),
                vector(1.0, 2.0, 5.0)
            );
        }

        Info<< "Mean " << dV.mean() << nl
            << "Median " << dV.median()
            << endl;

        dV.write("Distribution_vector_test", dV.normalised());
    }

    {
        // labelVector
        Distribution<labelVector> dLV(labelVector::one*10);

        label randomDistributionTestSize = 2000000;

        Info<< nl << "Distribution<labelVector>" << nl
            << "Sampling "
            << randomDistributionTestSize
            << " times from uniform distribution."
            << endl;

        for (label i = 0; i < randomDistributionTestSize; i++)
        {
            dLV.add
            (
                labelVector
                (
                    R.integer(-1000, 1000),
                    R.integer(-5000, 5000),
                    R.integer(-2000, 7000)
                )
            );
        }

        Info<< "Mean " << dLV.mean() << nl
            << "Median " << dLV.median()
            << endl;

        dLV.write("Distribution_labelVector_test", dLV.normalised());
    }

    {
        // tensor
        Distribution<tensor> dT(tensor::one*1e-2);

        label randomDistributionTestSize = 2000000;

        Info<< nl << "Distribution<tensor>" << nl
            << "Sampling "
            << randomDistributionTestSize
            << " times from uniform distribution."
            << endl;

        for (label i = 0; i < randomDistributionTestSize; i++)
        {
            dT.add(R.tensor01());
        }

        Info<< "Mean " << dT.mean() << nl
            << "Median " << dT.median()
            << endl;

        dT.write("Distribution_tensor_test", dT.normalised());
    }

    {
        // symmTensor
        Distribution<symmTensor> dSyT(symmTensor::one*1e-2);

        label randomDistributionTestSize = 2000000;

        Info<< nl << "Distribution<symmTensor>" << nl
            << "Sampling "
            << randomDistributionTestSize
            << " times from uniform distribution."
            << endl;

        for (label i = 0; i < randomDistributionTestSize; i++)
        {
            dSyT.add(R.symmTensor01());
        }

        Info<< "Mean " << dSyT.mean() << nl
            << "Median " << dSyT.median()
            << endl;

        dSyT.write("Distribution_symmTensor_test", dSyT.normalised());
    }

    {
        // sphericalTensor
        Distribution<sphericalTensor> dSpT(sphericalTensor::one*1e-2);

        label randomDistributionTestSize = 50000000;

        Info<< nl << "Distribution<sphericalTensor>" << nl
            << "Sampling "
            << randomDistributionTestSize
            << " times from uniform distribution."
            << endl;

        for (label i = 0; i < randomDistributionTestSize; i++)
        {
            dSpT.add(R.sphericalTensor01());
        }

        Info<< "Mean " << dSpT.mean() << nl
            << "Median " << dSpT.median()
            << endl;

        dSpT.write("Distribution_sphericalTensor_test", dSpT.normalised());
    }

    Info<< nl << "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
