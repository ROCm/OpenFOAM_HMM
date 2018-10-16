/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
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

Class
    Foam::functionObjects::stabilityBlendingFactor

Group
    grpFieldFunctionObjects

Description
    Calculates and outputs the stabilityBlendingFactor to be used by the
    local blended convection scheme. The output is a surface field weight
    between 0-1

    The weight of a blended scheme is given by a function of the blending
    factor, f:

    \f[
        weight = f scheme1 + (1 - f) scheme2
    \f]

    The factor is calculated based on six criteria:

    1) mesh non-orthogonality field
    2) magnitude of cell centres gradient
    3) convergence rate of residuals
    4) faceWeight
    5) skewness
    6) Co number

    The user can enable them individually.

    For option 1, the following relation is used:
        \f[
            fNon =
                min
                (
                    max
                    (
                        0.0,
                        (nonOrthogonality - maxNonOrthogonality)
                       /(minNonOrthogonality - maxNonOrthogonality)
                    ),
                    1.0
                )
        \f]

    For option 2, the following relation is used:
        \f[
            fMagGradCc =
                min
                (
                    max
                    (
                        0.0,
                        (magGradCC - maxGradCc)
                      / (minGradCc - maxGradCc)
                    ),
                    1.0
                )
        \f]
    Note that  magGradCC is equal to 3 for ortoghonal meshes

    For option 3 a PID control is used in order to control residual
    unbouded fluctuations for individual cells.

        \f[
            factor =
                P*residual
              + I*residualIntegral
              + D*residualDifferential
        \f]

        where P, I and D are user inputs

    The following relation is used:
    \f[
        fRes = (factor - meanRes)/(maxRes*meanRes);
    \f]

        where:
            meanRes = average(residual)
            maxRes is an user input


    fRes will blend more towards one as the cell residual is larger then
    the domain mean residuals.


    For option 4 , the following relation is used:

        \f[
            ffaceWeight = min
            (
                max
                (
                    0.0,
                    (minFaceWeight - faceWeights)
                  / (minFaceWeight - maxFaceWeight)
                ),
                1.0
            )
        \f]
    Note that faceWeights for a orthogonal mesh is 0.5.


    For option 5 , the following relation is used:

        \f[
            fskewness =
            min
            (
                max
                (
                    0.0,
                    (skewness    - maxSkewness)
                  / (minSkewness - maxSkewness)
                ),
                1.0
            )
        \f]


    For option 6 , the following relation is used:

    \f[
        fCoWeight = min(max((Co - Co1)/(Co2 - Co1), 0), 1)
    \f]
    where

    \vartable
        Co1 | Courant number below which scheme2 is used
        Co2 | Courant number above which scheme1 is used
    \endvartable

    The final factor is determined by:

    \f[
        f = max(fNon, fMagGradCc, fRes, ffaceWeight, fskewness, fCoWeight)
    \f]

    An indicator (volume) field, named blendedIndicator is generated if the log
    flag is on:
    - 1 represent scheme1 as active,
    - 0 represent scheme2 as active.

    Additional reporting is written to the standard output, providing
    statistics as to the number of cells used by each scheme.

Usage
    Example of function object specification to calculate the blending factor:
    \verbatim
    stabilityBlendingFactor1
    {
        type                stabilityBlendingFactor;
        libs                ("libfieldFunctionObjects.so");

        log                 true;
        writeToFile         false;

        switchNonOrtho      yes;
        switchGradCc        no;
        switchResiduals     yes;
        switchFaceWeight    no;
        switchSkewness      no;
        switchCo            no;

        maxNonOrthogonality 20;
        minNonOrthogonality 60;

        maxGradCc           3;
        minGradCc           4;

        maxFaceWeight       0.3;
        minFaceWeight       0.2;

        maxSkewness         2;
        minSkewness         3;

        maxResidual         10;

        result              UBlendingFactor;
        residual            initialResidual:p;
        P                   1.5;
        I                   0;
        D                   0.5;

        Co1                 1;
        Co2                 10;
        ...
        ...
        field               U;
    }

    Example of function object specification to calculate the residuals used
    by stabilityBlendingFactor. The following writes 'initialResidual:p'
    field

        residuals
        {
            type            residuals;
            libs            ("libutilityFunctionObjects.so");
            writeFields     true;
            writeControl    outputTime;
            fields          (p);
        }

    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description             | Required    | Default value
        type         | Type name: stabilityBlendingFactor | yes       |

        log          | Log to standard output  | no          | yes
        writeToFile  | Log to file             | no          | yes

        switchNonOrtho | non-orthogonal method | no | false
        switchGradCc | cell centre gradient method   | no | false
        switchResiduals | residual evolution method   | no | false
        switchFaceWeight | face weight method | no    | false
        switchSkewness | skewness method | no    | false
        switchCo       | Co blended |  no        | false

        maxNonOrthogonality| maximum non-orthogonal for scheme2 | no | 20
        minNonOrthogonality| minimum non-orthogonal for scheme1 | no | 60

        maxGradCc| maximum gradient for scheme2 | no | 2
        minGradCc| minimum gradient for scheme1 | no | 4

        maxResidual| maximum residual-mean ratio for scheme1 | no | 10
        P       |   proportional factor for PID     | no | 3
        I       |   integral factor for PID         | no | 0
        D       |   differential factor for PID     | no | 0.25

        maxFaceWeight | maximum face weight for scheme1 | no | 0.2
        minFaceWeight | minimum face weight for scheme2 | no | 0.3

        maxSkewness | maximum skewness for scheme2 | no | 2
        minSkewness | minimum skewness for scheme1 | no | 3


        faceWeight | Name of the faceWeight field | no | faceWeight
        skewness | Name of the skewness field | no | skewness
        nonOrthogonality | Name of the non-orthogonal field | no |
            nonOrthoAngle
        residual    | Name of the residual field | no | initialResidual:p
        U           | Name of the flux field for Co blended | no | U


        Co1 | Courant number below which scheme2 is used | no | 1
        Co2 | Courant number above which scheme1 is used | no | 10


        tolerance    | Tolerance for number of blended cells | no | 0.001
        field        | Name of field to evaluate | yes       |
        result    | Name of surface field to be used in the localBlended scheme
        | yes |
    \endtable

    The 'log' flag true write the number of cells on each scheme, plus it
    writes the field name "blendedIndicator".

    The 'result' entry is the field which is read by the localBlended scheme
    specified in fvSchemes. This name is determined by the localBlended class.


See also
    Foam::functionObjects::fieldExpression
    Foam::functionObjects::fvMeshFunctionObject
    Foam::functionObjects::writeFile

SourceFiles
    stabilityBlendingFactor.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_stabilityBlendingFactor_H
#define functionObjects_stabilityBlendingFactor_H

#include "fieldExpression.H"
#include "writeFile.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class stabilityBlendingFactor Declaration
\*---------------------------------------------------------------------------*/

class stabilityBlendingFactor
:
    public fieldExpression,
    public writeFile
{
    // Private member data

        //- Cell based blended indicator
        volScalarField indicator_;

        // Switches

            //- Switch for non-orthogonality
            Switch nonOrthogonality_;

            //- Switch for grad of cell centres
            Switch gradCc_;

            //- Switch for residuals
            Switch residuals_;

            //- Switch for face weight
            Switch faceWeight_;

            //- Switch for skewness
            Switch skewness_;

            //- Switch for Co
            Switch Co_;


        // Lower and upper limits

            //- Maximum non-orthogonality for fully scheme 2
            scalar maxNonOrthogonality_;

            //- Minimum non-orthogonality for fully scheme 1
            scalar minNonOrthogonality_;

            //- Maximum  gradcc for fully scheme 2
            scalar maxGradCc_;

            //- Minimum  gradcc for fully scheme 1
            scalar minGradCc_;

            //- Maximum ratio to average residual for scheme 2
            scalar maxResidual_;

            //- Minimum face weight for fully scheme 2
            scalar minFaceWeight_;

            //- Maximum face weight for fully scheme 1
            scalar maxFaceWeight_;

            //- Maximum skewness for fully scheme 2
            scalar maxSkewness_;

            //- Minimum skewness for fully scheme 1
            scalar minSkewness_;

            //- Maximum Co for fully scheme 2
            scalar Co1_;

            //- Minimum Co for fully scheme 1
            scalar Co2_;


        // File names

            //- Name of the non-orthogonalit field
            word nonOrthogonalityName_;

            //- Name of the face weight field
            word faceWeightName_;

            //- Name of the skewnes field
            word skewnessName_;

            //- Name of the residual field
            word residualName_;

            //- Name of the U used for Co based blended
            word UName_;


        //- Tolerance used when calculating the number of blended cells
        scalar tolerance_;


        //- Error fields
        scalarField error_;
        scalarField errorIntegral_;
        scalarField oldError_;
        scalarField oldErrorIntegral_;

        //- Proportional gain
        scalar P_;

        //- Integral gain
        scalar I_;

        //- Derivative gain
        scalar D_;


    // Private Member Functions

        //- Init fields
        bool init(bool first);

        //- Calculate statistics
        void calcStats(label&, label&, label&) const ;

        //- Calculate the blending factor field and return true if successful
        virtual bool calc();


protected:

    // Protected Member Functions

        //- Write the file header
        virtual void writeFileHeader(Ostream& os) const;


public:

    //- Runtime type information
    TypeName("stabilityBlendingFactor");


    // Constructors

        //- Construct from Time and dictionary
        stabilityBlendingFactor
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~stabilityBlendingFactor();


    // Member Functions

        //- Read the stabilityBlendingFactor data
        virtual bool read(const dictionary&);

        //- Write the stabilityBlendingFactor
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
