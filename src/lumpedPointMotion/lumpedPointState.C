/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 OpenCFD Ltd.
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

#include "lumpedPointState.H"
#include "demandDrivenData.H"
#include "EulerCoordinateRotation.H"
#include "unitConversion.H"

#include "ISstream.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::lumpedPointState::inputFormatType
>
Foam::lumpedPointState::formatNames
{
    { inputFormatType::PLAIN, "plain" },
    { inputFormatType::DICTIONARY, "dictionary" }
};


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//! \cond fileScope
static Foam::string getLineNoComment
(
    Foam::ISstream& is,
    const char comment = '#'
)
{
    Foam::string line;
    do
    {
        is.getLine(line);
    }
    while ((line.empty() || line[0] == comment) && is.good());

    return line;
}
//! \endcond


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::lumpedPointState::calcRotations() const
{
    rotationPtr_ = new tensorField(angles_.size());
    forAll(angles_, itemi)
    {
        rotationPtr_->operator[](itemi) = EulerCoordinateRotation
        (
            angles_[itemi],
            degrees_  // true=degrees, false=radians
        ).R();
    }
}


void Foam::lumpedPointState::readDict(const dictionary& dict)
{
    dict.lookup("points") >> points_;
    dict.lookup("angles") >> angles_;
    degrees_ = dict.lookupOrDefault("degrees", false);
    deleteDemandDrivenData(rotationPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::lumpedPointState::lumpedPointState()
:
    points_(0),
    angles_(0),
    degrees_(false),
    rotationPtr_(nullptr)
{}


Foam::lumpedPointState::lumpedPointState(const lumpedPointState& rhs)
:
    points_(rhs.points_),
    angles_(rhs.angles_),
    degrees_(rhs.degrees_),
    rotationPtr_(nullptr)
{}


Foam::lumpedPointState::lumpedPointState
(
    const pointField& pts
)
:
    points_(pts),
    angles_(points_.size(), Zero),
    degrees_(false),
    rotationPtr_(nullptr)
{}


Foam::lumpedPointState::lumpedPointState
(
    tmp<pointField>& pts
)
:
    points_(pts),
    angles_(points_.size(), Zero),
    degrees_(false),
    rotationPtr_(nullptr)
{}


Foam::lumpedPointState::lumpedPointState
(
    const dictionary& dict
)
:
    points_(0),
    angles_(0),
    degrees_(false),
    rotationPtr_(nullptr)
{
    readDict(dict);
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

Foam::lumpedPointState::~lumpedPointState()
{
    deleteDemandDrivenData(rotationPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::lumpedPointState::operator=(const lumpedPointState& rhs)
{
    points_  = rhs.points_;
    angles_  = rhs.angles_;
    degrees_ = rhs.degrees_;

    deleteDemandDrivenData(rotationPtr_);
}


void Foam::lumpedPointState::relax
(
    const scalar alpha,
    const lumpedPointState& prev
)
{
    points_ = prev.points_ + alpha*(points_ - prev.points_);

    scalar convert = 1.0;
    if (degrees_ != prev.degrees_)
    {
        if (prev.degrees_)
        {
            // Was degrees, now radians
            convert = degToRad();
        }
        else
        {
            // Was radians, now degrees
            convert = radToDeg();
        }
    }

    angles_ = convert*prev.angles_ + alpha*(angles_ - convert*prev.angles_);

    deleteDemandDrivenData(rotationPtr_);
}


bool Foam::lumpedPointState::readPlain(Istream& is)
{
    // Assume generic input stream so we can do line-based input
    ISstream& iss = dynamic_cast<ISstream&>(is);

    string line = getLineNoComment(iss);

    label count;
    {
        IStringStream isstr(line);
        isstr >> count;
    }

    points_.setSize(count);
    angles_.setSize(count);

    count = 0;
    forAll(points_, i)
    {
        iss.getLine(line);
        IStringStream isstr(line);

        isstr
            >> points_[count].x() >> points_[count].y() >> points_[count].z()
            >> angles_[count].x() >> angles_[count].y() >> angles_[count].z();

        ++count;
    }

    points_.setSize(count);
    angles_.setSize(count);

    degrees_ = false;
    deleteDemandDrivenData(rotationPtr_);

    return count;
}


bool Foam::lumpedPointState::readData(Istream& is)
{
    dictionary dict(is);
    readDict(dict);

    return points_.size();
}


bool Foam::lumpedPointState::writeData(Ostream& os) const
{
    writeDict(os);
    return true;
}


void Foam::lumpedPointState::writeDict(Ostream& os) const
{
    os.writeEntry("points", points_);
    os.writeEntry("angles", angles_);
    if (degrees_)
    {
        os.writeEntry("degrees", word("true"));
    }
}


void Foam::lumpedPointState::writePlain(Ostream& os) const
{
    os  <<"# input for OpenFOAM\n"
        <<"# N, points, angles\n"
        << points_.size() << "\n";

    forAll(points_, i)
    {
        const vector& pt = points_[i];

        os  << pt.x() << ' '
            << pt.y() << ' '
            << pt.z() << ' ';

        if (i < angles_.size())
        {
            os  << angles_[i].x() << ' '
                << angles_[i].y() << ' '
                << angles_[i].z() << '\n';
        }
        else
        {
            os  << "0 0 0\n";
        }
    }
}


bool Foam::lumpedPointState::readData
(
    const inputFormatType& fmt,
    const fileName& file
)
{
    bool ok = false;
    if (Pstream::master())
    {
        IFstream is(file);

        if (fmt == inputFormatType::PLAIN)
        {
            ok = this->readPlain(is);
        }
        else
        {
            ok = this->readData(is);
        }
    }

    if (Pstream::parRun())
    {
        // Scatter master data using communication scheme

        const List<Pstream::commsStruct>& comms =
        (
            (Pstream::nProcs() < Pstream::nProcsSimpleSum)
          ? Pstream::linearCommunication()
          : Pstream::treeCommunication()
        );

        // Get my communication order
        const Pstream::commsStruct& myComm = comms[Pstream::myProcNo()];

        // Receive from up
        if (myComm.above() != -1)
        {
            IPstream fromAbove
            (
                UPstream::commsTypes::scheduled,
                myComm.above(),
                0,
                Pstream::msgType(),
                Pstream::worldComm
            );

            fromAbove >> points_ >> angles_ >> degrees_;
        }

        // Send to downstairs neighbours
        forAllReverse(myComm.below(), belowI)
        {
            OPstream toBelow
            (
                UPstream::commsTypes::scheduled,
                myComm.below()[belowI],
                0,
                Pstream::msgType(),
                Pstream::worldComm
            );

            toBelow << points_ << angles_ << degrees_;
        }

        deleteDemandDrivenData(rotationPtr_);

        // MPI barrier
        Pstream::scatter(ok);
    }

    return ok;
}


// ************************************************************************* //
