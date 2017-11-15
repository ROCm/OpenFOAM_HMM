/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2017 OpenCFD Ltd.
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

#include "fieldAverageItem.H"
#include "IOstreams.H"
#include "dictionaryEntry.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fieldAverageItem::fieldAverageItem(Istream& is)
:
    active_(false),
    fieldName_("unknown"),
    mean_(false),
    meanFieldName_("unknown"),
    prime2Mean_(false),
    prime2MeanFieldName_("unknown"),
    base_(baseType::ITER),
    totalIter_(0),
    totalTime_(-1),
    window_(-1.0),
    windowName_(""),
    windowType_(windowType::NONE),

    windowTimes_(),
    windowFieldNames_(),
    allowRestart_(true)
{
    is >> *this;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::functionObjects::operator>>
(
    Istream& is,
    fieldAverageItem& faItem
)
{
    is.check(FUNCTION_NAME);

    const dictionaryEntry entry(dictionary::null, is);

    faItem.active_ = false;
    faItem.fieldName_ = entry.keyword();
    faItem.mean_ = readBool(entry.lookup("mean"));
    faItem.prime2Mean_ = readBool(entry.lookup("prime2Mean"));
    faItem.base_ = faItem.baseTypeNames_.lookup("base", entry);
    faItem.window_ = entry.lookupOrDefault<scalar>("window", -1.0);

    if (faItem.window_ > 0)
    {
        faItem.windowType_ =
            faItem.windowTypeNames_.lookup("windowType", entry);

        if (faItem.windowType_ != fieldAverageItem::windowType::NONE)
        {
            if
            (
                faItem.base_ == fieldAverageItem::baseType::ITER
             && label(faItem.window_) < 1
            )
            {
                FatalIOErrorInFunction(entry)
                    << "Window must be 1 or more for base type "
                    << faItem.baseTypeNames_[fieldAverageItem::baseType::ITER]
                    << exit(FatalIOError);
            }

            faItem.windowName_ = entry.lookupOrDefault<word>("windowName", "");

            if (faItem.windowType_ == fieldAverageItem::windowType::EXACT)
            {
                faItem.allowRestart_ = readBool(entry.lookup("allowRestart"));

                if (!faItem.allowRestart_)
                {
                    WarningInFunction
                        << faItem.windowTypeNames_[faItem.windowType_]
                        << " windowing for field " << faItem.fieldName_
                        << " will not write intermediate files and restart will"
                        << " not be possible.  To enable, please set"
                        << " 'allowRestart' to 'yes'"
                        << endl;
                }
            }
        }
        else
        {
            // Deactivate windowing
            faItem.window_ = -1;
        }
    }
    else
    {
        // Deactivate windowing
        faItem.window_ = -1;
    }

    faItem.meanFieldName_ = faItem.fieldName_ + fieldAverageItem::EXT_MEAN;
    faItem.prime2MeanFieldName_ =
        faItem.fieldName_ + fieldAverageItem::EXT_PRIME2MEAN;

    if ((faItem.window_ > 0) && (!faItem.windowName_.empty()))
    {
        faItem.meanFieldName_ =
            faItem.meanFieldName_ + "_" + faItem.windowName_;

        faItem.prime2MeanFieldName_ =
            faItem.prime2MeanFieldName_ + "_" + faItem.windowName_;
    }

    return is;
}


Foam::Ostream& Foam::functionObjects::operator<<
(
    Ostream& os,
    const fieldAverageItem& faItem
)
{
    os.check(FUNCTION_NAME);

    os.beginBlock(faItem.fieldName_);

    os.writeEntry("mean", faItem.mean_);
    os.writeEntry("prime2Mean", faItem.prime2Mean_);
    os.writeEntry("base", faItem.baseTypeNames_[faItem.base_]);

    if (faItem.window_ > 0)
    {
        os.writeEntry("window", faItem.window_);

        if (!faItem.windowName_.empty())
        {
            os.writeEntry("windowName", faItem.windowName_);
        }

        os.writeEntry
        (
            "windowType",
            faItem.windowTypeNames_[faItem.windowType_]
        );

        os.writeEntry("allowRestart", faItem.allowRestart_);
    }

    os.endBlock() << flush;

    os.check(FUNCTION_NAME);

    return os;
}


// ************************************************************************* //
