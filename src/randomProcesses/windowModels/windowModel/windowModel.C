#include "windowModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(windowModel, 0);
    defineRunTimeSelectionTable(windowModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::windowModel::windowModel(const dictionary& dict, const label nSamples)
:
    scalarField(nSamples),
    nOverlapSamples_(0),
    nWindow_(dict.lookupOrDefault("nWindow", -1))
{
    scalar prc = readScalar(dict.lookup("overlapPercent"));
    nOverlapSamples_ = floor(prc/scalar(100)*nSamples);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::windowModel::~windowModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::windowModel::nSamples() const
{
    return size();
}


Foam::label Foam::windowModel::nWindow() const
{
    return nWindow_;
}


Foam::label Foam::windowModel::nWindowsTotal(label nSamplesTotal) const
{
    const label nSamples = this->nSamples();

    return floor((nSamplesTotal - nSamples)/(nSamples - nOverlapSamples_)) + 1;
}


Foam::label Foam::windowModel::validate(const label nSamplesTotal)
{
    label nSamples = this->nSamples();

    if (nSamplesTotal < nSamples)
    {
        FatalErrorInFunction
            << "Block size N = " << nSamples
            << " is larger than total number of data points = " << nSamplesTotal
            << exit(FatalError);
    }

    label nWindowAvailable = nWindowsTotal(nSamplesTotal);

    if (nWindow_ == -1)
    {
        nWindow_ = nWindowAvailable;
    }

    if (nWindow_ > nWindowAvailable)
    {
        FatalErrorInFunction
            << "Number of data points calculated with " << nWindow_
            << " windows greater than the total number of data points"
            << nl
            << "    Block size, N = " << nSamples << nl
            << "    Total number of data points = " << nSamplesTotal << nl
            << "    Maximum number of windows = " << nWindowAvailable << nl
            << "    Requested number of windows = " << nWindow_
            << exit(FatalError);
    }

    label nRequiredSamples =
        nWindow_*nSamples - (nWindow_ - 1)*nOverlapSamples_;

    Info<< "Windowing:" << nl
        << "    Total samples              : " << nSamplesTotal << nl
        << "    Samples per window         : " << nSamples << nl
        << "    Number of windows          : " << nWindow_ << nl
        << "    Overlap size               : " << nOverlapSamples_ << nl
        << "    Required number of samples : " << nRequiredSamples
        << endl;

    return nRequiredSamples;
}


// ************************************************************************* //
