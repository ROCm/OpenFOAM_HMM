// As per "setRootCase.H", but used when command-arguments are non-mandatory

Foam::argList::nonMandatoryArgs();

#include "listOptions.H"

Foam::argList args(argc, argv, false); // With suppressed check for arguments
if (!args.checkRootCase())
{
    Foam::FatalError.exit();
}

#include "listOutput.H"

//// Any subsequent user code will need to perform any pending checks
//// on arguments. Eg,
//
//  if (!args.check(true, false))
//  {
//      FatalError.exit();
//  }
