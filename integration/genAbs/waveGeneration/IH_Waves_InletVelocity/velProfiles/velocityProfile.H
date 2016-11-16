/*---------------------------------------------------------------------------*\
   IH-Cantabria 2015 (http://www.ihcantabria.com/en/)
   IHFOAM 2015 (http://ihfoam.ihcantabria.com/) 

   Author(s):  Javier Lopez Lara (jav.lopez@unican.es)
               Gabriel Barajas   (barajasg@unican.es)
\*---------------------------------------------------------------------------*/

if ( waveType_ == "regular" )
{
    if ( waveTheory_ == "StokesI" )
    {
        #include "profileStokesI.H"
        #include "profileGENAB.H"
    }
    else if ( waveTheory_ == "StokesII" )
    {
        #include "profileStokesII.H"
        #include "profileGENAB.H"
    }
    else if ( waveTheory_ == "StokesV" )
    {
        #include "profileStokesV.H"
        #include "profileGENAB.H"
    }
    else if ( waveTheory_ == "Cnoidal" )
    {
        #include "profileCnoidal.H"
        #include "profileGENAB.H"
    }
}
else if ( waveType_ == "solitary" )
{
    #include "profileBoussinesq.H"
    #include "profileGENAB.H"
}
