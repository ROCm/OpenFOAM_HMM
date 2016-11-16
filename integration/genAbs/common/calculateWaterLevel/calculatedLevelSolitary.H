/*---------------------------------------------------------------------------*\
   IH-Cantabria 2015 (http://www.ihcantabria.com/en/)
   IHFOAM 2015 (http://ihfoam.ihcantabria.com/) 

   Author(s):  Javier Lopez Lara (jav.lopez@unican.es)
               Gabriel Barajas   (barajasg@unican.es)
\*---------------------------------------------------------------------------*/ 

   if ( waveTheory_ == "Boussinesq" )
   {
       forAll(calculatedLevel, it2)
       {
           calculatedLevel[it2] = RealwaterDepth_
               + BoussinesqFun :: eta
               (
                   waveHeight_, 
                   RealwaterDepth_,  
                   xGroup[it2], 
                   yGroup[it2], 
                   waveAngle, 
                   currTime, 
                   X0
               );
        }
   }
   else 
   {
        FatalError << "Wave theory not supported, use:\n" 
	    << "Boussinesq" << exit(FatalError);
   }
