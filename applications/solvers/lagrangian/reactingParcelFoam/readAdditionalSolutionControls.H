dictionary additional = mesh.solutionDict().subDict("additional");

// flag to activate solve transport for each specie (Y vector)
bool solveSpecies = additional.lookupOrDefault("solveSpecies", true);
