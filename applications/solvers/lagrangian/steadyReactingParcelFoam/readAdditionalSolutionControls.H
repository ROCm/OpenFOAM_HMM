dictionary additional = mesh.solutionDict().subDict("additional");

bool eWork = additional.lookupOrDefault("eWork", true);
bool hWork = additional.lookupOrDefault("hWork", true);

// flag to activate solve transport for each specie (Y vector)
bool solveSpecies = additional.lookupOrDefault("solveSpecies", true);
