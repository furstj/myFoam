dictionary lusgsDict = mesh.solutionDict().subDict("LUSGS");

label lusgsIntIters = 
    lusgsDict.lookupOrDefault<label>("internalIterations", 1);

scalar lusgsTolerance =
    lusgsDict.lookupOrDefault<scalar>("tolerance", 0.0);

scalar lusgsRelTol =
    lusgsDict.lookupOrDefault<scalar>("relTol", 0.0);

bool LTS = 
    lusgsDict.lookupOrDefault("localTimestepping", false);

bool turbOnFinalIterOnly =
    lusgsDict.lookupOrDefault<Switch>("turbOnFinalIterOnly", true);
