{
    volScalarField& rDeltaTau = trDeltaTau.ref();
    
    // Set the reciprocal time-step from the local Courant number
    rDeltaTau.ref() = max(
        1/dimensionedScalar("maxDeltaT", dimTime, maxDeltaT),
        fvc::surfaceSum(amaxSf)()() /
        ((2*maxCo)*mesh.V()) 
    );
    
    // Update tho boundary values of the reciprocal time-step
    rDeltaTau.correctBoundaryConditions();
    
    scalar rDeltaTauSmoothingCoeff
        (
            runTime.controlDict().lookupOrDefault<scalar>
            (
                "rDeltaTauSmoothingCoeff",
                0.02
            )
        );
    
    if (rDeltaTauSmoothingCoeff > 0.0) 
        fvc::smooth(rDeltaTau, rDeltaTauSmoothingCoeff);
    
    Info<< "Flow time scale min/max = "
        << gMin(1/rDeltaTau.primitiveField())
        << ", " << gMax(1/rDeltaTau.primitiveField()) << endl;
}
