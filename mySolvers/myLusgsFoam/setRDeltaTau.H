{
    volScalarField& rDeltaTau = trDeltaTau.ref();

	 volScalarField rDeltaTau0("rDeltaTau0", rDeltaTau);
    
    // Set the reciprocal time-step from the local Courant number
    rDeltaTau.ref() = max(
        1/dimensionedScalar("maxDeltaT", dimTime, maxDeltaT),
        fvc::surfaceSum(amaxSf)()() /
        ((2*maxCo)*mesh.V()) 
    );
    
    // Update the boundary values of the reciprocal time-step
    rDeltaTau.correctBoundaryConditions();

    Info<< "Flow time scale min/max = "
        << gMin(1/rDeltaTau.primitiveField())
        << ", " << gMax(1/rDeltaTau.primitiveField()) << endl;
    
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
    
    Info<< "Smoothed flow time scale min/max = "
        << gMin(1/rDeltaTau.primitiveField())
        << ", " << gMax(1/rDeltaTau.primitiveField()) << endl;

    scalar rDeltaTauDampingCoeff
        (
            runTime.controlDict().lookupOrDefault<scalar>
            (
                "rDeltaTauDampingCoeff",
                1.0
            )
        );

	  // Limit rate of change of time scale
     // - reduce as much as required
     // - only increase at a fraction of old time scale
     if
     (
         rDeltaTauDampingCoeff < 1.0
      && runTime.timeIndex() > runTime.startTimeIndex() + 1
     )
     {
         rDeltaTau =
             rDeltaTau0
            *max(rDeltaTau/rDeltaTau0, scalar(1) - rDeltaTauDampingCoeff);

         Info<< "Damped flow time scale min/max = "
             << gMin(1/rDeltaTau.primitiveField())
             << ", " << gMax(1/rDeltaTau.primitiveField()) << endl;
     }
}
