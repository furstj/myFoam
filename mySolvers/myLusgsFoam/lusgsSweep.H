{
    Field<scalar> D( mesh.V()/dt );
    if (LTS) {
      D *= dt * trDeltaTau();
    }

    scalar omega = 1.1;   // Relaxation factor 1<omega<=2

    const volScalarField kappa( thermo.Cp() / thermo.Cv() );
    const volScalarField a( sqrt( kappa / thermo.psi() ) );
    const volScalarField nuMax( "nuMax", max(
        4./3.0 * turbulence->nuEff(),
        kappa * turbulence->alphaEff()/rho
    ));

    // ALE mesh velocity + velocity due to MRF
    surfaceScalarField meshPhi( dbnsFlux.meshPhi() );
    MRF.makeAbsolute(meshPhi);

    Field<scalar> dRhoStar( dRho );
    Field<vector> dRhoUStar( dRhoU );
    Field<scalar> dRhoEStar( dRhoE );
    
    #include "forwardSweep.H"
    #include "backwardSweep.H"
}
