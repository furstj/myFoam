{
    volScalarField& he = thermo.he();

    surfaceScalarField mshPhi(
        IOobject("meshPhi", mesh.time().timeName(), mesh), 
        mesh, dimensionedScalar("0", dimVolume/dimTime, 0.0));

    MRF.makeAbsolute(mshPhi);

    fvScalarMatrix EEqn
    (
        fvm::div(phi, he)
      + (
            he.name() == "e"
          ? fvc::div(phi, volScalarField("Ekp", 0.5*magSqr(U) + p/rho))
          : fvc::div(phi, volScalarField("K", 0.5*magSqr(U)))
        )
      + thermophysicalTransport->divq(he)
      + fvc::div(mshPhi,p)
     ==
        fvModels.source(rho, he)
    );

    EEqn.relax();

    fvConstraints.constrain(EEqn);

    EEqn.solve();

    fvConstraints.constrain(he);

    thermo.correct();
}
