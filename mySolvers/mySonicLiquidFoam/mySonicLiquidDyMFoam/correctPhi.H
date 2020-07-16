#if (OPENFOAM >= 1912)
CorrectPhi
(
    U,
    phi,
    p,
    rho,
    psiField,
    dimensionedScalar("rAUf", dimTime, 1),
    divrhoU,
    pimple
);
#else
CorrectPhi
(
    U,
    phi,
    p,
    rho,
    psiField,
    dimensionedScalar("rAUf", dimTime, 1),
    divrhoU,
    pimple,
    true
);
#endif
