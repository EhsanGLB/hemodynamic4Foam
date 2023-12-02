{
    volScalarField DB("DB", kB*T/(3.0*3.14*turbulence->nu()*rho*dp) );
    volScalarField DT("DT", 0.26*(k/(2.0*k+kp))*(turbulence->nu()*C/T) );

    fvScalarMatrix CEqn
    ( 
          fvm::ddt(C)
        + fvm::div(phi, C)
        - fvm::laplacian(DB, C)
        - fvc::laplacian(DT, T)
    );

    CEqn.relax();
    CEqn.solve();
}

