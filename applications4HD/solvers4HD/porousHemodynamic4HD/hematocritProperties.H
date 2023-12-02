// Variables for the haemo-particle-transport model

Info << "Reading hematocritProperties\n" << endl;

IOdictionary hematocritProperties
(
    IOobject
    (
        "hematocritProperties", // name of the dictionary
        runTime.constant(), // location in the case - this one is in constant
        mesh, // needs the mesh object reference to do some voodoo - unimportant now
        IOobject::MUST_READ_IF_MODIFIED, // the file will be re-read if it gets modified during time stepping
        IOobject::NO_WRITE // read-only
    )
);


dimensionedScalar a(hematocritProperties.lookup("a"));
dimensionedScalar Kc(hematocritProperties.lookup("Kc"));
dimensionedScalar Kmu(hematocritProperties.lookup("Kmu"));
dimensionedScalar haemoSwitch(hematocritProperties.lookup("haemoSwitch"));
dimensionedScalar haemoSwitchTime(hematocritProperties.lookup("haemoSwitchTime"));
dimensionedScalar haemoSw = haemoSwitch;
dimensionedScalar linKc(hematocritProperties.lookup("linKc"));
dimensionedScalar nonlinKc(hematocritProperties.lookup("nonlinKc"));
dimensionedScalar nonlinKcRatio(hematocritProperties.lookup("nonlinKcRatio"));


if (linKc.value() == 1.0 && nonlinKc.value() == 0.0)
{
    Info << "linear Kc selected" << nl << Kmu << nl << Kc << nl << endl;
}
else if (linKc.value() == 0.0 && nonlinKc.value() == 1.0)
{
    Info << "non-linear Kc selected" << nl << Kmu << nl << nonlinKcRatio << nl << endl;
}
else
{
    Info << "Inconsistent Kc settings! If this is what you want, that's fine." << nl << "But it may not give useful results:" << nl << Kmu << nl << Kc << nl << linKc << nl << nonlinKc << nl << nonlinKcRatio << nl << endl;
}

