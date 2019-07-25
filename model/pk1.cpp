[PARAM] CL=1, VC=20, KA=1

[CMT] GUT CENT

[ODE]
dxdt_GUT = -KA*GUT;
dxdt_CENT = KA*GUT - (CL/VC)*CENT;

[TABLE] capture CP = CENT/VC;
