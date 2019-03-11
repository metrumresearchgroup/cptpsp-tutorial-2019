//Voriconazole PBPK model for a typical adult male

[PARAM] 
  //Tissue volumes (L); source: https://www.ncbi.nlm.nih.gov/pubmed/14506981
  Vad = 18.2 //adipose
  Vbo = 10.5 //bone
  Vbr = 1.45 //brain
  VguWall = 0.65  //gut wall
  VguLumen = 0.35 //gut lumen
  Vhe = 0.33 //heart
  Vki = 0.31 //kidneys
  Vli = 1.8 //liver
  Vlu = 0.5 //lungs
  Vmu = 29 //muscle
  Vsp = 0.15 //spleen
  Vbl = 5.6  //blood
  
  
  //Tissue blood flows (L/h); Cardiac output = 6.5 (L/min); source: https://www.ncbi.nlm.nih.gov/pubmed/14506981
  Qad = 0.05*6.5*60
  Qbo = 0.05*6.5*60
  Qbr = 0.12*6.5*60
  Qgu = 0.15*6.5*60 
  Qhe = 0.04*6.5*60
  Qki = 0.19*6.5*60
  Qmu = 0.17*6.5*60
  Qsp = 0.03*6.5*60
  Qha = 0.065*6.5*60  //hepatic artery
  Qlu = 6.5*60        //same as cardiac output
  
  //partition coefficients estimated by Poulin and Theil method https://jpharmsci.org/article/S0022-3549(16)30889-9/fulltext
  Kpad = 9.89  //adipose:plasma
  Kpbo = 7.91  //bone:plasma
  Kpbr = 7.35  //brain:plasma
  Kpgu = 5.82  //gut:plasma
  Kphe = 1.95  //heart:plasma
  Kpki = 2.9   //kidney:plasma
  Kpli = 4.66  //liver:plasma
  Kplu = 0.83  //lungs:plasma
  Kpmu = 2.94  //muscle:plasma; optimized
  Kpsp = 2.96  //spleen:plasma
  Kpre = 4     //calculated as average of non adipose Kps
  BP = 1       //blood:plasma ratio
  
  //other parameters
  WEIGHT = 73 //(kg)
  ka = 0.849  //absorption rate constant (/hr) 
  fup = 0.42  //fraction of unbound drug in plasma
  
  //in vitro hepatic clearance parameters http://dmd.aspetjournals.org/content/38/1/25.long
  fumic = 0.711 //fraction of unbound drug in microsomes
  MPPGL = 30.3  //adult mg microsomal protein per g liver (mg/g)
  VmaxH = 40    //adult hepatic Vmax (pmol/min/mg)
  KmH = 9.3     //adult hepatic Km (uM)


[CMT] 
  GUTLUMEN GUT ADIPOSE BRAIN HEART BONE 
  KIDNEY LIVER LUNG MUSCLE SPLEEN REST 
  ART VEN


[MAIN]
  //additional volume derivations
  double Vgu = VguWall + VguLumen;  //total gut volume
  double Vve = 0.705*Vbl; //venous blood
  double Var = 0.295*Vbl; //arterial blood
  double Vre = WEIGHT - (Vli+Vki+Vsp+Vhe+Vlu+Vbo+Vbr+Vmu+Vad+VguWall+Vbl); //volume of rest of the body compartment
  
  //additional blood flow derivation
  double Qli = Qgu + Qsp + Qha;
  double Qtot = Qli + Qki + Qbo + Qhe + Qmu + Qad + Qbr;
  double Qre = Qlu - Qtot;
  
  //intrinsic hepatic clearance calculation
  double scale_factor_H = MPPGL*Vli*1000; //hepatic scale factor (mg)
  double CLintHep = (VmaxH/KmH)*scale_factor_H*60*1e-6; //(L/hr)
  CLintHep = CLintHep/fumic; 
  
  //renal clearance https://link.springer.com/article/10.1007%2Fs40262-014-0181-y
  double CLrenal = 0.096; //(L/hr)


[ODE]
  //Calculation of tissue drug concentrations (mg/L)
  double Cadipose = ADIPOSE/Vad;
  double Cbone = BONE/Vbo;
  double Cbrain = BRAIN/Vbr; 
  double Cheart = HEART/Vhe; 
  double Ckidney = KIDNEY/Vki;
  double Cliver = LIVER/Vli; 
  double Clung = LUNG/Vlu; 
  double Cmuscle = MUSCLE/Vmu;
  double Cspleen = SPLEEN/Vsp;
  double Crest = REST/Vre;
  double Carterial = ART/Var;
  double Cvenous = VEN/Vve;
  double CgutLumen = GUTLUMEN/VguLumen;
  double Cgut = GUT/VguWall;
  
  //Free Concentration Calculations
  double Cliverfree = Cliver*fup; 
  double Ckidneyfree = Ckidney*fup;
  
  //ODEs
  dxdt_GUTLUMEN = -ka*GUTLUMEN;
  dxdt_GUT = ka*GUTLUMEN + Qgu*(Carterial - Cgut/(Kpgu/BP)); 
  dxdt_ADIPOSE = Qad*(Carterial - Cadipose/(Kpad/BP)); 
  dxdt_BRAIN = Qbr*(Carterial - Cbrain/(Kpbr/BP));
  dxdt_HEART = Qhe*(Carterial - Cheart/(Kphe/BP));
  dxdt_KIDNEY = Qki*(Carterial - Ckidney/(Kpki/BP)) - CLrenal*(Ckidneyfree/(Kpki/BP));
  dxdt_LIVER = Qgu*(Cgut/(Kpgu/BP)) + Qsp*(Cspleen/(Kpsp/BP)) + Qha*(Carterial) - Qli*(Cliver/(Kpli/BP)) - 
    CLintHep*(Cliverfree/(Kpli/BP)); 
  dxdt_LUNG = Qlu*(Cvenous - Clung/(Kplu/BP));
  dxdt_MUSCLE = Qmu*(Carterial - Cmuscle/(Kpmu/BP));
  dxdt_SPLEEN = Qsp*(Carterial - Cspleen/(Kpsp/BP));
  dxdt_BONE = Qbo*(Carterial - Cbone/(Kpbo/BP));
  dxdt_REST = Qre*(Carterial - Crest/(Kpre/BP));
  dxdt_VEN = Qad*(Cadipose/(Kpad/BP)) + Qbr*(Cbrain/(Kpbr/BP)) +
    Qhe*(Cheart/(Kphe/BP)) + Qki*(Ckidney/(Kpki/BP)) + Qli*(Cliver/(Kpli/BP)) + 
    Qmu*(Cmuscle/(Kpmu/BP)) + Qbo*(Cbone/(Kpbo/BP)) + Qre*(Crest/(Kpre/BP)) - Qlu*Cvenous;
  dxdt_ART = Qlu*(Clung/(Kplu/BP) - Carterial);


[TABLE]
  capture CP = Cvenous/BP;



