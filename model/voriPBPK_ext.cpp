//Voriconazole PBPK model for a typical adult male


$PROB voriPBPK


$PARAM 
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
Qlu = 6.5*60  //same as cardiac output

//partition coefficients estimated by Poulin and Theil method https://jpharmsci.org/article/S0022-3549(16)30889-9/fulltext
Kpad = 9.89  //adipose:plasma
Kpbo = 7.91  //bone:plasma
Kpbr = 7.35  //brain:plasma
Kpgu = 5.82  //gut:plasma
Kphe = 1.95  //heart:plasma
Kpki = 2.9  //kidney:plasma
Kpli = 4.66  //liver:plasma
Kplu = 0.83  //lungs:plasma
Kpmu = 2.94  //muscle:plasma; optimized
Kpsp = 2.96  //spleen:plasma
Kpre = 4 //calculated as average of non adipose Kps
BP = 1 //blood:plasma ratio; optimized

//other parameters
WEIGHT = 73 //(kg)
Ka = 0.849 //absorption rate constant(/hr) 
fup = 0.42 //fraction of unbound drug in plasma

//in vitro hepatic clearance parameters http://dmd.aspetjournals.org/content/38/1/25.long
fumic = 0.711 //fraction of unbound drug in microsomes
MPPGL = 30.3 //adult mg microsomal protein per g liver (mg/g)
VmaxH = 40 //adult hepatic Vmax (pmol/min/mg)
KmH = 9.3 //adult hepatic Km (uM)

//in vitro intestinal clearance parameters
MPPGI = 0 //adult mg microsomal protein per g intestine (mg/g)
VmaxG = 40 //adult intestinal Vmax (pmol/min/mg)
KmG = 9.3 //adult intestinal Km (uM)

//absorption model parameters
MW = 349.317  //(g/mol)
logP = 2.56  //log10 octanol oil:water partition coefficient; will be used as proxy for membrane affinity; preferably we will have phospholipid bilayer:water partition instead
S_lumen = 0.39*1000  //(mg/L) voriconazole intestinal lumen solubility https://www.ncbi.nlm.nih.gov/pubmed/24557773
L = 280  //(cm) small intestine length; from ICRP Publication 89
d = 2.5  //(cm) diameter of small intestine lumen
PF = 1.57  //3  //2.29  //average of 1.57 and 3; 1.57  //plicae circulare factor https://www.ncbi.nlm.nih.gov/pubmed/24694282
VF = 6.5  //villi factor
MF = 13  //microvilli factor
ITT = 3.32  //(h) small intestine transit time; https://www.ncbi.nlm.nih.gov/pubmed/25986421
A = 7440  //this and the rest of parameters are constants in the permeability calculation equation https://www.ncbi.nlm.nih.gov/pubmed/15267240
B = 1e7
alpha = 0.6
beta = 4.395
fabs = 1  //absorption factor to manipulate ka
fdis = 1  //disappearance from gut lumen factor to manipulate kd
fperm = 1  //permeability factor to manipulate Pm


$CMT 
GUTLUMEN GUTWALL GUT ADIPOSE BRAIN HEART BONE 
KIDNEY LIVER LUNG MUSCLE SPLEEN REST 
ART VEN


$MAIN
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

//intrinsic intestinal clearance calculation
double scale_factor_G = MPPGI*VguWall*1000; //intestinal scale factor (mg) 
double CLintGut = (VmaxG/KmG)*scale_factor_G*60*1e-6; //(L/hr)
CLintGut = CLintGut/fumic;

//renal clearance https://link.springer.com/article/10.1007%2Fs40262-014-0181-y
double CLrenal = 0.096; //(L/hr)

//absorption model parameters derivation
double SA_abs = M_PI*L*d*PF*VF*MF*1e-4;  //(m^2) mucosal absorption surface area https://www.ncbi.nlm.nih.gov/pubmed/24694282
double SA_basal = M_PI*L*d*PF*VF*1e-4;  //(m^2) basal membrane surface area https://www.ncbi.nlm.nih.gov/pubmed/24694282
double MA = pow(10,logP);  //membrane affinity
double MW_eff = MW - (3*17);  //effective molecular weight; voriconazole has 3 F atoms so we subtract 17 mass units per atom https://www.ncbi.nlm.nih.gov/pubmed/15267240
double Peff = fperm*A*((pow(MW_eff,(-alpha-beta))*MA)/(pow(MW_eff,(-alpha)) + B*pow(MW_eff,(-beta))*MA) * 1e-2 * 3600);  //(m/h) intestinal permeability 
double kd = fdis*Peff*SA_abs*1000/VguLumen;  //(h-1) rate constant for drug disappearing from lumen and into enterocytes
double ka = fabs*Peff*SA_basal*1000/VguWall; //(h-1) rate constant for drug absorption from enterocytes to gut circulation
double kt = 1/ITT;  //(h-1) intestinal transit rate constatnt


$ODE
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
double CgutWall = GUTWALL/VguWall;
double Cgut = GUT/VguWall;

//Free Concentration Calculations
double Cliverfree = Cliver*fup; 
double Ckidneyfree = Ckidney*fup;

//accounting for solubility
double f = 1;
if(CgutLumen > S_lumen){
  f = 0;
}

//ODEs
dxdt_GUTLUMEN = - kd*VguLumen*(f*CgutLumen + (1-f)*S_lumen) - kt*GUTLUMEN;
dxdt_GUTWALL = kd*VguLumen*(f*CgutLumen + (1-f)*S_lumen) - ka*GUTWALL - CLintGut*CgutWall;
dxdt_GUT = ka*GUTWALL + Qgu*(Carterial - Cgut/(Kpgu/BP)); 
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


$CAPTURE Cvenous
  
  
  