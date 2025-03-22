$PROB
# Model: 'Epidural sufentanil'

$PARAM @annotated
THETA1   : 66.2  : [L/day]
THETA2   : 44.0  : [L/day]
THETA3   : 9.35  : [L]
THETA4   : 504  : [L]
THETA5   : 1.99  : [1/day]
THETA6 : 10.9  : [1/day]
THETA7 : 0.162  : [1/day]

$OMEGA @annotated 
ECL: 0.0 : ETA

$MAIN
double CL = THETA1*exp(ECL);
double Q  = THETA2;
double V1 = THETA3;
double V2 = THETA4;
double KA12 = THETA5;
double KA14 = THETA6;
double KA41 = THETA7;
double K20 = CL/V1;
double K23 = Q/V1;
double K32 = Q/V2;
double ASS = 3/KA12;
  
$CMT @annotated
A1 : abs1
A2 : central
A3 : peripheral
A4 : abs2

$ODE

dxdt_A1 = KA12*ASS+KA14*ASS*exp(-KA41*SOLVERTIME)-KA12*A1-KA14*A1+KA41*A4;
dxdt_A2 = KA12*A1-(K20+K23)*A2+K32*A3;
dxdt_A3 = K23*A2-K32*A3;
dxdt_A4 = KA14*A1-KA41*A4;

$ERROR
double ABS=A1;
double CENT=A2/V1*1000;

$CAPTURE
ABS CENT
