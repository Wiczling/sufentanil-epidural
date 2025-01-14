$PROB
# Model: 'Epidural sufentanil'

$PARAM @annotated
THETA1   : 67.9  : [L/day]
THETA2   : 112  : [L/day]
THETA3   : 7.18  : [L]
THETA4   : 1280  : [L]
THETA5   : 1.07  : [1/day]
$OMEGA @annotated 
ECL: 0.0 : ETA

$MAIN
double CL = THETA1*exp(ECL);
double Q  = THETA2;
double V1 = THETA3;
double V2 = THETA4;
double KA12 = THETA5;
double K20 = CL/V1;
double K23 = Q/V1;
double K32 = Q/V2;
  
$CMT @annotated
A1 : abs1
A2 : central
A3 : peripheral

$ODE
dxdt_A1 = -KA12*A1;
dxdt_A2 = KA12*A1-(K20+K23)*A2+K32*A3;
dxdt_A3 = K23*A2-K32*A3;

$ERROR
double ABS=A1;
double CENT=A2/V1*1000;

$CAPTURE
ABS CENT
