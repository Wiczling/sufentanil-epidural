$PROBLEM sufentanil_epidural

$INPUT ID TIME AMT RATE DV CMT MDV EVID SEX AGE WT HEIGHT TINFS TINFR DOSES DOSER OUT NUM

$DATA ../../../../data/derived/nonmem_data_v02.csv  
IGNORE=I
IGNORE=(OUT.EQ.1)
IGNORE=(CMT.EQ.22)
IGNORE=(CMT.EQ.33)

$SUBROUTINES ADVAN13 TOL=6

$MODEL 
  COMP  = (ABS)
  COMP  = (CENTRAL)
  COMP  = (PERIPHERAL)
  COMP  = (ABS2)

$PK
;PK PARAMETERS
 CL = THETA(1)*EXP(ETA(1))
 V1 = THETA(2)*EXP(ETA(2))
 Q  = THETA(3)*EXP(ETA(3))
 V2 = THETA(4)*EXP(ETA(4))
 KA = THETA(5)*EXP(ETA(5))
 K14= THETA(6)*EXP(ETA(6))
 K41= THETA(7)*EXP(ETA(7))

;SECONDARY PARAMETERS 
 K20 = CL/V1;
 K23 = Q/V1;
 K32 = Q/V2;
 
$DES
 DADT(1) = -KA*A(1)-K14*A(1)+K41*A(4)
 DADT(2) =  KA*A(1)-K23*A(2)+K32*A(3)-K20*A(2)
 DADT(3) =          K23*A(2)-K32*A(3)
 DADT(4) =          K14*A(1)-K41*A(4)
 
$ERROR

IPRED = A(2)/V1*1000 ;ug/L->ng/L=pg/ml
Y=IPRED*(1+EPS(2))+EPS(1)  


$THETA  
; https://pubmed.ncbi.nlm.nih.gov/26105145/
; VC = 7.90 l, VT = 481 L, Cl = 45.3 L/h, and Q = 38.3 L/h
45.3 FIX    ;[L/h] CL
7.90 FIX    ;[L] V1
38.3 FIX    ;[L/h] Q
481 FIX     ;[L] V2
(0,1)       ;[1/h] KA
(0,0.1)    ;[1/h] K14
(0,0.01)     ;[1/h] K41

$OMEGA     
0.1  ;[P] ETA-CL
0 FIX ;[P] ETA-V1
0 FIX ;[P] ETA-Q
0 FIX ;[P] ETA-V2
0 FIX   ;[P] ETA-KA
0 FIX   ;[P] ETA-K14
0 FIX     ;[P] ETA-K41

$SIGMA
(0.01 FIX)  ;[P] ADD pg/ml
(0.4)       ;[P] PROP

$EST MAXEVAL=9999 METHOD=1 INTER SIGL=6 NSIG=2 PRINT=1 MSFO=./1.msf
$COV PRINT=E MATRIX=R

$TABLE NUM ID TIME EVID MDV CMT PRED IPRED CWRES WRES NPDE NOPRINT NOAPPEND ONEHEADERALL
FILE=sdtab1.tab
$TABLE NUM ID ETAS(1:LAST) KA CL V1 Q V2 K14 K41 NOPRINT NOAPPEND ONEHEADERALL
FILE=patab1.tab