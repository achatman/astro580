import math
import sys
# ***********************************************
#   old fortran   SUBROUTINE EOS (X,Y,TIN,VIN,P,E,PE,PV,PT,EV,ET)
# ***********************************************
#
def eos(X,Y,TIN,VIN):
# input X, Y, T, V (=1/rho, cgs)
# returns  P,E,PE,PV,PT,EV,ET
    R = 8.31434E7
    A = 7.56471E-15
    BK = 8.6170837E-5
    AVAGD = 6.02217E23
    AD3 = A/3.E0
    T3OUT = 1.665795163E-25
    T4OUT = 3.802592017E-28
    T2OUT = 5.347896E-35
    T5OUT = 6.614536E-34
#Ionization Potentials for Hydrogen and Helium
    XH = 13.595
    XHE = 24.581
    XHE2 = 54.403
#
    C1 = 4.009296e-9
    C2 = 1.00797
    C3 = 4.0026
    XM = 7.9
    CM = 0.7
    ZPZP = 0.12014
#
    PREC = 1.E-10
    ONHLF = 1.5
#
#
    V = VIN
    T = TIN
    Z = 1-X-Y
    FRE = 0.
    ENT = 0.
    PARGRM = X/C2 + Y/C3
    RMUC = 1.0/PARGRM
    RT=R*T
    TT4=T**4
    TK=1.0/(T*BK)
    SQT=math.sqrt(T)
# C1=ORIGINAL(C1(0.33334622))/R
    T1 = V*SQT**3*C1
    T2 = T2OUT
#
    if T > 2e3:
        T2 = math.exp(-XH*TK)
    T3 = T3OUT
    if T > 5e3:
        T3 = math.exp(-XHE*TK)
    T4 = T4OUT
    if  T > 10e3:
        T4 = math.exp(-XHE2*TK)
    T5 = T5OUT
    if T > 1.2e3:
        T5 = math.exp(-XM*TK)
#
    D=T1*T2
    B=4.0*T1*T3
    C=B*T1*T4
    DD=2.0*CM*T1*T5
    ZNA=Z*2.48e-3/24.969
    ZMG=Z*ZPZP/45.807
# CONVERGE ON ELECTRON DENSITY USING THE SAHA EQUATION.
    GES = (X+Y*0.5)/(1.0+Y/(4.0*C))
    if GES < X:
        GES = 0.5*(math.sqrt(D*(D+4.0*X))-D)
    if GES < 1.e-6*Z:
        GES = 1.e-6*Z
    XC2 = X/C2
    YC3 = Y/C3
# NEWTON METHOD FOR ELECTRON DENSITY.
    I=0
    while I <= 25:
        I=I+1
        T2 = C/GES+GES+B
        GEP = XC2*D/(GES+D)+YC3*(B+2.0*C/GES)/T2 + ZMG*DD/(GES+DD) + ZNA
        T1 = 1.0+XC2*D/(D+GES)**2+YC3/T2*(2.0*C/GES**2+(B+2.0*C/GES)*(1.0-C/GES**2)/T2)+ ZMG*DD/(GES+DD)**2
        DGES = (GEP-GES)/T1
        GES = GES+DGES
        if math.fabs(DGES)/GES < PREC:
            I = 27
    if I == 26:
# you've run out of attempts.  Stop before you get hurt too much
         print("No convergence in EOS T=", T, " V=", V,  "DGES/GES", math.fabs(DGES)/GES)
         sys.exit(1)
# otherwise you are done - lets compute some more fun stuff
#   ELECTRON PRESSURE
    PE=RT*GES/V
    TOTLN=PARGRM+GES
    XX=D/(GES+D)
    T2=GES+B+C/GES
    YY=B/T2
    ZZ=C/(GES*T2)
    WW=DD/(GES+DD)
#     DERIVATIVES OF THE SAHA EQUATION
    T1 = YC3*(B+2.0*C/GES)
    QC0 = 1.0+XC2*XX/(GES+D)+ZMG*WW/(GES+DD)+YC3/T2*(2.0*C/GES**2+(B+2.0*C/GES)*(1.0-C/GES**2)/T2)
    QC1 = XC2*(1.0-XX)/(GES+D)
    QC4 = ZMG*(1.0-WW)/(GES+DD)
    QC2 = (YC3-T1/T2)/T2
    QC3 = (YC3*2.0-T1/T2)/(GES*T2)
    QGV = (QC1*D+QC2*B+QC3*2.0*C+QC4*DD)/(QC0*V)
    QP1 = D*(1.5+XH*TK)/T
    QP2 = B*(1.5+XHE*TK)/T
    QP3 = C*(3.0+(XHE+XHE2)*TK)/T
    QP4 = DD*(1.5+XM*TK)/T
    QGT = (QC1*QP1+QC2*QP2+QC3*QP3+QC4*QP4)/QC0
#          ELECTRON PRESSURE DERIVATIVES.
    PET=1.0+QGT/GES
    PEV=-1.0+QGV/GES
#         PRESSURE DUE TO THE IDEAL GAS
    P=RT*TOTLN/V
    PT=P/T+RT*QGT/V
    PV=RT*QGV/V-P/V
    BP = R*TOTLN
    BPV = R*QGV
    BPT = R*QGT
#           ADD THE RADIATION PRESSURE
    P = P+AD3*TT4
    PT = PT+4.0*AD3*TT4/T
#          IONIZATION ENERGY
    EI = (R/BK)*(XH*XX*XC2+YC3*(XHE*YY+(XHE+XHE2)*ZZ)+ZMG*XM*WW + ZNA*5.524 )
#         TOTAL INTERNAL ENERGY
    E=1.5*RT*TOTLN+A*V*TT4+EI
    EV = T*PT-P
    QXT = ((1.0-XX)*QP1-XX*QGT)/(GES+D)
    DT2 = QGT*(1.0-C/GES**2)+QP2+QP3/GES
    QYT = (QP2-B*DT2/T2)/T2
    QZT = (QP3-C*QGT/GES-C*DT2/T2)/(T2*GES)
    QWT = ((1.0-WW)*QP4-WW*QGT)/(GES+DD)
    EIT = (R/BK)*(XH*QXT*XC2+YC3*(XHE*QYT+(XHE+XHE2)*QZT)+ZMG*XM*QWT)
    ET = 1.5*R*(TOTLN+T*QGT)+4.0*A*V*TT4/T+EIT
#
# Return derived EOS quantities
    return P,E,PE,PV,PT,EV,ET
#
#-----------
#
# Main program below - modify and loop as needed...
#

TIN = float(input("Enter T (in K) "))
RHOIN = float(input("Enter rho (g/cm3) "))
VIN = 1.0/RHOIN
X = float(input("Enter X "))
Y = float(input("Enter Y "))
# Call EOS function to return values of thermo variables given inputs
P,E,PE,PV,PT,EV,ET = eos(X,Y,TIN,VIN)
#

print(f"Pressure: P = {P} dynes")
print(f"Internal Energy: E = {E} ergs/g")
print(f"PE: {PE}", f"PV: {PV}", f'PT: {PT}', f'EV: {EV}', f'ET: {ET}')

print("chiT= ",PT*TIN/P)
print("chirho= ",-PV*VIN/P)
print("Cv= ",EV)

gamma1 = -(VIN/P)*PV
gamma2 = 1/(1 - (P/TIN)/PT)
gamma3 = gamma1 * (P/TIN)/PT + 1
print(f'Gamma_1: {gamma1}')
print(f'Gamma_2: {gamma2}')
print(f'Gamma_3: {gamma3}')

