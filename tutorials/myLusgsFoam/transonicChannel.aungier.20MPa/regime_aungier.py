#!/usr/bin/env python3

from math import *
from scipy.optimize import root

pTot = 20e6
Ttot = 400
M2i  = 0.675

def sign(s):
    if s>0: return 1
    if s<0: return -1
    return 0

class JANAF:
    # JANAF coeffs
    def __init__(self):
        self.a = [3.53101, -0.000123661, -5.02999e-07, 2.43531e-09, -1.40881e-12, -1046.98, 2.96747]
        self.R = 8.31446261815324 / 28.01348e-3

    def cp0(self, p, T):
        return self.R*((((self.a[4]*T + self.a[3])*T + self.a[2])*T + self.a[1])*T + self.a[0])
    
    def H0(self, p, T):
        return self.R*(((((self.a[4]/5*T + self.a[3]/4)*T + self.a[2]/3)*T + \
                         self.a[1]/2)*T + self.a[0])*T + self.a[5])
    
    def S0(self, p, T):
        return self.R*(((((self.a[4]/4*T + self.a[3]/3)*T + self.a[2]/2)*T + \
                         self.a[1])*T + self.a[0]*log(T)) + self.a[6])
    
    
        
class AungierRedlichKwong:
    def __init__(self, Pc = 3395800.0, Vc = 0.0031918,  Tc = 126.192, omega = 0.0372):
        self.R = 8.31446261815324 / 28.01348e-3
        self.Tc = Tc
        self.Vc = Vc
        self.Pc = Pc
        self.omega = omega
        self.n = 0.4986 + 1.1735*self.omega + 0.4754*self.omega**2
        self.a = 0.42747*(self.R*self.Tc)**2/self.Pc 
        self.b = 0.08664*self.R*self.Tc/self.Pc;
        self.c = self.R*self.Tc/(self.Pc + self.a/(self.Vc*(self.Vc + self.b))) + self.b - self.Vc;
        
    def p(self, rho, T):
        Tr = T/self.Tc
        alpha = self.a*Tr**(-self.n)
        v = 1.0/rho
        p = self.R*T/(v - self.b + self.c) - alpha/v/(v+self.b)
        return p
    
    def rho(self, p, T):
        return p/(self.Z(p,T)*self.R*T)

    def Z(self, p, T):
        Tr = T/self.Tc;
        RT = self.R*T;
        alpha = self.a*Tr**(-self.n)
        
        B = self.c*p/RT - 1;
        C = p*(-RT*self.b + alpha - self.b**2*p + self.b*self.c*p)/(RT)**2;
        D = alpha*p**2*(self.c - self.b)/(RT)**3;

        # Solve Z^3 + B*Z^2 + C*Z + D = 0, return maximal real root 
        root = 0;
        q = (3*C - B*B)/9;
        r = (-27*D + B*(9*C - 2*B*B))/54.0;
        disc = q**3 + r**2;
        term1 = B/3.0;
    
        if (disc > 0):
            # one root real, two are complex
            s = r + sqrt(disc);
            s = sign(s)*pow(abs(s), 1.0/3.0)
            t = r - sqrt(disc);
            t = sign(t)*pow(abs(t), 1.0/3.0)
            root = -term1 + s + t;
            
        elif (disc == 0):
            # All roots real, at least two are equal.
            r13 = sign(r)*pow(abs(r), 1.0/3.0)
            root = max(2*r13 - term1, -r13 - term1)
        else:
            # Only option left is that all roots are real and unequal (to get here, q < 0)
            dum1 = acos(r/sqrt(-q**3));
            r13 = 2.0*sqrt(-q);
            x1 = -term1 + r13*cos(dum1/3.0);
            x2 = -term1 + r13*cos((dum1 + 2.0*pi)/3.0);
            x3 = -term1 + r13*cos((dum1 + 4.0*pi)/3.0);
            root = max(x1, max(x2, x3));

        return root


    def dpdV(self, V, T, alpha):
        return (- self.R*T/(V - self.b + self.c)**2 + alpha/(V*(V + self.b)**2) + \
                alpha/(V**2*(V + self.b)))

    
    def dpdT(self, V, T, alpha):
        return self.R/(V - self.b + self.c) + alpha*self.n/(T*V*(V + self.b))

    
    def cpdep(self, p, T):
        Z = self.Z(p, T)
        V = Z*self.R*T/p
        Tr = T/self.Tc
        alpha = self.a*pow(Tr, -self.n)

        dVdT = - self.dpdT(V, T, alpha)/self.dpdV(V, T, alpha);
        cpdep = p*dVdT - self.R + \
            alpha/self.b*self.n*(self.n + 1)/T*log((V + self.b)/V) + \
            alpha*(self.n + 1)*dVdT/V/(V + self.b);

        return cpdep

    def cpmcv(self, p, T):
        Z = self.Z(p, T)
        V = Z*self.R*T/p
        Tr = T/self.Tc
        alpha = self.a*pow(Tr, -self.n)
        kappa= -1/(V*self.dpdV(V,T,alpha))
        beta = self.dpdT(V, T, alpha)*kappa
        return V*T*beta**2/kappa

    
    def Hdep(self, p, T):
         Z = self.Z(p, T)
         V = Z*self.R*T/p
         Tr = T/self.Tc
         Hdep = p*V - self.R*T - self.a/self.b*(self.n + 1)*Tr**(-self.n)*log((V + self.b)/V)
         return Hdep
     
     
    def Sdep(self, p, T):
        Tstd = 288.15
        Pstd = 101325
        Z = self.Z(p, T)
        V = Z*self.R*T/p
        Tr = T/self.Tc
        alpha = self.a*pow(Tr, -self.n)
        Vstd = self.R*Tstd/Pstd
        Sdep = -self.R*log(T/Tstd) + self.R*log((V - self.b + self.c)/Vstd) - \
            self.n*alpha/self.b/T*log((V + self.b)/V);
        return Sdep


class Gas:

    def __init__(self, eos, thermo):
        self.thermo_ = thermo
        self.eos_ = eos

    def H(self, p, T):
        return self.thermo_.H0(p,T) + self.eos_.Hdep(p,T)

    def S(self, p, T):
        return self.thermo_.S0(p,T) + self.eos_.Sdep(p,T)

    def cp(self, p, T):
        return self.thermo_.cp0(p,T) + self.eos_.cpdep(p,T)

    def cpmcv(self, p, T):
        return self.eos_.cpmcv(p,T)

    def pT_hs(self, h, S):

        def f(x):
            return [h - self.H(x[0],x[1]), S - self.S(x[0],x[1])]

        sol = root(f, [pTot, Ttot])
        return sol.x
        
    
    def SoundSpeed(self, p, T):
        v = 1/self.eos_.rho(p,T)
        p1 = (1 + 1.e-6)*p
        v1 = 1/self.eos_.rho(p1,T)
        beta_T = -1/v*(v1 - v)/(p1 - p)

        cp = self.cp(p,T)
        cv = cp - self.cpmcv(p,T)
        gamma_pv = cp/cv/beta_T/p
        z = self.eos_.Z(p,T)
        return sqrt(gamma_pv*z*self.eos_.R*T)

    
gas = Gas(AungierRedlichKwong(), JANAF())

H = gas.H(pTot, Ttot)
s = gas.S(pTot, Ttot)

u = 10.
for iter in range(20):
    h = H - 0.5*u**2
    p, T = gas.pT_hs(h, s)
    a = gas.SoundSpeed(p,T)
    Ma = u/a
    u = M2i*a

print(p)




    
