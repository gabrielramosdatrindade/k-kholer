import numpy as np
import matplotlib.pylab as plt
from numpy import *
from matplotlib.pylab import *
from prettytable import PrettyTable
import pandas as pd

#------------------------------------------------------------------------------
# Gráfico supersaturação x diâmetro úmido

D_wet1 = np.linspace(0.5*10**(-7),10**(-5),5000)
D_wet2 = np.linspace(1*10**(-7),10**(-5),5000)
D_wet3 = np.linspace(5*10**(-7),10**(-5),5000)


D1,D2,D3 = 0.5*10**(-7), 1*10**(-7), 5*10**(-7)  # diâmetro seco em m
k1,k2 = 0.61, 1.28 # parâmetro de higroscopicidade, k(NaCl) = 1.28 e k((NH4)2SO4) = 0.61
Ten_s = 0.073 # tensão superficial em J/m^2
T = 293 # temperatura em K
rho_w = 997.07 #densidade da água em kg/m^3
M_w = 18.01528 * 10**(-3) # massa molar da água em kg/mol
R = 8.314462 # constante universal dos gases perfeitos em J/(K*mol)
A = (4*Ten_s*M_w)/(R*T*rho_w)

def f():
    s1 = ((D_wet1**3 - D1**3)/(D_wet1**3 - (D1**3)*(1-k1)))*np.exp(A/D_wet1) #supersaturação p/ k1, D1
    s2 = ((D_wet1**3 - D1**3)/(D_wet1**3 - (D1**3)*(1-k2)))*np.exp(A/D_wet1) #supersaturação p/ k2, D1
    s3 = ((D_wet2**3 - D2**3)/(D_wet2**3 - (D2**3)*(1-k1)))*np.exp(A/D_wet2) #supersaturação p/ k1, D2
    s4 = ((D_wet2**3 - D2**3)/(D_wet2**3 - (D2**3)*(1-k2)))*np.exp(A/D_wet2) #supersaturação p/ k2, D2
    s5 = ((D_wet3**3 - D3**3)/(D_wet3**3 - (D3**3)*(1-k1)))*np.exp(A/D_wet3) #supersaturação p/ k1, D3
    s6 = ((D_wet3**3 - D3**3)/(D_wet3**3 - (D3**3)*(1-k2)))*np.exp(A/D_wet3) #supersaturação p/ k2, D3
        
    return (s1-1)*100, (s2-1)*100, (s3-1)*100, (s4-1)*100, (s5-1)*100, (s6-1)*100

S1,S2,S3,S4,S5,S6 = f()

df = pd.read_csv('sulfatodeamonio.csv', names = ["Wet diameter","Supersaturation"])

fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111, facecolor='w', axisbelow=True)
ax.plot(D_wet1/10**(-6), S1, 'k', alpha=1, lw=2, label ='D = 50 nm kappa = 0,61')
ax.plot(D_wet2/10**(-6), S3, 'b', alpha=1, lw=2, label ='D = 100 nm kappa = 0,61')
ax.plot(D_wet3/10**(-6), S5, 'r', alpha=1, lw=2, label ='D = 500 nm kappa = 0,61')
ax.scatter(df["Wet diameter"], df["Supersaturation"], label ='S&P')
ax.set_xlabel('Diâmetro úmido (nm)')
ax.set_ylabel('Supersaturação (%)')
ax.set_ylim(-0.3,0.5)
ax.set_xscale('log')
ax.set_title('Sulfato de amônio', fontsize = 15)


ax.grid(b=True, which='major', c='k', lw=1, ls='-')
ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2),fancybox=True, shadow=True, ncol=3)
# -----------------------------------------------------------------------------
# Tabela supersaturação crítica x valores de D,k
x = PrettyTable()

x.field_names = [ 'D e k','SC']

x.add_row(['D = 50 nm, kappa = 0,61', np.amax(S1)])
x.add_row(['D = 50 nm, kappa = 1,28', np.amax(S2)])
x.add_row(['D = 100 nm, kappa = 0,61', np.amax(S3)])
x.add_row(['D = 100 nm, kappa = 1,28', np.amax(S4)])
x.add_row(['D = 500 nm, kappa = 0,61', np.amax(S5)])
x.add_row(['D = 500 nm, kappa = 1,28', np.amax(S6)])

print(x)
# -----------------------------------------------------------------------------
# Tabela para a aproximação da máxima saturação

C1 = (16/3)*(1/(3*k1))**(1/2) # k1 = 0,61
C2 = (16/3)*(1/(3*k2))**(1/2) # k2 = 1,28
Z = ((Ten_s*M_w)/(R*T*rho_w))**(3/2)
D_c = (D1, D2, D3)
S_c1 = []
S_c2 = []

for i in D_c:
    S_c1.append((C1*Z*(1/i)**(3/2))*100)
    S_c2.append((C2*Z*(1/i)**(3/2))*100)

y = PrettyTable()

y.field_names = [ 'D e k','SC (aprox)']

y.add_row(['D = 50 nm, kappa = 0,61', S_c1[0]])
y.add_row(['D = 50 nm, kappa = 1,28', S_c2[0]])
y.add_row(['D = 200 nm, kappa = 0,61', S_c1[1]])
y.add_row(['D = 200 nm, kappa = 1,28', S_c2[1]])
y.add_row(['D = 350 nm, kappa = 0,61', S_c1[2]])
y.add_row(['D = 350 nm, kappa = 1,28', S_c2[2]])

print(y)
#------------------------------------------------------------------------------  







