import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.optimize import curve_fit
import pandas as pd


df = pd.read_csv('Curva1-SA.csv', names = ["Wet diameter","Supersaturation"])


### definindo a função

def func(D_wet1, D1, k1, A):
     return ((D_wet1**3 - D1**3)/(D_wet1**3 - (D1**3)*(1-k1)))*np.exp(A/D_wet1)
 
xdata = df['Wet diameter']
ydata = df['Supersaturation'] 

parametros_iniciais =[50*10**(-9),0.61,2.165690769263533e-09] # parâmetors iniciais
valores_estimados =[func(x,*parametros_iniciais ) for x in xdata] #valores estimados a partir dos parâmetros iniciais

### Ajuste

popt,pcov = curve_fit(func, xdata, ydata, parametros_iniciais) # (Documentação) Use non-linear least squares to fit a function, f, to data.

print (popt) #parâmetros
print (pcov) #covariância

### Gráficos

intervalo_x = np.exp(np.linspace(np.log(min(xdata)),np.log(max(xdata)),500))
dados_ajustados = [func(x, *popt) for x in intervalo_x]
fig1 = plt.figure(1)
ax=fig1.add_subplot(1,1,1)

### Três gráficos: Gráfico dos dados, dos parâmetros iniciais (chute) e do ajuste

ax.plot(xdata,ydata,linestyle='',marker='x', color='r',label="Dados")
ax.plot(xdata,valores_estimados,linestyle='',marker='x', color='b',label="Chute inicial")
ax.plot(intervalo_x,dados_ajustados,linestyle='-', marker = 'x', color='#900000',label="Ajuste com parâmetros ({0:0.2g},{1:0.2g},{2:0.2g})".format(*popt))

ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2),fancybox=True, shadow=True, ncol=3)
ax.set_ylabel("Supersaturação")
ax.set_xlabel("Diâmetro úmido (um)")

ax.set_xscale('log')

ax.grid()

### Matriz de covariância
#tab= [['{:.2g}'.format(j) for j in i] for i in pcov]
#the_table = plt.table(cellText=tab,
#                  colWidths = [0.2]*3,
#                  loc='upper right', bbox=[0.483, 0.35, 0.5, 0.25] )
#plt.text(250,65,'Covariância:',size=12)
#
#plt.show()
