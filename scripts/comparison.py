import matplotlib.pyplot as plt
import pandas as pd
import re


fileNameWSGG = '../res/qr_T_1000P_1XH2O_20XCO2_20WSGG.res'
fileNameSNB = '../res/qr_T_1000P_1XH2O_20XCO2_20SNB.res'
fileNameWSGG2 = re.split('\.|\/', fileNameWSGG)[-2]
fileNameSNB2 = re.split('\.|\/', fileNameSNB)[-2]

fileNameCassol = '../validation/Cassol2013/p2qr.csv'

qrWSGG = pd.read_csv(fileNameWSGG, delimiter='\t')
#qrSNB = pd.read_csv(fileNameSNB, delimiter='\t')
qrCassol = pd.read_csv(fileNameCassol)

plt.plot(qrWSGG['x [m]'], qrWSGG['qr [W/m2]']/1000, label="WSGG-FlameMaster")
#plt.plot(qrSNB['x [m]'], qrSNB['qr [W/m2]']/1000)
plt.plot(qrCassol['x'], qrCassol['qr'], marker = 's', color = 'k', fillstyle = 'none', linestyle ='none', label="Cassol 2013")
plt.xlabel("X[m]")
plt.ylabel("$\dot{q_r}$ [kW/m$^2$]")
plt.legend()
plt.savefig("../res/comparisonQr.png", dpi = 200, format = "png")
plt.show()

fileNameWSGG = '../res/dqr_T_1000P_1XH2O_20XCO2_20WSGG.res'
fileNameSNB = '../res/dqr_T_1000P_1XH2O_20XCO2_20SNB.res'
fileNameWSGG2 = re.split('\.|\/', fileNameWSGG)[-2]
fileNameSNB2 = re.split('\.|\/', fileNameSNB)[-2]

fileNameCassol = '../validation/Cassol2013/p2dqr.csv'

dqrWSGG = pd.read_csv(fileNameWSGG, delimiter='\t')
#dqrSNB = pd.read_csv(fileNameSNB, delimiter='\t')
dqrCassol = pd.read_csv(fileNameCassol)

plt.plot(dqrWSGG['x [m]'], -dqrWSGG['dqr [W/m3]']/1000, label="WSGG-FlameMaster")
#plt.plot(qrSNB['x [m]'], qrSNB['dqr [W/m3]']/1000)
plt.plot(dqrCassol['x'], dqrCassol['Curve1'], marker = 's', color = 'k', fillstyle = 'none', linestyle ='none', label="Cassol 2013")
plt.xlabel("X[m]")
plt.ylabel("$-\dot{dq_r}$ [kW/m$^3$]")
plt.legend()
plt.savefig("../res/comparisondQr.png", dpi = 200, format = "png")
plt.show()


