import time
import readThermoElement
import sessionRecord
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import numpy as np
import os
import pickle

class Canvas(FigureCanvas):

session = sessionRecord.seshRecord()
dir = "/Volumes/GoogleDrive/My Drive/LA-ICP-MS/External_User_Data/Ingersoll/RVI22-01/Data/RVI_22_01_SEQ1"
#dir = "/Users/jeremyhourigan/Dropbox/Ariolimax/LaserAblation/ThermoDAT/AOS26 DAT"
minStart = np.inf
file_names = [fn for fn in os.listdir(dir) if fn.endswith('dat')]
for i,f in enumerate(file_names):
    file=os.path.join(dir, f)
    print(f)
    s = readThermoElement.sampleData()
    s.parseDAT(file)
    startT = time.mktime(s.fileTimes["DAT"])
    if startT < minStart:
        minStart = startT
    absT = s.scans["scanTime"]+startT
    if i == 0:
        session.createSession(s.isotopes)
    for isotope in s.isotopes:
        s.masses[isotope].filterRaw(5E+6,1000,absT)
        session.appendToRecord(isotope, s.masses[isotope].filtered)

pickle.dump(session, open("../test/save.p", "wb"))

"""fig, axs = plt.subplots(len(s.isotopes),1)
ax = plt.gca()
for i, isotope in enumerate(s.isotopes):
    if session.record[isotope]["vals"] > 20:
        y0 = np.array(session.record[isotope]["pulse"])
        y1 = np.array(session.record[isotope]["analog"])
        x = session.record[isotope]["time"]
        y = y0/y1
        med = np.median(y)
        Q3 = np.percentile(y,75,interpolation='midpoint')
        Q1 = np.percentile(y,25,interpolation='midpoint')
        iqr = Q3-Q1
        inY = np.where(np.abs(y-med)/iqr <= 3, y, np.nan)
        inX = np.where(np.abs(y-med)/iqr <= 3, x, np.nan)
        outY = np.where(np.abs(y-med)/iqr > 3,y,np.nan)
        outX = np.where(np.abs(y-med)/iqr > 3,x,np.nan)
        axs[i].scatter(inX, inY, c = 'blue', s = 2)
        axs[i].scatter(outX, outY, c = 'red', s = 2)
        axs[i].set_yscale('linear')
axs[0].get_shared_x_axes().join(axs[0], *axs[1:])
axs[0].get_shared_y_axes().join(axs[0], *axs[1:])
plt.show()"""

fig = plt.figure()
ax = plt.axes(projection='3d')
if session.record["U238"]["vals"] > 20:
    y0 = np.array(session.record["U238"]["pulse"])
    y1 = np.array(session.record["U238"]["analog"])
    x = session.record["U238"]["time"] - minStart
    y = y0/y1
    med = np.median(y)
    Q3 = np.percentile(y,75,interpolation='midpoint')
    Q1 = np.percentile(y,25,interpolation='midpoint')
    iqr = Q3-Q1
    # print(iqr)
    # print(med)
    # inY = np.where(np.abs(y-med)/iqr <= 5, y, np.nan)
    # inX = np.where(np.abs(y-med)/iqr <= 5, x, np.nan)
    # outY = np.where(np.abs(y-med)/iqr > 5, y, np.nan)
    # outX = np.where(np.abs(y-med)/iqr > 5, x, np.nan)
    #plt.scatter(inX, inY, c = 'blue', s = 2)
    #plt.scatter(outX, outY, c = 'red', s = 2)
    ax.scatter3D(x, y0, y1 , s = 1)
plt.show()


