import csv
import time
from datetime import datetime
import numpy as np


def writeFIN2(dataDir, session = None, postProcessed = False):
    """
    Write chromatogram text file formatted as a ThermoFinnigan (*.FIN2) file

    :param dataDir: Import data directory object
    :param session: SessionFits object
    :param postProcessed: Boolean   True = use post-processed data in Session,
                                    False = use row count equivalents from import
    :return:
    """
    startTime ={}
    for d,files in dataDir.items():
        startTime[d] = time.localtime(datetime.now().timestamp())
        fileList =list(files.keys())
        s = dataDir[d][fileList[0]]['data']
        line6 = []
        header = ["Time"]
        ts = np.zeros ((len(s.scans["scanTime"]),s.metaData["masses"]+1))
        ts[:,0] = s.scans["scanTime"]
        for i, isotope in enumerate(s.isotopes):
            if postProcessed and session is not None:
                intensity = session.sessionFits[isotope].massTimeSeries
            else:
                intensity = s.masses[isotope].massTimeSeries
            fixedPrecision = 10 * intensity
            fixedPrecision = fixedPrecision.astype(int)
            fixedPrecision = fixedPrecision/10
            ts[:,i+1] = fixedPrecision
            header.append(f'{isotope}')
            line6.append("16")
        for file in files.keys():
            s = dataDir[d][file]["data"]
            with open(s.filePaths["FIN2"], 'w', newline='') as fin2:
                fin2writer = csv.writer(fin2, delimiter=',')
                fin2writer.writerow(["Finnigan MAT ELEMENT Raw Data"])
                timestamp = time.strftime('%A, %B %d, %Y, %H:%M:%S', s.fileTimes["DAT"])
                fin2.write(timestamp+'\r\n')
                fin2writer.writerow([s.FIN])
                fin2writer.writerow([f'{len(s.scans["scanTime"])}'])
                fin2writer.writerow([0])
                fin2writer.writerow(line6)
                fin2writer.writerow(["CPS"])
                fin2writer.writerow(header)
                fin2writer.writerows(ts)
                startTime[d] = min(startTime[d], s.fileTimes["DAT"])
        writeFIN(dataDir[d], startTime[d])

def writeFIN(seqFiles, startTime):
    """
    A single FIN file is generated for each "sequence" (directory).
    :param seqFiles: Object contain data from all files in an import directory
    :param startTime: datetime for the begnning of the session.
    :return:
    """
    fileList = list(seqFiles.keys())
    s = seqFiles[fileList[0]]['data']
    with open(s.filePaths["FIN"], 'w', newline='') as fin:
        finwriter = csv.writer(fin, delimiter=',')
        finwriter.writerow(["Finnigan MAT ELEMENT"])
        timestamp = time.strftime('%A, %B %d, %Y, %H:%M:%S', startTime)
        fin.write(timestamp+'\r\n')
        finwriter.writerow([s.filePaths["SEQ"]])
        finwriter.writerow([s.filePaths["MET"]])
        finwriter.writerow([s.filePaths["TPF"]])
        finwriter.writerow([s.filePaths["DAT"]])
        finwriter.writerow("")
        finwriter.writerow("")
        finwriter.writerow([s.metaData["masses"]])
        finwriter.writerow(s.isotopes)
        dwells = []
        for isotope,massData in s.masses.items():
            dwell = int(1000*massData.totaldwell)
            dwells = np.append(dwells, f'{dwell:d}')
        finwriter.writerow(dwells)
        for file in seqFiles.values():
            finwriter.writerow([f'{file["data"].name}.FIN2'])



