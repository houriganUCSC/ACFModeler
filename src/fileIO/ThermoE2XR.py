import csv
import struct
import time
import re
import numpy as np
import os
from datetime import datetime
#Project imports
from src.records.Session import Session, Mass, Sample


""" BEGIN INSTRUMENT CONSTANTS """
SAMPLE_NUM_SEPARATOR = "[-_]"  # Dash or underscore
MASS_NUMERIC_PRECISION = 0.5
PULSE_THRESHOLD = 4E+6
IGNORE_FARADAY = True

""" BEGIN INSTRUMENT CONSTANTS """
MAG_DAC_BITS = 18  # ELEMENT XR @ UCSC

""" BEGIN CONSTANTS FOR DAT PARSING"""
LEN_DAT_FILE_ID = 16

# DAT FILE MASKS
DAT_TYPE_MASK   = 0xF0000000  # TYPE NIBBLE IN SCAN ITEM
DAT_DATA_MASK   = 0x0FFFFFFF  # DATA BITS IN SCAN ITEM
DATA_FLAG_MASK  = 0x0F000000
DETECT_TYP_MASK = 0x00F00000
DATA_EXP_MASK   = 0x000F0000
DATA_BASE_MASK  = 0x0000FFFF

# DAT FILE OFFSETS
OFFSET_HEADER_START = 0x94
OFFSET_HEADER_LENGTH = 0xAC
OFFSET_TIMESTAMP = 0xB0

# DAT FILE KEYS
KEY_DWELL = 0x30000000  # MASS DWELL TIME
KEY_MAG = 0x20000000  # MAGNET MASS
KEY_MAGF = 0x40000000  # ACTUAL MASS
KEY_INTENSITY = 0x10000000  # INTENSITY
KEY_END_OF_MASS = 0x80000000  # END OF MASS
KEY_END_OF_SCAN = 0xF0000000  # END OF SCAN

KEY_PULSE = 0x00100000
KEY_ANALOG = 0x00000000
KEY_FARADAY = 0x00800000

EXP_SHIFT = 16

"""BEGIN CONSTANTS FOR INF PARSING"""
LEN_INF_FILE_ID = 64

# INF File Bit Masks
MASK_TYP = 0x00000000000000FF
MASK_PTR = 0x00000000FFFFFF00
MASK_TOK = 0x000000FF00000000
MASK_FLG = 0x00000F0000000000
MASK_LEN = 0xFFFFF00000000000

# INF File Data Keys
# Record identifiers for important stored data
KEY_TIMESTAMP = 0x84
KEY_DEADTIME = 0xB8
KEY_RUNS = 0x99
KEY_MASSES = 0x98
KEY_MASS_ID = 0xC2

# INF File Data Offsets
OFFSET_DATETIME = 0x84
OFFSET_FIELDS = 0x108
OFFSET_REGISTRY = 0x114

""" BEGIN SAMPLE RECORD CLASS """


class ThermoDAT(Sample):
    def __init__(self, datPath, session: Session):
        if os.path.exists(datPath):
            super().__init__()
            self.session = session
            self.session.isotopes = None

            # create unique sample primary key
            self.ID = self.session.unique
            self.session.unique += 1
            self.datPath = datPath

    def parseDAT(self):
        """
        High-level function that parses and extracts data from both *.dat and *.inf files.
        """
        _, basename = os.path.split(self.datPath)
        name = basename.strip(".dat")
        self.getDatMetaData()
        if os.path.exists(self.filePaths["INF"]):
            ThermoINF(self).parseINF()
            if self.session.method == None:
                self.session.method = {'isotopes': self.isotopes,
                                       'runs': self.metaData['runs'],
                                       'passes': self.metaData['passes'],
                                       'cycle': self.metaData['cycles'],
                                       'masses':self.metaData["masses"],
                                       'deadTime': self.metaData['deadTime']}

        elif self.session.method is not None:
            self.isotopes = self.session.method['isotopes']
            for key in ['runs', 'passes', 'cycle', 'masses', 'deadTime']:
                self.metaData[key] = self.session.method[key]

        else:
            self.session.isotopes = None
            print("No Isotopes")

        self.getDatScans()
        self.session.startTime = min(self.session.startTime, self.fileTimes["DAT"])
        self.session.isotopes = self.isotopes

    def getDatMetaData(self):
        """
        Extract Method (*.dat), Data (*.dat), and tune (*.tpf) file paths from *.dat file. Dat path is the original
        location (DAT0) of the dat file on the instrument computer.  The "DAT" entry in the filePaths dictionary is the
        current location of the dat file. Data are loaded into class data
        :param dat: file object for binary read (opened for binary read)
        :return:
        """
        self.filePaths["DAT"] = self.datPath
        self.seqDir, basename = os.path.split(self.datPath)
        head, seqName = os.path.split(self.seqDir)
        self.filePaths["SEQ"] = os.path.join(self.seqDir, seqName + ".seq")
        self.filePaths["FIN"] = self.filePaths["SEQ"].replace(".seq", ".FIN")
        self.FIN = seqName +".FIN"
        self.filePaths["INF"] = self.datPath.replace(".dat", ".inf")
        self.filePaths["FIN2"] = self.datPath.replace(".dat", ".FIN2")
        groupRegex = f"{SAMPLE_NUM_SEPARATOR}[0-9a-zA-Z]+.dat"
        self.group = re.split(groupRegex, basename)[0]
        self.name = re.split(".dat", basename)[0]
        numRegex = f"{self.group}{SAMPLE_NUM_SEPARATOR}"
        try:
            self.smpNum = int(re.split(numRegex,self.name)[1])
        except:
            self.smpNum = self.session.unique-1
        self.session.samples[self.name] = self

        with open(self.datPath, mode='rb') as dat:
            # Get file paths (*.dat, *.met, *.tpf) from DAT file
            pathOffset = 356
            paths = ["DAT0","MET","TPF"]
            offsets = [16, 0, 0]
            for i,p in enumerate(paths):
                # Read Dat File Path String
                dat.seek(pathOffset)
                readBytes = 2*(struct.unpack('<1L', dat.read(4))[0])
                pathOffset += 4
                dat.seek(pathOffset)
                path = dat.read(readBytes).decode("utf-16-le")
                self.filePaths[p] = path.rstrip('\x00')
                pathOffset += (readBytes + offsets[i])

            # Get Start time from DAT file
            dat.seek(OFFSET_TIMESTAMP)
            tmp = struct.unpack('<1L', dat.read(4))
            self.fileTimes["DAT"] = time.localtime(tmp[0])
        dat.close()

    def getDatScans(self, infRead = True):
        """
        Extracts raw data block based on entries in DatHdr.
        :param dat: file object for binary read (opened for binary read)
        :return:
        """
        with open(self.datPath, mode='rb') as dat:
            dat.seek(OFFSET_HEADER_START)
            scanStart = struct.unpack('<1L', dat.read(4))[0] + 4
            dat.seek(OFFSET_HEADER_LENGTH)
            length = (struct.unpack('<1L', dat.read(4))[0])
            scanBytes = 4 * length
            dat.seek(scanStart)
            datHdr = struct.unpack('<%dL' % length, dat.read(scanBytes))

            length = datHdr[1] - datHdr[0]
            nVals = int(length / 4)
            rows = len(datHdr)
            cols = nVals
            datScans = np.zeros((rows, cols), dtype=np.uint32)
            # print(f'Len: {length}, rows: {rows}, cols {cols}')
            for i, x in enumerate(datHdr):
                dat.seek(x)
                line = struct.unpack('<%dL' % nVals, dat.read(length))
                try:
                    datScans[i, :] = line
                except:
                    for j, val in enumerate(line):
                        try:
                            datScans[i, j] = np.uint32(val)
                        except:
                            print(f'{i},{j}: {type(val)}')

            ACF = np.array(datScans[:, 12] / 64)
            FCF = np.array(datScans[:, 34] >> 8)
            EDAC = np.array(datScans[:, 31])
            scanTime = np.array((datScans[:, 19] - datScans[0, 18]) / 1000)
            scanTime += time.mktime(self.fileTimes["DAT"])  #ScanTime needs to be absolute for this analysis
            sampleKeys = np.ones_like(scanTime)*self.ID

            # SELECT FIRST ROW TO MASK KEYS FOR POPULATING RAW DATA
            # dwell, magnet mass, actual mass, and truncated mass are pulled from 1st scan (i.e. not a time series)
            # for Intensities the column at the current index is selected (i.e. all scans)

            parseRow = datScans[0, :]
            massIdx = 0
            idx = 46
            channels = 0
            dwell = 0
            # String of isotope ID, used as key for dictionary
            magMasses = np.array([])
            actMasses = np.array([])
            pulse = np.array([])
            analog = np.array([])
            faraday = np.array([])

            for x in parseRow[idx:]:
                key = x & DAT_TYPE_MASK  # DAT_TYPE_MASK   = 0xF0000000
                dataBits = x & DAT_DATA_MASK  # DAT_DATA_MASK   = 0x0FFFFFFF
                if key == KEY_DWELL:
                    dwell = dataBits /1E+6
                    idx += 1
                elif key == KEY_MAG:
                    magMass = dataBits * 1.0 / (2.0 ** (MAG_DAC_BITS))
                    magMasses = np.append(magMasses, magMass)
                    idx += 1
                elif key == KEY_MAGF:
                    magMass = magMasses[-1]
                    actMass = 1 / (float(dataBits)) * magMass * EDAC[0] * 1000
                    actMasses = np.append(actMasses, actMass)
                    idx += 1
                    channels += 1
                elif key == KEY_INTENSITY:
                    detType = dataBits & DETECT_TYP_MASK
                    allScans = datScans[:, idx]
                    iExp = (allScans & DATA_EXP_MASK) >> 16 #EXP_SHIFT
                    iBase = allScans & DATA_BASE_MASK
                    iFlag = (allScans & DATA_FLAG_MASK) > 0
                    values = np.where(iFlag, iBase * float("nan"), iBase << iExp)
                    values = np.expand_dims(values, 1)
                    if detType == KEY_PULSE:
                        if channels == 1:
                            pulse = values
                        else:
                            pulse = np.append(pulse, values, axis=1)
                    elif detType == KEY_ANALOG:
                        if channels == 1:
                            analog = values
                        else:
                            analog = np.append(analog, values, axis=1)
                    elif detType == KEY_FARADAY:
                        if channels == 1:
                            faraday = values
                        else:
                            faraday = np.append(faraday, values, axis=1)
                    idx += 1
                elif key == KEY_END_OF_MASS:
                    if self.isotopes is None:
                        isotope = massIdx
                    else:
                        isotope = self.isotopes[massIdx]
                    # Calculate CHROM Data
                    pCross = self.session.pCross  # Cross-over to analog counts
                    analogCounts = (analog.T * ACF).T
                    # Filter pulse count array NaNs or P-greater-than-threshold values are replaced with corresponding ACF-scaled analog value
                    reported = np.where(pCross > pulse, pulse, analogCounts)
                    if not self.session.ignoreFaraday:
                        reported = np.where(np.isnan(reported), (faraday.T * FCF).T, reported)
                    timeSeries = np.nanmean(reported, axis=1)

                    if isotope not in self.session.masses.keys():
                        self.session.masses[isotope] = Mass(self.session)
                        self.session.masses[isotope].ACF = ACF
                        self.session.masses[isotope].pulse = pulse
                        self.session.masses[isotope].analog = analog
                        self.session.masses[isotope].faraday = faraday
                        self.session.masses[isotope].timeSeries = timeSeries
                        self.session.masses[isotope].chDwell = dwell
                        self.session.masses[isotope].aveMass = np.mean(actMasses)
                        truncMass = round(np.mean(actMasses) / MASS_NUMERIC_PRECISION) * MASS_NUMERIC_PRECISION
                        self.session.masses[isotope].truncMass = truncMass
                        self.session.masses[isotope].channels = channels
                        self.session.masses[isotope].totalDwell = channels * dwell
                        self.session.masses[isotope].magMasses = magMasses
                        self.session.masses[isotope].actMasses = actMasses
                    else:
                        self.session.masses[isotope].ACF = np.append(self.session.masses[isotope].ACF, ACF, axis=0)
                        self.session.masses[isotope].pulse = np.append(self.session.masses[isotope].pulse, pulse, axis=0)
                        self.session.masses[isotope].analog = np.append(self.session.masses[isotope].analog, analog, axis=0)
                        self.session.masses[isotope].faraday = np.append(self.session.masses[isotope].analog, faraday, axis=0)
                        self.session.masses[isotope].timeSeries = np.append(self.session.masses[isotope].timeSeries, timeSeries, axis=0)

                    channels = 0
                    massIdx += 1
                    idx += 1
                    pulse = np.array([[]])
                    analog = np.array([[]])
                    faraday = np.array([[]])
                    magMasses = np.array([])
                    actMasses = np.array([])

                elif key == KEY_END_OF_SCAN:
                    if self.session.FCF.size == 0:
                        self.session.FCF = FCF
                    else:
                        self.session.FCF = np.append(self.session.FCF, FCF)
                    if self.session.EDAC.size == 0:
                        self.session.EDAC = EDAC
                    else:
                        self.session.EDAC = np.append(self.session.EDAC, EDAC)
                    if self.session.scanTime.size == 0:
                        self.session.scanTime = scanTime
                    else:
                        self.session.scanTime = np.append(self.session.scanTime, scanTime)
                    if self.session.sampleKeys.size == 0:
                        self.session.sampleKeys = sampleKeys
                    else:
                        self.session.sampleKeys = np.append(self.session.sampleKeys, sampleKeys)

                    # print(f"End of Scan @{idx}")
                    if not infRead:
                        for key in self.session.masses.keys():
                            intMass = self.session.masses[key].truncMass
                        # TODO: Query
                    self.session.samples[self.name] = self
                    # IoLog.debug(f"{name} DAT file read complete")
        dat.close()



class ThermoINF:
    def __init__(self, datDataObject: ThermoDAT):
        self.dat = datDataObject

    def parseINF(self):
        """
        Partial parsing algorithm to extract, isotope name strings, deadtime, runs/pass, etc. from Thermo
        binary setup information file (*.inf).
        :param infPath:
        :return:
        """

        infPath = self.dat.filePaths["INF"]
        with open(infPath, mode='rb') as inf:
            infHdr = {}

            # Read time Inf file was generated
            inf.seek(OFFSET_DATETIME)
            infSecs = struct.unpack('<1l', inf.read(4))
            self.dat.fileTimes["INF"] = time.localtime(infSecs[0])

            # Read "table of contents"
            inf.seek(OFFSET_FIELDS)
            fields = ord(inf.read(1))
            inf.seek(OFFSET_REGISTRY)
            tmp = inf.read(fields * 8)
            infVals = struct.unpack('<%dQ' % fields, tmp)
            for x in infVals:
                key = (x & MASK_TOK) >> 32
                infHdr[key] = {
                    "type": (x & MASK_TYP),
                    "pointer": (x & MASK_PTR) >> 8,
                    "length": (x & MASK_LEN) >> 44,
                    "flag": (x & MASK_FLG) >> 40
                }

            # Get Deadtime
            tmp = infHdr.get(KEY_DEADTIME)
            inf.seek(tmp['pointer'])
            length = tmp['length']
            entries = length / 2
            dt = struct.unpack('<%dh' % entries, inf.read(length))[0]
            self.dat.metaData["deadTime"] = dt * 1.0E-9

            # Get runs, passes and cycles
            tmp = infHdr.get(KEY_RUNS)
            inf.seek(tmp['pointer'])
            length = tmp['length']
            entries = length / 2
            rp = struct.unpack('<%dh' % entries, inf.read(length))
            self.dat.metaData["runs"] = rp[0]
            self.dat.metaData["passes"] = rp[1]
            self.dat.metaData["cycles"] = rp[0] * rp[1]

            # Get number of masses in run table
            tmp = infHdr.get(KEY_MASSES)
            inf.seek(tmp['pointer'])
            length = tmp['length']
            entries = length / 2
            self.dat.metaData["masses"] = struct.unpack('<%dh' % entries, inf.read(length))[0]

            tmp = infHdr.get(KEY_MASS_ID)
            inf.seek(tmp['pointer'])
            length = tmp['length']
            entries = length / 8
            massIDs = struct.unpack('<%dQ' % entries, inf.read(length))
            isotopes = []
            for x in massIDs:
                inf.seek((x & MASK_PTR) >> 8)
                length = int((x & MASK_LEN) >> 44)
                masses = inf.read(length)
                masses = masses[10:40]
                massID = masses.decode("utf-16-le")
                massID = massID.rstrip('\x00')
                isotopes.append(massID)
            self.dat.isotopes = isotopes
        inf.close()


class ThermoFIN2:

    def __init__(self, session: Session):
        self.session = session
        self.dir = None

    def write(self, fin2Dir, finName, postProcessed: bool = False):
        startTime = time.localtime(datetime.now().timestamp())
        fNames = []
        finPath = os.path.join(fin2Dir, finName)
        for smpName, smpRecord in self.session.samples.items():
            fName = smpName +".FIN2"
            fPath = os.path.join(fin2Dir, fName)
            fNames.append(fName)
            mask = np.where(self.session.sampleKeys == smpRecord.ID)
            cycleTime = self.session.scanTime[mask]
            cycleTime = cycleTime - np.amin(cycleTime)
            cycleTime = cycleTime*10000
            cycleTime = cycleTime.astype(int)/10000
            chromData = np.expand_dims(cycleTime, axis=1)
            chromHdr = ['Time']
            line6 = []
            for massName, massObj in self.session.masses.items():
                if not postProcessed:
                    fixedPrecision = 100 * massObj.timeSeries[mask]
                else:
                    fixedPrecision = 100 * massObj.modeledTimeSeries[mask]
                fixedPrecision = fixedPrecision.astype(int)
                fixedPrecision = fixedPrecision / 100
                fixedPrecision = np.expand_dims(fixedPrecision,axis=1)
                chromData = np.append(chromData, fixedPrecision, axis=1)
                chromHdr.append(massName)
                line6.append("16")
            with open(fPath, 'w', newline='') as fin2:
                    fin2writer = csv.writer(fin2, delimiter=',')
                    fin2writer.writerow(["Finnigan MAT ELEMENT Raw Data"])  # Header 1:  File type description
                    timestamp = time.strftime('%A, %B %d, %Y %H:%M:%S', smpRecord.fileTimes["DAT"])
                    fin2.write(timestamp + '\r\n')  # Header 2:  timestamp
                    fin2writer.writerow([finName])  # Header 3:  FIN file name
                    fin2writer.writerow([f'{len(cycleTime)}'])  # Header 4:  n Scans
                    fin2writer.writerow([0])  # Header 5:  unknown
                    fin2writer.writerow(line6)  # Header 6:  unknown
                    fin2writer.writerow(["CPS"])  # Header 7:  units
                    fin2writer.writerow(chromHdr)  # Header 8:  column names
                    fin2writer.writerows(chromData)  # Data    :  2D data array
                    fin2.close()
        # Write FIN File
        with open(finPath, 'w', newline='') as fin:
            finwriter = csv.writer(fin, delimiter=',')
            finwriter.writerow(["Finnigan MAT ELEMENT"])  # Header 1:  File type description
            timestamp = time.strftime('%A, %B %d, %Y, %H:%M:%S', startTime)
            fin.write(timestamp + '\r\n')  # Header 2:  timestamp
            finwriter.writerow([smpRecord.filePaths["SEQ"]])  # Header 3:  SEQ Path
            finwriter.writerow([smpRecord.filePaths["MET"]])  # Header 4:  MET Path
            finwriter.writerow([smpRecord.filePaths["TPF"]])  # Header 5:  TPF Path
            finwriter.writerow([smpRecord.filePaths["DAT0"]])  # Header 6:  DAT Path
            finwriter.writerow("")  # Header 7:
            finwriter.writerow("")  # Header 8:
            finwriter.writerow([len(self.session.isotopes)])  # Header 9:  n masses
            finwriter.writerow(self.session.isotopes)  # Header 10:  isotope names
            dwells = []
            for isotope, massData in self.session.masses.items():
                dwell = int(1000 * massData.totalDwell)
                dwells = np.append(dwells, f'{dwell:d}')
            finwriter.writerow(dwells)  # Header 11:  dwell times
            for fileName in fNames:
                finwriter.writerow([fileName])  # Header 12:  FIN2 names
            fin.close()