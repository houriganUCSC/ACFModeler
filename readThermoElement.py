import struct
import time
import re
import numpy as np
import massRecord
import os


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

class sampleData:

    def __init__(self):
        self.name = ""
        self.group = ""
        self.isotopes = []
        self.filePaths = {}
        self.fileTimes = {}
        self.metaData = {
            "runs": 0,
            "passes": 0,
            "deadTime": 0,
            "masses": 0,
            "cycles": 0
        }
        self.masses = {}
        self.scans = {
            "ACF": [],
            "FCF": [],
            "EDAC": [],
            "scanTime": []
        }

    def constructPaths(self, datPath = None):
        """
        Function that reads in the *.dat file path. With os.path functions creates FIN, FIN2, INF file path
        Regex search determine the sample name (group in Iolite4 parlance). Results are stored in class data.
        :param datPath: Path-like
        :return:
        """
        self.filePaths["DAT"] = datPath
        #basename = os.path.basename(datPath)
        self.seqDir, basename = os.path.split(datPath)
        head, seqName = os.path.split(self.seqDir)
        self.filePaths["SEQ"]=os.path.join(self.seqDir, seqName + ".seq")
        self.filePaths["FIN"] = self.filePaths["SEQ"].replace(".seq", ".FIN")
        self.FIN = seqName +".FIN"
        self.filePaths["INF"] = datPath.replace(".dat", ".inf")
        self.filePaths["FIN2"] = datPath.replace(".dat", ".FIN2")
        groupRegex = f"{SAMPLE_NUM_SEPARATOR}[0-9a-zA-Z]+.dat"
        self.group = re.split(groupRegex, basename)[0]
        self.name = re.split(".dat", basename)[0]

    def getDatFilePaths(self, dat):
        """
        Extracts Method (*.dat), Data (*.dat), and tune (*.tpf) file paths from *.dat file. Dat path is the original
        location (DAT0) of the dat file on the instrument computer.  The "DAT" entry in the filePaths dictionary is the
        current location of the dat file. Data are loaded into class data
        :param dat: file object for binary read (opened for binary read)
        :return:
        """
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

    def getDatTimestamp(self, dat):
        """
        Get timestamp (analysis time) for *.dat file
        :param dat: file object for *dat file binary reads
        :return:
        """
        dat.seek(OFFSET_TIMESTAMP)
        tmp = struct.unpack('<1L', dat.read(4))
        self.fileTimes["DAT"]= time.localtime(tmp[0])

    def getInfData(self, debug=False):
        """
        Reads ThermoFinnigan setup file (*.inf) to extract runs/passes (cycles through run table), pulse counter
        deadtime, number of masses in run table and the strings for each mass (e.g. "Hg202").  Deadtime and mass strings
        are not recorded in the *.dat file.
        :param debug:
        :return:
        """
        inf = thermoINF()
        inf.parseINF(self.filePaths["INF"])
        self.metaData["runs"] = inf.runs
        self.metaData["passes"] = inf.passes
        self.metaData["deadTime"] = inf.deadTime
        self.metaData["masses"] = inf.masses
        self.metaData["cycles"] = inf.runs * inf.passes
        self.fileTimes["INF"] = inf.infTime
        self.isotopes = inf.isotopes
        if debug:
            print("Metadata updated in sample record")

    def createMassRecords(self, debug=False):
        """
        Creates instances of the MassRecord class for each isotope.  Isotope name strings derived from infReader.
        :param debug:
        :return:
        """
        for i,isotope in enumerate(self.isotopes):
            self.masses[isotope] = massRecord.massRecord()
        if debug:
            print(f"{i} masses uploaded into sample record")

    def getDatHdr(self, dat):
        """
        Extacts Dat file lookup table block which has pointers and lengths of cycle data
        :param dat: file object for binary read (opened for binary read)
        :return:
        """
        # Go to position and read pointer to header start
        dat.seek(OFFSET_HEADER_START)
        scanStart = struct.unpack('<1L', dat.read(4))[0] + 4
        dat.seek(OFFSET_HEADER_LENGTH)
        length = (struct.unpack('<1L', dat.read(4))[0])
        scanBytes = 4 * length
        dat.seek(scanStart)
        self.datHdr = struct.unpack('<%dL' % length, dat.read(scanBytes))


    def getDatScans(self, dat):
        """
        Extracts raw data block based on entries in DatHdr.
        :param dat: file object for binary read (opened for binary read)
        :return:
        """

        length = self.datHdr[1] - self.datHdr[0]
        nVals = int(length / 4)
        rows = len(self.datHdr)
        cols = nVals
        datScans = np.zeros((rows, cols), dtype=np.int32)
        for i,x in enumerate(self.datHdr):
            dat.seek(x)
            datScans[i, :] = struct.unpack('<%dL' % nVals, dat.read(length))
        self.scans["ACF"] = datScans[:, 12] / 64
        self.scans["FCF"] = datScans[:, 34] >> 8
        self.scans["FCF"] = datScans[:, 34] >> 8
        self.scans["EDAC"] = datScans[:, 31]
        self.scans["scanTime"] = (datScans[:, 19] - datScans[0, 18]) / 1000
        # SELECT FIRST ROW TO MASK KEYS FOR POPULATING RAW DATA
        # dwell, magnet mass, actual mass, and truncated mass are pulled from 1st scan (i.e. not a time series)
        # for Intensities the column at the current index is selected (i.e. all scans)
        parseRow = datScans[0, :]
        mass = 0
        idx = 46
        channels = 0
        dwell = 0
        ID = self.isotopes[mass]
        for x in parseRow[idx:]:
            key = x & DAT_TYPE_MASK  # DAT_TYPE_MASK   = 0xF0000000
            dataBits = x & DAT_DATA_MASK  # DAT_DATA_MASK   = 0x0FFFFFFF
            if key == KEY_DWELL:
                dwell = dataBits /1E+6
                idx += 1
            elif key == KEY_MAG:
                ID = self.isotopes[mass]
                self.masses[ID].magMass.append((dataBits * 1.0) / 2.0 ** (MAG_DAC_BITS))
                #print(f"{ID} magMass @{idx} is {self.masses[ID].magMass}")
                idx += 1
            elif key == KEY_MAGF:
                magMass = self.masses[ID].magMass[-1]
                edac = self.scans["EDAC"][0]
                actMass = 1 / (float(dataBits)) * magMass * edac * 1000
                self.masses[ID].actMass.append(actMass)
                #print(f"{ID} actMass @{idx} = {self.masses[ID].actMass}")
                idx += 1
                channels += 1
            elif key == KEY_INTENSITY:
                detType = dataBits & DETECT_TYP_MASK
                allScans = datScans[:, idx]
                iExp = (allScans & DATA_EXP_MASK) >> 16 #EXP_SHIFT
                iBase = allScans & DATA_BASE_MASK
                iFlag = (allScans & DATA_FLAG_MASK) > 0
                values = np.where(iFlag, iBase * float("nan"), iBase << iExp)
                if detType == KEY_PULSE:
                    detector = "pulse"
                elif detType == KEY_ANALOG:
                    detector = "analog"
                elif detType == KEY_FARADAY:
                    detector = "faraday"
                else:
                    detector = ""  # THROW EXCEPTION
                # np.append (self.masses[ID].pulse, values, axis = 0)
                dataToAppend = np.expand_dims(values, 1)
                if channels == 1:
                    self.masses[ID].raw[detector] = dataToAppend
                else:
                    self.masses[ID].raw[detector] = np.append(self.masses[ID].raw[detector], dataToAppend, axis=1)
                #print(self.masses[ID].raw["pulse"])
                idx += 1
            elif key == KEY_END_OF_MASS:
                # print(f"{ID} EndofMass @{idx} dwell = {dwell}, {channels} channels")
                self.masses[ID].dwell = dwell
                aveMass = np.mean(self.masses[ID].actMass)
                truncMass = round(aveMass / MASS_NUMERIC_PRECISION) * MASS_NUMERIC_PRECISION
                #print (f"Trunctated mass = {truncMass}")
                self.masses[ID].channels = channels
                self.masses[ID].totaldwell = channels*dwell
                channels = 0
                mass += 1
                idx += 1
            elif key == KEY_END_OF_SCAN:
                #print(f"End of Scan @{idx}")
                name = self.name
                # IoLog.debug(f"{name} DAT file read complete")

    def calculateMassIntensities(self):
        """
        For each isotope in the data file averages (nanmean) channels at each cycle to produce a single intensity
        for each cycle to report as chromatogram text.
        :return:
        """
        ACF = self.scans["ACF"]
        FCF = self.scans["FCF"]
        pMax = PULSE_THRESHOLD
        ignoreFaraday = IGNORE_FARADAY
        for isotope in self.masses:
            self.masses[isotope].processTimeSeries(pMax, ACF, FCF, ignoreFaraday)
            #print(isotope)

    def parseDAT(self, datPath):
        """
        High-level function that parses and extracts data from both *.dat and *.inf files.
        :param datPath: File-like string object
        :return:
        """
        self.constructPaths(datPath)
        with open(datPath, mode='rb') as dat:
            # self.getDatFilePaths(dat)
            self.getDatTimestamp(dat)
            self.getInfData()
            self.getDatFilePaths(dat)
            self.createMassRecords()
            self.getDatHdr(dat)
            self.getDatScans(dat)
            self.calculateMassIntensities()
        dat.close()

class thermoINF:
    def __init__(self):
        self.infTime = None
        self.deadtime = 0
        self.runs = 0
        self.passes = 0
        self.masses = 0
        self.isotopes = []
        self.infPath = ""

    def parseINF(self, infPath):
        """
        Partial parsing algorithm to extract, isotope name strings, deadtime, runs/pass, etc. from Thermo
        binary setup information file (*.inf).
        :param infPath:
        :return:
        """
        self.infPath = infPath
        with open(infPath, mode='rb') as inf:
            infHdr = {}
            inf.seek(OFFSET_DATETIME)
            infSecs = struct.unpack('<1l', inf.read(4))
            self.infTime = time.localtime(infSecs[0])

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
            tmp = infHdr.get(KEY_DEADTIME)
            inf.seek(tmp['pointer'])
            length = tmp['length']
            entries = length / 2
            self.deadTime = struct.unpack('<%dh' % entries, inf.read(length))[0]

            tmp = infHdr.get(KEY_RUNS)
            inf.seek(tmp['pointer'])
            length = tmp['length']
            entries = length / 2
            rp = struct.unpack('<%dh' % entries, inf.read(length))
            self.runs = rp[0],
            self.passes = rp[1]

            tmp = infHdr.get(KEY_MASSES)
            inf.seek(tmp['pointer'])
            length = tmp['length']
            entries = length / 2
            self.masses = struct.unpack('<%dh' % entries, inf.read(length))[0]

            tmp = infHdr.get(KEY_MASS_ID)
            inf.seek(tmp['pointer'])
            length = tmp['length']
            entries = length / 8
            massIDs = struct.unpack('<%dQ' % entries, inf.read(length))
            for x in massIDs:
                inf.seek((x & MASK_PTR) >> 8)
                length = int((x & MASK_LEN) >> 44)
                masses = inf.read(length)
                masses = masses[10:40]
                massID = masses.decode("utf-16-le")
                massID = massID.rstrip('\x00')
                self.isotopes.append(massID)
        inf.close()

