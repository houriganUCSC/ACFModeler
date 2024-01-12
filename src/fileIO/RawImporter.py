from src.records import Session
from src.fileIO.ThermoE2XR import ThermoDAT

class RawDataImporter:
    def __init__(self, fPath, session: Session):
        self.type = None
        self.filePath = fPath
        self.session = session

    def import(self, paths):
        pass


class ElementImporter(RawDataImporter):
    from src.fileIO.ThermoElement import SampleData
    def __init__(self, fPath, session: Session):
        self.type = 'Element DAT'
        super().__init__(fPath, session)
        self.smpData = SampleData()


    def import(self, paths):
        self.smpData.parseDAT(paths)


    def getStartTime(self):
        self.smpData.

class AgilentImporter(RawImporter):
    def __init__(self):
        self.type = 'Agilent'
        super.__init__()

    def import(self):
        pass

class PerkinElmerImporter(RawImporter):
    def __init__(self):
        self.type = 'Perkin-Elmer'
        super.__init__()

    def import(self):
        pass