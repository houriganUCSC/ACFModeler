class RawDataImporter:
    def __init__(self):
        self.type = None

    def import(self, paths):
        pass


class ElementImporter(RawImporter):
    from src.fileIO.ThermoElement import SampleData
    def __init__(self):
        self.type = 'Element'
        super.__init__()
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