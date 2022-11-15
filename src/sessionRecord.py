import numpy as np

class seshRecord:

    def __init__(self):
        self.record = {}

    def createSession(self, isotopes):
        for i in isotopes:
            self.record[i] = {
                "pulse":np.array([]),
                "analog":np.array([]),
                "time":np.array([]),
                "vals":0
            }
    def appendToRecord(self,isotope,filtered):
        for key in ["pulse", "analog", "time"]:
            self.record[isotope][key] = np.append(self.record[isotope][key],
                                                  filtered[key])
            self.record[isotope]["vals"] = len(self.record[isotope]["pulse"])
