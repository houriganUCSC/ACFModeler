class CrossCalHeaderCodes():
    def __init__(self):
        self.codes = {
            'Mass': {'hexCode':'Mass', "numFormat":'', "width": 75},
            'use ACF': {'hexCode': b"\x75\x73\x65\x20\xCE\xB1", "numFormat":'', "width": 50},
            'use Tau': {'hexCode':b"\x75\x73\x65\x20\xCF\x84", "numFormat":'', "width": 50},
            'ACF Type': {'hexCode':b"\xCF\x84\x20\x53\x6F\x75\x72\x63\x65", "numFormat":'', "width": 100},
            'Tau Type': {'hexCode':b"\xCE\xB1\x20\x53\x6F\x75\x72\x63\x65", "numFormat":'', "width": 100},
            'a1': {'hexCode':b"\x61\xE2\x82\x81", "numFormat":'', "width": 100},
            'a2': {'hexCode':b"\x61\xE2\x82\x82", "numFormat":'', "width": 100},
            'tau': {'hexCode':b"\xCF\x84", "numFormat":'', "width": 100},
            'int_a1': {'hexCode':b"\x69\x6E\x74\x20\x61\xE2\x82\x81", "numFormat":'', "width": 100},
            'int_a2': {'hexCode':b"\x69\x6E\x74\x20\x61\xE2\x82\x82", "numFormat":'', "width": 100},
            'int_tau': {'hexCode':b"\x69\x6E\x74\x20\xCF\x84", "numFormat":'', "width": 100},
            'se_a1': {'hexCode':b"\xCF\x83\x28\x61\xE2\x81\x81\x29", "numFormat":'', "width": 50},
            'se_a2': {'hexCode':b"\xCF\x83\x28\x61\xE2\x81\x82\x29", "numFormat":'', "width": 50},
            'se_tau': {'hexCode':b"\xCF\x83\x28\xCF\x84\x29", "numFormat":'', "width": 50},
            'rSqr': {'hexCode':b"\x72\xC2\xB2", "numFormat":'', "width": 50},
            "redChi2": {'hexCode':b"\xCF\x87\xE1\xB5\xA3\xC2\xB2", "numFormat":'', "width": 50}
            }

    def getHeaderLabel(self, key:str):
        hdr = self.codes[key]
        if isinstance(hdr, bytes):
            label = hdr.decode('UTF-8')
        else:
            label = hdr
        return label

    def getDefaultWidth(self, key:str):
        return self.codes[key]['width']

    def formatValue(self, key:str, val: float):
        valFormat =  self.codes[key]['numFormat']
        valText = f'{val:{valFormat}}'
        if 'G' in valFormat:
            valText = valText.replace("E-0", "E-").replace("E+0", "E+")
        return valText


class ModelHeaderCodes():
    def __init__(self):
        self.codes = {
            'ACF (a1)': b"\x41\x43\x46\x20\x28\x61\xE2\x82\x81\x29",
            'Drift (a2)': b"\x44\x72\x69\x66\x74\x20\x28\x61\xE2\x82\x82\x29",
            'Tau': b"\xCF\x84",
            }

    def getHeaderLabel(self, key:str):
        hdr = self.codes[key]['hexCode']
        if isinstance(hdr, bytes):
            label = hdr.decode('UTF-8')
        else:
            label = hdr
        return label


class IsotopeFitHeaderCodes():
    def __init__(self):
        self.codes = {
            'Mass': {'hexCode': 'Mass', "numFormat": '0.2f', "width": 75},
            'n Used': {'hexCode': 'n Used', "numFormat": 'd', "width": 75},
            'n': {'hexCode': 'n', "numFormat": 'd', "width": 75},
            'a1': {'hexCode': b"\x61\xE2\x82\x81", "numFormat": '.4G', "width": 75},
            'se_a1': {'hexCode': b"\xCF\x83\x28\x61\xE2\x81\x81\x29", "numFormat": '.4G', "width": 75},
            'a2': {'hexCode': b"\x61\xE2\x82\x82", "numFormat": '.4G', "width": 75},
            'se_a2': {'hexCode': b"\xCF\x83\x28\x61\xE2\x81\x82\x29", "numFormat": '.4G', "width": 75},
            'tau': {'hexCode': b"\xCF\x84", "numFormat": '.4G', "width": 75},
            'se_tau': {'hexCode': b"\xCF\x83\x28\xCF\x84\x29", "numFormat": '.4G', "width": 75},
            'rSqr': {'hexCode': b"\x72\xC2\xB2", "numFormat": '0.3f', "width": 50},
            "redChi2": {'hexCode': b"\xCF\x87\xE1\xB5\xA3\xC2\xB2", "numFormat": '0.2f', "width": 50},
            "ACF": {'hexCode': "ACF", "numFormat": '0.1f', "width": 75},
            "Drift": {'hexCode': "Drift", "numFormat": '0.1%', "width": 75}
            }


    def getHeaderLabel(self, key:str):
        hdr = self.codes[key]['hexCode']
        if isinstance(hdr, bytes):
            label = hdr.decode('UTF-8')
        else:
            label = hdr
        return label

    def getDefaultWidth(self, key:str):
        return self.codes[key]["width"]

    def formatValue(self, key:str, val: float):
        valFormat =  self.codes[key]['numFormat']
        valText = f'{val:{valFormat}}'
        if 'G' in valFormat:
            valText = valText.replace("E-0", "E-").replace("E+0", "E+")
        return valText