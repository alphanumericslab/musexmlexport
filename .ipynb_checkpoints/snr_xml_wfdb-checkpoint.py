import os, os.path, sys, numpy as np
import xml.parsers.expat
import re
import codecs
import base64
import array
import struct
import math
from scipy.io import loadmat
import warnings
import time
warnings.filterwarnings("ignore")
import pandas as pd
INDEPENDENT_LEADS = ("I", "II", "V1", "V2", "V3", "V4", "V5", "V6")


# Path To XML files
PathToXML='/labs/collab/12LeadECG/Emory/2022/XML_data'
# Path To WFDB files
PathToWFDB='/labs/collab/12LeadECG/Emory/2022/WFDB'
# path to where to save snr between xml and wfdb
path_to_snr = '/labs/collab/12LeadECG/Emory/2022/signal_to_noise_df.csv'

# Classes
###############################################################################
class XmlElementParser:
    """Abstract base class for a XML Parsing State. It contains methods for
    restoring the previous state and for tracking the character data between
    tags."""
    def __init__(self, old_State = None):
        self.__old_State = old_State
        self.__data_Text = ""

    def restoreState(self, context):
        """This method restores the previous state in the XML parser."""
        if self.__old_State:
            context.setState(self.__old_State)

    def clearData(self):
        """This method clears the character data that was collected during parsing"""
        self.__data_Text = ""

    def getData(self):
        """This method returns the character data that was collected during
        parsing and it strips any leading or trailing whitespace"""
        return self.__data_Text.strip()

    def start_element(self, name, attrs, context):
        print ("""abstract method, called at the start of an XML element""")
        sys.exit(0)

    def end_element(self, name, context):
        print ("""abstract method, called at the end of an XML element""")
        sys.exit(0)

    def char_data(self, data, context):
        """This method accumulates any character data"""
        self.__data_Text = self.__data_Text + data



class IdleParser(XmlElementParser):
    """State for handling the Idle condition"""
    def __init__(self):
        XmlElementParser.__init__(self)

    def start_element(self, name, attrs, context):
        if name == WaveformElementParser.Tag:
            context.setState(WaveformElementParser(self))
        if name == PatientDemoElementParser.Tag:
            context.setState(PatientDemoElementParser(self))
        if name == TestDemoElementParser.Tag:
            context.setState(TestDemoElementParser(self))
        if name == DiagnosisElementParser.Tag:
            context.setState(DiagnosisElementParser(self))


    def end_element(self, name, context):
        self.clearData()





class DiagnosisElementParser(XmlElementParser):
    """State for handling the Diagnosis element"""

    Tag = "Diagnosis"

    def __init__(self, old_State):
        XmlElementParser.__init__(self, old_State)

    def start_element(self, name, attrs, context):
        self.clearData()
        if name == DiagnosisStatementElementParser.Tag:
            context.setState(DiagnosisStatementElementParser(self))

    def end_element(self, name, context):
        if name == self.Tag:
            self.restoreState(context)


class DiagnosisStatementElementParser(XmlElementParser):
    """State for handling the DiagnosisStatement element"""

    Tag = "DiagnosisStatement"


    def __init__(self, old_State):
        XmlElementParser.__init__(self, old_State)

    def start_element(self, name, attrs, context):
        self.clearData()
        if name == StmtTextElementParser.Tag:
            context.setState(StmtTextElementParser(self))

    def end_element(self, name, context):
        if name == self.Tag:
            self.restoreState(context)



class StmtTextElementParser(XmlElementParser):
    """State for handling the StmtText element"""

    Tag = "StmtText"

    def __init__(self, old_State):
        XmlElementParser.__init__(self, old_State)

    def start_element(self, name, attrs, context):
        self.clearData()

    def end_element(self, name, context):
        if name == self.Tag:
            self.restoreState(context)
#            context.setStmtText(self.getData())
            context.addStmtText(self.getData())









class TestDemoElementParser(XmlElementParser):
    """State for handling the Demographics of test element"""

    Tag = "TestDemographics"

    def __init__(self, old_State):
        XmlElementParser.__init__(self, old_State)

    def start_element(self, name, attrs, context):
        self.clearData()
        if name == DataTypeElementParser.Tag:
            context.setState(DataTypeElementParser(self))
        elif name == AcquisitionTimeElementParser.Tag:
            context.setState(AcquisitionTimeElementParser(self))
        elif name == AcquisitionDateElementParser.Tag:
            context.setState(AcquisitionDateElementParser(self))

    def end_element(self, name, context):
        if name == self.Tag:
            self.restoreState(context)


class DataTypeElementParser(XmlElementParser):
    """State for handling the Data Type element"""

    Tag = "DataType"

    def __init__(self, old_State):
        XmlElementParser.__init__(self, old_State)

    def start_element(self, name, attrs, context):
        self.clearData()

    def end_element(self, name, context):
        if name == self.Tag:
            self.restoreState(context)
            context.setDataType(self.getData())


class AcquisitionTimeElementParser(XmlElementParser):
    """State for handling the AcquisitionTime element"""

    Tag = "AcquisitionTime"

    def __init__(self, old_State):
        XmlElementParser.__init__(self, old_State)

    def start_element(self, name, attrs, context):
        self.clearData()

    def end_element(self, name, context):
        if name == self.Tag:
            self.restoreState(context)
            context.setAcquisitionTime(self.getData())

class AcquisitionDateElementParser(XmlElementParser):
    """State for handling the AcquisitionDate element"""

    Tag = "AcquisitionDate"

    def __init__(self, old_State):
        XmlElementParser.__init__(self, old_State)

    def start_element(self, name, attrs, context):
        self.clearData()

    def end_element(self, name, context):
        if name == self.Tag:
            self.restoreState(context)
            context.setAcquisitionDate(self.getData())




class PatientDemoElementParser(XmlElementParser):
    """State for handling the Patient Demographics element"""

    Tag = "PatientDemographics"

    def __init__(self, old_State):
        XmlElementParser.__init__(self, old_State)

    def start_element(self, name, attrs, context):
        self.clearData()
        if name == PtIDElementParser.Tag:
            context.setState(PtIDElementParser(self))
        elif name == AgeElementParser.Tag:
            context.setState(AgeElementParser(self))
        elif name == GenderElementParser.Tag:
            context.setState(GenderElementParser(self))
        elif name == RaceElementParser.Tag:
            context.setState(RaceElementParser(self))

    def end_element(self, name, context):
        if name == self.Tag:
            self.restoreState(context)



class PtIDElementParser(XmlElementParser):
    """State for handling the Patient ID element"""

    Tag = "PatientID"

    def __init__(self, old_State):
        XmlElementParser.__init__(self, old_State)

    def start_element(self, name, attrs, context):
        self.clearData()

    def end_element(self, name, context):
        if name == self.Tag:
            self.restoreState(context)
            context.setPtID(self.getData())



class AgeElementParser(XmlElementParser):
    """State for handling the Age element"""

    Tag = "PatientAge"

    def __init__(self, old_State):
        XmlElementParser.__init__(self, old_State)

    def start_element(self, name, attrs, context):
        self.clearData()

    def end_element(self, name, context):
        if name == self.Tag:
            self.restoreState(context)
            context.setAge(self.getData())
#            print ("Patient Age is %s sps..." % (self.getData()))



class GenderElementParser(XmlElementParser):
    """State for handling the Gender element"""

    Tag = "Gender"

    def __init__(self, old_State):
        XmlElementParser.__init__(self, old_State)

    def start_element(self, name, attrs, context):
        self.clearData()

    def end_element(self, name, context):
        if name == self.Tag:
            self.restoreState(context)
            context.setGender(self.getData())


class RaceElementParser(XmlElementParser):
    """State for handling the Gender element"""

    Tag = "Race"

    def __init__(self, old_State):
        XmlElementParser.__init__(self, old_State)

    def start_element(self, name, attrs, context):
        self.clearData()

    def end_element(self, name, context):
        if name == self.Tag:
            self.restoreState(context)
            context.setRace(self.getData())






class WaveformElementParser(XmlElementParser):
    """State for handling the Waveform element"""

    Tag = "Waveform"

    def __init__(self, old_State):
        XmlElementParser.__init__(self, old_State)

    def start_element(self, name, attrs, context):
        self.clearData()
        if name == WaveformTypeElementParser.Tag:
            context.setState(WaveformTypeElementParser(self))
        elif name == LeadDataElementParser.Tag:
            context.setState(LeadDataElementParser(self))
        elif name == SampleBaseElementParser.Tag:
            context.setState(SampleBaseElementParser(self))

    def end_element(self, name, context):
        if name == self.Tag:
            self.restoreState(context)

class SampleBaseElementParser(XmlElementParser):
    """State for handling the SampleBase element"""

    Tag = "SampleBase"

    def __init__(self, old_State):
        XmlElementParser.__init__(self, old_State)

    def start_element(self, name, attrs, context):
        self.clearData()

    def end_element(self, name, context):
        if name == self.Tag:
            self.restoreState(context)
            if context.found_Rhythm:
                context.setSampleBase(self.getData())
#                print ("Sampling rate for rhythm is %s sps..." % (context.sample_Rate))

class LeadUnitsPerBitElementParser(XmlElementParser):
    """State for handling the LeadAmplitudeUnitsPerBit element"""

    Tag = "LeadAmplitudeUnitsPerBit"

    def __init__(self, old_State):
        XmlElementParser.__init__(self, old_State)

    def start_element(self, name, attrs, context):
        self.clearData()

    def end_element(self, name, context):
        if name == self.Tag:
            self.restoreState(context)
            context.setAdu(float(self.getData().strip()))

class LeadUnitsElementParser(XmlElementParser):
    """State for handling the LeadAmplitudeUnits element"""

    Tag = "LeadAmplitudeUnits"

    def __init__(self, old_State):
        XmlElementParser.__init__(self, old_State)

    def start_element(self, name, attrs, context):
        self.clearData()

    def end_element(self, name, context):
        if name == self.Tag:
            self.restoreState(context)
            context.setUnits(self.getData().strip())


class WaveformTypeElementParser(XmlElementParser):
    """State for handling the WaveformType element"""

    Tag = "WaveformType"

    def __init__(self, old_State):
        XmlElementParser.__init__(self, old_State)

    def start_element(self, name, attrs, context):
        self.clearData()

    def end_element(self, name, context):
        if name == self.Tag:
            self.restoreState(context)
            if self.getData().find("Rhythm") >= 0:
                context.setRhythmFound(1)
#                print ("ECG %s object found." % self.getData())
            else:
                context.setRhythmFound(0)

class LeadDataElementParser(XmlElementParser):
    """State for handling the LeadData element"""

    Tag = "LeadData"

    def __init__(self, old_State):
        XmlElementParser.__init__(self, old_State)

    def start_element(self, name, attrs, context):
        self.clearData()
        if name == LeadIdElementParser.Tag:
            context.setState(LeadIdElementParser(self))
        if name == WaveformDataElementParser.Tag:
            context.setState(WaveformDataElementParser(self))
        if name == LeadUnitsPerBitElementParser.Tag:
            context.setState(LeadUnitsPerBitElementParser(self))
        if name == LeadUnitsElementParser.Tag:
            context.setState(LeadUnitsElementParser(self))

    def end_element(self, name, context):
        if name == self.Tag:
            self.restoreState(context)

class LeadIdElementParser(XmlElementParser):
    """State for handling the LeadID element"""

    Tag = "LeadID"

    def __init__(self, old_State):
        XmlElementParser.__init__(self, old_State)

    def start_element(self, name, attrs, context):
        self.clearData()

    def end_element(self, name, context):
        if name == self.Tag:
            self.restoreState(context)
            if context.found_Rhythm:
#                sys.stdout.write("   Lead %2s found..." % self.getData())
                context.addLeadId(self.getData())

class WaveformDataElementParser(XmlElementParser):
    """State for handling the WaveformData element"""

    Tag = "WaveFormData"

    def __init__(self, old_State):
        XmlElementParser.__init__(self, old_State)

    def start_element(self, name, attrs, context):
        self.clearData()

    def end_element(self, name, context):
        if name == self.Tag:
            self.restoreState(context)
            if context.found_Rhythm:
#                print ("   Adding data for lead %2s." % context.lead_Id)
                context.addWaveformData(self.getData())






class MuseXmlParser:
    """This class is the parsing context in the object-oriented State pattern."""
    def __init__(self):
        self.ecg_Data = dict()
        self.ecg_Leads = list()
        self.__state = IdleParser()
        self.found_Rhythm = 0
        self.sample_Rate = 0
        self.adu_Gain = 1
        self.units = ""

        self.PtID = 0
        self.age = 0
        self.Gender = ""
        self.Race = ""

        self.DataType=""
        self.AcquisitionTime=""
        self.AcquisitionDate=""

        self.DiagnosisText=list()


    def setState(self, s):
        self.__state = s

    def getState(self):
        return self.__state

    def setSampleBase(self, text):
        if self.found_Rhythm:
            self.sample_Rate = int(text)

    def setAdu(self, gain):
        self.adu_Gain = gain

    def setUnits(self, units):
        self.units = units




    def setPtID(self, PtID):
        self.PtID = PtID
    def setAge(self, age):
        self.age = age
    def setGender(self, Gender):
        self.Gender = Gender
    def setRace(self, Race):
        self.Race = Race


    def setDataType(self, DataType):
        self.DataType = DataType
    def setAcquisitionTime(self, AcquisitionTime):
        self.AcquisitionTime = AcquisitionTime
    def setAcquisitionDate(self, AcquisitionDate):
        self.AcquisitionDate = AcquisitionDate

     
    
    def addStmtText(self,text):
#        self.setStmtText.append(text)
        self.DiagnosisText.append(text)


    def setRhythmFound(self, v):
        self.found_Rhythm = v

    def addLeadId(self, id):
        self.lead_Id = id

    def addWaveformData(self, text):
        self.ecg_Data[self.lead_Id] = base64.b64decode(text)
        self.ecg_Leads.append(self.lead_Id)




    def start_element(self, name, attrs):
        """This function trackes the start elements found in the XML file with a
        simple state machine"""
        self.__state.start_element(name, attrs, self)

    def end_element(self, name):
        self.__state.end_element(name, self)

    def char_data(self, data):
        self.__state.char_data(data, self)

    def makeZcg(self):
        """This function converts the data read from the XML file into a ZCG buffer
        suitable for storage in binary format."""
#        print(self.ecg_Leads)
        # All of the leads should have the same number of samples
        n = len(self.ecg_Data[self.ecg_Leads[0]])
        # We have 2 bytes per ECG sample, so make our buffer size n * DATAMUX
        self.zcg = array.array('d')
        # Verify that all of the independent leads are accounted for...
        for lead in INDEPENDENT_LEADS:
            if lead not in self.ecg_Leads:
                print ("Error! The XML file is missing data for lead ", lead)
                sys.exit(-1)

        # Append the data into our huge ZCG buffer in the correct order
        for t in range(0,n,2):
            for lead in self.ecg_Leads:
                sample = struct.unpack("h", self.ecg_Data[lead][t:t+2])
                self.zcg.append(sample[0])
       



    def read_xml_leads(self, file_Name,PathToWFDB):


        """This function writes the ZCG buffer to a header file"""
        std_Leads = set(INDEPENDENT_LEADS)
        header = ("I", "II", "III", "aVR", "aVL", "aVF", "V1", "V2", "V3", "V4", "V5", "V6")
        num_leads = len(header)

        extra_Leads = std_Leads.symmetric_difference(set(self.ecg_Leads))
        #print ("EXTRA LEADS: ", extra_Leads)


        # filenames
        
        header_file=file_Name + ".hea"
        mat_file = file_Name + ".mat"

        header_filename = os.path.join(PathToWFDB, header_file)
        mat_filename = os.path.join(PathToWFDB, mat_file)


        # Information needed to create the wfdb bin files
        record_byte_offset=0
        record_adcGain=1
        record_adcUnits=self.units
        record_baseline=0
        record_adcRes=16
        record_adcZero=0
        record_blockSize=0


        # units
        if self.units=="MICROVOLTS":
            units="1000/mv"
        else:
            units="1/"+self.units


        samples = dict()

        val_leads=[[] for i in range(12)]


        for i in range(0, len(self.zcg), len(self.ecg_Leads)):
            # The values in the ZCG buffer are stored in the same order
            # as the ecg leads are themselves...
            k = 0
            for lead in self.ecg_Leads:
                samples[lead] = self.zcg[i + k]
                k = k + 1
                # Output each sample, calculated and uncalcuated

            val_leads[0].append(int(samples["I"] * self.adu_Gain))
            val_leads[1].append(int(samples["II"] * self.adu_Gain))
            # II - I = III
            val_leads[2].append(int((samples["II"] - samples["I"]) * self.adu_Gain))
            # aVR = -(I + II)/2
            val_leads[3].append(int((-(samples["I"] + samples["II"])/2) * self.adu_Gain))
            # aVL = I - II/2
            val_leads[4].append(int((samples["I"] - samples["II"]/2) * self.adu_Gain))
            # aVF = II - I/2
            val_leads[5].append(int((samples["II"] - samples["I"]/2) * self.adu_Gain))
            # output the precordial leads
            val_leads[6].append(int(samples["V1"] * self.adu_Gain))
            val_leads[7].append(int(samples["V2"] * self.adu_Gain))
            val_leads[8].append(int(samples["V3"] * self.adu_Gain))
            val_leads[9].append(int(samples["V4"] * self.adu_Gain))
            val_leads[10].append(int(samples["V5"] * self.adu_Gain))
            val_leads[11].append(int(samples["V6"] * self.adu_Gain))


        # number of samples
        num_samples = len(val_leads[0])
        return(val_leads)




# Find XML files.
def find_xml_files(xml_directory):
    xml_files = list()
    for xml_f in sorted(os.listdir(xml_directory)):
        xml_file_path = os.path.join(xml_directory, xml_f) # Full path for label file
        if os.path.isfile(xml_file_path) and xml_f.lower().endswith('.xml') and not xml_f.lower().startswith('.'):
            xml_files.append(xml_file_path)

    return xml_files


# 3 handler functions for the expat parser
def start_element(name, attrs):
    g_Parser.start_element(name, attrs)
#    print ('Start element:', name, attrs)
def end_element(name):
    g_Parser.end_element(name)
#    print ('End element:', name)
def char_data(data):
#    print(data)
    g_Parser.char_data(data)


def read_encoding(xmlFileName):
    fid = open(xmlFileName, 'rb')
    if fid:
        pattern = "\\<\\?xml\\s+.*encoding=\"([\w-]+)\"\\?\\>"
        # print("pattern:", pattern)

        for bytes_data in fid.readlines():
            line = "".join(map(chr, bytes_data))
            #print("line:", line)

            result = re.match(pattern, line)
            if result:
#                print ("encoding is ", result.group(1))
                return result.group(1)

    return ""


def eval_check_sum(y):

    bit = 16
    check_sum = np.nansum(y)
    M = check_sum/pow(2,(bit-1))

    if M<0:
        check_sum = check_sum % pow(-2,(bit-1))

        if ((not check_sum) & (abs(M)<1)):
            check_sum = pow(-2,(bit-1))
        elif math.ceil(M) % 2:
            check_sum = pow(2,(bit-1)) + check_sum
    else:
        check_sum = check_sum % pow(2,(bit-1))
        if math.floor(M) % 2:
            check_sum = pow(-2,(bit-1))+check_sum
    return check_sum



PathToError=''

g_Parser = MuseXmlParser()

fileEncoding = None
p = xml.parsers.expat.ParserCreate()
p.StartElementHandler = start_element
p.EndElementHandler = end_element
p.CharacterDataHandler = char_data


xml_files = find_xml_files(PathToXML)
num_files = len(xml_files)

print("--------------------Done reading files--------------------")

start = time.time()
xmls = [".".join(f.split(".")[:-1]) for f in os.listdir(PathToXML)]
wfdbs = [".".join(f.split(".")[:-1]) for f in os.listdir(PathToWFDB)]


### Calculate SNR and save to the .csv
snr_vector = []
wfdb_vector = []
xml_vector = []
for i_file,xml_file_path in enumerate(xml_files):
    print("---------Iteration---------- ",i_file)
    try:
        xml_ready = os.path.basename(xml_file_path)[:-4]
        g_Parser = MuseXmlParser()
        p = xml.parsers.expat.ParserCreate()
        p.StartElementHandler = start_element
        p.EndElementHandler = end_element
        p.CharacterDataHandler = char_data
    
        # Automatically determine the encoding of XML file
        if not fileEncoding:
            fileEncoding = read_encoding(xml_file_path)
            if len(fileEncoding) == 0:
                print("ERROR: Cannot determine file encoding from XML file!")
                sys.exit(1)
    
        # Read the XML file and parse it
        f = codecs.open(xml_file_path, mode='r', encoding=fileEncoding)
        p.Parse(f.read())
        g_Parser.makeZcg()
    
        base_Name = os.path.splitext(os.path.basename(xml_file_path))[0]
    
        # read leads from xml
        leads = g_Parser.read_xml_leads(base_Name,PathToWFDB)
        leads = np.array(leads)
        
        # read mat
        wfdb_name = os.path.join(PathToWFDB, base_Name) + '.mat'
        xml_vector.append(xml_file_path)
        wfdb_vector.append(wfdb_name)
        if not os.path.exists(wfdb_name): # if wfdb doesn't exist
            snr = 0
            snr_vector.append(snr)
            continue
    
        mat_file = loadmat(wfdb_name)
        rec = mat_file['val']
    
    
        # take only lead I
        xml_leadI = leads[0,:]
        mat_leadI = rec[0,:]
    
        # Calculate the SNR 
        if len(xml_leadI) != len(mat_leadI): # if the length isn't the same snr 0
            snr = 0
            snr_vector.append(snr)
            continue
        
        dif = xml_leadI - mat_leadI # diference of signals, exptected SNR inf
        signal_power = np.mean(np.abs(xml_leadI)**2)
        noise_power = np.mean(np.abs(dif)**2)
        snr = 10 * np.log10(signal_power / noise_power)
        snr_vector.append(snr)
       
    
        print("Done.. {} / {}".format(i_file+1,num_files))
        # clean buff
        del p 


    except Exception as e:
        print("Error in file {} ".format(xml_file_path))
        pass
    #if i_file==20:
    #    break

print("--------------------Done converting files--------------------")

data = {'xml_path': xml_vector,
        'wfdb_path': wfdb_vector,
        'snr': snr_vector}
df = pd.DataFrame(data)
df.to_csv(path_to_snr)
