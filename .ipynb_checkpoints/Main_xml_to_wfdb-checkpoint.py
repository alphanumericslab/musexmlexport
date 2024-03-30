import os, os.path, sys, numpy as np
#import musexmlex_f
import xml.parsers.expat
import re
import codecs
import base64
import array
import struct
import math
from SNOMEDCT_mapper_master import extract_snomed
from SNOMEDCT_mapper_master import extract_MUSE
from scipy.io import savemat
import warnings
import time
warnings.filterwarnings("ignore")
import sys
from musexmlex import *
from helper_functions import *


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
def end_element(name):
    g_Parser.end_element(name)
def char_data(data):
    g_Parser.char_data(data)    




class CustomMuseXmlParser(MuseXmlParser):
    """This class is the parsing context in the object-oriented State pattern.
    It is based on the original MuseXmlParser class, with added function
    for demographics and wfdb (.mat, .hea) export"""
    def __init__(self):
        super().__init__()
        self.__state = IdleParser()
        self.PtID = 0
        self.age = 0
        self.Gender = ""
        self.Race = ""

        self.DataType=""
        self.AcquisitionTime=""
        self.AcquisitionDate=""

        self.DiagnosisText=list()
        self.__state = NewIdleParser()
    
    def setState(self, s):
        self.__state = s
        
    def getState(self):
        return self.__state


       
    def start_element(self, name, attrs):
        """This function trackes the start elements found in the XML file with a
        simple state machine"""
        self.__state.start_element(name, attrs, self)                
            
    def end_element(self, name):
        self.__state.end_element(name, self)                
            
    def char_data(self, data):
        self.__state.char_data(data, self)   

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
        self.DiagnosisText.append(text)
    

    def writeWFDB(self, file_Name,PathToWFDB):
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


        # write bin mat file
        #mdic = {"val": val_leads}
        #savemat(mat_filename, mdic)
        [record_fmt,record_byte_offset,record_baseline,record_adc_zero]=save_matfile(val_leads,num_leads,record_adcRes,mat_filename)

        # Extract snomed ct codes from diagnosis

        snomed_ct=extract_snomed.extract_snomed(str(self.DiagnosisText))
 
        # Extract Muse codes from diagnosis
        muse_codes=extract_MUSE.extract_MUSE(str(self.DiagnosisText))
        
        print(snomed_ct)

        # write header
        fd = open(header_filename, 'wt')
        if fd:
            self.AcquisitionDate = self.AcquisitionDate.replace('-', '/')
            month, day, year = self.AcquisitionDate.split('/')
            # Switch month and day
            self.AcquisitionDate= f"{day}/{month}/{year}"
            fd.write("{} {} {} {} {} {}\n".format(file_Name,len(header),self.sample_Rate,num_samples,self.AcquisitionTime,self.AcquisitionDate))
            for ii,leads_i in enumerate(header):
                adcResNoffset=str(record_adcRes) + "+" + str(record_byte_offset[ii])
                fd.write("{} {} {} {} {} {} {} {} {}\n".format(mat_file,adcResNoffset,units,record_adcRes,record_adc_zero[ii],val_leads[ii][0],eval_check_sum(val_leads[ii]),0,leads_i))
 

            #fd.write("#Age : {}\n".format(self.age))
            #fd.write("#Sex : {}\n".format(self.Gender))
            #fd.write("#Race : {}\n".format(self.Race))
            #fd.write("#PtID : {}\n".format(self.PtID))
            #fd.write("#Dx : {}\n".format(str(snomed_ct)[1:-1]))
            #fd.write("#DxMuseCodes : {}\n".format(str(muse_codes)[1:-1]))
            #fd.write("#DxComment : {}\n".format(",".join(self.DiagnosisText)))


if __name__ == '__main__':
    
    INDEPENDENT_LEADS = ("I", "II", "V1", "V2", "V3", "V4", "V5", "V6")
    # Path To XML files
    PathToXML = '/Users/zkoscov/Desktop/work/EMR/XML_data'
    # Path To WFDB files
    PathToWFDB = '/Users/zkoscov/Desktop/work/EMR/wfdb_data'
    PathToError=''
    
    xml_files = find_xml_files(PathToXML)
    num_files = len(xml_files)
    
    print("--------------------Done reading files--------------------")
    
    start = time.time()
    xmls = [".".join(f.split(".")[:-1]) for f in os.listdir(PathToXML)]
    wfdbs = [".".join(f.split(".")[:-1]) for f in os.listdir(PathToWFDB)]
    
    
    #g_Parser = MuseXmlParser()
    g_Parser = CustomMuseXmlParser()
    
    fileEncoding = None
    p = xml.parsers.expat.ParserCreate()
    p.StartElementHandler = start_element
    p.EndElementHandler = end_element
    p.CharacterDataHandler = char_data
    
    
    for i_file,xml_file_path in enumerate(xml_files):
        print("---------Iteration---------- ",i_file)
        try:
            xml_ready = os.path.basename(xml_file_path)[:-4]
            #if xml_ready in wfdbs:
            #    continue
            # XML parser
            g_Parser = CustomMuseXmlParser()
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
        
            # Write wfdb: .mat and .hea files
            g_Parser.writeWFDB(base_Name,PathToWFDB)
        
            # Write CSV
            #g_Parser.writeCSV(base_Name,PathToCSV)
        
            print("Done.. {} / {}".format(i_file+1,num_files))
            # clean buff
            del p 
    
    
        except Exception as e:
            print("Error in file {} ".format(xml_file_path))
            pass
    
    print("--------------------Done converting files--------------------")
    
    
    
    
    
    
