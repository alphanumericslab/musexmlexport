import os, os.path, sys, numpy as np
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
from helper_functions import *
warnings.filterwarnings("ignore")
import sys
sys.path.append('/Users/zkoscov/Desktop/work/musexmlexport-master')
from musexmlex import *



class NewIdleParser(XmlElementParser):
    """State for handling the Idle conditio
     modified from orifinal musexml version (added patient demographics for header) """
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


def save_matfile(data,num_channels,record_adcRes,mat_filename):


# Function modified from The native Python waveform-database (WFDB) package. A library of tools for reading, writing, and processing WFDB signals and annotations.
# The development version is hosted at: https://github.com/MIT-LCP/wfdb-python. 
# The development repository is hosted at: https://github.com/MIT-LCP/wfdb-python

# The package is to be expanded with physiological signal-processing tools, and general improvements. Development is made for Python 3.6+ only.

    record_name_out = mat_filename
    record_n_sig = num_channels
    record_p_signal = np.array(data).transpose() 
    record_adc_gain = 1
    record_baseline = num_channels * [0]
    record_adc_zero = num_channels * [0]
    record_fmt =  record_n_sig * str(record_adcRes)
    record_adc_res=record_n_sig * [record_adcRes]




    # Some variables describing the format of the .mat file
    field_version = 256         # 0x0100 or 256
    endian_indicator = b'IM'    # little endian
    master_type = 14            # matrix
    sub1_type = 6               # UINT32
    sub2_type = 5               # INT32
    sub3_type = 1               # INT8
    sub1_class = 6              # double precision array

    sampto = None
    sampfrom=0

    bytes_per_element = 2
    sub4_type = 3       # MAT16
    mut_type = '<i2'    # np.int16
    wfdb_type = '16'    # Align with byte boundary (16)
    offset = 0          # Offset between sample values and the raw
                            # byte/word values as interpreted by Matlab/Octave




    bytes_per_element = 1
    for i in range(record_n_sig):
        if (record_adc_res[i] > 0):
            if (record_adc_res[i] > 16):
                bytes_per_element = 4
            elif (record_adc_res[i] > 8) and (bytes_per_element < 2):
                bytes_per_element = 2
        else:
            # adc_res not specified.. try to guess from format
            if (record_fmt[i] == '24') or (record_fmt[i] == '32'):
                bytes_per_element = 4
            elif (record_fmt[i] != '80') and (bytes_per_element < 2):
                bytes_per_element = 2

    if (bytes_per_element == 1):
        sub4_type = 2       # MAT8
        out_type = '<u1'    # np.uint8
        wfdb_type = '80'    # Offset binary form (80)
        offset = 128        # Offset between sample values and the raw
                            # byte/word values as interpreted by Matlab/Octave
    elif (bytes_per_element == 2):
        sub4_type = 3       # MAT16
        out_type = '<i2'    # np.int16
        wfdb_type = '16'    # Align with byte boundary (16)
        offset = 0          # Offset between sample values and the raw
                            # byte/word values as interpreted by Matlab/Octave
    else:
        sub4_type = 5       # MAT32
        out_type = '<i4'    # np.int32
        wfdb_type = '32'    # Align with byte boundary (32)
        offset = 0          # Offset between sample values and the raw
                            # byte/word values as interpreted by Matlab/Octave






    # Ensure the signal size does not exceed the 2^31 byte limit
    max_length = int((2**31) / bytes_per_element / record_n_sig)
    if sampto is None:
        sampto = record_p_signal.shape[0]
    desired_length = sampto - sampfrom
    # Snip record
    if desired_length > max_length:
        raise Exception("Can't write .mat file: data size exceeds 2GB limit")


    # Bytes of actual data
    bytes_of_data = bytes_per_element * record_n_sig * desired_length
    # This is the remaining number of bytes that don't fit into integer
    # multiple of 8: i.e. if 18 bytes, bytes_remain = 2, from 17 to 18
    bytes_remain = bytes_of_data % 8


    # master_bytes = (8 + 8) + (8 + 8) + (8 + 8) + (8 + bytes_of_data) + padding
    # Must be integer multiple 8
    if bytes_remain == 0:
        master_bytes = bytes_of_data + 56
    else:
        master_bytes = bytes_of_data + 64 - (bytes_remain)




    # Start writing the file
    output_file = record_name_out 
    with open(output_file, 'wb') as f:
        # Descriptive text (124 bytes)
        f.write(struct.pack('<124s', b'MATLAB 4.0'))
        # Version (2 bytes)
        f.write(struct.pack('<H', field_version))
        # Endian indicator (2 bytes)
        f.write(struct.pack('<2s', endian_indicator))

        # Master tag data type (4 bytes)
        f.write(struct.pack('<I', master_type))
        # Master tag number of bytes (4 bytes)
        # Number of bytes of data element
        #     = (8 + 8) + (8 + 8) + (8 + 8) + (8 + bytes_of_data)
        #     = 56 + bytes_of_data
        f.write(struct.pack('<I', master_bytes))

        # Matrix data has 4 subelements (5 if imaginary):
        #     Array flags, dimensions array, array name, real part
        # Each subelement has its own subtag, and subdata

        # Subelement 1: Array flags
        # Subtag 1: data type (4 bytes)
        f.write(struct.pack('<I', sub1_type))
        # Subtag 1: number of bytes (4 bytes)
        f.write(struct.pack('<I', 8))
        # Value class indication the MATLAB data type (8 bytes)
        f.write(struct.pack('<Q', sub1_class))

        # Subelement 2: Rows and columns
        # Subtag 2: data type (4 bytes)
        f.write(struct.pack('<I', sub2_type))
        # Subtag 2: number of bytes (4 bytes)
        f.write(struct.pack('<I', 8))
        # Number of signals (4 bytes)
        f.write(struct.pack('<I', record_n_sig))
        # Number of rows (4 bytes)
        f.write(struct.pack('<I', desired_length))

        # Subelement 3: Array name
        # Subtag 3: data type (4 bytes)
        f.write(struct.pack('<I', sub3_type))
        # Subtag 3: number of bytes (4 bytes)
        f.write(struct.pack('<I', 3))
        # Subtag 3: name of the array (8 bytes)
        f.write(struct.pack('<8s', b'val'))

        # Subelement 4: Signal data
        # Subtag 4: data type (4 bytes)
        f.write(struct.pack('<I', sub4_type))
        # Subtag 4: number of bytes (4 bytes)
        f.write(struct.pack('<I', bytes_of_data))
        # Total size of everything before actual data:
        #     128 byte header
        #     + 8 byte master tag
        #     + 56 byte subelements (48 byte default + 8 byte name)
        #     = 192

        # Copy the selected data into the .mat file
        out_data = record_p_signal * record_adc_gain + record_baseline - record_adc_zero
        # Cast the data to the correct type base on the bytes_per_element
        out_data = np.around(out_data).astype(out_type)
        # out_data should be [r1c1, r1c2, r2c1, r2c2, etc.]
        out_data = out_data.flatten()
        out_fmt = '<%sh' % len(out_data)
        f.write(struct.pack(out_fmt, *out_data))




    # Modify the record file to reflect the new data
        record_fmt = num_channels * [wfdb_type]
        record_byte_offset = num_channels * [192]
        record_baseline = [b - record_adc_zero[i] for i,b in enumerate(record_baseline)]
        record_adc_zero = num_channels * [0]

        return [record_fmt,record_byte_offset,record_baseline,record_adc_zero]




