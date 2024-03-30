import numpy as np
import os
import pandas as pd
import glob
import xml.etree.ElementTree as ET
import tqdm

# Find XML files.
def find_xml_files(xml_directory):
    xml_files = list()
    for xml_f in sorted(os.listdir(xml_directory)):
        xml_file_path = os.path.join(xml_directory, xml_f) # Full path for label file
        if os.path.isfile(xml_file_path) and xml_f.lower().endswith('.xml') and not xml_f.lower().startswith('.'):
            xml_files.append(xml_file_path)

    return xml_files


if __name__ == '__main__':
    PathToXML='/labs/collab/12LeadECG/Emory/2022/XML_data'
    PatToOutputCSV = '/labs/collab/12LeadECG/Emory/2022/num_bit_error_files.csv'
    xml_files = find_xml_files(PathToXML)
    num_files = len(xml_files)
    print("--------------------Done reading files--------------------")
    
    incorrect_bits = []
    for i_file, xml_file in enumerate(xml_files):
        fileid = os.path.basename(xml_file)[:-4]
        tree = ET.parse(xml_file)
        root = tree.getroot()
        for parent in root:
            if parent.tag=='Waveform' or parent.tag=='WaveForm':
                for child in parent.iter():
                    for c in child:
                        if c.tag == 'LeadData':
                            for l in c:
                                if l.tag=='LeadHighLimit':
                                    if l.text != '32767':
                                        print('different number of bits')
                                        incorrect_bits.append(xml_file)
                                        
                                    
        print(i_file, " done")
        if i_file == 10:
            break
    data = {'file_names':incorrect_bits}
    df = pd.DataFrame(data)
    df.to_csv(PatToOutputCSV)
    