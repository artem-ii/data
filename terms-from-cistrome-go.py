import numpy as np
import pandas as pd
import os
import sys
import fileinput

if sys.version_info[0] < 3:
    from StringIO import StringIO
else:
    from io import StringIO

files = input("List of files: ")
files_list = files.split()

flag = False
filedata = str()
for f in files_list:
    filedata = str()
    print("Processing" + f)
    with open(f, "r") as file:
        for line in file:
            if "# Term" in line:
                flag = True
            if flag == True:
                filedata += line
        filedata = StringIO(filedata)
        df = pd.read_csv(filedata, sep = "\t")
        print(df.head())
        sig = df['FDR'] <= 0.05
        sigterms = df[sig]
        sigterms = sigterms.iloc[:,0]
        sigterms_list = list(sigterms)
        print('There are ' + str(len(sigterms)) +
         ' significantly changed terms')
        terms = str()
        for i in sigterms_list:
            splt = i.split("(G")
            if len(splt) > 1:
                term = splt[1]
                term = "G" + term[:-1] + "\n"
            else:
                term = i + "\n"
            print(term)
            terms += term
    name = str(f).split("/")[-1]
    name = str(name)
    with open(str("terms-" + str(name)), 'w') as outfile:
        outfile.write(terms)
        flag = False
    os.system("cat terms-* > merged-terms.txt")
    #print(sigterms_list)
