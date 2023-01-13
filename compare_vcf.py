import pandas as pd
import sys
import csv

args = sys.argv
str_args1 = str(args[1])
str_args2 = str(args[2])
str_args3 = str(args[3])
str_args4 = str(args[4])


def vcfparse(filename):
    df = pd.read_table(filename,comment="#",header=None)
    data = df.iloc[:,[2,9,10,11,12,13,14,15,16]]
    for i in range(1,len(data.columns)):
        data.iloc[:,i] = data.iloc[:,i].str[:3]
    data_=data.set_index(2)
    return data_


def vcfcompare(data_, data_1):
    i_list = data_.index.tolist()
    i_list1 = data_1.index.tolist()
    i_data1_2 = set(i_list)&set(i_list1)
    i_data1_2 = list(i_data1_2)
    share_allele1 = data_1.loc[i_data1_2]
    share_allele = data_.loc[i_data1_2]
    equal_loci =[]
    for loci in i_data1_2:
        a = data_1.loc[loci]
        b = data_.loc[loci]
        if a.equals(b) == True:
            equal_loci.append(loci)
    return equal_loci


data_= vcfparse(str_args1)
data_1 = vcfparse(str_args2)
data_2 = vcfparse(str_args3)
data_3 = vcfparse(str_args4)

pair12 = vcfcompare(data_,data_1)
pair13 = vcfcompare(data_,data_2)
pair14 = vcfcompare(data_,data_3)
pair23 = vcfcompare(data_1,data_2)
pair24 = vcfcompare(data_1,data_3)
pair34 = vcfcompare(data_2,data_3)

share_loci12 = data_.loc[pair12]
share_loci34 = data_2.loc[pair34]
pair1234 = vcfcompare(share_loci12,share_loci34)

length = []

length.append(len(pair12))
length.append(len(pair13))
length.append(len(pair14))
length.append(len(pair23))
length.append(len(pair24))
length.append(len(pair34))
length.append(len(pair1234))


with open('length.csv', 'w') as f:
    writer = csv.writer(f)
    writer.writerow(length)



str_ = '\n'.join(pair1234)
with open("share_loci.txt", 'wt') as f:
    f.write(str_)


