import random
import time

import pickle

import numpy as np

from functools import reduce

import pandas as pd
#import homelette as hm

#from modeller import *
#from modeller.automodel import *    # Load the automodel class

from Bio.Align.Applications import ClustalwCommandline

from Bio import AlignIO

from io import StringIO


import os



import re
 
 
def replace_substring(test_str, sb, ss):
    # Replacing all occurrences of substring s1 with s2
    test_str = re.sub(sb, ss, test_str)
    return test_str
 

    

def replacer(s, newstring, index, nofail=False):
            # raise an error if index is outside of the string
        if not nofail and index not in range(len(s)):
            raise ValueError("index outside given string")

            # if not erroring, but the index is still not in the correct range..
        if index < 0:  # add it to the beginning
            return newstring + s
        if index > len(s):  # add it to the end
            return s + newstring

            # insert the new string between "slices" of the original
        return s[:index] + newstring + s[index + 1:]
    
    
    
def find_all(a_str, sub):
    start = 0
    while True:
        start = a_str.find(sub, start)
        if start == -1: return
        yield start
        start += len(sub) # use start += 1 to find overlapping matches




def check_sets(list1, list2):
    set1 = set(list1)
    set2 = set(list2)
 
    if set1 == set2:
        return True
    else:
        return False
 

track10 = ['0']
track2 = []
track2.append(track10)



main = 0
Spike_seq = 'MFIFLLFLTLTSGSDLDRCTTFDDVQAPNYTQHTSSMRGVYYPDEIFRSDTLYLTQDLFLPFYSNVTGFHTINHTFDNPVIPFKDGIYFAATEKSNVVRGWVFGSTMNNKSQSVIIINNSTNVVIRACNFELCDNPFFAVSKPMGTQTHTMIFDNAFNCTFEYISDAFSLDVSEKSGNFKHLREFVFKNKDGFLYVYKGYQPIDVVRDLPSGFNTLKPIFKLPLGINITNFRAILTAFSPAQDTWGTSAAAYFVGYLKPTTFMLKYDENGTITDAVDCSQNPLAELKCSVKSFEIDKGIYQTSNFRVVPSGDVVRFPNITNLCPFGEVFNATKFPSVYAWERKKISNCVADYSVLYNSTFFSTFKCYGVSATKLNDLCFSNVYADSFVVKGDDVRQIAPGQTGVIADYNYKLPDDFMGCVLAWNTRNIDATSTGNYNYKYRYLRHGKLRPFERDISNVPFSPDGKPCTPPALNCYWPLNDYGFYTTTGIGYQPYRVVVLSFELLNAPATVCGPKLSTDLIKNQCVNFNFNGLTGTGVLTPSSKRFQPFQQFGRDVSDFTDSVRDPKTSEILDISPCSFGGVSVITPGTNASSEVAVLYQDVNCTDVSTAIHADQLTPAWRIYSTGNNVFQTQAGCLIGAEHVDTSYECDIPIGAGICASYHTVSLLRSTSQKSIVAYTMSLGADSSIAYSNNTIAIPTNFSISITTEVMPVSMAKTSVDCNMYICGDSTECANLLLQYGSFCTQLNRALSGIAAEQDRNTREVFAQVKQMYKTPTLKYFGGFNFSQILPDPLKPTKRSFIEDLLFNKVTLADAGFMKQYGECLGDINARDLICAQKFNGLTVLPPLLTDDMIAAYTAALVSGTATAGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKQIANQFNKAISQIQESLTTTSTALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQAAPHGVVFLHVTYVPSQERNFTTAPAICHEGKAYFPREGVFVFNGTSWFITQRNFFSPQIITTDNTFVSGNCDVVIGIINNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYVWLGFIAGLIAIVMVTILLCCMTSCCSCLKGACSCGSCCKFDEDDSEPVLKGVKLHYT'

#h1ag31 = 'MFVLLVAYALLHIAGCQTTNGLNTSYSVCNGCVGYSENVFAVESGGYIPSDFAFNNWFLLTNTSSVVDGVVRSFQPLLLNCLWSVSGLRFTTGFVYFNGTGRGDCKGFSSDVLSDVIRYNLNFEENLRRGTILFKTSYGVVVFYCTNNTLVSGDAHIPFGTVLGNFYCFVNTTIGNETTSAFVGALPKTVREFVISRTGHFYINGYRYFTLGNVEAVNFNVTTAETTDFCTVALASYADVLVNVSQTSIANIIYCNSVINRLRCDQLSFDVPDGFYSTSPIQSVELPVSIVSLPVYHKHTFIVLYVDFKPQSGGGKCFNCYPAGVNITLANFNETKGPLCVDTSHFTTKYVAVYANVGRWSASINTGNCPFSFGKVNNFVKFGSVCFSLKDIPGGCAMPIVANWAYSKYYTIGSLYVSWSDGDGITGVPQPVEGVSSFMNVTLDKCTKYNIYDVSGVGVIRVSNDTFLNGITYTSTSGNLLGFKDVTKGTIYSITPCNPPDQLVVYQQAVVGAMLSENFTSYGFSNVVELPKFFYASNGTYNCTDAVLTYSSFGVCADGSIIAVQPRNVSYDSVSAIVTANLSIPSNWTTSVQVEYLQITSTPIVVDCSTYVCNGNVRCVELLKQYTSACKTIEDALRNSAMLESADVSEMLTFDKKAFTLANVSSFGDYNLSSVIPSLPRSGSRVAGRSAIEDILFSKLVTSGLGTVDADYKKCTKGLSIADLACAQYYNGIMVLPGVADAERMAMYTGSLIGGIALGGLTSAASIPFSLAIQSRLNYVALQTDVLQENQKILAASFNKAMTNIVDAFTGVNDAITQTSQALQTVATALNKIQDVVNQQGNSLNHLTSQLRQNFQAISSSIQAIYDRLDIIQADQQVDRLITGRLAALNVFVSHTLTKYTEVRASRQLAQQKVNECVKSQSKRYGFCGNGTHIFSLVNAAPEGLVFLHTVLLPTQYKDVEAWSGLCVDGRNGYVLRQPNLALYKEGNYYRITSRIMFEPRIPTIADFVQIENCNVTFVNISRSELQTIVPEYIDVNKTLQELSYKLPNYTVPDLVVEQYNQTILNLTSEISTLENKSAELNYTVQKLQTLIDNINSTLVDLKWLNRVETYIKWPW'


print(len(Spike_seq))
delt = 0
cov_2 = 'MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT'

while main < 8000000000:
    try: 
        step = 0
        length = 0
        s2 = Spike_seq
        track1 = []
        while step < 100:


            if length == 0:

                lp1 = ['D','P','D','D','P','I']


                st1 = random.choice(lp1)
                #print(st1)



                if st1 == 'P':

                    dfp = pd.read_csv('SNPs_weights.csv')

                    weighted_sample = dfp.sample(n=100, weights="Weights", random_state=np.random.RandomState())
                    #weighted_sample = weighted_sample.sort_values('Frequency', ascending=False)
                    #print(weighted_sample)
                    con1 = weighted_sample['Org'].tolist()
                    con2 = weighted_sample['Mut'].tolist()


                    for ci in con1:
                        aw = list(find_all(s2, ci))

                        if len(aw) > 0:

                            ind1 = random.choice(aw)

                            indz1 = con1.index(ci)


                            result_seq = replacer(s2, con2[indz1] , ind1)
                            break



                if st1 == 'I':

                    dfi = pd.read_csv('INDEL_1_weights.csv')

                    weighted_sample = dfi.sample(n=500, weights="Weights", random_state=np.random.RandomState())


                    con1 = weighted_sample['Org'].tolist()
                    con2 = weighted_sample['Mut'].tolist()

                    #print('insertion')
                    for ci in con1:

                        uj = replace_substring(ci, '-', '')
                        #print(uj)
                        #uj = str(ci).replace('-', '')
                        #print(uj)
                        aw = list(find_all(s2, uj))
                        #print(aw)
                        #print(aw)

                        if len(aw) > 0:

                            ind1 = random.choice(aw)

                            indz1 = con1.index(ci)


                            #print(con2[indz1])
                            #print(con1[indz1])
                            #print(s2)
                            #print(len(s2))
                            result_seq = replacer(s2, con2[indz1] , ind1)
                            #print(result_seq)
                            #print(len(result_seq))
                            #print('break')
                            #time.sleep(1000000)
                            #result_seq = result_seq.replace('-', '')

                            result_seq = replace_substring(result_seq, '-', '')
                            #print(result_seq)
                            #print(len(result_seq))
                            #time.sleep(100000)
                            trh = str(con1[indz1]) + '_' + str(con2[indz1])
                            track1.append(trh)
                            break


                if st1 == 'D':
                    #print('D')
                    dfd = pd.read_csv('INDEL_2_weights.csv')

                    weighted_sample = dfd.sample(n=500, weights="Weights", random_state=np.random.RandomState())


                    con1 = weighted_sample['Mut'].tolist()
                    con2 = weighted_sample['Org'].tolist()

                    #print(len(s2))
                    #print(s2)
                    
                    for ci in con1:
                        #print(ci)
                        
                        uj = replace_substring(ci, '-', '')
                        aw = list(find_all(s2, uj))

                        if len(aw) > 0:

                            ind1 = random.choice(aw)

                            indz1 = con1.index(ci)

                            #print(con2[indz1])
                            #print(con1[indz1])
                            
                            result_seq = replacer(s2, con2[indz1] , ind1)
                            result_seq = s2.replace(con1[indz1], con2[indz1])
                            #print(result_seq)
                            #print(len(result_seq))
                            #time.sleep(1000000)

                            #result_seq = result_seq.replace('-', '')
                            result_seq = replace_substring(result_seq, '-', '')
                            #print(result_seq)
                            
                            trh = str(con1[indz1]) + '_' + str(con2[indz1])
                            track1.append(trh)
                            break

                #print(trh)
                s2 = result_seq
                #print(len(s2))
                #print(s2)
                #time.sleep(10000)
                #print(len(s2))
                if (len(s2) > 1269) and (len(s2) < 1275):
                    if length == 0:
                        length = 11 
                        fZ = 1
                        


                        
                        
                        for tr1 in track2:
                            #print(tr1)
                            if check_sets(tr1, track1):
                                fZ = 0
                                print('same')
                                break
                        
                        #print('not same')        
                        if fZ == 1:
                            track2.append(track1)
                            print(len(track1))
                            print(reduce(lambda count, l: count + len(l), track2, 0))
                            print('yup')
                            lk = open('only_sequences_CoVs.txt', 'a+')
                            
                            lk.write('>' +str(main) + '_' + str(step) + '____'+ '\t' + str(track1) +'\n')
                            lk.write(s2 + '\n')
                            lk.close
                    

                            break






            step += 1
        
    except:
        pass   
            
    
    
    print(main)
    main += 1
        




