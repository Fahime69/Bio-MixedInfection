import os
import sys
import numpy as np
import io
import math
import collections
############################## Split PE-family ####################
Amb_low=[] #location of PE-family
Amb_high=[]
with open('decor.gff') as check:
     for line in check:
             info = line.split('\t')
             t=info[8].split(';')
             p=t[1].split('=')
             name=p[1].split(' ')
             if str(name[0])=='PPE' or str(name[0])=='PE' or str(name[0])=='PE-PGRS':
                 
                  high=info[4]
                  low=info[3]
                  Amb_low.append(int(low))
                  Amb_high.append(int(high)) 

############################# Read final-drug file #############################
gene=[]
position=[]
mutation_code=[]
drug=[]

with open ('drugs_final2.txt') as check:
    for line in check:
        if line[0] != '#':
            
            subline = line.split('\t')
            pos = subline[0].split(':')
            position.append(pos[0])
            drug.append(subline[2].rstrip())
            mutation_code.append(pos[1])
            annotation= subline[1].split(',')
            if len(annotation)>2 :
               
                t=annotation[3].split(' ')
                gene.append(t[len(t)-1])
            else:
                t=annotation[0].split('*')
                gene.append(t[1])
               
te = np.array(gene)
uniq_gene, indices = np.unique(te, return_index=True)  
gene_number = len(uniq_gene)

#################################### type of drug columns ##################################
# O 	Ofloxacin 	=column[0]
# L 	linezolid 	=column[1]
# I 	Isoniazid 	=column[2]
# S 	streptomycin=column[3]
# E 	ethambutol 	=column[4]
# C 	cycloserine =column[5] 	
# R 	rifampicin 	=column[6]
# M 	ethionamide =column[7]	
# B 	bedaquiline =column[8]
# F 	clofazamine =column[9]
# K 	kanamycin 	=column[10]
# A 	capreomycin =column[11]	
# P 	pyrazinamide =column[12]	
# F 	clofazamine =column[13]	
# W 	beta-lactams =column[14]	
# Q 	PAS (4-aminosalicylic acid) =column[15]	
# D 	common to aminoglycosides 	=column[16]
# T 	common to isoniazid and ethionamide =column[17] 	
# X 	rifampin compensatory 	=column[18]
# Y 	isoniazid compensatory 	=column[19]
# Z 	ethambutol compensatory =column[20]	
# G 	isoniazid secondary compensatory =column[21]
################################## type of rows ############################
#row[1] = 'GyrB'
#row[2] = 'EmbB'
#row[3] =  Gid'
#row[4] = 'EthA'
#row[5] = 'GyrA'
#row[6] = 'KasA'
#row[7] = 'InhA'
#row[8] = 'PncA'
#row[9] = 'RpoB'
#row[10] = RpsL'
#row[11] = 'KatG'
drug_ind = {}
gene_ind = {}
drug_ind = {'A':0,'B':1,'C':2,'D':3,'E':4,'F':5,'G':6,'I':7,'K':8,'L':9,'M':10,'N':11,'O':12,'Q':13,'R':14,'S':15,'T':16,'X':17,'W':18,'Y':19,'Z':20}
gene_ind= {'GyrB':0,'EmbB':1,'Gid':2,'EthA':3,'GyrA':4,'KasA':5,'InhA':6,'PncA':7,'RpoB':8,'RpsL':9,'KatG':10,'':11}
al_ind = {'A':0,'C':1,'G':2,'T':3}

file0 = open('all_information.txt','w')
file0.write('Filename \t Suseptible \t Resistance \t Het-resistance \t Amb-resistance \t Drug-resistance '+'\n')
###########################   Read vcf files ##################################
with open('Alldir-com.txt') as allfile:

  for dir in allfile:
     vcfname1=dir.split('/')
     aa=vcfname1[len(vcfname1)-2]
     bb=vcfname1[len(vcfname1)-1]
     br=bb.split('.')
     vcfname=aa+br[0]
     AF1=[]
     AF2=[]
     mutatiom = []
     vcfmutation = []
     drug_counter_all = [[0 for x in range(21)] for y in range(12)]
     dca =[[0 for x in range(21)] for y in range(12)]
     with open(dir.rstrip()) as vcffile:
          num_amb=0 #number of ambiguous base
          num_all_line_vcf=0 #number of all vcf line
          num_lowcov=0 #number of low coverage base
          num_del=0 #number of delete base
          num_pass=0 #number of pass base
          num_indel=0 #number of base which len(ref) != len(alt) != 1
          num_mixed=0 #number of base mixed filter (amb,lowcov)(amb,del) ...
          num_PE_family=0# number of ambiguous base which is member of PE family
          num_hetero=0 #number of base which have heterozigous > 0.1
          
          line_count = 0 # number of lines from vcf file which allele frequency and heterozygous are computed
          num_of_Hetero_resistance = 0
          num_amb_res = 0
          counter_res = 0
          for line in vcffile: 
            if line[0] !='#':
                 num_all_line_vcf=num_all_line_vcf+1
                 base=line.split('\t')
                 filt=base[6]
                 filter=filt.split(';')
                 Poss=int(base[1])
               
                 if len(filter)==1:
                    if filter[0] == 'LowCov':
                         num_lowcov=num_lowcov+1
                    if filter[0] == 'Del':
                         num_del=num_del+1
                    if filter[0] == 'PASS' or filter[0] == 'Amb':

                        if filter[0] == 'Amb':
                            num_amb=num_amb+1
                       
                        if filter[0]== 'PASS':
                               num_pass=num_pass+1
############################################ check Ambiguous bases from PE- family ####################################
                                               
                        flag=0 # flag for PE-family
                        i=0
                        while i<len(Amb_high):
                            
                                lo=Amb_low[i]
                                hi=Amb_high[i]
                                i=i+1    
                                if Poss>lo and Poss<hi:
                                   num_PE_family=num_PE_family+1
                                   flag=1
############################################## Compute Allele Frequency ###################################
                        if len(base[3]) != len(base[4]): # count bases which have insert or delete
                                     num_indel= num_indel+1
                        
                        if len(base[3]) == len(base[4]) == 1 and flag ==0:
                            line_count +=1
                            info = base[7].split(';')
                            dp=info[0].split('=')
                            alt=base[4]
                            ref=base[3]
                            if alt == '.':
                                alt = ref
                                								
                            t=info[5].split('=')
                            bc=t[1].split(',')
                            d = int(bc[0])+int(bc[1])+int(bc[2])+int(bc[3])
                            if d != 0:
                                A=(float(float(bc[0])/float(d)))
                                C=(float(float(bc[1])/float(d)))
                                G=(float(float(bc[2])/float(d)))
                                T=(float(float(bc[3])/float(d)))
                                fr_al = [A,C,G,T]
                                
                                indices = [wq for wq, x in enumerate(position) if int(x) == Poss]
                                for ind in indices:
                                    temp1 = mutation_code[ind]
                                    temp2 = temp1.split('->')
                                    ref_D = temp2[0]
                                    alt_D = temp2[1]
                                    flag2 = 0
                                   
                                    
                                    
                                    if alt == alt_D:
                                        counter_res +=1
                                    if ref_D == ref and len(alt_D)==1 and len(ref_D) ==1:
                                        ref_ind = al_ind[ref_D]
                                        alt_ind = al_ind[alt_D]
                                        print fr_al[ref_ind],fr_al[alt_ind]
                                        if fr_al[ref_ind] == float(0.5) and fr_al[alt_ind] == float(0.5):
                                            num_of_Hetero_resistance +=1
                                            flag2 = 1											
                                        drug_type = drug[ind]
                                        gene_type = gene[ind]
                                        mutatiom.append(temp1)
                                        
                                        tempo = str(ref)+'->'+str(alt)
                                        vcfmutation.append(tempo)
                                        if flag2 == 1:

                                            drug_counter_all [gene_ind[gene_type]][drug_ind[drug_type]] = drug_counter_all[gene_ind[gene_type]][drug_ind[drug_type]] +1
                                        if flag2 == 0:
                                            num_amb_res += 1 
                                        AF1.append(fr_al[ref_ind])
                                        AF2.append(fr_al[alt_ind])
                                        if alt == alt_D:
                                            u1= gene_ind.get(gene_type,None)
                                            u2 = drug_ind.get(drug_type,None)
                                            if u1 != None and u2 != None:
                                                dca[gene_ind[gene_type]][drug_ind[drug_type]] = dca[gene_ind[gene_type]][drug_ind[drug_type]] + 1
                                                print gene_ind[gene_type] , drug_ind[drug_type]
     file1 = open(str(vcfname) + '-DrugInfo.txt','w')
     for q in range(10):
         for p in range(21):
                file1.write(str(drug_counter_all[q][p]) + '\t')
         file1.write('\n')
     file1.close()     
     file2 = open(str(vcfname) + '-info.txt','w')
     file2.write( str(line_count) + 'Number of line' + '\n'+\
                  str(num_of_Hetero_resistance) + 'Number of hetero resistance bases'+'\n'+\
                  str(num_amb_res) + 'Number of ambiguous bases')
     file2.close()
     file3 = open(str(vcfname) + '-Freqency.txt','w')
     for km in range(len(AF1)):
         file3.write('Mutation:'+str(mutatiom[km])+'\t In the sample: \t'+str(vcfmutation[km])+'\t Frequency: \t'+str(AF1[km]) +'->'+ str(AF2[km])+'\n')
                     
     file3.close()
     a1 =[sum([row[i] for row in dca]) for i in range(0,len(dca[0]))]
     
     file0.write(str(vcfname)+'\t'+str(len(AF1)) + '\t' + str(counter_res)+'\t'+ str(num_of_Hetero_resistance)+'\t'+str(num_amb_res)+'\t'+str(a1[0])+'-'+str(a1[1])+'-'+str(a1[2])+'-'+str(a1[3])+'-'+str(a1[4])+'-'+str(a1[5])+'-'+str(a1[6])+'-'+str(a1[7])+'-'+str(a1[8])+'-'+str(a1[9])+'-'+str(a1[10])+'-'+str(a1[11])+'-'+str(a1[12])+'-'+str(a1[13])+'-'+str(a1[14])+'-'+str(a1[15])+'-'+str(a1[16])+'-'+str(a1[17])+'-'+str(a1[18])+'-'+str(a1[19])+'-'+str(a1[20])+'\n')
file0.close()       