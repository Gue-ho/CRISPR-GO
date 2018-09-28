import os
from CRISPR_GO import *

os.chdir('/home/baelab/Desktop/CRISPR-GO/20170905')

f=open('clinvar_pathogenic_dbSNP_NSN.txt').readlines()
fbio=open('mart_export.txt').readlines()
fw=open('clinvar_pathogenic_NSN_supple.txt','w')
fw.write('dbSNP #\tGenotype\tProtospacer and PAM sequence\tAssociated genetic disease\n')

biomart={}

amino = {'GAC': 'Asp', 'TGG': 'Trp', 'CCT': 'Pro', 'CTT': 'Leu', 'ATT': 'Ile', 'ATC': 'Ile', 'AAA': 'Lys', 'CCA': 'Pro', 'GTT': 'Val', 'AGG': 'Arg', 'CTC': 'Leu', 'CAT': 'His', 'GGT': 'Gly', 'CGT': 'Arg', 'GTA': 'Val', 'TTT': 'Phe', 'ATG': 'Met', 'TGA': 'Ter', 'CAG': 'Gln', 'GCG': 'Ala', 'CAA': 'Gln', 'GGC': 'Gly', 'TTA': 'Leu', 'CTG': 'Leu', 'TAT': 'Tyr', 'TCC': 'Ser', 'TAA': 'Ter', 'CGG': 'Arg', 'ATA': 'Ile', 'CGA': 'Arg', 'TCG': 'Ser', 'TAG': 'Ter', 'GGA': 'Gly', 'AAG': 'Lys', 'TGC': 'Cys', 'TCA': 'Ser', 'GCC': 'Ala', 'AAC': 'Asn', 'AGA': 'Arg', 'ACA': 'Thr', 'GCT': 'Ala', 'GTG': 'Val', 'GTC': 'Val', 'GAA': 'Glu', 'AGC': 'Ser', 'GAT': 'Asp', 'CGC': 'Arg', 'GAG': 'Glu', 'TCT': 'Ser', 'TAC': 'Tyr', 'GGG': 'Gly', 'ACC': 'Thr', 'TGT': 'Cys', 'ACT': 'Thr', 'CCG': 'Pro', 'CTA': 'Leu', 'CAC': 'His', 'ACG': 'Thr', 'AAT': 'Asn', 'GCA': 'Ala', 'TTG': 'Leu', 'AGT': 'Ser', 'TTC': 'Phe', 'CCC': 'Pro'}

for i in fbio[1:]:
    try:
        test=biomart[i.split()[0]]
    except:
        biomart[i.split()[0]]=i.split()[1]

pam_list=['GG','AG','GA','GC','GT','GAN','AA'] #SpCas9 frist!!!

cnt_dic={'disease':0, 'per_all':0, 'all':0, 'aa_per_all':0}
for i in pam_list+['sp']:
    cnt_dic[i]=0
    cnt_dic[i+'_per']=0
    cnt_dic[i+'_aa_per']=0

err_dic={'nogene':0, 'err':0}

repaired_codon={'TGG':0, 'CGA':0, 'CAG':0, 'CAA':0}

for i in f:
    cnt_dic['disease']+=1
    per=0
    seq_list=i.split()[8:]
    site=-1
    gene=0
    rb=''
    gRNA_list=[]

    info_dic=info_extractor(i.split()[7])
    
    gene=info_dic['GENEINFO']
    gene=gene[:gene.find(':')]

    try: info_dic['GENEINFO']
    except: 
        err_dic['err']+=1
        continue

    try: strand=biomart[gene]
    except:
        err_dic['nogene']+=1
        continue
    strand=1
    potential_list=[[],[],[],[],[]] #popam, repopam, popam_1, repopam_1, stop
    
    for x in range(int(len(seq_list)/2)):
        seq=seq_list[2*x]
        codon=seq_list[2*x+1]
        for stop in ['TAA','TAG','TGA']:
            n=0
            stop_find = stop_finder(seq, strand, stop)
            for po in stop_find[1:]:
                potential_list[n].append(stop_find[0]+','+po+','+codon)
                n+=1
    
    if potential_list==[[],[],[],[],[]]: continue

    pam_c=0
    repam_c=0
    per_c=0
    aa_per_c=0
    sp_per_c=0
    sp_aa_per_c=0
    repaired_c=[0,0,0,0] #TGG, CGA, CAG, CAA

    for pam in pam_list:
        if pam=='GAN':
            (cnt_dic, pam_c, repam_c, per_c, aa_per_c, gRNA_list, repaired_c, sp_per_c, sp_aa_per_c) = pam_finder(potential_list[2], potential_list[3], 'GAN', cnt_dic, pam_c, repam_c, per_c, aa_per_c, i.split()[3], gRNA_list, repaired_c, sp_per_c, sp_aa_per_c)
        else: (cnt_dic, pam_c, repam_c, per_c, aa_per_c, gRNA_list, repaired_c, sp_per_c, sp_aa_per_c) = pam_finder(potential_list[0], potential_list[1], pam, cnt_dic, pam_c, repam_c, per_c, aa_per_c, i.split()[3], gRNA_list, repaired_c, sp_per_c, sp_aa_per_c)

    if repaired_c[0]>=1: repaired_codon['TGG']+=1
    if repaired_c[1]>=1: repaired_codon['CGA']+=1
    if repaired_c[2]>=1: repaired_codon['CAG']+=1
    if repaired_c[3]>=1: repaired_codon['CAA']+=1
    if per_c==1: cnt_dic['per_all']+=1
    if pam_c==1 or repam_c==1: cnt_dic['all']+=1
    if aa_per_c==1: cnt_dic['aa_per_all']+=1
    
    mut_symbol=info_dic['CLNHGVS']
    mut_symbol=mut_symbol[mut_symbol.find(':')+1:]

    disease_name=info_dic['CLNDBN']
    if disease_name.find('|')!=-1: disease_name=disease_name[:disease_name.find('|')]

    fw.write(info_dic['RS']+'\t'+gene+';'+mut_symbol+'\t')
    for seq in gRNA_list:
        if seq==gRNA_list[-1]:
            fw.write(seq)
            break
        fw.write(seq+',')
    fw.write('\t'+disease_name+'\n')

fw.close()

print('All : '+str(cnt_dic['disease']))
print('correct : '+str(cnt_dic['all']))
print('per cor : '+str(cnt_dic['per_all']))
print('aa per cor : '+str(cnt_dic['aa_per_all']))

for pam in pam_list+['sp']:
    print(pam+' : '+str(cnt_dic[pam]))
    print(pam+'_per : '+str(cnt_dic[pam+'_per']))
    print(pam+'_aa_per : '+str(cnt_dic[pam+'_aa_per']))

print(repaired_codon)
