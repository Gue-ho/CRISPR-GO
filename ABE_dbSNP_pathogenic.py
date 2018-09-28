import os, time, requests
from CRISPR_GO import *

chro_list=['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y']

amino = {'GAC': 'Asp', 'TGG': 'Trp', 'CCT': 'Pro', 'CTT': 'Leu', 'ATT': 'Ile', 'ATC': 'Ile', 'AAA': 'Lys', 'CCA': 'Pro', 'GTT': 'Val', 'AGG': 'Arg', 'CTC': 'Leu', 'CAT': 'His', 'GGT': 'Gly', 'CGT': 'Arg', 'GTA': 'Val', 'TTT': 'Phe', 'ATG': 'Met', 'TGA': 'Ter', 'CAG': 'Gln', 'GCG': 'Ala', 'CAA': 'Gln', 'GGC': 'Gly', 'TTA': 'Leu', 'CTG': 'Leu', 'TAT': 'Tyr', 'TCC': 'Ser', 'TAA': 'Ter', 'CGG': 'Arg', 'ATA': 'Ile', 'CGA': 'Arg', 'TCG': 'Ser', 'TAG': 'Ter', 'GGA': 'Gly', 'AAG': 'Lys', 'TGC': 'Cys', 'TCA': 'Ser', 'GCC': 'Ala', 'AAC': 'Asn', 'AGA': 'Arg', 'ACA': 'Thr', 'GCT': 'Ala', 'GTG': 'Val', 'GTC': 'Val', 'GAA': 'Glu', 'AGC': 'Ser', 'GAT': 'Asp', 'CGC': 'Arg', 'GAG': 'Glu', 'TCT': 'Ser', 'TAC': 'Tyr', 'GGG': 'Gly', 'ACC': 'Thr', 'TGT': 'Cys', 'ACT': 'Thr', 'CCG': 'Pro', 'CTA': 'Leu', 'CAC': 'His', 'ACG': 'Thr', 'AAT': 'Asn', 'GCA': 'Ala', 'TTG': 'Leu', 'AGT': 'Ser', 'TTC': 'Phe', 'CCC': 'Pro'}

os.chdir('/home/baelab/Desktop/CRISPR-GO/20170905')

f=open('clinvar_20170905.vcf').readlines()

chro='0'
chro_v=-1

cnt_dic={'all':0, 'indel':0, 'silent':0, 'NSN':0, 'NSM':0}

f_nsn=open('clinvar_pathogenic_dbSNP_NSN.txt','w')
f_nsm=open('clinvar_pathogenic_dbSNP_NSM.txt','w')
fe=open('clinavr_pathogenic_dbSNP_error.txt','w')

flen=len(f)
line=0
stand=1
t=time.time()

err_cnt=0
vvv=0

for i in f:

    line+=1

    if i[0]=='#': continue

    isp=i.split()
    
    if isp[0]=='MT': continue
    if isp[0]!=chro:
        chro_v+=1
        chro=chro_list[chro_v]
        fc=open('/home/baelab/Desktop/CRISPR-GO/dbSNP/info_chr_result.txt').readlines()
        snp_dic={}
        vvv=1
        for x in fc:
            xsp=x.split()
            snp_dic[xsp[1]]=[xsp[2],xsp[3],xsp[4],xsp[5],xsp[6][1:],xsp[7]]
    rsid=isp[2][2:]
        
    try: snp_dic[rsid]
    except: continue
        

    wt_nt=isp[3]
    mut_nt=isp[4].split(',')
    try: 
        mut_nt.index('N')
        mut_nt.remove('N')
    except: pass
    
    cnt_dic['all']+=len(mut_nt)

    if len(wt_nt)!=len(mut_nt[0]): 
        cnt_dic['indel']+=len(mut_nt)
        continue
    
    if isp[7].find('NSN')==-1 and isp[7].find('NSM')==-1:
        cnt_dic['silent']+=len(mut_nt)
        continue
    
    try:
        snp_info=snp_dic[rsid]
        snp_ref_nt=snp_info[3][snp_info[3].find('&')-1:snp_info[3].find('&')]
    except:
        continue
        print(line)
        v=0
        while v<5:
            try:
                req=requests.get(url='https://www.ncbi.nlm.nih.gov/snp/?term=rs'+rsid)
                break
            except:
                v+=1
                print('try request url '+str(v))
        html=req.text.split('\n')
        print(rsid)
        for x in html:
            if x.find('Chromosome:')!=-1:
                st=x.find('class="snp_flanks">')+len('class="snp_flanks">')
                seq5=x[st:st+25].upper()
                seq3=x[x.find('</span>',st)+7:x.find('</span>',st)+32].upper()
            if x.find('HGVS')!=-1:
                cv=1
                while cv<15:
                    CDS=x[x.find('c.')+2:x.find('c.')+2+cv]
                    try:
                        int(CDS)
                        cv+=1
                        continue
                    except:
                        CDS=x[x.find('c.')+2:x.find('c.')+1+cv]
                        snp_ref_nt=x[x.find('c.')+2+cv:x.find('c.')+3+cv]
                        break
                break

        if str(CDS).find('*')!=-1 or CDS=='':
            cv=1
            while cv<15:
                CDS=x[x.find('c.')+3:x.find('c.')+3+cv].replace('*','')
                try:
                    int(CDS)
                    cv+=1
                    continue
                except:
                    CDS=x[x.find('c.')+3:x.find('c.')+2+cv].replace('*','')
                    snp_ref_nt=x[x.find('c.')+2+cv:x.find('c.')+3+cv]
                    break
        
        snp_info=[seq5, seq3, CDS]
        
    if len(snp_info[0])<25 or len(snp_info[1])<25:
        pos=int(snp_info[4])
        fchro=open('/home/baelab/Desktop/Genome/Human_GRCh38/Homo_sapiens.GRCh38.dna.chromosome.'+chro+'.fa').readlines()
        row=int(pos/60)+1
        #cseq=(fchro[row-1]+fchro[row]+fchro[row+1]).replace('\n','')
        #rpos=pos-(row-1)*60-1
        cseq=''
        strand=1
        for x in range(row-50,row+50):
            cseq+=fchro[x].replace('\n','')
        rpos=cseq.find(snp_info[0])+len(snp_info[0])
        if rpos==-1+len(snp_info[0]): 
            cseq=cseq.translate(cseq.maketrans('ATGC','TACG'))[::-1]
            rpos=cseq.find(snp_info[0])+len(snp_info[0])
            strand=-1
        seq5=cseq[rpos-25:rpos].upper()
        seq3=cseq[rpos+1:rpos+26].upper()
        if seq5.find(snp_info[0])==-1:
            v=0
            while v<5:
                v+=1
                try:
                    req=requests.get('https://genome.ucsc.edu/cgi-bin/das/hg19/dna?segment=chr'+chro+':'+str(pos-25)+','+str(pos+25))   
                    break
                except:
                    continue
                
            html=req.text.split('\n')
            cseq=''
            for x in html:
                if x.find('<')!=-1: continue
                print('a')
                cseq+=x.replace(' ','')
            print(html)
            print(cseq)
            cseq=cseq.upper()
            if cseq.find(snp_info[0])==-1:
                cseq=cseq.translate(cseq.maketrans('ATGC','TACG'))[::-1]
            seq5=cseq[:25].upper()
            seq3=cseq[26:].upper()
            if seq5.find(snp_info[0])==-1:
                print(snp_info)
                input()
        if seq5[-len(snp_info[0]):]!=snp_info[0]:
            print(snp_info)
            print(seq5)
            print(seq5[:len(snp_info[0])])
            input()
        snp_info=[seq5, seq3, CDS]

    try:
        CDS=int(snp_info[2])%3
    except:
        if len(mut_nt)==1:
            if isp[7].find('NSN'): cnt_dic['NSN']+=1
            if isp[7].find('NSM'): cnt_dic['NSM']+=1
            continue
        else:
            fe.write(i)
            continue

    if CDS==0: CDS=3
    CDS-=1
    
    wt_codon=(snp_info[0][-25:]+wt_nt+snp_info[1][:25])[25-CDS:28-CDS].upper()
    try: wt_aa=amino[wt_codon]
    except:
        print(rsid)
        print(snp_info)

    for nt in mut_nt:
        if snp_ref_nt!=wt_nt: nt=nt.translate(nt.maketrans('ATGC','TACG'))[::-1]
        seq=snp_info[0][-25:]+nt+snp_info[1][:25]
        codon=seq[25-CDS:28-CDS].upper()
        try: aa=amino[codon]
        except:
            if len(mut_nt)==1:
                if isp[7].find('NSM')!=-1: cnt_dic['NSM']+=1
                elif isp[7].find('NSN')!=-1: cnt_dic['NSN']+=1
                else: cnt_dic['silent']+=1
            fe.write(i+'\t'+snp_info[0]+'\t'+snp_info[1]+'\t'+str(snp_info[2])+'\n')
            continue
        if aa=='Ter':
            cnt_dic['NSN']+=1
            f_nsn.write(i.replace('\n','')+'\t'+seq+'\t'+str(CDS)+'\n')
        elif wt_aa!=aa: 
            cnt_dic['NSM']+=1
            f_nsm.write(i.replace('\n','')+'\t'+seq+'\t'+str(CDS)+'\n')
        elif wt_aa==aa: cnt_dic['silent']+=1

    if line*1000/flen > stand:
        print(str(round(0.1*stand,2))+' %\t'+str(round(time.time()-t,2)))
        stand+=1

f_nsn.close()
f_nsm.close()
fe.close()
for i in ['all','indel','silent','NSM','NSN']:
    print(i+'\t:\t'+str(cnt_dic[i]))    
        
