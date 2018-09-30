ABE_window={'TAA':[39,44,4,10,38,43,5,11], 'TAG':[38,44,4,10,37,43,5,11], 'TGA':[39,45,4,10,38,44,5,11], 'TAA':[39,45,5,10,38,44,6,11], 'CTA':[39,45,5,11,38,44,6,12], 'TCA':[39,45,4,10,38,44,5,11]}

amino = {'GAC': 'Asp', 'TGG': 'Trp', 'CCT': 'Pro', 'CTT': 'Leu', 'ATT': 'Ile', 'ATC': 'Ile', 'AAA': 'Lys', 'CCA': 'Pro', 'GTT': 'Val', 'AGG': 'Arg', 'CTC': 'Leu', 'CAT': 'His', 'GGT': 'Gly', 'CGT': 'Arg', 'GTA': 'Val', 'TTT': 'Phe', 'ATG': 'Met', 'TGA': 'Ter', 'CAG': 'Gln', 'GCG': 'Ala', 'CAA': 'Gln', 'GGC': 'Gly', 'TTA': 'Leu', 'CTG': 'Leu', 'TAT': 'Tyr', 'TCC': 'Ser', 'TAA': 'Ter', 'CGG': 'Arg', 'ATA': 'Ile', 'CGA': 'Arg', 'TCG': 'Ser', 'TAG': 'Ter', 'GGA': 'Gly', 'AAG': 'Lys', 'TGC': 'Cys', 'TCA': 'Ser', 'GCC': 'Ala', 'AAC': 'Asn', 'AGA': 'Arg', 'ACA': 'Thr', 'GCT': 'Ala', 'GTG': 'Val', 'GTC': 'Val', 'GAA': 'Glu', 'AGC': 'Ser', 'GAT': 'Asp', 'CGC': 'Arg', 'GAG': 'Glu', 'TCT': 'Ser', 'TAC': 'Tyr', 'GGG': 'Gly', 'ACC': 'Thr', 'TGT': 'Cys', 'ACT': 'Thr', 'CCG': 'Pro', 'CTA': 'Leu', 'CAC': 'His', 'ACG': 'Thr', 'AAT': 'Asn', 'GCA': 'Ala', 'TTG': 'Leu', 'AGT': 'Ser', 'TTC': 'Phe', 'CCC': 'Pro'}

def rev_comp (seq):
    return seq.translate(seq.maketrans('ATGC','TACG'))[::-1]

def info_extractor (info):
    info_dic={}
    for i in info.split(';'):
        info_dic[i[:i.find('=')]]=i[i.find('=')+1:]
    return info_dic

def stop_finder (seq, strand, stop):
    if strand==-1: stop=rev_comp(stop)
    if seq.find(stop)==-1: return []
    site=seq[23:28].find(stop)
    if site==-1: return []
    l=ABE_window[stop]
    return [seq,seq[l[0]+site:l[1]+site], seq[l[2]+site:l[3]+site], seq[l[4]+site:l[5]+site], seq[l[6]+site:l[7]+site], stop]
    
def pam_finder (seq_list, reverse_list, pam, cnt_dic, pam_c, repam_c, per_c, aa_per_c, before_info, gRNA_list, repaired_c, sp_per_c, sp_aa_per_c):
    pam_strand=0
    pam_reverse=0
    per=0
    aa_per=0
    sp_c=0
    true_pam=0
    tgg=0
    cag=0
    cga=0
    caa=0
    if pam=='GAN':
        true_pam='GAN'
        pam='GA'
    if pam_c==1 or repam_c==1:
        sp_c=1
        if per_c==1:
            sp_per_c=1
    re_pam=rev_comp(pam)
    for seq in seq_list:
        fullseq=seq.split(',')[0]
        codon=int(seq.split(',')[2])
        seq=seq.split(',')[1]
        codon_seq='TGG'
        if seq.find(pam)!=-1:
            ed=seq.find(pam)+len(pam)
            full_ed=fullseq.find(seq)+ed
            pam_strand=1
            pam_c=1
            gRNA_list.append(fullseq[full_ed-23:full_ed])
            repaired_c[0]+=1
            wt_seq=fullseq[:25]+before_info+fullseq[26:]
            wt_codon_seq=wt_seq[25-codon:28-codon]
            
            if fullseq[25-codon:28-codon]!='TAA' and before_info=='G' and fullseq[25]=='A': per=1
            if amino[wt_codon_seq]==amino[codon_seq]: aa_per=1

    for seq in reverse_list:
        fullseq=seq.split(',')[0]
        codon=int(seq.split(',')[2])
        codon_seq=fullseq[25-codon:28-codon]
        seq=seq.split(',')[1]
        codon_seq='C'+codon_seq[1:]
        if seq.find(re_pam)!=-1:
            st=seq.find(re_pam)
            full_st=fullseq.find(seq)+st
            pam_reverse=1
            repam_c=1
            if codon_seq=='CAG': repaired_c[1]+=1
            elif codon_seq=='CGA': repaired_c[2]+=1
            elif codon_seq=='CAA': repaired_c[3]+=1
            gRNA_list.append(rev_comp(fullseq[full_st:full_st+23]))
            if before_info=='C' and fullseq[25]=='T':per=1
            wt_seq=fullseq[:25]+before_info+fullseq[26:]
            wt_codon_seq=wt_seq[25-codon:28-codon]
            if amino[wt_codon_seq]==amino[codon_seq]: aa_per=1
    
    if true_pam=='GAN': pam='GAN'
    if pam_strand==1 or pam_reverse==1:
        cnt_dic[pam]+=1
        reparied=1
        if pam in ['GG', 'AG'] and sp_c==0:
            cnt_dic['sp']+=1
        if per==1:
            per_c=1
            cnt_dic[pam+'_per']+=1
            if pam in ['GG', 'AG'] and sp_per_c==0:
                sp_per_c=1
                cnt_dic['sp_per']+=1
        if aa_per==1:
            aa_per_c=1
            cnt_dic[pam+'_aa_per']+=1
            if pam in ['GG', 'AG'] and sp_aa_per_c==0:
                cnt_dic['sp_aa_per']+=1
                sp_aa_per_c=1
    return cnt_dic, pam_c, repam_c, per_c, aa_per_c, gRNA_list, repaired_c, sp_per_c, sp_aa_per_c

def int_extract (a):
    res=''
    for i in a:
        try:
            int(a)
            res+=a
        except:
            continue
    return res
