def fasta_bind_site(c, seq1, seq2):
    a = []
    for num in range(len(c)):
        m = -1
        k1 = 0
        while m != c[num]:
            if seq1[k1] != '-':
                m += 1
            k1 += 1
        k1 = k1 - 1
        m = 0
        k2 = -1
        if seq2[k1] != '-':
            while m != k1+1:
                if seq2[m] != '-': k2 += 1
                m += 1
        elif seq2[k1] == '-': k2 = -1

        a.append(k2)
    return(a)
