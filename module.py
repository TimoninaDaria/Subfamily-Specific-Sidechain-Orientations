def names_of_prot(fasta):
    fassta = fasta.read().split(">")
    names_of_prot = []
    k = 0
    for i in fassta:
        arr = i.split('\n')
        names_of_prot.append(arr[0])
        arr[0] = ''
        fassta[k] = ''
        for n in arr:
            fassta[k] += n
        k += 1
    del fassta[0]
    del names_of_prot[0]
    fasta.seek(0)
    return names_of_prot
  
def seq_from_fasta(fasta):
    fassta = fasta.read().split(">")
    k = 0
    for i in fassta:
        arr = i.split('\n')
        arr[0] = ''
        fassta[k] = ''
        for n in arr:
            fassta[k] += n
        k += 1
    del fassta[0]
    fasta.seek(0)
    return fassta



