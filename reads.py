
import sys

#cree 2 listes comportant les k-mers avec leur nombres dans l'autre liste
def create_k_mer(fichier,len_k_mer=30):
    with open(fichier,'r') as input_f:
        ind=0
        dataset=[]
        nb_dataset=[]
        k_mer=""
        for line in input_f:
            line.strip()
            if not(ind==0):
                for i in range(len(line)-len_k_mer-1):
                    k_mer=""
                    for j in range(len_k_mer):
                        k_mer+=line[i+j]
                    #si le k-mer est déjà dans le dataset, incrémente son compte de 1
                    if dataset.__contains__(k_mer):
                        nb_dataset[dataset.index(k_mer)]+=1
                    #si le k_mer n'est pas dans le dataset, le rajoute au dataset et ajoute son compte 
                    else:
                        dataset.append(k_mer)
                        nb_dataset.append(1)
            ind+=1
        return (dataset,nb_dataset)

#cree un dictionnaire avec comme clefs les k-mers et pour valeurs leur nombre d'occurences et les sequences contenant ces k-mers
def create_k_mer_Dico(fichier,len_k_mer=30):
    with open(fichier,'r') as input_f:
        ind=0
        dico={}
        k_mer=""
        for line in input_f.readlines():
            if (ind%10000==0):
                print('Ligne n°'+str(ind))
            if ">" in line:
                header=line.strip()[1:]
            for i in range(len(line)-len_k_mer-1):
                k_mer=""
                for j in range(len_k_mer):
                    k_mer+=line[i+j]
                #si le k-mer est déjà dans le dataset, incrémente son compte de 1
                if k_mer in dico:
                    dico[k_mer][0]+=1
                    dico[k_mer][1].append(header)
                #si le k_mer n'est pas dans le dataset, le rajoute au dataset et ajoute son compte 
                else:
                    dico[k_mer]=[1,[header]]
            ind+=1
        return dico

def mal_sequences(dictio):
    results=set()
    for key,values in dictio.items():
        if (values[0]<=2):
            results.update(values[1])
    return results

def main():
    len_k_mer=30
    fichier="reads.fasta"
    #REPO="C:\\Users\\mivip\\OneDrive\\Bureau\\Travail_cours\\2A\\algo2A\\projet\\"+fichier
    """
    DATA=create_k_mer(REPO,len_k_mer)
    print(len(DATA[0]),len(DATA[1]))
    print(DATA[0][0],DATA[1][0])
    """
    Dictio=create_k_mer_Dico(fichier,len_k_mer)
    print(len(Dictio.keys()),len(Dictio.values()))
    results=mal_sequences(Dictio)
    print(len(results))
    print(results)

main()
sys.exit(0)
