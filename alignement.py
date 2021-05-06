#!/usr/bin/env python3

# -*- coding:utf-8 -*-



from Bio.Seq import Seq

from Bio import SeqIO

import numpy as np

from Bio import Entrez

from Bio.Blast import NCBIWWW

from Bio.Blast import NCBIXML

from Bio.SeqUtils import GC





mail = input ("Entrez votre mail please : \n")

Entrez.email = mail



terme = input("\n" + "Saisissez votre requete (par exemple : arabidopsis thaliana [organism] AND SEX1[gene] AND biomol_mRNA [properties] : \n")



ma_req = Entrez.esearch(db="nucleotide", term= terme)

mon_res = Entrez.read(ma_req)



#print(mon_res["Count"])

#print(mon_res["IdList"])

ma_req.close()

list_id = mon_res["IdList"]

id_seq = list_id[0]

fic_seq = Entrez.efetch(db="nucleotide", id=id_seq, rettype="gb")

ma_seq = SeqIO.read(fic_seq,format = "gb")





with open("code.fasta", "w") as code:

	code.write(ma_seq.format("fasta"))





print("\n" + "La sequence est disponible dans le repertoire courant dans le fichier code.fasta"+"\n")

print("\n" + "Voici un apercu de ce fichier :" +"\n")

maSeq = SeqIO.read("code.fasta", "fasta")

print(maSeq)

GC= GC(maSeq.seq)

print("\n"+ "le taux de GC de cette sequence est de : "+ (str(GC)) + " %" +"\n")

fic_seq.close()



inFile = open('code.fasta','r')



maSeq = SeqIO.read('code.fasta','fasta')

result = NCBIWWW.qblast('blastn', 'nt', maSeq.format('fasta'))



blast_record = NCBIXML.read(result)



liste = []

for alignment in blast_record.alignments:

	for hsp in alignment.hsps:

		if hsp.expect<1e-37:

			#if alignment.accession not in liste :

			liste.append(alignment.accession)

			



liste_unique = []

for i in liste:

	if i not in liste_unique:

		liste_unique.append(i)

		



var = liste_unique[0]

print("la sequence la plus similaire pour la sequence etudier possede le numero d'accession suivant : " + liste_unique[0] + " cela correspond à la sequence qui vous ai presente ci-dessous :" + "\n")



ma_req2 = Entrez.esearch(db="nucleotide", term= var)

mon_res2 = Entrez.read(ma_req2)



#print(mon_res2["Count"])

#print(mon_res2["IdList"])

ma_req2.close()

list_id2 = mon_res2["IdList"]

id_seq2 = list_id2[0]

fic_seq2 = Entrez.efetch(db="nucleotide", id=id_seq2, rettype="gb")

ma_seq2 = SeqIO.read(fic_seq2,format = "gb")





with open("seqLaPlusSimilaire.fasta", "w") as seqLaPlusSimilaire:

	seqLaPlusSimilaire.write(ma_seq.format("fasta"))





print("\n" + "Voici un apercu de ce fichier :" +"\n")

maSeq2 = SeqIO.read("seqLaPlusSimilaire.fasta", "fasta")
print("\n" + "La sequence la plus similaire est disponible dans le repertoire courant dans le fichier seqLaPlusSimilaire.fasta"+"\n")

print(maSeq2)

fic_seq.close()



result.close()





maSequence = SeqIO.read('code.fasta', "fasta")

maSequence2 = SeqIO.read('maseqmuter.fasta', "fasta")



#on recuperer uniquement les sequences nucleotidique de nos fichier fasta

x=maSequence.seq

y=maSequence2.seq





def alignement(x, y, match = 2, mismatch = 2, gap = 1):

    long_x = len(x)

    long_y = len(y)
    

    #Initialisation des matrices
    

    M = np.zeros((long_x + 1, long_y + 1)) # Création d'une table avec long_x + 1 lignes et long_y + 1 colonnes

    M[:,0] = np.linspace(0, -long_x, long_x + 1) # Permet d’obtenir un tableau à 1D allant d’une valeur de départ à une valeur de fin avec un nombre donné d’éléments.

    M[0,:] = np.linspace(0, -long_y, long_y + 1) # On va de 0 à -long_y en long_y + 1 valeur



    S = np.zeros((long_x + 1, long_y + 1))

    S[:,0] = 3

    S[0,:] = 4



    # Iic on va calculer les differents score pour choisir le max et affectation dans matrice

    score = np.zeros(3) # Table avec 3 zéros



    for i in range(long_x):

        for j in range(long_y):

            if x[i] == y[j]:

                score[0] = M[i,j] + match

            else:

                score[0] = M[i,j] - mismatch #A la place d'additionner on inverse les signes pour avoir que des chiffres positifs en argument de notre fonction alignement

            score[1] = M[i,j+1] - gap #insertion

            score[2] = M[i+1,j] - gap #deletion

            scoremax = np.max(score)

            M[i+1,j+1] = scoremax #on vient mettre le score max de nos 3 scores calculé dans notre matrice M

            if score[0] == scoremax: #on vient ici en fonction de si la valeur max vient de la premiere valeur de notre liste score, deuxieme ou troisieme afin de retracer l'origine du score max de la case, on va en fonction affecté une valeur ou si ca vient de match ou mismatch, insertion ou deletion, et on vient mettre ce score dans notre seconde matrice S.
            

                S[i+1,j+1] += 2  #match ou mismatch

            if score[1] == scoremax:

                S[i+1,j+1] += 3  #insertion

            if score[2] == scoremax:

                S[i+1,j+1] += 4  #deletion
                

   # Recuperation de l'alignement optimal

    i = long_x

    j = long_y

    ntx = []

    nty = []

    while i > 0 or j > 0:
    
#les chiffre entre crocher qui vont suivre font referance à la matrice S à qui on a affecter des valeurs en fonction de l'origine du score à une certaines position de la matrice M. Dans le cas où ce score peut provenir de deux calculs ou trois different, alors on additione les chiffres mais ces derniers ont été choisit de sorte à ce que chacun traduise un cas unique. Si c'est 2 c'est originaire d'une fleche en diagonal, si c'est 5 c'est soit fleche diagonale ou horizontale, etc.. 9 c'est dans le cas ou le score max est trouvé grace aux trois calcule.


        if S[i,j] in [2, 5, 6, 9]: # Pour la valeur P[i,j] = 2, 5, 6 ou 9

            ntx.append(x[i-1]) # On ajoute dans la liste ntx le nucléotide à la position i-1 (car on a decale de 1 dans une matrice)

            nty.append(y[j-1]) # On ajoute dans la liste nty le nucléotide à la position j-1

            i -= 1 # on decremante pour passer au nt suivant

            j -= 1 

        elif S[i,j] in [3,7]:

            ntx.append(x[i-1])

            nty.append('-')

            i -= 1

        elif S[i,j] in [4]:

            ntx.append('-')

            nty.append(y[j-1])

            j -= 1

    # on reverse pour recuperer le sequence dans l'odre en partant du debut ( vu qu'on est partie de la fin)

    ntx = ''.join(ntx)[::-1]

    nty = ''.join(nty)[::-1]

    return '\n'.join([ntx, nty])







var = open('resultats.txt', 'w')

var.write(alignement(x, y))

print("\n"+"Le fichier resultats est pret, vous pouvez le consulter dans votre repertoir courant.")


# si vous souhaitez lancer le programme sur des sequences plus petites afin d'avoir un résultat plus lisible pour juger de la fonctionnalité du code vous pouvez decomentez les lignes ci dessous afin de generer des sequences:

# x = np.random.choice(['A', 'T', 'G', 'C'], 50)

# y = np.random.choice(['A', 'T', 'G', 'C'], 50)
