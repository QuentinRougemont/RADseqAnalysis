#!/usr/bin/python
#-*- coding: utf-8 -*-
#######Script python p
#script by G. Lasalle! : Gille.lassalle@rennes.inra.fr

#Pre : repertoire de read1 + rep de read2
#post : lancement auto de clonefilter

#Input : chemin du repertoire des read1 et read2 
#Output :lance clone filter pour tous les couples de fichiers
#date 11 mars 2015

import os
import argparse

#import csv
#import operator

parser = argparse.ArgumentParser()
parser.add_argument('-i', action='store', dest='inputFile', help='Pair Read-files Directory')
#parser.add_argument('-o', action='store', dest='OutDir', help='results dir')


parser.add_argument('--version', action='version', version='%(prog)s 0.1')
results = parser.parse_args()

print 'repertoire de reads		=', results.inputFile
#print 'Out Directory map&fst	=', results.OutDir


##ON Parcours le Dir des fichiers sortis de ProcessRadtags

chemin=os.path.abspath(results.inputFile)  # le chemin absolu
consensusFile=os.path.basename(results.inputFile)  #le nom du fichier
#outFilesDir=os.path.abspath(results.OutDir)

print chemin
print consensusFile



dirs = os.listdir(chemin) #on regarde le contenu du rep
print dirs[0]
lenght=len(dirs)
print "nombre de fichiers théorique : "+str(lenght)
numFile=0
for file1 in dirs:
    if os.path.isfile(os.path.join(chemin, file1)): #traitement des dir et . .. one ne s'occupe que des fichiers  
        fichierR1=file1.split("_")
        #print fichierR1[3] 
        if fichierR1[5] =="1.fq.gz":
            print fichierR1[5]
            file2=str(fichierR1[0])+"_"+str(fichierR1[1])+"_"+str(fichierR1[2])+"_"+str(fichierR1[3])+"_"+str(fichierR1[4])+"_2.fq.gz"
            print file2
            if os.path.isfile(os.path.join(chemin, file2)):
                print "fichier R2 ok !! "+file2            
                commande1="clone_filter -1 "+chemin+"/"+file1+" -2 "+chemin+"/"+file2+" -o "+chemin+"/ -i gzfastq"
                print commande1
                os.system(commande1)
                commande2="rm "+chemin+"/"+file1
                commande3="rm "+chemin+"/"+file2
                print "nettoyage "+file1
                os.system(commande2)
                print "nettoyage "+file2
                os.system(commande3)
