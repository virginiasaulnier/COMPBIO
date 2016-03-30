import pdb
import os
from Bio import Entrez
import time 
import urllib2
import re
print time.time()
Entrez.email = "jeffng29@gmail.com"
os.chdir ('F:\COMP 383 Final Project')
##needs to be changed to Version 2, not Version 3. Version 2 is used in HapMap
reader = open('HumanWG-6_V3_0_R3_11282955_A.txt','r')
##gene ids
reader1 = open ('idsTable.txt','r')

i=0
#Moves to the correct line in the HumanWG file, skips comment lines
while i<1:
    line = reader.readline()
    if line.startswith('[Probes]'):
        i+=1

i=0

line= reader.readline()
line1 = reader1.readline()
#Two dictionaries to hold ENSG/chromo information from idsTable
ilmntoENSG = {}
ilmnchromo = {}
#Extracts data from idsTable, placing the information into the various dictionaries
#tests time
print time.time()
#pulls data from table
while i < 1:
    line1 = reader1.readline()
    if line1 == '':
        i+=1
    else:
        linelist1= line1.strip().split()
        ilmntoENSG[linelist1[1]] = ''
        ilmnchromo[linelist1[1]]=''
        if linelist1[2]!='NA':
            ilmnchromo[linelist1[1]] = linelist1[3]
            if linelist1[4].startswith('"'):
                index = len(linelist1) -1
                actualid = linelist1[index].strip('"')
                ilmntoENSG[linelist1[1]] = actualid
            else:
                ilmntoENSG[linelist1[1]]=linelist1[4]
#testers 
i=0
hi=0
hihi=0
while i<1:
    d=0
    name=''
    
    line = reader.readline()
    if line == '':
        i+=1
    else:
        linelist = line.strip().split()
        x =0
        j=0
        xlist = linelist[0:25]
        #looks for chromosome position data, marked by a '-', returns index j in xlist
        while x < len(xlist):
            if '-' in xlist[j]:
                if xlist[j][0].isdigit():
                    break
            x+=1
            j+=1
        #gets the ILMN id, make sure it matches with one from idTables
        for y in linelist:
            if y.startswith('ILMN'):
                if y in ilmntoENSG:
                    name = y
                    d+=1
        x=0
      
        
        #If there is known position/chromosome data, if not known (j=len(xlist)), will search ncbi to look for one
        if d==1 and j!= len(xlist):
            chromo = linelist[j]
            counter = str(linelist[j]).count(':')
            if counter >=1:
                x=0
                chromopos=''
                while x <= counter:
                    chromopos = chromopos + chromo.split(':')[x] + ' '
                    x+=1
                print 'chr' + ilmnchromo[name]+ ' '+chromopos + ilmntoENSG[name]
                hihi+=1
            else:
                chromopos1 = chromo.split('-')[0]
                chromopos2 = chromo.split('-')[1]
                hihi+=1
                print 'chr'+ ilmnchromo[name] + ' ' +chromopos1 + ' ' + chromopos2 + ' ' + ilmntoENSG[name]
        elif name!='':
            Entrez.email = 'jeffng29@gmail.com'
            ilmn = name
            handle = Entrez.esearch('geoprofiles', ilmn)
            record=Entrez.read(handle)
            id = record['IdList'][1]
            req = urllib2.Request('http://www.ncbi.nlm.nih.gov/geoprofiles/?term=' + id)
            response = urllib2.urlopen(req)
            page=response.read()
            x = page.find('http://www.ncbi.nlm.nih.gov/geoprofiles?Db=gene&amp;DbFrom=geoprofiles&amp;Cmd=Link&amp;LinkName=geoprofiles_gene&amp;')
            if x != -1:
                rest = page[x:]
                end = rest.find ('">')
                url = rest[0:end]
                nameend = rest.find('</a>')
                namestart = end+2
                name = rest[namestart:nameend]
                handle1 = Entrez.esearch('gene',name)
                record=Entrez.read(handle1)

                req = urllib2.Request(url)
                response = urllib2.urlopen(req)
                page= response.read()
                x = page.find('<select name="EntrezSystem2.PEntrez.Gene.Gene_ResultsPanel.Gene_RVFullReport.Gene_GenomicProductsP.accessionList" sid="1" id="accessionList">')

                rest = page[x:]
                start = rest.find ('<option value="')
                rest1 = rest[start:]
                end = rest1.find ('">')

                idandpos = rest1[0:end]
                listid = idandpos.split('_')
                idtag = listid[0][(len(listid[0]))-2:] + '_'+listid[1]
                pos1 = listid[2]
                pos2 = listid[3]

                chromostart = rest1.find(idtag+ ' ')
                length = len(idtag)+1
                chromostart = chromostart + length
                start= rest1[chromostart:]

                chromo = rest1[chromostart:(chromostart+13)].strip()
                print chromo + ' ' + pos1 + ' '+pos2
                hihi+=1

            else:
                hi+=1
                print hi
print time.time()
print hihi
print hi

    
