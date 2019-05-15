import re
import os
import subprocess
import numpy as np
import pandas as pd
import Bio.Align.Applications as ba
#from Bio.Align.Applications import ClustalwCommandline
from Bio import Phylo
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import AlignInfo
#import Bio.Align.Applications
import tkinter as tk
from tkinter import filedialog
import time

class Utils():

        # Direccion del ejecutable clustalW, debería leerse desde un archivo config.
        clustalw_exe = r"C:\Program Files (x86)\ClustalW2\clustalw2.exe"
        df = pd.DataFrame({})
        consensus = ''
        alignment = ''
        
        def __init__(self):
                pass

        def select_file(self):
                # para quitar la ventana que aparece
                root = tk.Tk()
                root.withdraw()
                
                # dialogo de seleccion de archivo
                file_path = filedialog.askopenfilename()
                return file_path
        
        def print_secuences(self,secs):
                for seq_record in secs:
                        print(seq_record.id)
                        print(seq_record.seq)
                        print(len(seq_record))

        def exec_clustalW(self,f):
                start = time.time()
                clustalw_cline = ba.ClustalwCommandline(self.clustalw_exe, infile=f)
                stdout, stderr = clustalw_cline()
                end = time.time()
                timeCost = end - start
                print("El alineamiento con CLustalW ha tardado {0:.2f} segundos.".format(timeCost))
                
        def exec_Muscle(self,f):
                start = time.time()
                alignFile = f.replace('.fasta', '.aln')
                muscle_cline = ba.MuscleCommandline(input=f, out=alignFile)
                stdout, stderr = muscle_cline()
                #stdout, stderr = clinetree()
                treeFile = f.replace('.fasta', '.dnd')
                print("Arbol filogenetico")
                clinetree = "muscle -maketree -in " + alignFile + " -out " + treeFile + " -cluster neighborjoining"
                subprocess.call(clinetree)
                end = time.time()
                timeCost = end - start
                print("El alineamiento con Muscle ha tardado {0:.2f} segundos.".format(timeCost))

        def print_phylo_tree(self,f):
                tree = Phylo.read(f, "newick")
                Phylo.draw_ascii(tree)
                
        def get_alignment(self,f,format):
                alignment = AlignIO.read(f, format)
                self.alignment = alignment
                return alignment

        def print_alignment(self, alignment):
                print("Longitud total del alineamiento %i" % alignment.get_alignment_length())
                for record in alignment:
                        print("%s - %s" % (record.seq, record.id))

        def run_alignment(self,tool):
                # Alineamiento
                fastaFile = self.select_file()
                if (tool == "clustalw"):
                        self.exec_clustalW(fastaFile)
                        alignFormat = 'clustal'
                elif (tool == "muscle"):
                        self.exec_Muscle(fastaFile)
                        alignFormat = 'fasta'
                else: 
                        raise Exception('Herramienta de alineamiento no reconocida')
                alignFile = fastaFile.replace('.fasta','.aln')

                alignment = self.get_alignment(alignFile, alignFormat)
                self.print_alignment(alignment)
                # Arbol filogenetico
                treeFile = fastaFile.replace('.fasta', '.dnd')
                self.print_phylo_tree(treeFile)

                # Consensus secuence  
                summary_align = AlignInfo.SummaryInfo(self.alignment)
                self.consensus = summary_align.dumb_consensus(threshold=0.7,require_multiple=1,ambiguous='?')

                return alignment
        
        def get_gene_from_description(self,desc):
                geneField = re.findall('GN=\S+',desc)
                if (len(geneField) == 1):
                        gen = geneField[0].split('=')[1]
                        return gen
                else:
                        return 'unknown_gene'

        
        def load_kinases(self):
                f1 = self.select_file()
                healthySeqRecords = SeqIO.parse(f1,'fasta')
                rows_list =[]
                for seq in healthySeqRecords:
                        dict1 = {}
                        dict1.update({'id' : seq.id,
                                      'state' : 'healthy',
                                      'length': len(seq),
                                      'gene' : self.get_gene_from_description(seq.description)
                        })  
                        rows_list.append(dict1)
                f2 = self.select_file()
                patSeqRecords = SeqIO.parse(f2, 'fasta')
                for seq in patSeqRecords:
                        dict1 = {}
                        dict1.update({'id' : seq.id, 
                                      'state': 'pathologic',
                                      'length': len(seq),
                                      'gene': self.get_gene_from_description(seq.description)
                        })
                        rows_list.append(dict1)
                self.df = pd.DataFrame(rows_list)
                self.df = self.df.set_index('id')

                print(self.df.to_string())
                filenames = [f1,f2]
                fastaout = os.path.dirname(os.path.abspath(__file__)) + '\\allkin.fasta'
                with open(fastaout, 'w') as outfile:
                        for fname in filenames:
                                with open(fname) as infile:
                                        outfile.write(infile.read())

        def setnA(self):
                # Columna nueva y valor por defecto
                self.df['nA'] = -1
                # Recorrer los alineamientos y establecer el numero aminoacidos
                for alignSec in self.alignment:
                        id = alignSec.id
                        nA = alignSec.seq.count('A')
                        self.df.at[id, 'nA'] = nA

        def setnB(self):
                # Columna nueva y valor por defecto
                self.df['nB'] = -1
                # Recorrer los alineamientos y establecer el numero aminoacidos
                for alignSec in self.alignment:
                        id = alignSec.id
                        nB = alignSec.seq.count('B')
                        self.df.at[id, 'nB'] = nB

        def setnC(self):
                # Columna nueva y valor por defecto
                self.df['nC'] = -1
                # Recorrer los alineamientos y establecer el numero aminoacidos
                for alignSec in self.alignment:
                        id = alignSec.id
                        nC = alignSec.seq.count('C')
                        self.df.at[id, 'nC'] = nC

        def setnD(self):
                # Columna nueva y valor por defecto
                self.df['nD'] = -1
                # Recorrer los alineamientos y establecer el numero aminoacidos
                for alignSec in self.alignment:
                        id = alignSec.id
                        nD = alignSec.seq.count('D')
                        self.df.at[id, 'nD'] = nD

        def setnE(self):
                # Columna nueva y valor por defecto
                self.df['nE'] = -1
                # Recorrer los alineamientos y establecer el numero aminoacidos
                for alignSec in self.alignment:
                        id = alignSec.id
                        nE = alignSec.seq.count('E')
                        self.df.at[id, 'nE'] = nE


        def setnF(self):
                # Columna nueva y valor por defecto
                self.df['nF'] = -1
                # Recorrer los alineamientos y establecer el numero aminoacidos
                for alignSec in self.alignment:
                        id = alignSec.id
                        nF = alignSec.seq.count('F')
                        self.df.at[id, 'nF'] = nF

        def setnG(self):
                # Columna nueva y valor por defecto
                self.df['nG'] = -1
                # Recorrer los alineamientos y establecer el numero aminoacidos
                for alignSec in self.alignment:
                        id = alignSec.id
                        nG = alignSec.seq.count('G')
                        self.df.at[id, 'nG'] = nG

        def setnH(self):
                # Columna nueva y valor por defecto
                self.df['nH'] = -1
                # Recorrer los alineamientos y establecer el numero aminoacidos
                for alignSec in self.alignment:
                        id = alignSec.id
                        nH = alignSec.seq.count('H')
                        self.df.at[id, 'nH'] = nH

        def setnI(self):
                # Columna nueva y valor por defecto
                self.df['nI'] = -1
                # Recorrer los alineamientos y establecer el numero aminoacidos
                for alignSec in self.alignment:
                        id = alignSec.id
                        nI = alignSec.seq.count('I')
                        self.df.at[id, 'nI'] = nI

        def setnJ(self):
                # Columna nueva y valor por defecto
                self.df['nJ'] = -1
                # Recorrer los alineamientos y establecer el numero aminoacidos
                for alignSec in self.alignment:
                        id = alignSec.id
                        nJ = alignSec.seq.count('J')
                        self.df.at[id, 'nJ'] = nJ

        def setnK(self):
                # Columna nueva y valor por defecto
                self.df['nK'] = -1
                # Recorrer los alineamientos y establecer el numero aminoacidos
                for alignSec in self.alignment:
                        id = alignSec.id
                        nK = alignSec.seq.count('K')
                        self.df.at[id, 'nK'] = nK

        def setnL(self):
                # Columna nueva y valor por defecto
                self.df['nL'] = -1
                # Recorrer los alineamientos y establecer el numero aminoacidos
                for alignSec in self.alignment:
                        id = alignSec.id
                        nL = alignSec.seq.count('L')
                        self.df.at[id, 'nL'] = nL

        def setnM(self):
                # Columna nueva y valor por defecto
                self.df['nM'] = -1
                # Recorrer los alineamientos y establecer el numero aminoacidos
                for alignSec in self.alignment:
                        id = alignSec.id
                        nM = alignSec.seq.count('M')
                        self.df.at[id, 'nM'] = nM
        
        def setnN(self):
                # Columna nueva y valor por defecto
                self.df['nN'] = -1
                # Recorrer los alineamientos y establecer el numero aminoacidos
                for alignSec in self.alignment:
                        id = alignSec.id
                        nN = alignSec.seq.count('N')
                        self.df.at[id, 'nN'] = nN

        def setnO(self):
                # Columna nueva y valor por defecto
                self.df['nO'] = -1
                # Recorrer los alineamientos y establecer el numero aminoacidos
                for alignSec in self.alignment:
                        id = alignSec.id
                        nO = alignSec.seq.count('O')
                        self.df.at[id, 'nO'] = nO

        def setnP(self):
                # Columna nueva y valor por defecto
                self.df['nP'] = -1
                # Recorrer los alineamientos y establecer el numero aminoacidos
                for alignSec in self.alignment:
                        id = alignSec.id
                        nP = alignSec.seq.count('P')
                        self.df.at[id, 'nP'] = nP

        def setnQ(self):
                # Columna nueva y valor por defecto
                self.df['nQ'] = -1
                # Recorrer los alineamientos y establecer el numero aminoacidos
                for alignSec in self.alignment:
                        id = alignSec.id
                        nQ = alignSec.seq.count('Q')
                        self.df.at[id, 'nQ'] = nQ

        def setnR(self):
                # Columna nueva y valor por defecto
                self.df['nR'] = -1
                # Recorrer los alineamientos y establecer el numero aminoacidos
                for alignSec in self.alignment:
                        id = alignSec.id
                        nR = alignSec.seq.count('R')
                        self.df.at[id, 'nR'] = nR

        def setnS(self):
                # Columna nueva y valor por defecto
                self.df['nS'] = -1
                # Recorrer los alineamientos y establecer el numero aminoacidos
                for alignSec in self.alignment:
                        id = alignSec.id
                        nS = alignSec.seq.count('S')
                        self.df.at[id, 'nS'] = nS

        def setnT(self):
                # Columna nueva y valor por defecto
                self.df['nT'] = -1
                # Recorrer los alineamientos y establecer el numero aminoacidos
                for alignSec in self.alignment:
                        id = alignSec.id
                        nT = alignSec.seq.count('T')
                        self.df.at[id, 'nT'] = nT

        def setnU(self):
                # Columna nueva y valor por defecto
                self.df['nU'] = -1
                # Recorrer los alineamientos y establecer el numero aminoacidos
                for alignSec in self.alignment:
                        id = alignSec.id
                        nU = alignSec.seq.count('U')
                        self.df.at[id, 'nU'] = nU

        def setnV(self):
                # Columna nueva y valor por defecto
                self.df['nV'] = -1
                # Recorrer los alineamientos y establecer el numero aminoacidos
                for alignSec in self.alignment:
                        id = alignSec.id
                        nV = alignSec.seq.count('V')
                        self.df.at[id, 'nV'] = nV

        def setnW(self):
                # Columna nueva y valor por defecto
                self.df['nW'] = -1
                # Recorrer los alineamientos y establecer el numero aminoacidos
                for alignSec in self.alignment:
                        id = alignSec.id
                        nW = alignSec.seq.count('W')
                        self.df.at[id, 'nW'] = nW

        def setnX(self):
                # Columna nueva y valor por defecto
                self.df['nX'] = -1
                # Recorrer los alineamientos y establecer el numero aminoacidos
                for alignSec in self.alignment:
                        id = alignSec.id
                        nX = alignSec.seq.count('X')
                        self.df.at[id, 'nX'] = nX

        def setnY(self):
                # Columna nueva y valor por defecto
                self.df['nY'] = -1
                # Recorrer los alineamientos y establecer el numero aminoacidos
                for alignSec in self.alignment:
                        id = alignSec.id
                        nY = alignSec.seq.count('Y')
                        self.df.at[id, 'nY'] = nY

        def setnZ(self):
                # Columna nueva y valor por defecto
                self.df['nZ'] = -1
                # Recorrer los alineamientos y establecer el numero aminoacidos
                for alignSec in self.alignment:
                        id = alignSec.id
                        nZ = alignSec.seq.count('Z')
                        self.df.at[id, 'nZ'] = nZ

        def setNGap(self):
                # Columna nueva y valor por defecto
                self.df['nGap'] = -1
                # Recorrer los alineamientos y establecer el numero de '-'
                for alignSec in self.alignment:
                        id = alignSec.id
                        nGap = alignSec.seq.count('-')
                        self.df.at[id, 'nGap'] = nGap

        def setNGaps(self):
                # Columna nueva y valor por defecto
                self.df['nGaps'] = -1
                # Recorrer los alineamientos y establecer el numero conjuntos de Gaps(grupos de - entre aminoacidos)
                for alignSec in self.alignment:
                        id = alignSec.id
                        nGaps = len(re.findall('-+', str(alignSec.seq).strip('-')))
                        self.df.at[id, 'nGaps'] = nGaps

        def setConsensusPercentage(self):
                # Columna nueva y valor por defecto
                self.df['consensus'] = -1.0
                consensusStr = str(self.consensus)
                n = len(consensusStr)
                for alignSec in self.alignment:
                        id = alignSec.id
                        seqStr = str(alignSec.seq)                       
                        x = 0
                        for i in range(0,n):
                                if (seqStr[i] == consensusStr[i]):
                                        x+=1
                        score = (x / n) * 100
                        self.df.at[id,'consensus'] = score

        


        

        
        def populate_features(self):
                self.setnA()
                self.setnB()
                self.setnC()
                self.setnD()
                self.setnE()
                self.setnF()
                self.setnG()
                self.setnH()
                self.setnI()
                self.setnJ()
                self.setnK()
                self.setnL()
                self.setnM()
                self.setnN()
                self.setnO()
                self.setnP()
                self.setnQ()
                self.setnR()
                self.setnS()
                self.setnT()
                self.setnU()
                self.setnV()
                self.setnX()
                self.setnY()
                self.setnW()
                self.setnZ()
                self.setNGap()
                self.setNGaps()
                self.setConsensusPercentage()
                print(self.df.to_string())
                
        
        


       
class AlignSequence():
        
        def __init__(self,seq):
                self.recSeq = seq
        


if __name__ == '__main__':
        
        util = Utils()
        # Cargar kinasas no patológicas
        util.load_kinases()
        print("")

        # Ejecucion del alineamiento
        print("Ejecucion del alineamiento") 
        aligmentTool = 'clustalw'
        als = util.run_alignment(aligmentTool)

        # Crear las features

        util.populate_features()

        # Secuencias consenso
        print("Secuencia consenso")
        print(util.consensus)


        # Sacar una secuencia del alineamiento
        print("pruebas con una secuencia alineada")
        a = als[5]
        print(a.id)
        print(a.seq)
        print(a.seq.count('A'))
        print(a.seq.count('B'))
        print(a.seq.count('C'))
        print(a.seq.count('D'))
        print(a.seq.count('E'))
        print(a.seq.count('F'))
        print(a.seq.count('G'))
        print(a.seq.count('U'))
        print(a.seq.count('V'))
        print(a.seq.count('W'))
        print(a.seq.count('X'))
        print(a.seq.count('Y'))
        print(a.seq.count('Z'))
        print(a.seq.count('-'))

        
                

        

