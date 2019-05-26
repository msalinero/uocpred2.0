from sklearn.ensemble import RandomForestClassifier
from sklearn.tree import DecisionTreeClassifier, DecisionTreeRegressor
from numpy.random import permutation
from numpy import array_split, concatenate
from sklearn.metrics import mean_squared_error, roc_curve
import sys
import re
import os
import subprocess
import numpy as np
import pandas as pd
import Bio.Align.Applications as ba
from Bio import Phylo
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import AlignInfo
import tkinter as tk
from tkinter import filedialog
from tkinter.ttk import *
from tkinter import scrolledtext
from tkinter import *
import time
import configparser
from sklearn.metrics import accuracy_score, confusion_matrix
from sklearn.metrics import f1_score, balanced_accuracy_score
import matplotlib.pyplot as plt
import datetime


class Utils():
        
        def __init__(self, df=pd.DataFrame({}), consensus='',
         alignment='', alignment_tool='', summary_align='', 
                     fastaout='', clustalw_exe='', nfold=2, criterion='', outputPath='', currpath ='', tmstmp=''):
                
                # Directorio actual del script
                self.currpath = os.path.dirname(os.path.abspath(__file__))
                # Directorio para los archivos de salida
                self.outputPath = self.currpath + '\\OutputFiles\\'
                # TimeStamp para identificar los archivos de esta ejecucion
                self.tmstmp = datetime.datetime.now().strftime("%d%m%Y-%H%M%S")

                # Creamos el directorio de archivos de salida si no existe
                if not os.path.exists(self.outputPath):
                        os.makedirs(self.outputPath)

                # Leemos parámetros del archivo de configuracion
                try :
                        config = configparser.ConfigParser()
                        config.read(self.currpath + '\\config.file')
                        self.clustalw_exe = config['DEFAULT']['muscle_exe']
                        self.nfolds = int(config['ML']['kfolds'])
                        self.criterion = config['ML']['criterion']

                except:
                        print("Problema con el archivo de configuración")
        
        
        def select_fasta_file(self,message):
                # para quitar la ventana 
                root = tk.Tk()
                root.withdraw()
                
                # dialogo de seleccion de archivo
                file_path = filedialog.askopenfilename(
                        initialdir = "./",title = message,
                        filetypes = (("FASTA files","*.fasta"),("all files","*.*")))
                root.destroy()
                root.quit()
                return file_path

        def load_kinases(self):

                # Lista que va a contener todos los diccionarios
                rows_list = []

                # Archivo conteniendo las quinasas benignas
                print("Carga de las quinasas benignas\n")
                f1 = self.select_fasta_file(
                    "Selecciona el archivo FASTA con las quinasas benignas")
                # Leemos el fichero y creamos la lista de SeqRecord
                healthySeqRecords = SeqIO.parse(f1, 'fasta')     
                # Iteramos y vamos generando las entradas del 
                # diccionario de cada secuencia           
                numHealth = 0
                # Archivo que recoge todas las quinasas
                self.fastaout = self.outputPath + self.tmstmp + 'allkin.fasta'
                outfile = open(self.fastaout, 'w')
                for seq in healthySeqRecords:
                        numHealth += 1
                        dict1 = {}
                        dict1.update({'id': seq.id,
                                      'state': 'healthy',
                                      'length': len(seq),
                                      'gene': self.get_gene_from_description(seq.description)
                                      })
                        rows_list.append(dict1)
                        outfile.write(">" + seq.id + "\n")
                        outfile.write(str(seq.seq) + "\n")
                
                print("Cargadas {} quinasas benignas\n".format(numHealth))

                # Archivo conteniendo las quinasas patologicas
                print("Carga de las quinasas patologicas\n")
                f2 = self.select_fasta_file(
                    "Selecciona el archivo FASTA con las quinasas patológicas")
                # Leemos el fichero y creamos la lista de SeqRecord
                patSeqRecords = SeqIO.parse(f2, 'fasta')
                # Iteramos y vamos generando las entradas del
                # diccionario de cada secuencia
                numPats = 0
                for seq in patSeqRecords:
                        numPats +=1
                        dict1 = {}
                        dict1.update({'id': seq.id,
                                      'state': 'pathologic',
                                      'length': len(seq),
                                      'gene': self.get_gene_from_description(seq.description)
                                      })
                        rows_list.append(dict1)
                        outfile.write(">" + seq.id + "\n")
                        outfile.write(str(seq.seq) + "\n")

                print("Cargadas {} quinasas patológicas\n".format(numPats))   

                outfile.close()

                # Creamos el DataFrame con la totalidad de las filas
                # y las columnas hasta ahora creadas                
                self.df = pd.DataFrame(rows_list)
                self.df = self.df.set_index('id')

                print("Cargadas un total de {} quinasas\n".format(len(self.df.index)))

                # Volcamos todas las entradas de los dos archivos en un mismo
                # archivo para el alineamiento
                filenames = [f1, f2]
                
                print("Generando el archivo FASTA con la totalidad de las quinasas")
                #self.fastaout = self.outputPath + self.tmstmp + 'allkin.fasta'
                #with open(self.fastaout, 'w') as outfile:
                #        for seq in healthySeqRecords:
                #                outfile.write("Hello World")
                #                outfile.write(str(seq.id))
                #                outfile.write(str(seq.seq))
                #        for seq in patSeqRecords:
                #                outfile.write(str(seq.id))
                #                outfile.write(str(seq.seq))
                #with open(self.fastaout, 'w') as outfile:
                #        for fname in filenames:
                #                with open(fname) as infile:
                #                        outfile.write(infile.read())

        
        def exec_clustalW(self,f):
                # Tiempo de inicio
                start = time.time()
                # Linea de comando para ejecutar el alineamiento ClustalW
                clustalw_cline = ba.ClustalwCommandline(self.clustalw_exe, infile=f)
                # Ejecutamos la linea de ClustalW
                stdout, stderr = clustalw_cline()
                # Final de ejecucion
                end = time.time()
                timeCost = end - start
                print("El alineamiento con CLustalW ha tardado {0:.2f} segundos.\n".format(timeCost))
                
        def exec_Muscle(self,f):
                # Tiempo de inicio
                start = time.time()
                alignFile = f.replace('.fasta', '.aln')
                # Comando de ejecucion de muscle
                muscle_cline = ba.MuscleCommandline(input=f, out=alignFile)
                # Ejecucion del comando
                stdout, stderr = muscle_cline()
                treeFile = f.replace('.fasta', '.dnd')
                print("Arbol filogenetico")
                # Comando para generar el arbol filogenético
                clinetree = "muscle -maketree -in " + alignFile + " -out " + treeFile + " -cluster neighborjoining"
                # Ejecucion del comando
                subprocess.call(clinetree)
                # Final de ejecucion
                end = time.time()
                timeCost = end - start
                print("El alineamiento con Muscle ha tardado {0:.2f} segundos.\n".format(timeCost))

        def print_phylo_tree(self,f):
                tree = Phylo.read(f, "newick")
                Phylo.draw(tree)
                

        def get_alignment(self,f,format):
                alignment = AlignIO.read(f, format)
                self.alignment = alignment
                return alignment

        def draw_alignment(self, alignfile):
                alignWindow = tk.Tk()
                alignWindow.title("Alineamiento")
                # Crear el area de texto
                txt = scrolledtext.ScrolledText(alignWindow, width=100, height=40)
                txt.grid(column=0, row=0)
                with open(alignfile, 'r') as f:
                        txt.insert('1.0',f.read())
                # Crear el boton
                def clicked():
                        alignWindow.destroy()
                        alignWindow.quit()
                btn = Button(alignWindow, text="OK", command=clicked)
                btn.grid(column=0, row=2)
                alignWindow.mainloop()

         
        def select_tool(self):
                toolWindow = tk.Tk()
                toolWindow.title("Herramienta de alineamiento")
                toolWindow.geometry('250x200')
                # Crear el label
                lbl = Label(toolWindow, text = "Selección del motor de alineamiento")
                lbl.grid(column=0, row=0)
                # Crear el combobox
                combo = Combobox(toolWindow, state="readonly")
                combo['values'] = ("ClustalW", "Muscle")
                combo.current(0)
                combo.grid(column=0, row=1)
                # Crear el boton
                def clicked():
                        self.alignment_tool = combo.get()
                        toolWindow.destroy()
                        toolWindow.quit()

                btn = Button(toolWindow, text="OK", command=clicked)
                btn.grid(column=0, row=2)

                toolWindow.mainloop()
                

        def run_alignment(self):
                # Alineamiento
                fastaFile = self.fastaout
                self.select_tool()
                # Barra de progreso
                print("Ejecucion del alineamiento con " + self.alignment_tool + "\n")
                                              
                if (self.alignment_tool == "ClustalW"):
                        self.exec_clustalW(fastaFile)
                        alignFormat = 'clustal'
                elif (self.alignment_tool == "Muscle"):
                        self.exec_Muscle(fastaFile)
                        alignFormat = 'fasta'
                else: 
                        raise Exception('Herramienta de alineamiento no reconocida')
                
                alignFile = fastaFile.replace('.fasta','.aln')                

                alignment = self.get_alignment(alignFile, alignFormat)
                self.draw_alignment(alignFile)
                
                # Arbol filogenetico
                treeFile = fastaFile.replace('.fasta', '.dnd')
                self.print_phylo_tree(treeFile)
                
                

                # Consensus secuence  
                self.summary_align = AlignInfo.SummaryInfo(self.alignment)
                self.consensus = self.summary_align.dumb_consensus(threshold=0.7,require_multiple=1,ambiguous='?')

                return alignment
        
        
        ## Investigar
        def get_activation_loop(self):
                pass

        ######################################
        #                                              
        # Funciones generadoras de Features
        # 
        ######################################

        # Si el campo GN='' del archivo fasta tiene la informacion del GEN, fijamos el gen
        def get_gene_from_description(self, desc):
                geneField = re.findall('GN=\S+', desc)
                if (len(geneField) == 1):
                        gen = geneField[0].split('=')[1]
                        return gen
                else:
                        return 'unknown_gene'


        def setnA(self):
                # Columna nueva y valor por defecto
                self.df['nA'] = -1.0
                # Recorrer los alineamientos y fijar la frecuencia de determinado Aminoacido
                for alignSec in self.alignment:
                        id = alignSec.id
                        nA = alignSec.seq.count('A') / self.df.at[id, 'length'] 
                        self.df.at[id, 'nA'] = nA

        def setnC(self):
                # Columna nueva y valor por defecto
                self.df['nC'] = -1.0
                # Recorrer los alineamientos y fijar la frecuencia de determinado Aminoacido
                for alignSec in self.alignment:
                        id = alignSec.id 
                        nC = alignSec.seq.count('C') / self.df.at[id, 'length']
                        self.df.at[id, 'nC'] = nC

        def setnD(self):
                # Columna nueva y valor por defecto
                self.df['nD'] = -1.0
                # Recorrer los alineamientos y fijar la frecuencia de determinado Aminoacido
                for alignSec in self.alignment:
                        id = alignSec.id
                        nD = alignSec.seq.count('D') / self.df.at[id, 'length'] 
                        self.df.at[id, 'nD'] = nD

        def setnE(self):
                # Columna nueva y valor por defecto
                self.df['nE'] = -1.0
                # Recorrer los alineamientos y fijar la frecuencia de determinado Aminoacido
                for alignSec in self.alignment:
                        id = alignSec.id
                        nE = alignSec.seq.count('E') / self.df.at[id, 'length'] 
                        self.df.at[id, 'nE'] = nE


        def setnF(self):
                # Columna nueva y valor por defecto
                self.df['nF'] = -1.0
                # Recorrer los alineamientos y fijar la frecuencia de determinado Aminoacido
                for alignSec in self.alignment:
                        id = alignSec.id
                        nF = alignSec.seq.count('F') / self.df.at[id, 'length']
                        self.df.at[id, 'nF'] = nF

        def setnG(self):
                # Columna nueva y valor por defecto
                self.df['nG'] = -1.0
                # Recorrer los alineamientos y fijar la frecuencia de determinado Aminoacido
                for alignSec in self.alignment:
                        id = alignSec.id
                        nG = alignSec.seq.count('G') / self.df.at[id, 'length'] 
                        self.df.at[id, 'nG'] = nG

        def setnH(self):
                # Columna nueva y valor por defecto
                self.df['nH'] = -1.0
                # Recorrer los alineamientos y fijar la frecuencia de determinado Aminoacido
                for alignSec in self.alignment:
                        id = alignSec.id
                        nH = alignSec.seq.count('H') / self.df.at[id, 'length'] 
                        self.df.at[id, 'nH'] = nH

        def setnI(self):
                # Columna nueva y valor por defecto
                self.df['nI'] = -1.0
                # Recorrer los alineamientos y fijar la frecuencia de determinado Aminoacido
                for alignSec in self.alignment:
                        id = alignSec.id
                        nI = alignSec.seq.count('I') / self.df.at[id, 'length']
                        self.df.at[id, 'nI'] = nI

        def setnK(self):
                # Columna nueva y valor por defecto
                self.df['nK'] = -1.0
                # Recorrer los alineamientos y fijar la frecuencia de determinado Aminoacido
                for alignSec in self.alignment:
                        id = alignSec.id
                        nK = alignSec.seq.count('K') / self.df.at[id, 'length']
                        self.df.at[id, 'nK'] = nK

        def setnL(self):
                # Columna nueva y valor por defecto
                self.df['nL'] = -1.0
                # Recorrer los alineamientos y fijar la frecuencia de determinado Aminoacido
                for alignSec in self.alignment:
                        id = alignSec.id
                        nL = alignSec.seq.count('L') / self.df.at[id, 'length'] 
                        self.df.at[id, 'nL'] = nL

        def setnM(self):
                # Columna nueva y valor por defecto
                self.df['nM'] = -1.0
                # Recorrer los alineamientos y fijar la frecuencia de determinado Aminoacido
                for alignSec in self.alignment:
                        id = alignSec.id
                        nM = alignSec.seq.count('M') / self.df.at[id, 'length'] 
                        self.df.at[id, 'nM'] = nM
        
        def setnN(self):
                # Columna nueva y valor por defecto
                self.df['nN'] = -1.0
                # Recorrer los alineamientos y fijar la frecuencia de determinado Aminoacido
                for alignSec in self.alignment:
                        id = alignSec.id
                        nN = alignSec.seq.count('N') / self.df.at[id, 'length'] 
                        self.df.at[id, 'nN'] = nN

        def setnP(self):
                # Columna nueva y valor por defecto
                self.df['nP'] = -1.0
                # Recorrer los alineamientos y fijar la frecuencia de determinado Aminoacido
                for alignSec in self.alignment:
                        id = alignSec.id
                        nP = alignSec.seq.count('P') / self.df.at[id, 'length'] 
                        self.df.at[id, 'nP'] = nP

        def setnQ(self):
                # Columna nueva y valor por defecto
                self.df['nQ'] = -1.0
                # Recorrer los alineamientos y fijar la frecuencia de determinado Aminoacido
                for alignSec in self.alignment:
                        id = alignSec.id
                        nQ = alignSec.seq.count('Q') / self.df.at[id, 'length'] 
                        self.df.at[id, 'nQ'] = nQ

        def setnR(self):
                # Columna nueva y valor por defecto
                self.df['nR'] = -1.0
                # Recorrer los alineamientos y fijar la frecuencia de determinado Aminoacido
                for alignSec in self.alignment:
                        id = alignSec.id
                        nR = alignSec.seq.count('R') / self.df.at[id, 'length'] 
                        self.df.at[id, 'nR'] = nR

        def setnS(self):
                # Columna nueva y valor por defecto
                self.df['nS'] = -1.0
                # Recorrer los alineamientos y fijar la frecuencia de determinado Aminoacido
                for alignSec in self.alignment:
                        id = alignSec.id
                        nS = alignSec.seq.count('S') / self.df.at[id, 'length']
                        self.df.at[id, 'nS'] = nS

        def setnT(self):
                # Columna nueva y valor por defecto
                self.df['nT'] = -1.0
                # Recorrer los alineamientos y fijar la frecuencia de determinado Aminoacido
                for alignSec in self.alignment:
                        id = alignSec.id
                        nT = alignSec.seq.count('T') / self.df.at[id, 'length'] 
                        self.df.at[id, 'nT'] = nT

        def setnV(self):
                # Columna nueva y valor por defecto
                self.df['nV'] = -1.0
                # Recorrer los alineamientos y fijar la frecuencia de determinado Aminoacido
                for alignSec in self.alignment:
                        id = alignSec.id
                        nV = alignSec.seq.count('V') / self.df.at[id, 'length']
                        self.df.at[id, 'nV'] = nV

        def setnW(self):
                # Columna nueva y valor por defecto
                self.df['nW'] = -1.0
                # Recorrer los alineamientos y fijar la frecuencia de determinado Aminoacido
                for alignSec in self.alignment:
                        id = alignSec.id
                        nW = alignSec.seq.count('W') / self.df.at[id, 'length'] 
                        self.df.at[id, 'nW'] = nW

        def setnX(self):
                # Columna nueva y valor por defecto
                self.df['nX'] = -1.0
                # Recorrer los alineamientos y fijar la frecuencia de determinado Aminoacido
                for alignSec in self.alignment:
                        id = alignSec.id
                        nX = alignSec.seq.count('X') / self.df.at[id, 'length'] 
                        self.df.at[id, 'nX'] = nX

        def setnY(self):
                # Columna nueva y valor por defecto
                self.df['nY'] = -1.0
                # Recorrer los alineamientos y fijar la frecuencia de determinado Aminoacido
                for alignSec in self.alignment:
                        id = alignSec.id
                        nY = alignSec.seq.count('Y') / self.df.at[id, 'length'] 
                        self.df.at[id, 'nY'] = nY

        def setNGap(self):
                # Columna nueva y valor por defecto
                self.df['nGap'] = -1
                # Recorrer los alineamientos y establecer la proporcion de gaps
                for alignSec in self.alignment:
                        id = alignSec.id
                        nGap = alignSec.seq.count('-') / self.df.at[id, 'length'] 
                        self.df.at[id, 'nGap'] = nGap

        def setNGaps(self):
                # Columna nueva y valor por defecto
                self.df['nGaps'] = -1
                # Recorrer los alineamientos y establecer el numero de conjuntos de Gaps (grupos de - entre aminoacidos)
                for alignSec in self.alignment:
                        id = alignSec.id
                        nGaps = len(re.findall('-+', str(alignSec.seq).strip('-')))
                        self.df.at[id, 'nGaps'] = nGaps

        def setConsensusPercentage(self):
                # Columna nueva y valor por defecto
                self.df['consensus'] = -1.0
                consensusStr = str(self.consensus)
                n = len(consensusStr)
                # Porcentage de coincidencia de la secuencia alineada con la secuencia consenso
                for alignSec in self.alignment:
                        id = alignSec.id
                        seqStr = str(alignSec.seq)                       
                        x = 0
                        for i in range(0,n):
                                if (seqStr[i] == consensusStr[i]):
                                        x+=1
                        score = (x / n) * 100
                        self.df.at[id,'consensus'] = score

        def setHydrophobicity_Group1(self):
                # Columna nueva y valor por defecto
                self.df['hydroG1'] = -1.0
                # Grupo 1 AA segun hidrofobia R, K, E, D, Q, N
                for alignSec in self.alignment:
                        id = alignSec.id
                        seq = str(alignSec.seq)
                        n = len(re.findall('[RKEDQN]', seq)) / self.df.at[id, 'length']
                        self.df.at[id, 'hydroG1'] = n

        def setHydrophobicity_Group2(self):
                # Columna nueva y valor por defecto
                self.df['hydroG2'] = -1.0
                # Grupo 2 AA segun hidrofobia G, A, S, T, P, H, Y
                for alignSec in self.alignment:
                        id = alignSec.id
                        seq = str(alignSec.seq)
                        n = len(re.findall('[GASTPHY]', seq)
                                ) / self.df.at[id, 'length']
                        self.df.at[id, 'hydroG2'] = n

        def setHydrophobicity_Group3(self):
                # Columna nueva y valor por defecto
                self.df['hydroG3'] = -1.0
                # Grupo 3 AA segun hidrofobia C, L, V, I, M, F, W
                for alignSec in self.alignment:
                        id = alignSec.id
                        seq = str(alignSec.seq)
                        n = len(re.findall('[CLVIMFW]', seq)
                                ) / self.df.at[id, 'length']
                        self.df.at[id, 'hydroG3'] = n
                
        def setVanderWaals_Group1(self):
                # Columna nueva y valor por defecto
                self.df['vdwG1'] = -1.0
                # Grupo 1 AA segun van der Waals Volume G, A, S, T, P, D, C
                for alignSec in self.alignment:
                        id = alignSec.id
                        seq = str(alignSec.seq)
                        n = len(re.findall('[GASTPDC]', seq)
                                ) / self.df.at[id, 'length']
                        self.df.at[id, 'vdwG1'] = n

        def setVanderWaals_Group2(self):
                # Columna nueva y valor por defecto
                self.df['vdwG2'] = -1.0
                # Grupo 1 AA segun van der Waals Volume N, V, E, Q, I, L
                for alignSec in self.alignment:
                        id = alignSec.id
                        seq = str(alignSec.seq)
                        n = len(re.findall('[NVEQIL]', seq)
                                ) / self.df.at[id, 'length']
                        self.df.at[id, 'vdwG2'] = n

        def setVanderWaals_Group3(self):
                # Columna nueva y valor por defecto
                self.df['vdwG3'] = -1.0
                # Grupo 1 AA segun van der Waals Volume M, H, K, F, R, Y, W
                for alignSec in self.alignment:
                        id = alignSec.id
                        seq = str(alignSec.seq)
                        n = len(re.findall('[MHKFRYW]', seq)
                                ) / self.df.at[id, 'length']
                        self.df.at[id, 'vdwG3'] = n

        def setPolarity_Group1(self):
                # Columna nueva y valor por defecto
                self.df['PG1'] = -1.0
                # Grupo 1 AA segun Polarity L, I, F, W, C, M, V, Y
                for alignSec in self.alignment:
                        id = alignSec.id
                        seq = str(alignSec.seq)
                        n = len(re.findall('[LIFWCMVY]', seq)
                                ) / self.df.at[id, 'length']
                        self.df.at[id, 'PG1'] = n

        def setPolarity_Group2(self):
                # Columna nueva y valor por defecto
                self.df['PG2'] = -1.0
                # Grupo 1 AA segun Polarity P, A, T, G, S
                for alignSec in self.alignment:
                        id = alignSec.id
                        seq = str(alignSec.seq)
                        n = len(re.findall('[PATGS]', seq)
                                ) / self.df.at[id, 'length']
                        self.df.at[id, 'PG2'] = n

        def setPolarity_Group3(self):
                # Columna nueva y valor por defecto
                self.df['PG3'] = -1.0
                # Grupo 1 AA segun Polarity H, Q, R, K, N, E, D
                for alignSec in self.alignment:
                        id = alignSec.id
                        seq = str(alignSec.seq)
                        n = len(re.findall('[HQRKNED]', seq)
                                ) / self.df.at[id, 'length']
                        self.df.at[id, 'PG3'] = n

        def setPolarizability_Group1(self):
                # Columna nueva y valor por defecto
                self.df['PolG1'] = -1.0
                # Grupo 1 AA segun Polarizability G, A, S, D, T
                for alignSec in self.alignment:
                        id = alignSec.id
                        seq = str(alignSec.seq)
                        n = len(re.findall('[GASDT]', seq)
                                ) / self.df.at[id, 'length']
                        self.df.at[id, 'PolG1'] = n
        
        def setPolarizability_Group2(self):
                # Columna nueva y valor por defecto
                self.df['PolG2'] = -1.0
                # Grupo 2 AA segun Polarizability C, P, N, V, E, Q, I, L
                for alignSec in self.alignment:
                        id = alignSec.id
                        seq = str(alignSec.seq)
                        n = len(re.findall('[CPNVEQIL]', seq)
                                ) / self.df.at[id, 'length']
                        self.df.at[id, 'PolG2'] = n

        def setPolarizability_Group3(self):
                # Columna nueva y valor por defecto
                self.df['PolG3'] = -1.0
                # Grupo 3 AA segun Polarizability K, M, H, F, R, Y, W
                for alignSec in self.alignment:
                        id = alignSec.id
                        seq = str(alignSec.seq)
                        n = len(re.findall('[KMHFRYW]', seq)
                                ) / self.df.at[id, 'length']
                        self.df.at[id, 'PolG3'] = n

        def setCharge_Group1(self):
                # Columna nueva y valor por defecto
                self.df['ChargeG1'] = -1.0
                # Grupo 3 AA segun carga K, R
                for alignSec in self.alignment:
                        id = alignSec.id
                        seq = str(alignSec.seq)
                        n = len(re.findall('[KR]', seq)
                                ) / self.df.at[id, 'length']
                        self.df.at[id, 'ChargeG1'] = n

        def setCharge_Group2(self):
                # Columna nueva y valor por defecto
                self.df['ChargeG2'] = -1.0
                # Grupo 3 AA segun carga A, N, C, Q, G, H, I, L, M, F, P, S, T, W, Y, V
                for alignSec in self.alignment:
                        id = alignSec.id
                        seq = str(alignSec.seq)
                        n = len(re.findall('[ANCQGHILMFPSTWYV]', seq)
                                ) / self.df.at[id, 'length']
                        self.df.at[id, 'ChargeG2'] = n

        def setCharge_Group3(self):
                # Columna nueva y valor por defecto
                self.df['ChargeG3'] = -1.0
                # Grupo 3 AA segun carga D, E
                for alignSec in self.alignment:
                        id = alignSec.id
                        seq = str(alignSec.seq)
                        n = len(re.findall('[DE]', seq)
                                ) / self.df.at[id, 'length']
                        self.df.at[id, 'ChargeG3'] = n

        def setSecondaryStructure_Group1(self):
                # Columna nueva y valor por defecto
                self.df['SEG1'] = -1.0
                # Grupo 3 AA segun estructura secundaria E, A, L, M, Q, K, R, H
                for alignSec in self.alignment:
                        id = alignSec.id
                        seq = str(alignSec.seq)
                        n = len(re.findall('[EALMQKRH]', seq)
                                ) / self.df.at[id, 'length']
                        self.df.at[id, 'SEG1'] = n

        def setSecondaryStructure_Group2(self):
                # Columna nueva y valor por defecto
                self.df['SEG2'] = -1.0
                # Grupo 3 AA segun estructura secundaria V, I, Y, C, W, F, T
                for alignSec in self.alignment:
                        id = alignSec.id
                        seq = str(alignSec.seq)
                        n = len(re.findall('[VIYCWFT]', seq)
                                ) / self.df.at[id, 'length']
                        self.df.at[id, 'SEG2'] = n

        def setSecondaryStructure_Group3(self):
                # Columna nueva y valor por defecto
                self.df['SEG3'] = -1.0
                # Grupo 3 AA segun estructura secundaria D, E
                for alignSec in self.alignment:
                        id = alignSec.id
                        seq = str(alignSec.seq)
                        n = len(re.findall('[GNPSD]', seq)
                                ) / self.df.at[id, 'length']
                        self.df.at[id, 'SEG3'] = n

        def setSolventAccessibility_Group1(self):
                # Columna nueva y valor por defecto
                self.df['SaG1'] = -1.0
                # Grupo 3 AA segun estructura secundaria D, E
                for alignSec in self.alignment:
                        id = alignSec.id
                        seq = str(alignSec.seq)
                        n = len(re.findall('[ALFCGIVW]', seq)
                                ) / self.df.at[id, 'length']
                        self.df.at[id, 'SaG1'] = n

        def setSolventAccessibility_Group2(self):
                # Columna nueva y valor por defecto
                self.df['SaG2'] = -1.0
                # Grupo 3 AA segun estructura secundaria D, E
                for alignSec in self.alignment:
                        id = alignSec.id
                        seq = str(alignSec.seq)
                        n = len(re.findall('[RKQEND]', seq)
                                ) / self.df.at[id, 'length']
                        self.df.at[id, 'SaG2'] = n

        def setSolventAccessibility_Group3(self):
                # Columna nueva y valor por defecto
                self.df['SaG3'] = -1.0
                # Grupo 3 AA segun estructura secundaria D, E
                for alignSec in self.alignment:
                        id = alignSec.id
                        seq = str(alignSec.seq)
                        n = len(re.findall('[MSPTHY]', seq)
                                ) / self.df.at[id, 'length']
                        self.df.at[id, 'SaG3'] = n
             
        def populate_features(self):
                print("Generando las features\n")
                self.setnA()
                self.setnC()
                self.setnD()
                self.setnE()
                self.setnF()
                self.setnG()
                self.setnH()
                self.setnI()
                self.setnK()
                self.setnL()
                self.setnM()
                self.setnN()
                self.setnP()
                self.setnQ()
                self.setnR()
                self.setnS()
                self.setnT()
                self.setnV()
                self.setnX()
                self.setnY()
                self.setnW()
                self.setNGap()
                self.setNGaps()
                self.setConsensusPercentage()
                self.setHydrophobicity_Group1()
                self.setHydrophobicity_Group2()
                self.setHydrophobicity_Group3()
                self.setVanderWaals_Group1()
                self.setVanderWaals_Group2()
                self.setVanderWaals_Group3()
                self.setPolarity_Group1()
                self.setPolarity_Group2()
                self.setPolarity_Group3()
                self.setPolarizability_Group1()
                self.setPolarizability_Group2()
                self.setPolarizability_Group3()
                self.setCharge_Group1()
                self.setCharge_Group2()
                self.setCharge_Group3()
                self.setSecondaryStructure_Group1()
                self.setSecondaryStructure_Group2()
                self.setSecondaryStructure_Group3()
                self.setSolventAccessibility_Group1()
                self.setSolventAccessibility_Group2()
                self.setSolventAccessibility_Group3()
                excelFile = self.outputPath + self.tmstmp + "dataframe.xlsx"
                print("Generando el excel file con el dataframe\n")
                self.df.to_excel(excelFile)
                
                
       

class ProteinProblem(object):

        def __init__(self,data):

                try:
                        config = configparser.ConfigParser()
                        config.read(os.path.dirname(os.path.abspath(__file__)) + '\\config.file')
                        
                        self.kfolds = int(config['ML']['kfolds'])
                        self.criterion = config['ML']['criterion']

                except:
                        print("Problema con el archivo de configuración")

                self.data_frame = data

                # Factorizamos las columnas que no son numericas(sabemos que las numéricas no son categoricas)                 
                cols = self.data_frame.columns
                numeric_cols = self.data_frame._get_numeric_data().columns
                nc = list(set(cols) - set(numeric_cols))

                for k in nc:
                        self.data_frame[k], _ = pd.factorize(self.data_frame[k])
                
                # Categorias de nuestro target 'state'
                categories = sorted(pd.Categorical(self.data_frame['state']).categories)
                
                self.classes = np.array(categories)
                self.features = self.data_frame.columns[self.data_frame.columns != 'state']
        
        @staticmethod
        def __factorize(data):

                y, _ = pd.factorize(pd.Categorical(data['state']), sort=True)
                return y

        def train(self, X, Y):
    
                raise NotImplementedError
                
        def validation_data(self):
                folds = self.kfolds
                df = self.data_frame
                response = []
                # Comprobar que hay mas observaciones que folds
                assert len(df.index) > folds
                # Generar particiones
                perms = array_split(permutation(len(df.index)), folds)
                

                for i in range(folds):
                        train_idxs = list(range(folds))
                        train_idxs.pop(i)
                        train = []
                        for idx in train_idxs:
                                train.append(perms[idx])
                        train = concatenate(train)

                        test_idx = perms[i]
                
                        # Observaciones para training
                        training = df.iloc[train]
                        # Observaciones para test
                        test_data = df.iloc[test_idx]

                        # Classes del conjunto de training
                        y = self.__factorize(training)
                        # Entrenar el modelo
                        classifier = self.train(training[self.features], y)
                        # Predecimos para el conjunto test
                        predictions = classifier.predict(test_data[self.features])
                        # Resultados esperados para el conjunto test
                        expected = self.__factorize(test_data)
                        response.append([predictions, expected])
                return response


class ProteinClassifier(ProteinProblem):

        def validate(self):

                confusion_matrices = []
                ncfs = 0
                avgAcc = 0.0
                for test, training in self.validation_data():

                        # Pintar confusion matrix
                        cm = self.confusion_matrix(training, test)
                        print(cm)
                        print()
                        ncfs += 1
                        cm1 = confusion_matrix(training, test, labels=[0, 1])
                        # Calculamos la precision de esta confusion matrix
                        accuracy = accuracy_score(training,test)
                        print("Accuracy:{0:.2f}\n".format(accuracy))
                        # Balance accuracy
                        print("Balance accuracy score:{0:.2f}\n".format(
                                balanced_accuracy_score(training, test)))
                        # F-Score
                        print("F-score:{0:.2f}\n".format(
                            f1_score(training, test)))

                        #fpr, tpr, threshold = metrics.roc_curve(y_test, preds)

                        avgAcc += accuracy
                        
                        confusion_matrices.append(cm1)                        

                avgAcc = avgAcc / ncfs
                print('Precision media:{0:.2f}'.format(avgAcc))

                return confusion_matrices

        @staticmethod
        def confusion_matrix(train, test):
                return pd.crosstab(train, test, rownames=['actual'], colnames=['preds'])
        
                         



class ProteinForest(ProteinClassifier):


        def train(self, X, Y):

                classifier = RandomForestClassifier(criterion=self.criterion, n_jobs=2, n_estimators=100)
                classifier = classifier.fit(X, Y)
                return classifier


if __name__ == '__main__':
        
        util = Utils()
        # Cargar quinasas patológicas y patologicas
        util.load_kinases()
        
        # Ejecucion del alineamiento
        als = util.run_alignment()

        # Crear las features
        util.populate_features()

        

        #Pruebas ML
        rf = ProteinForest(util.df)
        cfs = rf.validate()
        
        
                

        

