from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.svm import SVC
from numpy.random import permutation
from numpy import array_split, concatenate
from sklearn.metrics import mean_squared_error, roc_curve, auc
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
import itertools
from scipy import interp
from sklearn.decomposition import PCA
import seaborn as sns

class Utils():
        
        def __init__(self):
                
                self.df = pd.DataFrame({})
                self.consensus = ''
                self.alignment = ''
                self.alignment_tool = ''
                self.summary_align = ''
                self.fastaout = '' 
                self.ncomp = 0

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
                        self.doPCA = config['DEFAULT']['PCA']
                        
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

        def load_proteins(self):

                # Lista que va a contener todos los diccionarios
                rows_list = []

                # Archivo conteniendo las proteinas no patologicas
                print("Carga de las proteinas no patologicas\n")
                f1 = self.select_fasta_file(
                    "Selecciona el archivo FASTA con las proteinas no patologicas")
                # Leemos el fichero y creamos la lista de SeqRecord
                healthySeqRecords = SeqIO.parse(f1, 'fasta')     
                # Iteramos y vamos generando las entradas del 
                # diccionario de cada secuencia           
                numHealth = 0
                # Archivo que recoge todas las proteinas
                self.fastaout = self.outputPath + self.tmstmp + 'AllProt.fasta'
                outfile = open(self.fastaout, 'w')
                for seq in healthySeqRecords:
                        numHealth += 1
                        dict1 = {}

                        # Generar features preliminares
                        dict1.update({'id': seq.id,
                                      'state': 'healthy',
                                      'length': len(seq),
                                      'gene': self.get_gene_from_description(seq.description)
                                      })
                        
                        # Generar features de proporcion de aminoacidos
                        self.generate_aminoacid_features(str(seq.seq), dict1)
                        
                        # Generar features de proporcion de dipeptidos
                        self.generate_dipeptides_features(str(seq.seq), dict1)
                        
                        rows_list.append(dict1)

                        outfile.write(">" + seq.id + "\n")
                        outfile.write(str(seq.seq) + "\n")
                
                print("Cargadas {} proteinas no patologicas\n".format(numHealth))

                # Archivo conteniendo las proteinas patologicas
                print("Carga de las proteinas patologicas\n")
                f2 = self.select_fasta_file(
                    "Selecciona el archivo FASTA con las proteinas patológicas")
                # Leemos el fichero y creamos la lista de SeqRecord
                patSeqRecords = SeqIO.parse(f2, 'fasta')
                # Iteramos y vamos generando las entradas del
                # diccionario de cada secuencia
                numPats = 0
                for seq in patSeqRecords:
                        numPats +=1
                        dict1 = {}
                        # Generar features preliminares
                        dict1.update({'id': seq.id,
                                      'state': 'pathologic',
                                      'length': len(seq),
                                      'gene': self.get_gene_from_description(seq.description)
                                      })

                        # Generar features de proporcion de aminoacidos
                        self.generate_aminoacid_features(str(seq.seq), dict1)

                        # Generar features de proporcion de dipeptidos
                        self.generate_dipeptides_features(str(seq.seq), dict1)

                        rows_list.append(dict1)

                        outfile.write(">" + seq.id + "\n")
                        outfile.write(str(seq.seq) + "\n")
                
                print("Cargadas {} proteinas patológicas\n".format(numPats))   

                outfile.close()               

                # Creamos el DataFrame siendo cada fila un diccionario 
                # y el valor de cada clave el valor de esa columna
                # determinada de cada 
                self.df = pd.DataFrame(rows_list)
                self.df = self.df.set_index('id')

                print("La totalidad de {} proteinas se han volcado en el archivo {}.\n".format(len(self.df.index),
                                                                                          self.fastaout))

               
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
                Phylo.draw(tree, do_show=True, show_confidence=True)
                

        def get_alignment(self,f,format):
                alignment = AlignIO.read(f, format)
                self.alignment = alignment
                return alignment

        def draw_alignment(self, alignfile):
                alignWindow = tk.Tk()
                alignWindow.title("Alineamiento")
                alignWindow.resizable(False, False)
                # Crear el area de texto
                txt = scrolledtext.ScrolledText(alignWindow, width=100, height=40)
                txt.pack()
                with open(alignfile, 'r') as f:
                        txt.insert('1.0',f.read())
                # Crear el boton
                def clicked():
                        alignWindow.destroy()
                        alignWindow.quit()
                btn = Button(alignWindow, text="OK", command=clicked)
                btn.pack()
                alignWindow.mainloop()

         
        def select_tool(self):
                toolWindow = tk.Tk()
                toolWindow.title("Herramienta de alineamiento")
                toolWindow.geometry('250x200')
                # Crear el label
                lbl = Label(toolWindow, text = "Selecciona el motor de alineamiento")
                lbl.pack()
                # Crear el combobox
                combo = Combobox(toolWindow, state="readonly")
                combo['values'] = ("ClustalW", "Muscle")
                combo.current(0)
                combo.pack()
                # Crear el boton
                def clicked():
                        self.alignment_tool = combo.get()
                        toolWindow.destroy()
                        toolWindow.quit()

                btn = Button(toolWindow, text="OK", command=clicked)
                btn.pack()

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
                print("Representacion del arbol filogenetico\n")
                treeFile = fastaFile.replace('.fasta', '.dnd')
                self.print_phylo_tree(treeFile)
                
                

                # Consensus secuence  
                self.summary_align = AlignInfo.SummaryInfo(self.alignment)
                self.consensus = self.summary_align.dumb_consensus(threshold=0.7,require_multiple=1,ambiguous='?')

                return alignment
        
                
        ######################################
        #                                              
        # Funciones generadoras de Features
        # 
        ######################################

        # Proporcion de aminoacidos
        def generate_aminoacid_features(self, secuence, dict):

                # Totalidad de aminoacidos
                amino_acids = 'ACDEFGHIKLMNPQRSTVWY'

                for a in amino_acids:
                        dict[a] = secuence.count(a) / len(secuence)

        # Proporcion de dipeptidos
        def generate_dipeptides_features(self, secuence, dict):

                # Totalidad de aminoacidos
                amino_acids = 'ACDEFGHIKLMNPQRSTVWY'

                dipeptides = itertools.product(amino_acids, repeat=2)
                for dipep in dipeptides:
                        d = "".join(dipep)
                        dict[d] = secuence.count(d) / (len(secuence)-1)


        # Si el campo GN='' del archivo fasta tiene la informacion del GEN, fijamos el gen
        def get_gene_from_description(self, desc):
                geneField = re.findall('GN=\S+', desc)
                if (len(geneField) == 1):
                        gen = geneField[0].split('=')[1]
                        return gen
                else:
                        return 'unknown_gene'
        
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
                # Generacion de features adicionales
                print("Generando las features\n")
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

                # Exportacion del dataframe y exportacion de EDA
                
                excelDf = self.outputPath + self.tmstmp + "Dataframe.xlsx"

                self.df.to_excel(excelDf)
                
                print("Dataframe exportado al archivo {}\n".format(excelDf))

                excelSummary = self.outputPath + self.tmstmp + "Summary.xlsx"

                self.df.describe(include='all').to_excel(excelSummary)
                
                print("Exploracion de los datos exportado al archivo {}\n".format(excelSummary))
                

                # PCA
                if(self.doPCA=='True'):
                        self.runPca()

                        

        def runPca(self):

                # Factorizamos las columnas que no son numericas(sabemos que las numéricas no son categoricas)
                cols = self.df.columns
                numeric_cols = self.df._get_numeric_data().columns
                nc = list(set(cols) - set(numeric_cols))

                for k in nc:
                        self.df[k], _ = pd.factorize(
                        self.df[k])
                

                dfAux = self.df.loc[:, self.df.columns != 'state']
                pca = PCA().fit(dfAux)
                plt.plot(np.cumsum(pca.explained_variance_ratio_))
                plt.xlabel('Número de componentes')
                plt.ylabel('Varianza explicada acumulada')
                plt.show(block=False)

                self.set_ncomponents()

                pca = PCA(n_components = self.ncomp)
                pca.fit(dfAux) 
                dfPCA = pd.DataFrame(pca.transform(dfAux), columns=[
                        'PCA%i' % i for i in range(self.ncomp)],index = dfAux.index)
                dfPCA = dfPCA.join(
                        self.df.loc[:, self.df.columns == 'state'], on='id')

                excelPCA = self.outputPath + self.tmstmp + "DfPCA.xlsx"
                dfPCA.to_excel(excelPCA)
                print("Dataframe con el PCA efectuado exportado al archivo {}\n".format(excelPCA))
                
                self.df = dfPCA                


        def set_ncomponents(self):
                toolWindow = tk.Tk()
                toolWindow.title("Numero de componentes")
                toolWindow.geometry('250x200')
                # Crear el label
                lbl = Label(
                    toolWindow, text="Introduce el numero de componentes principales")
                lbl.pack()
                # Crear el combobox
                ent = tk.Entry(toolWindow)
                ent.pack()
                # Crear el boton
                def clicked():
                        plt.close('all')
                        self.ncomp = int(ent.get())
                        toolWindow.destroy()
                        toolWindow.quit()
                        

                btn = Button(toolWindow, text="OK", command=clicked)
                btn.pack()

                toolWindow.mainloop()

               

class ProteinProblem(object):

        def __init__(self,data,tmstmp,outputPath):

                self.data_frame = data
                self.tmstmp = tmstmp
                self.outputPath = outputPath
                self.validationFile = self.outputPath + self.tmstmp + self.name +'Validation.out'

                try:
                        

                        config = configparser.ConfigParser()
                        config.read(os.path.dirname(os.path.abspath(__file__)) + '\\config.file')
                        # Numero de folds
                        self.kfolds = int(config['ML']['kfolds'])
                        

                except:
                        print("Problema con el archivo de configuración")


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
                
                #Roc Curve
                tprs = []
                aucs = []
                mean_fpr = np.linspace(0, 1, 100)
                

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

                        # Clases del conjunto de training
                        y = self.__factorize(training)
                        # Entrenar el modelo
                        classifier = self.train(training[self.features], y)                        
                        # Predecimos para el conjunto test
                        predictions = classifier.predict(test_data[self.features])
                        # Resultados esperados para el conjunto test
                        expected = self.__factorize(test_data)
                        #Roc Curve
                        probas_ = classifier.predict_proba(test_data[self.features])
                        fpr, tpr, thresholds = roc_curve(expected, probas_[:, 1])
                        tprs.append(interp(mean_fpr, fpr, tpr))
                        tprs[-1][0] = 0.0
                        roc_auc = auc(fpr, tpr)
                        aucs.append(roc_auc)
                        plt.plot(fpr, tpr, lw=1, alpha=0.3,
                        label='ROC fold %d (AUC = %0.2f)' % (i, roc_auc))

                        response.append([predictions, expected])

                plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r',
                         label='Chance', alpha=.8)


                mean_tpr = np.mean(tprs, axis=0)
                mean_tpr[-1] = 1.0
                mean_auc = auc(mean_fpr, mean_tpr)
                std_auc = np.std(aucs)
                plt.plot(mean_fpr, mean_tpr, color='b',
                        label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc),
                        lw=2, alpha=.8)

                std_tpr = np.std(tprs, axis=0)
                tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
                tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
                plt.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2,
                                label=r'$\pm$ 1 std. dev.')

        
                plt.xlim([-0.05, 1.05])
                plt.ylim([-0.05, 1.05])
                plt.xlabel('False Positive Rate')
                plt.ylabel('True Positive Rate')
                plt.title('Curvas ROC de {}-fold Cross Validation para {}'.format(self.kfolds, self.name))
                plt.legend(loc="lower right")

                plotFile = self.outputPath + self.tmstmp + self.name + "KCross.png"
                plt.savefig(plotFile, bbox_inches='tight')

                plt.show(block=False)

                return response


class ProteinClassifier(ProteinProblem):

        def validate(self):

                print("Ejecutando {}-folds Cross Validation para {}\n".format(self.kfolds,self.name))

                original = sys.stdout

                sys.stdout = open(self.validationFile, 'w')

                print("{}-folds Cross Validation para {}.\n".format(self.kfolds,self.name))

                confusion_matrices = []
                ncfs = 0
                avgAcc = 0.0
                avgs = []
                
                for test, training in self.validation_data():

                        print("Fold {}\n".format(ncfs))

                        # Pintar confusion matrix
                        cm = self.confusion_matrix(training, test)
                        print(cm)
                        print()
                        
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

                        avgAcc += accuracy
                        
                        confusion_matrices.append(cm1)

                        fold_name = "{} fold {}".format(self.name,ncfs)

                        avgs.append([fold_name,accuracy])

                        ncfs += 1

                avgAcc = avgAcc / ncfs
                print('Precision media:{0:.2f}'.format(avgAcc))

                avgs.append(["Media {}".format(self.name), avgAcc])

                sys.stdout = original

                self.draw_validation()

                return avgs

        def draw_validation(self):
                validatioWindow = tk.Tk()
                validatioWindow.title("Validation")
                validatioWindow.resizable(False, False)
                # Crear el area de texto
                txt = scrolledtext.ScrolledText(
                    validatioWindow, width=100, height=40)
                txt.pack()
                with open(self.validationFile, 'r') as f:
                        txt.insert('1.0', f.read())
                # Crear el boton
                def clicked():
                        plt.close('all')
                        validatioWindow.destroy()
                        validatioWindow.quit()
                        
                btn = Button(validatioWindow, text="OK", command=clicked)
                btn.pack()
                validatioWindow.mainloop()

        @staticmethod
        def confusion_matrix(train, test):
                return pd.crosstab(train, test, rownames=['actual'], colnames=['preds'])
        
                         
class ProteinForest(ProteinClassifier):

        def __init__(self,data,tmstmp,outputPath):

                self.name = 'RandomForest'

                super().__init__(data, tmstmp, outputPath)

                # Hiperparametros random forest
                config = configparser.ConfigParser()
                config.read(os.path.dirname(os.path.abspath(__file__)) + '\\config.file')
                self.criterion = config['RANDOMFOREST']['criterion']
                self.n_estimators = int(config['RANDOMFOREST']['n_estimators'])
                self.n_jobs = int(config['RANDOMFOREST']['n_jobs'])

                


        def train(self, X, Y):

                classifier = RandomForestClassifier(criterion=self.criterion,
                                                    n_jobs=self.n_jobs, n_estimators=self.n_estimators)
                classifier = classifier.fit(X, Y)
                return classifier


class ProteinSVM(ProteinClassifier):

        def __init__(self, data, tmstmp, outputPath):

                self.name = 'SVC'

                super().__init__(data, tmstmp, outputPath)
                # Hiperparametros SVC
                config = configparser.ConfigParser()
                config.read(os.path.dirname(os.path.abspath(__file__)) + '\\config.file')
                self.kernel = config['SVC']['kernel']

               

        def train(self, X, Y):

                classifier = SVC(gamma='auto', probability=True, kernel = self.kernel)
                classifier = classifier.fit(X, Y)
                return classifier


class ProteinAdaBoost(ProteinClassifier):

        def __init__(self, data, tmstmp, outputPath):

                self.name = 'AdaBoost'

                super().__init__(data, tmstmp, outputPath)

                # Hiperparametros AdaBoost
                config = configparser.ConfigParser()
                config.read(os.path.dirname(
                    os.path.abspath(__file__)) + '\\config.file')          

        def train(self, X, Y):

                classifier = AdaBoostClassifier()
                classifier = classifier.fit(X, Y)
                return classifier


class ProteinKNN(ProteinClassifier):

        def __init__(self, data, tmstmp, outputPath):

                self.name = 'KNN'

                super().__init__(data, tmstmp, outputPath)

                # Hiperparametros KNN
                config = configparser.ConfigParser()
                config.read(os.path.dirname(
                    os.path.abspath(__file__)) + '\\config.file')
                self.neighbors = int(config['KNN']['neighbors'])
                

        def train(self, X, Y):

                classifier = KNeighborsClassifier(n_neighbors=self.neighbors)
                classifier = classifier.fit(X, Y)
                return classifier


class Predictor(object):

        def __init__(self):

                # Hiperparametros AdaBoost
                config = configparser.ConfigParser()
                config.read(os.path.dirname(
                    os.path.abspath(__file__)) + '\\config.file')

                self.doSVC = config['SVC']['doSVC']
                self.doRandomForest = config['RANDOMFOREST']['doRF']
                self.doAdaBoost = config['ADABOOST']['doAda']
                self.doKNN = config['KNN']['doKNN']

       
        def pipeline(self):

                util = Utils()
                # Cargar proteinas no patologicas y patologicas
                util.load_proteins()

                # Ejecucion del alineamiento
                als = util.run_alignment()

                # Crear las features
                util.populate_features()

                # Algoritmos ML
                algoritmos_mls = []
                if (self.doRandomForest == 'True'):
                        algoritmos_mls.append(ProteinForest(util.df, util.tmstmp, util.outputPath))
                if (self.doSVC == 'True'):
                        algoritmos_mls.append(ProteinSVM(util.df, util.tmstmp, util.outputPath))
                if (self.doAdaBoost == 'True'):
                        algoritmos_mls.append(ProteinAdaBoost(util.df, util.tmstmp, util.outputPath))
                if (self.doKNN == 'True'):
                        algoritmos_mls.append(ProteinKNN(util.df, util.tmstmp, util.outputPath))


                # Validar los algoritmos
                log_cols = ["Clasificador", "Precision"]
                log = pd.DataFrame(columns=log_cols)
                

                for ml in algoritmos_mls:
                        for accs in ml.validate():
                                log_entry = pd.DataFrame([accs], columns=log_cols)
                                log = log.append(log_entry)
                
                # Tabla de resultados
                sns.set_color_codes("muted")
                clrs = ['r' if ('Media' in y) else 'b' for y in log["Clasificador"]]
                sns.barplot(x="Precision", y="Clasificador",
                            data=log, palette =clrs)

                plt.xlabel('Accuracy %')
                plt.title('Accuracy Clasificadores')
                plotFile = util.outputPath + util.tmstmp + "PerformanceComp.png"
                plt.savefig(plotFile, bbox_inches='tight')
                plt.show()

                #Exportar los resultados
                excelLog = util.outputPath + util.tmstmp + "AlgLog.xlsx"

                log = log.set_index('Clasificador')
                log.to_excel(excelLog)
   


                print("Fin de la ejecucion\n")






if __name__ == '__main__':
        
        pred = Predictor()
        pred.pipeline()

       
        
        
                

        

