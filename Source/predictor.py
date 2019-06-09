from sklearn.ensemble import RandomForestClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.svm import SVC
from sklearn.model_selection import GridSearchCV
from sklearn.preprocessing import MinMaxScaler
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
from sklearn.metrics import f1_score, precision_score, recall_score
import matplotlib.pyplot as plt
import datetime
import itertools
from scipy import interp
from sklearn.decomposition import PCA
import seaborn as sns
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

pd.options.mode.chained_assignment = None  # default='warn'

# Clase que recopila y prepara los datos
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
                        self.clustalw_exe = config['DEFAULT']['clustalw_exe']
                        self.doPCA = config['ML']['PCA']
                        
                except:
                        print("Problema con el archivo de configuración.")
        
        
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

                        # Generar Transition features
                        self.generate_P_transition_features(
                            str(seq.seq), dict1)
                        self.generate_Pol_transition_features(
                            str(seq.seq), dict1)
                        self.generate_SE_transition_features(
                            str(seq.seq), dict1)
                        self.generate_Sa_transition_features(
                            str(seq.seq), dict1)
                        self.generate_charge_transition_features(
                            str(seq.seq), dict1)
                        self.generate_hydro_transition_features(
                            str(seq.seq), dict1)
                        self.generate_vdw_transition_features(
                            str(seq.seq), dict1)
                        
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

                        # Generar Transition features
                        self.generate_P_transition_features(
                                str(seq.seq), dict1)
                        self.generate_Pol_transition_features(
                                str(seq.seq), dict1)
                        self.generate_SE_transition_features(
                                str(seq.seq), dict1)
                        self.generate_Sa_transition_features(
                                str(seq.seq), dict1)
                        self.generate_charge_transition_features(
                                str(seq.seq), dict1)
                        self.generate_hydro_transition_features(
                                str(seq.seq), dict1)
                        self.generate_vdw_transition_features(
                                str(seq.seq), dict1)

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

        # Representar el arbol filogenetico
        def print_phylo_tree(self,f):
                tree = Phylo.read(f, "newick")
                Phylo.draw(tree, do_show=True, show_confidence=True)
                

        def get_alignment(self,f,format):
                alignment = AlignIO.read(f, format)
                self.alignment = alignment
                return alignment

        # Ventana con los alineamientos desde el archivo generado
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

        # Dialogo de seleccion del alineamiento
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
                
        # Ejecutar el alineamiento seleccionado
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
        
        # Propocion de gaps en el alineamiento
        def setNGap(self):
                # Columna nueva y valor por defecto
                self.df['nGap'] = -1.0
                # Recorrer los alineamientos y establecer la proporcion de gaps
                for alignSec in self.alignment:
                        id = alignSec.id
                        nGap = alignSec.seq.count('-') / len(alignSec.seq) 
                        self.df.at[id, 'nGap'] = nGap

        # Grupos de gaps en la secuencia alineada
        def setNGaps(self):
                # Columna nueva y valor por defecto
                self.df['nGaps'] = -1
                # Recorrer los alineamientos y establecer el numero de conjuntos de Gaps (grupos de - entre aminoacidos)
                for alignSec in self.alignment:
                        id = alignSec.id
                        nGaps = len(re.findall('-+', str(alignSec.seq).strip('-')))
                        self.df.at[id, 'nGaps'] = nGaps

        # Proporcion de coincidencia con la secuencia consenso
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

        # Proporcion de aminoacidos de este grupo para la propiedad
        def setHydrophobicity_Group1(self):
                # Columna nueva y valor por defecto
                self.df['hydroG1'] = -1.0
                # Grupo 1 AA segun hidrofobia R, K, E, D, Q, N
                for alignSec in self.alignment:
                        id = alignSec.id
                        seq = str(alignSec.seq)
                        n = len(re.findall('[RKEDQN]', seq)) / self.df.at[id, 'length']
                        self.df.at[id, 'hydroG1'] = n

        # Proporcion de aminoacidos de este grupo para la propiedad
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

        # Proporcion de aminoacidos de este grupo para la propiedad
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

        # Transition features para esta propiedad
        def generate_hydro_transition_features(self, secuence, dict):
                groups = ['','[RKEDQN]','[GASTPHY]','[CLVIMFW]']

                for r,s in [[1,2],[1,3],[2,3]]:
                        grs = re.findall(groups[r]+groups[s], secuence)
                        gsr = re.findall(groups[s]+groups[r], secuence)
                        column_name = 'hydro.Tr' + str(r) + str(s) + str(s) + str(r)
                        dict[column_name] = (len(grs) + len(gsr))/ (len(secuence) - 1)

        # Proporcion de aminoacidos de este grupo para la propiedad
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

        # Proporcion de aminoacidos de este grupo para la propiedad
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

        # Proporcion de aminoacidos de este grupo para la propiedad
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

        # Transition features para esta propiedad
        def generate_vdw_transition_features(self, secuence, dict):
                groups = ['', '[GASTPDC]', '[NVEQIL]', '[MHKFRYW]']

                for r, s in [[1, 2], [1, 3], [2, 3]]:
                        grs = re.findall(groups[r]+groups[s], secuence)
                        gsr = re.findall(groups[s]+groups[r], secuence)
                        column_name = 'vdw.Tr' + \
                            str(r) + str(s) + str(s) + str(r)
                        dict[column_name] = (
                            len(grs) + len(gsr)) / (len(secuence) - 1)

        # Proporcion de aminoacidos de este grupo para la propiedad
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


        # Proporcion de aminoacidos de este grupo para la propiedad
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

        # Proporcion de aminoacidos de este grupo para la propiedad
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

        # Transition features para esta propiedad
        def generate_P_transition_features(self, secuence, dict):
                groups = ['', '[LIFWCMVY]', '[PATGS]', '[HQRKNED]']

                for r, s in [[1, 2], [1, 3], [2, 3]]:
                        grs = re.findall(groups[r]+groups[s], secuence)
                        gsr = re.findall(groups[s]+groups[r], secuence)
                        column_name = 'P.Tr' + \
                            str(r) + str(s) + str(s) + str(r)
                        dict[column_name] = (
                            len(grs) + len(gsr)) / (len(secuence) - 1)

        # Proporcion de aminoacidos de este grupo para la propiedad
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
        
        # Proporcion de aminoacidos de este grupo para la propiedad
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

        # Proporcion de aminoacidos de este grupo para la propiedad
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

        # Transition features para esta propiedad
        def generate_Pol_transition_features(self, secuence, dict):
                groups = ['', '[GASDT]', '[CPNVEQIL]', '[KMHFRYW]']

                for r, s in [[1, 2], [1, 3], [2, 3]]:
                        grs = re.findall(groups[r]+groups[s], secuence)
                        gsr = re.findall(groups[s]+groups[r], secuence)
                        column_name = 'Pol.Tr' + \
                            str(r) + str(s) + str(s) + str(r)
                        dict[column_name] = (
                            len(grs) + len(gsr)) / (len(secuence) - 1)

        # Proporcion de aminoacidos de este grupo para la propiedad
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

        # Proporcion de aminoacidos de este grupo para la propiedad
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

        # Proporcion de aminoacidos de este grupo para la propiedad
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

        # Transition features
        def generate_charge_transition_features(self, secuence, dict):
                groups = ['', '[KR]', '[ANCQGHILMFPSTWYV]', '[DE]']

                for r, s in [[1, 2], [1, 3], [2, 3]]:
                        grs = re.findall(groups[r]+groups[s], secuence)
                        gsr = re.findall(groups[s]+groups[r], secuence)
                        column_name = 'Charge.Tr' + \
                            str(r) + str(s) + str(s) + str(r)
                        dict[column_name] = (
                            len(grs) + len(gsr)) / (len(secuence) - 1)

        # Proporcion de aminoacidos de este grupo para la propiedad
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

        # Proporcion de aminoacidos de este grupo para la propiedad
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

        # Proporcion de aminoacidos de este grupo para la propiedad
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

        # Transition features para esta propiedad
        def generate_SE_transition_features(self, secuence, dict):
                groups = ['', '[EALMQKRH]', '[VIYCWFT]', '[GNPSD]']

                for r, s in [[1, 2], [1, 3], [2, 3]]:
                        grs = re.findall(groups[r]+groups[s], secuence)
                        gsr = re.findall(groups[s]+groups[r], secuence)
                        column_name = 'SE.Tr' + \
                            str(r) + str(s) + str(s) + str(r)
                        dict[column_name] = (
                            len(grs) + len(gsr)) / (len(secuence) - 1)

        # Proporcion de aminoacidos de este grupo para la propiedad
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

        # Proporcion de aminoacidos de este grupo para la propiedad
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

        # Proporcion de aminoacidos de este grupo para la propiedad
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

        # Transition features para esta propiedad
        def generate_Sa_transition_features(self, secuence, dict):
                groups = ['', '[ALFCGIVW]', '[RKQEND]', '[MSPTHY]']

                for r, s in [[1, 2], [1, 3], [2, 3]]:
                        grs = re.findall(groups[r]+groups[s], secuence)
                        gsr = re.findall(groups[s]+groups[r], secuence)
                        column_name = 'Sa.Tr' + \
                            str(r) + str(s) + str(s) + str(r)
                        dict[column_name] = (
                            len(grs) + len(gsr)) / (len(secuence) - 1)
             
        # Generar el resto de features y exportar los datos
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

                shape = self.df.shape
                ncols = shape[1]
                nsamples = shape[0]

                print("El dataframe generado tiene unas dimensiones de {} muestras y {} columnas\n".format(nsamples,ncols))
                
                print("Dataframe exportado al archivo {}\n".format(excelDf))

                excelSummary = self.outputPath + self.tmstmp + "Summary.xlsx"

                self.df.describe(include='all').to_excel(excelSummary)
                
                print("Exploracion de los datos exportado al archivo {}\n".format(excelSummary))
                

                

                        

        def runPca(self):

                # Factorizamos las columnas que no son numericas(sabemos que las numéricas no son categoricas)
                cols = self.df.columns
                numeric_cols = self.df._get_numeric_data().columns
                nc = list(set(cols) - set(numeric_cols))

                for k in nc:
                        self.df[k], _ = pd.factorize(
                        self.df[k])
                

                dfAux = self.df.loc[:, self.df.columns != 'state']

                scaler = MinMaxScaler()
                dfAux[dfAux.columns] = scaler.fit_transform(dfAux[dfAux.columns])

                excelScale = self.outputPath + self.tmstmp + "DfScale.xlsx"
                dfAux.to_excel(excelScale)

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

        # Dialogo para fijar el numero de componentes principales
        def set_ncomponents(self):
                toolWindow = tk.Tk()
                toolWindow.title("Número de componentes")
                toolWindow.geometry('250x200')
                # Crear el label
                lbl = Label(
                    toolWindow, text="Número de componentes principales")
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

        # Pipeline para generar los datos para los algoritmos de ML
        def prepareData(self):

                 # Cargar proteinas no patologicas y patologicas
                self.load_proteins()

                # Ejecucion del alineamiento
                self.run_alignment()

                # Crear las features
                self.populate_features()

                # PCA
                if(self.doPCA == 'True'):
                        self.runPca()
               
# Clase padre de los clasificadores
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
                        self.target_metric = config['ML']['target_metric']
                        

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
        
        # Factorizar el label objetivo (0=healthy, 1=pathologic)
        @staticmethod
        def __factorize(data):

                y, _ = pd.factorize(pd.Categorical(data['state']), sort=True)
                return y

        def train(self, X, Y):
    
                raise NotImplementedError
                
        # K-fold Cross Validation y generar curvas ROC
        def validation_data(self):
                folds = self.kfolds
                df = self.data_frame
                response = []
                # Comprobar que hay mas observaciones que folds
                assert len(df.index) > folds
                # Generar particiones
                perms = array_split(permutation(len(df.index)), folds)
                
                # Roc Curve
                tprs = []
                aucs = []
                mean_fpr = np.linspace(0, 1, 100)

                # Hiperparámetros usados
                dfGrid = pd.DataFrame({})
                rows_list = []
                         

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
                        # Imprimimos los parametros elegidos por grid search
                        print("Parametros elegidos para fold {}\n".format(i))
                        print(classifier.best_params_ )
                        print("")
                        dict1 = {}
                        dict1.update(classifier.best_params_)
                        rows_list.append(dict1)
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

                # Grafica que recoge las curvas ROC
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

                plotFile = self.outputPath + self.tmstmp + self.name + "RocCuves.png"
                plt.savefig(plotFile, bbox_inches='tight')

                plt.show(block=False)

                # Exportacion de grid search
                dfGrid = pd.DataFrame(rows_list)
                gridFile = self.outputPath + self.tmstmp + self.name + "GridSearch.xlsx"
                dfGrid.to_excel(gridFile)

                return response

# Clase intermedia
class ProteinClassifier(ProteinProblem):

        def validate(self):

                print("Ejecutando {}-folds Cross Validation para {}\n".format(self.kfolds,self.name))

                original = sys.stdout

                sys.stdout = open(self.validationFile, 'w')

                print("{}-folds Cross Validation para {}.\n".format(self.kfolds,self.name))

                confusion_matrices = []
                ncfs = 0
                mean_accuracy = 0.0
                mean_precision = 0.0
                mean_recall = 0.0
                mean_f1 = 0.0
                metrics = []
                
                for test, training in self.validation_data():

                        print("Fold {}\n".format(ncfs))

                        # Pintar confusion matrix
                        cm = self.confusion_matrix(training, test)
                        print(cm)
                        print()
                        
                        cm1 = confusion_matrix(training, test, labels=[0, 1])

                        # Calculamos la precision de esta confusion matrix
                        accuracy = accuracy_score(training,test)
                        mean_accuracy = mean_accuracy + accuracy
                        print("Accuracy:{0:.2f}\n".format(accuracy))
                        # Precision
                        precision = precision_score(training, test)
                        mean_precision = mean_precision + precision
                        print("Precision score:{0:.2f}\n".format(precision))
                        # Recall Score
                        recall = recall_score(training, test)
                        mean_recall = mean_recall + recall
                        print("Recall score:{0:.2f}\n".format(recall))
                        # F-Score
                        f1 = f1_score(training, test)
                        mean_f1 = mean_f1 + f1
                        print("F-score:{0:.2f}\n".format(f1))

                        confusion_matrices.append(cm1)

                        fold_name = "{} fold {}".format(self.name,ncfs)

                        metrics.append(
                            [fold_name, accuracy, precision, recall, f1])

                        ncfs += 1

                mean_accuracy = mean_accuracy / ncfs
                print('Accuracy media:{0:.2f}\n'.format(mean_accuracy))

                mean_precision = mean_precision / ncfs
                print('Precision media:{0:.2f}\n'.format(mean_precision))

                mean_recall = mean_recall / ncfs
                print('Recall media:{0:.2f}\n'.format(mean_recall))

                mean_f1 = mean_f1 / ncfs
                print('F1 media:{0:2f}\n'.format(mean_f1))

                metrics.append(["Media {}".format(self.name), mean_accuracy,mean_precision,mean_recall,mean_f1])

                sys.stdout = original

                self.draw_validation()

                return metrics


        # Pinta los hiperparametros, confusion matrix y metricas
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
        
# Implementacion de la clase random forest
class ProteinForest(ProteinClassifier):

        def __init__(self,data,tmstmp,outputPath):

                self.name = 'RandomForest'

                super().__init__(data, tmstmp, outputPath)

                # Hiperparametros random forest
                config = configparser.ConfigParser()
                config.read(os.path.dirname(os.path.abspath(__file__)) + '\\config.file')
                self.n_estimators = eval(config['RANDOMFOREST']['n_estimators'])
                self.max_features = eval(config['RANDOMFOREST']['max_features'])
                self.max_depth = eval(config['RANDOMFOREST']['max_depth'])
                self.criterion = eval(config['RANDOMFOREST']['criterion'])
                
                


        def train(self, X, Y):

                # Hacemos grid search para determinar los mejores parámetros
                param_grid = { 
                        'n_estimators': self.n_estimators,
                        'max_features': self.max_features,
                        'max_depth' : self.max_depth,
                        'criterion' : self.criterion
                }

                classifier = RandomForestClassifier()
                CV_rfc = GridSearchCV(
                    estimator=classifier, param_grid=param_grid, cv=3, scoring = self.target_metric)
                
                classifier = CV_rfc.fit(X, Y)
                return classifier

# Clase que implementa support vector machine
class ProteinSVM(ProteinClassifier):

        def __init__(self, data, tmstmp, outputPath):

                # Hiperparametros SVC
                config = configparser.ConfigParser()
                config.read(os.path.dirname(
                    os.path.abspath(__file__)) + '\\config.file')
                self.kernel = eval(config['SVC']['kernel'])
                self.C = eval(config['SVC']['C'])
                self.gamma = eval(config['SVC']['gamma'])
                self.decision_function_shape = eval(
                    config['SVC']['decision_function_shape'])
                self.shrinking = eval(config['SVC']['shrinking'])


                self.name = 'SVM'

                super().__init__(data, tmstmp, outputPath)
                

                

               
        # Crea y entrena para los datos proporcionados
        def train(self, X, Y):

                # Hacemos grid search para determinar los mejores parámetros
                param_grid = {
                    'kernel': self.kernel,
                    'C': self.C,
                    'gamma': self.gamma,
                    'decision_function_shape': self.decision_function_shape,
                    'shrinking': self.shrinking,
                }

                classifier = SVC(probability = True)

                CV_svm = GridSearchCV(
                    estimator=classifier, param_grid=param_grid, cv=3, scoring= self.target_metric)
                classifier = CV_svm.fit(X, Y)
                return classifier

# Clase que implementa decision tree
class ProteinTree(ProteinClassifier):

        def __init__(self, data, tmstmp, outputPath):

                self.name = 'DecTree'

                super().__init__(data, tmstmp, outputPath)

                # Hiperparametros Tree
                config = configparser.ConfigParser()
                config.read(os.path.dirname(
                    os.path.abspath(__file__)) + '\\config.file')
                self.criterion = eval(config['TREE']['criterion'])
                self.min_samples_leaf = eval(config['TREE']['min_samples_leaf'])
                self.max_depth = eval(config['TREE']['max_depth'])

        def train(self, X, Y):

                param_grid = {
                    'criterion': self.criterion,
                    'min_samples_leaf': self.min_samples_leaf,
                    'max_depth': self.max_depth,
                }

                classifier = DecisionTreeClassifier()
                CV_tree = GridSearchCV(
                    estimator=classifier, param_grid=param_grid, cv=3, scoring = self.target_metric)
                classifier = CV_tree.fit(X, Y)
                return classifier

# Clase que implementa KNN
class ProteinKNN(ProteinClassifier):

        def __init__(self, data, tmstmp, outputPath):

                # Hiperparametros KNN
                config = configparser.ConfigParser()
                config.read(os.path.dirname(
                    os.path.abspath(__file__)) + '\\config.file')
                self.n_neighbors = eval(config['KNN']['n_neighbors'])

                self.name = 'KNN' 

                super().__init__(data, tmstmp, outputPath)

                
                

        def train(self, X, Y):

                param_grid = {
                    'n_neighbors': self.n_neighbors
                }

                classifier = KNeighborsClassifier()
                CV_KNN = GridSearchCV(
                    estimator=classifier, param_grid=param_grid, cv=3, scoring = self.target_metric)
                classifier = CV_KNN.fit(X, Y)
                return classifier

# Clase que implementa el pipeline principal
class Predictor(object):

        def __init__(self):

                # Hiperparametros
                config = configparser.ConfigParser()
                config.read(os.path.dirname(
                    os.path.abspath(__file__)) + '\\config.file')

                self.doSVC = config['SVC']['doSVC']
                self.doRandomForest = config['RANDOMFOREST']['doRF']
                self.doTree = config['TREE']['doTree']
                self.doKNN = config['KNN']['doKNN']
                self.target_metric = config['ML']['target_metric']

        # Ejecuta los clasificadores y exporta los resultados
        def runClassifiers(self, ut):

                # Algoritmos ML
                algoritmos_mls = []
                if (self.doRandomForest == 'True'):
                        algoritmos_mls.append(ProteinForest(
                            ut.df, ut.tmstmp, ut.outputPath))
                if (self.doSVC == 'True'):
                        algoritmos_mls.append(ProteinSVM(
                            ut.df, ut.tmstmp, ut.outputPath))
                if (self.doKNN == 'True'):
                        algoritmos_mls.append(ProteinKNN(
                            ut.df, ut.tmstmp, ut.outputPath))
                if (self.doTree == 'True'):
                        algoritmos_mls.append(ProteinTree(
                            ut.df, ut.tmstmp, ut.outputPath))

                # Validar los algoritmos
                log_cols = ["Clasificador", "accuracy",
                            "precision", "recall", "f1"]
                log = pd.DataFrame(columns=log_cols)

                for ml in algoritmos_mls:
                        for metrics in ml.validate():
                                log_entry = pd.DataFrame(
                                    [metrics], columns=log_cols)
                                log = log.append(log_entry)

                # Tabla de resultados
                sns.set_color_codes("muted")
                clrs = ['g' if ('Media' in y)
                        else 'b' for y in log["Clasificador"]]
                sns.barplot(x=self.target_metric, y="Clasificador",
                            data=log, palette=clrs)

                plt.xlabel(self.target_metric)
                plt.title('{} Clasificadores'.format(self.target_metric))
                plotFile = ut.outputPath + ut.tmstmp + "PerformanceComp.png"
                plt.savefig(plotFile, bbox_inches='tight')
                plt.show()

                #Exportar los resultados
                excelLog = ut.outputPath + ut.tmstmp + "AlgLog.xlsx"

                log = log.set_index('Clasificador')
                log.to_excel(excelLog)



       
        def pipeline(self):

                util = Utils()

                # Preparar los datos
                util.prepareData()

                # Ejecutar clasificadores
                self.runClassifiers(util)

                print("Fin de la ejecución\n")






if __name__ == '__main__':
        
        pred = Predictor()
        pred.pipeline()

       
        
        
                

        

