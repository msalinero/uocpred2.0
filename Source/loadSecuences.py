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

                return alignment
        
        
        def load_kinases(self):
                f1 = self.select_file()
                healthySeqRecords = SeqIO.parse(f1,'fasta')
                rows_list =[]
                for seq in healthySeqRecords:
                        dict1 = {}
                        dict1.update({'id' : seq.id,
                                      'state' : 'healthy',
                                      'length': len(seq)
                        })  
                        rows_list.append(dict1)
                f2 = self.select_file()
                patSeqRecords = SeqIO.parse(f2, 'fasta')
                for seq in patSeqRecords:
                        dict1 = {}
                        dict1.update({'id' : seq.id, 
                                      'state': 'pathologic',
                                      'length': len(seq)
                        })
                        rows_list.append(dict1)
                self.df = pd.DataFrame(rows_list)
                self.df = self.df.set_index('id',drop = False)
                print(self.df)
                filenames = [f1,f2]
                fastaout = os.path.dirname(os.path.abspath(__file__)) + '\\allkin.fasta'
                with open(fastaout, 'w') as outfile:
                        for fname in filenames:
                                with open(fname) as infile:
                                        outfile.write(infile.read())

        def populate_features(self):
                pass
        
        def setnA(self):
                pass


       
class AlignSequence():
        
        def __init__(self,seq):
                self.recSeq = seq
        


if __name__ == '__main__':
        
        util = Utils()
        # Cargar kinasas no patológicas
        util.load_kinases()
        print("")

        # Ejecucion con ClustalW
         
        aligmentTool = 'clustalw'
        als = util.run_alignment(aligmentTool)

        # Sacar una secuencia del alineamiento
        a = als[0]
        aminoacids = set(a.seq)
        for aminoacid in aminoacids:
                print(aminoacid)
        print(len(aminoacids))
        print(a.seq.count('T'))
        print(len(a))
        print(len(a.seq))
         
        # Consensus secuence  
        summary_align = AlignInfo.SummaryInfo(als)
        consensus = summary_align.dumb_consensus(threshold=0.6,require_multiple=1,ambiguous='?')
        print(consensus)
        print(len(consensus))
        print((consensus.count('?') / len(consensus)) * 100)


        

