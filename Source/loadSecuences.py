from Bio.Align.Applications import ClustalwCommandline
from Bio import SeqIO
from Bio import AlignIO
import Bio.Align.Applications
import tkinter as tk
from tkinter import filedialog

def get_file():
    root = tk.Tk()
    root.withdraw()
    
    file_path = filedialog.askopenfilename()
    return file_path


#def get_secuences():
#    file_path = get_file()
#    return SeqIO.parse(file_path, "fasta")

#def print_secuences(secs):
#    for seq_record in secs:
#        print(seq_record.id)
#        print(seq_record.seq)
#        print(len(seq_record))


#print("Secuencias e isoformas")
#benignSecuences = get_secuences()
#print_secuences(benignSecuences)

#print("Secuencias con GKD")
#GKDSecuences = get_secuences()
#print_secuences(GKDSecuences)

 
### Alineamiento con ClustalW

#help(ClustalwCommandline)

f = get_file()
cline = ClustalwCommandline("clustalw2", infile=f)
clustalw_exe = r"C:\Program Files (x86)\ClustalW2\clustalw2.exe"
clustalw_cline = ClustalwCommandline(clustalw_exe, infile=f)
stdout, stderr = clustalw_cline()
end = time.time()
print(end - start)
align = AlignIO.read("allhumankinases2.aln", "clustal")
print(align)
