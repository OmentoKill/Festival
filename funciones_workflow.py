from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import gc_fraction
import re
import os
import matplotlib.pyplot as plt
from collections import Counter
import numpy as np


#encargado de la realización de api y descarga de archivos .fasta
def descargar_fasta(gen, especie, email, folder="./fastas", prefix="secuencia"):
    if not os.path.exists(folder):
        os.makedirs(folder)
    n = 1
    while True:
        filename = os.path.join(folder, f"{prefix}[{n}].fasta")
        if not os.path.exists(filename):
            break
        n += 1

    Entrez.email = email
    query = f"{gen}[Gene] AND {especie}[Organism]"
    handle = Entrez.esearch(db="nucleotide", term=query, retmax=1)
    record = Entrez.read(handle)
    handle.close()
    if not record['IdList']:
        print("No se encontró ningún resultado.")
        return None
    seq_id = record['IdList'][0]
    handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="fasta", retmode="text")
    with open(filename, "w") as out_f:
        out_f.write(handle.read())
    handle.close()
    print(f"Secuencia descargada en: {filename}")
    return filename

def limpiar_secuencia(seq):
    # Solo A, T, G, C
    return ''.join([b for b in seq.upper() if b in 'ATGC'])

def fusionar_fasta(fasta_files, output_fasta="secuencia_fusionada.fasta"):
    sec_total = ''
    ids = []
    longitudes = []
    for ff in fasta_files:
        record = SeqIO.read(ff, "fasta")
        sec = str(record.seq)
        sec_total += sec
        ids.append(record.id)
        longitudes.append(len(limpiar_secuencia(sec)))
    record_fusion = SeqRecord(Seq(sec_total), id="_".join(ids), description="Fusión de secuencias FASTA")
    SeqIO.write(record_fusion, output_fasta, "fasta")
    print(f"Secuencia fusionada guardada en: {output_fasta}")
    return sec_total, longitudes

def traducir_a_proteina(secuencia):
    """
    Traduce una secuencia de ADN a proteína, manejando correctamente secuencias
    cuya longitud no es un múltiplo de 3 para evitar advertencias.
    """
    resto = len(secuencia) % 3
    if resto != 0:
        secuencia_recortada = secuencia[:-resto]
    else:
        secuencia_recortada = secuencia

    return str(Seq(secuencia_recortada).translate(to_stop=True))

def contenido_gc(secuencia):
    return gc_fraction(secuencia) * 100

def secuencia_complementaria(secuencia):
    return str(Seq(secuencia).complement())

def buscar_patron(secuencia, patron):
    return [(m.start(), m.group()) for m in re.finditer(patron, secuencia)]

