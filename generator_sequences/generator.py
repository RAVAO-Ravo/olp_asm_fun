#!/bin/python3
#-*- coding:utf-8 -*-

import random
from argparse import ArgumentParser, MetavarTypeHelpFormatter
from typing import List

def generate_genome(length_genome: int) -> str:
    """
    Génère une séquence d'ADN aléatoire de la longueur spécifiée.

    :param length_genome: La longueur de la séquence d'ADN générée.
    :return: La séquence d'ADN générée.
    """
    # Définition des bases nucléotidiques à utiliser 
    bases = ['A', 'C', 'G', 'T']
    
    # Générer le génome
    return ''.join([random.choice(bases) for _ in range(length_genome)])

def generate_reads(genome: str, n_reads: int, min_read_length: int, max_read_length: int) -> List[str]:
    """
    Génère une liste de reads à partir d'une séquence d'ADN donnée.

    :param genome: La séquence d'ADN à partir de laquelle générer les reads.
    :param n_reads: Le nombre de reads à générer.
    :param min_read_length: La longueur minimale d'un read.
    :param max_read_length: La longueur maximale d'un read.
    :return: Une liste de reads générés.
    """
    # Initialiser la liste des reads
    reads = ['']*n_reads

    # Pour le nombre de reads indiqué
    for i in range(n_reads):

        # Sélectionner une longueur de read entre les deux bornes
        read_length = random.randint(min_read_length, max_read_length)
        
        # Sélectionner un point de départ dans le génome 
        start_pos = random.randint(0, len(genome) - read_length)
        
        # Générer le read
        read = genome[start_pos : start_pos + read_length]
        
        # Ajouter le read à la liste de reads
        reads[i] = read

    return reads

def save_genome(genome: str, filename: str) -> None:
    """
    Enregistre la séquence d'ADN générée au format FASTA.

    :param genome: La séquence d'ADN à enregistrer.
    :param filename: Le nom du fichier de sortie.
    """
    # Ouvrir le fichier qui contiendra le génome
    with open(file=filename, mode='w') as fasta:
        fasta.write(f">genome\n{genome}\n")

def save_reads(reads: List[str], filename: str) -> None:
    """
    Enregistre la liste de reads générés au format FASTQ.

    :param reads: La liste de reads à enregistrer.
    :param filename: Le nom du fichier de sortie.
    """
    # Ouvrir le fichier qui contiendra les reads
    with open(file=filename, mode='w') as fastq:
        for i, read in enumerate(iterable=reads, start=1):
            fastq.write(f"@Read_{i}\n{read}\n+\n{'~' * len(read)}\n")

if __name__ == "__main__":
    # Parseur d'arguments
    parser = ArgumentParser(description="Génère une séquence d'ADN aléatoire et des lectures.", 
                            formatter_class=lambda prog : MetavarTypeHelpFormatter(prog=prog, max_help_position=50, width=500))
    parser.add_argument("-G", "--length_genome", type=int, default=300, help="Longueur de la séquence d'ADN.")
    parser.add_argument("-m", "--min_read_length", type=int, default=40, help="Longueur minimale des lectures générées.")
    parser.add_argument("-M", "--max_read_length", type=int, default=80, help="Longueur maximale des lectures générées.")
    parser.add_argument("-n", "--n_reads", type=int, default=10000, help="Nombre de lectures à générer.")
    parser.add_argument("-f", "--fasta", type=str, default="genome.fasta", help="Nom du fichier fasta.")
    parser.add_argument("-q", "--fastq", type=str, default="reads.fastq", help="Nom du fichier fastq.")
    args = parser.parse_args()

    # Paramètres
    length_genome = args.length_genome
    min_read_length = args.min_read_length
    max_read_length = args.max_read_length
    n_reads = args.n_reads
    fasta = args.fasta
    fastq = args.fastq

    # Définir la seed
    random.seed(42)

    # Générer la séquence d'ADN
    genome = generate_genome(length_genome=length_genome)

    # Générer les reads à partir de la séquence d'ADN
    reads = generate_reads(genome=genome, n_reads=n_reads, min_read_length=min_read_length, max_read_length=max_read_length)

    # Sauvegarder la séquence d'ADN au format FASTA
    save_genome(genome=genome, filename=fasta)

    # Sauvegarder les reads au format FASTQ
    save_reads(reads=reads, filename=fastq)