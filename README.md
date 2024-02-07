# olp_asm_fun

## Description

olp_asm_fun est un projet éducatif visant à comprendre les processus d'assemblage de génomes. Ce projet n'est pas destiné à un usage sérieux. Écrit en C++, il sert de programme d'assemblage utilisant un graphe de chevauchement et des *k*-mers si spécifiés.

## Structure du dépôt 

La structure du dépôt comprend trois répertoires principaux :

- `generator_sequences`: Contient un générateur de génome et de lectures codé en Python3.
- `include`: Contient des fichiers d'en-tête (.hpp).
- `src`: Contient des fichiers source (.cpp).

Le projet inclut également un Makefile pour générer l'exécutable "olp_asm" et des fichiers de test.

## Installation et dépendances

Le projet nécessite un compilateur C++ et un interpréteur Python3. Pour commencer, clonez le dépôt :

```bash
git clone https://github.com/RAVAO-Ravo/olp_asm_fun.git
```

Une fois le dépot récupéré, allez dans le répertoire, et exécutez la commande `make`.

## Utilisation

olp_asm_fun propose un ensemble d'options en ligne de commande pour assembler les séquences :

```bash
olp_asm -q <nom_fichier_fastq> [-k <longueur_kmers>] [-s <seuil>] -f <nom_fichier_fasta> [-m <longueur_minimale>] [-h]
```

- `-q <nom_fichier_fastq>` : Spécifie le nom du fichier FastQ à utiliser.
- `-k <longueur_kmers>` : Définit la longueur des *k*-mers à utiliser. (Optionnel, valeur par défaut : -1)
- `-s <seuil>` : Définit le score de chevauchement minimum pour conserver un nœud dans le graphe. (Optionnel, valeur par défaut : 10)
- `-f <nom_fichier_fasta>` : Spécifie le nom du fichier Fasta pour stocker les contigs.
- `-m <longueur_minimale>` : Définit la longueur minimale d'un contig à conserver. (Optionnel, valeur par défaut : 0)
- `-h` : Affiche ce message d'aide. (Optionnel)

Le projet propose également un générateur de séquences Python3 avec les options suivantes :

```bash
python3 generate_sequences.py -G <longueur_genome> [-m <longueur_min_lecture>] [-M <longueur_max_lecture>] [-n <nb_lectures>] [-f <nom_fichier_fasta>] [-q <nom_fichier_fastq>]
```

- `-G <longueur_genome>` : Définit la longueur de la séquence d'ADN à générer.
- `-m <longueur_min_lecture>` : Définit la longueur minimale des lectures générées. (Optionnel, valeur par défaut : 40)
- `-M <longueur_max_lecture>` : Définit la longueur maximale des lectures générées. (Optionnel, valeur par défaut : 80)
- `-n <nb_lectures>` : Définit le nombre de lectures à générer. (Optionnel, valeur par défaut : 10000)
- `-f <nom_fichier_fasta>` : Spécifie le nom du fichier Fasta pour stocker la séquence d'ADN générée. (Optionnel, valeur par défaut : "genome.fasta")
- `-q <nom_fichier_fastq>` : Spécifie le nom du fichier FastQ pour stocker les lectures générées. (Optionnel, valeur par défaut : "reads.fastq")


## Fonctionnement

Le fonctionnement de olp_asm_fun peut être décrit en plusieurs étapes :

### 1. Récupération des Séquences
Le programme commence par récupérer les séquences à partir d'un fichier FastQ fourni en entrée.

### 2. Création de *k*-mers (Optionnel)
Si spécifié, le programme crée des *k*-mers à partir des séquences.

### 3. Construction du Graphe de Chevauchement
En utilisant les séquences (et éventuellement les *k*-mers), le programme construit un graphe de chevauchement où les nœuds représentent des séquences et les arêtes représentent les chevauchements entre ces séquences.

### 4. Nettoyage du Graphe
Le graphe est nettoyé en retirant les nœuds ayant des chevauchements de score insuffisant. Cela permet de réduire le bruit et d'améliorer la précision de l'assemblage.

### 5. Assemblage des Contigs
Les contigs sont assemblés à partir des chemins dans le graphe. Le programme privilégie les chemins les plus longs dans le graphe, car ils représentent les séquences les plus probablement correctes, et termine par les chemins les plus courts.

### 6. Vérification et Validation
Une vérification est effectuée pour s'assurer que les nœuds précédemment retirés sont contenus dans les contigs. De plus, les contigs les plus longs sont vérifiés pour contenir les contigs plus courts, ce qui renforce la cohérence de l'assemblage.

### 7. Enregistrement des Contigs
Enfin, les contigs résultants sont enregistrés dans un fichier au format Fasta pour une utilisation ultérieure.

Cette approche permet de reconstruire un génome approximatif à partir des données de séquençage brutes, en utilisant des techniques de graphes pour résoudre les chevauchements entre les séquences et assembler les régions de manière cohérente.

## Licence

Ce projet est sous licence [Attribution-ShareAlike 4.0 International (CC BY-SA 4.0)](https://creativecommons.org/licenses/by-sa/4.0/).

Vous êtes libre de :

- **Partager** : copier et redistribuer le matériel sous n'importe quel format ou médium.
- **Adapter** : remixer, transformer et construire à partir du matériel.

Selon les conditions suivantes :

- **Attribution** : Vous devez donner le crédit approprié, fournir un lien vers la licence et indiquer si des modifications ont été faites. Vous devez le faire d'une manière raisonnable, mais d'une manière qui n'implique pas que l'auteur vous approuve ou approuve votre utilisation du matériel.
- **ShareAlike** : Si vous remixez, transformez ou construisez à partir du matériel, vous devez distribuer vos contributions sous la même licence que l'original.

Veuillez consulter le fichier [LICENSE](LICENSE) pour plus de détails.