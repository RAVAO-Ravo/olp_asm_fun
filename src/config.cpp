#include "../include/config.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_set>

std::vector<std::string> read_fastq(const std::string& filename) {
	// Ouvrir le fichier FASTQ en mode lecture
	std::ifstream fastq(filename);

	// Vérifier si le fichier est ouvert avec succès
	if (!fastq.is_open()) {
		std::cerr << "Erreur d'ouverture du fichier FASTQ." << std::endl;
		// Retourner un vecteur vide en cas d'erreur
		return std::vector<std::string>();
	}

	// Lire toutes les lignes du fichier FASTQ
	std::vector<std::string> lines;
	std::string line;
	while (std::getline(fastq, line)) {
		lines.push_back(line);
	}

	// Initialiser l'itérateur de comptage
	int cpt = 1;

	// Extraire les séquences
	std::vector<std::string> sequences;
	for (size_t idx = 1; idx < lines.size(); idx += 4) {
		sequences.push_back(lines[idx]);

		// Affichage de la progression
		std::cout << "\rNombre de reads récupérés : [" << cpt << "]" << std::flush;
		cpt++;
	}

	// Faire le saut de ligne
	std::cout << std::endl;

	// Fermer le fichier FASTQ
	fastq.close();

	// Retourner le vecteur des séquences
	return sequences;
}

std::vector<std::string> kmerisation(const std::vector<std::string>& sequences, int k) {
	// Fonction locale pour générer les k-mers à partir d'une séquence
	auto generate_kmers = [](const std::string& sequence, int k) -> std::vector<std::string> {
		std::vector<std::string> kmers;
		for (size_t i = 0; i <= sequence.length() - k; ++i) {
			kmers.push_back(sequence.substr(i, k));
		}
		return kmers;
	};

	// Initialiser l'itérateur de comptage
	int cpt = 0;

	// Utilise une approche de compréhension de liste pour générer tous les k-mers uniques
	std::unordered_set<std::string> kmer_set;
	for (const auto& sequence : sequences) {
		auto kmers = generate_kmers(sequence, k);
		kmer_set.insert(kmers.begin(), kmers.end());

		// Affichage de la progression
		cpt += kmers.size();
		std::cout << "\rNombre de " << k << "-mers crées : [" << cpt << "]" << std::flush;
	}
	
	// Faire le saut de ligne
	std::cout << std::endl;

	// Convertit l'ensemble de k-mers en liste
	std::vector<std::string> kmer_list(kmer_set.begin(), kmer_set.end());

	// Affichage du nombre de k-mers crées
	std::cout << "\rNombre de " << k << "-mers uniques : [" << kmer_list.size() << "]" << std::endl;

	// Retourne la liste de k-mers uniques
	return kmer_list;
}

int compute_overlap(const std::string& seq1, const std::string& seq2) {
	// Fusionne les deux séquences avec un caractère spécial '$'
	std::string seq_merged = seq2 + '$' + seq1;

	// Obtient la longueur des deux séquences
	int length_seq1 = seq1.length();
	int length_seq2 = seq2.length();

	// Vérifie si l'une des séquences est vide, dans ce cas, le chevauchement est nul
	if (length_seq1 == 0 || length_seq2 == 0) {
		return 0;
	}

	// Initialise les indices pour le chevauchement
	int n = length_seq1 + length_seq2 + 1;
	int end_suffix = n - 1;
	int end_prefix = length_seq2 - 1;

	// Parcourt les séquences pour trouver le chevauchement maximal
	while (end_prefix >= 0) {
		if (seq_merged[end_prefix] == seq_merged[end_suffix]) {
			end_suffix--;
			end_prefix--;
		} else {
			// En cas de mismatch, ajuste les indices
			if (end_suffix != n - 1) {
				end_suffix = n - 1;
			} else {
				end_prefix--;
			}
		}
	}

	// Calcule et retourne la longueur du chevauchement maximal
	return n - end_suffix - 1;
}

std::vector<std::string> calculate_overlap(const std::string& seq1, const std::vector<std::string>& sequences) {
	// Initialise les variables pour le score et le meilleur score
	int score = 0;
	int best_score = 0;

	// Initialise la liste pour stocker le meilleur chevauchement et le score associé
	std::vector<std::string> result_tuple{"", "0"};
	std::string seq2 = "";

	// Parcourt toutes les séquences dans la liste
	for (const auto& seq2 : sequences) {
		// Vérifie que la séquence n'est pas la même que la séquence de référence
		if (seq1 != seq2) {
			// Calcule le score de chevauchement avec la séquence de référence
			score = compute_overlap(seq1, seq2);

			// Met à jour le meilleur score et le meilleur chevauchement si le score actuel est supérieur
			if (score > best_score) {
				best_score = score;
				result_tuple[0] = seq2;
				result_tuple[1] = std::to_string(best_score);
			}
		}
	}

	// Retourne la liste contenant le meilleur chevauchement et le score associé
	return result_tuple;
}

std::string concat_sequences(const std::string& seq1, const std::string& seq2, int score) {
	// Concaténer la séquence 1 avec la séquence 2 à partir de l'indice spécifié par le score
	std::string result_sequence = seq1 + seq2.substr(score);

	// Retourner la séquence résultante
	return result_sequence;
}