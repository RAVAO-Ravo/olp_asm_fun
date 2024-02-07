#include "../include/config.hpp"
#include "../include/OverlapAssembler.hpp"
#include <iostream>
#include <vector>
#include <unordered_map>
#include <string>
#include <fstream>

OverlapAssembler::OverlapAssembler(const std::string& filename, int k) {
	// Vérifie si la longueur des k-mers est spécifiée
	if (k != -1) {
		// k-mérise les séquences à partir du fichier FASTQ
		sequences = kmerisation(read_fastq(filename), k);
	} else {
		// Utilise les séquences brutes du fichier FASTQ
		sequences = read_fastq(filename);
	}

	// Initialise les autres attributs de l'assembleur
	this->k = k;
	this->overlap_graph = {};
	this->contigs = {};
	this->trash = {};
}

void OverlapAssembler::MakeGraph() {
	// Récupère les séquences à partir de l'assembleur
	std::vector<std::string> sequences = this->sequences;

	// Initialiser l'itérateur et le total de la barre de progression
	int cpt = 1;
	const int total = sequences.size();

	// Créer le graphe de chevauchement
	for (const std::string& seq1 : sequences) {
		this->overlap_graph[seq1] = calculate_overlap(seq1, sequences);

		// Affichage de la progression
		std::cout << "\rNombre de nœuds crées : [" << cpt << "/" << total << "]" << std::flush;
		cpt++;
	}

	// Faire le saut de ligne
	std::cout << std::endl;
}

void OverlapAssembler::CleanGraph(int seuil) {
	// Récupère le graphe de chevauchement à partir de l'assembleur
	std::unordered_map<std::string, std::vector<std::string>> overlap_graph = this->overlap_graph;

	// Vérifie si la longueur des k-mers est spécifiée et ajuste le seuil si nécessaire
	if (this->k != -1 && this->k < seuil) {
		seuil = this->k - static_cast<int>(0.2 * this->k);
		std::cout << "\nLe seuil a été changé en " << seuil << ", car la valeur était supérieure à la taille des k-mers.\n" << std::endl;
	}

	// Initialiser l'itérateur et le total de la barre de progression
	int cpt = 1;
	const int total = overlap_graph.size();

	// Identifie les séquences à supprimer (celles avec un score inférieur au seuil)
	std::vector<std::string> trash;
	for (const auto& item : overlap_graph) {
		const std::string& seq1 = item.first;
		int score = std::stoi(item.second[1]);

		if (score < seuil) {
			trash.push_back(seq1);
		}

		// Affichage de la progression
		std::cout << "\rNombre de nœuds traités : [" << cpt << "/" << total << "]" << std::flush;
		cpt++;
	}

	// Faire le saut de ligne
	std::cout << std::endl;

	// Supprime les séquences "inutiles" du graphe de chevauchement
	for (const auto& sequence : trash) {
		overlap_graph.erase(sequence);
	}

	// Affichage du nombre de nœuds restants et des nœuds rejetés
	std::cout << "Nombre de nœuds restants : [" << overlap_graph.size() << "]" << std::endl;
	std::cout << "Nombre de nœuds rejetés : [" << trash.size() << "]" << std::endl;

	// Met à jour l'attribut overlap_graph de l'assembleur
	this->overlap_graph = overlap_graph;

	// Stocke les séquences "inutiles" dans l'attribut trash de l'assembleur
	this->trash = trash;
}

std::vector<std::string> OverlapAssembler::FindBestNode(const std::unordered_map<std::string, std::vector<std::string>>& overlap_graph) {
	// Fonction locale pour calculer la longueur du chemin à partir d'un nœud
	auto GetLengthPath = [&](const std::string& node) -> int {
		std::vector<std::string> already_done;

		// Vérifie si le nœud est présent dans le graphe et n'a pas de successeur
		if (overlap_graph.find(node) == overlap_graph.end() || overlap_graph.at(node)[0] == "") {
			return 0;
		}

		// Initialise la liste des nœuds déjà visités avec le nœud actuel
		already_done.push_back(node);

		// Initialise la longueur du chemin
		int length = 1;
		std::string next_node = overlap_graph.at(node)[0];

		// Parcourt les nœuds suivants dans le chemin jusqu'à trouver un nœud sans successeur ou un nœud déjà visité
		while (overlap_graph.find(next_node) != overlap_graph.end() && next_node != "" && std::find(already_done.begin(), already_done.end(), next_node) == already_done.end()) {
			already_done.push_back(next_node);
			length++;
			next_node = overlap_graph.at(next_node)[0];
		}

		// Retourner la longueur du chemin
		return length;
	};

	// Initialise les variables pour stocker le meilleur nœud
	int best_length = 0;
	std::vector<std::string> best_node{"", "", "0"};

	// Parcourt tous les nœuds dans le graphe
	for (const auto& item : overlap_graph) {
		const std::string& seq1 = item.first;
		const std::string& seq2 = item.second[0];

		// Calcule la longueur du chemin pour le nœud actuel
		int current_length = GetLengthPath(seq1);

		// Met à jour le meilleur nœud si la longueur actuelle est supérieure
		if (best_length < current_length) {
			best_length = current_length;
			best_node[0] = seq1;
			best_node[1] = seq2;
			best_node[2] = item.second[1];
		}
	}

	// Retourne le meilleur nœud trouvé
	return best_node;
};

void OverlapAssembler::AssembleContigs() {
	// Récupère le graphe de chevauchement à partir de l'assembleur
	std::unordered_map<std::string, std::vector<std::string>> overlap_graph = this->overlap_graph;

	// Initialise la liste des contigs
	std::vector<std::string> contigs;

	// Trouve le meilleur nœud de départ
	std::vector<std::string> best_node = FindBestNode(overlap_graph);

	// Initialise le premier contig
	std::string contig = concat_sequences(best_node[0], best_node[1], std::stoi(best_node[2]));

	// Initialise le total de la barre de progression
	const int total = overlap_graph.size();

	// Boucle jusqu'à ce que tous les nœuds soient utilisés
	while (!overlap_graph.empty()) {
		// Vérifie si le nœud suivant est présent dans le graphe
		if (overlap_graph.find(best_node[1]) != overlap_graph.end()) {
			// Supprime le nœud courant du graphe
			overlap_graph.erase(best_node[0]);

			// Met à jour les nœuds courant et suivant, ainsi que le score
			best_node = {best_node[1], overlap_graph.at(best_node[1])[0], overlap_graph.at(best_node[1])[1]};

			// Concatène la séquence du nœud suivant au contig actuel
			contig = concat_sequences(contig, best_node[1], std::stoi(best_node[2]));
		} 
		else {
			// Ajoute le contig actuel à la liste des contigs
			contigs.push_back(contig);

			// Supprime le nœud courant du graphe
			overlap_graph.erase(best_node[0]);

			// Trouve le meilleur nœud de départ pour le prochain contig
			best_node = FindBestNode(overlap_graph);

			// Initialise un nouveau contig avec le meilleur nœud trouvé
			contig = concat_sequences(best_node[0], best_node[1], std::stoi(best_node[2]));
		}

		// Affichage de la progression
		std::cout << "\rNombre de nœuds traités : [" << total - overlap_graph.size() << "/" << total << "]" << std::flush;
	}

	// Faire le saut de ligne
	std::cout << std::endl;

	// Affichage du nombre de reads obtenus
	std::cout << "Nombre de contigs générés: [" << contigs.size() << "]" << std::endl;

	// Stocke les contigs résultants dans l'attribut contigs de l'assembleur
	this->contigs = contigs;
}

void OverlapAssembler::remove_contained_sequences() {
	// Récupérer les contigs
	std::vector<std::string> contigs = this->contigs;
	contigs.insert(contigs.end(), this->trash.begin(), this->trash.end());

	// Initialiser la liste des index à retirer
	std::vector<int> contained_sequences_index;

	// Initialiser l'itérateur et le total de la barre de progression
	int cpt = 1;
	const int total = contigs.size();

	// Faire les comparaisons
	for (size_t i = 0; i < contigs.size(); ++i) {
		for (size_t j = 0; j < contigs.size(); ++j) {
			// Si le contig1 est plus grand que le contig2, passer à l'itération suivante
			if (contigs[i].length() > contigs[j].length()) {
				continue;
			}
			// Si le contig1 est plus petit ou égal au contig2
			else if (contigs[i].length() <= contigs[j].length()) {
				// Vérifier si le contig1 est contenu dans le contig2
				if (i != j && contigs[j].find(contigs[i]) != std::string::npos) {
					contained_sequences_index.push_back(i);
					break;
				}
			}
		}

		// Affichage de la progression
		std::cout << "\rNombre de contigs vérifiés : [" << cpt << "/" << total << "]" << std::flush;
		cpt++;
	}

	// Faire le saut de ligne
	std::cout << std::endl;

	// Retirer les contigs contenus dans d'autres
	for (auto it = contained_sequences_index.rbegin(); it != contained_sequences_index.rend(); ++it) {
		contigs.erase(contigs.begin() + *it);
	}

	// Affichage du nombre de contigs restants
	std::cout << "Nombre de contigs restant : [" << contigs.size() << "]" << std::endl;

	// Mettre à jour la liste des contigs
	this->contigs = contigs;
}

void OverlapAssembler::SaveContigs(const std::string& filename, int min_length) {
	// Combine les contigs et la corbeille
	std::vector<std::string> contigs = this->contigs;

	// Ouvre le fichier FASTA en mode écriture
	std::ofstream fasta(filename);

	if (!fasta.is_open()) {
		std::cerr << "Erreur lors de l'ouverture du fichier : " << filename << std::endl;
		return;
	}

	// Initialiser l'itérateur et le total de la barre de progression
	int cpt = 1;
	const int total = contigs.size();

	// Parcourt tous les contigs
	for (const std::string& contig : contigs) {
		// Vérifie si la longueur du contig est supérieure à la longueur minimale
		int contig_length = contig.length();
		if (min_length <= contig_length) {
			// Génère une ligne au format FASTA pour le contig et l'écrit dans le fichier
			fasta << ">contig" << cpt << '\n' << contig << '\n' << '\n';

			// Affichage de la progression
			std::cout << "\rNombre de contigs sauvegardés : [" << cpt << "/" << total << "]" << std::flush;
			cpt++;
		}
	}

	// Faire le saut de ligne
	std::cout << std::endl;

	// Ferme le fichier
	fasta.close();
}