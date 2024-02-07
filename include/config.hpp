#ifndef __CONFIG__
	#define __CONFIG__

	#include <vector>
	#include <string>
	#include <unordered_set>
	
	/**
	 * @brief Lit un fichier FASTQ et extrait les séquences.
	 *
	 * @param filename Le nom du fichier FASTQ.
	 * 
	 * @return Un vecteur des séquences extraites.
	**/
	std::vector<std::string> read_fastq(const std::string& filename);

	/**
	 * @brief Effectue la k-mérisation des séquences en générant tous les k-mers uniques.
	 *
	 * @param sequences Une liste de séquences à k-mériser.
	 * @param k La longueur des k-mers.
	 * 
	 * @return Une liste de tous les k-mers uniques générés à partir des séquences.
	**/
	std::vector<std::string> kmerisation(const std::vector<std::string>& sequences, int k);

	/**
	 * @brief Concatène deux séquences en utilisant un score pour déterminer le point de départ de la deuxième séquence.
	 *
	 * @param seq1 La première séquence.
	 * @param seq2 La deuxième séquence.
	 * @param score Le score pour déterminer le point de départ de la deuxième séquence.
	 *
	 * @return La séquence résultante après la concaténation.
	**/
	std::string concat_sequences(const std::string& seq1, const std::string& seq2, int score);

	/**
	 * @brief Calcule la longueur du chevauchement maximal entre deux séquences.
	 *
	 * @param seq1 La première séquence.
	 * @param seq2 La deuxième séquence.
	 * 
	 * @return La longueur du chevauchement maximal entre les deux séquences.
	**/
	int compute_overlap(const std::string& seq1, const std::string& seq2);

	/**
	 * @brief Calcule le meilleur chevauchement et le score associé avec une séquence parmi une liste de séquences.
	 *
	 * @param seq1 La séquence de référence.
	 * @param sequences Une liste de séquences avec lesquelles comparer la séquence de référence.
	 * 
	 * @return Une liste contenant le meilleur chevauchement et le score associé.
	**/
	std::vector<std::string> calculate_overlap(const std::string& seq1, const std::vector<std::string>& sequences);

#endif