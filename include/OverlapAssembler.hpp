#ifndef __OVERLAPASSEMBLER__
	#define __OVERLAPASSEMBLER__

	#include "config.hpp"
	#include <vector>
	#include <unordered_map>
	#include <string>

	class OverlapAssembler {
		private:
			int k;
			std::vector<std::string> sequences{};
			std::unordered_map<std::string, std::vector<std::string>> overlap_graph{};
			std::vector<std::string> trash{};
			std::vector<std::string> contigs{};

		public:
			/**
			 * @brief Initialise l'assembleur avec un fichier FASTQ et une longueur de k-mers optionnelle.
			 *
			 * @param filename Le nom du fichier FASTQ.
			 * @param k La longueur des k-mers à utiliser, si spécifiée.
			 * 
			 * @return Une instance BrutFoceAssembler.
			**/
			OverlapAssembler(const std::string& filename, int k = -1);

			/**
			 * @brief Crée le graphe de chevauchement à partir des séquences stockées dans l'assembleur.
			**/
			void MakeGraph();

			/**
			 * @brief Nettoie le graphe de chevauchement en supprimant les séquences ayant un score de chevauchement inférieur au seuil.
			 *
			 * @param seuil Le seuil à partir duquel les séquences sont considérées comme "inutiles" et sont supprimées.
			**/
			void CleanGraph(int seuil);

			/**
			 * @brief Recherche le meilleur nœud dans le graphe de chevauchement basé sur la longueur du chemin.
			 *
			 * @param overlap_graph Le graphe de chevauchement.
			 * 
			 * @return Une liste représentant le meilleur nœud trouvé, avec [seq1, seq2, score].
			**/
			std::vector<std::string> FindBestNode(const std::unordered_map<std::string, std::vector<std::string>>& overlap_graph);

			/**
			 * @brief Assemble les contigs à partir du graphe de chevauchement.
			**/
			void AssembleContigs();

			/**
			 * @brief les séquences contenues dans d'autres séquences de la liste donnée.
			**/
			void remove_contained_sequences();

			/**
			 * @brief Enregistre les contigs dans un fichier au format FASTA, en excluant ceux en dessous d'une longueur minimale.
			 *
			 * @param filename Le nom du fichier dans lequel enregistrer les contigs.
			 * 
			 * @param min_length La longueur minimale des contigs à conserver.
			**/
			void SaveContigs(const std::string& filename, int min_length);

	};

#endif