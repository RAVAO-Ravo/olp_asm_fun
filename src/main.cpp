#include "../include/config.hpp"
#include "../include/OverlapAssembler.hpp"
#include "../include/cxxopts.hpp"
#include <iostream>
#include <string>
#include <chrono>

int main(int argc, char* argv[]) {
	// Parseur d'arguments
	cxxopts::Options options("olp_asm", "\nRéalise un assemblage (par un graphe de chevauchement) à partir d'un fichier FastQ.\n");
	options.add_options()
		("q,fastq", "Nom du fichier fastq à utiliser.", cxxopts::value<std::string>())
		("k,kmers_length", "Longueur des k-mers à utiliser.", cxxopts::value<int>()->default_value("-1"))
		("s,seuil", "Le score de chevauchement minimum pour garder un nœud dans le graphe.", cxxopts::value<int>()->default_value("10"))
		("f,fasta", "Nom du fichier fasta qui contiendra les contigs.", cxxopts::value<std::string>())
		("m,min_length", "Longueur minimum d'un contig pour être garder.", cxxopts::value<int>()->default_value("0"))
		("h,help", "Affiche l'aide.");
	auto result = options.parse(argc, argv);

	// Affiche l'aide en cas de demande
	if (result.count("help")) {
		std::cout << options.help() << std::endl;
		return 0;
	}
	
	// Récupère les paramètres
	std::string fastq = result["fastq"].as<std::string>();
	int kmers_length = result["kmers_length"].as<int>();
	int seuil = result["seuil"].as<int>();
	std::string fasta = result["fasta"].as<std::string>();
	int min_length = result["min_length"].as<int>();
	
	// Démmarrage des traitements
	std::cout << "\n--- DÉBUT ---" << std::endl;

	// Récupération des séquences à utiliser
	std::cout << "\n- Récupération des séquences -" << std::endl;
	auto start_time = std::chrono::high_resolution_clock::now();
	OverlapAssembler assembler = OverlapAssembler(fastq, kmers_length);
	auto end_time = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
	std::cout << "Temps d'exécution : " << duration.count() << " seconds" << std::endl;

	// Création du graphe de chevauchement
	start_time = std::chrono::high_resolution_clock::now();
	std::cout << "\n- Création du graphe de chevauchement -" << std::endl;
	assembler.MakeGraph();
	end_time = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
	std::cout << "Temps d'exécution : " << duration.count() << " seconds" << std::endl;

	// Nettoyage du graphe
	start_time = std::chrono::high_resolution_clock::now();
	std::cout << "\n- Nettoyage du graphe (seuil = " << seuil << ") -" << std::endl;
	assembler.CleanGraph(seuil);
	end_time = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
	std::cout << "Temps d'exécution : " << duration.count() << " seconds" << std::endl;

	// Assemblage des contigs
	start_time = std::chrono::high_resolution_clock::now();
	std::cout << "\n- Assemblage des contigs -" << std::endl;
	assembler.AssembleContigs();
	end_time = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
	std::cout << "Temps d'exécution : " << duration.count() << " seconds" << std::endl;

	// Vérifier si des séquences sont contenues dans d'autres
	start_time = std::chrono::high_resolution_clock::now();
	std::cout << "\n- Retrait des contigs contenus -" << std::endl;
	assembler.remove_contained_sequences();
	end_time = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
	std::cout << "Temps d'exécution : " << duration.count() << " seconds" << std::endl;

	// Sauvegarder les contigs obtenus
	start_time = std::chrono::high_resolution_clock::now();
	std::cout << "\n- Sauvegarde des contigs ⩾ " << min_length << " -" << std::endl;
	assembler.SaveContigs(fasta, min_length);
	end_time = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
	std::cout << "Temps d'exécution : " << duration.count() << " seconds" << std::endl;

	// Fin des traitements
	std::cout << "\n--- FIN ---\n" << std::endl;

	return 0;
}