#include <iostream>
#include <ctime>
#include <cmath>
#include <cstdio>
#include <cstddef>
#include <algorithm>
#include <iomanip>
#include <unistd.h>

#include "integer_ga.h"
#include "evaluate.hpp"
#include "SSA/option.h"
#include "global.hpp"

using namespace genetic_algorithm;

extern Graph g_g;
extern HyperGraph g_hg, g_hg_2;

int main(int argc, char **argv)
{
	srand(time(NULL));

	OptionParser op(argc, argv);
	if (!op.validCheck())
	{
		printf("Parameters error, please check the readme.txt file for correct format!\n");
		return -1;
	}
	char *inFile = op.getPara("-i");
	if (inFile == NULL)
	{
		printf("Please input inFile name!\n");
		return -1;
	}

	char *outFile = op.getPara("-o");
	if (outFile == NULL)
	{
		printf("Output outFile name was not given.\n");
	}

	char *model = op.getPara("-m");
	if (model == NULL)
		model = (char *)"LT";

	if (strcmp(model, "LT") == 0)
	{
		g_g.readGraphLT(inFile);
	}
	else if (strcmp(model, "IC") == 0)
	{
		g_g.readGraphIC(inFile);
	}
	else
	{
		printf("Incorrect model option!");
		return -1;
	}

	int n = g_g.getSize();
	// std::cout << "edge: " << g_g.getEdge() << std::endl;
	g_hg.init(n);

	char *tmp;
	int t = 1;
	// todo thread
	// todo sample size maybe empirical parameter

	int k = 50;
	tmp = op.getPara("-k");
	if (tmp != NULL)
	{
		k = atoi(tmp);
	}

	int mo = 0;
	if (strcmp(model, "IC") == 0)
	{
		mo = 1;
	}

	// genetic paras
	tmp = op.getPara("-psize");
	int population_size = 800;
	if (tmp != NULL)
	{
		population_size = atoi(tmp);
	}

	tmp = op.getPara("-crate");
	double crossover_rate = 1.0;
	if (tmp != NULL)
	{
		crossover_rate = atof(tmp);
	}

	tmp = op.getPara("-rrate");
	double repair_rate = 1.0;
	if (tmp != NULL)
	{
		repair_rate = atof(tmp);
	}

	tmp = op.getPara("-maxgen");
	int max_generation = 1000;
	if (tmp != NULL)
	{
		max_generation = atoi(tmp);
	}

	tmp = op.getPara("-rpratio");
	double repair_ratio = 0.2;
	if (tmp != NULL)
	{
		repair_ratio = atof(tmp);
	}

	tmp = op.getPara("-numrf");
	if (tmp != NULL)
	{
		num_refresh = atol(tmp);
	}

	tmp = op.getPara("-amode");
	int max_age = 20;
	if (tmp != NULL)
	{
		max_age = atoi(tmp);
	}

	// Node idx starts from 1
	// One ahromosome takes n ioci
	clock_t start = clock();
	IMP imFunction = IMP(n, t, mo);

	extern vector<int> g_nodeDegree;
	extern vector<bool> g_edgeMark;
	g_nodeDegree.resize(n + 1, 0);
	g_nodeDegree_cpy.resize(n + 1, 0);
	g_edgeMark.resize(5000, false);
	g_edgeList3.reserve(num_refresh);

	size_t chrom_len = k;
	IntegerGA GA(chrom_len, imFunction, n);

	GA.population_size(population_size);
	GA.crossover_rate(crossover_rate);
	GA.repair_ratio(repair_ratio);
	GA.rp_rate(repair_rate);
	GA.max_age(max_age);
	GA.max_gen(max_generation);
	GA.crossover_method(IntegerGA::CrossoverMethod::set_based);
	GA.selection_method(IntegerGA::SogaSelection::tournament);

	GA.run();
	
	float run_time = (float)(clock() - start) / CLOCKS_PER_SEC;
	float run_memory = getCurrentMemoryUsage();
	cout << "Time: " << run_time << "s" << endl;
	cout << "Memory: " << run_memory << " MB" << endl;

	auto sols = GA.population();
	size_t num_evals = GA.num_fitness_evals();
	std::cout << "EvalTimes: " << num_evals << std::endl;

	std::map<string, double> ht;
	for (const auto &sol : sols)
	{
		string s;
		for (const auto &gene : sol.chromosome)
		{
			// cout << gene + 1 << ", ";
			s.append(to_string(gene + 1) + " ");
		}
		// cout << ")" << endl;

		if (ht.contains(s))
		{
			++ht[s];
		}
		else
		{
			ht[s] = 0.0;
		}
	}

	cout << "Best: ";
	std::vector<pair<string, double>> tmp_v(ht.begin(), ht.end());
	std::sort(tmp_v.begin(), tmp_v.end(), [](std::pair<string, double> &a, std::pair<string, double> &b)
			  { return a.second > b.second; });

	string p_best = tmp_v[0].first;
	double best_pct = tmp_v[0].second / static_cast<double>(sols.size());
	cout << p_best << endl
		 << best_pct << "%" << endl;
	
	/* file stuff */
	std::ofstream file;
	file.open("test.seeds");
	file << p_best;
	file.close();
	if (outFile == NULL)
		return 0;

	// open a new file
	size_t new_order = 1;
	char cbuffer[100];
	char append[10];
	strcpy(cbuffer, outFile);
	while(access(cbuffer, F_OK) != -1)
	{
		strcpy(cbuffer, outFile);
		sprintf(append, "(%ld)", new_order);
		strcat(cbuffer, append);
		++new_order;
	}

	file.open(cbuffer);

	file << "[INPUT]: " << inFile << '\n';
	file << "[MODEL]: " << model << '\n';
	file << "[k]: " << k << '\n';
	file << "[TIME]: " << run_time << '\n';
	file << "[MEMORY]: " << run_memory << '\n';
	file << "[GENERATION]: " << GA.generation_cntr() << '\n';
	file << "[POP_SIZE]: " << GA.population_size() << '\n';
	file << "[EVAL_TIME]: " << GA.num_fitness_evals() << '\n';
	file << "[BEST]:" << '\n';
	file << p_best << "\n" << best_pct << "%" << '\n';
	file << "[POPULATION]:" << '\n';
	for (const auto &sol : sols)
	{
		file << sol.fitness[0] << ", " << sol.fitness[1] << " : ";
		for (const auto &gene : sol.chromosome)
		{
			file << gene + 1 << " ";
		}
		file << '\n';
	}
	file.close();

	return 0;
}