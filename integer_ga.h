/*
 *  MIT License
 *
 *  Copyright (c) 2021 Kriszti�n Rug�si
 *
 *  Permission is hereby granted, free of charge, to any person obtaining a copy
 *  of this softwareand associated documentation files (the "Software"), to deal
 *  in the Software without restriction, including without limitation the rights
 *  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 *  copies of the Software, and to permit persons to whom the Software is
 *  furnished to do so, subject to the following conditions:
 *
 *  The above copyright noticeand this permission notice shall be included in all
 *  copies or substantial portions of the Software.
 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 *  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 *  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 *  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 *  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 *  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 *  SOFTWARE.
 */

/**
 * This file contains the integer coded genetic algorithm class.
 *
 * @file integer_ga.h
 */

#ifndef GA_INTEGER_GA_H
#define GA_INTEGER_GA_H

#include <cstddef>
#include <algorithm>
#include <malloc.h>

#include "GA_src/base_ga.h"
#include "evaluate.hpp"
#include "global.hpp"

namespace genetic_algorithm
{
    /**
     * Integer coded GA. \n
     * Same as @ref BinaryGA, but the genes of the chromosomes can be any integer on [0, base], not just 0 or 1. \n
     * It also uses a slightly different mutation function with swaps and inversions.
     */
    class IntegerGA : public GA<int>
    {
    public:
        /**
         * Possible crossover operators that can be used in the IntegerGA. \n
         * These crossover method are the same as the ones used in the binary coded algorithm. \n
         * Set the crossover method used in the algorithm with @ref crossover_method. \n
         * The function used for the crossovers with the custom method can be set with @ref setCrossoverFunction.
         */
        enum class CrossoverMethod
        {
            set_based,
            custom /**< Custom crossover operator defined by the user. @see setCrossoverMethod */
        };

        /**
         * Basic contructor for the IntegerGA.
         *
         * @param chrom_len The number of genes in each chromosome.
         * @param fitness_function The fitness function used in the algorithm.
         * @param base The number of values a gene can take. Must be > 1. If 2, same as the @ref BinaryGA.
         */
        IntegerGA(size_t chrom_len, IMP im_function, size_t base);
        // ~IntegerGA();

        /*
         * self-defined run function
         */
        CandidateVec run();

        /**
         * Sets the crossover function used in the algorithm to @f.
         * @see CrossoverMethod
         *
         * @param method The crossover function to use.
         */
        void crossover_method(crossoverFunction_t f);

        /**
         * Sets the crossover method used in the algorithm to @p method.
         * @see CrossoverMethod
         *
         * @param method The crossover method to use.
         */
        void crossover_method(CrossoverMethod method);
        [[nodiscard]] CrossoverMethod crossover_method() const;

        /**
         * Sets the number of crossover points used in the crossovers to @p n if the n_point crossover method selected. \n
         * The number of crossover points must be at least 1.
         * @see crossover_method @see CrossoverMethod
         *
         * @param n The number of crossover points.
         */
        void num_crossover_points(size_t n);
        [[nodiscard]] size_t num_crossover_points() const;

        /**
         * Sets the number of values a gene can take to @p base. \n
         * The value of the base must be at least 2, and the GA is essentially the same as the
         * BinaryGA if the base is set to 2.
         *
         * @param base The number of values a gene can be.
         */
        void base(size_t base);
        [[nodiscard]] size_t base() const;

        void rp_rate(double rp_rate);
        [[nodiscard]] double rp_rate() const;

        void repair_ratio(double repair_ratio);
        [[nodiscard]] double repair_ratio() const;

        void max_age(int max_age);
        [[nodiscard]] int max_age() const;

    private:
        /* Find all Pareto fronts in the population and also assign the nondomination ranks of the candidates (assuming fitness maximization). */
        static std::vector<std::vector<size_t>> nonDominatedSort(Population &pop);

        /* Calculate the crowding distances of the candidates in each pareto front in pfronts of the population. */
        static void calcCrowdingDistances(Population &pop, std::vector<std::vector<size_t>> &pfronts);

        /* Returns true if lhs is better than rhs. */
        static bool crowdedCompare(const Candidate &lhs, const Candidate &rhs);

        static Candidate nsga2Select(const Population &pop);

        /* Create the population of the next generation from the old population and the children. */
        Population updateNsga2Population(Population &old_pop, CandidateVec &children) const;

        CrossoverMethod crossover_method_ = CrossoverMethod::set_based;
        size_t num_crossover_points_ = 3;
        size_t base_ = 4;

        double rp_rate_ = 0.75;
        double repair_ratio_ = 0.2;
        int max_age_ = 20;

        Candidate generateCandidate() const override;
        CandidatePair crossover(const Candidate &parent1, const Candidate &parent2) const override;
        inline void repair(Population &pop, double reservation);
        inline void mutate(Candidate &child) const;

        static CandidatePair setCrossover(const Candidate &parent1, const Candidate &parent2, double pc);
        static Candidate setRepair(const Candidate &parent1, const Candidate &parent2);

        inline void evaluate_older(Population &pop);
        inline void evaluate_new(Population &pop);

        IMP im_function_;

        vector<double> degree;
        vector<int> seeds;
        /* Candiate pool stuff */
    };

} // namespace genetic_algorithm

/* IMPLEMENTATION */

#include <algorithm>
#include <vector>
#include <unordered_set>
#include <utility>
#include <stdexcept>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <map>
#include <unistd.h>
#include "GA_src/rng.h"

namespace genetic_algorithm
{
    inline IntegerGA::IntegerGA(size_t chrom_len, IMP im_function, size_t base)
        : GA(chrom_len, im_function), base_(base), im_function_(im_function)
    {
        if (base < 2)
            throw std::invalid_argument("The base must be at least 2.");

        if (chrom_len == 0)
        {
            throw std::invalid_argument("The chromosome length must be at least 1.");
        }
        if (fitnessFunction == nullptr)
        {
            throw std::invalid_argument("The fitness function is a nullptr.");
        }
    }

    inline IntegerGA::CandidateVec IntegerGA::run()
    {
        using namespace std;

        // std::ofstream file;
        // stringstream ss;
        // ss << unsigned(getpid());
        // string filename = "his/";
        // filename.append(ss.str());
        // filename.append("his.txt");
        // file.open(filename.c_str());

        // cout << filename << endl;

        /* Init. */
        num_fitness_evals_ = 0;
        population_.clear();

        degree.resize(chrom_len_ + 1, 0);
        seeds.resize(chrom_len_ + 1, 0);

        /* Create and evaluate the initial population. */
        population_.reserve(population_size_ + 1);
        refresh(im_function_.getMo());
        count_edge_3();
        while (population_.size() < population_size_)
        {
            vector<int> chrom;
            chrom.reserve(chrom_len_);
            while (chrom.size() < chrom_len_)
            {
                size_t e_idx = rng::randomIdx(g_edgeList3.size());
                int nodeIdx = g_edgeList3[e_idx][rng::randomIdx(g_edgeList3[e_idx].size() - 1) + 1] - 1;
                if (!binary_search(chrom.begin(), chrom.end(), nodeIdx))
                {
                    chrom.push_back(nodeIdx);
                    sort(chrom.begin(), chrom.end());
                }
            }
            Candidate c(chrom);
            c.fitness = {0, 0};
            c.is_evaluated = false;
            population_.push_back(c);
        }
        seeds.clear();
        buildSeedSet(g_hg, seeds, g_g.getSize(), chrom_len_, degree);
        repair(population_, 0.0);

        evaluate_new(population_);

        refresh(im_function_.getMo());

        // cout << "rpatio: " << repair_ratio_ << endl;
        // cout << "maxage: " << max_age_ << endl;
        // cout << "rfnum: " << num_refresh << endl;
        // cout << "psize: " << population_size_ << endl;
        // cout << "maxgen: " << max_gen_ << endl;
        // cout << "rprate: " << rp_rate_ << endl;
        cout << repair_ratio_ << "\t" << max_age_ << "\t" << num_refresh << "\t" << population_size_ << "\t" << max_gen_ << "\t" << rp_rate_ << endl;

        // printf("gen: %2.d\r", 0);
        // fflush(stdout);

        /* Other generations. */
        generation_cntr_ = 1;
        size_t num_children = population_size_ + population_size_ % 2;
        vector<CandidatePair> parent_pairs(num_children / 2);
        vector<Candidate> children;
        children.reserve(num_children * 2 + 1);

        // appearence map for early stop
        map<string, double> ht;
        map<string, vector<size_t>> ht_idx;

        while (generation_cntr_ < max_gen_)
        {
            /* Selections. */
            generate(execution::par_unseq, parent_pairs.begin(), parent_pairs.end(),
                     [this]() -> CandidatePair
                     {
                         return move(make_pair(nsga2Select(population_), nsga2Select(population_)));
                     });

            /* Crossovers. */
            for_each(execution::par_unseq, parent_pairs.begin(), parent_pairs.end(),
                     [this](CandidatePair &p) -> void
                     {
                         p = setCrossover(p.first, p.second, crossover_rate_);
                     });

            children.clear();
            for (size_t i = 0; i < parent_pairs.size(); i++)
            {
                children.push_back(move(parent_pairs[i].first));
                children.push_back(move(parent_pairs[i].second));
            }

            seeds.clear();
            buildSeedSet(g_hg, seeds, g_g.getSize(), chrom_len_, degree);

            // Candidate best_c(seeds);
            // best_c.is_evaluated = false;
            // size_t init_size = children.size();
            // children.resize(init_size * 2);
            // for (size_t i = 0; i < init_size; i++)
            // {
            //     // printf("%d %o\t",i , &(pop[i + init_size].chromosome));
            //     children[i + init_size] = setRepair(children[i], best_c);
            // }

            /* Apply repair function to the children if set. */
            repair(children, 1.0);

            /* Overwrite the current population with the children. */
            evaluate_new(children);
            population_ = updateNsga2Population(population_, children);

            // early stop
            // file << "g:" << generation_cntr_ << endl;
            if (generation_cntr_ % 1 == 0)
            {
                ht.clear();
                ht_idx.clear();
                for (size_t i = 0; i < population_.size(); i++)
                {
                    string s;
                    for (const auto &gene : population_[i].chromosome)
                    {
                        s.append(to_string(gene + 1) + " ");
                    }
                    if (ht.contains(s))
                    {
                        ++ht[s];
                        ht_idx[s].push_back(i);
                    }
                    else
                    {
                        ht[s] = 0.0;
                        ht_idx[s].reserve(480);
                    }
                }

                vector<pair<string, double>> ht_v(ht.begin(), ht.end());
                partial_sort(ht_v.begin(), ht_v.begin() + min((size_t)ht_v.size(), (size_t)5), ht_v.end(), [](pair<string, double> &a, pair<string, double> &b)
                             { return a.second > b.second; });
                pair<string, double> max_set = ht_v[0];

                // cout << generation_cntr_ << " =========== " << endl;
                // for (size_t i = 0; i < 5; i++)
                // {
                //     cout << ht_v[i].first << " : " << ht_v[i].second / static_cast<double>(population_size_) << endl;
                // }
                // file << max_set.first << " : " << max_set.second << endl;
                double convergence_ratio = max_set.second / static_cast<double>(population_size_);
                crossover_rate_ = 1 - convergence_ratio;
                // rp_rate_ = 1 - sqrt(max_set.second / static_cast<double>(population_size_));
                // repair_ratio_ =  max(min(sqrt(max_set.second / static_cast<double>(population_size_)), 0.1), 0.3);

                vector<pair<string, double>>().swap(ht_v); // memory leaking?
                if (convergence_ratio > 0.6)
                {
                    cout << "Convergence: " << convergence_ratio << "%" << endl;
                    break;
                }

                // printf("gen: %2.ld   %.4lf\r", generation_cntr_, max_set.second / static_cast<double>(population_size_));
                // fflush(stdout);
            }

            refresh(im_function_.getMo());
            evaluate_older(population_);

            generation_cntr_++;
        }
        cout << "Generation: " << generation_cntr_ << endl;

        map<string, double>().swap(ht);
        vector<CandidatePair>().swap(parent_pairs);
        vector<Candidate>().swap(children);
        vector<double>().swap(degree);
        vector<int>().swap(seeds);
        malloc_trim(0);

        // file.close();

        return solutions_;
    }

    // ~IntegerGA()
    // {
    //     /* release memory */
    // }

    inline void IntegerGA::repair(Population &pop, double reservation)
    {
        if (rng::randomReal(0.0, 1.0) < rp_rate_ && repair_ratio_ > 0.0)
        {
            size_t tmp = chrom_len_ * repair_ratio_;
            if (tmp == 0)
                tmp = 1;
            vector<size_t> chosen(tmp);
            for (size_t i = 0; i < tmp; i++)
            {
                chosen[i] = rng::randomIdx(chrom_len_);
            }

            size_t init_size = pop.size();
            if (rng::randomReal(0.0, 1.0) < reservation)
            {
                pop.resize(init_size * 2);
                for (size_t i = 0; i < init_size; i++)
                {
                    // printf("%d %o\t",i , &(pop[i + init_size].chromosome));
                    pop[i + init_size] = Candidate(pop[i].chromosome);
                }
            }

            for (size_t i = 0; i < init_size; ++i)
            {
                for (size_t j = 0; j < tmp; j++)
                {
                    int nodeIdx = seeds[rng::randomIdx(seeds.size())] - 1;
                    if (!binary_search(pop[i].chromosome.begin(), pop[i].chromosome.end(), nodeIdx))
                    {
                        pop[i].chromosome[chosen[j]] = nodeIdx;
                        sort(pop[i].chromosome.begin(), pop[i].chromosome.end());
                    }
                }
                pop[i].is_evaluated = false;
            }
        }
    }

    inline void IntegerGA::crossover_method(crossoverFunction_t f)
    {
        if (f == nullptr)
            throw std::invalid_argument("The function used for the crossovers can't be a nullptr.");

        crossover_method_ = CrossoverMethod::custom;
        customCrossover = f;
    }

    inline void IntegerGA::crossover_method(CrossoverMethod method)
    {
        if (static_cast<size_t>(method) > 4)
            throw std::invalid_argument("Invalid crossover method selected.");

        crossover_method_ = method;
    }

    inline IntegerGA::CrossoverMethod IntegerGA::crossover_method() const
    {
        return crossover_method_;
    }

    inline void IntegerGA::base(size_t base)
    {
        if (base < 2)
            throw std::invalid_argument("The base must be at least 2.");

        base_ = base;
    }

    inline size_t IntegerGA::base() const
    {
        return base_;
    }

    inline void IntegerGA::rp_rate(double rp_rate)
    {
        rp_rate_ = rp_rate;
    }
    inline double IntegerGA::rp_rate() const
    {
        return rp_rate_;
    }

    inline void IntegerGA::repair_ratio(double repair_ratio)
    {
        repair_ratio_ = repair_ratio;
    }
    inline double IntegerGA::repair_ratio() const
    {
        return repair_ratio_;
    }

    inline void IntegerGA::max_age(int max_age)
    {
        max_age_ = max_age;
    }
    inline int IntegerGA::max_age() const
    {
        return max_age_;
    }

    inline IntegerGA::Candidate IntegerGA::generateCandidate() const
    {
        assert(chrom_len_ > 0);
        assert(base_ > 1);

        Candidate sol;
        sol.chromosome.reserve(chrom_len_);
        for (size_t i = 0; i < chrom_len_; i++)
        {
            sol.chromosome.push_back(rng::randomInt(size_t{0}, base_ - 1));
        }

        std::sort(sol.chromosome.begin(), sol.chromosome.end());
        return sol;
    }

    inline IntegerGA::CandidatePair IntegerGA::crossover(const Candidate &parent1, const Candidate &parent2) const
    {
        /* Edge case. No point in performing the mutations if the parents are the same. */
        if (parent1 == parent2)
            return std::make_pair(parent1, parent2);

        return setCrossover(parent1, parent2, crossover_rate_);
    }

    inline IntegerGA::CandidatePair IntegerGA::setCrossover(const Candidate &parent1, const Candidate &parent2, double pc)
    {
        using namespace std;
        assert(parent1.chromosome.size() == parent2.chromosome.size());
        assert(0.0 <= pc && pc <= 1.0);
        Candidate child1(parent1), child2(parent2);

        /* Perform crossover with pc probability. */
        if (rng::randomReal() <= pc)
        {
            size_t len = child1.chromosome.size();
            // find unique loci and shared loci
            // random assign unique loci to vectors
            // merge assigned loci and shared loci together

            vector<int> uniqueList, overlappingList;
            uniqueList.reserve(len);
            overlappingList.reserve(len);

            size_t i = 0, j = 0, k;

            while (i < len && j < len)
            {
                while (i < len && parent1.chromosome[i] < parent2.chromosome[j])
                {
                    uniqueList.push_back(parent1.chromosome[i]);
                    ++i;
                }
                while (j < len && (i == len || parent1.chromosome[i] > parent2.chromosome[j]))
                {
                    uniqueList.push_back(parent2.chromosome[j]);
                    ++j;
                }
                if (j < len && i < len && parent1.chromosome[i] == parent2.chromosome[j])
                {
                    overlappingList.push_back(parent2.chromosome[j]);
                    ++i, ++j;
                }
            }

            while (i < len)
            {
                uniqueList.push_back(parent1.chromosome[i]);
                ++i;
            }

            while (j < len)
            {
                uniqueList.push_back(parent2.chromosome[j]);
                ++j;
            }

            vector<int> v1, v2;
            v1.reserve(len);
            v2.reserve(len);
            i = len - overlappingList.size(), j = i, k = 0;
            while (i != 0 && j != 0)
            {
                if (rng::randomBool())
                {
                    v1.push_back(uniqueList[k]);
                    --i;
                }
                else
                {
                    v2.push_back(uniqueList[k]);
                    --j;
                }
                k++;
            }

            while (i != 0)
            {
                v1.push_back(uniqueList[k]);
                --i;
                k++;
            }

            while (j != 0)
            {
                v2.push_back(uniqueList[k]);
                --j;
                k++;
            }

            merge(v1.begin(), v1.end(), overlappingList.begin(), overlappingList.end(), child1.chromosome.begin());
            merge(v2.begin(), v2.end(), overlappingList.begin(), overlappingList.end(), child2.chromosome.begin());

            if (child1 != parent1)
            {
                child1.is_evaluated = false;
                child2.is_evaluated = false;
            }

            vector<int>().swap(v1);
            vector<int>().swap(v2);
            vector<int>().swap(uniqueList);
            vector<int>().swap(overlappingList);
        }

        return make_pair(child1, child2);
    }

    inline IntegerGA::Candidate IntegerGA::setRepair(const Candidate &parent1, const Candidate &parent2)
    {
        using namespace std;
        assert(parent1.chromosome.size() == parent2.chromosome.size());
        Candidate child1(parent1);

        /* Perform crossover with pc probability. */
        size_t len = child1.chromosome.size();
        // find unique loci and shared loci
        // random assign unique loci to vectors
        // merge assigned loci and shared loci together

        vector<int> uniqueList, overlappingList;
        uniqueList.reserve(len);
        overlappingList.reserve(len);

        size_t i = 0, j = 0, k;

        while (i < len && j < len)
        {
            while (i < len && parent1.chromosome[i] < parent2.chromosome[j])
            {
                uniqueList.push_back(parent1.chromosome[i]);
                ++i;
            }
            while (j < len && (i == len || parent1.chromosome[i] > parent2.chromosome[j]))
            {
                uniqueList.push_back(parent2.chromosome[j]);
                ++j;
            }
            if (j < len && i < len && parent1.chromosome[i] == parent2.chromosome[j])
            {
                overlappingList.push_back(parent2.chromosome[j]);
                ++i, ++j;
            }
        }

        while (i < len)
        {
            uniqueList.push_back(parent1.chromosome[i]);
            ++i;
        }

        while (j < len)
        {
            uniqueList.push_back(parent2.chromosome[j]);
            ++j;
        }

        vector<int> v1;
        v1.reserve(len);
        i = len - overlappingList.size(), j = i, k = 0;
        while (i != 0 && j != 0)
        {
            if (rng::randomBool())
            {
                v1.push_back(uniqueList[k]);
                --i;
            }
            else
            {
                --j;
            }
            k++;
        }

        while (i != 0)
        {
            v1.push_back(uniqueList[k]);
            --i;
            k++;
        }

        merge(v1.begin(), v1.end(), overlappingList.begin(), overlappingList.end(), child1.chromosome.begin());

        if (child1 != parent1)
        {
            child1.is_evaluated = false;
        }

        vector<int>().swap(v1);
        vector<int>().swap(uniqueList);
        vector<int>().swap(overlappingList);

        return child1;
    }

    inline void IntegerGA::evaluate_older(Population &pop)
    {
        assert(fitnessFunction != nullptr);

        std::for_each(std::execution::par_unseq, pop.begin(), pop.end(),
                      [this](Candidate &sol)
                      {
                          double f1 = sol.fitness[0];
                          double f2 = sol.fitness[1];
                          sol.fitness = im_function_(sol.chromosome);

                          double age = f2 + 1;
                          sol.fitness[0] = 1.0 / age * sol.fitness[0] + (age - 1) / age * f1;
                          if (age > max_age_)
                          {
                              sol.fitness[1] = 1.0;
                          }
                          else
                          {
                              sol.fitness[1] = age;
                          }

                          sol.is_evaluated = true;

                          num_fitness_evals_++;
                      });
    }

    inline void IntegerGA::evaluate_new(Population &pop)
    {
        assert(fitnessFunction != nullptr);

        std::for_each(std::execution::par_unseq, pop.begin(), pop.end(),
                      [this](Candidate &sol)
                      {
                          if (!sol.is_evaluated)
                          {
                              sol.fitness = im_function_(sol.chromosome);
                              sol.fitness[1] = 1.0;
                              num_fitness_evals_++;

                              sol.is_evaluated = true;
                          }
                      });
    }

    inline void IntegerGA::mutate(Candidate &child) const
    {
    }

    inline std::vector<std::vector<size_t>> IntegerGA::nonDominatedSort(Population &pop)
    {
        using namespace std;

        /* Calc the number of candidates which dominate each candidate, and the indices of the candidates it dominates. */
        vector<size_t> dom_count(pop.size(), 0);
        vector<vector<bool>> dom_list(pop.size());
        for (size_t i = 0; i < pop.size(); i++)
        {
            dom_list[i].resize(pop.size(), false);
        }

        for (size_t i = 0; i < pop.size(); i++)
        {
            for (size_t j = 0; j < i; j++)
            {
                if (detail::paretoCompare(pop[j].fitness, pop[i].fitness))
                {
                    dom_count[j]++;
                    dom_list[i][j] = true;
                }
                else if (detail::paretoCompare(pop[i].fitness, pop[j].fitness))
                {
                    dom_count[i]++;
                    dom_list[j][i] = true;
                }
            }
        }

        /* Find the indices of all non-dominated candidates (first/best pareto front). */
        vector<size_t> front;
        for (size_t i = 0; i < pop.size(); i++)
        {
            if (dom_count[i] == 0)
            {
                front.push_back(i);
                pop[i].rank = 0;
            }
        }
        /* Find all the other pareto fronts. */
        vector<vector<size_t>> pareto_fronts;
        size_t front_idx = 1;
        while (!front.empty())
        {
            /* "Remove" the current front and find the next one. */
            vector<size_t> next_front;
            for (const auto &i : front)
            {
                for (size_t j = 0; j < pop.size(); j++)
                {
                    /* j belongs to the next front if it's domination count will become 0. */
                    if (dom_list[i][j] && --dom_count[j] == 0)
                    {
                        next_front.push_back(j);
                        pop[j].rank = front_idx;
                    }
                }
            }
            pareto_fronts.push_back(front);
            front = next_front;
            front_idx++;
        }

        for (size_t i = 0; i < pop.size(); i++)
        {
            vector<bool>().swap(dom_list[i]);
        }
        vector<vector<bool>>().swap(dom_list);

        return pareto_fronts;
    }

    inline void IntegerGA::calcCrowdingDistances(Population &pop, std::vector<std::vector<size_t>> &pfronts)
    {
        using namespace std;
        assert(!pop.empty());

        for (const auto &pfront : pfronts)
        {
            for (const auto &idx : pfront)
            {
                pop[idx].distance = 0.0;
            }
        }

        for_each(execution::par_unseq, pfronts.begin(), pfronts.end(),
                 [&pop](vector<size_t> &pfront)
                 {
                     /* Calc the distances in each fitness dimension. */
                     for (size_t d = 0; d < pop[0].fitness.size(); d++)
                     {
                         sort(pfront.begin(), pfront.end(),
                              [&pop, &d](size_t lidx, size_t ridx)
                              {
                                  return pop[lidx].fitness[d] < pop[ridx].fitness[d];
                              });

                         /* Calc the crowding distance for each solution. */
                         double finterval = pop[pfront.back()].fitness[d] - pop[pfront.front()].fitness[d];
                         finterval = max(finterval, 1E-6);

                         pop[pfront.front()].distance = numeric_limits<double>::infinity();
                         pop[pfront.back()].distance = numeric_limits<double>::infinity();
                         for (size_t i = 1; i < pfront.size() - 1; i++)
                         {
                             pop[pfront[i]].distance += (pop[pfront[i + 1]].fitness[d] - pop[pfront[i - 1]].fitness[d]) / finterval;
                         }
                     }
                 });
    }

    inline bool IntegerGA::crowdedCompare(const Candidate &lhs, const Candidate &rhs)
    {
        if (rhs.rank > lhs.rank)
            return true;
        else if (lhs.rank == rhs.rank)
            return lhs.distance > rhs.distance;
        else
            return false;
    }

    inline IntegerGA::Candidate IntegerGA::nsga2Select(const Population &pop)
    {
        assert(!pop.empty());

        size_t idx1 = rng::randomIdx(pop.size());
        size_t idx2 = rng::randomIdx(pop.size());

        return crowdedCompare(pop[idx1], pop[idx2]) ? pop[idx1] : pop[idx2];
    }

    inline IntegerGA::Population IntegerGA::updateNsga2Population(Population &old_pop, CandidateVec &children) const
    {
        using namespace std;
        assert(old_pop.size() == population_size_);
        assert(!children.empty());
        assert(all_of(old_pop.begin(), old_pop.end(), [](const Candidate &sol)
                      { return sol.is_evaluated; }));
        assert(all_of(children.begin(), children.end(), [](const Candidate &sol)
                      { return sol.is_evaluated; }));

        Population new_pop;
        new_pop.reserve(population_size_);

        old_pop.insert(old_pop.end(), make_move_iterator(children.begin()), make_move_iterator(children.end()));
        vector<vector<size_t>> pareto_fronts = nonDominatedSort(old_pop);
        calcCrowdingDistances(old_pop, pareto_fronts);

        /* Add entire fronts while possible. */
        size_t front_idx = 0;
        while (new_pop.size() + pareto_fronts[front_idx].size() <= population_size_)
        {
            for (const auto &idx : pareto_fronts[front_idx])
            {
                new_pop.push_back(move(old_pop[idx]));
            }
            front_idx++;
        }

        /* Add the remaining candidates from the partial front if there is one. */
        if (new_pop.size() != population_size_)
        {
            vector<size_t> added_indices(population_size_ - new_pop.size()); /* For updating the crowding distances in this front. */
            iota(added_indices.begin(), added_indices.end(), new_pop.size());

            vector<size_t> partial_front = pareto_fronts[front_idx];

            sort(partial_front.begin(), partial_front.end(),
                 [&old_pop](size_t lidx, size_t ridx)
                 {
                     return crowdedCompare(old_pop[lidx], old_pop[ridx]);
                 });

            for (const auto &idx : partial_front)
            {
                new_pop.push_back(move(old_pop[idx]));
                if (new_pop.size() == population_size_)
                    break;
            }

            vector<vector<size_t>> temp = {added_indices};
            calcCrowdingDistances(new_pop, temp);
        }

        return new_pop;
    }

} // namespace genetic_algorithm

#endif // !GA_INTEGER_GA_H