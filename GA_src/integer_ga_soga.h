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

#include "base_ga.h"
#include "../ga_operators_copy.hpp"

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

    private:
        CrossoverMethod crossover_method_ = CrossoverMethod::set_based;
        size_t num_crossover_points_ = 3;
        size_t base_ = 4;
        double swap_rate_ = 0.1;
        double inversion_rate_ = 0.1;
        double child_rate_ = 0.0;

        Candidate generateCandidate() const override;
        CandidatePair crossover(const Candidate &parent1, const Candidate &parent2) const override;
        inline void repair(Population &pop);
        inline void mutate(Candidate &child) const;
        inline Population generateInitialPopulation() const;
        inline bool stopCondition() const;

        static CandidatePair setCrossover(const Candidate &parent1, const Candidate &parent2, double pc);

        inline void evaluate_new(Population &pop);

        IMP im_function_;

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

#include "rng.h"

namespace genetic_algorithm
{
    inline IntegerGA::IntegerGA(size_t chrom_len, IMP im_function, size_t base)
        : GA(chrom_len, im_function), base_(base), im_function_(im_function)
    {
        if (base < 2)
            throw std::invalid_argument("The base must be at least 2.");
    }

    inline IntegerGA::CandidateVec IntegerGA::run()
    {
        using namespace std;

        /* Init. */
        num_fitness_evals_ = 0;
        solutions_.clear();
        population_.clear();
        soga_history_.clear();
        soga_history_.reserve(max_gen_);

        for (vector<int>::const_iterator it = g_hg.getNodeMapVec().cbegin(); it != g_hg.getNodeMapVec().cend(); it++)
            g_nodeDegree[*it - 1] = g_hg.getNumEdgeNode(*it);

        population_.reserve(population_size_);
        for (size_t i = 0; i < std::min(population_size_, initial_population_preset_.size()); i++)
        {
            population_.push_back(initial_population_preset_[i]);
        }
        while (population_.size() < population_size_)
        {
            vector<int> chrom;
            chrom.reserve(chrom_len_);
            while (chrom.size() < chrom_len_)
            {
                size_t e_idx = rng::randomIdx(g_hg.getNumEdge());
                int nodeIdx = g_hg.getEdge(e_idx)[rng::randomIdx(g_hg.getEdge(e_idx).size())] - 1;
                if (!binary_search(chrom.begin(), chrom.end(), nodeIdx))
                {
                    chrom.push_back(nodeIdx);
                    sort(chrom.begin(), chrom.end());
                }
            }
            Candidate c(chrom);
            c.fitness = {0};
            c.is_evaluated = false;
            population_.push_back(c);
        }

        evaluate_new(population_);
        updateStats(population_);
        /* Other generations. */
        generation_cntr_ = 1;
        // printf("gen : %2.ld\r", generation_cntr_);
        // fflush(stdout);
        size_t num_children = population_size_ + population_size_ % 2;
        vector<CandidatePair> parent_pairs(num_children / 2);
        vector<Candidate> children;
        children.reserve(num_children);

        // appearence map for early stop
        std::map<string, double> ht;

        while (generation_cntr_ < max_gen_)
        {
            /* Selections. */
            generate(execution::par_unseq, parent_pairs.begin(), parent_pairs.end(),
                     [this]() -> CandidatePair
                     {
                         return move(make_pair(select(population_), select(population_)));
                     });

            /* Crossovers. */
            for_each(execution::par_unseq, parent_pairs.begin(), parent_pairs.end(),
                     [this](CandidatePair &p) -> void
                     {
                         p = crossover(p.first, p.second);
                     });

            children.clear();
            for (size_t i = 0; i < parent_pairs.size(); i++)
            {
                children.push_back(move(parent_pairs[i].first));
                children.push_back(move(parent_pairs[i].second));
            }

            /* Apply repair function to the children if set. */
            repair(children);

            /* Overwrite the current population with the children. */
            evaluate_new(children);
            // population_ = updatePopulation(population_, children);
            std::partial_sort(population_.begin(), population_.begin() + 1, population_.end(),
                              [](const Candidate &lhs, const Candidate &rhs)
                              {
                                  return lhs.fitness[0] > rhs.fitness[0];
                              });
            std::partial_sort(children.begin(), children.begin() + 1, children.end(),
                              [](const Candidate &lhs, const Candidate &rhs)
                              {
                                  return lhs.fitness[0] < rhs.fitness[0];
                              });
            children[0] = population_[0];
            population_ = children;

            generation_cntr_++;
            updateStats(population_);
            // printf("gen : %2.ld\r", generation_cntr_);
            // fflush(stdout);

            // early stop
            cout << "gen : " << generation_cntr_ << " : " << soga_history_.fitness_max.back() << endl;
        }
        printf("\npop size = %ld\n", population_.size());
        cout << generation_cntr_ << " generation " << endl;

        std::map<string, double>().swap(ht);
        vector<CandidatePair>().swap(parent_pairs);
        vector<Candidate>().swap(children);
        malloc_trim(0);

        return solutions_;
    }

    // ~IntegerGA()
    // {
    //     /* release memory */
    // }

    inline void IntegerGA::repair(Population &pop)
    {
        size_t tmp = chrom_len_ / 2;
        if (tmp == 0)
            tmp = 1;
        vector<size_t> chosen(tmp);
        for (size_t i = 0; i < tmp; i++)
        {
            chosen[i] = rng::randomIdx(chrom_len_);
        }

        size_t init_size = pop.size();
        // if (rng::randomReal(0.0, 1.0) > 0.75)
        // {
        //     pop.resize(init_size * 2);
        //     for (size_t i = 0; i < init_size; i++)
        //     {
        //         // printf("%d %o\t",i , &(pop[i + init_size].chromosome));
        //         pop[i + init_size] = Candidate(pop[i].chromosome);
        //         pop[i + init_size].fitness = vector<double>(1, 0);
        //     }
        // }
        for (size_t i = 0; i < init_size; ++i)
        {
            for (size_t j = 0; j < tmp; j++)
            {
                size_t e_idx = rng::randomIdx(g_hg.getNumEdge());
                int nodeIdx = g_hg.getEdge(e_idx)[rng::randomIdx(g_hg.getEdge(e_idx).size())] - 1;
                if (!binary_search(pop[i].chromosome.begin(), pop[i].chromosome.end(), nodeIdx))
                {
                    pop[i].chromosome[chosen[j]] = nodeIdx;
                    sort(pop[i].chromosome.begin(), pop[i].chromosome.end());
                }
            }
            pop[i].is_evaluated = false;
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

        switch (crossover_method_)
        {
        case CrossoverMethod::set_based:
            return setCrossover(parent1, parent2, crossover_rate_);
        case CrossoverMethod::custom:
            return customCrossover(parent1, parent2, crossover_rate_);
        default:
            assert(false); /* Invalid crossover method. Shouldn't get here. */
            std::abort();
        }
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
                while (j < len && parent1.chromosome[i] > parent2.chromosome[j])
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

    inline void IntegerGA::evaluate_new(Population &pop)
    {
        assert(fitnessFunction != nullptr);

        std::for_each(std::execution::par_unseq, pop.begin(), pop.end(),
                      [this](Candidate &sol)
                      {
                          sol.fitness[0] = im_function_(sol.chromosome)[0];
                          sol.is_evaluated = true;
                      });
        if (pop.size() > 1600)
            cout << "invalid " << endl;
    }

    inline IntegerGA::Population IntegerGA::generateInitialPopulation() const
    {
        assert(population_size_ > 0);

        Population pop;
        pop.reserve(population_size_);

        for (size_t i = 0; i < std::min(population_size_, initial_population_preset_.size()); i++)
        {
            pop.push_back(initial_population_preset_[i]);
        }
        while (pop.size() < population_size_)
        {
            Candidate c = generateCandidate();
            c.fitness = {0, 0};
            pop.push_back(c);
        }

        return pop;
    }

    inline bool IntegerGA::stopCondition() const
    {
        /* Always stop when reaching max_gen regardless of stop condition. */
        if (generation_cntr_ >= max_gen_ - 1)
            return true;

        /* Early-stop conditions. */
        double metric_now, metric_old;
        switch (stop_condition_)
        {
        case StopCondition::max_gen:
            /* Already checked above. */
            return false;

        case StopCondition::fitness_value:
            return std::any_of(population_.begin(), population_.end(),
                               [this](const Candidate &sol)
                               {
                                   return detail::paretoCompare(fitness_reference_, sol.fitness);
                               });

        case StopCondition::fitness_evals:
            return num_fitness_evals_ >= max_fitness_evals_;

        case StopCondition::fitness_mean_stall:
            if (generation_cntr_ >= stall_gen_count_)
            {
                metric_now = soga_history_.fitness_mean[generation_cntr_];
                metric_old = soga_history_.fitness_mean[generation_cntr_ - stall_gen_count_];

                return (metric_now - metric_old) < stall_threshold_;
            }
            else
                return false;

        case StopCondition::fitness_best_stall:
            if (generation_cntr_ >= stall_gen_count_)
            {
                metric_now = soga_history_.fitness_mean[generation_cntr_];
                metric_old = soga_history_.fitness_mean[generation_cntr_ - stall_gen_count_];

                return (metric_now - metric_old) < stall_threshold_;
            }
            else
                return false;

        default:
            assert(false); /* Invalid stop condition. Shouldn't get here. */
            std::abort();
        }
    }

    inline void IntegerGA::mutate(Candidate &child) const
    {
    }

} // namespace genetic_algorithm

#endif // !GA_INTEGER_GA_H