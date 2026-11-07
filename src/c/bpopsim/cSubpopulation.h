#pragma once

#include "common.h"

namespace bpopsim {

  // Standard library names used unqualified throughout the bpopsim headers.
  using std::map;
  using std::pair;
  using std::set;
  using std::string;
  using std::vector;

  // Predef
  class SimulationParameters;

  // Defined
  struct cGenotype
  {
    uint32_t unique_node_id;
    string name;
    double fitness;
    double this_mutation_fitness_effect;  // Size of the mutation that created this genotype
    uint32_t total_mutation_count;        // total number of mutations
    vector<uint32_t> mutation_counts;     // number of mutations in each category
    bool marked_for_deletion;             // do individuals of this type still exist?

    cGenotype(SimulationParameters& in_simulation_parameters, uint32_t in_node_id);

    // Constructor for children, initializes wrt ancestor correctly.
    cGenotype(SimulationParameters& in_simulation_parameters, uint32_t in_node_id, cGenotype& in_ancestor);

    // Applies one (random) mutation to genotype.
    bool AddOneMutation(SimulationParameters& in_simulation_parameters, RandomGenerator& rng);
  };

  struct cFrequencySlice
  {
    cFrequencySlice(uint32_t node_id, double in_low, double in_high)
    : unique_node_id(node_id)
    , low(in_low)
    , high(in_high)
    {}
    uint32_t unique_node_id;
    double low;
    double high;
  };

  struct cSortByLow
  {
    bool operator()(const cFrequencySlice& lhs, const cFrequencySlice& rhs) const { return lhs.low < rhs.low; }
  };

  class cSubpopulation
  {

  private:
    double m_number{0};
    char m_marker{'n'};
    tree<cGenotype>::iterator m_genotype;

  public:
    cSubpopulation() = default;
    virtual ~cSubpopulation() = default;

    double GetFitness() { return m_genotype->fitness; }
    string GetName() { return m_genotype->name; }
    tree<cGenotype>::iterator GetGenotypeIter() { return m_genotype; }
    cGenotype& GetGenotype() { return *m_genotype; }
    uint32_t GetNode_id() { return m_genotype->unique_node_id; }
    double GetNumber() const { return m_number; }
    double GetIntegralNumber() const { return floor(m_number); }
    char GetMarker() const { return m_marker; }
    uint32_t GetTotalMutationCount() { return m_genotype->total_mutation_count; }
    uint32_t GetMutationCountInCategory(size_t in_category) { return m_genotype->mutation_counts[in_category]; }

    void SetGenotype(tree<cGenotype>::iterator location) { m_genotype = location; }
    void SetNumber(const double in_number) { m_number = in_number; }
    void SetMarker(const char in_marker) { m_marker = in_marker; }
    void Transfer(double success_prob, RandomGenerator& randomgenerator);
    virtual bool CreateDescendant(RandomGenerator& randomgenerator, cSubpopulation& ancestor,
                                  SimulationParameters& simulation_parameters, tree<cGenotype>& in_tree,
                                  uint32_t node_id, bool verbose = false);

    virtual void AddToTree(tree<cGenotype>& in_tree, tree<cGenotype>::iterator parent, cGenotype child);
  };

}  // namespace bpopsim
