#ifndef cSubpopulation_h
#define cSubpopulation_h

#include "common.h"

namespace bpopsim {

  // Predef
  class SimulationParameters;
  
  // Defined
  struct cGenotype 
  {
    uint32_t  unique_node_id;
    string    name;
    double    fitness;
    double this_mutation_fitness_effect; // Size of the mutation that created this genotype
    uint32_t  total_mutation_count;   // total number of mutations
    vector<uint32_t> mutation_counts; // number of mutations in each category
    bool      marked_for_deletion;      // do individuals of this type still exist?
    
    cGenotype(
              SimulationParameters& in_simulation_parameters, 
              uint32_t in_node_id
              );
    
    // Constructor for children, initializes wrt ancestor correctly.
    cGenotype(
              SimulationParameters& in_simulation_parameters, 
              uint32_t in_node_id,
              cGenotype& in_ancestor
              );
    
    // Applies one (random) mutation to genotype.
    bool AddOneMutation(
                        SimulationParameters& in_simulation_parameters,
                        gsl_rng * rng
                        );
    
  };

  struct cFrequencySlice {
    cFrequencySlice(uint32_t node_id, double in_low, double in_high) : unique_node_id(node_id), low(in_low), high(in_high) {}
    uint32_t unique_node_id;
    double low;
    double high;
  };

  struct cSortByLow {
    bool operator () (const cFrequencySlice & lhs , const cFrequencySlice & rhs) const {
      return lhs.low < rhs.low;
    }
  };

  class cSubpopulation {

  private:
    double m_number;
    char   m_marker;
    tree<cGenotype>::iterator m_genotype;

  public:
    cSubpopulation() 
    : m_number(0)
    , m_marker('n') 
    { };
    
    cSubpopulation(const cSubpopulation& in); // Copy constructor
    
    virtual ~cSubpopulation() { ; }; 

    const double GetFitness() { return m_genotype->fitness; }
    const string GetName() { return m_genotype->name; }
    tree<cGenotype>::iterator GetGenotypeIter() { return m_genotype; }
    cGenotype& GetGenotype() { return *m_genotype; }
    const uint32_t GetNode_id() { return m_genotype->unique_node_id; }
    const double GetNumber() { return m_number; }
    const double GetIntegralNumber() { return floor(m_number); }
    const char GetMarker() { return m_marker; }
    const uint32_t GetTotalMutationCount() { return m_genotype->total_mutation_count; }
    const uint32_t GetMutationCountInCategory(size_t in_category) { return m_genotype->mutation_counts[in_category]; }
    
    void SetGenotype(tree<cGenotype>::iterator location) { m_genotype = location; }
    void SetNumber(const double in_number) { m_number = in_number; }
    void SetMarker(const char in_marker) { m_marker = in_marker; } 
    void Transfer(double success_prob, gsl_rng * randomgenerator);
    double MutantFitness(double in_fitness, 
                         double in_average_mutation_s, 
                         const string& in_type_of_mutation, 
                         gsl_rng * randomgenerator);  
    virtual bool CreateDescendant(gsl_rng * randomgenerator,
                                  cSubpopulation &ancestor, 
                                  SimulationParameters& simulation_parameters, 
                                  tree<cGenotype>& in_tree, 
                                  uint32_t node_id);
    
    virtual void AddToTree(tree<cGenotype>& in_tree,
                           tree<cGenotype>::iterator parent,
                           cGenotype child);
  };
  
}
  
#endif
