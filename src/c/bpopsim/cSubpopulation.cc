#include "common.h"
#include "cSubpopulation.h"
#include "cPopulation.h"

using namespace std;

namespace bpopsim {

cGenotype::cGenotype(
          SimulationParameters& in_simulation_parameters, 
          uint32_t in_node_id
          ) 
: unique_node_id(in_node_id)
, name("")
, fitness(1.0)
, this_mutation_fitness_effect(0.0)
, total_mutation_count(0)
, marked_for_deletion(false)
{ 
  // set mutation counts vector to be as large as the number of categories
  mutation_counts.resize(in_simulation_parameters.mutation_rates_per_division.size(), 0); 
};
  
cGenotype::cGenotype(
                     SimulationParameters& in_simulation_parameters, 
                     uint32_t in_node_id,
                     cGenotype& in_ancestor
                     ) 
: unique_node_id(in_node_id)
, name("")
, fitness(in_ancestor.fitness)
, this_mutation_fitness_effect(0.0)
, total_mutation_count(in_ancestor.total_mutation_count)
, mutation_counts(in_ancestor.mutation_counts)
, marked_for_deletion(false)
{ 
  // set mutation counts vector to be as large as the number of categories
  mutation_counts.resize(in_simulation_parameters.mutation_rates_per_division.size(), 0); 
};
  
bool cGenotype::AddOneMutation(
                               SimulationParameters& simulation_parameters,
                               gsl_rng * rng
                               )
{
  // @JEB: more of this could be moved to SimulationParameters object
  
  // override fitness effects of first mutations
  if (this->total_mutation_count < simulation_parameters.first_mutation_fitness_effects.size()) {
    this->this_mutation_fitness_effect = simulation_parameters.first_mutation_fitness_effects[this->total_mutation_count];
    
    // Note: we increment only the total number of mutations, not the number in any given category
  }
  else { // default behavior
  
    uint32_t this_mutation_category = simulation_parameters.DetermineMutationCategory(rng);
    double average_mutation_fitness_effect = simulation_parameters.mutation_fitness_effects[this_mutation_category];
    
    if(simulation_parameters.mutation_fitness_effect_model=="e")
    {
      this->this_mutation_fitness_effect = gsl_ran_exponential(rng,average_mutation_fitness_effect);
    }
      
    if (simulation_parameters.mutation_fitness_effect_model=="u")
    {
      this->this_mutation_fitness_effect = average_mutation_fitness_effect;
    }
    
    // One-time mutations
    if (simulation_parameters.mutation_fitness_effect_model=="o")
    {
      if (this->mutation_counts[this_mutation_category] < 1) {
        this->this_mutation_fitness_effect = average_mutation_fitness_effect;
      } else {
        return false; // short circuit = do not add population
      }
    }
    
    this->mutation_counts[this_mutation_category]++;
  }
  
  // Update our information
  this->total_mutation_count++;
  this->fitness += this->this_mutation_fitness_effect;
  
  return true;
}
  
/* Copy constructor */
cSubpopulation::cSubpopulation(const cSubpopulation& in) :
  m_number(0),
  m_marker(0),
  m_genotype(0) { 
  m_number = in.m_number;
  m_marker = in.m_marker;
  m_genotype = in.m_genotype;
};

bool cSubpopulation::CreateDescendant(  gsl_rng * randomgenerator,
                                        cSubpopulation &ancestor, 
                                        SimulationParameters& simulation_parameters,
                                        tree<cGenotype>& in_tree, 
                                        uint32_t node_id
                                      ) 
{	
  if (g_verbose) cout << "Creating descendant with node_id: " << node_id << endl;
	
	cGenotype new_genotype(simulation_parameters, node_id, ancestor.GetGenotype());
  bool success_adding = new_genotype.AddOneMutation(simulation_parameters, randomgenerator);
  if (!success_adding) return false;
  
  tree<cGenotype>::iterator_base new_geno_it;
  // There is only one new one...
  SetNumber(1.0);
  
  // ...taken from the ancestor.
  ancestor.SetNumber(ancestor.GetNumber()-1.0);
  
  SetMarker(ancestor.GetMarker());
  
  new_geno_it = in_tree.append_child(ancestor.m_genotype, new_genotype); 
	SetGenotype(new_geno_it);
  return true;
}

void cSubpopulation::AddToTree(tree<cGenotype> &in_tree, 
                               tree<cGenotype>::iterator parent, 
                               cGenotype child) {
  in_tree.append_child(parent, child);
}

void cSubpopulation::Transfer(double success_prob, 
                              gsl_rng * randomgenerator) {
		
	uint32_t random_gsl_int = gsl_ran_binomial(randomgenerator, 
                                             success_prob, 
                                             GetIntegralNumber());
   
	m_number = random_gsl_int;
}

}
