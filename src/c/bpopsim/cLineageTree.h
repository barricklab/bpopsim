#ifndef cLineageTree_h
#define cLineageTree_h

//#include "cPopulation.h"

class cPopulation;
#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip> 

class cLineageTree
{
   private:
     std::vector<long double> m_tree;
     std::vector<int> m_lineages;
     std::vector<long double> m_frequencies;
   
   public:
     cLineageTree() {;};
     
   virtual ~cLineageTree() {;};

   void AddNode(long double in_mutation){ m_tree.push_back(in_mutation);}
   void SetLineage(int in_previous){m_lineages.push_back(in_previous);}
   void ClearRuns();
   void PrintTree(const std::string& output_file_name);
   void CalculateFrequencies(cPopulation& in, const std::string& output_file_name);     
      

   std::vector<long double> GetRuns() { return m_tree;}
   std::vector<int> GetLineages() { return m_lineages;}
   int GetSizeOfTree() {return m_tree.size();}
   
   
}; 
#endif
