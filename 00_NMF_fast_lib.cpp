#include <Rcpp.h>
#include <string>
#include <set>
#include <vector>

// [[Rcpp::export]]
Rcpp::IntegerMatrix calculateOverlapOptimized(Rcpp::List input_list) {
  unsigned int nrow = input_list.size();
  std::vector<int> results(nrow * nrow, 0); // 1D vector to store results
  
  // Create a vector of sets to store the gene sets
  std::vector<std::set<std::string>> gene_sets(nrow);
  
  // Convert CharacterVectors to sets before calculating Jaccard
  for (int i = 0; i < nrow; ++i) {
    Rcpp::CharacterVector genes_i = input_list[i];
    gene_sets[i] = std::set<std::string>(genes_i.begin(), genes_i.end());
  }
  
  // Calculate the Jaccard index (or intersection size)
  for (int i = 0; i < nrow; ++i) {
    for (int j = 0; j < nrow; ++j) {
      if (i <= j) {
        // Compute the intersection size
        std::vector<std::string> intersection;
        std::set_intersection(gene_sets[i].begin(), gene_sets[i].end(),
                              gene_sets[j].begin(), gene_sets[j].end(),
                              std::back_inserter(intersection));
        
        int intersection_size = intersection.size();
        // Store the size of the intersection in matrix J
        results[i * nrow + j] = intersection_size;
        results[j * nrow + i] = intersection_size;  // Exploit symmetry
      }
    }
  }
  
  Rcpp::IntegerMatrix J(nrow, nrow);
  std::copy(results.begin(), results.end(), J.begin());
  return J;
}