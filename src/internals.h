#ifndef OUTBREAKER2_INTERNALS_H
#define OUTBREAKER2_INTERNALS_H

#include <Rcpp.h>

std::vector<int> cpp_are_possible_ancestors(Rcpp::IntegerVector t_inf, size_t i);

size_t cpp_sample1(Rcpp::IntegerVector x);

size_t cpp_pick_possible_ancestor(Rcpp::IntegerVector t_inf, size_t i);

Rcpp::IntegerVector cpp_find_descendents(Rcpp::IntegerVector alpha, size_t i);

Rcpp::IntegerVector cpp_find_local_cases(Rcpp::IntegerVector alpha, size_t i);

Rcpp::List cpp_swap_cases(Rcpp::List param, size_t i);

size_t cpp_get_n_mutations(Rcpp::List data, size_t i, size_t j);

Rcpp::List cpp_lookup_sequenced_ancestor(Rcpp::List data, Rcpp::List param, size_t i);

void lookup_sequenced_ancestor(Rcpp::IntegerVector alpha, Rcpp::IntegerVector kappa, 
			       Rcpp::LogicalVector has_dna, size_t i, 
			       size_t *out_alpha, 
			       size_t *out_n_generations, 
			       bool *found_sequenced_ancestor
			       );

bool is_between_ward(Rcpp::NumericMatrix ward_matrix,
		     Rcpp::IntegerVector t_inf,
		     Rcpp::IntegerVector t_onw,
		     Rcpp::IntegerVector alpha,
		     int C_ind,
		     size_t j);

#endif
