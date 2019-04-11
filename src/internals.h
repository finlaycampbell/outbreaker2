#ifndef OUTBREAKER2_INTERNALS_H
#define OUTBREAKER2_INTERNALS_H

#include <Rcpp.h>
#include <RcppEigen.h>

std::vector<int> cpp_are_possible_ancestors(Rcpp::IntegerVector t_inf, size_t i);

size_t cpp_sample1(Rcpp::IntegerVector x);

size_t cpp_pick_possible_ancestor(Rcpp::IntegerVector t_inf, size_t i);

Rcpp::IntegerVector cpp_find_descendents(Rcpp::IntegerVector alpha, size_t i);

Rcpp::IntegerVector cpp_find_local_cases(Rcpp::IntegerVector alpha, size_t i);

Rcpp::List cpp_swap_cases(Rcpp::List param, size_t i, bool swap_place);

size_t cpp_get_n_mutations(Rcpp::List data, size_t i, size_t j);

Rcpp::List cpp_lookup_sequenced_ancestor(Rcpp::List data, Rcpp::List param, size_t i);

void lookup_sequenced_ancestor(Rcpp::IntegerVector alpha, Rcpp::IntegerVector kappa, 
			       Rcpp::LogicalVector has_dna, size_t i, 
			       size_t *out_alpha, 
			       size_t *out_n_generations, 
			       bool *found_sequenced_ancestor
			       );

bool is_between_place(Rcpp::NumericMatrix place_matrix,
		     Rcpp::IntegerVector t_inf,
		     Rcpp::IntegerVector t_onw,
		     Rcpp::IntegerVector alpha,
		     int C_ind,
		     size_t j);

Rcpp::NumericVector get_transition_mat(Rcpp::NumericMatrix p_trans,
				       Rcpp::NumericVector p_place,
				       double eps,
				       double tau,
				       double N_place_unobserved,
				       int max_gamma);

Rcpp::NumericVector get_transition_mat_1(Rcpp::NumericMatrix p_trans,
					 Rcpp::NumericVector p_place,
					 double eps,
					 double N_place_unobserved);

Eigen::MatrixXd calc_trans(Rcpp::NumericMatrix p_trans,
			   Rcpp::NumericVector p_place,
			   double prob,
			   double N_place_unobserved);

Rcpp::NumericMatrix t_inf_change(Rcpp::List data,
				 Rcpp::IntegerVector alpha,
				 Rcpp::IntegerVector kappa,
				 size_t p,
				 size_t t_inf_1,
				 size_t t_inf_2);

Rcpp::NumericMatrix alpha_change(Rcpp::List data,
				 size_t p,
				 size_t kappa,
				 size_t t_inf,
				 size_t alpha_1,
				 size_t alpha_2);

Rcpp::NumericMatrix local_n_contacts(Rcpp::List data,
				     Rcpp::List param,
				     Rcpp::IntegerVector p);

Rcpp::NumericMatrix kappa_change(Rcpp::List data,
				 Rcpp::List param,
				 size_t p,
				 size_t t_inf,
				 size_t alpha,
				 size_t kappa1,
				 size_t kappa2);

Rcpp::NumericMatrix swap_cases_change(Rcpp::List data,
				      Rcpp::List param,
				      Rcpp::List new_param,
				      size_t i,
				      Rcpp::IntegerVector alpha,
				      Rcpp::IntegerVector t_inf,
				      Rcpp::IntegerVector kappa,
				      Rcpp::IntegerVector local_cases,
				      size_t n_mat);

Rcpp::NumericMatrix cpp_find_ancestors(Rcpp::IntegerVector alpha,
				       Rcpp::NumericMatrix ancestors,
				       SEXP i = R_NilValue);

Rcpp::IntegerVector cpp_find_mrca(size_t i,
				  size_t j,
				  Rcpp::NumericMatrix ancestors);

Rcpp::NumericMatrix update_mrca(Rcpp::NumericMatrix combn,
				Rcpp::NumericMatrix ancestors);

#endif
