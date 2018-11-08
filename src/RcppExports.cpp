// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// cpp_are_possible_ancestors
std::vector<int> cpp_are_possible_ancestors(Rcpp::IntegerVector t_inf, size_t i);
RcppExport SEXP _outbreaker2_cpp_are_possible_ancestors(SEXP t_infSEXP, SEXP iSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type t_inf(t_infSEXP);
    Rcpp::traits::input_parameter< size_t >::type i(iSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_are_possible_ancestors(t_inf, i));
    return rcpp_result_gen;
END_RCPP
}
// cpp_sample1
size_t cpp_sample1(std::vector<int> x);
RcppExport SEXP _outbreaker2_cpp_sample1(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<int> >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_sample1(x));
    return rcpp_result_gen;
END_RCPP
}
// cpp_pick_possible_ancestor
size_t cpp_pick_possible_ancestor(Rcpp::IntegerVector t_inf, size_t i);
RcppExport SEXP _outbreaker2_cpp_pick_possible_ancestor(SEXP t_infSEXP, SEXP iSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type t_inf(t_infSEXP);
    Rcpp::traits::input_parameter< size_t >::type i(iSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_pick_possible_ancestor(t_inf, i));
    return rcpp_result_gen;
END_RCPP
}
// cpp_find_descendents
Rcpp::IntegerVector cpp_find_descendents(Rcpp::IntegerVector alpha, size_t i);
RcppExport SEXP _outbreaker2_cpp_find_descendents(SEXP alphaSEXP, SEXP iSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< size_t >::type i(iSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_find_descendents(alpha, i));
    return rcpp_result_gen;
END_RCPP
}
// cpp_find_local_cases
Rcpp::IntegerVector cpp_find_local_cases(Rcpp::IntegerVector alpha, size_t i);
RcppExport SEXP _outbreaker2_cpp_find_local_cases(SEXP alphaSEXP, SEXP iSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< size_t >::type i(iSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_find_local_cases(alpha, i));
    return rcpp_result_gen;
END_RCPP
}
// cpp_swap_cases
Rcpp::List cpp_swap_cases(Rcpp::List param, size_t i, bool swap_ward);
RcppExport SEXP _outbreaker2_cpp_swap_cases(SEXP paramSEXP, SEXP iSEXP, SEXP swap_wardSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type param(paramSEXP);
    Rcpp::traits::input_parameter< size_t >::type i(iSEXP);
    Rcpp::traits::input_parameter< bool >::type swap_ward(swap_wardSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_swap_cases(param, i, swap_ward));
    return rcpp_result_gen;
END_RCPP
}
// cpp_get_n_mutations
size_t cpp_get_n_mutations(Rcpp::List data, size_t i, size_t j);
RcppExport SEXP _outbreaker2_cpp_get_n_mutations(SEXP dataSEXP, SEXP iSEXP, SEXP jSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< size_t >::type i(iSEXP);
    Rcpp::traits::input_parameter< size_t >::type j(jSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_get_n_mutations(data, i, j));
    return rcpp_result_gen;
END_RCPP
}
// cpp_lookup_sequenced_ancestor
Rcpp::List cpp_lookup_sequenced_ancestor(Rcpp::List data, Rcpp::List param, size_t i);
RcppExport SEXP _outbreaker2_cpp_lookup_sequenced_ancestor(SEXP dataSEXP, SEXP paramSEXP, SEXP iSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type param(paramSEXP);
    Rcpp::traits::input_parameter< size_t >::type i(iSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_lookup_sequenced_ancestor(data, param, i));
    return rcpp_result_gen;
END_RCPP
}
// is_between_ward
bool is_between_ward(Rcpp::NumericMatrix ward_matrix, Rcpp::IntegerVector t_inf, Rcpp::IntegerVector t_onw, Rcpp::IntegerVector alpha, int C_ind, size_t j);
RcppExport SEXP _outbreaker2_is_between_ward(SEXP ward_matrixSEXP, SEXP t_infSEXP, SEXP t_onwSEXP, SEXP alphaSEXP, SEXP C_indSEXP, SEXP jSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type ward_matrix(ward_matrixSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type t_inf(t_infSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type t_onw(t_onwSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< int >::type C_ind(C_indSEXP);
    Rcpp::traits::input_parameter< size_t >::type j(jSEXP);
    rcpp_result_gen = Rcpp::wrap(is_between_ward(ward_matrix, t_inf, t_onw, alpha, C_ind, j));
    return rcpp_result_gen;
END_RCPP
}
// get_ward_p
Rcpp::NumericVector get_ward_p(Rcpp::NumericVector p_ward, double eps, double tau, int max_gamma);
RcppExport SEXP _outbreaker2_get_ward_p(SEXP p_wardSEXP, SEXP epsSEXP, SEXP tauSEXP, SEXP max_gammaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type p_ward(p_wardSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< int >::type max_gamma(max_gammaSEXP);
    rcpp_result_gen = Rcpp::wrap(get_ward_p(p_ward, eps, tau, max_gamma));
    return rcpp_result_gen;
END_RCPP
}
// cpp_ll_genetic
double cpp_ll_genetic(Rcpp::List data, Rcpp::List param, SEXP i, Rcpp::RObject custom_function);
RcppExport SEXP _outbreaker2_cpp_ll_genetic(SEXP dataSEXP, SEXP paramSEXP, SEXP iSEXP, SEXP custom_functionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type param(paramSEXP);
    Rcpp::traits::input_parameter< SEXP >::type i(iSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type custom_function(custom_functionSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_ll_genetic(data, param, i, custom_function));
    return rcpp_result_gen;
END_RCPP
}
// cpp_ll_timing_infections
double cpp_ll_timing_infections(Rcpp::List data, Rcpp::List param, SEXP i, Rcpp::RObject custom_function);
RcppExport SEXP _outbreaker2_cpp_ll_timing_infections(SEXP dataSEXP, SEXP paramSEXP, SEXP iSEXP, SEXP custom_functionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type param(paramSEXP);
    Rcpp::traits::input_parameter< SEXP >::type i(iSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type custom_function(custom_functionSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_ll_timing_infections(data, param, i, custom_function));
    return rcpp_result_gen;
END_RCPP
}
// cpp_ll_timing_sampling
double cpp_ll_timing_sampling(Rcpp::List data, Rcpp::List param, SEXP i, Rcpp::RObject custom_function);
RcppExport SEXP _outbreaker2_cpp_ll_timing_sampling(SEXP dataSEXP, SEXP paramSEXP, SEXP iSEXP, SEXP custom_functionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type param(paramSEXP);
    Rcpp::traits::input_parameter< SEXP >::type i(iSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type custom_function(custom_functionSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_ll_timing_sampling(data, param, i, custom_function));
    return rcpp_result_gen;
END_RCPP
}
// cpp_ll_reporting
double cpp_ll_reporting(Rcpp::List data, Rcpp::List param, SEXP i, Rcpp::RObject custom_function);
RcppExport SEXP _outbreaker2_cpp_ll_reporting(SEXP dataSEXP, SEXP paramSEXP, SEXP iSEXP, SEXP custom_functionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type param(paramSEXP);
    Rcpp::traits::input_parameter< SEXP >::type i(iSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type custom_function(custom_functionSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_ll_reporting(data, param, i, custom_function));
    return rcpp_result_gen;
END_RCPP
}
// cpp_ll_contact
double cpp_ll_contact(Rcpp::List data, Rcpp::List param, SEXP i, Rcpp::RObject custom_function);
RcppExport SEXP _outbreaker2_cpp_ll_contact(SEXP dataSEXP, SEXP paramSEXP, SEXP iSEXP, SEXP custom_functionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type param(paramSEXP);
    Rcpp::traits::input_parameter< SEXP >::type i(iSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type custom_function(custom_functionSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_ll_contact(data, param, i, custom_function));
    return rcpp_result_gen;
END_RCPP
}
// cpp_ll_timing
double cpp_ll_timing(Rcpp::List data, Rcpp::List param, SEXP i, Rcpp::RObject custom_functions);
RcppExport SEXP _outbreaker2_cpp_ll_timing(SEXP dataSEXP, SEXP paramSEXP, SEXP iSEXP, SEXP custom_functionsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type param(paramSEXP);
    Rcpp::traits::input_parameter< SEXP >::type i(iSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type custom_functions(custom_functionsSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_ll_timing(data, param, i, custom_functions));
    return rcpp_result_gen;
END_RCPP
}
// cpp_ll_all
double cpp_ll_all(Rcpp::List data, Rcpp::List param, SEXP i, Rcpp::RObject custom_functions);
RcppExport SEXP _outbreaker2_cpp_ll_all(SEXP dataSEXP, SEXP paramSEXP, SEXP iSEXP, SEXP custom_functionsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type param(paramSEXP);
    Rcpp::traits::input_parameter< SEXP >::type i(iSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type custom_functions(custom_functionsSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_ll_all(data, param, i, custom_functions));
    return rcpp_result_gen;
END_RCPP
}
// cpp_move_mu
Rcpp::List cpp_move_mu(Rcpp::List param, Rcpp::List data, Rcpp::List config, Rcpp::RObject custom_ll, Rcpp::RObject custom_prior);
RcppExport SEXP _outbreaker2_cpp_move_mu(SEXP paramSEXP, SEXP dataSEXP, SEXP configSEXP, SEXP custom_llSEXP, SEXP custom_priorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type param(paramSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type config(configSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type custom_ll(custom_llSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type custom_prior(custom_priorSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_move_mu(param, data, config, custom_ll, custom_prior));
    return rcpp_result_gen;
END_RCPP
}
// cpp_move_pi
Rcpp::List cpp_move_pi(Rcpp::List param, Rcpp::List data, Rcpp::List config, Rcpp::RObject custom_ll, Rcpp::RObject custom_prior);
RcppExport SEXP _outbreaker2_cpp_move_pi(SEXP paramSEXP, SEXP dataSEXP, SEXP configSEXP, SEXP custom_llSEXP, SEXP custom_priorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type param(paramSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type config(configSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type custom_ll(custom_llSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type custom_prior(custom_priorSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_move_pi(param, data, config, custom_ll, custom_prior));
    return rcpp_result_gen;
END_RCPP
}
// cpp_move_tau
Rcpp::List cpp_move_tau(Rcpp::List param, Rcpp::List data, Rcpp::List config, Rcpp::RObject custom_ll, Rcpp::RObject custom_prior);
RcppExport SEXP _outbreaker2_cpp_move_tau(SEXP paramSEXP, SEXP dataSEXP, SEXP configSEXP, SEXP custom_llSEXP, SEXP custom_priorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type param(paramSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type config(configSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type custom_ll(custom_llSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type custom_prior(custom_priorSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_move_tau(param, data, config, custom_ll, custom_prior));
    return rcpp_result_gen;
END_RCPP
}
// cpp_move_eps
Rcpp::List cpp_move_eps(Rcpp::List param, Rcpp::List data, Rcpp::List config, Rcpp::RObject custom_ll, Rcpp::RObject custom_prior);
RcppExport SEXP _outbreaker2_cpp_move_eps(SEXP paramSEXP, SEXP dataSEXP, SEXP configSEXP, SEXP custom_llSEXP, SEXP custom_priorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type param(paramSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type config(configSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type custom_ll(custom_llSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type custom_prior(custom_priorSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_move_eps(param, data, config, custom_ll, custom_prior));
    return rcpp_result_gen;
END_RCPP
}
// cpp_move_lambda
Rcpp::List cpp_move_lambda(Rcpp::List param, Rcpp::List data, Rcpp::List config, Rcpp::RObject custom_ll, Rcpp::RObject custom_prior);
RcppExport SEXP _outbreaker2_cpp_move_lambda(SEXP paramSEXP, SEXP dataSEXP, SEXP configSEXP, SEXP custom_llSEXP, SEXP custom_priorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type param(paramSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type config(configSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type custom_ll(custom_llSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type custom_prior(custom_priorSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_move_lambda(param, data, config, custom_ll, custom_prior));
    return rcpp_result_gen;
END_RCPP
}
// cpp_move_t_inf
Rcpp::List cpp_move_t_inf(Rcpp::List param, Rcpp::List data, Rcpp::RObject list_custom_ll);
RcppExport SEXP _outbreaker2_cpp_move_t_inf(SEXP paramSEXP, SEXP dataSEXP, SEXP list_custom_llSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type param(paramSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type list_custom_ll(list_custom_llSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_move_t_inf(param, data, list_custom_ll));
    return rcpp_result_gen;
END_RCPP
}
// cpp_move_alpha
Rcpp::List cpp_move_alpha(Rcpp::List param, Rcpp::List data, Rcpp::RObject list_custom_ll);
RcppExport SEXP _outbreaker2_cpp_move_alpha(SEXP paramSEXP, SEXP dataSEXP, SEXP list_custom_llSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type param(paramSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type list_custom_ll(list_custom_llSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_move_alpha(param, data, list_custom_ll));
    return rcpp_result_gen;
END_RCPP
}
// cpp_move_model
Rcpp::List cpp_move_model(Rcpp::List param, Rcpp::List data, Rcpp::List config, Rcpp::RObject list_custom_ll);
RcppExport SEXP _outbreaker2_cpp_move_model(SEXP paramSEXP, SEXP dataSEXP, SEXP configSEXP, SEXP list_custom_llSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type param(paramSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type config(configSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type list_custom_ll(list_custom_llSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_move_model(param, data, config, list_custom_ll));
    return rcpp_result_gen;
END_RCPP
}
// cpp_move_joint
Rcpp::List cpp_move_joint(Rcpp::List param, Rcpp::List data, Rcpp::List config, Rcpp::RObject list_custom_ll);
RcppExport SEXP _outbreaker2_cpp_move_joint(SEXP paramSEXP, SEXP dataSEXP, SEXP configSEXP, SEXP list_custom_llSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type param(paramSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type config(configSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type list_custom_ll(list_custom_llSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_move_joint(param, data, config, list_custom_ll));
    return rcpp_result_gen;
END_RCPP
}
// cpp_move_swap_cases
Rcpp::List cpp_move_swap_cases(Rcpp::List param, Rcpp::List data, Rcpp::RObject list_custom_ll);
RcppExport SEXP _outbreaker2_cpp_move_swap_cases(SEXP paramSEXP, SEXP dataSEXP, SEXP list_custom_llSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type param(paramSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type list_custom_ll(list_custom_llSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_move_swap_cases(param, data, list_custom_ll));
    return rcpp_result_gen;
END_RCPP
}
// cpp_move_kappa
Rcpp::List cpp_move_kappa(Rcpp::List param, Rcpp::List data, Rcpp::List config, Rcpp::RObject list_custom_ll);
RcppExport SEXP _outbreaker2_cpp_move_kappa(SEXP paramSEXP, SEXP dataSEXP, SEXP configSEXP, SEXP list_custom_llSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type param(paramSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type config(configSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type list_custom_ll(list_custom_llSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_move_kappa(param, data, config, list_custom_ll));
    return rcpp_result_gen;
END_RCPP
}
// cpp_prior_mu
double cpp_prior_mu(Rcpp::List param, Rcpp::List config, Rcpp::RObject custom_function);
RcppExport SEXP _outbreaker2_cpp_prior_mu(SEXP paramSEXP, SEXP configSEXP, SEXP custom_functionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type param(paramSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type config(configSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type custom_function(custom_functionSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_prior_mu(param, config, custom_function));
    return rcpp_result_gen;
END_RCPP
}
// cpp_prior_pi
double cpp_prior_pi(Rcpp::List param, Rcpp::List config, Rcpp::RObject custom_function);
RcppExport SEXP _outbreaker2_cpp_prior_pi(SEXP paramSEXP, SEXP configSEXP, SEXP custom_functionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type param(paramSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type config(configSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type custom_function(custom_functionSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_prior_pi(param, config, custom_function));
    return rcpp_result_gen;
END_RCPP
}
// cpp_prior_tau
double cpp_prior_tau(Rcpp::List param, Rcpp::List config, Rcpp::RObject custom_function);
RcppExport SEXP _outbreaker2_cpp_prior_tau(SEXP paramSEXP, SEXP configSEXP, SEXP custom_functionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type param(paramSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type config(configSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type custom_function(custom_functionSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_prior_tau(param, config, custom_function));
    return rcpp_result_gen;
END_RCPP
}
// cpp_prior_eps
double cpp_prior_eps(Rcpp::List param, Rcpp::List config, Rcpp::RObject custom_function);
RcppExport SEXP _outbreaker2_cpp_prior_eps(SEXP paramSEXP, SEXP configSEXP, SEXP custom_functionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type param(paramSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type config(configSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type custom_function(custom_functionSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_prior_eps(param, config, custom_function));
    return rcpp_result_gen;
END_RCPP
}
// cpp_prior_lambda
double cpp_prior_lambda(Rcpp::List param, Rcpp::List config, Rcpp::RObject custom_function);
RcppExport SEXP _outbreaker2_cpp_prior_lambda(SEXP paramSEXP, SEXP configSEXP, SEXP custom_functionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type param(paramSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type config(configSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type custom_function(custom_functionSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_prior_lambda(param, config, custom_function));
    return rcpp_result_gen;
END_RCPP
}
// cpp_prior_all
double cpp_prior_all(Rcpp::List param, Rcpp::List config, Rcpp::RObject custom_functions);
RcppExport SEXP _outbreaker2_cpp_prior_all(SEXP paramSEXP, SEXP configSEXP, SEXP custom_functionsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type param(paramSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type config(configSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type custom_functions(custom_functionsSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_prior_all(param, config, custom_functions));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_outbreaker2_cpp_are_possible_ancestors", (DL_FUNC) &_outbreaker2_cpp_are_possible_ancestors, 2},
    {"_outbreaker2_cpp_sample1", (DL_FUNC) &_outbreaker2_cpp_sample1, 1},
    {"_outbreaker2_cpp_pick_possible_ancestor", (DL_FUNC) &_outbreaker2_cpp_pick_possible_ancestor, 2},
    {"_outbreaker2_cpp_find_descendents", (DL_FUNC) &_outbreaker2_cpp_find_descendents, 2},
    {"_outbreaker2_cpp_find_local_cases", (DL_FUNC) &_outbreaker2_cpp_find_local_cases, 2},
    {"_outbreaker2_cpp_swap_cases", (DL_FUNC) &_outbreaker2_cpp_swap_cases, 3},
    {"_outbreaker2_cpp_get_n_mutations", (DL_FUNC) &_outbreaker2_cpp_get_n_mutations, 3},
    {"_outbreaker2_cpp_lookup_sequenced_ancestor", (DL_FUNC) &_outbreaker2_cpp_lookup_sequenced_ancestor, 3},
    {"_outbreaker2_is_between_ward", (DL_FUNC) &_outbreaker2_is_between_ward, 6},
    {"_outbreaker2_get_ward_p", (DL_FUNC) &_outbreaker2_get_ward_p, 4},
    {"_outbreaker2_cpp_ll_genetic", (DL_FUNC) &_outbreaker2_cpp_ll_genetic, 4},
    {"_outbreaker2_cpp_ll_timing_infections", (DL_FUNC) &_outbreaker2_cpp_ll_timing_infections, 4},
    {"_outbreaker2_cpp_ll_timing_sampling", (DL_FUNC) &_outbreaker2_cpp_ll_timing_sampling, 4},
    {"_outbreaker2_cpp_ll_reporting", (DL_FUNC) &_outbreaker2_cpp_ll_reporting, 4},
    {"_outbreaker2_cpp_ll_contact", (DL_FUNC) &_outbreaker2_cpp_ll_contact, 4},
    {"_outbreaker2_cpp_ll_timing", (DL_FUNC) &_outbreaker2_cpp_ll_timing, 4},
    {"_outbreaker2_cpp_ll_all", (DL_FUNC) &_outbreaker2_cpp_ll_all, 4},
    {"_outbreaker2_cpp_move_mu", (DL_FUNC) &_outbreaker2_cpp_move_mu, 5},
    {"_outbreaker2_cpp_move_pi", (DL_FUNC) &_outbreaker2_cpp_move_pi, 5},
    {"_outbreaker2_cpp_move_tau", (DL_FUNC) &_outbreaker2_cpp_move_tau, 5},
    {"_outbreaker2_cpp_move_eps", (DL_FUNC) &_outbreaker2_cpp_move_eps, 5},
    {"_outbreaker2_cpp_move_lambda", (DL_FUNC) &_outbreaker2_cpp_move_lambda, 5},
    {"_outbreaker2_cpp_move_t_inf", (DL_FUNC) &_outbreaker2_cpp_move_t_inf, 3},
    {"_outbreaker2_cpp_move_alpha", (DL_FUNC) &_outbreaker2_cpp_move_alpha, 3},
    {"_outbreaker2_cpp_move_model", (DL_FUNC) &_outbreaker2_cpp_move_model, 4},
    {"_outbreaker2_cpp_move_joint", (DL_FUNC) &_outbreaker2_cpp_move_joint, 4},
    {"_outbreaker2_cpp_move_swap_cases", (DL_FUNC) &_outbreaker2_cpp_move_swap_cases, 3},
    {"_outbreaker2_cpp_move_kappa", (DL_FUNC) &_outbreaker2_cpp_move_kappa, 4},
    {"_outbreaker2_cpp_prior_mu", (DL_FUNC) &_outbreaker2_cpp_prior_mu, 3},
    {"_outbreaker2_cpp_prior_pi", (DL_FUNC) &_outbreaker2_cpp_prior_pi, 3},
    {"_outbreaker2_cpp_prior_tau", (DL_FUNC) &_outbreaker2_cpp_prior_tau, 3},
    {"_outbreaker2_cpp_prior_eps", (DL_FUNC) &_outbreaker2_cpp_prior_eps, 3},
    {"_outbreaker2_cpp_prior_lambda", (DL_FUNC) &_outbreaker2_cpp_prior_lambda, 3},
    {"_outbreaker2_cpp_prior_all", (DL_FUNC) &_outbreaker2_cpp_prior_all, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_outbreaker2(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
