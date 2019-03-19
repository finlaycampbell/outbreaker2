#include <Rcpp.h>
#include <Rmath.h>
#include "internals.h"
#include "likelihoods.h"
#include "priors.h"



// IMPORTANT: ON INDEXING VECTORS AND ANCESTRIES

// Most of the functions implemented here are susceptible to be called from R
// via Rcpp, and are therefore treated as interfaces. This causes a number of
// headaches when using indices of cases defined in R (1:N) to refer to elements
// in Rcpp / Cpp vectors (0:N-1). By convention, we store all data on the
// original scale (1:N), and modify indices whenever accessing elements of
// vectors. In other words, in an expression like 'alpha[j]', 'j' should always
// be on the internal scale (0:N-1).

// In all these functions, 'SEXP i' is an optional vector of case indices, on
// the 1:N scale.






// ---------------------------

// Movement of the mutation rate 'mu' is done using a dumb normal proposal. This
// is satisfying for now - we only reject a few non-sensical values outside the
// range [0;1]. The SD of the proposal (implicitely contained in rand$mu.rnorm1,
// but really provided through 'config', seems fine as the range of real values
// will never change much. Probably not much point in using auto-tuning here.

// [[Rcpp::export(rng = true)]]
Rcpp::List cpp_move_mu(Rcpp::List param, Rcpp::List data, Rcpp::List config,
		       Rcpp::RObject custom_ll = R_NilValue,
		       Rcpp::RObject custom_prior = R_NilValue) {

  // deep copy here for now, ultimately should be an arg.
  Rcpp::List new_param = clone(param);
  Rcpp::NumericVector mu = param["mu"];
  Rcpp::NumericVector new_mu = new_param["mu"];

  double sd_mu = static_cast<double>(config["sd_mu"]);

  double old_logpost = 0.0, new_logpost = 0.0, p_accept = 0.0;

  // proposal (normal distribution with SD: config$sd_mu)

  new_mu[0] += R::rnorm(0.0, sd_mu); // new proposed value


  // automatic rejection of negative mu

  if (new_mu[0] < 0.0) {
    return param;
  }


  // compute likelihoods
  old_logpost = cpp_ll_genetic(data, param, R_NilValue, custom_ll);
  new_logpost = cpp_ll_genetic(data, new_param, R_NilValue, custom_ll);


  // compute priors

  old_logpost += cpp_prior_mu(param, config, custom_prior);
  new_logpost += cpp_prior_mu(new_param, config, custom_prior);


  // acceptance term

  p_accept = exp(new_logpost - old_logpost);


  // acceptance: the new value is already in mu, so we only act if the move is
  // rejected, in which case we restore the previous ('old') value

  if (p_accept < unif_rand()) { // reject new values
    return param;
  }

  return new_param;
}






// ---------------------------

// movement of the Reporting probability 'pi' is done using a dumb normal
// proposal. This is satisfying for now - we only reject a few non-sensical
// values outside the range [0;1]. The SD of the proposal (implicitely contained
// in rand$pi.rnorm1, but really provided through 'config', seems fine as the
// range of real values will never change much. Probably not much point in using
// auto-tuning here.

// [[Rcpp::export(rng = true)]]
Rcpp::List cpp_move_pi(Rcpp::List param, Rcpp::List data, Rcpp::List config,
		       Rcpp::RObject custom_ll = R_NilValue,
		       Rcpp::RObject custom_prior = R_NilValue) {

  // deep copy here for now, ultimately should be an arg.

  Rcpp::List new_param = clone(param);
  Rcpp::NumericVector pi = param["pi"]; // these are just pointers
  Rcpp::NumericVector new_pi = new_param["pi"]; // these are just pointers

  double sd_pi = static_cast<double>(config["sd_pi"]);

  double old_logpost = 0.0, new_logpost = 0.0, p_accept = 0.0;


  // proposal (normal distribution with SD: config$sd_pi)

  new_pi[0] += R::rnorm(0.0, sd_pi); // new proposed value


  // automatic rejection of pi outside [0;1]

  if (new_pi[0] < 0.0 || new_pi[0] > 1.0) {
    return param;
  }


  // compute likelihoods
  old_logpost = cpp_ll_reporting(data, param, R_NilValue, custom_ll);
  new_logpost = cpp_ll_reporting(data, new_param, R_NilValue, custom_ll);


  // compute priors

  old_logpost += cpp_prior_pi(param, config, custom_prior);
  new_logpost += cpp_prior_pi(new_param, config, custom_prior);


  // acceptance term

  p_accept = exp(new_logpost - old_logpost);


  // acceptance: the new value is already in pi, so we only act if the move is
  // rejected, in which case we restore the previous ('old') value

  if (p_accept < unif_rand()) { // reject new values
    return param;
  }

  return new_param;
}






// ---------------------------

// movement of the Reporting probability 'pi' is done using a dumb normal
// proposal. This is satisfying for now - we only reject a few non-sensical
// values outside the range [0;1]. The SD of the proposal (implicitely contained
// in rand$pi.rnorm1, but really provided through 'config', seems fine as the
// range of real values will never change much. Probably not much point in using
// auto-tuning here.

// [[Rcpp::export(rng = true)]]
Rcpp::List cpp_move_tau(Rcpp::List param, Rcpp::List data, Rcpp::List config,
			Rcpp::RObject custom_ll = R_NilValue,
			Rcpp::RObject custom_prior = R_NilValue) {

  // deep copy here for now, ultimately should be an arg.

  Rcpp::List new_param = clone(param);
  Rcpp::NumericVector tau = param["tau"]; // these are just pointers
  Rcpp::NumericVector new_tau = new_param["tau"]; // these are just pointers

  int max_gamma = static_cast<int>(config["max_kappa"]);
  double eps = param["eps"];
  Rcpp::NumericVector p_ward = data["p_ward"];
  
  double sd_pi = static_cast<double>(config["sd_pi"]);

  double old_logpost = 0.0, new_logpost = 0.0, p_accept = 0.0;


  // proposal (normal distribution with SD: config$sd_pi)

  new_tau[0] += R::rnorm(0.0, sd_pi); // new proposed value


  // automatic rejection of tau outside [0;1]

  if (new_tau[0] < 0.0 || new_tau[0] > 1.0) {
    return param;
  }

  // update the transition probabilities between wards
  new_param["ward_mat"] = get_ward_p(p_ward, eps, new_tau[0], max_gamma);
  
  // compute likelihoods
  old_logpost = cpp_ll_contact(data, param, R_NilValue, custom_ll);
  new_logpost = cpp_ll_contact(data, new_param, R_NilValue, custom_ll);


  // compute priors

  old_logpost += cpp_prior_tau(param, config, custom_prior);
  new_logpost += cpp_prior_tau(new_param, config, custom_prior);


  // acceptance term

  p_accept = exp(new_logpost - old_logpost);

  // acceptance: the new value is already in tau, so we only act if the move is
  // rejected, in which case we restore the previous ('old') value

  if (p_accept < unif_rand()) { // reject new values
    return param;
  } 

  return new_param;
}






// ---------------------------

// movement of the contact reporting coverage 'eps' is done using a dumb normal
// proposal. This is satisfying for now - we only reject a few non-sensical
// values outside the range [0;1]. The SD of the proposal is provided through
// 'config'; this seems fine as the range of real values will never change
// much. Probably not much point in using auto-tuning here.

// [[Rcpp::export(rng = true)]]
Rcpp::List cpp_move_eps(Rcpp::List param, Rcpp::List data, Rcpp::List config,
			Rcpp::RObject custom_ll = R_NilValue,
			Rcpp::RObject custom_prior = R_NilValue) {

  // deep copy here for now, ultimately should be an arg.

  Rcpp::List new_param = clone(param);
  Rcpp::NumericVector eps = param["eps"]; // these are just pointers
  Rcpp::NumericVector new_eps = new_param["eps"]; // these are just pointers
  Rcpp::LogicalVector move_eps = config["move_eps"]; // these are just pointers
  
  int max_gamma = static_cast<int>(config["max_kappa"]);
  double tau = param["tau"];
  Rcpp::NumericVector p_ward = data["p_ward"];
  Rcpp::NumericMatrix ward_matrix = Rcpp::as<Rcpp::NumericMatrix>(data["ward_matrix"]);

  //  double sd_eps = static_cast<double>(config["sd_eps"]);
  Rcpp::NumericVector sd_eps = config["sd_eps"];

  double old_logpost = 0.0, new_logpost = 0.0, p_accept = 0.0;

  for (size_t i = 0; i < eps.size(); i++) {

    if(move_eps[i]) {
      
      // proposal (normal distribution with SD: config$sd_eps)

      new_eps[i] += R::rnorm(0.0, sd_eps[i]); // new proposed value

      // automatic rejection of eps outside [0;1]

      if (new_eps[i] < 0.0 || new_eps[i] > 1.0) {
	return param;
      }

      if(ward_matrix.nrow() > 0) {
	new_param["ward_mat"] = get_ward_p(p_ward, new_eps[i], tau, max_gamma);
	new_param["ward_mat_1"] = get_ward_p(p_ward, new_eps[i], tau, 1);
      }

      // compute likelihoods
      old_logpost = cpp_ll_contact(data, param, R_NilValue, custom_ll);
      new_logpost = cpp_ll_contact(data, new_param, R_NilValue, custom_ll);

      // compute priors
      old_logpost += cpp_prior_eps(param, config, custom_prior);
      new_logpost += cpp_prior_eps(new_param, config, custom_prior);

      // acceptance term
      p_accept = exp(new_logpost - old_logpost);

      // acceptance: the new value is already in eps, so we only act if the move is
      // rejected, in which case we restore the previous ('old') value

      if (p_accept < unif_rand()) { // reject new values
	new_eps[i] = eps[i];
      }

    } else {

      new_eps[i] = eps[i];
	
    }

  }
  
  return new_param;
  
}






// ---------------------------

// movement of the contact sensitivity 'eta' is done using a dumb normal
// proposal. This is satisfying for now - we only reject a few non-sensical
// values outside the range [0;1]. The SD of the proposal is provided through
// 'config'; this seems fine as the range of real values will never change
// much. Probably not much point in using auto-tuning here.

// [[Rcpp::export(rng = true)]]
Rcpp::List cpp_move_eta(Rcpp::List param, Rcpp::List data, Rcpp::List config,
			Rcpp::RObject custom_ll = R_NilValue,
			Rcpp::RObject custom_prior = R_NilValue) {

  // deep copy here for now, ultimately should be an arg.

  Rcpp::List new_param = clone(param);
  Rcpp::NumericVector eta = param["eta"]; // these are just pointers
  Rcpp::NumericVector new_eta = new_param["eta"]; // these are just pointers
  Rcpp::LogicalVector move_eta = config["move_eta"]; // these are just pointers

  //  double sd_eta = static_cast<double>(config["sd_eta"]);
  Rcpp::NumericVector sd_eta = config["sd_eta"];
  
  double old_logpost = 0.0, new_logpost = 0.0, p_accept = 0.0;

  for (size_t i = 0; i < eta.size(); i++) {

    if(move_eta[i]) {
      
      // proposal (normal distribution with SD: config$sd_eta)

      new_eta[i] += R::rnorm(0.0, sd_eta[i]); // new proposed value

      // automatic rejection of eta outside [0;1]

      if (new_eta[i] < 0.0 || new_eta[i] > 1.0) {
	return param;
      }

      // compute likelihoods
      old_logpost = cpp_ll_contact(data, param, R_NilValue, custom_ll);
      new_logpost = cpp_ll_contact(data, new_param, R_NilValue, custom_ll);

      // compute priors
      old_logpost += cpp_prior_eta(param, config, custom_prior);
      new_logpost += cpp_prior_eta(new_param, config, custom_prior);

      // acceptance term
      p_accept = exp(new_logpost - old_logpost);

      // acceptance: the new value is already in eta, so we only act if the move is
      // rejected, in which case we restore the previous ('old') value

      if (p_accept < unif_rand()) { // reject new values
	new_eta[i] = eta[i];
      }

    } else {

      new_eta[i] = eta[i];
	
    }

  }
  
  return new_param;
  
}






// ---------------------------

// movement of the non-infectious contact rate 'eps' is done using a dumb
// normal proposal. This is satisfying for now - we only reject a few
// non-sensical values outside the range [0;1]. The SD of the proposal is
// provided through 'config'; this seems fine as the range of real values will
// never change much. Probably not much point in using auto-tuning here.

// [[Rcpp::export(rng = true)]]
Rcpp::List cpp_move_lambda(Rcpp::List param, Rcpp::List data, Rcpp::List config,
			   Rcpp::RObject custom_ll = R_NilValue,
			   Rcpp::RObject custom_prior = R_NilValue) {

  // deep copy here for now, ultimately should be an arg.

  Rcpp::List new_param = clone(param);
  Rcpp::NumericVector lambda = param["lambda"]; // these are just pointers
  Rcpp::NumericVector new_lambda = new_param["lambda"]; // these are just pointers
  Rcpp::LogicalVector move_lambda = config["move_lambda"]; // these are just pointers

  //  double sd_lambda = static_cast<double>(config["sd_lambda"]);
  Rcpp::NumericVector sd_lambda = config["sd_lambda"];
  
  double old_logpost = 0.0, new_logpost = 0.0, p_accept = 0.0;


  for (size_t i = 0; i < lambda.size(); i++) {

    if(move_lambda[i]) {
      
      // proposal (normal distribution with SD: config$sd_lambda)

      new_lambda[i] += R::rnorm(0.0, sd_lambda[i]); // new proposed value


      // automatic rejection of lambda outside [0;1]

      if (new_lambda[i] < 0.0 || new_lambda[i] > 1.0) {
	return param;
      }


      // compute likelihoods
      old_logpost = cpp_ll_contact(data, param, R_NilValue, custom_ll);
      new_logpost = cpp_ll_contact(data, new_param, R_NilValue, custom_ll);


      // compute priors

      old_logpost += cpp_prior_lambda(param, config, custom_prior);
      new_logpost += cpp_prior_lambda(new_param, config, custom_prior);


      // acceptance term

      p_accept = exp(new_logpost - old_logpost);


      // acceptance: the new value is already in lambda, so we only act if the move is
      // rejected, in which case we restore the previous ('old') value

      if (p_accept < unif_rand()) { // reject new values
	new_lambda[i] = lambda[i];
      }

    } else {

      new_lambda[i] = lambda[i];
      
    }
  
  }

  return new_param;
}






// ---------------------------

// Movement of infection dates are +/- 1 from current states. These movements
// are currently vectorised, i.e. a bunch of dates are proposed all together;
// this may not be sustainable for larger datasets. The non-vectorised option
// will be slower and speed-up with C/C++ will be more substantial then.

// This version differs from the initial R implementation in several points:

// 1. all cases are moved
// 2. cases are moved one by one
// 3. movement for each case is +/- 1 time unit

// Notes

// - when computing the timing log-likelihood, the descendents of each
// case are also affected.

// - we generate a new vector 'new_t_inf', which will replace the
// previous pointer defining param["t_inf"].

// [[Rcpp::export(rng = true)]]
Rcpp::List cpp_move_t_inf(Rcpp::List param, Rcpp::List data,
			  Rcpp::RObject list_custom_ll = R_NilValue) {

  // deep copy here for now, ultimately should be an arg.
  Rcpp::List new_param = clone(param);
  Rcpp::IntegerVector t_inf = param["t_inf"];
  Rcpp::IntegerVector kappa = param["kappa"];
  Rcpp::IntegerVector new_t_inf = new_param["t_inf"];
  Rcpp::NumericMatrix n_contacts = param["n_contacts"];
  Rcpp::NumericMatrix old_n_contacts = clone(n_contacts);
  Rcpp::NumericMatrix new_n_contacts = clone(n_contacts);
  Rcpp::IntegerVector alpha = param["alpha"];
  Rcpp::IntegerVector local_cases;
  Rcpp::NumericMatrix diff_n_contacts;
  Rcpp::List list_functions;
  size_t N = static_cast<size_t>(data["N"]);

  double old_loglike = 0.0, new_loglike = 0.0, p_accept = 0.0;
  double old_loc_loglike = 0.0, new_loc_loglike = 0.0, p_loc_accept = 0.0;

  for (size_t i = 0; i < N; i++) {

    local_cases = cpp_find_descendents(param["alpha"], i+1);

    // loglike with current value
    // old_loglike = cpp_ll_timing(data, param, R_NilValue);
    old_loc_loglike = cpp_ll_timing(data, param, i+1, list_custom_ll); // term for case 'i' with offset

    // term descendents of 'i'
    if (local_cases.size() > 0) {
      old_loc_loglike += cpp_ll_timing(data, param, local_cases, list_custom_ll);
    }

    // the contact likelihood isn't subsetted by i, so only add once
    // we have to pass a single function, not list of functions, to ll_contact
    if(list_custom_ll != R_NilValue) {
      Rcpp::List list_functions = Rcpp::as<Rcpp::List>(list_custom_ll);
      old_loc_loglike += cpp_ll_contact(data, param, i, list_functions["contact"]);
    } else {
      old_loc_loglike += cpp_ll_contact(data, param, i);
    }
    
    // proposal (+/- 1)
    new_t_inf[i] += unif_rand() > 0.5 ? 1 : -1; // new proposed value

    // calculate changes in contact types upon changing infection time
    if(n_contacts.nrow() > 0) {
      diff_n_contacts = t_inf_change(data, alpha, kappa, i+1, t_inf[i], new_t_inf[i]);
      for(size_t j = 0; j < n_contacts.nrow(); j++) {
	for(size_t k = 0; k < 6; k++) {
	  new_n_contacts(j, k) = old_n_contacts(j, k) + diff_n_contacts(j, k);
	}
      }
    }

    new_param["n_contacts"] = clone(new_n_contacts);
    
    // loglike with new value
    // new_loglike = cpp_ll_timing(data, new_param, R_NilValue);
    new_loc_loglike = cpp_ll_timing(data, new_param, i+1, list_custom_ll); // term for case 'i' with offset

    // term descendents of 'i'
    if (local_cases.size() > 0) {
      new_loc_loglike += cpp_ll_timing(data, new_param, local_cases, list_custom_ll);
    }

    // the contact likelihood isn't subsetted by i, so only add once
    // we have to pass a single function, not list of functions, to ll_contact
    if(list_custom_ll != R_NilValue) {
      list_functions = Rcpp::as<Rcpp::List>(list_custom_ll);
      new_loc_loglike += cpp_ll_contact(data, new_param, i, list_functions["contact"]);
    } else {
      new_loc_loglike += cpp_ll_contact(data, new_param, i);
    }

    // acceptance term
    // p_accept = exp(new_loglike - old_loglike);
    p_loc_accept = exp(new_loc_loglike - old_loc_loglike);

    //    std::cout << old_loc_loglike << " | " << new_loc_loglike << std::endl;
    
    // acceptance: the new value is already in t_inf, so we only act if the move
    // is rejected, in which case we restore the previous ('old') value

    if (p_loc_accept < unif_rand()) { // reject new values
      new_t_inf[i] = t_inf[i];
      new_n_contacts = clone(old_n_contacts);
      new_param["n_contacts"] = clone(old_n_contacts);
      //      std::cout << "reject" << std::endl;
    } else {
      t_inf[i] = new_t_inf[i];
      old_n_contacts = clone(new_n_contacts);
      param["n_contacts"] = clone(new_n_contacts);
      //      std::cout << "accept" << std::endl;
    }
  }

  new_param["n_contacts"] = clone(new_n_contacts);
  return new_param;

}






// ---------------------------

// Movement of ancestries ('alpha') is not vectorised, movements are made one
// case at a time. This procedure is simply about picking an infector at random
// amongst cases preceeding the case considered. The original version in
// 'outbreaker' used to move simultaneously 'alpha', 'kappa' and 't_inf', but
// current implementation is simpler and seems to mix at least as well. Proper
// movement of 'alpha' needs this procedure as well as a swapping procedure
// (swaps are not possible through move.alpha only); in all instances, 'alpha'
// is on the scale 1:N.

// [[Rcpp::export(rng = true)]]
Rcpp::List cpp_move_alpha(Rcpp::List param, Rcpp::List data,
			  Rcpp::RObject list_custom_ll = R_NilValue) {
  Rcpp::List new_param = clone(param);
  Rcpp::IntegerVector alpha = param["alpha"]; // pointer to param$alpha
  Rcpp::IntegerVector t_inf = param["t_inf"]; // pointer to param$t_inf
  Rcpp::IntegerVector kappa = param["kappa"]; // pointer to param$kappa
  Rcpp::IntegerVector new_alpha = new_param["alpha"];
  Rcpp::NumericMatrix n_contacts = param["n_contacts"];
  Rcpp::NumericMatrix new_n_contacts = new_param["n_contacts"];
  
  // Deep copy to prevent pointer exchange
  // Maybe would be quicker to apply alpha_change instead of constant cloning
  Rcpp::NumericMatrix diff_n_contacts;
  
  size_t N = static_cast<size_t>(data["N"]);

  double old_loglike = 0.0, new_loglike = 0.0, p_accept = 0.0;

  //  std::cout << cpp_ll_all(data, new_param, R_NilValue, list_custom_ll)  << std::endl;
  
  for (size_t i = 0; i < N; i++) {

    // only non-NA ancestries are moved, if there is at least 1 choice
    if (alpha[i] != NA_INTEGER && sum(t_inf < t_inf[i]) > 0) {

      // loglike with current value
      // old_loglike = cpp_ll_all(data, param, R_NilValue);
      old_loglike = cpp_ll_all(data, param, i+1, list_custom_ll); // offset

      // proposal (+/- 1)
      new_alpha[i] = cpp_pick_possible_ancestor(t_inf, i+1); // new proposed value (on scale 1 ... N)

      // calculate changes in contact types upon changing infector
      if(n_contacts.nrow() > 0) {
	diff_n_contacts = alpha_change(data, i + 1, kappa[i], t_inf[i], alpha[i], new_alpha[i]);
	for(size_t j = 0; j < n_contacts.nrow(); j++) {
	  for(size_t k = 0; k < 6; k++) {
	    new_n_contacts(j, k) = n_contacts(j, k) + diff_n_contacts(j, k);
	  }
	}
      }

      new_param["n_contacts"] = clone(new_n_contacts);
      
      // loglike with current value
      new_loglike = cpp_ll_all(data, new_param, i+1, list_custom_ll);

      // acceptance term
      p_accept = exp(new_loglike - old_loglike);

      // which case we restore the previous ('old') value
      if (p_accept < unif_rand()) { // reject new values
	new_alpha[i] = alpha[i];
	new_param["n_contacts"] = clone(n_contacts);
	new_n_contacts = clone(n_contacts);
      } else {
	param["n_contacts"] = clone(new_n_contacts);
	alpha[i] = new_alpha[i];
	n_contacts = clone(new_n_contacts);
      }
    }
  }
  
  return new_param;
  
}





// ---------------------------

// Movement of ancestries ('alpha') is not vectorised, movements are made one
// case at a time. This procedure is simply about picking an infector at random
// amongst cases preceeding the case considered. The original version in
// 'outbreaker' used to move simultaneously 'alpha', 'kappa' and 't_inf', but
// current implementation is simpler and seems to mix at least as well. Proper
// movement of 'alpha' needs this procedure as well as a swapping procedure
// (swaps are not possible through move.alpha only); in all instances, 'alpha'
// is on the scale 1:N.

// [[Rcpp::export(rng = true)]]
Rcpp::List cpp_move_model(Rcpp::List param, Rcpp::List data, Rcpp::List config,
			  Rcpp::RObject list_custom_ll = R_NilValue) {
  Rcpp::List new_param = clone(param);

  // Only propose a model swap prop_model_move % of the time
  double prop_model_move = config["prop_model_move"];
  if(prop_model_move < unif_rand()) {
    return(new_param);
  }
  
  Rcpp::IntegerVector alpha = param["alpha"]; // pointer to param$alpha
  Rcpp::IntegerVector t_inf = param["t_inf"]; // pointer to param$t_inf
  Rcpp::IntegerVector t_onw = param["t_onw"]; // pointer to param$t_onw
  Rcpp::IntegerVector kappa = param["kappa"]; // pointer to param$kappa
  Rcpp::IntegerVector ward = param["ward"]; // pointer to param$kappa

  size_t N = static_cast<size_t>(data["N"]);
  size_t K = static_cast<size_t>(config["max_kappa"]);
  
  Rcpp::IntegerVector new_t_onw = new_param["t_onw"];
  Rcpp::IntegerVector new_kappa = new_param["kappa"];
  Rcpp::IntegerVector new_alpha = new_param["alpha"];
  Rcpp::IntegerVector new_ward = new_param["ward"];

  double ncol = data["ward_ncol"];
  int N_ward = data["N_ward"];
  
  double old_loglike = 0.0, new_loglike = 0.0, p_accept = 0.0;
    
  for (size_t i = 0; i < N; i++) {

    // only non-NA ancestries are moved, if there is at least 1 choice
    if (alpha[i] != NA_INTEGER && sum(t_inf < t_inf[i]) > 0) {

      // loglike with current value
      // old_loglike = cpp_ll_all(data, param, R_NilValue);
      old_loglike = cpp_ll_all(data, param, i+1, list_custom_ll); // offset

      // Pick a new sampled ancestor
      new_alpha[i] = cpp_pick_possible_ancestor(t_inf, i+1); // new proposed value (on scale 1 ... N)

      // Make sure m1 <-> m2 are proposed 50% of the time
      bool m1_to_m2 = (unif_rand() > 0.5) ? true : false;
      
      // jump from Model1 to Model2
      if(kappa[i] == 1 && m1_to_m2) {

	// In Model2 (kappa > 1), infection times must be at least two days apart
	if(t_inf[i] - t_inf[new_alpha[i] - 1] < 2) {
	  p_accept = 0.0;
	} else {

	  // Propose t_onw between the infection times of the two sampled cases
	  Rcpp::IntegerVector t_seq = Rcpp::seq(t_inf[new_alpha[i] - 1] + 1, t_inf[i] - 1);
	  new_t_onw[i] = Rcpp::as<int>(Rcpp::sample(t_seq, 1, true));

	  // Propose kappa between the 2 and max_kappa
	  Rcpp::IntegerVector kappa_seq = Rcpp::seq(2, K);
	  new_kappa[i] = Rcpp::as<int>(Rcpp::sample(kappa_seq, 1, true));

	  // Propose a new ward
	  Rcpp::IntegerVector ward_seq = Rcpp::seq(1, N_ward);
	  new_ward[i] = Rcpp::as<int>(Rcpp::sample(ward_seq, 1, true));
	  
	  // loglike with current value
	  new_loglike = cpp_ll_all(data, new_param, i+1, list_custom_ll);

	  // The ratio of proposal distributions is length((t_inf_i + 1):(t_inf(alpha_i)-1))
	  if(new_loglike == R_NegInf) {
	    p_accept = 0.0;
	  } else {
	    p_accept = exp(new_loglike - old_loglike)*t_seq.size()/ncol;
	  }
	}
	// Jump from Model2 to Model1
      } else if(kappa[i] > 1 && !m1_to_m2) {

	// No longer inferring t_onw - set to -1000 (for some reason NA_INTEGER causes segfault)
	new_t_onw[i] = -1000;
	new_kappa[i] = 1;
	new_ward[i] = 0;

	// loglike with current value
	new_loglike = cpp_ll_all(data, new_param, i+1, list_custom_ll);
	
	// if the new likelihood state is impossible, reject it - otherwise it
	// might get accepted because the ratio of proposal distribution says so 
	if(new_loglike == R_NegInf) {
	  p_accept = 0.0;
	} else {
	// The proposal probability of going from state M2 -> M1 is 1 (there is
	// only possible Model 1, in that we remove t_onw and kappa_beta_i); the
	// proposal probability of going from M1 -> M2 is 1/(range of possible
	// t_onw * range of possible kappa_beta_i) - we therefore need adjust
	// for this - however the prior on these infection times is 1/(Tend -
	// T0), which goes on the numerator
	  double tdiff = t_inf[i] - t_inf[new_alpha[i] - 1] - 1;
	  if(tdiff == 0.0) {
	    p_accept = 0;
	  } else {
	    p_accept = exp(new_loglike - old_loglike)*ncol/tdiff;
	  }
	}
      } else {
	p_accept = 0.0;
      }

      // which case we restore the previous ('old') value
      if (p_accept < unif_rand()) { // reject new values
	new_alpha[i] = alpha[i];
	new_t_onw[i] = t_onw[i];
	new_kappa[i] = kappa[i];
	new_ward[i] = ward[i];
      } else {
	//	std::cout << "accept" << std::endl;
      }
    }
  }

  return new_param;

}







// ---------------------------

// Movement of ancestries ('alpha') is not vectorised, movements are made one
// case at a time. This procedure is simply about picking an infector at random
// amongst cases preceeding the case considered. The original version in
// 'outbreaker' used to move simultaneously 'alpha', 'kappa' and 't_inf', but
// current implementation is simpler and seems to mix at least as well. Proper
// movement of 'alpha' needs this procedure as well as a swapping procedure
// (swaps are not possible through move.alpha only); in all instances, 'alpha'
// is on the scale 1:N.

// [[Rcpp::export(rng = true)]]
Rcpp::List cpp_move_joint(Rcpp::List param, Rcpp::List data, Rcpp::List config,
			  Rcpp::RObject list_custom_ll = R_NilValue) {
  Rcpp::List new_param = clone(param);

  Rcpp::IntegerVector alpha = param["alpha"]; // pointer to param$alpha
  Rcpp::IntegerVector t_inf = param["t_inf"]; // pointer to param$t_inf
  Rcpp::IntegerVector t_onw = param["t_onw"];
  Rcpp::IntegerVector kappa = param["kappa"];
  Rcpp::IntegerVector ward = param["ward"];

  size_t N = static_cast<size_t>(data["N"]);
  size_t K = static_cast<size_t>(config["max_kappa"]);

  int N_ward = data["N_ward"];
  
  int jump;

  double sd_t_onw = static_cast<double>(config["sd_t_onw"]);

  Rcpp::IntegerVector new_t_onw = new_param["t_onw"];
  Rcpp::IntegerVector new_kappa = new_param["kappa"];
  Rcpp::IntegerVector new_alpha = new_param["alpha"];
  Rcpp::IntegerVector new_ward = new_param["ward"];
  
  double old_loglike = 0.0, new_loglike = 0.0, p_accept = 0.0;
    
  for (size_t i = 0; i < N; i++) {

    // only non-NA ancestries are moved, if there is at least 1 choice
    if (alpha[i] != NA_INTEGER && sum(t_inf < t_inf[i]) > 0) {

      // loglike with current value
      // old_loglike = cpp_ll_all(data, param, R_NilValue);
      old_loglike = cpp_ll_all(data, param, i+1, list_custom_ll); // offset

      // Within Model 1 move; just an alpha move (t_onw and kappa are not inferred)
      if(kappa[i] == 1) {

	// proposal (+/- 1)
	new_alpha[i] = cpp_pick_possible_ancestor(t_inf, i+1); // new proposed value (on scale 1 ... N)
	
	// loglike with current value
	new_loglike = cpp_ll_all(data, new_param, i+1, list_custom_ll);

	// acceptance term
	p_accept = exp(new_loglike - old_loglike);

	// Within Model 2 move; move alpha, kappa and t_onw
      } else if(kappa[i] > 1) {

	jump = (unif_rand() > 0.5) ? 1 : -1;
	new_kappa[i] = new_kappa[i] + jump;
	
	// A move back to kappa = 1 is impossible; that would represent a move to Model 1
	if (new_kappa[i] < 2 || new_kappa[i] > K) {
	  p_accept = 0.0;
	} else {
	  // Only propose a new ancestry prop_alpha_move % of the time - this
	  // allows t_onw to mix independently of alpha at times
	  double prop_alpha_move = config["prop_alpha_move"];
	  if(prop_alpha_move < unif_rand()) {
	    // proposal (+/- 1)
	    new_alpha[i] = cpp_pick_possible_ancestor(t_inf, i+1); // new proposed value (on scale 1 ... N)
	  }
	  // Normal MH move from previous state
	  new_t_onw[i] += std::round(R::rnorm(0.0, sd_t_onw));
	  // Swap ward half the time
	  if(unif_rand() > 0.5) {
	    Rcpp::IntegerVector ward_seq = Rcpp::seq(1, N_ward);
	    new_ward[i] = Rcpp::as<int>(Rcpp::sample(ward_seq, 1, true));
	  }
	  // loglike with current value
	  new_loglike = cpp_ll_all(data, new_param, i+1, list_custom_ll);

	  // acceptance term
	  p_accept = exp(new_loglike - old_loglike);
	}
      }

      // Reject new values if not accepted
      if (p_accept < unif_rand()) { // reject new values
	new_alpha[i] = alpha[i];
	new_t_onw[i] = t_onw[i];
	new_kappa[i] = kappa[i];
	new_ward[i] = ward[i];
      }
    }
  }

  return new_param;

}







// ---------------------------

// The basic movement of ancestries (picking an ancestor at random amongst in
// previous cases) makes swaps of ancestries (A->B) to (B->A) very
// difficult. This function addresses the issue. It is computer-intensive, but
// likely a determining factor for faster mixing. Unlike previous versions in
// the original 'outbreaker' package, all cases are 'moved' here. A swap is
// defined as:

// x -> a -> b  becomes a -> x -> b

// Obviously cases are moved one at a time. We need to use local likelihood
// changes for this move to scale well with outbreak size. The complicated bit
// is that the move impacts all descendents from 'a' as well as 'x'.

// [[Rcpp::export(rng = true)]]
Rcpp::List cpp_move_swap_cases(Rcpp::List param, Rcpp::List data,
			       Rcpp::RObject list_custom_ll = R_NilValue) {

  Rcpp::List new_param = clone(param);
  Rcpp::IntegerVector alpha = param["alpha"]; // pointer to param$alpha
  Rcpp::IntegerVector t_inf = param["t_inf"]; // pointer to param$t_inf
  Rcpp::IntegerVector kappa = param["kappa"]; // pointer to param$kappa
  Rcpp::IntegerVector t_onw = param["t_onw"]; // pointer to param$t_inf
  Rcpp::IntegerVector ward = param["ward"]; // pointer to param$ward
  Rcpp::NumericMatrix n_contacts = param["n_contacts"];
  Rcpp::NumericMatrix new_n_contacts = clone(n_contacts);
  Rcpp::NumericMatrix old_n_contacts = clone(n_contacts);
  Rcpp::NumericMatrix diff_n_contacts;

  Rcpp::List swapinfo; // contains alpha and t_inf
  Rcpp::IntegerVector local_cases;

  bool swap_ward = Rcpp::as<bool>(data["swap_ward"]);
  size_t N = static_cast<size_t>(data["N"]);

  double old_loglike = 0.0, new_loglike = 0.0, p_accept = 0.0;

  Rcpp::List ctd_timed_matrix_list = Rcpp::as<Rcpp::List>(data["ctd_timed_matrix"]);
  size_t n_mat = ctd_timed_matrix_list.size();
  
  for (size_t i = 0; i < N; i++) {

    // only non-NA ancestries are moved, if there is at least 1 choice
    if (alpha[i] != NA_INTEGER && sum(t_inf < t_inf[i]) > 0) {

      // The local likelihood is defined as the likelihood computed for the
      // cases affected by the swap; these include:

      // - 'i'
      // - the descendents of 'i'
      // - 'alpha[i]'
      // - the descendents of 'alpha[i]' (other than 'i')

      local_cases = cpp_find_local_cases(param["alpha"], i+1);


      // loglike with current parameters

      old_loglike = cpp_ll_all(data, param, local_cases, list_custom_ll); // offset

      // proposal: swap case 'i' and its ancestor

      swapinfo = cpp_swap_cases(param, i+1, swap_ward);
      new_param["alpha"] = swapinfo["alpha"];
      new_param["t_inf"] = swapinfo["t_inf"];
      new_param["t_onw"] = swapinfo["t_onw"];
      new_param["kappa"] = swapinfo["kappa"];
      new_param["ward"] = swapinfo["ward"];

      // calculate changes in contact types upon swapping cases
      if(n_contacts.nrow() > 0) {
	//	old_local_n_contacts = local_n_contacts(data, param, local_cases);
	//	new_local_n_contacts = local_n_contacts(data, new_param, local_cases);

	diff_n_contacts = swap_cases_change(data,
					    param,
					    new_param,
					    i + 1,
					    alpha,
					    t_inf,
					    kappa,
					    local_cases,
					    n_mat);

	for(size_t j = 0; j < n_contacts.nrow(); j++) {
	  for(size_t k = 0; k < 6; k++) {
	    new_n_contacts(j, k) = old_n_contacts(j, k) +
	      diff_n_contacts(j, k);
	  }
	}
      }

      new_param["n_contacts"] = clone(new_n_contacts);
      
      // loglike with new parameters

      new_loglike = cpp_ll_all(data, new_param, local_cases, list_custom_ll);


      // acceptance term

      p_accept = exp(new_loglike - old_loglike);
	          
      // acceptance: change param only if new values is accepted

      if (p_accept >= unif_rand()) { // accept new parameters
	param["alpha"] = new_param["alpha"];
	param["t_inf"] = new_param["t_inf"];
	param["t_onw"] = new_param["t_onw"];
	param["kappa"] = new_param["kappa"];
	param["ward"] = new_param["ward"];
	old_n_contacts = clone(new_n_contacts);
	param["n_contacts"] = clone(new_n_contacts);
      }
    }
  }

  param["n_contacts"] = clone(old_n_contacts);
  
  return param;
}






// ---------------------------


// Movement of the number of generations on transmission chains ('kappa') is
// done for one ancestry at a time. As for infection times ('t_inf') we use a
// dumb, symmetric +/- 1 proposal. But because values are typically in a short
// range (e.g. [1-3]) we probably propose more dumb values here. We may
// eventually want to bounce back or use and correct for assymetric proposals.

// [[Rcpp::export(rng = true)]]
Rcpp::List cpp_move_kappa(Rcpp::List param, Rcpp::List data, Rcpp::List config,
			  Rcpp::RObject list_custom_ll = R_NilValue) {
  Rcpp::List new_param = clone(param);
  Rcpp::IntegerVector alpha = param["alpha"]; // pointer to param$alpha
  Rcpp::IntegerVector kappa = param["kappa"]; // pointer to param$kappa
  Rcpp::IntegerVector t_inf = param["t_inf"]; // pointer to param$t_inf
  Rcpp::IntegerVector new_kappa = new_param["kappa"];
  Rcpp::NumericMatrix n_contacts = param["n_contacts"];
  Rcpp::NumericMatrix new_n_contacts = new_param["n_contacts"];
  //  Rcpp::NumericMatrix new_n_contacts = clone(n_contacts);
  Rcpp::NumericMatrix diff_n_contacts;

  
  size_t N = static_cast<size_t>(data["N"]);
  size_t K = static_cast<size_t>(config["max_kappa"]);
  size_t jump;

  double old_loglike = 0.0, new_loglike = 0.0, p_accept = 0.0;

  for (size_t i = 0; i < N; i++) {

    // only non-NA ancestries are moved
    if (alpha[i] != NA_INTEGER) {

      // propose new kappa
      jump = (unif_rand() > 0.5) ? 1 : -1;
      new_kappa[i] = kappa[i] + jump;

      // only look into this move if new kappa is positive and smaller than the
      // maximum value; if not, remember to reset the value of new_kappa to that
      // of kappa, otherwise we implicitely accept stupid moves automatically

      if (new_kappa[i] < 1 || new_kappa[i] > K) {
	new_kappa[i] = kappa[i];
      } else {

	// loglike with current parameters
	old_loglike = cpp_ll_all(data, param, i+1, list_custom_ll);

	// calculate changes in contact types upon changing infection time
	if(n_contacts.nrow() > 0) {
	  diff_n_contacts = kappa_change(data, param, i+1, t_inf[i], alpha[i], kappa[i], new_kappa[i]);
	  for(size_t j = 0; j < n_contacts.nrow(); j++) {
	    for(size_t k = 0; k < 6; k++) {
	      new_n_contacts(j, k) = n_contacts(j, k) + diff_n_contacts(j, k);
	    }
	  }
	}

	new_param["n_contacts"] = clone(new_n_contacts);
	
	// loglike with new parameters
	new_loglike = cpp_ll_all(data, new_param, i+1, list_custom_ll);

	// acceptance term
	p_accept = exp(new_loglike - old_loglike);
	
	// acceptance: change param only if new values is accepted
	if (p_accept >= unif_rand()) { // accept new parameters
	  //	  	  	   Rprintf("\naccepting kappa:%d  (p: %f  old ll:  %f  new ll: %f",
	  //	   	   new_kappa[i], p_accept, old_loglike, new_loglike);
	  kappa[i] = new_kappa[i];
	  n_contacts = clone(new_n_contacts);
	  //param["kappa"] = new_kappa;
	} else {
	  new_kappa[i] = kappa[i];
	  new_n_contacts = clone(n_contacts);
	}
      }
    }

  }

  param["n_contacts"] = clone(n_contacts);
  
  return param;
}
