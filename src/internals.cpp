
#include "internals.h"
#include <RcppEigen.h>


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

//   This function returns a vector of indices of cases which could be infector
//   of 'i' (i.e., their infection dates preceed that of 'i'). Only tricky bit
//   here is keep in mind that 't_inf' is indexed from 0 to N-1, while 'i' and
//   'alpha' (ancestors) are values from 1 to N.

//   Original R code:

// are.possible.alpha <- function(t_inf, i) {
//     if (length(i)>1) {
//         stop("i has a length > 1")
//     }
//     if (any(t_inf[i]==min(t_inf))) {
//         return(NA)
//     }
//     return(which(t_inf < t_inf[i[1]]))
// }

// [[Rcpp::export()]]
std::vector<int> cpp_are_possible_ancestors(Rcpp::IntegerVector t_inf, size_t i) {
  size_t n = t_inf.size();
  std::vector<int> out;
  out.reserve(n);
  for (size_t j = 0; j < n; j++) {
    if (t_inf[j] < t_inf[i-1]) { // offset
      out.push_back(j+1); // +1 needed for range 1 ... N
    }
  }
  return out;
}





// ---------------------------

//  This function samples a single value from a vector of integers.

// [[Rcpp::export()]]
size_t cpp_sample1(std::vector<int> x) {
  if (x.size() < 1) {
    Rcpp::Rcerr << "Trying to sample from empty vector" << std::endl;
    Rcpp::stop("Trying to sample from empty vector");
  }

  return x[unif_rand() * x.size()];
}




// ---------------------------

//    This function choose a possible infector for case 'i'; 'i' is on the scale
//    1:N

//    Original R version:

// .choose.possible.alpha <- function(t_inf, i) {
//     return(sample(are.possible.alpha(t_inf=t_inf, i=i), 1))
// }

// [[Rcpp::export()]]
size_t cpp_pick_possible_ancestor(Rcpp::IntegerVector t_inf, size_t i) {
  return cpp_sample1(cpp_are_possible_ancestors(t_inf, i));
}






// ---------------------------

// This function returns the descendents of a given case 'i' in the current
// ancestries; 'i' is on the scale 1:N. The output is also on the scale 1:N.

// Original R version:

// find.descendents <- function(param, i) {
//   ## find descendents
//     which(param.current$alpha==i)
//  }

// [[Rcpp::export()]]
Rcpp::IntegerVector cpp_find_descendents(Rcpp::IntegerVector alpha, size_t i) {
  size_t counter = 0, n = 0;

  // determine size of output vector and create it
  for (size_t j = 0; j < alpha.size(); j++) {
    if (alpha[j] == i) n++;
  }
  
  Rcpp::IntegerVector out(n);

  // fill in output vector
  for (size_t j = 0; j < alpha.size(); j++) {
    if (alpha[j] == i) {
      out[counter++] = j + 1; // offset
    }
  }
  return out;
}







// ---------------------------

// This function returns a vector of indices of cases which are 'local' to a
// case 'i'. Locality is defined as the following set of cases:

// - 'i'
// - the descendents of 'i'
// - 'alpha[i-1]'
// - the descendents of 'alpha[i]' (excluding 'i')

// where 'alpha' is a IntegerVector storing ancestries. Note that 'i' and
// 'alpha' are on the scale 1:N. 

// [[Rcpp::export()]]
Rcpp::IntegerVector cpp_find_local_cases(Rcpp::IntegerVector alpha, size_t i) {
  // determine descendents of 'i':
  Rcpp::IntegerVector desc_i = cpp_find_descendents(alpha, i);
  size_t n = desc_i.size() + 1; // +1 is to count 'i' itself
  
  // determine descendents of 'alpha[i]':
  Rcpp::IntegerVector desc_alpha_i = cpp_find_descendents(alpha,
							  (size_t) alpha[i-1]);
  if (alpha[i-1] != NA_INTEGER) {
    n += desc_alpha_i.size();
  }

  // create output
  Rcpp::IntegerVector out(n);
  size_t counter = 0;

  // 'i'
  out[counter++] = i;

  // 'descendents of 'i'
  for (size_t j = 0; j < desc_i.size(); j++) {
    out[counter++] = desc_i[j];
  }
  
  if (alpha[i-1] != NA_INTEGER) {
    // alpha[i-1] ...
    out[counter++] = alpha[i-1];
    
    // ... and its descendents
    for (size_t j = 0; j < desc_alpha_i.size(); j++) {
      if ( desc_alpha_i[j] != i) {
	out[counter++] = desc_alpha_i[j];
      }
    }
  }

  return out;
}








// ---------------------------

// This function swaps cases in a transmission tree. The focus case is 'i', and
// is swapped with its ancestor 'x=alpha[i-1]'. In other words the change is
// from: x -> i to i -> x
// Involved changes are:

// - descendents of 'i' become descendents of 'x'
// - descendents of 'x' become descendents of 'i'
// - the infector if 'i' becomes the infector of 'x' (i.e. alpha[x-1])
// - the infector if 'x' becomes 'i'
// - infection time of 'i' becomes that of 'x'
// - infection time of 'x' becomes that of 'i'

// Note on indexing: 'i', 'x', and values of alpha are on the scale 1:N. The
// function's output is a list with new alpha and t_inf.

// Note on forbidden swaps: two types of swaps are excluded:
// - 'i' is imported, so that 'alpha[i-1]' is NA_INTEGER
// - 'x' is imported, so that 'alpha[x-1]' is NA_INTEGER

// [[Rcpp::export()]]
Rcpp::List cpp_swap_cases(Rcpp::List param, size_t i, bool swap_place) {
  
  Rcpp::IntegerVector alpha_in = param["alpha"];
  Rcpp::IntegerVector t_inf_in = param["t_inf"];
  Rcpp::IntegerVector t_onw_in = param["t_onw"];
  Rcpp::IntegerVector kappa_in = param["kappa"];
  
  Rcpp::IntegerVector alpha_out = clone(alpha_in);
  Rcpp::IntegerVector t_inf_out = clone(t_inf_in);
  Rcpp::IntegerVector t_onw_out = clone(t_onw_in);
  Rcpp::IntegerVector kappa_out = clone(kappa_in);
  
  Rcpp::List out;
  out["alpha"] = alpha_out;
  out["t_inf"] = t_inf_out;
  out["t_onw"] = t_onw_out;
  out["kappa"] = kappa_out;

  size_t N = alpha_in.size();
  
  // escape if the case is imported, i.e. alpha[i-1] is NA
  
  if (alpha_in[i-1] == NA_INTEGER) {
    return out;
  }


  // escape if ancestor of the case is imported, i.e. alpha[x-1] is NA
  // if we want to swap imports, we can set swap_place to TRUE
  
  size_t x = (size_t) alpha_in[i-1];
  if (alpha_in[x-1] == NA_INTEGER) {
   return out;
  }
  
 
  // replace ancestries:
  // - descendents of 'i' become descendents of 'x'
  // - descendents of 'x' become descendents of 'i'

  // do this 3/4 of the time; the other 1/4  keep downstream ancestries the same
  // this improves mixing
  
  if(unif_rand() > 0.25) {
    for (size_t j = 0; j < N; j++) {
      if (alpha_in[j] == i) {
	alpha_out[j] = x;
      } else if (alpha_in[j] == x) {
	alpha_out[j] = i;
      }
    }
  }

  // the ancestor of 'i' becomes an ancestor of 'x'

  alpha_out[i-1] = alpha_in[x-1];

  // 'i' is now the ancestor of 'x'
  alpha_out[x-1] = i;
  
  // swap infections times of 'i' and 'x'
  t_inf_out[i-1] =   t_inf_in[x-1];
  t_inf_out[x-1] =   t_inf_in[i-1];

  // swap_place is a variable used to fix the imputed place, gamma and infection
  // time of a case when initiating the tree with eps = 1.0 - otherwise kappa
  // will be changed even when move_kappa = FALSE and the tree will be messed up
  if(swap_place) {

    t_onw_out[i-1] =   t_onw_in[x-1];
    t_onw_out[x-1] =   t_onw_in[i-1];
  
    kappa_out[i-1] =   kappa_in[x-1];
    kappa_out[x-1] =   kappa_in[i-1];
    
  }
  
  return out;
  
}






// ---------------------------

// This function returns the number of mutations between two cases from a 'data'
// object. It uses the indexing of cases in the DNA matrix to ensure
// correspondance between cases and their sequences (not all cases may have a
// sequence). 

// i and j are indices of cases on the scale 1:N; note that the vectors and
// matrices are indexed on 0:(N-1).

// [[Rcpp::export()]]
size_t cpp_get_n_mutations(Rcpp::List data, size_t i, size_t j) {
  Rcpp::LogicalVector has_dna = data["has_dna"];
  Rcpp::IntegerVector id_in_dna = data["id_in_dna"];
  Rcpp::IntegerMatrix D = data["D"];

  
  // Ideally we should return NA_integer here, but then the type of the function
  // should be a Rcpp::IntegerVector, which would complicate things. The second
  // best thing we can do here really is to issue an error when trying to get
  // number of mutations between cases with missing sequences.
  
  if (!(has_dna[i-1] && has_dna[j-1])) {
    Rcpp::stop("Trying to get genetic distances between missing sequences.");
  } 

  size_t out = D(id_in_dna[i-1] - 1, id_in_dna[j-1] - 1);

  return out;
  
}






// ---------------------------

// This function looks up a transmission chain to find the most recent ancestor
// with a sequence, for a given case 'i'. It stops at two conditions: i) it
// finds a sequenced ancestor, or ii) the current ancestor is 'NA'. It returns a
// List with three values: i) the index of the most recent ancestor (on the
// scale 1:N), ii) the total number of generations between this case and 'i',
// and iii) a logical 'found_sequenced_ancestor'. If the latter is FALSE, then
// previous values are 'NA_INTEGER'.

// This is the exported interface. It calls upon a non-exported function
// (lookup_sequenced_ancestor) which does not make memory allocation for the
// output, but instead modifies one of its arguments. This trade-off pays as it
// allows for unit testing via the interface, but remains quite fast as the
// non-exported function can be used internally. "i" is indexed on 1:N.

// [[Rcpp::export()]]
Rcpp::List cpp_lookup_sequenced_ancestor(Rcpp::List data, Rcpp::List param, size_t i) {
 
  Rcpp::IntegerVector alpha = param["alpha"];
  Rcpp::IntegerVector kappa = param["kappa"];
  Rcpp::LogicalVector has_dna = data["has_dna"];
  
  Rcpp::List out;
  Rcpp::IntegerVector out_ances(1);
  Rcpp::IntegerVector out_n_generations(1);
  Rcpp::LogicalVector out_found_sequenced_ancestor(1);
  out["alpha"] = out_ances;
  out["n_generations"] = out_n_generations;
  out["found_sequenced_ancestor"] = out_found_sequenced_ancestor;
  
  size_t ances[1];
  size_t n_generations[1];
  bool found_sequenced_ancestor[1];

  ances[0] = NA_INTEGER;
  n_generations[0] = NA_INTEGER;
  found_sequenced_ancestor[0] = false;

  // This function modifies its last argument
  lookup_sequenced_ancestor(alpha, kappa, has_dna, i, // inputs
			    ances, n_generations, 
			    found_sequenced_ancestor); // outputs


  out_ances[0] = static_cast<int>(ances[0]);
  out_n_generations[0] =  static_cast<int>(n_generations[0]); 
  out_found_sequenced_ancestor[0] = found_sequenced_ancestor[0];

  return out;
}






// ---------------------------

// This function is the internal version of cpp_lookup_sequenced_ancestor. It is
// not meant to be called by users, only by internal procedures, as it modifies
// the content of its last argument rather than creating a new object, which is
// obviously dangerous. Only use it carefully if you handled the creating of its
// last argument 'out'. 'out_' are technically outputs with three components:
// "ances" (IntegerVector of size 1), "n_generations" (same), and
// "found_sequenced_ancestor" (LogicalVector of length 1). "i" is indexed on
// 1:N.

void lookup_sequenced_ancestor(Rcpp::IntegerVector alpha, Rcpp::IntegerVector kappa, 
			       Rcpp::LogicalVector has_dna, size_t i, 
			       size_t *out_alpha, 
			       size_t *out_n_generations, 
			       bool *out_found_sequenced_ancestor
			       ) {

  if (!has_dna[i - 1] || alpha[i - 1] == NA_INTEGER) {
    return;
  }
  
 
  size_t current_case = i; // this one is indexed on 1:N
  size_t n_generations = kappa[current_case - 1];
  bool ances_has_dna = has_dna[alpha[current_case - 1] - 1]; // offset for indexing vectors
    

  // look recursively for ancestor with sequence if needed
	    
  while (!ances_has_dna && (alpha[current_case - 1] != NA_INTEGER)) {
    current_case = alpha[current_case - 1]; // 1 step back up the transmission chain
    ances_has_dna = (alpha[current_case - 1] != NA_INTEGER) && // need to test for NA *first*
      has_dna[alpha[current_case - 1] - 1]; // offset for indexing vectors
    n_generations += kappa[current_case - 1];
  }


  // change outputs as needed

  
  if (ances_has_dna) {
      out_alpha[0] = alpha[current_case - 1];
      out_n_generations[0] = n_generations;
      out_found_sequenced_ancestor[0] = true;
  } else {
    out_alpha[0] = NA_INTEGER;
    out_n_generations[0] = NA_INTEGER;
    out_found_sequenced_ancestor[0] = false;
  }
  
}



// ---------------------------

// This function returns a boolean indicating if an unobserved case moved
// between places or not - this boils down to determining if the place at the time
// of infection (t_inf) is the same as the place of the infector at the time of
// onplace infection (t_onw)

// - 'i'
// - the descendents of 'i'
// - 'alpha[i-1]'
// - the descendents of 'alpha[i]' (excluding 'i')

// where 'alpha' is a IntegerVector storing ancestries. Note that 'i' and
// 'alpha' are on the scale 1:N. 

// [[Rcpp::export()]]
bool is_between_place(Rcpp::NumericMatrix place_matrix, Rcpp::IntegerVector t_inf,
		     Rcpp::IntegerVector t_onw, Rcpp::IntegerVector alpha,
		     int C_ind, size_t j) {
  
  int ind1 = t_inf[j] + C_ind;
  int ind2 = t_onw[j] + C_ind;
  int place1 = place_matrix(j, ind1);
  int place2 = place_matrix(alpha[j] - 1, ind2);
  bool out = (place1 != place2 &&
	      place1 != 0 &&
	      place2 != 0 &&
	      ind1 >= 0 &&
	      ind2 >= 0 &&
	      ind1 < place_matrix.ncol() &&
	      ind2 < place_matrix.ncol());

  return out;
}






// ---------------------------

// This function returns a 3D matrix (well, just a vector) of transition
// probabilities from location k to location l for all 1:max_kappa generations. It is
// indexed mat[k, l, gamma], ie mat[N*N*(gamma - 1) + N*l + k], where N is the
// number of unique locations

// loc_mat is a matrix x[i,j] of naive transition probabilities of moving from i to j

// int_mat is a matrix x[i,j] of transition probabilities of moving from i
// to j, summed over all intermediate places that are not i or j

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export()]]
Rcpp::NumericVector get_transition_mat(Rcpp::NumericMatrix p_trans,
				       Rcpp::NumericVector p_place,
				       Rcpp::NumericVector p_place_adj,
				       double eps,
				       double tau,
				       double prop_place_unobserved,
				       int max_kappa) {

  // number of unique places (add one to account for unknown place)
  int N = p_trans.nrow() + 1;

  // output 3D matrix (ie vector) of transition probabilities
  Rcpp::NumericVector out(N*N*(max_kappa));
  
    // integrate prob into transition matrix
  Rcpp::NumericMatrix Reps_trans = p_trans * (1 - eps);
  Rcpp::NumericMatrix Rtau_trans = p_trans * (1 - tau);
  Reps_trans.fill_diag(eps);
  Rtau_trans.fill_diag(tau);

  // this represents the average transition probability from unknown - unkown
  // (not unobservd to unobserved) (i.e. integrating over 1-eps and eps) - this
  // simplifies to 1/N_u - this is just a simplifcation so that we don't have to
  // calculate the transition values indivdiually for all unobserved places - we
  // assume they have equal prior probabilities and therefore we can calculate
  // the average value analytically - normally this cell would just contain eps
  // (because thats the probability of i -> i) but actually this represents the
  // average transition probability between all unknown places, most of which
  // are not the same place
  if(prop_place_unobserved > 0) {
    Reps_trans(N-2, N-2) = eps + (1 - eps)*(1 - prop_place_unobserved);
    Rtau_trans(N-2, N-2) = tau + (1 - tau)*(1 - prop_place_unobserved);
  }

  // convert to Eigen for matrix multiplication
  // we need to copy not map so we can manipulate it later
  Eigen::MatrixXd eps_trans(Rcpp::as<Eigen::MatrixXd>(Reps_trans));
  Eigen::MatrixXd tau_trans(Rcpp::as<Eigen::MatrixXd>(Rtau_trans));

  // prior probability of being in a given place
  Eigen::VectorXd ep_place(Rcpp::as<Eigen::VectorXd>(p_place));
  Eigen::VectorXd ep_place_adj(Rcpp::as<Eigen::VectorXd>(p_place_adj));

  // add marginal probabilities for place = 0
  Eigen::MatrixXd marg = get_marginal_trans(eps_trans, ep_place, ep_place_adj);

  // the value that we save has to be readjusted by N_place_unobserved so that
  // we get a probability *to a single place* - but for the iterative
  // integrations we need to keep the complete value
  
  // assign to out (convert to vector then fill out)
  // this is the transition probability when kappa = 1
  std::vector<double> vec_1(marg.data(), marg.data() + marg.size());
  std::copy(vec_1.begin(), vec_1.end(), out.begin());
  int ind = vec_1.size();
  
  // transition matrix in one unobserved generation
  const Eigen::MatrixXd epstau_trans = eps_trans * tau_trans;

  // trans is the matrix to be iteratively updated - it starts with the
  // probabilities after a single observed generation (i.e. a single eps move)
  Eigen::MatrixXd gen_n = eps_trans;
  
  for (size_t kappa = 2; kappa <= max_kappa; kappa++) {

    // matrix multiplication
    gen_n = gen_n * epstau_trans;

    // add marginal probabilities for place = 0
    marg = get_marginal_trans(gen_n, ep_place, ep_place_adj);

    // assign to out
    std::vector<double> vec_n(marg.data(), marg.data() + marg.size());
    std::copy(vec_n.begin(), vec_n.end(), out.begin() + ind);
    
    // update vector index
    ind += vec_n.size();

  }

  return out;
  
}










// ---------------------------

// This function will compute the rowmeans and colmeans weighted by the prior
// probability of being in a given place - these represent the transition
// probabilities to and from an 'unknown' place, which represent the integration
// over all possible routes

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export()]]
Eigen::MatrixXd get_marginal_trans(Eigen::MatrixXd p_trans,
				   Eigen::VectorXd p_place,
				   Eigen::VectorXd p_place_adj) {

  Eigen::VectorXd cmean = p_place_adj.transpose() * p_trans;

  Eigen::VectorXd rmean = p_trans * p_place;
  
  // bigger matrix with extra row and column for place = 0
  Eigen::MatrixXd out(p_trans.rows()+1, p_trans.cols()+1);
  out.bottomRightCorner(p_trans.rows(), p_trans.cols()) = p_trans;

  // fill with row mean and col means
  out.block(1, 0, out.rows()-1, 1) = rmean;
  out.block(0, 1, 1, out.rows()-1) = cmean.transpose();

  // fill (0, 0) with global mean (rowmeans weighted by prior)
  //  out(0,0) = (rmean.array() * p_place.array()).sum();
  out(0,0) = rmean.transpose() * p_place_adj;

  return out;
  
}










// ---------------------------
// Find the number of contacts between cases on their day of infection

// [[Rcpp::export()]]
Rcpp::NumericMatrix t_inf_change(Rcpp::List data,
				 Rcpp::IntegerVector alpha,
				 Rcpp::IntegerVector kappa,
				 size_t p,
				 size_t t_inf_1,
				 size_t t_inf_2) {

  Rcpp::List ctd_timed_matrix_list = Rcpp::as<Rcpp::List>(data["ctd_timed_matrix"]);
  Rcpp::NumericMatrix n_contacts(ctd_timed_matrix_list.size(), 6);

  // column to access in contact data
  int ind1;
  int ind2;
  int place1;
  int place2;
  int place3;
  int place4;
  bool place_changed;
  bool match2;
  bool match3;
  bool ind_error1;
  bool ind_error2;
  int t_inf;
  
  int C_ind = static_cast<int>(data["C_ind"]);
  size_t N = static_cast<size_t>(data["N"]);
  int mat_size;
    
  for (size_t i = 0; i < ctd_timed_matrix_list.size(); i++) {
    
    Rcpp::NumericMatrix contacts_timed = Rcpp::as<Rcpp::NumericMatrix>(ctd_timed_matrix_list[i]);

    mat_size = contacts_timed.ncol();
    
    ind1 = t_inf_1 + C_ind;
    ind2 = t_inf_2 + C_ind;

    ind_error1 = ind1 < 0 || ind1 >= mat_size;
    ind_error2 = ind2 < 0 || ind2 >= mat_size;

    // avoid negative indexing
    if(ind_error1) {
      place1 = 0;
    } else {
      place1 = static_cast<int>(contacts_timed(p - 1, ind1));
    }
    
    // avoid negative indexing
    if(ind_error2) {
      place2 = 0;
    } else {
      place2 = static_cast<int>(contacts_timed(p - 1, ind2));
    }

    // First check if the place has even changed
    place_changed = place1 != 0 &&
      place2 != 0 &&
      place1 != place2;

    for(size_t k = 0; k < N; k++) {
      // skip comparison with self
      if(p - 1 != k) {

	if(ind_error1) {
	  match2 = false;
	} else {
	  // Are they in the same place on the day of infection?
	  place3 = static_cast<int>(contacts_timed(k, ind1));
	  match2 = place1 != 0 &&
	    place1 == place3;
	}

	if(ind_error2) {
	  match3 = false;
	} else {
	  // Are they in the same place on the day of infection?
	  place4 = static_cast<int>(contacts_timed(k, ind2));
	  match3 = place2 != 0 &&
	    place2 == place4;
	}

	// If there is a contact in the old proposal and not the new proposal
	if(match2 && !match3) {
	  // when considering the infector - only if kappa = 1
	  if(alpha(p - 1) - 1 == k && kappa(p - 1) == 1) {
	    // we have lost a true pos
	    n_contacts(i, 0) += -1;
	    // we have gained a false neg
	    n_contacts(i, 3) += 1;
	    // when considering a non-infector
	  } else if(alpha(p - 1) != NA_INTEGER) {
	    // we have lost a false pos
	    n_contacts(i, 1) += -1;
	    // we have gained a true neg
	    n_contacts(i, 2) += 1;
	  }
	  // If there is a contact in the new proposal and not the old proposal
	} else if(!match2 && match3) {
	  // when considering the infector - only if kappa = 1
	  if(alpha(p - 1) - 1 == k && kappa(p - 1) == 1) {
	    // we have gained a true pos
	    n_contacts(i, 0) += 1;
	    // we have lost a false neg
	    n_contacts(i, 3) += -1;
	    // when considering a non-infector
	  } else if(alpha(p - 1) != NA_INTEGER) {
	    // we have gained a false pos
	    n_contacts(i, 1) += 1;
	    // we have lost a true neg
	    n_contacts(i, 2) += -1;
	  }
	}
      }
    }
    
  }

  return(n_contacts);
  
}






// ---------------------------
// Find how the types of contacts change with an alpha change

// [[Rcpp::export()]]
Rcpp::NumericMatrix alpha_change(Rcpp::List data,
				 size_t p,
				 size_t kappa,
				 size_t t_inf,
				 size_t alpha_1,
				 size_t alpha_2) {

  Rcpp::List ctd_timed_matrix_list = Rcpp::as<Rcpp::List>(data["ctd_timed_matrix"]);
  Rcpp::NumericMatrix n_contacts(ctd_timed_matrix_list.size(), 6);

  //  size_t true_pos = 0;
  //  size_t false_pos = 0;
  //  size_t true_neg = 0;
  //  size_t false_neg = 0;
  //  size_t imports = 0;
  //  size_t unobsv_case = 0;

  int ind;
  int place1;
  int place2;
  int place3;
  bool match1;
  bool match2;
  
  int C_ind = static_cast<int>(data["C_ind"]);
  size_t N = static_cast<size_t>(data["N"]);
    
  for (size_t i = 0; i < ctd_timed_matrix_list.size(); i++) {
    
    Rcpp::NumericMatrix contacts_timed = Rcpp::as<Rcpp::NumericMatrix>(ctd_timed_matrix_list[i]);

    ind = t_inf + C_ind;

    // avoid negative indexing
    if(ind < 0 || ind >= contacts_timed.ncol()) {
      place1 = 0;
      place2 = 0;
      place3 = 0;
    } else {
      place1 = static_cast<int>(contacts_timed(p - 1, ind));
      place2 = static_cast<int>(contacts_timed(alpha_1 - 1, ind));
      place3 = static_cast<int>(contacts_timed(alpha_2 - 1, ind));
    }

    // First check if the place has even changed
    // Does the tpair have a contact in the old parameter state
    match1 = place1 != 0 &&
      place1 == place2;

    // Does the tpair have a contact in the new parameter state
    match2 = place1 != 0 &&
      place1 == place3; 

    // This currently does not deal with NAs, but imports are skipped in move_alpha anyways
    // in this manner so it doesn't matter

    // If kappa > 1, it basically means all cases are non-transmission pairs so moving alpha
    // doesn't change anything
    if(kappa == 1) {
      if(match1 && !match2) {
	// lost a true pos
	n_contacts(i, 0) += -1;
	// lost a true neg
	n_contacts(i, 2) += -1;
	// gained a false pos
	n_contacts(i, 1) += 1;
	// gained a false neg
	n_contacts(i, 3) += 1;
      } else if(!match1 && match2) {
	// gained a true pos
	n_contacts(i, 0) += 1;
	// lost a palse pos
	n_contacts(i, 1) += -1;
	// gained a true neg
	n_contacts(i, 2) += 1;
	// lost a false neg
	n_contacts(i, 3) += -1;
      }
    }
  }

  return(n_contacts);
  
}









// ---------------------------
// Find the number of contacts between cases on their day of infection

// [[Rcpp::export()]]
Rcpp::NumericMatrix local_n_contacts(Rcpp::List data,
				     Rcpp::List param,
				     Rcpp::IntegerVector p) {

  Rcpp::List ctd_timed_matrix_list = Rcpp::as<Rcpp::List>(data["ctd_timed_matrix"]);
  Rcpp::NumericMatrix n_contacts(ctd_timed_matrix_list.size(), 6);

  Rcpp::IntegerVector alpha = param["alpha"];
  Rcpp::IntegerVector t_inf = param["t_inf"];
  Rcpp::IntegerVector kappa = param["kappa"];

  // column to access in contact data
  int ind;
  int j;
  int place1;
  int place2;
  bool match;
  bool ind_correct;

  int C_ind = static_cast<int>(data["C_ind"]);
  size_t N = static_cast<size_t>(data["N"]);
  int mat_size;
  int alph;
  bool unobs_case;

  for (size_t i = 0; i < ctd_timed_matrix_list.size(); i++) {

    Rcpp::NumericMatrix contacts_timed = Rcpp::as<Rcpp::NumericMatrix>(ctd_timed_matrix_list[i]);

    mat_size = contacts_timed.ncol();
    
    for (size_t l = 0; l < p.size(); l++) {

      // Count only for IDs in p
      j = p[l] - 1;

      // Make sure index is within contact matrix
      ind = t_inf[j] + C_ind;
      ind_correct = ind >= 0 && ind < mat_size;

      alph = alpha[j] - 1;
      unobs_case = kappa[j] > 1;
      
      for(size_t k = 0; k < N; k++) {
      	// skip comparison with self
      	if(j != k) {

      	  if(ind_correct) {
      	    place1 = static_cast<int>(contacts_timed(j, ind));
      	    place2 = static_cast<int>(contacts_timed(k, ind));
      	    match = place1 != 0 && place1 == place2;
      	  } else {
	    match = false;
	  }

      	  // if there is a contact
      	  if(match) {
      	    // when considering the infector
      	    if(alph == k) {
      	      // we have an unobserved case - like an unconnected case
      	      // therefore this is a false positive
      	      if(unobs_case) {
      		n_contacts(i, 1) += 1;
      	      } else {
      		// we have gained a true pos
      		n_contacts(i, 0) += 1;
      	      }
      	    } else {
      	      // we have gained a false pos
      	      n_contacts(i, 1) += 1;
      	    }
      	    // if there is no contact
      	  } else {
      	    // when considering the infector
      	    if(alph == k) {
      	      // we have an unobserved case - like an unconnected case
      	      // therefore this is a true negative
      	      if(unobs_case) {
      		n_contacts(i, 2) += 1;
      	      } else {
      		// we have gained a false negative
      		n_contacts(i, 3) += 1;
      	      }
      	    } else {
      	      // we have gained a true negative
      	      n_contacts(i, 2) += 1;
      	    }
      	  }
      	}
      }
      
    }
  }

  return(n_contacts);
  
}








// ---------------------------
// Find the change in number of contacts upon moving kappa

// [[Rcpp::export()]]
Rcpp::NumericMatrix kappa_change(Rcpp::List data,
				 Rcpp::List param,
				 size_t p,
				 size_t t_inf,
				 size_t alpha,
				 size_t kappa1,
				 size_t kappa2) {

  Rcpp::List ctd_timed_matrix_list = Rcpp::as<Rcpp::List>(data["ctd_timed_matrix"]);
  Rcpp::NumericMatrix n_contacts(ctd_timed_matrix_list.size(), 6);

  //  size_t true_pos = 0;
  //  size_t false_pos = 0;
  //  size_t true_neg = 0;
  //  size_t false_neg = 0;
  //  size_t imports = 0;
  //  size_t unobsv_case = 0;

  int ind;
  int place1;
  int place2;
  bool match;
  
  int C_ind = static_cast<int>(data["C_ind"]);
  size_t N = static_cast<size_t>(data["N"]);
  int mat_size;
    
  for (size_t i = 0; i < ctd_timed_matrix_list.size(); i++) {
    
    Rcpp::NumericMatrix contacts_timed = Rcpp::as<Rcpp::NumericMatrix>(ctd_timed_matrix_list[i]);

    mat_size = contacts_timed.ncol();
    
    ind = t_inf + C_ind;

    match = ind >= 0 && ind < mat_size;

    // Avoid indexing by negative
    // Does the tpair have a contact in the old parameter state
    if(match) {
      place1 = static_cast<int>(contacts_timed(p - 1, ind));
      place2 = static_cast<int>(contacts_timed(alpha - 1, ind));
      match = place1 != 0 && place1 == place2;
    }
    
    // Increasing kappa from 1 -> >1 changes a true positive to a false positive
    // (i.e. they are going from a transmission pair to a non-transmission
    // pair), and it goes from a false negative to a true negative

    // If they are not a transmission pair, it doesn't change anything - kappa
    // only refers to the relationship between tpairs

    // From >1 -> 1 changes from a ntp to a tp. So either a true negative to a
    // false negative, or a false positive to a true positive
    
    // // First check if they are in the same place


    // match = place1 != 0 &&
    //   place1 == place2 &&
    //   ind >= 0 &&
    //   ind < contacts_timed.ncol();

    // a contact with kappa > 1 converts a true positive to a false positive
    if(kappa1 == 1 && kappa2 > 1) {
      // if they previously matched
      if(match) {
	// gained a false positive
	n_contacts(i, 1) += 1;
	// lost a true pos
	n_contacts(i, 0) += -1;
      } else {
	// gained a true negative
	n_contacts(i, 2) += 1;
	// lost a false negative
	n_contacts(i, 3) += -1;
      }
      // Removing an unobserved case
    } else if(kappa1 > 1 && kappa2 == 1) {
      // if they match
      if(match) {
	// gained a true positive
	n_contacts(i, 0) += 1;
	// lost a false pos
	n_contacts(i, 1) += -1;
      } else {
	// lost a true negative
	n_contacts(i, 2) += -1;
	// gained a false negative
	n_contacts(i, 3) += 1;
      }
    }
  }

  return(n_contacts);
  
}






// ---------------------------
// Find the change in number of contacts upon swapping cases
// This involves re-calculating all contacts for i and alpha[i], and the
// infectious contacts for the remaining local contacts

// [[Rcpp::export()]]
Rcpp::NumericMatrix swap_cases_change(Rcpp::List data,
				      Rcpp::List param,
				      Rcpp::List new_param,
				      size_t i,
				      Rcpp::IntegerVector alpha,
				      Rcpp::IntegerVector t_inf,
				      Rcpp::IntegerVector kappa,
				      Rcpp::IntegerVector local_cases,
				      size_t n_mat) {

  Rcpp::IntegerVector new_alpha = new_param["alpha"]; // pointer to param$alpha
  
  Rcpp::NumericMatrix diff_n_contacts(n_mat, 6);
  Rcpp::NumericMatrix tmp_n_contacts;

  Rcpp::IntegerVector i_alpha = Rcpp::IntegerVector::create(i, alpha[i - 1]);

  int id;

  Rcpp::NumericMatrix old_local_n_contacts = local_n_contacts(data, param, local_cases);
  Rcpp::NumericMatrix new_local_n_contacts = local_n_contacts(data, new_param, local_cases);

  for(size_t k = 0; k < n_mat; k++) {
    for(size_t l = 0; l < 6; l++) {
      diff_n_contacts(k, l) = new_local_n_contacts(k, l) - old_local_n_contacts(k, l);
    }
  }

  // If we can get this working it gives a ~20% speed improvement
  // for(size_t j = 0; j < local_cases.size(); j++) {

  //   id = local_cases[j];
    
  //   if(id != i && id != alpha[i - 1]) {

  //     tmp_n_contacts = alpha_change(data,
  // 				    id,
  // 				    kappa[id - 1],
  // 				    t_inf[id - 1],
  // 				    alpha[id - 1],
  // 				    new_alpha[id - 1]);

  //     //      Rcpp::Rcout << tmp_n_contacts << std::endl;

  //     for(size_t k = 0; k < n_mat; k++) {
  // 	for(size_t l = 0; l < 6; l++) {
  // 	  diff_n_contacts(k, l) = diff_n_contacts(k, l) + tmp_n_contacts(k, l);
  // 	}
  //     }
      
  //   }
  // }

  return(diff_n_contacts);

}







// ---------------------------

// This function returns a matrix of ancestries. Each row represents the
// ancestors of one case, ordered from the root to leaf. Indexing by i will only
// update the ancestries of those cases.

// [[Rcpp::export()]]
Rcpp::NumericMatrix cpp_find_ancestors(Rcpp::IntegerVector alpha,
				       Rcpp::NumericMatrix ancestors,
				       SEXP i) {

  if(ancestors.nrow() == 0) return(ancestors);

  Rcpp::NumericMatrix out = clone(ancestors);
  
  Rcpp::IntegerVector vec_i;
  size_t depth, infector, node, leaf;
  Rcpp::IntegerVector::iterator it;
  int start;
  
  size_t N = alpha.size();
  
  if(i == R_NilValue) {
    // evaluate for all cases
    vec_i = Rcpp::seq(1, N);
  } else {
    // only the cases listed in 'i' are retained
    // pass via tmp to so you don't run into declaration issues with vec_i
    Rcpp::IntegerVector tmp(i);
    vec_i = tmp;
  }

  // zero vector for re-setting ancestor row
  Rcpp::IntegerVector zero_vec(out.ncol());
  
  for (size_t j = 0; j < vec_i.size(); j++) {
	
    leaf = vec_i[j];
  
    Rcpp::IntegerVector vec(N);
    depth = 0;
    infector = alpha[leaf-1];
    node = leaf;
    
    while(infector != NA_INTEGER) {
      vec(depth) = infector;
      depth += 1;
      infector = alpha[infector-1];
    }

    // reverse the order so the root is listed first (this makes identification
    // of the MRCA much easier)
    it = std::find(vec.begin(), vec.end(), 0);
    start = std::distance(vec.begin(), it);
    
    out.row(leaf-1) = zero_vec;
    
    out(leaf-1, start) = leaf;
    start -= 1;
    for(int k = start; k >= 0; k--) {
      out(leaf-1, start-k) = vec[k];
    }
    
  }
  
  return out;
  
}





// Returns a vector with the two nodes downstream of the MRCA - these are the
// nodes that we need to calculate divergence times. If one of the leaves is the
// MRCA, return 0 as that node ID. If there is no MRCA, return (0, 0).

// [[Rcpp::export()]]
Rcpp::IntegerVector cpp_find_mrca(size_t i,
				  size_t j,
				  Rcpp::NumericMatrix ancestors) {

  size_t alpha_1, alpha_2;
  
  for(size_t k = 0; k < ancestors.ncol(); k++) {
    alpha_1 = ancestors(i-1,k);
    alpha_2 = ancestors(j-1,k);
    if(alpha_1 != alpha_2) {
      return Rcpp::IntegerVector::create(static_cast<int>(alpha_1),
					 static_cast<int>(alpha_2));
    }
  }

  return Rcpp::IntegerVector::create(0, 0);
  
}





// Returns the output from cpp_find_mrca across all combinations of cases
// provided in combn.

// [[Rcpp::export()]]
Rcpp::NumericMatrix update_mrca(Rcpp::NumericMatrix combn,
				Rcpp::NumericMatrix ancestors) {

  size_t N = combn.nrow();

  if(N == 0) return(combn);
  
  Rcpp::NumericMatrix out(N, 2);

  for(size_t i = 0; i < N; i++) {
    out.row(i) = cpp_find_mrca(combn(i,0), combn(i,1), ancestors);
  }

  return out;
  
}
