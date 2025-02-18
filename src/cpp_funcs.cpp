#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double ddirmultinomial_cpp(NumericVector x, double n, NumericVector alpha, bool log) {

  int R = alpha.size();
  double sum_alpha = sum(alpha);
  double temp = 0.0;

  for (int r = 0; r < R; ++r) {
    temp += lgamma(x[r] + alpha[r]) - lgamma(alpha[r]) - lgamma(x[r] + 1.0);
  }

  double log_prob = lgamma(sum_alpha) + lgamma(n + 1.0)- lgamma(n + sum_alpha) + temp;

  if (log) {
    return log_prob;
  } else {
    return exp(log_prob);
  }


}

// [[Rcpp::export]]
NumericMatrix compute_probability_cpp(NumericMatrix Y, NumericVector omega_M,
                                      NumericMatrix q_star_1_J,
                                      NumericVector gamma_1_J_star) {

  // load user-defined R function for density of DirMulti
  //Function ddirmultinomial("ddirmultinomial");

  int C = Y.ncol(); // Number of columns in Y (C)
  int J = q_star_1_J.nrow(); // Number of rows in mu q_star (J)

  // Create a matrix to store the result (C x J)
  NumericMatrix probabilities(C, J);

  // Loop over each column (1:C)
  for (int c = 0; c < C; ++c) {

    // Initialize a vector to store log probabilities (LP) for the current cell
    NumericVector LP(J);

    // Loop over each cluster j (1:J)
    for (int j = 0; j < J; ++j) {

      NumericVector alpha = q_star_1_J(j, _) * gamma_1_J_star[j];
      double n = sum(Y(_, c));
      NumericVector Y_c = Y(_, c);
      //double log_prob = as<double>(ddirmultinomial(Y_c, n, alpha, true));
      double log_like = ddirmultinomial_cpp(Y_c, n, alpha, true);

      LP[j] = log_like + log(omega_M[j]);

    }

    // Compute the normalizing constant (nc) as the negative max(LP)
    double nc = -max(LP);

    // Compute the probabilities (P) and normalize
    NumericVector P = exp(LP + nc) / sum(exp(LP + nc));

    // Store probabilities
    for (int j = 0; j < J; ++j) {
      probabilities(c, j) = P[j];
    }
  }

  return probabilities;
}

// [[Rcpp::export]]
int draw_cat(NumericVector prob){

  // number of categories
  int k = prob.size();

  // the function returns a one-hot encoding
  // set first argument to 1 for a Categorical distribution
  IntegerVector outcome(k);
  rmultinom(1, prob.begin(), k, outcome.begin());

  // return the label
  // we need to target the label '1' (others are all 0)


  // if use match function, the first element should be a vector,
  // and it returns type IntegerVector. No need to increment the index by 1.
  //IntegerVector target(1,1);
  //IntegerVector ix = match(target, Z);

  int ix = which_max(outcome) + 1;

  return ix;

}


// [[Rcpp::export]]
IntegerVector Z_sample_cpp(NumericMatrix Prob){

  int C = Prob.nrow();

  IntegerVector Z(C);

  for (int c = 0; c < C; ++c){
    NumericVector prob_c = Prob.row(c);
    Z[c] = draw_cat(prob_c);
  }

  return Z;

}


// [[Rcpp::export]]
double q_star_logprob_cpp(NumericMatrix Y_sub,
                          NumericVector q_j_star,
                          double gamma_j_star,
                          NumericVector alpha_h){

  int R = Y_sub.nrow();
  int C_j = Y_sub.ncol();

  double log_prob1 = 0.0;

  // from likelihood
  for (int c = 0; c < C_j; ++c) {
    for (int r = 0; r < R; ++r){
      log_prob1 += lgamma(Y_sub(r, c) + q_j_star[r] * gamma_j_star) - lgamma(q_j_star[r] * gamma_j_star) -lgamma(Y_sub(r, c) + 1.0);
    }
  }


  // from prior
  double log_prob2 = sum((alpha_h - 1.0) * log(q_j_star));

  double log_prob = log_prob1 + log_prob2;

  return log_prob;

}


// [[Rcpp::export]]
double gamma_logprob_cpp(NumericMatrix Y_sub,
                         NumericVector q_j_star,
                         double gamma_j_star,
                         double a_gamma,
                         double b_gamma,
                         double lb_gamma){

  int R = Y_sub.nrow();
  int C_j = Y_sub.ncol();

  double log_prob1 = 0.0;

  // from likelihood
  for (int c = 0; c < C_j; ++c) {

    // Total counts N
    double N = 0.0;

    for (int r = 0; r < R; ++r){
      N += Y_sub(r, c);
      log_prob1 += lgamma(Y_sub(r, c) + q_j_star[r] * gamma_j_star) - lgamma(q_j_star[r] * gamma_j_star) -lgamma(Y_sub(r, c) + 1.0);
    }

    log_prob1 += lgamma(gamma_j_star) + lgamma(N + 1.0) - lgamma(N + gamma_j_star);

  }


  // from prior
  double log_prob2 = (a_gamma - 1.0) * log(gamma_j_star) - b_gamma * gamma_j_star;

  double log_prob = log_prob1 + log_prob2;

  if (gamma_j_star < lb_gamma) {
    log_prob = R_NegInf;
  }

  return log_prob;

}


