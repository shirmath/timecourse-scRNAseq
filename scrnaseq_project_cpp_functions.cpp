#include <RcppArmadillo.h> 

// [[Rcpp::depends(RcppArmadillo)]] 
using namespace Rcpp; 
using namespace arma; 

// [[Rcpp::export]] 
double obj_function(
    const Rcpp::List & data  , // List(Y, X, O, w)
    const Rcpp::List & params,
    const double scale = 1) {
  
  // Conversion from R to arma data structures
  
  const arma::cube & Y = Rcpp::as<arma::cube>(data["Y"]); // responses (m,J,n)
  /*
   const arma::mat & X = Rcpp::as<arma::mat>(data["X"]); // covariates (n,d)
   const arma::mat & O = Rcpp::as<arma::mat>(data["O"]); // offsets (n,J)
   const arma::vec & w = Rcpp::as<arma::vec>(data["w"]); // weights (n) 
   */ 
  const auto mu = Rcpp::as<arma::mat>(params["mu"]);     // true mean (p)
  const auto init_M = Rcpp::as<arma::cube>(params["M"]);     // variational approximation mean (m,J,n)
  const auto init_S = Rcpp::as<arma::cube>(params["S"]);     // variational approximation variances (m,J,n)
  const auto Sigma_Z = Rcpp::as<arma::mat>(params["Sigma_Z"]); // covinv for particular VAR instance (p,p) 
  // const auto Sigma = Rcpp::as<arma::mat>(params["Sigma"]); // covinv of noise in VAR process (p,p)
  const auto A = Rcpp::as<arma::mat>(params["A"]); // VAR(1) coefficient matrix (p,p)

  // set up useful data structures and base quantities
  arma::mat Sigma = Sigma_Z - A * Sigma_Z * trans(A);
  arma::mat mu_rep = repelem(mu, Y.n_rows, 1);
  arma::cube M_and_S_half = init_M + 0.5 * init_S; 
  arma::mat Omega_Z = inv(Sigma_Z);
  arma::mat Omega = inv(Sigma);
  
  arma::mat M1 = init_M.row_as_mat(0); 
  arma::cube M1_M1t_S1 = cube(Y.n_cols, Y.n_cols, Y.n_slices);
  arma::cube Mt_Mt1 = cube(Y.n_rows-1, Y.n_rows-1, Y.n_slices);
  
  for (int i=0; i != Y.n_slices; ++i) {
    M1_M1t_S1.slice(i) = trans(M1.row(i)) * M1.row(i) + diagmat(init_S.slice(i).row(0));
  }
  
  arma::mat M1t_M1_S1_sum = sum(M1_M1t_S1,2);
  
  arma::mat M_temp = mat(init_M.n_cols, init_M.n_rows);
  arma::mat S_temp = mat(init_S.n_cols, init_S.n_cols);
  arma::mat T4_temp = mat(init_S.n_cols, init_S.n_cols);
  
  double term_4 = 0;
  
  for (int i=0; i!= Y.n_slices; ++i) {
    M_temp = trans(init_M.slice(i));
    S_temp = init_S.slice(i);
    for (int t=0; t != Y.n_rows-1; ++t) {
      
      T4_temp = (M_temp.col(t+1) - A*M_temp.col(t)) * trans((M_temp.col(t+1) - A*M_temp.col(t)))  + diagmat(S_temp.row(t+1)) + A * diagmat(S_temp.row(t)) * trans(A);
      term_4 = term_4 - 0.5*trace(T4_temp * Omega);
    }
  }

  // Compute relevant quantities
  double term_1 = accu(Y % (init_M.each_slice() + mu_rep) - exp(M_and_S_half.each_slice() + mu_rep));
  double term_2 = 0.5*(Y.n_cols * Y.n_rows * Y.n_slices - Y.n_slices*log_det(Sigma_Z).real() - Y.n_slices*(Y.n_rows - 1)*log_det(Sigma).real() + accu(log(init_S)));
  double term_3 = -0.5*trace(M1t_M1_S1_sum * Omega_Z);
  
  
  // print statements for debugging
  /*
   cout << "Term 1: " << term_1 << "\n";
   cout << "Term 2: " << term_2 << "\n";
   cout << "Term 3: " << term_3 << "\n";
   cout << "Term 4: " << term_4 << "\n";
   */
  
  return (scale*(term_1 + term_2 + term_3 + term_4));
}

// [[Rcpp::export]] 
arma::mat mu_grad(
    const arma::vec & mu_init, // value of mu at which to evaluate gradient
    const Rcpp::List & data  , // List(Y, X, O, w)
    const Rcpp::List & params,
    const double scale = 1) {
  
  // Conversion from R to arma data structures
  
  const arma::cube & Y = Rcpp::as<arma::cube>(data["Y"]); // responses (m,J,n)
  /*
   const arma::mat & X = Rcpp::as<arma::mat>(data["X"]); // covariates (n,d)
   const arma::mat & O = Rcpp::as<arma::mat>(data["O"]); // offsets (n,J)
   const arma::vec & w = Rcpp::as<arma::vec>(data["w"]); // weights (n) 
   */ 
  
  //arma::mat mu = Rcpp::as<arma::mat>(mu_input);     // true mean (p)
  const auto mu = trans(mu_init);
  const auto init_M = Rcpp::as<arma::cube>(params["M"]);     // variational approximation mean (m,J,n)
  const auto init_S = Rcpp::as<arma::cube>(params["S"]);     // variational approximation variances (m,J,n)

  // set up useful data structures and base quantities
  arma::mat mu_rep = repelem(mu, Y.n_rows, 1);
  arma::cube M_and_S_half = init_M + 0.5 * init_S; 
  
  //compute gradient
  arma::mat grad = scale*sum(sum((Y - exp(M_and_S_half.each_slice() + mu_rep)), 2),0);
  
  return (grad);
}

// [[Rcpp::export]] 
arma::cube M_grad(
    const arma::vec & M_vec, // value of M at which to evaluate gradient
    const Rcpp::List & data  , // List(Y, X, O, w)
    const Rcpp::List & params,
    const double scale = 1) {
  
  // Conversion from R to arma data structures
  
  const arma::cube & Y = Rcpp::as<arma::cube>(data["Y"]); // responses (m,J,n)
  int m = Y.n_rows, J = Y.n_cols, n = Y.n_slices;
  /*
   const arma::mat & X = Rcpp::as<arma::mat>(data["X"]); // covariates (n,d)
   const arma::mat & O = Rcpp::as<arma::mat>(data["O"]); // offsets (n,J)
   const arma::vec & w = Rcpp::as<arma::vec>(data["w"]); // weights (n) 
   */ 
  arma::cube init_M(m, J, n);
  std::copy(M_vec.begin(), M_vec.end(), init_M.begin()); 
  const arma::mat mu = Rcpp::as<arma::mat>(params["mu"]);     // true mean (p)
  const auto init_S = Rcpp::as<arma::cube>(params["S"]);     // variational approximation variances (m,J,n)
  const auto Sigma_Z = Rcpp::as<arma::mat>(params["Sigma_Z"]); // covinv for particular VAR instance (p,p) 
  const auto A = Rcpp::as<arma::mat>(params["A"]); // VAR(1) coefficient matrix (p,p)
  
  // set up useful data structures and base quantities
  arma::mat Sigma = Sigma_Z - A * Sigma_Z * trans(A);
  arma::mat mu_rep = repelem(mu, Y.n_rows, 1);
  arma::cube M_and_S_half = init_M + 0.5 * init_S; 
  arma::mat Omega_Z = inv(Sigma_Z);
  arma::mat Omega = inv(Sigma);
  
  //compute gradient
  arma::cube grad = arma::cube(m, J, n);
  
  for (int t=0; t != m; ++t) {
    //compute gradient value
    mat temp_grad(n, J);
    temp_grad = Y.row_as_mat(t) - exp(M_and_S_half.each_slice() + mu_rep).row_as_mat(t);
    
    if (t == 0) {
      temp_grad = temp_grad - init_M.row_as_mat(0) * Omega_Z;
      temp_grad = temp_grad - (init_M.row_as_mat(0)*trans(A) - init_M.row_as_mat(1))*Omega*A;
    }
    
    
    if (0 < t & t < m-1) {
      temp_grad = temp_grad - (init_M.row_as_mat(t)*trans(A) - init_M.row_as_mat(t+1))*Omega*A;
      temp_grad = temp_grad - (init_M.row_as_mat(t) - init_M.row_as_mat(t-1)*trans(A))*Omega;
    }
    
    
    if (t == m-1) {
      temp_grad = temp_grad - (init_M.row_as_mat(t) - init_M.row_as_mat(t-1)*trans(A))*Omega;
    }
    
    //copy to object to be returned
    for(int j=0; j!= J; ++j) {
      for(int i=0; i != n; ++i) {
        grad(t,j,i) = temp_grad(i,j);
      }
    }
  }
  
  return (scale*grad);
}

// [[Rcpp::export]] 
arma::cube S_grad(
    const arma::vec & S_vec, // value of M at which to evaluate gradient
    const Rcpp::List & data  , // List(Y, X, O, w)
    const Rcpp::List & params,
    const double scale = 1) {
  
  // Conversion from R to arma data structures
  
  const arma::cube & Y = Rcpp::as<arma::cube>(data["Y"]); // responses (m,J,n)
  int m = Y.n_rows, J = Y.n_cols, n = Y.n_slices;
  /*
   const arma::mat & X = Rcpp::as<arma::mat>(data["X"]); // covariates (n,d)
   const arma::mat & O = Rcpp::as<arma::mat>(data["O"]); // offsets (n,J)
   const arma::vec & w = Rcpp::as<arma::vec>(data["w"]); // weights (n) 
   */ 
  arma::cube init_S(m, J, n);
  std::copy(S_vec.begin(), S_vec.end(), init_S.begin()); 
  const arma::mat mu = Rcpp::as<arma::mat>(params["mu"]);     // true mean (p)
  const auto init_M = Rcpp::as<arma::cube>(params["M"]);     // variational approximation variances (m,J,n)
  const auto Sigma_Z = Rcpp::as<arma::mat>(params["Sigma_Z"]); // covinv for particular VAR instance (p,p) 
  const auto A = Rcpp::as<arma::mat>(params["A"]); // VAR(1) coefficient matrix (p,p)
  
  // set up useful data structures and base quantities
  arma::mat Sigma = Sigma_Z - A * Sigma_Z * trans(A);
  arma::mat mu_rep = repelem(mu, Y.n_rows, 1);
  arma::cube M_and_S_half = init_M + 0.5 * init_S; 
  arma::mat Omega_Z = inv(Sigma_Z);
  arma::mat Omega = inv(Sigma);
  arma::mat At_Omega_A = trans(A) * Omega * A;
  //compute gradient
  arma::cube grad = arma::cube(m, J, n);
  
  for (int t=0; t != m; ++t) {
    //compute gradient value
    mat temp_grad(n, J);
    temp_grad = - 0.5*exp(M_and_S_half.each_slice() + mu_rep).row_as_mat(t) + 0.5*pow(init_S.row_as_mat(t), -1);
    
    if (t == 0) {
      temp_grad = temp_grad - 0.5*mat(n,1,fill::ones)*trans(Omega_Z.diag() + At_Omega_A.diag());
    }
    
    
    if (0 < t & t < m-1) {
      temp_grad = temp_grad - 0.5*mat(n,1,fill::ones)*trans(Omega.diag() + At_Omega_A.diag());
    }
    
    
    if (t == m-1) {
      temp_grad = temp_grad - 0.5*mat(n,1,fill::ones)*trans(Omega.diag());
    }
    
    //copy to object to be returned
    for(int j=0; j!= J; ++j) {
      for(int i=0; i != n; ++i) {
        grad(t,j,i) = temp_grad(i,j);
      }
    }
  }
  
  return (scale*grad);
}

// [[Rcpp::export]] 
arma::vec A_grad(
    const arma::vec & A_vec, // value of A at which to evaluate gradient
    const Rcpp::List & data  , // List(Y, X, O, w)
    const Rcpp::List & params,
    const bool & elbo_2 = false,
    const double scale = 1) {
  
  // Conversion from R to arma data structures
  
  const arma::cube & Y = Rcpp::as<arma::cube>(data["Y"]); // responses (m,J,n)
  int m = Y.n_rows, J = Y.n_cols, n = Y.n_slices;
  /*
   const arma::mat & X = Rcpp::as<arma::mat>(data["X"]); // covariates (n,d)
   const arma::mat & O = Rcpp::as<arma::mat>(data["O"]); // offsets (n,J)
   const arma::vec & w = Rcpp::as<arma::vec>(data["w"]); // weights (n) 
   */ 
  arma::mat A(J, J);
  std::copy(A_vec.begin(), A_vec.end(), A.begin()); 
  const arma::mat mu = Rcpp::as<arma::mat>(params["mu"]);     // true mean (p)
  const auto init_M = Rcpp::as<arma::cube>(params["M"]);     // variational approximation means (m,J,n)
  const auto init_S = Rcpp::as<arma::cube>(params["S"]);  // variational approximation variances (m,J,n)
  arma::mat grad = arma::mat(J, J);
  
  if (elbo_2) {
    const auto Sigma = Rcpp::as<arma::mat>(params["Sigma"]); // covinv for VAR process (p,p)
    arma::mat t1 = arma::mat(J,J, fill::zeros);
    
    for (int t=0; t < m-1; t++) {
      t1 = t1 + trans(init_M.row_as_mat(t+1)) * init_M.row_as_mat(t);
      t1 = t1 - A * (trans(init_M.row_as_mat(t)) * init_M.row_as_mat(t) + diagmat(sum(init_S.row_as_mat(t), 0)));
    }
    
    grad = inv(Sigma) * t1;
    
  } else {
    const auto Sigma_Z = Rcpp::as<arma::mat>(params["Sigma_Z"]); // covinv for particular VAR instance (p,p) 
    
    // set up useful data structures and base quantities
    arma::mat Sigma = Sigma_Z - A * Sigma_Z * trans(A);
    arma::mat Omega_Z = inv(Sigma_Z);
    arma::mat Omega = inv(Sigma);
    arma::mat Omega_A_Sigma_Z = Omega * A * Sigma_Z;
    
    arma::mat t2 = arma::mat(J, J, fill::zeros);
    arma::mat quad_term = arma::mat(J, J, fill::zeros);
    
    for (int t = 0; t < m-1; ++t) {
      quad_term = init_M.row_as_mat(t+1) - init_M.row_as_mat(t) * trans(A);
      t2 = t2 + trans(quad_term) * (quad_term * Omega_A_Sigma_Z - init_M.row_as_mat(t));
      t2 = t2 + diagmat(sum(init_S.row_as_mat(t+1), 0)) * Omega_A_Sigma_Z + A * diagmat(sum(init_S.row_as_mat(t), 0)) * (eye(J, J) + trans(A) * Omega_A_Sigma_Z);
    }
    
    //compute gradient
    grad = n * (m-1) * Omega_A_Sigma_Z - Omega * t2;
  }
  
  return (scale*vectorise(grad));
}

// [[Rcpp::export]] 
arma::mat Sigma_Z_grad(
    const arma::vec & Sigma_Z_vec, // value of M at which to evaluate gradient
    const Rcpp::List & data  , // List(Y, X, O, w)
    const Rcpp::List & params,
    const double scale = 1) {
  
  // Conversion from R to arma data structures
  
  const arma::cube & Y = Rcpp::as<arma::cube>(data["Y"]); // responses (m,J,n)
  int m = Y.n_rows, J = Y.n_cols, n = Y.n_slices;
  /*
   const arma::mat & X = Rcpp::as<arma::mat>(data["X"]); // covariates (n,d)
   const arma::mat & O = Rcpp::as<arma::mat>(data["O"]); // offsets (n,J)
   const arma::vec & w = Rcpp::as<arma::vec>(data["w"]); // weights (n) 
   */ 
  arma::mat Sigma_Z(J, J);
  std::copy(Sigma_Z_vec.begin(), Sigma_Z_vec.end(), Sigma_Z.begin()); 
  const arma::mat mu = Rcpp::as<arma::mat>(params["mu"]);     // true mean (p)
  const auto init_M = Rcpp::as<arma::cube>(params["M"]);     // variational approximation variances (m,J,n)
  const auto init_S = Rcpp::as<arma::cube>(params["S"]); // covinv for particular VAR instance (p,p) 
  const auto A = Rcpp::as<arma::mat>(params["A"]); // VAR(1) coefficient matrix (p,p)
  
  // set up useful data structures and base quantities
  arma::mat Sigma = Sigma_Z - A * Sigma_Z * trans(A);
  arma::mat Omega_Z = inv(Sigma_Z);
  arma::mat Omega = inv(Sigma);
  arma::mat At_Omega_A = trans(A) * Omega * A;
  
  arma::mat S1_sum = diagmat(sum(init_S.row_as_mat(0), 0));
  arma::mat B = arma::mat(J, J, fill::zeros);
  arma::mat quad_term = arma::mat(J, J, fill::zeros);
  
  
  
  for (int t = 0; t < m-1; ++t) {
    quad_term = init_M.row_as_mat(t+1) - init_M.row_as_mat(t) * trans(A);
    B = B + trans(quad_term) * quad_term + diagmat(sum(init_S.row_as_mat(t+1), 0)) + A * diagmat(sum(init_S.row_as_mat(t), 0)) * trans(A);
  }
  
  
  //compute gradient
  arma::mat grad = arma::mat(J, J);
  arma::mat mid_term = Omega * B * Omega;
  
  
  arma::mat t1 = -0.5 * n * (Omega_Z + (m-1)*(Omega - trans(A) * Omega * A));
  arma::mat t2 = grad + 0.5 * Omega_Z * (trans(init_M.row_as_mat(0)) * init_M.row_as_mat(0) + S1_sum) * Omega_Z;
  arma::mat t3 = grad + 0.5 * mid_term;
  arma::mat t4 = grad - 0.5 * trans(A) * mid_term * A;
  
  grad = t1 + t2 + t3 + t4;
  
  return (scale*grad);
}

// [[Rcpp::export]] 
arma::vec Sigma_Z_constraint(
    const arma::vec & Sigma_Z_vec, // value of M at which to evaluate gradient
    const Rcpp::List & data, // List(Y, X, O, w)
    const Rcpp::List & params,
    const double scale) {
  
  const arma::cube & Y = Rcpp::as<arma::cube>(data["Y"]); // responses (m,J,n)
  int J = Y.n_cols;
  
  // initalize Sigma_Z matrix
  arma::mat Sigma_Z(J, J);
  std::copy(Sigma_Z_vec.begin(), Sigma_Z_vec.end(), Sigma_Z.begin()); 
  
  vec diag_vals(J);
  
  
  // compute diagonal entries of D matrix in LDL decomposition
  if (Sigma_Z.is_sympd()) {
    mat L;
    uvec P;
    chol(L, P, Sigma_Z, "lower", "vector");
    diag_vals = pow(L.diag(),2);
  } else {
    diag_vals.fill(-datum::inf);
  }
  
  return (-diag_vals);
}


// [[Rcpp::export]] 
arma::mat Sigma_Z_constraint_jac(
    const arma::vec & Sigma_Z_vec, // value of M at which to evaluate gradient
    const Rcpp::List & data, // List(Y, X, O, w)
    const Rcpp::List & params,
    const double scale) {
  
  const arma::cube & Y = Rcpp::as<arma::cube>(data["Y"]); // responses (m,J,n)
  int J = Y.n_cols;
  
  // initalize Sigma_Z matrix
  arma::mat Sigma_Z(J, J);
  std::copy(Sigma_Z_vec.begin(), Sigma_Z_vec.end(), Sigma_Z.begin()); 
  
  mat jac_res(J, J*J);
  cube jac(J, J, J, fill::zeros);
  jac(0,0,0) = 1;
  
  // compute jacobian for all diagonals beyond the first one
  if (Sigma_Z.is_sympd()) {
    for (int i=1; i < J; i++) {
      jac(i,i,i) = 1;
      mat Z_inv = Sigma_Z.submat(0,i-1,0,i-1).i();
      vec y = Sigma_Z.col(i).head(i);
      vec Z_inv_y = Z_inv * y;
      
      for (int j=0; j < i+1; j++) {
        for (int k = 0; k < i; k++) {
          if (j == i) {
            jac(i,j,k) = -2*Z_inv_y(k);
            jac(i,k,j) = -2*Z_inv_y(k);
          } else {
            jac(i, j, k) = -Z_inv_y(j) * Z_inv_y(k);
            jac(i, k, j) = -Z_inv_y(j) * Z_inv_y(k);
          }
          
        }
      }
    }
    
  } else {
    jac_res.fill(-datum::inf);
  }
  
  for (int l=0; l < J; l++) {
    jac_res.row(l) = vectorise(jac.slice(l)).t();
  }
  
  return (-jac_res);
}


// FUNCTIONS FOR NEW ELBO VI OPTIMIZATION SCHEME
// [[Rcpp::export]] 
double obj_function2(
    const Rcpp::List & data  , // List(Y, X, O, w)
    const Rcpp::List & params,
    const double scale = 1) {
  
  // Conversion from R to arma data structures
  
  const arma::cube & Y = Rcpp::as<arma::cube>(data["Y"]); // responses (m,J,n)
  const auto mu = Rcpp::as<arma::mat>(params["mu"]);     // true mean (J)
  const auto init_M = Rcpp::as<arma::cube>(params["M"]);     // variational approximation mean (m,J,n)
  const auto init_S = Rcpp::as<arma::cube>(params["S"]);     // variational approximation variances (m,J,n)
  // const auto Sigma_Z = Rcpp::as<arma::mat>(params["Sigma_Z"]); // covinv for particular VAR instance (p,p) 
  const auto Sigma = Rcpp::as<arma::mat>(params["Sigma"]); // covinv of noise in VAR process (p,p)
  const auto A = Rcpp::as<arma::mat>(params["A"]); // VAR(1) coefficient matrix (p,p)
  
  // set up useful data structures and base quantities
  arma::mat mu_rep = repelem(mu, Y.n_rows, 1);
  arma::cube M_and_S_half = init_M + 0.5 * init_S; 
  arma::mat Omega = inv(Sigma);
  
  arma::mat M1 = init_M.row_as_mat(0); 
  arma::cube M1_M1t_S1 = cube(Y.n_cols, Y.n_cols, Y.n_slices);
  arma::cube Mt_Mt1 = cube(Y.n_rows-1, Y.n_rows-1, Y.n_slices);
  
  for (int i=0; i != Y.n_slices; ++i) {
    M1_M1t_S1.slice(i) = trans(M1.row(i)) * M1.row(i) + diagmat(init_S.slice(i).row(0));
  }
  
  arma::mat M1t_M1_S1_sum = sum(M1_M1t_S1,2);
  
  arma::mat M_temp = mat(init_M.n_cols, init_M.n_rows);
  arma::mat S_temp = mat(init_S.n_cols, init_S.n_cols);
  arma::mat T4_temp = mat(init_S.n_cols, init_S.n_cols);
  
  double term_4 = 0;
  
  for (int i=0; i!= Y.n_slices; ++i) {
    M_temp = trans(init_M.slice(i));
    S_temp = init_S.slice(i);
    for (int t=0; t != Y.n_rows-1; ++t) {
      
      T4_temp = (M_temp.col(t+1) - A*M_temp.col(t)) * trans((M_temp.col(t+1) - A*M_temp.col(t)))  + diagmat(S_temp.row(t+1)) + A * diagmat(S_temp.row(t)) * trans(A);
      term_4 = term_4 - 0.5*trace(T4_temp * Omega);
    }
  }
  
  // Compute relevant quantities
  double term_1 = accu(Y % (init_M.each_slice() + mu_rep) - exp(M_and_S_half.each_slice() + mu_rep));
  double term_2 = 0.5*(Y.n_cols * (Y.n_rows - 1) * Y.n_slices - Y.n_slices*(Y.n_rows - 1)*log_det(Sigma).real() + accu(log(init_S)) - accu(log(init_S.row_as_mat(0))));
  
  // print statements for debugging
  /*
   cout << "Term 1: " << term_1 << "\n";
   cout << "Term 2: " << term_2 << "\n";
   cout << "Term 3: " << term_3 << "\n";
   cout << "Term 4: " << term_4 << "\n";
   */
  
  return (scale*(term_1 + term_2 + term_4));
}

// [[Rcpp::export]] 
double obj_function2_cov(
    const Rcpp::List & data  , // List(Y, X, O, w)
    const Rcpp::List & params,
    const double scale = 1) {
  
  // Conversion from R to arma data structures
  const arma::cube & Y = Rcpp::as<arma::cube>(data["Y"]); // responses (m,J,n)
  const arma::cube & X = Rcpp::as<arma::cube>(data["X"]); // responses (m,p,n)
  const arma::mat & O = Rcpp::as<arma::mat>(data["O"]); // offsets (m, n)
  const auto init_M = Rcpp::as<arma::cube>(params["M"]);     // variational approximation mean (m,J,n)
  const auto init_S = Rcpp::as<arma::cube>(params["S"]);     // variational approximation variances (m,J,n)
  const auto Sigma = Rcpp::as<arma::mat>(params["Sigma"]); // cov of noise in VAR process (J,J)
  const auto A = Rcpp::as<arma::mat>(params["A"]); // VAR(1) coefficient matrix (J,J)
  const auto beta = Rcpp::as<arma::mat>(params["Beta"]);     // coefficient matrix for all categories (p+1 x J)
  
  
  // set up useful data structures and base quantities
  arma::cube M_and_S_half = init_M + 0.5 * init_S; 
  arma::mat Omega = inv(Sigma);
  arma::cube mu_cube = cube(size(Y));
  // populate mu cube appropriately by looping over all covariates and computing their specific mean
  for (int i=0; i != Y.n_slices; ++i) {
    mu_cube.slice(i) = X.slice(i) * beta.rows(1,beta.n_rows-1) + repelem(beta.row(0), Y.n_rows, 1);
  }
  // add offset to mu_cube value
  for (int i=0; i != Y.n_cols; ++i) {
    mu_cube.col(i) = mu_cube.col_as_mat(i) + O;
  }

  arma::mat M1 = init_M.row_as_mat(0); 
  arma::cube M1_M1t_S1 = cube(Y.n_cols, Y.n_cols, Y.n_slices);
  arma::cube Mt_Mt1 = cube(Y.n_rows-1, Y.n_rows-1, Y.n_slices);
  
  for (int i=0; i != Y.n_slices; ++i) {
    M1_M1t_S1.slice(i) = trans(M1.row(i)) * M1.row(i) + diagmat(init_S.slice(i).row(0));
  }
  
  arma::mat M1t_M1_S1_sum = sum(M1_M1t_S1,2);
  
  arma::mat M_temp = mat(init_M.n_cols, init_M.n_rows);
  arma::mat S_temp = mat(init_S.n_cols, init_S.n_cols);
  arma::mat T3_temp = mat(init_S.n_cols, init_S.n_cols);
  
  double term_3 = 0;
  
  for (int i=0; i!= Y.n_slices; ++i) {
    M_temp = trans(init_M.slice(i));
    S_temp = init_S.slice(i);
    for (int t=0; t != Y.n_rows-1; ++t) {
      
      T3_temp = (M_temp.col(t+1) - A*M_temp.col(t)) * trans((M_temp.col(t+1) - A*M_temp.col(t)))  + diagmat(S_temp.row(t+1)) + A * diagmat(S_temp.row(t)) * trans(A);
      term_3 = term_3 - 0.5*trace(T3_temp * Omega);
    }
  }
  
  // Compute relevant quantities
  double term_1 = accu(Y % (init_M + mu_cube) - exp(M_and_S_half + mu_cube));
  double term_2 = 0.5*(Y.n_cols * (Y.n_rows - 1) * Y.n_slices - Y.n_slices*(Y.n_rows - 1)*log_det(Sigma).real() + accu(log(init_S)) - accu(log(init_S.row_as_mat(0))));
  
  // print statements for debugging
  /*
   cout << "Term 1: " << term_1 << "\n";
   cout << "Term 2: " << term_2 << "\n";
   cout << "Term 3: " << term_3 << "\n";
   */
  
  
  
  return (scale*(term_1 + term_2 + term_3));
}

// [[Rcpp::export]]
arma::mat beta_grad2(
    const arma::vec & beta_vec, // value of beta (not including row for intercepts!) at which to evaluate gradient (as vector)
    const Rcpp::List & data  , // List(Y, X, O)
    const Rcpp::List & params,
    const double scale = 1
) {
  // Conversion from R to arma data structures
  const arma::cube & Y = Rcpp::as<arma::cube>(data["Y"]); // responses (m,J,n)
  const arma::cube & X = Rcpp::as<arma::cube>(data["X"]); // covariates (m,p,n)
  const arma::mat & O = Rcpp::as<arma::mat>(data["O"]); // offsets (m,n)
  const auto init_M = Rcpp::as<arma::cube>(params["M"]);     // variational approximation means (m,J,n)
  const auto init_S = Rcpp::as<arma::cube>(params["S"]);     // variational approximation variances (m,J,n)
  const arma::vec & beta_0 = Rcpp::as<arma::mat>(params["Beta"]).row(0).as_col(); //intercept parameters (J) obtained as first row of beta in supplied params list
  int J = Y.n_cols, p = X.n_cols;
  
  //initialize beta matrix for evaluating gradient at input value of beta
  arma::mat beta(p, J);
  std::copy(beta_vec.begin(), beta_vec.end(), beta.begin());

  //make matrix to store gradient to be returned
  arma::mat grad = arma::mat(p, J);
  
  //set up useful quantities for computations
  arma::cube M_and_S_half = init_M + 0.5 * init_S;
  arma::cube xbeta = X.each_slice() * beta;
  
  //compute gradient
  for (int j=0; j != J; ++j) {
    for (int k=0; k != p; ++k) {
      grad(k,j) = accu(X.col_as_mat(k) % (Y.col_as_mat(j) - exp(xbeta.col_as_mat(j) + M_and_S_half.col_as_mat(j) + beta_0(j) + O)));
    }
  }
  
  return(scale*grad);
}

// [[Rcpp::export]] 
arma::cube M_grad2(
    const arma::vec & M_vec, // value of M at which to evaluate gradient
    const Rcpp::List & data  , // List(Y, X, O)
    const Rcpp::List & params,
    const double scale = 1) {
  
  // Conversion from R to arma data structures
  const arma::cube & Y = Rcpp::as<arma::cube>(data["Y"]); // responses (m,J,n)
  int m = Y.n_rows, J = Y.n_cols, n = Y.n_slices;
  /*
   const arma::mat & X = Rcpp::as<arma::mat>(data["X"]); // covariates (n,d)
   const arma::mat & O = Rcpp::as<arma::mat>(data["O"]); // offsets (n,J)
   const arma::vec & w = Rcpp::as<arma::vec>(data["w"]); // weights (n) 
   */ 
  arma::cube init_M(m, J, n);
  std::copy(M_vec.begin(), M_vec.end(), init_M.begin()); 
  const arma::mat mu = Rcpp::as<arma::mat>(params["mu"]);     // true mean (p)
  const auto init_S = Rcpp::as<arma::cube>(params["S"]);     // variational approximation variances (m,J,n)
  const auto Sigma = Rcpp::as<arma::mat>(params["Sigma"]); // covinv for particular VAR instance (p,p) 
  const auto A = Rcpp::as<arma::mat>(params["A"]); // VAR(1) coefficient matrix (p,p)
  
  // set up useful data structures and base quantities
  arma::mat mu_rep = repelem(mu, Y.n_rows, 1);
  arma::cube M_and_S_half = init_M + 0.5 * init_S; 
  arma::mat Omega = inv(Sigma);
  
  //compute gradient
  arma::cube grad = arma::cube(m, J, n);
  
  for (int t=0; t != m; ++t) {
    //compute gradient value
    mat temp_grad(n, J);
    temp_grad = Y.row_as_mat(t) - exp(M_and_S_half.each_slice() + mu_rep).row_as_mat(t);
    
    if (t == 0) {
      temp_grad = temp_grad - (init_M.row_as_mat(0)*trans(A) - init_M.row_as_mat(1))*Omega*A;
    }
    
    
    if (0 < t & t < m-1) {
      temp_grad = temp_grad - (init_M.row_as_mat(t)*trans(A) - init_M.row_as_mat(t+1))*Omega*A;
      temp_grad = temp_grad - (init_M.row_as_mat(t) - init_M.row_as_mat(t-1)*trans(A))*Omega;
    }
    
    
    if (t == m-1) {
      temp_grad = temp_grad - (init_M.row_as_mat(t) - init_M.row_as_mat(t-1)*trans(A))*Omega;
    }
    
    //copy to object to be returned
    for(int j=0; j!= J; ++j) {
      for(int i=0; i != n; ++i) {
        grad(t,j,i) = temp_grad(i,j);
      }
    }
  }
  
  return (scale*grad);
}

// [[Rcpp::export]] 
arma::cube M_grad2_cov(
    const arma::vec & M_vec, // value of M at which to evaluate gradient
    const Rcpp::List & data  , // List(Y, X, O)
    const Rcpp::List & params,
    const double scale = 1) {
  
  // Conversion from R to arma data structures
  const arma::cube & Y = Rcpp::as<arma::cube>(data["Y"]); // responses (m,J,n)
  int m = Y.n_rows, J = Y.n_cols, n = Y.n_slices;
  const arma::cube & X = Rcpp::as<arma::cube>(data["X"]); // covariates (m, p, n)
  const arma::mat & O = Rcpp::as<arma::mat>(data["O"]); // offsets (m,n)
  arma::cube init_M(m, J, n);
  std::copy(M_vec.begin(), M_vec.end(), init_M.begin()); 
  const arma::mat beta = Rcpp::as<arma::mat>(params["Beta"]);     // coefficient matrix (p, J)
  const auto init_S = Rcpp::as<arma::cube>(params["S"]);     // variational approximation variances (m,J,n)
  const auto Sigma = Rcpp::as<arma::mat>(params["Sigma"]); // covinv for particular VAR instance (J,J) 
  const auto A = Rcpp::as<arma::mat>(params["A"]); // VAR(1) coefficient matrix (J,J)
  
  // set up useful data structures and base quantities
  arma::cube mu_cube = cube(size(Y));
  // populate mu cube appropriately by looping over all covariates and computing their specific mean
  for (int i=0; i != Y.n_slices; ++i) {
    mu_cube.slice(i) = X.slice(i) * beta.rows(1,beta.n_rows-1) + repelem(beta.row(0), Y.n_rows, 1);
  }
  // add offset to mu_cube value
  for (int i=0; i != Y.n_cols; ++i) {
    mu_cube.col(i) = mu_cube.col_as_mat(i) + O;
  }
  arma::cube M_and_S_half = init_M + 0.5 * init_S; 
  arma::mat Omega = inv(Sigma);
  
  //compute gradient
  arma::cube grad = arma::cube(m, J, n);
  
  for (int t=0; t != m; ++t) {
    //compute gradient value
    mat temp_grad(n, J);
    temp_grad = Y.row_as_mat(t) - exp(M_and_S_half + mu_cube).row_as_mat(t);
    
    if (t == 0) {
      temp_grad = temp_grad - (init_M.row_as_mat(0)*trans(A) - init_M.row_as_mat(1))*Omega*A;
    }
    
    
    if (0 < t & t < m-1) {
      temp_grad = temp_grad - (init_M.row_as_mat(t)*trans(A) - init_M.row_as_mat(t+1))*Omega*A;
      temp_grad = temp_grad - (init_M.row_as_mat(t) - init_M.row_as_mat(t-1)*trans(A))*Omega;
    }
    
    
    if (t == m-1) {
      temp_grad = temp_grad - (init_M.row_as_mat(t) - init_M.row_as_mat(t-1)*trans(A))*Omega;
    }
    
    //copy to object to be returned
    for(int j=0; j!= J; ++j) {
      for(int i=0; i != n; ++i) {
        grad(t,j,i) = temp_grad(i,j);
      }
    }
  }
  
  return (scale*grad);
}

// [[Rcpp::export]] 
arma::cube S_grad2(
    const arma::vec & S_vec, // value of M at which to evaluate gradient
    const Rcpp::List & data  , // List(Y, X, O, w)
    const Rcpp::List & params,
    const double scale = 1) {
  
  // Conversion from R to arma data structures
  
  const arma::cube & Y = Rcpp::as<arma::cube>(data["Y"]); // responses (m,J,n)
  int m = Y.n_rows, J = Y.n_cols, n = Y.n_slices;
  /*
   const arma::mat & X = Rcpp::as<arma::mat>(data["X"]); // covariates (n,d)
   const arma::mat & O = Rcpp::as<arma::mat>(data["O"]); // offsets (n,J)
   const arma::vec & w = Rcpp::as<arma::vec>(data["w"]); // weights (n) 
   */ 
  arma::cube init_S(m, J, n);
  std::copy(S_vec.begin(), S_vec.end(), init_S.begin()); 
  const arma::mat mu = Rcpp::as<arma::mat>(params["mu"]);     // true mean (p)
  const auto init_M = Rcpp::as<arma::cube>(params["M"]);     // variational approximation means (m,J,n)
  const auto Sigma = Rcpp::as<arma::mat>(params["Sigma"]); // covinv for particular VAR instance (p,p) 
  const auto A = Rcpp::as<arma::mat>(params["A"]); // VAR(1) coefficient matrix (p,p)
  
  // set up useful data structures and base quantities
  arma::mat mu_rep = repelem(mu, Y.n_rows, 1);
  arma::cube M_and_S_half = init_M + 0.5 * init_S; 
  arma::mat Omega = inv(Sigma);
  arma::mat At_Omega_A = trans(A) * Omega * A;
  //compute gradient
  arma::cube grad = arma::cube(m, J, n);
  
  for (int t=0; t != m; ++t) {
    //compute gradient value
    mat temp_grad(n, J);
    temp_grad = - 0.5*exp(M_and_S_half.each_slice() + mu_rep).row_as_mat(t) + 0.5*pow(init_S.row_as_mat(t), -1);
    
    if (t == 0) {
      temp_grad = temp_grad - 0.5*mat(n,1,fill::ones)*trans(At_Omega_A.diag()) - 0.5*pow(init_S.row_as_mat(t), -1);
    }
    
    
    if (0 < t & t < m-1) {
      temp_grad = temp_grad - 0.5*mat(n,1,fill::ones)*trans(Omega.diag() + At_Omega_A.diag());
    }
    
    
    if (t == m-1) {
      temp_grad = temp_grad - 0.5*mat(n,1,fill::ones)*trans(Omega.diag());
    }
    
    //copy to object to be returned
    for(int j=0; j!= J; ++j) {
      for(int i=0; i != n; ++i) {
        grad(t,j,i) = temp_grad(i,j);
      }
    }
  }
  
  return (scale*grad);
}

// [[Rcpp::export]] 
arma::cube S_grad2_cov(
    const arma::vec & S_vec, // value of M at which to evaluate gradient
    const Rcpp::List & data  , // List(Y, X, O, w)
    const Rcpp::List & params,
    const double scale = 1) {
  
  // Conversion from R to arma data structures
  
  const arma::cube & Y = Rcpp::as<arma::cube>(data["Y"]); // responses (m,J,n)
  int m = Y.n_rows, J = Y.n_cols, n = Y.n_slices;
  const arma::cube & X = Rcpp::as<arma::cube>(data["X"]); // covariates (m, p, n)
  const arma::mat & O = Rcpp::as<arma::mat>(data["O"]); // offsets (m,n)
  
  arma::cube init_S(m, J, n);
  std::copy(S_vec.begin(), S_vec.end(), init_S.begin()); 
  const arma::mat beta = Rcpp::as<arma::mat>(params["Beta"]);     // coefficient matrix (p, J)
  const auto init_M = Rcpp::as<arma::cube>(params["M"]);     // variational approximation means (m,J,n)
  const auto Sigma = Rcpp::as<arma::mat>(params["Sigma"]); // covinv for particular VAR instance (p,p) 
  const auto A = Rcpp::as<arma::mat>(params["A"]); // VAR(1) coefficient matrix (p,p)
  
  // set up useful data structures and base quantities
  // set up useful data structures and base quantities
  arma::cube mu_cube = cube(size(Y));
  // populate mu cube appropriately by looping over all covariates and computing their specific mean
  for (int i=0; i != Y.n_slices; ++i) {
    mu_cube.slice(i) = X.slice(i) * beta.rows(1,beta.n_rows-1) + repelem(beta.row(0), Y.n_rows, 1);
  }
  // add offset to mu_cube value
  for (int i=0; i != Y.n_cols; ++i) {
    mu_cube.col(i) = mu_cube.col_as_mat(i) + O;
  }
  arma::cube M_and_S_half = init_M + 0.5 * init_S; 
  arma::mat Omega = inv(Sigma);
  arma::mat At_Omega_A = trans(A) * Omega * A;
  //compute gradient
  arma::cube grad = arma::cube(m, J, n);
  
  for (int t=0; t != m; ++t) {
    //compute gradient value
    mat temp_grad(n, J);
    temp_grad = - 0.5*exp(M_and_S_half + mu_cube).row_as_mat(t) + 0.5*pow(init_S.row_as_mat(t), -1);
    
    if (t == 0) {
      temp_grad = temp_grad - 0.5*mat(n,1,fill::ones)*trans(At_Omega_A.diag()) - 0.5*pow(init_S.row_as_mat(t), -1);
    }
    
    
    if (0 < t & t < m-1) {
      temp_grad = temp_grad - 0.5*mat(n,1,fill::ones)*trans(Omega.diag() + At_Omega_A.diag());
    }
    
    
    if (t == m-1) {
      temp_grad = temp_grad - 0.5*mat(n,1,fill::ones)*trans(Omega.diag());
    }
    
    //copy to object to be returned
    for(int j=0; j!= J; ++j) {
      for(int i=0; i != n; ++i) {
        grad(t,j,i) = temp_grad(i,j);
      }
    }
  }
  
  return (scale*grad);
}