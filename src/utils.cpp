#include <RcppArmadillo.h>

//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::plugins(cpp11)]]

// calculate gene size-----------------------------------------------------
unsigned int cal_gene_size(const std::vector<unsigned int> block_start_vec,
                           const std::vector<unsigned int> block_end_vec) {

  unsigned int res = 0;
  unsigned int temp;
  std::vector<unsigned int>::const_iterator it;
  for (it = block_start_vec.begin(); it != block_start_vec.end(); ++it) {
    temp = (*std::next(block_end_vec.begin(), std::distance(block_start_vec.begin(), it))) - (*it) + 1;
    res = res + temp;
  }
  return res;
}

// create sequence from start and end vectors----------------------------------
std::vector<unsigned int> seq_v(const std::vector<unsigned int>& x,
                                const std::vector<unsigned int>& y) {
  std::vector<unsigned int> res;
  std::vector<unsigned int>::const_iterator x_it, y_it;
  unsigned int i;
  for (x_it = x.begin(); x_it != x.end(); ++x_it) {
    y_it = std::next(y.begin(), std::distance(x.begin(), x_it));
    for (i = (*x_it); i <= (*y_it); ++i) {
      res.push_back(i);
    }
  }
  return res;
}

// get backward subset of a read type------------------------------------------------------
std::vector<std::string> get_backward_subset(const std::string& s) {

  // search for position of '-' in s
  std::vector<std::string::size_type> index;
  std::string::const_iterator s_it;
  for (s_it = s.begin(); s_it != s.end(); ++s_it) {
    if ((*s_it) == '-') {
      index.push_back(std::distance(s.begin(), s_it));
    }
  }

  // get backward substring from s
  std::vector<std::string> res(index.size() + 1);
  std::vector<std::string>::iterator res_it;
  std::vector<std::string::size_type>::iterator index_it = index.begin();
  res[0] = s;
  for (res_it = std::next(res.begin(), 1); res_it != res.end(); ++res_it) {
    (*res_it) = s.substr((*index_it) + 1, s.size() - (*index_it) - 1);
    index_it = std::next(index_it, 1);
  }
  return res;
}

// get forward subset of a read type----------------------------------------------------------
std::vector<std::string> get_forward_subset(const std::string& s) {

  // search for position of '-' in s
  std::vector<std::string::size_type> index;
  std::string::const_iterator s_it;
  for (s_it = s.begin(); s_it != s.end() ; ++s_it) {
    if ((*s_it) == '-') {
      index.push_back(std::distance(s.begin(), s_it));
    }
  }

  // get forward substring from s
  std::vector<std::string> res(index.size() + 1);
  std::vector<std::string>::iterator res_it;
  std::vector<std::string::size_type>::iterator index_it = index.begin();
  for (res_it = res.begin(); res_it != std::prev(res.end(), 1); ++res_it) {
    (*res_it) = s.substr(0, (*index_it));
    index_it = std::next(index_it, 1);
  }
  res.back() = s;
  return res;
}

// check if two read type can connect together---------------------------------------------
bool check_connection(const std::string& s_small,
                      const std::string& s_big) {

  std::vector<std::string> s1_backward = get_backward_subset(s_small);
  std::vector<std::string> s2_forward = get_forward_subset(s_big);
  std::sort(s1_backward.begin(), s1_backward.end());
  std::sort(s2_forward.begin(), s2_forward.end());
  std::vector<std::string> res(s1_backward.size() + s2_forward.size(), "");
  std::vector<std::string>::iterator it;
  it = std::set_intersection(s1_backward.begin(), s1_backward.end(), s2_forward.begin(), s2_forward.end(), res.begin());
  if (res[0] != "") {
    return true;
  }
  else {
    return false;
  }
}

// create graph matrix----------------------------------------------------------------------
arma::imat create_graph_matrix(const std::vector<std::string>& g_read_type_vec) {

  arma::imat res(g_read_type_vec.size(), g_read_type_vec.size(), arma::fill::zeros);
  std::vector<std::string>::const_iterator i, j;
  for (i = g_read_type_vec.begin(); i != g_read_type_vec.end(); ++i) {
    std::vector<std::string> connection_nodes_reminder;
    for (j = std::next(i, 1); j != g_read_type_vec.end(); ++j) {
      if (check_connection(*i, *j)) {
        bool valid = true;
        std::vector<std::string>::iterator k;
        for (k = connection_nodes_reminder.begin(); k != connection_nodes_reminder.end(); ++k) {
          if (check_connection(*k, *j)) {
            valid = false;
            break;
          }
        }
        if (valid) {
          res(std::distance(g_read_type_vec.begin(), i), std::distance(g_read_type_vec.begin(), j)) = 1;
          connection_nodes_reminder.push_back(*j);
        }
      }
    }
  }
  return res;
}

// create read type group for one gene----------------------------------------------------------
arma::uvec create_read_type_group_g(const std::vector<std::string>& read_type,
                                    const arma::imat& graph_matrix,
                                    const double min_reduce) {

  // containers for result saving
  arma::uvec read_type_group(read_type.size(), arma::fill::zeros);
  arma::uvec read_type_group_default(read_type.size(), arma::fill::ones);
  // remember grouped read types
  std::vector<std::string> read_type_grouped;
  unsigned int group_id, idx, idx_curr, idx_next;
  arma::uvec idx_next_vec;
  group_id = 1;

  // loop for each read type, to search read types that can be grouped with it
  for (idx = 0; idx < read_type.size(); ++idx) {
    std::vector<std::string>::iterator find_grouped_it;
    find_grouped_it = std::find(read_type_grouped.begin(), read_type_grouped.end(), read_type[idx]);
    if (find_grouped_it == read_type_grouped.end()) {
      // save index of read types that can be grouped together
      arma::uvec read_type_group_idx;
      arma::uvec new_idx(1);
      new_idx[0] = idx;
      read_type_grouped.push_back(read_type[idx]);
      read_type_group_idx.insert_rows(read_type_group_idx.n_rows, new_idx);
      idx_curr = idx;
      // keep going to find all read types that can be grouped with read_type_vec[idx]
      // check if out degree equals 1
      while (arma::sum(graph_matrix.row(idx_curr)) == 1) {
        idx_next_vec = arma::find(graph_matrix.row(idx_curr) == 1);
        idx_next = idx_next_vec(0);
        // check if in degree equals 1
        if (arma::sum(graph_matrix.col(idx_next)) == 1) {
          new_idx[0] = idx_next;
          read_type_group_idx.insert_rows(read_type_group_idx.n_rows, new_idx);
          read_type_grouped.push_back(read_type[idx_next]);
          idx_curr = idx_next;
        }
        else {
          break;
        }
      }
      read_type_group(read_type_group_idx).fill(group_id);
      group_id++;
    }
  }
  if (((double)read_type.size() - (double)group_id + 1) / (double)read_type.size() < min_reduce){
    return read_type_group_default;
  }
  else{
    return read_type_group;
  }
}

// calculate initial value----------------------------------------------------------------------------
Rcpp::List calculate_initial_val(Rcpp::NumericMatrix g_data_matrix,
                                 Rcpp::NumericVector lib_size_norm,
                                 Rcpp::IntegerVector group_sample_num) {

  int group0_sample_num = group_sample_num[0];
  int group1_sample_num = group_sample_num[1];

  // create theta matrix
  Rcpp::NumericMatrix g_theta_matrix(g_data_matrix.nrow(), g_data_matrix.ncol());
  for (unsigned int i = 0; i < g_data_matrix.nrow(); ++i) {
    g_theta_matrix.row(i) = (g_data_matrix.row(i) + 0.01) / lib_size_norm;
  }

  // calculate initial value H0 U H1
  Rcpp::NumericMatrix g_theta_matrix_0 = g_theta_matrix(Rcpp::_, Rcpp::Range(0, group0_sample_num - 1));
  Rcpp::NumericMatrix g_theta_matrix_1 = g_theta_matrix(Rcpp::_, Rcpp::Range(group0_sample_num, g_theta_matrix.ncol() - 1));
  Rcpp::NumericVector init_val_theta_0(g_theta_matrix.nrow());
  Rcpp::NumericVector init_val_theta_1(g_theta_matrix.nrow());
  Rcpp::NumericVector theta_var_0(g_theta_matrix.nrow());
  Rcpp::NumericVector theta_var_1(g_theta_matrix.nrow());
  for (unsigned int i = 0; i < g_theta_matrix.nrow(); ++i) {
    init_val_theta_0[i] = Rcpp::mean(g_theta_matrix_0.row(i));
    init_val_theta_1[i] = Rcpp::mean(g_theta_matrix_1.row(i));
    theta_var_0[i] = Rcpp::var(g_theta_matrix_0.row(i));
    theta_var_1[i] = Rcpp::var(g_theta_matrix_1.row(i));
  }
  Rcpp::NumericVector init_val_phi = (theta_var_0 / Rcpp::pow(init_val_theta_0, 2) +
                                      theta_var_1 / Rcpp::pow(init_val_theta_1, 2)) / 2;

  // calculate initial value H0
  Rcpp::NumericVector init_val_theta_H0(g_theta_matrix.nrow());
  Rcpp::NumericVector theta_var_H0(g_theta_matrix.nrow());
  for (unsigned int i = 0; i < g_theta_matrix.nrow(); ++i) {
    init_val_theta_H0[i] = Rcpp::mean(g_theta_matrix.row(i));
    theta_var_H0[i] = Rcpp::var(g_theta_matrix.row(i));
  }
  Rcpp::NumericVector init_val_phi_H0 = theta_var_H0 / Rcpp::pow(init_val_theta_H0, 2);

  // result
  Rcpp::List res = Rcpp::List::create(Rcpp::_["init_val_theta_0"] = init_val_theta_0,
                                      Rcpp::_["init_val_theta_1"] = init_val_theta_1,
                                      Rcpp::_["init_val_phi"] = init_val_phi,
                                      Rcpp::_["init_val_theta_H0"] = init_val_theta_H0,
                                      Rcpp::_["init_val_phi_H0"] = init_val_phi_H0);
  return res;
}

// overload "+" for matrix + matrix--------------------------------------------------------------------
Rcpp::NumericMatrix operator+(Rcpp::NumericMatrix x, Rcpp::NumericMatrix y) {
  Rcpp::NumericMatrix res(x.nrow(), x.ncol());
  for (unsigned int i = 0; i < x.nrow(); ++i) {
    res.row(i) = x.row(i) + y.row(i);
  }
  return res;
}

// overload "*" for matrix * matrix------------------------------------------------------------------
Rcpp::NumericMatrix operator*(Rcpp::NumericMatrix x, Rcpp::NumericMatrix y) {
  Rcpp::NumericMatrix res(x.nrow(), x.ncol());
  for (unsigned int i = 0; i < x.nrow(); ++i) {
    res.row(i) = x.row(i) * y.row(i);
  }
  return res;
}

// overload "/" for matrix / matrix-----------------------------------------------------
Rcpp::NumericMatrix operator/(Rcpp::NumericMatrix x, Rcpp::NumericMatrix y) {
  Rcpp::NumericMatrix res(x.nrow(), x.ncol());
  for (unsigned int i = 0; i < x.nrow(); ++i) {
    res.row(i) = x.row(i) / y.row(i);
  }
  return res;
}

// log for matrix------------------------------------------------------------------------
Rcpp::NumericMatrix log_mat(Rcpp::NumericMatrix x) {
  Rcpp::NumericMatrix res(x.nrow(), x.ncol());
  for (unsigned int i = 0; i < x.nrow(); ++i) {
    res.row(i) = Rcpp::log(x.row(i));
  }
  return res;
}

// lgamma for matrix-----------------------------------------------------------------------------
Rcpp::NumericMatrix lgamma_mat(Rcpp::NumericMatrix x) {
  Rcpp::NumericMatrix res(x.nrow(), x.ncol());
  for (unsigned int i = 0; i < x.nrow(); ++i) {
    res.row(i) = Rcpp::lgamma(x.row(i));
  }
  return res;
}

// digamma for matrix-------------------------------------------------------------------
Rcpp::NumericMatrix digamma_mat(Rcpp::NumericMatrix x) {
  Rcpp::NumericMatrix res(x.nrow(), x.ncol());
  for (unsigned int i = 0; i < x.nrow(); ++i) {
    res.row(i) = Rcpp::digamma(x.row(i));
  }
  return res;
}

// sum for matrix---------------------------------------------------------------------
double sum_mat(Rcpp::NumericMatrix x) {
  double res = 0;
  for (unsigned int i = 0; i < x.nrow(); ++i) {
    res = res + Rcpp::sum(x.row(i));
  }
  return res;
}

// log likelihood function------------------------------------------------------------------------
double ll(Rcpp::NumericMatrix g_data_matrix,
          Rcpp::NumericVector lib_size_norm,
          Rcpp::IntegerVector group_sample_num,
          Rcpp::NumericVector theta0,
          Rcpp::NumericVector theta1,
          Rcpp::NumericVector phi) {

  int group0_sample_num = group_sample_num[0];
  int group1_sample_num = group_sample_num[1];

  // prepare matrix
  Rcpp::NumericMatrix theta_matrix(g_data_matrix.nrow(), g_data_matrix.ncol());
  Rcpp::NumericMatrix phi_matrix(g_data_matrix.nrow(), g_data_matrix.ncol());
  Rcpp::NumericMatrix phi_matrix_inv(g_data_matrix.nrow(), g_data_matrix.ncol());
  Rcpp::NumericMatrix lib_size_norm_matrix(g_data_matrix.nrow(), g_data_matrix.ncol());
  for (unsigned int i = 0; i < group0_sample_num; ++i) {
    theta_matrix.column(i) = theta0;
    phi_matrix.column(i) = phi;
    phi_matrix_inv.column(i) = Rcpp::pow(phi, -1);
  }
  for (unsigned int i = group0_sample_num; i < g_data_matrix.ncol(); ++i) {
    theta_matrix.column(i) = theta1;
    phi_matrix.column(i) = phi;
    phi_matrix_inv.column(i) = Rcpp::pow(phi, -1);
  }
  for (unsigned int i = 0; i < g_data_matrix.nrow(); ++i) {
    lib_size_norm_matrix.row(i) = lib_size_norm;
  }

  // calculate ll
  double res;
  res = sum_mat(lgamma_mat(g_data_matrix + phi_matrix_inv)) - sum_mat(lgamma_mat(g_data_matrix + 1)) -
        (group0_sample_num + group1_sample_num) * Rcpp::sum(Rcpp::lgamma(Rcpp::pow(phi, -1))) -
        sum_mat(phi_matrix_inv * log_mat((phi_matrix * lib_size_norm_matrix * theta_matrix) + 1)) +
        sum_mat(g_data_matrix * log_mat((lib_size_norm_matrix * theta_matrix) /
                                        ((lib_size_norm_matrix * theta_matrix) + phi_matrix_inv)));
  return res;
}

// log likelihood function H0-----------------------------------------------------------------------
double ll_H0(Rcpp::NumericMatrix g_data_matrix,
             Rcpp::NumericVector lib_size_norm,
             Rcpp::IntegerVector group_sample_num,
             Rcpp::NumericVector theta,
             Rcpp::NumericVector phi) {

  Rcpp::NumericMatrix theta_matrix(g_data_matrix.nrow(), g_data_matrix.ncol());
  Rcpp::NumericMatrix phi_matrix(g_data_matrix.nrow(), g_data_matrix.ncol());
  Rcpp::NumericMatrix phi_matrix_inv(g_data_matrix.nrow(), g_data_matrix.ncol());
  Rcpp::NumericMatrix lib_size_norm_matrix(g_data_matrix.nrow(), g_data_matrix.ncol());
  for (unsigned int i = 0; i < g_data_matrix.ncol(); ++i) {
    theta_matrix.column(i) = theta;
    phi_matrix.column(i) = phi;
    phi_matrix_inv.column(i) = Rcpp::pow(phi, -1);
  }
  for (unsigned int i = 0; i < g_data_matrix.nrow(); ++i) {
    lib_size_norm_matrix.row(i) = lib_size_norm;
  }

  // calculate ll_H0
  double res;
  res = sum_mat(lgamma_mat(g_data_matrix + phi_matrix_inv)) - sum_mat(lgamma_mat(g_data_matrix + 1)) -
        (group_sample_num[0] + group_sample_num[1]) * Rcpp::sum(Rcpp::lgamma(Rcpp::pow(phi, -1))) -
        sum_mat(phi_matrix_inv * log_mat((phi_matrix * lib_size_norm_matrix * theta_matrix) + 1)) +
        sum_mat(g_data_matrix * log_mat((lib_size_norm_matrix * theta_matrix) /
                                          ((lib_size_norm_matrix * theta_matrix) + phi_matrix_inv)));
  return res;
}

// estimate parameter-------------------------------------------------------------------------
Rcpp::List estimate_par(Rcpp::NumericMatrix g_data_matrix,
                        Rcpp::NumericVector lib_size_norm,
                        Rcpp::IntegerVector group_sample_num,
                        Rcpp::NumericVector init_val_theta0,
                        Rcpp::NumericVector init_val_theta1,
                        Rcpp::NumericVector init_val_phi,
                        double tol,
                        int times) {

  int group0_sample_num = group_sample_num[0];
  int group1_sample_num = group_sample_num[1];
  Rcpp::NumericVector theta0_old, theta1_old, theta0_new, theta1_new, phi_old, phi_new;
  theta0_old = Rcpp::clone(init_val_theta0);
  theta1_old = Rcpp::clone(init_val_theta1);
  phi_old = Rcpp::clone(init_val_phi);
  theta0_new = Rcpp::clone(theta0_old);
  theta1_new = Rcpp::clone(theta1_old);
  phi_new = Rcpp::clone(phi_old);
  Rcpp::NumericMatrix lib_size_norm_matrix(g_data_matrix.nrow(), g_data_matrix.ncol());
  for (unsigned int i = 0; i < g_data_matrix.nrow(); ++i) {
    lib_size_norm_matrix.row(i) = lib_size_norm;
  }

  // estimate
  Rcpp::NumericMatrix theta_matrix(g_data_matrix.nrow(), g_data_matrix.ncol());
  Rcpp::NumericMatrix phi_matrix(g_data_matrix.nrow(), g_data_matrix.ncol());
  Rcpp::NumericMatrix phi_matrix_inv(g_data_matrix.nrow(), g_data_matrix.ncol());
  Rcpp::NumericMatrix C_matrix(g_data_matrix.nrow(), g_data_matrix.ncol());
  Rcpp::NumericMatrix C0_matrix(g_data_matrix.nrow(), group0_sample_num);
  Rcpp::NumericMatrix C1_matrix(g_data_matrix.nrow(), group1_sample_num);
  Rcpp::NumericVector C0(g_data_matrix.nrow()), C1(g_data_matrix.nrow());
  Rcpp::NumericMatrix B_matrix(g_data_matrix.nrow(), g_data_matrix.ncol());
  Rcpp::NumericVector B(g_data_matrix.nrow());
  Rcpp::NumericVector phi_before, f_phi_before, ff_phi_before;
  Rcpp::NumericVector phi_after, f_phi_after, ff_phi_after;
  for (int run = 1; run <= times; ++run) {

    // prepare matrix
    for (unsigned int i = 0; i < group0_sample_num; ++i) {
      theta_matrix.column(i) = theta0_old;
      phi_matrix.column(i) = phi_old;
      phi_matrix_inv.column(i) = Rcpp::pow(phi_old, -1);
    }
    for (unsigned int i = group0_sample_num; i < g_data_matrix.ncol(); ++i) {
      theta_matrix.column(i) = theta1_old;
      phi_matrix.column(i) = phi_old;
      phi_matrix_inv.column(i) = Rcpp::pow(phi_old, -1);
    }
    C_matrix = ((g_data_matrix + phi_matrix_inv) * theta_matrix) /
               ((lib_size_norm_matrix * theta_matrix) + phi_matrix_inv);
    B_matrix = digamma_mat(g_data_matrix + phi_matrix_inv) +
               log_mat(theta_matrix / ((lib_size_norm_matrix * theta_matrix) + phi_matrix_inv));
    C0_matrix = C_matrix(Rcpp::_, Rcpp::Range(0, group0_sample_num - 1));
    C1_matrix = C_matrix(Rcpp::_, Rcpp::Range(group0_sample_num, g_data_matrix.ncol() - 1));
    for (unsigned int i = 0; i < g_data_matrix.nrow(); ++i) {
      C0[i] = Rcpp::sum(C0_matrix.row(i));
      C1[i] = Rcpp::sum(C1_matrix.row(i));
      B[i] = Rcpp::sum(B_matrix.row(i));
    }

    // update new theta
    theta0_new = C0 / group0_sample_num;
    theta1_new = C1 / group1_sample_num;

    // update new phi
    phi_before = Rcpp::clone(init_val_phi);
    for (int run = 1; run <= times; ++run) {
      f_phi_before = (group0_sample_num + group1_sample_num) * Rcpp::log(phi_before) -
                     (group0_sample_num + group1_sample_num) +
                     (group0_sample_num + group1_sample_num) * Rcpp::digamma(Rcpp::pow(phi_before, -1)) +
                     group0_sample_num * Rcpp::log(theta0_new) + group1_sample_num * Rcpp::log(theta1_new) -
                     B + Rcpp::pow(theta0_new, -1) * C0 + Rcpp::pow(theta1_new, -1) * C1;
      ff_phi_before = (group0_sample_num + group1_sample_num) * Rcpp::pow(phi_before, -1) -
                      (group0_sample_num + group1_sample_num) * Rcpp::trigamma(Rcpp::pow(phi_before, -1)) *
                      Rcpp::pow(phi_before, -2);
      phi_after = phi_before - f_phi_before / ff_phi_before;
      f_phi_after = (group0_sample_num + group1_sample_num) * Rcpp::log(phi_after) -
                    (group0_sample_num + group1_sample_num) +
                    (group0_sample_num + group1_sample_num) * Rcpp::digamma(Rcpp::pow(phi_after, -1)) +
                    group0_sample_num * Rcpp::log(theta0_new) + group1_sample_num * Rcpp::log(theta1_new) -
                    B + Rcpp::pow(theta0_new, -1) * C0 + Rcpp::pow(theta1_new, -1) * C1;
      double grad_norm = std::sqrt(Rcpp::sum(Rcpp::pow(f_phi_after, 2)));
      phi_before = Rcpp::clone(phi_after);
      if (grad_norm < tol) {
        break;
      }
    }
    phi_new = Rcpp::clone(phi_before);
    double ll_error = std::abs((ll(g_data_matrix, lib_size_norm, group_sample_num, theta0_new, theta1_new, phi_new) -
                               ll(g_data_matrix, lib_size_norm, group_sample_num, theta0_old, theta1_old, phi_old)) /
                               ll(g_data_matrix, lib_size_norm, group_sample_num, theta0_old, theta1_old, phi_old));
    theta0_old = Rcpp::clone(theta0_new);
    theta1_old = Rcpp::clone(theta1_new);
    phi_old = Rcpp::clone(phi_new);
    if (ll_error < tol) {
      break;
    }
  }

  // result
  Rcpp::List res = Rcpp::List::create(Rcpp::_["theta0"] = theta0_new,
                                      Rcpp::_["theta1"] = theta1_new,
                                      Rcpp::_["phi"] = phi_new);
  return res;
}

// estimate parameter H0---------------------------------------------------------------------------------
Rcpp::List estimate_par_H0(Rcpp::NumericMatrix g_data_matrix,
                           Rcpp::NumericVector lib_size_norm,
                           Rcpp::IntegerVector group_sample_num,
                           Rcpp::NumericVector init_val_theta_H0,
                           Rcpp::NumericVector init_val_phi_H0,
                           double tol,
                           int times) {

  int group0_sample_num = group_sample_num[0];
  int group1_sample_num = group_sample_num[1];
  Rcpp::NumericVector theta_old, theta_new, phi_old, phi_new;
  theta_old = Rcpp::clone(init_val_theta_H0);
  phi_old = Rcpp::clone(init_val_phi_H0);
  theta_new = Rcpp::clone(theta_old);
  phi_new = Rcpp::clone(phi_old);
  Rcpp::NumericMatrix lib_size_norm_matrix(g_data_matrix.nrow(), g_data_matrix.ncol());
  for (unsigned int i = 0; i < g_data_matrix.nrow(); ++i) {
    lib_size_norm_matrix.row(i) = lib_size_norm;
  }

  // estimate
  Rcpp::NumericMatrix theta_matrix(g_data_matrix.nrow(), g_data_matrix.ncol());
  Rcpp::NumericMatrix phi_matrix(g_data_matrix.nrow(), g_data_matrix.ncol());
  Rcpp::NumericMatrix phi_matrix_inv(g_data_matrix.nrow(), g_data_matrix.ncol());
  Rcpp::NumericMatrix C_matrix(g_data_matrix.nrow(), g_data_matrix.ncol());
  Rcpp::NumericVector C(g_data_matrix.nrow());
  Rcpp::NumericMatrix B_matrix(g_data_matrix.nrow(), g_data_matrix.ncol());
  Rcpp::NumericVector B(g_data_matrix.nrow());
  Rcpp::NumericVector phi_before, f_phi_before, ff_phi_before;
  Rcpp::NumericVector phi_after, f_phi_after, ff_phi_after;
  for (int run = 1; run <= times; ++run) {
    // prepare matrix
    for (unsigned int i = 0; i < g_data_matrix.ncol(); ++i) {
      theta_matrix.column(i) = theta_old;
      phi_matrix.column(i) = phi_old;
      phi_matrix_inv.column(i) = Rcpp::pow(phi_old, -1);
    }
    C_matrix = ((g_data_matrix + phi_matrix_inv) * theta_matrix) / ((lib_size_norm_matrix * theta_matrix) + phi_matrix_inv);
    B_matrix = digamma_mat(g_data_matrix + phi_matrix_inv) +
               log_mat(theta_matrix / ((lib_size_norm_matrix * theta_matrix) + phi_matrix_inv));
    for (unsigned int i = 0; i < g_data_matrix.nrow(); ++i) {
      C[i] = Rcpp::sum(C_matrix.row(i));
      B[i] = Rcpp::sum(B_matrix.row(i));
    }

    // update new theta
    theta_new = C / (group0_sample_num + group1_sample_num);

    // update new phi
    phi_before = Rcpp::clone(init_val_phi_H0);
    for (int run = 1; run <= times; ++run) {
      f_phi_before = (group0_sample_num + group1_sample_num) * Rcpp::log(phi_before) -
                     (group0_sample_num + group1_sample_num) +
                     (group0_sample_num + group1_sample_num) * Rcpp::digamma(Rcpp::pow(phi_before, -1)) +
                     (group0_sample_num + group1_sample_num) * Rcpp::log(theta_new) -
                     B + Rcpp::pow(theta_new, -1) * C;
      ff_phi_before = (group0_sample_num + group1_sample_num) * Rcpp::pow(phi_before, -1) -
                      (group0_sample_num + group1_sample_num) * Rcpp::trigamma(Rcpp::pow(phi_before, -1)) *
                      Rcpp::pow(phi_before, -2);
      phi_after = phi_before - f_phi_before / ff_phi_before;
      f_phi_after = (group0_sample_num + group1_sample_num) * Rcpp::log(phi_after) -
                    (group0_sample_num + group1_sample_num) +
                    (group0_sample_num + group1_sample_num) * Rcpp::digamma(Rcpp::pow(phi_after, -1)) +
                    (group0_sample_num + group1_sample_num) * Rcpp::log(theta_new) -
                    B + Rcpp::pow(theta_new, -1) * C;
      double grad_norm = std::sqrt(Rcpp::sum(Rcpp::pow(f_phi_after, 2)));
      phi_before = Rcpp::clone(phi_after);
      if (grad_norm < tol) {
        break;
      }
    }
    phi_new = Rcpp::clone(phi_before);
    double ll_error = std::abs((ll_H0(g_data_matrix, lib_size_norm, group_sample_num, theta_new, phi_new) -
                               ll_H0(g_data_matrix, lib_size_norm, group_sample_num, theta_old, phi_old)) /
                               ll_H0(g_data_matrix, lib_size_norm, group_sample_num, theta_old, phi_old));
    theta_old = Rcpp::clone(theta_new);
    phi_old = Rcpp::clone(phi_new);
    if (ll_error < tol) {
      break;
    }
  }

  // result
  Rcpp::List res = Rcpp::List::create(Rcpp::_["theta"] = theta_new,
                                      Rcpp::_["phi"] = phi_new);
  return res;
}
/*_____________________________________________________________________________________________________________________
_____________________________________________________________________________________________________________________*/

// create block annotation---------------------------------------
//[[Rcpp::export]]
Rcpp::List create_block_cpp(const Rcpp::DataFrame ann) {

  // output result
  Rcpp::List block_lookup_list;
  Rcpp::List gene_range_lookup_list;
  Rcpp::List gene_size_lookup_list;

  // get annotation file columns
  std::vector<std::string> ann_chr_vec = ann["V1"];
  std::vector<unsigned int> ann_start_vec = ann["V4"];
  std::vector<unsigned int> ann_end_vec = ann["V5"];
  std::vector<std::string> ann_strand_vec = ann["V7"];
  std::vector<std::string> ann_gene_vec = ann["V9"];

  // get unique chr names
  std::vector<std::string> unique_chr_vec(ann_chr_vec);
  std::sort(unique_chr_vec.begin(), unique_chr_vec.end());
  std::vector<std::string>::iterator unique_chr_vec_it;
  unique_chr_vec_it = std::unique(unique_chr_vec.begin(), unique_chr_vec.end());
  unique_chr_vec.erase(unique_chr_vec_it, unique_chr_vec.end());

  // create block for each chr
  for (unique_chr_vec_it = unique_chr_vec.begin(); unique_chr_vec_it != unique_chr_vec.end(); ++unique_chr_vec_it) {

    Rcpp::List g_block_list;
    std::vector<std::string> g_range_gene_vec;
    std::vector<unsigned int> g_range_start_vec, g_range_end_vec;
    // find chr first and last position in annotation file
    std::vector<std::string>::iterator chr_first, chr_last;
    unsigned int chr_first_d, chr_last_d;
    chr_first = std::find_first_of(ann_chr_vec.begin(), ann_chr_vec.end(),
                                   unique_chr_vec_it, std::next(unique_chr_vec_it, 1));
    chr_last = std::find_end(ann_chr_vec.begin(), ann_chr_vec.end(),
                             unique_chr_vec_it, std::next(unique_chr_vec_it, 1));
    chr_first_d = std::distance(ann_chr_vec.begin(), chr_first);
    chr_last_d = std::distance(ann_chr_vec.begin(), chr_last);

    // get unique gene names
    std::vector<std::string> unique_gene_vec(std::next(ann_gene_vec.begin(), chr_first_d),
                                             std::next(ann_gene_vec.begin(), chr_last_d + 1));
    std::sort(unique_gene_vec.begin(), unique_gene_vec.end());
    std::vector<std::string>::iterator unique_gene_vec_it;
    unique_gene_vec_it = std::unique(unique_gene_vec.begin(), unique_gene_vec.end());
    unique_gene_vec.erase(unique_gene_vec_it, unique_gene_vec.end());

    // create block for each gene g
    for (unique_gene_vec_it = unique_gene_vec.begin(); unique_gene_vec_it != unique_gene_vec.end(); ++unique_gene_vec_it) {
      // find gene first and last position in annotation file relative to chr
      std::vector<std::string>::iterator g_first, g_last;
      unsigned int g_first_d, g_last_d;
      g_first = std::find_first_of(std::next(ann_gene_vec.begin(), chr_first_d), std::next(ann_gene_vec.begin(), chr_last_d + 1),
                                   unique_gene_vec_it, std::next(unique_gene_vec_it, 1));
      g_last = std::find_end(std::next(ann_gene_vec.begin(), chr_first_d), std::next(ann_gene_vec.begin(), chr_last_d + 1),
                             unique_gene_vec_it, std::next(unique_gene_vec_it, 1));
      g_first_d = std::distance(std::next(ann_gene_vec.begin(), chr_first_d), g_first);
      g_last_d = std::distance(std::next(ann_gene_vec.begin(), chr_first_d), g_last);

      // get gene g chromosome
      std::string g_chr = *(std::next(ann_chr_vec.begin(), chr_first_d + g_first_d));

      // get gene g strand
      std::string g_strand = *(std::next(ann_strand_vec.begin(), chr_first_d + g_first_d));

      // get gene g unique start and end positions
      std::vector<unsigned int> g_start_vec(std::next(ann_start_vec.begin(), chr_first_d + g_first_d),
                                            std::next(ann_start_vec.begin(), chr_first_d + g_last_d + 1));
      std::vector<unsigned int> g_end_vec(std::next(ann_end_vec.begin(), chr_first_d + g_first_d),
                                          std::next(ann_end_vec.begin(), chr_first_d + g_last_d + 1));
      std::vector<unsigned int> g_start_unique_vec(g_start_vec), g_end_unique_vec(g_end_vec);
      std::vector<unsigned int>::iterator g_start_unique_vec_it, g_end_unique_vec_it;
      std::sort(g_start_unique_vec.begin(), g_start_unique_vec.end());
      g_start_unique_vec_it = std::unique(g_start_unique_vec.begin(), g_start_unique_vec.end());
      g_start_unique_vec.erase(g_start_unique_vec_it, g_start_unique_vec.end());
      std::sort(g_end_unique_vec.begin(), g_end_unique_vec.end());
      g_end_unique_vec_it = std::unique(g_end_unique_vec.begin(), g_end_unique_vec.end());
      g_end_unique_vec.erase(g_end_unique_vec_it, g_end_unique_vec.end());

      // get gene range
      unsigned int g_range_start, g_range_end;
      g_range_start = g_start_unique_vec.front();
      g_range_end = g_end_unique_vec.back();

      // concatenate start and end positions in a multiset
      std::multiset<std::pair<unsigned int, bool> > g_start_end_vec;
      for (g_start_unique_vec_it = g_start_unique_vec.begin(); g_start_unique_vec_it != g_start_unique_vec.end(); ++g_start_unique_vec_it) {
        g_start_end_vec.insert(std::pair<unsigned int, bool>(*g_start_unique_vec_it, 0));
      }
      for (g_end_unique_vec_it = g_end_unique_vec.begin(); g_end_unique_vec_it != g_end_unique_vec.end(); ++ g_end_unique_vec_it) {
        g_start_end_vec.insert(std::pair<unsigned int, bool>(*g_end_unique_vec_it, 1));
      }

      // create block
      unsigned int count = 1;
      unsigned int block_pos_number = 0;
      std::vector<unsigned int> block_start_vec, block_end_vec;
      std::vector<std::string> block_id_vec;
      std::multiset<std::pair<unsigned int, bool> >::iterator g_start_end_vec_it1, g_start_end_vec_it2;
      for (g_start_end_vec_it1 = g_start_end_vec.begin(); g_start_end_vec_it1 != std::prev(g_start_end_vec.end(), 1); ++g_start_end_vec_it1) {
        g_start_end_vec_it2 = std::next(g_start_end_vec_it1, 1);
        // get potential block start and end
        unsigned int block_start, block_end;
        if ((*g_start_end_vec_it1).second == 0) {
          block_start = (*g_start_end_vec_it1).first;
        }
        else if ((*g_start_end_vec_it1).second == 1) {
          block_start = (*g_start_end_vec_it1).first + 1;
        }
        if ((*g_start_end_vec_it2).second == 0) {
          block_end = (*g_start_end_vec_it2).first - 1;
        }
        else if ((*g_start_end_vec_it2).second == 1) {
          block_end = (*g_start_end_vec_it2).first;
        }

        // check if the block start and end is valid
        if (block_start <= block_end) {
          std::vector<unsigned int>::iterator g_start_vec_it, g_end_vec_it;
          bool block_start_end_valid = false;
          for (g_start_vec_it = g_start_vec.begin(); g_start_vec_it != g_start_vec.end(); ++g_start_vec_it) {
            g_end_vec_it = std::next(g_end_vec.begin(), std::distance(g_start_vec.begin(), g_start_vec_it));
            if ((block_start >= (*g_start_vec_it)) && (block_end <= (*g_end_vec_it))) {
              block_start_end_valid = true;
              break;
            }
          }
          if (block_start_end_valid == true) {
            block_start_vec.push_back(block_start);
            block_end_vec.push_back(block_end);
            block_id_vec.push_back(std::to_string(count));
            block_pos_number = block_pos_number + block_end - block_start + 1;
            count = count + 1;
          }
        }
      }

      // modify block id to the same length string
      unsigned int block_id_pad = (block_id_vec.back()).length();
      std::vector<std::string>::iterator block_id_vec_it;
      for(block_id_vec_it = block_id_vec.begin(); block_id_vec_it != block_id_vec.end(); ++block_id_vec_it) {
        (*block_id_vec_it).insert((*block_id_vec_it).begin(), block_id_pad - (*block_id_vec_it).size(), '0');
      }

      // check strand
      if (g_strand == "-") {
        std::reverse(block_id_vec.begin(), block_id_vec.end());
      }

      // save gene result
      Rcpp::DataFrame g_block_df = Rcpp::DataFrame::create(
        Rcpp::_["block_start"] = block_start_vec,
        Rcpp::_["block_end"] = block_end_vec,
        Rcpp::_["block_id"] = block_id_vec,
        Rcpp::_["stringsAsFactors"] = false);
      g_block_list[(*unique_gene_vec_it)] = g_block_df;
      g_range_gene_vec.push_back(*unique_gene_vec_it);
      g_range_start_vec.push_back(g_range_start);
      g_range_end_vec.push_back(g_range_end);
      gene_size_lookup_list[(*unique_gene_vec_it)] = cal_gene_size(block_start_vec, block_end_vec);
    }

    // save chr result
    block_lookup_list[(*unique_chr_vec_it)] = g_block_list;
    Rcpp::DataFrame g_range_df = Rcpp::DataFrame::create(
      Rcpp::_["gene_id"] = g_range_gene_vec,
      Rcpp::_["gene_start"] = g_range_start_vec,
      Rcpp::_["gene_end"] = g_range_end_vec,
      Rcpp::_["stringsAsFactors"] = false);
    gene_range_lookup_list[(*unique_chr_vec_it)] = g_range_df;
  }

  return Rcpp::List::create(
    Rcpp::_["block"] = block_lookup_list,
    Rcpp::_["range"] = gene_range_lookup_list,
    Rcpp::_["size"] = gene_size_lookup_list);
}


// summarize read into read type-----------------------------------
//[[Rcpp::export]]
Rcpp::DataFrame create_read_type_cpp(const std::string& input_sam_path,
                                     const Rcpp::List block_ann,
                                     const Rcpp::List gene_range,
                                       int min_overlap) {

  // output result
  std::vector<std::string> read_id_vec, read_type_vec, read_gene_vec;

  // open file for reading
  std::fstream input_sam;
  input_sam.open(input_sam_path.c_str(), std::fstream::in);

  // create read type for each read
  if (input_sam.is_open()) {

    // used for saving read information for a read
    std::string input_sam_line, read_id, read_chr, read_cigar;
    unsigned int read_start;

    // loop for each read
    while (std::getline(input_sam, input_sam_line)) {

      try {
        // get read information
        std::stringstream input_sam_line_ss(input_sam_line);
        std::getline(input_sam_line_ss, read_id, '\t');
        input_sam_line_ss.ignore(input_sam_line.size(), '\t');
        std::getline(input_sam_line_ss, read_chr, '\t');
        input_sam_line_ss >> read_start;
        input_sam_line_ss.ignore();
        input_sam_line_ss.ignore(input_sam_line.size(), '\t');
        std::getline(input_sam_line_ss, read_cigar, '\t');

        // extract chromosome number
        std::string::size_type found = read_chr.find("chr");
        if (found != std::string::npos) {
          read_chr = read_chr.substr(found + 3);
        }

        // extract cigar flag char and number
        std::vector<char> read_cigar_char_vec;
        std::vector<int> read_cigar_num_vec;
        std::string::size_type found_start = 0;
        found = read_cigar.find_first_of("MIDNSHP");
        while (found != std::string::npos) {
          read_cigar_char_vec.push_back(read_cigar[found]);
          read_cigar_num_vec.push_back(std::stoi(read_cigar.substr(found_start, found - found_start)));
          found_start = found + 1;
          found = read_cigar.find_first_of("MIDNSHP", found_start);
        }

        // create read position sequence
        std::vector<unsigned int> read_seq_start_vec, read_seq_end_vec;
        std::vector<char>::iterator read_cigar_char_vec_it;
        unsigned int read_seq_start_last = read_start;
        unsigned int read_end;
        for(read_cigar_char_vec_it = read_cigar_char_vec.begin(); read_cigar_char_vec_it != read_cigar_char_vec.end();
        ++read_cigar_char_vec_it) {
          if (((*read_cigar_char_vec_it) == 'M') | ((*read_cigar_char_vec_it) == 'D')) {
            read_seq_start_vec.push_back(read_seq_start_last);
            read_seq_end_vec.push_back(read_seq_start_last +
              (*std::next(read_cigar_num_vec.begin(),
                          std::distance(read_cigar_char_vec.begin(), read_cigar_char_vec_it))) - 1);
            read_seq_start_last = read_seq_start_last + (*std::next(read_cigar_num_vec.begin(),
                                                                    std::distance(read_cigar_char_vec.begin(), read_cigar_char_vec_it)));
          }
          else if ((*read_cigar_char_vec_it) == 'N') {
            read_seq_start_last = read_seq_start_last + (*std::next(read_cigar_num_vec.begin(),
                                                                    std::distance(read_cigar_char_vec.begin(), read_cigar_char_vec_it)));
          }
        }
        std::vector<unsigned int> read_seq_vec;
        read_seq_vec = seq_v(read_seq_start_vec, read_seq_end_vec);
        read_end = read_seq_vec.back();

        // get potential gene names
        Rcpp::DataFrame chr_gene_range_df = Rcpp::as<Rcpp::DataFrame>(gene_range[read_chr]);
        Rcpp::CharacterVector chr_gene_range_id_vec = chr_gene_range_df["gene_id"];
        Rcpp::IntegerVector chr_gene_range_start_vec = chr_gene_range_df["gene_start"];
        Rcpp::IntegerVector chr_gene_range_end_vec = chr_gene_range_df["gene_end"];
        std::vector<std::string> potential_gene_id;
        potential_gene_id = Rcpp::as<std::vector<std::string> >(chr_gene_range_id_vec[(chr_gene_range_start_vec <= read_start) &
          (chr_gene_range_end_vec >= read_end)]);

        // check each potential gene
        std::vector<std::string>::iterator potential_gene_id_it;
        for (potential_gene_id_it = potential_gene_id.begin(); potential_gene_id_it != potential_gene_id.end(); ++ potential_gene_id_it) {

          // get block annotation for the gene
          Rcpp::List chr_block_ls = block_ann[read_chr];
          Rcpp::DataFrame chr_gene_block_df = Rcpp::as<Rcpp::DataFrame>(chr_block_ls[(*potential_gene_id_it)]);
          Rcpp::IntegerVector chr_gene_block_start = chr_gene_block_df["block_start"];
          Rcpp::IntegerVector chr_gene_block_end = chr_gene_block_df["block_end"];
          Rcpp::CharacterVector chr_gene_block_id = chr_gene_block_df["block_id"];

          // check each read position to get their block id
          std::vector<unsigned int>::iterator read_seq_vec_it;
          Rcpp::CharacterVector read_block_id_temp;
          std::vector<std::string> read_block_id_vec;
          for (read_seq_vec_it = read_seq_vec.begin(); read_seq_vec_it != read_seq_vec.end(); ++read_seq_vec_it) {
            read_block_id_temp = chr_gene_block_id[(chr_gene_block_start <= (*read_seq_vec_it)) & (chr_gene_block_end >= (*read_seq_vec_it))];
            if(read_block_id_temp.size() == 1) {
              read_block_id_vec.push_back(Rcpp::as<std::string>(read_block_id_temp[0]));
            }
            else if(read_block_id_temp.size() > 1) {
              Rcpp::Rcout << "block overlap error" << std::endl;
            }
          }

          if (read_block_id_vec.size() >= min_overlap) {
            // get unique block id
            std::vector<std::string>::iterator read_block_id_vec_it;
            std::sort(read_block_id_vec.begin(), read_block_id_vec.end());
            read_block_id_vec_it = std::unique(read_block_id_vec.begin(), read_block_id_vec.end());
            read_block_id_vec.erase(read_block_id_vec_it, read_block_id_vec.end());

            // concatenate block id
            std::string read_type = read_block_id_vec[0];
            for (read_block_id_vec_it = read_block_id_vec.begin() + 1; read_block_id_vec_it != read_block_id_vec.end();
            ++read_block_id_vec_it) {
              read_type = read_type + "-" + (*read_block_id_vec_it);
            }

            // save read type
            read_id_vec.push_back(read_id);
            read_gene_vec.push_back(*potential_gene_id_it);
            read_type_vec.push_back(read_type);
          }
        }
      }
      catch(...) {
        Rcpp::Rcout << "summarize failure" << std::endl;
        Rcpp::Rcout << "read_id: " << read_id << " | " << "read_chr: " << read_chr << std::endl;
        Rcpp::Rcout << "==============================================================" << std::endl;
      }
    }
  }
  else if (!input_sam.is_open()) {
    Rcpp::Rcout << "file open failed" << '\n';
  }

  // close input file
  input_sam.close();

  // output
  Rcpp::DataFrame res = Rcpp::DataFrame::create(Rcpp::_["read_id"] = read_id_vec,
                                                Rcpp::_["read_gene"] = read_gene_vec,
                                                Rcpp::_["read_type"] = read_type_vec,
                                                Rcpp::_["stringsAsFactors"] = false);
  return res;
}


// create read type group---------------------------------------------------
// read_gene_vec and read_type_vec must be sorted in increasing order
//[[Rcpp::export]]
Rcpp::DataFrame create_read_type_group_cpp(const std::vector<std::string>& read_gene_unique_vec,
                                           const std::vector<std::string>& read_gene_vec,
                                           const std::vector<std::string>& read_type_vec,
                                           double min_reduce) {
  // containers for result saving
  std::vector<std::string> res_gene_vec, res_read_type_vec;
  arma::uvec res_read_type_group_vec, read_type_group_vec;

  // create read type group for each gene
  std::vector<std::string>::const_iterator read_gene_unique_vec_it;
  for (read_gene_unique_vec_it = read_gene_unique_vec.begin(); read_gene_unique_vec_it != read_gene_unique_vec.end();
  ++ read_gene_unique_vec_it) {
    // find gene range
    unsigned int g_first, g_last;
    std::vector<std::string>::const_iterator g_first_it, g_last_it;
    g_first_it = std::find_first_of(read_gene_vec.begin(), read_gene_vec.end(), read_gene_unique_vec_it,
                                    std::next(read_gene_unique_vec_it, 1));
    g_last_it = std::find_end(read_gene_vec.begin(), read_gene_vec.end(), read_gene_unique_vec_it,
                              std::next(read_gene_unique_vec_it, 1));
    g_first = std::distance(read_gene_vec.begin(), g_first_it);
    g_last = std::distance(read_gene_vec.begin(), g_last_it);

    // create connection graph matrix
    std::vector<std::string> g_read_type_vec(std::next(read_type_vec.begin(), g_first),
                                             std::next(read_type_vec.begin(), g_last + 1));
    arma::imat g_graph_matrix = create_graph_matrix(g_read_type_vec);

    // create read type group
    read_type_group_vec = create_read_type_group_g(g_read_type_vec, g_graph_matrix, min_reduce);

    // save result
    res_gene_vec.insert(res_gene_vec.end(), g_first_it, std::next(g_last_it, 1));
    res_read_type_vec.insert(res_read_type_vec.end(), std::next(read_type_vec.begin(), g_first),
                             std::next(read_type_vec.begin(), g_last + 1));
    res_read_type_group_vec.insert_rows(res_read_type_group_vec.n_rows, read_type_group_vec);
  }

  // return result
  Rcpp::DataFrame res = Rcpp::DataFrame::create(Rcpp::_["read_gene"] = res_gene_vec,
                                                Rcpp::_["read_type"] = res_read_type_vec,
                                                Rcpp::_["read_type_group"] = res_read_type_group_vec,
                                                Rcpp::_["stringsAsFactors"] = false);
  return res;
}


// fit SGNB model and calculate test result----------------------------------------------------------
//[[Rcpp::export]]
Rcpp::DataFrame fit_SGNB_cpp (const std::vector<std::string>& read_gene_unique_vec,
                              const std::vector<std::string>& read_gene_vec,
                              Rcpp::NumericMatrix data_matrix,
                              Rcpp::NumericVector lib_size_norm,
                              Rcpp::IntegerVector group_sample_num,
                              Rcpp::List gene_size_ls,
                              double tol,
                              int times) {

  // containers for saving result
  Rcpp::NumericVector p_value_vec, gene_expr_0, gene_expr_1, lr_ts_vec, df_vec;

  // fit SGNB for each gene
  std::vector<std::string>::const_iterator read_gene_unique_vec_it;
  for (read_gene_unique_vec_it = read_gene_unique_vec.begin(); read_gene_unique_vec_it != read_gene_unique_vec.end();
       ++read_gene_unique_vec_it) {

    try {

      unsigned int g_first, g_last;

      // find gene range
      std::vector<std::string>::const_iterator g_first_it, g_last_it;
      g_first_it = std::find_first_of(read_gene_vec.begin(), read_gene_vec.end(),
                                      read_gene_unique_vec_it, std::next(read_gene_unique_vec_it, 1));
      g_last_it = std::find_end(read_gene_vec.begin(), read_gene_vec.end(),
                                read_gene_unique_vec_it, std::next(read_gene_unique_vec_it, 1));
      g_first = std::distance(read_gene_vec.begin(), g_first_it);
      g_last = std::distance(read_gene_vec.begin(), g_last_it);

      // get gene count data
      Rcpp::NumericMatrix g_data_matrix = data_matrix(Rcpp::Range(g_first, g_last), Rcpp::_);

      // calculate initial value
      // H0 and H1
      Rcpp::List init_val_ls = calculate_initial_val(g_data_matrix, lib_size_norm, group_sample_num);
      Rcpp::NumericVector init_val_theta_0 = init_val_ls["init_val_theta_0"];
      Rcpp::NumericVector init_val_theta_1 = init_val_ls["init_val_theta_1"];
      Rcpp::NumericVector init_val_phi = init_val_ls["init_val_phi"];
      Rcpp::NumericVector init_val_theta_H0 = init_val_ls["init_val_theta_H0"];
      Rcpp::NumericVector init_val_phi_H0 = init_val_ls["init_val_phi_H0"];

      // estimate parameter
      Rcpp::List fitted_val_ls = estimate_par(g_data_matrix, lib_size_norm, group_sample_num, init_val_theta_0,
                                              init_val_theta_1, init_val_phi, tol, times);
      Rcpp::List fitted_val_H0_ls = estimate_par_H0(g_data_matrix, lib_size_norm, group_sample_num, init_val_theta_H0,
                                                    init_val_phi_H0, tol, times);

      // result
      Rcpp::NumericVector theta0, theta1, phi, theta_H0, phi_H0;
      theta0 = fitted_val_ls["theta0"];
      theta1 = fitted_val_ls["theta1"];
      phi = fitted_val_ls["phi"];
      theta_H0 = fitted_val_H0_ls["theta"];
      phi_H0 = fitted_val_H0_ls["phi"];
      Rcpp::NumericVector lr_ts(1);
      lr_ts[0] = -2 * (ll_H0(g_data_matrix, lib_size_norm, group_sample_num, theta_H0, phi_H0) -
        ll(g_data_matrix, lib_size_norm, group_sample_num, theta0, theta1, phi));
      int df = g_data_matrix.nrow();
      Rcpp::NumericVector p_value = Rcpp::pchisq(lr_ts, df, false);
      unsigned int gene_size = gene_size_ls[(*read_gene_unique_vec_it)];
      gene_expr_0.push_back(Rcpp::sum(theta0) / gene_size);
      gene_expr_1.push_back(Rcpp::sum(theta1) / gene_size);
      p_value_vec.push_back(p_value[0]);
      lr_ts_vec.push_back(lr_ts[0]);
      df_vec.push_back(df);
    }
    catch (...) {
      Rcpp::Rcout << "error gene : " << *read_gene_unique_vec_it << '\n';
      p_value_vec.push_back(NA_REAL);
      gene_expr_0.push_back(NA_REAL);
      gene_expr_1.push_back(NA_REAL);
      lr_ts_vec.push_back(NA_REAL);
      df_vec.push_back(NA_REAL);
    }
  }

  Rcpp::DataFrame res = Rcpp::DataFrame::create(Rcpp::_["gene_id"] = read_gene_unique_vec,
                                                Rcpp::_["p_value"] = p_value_vec,
                                                Rcpp::_["lr_ts"] = lr_ts_vec,
                                                Rcpp::_["df"] = df_vec,
                                                Rcpp::_["gene_expr_0"] = gene_expr_0,
                                                Rcpp::_["gene_expr_1"] = gene_expr_1,
                                                Rcpp::_["stringsAsFactors"] = false);
  return res;
}
