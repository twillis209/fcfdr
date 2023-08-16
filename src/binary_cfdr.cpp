#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec arma_pmax(arma::vec x, arma::vec y) {
  arma::vec p(x.size());
  
  for(int i = 0; i < x.size(); ++i) {
    p(i) = std::max(x(i), y(i));
  }

  return p;
}

// [[Rcpp::export]]
arma::vec arma_pmin(arma::vec x, arma::vec y) {
  arma::vec p(x.size());
  
  for(int i = 0; i < x.size(); ++i) {
    p(i) = std::min(x(i), y(i));
  }

  return p;
}

// [[Rcpp::export]]
arma::vec approxExtrap_rcpp(arma::vec x, arma::vec y, arma::vec xout) {
  x = x(arma::sort_index(x));
  y = y(arma::sort_index(x));

  arma::vec yout(xout.size());

  double x_min = x(0);
  double x_max = x(x.size()-1);

  for(int i = 0; i < xout.size(); ++i) {
    if(xout(i) > x_max) {
      yout(i) = (y(y.size()-1) - y(y.size()-2))/(x(x.size()-1) - x(x.size()-2))*(xout(i)-x(x.size()-2))+y(y.size()-2);
    } else if(xout(i) < x_min) {
      yout(i) = (y(1) - y(0))/(x(1) - x(0))*(xout(i)-x(0))+y(0);
    } else if(xout(i) == x_min) {
      yout(i) = y(0);
    } else if(xout(i) == x_max) {
      yout(i) = y(y.size()-1);
    } else {
      auto upper = std::lower_bound(x.begin(), x.end(), xout(i));

      arma::uword x1_idx = std::distance(x.begin(), upper);

      double x1 = x(x1_idx);
      double x0 = x(x1_idx-1);
      double y1 = y(x1_idx);
      double y0 = y(x1_idx-1);

      yout(i) = y0+((y1-y0)/(x1-x0))*(xout(i)-x0);
    }
  }

  return yout;
}

//' @title Linear interpolation function
//'
//' @description This function carries out linear interpolation, reproducing a subset of the behaviour of stats::approx in R.
//'
//' @details Linear interpolation is performed with f(x)=f(x0)+[(f(x1)-f(x0))/(x1-x0)]*(x-x0). x is assumed to be sorted; if sorted in descending order, the order of elements in x and y is reversed. If a value of xout lies outside the range [min(x), max(x)], the value is interpolated using the nearest extremum (this corresponds to the behaviour of stats::approx with method=2).
//'
//' @param x x coordinates of points to be interpolated
//' @param y y coordinates of points to be interpolated
//' @param xout x coordinates at which to interpolate
//'
//' @return interpolated values 
//' 
//' @author Tom Willis
// [[Rcpp::export]]
arma::vec approxfun_rcpp(arma::vec x, arma::vec y, arma::vec xout) {
  x = x(arma::sort_index(x));
  y = y(arma::sort_index(x));

  arma::vec yout(xout.size());

  for(int i = 0; i < xout.size(); ++i) {

    auto upper = std::lower_bound(x.begin(), x.end(), xout(i));

    if(upper == x.begin()) {
      yout(i) = y(0);
    } else if(upper == x.end()) {
      yout(i) = y(y.size()-1);
    } else {
      arma::uword x1_idx = std::distance(x.begin(), upper);

      double x1 = x(x1_idx);
      double x0 = x(x1_idx-1);
      double y1 = y(x1_idx);
      double y0 = y(x1_idx-1);

      yout(i) = y0+((y1-y0)/(x1-x0))*(xout(i)-x0);
    }
  }

  return yout;
}

// [[Rcpp::export]]
arma::vec per_group_binary_cfdr_rcpp(arma::vec p_loo, arma::vec q_loo, arma::vec ps, arma::vec qs, arma::vec x, bool verbose = false) {
    double q0 = static_cast<double>(arma::conv_to<arma::uvec>::from(arma::intersect(find(q_loo == 1), find(p_loo > 0.5))).size()) / static_cast<double>(arma::conv_to<arma::uvec>::from(find(p_loo > 0.5)).size());

    double mult = static_cast<double>(arma::conv_to<arma::uvec>::from(arma::intersect(find(q_loo == 0), find(p_loo > 0.5))).size()) / static_cast<double>(arma::conv_to<arma::uvec>::from(arma::intersect(find(q_loo == 1), find(p_loo > 0.5))).size());

    arma::vec q0_sol(ps.size());
    arma::vec q1_sol(ps.size());

    for (int i = 0; i < ps.size(); ++i) {
        double p = ps(i);
        q0_sol(i) = std::max((int) arma::conv_to<arma::uvec>::from(arma::intersect(find(p_loo <= p), find(q_loo == 0))).size(), 1);
        q1_sol(i) = std::max((int) arma::conv_to<arma::uvec>::from(arma::intersect(find(p_loo <= p), find(q_loo == 1))).size(), 1);
    }

    arma::vec sol(qs.size());

    for (int i = 0; i < qs.size(); ++i) {
        if (qs(i) == 0) {
            sol(i) = mult * ps(i) / q0_sol(i);
        } else {
            sol(i) = (1. / mult) * ps(i) / q1_sol(i);
        }
    }

    arma::vec y(x.size());

    for(int i = 0; i < x.size(); ++i) {
      y(i) = (double) x(i) / (double) std::max((int) arma::conv_to<arma::uvec>::from(arma::intersect(arma::find(p_loo <= x(i)), arma::find(q_loo == 0))).size(), 1);
    }

    arma::vec extr_x = sol.elem(arma::find_unique(sol));
    arma::vec extr_y = approxExtrap_rcpp(y, x, extr_x);

    arma::vec invg0_y = arma_pmax(arma_pmin(extr_y, arma::vec(extr_y.size()).fill(1.)), arma::vec(extr_y.size()).fill(0.));

    arma::vec y1(x.size());

    for(int i = 0; i < x.size(); ++i) {
      y1(i) = (double) x(i) / (double) std::max((int) arma::conv_to<arma::uvec>::from(arma::intersect(arma::find(p_loo <= x(i)), arma::find(q_loo == 1))).size(), 1);
    }

    arma::vec extr1_y = approxExtrap_rcpp(y1, x, extr_x);
    arma::vec invg1_y = arma_pmax(arma_pmin(extr1_y, arma::vec(arma::size(extr1_y)).fill(1.)), arma::vec(arma::size(extr1_y)).fill(0.));

    arma::vec p1 = arma::zeros<arma::vec>(x.size());
    arma::vec p0 = arma::zeros<arma::vec>(x.size());

    for (int i = 0; i < qs.size(); ++i) {
        if (qs(i) == 0) {
          p1(i) = (double) approxfun_rcpp(extr_x, invg1_y, arma::vec(1).fill(sol(i)))(0);
        } else {
            p1(i) = ps(i);
        }

        if (qs(i) == 1) {
          p0(i) = (double) approxfun_rcpp(extr_x, invg0_y, arma::vec(1).fill(sol(i)))(0);
        } else {
            p0(i) = ps(i);
        }
    }

    arma::vec v = p0 * (1 - q0) + p1 * q0;

    if(verbose) {
      Rcout << "q0: " << q0 << "\n";
      Rcout << "mult: " << mult << "\n";
      Rcout << "q0_sol:" << "\n";
      for(int i = 0; i < q0_sol.size(); ++i) {
        Rcout << q0_sol(i) << "\n";
      }
      
      Rcout << "q1_sol:" << "\n";
      for(int i = 0; i < q1_sol.size(); ++i) {
        Rcout << q1_sol(i) << "\n";
      }
       
      Rcout << "sol:" << "\n";
      for(int i = 0; i < sol.size(); ++i) {
        Rcout << sol(i) << "\n";
      }
      
      Rcout << "y:" << "\n";
      for(int i = 0; i < y.size(); ++i) {
        Rcout << y(i) << "\n";
      }
      
      Rcout << "extr_x:" << "\n";
      for(int i = 0; i < extr_x.size(); ++i) {
        Rcout << extr_x(i) << "\n";
      }

      Rcout << "invg0_y:" << "\n";
      for(int i = 0; i < invg0_y.size(); ++i) {
        Rcout << invg0_y(i) << "\n";
      }
      
      Rcout << "y1:" << "\n";
      for(int i = 0; i < y1.size(); ++i) {
        Rcout << y1(i) << "\n";
      }

      // This differs with the R code
      Rcout << "extr1_y:" << "\n";
      for(int i = 0; i < extr1_y.size(); ++i) {
        Rcout << extr1_y(i) << "\n";
      }
      
      Rcout << "invg1_y:" << "\n";
      for(int i = 0; i < invg1_y.size(); ++i) {
        Rcout << invg1_y(i) << "\n";
      }
      
      Rcout << "p1:" << "\n";
      for(int i = 0; i < p1.size(); ++i) {
        Rcout << p1(i) << "\n";
      }
      
      Rcout << "p0:" << "\n";
      for(int i = 0; i < p0.size(); ++i) {
        Rcout << p0(i) << "\n";
      }
      
      Rcout << "v:" << "\n";
      for(int i = 0; i < v.size(); ++i) {
        Rcout << v(i) << "\n";
      }

    }
    
    return Rcpp::NumericVector(v.begin(), v.end());
}

