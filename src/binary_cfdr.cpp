// force return as a vector rather than single column matrix
#define RCPP_ARMADILLO_RETURN_ANYVEC_AS_VECTOR

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
  arma::vec xc = x(arma::sort_index(x));
  arma::vec yc = y(arma::sort_index(x));

  arma::vec yout(xout.size());

  double x_min = xc(0);
  double x_max = xc(xc.size()-1);

  for(int i = 0; i < xout.size(); ++i) {
    if(xout(i) > x_max) {
      yout(i) = (yc(yc.size()-1) - yc(yc.size()-2))/(xc(xc.size()-1) - xc(xc.size()-2))*(xout(i)-xc(xc.size()-2))+yc(yc.size()-2);
    } else if(xout(i) < x_min) {
      yout(i) = (yc(1) - yc(0))/(xc(1) - xc(0))*(xout(i)-xc(0))+yc(0);
    } else if(xout(i) == x_min) {
      yout(i) = yc(0);
    } else if(xout(i) == x_max) {
      yout(i) = yc(yc.size()-1);
    } else {
      auto upper = std::lower_bound(xc.begin(), xc.end(), xout(i));

      arma::uword x1_idx = std::distance(xc.begin(), upper);

      double x1 = xc(x1_idx);
      double x0 = xc(x1_idx-1);
      double y1 = yc(x1_idx);
      double y0 = yc(x1_idx-1);

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
  arma::vec xc = x(arma::sort_index(x));
  arma::vec yc = y(arma::sort_index(x));

  arma::vec yout(xout.size());

  for(int i = 0; i < xout.size(); ++i) {

    auto upper = std::lower_bound(xc.begin(), xc.end(), xout(i));

    if(upper == xc.begin()) {
      yout(i) = yc(0);
    } else if(upper == xc.end()) {
      yout(i) = yc(yc.size()-1);
    } else {
      arma::uword x1_idx = std::distance(xc.begin(), upper);

      double x1 = xc(x1_idx);
      double x0 = xc(x1_idx-1);
      double y1 = yc(x1_idx);
      double y0 = yc(x1_idx-1);

      yout(i) = y0+((y1-y0)/(x1-x0))*(xout(i)-x0);
    }
  }

  return yout;
}

// [[Rcpp::export]]
arma::vec per_group_binary_cfdr_rcpp(arma::vec p_loo, arma::vec q_loo, arma::vec ps, arma::vec qs, arma::vec x) {
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

    arma::vec p1 = arma::vec(qs.size()).fill(R_NaN);
    arma::vec p0 = arma::vec(qs.size()).fill(R_NaN);

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

    return v;
}

