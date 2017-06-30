#include "mcPAFit_header.h"
// Calculate log posterior by processing the timeline
// Assumption: one edge at a time
// [[Rcpp::export(".estimate_alpha_core")]]
SEXP estimate_alpha_core(SEXP a,
                         SEXP b,
                         SEXP c,
                         SEXP d) {
    //std::cout << "Reached here 0\n";
    NumericVector Sum_m_k(a);
    NumericMatrix n_tk(b);
    NumericVector m_t(c);
    NumericVector center_k(d);
    double first_term = 0;
    //std::cout << center_k.size() << "\n";
    for (unsigned long k = 0; k < center_k.size(); ++k) {
        first_term += Sum_m_k.at(k) * log(center_k.at(k));
    }
    //std::cout << "Reached here 1\n";
    unsigned long K = center_k.size();
    unsigned long T = n_tk.nrow();
    auto f =[&] (double alpha) {
        double second_term = 0;
        #pragma omp parallel for \
        default(shared)        \
        reduction(+:second_term)
        for (unsigned long t = 0; t < T; ++t) {
            double upper = 0;
            double lower = 0;
            for (unsigned long k = 0; k < K; ++k) {
                lower += n_tk.at(t,k) * pow(center_k.at(k),alpha);
                upper += m_t.at(t) * n_tk.at(t,k) * pow(center_k.at(k),alpha) * log(center_k.at(k));
            }
            second_term += upper / lower;
        }
        return(first_term - second_term);
    };
    double alpha = my_zeroin(-2, 2, f, DBL_EPSILON, 500);
    //std::cout << "Reached here 2\n";
    // calculate standard deviation of alpha
    double var   = 0; // variance
    #pragma omp parallel for \
    default(shared)        \
    reduction(+:var)
      for (unsigned long t = 0; t < T; ++t) {
        double u_u = 0;
        double v_v = 0;
        double u_v = 0;
        for (unsigned long k = 0; k < K; ++k) {
            u_u += n_tk.at(t,k) * pow(center_k.at(k),alpha) * pow(log(center_k.at(k)),2);
            u_v += n_tk.at(t,k) * pow(center_k.at(k),alpha) * log(center_k.at(k));
            v_v += n_tk.at(t,k) * pow(center_k.at(k),alpha);
        }
        var += m_t.at(t) * (u_u * v_v - pow(u_v,2)) / pow(v_v,2);
      }
      
    return Rcpp::List::create(Rcpp::Named("alpha")    = alpha,
                              Rcpp::Named("variance") = 1/ var);
}

double my_zeroin(double ax,double bx,std::function <double (double)> f,double tol,long max_iter)    /* An estimate to the root	*/
{
  double a,b,c;				/* Abscissae, descr. see above	*/
  double fa;				/* f(a)				*/
  double fb;				/* f(b)				*/
  double fc;				/* f(c)				*/
  
  a = ax;  b = bx;  fa = f(a);  fb = f(b);
  c = a;   fc = fa;
  long count = 0;
  for(;;)		/* Main iteration loop	*/
  { count++;
    //printf("%ld ",count);
    if (count > max_iter)
      return b;
    double prev_step = b-a;		/* Distance from the last but one*/
    /* to the last approximation	*/
    double tol_act;			/* Actual tolerance		*/
    double p;      			/* Interpolation step is calcu- */
    double q;      			/* lated in the form p/q; divi- */
    /* sion operations is delayed   */
    /* until the last moment	*/
    double new_step;      		/* Step at this iteration       */
    
    if( fabs(fc) < fabs(fb) )
    {                         		/* Swap data for b to be the 	*/
    a = b;  b = c;  c = a;          /* best approximation		*/
    fa=fb;  fb=fc;  fc=fa;
    }
    tol_act = 2*LDBL_EPSILON*fabs(b) + tol/2;
    new_step = (c-b)/2;
    
    if( fabs(new_step) <= tol_act || fb == (double)0 )
      return b;				/* Acceptable approx. is found	*/
    
    /* Decide if the interpolation can be tried	*/
    if( fabs(prev_step) >= tol_act	/* If prev_step was large enough*/
    && fabs(fa) > fabs(fb) )	/* and was in true direction,	*/
    {					/* Interpolatiom may be tried	*/
    double t1,cb,t2;
      cb = c-b;
      if( a==c )			/* If we have only two distinct	*/
      {				/* points linear interpolation 	*/
    t1 = fb/fa;			/* can only be applied		*/
    p = cb*t1;
    q = 1.0 - t1;
      }
      else				/* Quadric inverse interpolation*/
      {
        q = fa/fc;  t1 = fb/fc;  t2 = fb/fa;
        p = t2 * ( cb*q*(q-t1) - (b-a)*(t1-1.0) );
        q = (q-1.0) * (t1-1.0) * (t2-1.0);
      }
      if( p>(double)0 )		/* p was calculated with the op-*/
    q = -q;			/* posite sign; make p positive	*/
    else				/* and assign possible minus to	*/
    p = -p;			/* q				*/
    
    if( p < (0.75*cb*q-fabs(tol_act*q)/2)	/* If b+p/q falls in [b,c]*/
    && p < fabs(prev_step*q/2) )	/* and isn't too large	*/
    new_step = p/q;			/* it is accepted	*/
    /* If p/q is too large then the	*/
    /* bissection procedure can 	*/
    /* reduce [b,c] range to more	*/
    /* extent			*/
    }
    
    if( fabs(new_step) < tol_act )	/* Adjust the step to be not less*/
    {
      if ( new_step > (double)0 ) {	/* than tolerance		*/
    new_step = tol_act;
      } else {
        new_step = -tol_act;
      }
    }
    a = b;  fa = fb;			/* Save the previous approx.	*/
    b += new_step;  fb = f(b);	/* Do step to a new approxim.	*/
    if( (fb > 0 && fc > 0) || (fb < 0 && fc < 0) )
    {                 			/* Adjust c for it to have a sign*/
    c = a;  fc = fa;                  /* opposite to that of b	*/
    }
  }
}




