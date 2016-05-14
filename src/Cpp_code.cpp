//// Cpp functions 2015-3-11 Thong Pham
#include <Rcpp.h>
#include <iostream>
#include <vector>
#include <string>
#include <functional>
#include <fstream>
#include <cmath>
#include <random>
#include <cstdlib>
#include <stdio.h>
#include <float.h>
#include <omp.h>
using namespace Rcpp;
// [[Rcpp::plugins("cpp11")]]

//// //////////////////////////////////////////////////////////////////////////////////////////////////////
struct edge{
  long from_id;  // id of from_node
  long to_id;    // id of to_node
};

struct node_array {
  std::vector<long>* id_vec;
  long         N;      // true number of distinc ID
  long         max_id; // the maximum ID, equal to id_vec.size() - 1
};

struct degree_array {
  std::vector<std::vector <long> *>  degree;
  std::vector<std::vector <long> *>  out_deg;
  // note that in degree_array, the time start from time 1 to T
};

struct timeline {
  std::vector<std::vector<edge> *> stamp;
  long               T;        // The final time-step in the sampler, note that the timeline start from 0 to T.
  // This mean that the total number of time steps is T + 1 <= M
};

struct PA_function {
  std::vector<double> * A;
};

double log_poisson(long , double );
double my_zeroin(double, double, std::function <double (double)>, double, long);
double log_prime(long);
//double log(double );
std::vector<edge> * input(std::string  );
node_array get_id_node(std::vector<edge> *);
degree_array get_degree(std::vector<edge> *,node_array,timeline);

PA_function read_PA_function(long G, std::string f_namename) ;
std::vector<double> * read_f(long N, std::string f_name);

long max_deg(degree_array, long );
PA_function init_PA_function(long, std::vector<double> &);
double cal_pos(const PA_function &PA, const std::vector<double> *f, const degree_array &D,const  std::vector<double> & Z,
               double lambda1, double lambda2, double s, double scale_f, long T, std::vector<long> &bin,
               std::vector<double> & center_k, std::vector<int> &active);

int swap(const PA_function &PA, const std::vector<double> *f, degree_array &D, const  std::vector<double> & Z,
         timeline &order, const node_array & V, long pos, const double u, double lambda2,std::vector<long > & bin);

int even_swap(const PA_function &, const std::vector<double> * , degree_array &, const  std::vector<double> &,
              timeline &, const node_array &, double &,const double, const std::vector<double> &, long, double, std::vector<long> &bin);

int odd_swap(const PA_function &, const std::vector<double> *, degree_array &, const  std::vector<double> &,
             timeline &, const node_array &, double &,const double, const std::vector<double> &, long, double, std::vector<long> &bin);


int update_f(PA_function &PA, std::vector<double> * &f, const degree_array &D, const std::vector <double> &Z,
             const double &s,  std::default_random_engine &gen,
             long T,std::vector<long> &bin,std::vector<double> &center_k, std::vector<int> &active);

int update_Z(const PA_function &PA, const std::vector<double> * f, const degree_array &D, std::vector <double> &Z,
             std::default_random_engine &gen, long T, std::vector<long> &bin, std::vector<int> &active) ;

int update_s(const std::vector<double> * , const PA_function &PA, std::vector<double> &center_k,
             double &, double, double, double, double ,double, std::vector<int> &active);

int para_gamma(double mean, double var, double &shape, double &rate);
int random_scaling_f(std::vector<double> * &f, const double &s, double & scale_f, std::default_random_engine &gen, std::vector<int> &active);
int print_deg(degree_array,long T) ;
double log_gamma(double );
double log_gamma(long );
double log_gamma_density(double , double , double );
int binning(long deg_max, long &G, std::vector<long> &bin_vector,std::vector<double> &center_k);
//// //////////////////////////////////////////////////////////////////////////////////////////////////////


// Important note about implementation:
// timeline: useable index from 0 to T
// degree_array: useable index from 1 to T + 1
// degree_array(T+1) will be the final degree
// to calculate degree array: loop on timeline from 0 to T, at each time t fill in the degree vector at t + 1
// to calculate posterior: loop on degree_array from 1 to T, at each time t compare degree t + 1 and t
// Z: index from 1 to T (length T+1)
// T should be maximum M - 1 (M: number of edges)
// the swap will go from index 0 to T

// [[Rcpp::export(".mcmc_running")]]
int mcmc_running(std::vector< std::string > argv) {
  //get the matrix of edges

  //std::cout << " wtf ";
  std::vector<edge> * E = input(argv[0]);
  node_array V = get_id_node(E);
  //get the list of nodes
  std::random_device rd;
  std::default_random_engine gen(3);
  std::uniform_real_distribution<double> uni(0.0,1.0);
  //std::uniform_int_distribution<int> uni_int(0,normal.order - 1);
  std::normal_distribution<double> normal_s(0.0,std::stod(argv[5]));      // proposal distribution for s


  //for (long i = 0; i < E->size(); ++i)
  //   std::cout<<E->at(i).from_id<<" "<<E->at(i).to_id <<std::endl;

  timeline normal_order = {
    std::vector<std::vector<edge>*>(E->size()),
    E->size() - 1
  };
  long T = E->size() - 1;
  printf("%d \n",T);
  //fflush(stdout);
  for (long i = 0; i <= T; ++ i) {
    normal_order.stamp.at(i) =  new std::vector<edge>;
    normal_order.stamp.at(i)->reserve(E->size());
    normal_order.stamp.at(i)->push_back(E->at(i));
  }
  //std::cout << normal_order.stamp.at(1)->size();
  //for (long i = 0; i < E->size(); ++i)
  //    std::cout<<normal_order.stamp.at(i)->at(0).from_id<<" " << normal_order.stamp.at(i)->at(0).to_id<<std::endl;
  //get the degree array
  degree_array D = get_degree(E,V,normal_order);

  // auxilary variable Z
  std::vector<double> Z(T + 1,2);


  double lambda1 = std::stod(argv[6]);
  double lambda2 = std::stod(argv[6]);
  double s       = std::stod(argv[8]); // starting value for prior of f


  long   G       = std::stol(argv[7]);
  long   burn_in = std::stol(argv[1]); //number of burn-in;
  long   M       = std::stol(argv[2]); //number of needed sample;
  long   skip    = std::stol(argv[3]); //number of skip;
  double swap_accept_percentage = 0;

  double h_s_shape = std::stod(argv[9]);
  double h_s_rate  = std::stod(argv[10]);

  long   s_accept     = 0;
  long   B            = std::stol(argv[4]); //number of combination of even and odd sweep inside an iteration
  // so the total iteration is burn_in + M*skip
  // log posterior

  long deg_max = max_deg(D,T);
  if (G >= deg_max + 1)
      G = deg_max + 1;  
  std::vector<double>  log_pos(0,burn_in + M*skip + 1);
  std::vector<long>    bin(deg_max + 1,-1);
  std::vector<double>  center_k(G,-1);
  printf("Deg max: %d\n",deg_max);
  printf("s: %f\n",s);
  //return 1;
  binning(deg_max,G,bin,center_k);
  //for (long i = 0; i < bin.size();++i)
  //    std::cout<<bin.at(i)<<" ";
  //std::cout<<"size bin:"<<bin.size()<<"\n";
  //for (long i = 0; i < center_k.size();++i)
  //   std::cout<<center_k.at(i)<<" ";
  //printf("reach here");
  long not_use_A = std::stol(argv[11]);
  long not_use_f = std::stol(argv[13]);
  //printf("Cant reach here");
  double scale_f = 1;
  // the PA function
  PA_function A;

  if (not_use_A == 1)
    A  = init_PA_function(G,center_k);
  else A = read_PA_function(G,argv[12]);
  //should find the active set: set of nodes with in-degree >= 1: we only calculate posterior and log-likelihood inside this set
  std::vector<int> active;
  active.reserve(V.N);
  for (int j = 0; j < V.N; ++j)
    if (D.degree.at(T+1)->at(j) > 0)
      active.push_back(j);

    printf("Size of active: %d\n",active.size());

    // fitness
    std::vector<double> * f ;
    if (not_use_f == 1)
      f = new std::vector<double>(V.N,1);
    else f = read_f(V.N, argv[14]);
   
    //print_deg(D,T);

    //std::cout<<"\n";
    //for (long i = 0; i < A.A->size();++i)
    //   std::cout<<A.A->at(i)<<" ";
    //"file_no_" + std::to_string(k) + ".vtk"
    std::string string_file = argv[15];

    std::ofstream f_out;
    f_out.open("f_out_" + string_file +  ".binary",std::ios::out | std::ios::binary);

    std::ofstream Z_out;
    Z_out.open("z_out_" + string_file +  ".binary", std::ios::out | std::ios::binary);

    std::ofstream order_out;
    order_out.open("order_out_" + string_file + ".binary",std::ios::out | std::ios::binary);

    std::ofstream s_out;
    s_out.open("s_out_" + string_file + ".binary",std::ios::out | std::ios::binary);

    std::ofstream log_pos_out;
    log_pos_out.open("log_pos_" + string_file + ".binary",std::ios::out | std::ios::binary);

    log_pos.push_back(cal_pos(A,f,D,Z,lambda1,lambda2,s, scale_f, T,bin,center_k,active));
    log_pos_out.write((char *) &log_pos.back(),sizeof(double));


    for (long i = 1; i <= burn_in + M*skip; i++) {
      if (i % 1000 == 0)
        printf("%d ",i);

      for (long j = 0; j < B; ++j) {
        std::vector<double> u_array_odd((long)ceil((T - 1.0)/2.0),0);
        for (long k = 0; k < u_array_odd.size(); k++)
          u_array_odd.at(k) = uni(gen);
        odd_swap(A,f,D,Z,normal_order,V, swap_accept_percentage,(i - 1) * B + j,u_array_odd, T,lambda2,bin);
        std::vector<double> u_array_even((long)floor((T - 1.0)/2.0),0);
        for (long k = 0; k < u_array_even.size(); k++)
            u_array_even.at(k) = uni(gen);
        even_swap(A,f,D,Z,normal_order, V ,swap_accept_percentage,(i - 1) * B + j,u_array_even, T, lambda2,bin);
      }

      //std::cout<<i<<" ";
      if (i % 1 == 0) {
        //  Update f
        update_f(A,f,D,Z,s,gen,T,bin,center_k, active);
        // random scaling for f
        random_scaling_f(f, s,scale_f, gen, active);
        update_Z(A,f,D,Z,gen,T,bin,active);
      }
      double r_u, r_normal;
      for (long kk = 0; kk < 1; ++ kk) {

      // Update s: hyper parameter of f
          r_u       = uni(gen);
          r_normal  = normal_s(gen);
          s_accept += update_s(f,A, center_k, s, scale_f, r_u,r_normal,h_s_shape,h_s_rate, active);
      }


      log_pos.push_back(cal_pos(A,f,D,Z,lambda1,lambda2,s, scale_f, T,bin,center_k,active));
      //fwrite(&log_pos.back(),1, sizeof(double),log_pos_out);

      if ((i > burn_in) && (i % skip == 0)) {
        f_out.write( (char *) &f->at(0), f->size()*sizeof(double));
        //Z_out.write( (char *) &Z[0], Z.size()*sizeof(double));
        for (long tt = 0; tt <= normal_order.T; ++tt) {
          order_out.write( (char *) &normal_order.stamp.at(tt)->at(0), sizeof(struct edge));
        }
        s_out.write( (char *) &s, sizeof(double));
        //std::cout << "finish write\n";
      }
    }
    log_pos_out.write( (char *) &log_pos.at(0), log_pos.size()*sizeof(double));

    //std::cout << "Done main loop\n";

    //std::cout << "\n swap accept percentage: " << swap_accept_percentage << "\n";


    //std::cout << "s accept:" << s_accept << " \n";

    //std::cout<<"size of A:" << A.A->size()<<"\n";
    //std::cout<<"size of f:" << f->size() << "\n";
    //std::cout<<"number of edges:"<< normal_order.T<< "\n";
    //std::cout<<"number of samples:" << M <<"\n";

    log_pos_out.close();

    f_out.close();
    Z_out.close();
    order_out.close();
    s_out.close();

    delete E;
    delete V.id_vec;
    for (long t = 0; t < normal_order.stamp.size(); t++) {
      delete normal_order.stamp.at(t);
    }
    for (long t = 1; t <= normal_order.T + 1; t++) {
      delete D.degree.at(t);
      delete D.out_deg.at(t);
    }

    delete A.A;
    delete f;
    //std::cout << "here";
    return 0;
}

std::vector<edge> * input(std::string  f_name) {
  std::ifstream f_in;
  f_in.open(f_name);
  // printf("%s \n",f_name);
  //fflush(stdout);
  long a;
  long b;
  std::vector<edge> * result = new std::vector<edge>;
  while (f_in >> a) {
    f_in >> b;
    edge temp = {a,b};
    //    printf("%d %d\n",a,b);
    //fflush(stdout);
    result->push_back(temp);
  }
  f_in.close();
  return result;
}

node_array get_id_node(std::vector<edge> *E){
  long max = 0; // the maximum id in the network
  long temp;
  //  printf("%d \n",E->size());
  //fflush(stdout);
  for (unsigned long i = 0; i < E->size(); ++i) {
    if ((temp = E->at(i).from_id) > max)
      max = temp;
    if ((temp = E->at(i).to_id) > max)
      max = temp;
  }
  node_array result;
  result.id_vec = new std::vector<long>(max + 1,-1);
  long count = -1;

  // the new id will start from 0
  for (unsigned long i = 0; i < E->size(); ++i) {
    if ((result.id_vec->at(E->at(i).from_id) == -1)) {
      count++;
      result.id_vec->at(E->at(i).from_id) = count;
    }
    if ((result.id_vec->at(E->at(i).to_id) == -1)) {
      count++;
      result.id_vec->at(E->at(i).to_id) = count;
    }
  }
  //the true number of distinc IDs
  result.N      = count + 1;
  result.max_id = max;
  //printf("%d %d\n",result.N,result.max_id);
  //fflush(stdout);
  return result;
}

degree_array get_degree(std::vector<edge> * E, node_array V,timeline time_order) {
  long T = time_order.T; //number of current time step in the time_order
  long N = V.N;          // number of nodes
  unsigned long M = E->size();    //number of edges
  //  printf("%d %d\n",T,N);
  //fflush(stdout);
  degree_array D = {
    std::vector<std::vector<long>*>(T + 2), // index start from 1
    std::vector<std::vector<long>*>(T + 2) // index start from 1
  };
  for (long t = 1; t <= T + 1; t++) {
    D.degree.at(t)  =  new std::vector<long>(N,-1);
    D.out_deg.at(t) =  new std::vector<long>(N,-1);
  }
  for (long t = 0; t <= T; ++t) {
    unsigned long j = 0;
    if (t > 0)
      for (unsigned long i = 0; i <  D.degree.at(t)->size(); ++i) {
        D.degree.at(t+1)->at(i)  = D.degree.at(t)->at(i);
        D.out_deg.at(t+1)->at(i) = D.out_deg.at(t)->at(i);
      }
      while ((j < M) && (j < time_order.stamp.at(t)->size())) {
        //std::cout<<"must be here\n";
        long to_node   = V.id_vec->at(time_order.stamp.at(t)->at(j).to_id);
        //std::cout<<"must be here!!\n";
        //std::cout << to_node << " ";
        long from_node = V.id_vec->at(time_order.stamp.at(t)->at(j).from_id);
        //std::cout << from_node << " \n";
        if (D.out_deg.at(t + 1)->at(from_node) == -1)
          D.out_deg.at(t + 1)->at(from_node) = 0;
        if (D.out_deg.at(t + 1)->at(to_node) == -1)
          D.out_deg.at(t + 1)->at(to_node) = 0;
        if (D.degree.at(t + 1)->at(to_node) == -1)
          D.degree.at(t + 1)->at(to_node) = 0;
        if (D.degree.at(t + 1)->at(from_node) == -1)
          D.degree.at(t + 1)->at(from_node) = 0;
        D.degree.at(t + 1)->at(to_node)++;
        D.out_deg.at(t + 1)->at(from_node)++;
        ++j;
      }
  }
  return D;
}

long max_deg(degree_array D, long T) {
  //   printf("%d\n",D.degree.size());
  // fflush(stdout);
  long N      = D.degree.at(1)->size() - 1;
  //    printf("shit\n");
  //fflush(stdout);
  long result = 0;
  long temp;
  for (long i = 0; i < N; ++i)
    if ((temp = D.degree.at(T+1)->at(i)) > result)
      result = temp;
    return result;
}

PA_function init_PA_function(long G, std::vector<double> & center_k) {
  //std::cout<<"how how";
  PA_function result;
  G = center_k.size();
  result.A = new std::vector<double>(G,1);
  //std::cout<<"why not here";
  for (unsigned long i = 0; i < result.A->size(); ++i)
    result.A->at(i) = center_k.at(i);
  return result;
}

PA_function read_PA_function(long G, std::string f_name) {
  PA_function result;
  result.A = new std::vector<double>(G,1);
  std::ifstream f_in;
  f_in.open(f_name);
  double a;
  long   i = 0;
  while (f_in >> a) {
    if (a != 0)
      result.A->at(i) = a;
    i++;
  }
  f_in.close();
  return result;
}

std::vector<double> * read_f(long N, std::string f_name) {
  std::vector<double> * result = new std::vector<double>(N,1);
  std::ifstream f_in;
  f_in.open(f_name);
  double a;
  long   i = 0;
  while (f_in >> a) {
    if (a != 0)
      result->at(i) = a;
    i++;
  }
  f_in.close();
  return result;
}

//calculate log posterior
double cal_pos(const PA_function &PA, const std::vector<double> *f, const degree_array &D,const  std::vector<double> & Z,
               double lambda1, double lambda2, double s, double scale_f, long T, std::vector<long> &bin,
               std::vector<double> &center_k, std::vector<int> &active){
  long N = f->size();   // number of nodes. remember the index start from 0;
  double result = 0;
  #pragma omp parallel for \
  default(shared)        \
    reduction (+:result)
    for (long t = 1; t <= T; ++t) {
      double norm = 0;
      long   m_t = 0; // total new edges at time t
      long   n_t = 0; // total new nodes at time t
      for (long i = 0; i < N; ++i)
          if (f->at(i) > 0) {  
              if (D.degree.at(t)->at(i) != -1){
              long m = D.degree.at(t+1)->at(i) - D.degree.at(t)->at(i);
              long temp = D.degree.at(t)->at(i);
              //long temp2 = bin.at(temp);
              result += m * (log(PA.A->at(bin[temp])) + log(f->at(i)));
              result += -log_prime(m);
              m_t    += m;
              norm   += PA.A->at(bin[temp]) * f->at(i);
              }
            else if (D.degree.at(t+1)->at(i) != - 1) {
               ++n_t;
               long m = D.degree.at(t+1)->at(i);
               result += m * (log(PA.A->at(0)) + log(f->at(i)));
              result += -log_prime(m);
              m_t    += m;
              norm   += PA.A->at(0) * f->at(i);
            }
          }
        //prior of m_t and n_t
        result += m_t * log(lambda1) - log_prime(m_t) - lambda1;
        result += n_t * log(lambda2) - log_prime(n_t) - lambda2;
        result += -log_prime(m_t) - log_gamma(m_t) - norm * Z[t];
        if ((m_t > 0) && (Z[t] >= 0))
          result += (m_t - 1) * log(Z[t]);
    }
    //log probability of m_0, n_0
    long m_0 = 0;
  long n_0 = 0;
  for (long i = 0; i < N; ++i)
    if (D.degree.at(1)->at(i) != -1) {
      ++n_0;
      m_0 += D.degree.at(1)->at(i);
    }
    result += m_0 * log(lambda1) - log_prime(m_0) - lambda1;
    result += n_0 * log(lambda2) - log_prime(n_0) - lambda2;
    // log probability of G_0
    // 1. nodes are distinguiable, edges are indistinguiable
    // 2. allows dupplications of nodes
    // so the probability is (1/n_0)^m_0
    result -= m_0*log(n_0);

    //prior of f: gamma prior
     #pragma omp parallel for \
      default(shared)      \
      reduction (+:result)
      for (long i = 0; i < N; ++ i){
          double shape = s;
          double rate  = s;
          result += shape * log(rate) - log_gamma(shape) + (shape - 1) * log(f->at(i))- rate * f->at(i);
      }
      //prior of b and s
      //result += 2 * log(s);
      return result;
}

//perform one adjacency switch by processing the timeline
int swap(const PA_function &PA, const std::vector<double> *f, degree_array &D, const  std::vector<double> & Z,
         timeline &order, const node_array & V, long pos, const double u, double lambda2, std::vector<long> & bin) {
  double  delta_log   = 0;   // difference of log posterior between after and before the swap
  long    N = f->size();
  long    n_1_b4      = 0; // number of new nodes at time t before the swap
  long    n_2_b4      = 0; // number of new nodes at time t + 1 before the swap
  long    n_1_after   = 0; // number of new nodes at time t after the swap
  long    n_2_after   = 0; // number of new nodes at time t + 1 after the swap
  long    m_t_1_b4    = 0; // total edges before swap between time t and t + 1
  long    m_t_1_after = 0; // total edges after swap between time t+1 and t+2
  long    m_t_2_b4    = 0; // total edges before swap between time t+1 and t+2
  long    m_t_2_after = 0; // total edges after swap between time t+1 and t+2
  //std::cout<<"go in swap ";
  std::vector<long> * deg_in_temp;
  if (pos != 0)
    deg_in_temp = new std::vector<long>(D.degree.at(pos)->begin(),D.degree.at(pos)->end());
  else
    deg_in_temp = new std::vector<long> (N,-1);

  std::vector<int>  is_appeared_pos1(N,0);
  std::vector<int>  is_appeared_pos(N,0);
  std::vector<int>  subtract_pos1(N,0);
  std::vector<int>  subtract_pos(N,0);

  for (unsigned long m = 0; m < order.stamp.at(pos + 1)->size(); ++ m) {
    long id_from = V.id_vec->at(order.stamp.at(pos + 1)->at(m).from_id);
    long id_to   = V.id_vec->at(order.stamp.at(pos + 1)->at(m).to_id);

    if (D.degree.at(pos + 1)->at(id_to) != - 1) {
      delta_log -=  (log(f->at(id_to)) + log(PA.A->at(bin[D.degree.at(pos + 1)->at(id_to)])));
      // delta_log +=  Z.at(pos + 1)*f->at(id_to)*PA.A->at(bin[D.degree.at(pos + 1)->at(id_to)]);

    }
    if ((pos != 0) && (D.degree.at(pos)->at(id_to) != -1)){
      delta_log +=  log(f->at(id_to))+ log(PA.A->at(bin[D.degree.at(pos)->at(id_to)]));
      //delta_log -=  (Z.at(pos)*f->at(id_to)*PA.A->at(bin[D.degree.at(pos)->at(id_to)]));

    }

    if (deg_in_temp->at(id_to) != - 1)
      deg_in_temp->at(id_to)++;
    else deg_in_temp->at(id_to) = 1;
    if (deg_in_temp->at(id_from) == - 1)
      deg_in_temp->at(id_from) = 0;
    m_t_1_after++; m_t_2_b4++;

    if (is_appeared_pos1.at(id_to) == 0) {
      if (pos != 0) {
        if (D.degree.at(pos)->at(id_to) == - 1)
          n_1_after++;
      }
      else n_1_after++;
      if (D.degree.at(pos + 1)->at(id_to) != - 1) {
        delta_log += Z.at(pos+1)*f->at(id_to)*PA.A->at(bin[D.degree.at(pos + 1)->at(id_to)]);
        if ((pos != 0) && (D.degree.at(pos)->at(id_to) != -1))
          delta_log -= Z.at(pos+1)*f->at(id_to)*PA.A->at(bin[D.degree.at(pos)->at(id_to)]);
      }
      if (D.degree.at(pos + 1)->at(id_to) == - 1)
        n_2_b4++;
      is_appeared_pos1.at(id_to) = 1;
    }
    if (is_appeared_pos1.at(id_from) == 0) {
      if (pos != 0) {
        if (D.degree.at(pos)->at(id_from) == - 1)
          n_1_after++;
      } else n_1_after++;
      if (D.degree.at(pos + 1)->at(id_from) == - 1)
        n_2_b4++;
      is_appeared_pos1.at(id_from) = 1;
    }

  }
  for (unsigned long m = 0; m < order.stamp.at(pos)->size(); ++ m) {
    long id_from = V.id_vec->at(order.stamp.at(pos)->at(m).from_id);
    long id_to   = V.id_vec->at(order.stamp.at(pos)->at(m).to_id);
    if (deg_in_temp->at(id_to) != - 1) {
      delta_log += log(f->at(id_to)) + log(PA.A->at(bin[deg_in_temp->at(id_to)]));

    }
    if ((pos != 0)&&(D.degree.at(pos)->at(id_to) != -1)){
      delta_log -= (log(f->at(id_to)) + log(PA.A->at(bin[D.degree.at(pos)->at(id_to)])));

    }

    m_t_1_b4++;
    m_t_2_after++;
    if (is_appeared_pos.at(id_to) == 0) {
      if (pos != 0) {
        if (D.degree.at(pos)->at(id_to) == - 1)
          n_1_b4++;
      } else n_1_b4++;
      if (deg_in_temp->at(id_to) == - 1)
        n_2_after++;
      if (deg_in_temp->at(id_to) != - 1) {
        delta_log -= Z.at(pos+1)*f->at(id_to)*PA.A->at(bin[deg_in_temp->at(id_to)]);
        if ((pos != 0) && (is_appeared_pos1.at(id_to) == 1)) {
          delta_log -= Z.at(pos+1)*f->at(id_to)*PA.A->at(bin[D.degree.at(pos + 1)->at(id_to)]);

          if (D.degree.at(pos)->at(id_to) != - 1) {
            delta_log += Z.at(pos+1)*f->at(id_to)*PA.A->at(bin[D.degree.at(pos)->at(id_to)]);
          }
        }
        if ((is_appeared_pos1.at(id_to) == 0) && (D.degree.at(pos + 1)->at(id_to) != - 1))
          delta_log += Z.at(pos+1)*f->at(id_to)*PA.A->at(bin[D.degree.at(pos + 1)->at(id_to)]);
      }
      is_appeared_pos.at(id_to) = 1;
    }
    if (is_appeared_pos.at(id_from) == 0) {
      if (pos != 0) {
        if (D.degree.at(pos)->at(id_from) == - 1)
          n_1_b4++;
      } else n_1_b4++;
      if (deg_in_temp->at(id_from) == - 1)
        n_2_b4++;
      is_appeared_pos.at(id_from) = 1;
    }
  }

  // can convert this part to be timeline-based

  //for (unsigned long i = 0; i < N; ++i) {
  //    if (D.degree.at(pos + 1)->at(i) != - 1)
  //        delta_log += Z.at(pos+1)*f->at(i)*PA.A->at(bin[D.degree.at(pos + 1)->at(i)]);
  //    if (deg_in_temp->at(i) != - 1)
  //        delta_log -= Z.at(pos+1)*f->at(i)*PA.A->at(bin[deg_in_temp->at(i)]);
  //}

  //further calculation of the difference of log posterior
  delta_log += log(Z.at(pos + 1))*(m_t_2_after - m_t_2_b4);
  if (pos != 0)
    delta_log += log(Z.at(pos))*(m_t_1_after - m_t_1_b4);
  delta_log += log_poisson(n_1_after,lambda2) + log_poisson(n_2_after,lambda2) -
    log_poisson(n_1_b4,lambda2) - log_poisson(n_2_b4,lambda2);
  // if this is the initial time
  if (pos == 0)
    delta_log -= m_t_1_after * log(n_1_after) + m_t_1_b4 * log(n_1_b4); //the probability of the initial graph
  // MH step
  if (delta_log >= log(u)) { //accept
    delete D.degree.at(pos + 1);
    D.degree.at(pos + 1) = deg_in_temp;
    std::vector<edge> * temp;
    temp = order.stamp.at(pos);
    order.stamp.at(pos)     = order.stamp.at(pos + 1);
    order.stamp.at(pos + 1) = temp;
    return 1;
  }
  else { //reject
    delete deg_in_temp;
    return  0;
  }
}

int odd_swap(const PA_function &PA, const std::vector<double> *f, degree_array &D, const  std::vector<double> & Z,timeline &order, const node_array & V,
             double &percentage,const double num, const std::vector<double> &u_array, long T,double lambda2, std::vector<long> &bin) {

  //std::cout<<"enter odd swap ";
  double temp = 0;

  std::vector<long> position(u_array.size(),0);
  long count = 0;

  for (long i = 0; i + 2 <= T; i = i + 2) {
    position[count] = i;
    count++;
  }
  //std::cout << "ok here";
  //std::cout<<"num of proc"<<omp_get_num_procs();

#pragma omp parallel for \
  default(shared)        \
    reduction(+:temp)
    for (unsigned long i = 0; i < u_array.size(); ++i) {

      temp += swap(PA,f,D,Z,order, V, position[i],u_array.at(i),lambda2,bin);
    }
    percentage = (percentage * num * u_array.size() + temp)/ (num * u_array.size() + u_array.size());
  //delta_log = temp;
  //std::cout<<"go out odd swap";
  return 0;
}

int even_swap(const PA_function &PA, const std::vector<double> * f, degree_array &D, const  std::vector<double> & Z,
              timeline &order, const node_array &V, double &percentage, const double num,
              const std::vector<double> &u_array, long T, double lambda2, std::vector<long> & bin) {

  double temp = 0;
  long count  = 0;

  std::vector<long> position(u_array.size(),0);
  for (long i = 1; i + 2 <= T; i = i + 2) {
    position[count] = i;
    count++;
  }
#pragma omp parallel for \
  default(shared)        \
    reduction(+:temp)
    for (unsigned long i = 0; i < u_array.size(); i++){
      temp += swap(PA,f,D,Z,order, V, position[i],u_array.at(i),lambda2,bin);
      //temp += delta_temp;
    }
    percentage = (percentage * num * u_array.size() + temp)/ (num * u_array.size() + u_array.size());
  //delta_log = temp;
  return 0;
}

//Gibbs update for f
int update_f(PA_function &PA, std::vector<double> * &f, const degree_array &D, const std::vector <double> &Z,
             const double &s, std::default_random_engine &gen, long T,
             std::vector<long> & bin, std::vector<double> &center_k, std::vector<int> &active){

  long N = f->size(); // number of nodes. remember the index start from 0;

  std::uniform_int_distribution<int> uni(0,N - 1);
  for (int i = 0; i < N; ++ i){
    double shape_f_single = s;
    double rate_f_single  = s;
    #pragma omp parallel for \
    default(shared)      \
      reduction (+:shape_f_single,rate_f_single)
      for (long t = 1; t <= T; ++t) {
        if (D.degree.at(t)->at(i) != -1){
          long num_new_edge = D.degree.at(t+1)->at(i) - D.degree.at(t)->at(i);
          shape_f_single += num_new_edge;
          rate_f_single  += PA.A->at(bin[D.degree.at(t)->at(i)]) * Z.at(t);
        }
        else if (D.degree.at(t + 1)->at(i) != -1) {
            long num_new_edge = D.degree.at(t+1)->at(i);
            shape_f_single += num_new_edge;
            rate_f_single  += PA.A->at(0) * Z.at(t);
        }
      }
      std::gamma_distribution<double> gamma_dist(shape_f_single,1/rate_f_single);
      f->at(i) = gamma_dist(gen);
  }
  return 0;
}

// Random rescaling for f
int random_scaling_f(std::vector<double> * &f, const double &s,double & scale_f, std::default_random_engine &gen, std::vector<int> &active) {
  long       N = f->size();
  double shape = s*N;
  double rate  = 1;
  double sum   = 0;
  for (long i = 0; i < N; ++i)
    sum += f->at(i);
  std::gamma_distribution<double> gamma_dist(shape,1/rate);
  scale_f = gamma_dist(gen);
  for (long i = 0; i < N; ++i)
    f->at(i) = scale_f * f->at(i)/sum;
  return 0;
}

//sample the auxilary variable Z
int update_Z(const PA_function &PA, const std::vector<double> * f, const degree_array &D, std::vector <double> &Z,
             std::default_random_engine &gen, long T, std::vector<long> &bin, std::vector<int> &active) {
  long N = f->size();     // number of nodes. remember the index start from 0;

  std::vector<double> shape(T + 1,0);
  std::vector<double> rate(T + 1,0);
#pragma omp parallel for
  for (long t = 1; t <= T; ++t) {
    double norm = 0;
    long   m_t  = 0;           // total edges at the time t
    for (long i = 0; i < N; ++i) {
      if (D.degree.at(t)->at(i) != -1){
        m_t   += D.degree.at(t+1)->at(i) - D.degree.at(t)->at(i);
        norm  += PA.A->at(bin[D.degree.at(t)->at(i)])*f->at(i);
      }

      else if (D.degree.at(t + 1)->at(i) != -1) {
        m_t  += D.degree.at(t+1)->at(i);
        norm += PA.A->at(0)*f->at(i);
      }
    }
    shape.at(t) = m_t;
    rate.at(t)  = norm;
  }
  for (long t = 1; t <= T; ++t)
    if ((shape.at(t) > 0) && (rate.at(t) > 0)){
      //std::cout<<shape[t]<<" "<<rate[t]<<" ";
      std::gamma_distribution<double> gamma_dist(shape.at(t),1/rate.at(t));
      Z.at(t) = gamma_dist(gen);
      //std::cout<<Z[t]<< " \n";
    }
    else if (shape.at(t) == 0) {
      Z.at(t) = 0;
    }
    return 0;
}


int update_s(const std::vector<double> * f, const PA_function &PA, std::vector<double> &center_k,
             double &s, double scale_f, double u, double normal, double h_s_shape, double h_s_rate,
              std::vector<int> &active){
  //do nothing
  double candidate = exp(log(s) + normal);
  double delta_log = 0;

  //prior of f
  for (unsigned long i = 0; i <= active.size() - 1; ++ i) {
    double shape   = s;
    double rate    = s;
    double shape_n = candidate;
    double rate_n  = candidate;

    delta_log += shape_n*log(rate_n)- log_gamma(shape_n) -
      shape * log(rate) + log_gamma(shape) +
      (shape_n - shape)*log(f->at(active.at(i))) -
      (rate_n - rate)* f->at(active.at(i));

  }
  //prior of s
  //delta_log += log_gamma_density(candidate,h_s_shape,h_s_rate) - log_gamma_density(s,h_s_shape,h_s_rate);
  //delta_log += 2 * (log(candidate) - log(s));
  if (log(u) <= delta_log) {
    //accept
    s = candidate;
    return 1;
  }
  else {

    return 0;
  }

  return 0;
}

double log_prime(long n) {
  if (n <0)
    printf("Wrong");
  else if (n==0)
    return 0;
  else {
    double result = 0;
    for (long i = 1; i <= n; ++i)
      result += log(i);
    return result;
  }
  return 0;
}

int para_gamma(double mean, double var, double &shape, double &rate){
  rate  = mean/var;
  shape = mean*rate;
  return 0;
}

int binning(long deg_max, long & G, std::vector<long> & bin, std::vector<double> &center_k) {
  if (G > deg_max + 1) {
    bin.resize(deg_max + 1);
    center_k.resize(deg_max + 1);
    G = deg_max + 1;
    for (long i = 0; i < G; ++ i) {
      bin.at(i) = i;
      center_k.at(i) = i;
    }
    return 0;
  }
  else {
    auto f = [deg_max, G] (double x) {
      double value = deg_max + 1.0;
      for (long i = 0; i < G; ++i)
        value -= floor(pow(x,i));
      return(value);
    };
    double b = my_zeroin(1+pow(10,-10), deg_max + G + 1.1, f, pow(10,-15),1000);
    long count = 0;
    unsigned long pos   = 0;
    for (long i = 0; i < G; ++i) {
      long start = pos;
      long interval_length = (long) floor(pow(b,i));
      for (long j = 0; j < interval_length; ++j) {
        bin.at(pos) = count;
        pos++;
      }
      long end = pos - 1;
      if (end == start)
        center_k.at(i) = start;
      else if (start == 0)
        center_k.at(i)= (start + end) / 2.0;
      else
        center_k.at(i) = start * sqrt((1.0 * end/start));
      count++;
    }
    if (pos == bin.size() - 1)
      bin.at(pos) = bin.at(pos - 1);
    center_k.at(0) = 1;
  }
  return 0;
}

int print_deg(degree_array D, long T) {
  printf("Degree-------\n");
  for (long i = 0; i < T; ++i) {
    for (unsigned long j = 0; j < D.degree.at(0)->size(); ++ j)
      printf("%d ",D.degree.at(i)->at(j));
    printf("     ");
    for (unsigned long j = 0; j < D.degree.at(0)->size(); ++ j)
      printf("%d ",D.out_deg.at(i)->at(j));
    printf("\n");
  }
  return 0;
}

double log_gamma(double x) {
  if (x > 0)
    return lgamma(x);
  else
    return 0;
}

double log_gamma(long x) {
  if (x > 0)
    return lgamma(x);
  else
    return 0;
}

double log_gamma_density(double x, double shape, double rate){
  if (x > 0)
    return shape*log(rate) - lgamma(shape) + (shape - 1)* log(x) - rate * x;
  else
    return 0;
}

double log_poisson(long n, double lambda){
  return(n*log(lambda) - lambda - log_prime(n));
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
    register double t1,cb,t2;
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
    
    if( fabs(new_step) < tol_act ) {	/* Adjust the step to be not less*/
    if( new_step > (double)0 )	/* than tolerance		*/
    new_step = tol_act;
    else
      new_step = -tol_act;
    }
    a = b;  fa = fb;			/* Save the previous approx.	*/
    b += new_step;  fb = f(b);	/* Do step to a new approxim.	*/
    if( (fb > 0 && fc > 0) || (fb < 0 && fc < 0) )
    {                 			/* Adjust c for it to have a sign*/
    c = a;  fc = fa;                  /* opposite to that of b	*/
    }
  }
}



