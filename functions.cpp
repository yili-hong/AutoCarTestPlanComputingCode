#define ARMA_USE_BLAS
#define ARMA_USE_LAPACK

#include <RcppArmadillo.h>
#include <iostream>
#include <Rmath.h>
#include <Rcpp.h>
#include <chrono>

using namespace Rcpp;
using namespace std;
using namespace arma;

//[[Rcpp::depends(RcppArmadillo)]]
 
double psev(double z)
{
  double res;
  res=1-exp(-exp(z));
  return res;
}

 
vec psev_vec(const vec& z)
{
  vec res;
  res=1-exp(-exp(z));
  return res;  
}


vec pow_vec(double a, const vec& b)
{
  int n=b.n_elem;
  vec res=b;
  int i;
  
  for(i=0; i<n; i++)
  {
    res(i)=pow(a, b(i));
  }
  
  return(res);
}



double dsev(double z)
{
  double res;
  res=exp(z-exp(z));
  return(res);
}

vec dsev_vec(const vec& z)
{
  vec res;
  res=exp(z-exp(z));
  return(res);
}

double pnor(double z)
{
  double res;
  res=R::pnorm(z,0.0,1.0,1,0);
  return res;
}

vec pnor_vec(const vec& z)
{
  vec res;
  res=normcdf(z);
  return res;
}

double dnor(double z)
{
  double res;
  res=R::dnorm(z,0.0,1.0,0);
  return res;
}


vec dnor_vec(const vec& z)
{
  vec res;
  res=normpdf(z);
  return res;
}

vec tmpOneFun(const vec& zz)
{
   vec res=zz;
   return(res.ones());
}
              
vec tmpIdFun(const vec& zz)
{
   return(zz);
}

vec tmpLogFun(const vec& zz)
{
   vec res=log(zz);
   return(res);
}
       
vec tmpExpFun(const vec& zz)
{
   vec res=exp(zz);
   return(res);
}

vec sev_dphiphi(const vec& zz)
{
   vec res=1-exp(zz);
   return(res);
}

vec nor_dphiphi(const vec& zz)
{
   vec res=(-1)*zz;
   return(res);
}

class distFuns
{
   public:
       vec (*Phi)(const vec&);   
       vec (*phi)(const vec&);        
       vec (*xtran)(const vec&);
       vec (*ixtran)(const vec&);
       vec (*dent)(const vec&);
       vec (*dphiphi)(const vec&);
       String dtype;
                         
       distFuns(String dd)
       {
         dtype=dd;
         if(dd=="weibull")   //Weibull distribution
         {
           Phi=&psev_vec;
           phi=&dsev_vec;
           xtran=&tmpLogFun;
           ixtran=&tmpExpFun;
           dent=&tmpIdFun;
           dphiphi=&sev_dphiphi;
         }
         
         if(dd=="lognormal")   //lognormal distribution
         {
           Phi=&pnor_vec;
           phi=&dnor_vec;
           xtran=&tmpLogFun;
           ixtran=&tmpExpFun;
           dent=&tmpIdFun;
           dphiphi=&nor_dphiphi;
         }         
         
       } 
        
};


/*  
double test_fun(double tt)
{
   double res=exp((log(tt)-1.0)/2.0);
   return(res);
}
*/


class gibbsModel
{
   public:
       List dat_list;
       double (*log_post_ratio)(double, int, const vec&, const List&);   
};

class gibbsOutput
{
   public:
      mat vmat;
      vec arate;
      vec th;
};

class gibsAdaptOutput
{
   public:
      int adap;
      vec scale;
};

gibsAdaptOutput gibbs_adap_fun(const vec& scale0, const vec& acc, int no_iter)
{
  double tmp, tmp1, dd;
  int i, ind;
  ind=0;
  int p=acc.n_elem;
  
  vec scale=scale0;
  
  tmp=(double) no_iter;
  tmp1=1/tmp;
  
  int adap=0;

  if(tmp1>.02)
  {
    dd=.02;
  }else{
         dd=tmp1;
        }


  for(i=0; i<p;i++)
  {
    if(acc[i]<=.35 && acc[i]>=.15)
    {
     ind++;
    }
  }

  if(ind==p)
  {
    adap=1;
    
  }else{

     for(i=0;i<p;i++)
     {
      if(acc[i]>.35)
      {
        scale[i]+=dd;
      }
      if(acc[i]<.15)
      {
        scale[i]-=dd;
      }
     }
    }
  for(i=0;i<p;i++)
  {
   if(scale[i]<=0)
   {
    scale[i]/=2.0;
   }
  }
  
  gibsAdaptOutput res;
  res.adap=adap;
  res.scale=scale;
  
  return res;
}


/******************************************************************************/
gibbsOutput Gibbs(gibbsModel gibbs_model, const vec& th0, const vec& scale, int mm)
{
  int p, m, i,j;
  double f1, uu, ww, u;
  double yy;
  
  int pp=th0.n_elem;
  mat vmat(mm, pp, fill::zeros);
  vec th=th0;
  vec arate=zeros(pp);
  
  u=0.0;
  p=pp;  //dimension of unknown parameter vector
  m=mm;  //number of mcmc iter


  GetRNGstate();

  for (i=0;i<m; i++)
  {
     for (j=0;j<p;j++)
     {
       ww=norm_rand();//rnorm(0.0,1.0);
       yy=ww * scale[j];
       
       f1=gibbs_model.log_post_ratio(yy, j, th, gibbs_model.dat_list);
       ww=exp(f1);

       uu=unif_rand();//runif(0.0,1.0);

       //Rprintf("%u, %lf, %lf, %lf, %lf,%lf,%lf \n",j, scale(j), th0(j), yy, f1, ww, uu);
       //Rprintf("%u,%lf,%lf,%lf \n",j, ww,f1,f0);

      if(uu<ww)
      {
        u=1.0;
        th[j]+=yy;
      }else{
             u=0.0;
            }

      //Rprintf("%u, %lf, %lf, %lf, %lf, %lf, %lf, %lf \n",j, scale(j), th0(j), yy, f1, ww, uu, u);
            
      vmat(i,j)=th[j];   //mcmc matrix
      arate[j]+=u;
     }
    }
  PutRNGstate();

  for(i=0;i<p;i++)
  {
    arate[i]/=m;
  }
  
  
  gibbsOutput res;
  res.vmat=vmat;
  res.th=th;
  res.arate=arate;

  return res;
}




/******************************************************************************/
List Gibbs_adaptive(gibbsModel gibbs_model, const vec& th0, const vec& scale0, int mm, int mm1, int max_iter=400)
{

  int adaptive, no_iter, adap;
  adaptive=1;
  adap=0;
  no_iter=0;
  
  int pp=th0.n_elem;

  vec scale=scale0; 
  vec th=th0; 
  vec arate=zeros(pp);
  
  while((adaptive==1)&& no_iter<max_iter)
  {
    no_iter++;
    gibbsOutput gres1=Gibbs(gibbs_model, th, scale, mm1);
    th=gres1.th;
    arate=gres1.arate;
    
    gibsAdaptOutput ares=gibbs_adap_fun(scale, arate, no_iter);
    adap=ares.adap;
    scale=ares.scale;
    
    Rcout<<"Adaptation iter: "<< no_iter<<", arate range: ("<< arate.min()<<", "<<arate.max()<< "), scale range: ("<<scale.min()<<", "<<scale.max()<<")."<<endl;
    
    if(adap==1)
    {
      adaptive=0;
      
      Rcout<<"Adaptation done after "<< no_iter<<" iterations."<<endl;
    }
  }

  if(no_iter==max_iter)
  {
    Rcout<<"max iter reached: "<< max_iter<<endl;
  }
  
  Rcout<<"Gibbs started. It may take a while."<<endl;
  
  gibbsOutput gres=Gibbs(gibbs_model, th, scale, mm);
  arate=gres.arate;
  mat vmat=gres.vmat;
  
  Rcout<<"Gibbs done."<<endl;

  List res=List::create(Named("vmat")=vmat, Named("arate")=arate, Named("scale")=scale);

  return res;
}

/*NUTS general code------------------------------------------------*/

class nutsModel
{
   public:
       //List dat_list;
       double (*fn_val)(const vec&, const List&);   
       vec (*gn_val)(const vec&, const List&);        
       double (*joint_fn_val)(const vec&, const vec&, const vec&, const List&);
};


vec leapfrog_step_cpp(nutsModel nuts_model, const List& dat_list, const vec& vec_theta, const vec& vec_r, double eps, const vec& M_diag)
{   
    //Rcout<<"leap frog eps="<<eps<<endl;
    
    int int_dim = vec_theta.n_elem;
    vec vec_r_tilde(int_dim);
    vec vec_theta_tilde(int_dim);
    vec ret_vec(2 * int_dim);
    // r_tilde <- r + 0.5 * eps * grad_f(theta=theta, dat=dat)
    vec_r_tilde = vec_r + 0.5 * eps * nuts_model.gn_val(vec_theta, dat_list);
    //vec_r_tilde = vec_r + 0.5 * eps * nuts_model.gn_val(vec_theta);    
    // theta_tilde <- theta + eps * r_tilde/M_diag
    vec_theta_tilde = vec_theta + eps * vec_r_tilde/M_diag;
    // r_tilde <- r_tilde + 0.5 * eps * grad_f(theta=theta_tilde, dat=dat)
    vec_r_tilde = vec_r_tilde + 0.5 * eps * nuts_model.gn_val(vec_theta_tilde, dat_list);
    //vec_r_tilde = vec_r_tilde + 0.5 * eps * nuts_model.gn_val(vec_theta_tilde);
    // res=list(theta = theta_tilde, r = r_tilde)
    ret_vec.subvec(0,int_dim-1) = vec_theta_tilde;
    ret_vec.subvec(int_dim,2*int_dim-1) = vec_r_tilde;
    return(ret_vec);
}



double find_reasonable_epsilon_cpp(nutsModel nuts_model, const List& dat_list, const vec& vec_theta, const vec& vec_M_diag, double num_eps = 1.0)
{
    int int_dim = vec_theta.n_elem;
    vec vec_r(int_dim);
    vec_r = randn(int_dim)%sqrt(vec_M_diag);
    
    vec proposed = leapfrog_step_cpp(nuts_model, dat_list, vec_theta, vec_r, num_eps, vec_M_diag);
    double log_ratio =nuts_model.joint_fn_val(proposed.subvec(0,int_dim-1), proposed.subvec(int_dim,2*int_dim-1), vec_M_diag, dat_list) - nuts_model.joint_fn_val(vec_theta, vec_r, vec_M_diag, dat_list);
    int int_alp = -1;
    if(isnan(log_ratio))
    {
        int_alp = -1;
    } else if(exp(log_ratio) > 0.5)
    {
        int_alp = 1;
    }
    int count = 1;
    while (int_alp * log_ratio > (-int_alp)*log(2) || isnan(log_ratio))
    {
        num_eps = pow(2.0,int_alp) * num_eps;
        proposed = leapfrog_step_cpp(nuts_model, dat_list, vec_theta, vec_r, num_eps, vec_M_diag);
        log_ratio =nuts_model.joint_fn_val(proposed.subvec(0,int_dim-1), proposed.subvec(int_dim,2*int_dim-1), vec_M_diag, dat_list) - nuts_model.joint_fn_val(vec_theta, vec_r, vec_M_diag, dat_list);
        count = count +1;
        if(count > 100)
        {
            Rcout<<"Could not find reasonable epsilon in 100 iterations!"<<endl;
            break;
        }
    }
    Rcout<<"Reasonable epsilon found after "<<count<<" steps.   eps="<<num_eps<<endl;
    return num_eps;
}

class hmctree
{
    public:
    //list(theta_minus=theta, theta_plus=theta, theta=theta, r_minus=r, r_plus=r, s=s, n=n, alpha=alpha, n_alpha=1)
        vec theta_minus;
        vec theta_plus;
        vec theta;
        vec r_minus;
        vec r_plus;
        int s;
        int n;
        double alpha;
        int n_alpha;    
        hmctree(int int_ndim)
        {
            n = int_ndim;
            theta_minus.set_size(n);
            theta_plus.set_size(n);
            theta.set_size(n);
            r_minus.set_size(n);
            r_plus.set_size(n);
            n_alpha = 1;
        }
        hmctree()
        {
            n_alpha = 1;
        }
};


bool check_NUTS_cpp(int int_s, const vec& vec_theta_plus, const vec& vec_theta_minus, const vec& vec_r_plus, const vec& vec_r_minus)
{
    bool condition1 = dot(vec_theta_plus - vec_theta_minus, vec_r_minus) >=0;
    // crossprod(vec_theta_plus - theta_minus, r_minus) >= 0
    bool condition2 = dot(vec_theta_plus - vec_theta_minus, vec_r_plus) >=0;
    // crossprod(vec_theta_plus - theta_minus, r_plus) >= 0
    bool res=(int_s & condition1 & condition2);
    if(isnan(int_s))
    {
        return false;
    }
    return(res);
}

hmctree build_tree_cpp(nutsModel nuts_model, const List& dat_list, vec vec_theta, vec vec_r, double num_u, int int_v, int int_j, double num_eps, vec vec_theta0, vec vec_r0, vec vec_M_diag, int Delta_max = 1000)
{
    int int_dim = vec_theta.n_elem;
    vec proposed(2*int_dim);
    //proposed = (theta,rtilde)
    double log_prob,log_prob0,alpha;
    int int_n, int_s, n_alpha;
    vec vec_theta_minus;
    vec vec_r_minus;
    vec vec_theta_plus;
    vec vec_r_plus;
    hmctree res(int_dim);
    hmctree obj0(int_dim);
    hmctree obj1(int_dim);
    if(int_j == 0)
    {
        proposed = leapfrog_step_cpp(nuts_model, dat_list, vec_theta, vec_r, int_v*num_eps , vec_M_diag);
        vec_theta = proposed.subvec(0,int_dim-1);
        vec_r = proposed.subvec(int_dim, 2*int_dim-1);
        log_prob = nuts_model.joint_fn_val(vec_theta, vec_r, vec_M_diag, dat_list);
        log_prob0 = nuts_model.joint_fn_val(vec_theta0, vec_r0, vec_M_diag, dat_list);
        int_n = (log(num_u) <= log_prob);
        int_s = (log(num_u) < (Delta_max + log_prob));
        alpha = min(1.0, exp(log_prob - log_prob0));
        // if(isnan(alpha))
        // {
        //     Rcout<<"we have a nan alpha!"<<endl;
        //     terminate();
        // }
        if(isnan(int_s))
        {
            Rcout<<"nan s produced"<<endl;
            int_s = 0;
        }
        if(isnan(int_n))
        {
            Rcout<<"nan n produced"<<endl;
            int_n = 0;
        }

        // res=list(theta_minus=theta, theta_plus=theta, theta=theta, r_minus=r, r_plus=r, s=s, n=n, alpha=alpha, n_alpha=1)
        res.theta_minus = vec_theta;
        res.theta_plus = vec_theta;
        res.theta = vec_theta;
        res.r_minus = vec_r;
        res.r_plus = vec_r;
        res.s = int_s;
        res.n = int_n;
        res.alpha = alpha;
        res.n_alpha = 1;
        return(res);
    } else
    {
        obj0 = build_tree_cpp(nuts_model, dat_list, vec_theta, vec_r, num_u, int_v, int_j-1, num_eps, vec_theta0, vec_r0,vec_M_diag);
        vec_theta_minus = obj0.theta_minus;
        vec_r_minus = obj0.r_minus;
        vec_theta_plus = obj0.theta_plus;
        vec_r_plus = obj0.r_plus;
        vec_theta = obj0.theta;
        if(obj0.s ==1)
        {
            if(int_v == -1)
            {
                obj1 = build_tree_cpp(nuts_model, dat_list, obj0.theta_minus, obj0.r_minus, num_u, int_v, int_j-1, num_eps, vec_theta0, vec_r0, vec_M_diag);
                vec_theta_minus = obj1.theta_minus;
                vec_r_minus = obj1.r_minus;
            } else
            {
                obj1 = build_tree_cpp(nuts_model, dat_list, obj0.theta_plus, obj0.r_plus, num_u, int_v, int_j-1, num_eps, vec_theta0, vec_r0,  vec_M_diag);
                vec_theta_plus = obj1.theta_plus;
                vec_r_plus = obj1.r_plus;
            }
            int_n = obj0.n + obj1.n;
            if(int_n != 0)
            {
                double num_prob =(double) obj1.n / (double) int_n;
                // updating theta
                // vec vec_unif(1000);
                // vec_unif.load("vec_unif.txt");
                // Rcout<<"int_j = "<<int_j<<endl;
                // Rcout<<"for debug, using unif.txt, unif = "<<vec_unif(int_j)<<endl;
                if(randu() < num_prob)
                {
                    vec_theta = obj1.theta;
                }
            }
            int_s = check_NUTS_cpp(obj1.s, vec_theta_plus, vec_theta_minus, vec_r_plus, vec_r_minus);
            alpha = obj0.alpha + obj1.alpha;
            n_alpha = obj0.n_alpha + obj1.n_alpha;
        } else
        {
            int_n = obj0.n;
            int_s = obj0.s;
            alpha = obj0.alpha;
            n_alpha = obj0.n_alpha;
        }
        if(isnan(int_s))
        {
            int_s = 0;
        }
        if(isnan(int_n))
        {
            int_n = 0;
        }

        res.theta_minus = vec_theta_minus;
        res.theta_plus = vec_theta_plus;
        res.theta = vec_theta;
        res.r_minus = vec_r_minus;
        res.r_plus = vec_r_plus;
        res.s = int_s;
        res.n = int_n;
        res.alpha = alpha;
        res.n_alpha = n_alpha;
        return(res);
    }
}


// pars = list(eps = eps, eps_bar = eps_bar, H = H, mu = mu, M_adapt = M_adapt, M_diag = M_diag))
class parlist
{
    public:
        double eps;
        double eps_bar;
        double hmcmu;
        double H;
        int M_adapt;
        vec M_diag;
        parlist(int int_dim)
        {
            M_diag.set_size(int_dim);
            M_diag.fill(1.0);
            M_adapt = 50;
            eps = 1.0;
            hmcmu = log(10.0*eps);
            H = 0.0;
            eps_bar = 1.0;
        }
        parlist()
        {
            M_adapt = 50;
            eps = 1.0;
            hmcmu = log(10.0*eps);
            H = 0.0;
            eps_bar = 1.0;
            M_diag.set_size(10);
            M_diag.fill(datum::nan);
        }
};

class thetaparlist
{
    public:
        parlist pl;
        vec theta;
        thetaparlist(int int_dim):
        pl(int_dim)
        {
        }
        thetaparlist():
        pl()
        {
        }
};

void adapt(double &num_H, double &log_eps, double &num_epsbar, double &num_eps, const int& int_iter, const double& num_t0, const double &num_delta, const double &num_gamma, const double &num_kappa, const double &num_hmcmu, const int &int_alpha, const int &n_alpha)
{
        num_H = (1.0 - 1.0/(double) (int_iter + 1 + num_t0))*num_H + 
            1.0/(double) (int_iter + 1 + num_t0) * (num_delta - int_alpha /(double) n_alpha);
        log_eps = num_hmcmu - sqrt((double) (int_iter + 1))/(double) num_gamma * num_H;
        num_epsbar = exp(pow((double) (int_iter+1),-1.0*num_kappa) * log_eps + 
            (1.0 - pow((double) (int_iter+1),(-1.0 * num_kappa))) * log(num_epsbar));
        num_eps = exp(log_eps);
}

thetaparlist NUTS_one_step_cpp(nutsModel nuts_model, const List& dat_list, vec vec_theta, int int_iter, parlist pl,
    double num_delta = 0.5, int max_treedepth = 10, double num_eps = 1)
{
    int int_dim = vec_theta.n_elem;
    double num_kappa = 0.75;
    double num_t0 = 10;
    double num_gamma = 0.05;
    int int_M_adapt = pl.M_adapt;
    vec vec_M_diag(int_dim);
    double num_hmcmu;
    double num_H;
    double num_epsbar;
    if(pl.M_diag.has_nan())
    {
        vec_M_diag.fill(1.0);
    } else{
        vec_M_diag = pl.M_diag;
    }

    if(int_iter == 0)
    {
        num_eps = find_reasonable_epsilon_cpp(nuts_model, dat_list, vec_theta, vec_M_diag, num_eps);
        num_hmcmu = log(10.0*num_eps);
        num_H = 0.0;
        num_epsbar = 1.0;
    } else{
        num_eps = pl.eps;
        num_epsbar = pl.eps_bar;
        num_H = pl.H;
        num_hmcmu = pl.hmcmu;
    }
    // r0 <- rnorm(length(theta), 0, sqrt(M_diag))
    // u <- runif(1, 0, exp(joint_log_density(theta, r0, f, M_diag, dat)))
    vec vec_r0 = randn(int_dim)%sqrt(vec_M_diag);
    // Rcout<<"for debug, use constant vec_r0"<<endl;
    // vec vec_r0(int_dim);
    // vec_r0.load("vec_r0.txt");
    double num_u = exp(nuts_model.joint_fn_val(vec_theta, vec_r0, vec_M_diag, dat_list));
    num_u = randu()*num_u;
    if(isnan(num_u))
    {
        Rcout<<"NUTS: sampled slice u is NaN"<<endl;
        num_u = randu()*1e5;
    }

    vec vec_theta_minus = vec_theta;
    vec vec_theta_plus = vec_theta;
    vec vec_r_minus = vec_r0;
    vec vec_r_plus = vec_r0;
    int int_j=0;
    int int_n=1;
    int int_s=1;

    if(int_iter > int_M_adapt-1)
    {
        num_eps = 0.9*num_epsbar + randu()*0.2*num_epsbar;// runif(1, 0.9*eps_bar, 1.1*eps_bar)
    }

    hmctree temp;
    while(int_s == 1)
    {
        int int_direction = randi(distr_param(0,1));
        if (int_direction==0)
        {
            int_direction = -1;
        }// sample(c(-1, 1), 1)
        // Rcout<<"for debug set int_direction = "<<-1<<endl;
        // int_direction = -1;
        if(int_direction == -1)
        {
            temp = build_tree_cpp(nuts_model, dat_list, vec_theta_minus, vec_r_minus, num_u, int_direction, int_j, num_eps, vec_theta, vec_r0, vec_M_diag);
            vec_theta_minus = temp.theta_minus;
            vec_r_minus = temp.r_minus;
        } else{
            temp = build_tree_cpp(nuts_model, dat_list, vec_theta_plus, vec_r_plus, num_u, int_direction, int_j, num_eps, vec_theta, vec_r0, vec_M_diag);
            vec_theta_plus = temp.theta_plus;
            vec_r_plus = temp.r_plus;
        }
        if(isnan(temp.s))
        {
            temp.s = 0;
        }

        if(temp.s == 1)
        {
            if(randu() < min(1.0, (double) temp.n/(double) int_n))
            {
                vec_theta = temp.theta;
            }
            // {
            //     Rcout<<"ratio = "<< (double) temp.n/(double) int_n<<endl;
            //     Rcout<<"temp.n = "<<temp.n<<"\t int_n = "<<int_n<<endl;
            //     Rcout<<"new theta not accepted!"<<endl;
            // }
            
        }
        int_n = int_n + temp.n;
        int_s = check_NUTS_cpp(temp.s, vec_theta_plus, vec_theta_minus, vec_r_plus, vec_r_minus);
        int_j = int_j + 1;
        if(int_j > max_treedepth)
        {
            // warning("NUTS: Reached max tree depth")
            // Rcout<<"NUTS: Reached max tree depth, int_iter = "<<int_iter<<endl;
            break;
        }
    }
    double log_eps;
    if(int_iter <= int_M_adapt-1)
    {
        
        adapt(num_H, log_eps, num_epsbar, num_eps, int_iter, num_t0, num_delta, num_gamma, num_kappa, num_hmcmu, temp.alpha, temp.n_alpha);
        
    } else{
        num_eps = num_epsbar;
    }
    
    //Rcout<<"iter = "<<int_iter<<"  num_eps = "<<num_eps<<" "<<temp.alpha<<" "<<temp.n_alpha<<endl; 
    
    thetaparlist res;
    parlist respl;
    respl.eps = num_eps;
    respl.eps_bar = num_epsbar;
    respl.H = num_H;
    respl.hmcmu = num_hmcmu;
    respl.M_adapt = int_M_adapt;
    respl.M_diag = vec_M_diag;

    res.theta = vec_theta;
    res.pl = respl;
    // res=list(theta = theta, pars = list(eps = eps, eps_bar = eps_bar, H = H, mu = mu, M_adapt = M_adapt, M_diag = M_diag))
    return res;
}

// NUTS <- function(theta, f, grad_f, dat, n_iter, M_diag = NULL, M_adapt = 50, delta = 0.5, max_treedepth = 10, eps = 1, verbose = TRUE)

mat Nuts_cpp(nutsModel nuts_model, const List& dat_list, vec vec_theta, int int_iter, vec vec_M_diag, int int_M_adapt = 50, double num_delta = 0.5, int int_max_treedepth = 10, double num_eps = 1)
{
    int int_ndim = vec_theta.n_elem;
    mat mat_theta_trace(int_ndim, int_iter, fill::ones);
    //par_list <- list(M_adapt = M_adapt)
    parlist pl(int_ndim);
    pl.M_diag = vec_M_diag;
    for(int i = 0; i<int_iter; i++)
    {
        //nuts <- NUTS_one_step(theta, iter, f, grad_f, par_list, delta = delta, max_treedepth = max_treedepth, eps = eps, verbose = verbose, dat)
        thetaparlist nuts = NUTS_one_step_cpp(nuts_model, dat_list, vec_theta, i, pl, num_delta, int_max_treedepth, num_eps);
        vec_theta = nuts.theta;
        mat_theta_trace.col(i) = nuts.theta;
        pl = nuts.pl;
        
        
        if((i==0)| (i%1000 == 999))
        {
            Rcout<<"iter = "<<i + 1<<" finished."<<" eps = "<<nuts.pl.eps<<endl;
        }
    }
    return mat_theta_trace.t();
}      

/*------------------------------------------------*/
//customized code for specific models
           
double my_fn(const vec& vec_theta, const List& dat_list)
{
   double num_val;
   double num_mu=dat_list["num_mu"];
   vec vec_w=dat_list["vec_w"];
   mat mat_eps_inv=dat_list["mat_eps_inv"]; 
   mat mat_alp_inv=dat_list["mat_alp_inv"];
   
   num_val = as_scalar(vec_theta.t() * mat_alp_inv * vec_theta);
   num_val += as_scalar((vec_w - num_mu - vec_theta).t() * mat_eps_inv * (vec_w - num_mu - vec_theta));
   return num_val * (-0.5);
}   

vec my_gn(const vec& vec_theta, const List& dat_list)
{
   vec vec_val;
   double num_mu=dat_list["num_mu"];
   vec vec_w=dat_list["vec_w"];
   mat mat_eps_inv=dat_list["mat_eps_inv"]; 
   mat mat_alp_inv=dat_list["mat_alp_inv"];
   
   vec_val = 2.0 * mat_alp_inv * vec_theta;
   vec_val = vec_val - 2.0*mat_eps_inv * (vec_w - num_mu - vec_theta);
   return vec_val*(-0.5);
}

double normal_fn(const vec& vec_theta, const List& dat_list)
{
   double num_val;

   vec mu=dat_list["mu"];
   mat sigma_inv=dat_list["sigma_inv"]; 
   
   num_val = as_scalar((vec_theta-mu).t() * sigma_inv * (vec_theta-mu));
   return num_val * (-0.5);
} 

 
vec normal_gn(const vec& vec_theta, const List& dat_list)
{
   vec vec_val;

   vec mu=dat_list["mu"];
   mat sigma_inv=dat_list["sigma_inv"]; 
   
   vec_val = sigma_inv * (vec_theta-mu);
   return vec_val*(-1.0);
}
 

double loglik_ld_val(distFuns dfuns, const mat& Ymat, const vec& mu_vec, const vec& sigma_vec)
{	
	int nn=Ymat.n_rows;
	vec ll(nn, fill::zeros);
  vec delta(nn, fill::zeros);
  
  vec tt=Ymat.col(2); //respU
  vec censor_cate=Ymat.col(3);
  vec wts_col=Ymat.col(4); 
  vec trun_cate=Ymat.col(7);
  
  vec zL=(dfuns.xtran(Ymat.col(1))-mu_vec)/sigma_vec;
  vec zU=(dfuns.xtran(Ymat.col(2))-mu_vec)/sigma_vec;
      
  vec dd=dfuns.phi(zU)/(sigma_vec % dfuns.dent(tt));
    
  vec pzU=dfuns.Phi(zU);
  vec pzL=dfuns.Phi(zL);

  uvec idx = find(censor_cate==1.0); //1 is none for censoring type
  delta.elem(idx).fill(1.0);
  
  idx = find(censor_cate==3.0); //left
  pzL.elem(idx).fill(0.0);

  idx = find(censor_cate==2.0); //right
  pzU.elem(idx).fill(1.0);
    
  vec ppc=pzU-pzL;
  
  //
  vec ztL=(dfuns.xtran(Ymat.col(5))-mu_vec)/sigma_vec;
  vec ztR=(dfuns.xtran(Ymat.col(6))-mu_vec)/sigma_vec;  

  vec pztL=dfuns.Phi(ztL);
  vec pztR=dfuns.Phi(ztR);
  
  idx = find((trun_cate==2.0) || (trun_cate==1.0)); // right or none
  pztL.elem(idx).fill(0.0);

  idx = find((trun_cate==3.0) || (trun_cate==1.0)); // left or none
  pztR.elem(idx).fill(1.0);
  
  vec ppt=pztR-pztL;
  
  vec term1=log(dd)-log(ppt);     //exact
  vec term2=log(ppc)-log(ppt);    //censored
 
  idx = find(delta==1.0);
  ll.elem(idx)=term1.elem(idx);
  
  idx = find(delta!=1.0);
  ll.elem(idx)=term2.elem(idx);
  
  //replace -inf values
  idx = find(ll== (-1)*datum::inf);
  vec rep_m_inf(nn);
  rep_m_inf.fill(-743.7469);
  ll.elem(idx)=rep_m_inf.elem(idx);
  
  double res=as_scalar(sum(ll%wts_col));
  
  //Rcout<<"beta0 is: "<<mu_vec(0)<<". loglik is: "<<sum(ll)<<endl;
  
  return(res);
}

    
vec loglik_ld_der(distFuns dfuns, const mat& Ymat, const vec& mu_vec, const vec& sigma_vec)	
{
	int nn=Ymat.n_rows;
	vec res(nn+1);

	vec vone(nn, fill::ones);
	mat der(nn, 2, fill::zeros);
  vec tres(nn, fill::zeros);
    
  vec delta(nn, fill::zeros);
  //vec tt=Ymat.col(2); //respU
  vec censor_cate=Ymat.col(3);
  vec wts_col=Ymat.col(4); 
  vec trun_cate=Ymat.col(7);
  
  vec zL=(dfuns.xtran(Ymat.col(1))-mu_vec)/sigma_vec;
  vec zU=(dfuns.xtran(Ymat.col(2))-mu_vec)/sigma_vec;
  
  //Rcout<<Ymat.row(0)<<endl;

  vec dd=(-1)*dfuns.dphiphi(zU)/sigma_vec;     //wrt mu
  vec s_dd=(-1)*dfuns.dphiphi(zU)%zU;   //wrt log(sigma)
  
  //mat debug1 = join_horiz(zU, mu_vec, sigma_vec);
  //Rcout<< debug1.row(0) <<endl;
  
      
  vec pzU=dfuns.Phi(zU);
  vec pzL=dfuns.Phi(zL);
  
  vec dzU=(-1)*dfuns.phi(zU)/sigma_vec;
  vec dzL=(-1)*dfuns.phi(zL)/sigma_vec;

  vec s_dzU=(-1)*dfuns.phi(zU)%zU;
  vec s_dzL=(-1)*dfuns.phi(zL)%zL;
  
  uvec idx = find(censor_cate==1.0); //1 is none for censoring type
  delta.elem(idx).fill(1.0);

  idx = find(censor_cate==2.0); //right
  pzU.elem(idx).fill(1.0);
  dzU.elem(idx).fill(0.0);
  s_dzU.elem(idx).fill(0.0);
  
  idx = find(censor_cate==3.0); //left
  pzL.elem(idx).fill(0.0);
  dzL.elem(idx).fill(0.0);
  s_dzL.elem(idx).fill(0.0);
  
  vec ppc=pzU-pzL;
  vec dpc=dzU-dzL;    
  vec s_dpc=s_dzU-s_dzL; 
      
  /*-----------------------------*/
  vec ztL=(dfuns.xtran(Ymat.col(5))-mu_vec)/sigma_vec;
  vec ztR=(dfuns.xtran(Ymat.col(6))-mu_vec)/sigma_vec;  
  
  vec pztL=dfuns.Phi(ztL);
  vec pztR=dfuns.Phi(ztR);
  
  vec dztL=(-1)*dfuns.phi(ztL)/sigma_vec;
  vec dztR=(-1)*dfuns.phi(ztR)/sigma_vec;  

  vec s_dztL=(-1)*dfuns.phi(ztL)%ztL;
  vec s_dztR=(-1)*dfuns.phi(ztR)%ztR;    

  idx = find((trun_cate==2.0) || (trun_cate==1.0)); // right or none
  pztL.elem(idx).fill(0.0);
  dztL.elem(idx).fill(0.0);
  s_dztL.elem(idx).fill(0.0);
      
  idx = find((trun_cate==3.0) || (trun_cate==1.0)); // left or none
  pztR.elem(idx).fill(1.0);
  dztR.elem(idx).fill(0.0);
  s_dztR.elem(idx).fill(0.0);

  vec ppt=pztR-pztL;
  vec dpt=dztR-dztL;
  vec s_dpt=s_dztR-s_dztL;
  
  //for mu  
  vec term1=dd-(dpt/ppt);     //exact
  vec term2=dpc/ppc-(dpt/ppt);    //censored

  //mat debug =join_horiz(join_horiz(dd, dpc, ppc), dpt, ppt);
  //Rcout<< debug.row(0) <<endl;

  idx = find(delta==1.0);
  tres.elem(idx)=term1.elem(idx);

  idx = find(delta!=1.0);
  tres.elem(idx)=term2.elem(idx);
  
  der.col(0)=tres;
  
  //for log(sigma)
  term1=(-1)*vone+s_dd-(s_dpt/ppt);     //exact
  term2=s_dpc/ppc-(s_dpt/ppt);    //censored

  idx = find(delta==1.0);
  tres.elem(idx)=term1.elem(idx);

  idx = find(delta!=1.0);
  tres.elem(idx)=term2.elem(idx);  
  
  der.col(1)=tres;
  
  //derivatives
  res.subvec(0, nn-1)=der.col(0)%wts_col;
  res(nn)=as_scalar(sum(der.col(1)%wts_col));
  
  //Rcout<<"res der"<<endl;
  //Rcout<<der.col(0)<<endl;
  

	return(res);
}	
	

//general likelihood function

double gen_loglik_lifedata(const vec& pars, const List& dat)
{  
  mat Ymat=dat["Ymat"]; 
  String dist_name=dat["dist_name"]; 
  String model_name=dat["model_name"];
  
  //Rcout<<"pars: \n"<< pars<<endl;
  //Rcout<<"Ymat\n"<< Ymat<<endl;
   
  distFuns dfuns(dist_name);
  
  //Rcout<<pars<<endl;
  //Rcout<<dfuns.Phi(pars)<<endl;
  
  int nn=Ymat.n_rows;
  vec mu_vec(nn);
  vec sigma_vec(nn);
  double res=0.0;
  
  //ld_fit
  if(model_name=="ld_fit")
  {
    double mu=pars(0);
    double sigma=exp(pars(1));
    mu_vec.fill(mu);
    sigma_vec.fill(sigma);
  
    res=loglik_ld_val(dfuns, Ymat, mu_vec, sigma_vec);  
  }
  
  // ldme_fit
  if(model_name=="ldme_fit")
  {
    mat xmat=dat["xmat"]; 
    mat zmat=dat["zmat"]; 
    
    int n_fixed=xmat.n_cols;
    int n_rand=zmat.n_cols;
    
    vec beta_vec=pars.subvec(0, n_fixed-1);
    vec w_vec=pars.subvec(n_fixed, n_fixed+n_rand-1);
    
    double sigma=exp(pars(n_fixed+n_rand-1+1));
    double t_sigma_w=pars(n_fixed+n_rand-1+2);
    double sigma_w=exp(t_sigma_w);
        
    mu_vec=xmat*beta_vec+zmat*w_vec; 
    sigma_vec.fill(sigma);
    
    double term1=loglik_ld_val(dfuns, Ymat, mu_vec, sigma_vec); 
  
    double term2=(-1)*n_rand*log(sigma_w)-0.5*(as_scalar(sum(w_vec%w_vec)))/(sigma_w*sigma_w);
    double term3=(-2)*t_sigma_w-1/(sigma_w*sigma_w);
    
    res=term1+term2+term3;
  
  }  

  /*
  if(model_name=="ldsp0_fit")
  {
    mat xmat=dat["xmat"]; 
    mat zmat=dat["zmat"];
    mat dmat=dat["dmat"]; 
    
    int n_fixed=xmat.n_cols;
    int n_rand=zmat.n_cols;
    
    vec beta_vec=pars.subvec(0, n_fixed-1);
    vec w_vec=pars.subvec(n_fixed, n_fixed+n_rand-1);
    
    double sigma=exp(pars(n_fixed+n_rand-1+1));
    double t_sigma_w=pars(n_fixed+n_rand-1+2);
    double sigma_w=exp(t_sigma_w);
    double t_nu=pars(n_fixed+n_rand-1+3);
    double nu=exp(t_nu);
    
    double t_kappa=pars(n_fixed+n_rand-1+4);
    double kappa=2*exp(t_kappa)/(1+exp(t_kappa));
            
    mu_vec=xmat*beta_vec+zmat*w_vec; 
    sigma_vec.fill(sigma);
    
    double term1=loglik_ld_val(dfuns, Ymat, mu_vec, sigma_vec); 
  
    //Rcout << "nu="<<nu<<" kappa="<< kappa << endl;
    
    mat bb=pow(dmat/nu, kappa);
    mat omat=exp((-1)*bb);
    omat.diag().fill(1+0.000001);
    
    //Rcout << omat.is_symmetric() << endl;
    
    mat omat_inv=inv(omat);
  
    double bb1=log(det(omat));
    
    double term2=(-1)*n_rand*log(sigma_w)-0.5*bb1-0.5*(as_scalar(w_vec.t()*omat_inv*w_vec))/(sigma_w*sigma_w);
    
    double term3=(-2)*t_sigma_w-1/(sigma_w*sigma_w)-3*t_nu-1/nu;
    
    double term4=t_kappa-2*log(1+exp(t_kappa));
    
    res=term1+term2+term3+term4;
    
    //Rcout<<"term1="<<term1<<" term2="<<term2<<" term3="<<term3<<endl;
  
  }  
  */

  if(model_name=="ldsp_fit")
  {
    mat xmat=dat["xmat"]; 
    mat zmat=dat["zmat"];
    mat dmat=dat["dmat"];
    mat cmat=dat["cmat"];
    vec ig_par=dat["ig_par"]; 
    vec beta_par=dat["beta_par"]; 
    
    int n_fixed=xmat.n_cols;
    int n_rand=zmat.n_cols-1;    //sum to zero constraint.
    
    vec beta_vec=pars.subvec(0, n_fixed-1);
    vec w_vec=pars.subvec(n_fixed, n_fixed+n_rand-1);  
    vec w_vec0=cmat*w_vec;
    
    double sigma=exp(pars(n_fixed+n_rand-1+1));
    double t_sigma_w=pars(n_fixed+n_rand-1+2);
    double sigma_w=exp(t_sigma_w);
    double t_nu=pars(n_fixed+n_rand-1+3);
    double nu=exp(t_nu);
    
    double t_kappa=pars(n_fixed+n_rand-1+4);
    double kappa=2*exp(t_kappa)/(1+exp(t_kappa));
            
    mu_vec=xmat*beta_vec+zmat*w_vec0; 
    sigma_vec.fill(sigma);
    
    double term1=loglik_ld_val(dfuns, Ymat, mu_vec, sigma_vec); 
  
    mat bb=pow(dmat/nu, kappa);
    mat omat0=exp((-1)*bb);
    omat0.diag().fill(1+0.000001);
    
    mat omat=cmat.t()*omat0*cmat;
    mat omat_inv=inv(omat);
    
    //Rcout << "nu="<<nu<<" kappa="<< kappa << endl;
    //Rcout << omat.is_symmetric() << endl;
    
    //Rcout<<omat.submat(0, 0, 10, 10)<<endl;
  
    double bb1=log(det(omat));
    
    //Rcout << "bb1="<<bb1<<endl;
    
    double term2=(-1.0)*n_rand*log(sigma_w)-0.5*bb1-0.5*(as_scalar(w_vec.t()*omat_inv*w_vec))/(sigma_w*sigma_w);
    
    double ig_a=ig_par(0);
    double ig_b=ig_par(1);
    double term3=(-2.0)*t_sigma_w-1.0/(sigma_w*sigma_w)-ig_a*t_nu-ig_b/nu;

    double beta_a=beta_par(0);
    double beta_b=beta_par(1);    
    double term4=beta_a*t_kappa-(beta_a+beta_b)*log(1.0+exp(t_kappa));
    
    res=term1+term2+term3+term4;
    
    //Rcout<<"term1="<<term1<<" term2="<<term2<<" term3="<<term3<<endl;
    //Rcout<<"loglik="<<res<<endl;
  
  }  


  
  return(res);

}

//derivatives

vec gen_loglik_lifedata_der(const vec& pars, const List& dat)
{
  int np=pars.n_elem;
  vec res(np, fill::zeros);

  mat Ymat=dat["Ymat"]; 
  
  String dist_name=dat["dist_name"]; 
  String model_name=dat["model_name"];
     
  distFuns dfuns(dist_name);
  
  int nn=Ymat.n_rows; 
  vec mu_vec(nn);
  vec sigma_vec(nn);
  
  
  if(model_name=="ld_fit")
  {
    double mu=pars(0);
    double sigma=exp(pars(1));
    mu_vec.fill(mu);
    sigma_vec.fill(sigma);
  
    vec tmp_der=loglik_ld_der(dfuns, Ymat, mu_vec, sigma_vec); 
    
    vec mu_der=tmp_der.subvec(0, nn-1); 
  
    res(0)=sum(mu_der);
    res(1)=tmp_der(nn);
  
  }

  //ldme_fit
  if(model_name=="ldme_fit")
  {
    mat xmat=dat["xmat"]; 
    mat zmat=dat["zmat"];
     
    int n_fixed=xmat.n_cols;
    int n_rand=zmat.n_cols;
    
    vec beta_vec=pars.subvec(0, n_fixed-1);
    vec w_vec=pars.subvec(n_fixed, n_fixed+n_rand-1);
    
    double sigma=exp(pars(n_fixed+n_rand-1+1));
    double t_sigma_w=pars(n_fixed+n_rand-1+2);
    double sigma_w=exp(t_sigma_w);
        
    mu_vec=xmat*beta_vec+zmat*w_vec; 
    sigma_vec.fill(sigma);
    
    vec tmp_der=loglik_ld_der(dfuns, Ymat, mu_vec, sigma_vec); 
    
    vec mu_der=tmp_der.subvec(0, nn-1);
    
    //Rcout<<n_fixed<<endl;
    
    for(int i=0; i<n_fixed; i++)
    {
      //Rcout<<"i="<<i<<endl;
      res(i)=as_scalar(sum(mu_der%xmat.col(i)));
    }
    
    for(int i=0; i<n_rand; i++)
    {
      res(n_fixed+i)=as_scalar(sum(mu_der%zmat.col(i)))-w_vec(i)/(sigma_w*sigma_w);
    }
       
    res(n_fixed+n_rand-1+1)=tmp_der(nn);
    res(n_fixed+n_rand-1+2)=(-1)*n_rand+(as_scalar(sum(w_vec%w_vec)))/(sigma_w*sigma_w)-2+2/(sigma_w*sigma_w);
  
  }


  //ldsp_fit
  if(model_name=="ldsp_fit")
  {
    mat xmat=dat["xmat"]; 
    mat zmat=dat["zmat"];
    mat dmat=dat["dmat"];
    mat cmat=dat["cmat"]; 
    vec ig_par=dat["ig_par"]; 
    vec beta_par=dat["beta_par"];    
    
    int n_fixed=xmat.n_cols;
    int n_rand=zmat.n_cols-1;
    
    vec beta_vec=pars.subvec(0, n_fixed-1);
    vec w_vec=pars.subvec(n_fixed, n_fixed+n_rand-1);
    vec w_vec0=cmat*w_vec;
    
    double sigma=exp(pars(n_fixed+n_rand-1+1));
    double t_sigma_w=pars(n_fixed+n_rand-1+2);
    double sigma_w=exp(t_sigma_w);
    double t_nu=pars(n_fixed+n_rand-1+3);
    double nu=exp(t_nu);
    double t_kappa=pars(n_fixed+n_rand-1+4);
    double kappa=2*exp(t_kappa)/(1+exp(t_kappa));
        
    mu_vec=xmat*beta_vec+zmat*w_vec0; 
    sigma_vec.fill(sigma);
    
    vec tmp_der=loglik_ld_der(dfuns, Ymat, mu_vec, sigma_vec); 
    
    vec mu_der=tmp_der.subvec(0, nn-1);
        
    for(int i=0; i<n_fixed; i++)
    {
      //Rcout<<"i="<<i<<endl;
      res(i)=as_scalar(sum(mu_der%xmat.col(i)));
    }
    
    mat ff=dmat/nu;
    mat bb=pow(ff, kappa);
     
    mat omat0=exp((-1)*bb);
    omat0.diag().fill(1+0.000001);
    
    mat omat=cmat.t()*omat0*cmat;
    mat omat_inv=inv(omat);
  
    //double bb1=log(det(omat));
    vec bb3=omat_inv*w_vec;
       
    for(int i=0; i<n_rand; i++)
    {
      res(n_fixed+i)=as_scalar(sum(mu_der%zmat.col(i)))-bb3(i)/(sigma_w*sigma_w);
    }
       
    //sigma
    res(n_fixed+n_rand-1+1)=tmp_der(nn);
    
    double bb4=as_scalar(w_vec.t()*omat_inv*w_vec);
    
    //sigma_w    
    res(n_fixed+n_rand-1+2)=(-1)*n_rand+bb4/(sigma_w*sigma_w)-2+2/(sigma_w*sigma_w);
    
    //nu
    mat bb5=kappa*(omat0%bb);
    double cc1=(-0.5)*trace(omat_inv*cmat.t()*bb5*cmat);
    double cc2=0.5*(1/(sigma_w*sigma_w))*as_scalar(bb3.t()*cmat.t()*bb5*cmat*bb3);
    
    double ig_a=ig_par(0);
    double ig_b=ig_par(1);
    
    res(n_fixed+n_rand-1+3)=cc1+cc2-ig_a+ig_b/nu;
    
    //kappa
    mat ldmat=log(dmat);
    ldmat.diag().zeros();
    bb5=(-1)*kappa*(1-kappa/2)*(omat0%bb%(ldmat-log(nu)));
    
    cc1=(-0.5)*trace(omat_inv*cmat.t()*bb5*cmat);
    cc2=0.5*(1/(sigma_w*sigma_w))*as_scalar(bb3.t()*cmat.t()*bb5*cmat*bb3);
    
    double beta_a=beta_par(0);
    double beta_b=beta_par(1);
    res(n_fixed+n_rand-1+4)=cc1+cc2+beta_a-(beta_a+beta_b)*kappa/2.0;
    
  }


  return res;    
}


double gen_loglik_lifedata_joint(const vec& vec_theta, const vec& vec_r, const vec& M_diag, const List& dat_list)
{
   double num_val;
   num_val = gen_loglik_lifedata(vec_theta, dat_list)-0.5 * sum(vec_r % vec_r / M_diag);
   return num_val;
}


//[[Rcpp::export]]       
mat Nuts_lifedata_posterior(vec vec_theta, int int_iter, List dat_list, List nuts_ctrl_list)
{
  nutsModel nuts_model;
  
  nuts_model.fn_val=&gen_loglik_lifedata;  
  nuts_model.gn_val=&gen_loglik_lifedata_der;
  nuts_model.joint_fn_val=&gen_loglik_lifedata_joint;
    
  vec vec_M_diag=dat_list["vec_M_diag"];
  //read control parameters   
  int int_M_adapt=nuts_ctrl_list["int_M_adapt"];
  double num_delta=nuts_ctrl_list["num_delta"];
  int int_max_treedepth=nuts_ctrl_list["int_max_treedepth"];
  double num_eps=nuts_ctrl_list["num_eps"]; 
      
  return Nuts_cpp(nuts_model, dat_list, vec_theta, int_iter, vec_M_diag, int_M_adapt, num_delta, int_max_treedepth, num_eps);
   

}



/*-----------------------------------------------------------------*/
//customized code for specific models for Gibbs

double my_log_post_ratio(double yy, int j, const vec& th, const List& dat_list)
{
  double res, aa1, aa2;
  int i;
  
  mat dat=dat_list["dat"]; 
  mat start_end_mat=dat_list["start_end_mat"];
  vec theta=dat_list["theta"];
  
  double mu=theta[0];
  double sigma=theta[1];;
  double sigma_gamma=theta[2];
  
  aa1=0.0;
  aa2=0.0;
  
  double gamma0=th[j];
  double gamma1=th[j]+yy;
  
  for(i=start_end_mat(0, j); i<=start_end_mat(1, j); i++)
  {
     double tt=dat(i,0);
     int delta=dat(i,1);
       
     double zz0=(log(tt)-(mu+gamma0))/sigma;
     double zz1=(log(tt)-(mu+gamma1))/sigma;
     
     if(delta==1)
     {
        aa1+=log(dsev(zz1))-log(dsev(zz0));        
     }else{
            aa1+=log(1-psev(zz1))-log(1-psev(zz0));
          }
  }
  
  aa2=-0.5*(gamma1*gamma1-gamma0*gamma0)/(sigma_gamma*sigma_gamma);
 
  res=aa1+aa2;
  
  return res;
}

 
//[[Rcpp::export]]       
List my_Gibbs_rand_par(vec vec_theta, List dat_list, vec vec_scale, int int_iter)
{
  gibbsModel gibbs_model;
  gibbs_model.dat_list=dat_list;
  gibbs_model.log_post_ratio=&my_log_post_ratio;
      
  List res=Gibbs_adaptive(gibbs_model, vec_theta, vec_scale, int_iter, 500, 400);
  return res;
}

//Bayesian log post ratio
double my_bayes_log_post_ratio(double yy, int j, const vec& th, const List& dat_list)
{
  double res, aa1, aa2, zz, zz1, tt, dyy, gamma, gamma1;
  int i, delta;
  
  mat dat=dat_list["dat"]; 
  mat se_mat=dat_list["start_end_mat"];
  mat pmat=dat_list["pmat"];
  
  int nn=dat.n_rows;
  int pp=se_mat.n_cols;
  //int np=pmat.n_rows;
  
  //Rcout<<"nn="<<nn<<"  pp="<<pp<<"  np="<<np<<endl;  
  
  vec th1=th;
  th1[j]+=yy;
  
  double mu=th[pp+0];
  double sigma=exp(th[pp+1]);
  double sigma_g=exp(th[pp+2]);
  
  double mu1=th1[pp+0];
  double sigma1=exp(th1[pp+1]);
  double sigma_g1=exp(th1[pp+2]);

  aa1=0.0;
  aa2=0.0;
  
  if(j<pp)
  {
   gamma=th[j];
   gamma1=th1[j];
  
   for(i=se_mat(0, j); i<=se_mat(1, j); i++)
   {
     tt=dat(i,0);
     delta=dat(i,1);
       
     zz=(log(tt)-(mu+gamma))/sigma;
     zz1=(log(tt)-(mu1+gamma1))/sigma1;
     
     if(delta==1)
     {
        aa1+=log(dsev(zz1))-log(dsev(zz));        
     }else{
            aa1+=log(1-psev(zz1))-log(1-psev(zz));
          }
   }
  
   aa2=-0.5*(gamma1*gamma1-gamma*gamma)/(sigma_g*sigma_g);
  
  }
  
  if( (j==pp)|(j==(pp+1)))
  {
   dyy=log(sigma1)-log(sigma);
   
   for(i=0; i<nn; i++)
   {
     tt=dat(i,0);
     delta=dat(i,1);
       
     zz=(log(tt)-mu)/sigma;
     zz1=(log(tt)-mu1)/sigma1;
     
     if(delta==1)
     {
        aa1+=log(dsev(zz1))-log(dsev(zz))-dyy;        
     }else{
            aa1+=log(1-psev(zz1))-log(1-psev(zz));
          }
   }
   
   if(j==pp)
   {
     aa2=-0.5*(pow(mu1-pmat(0,0), 2)-pow(mu-pmat(0,0), 2))/pow(pmat(0,1), 2);
   }

   if(j==(pp+1))
   {
     aa2=-0.5*(pow(log(sigma1)-pmat(1,0), 2)-pow(log(sigma)-pmat(1,0), 2))/pow(pmat(1,1), 2);
   }
   
  }

  if(j==(pp+2))
  {
    dyy=log(sigma_g1)-log(sigma_g);
   
    for(i=0; i<pp; i++)
    {
      gamma=th[i];
      aa1+=-0.5*gamma*gamma*(1/(sigma_g1*sigma_g1)-1/(sigma_g*sigma_g))-dyy;        
     
    }

     aa2=-0.5*(pow(log(sigma_g1)-pmat(2,0), 2)-pow(log(sigma_g)-pmat(2,0), 2))/pow(pmat(2,1), 2);
  }
 
 
  res=aa1+aa2;
 
  return res;
}

 
//[[Rcpp::export]]       
List my_Bayes_Gibbs_rand_par(vec vec_theta, List dat_list, vec vec_scale, int int_iter)
{
  gibbsModel gibbs_model;
  gibbs_model.dat_list=dat_list;
  gibbs_model.log_post_ratio=&my_bayes_log_post_ratio;
      
  List res=Gibbs_adaptive(gibbs_model, vec_theta, vec_scale, int_iter, 500, 400);
  return res;
}  


double bayes_rand_par_igamma_log_post_ratio(double yy, int j, const vec& th, const List& dat_list)
{
  double res, aa1, aa2, zz, zz1, tt, dyy, gamma, gamma1;
  int i, delta;
  
  mat dat=dat_list["dat"]; 
  mat se_mat=dat_list["start_end_mat"];
  
  int nn=dat.n_rows;
  int pp=se_mat.n_cols;
  //int np=pmat.n_rows;
  
  //Rcout<<"nn="<<nn<<"  pp="<<pp<<"  np="<<np<<endl;  
  
  vec th1=th;
  th1[j]+=yy;
  
  double mu=th[pp+0];
  double sigma=exp(th[pp+1]);
  double sigma_g=exp(th[pp+2]);
  
  double mu1=th1[pp+0];
  double sigma1=exp(th1[pp+1]);
  double sigma_g1=exp(th1[pp+2]);

  aa1=0.0;
  aa2=0.0;
  
  if(j<pp)
  {
   gamma=th[j];
   gamma1=th1[j];
  
   for(i=se_mat(0, j); i<=se_mat(1, j); i++)
   {
     tt=dat(i,0);
     delta=dat(i,1);
       
     zz=(log(tt)-(mu+gamma))/sigma;
     zz1=(log(tt)-(mu1+gamma1))/sigma1;
     
     if(delta==1)
     {
        aa1+=log(dsev(zz1))-log(dsev(zz));        
     }else{
            aa1+=log(1-psev(zz1))-log(1-psev(zz));
          }
   }
  
   aa2=-0.5*(gamma1*gamma1-gamma*gamma)/(sigma_g*sigma_g);
  
  }
  
  if( (j==pp)|(j==(pp+1)))
  {
   dyy=log(sigma1)-log(sigma);
   
   for(i=0; i<nn; i++)
   {
     tt=dat(i,0);
     delta=dat(i,1);
       
     zz=(log(tt)-mu)/sigma;
     zz1=(log(tt)-mu1)/sigma1;
     
     if(delta==1)
     {
        aa1+=log(dsev(zz1))-log(dsev(zz))-dyy;        
     }else{
            aa1+=log(1-psev(zz1))-log(1-psev(zz));
          }
   }
   
   if(j==pp) //mu
   {
     aa2=0.0;
   }

   if(j==(pp+1))   //sigma
   {
     aa2=0.0;
   }
   
  }

  if(j==(pp+2))   //sigma_r
  {
    dyy=log(sigma_g1)-log(sigma_g);
   
    for(i=0; i<pp; i++)
    {
      gamma=th[i];
      aa1+=-0.5*gamma*gamma*(1/(sigma_g1*sigma_g1)-1/(sigma_g*sigma_g))-dyy;        
     
    }

     aa2=-2.0*log(sigma_g1)-1.0/(sigma_g1*sigma_g1)+2.0*log(sigma_g)+1.0/(sigma_g*sigma_g);
  }
 
 
  res=aa1+aa2;
 
  return res;
}



//[[Rcpp::export]]       
List Bayes_Gibbs_rand_par_IGamma(vec vec_theta, List dat_list, vec vec_scale, int int_iter)
{
  gibbsModel gibbs_model;
  gibbs_model.dat_list=dat_list;
  gibbs_model.log_post_ratio=&bayes_rand_par_igamma_log_post_ratio;
      
  List res=Gibbs_adaptive(gibbs_model, vec_theta, vec_scale, int_iter, 500, 400);
  return res;
}

double log_post_ratio_lifedata(double yy, int j, const vec& th, const List& dat_list)
{
  double res, aa0, aa1;
  
  vec th1=th;
  th1(j)=th(j)+yy;
  
  aa0=gen_loglik_lifedata(th, dat_list);
  aa1=gen_loglik_lifedata(th1, dat_list);
  
  res=aa1-aa0;
  
  //Rcout<<"j="<<j<<" yy="<<yy<<" aa1="<<aa1<<"   aa0="<<aa0<<endl;
    
  return res;
}

mat cov_mat_prep(const vec& val, int n_gamma)  
{
  mat Dmat(n_gamma, n_gamma, fill::zeros);
  mat rhomat(n_gamma, n_gamma, fill::ones);
  
  int nn=val.n_elem;
  
  vec a1=exp(val.subvec(0,n_gamma-1));
  
  for(int i=0; i<n_gamma; i++)
  {
    Dmat(i,i)=a1(i);
  } 
  
  if(n_gamma>1)
  {
    vec tmp=val.subvec(n_gamma, nn-1);
    vec a2=2*exp(tmp)/(1+exp(tmp))-1;
    
    int i, j, cc;
    
    cc=0;
    
    //Rcout<<"a2="<<a2<<endl;
    
    for(j=0; j<n_gamma-1; j++)
    {
       for(i=j+1; i<n_gamma; i++)
       {
         //Rcout<<"j="<<j<<" i="<<i<<endl;
              
         double tmp1=a2(cc+i-j-1);
         rhomat(i,j)=tmp1;
         rhomat(j, i)=tmp1;
         cc=cc+(n_gamma-j-1);
       }
    }
    
    //Rcout<<"Dmat="<<Dmat<<" rhomat="<<rhomat<<endl;
    
  }
  
  mat Sigma_mat=Dmat*rhomat*Dmat;

  return(Sigma_mat);
}

vec Dfun_loglogis(const vec& tt, const vec& alpha_vec, const vec& mu_vec, const mat& gamma_mat)
{
  vec t_alpha=exp(alpha_vec);
  //vec a1=(-1)*t_alpha(0)*exp(gamma_mat.col(0));
  //vec a2=(log(tt)-mu_vec-gamma_mat.col(1))/t_alpha(1);
  double a1=(-1)*t_alpha(0);
  vec a2=(log(tt)-mu_vec-gamma_mat.col(0))/(t_alpha(1)*exp(gamma_mat.col(1)));
  
  vec res=a1/(1.0+exp(-a2)); 
  return(res); 
}

vec Dfun_loglogis1(const vec& tt, const vec& alpha_vec, const vec& mu_vec, const mat& gamma_mat)
{
  vec t_alpha=exp(alpha_vec);
  //vec a1=(-1)*t_alpha(0)*exp(gamma_mat.col(0));
  //vec a2=(log(tt)-mu_vec)/t_alpha(1);
  double a1=(-1)*t_alpha(0);
  vec a2=(log(tt)-mu_vec-gamma_mat.col(0))/t_alpha(1);
  
  vec res=a1/(1.0+exp(-a2)); 
  return(res); 
}

vec Dfun_lnrate(const vec& tt, const vec& alpha_vec, const vec& mu_vec, const mat& gamma_mat)
{
  vec res=alpha_vec(0)+exp(mu_vec+gamma_mat.col(0))%tt;
  
  //Rcout<<"alpha="<<alpha_vec<<endl;
  
  return(res); 
}

vec Dfun_concave(const vec& tt, const vec& alpha_vec, const vec& mu_vec, const mat& gamma_mat)
{
  vec t_alpha=exp(alpha_vec);
  vec a1=(-1)*t_alpha(0)*exp(gamma_mat.col(0));
  vec a2=log(tt)-mu_vec-gamma_mat.col(1);  
  vec res=a1%(1.0-exp(-exp(a2))); 
  
  //Rcout<<"alpha="<<t_alpha<<endl;

  return(res); 
}


class rmdtFuns
{
   public:
       vec (*Dfun)(const vec&, const vec&, const vec&, const mat&);
       int n_alpha;       //no. of fixed pars, minus no. in betavec
       int n_gamma;       //no. of random pars
       int trend;         //1 for increasing, -1 for decreasing paths.
       String mname;      //model name
                         
       rmdtFuns(String dd)
       {
         mname=dd;
         if(dd=="rmdt_loglogis")   //rmdt.loglogis 
         {
           Dfun=&Dfun_loglogis;
           n_alpha=2;
           n_gamma=2;
           trend= -1;

         }       

         if(dd=="rmdt_loglogis1")   //rmdt.loglogis1 
         {
           Dfun=&Dfun_loglogis1;
           n_alpha=2;
           n_gamma=1;
           trend= -1;

         } 

         if(dd=="rmdt_lnrate")   //rmdt.lnrate 
         {
           Dfun=&Dfun_lnrate;
           n_alpha=1;
           n_gamma=1;
           trend= 1;

         }                   

         if(dd=="rmdt_concave")   //rmdt.concave 
         {
           Dfun=&Dfun_concave;
           n_alpha=1;
           n_gamma=2;
           trend= -1;

         }   
         
       } 
        
};



vec gompertz_fun(const vec& tt, double alpha, double kappa, double beta)
{
   vec res;
   vec tmp=pow_vec(beta, tt);
   vec tmp1=pow_vec(kappa, tmp);
   
   res=alpha*tmp1-alpha*kappa;

   return(res);

}

vec gompertz_der_fun(const vec& tt, double alpha, double kappa, double beta)
{
   vec res;
   vec tmp=pow_vec(beta, tt);
   vec tmp1=pow_vec(kappa, tmp);
   vec tmp2=tmp%tmp1;
   
   res=alpha*tmp2*log(kappa)*log(beta);
   
   return(res);
}


vec wei_altp_fun(const vec& tt, double alpha, double mu, double sigma)
{
   vec zz=(log(tt)-mu)/sigma;
   vec res=alpha*psev_vec(zz);

   return(res);
}

vec wei_altp_der_fun(const vec& tt, double alpha, double mu, double sigma)
{
   vec zz=(log(tt)-mu)/sigma;
   vec tmp=1/(sigma*tt);
   vec res=alpha*dsev_vec(zz)%tmp;
   
   return(res);
}

double loglik_recurr_gomp(const vec& pars, const List& dat)
{  
  mat ne_mat=dat["Non_event_mat"];
  mat emat=dat["Event_mat"];
  
  double alpha=exp(pars(0));
  double beta=exp(pars(1))/(1+exp(pars(1))); 
  double kappa=exp(pars(2))/(1+exp(pars(2)));
   
  vec aa1=gompertz_der_fun(emat.col(2), alpha, kappa, beta);
  double pp1=(-1)*sum(log(aa1%emat.col(3)%emat.col(4)));
   
  vec bb1=gompertz_fun(ne_mat.col(1), alpha, kappa, beta);
  vec bb2=gompertz_fun(ne_mat.col(2), alpha, kappa, beta);

  vec dd=bb2-bb1; 
  double mdd=max(dd);
   
  uvec idx = find(dd==0.0); 
  dd.elem(idx).fill(mdd);
   
  double pp2=sum(ne_mat.col(3)%dd%ne_mat.col(4));
   
  double res=(-1)*(pp1+pp2);
   
  return(res);
   
}    


double loglik_recurr_wei_altp(const vec& pars, const List& dat)
{  
  mat ne_mat=dat["Non_event_mat"];
  mat emat=dat["Event_mat"];
  vec p_inf=dat["prior_inf"];
  
  double alpha=exp(pars(0));
  double mu=pars(1); 
  double sigma=exp(pars(2));
  
  double ww0=(pars(0)-p_inf(0))/p_inf(1);
  double ww1=(pars(1)-p_inf(2))/p_inf(3);
  double ww2=(pars(2)-p_inf(4))/p_inf(5);  
  
  double cc0=dnor(ww0)/p_inf(1);
  double cc1=dnor(ww1)/p_inf(3);
  double cc2=dnor(ww2)/p_inf(5);  
  double pp3=log(cc0)+log(cc1)+log(cc2);
   
  vec aa1=wei_altp_der_fun(emat.col(2), alpha, mu, sigma);
  double pp1=sum(log(aa1%emat.col(3)%emat.col(4)));
   
  vec bb1=wei_altp_fun(ne_mat.col(1), alpha, mu, sigma);
  vec bb2=wei_altp_fun(ne_mat.col(2), alpha, mu, sigma);
   
  vec dd=bb2-bb1; 
  double mdd=max(dd);
   
  uvec idx = find(dd==0.0); 
  dd.elem(idx).fill(mdd);
   
  double pp2=(-1)*sum(ne_mat.col(3)%dd%ne_mat.col(4));
   
  double res=pp1+pp2+pp3;
  
  //Rcout<<pars<<endl;
  //Rcout<<res<<endl;
   
  return(res);
   
}    


double gen_loglik_rmdt(const vec& pars, const List& dat)
{  
  mat Ymat=dat["Ymat"]; 
  String model_name=dat["model_name"];
  
  //Rcout<<"pars: \n"<< pars<<endl;
  //Rcout<<"Ymat\n"<< Ymat<<endl;
   
  rmdtFuns dfuns(model_name);
  
  //Rcout<<pars<<endl;
  
  int nn=Ymat.n_rows;
  vec mu_vec(nn);
  //vec sigma_vec(nn);
  double res=0.0;
  
  //Rcout << "ok"<<endl;
 
  //if((model_name=="rmdt_loglogis") | (model_name=="rmdt_loglogis1") | (model_name=="rmdt_lnrate"))
  //{
    mat xmat=dat["xmat"]; 
    mat zmat=dat["zmat"];
    mat cmat=dat["cmat"];
    
    double hcau_par=dat["hcau_par"]; 
    int n_alpha=dfuns.n_alpha; 
    int n_gamma=dfuns.n_gamma;
    int n_pars=pars.n_elem;
    
    //Rcout << "n_alpha="<<n_alpha<<endl;
        
    int n_fixed=xmat.n_cols;
    int n_rand=zmat.n_cols;     
    int n_randc=cmat.n_cols;     

    vec alpha_vec=pars.subvec(0, n_alpha-1);  
    //alpha_vec=exp(alpha_vec);
    //Rcout << "alpha="<<alpha_vec<<endl;
          
    vec beta_vec=pars.subvec(n_alpha, n_alpha+n_fixed-1);
    
    //Rcout << "beta="<<beta_vec<<endl;
    
    double sigma=pars(n_alpha+n_fixed);
    sigma=exp(sigma);    
    //Rcout << "sigma="<<sigma<<endl;

    int idx=n_alpha+n_fixed+1;
    vec w_vec=pars.subvec(idx, idx+n_gamma*n_randc-1); 
    
    mat w_matc=w_vec;
    w_matc.reshape(n_gamma, n_randc);  //gamma_vec

    mat w_mat=w_matc*cmat.t();
      
    int idx1=idx+n_gamma*n_randc;
    vec sigmak=pars.subvec(idx1, idx1+n_gamma-1);
    sigmak=exp(sigmak);
    
    vec val=pars.subvec(idx1, n_pars-1);
    mat Sigma_mat=cov_mat_prep(val, n_gamma);
         
    //Rcout <<"ok"<< Sigma_mat<<endl;
      
    mu_vec=xmat*beta_vec;
    mat gamma_mat=zmat*w_mat.t(); 
    
    vec mm=dfuns.Dfun(Ymat.col(1), alpha_vec, mu_vec, gamma_mat); 
    vec zz=(Ymat.col(2)-mm)/sigma;
    
    //Rcout <<"ok "<< size(gamma_mat)<<endl;
    
    double term1=as_scalar(sum(log(dnor_vec(zz))-log(sigma)));
    
    mat smat_inv=inv(Sigma_mat);
    mat mat1=w_mat.t()*smat_inv*w_mat;
    vec bb=mat1.diag();
    
    //Rcout<<omat.submat(0, 0, 10, 10)<<endl;
  
    double bb1=log(det(Sigma_mat));
    
    //Rcout << "bb1="<<bb1<<endl;
    
    double term2=(-0.5)*n_rand*bb1-0.5*as_scalar(sum(bb));
        
    double term3=as_scalar(sum(log(sigmak)-log(hcau_par*hcau_par+sigmak%sigmak)));

    double term4=0.0;
    
    if(n_gamma>1)
    {
      vec rho_vec=pars.subvec(idx1+n_gamma, n_pars-1);
      rho_vec=2*exp(rho_vec)/(1+exp(rho_vec))-1;
   
      term4=as_scalar(sum(log(1-rho_vec%rho_vec)));
    }
    
    res=term1+term2+term3+term4;
    
    //Rcout<<"term1="<<term1<<" term2="<<term2<<" term3="<<term3<<" term4="<<term4<<endl;
    //Rcout<<"loglik="<<res<<endl;
  //} 
    
  return(res);

}


double log_post_ratio_rmdt(double yy, int j, const vec& th, const List& dat_list)
{
  double res, aa0, aa1;
  
  vec th1=th;
  th1(j)=th(j)+yy;
  
  aa0=gen_loglik_rmdt(th, dat_list);
  aa1=gen_loglik_rmdt(th1, dat_list);
  
  res=aa1-aa0;
  
  //Rcout<<"j="<<j<<" yy="<<yy<<" aa1="<<aa1<<"   aa0="<<aa0<<endl;
    
  return res;
}


double log_post_ratio_gomp(double yy, int j, const vec& th, const List& dat_list)
{
  double res, aa0, aa1;
  
  vec th1=th;
  th1(j)=th(j)+yy;
  
  aa0=loglik_recurr_gomp(th, dat_list);
  aa1=loglik_recurr_gomp(th1, dat_list);
  
  res=aa1-aa0;
  
  //Rcout<<"j="<<j<<" yy="<<yy<<" aa1="<<aa1<<"   aa0="<<aa0<<endl;
    
  return res;
}

double log_post_ratio_wei_altp(double yy, int j, const vec& th, const List& dat_list)
{
  double res, aa0, aa1;
  
  vec th1=th;
  th1(j)=th(j)+yy;
  
  aa0=loglik_recurr_wei_altp(th, dat_list);
  aa1=loglik_recurr_wei_altp(th1, dat_list);
  
  res=aa1-aa0;
  
  //Rcout<<"j="<<j<<" yy="<<yy<<" aa1="<<aa1<<"   aa0="<<aa0<<endl;
    
  return res;
}

 
//[[Rcpp::export]]       
List Gibbs_lifedata(vec vec_theta, List dat_list, vec vec_scale, int int_iter)
{
  gibbsModel gibbs_model;
  gibbs_model.dat_list=dat_list;
  gibbs_model.log_post_ratio=&log_post_ratio_lifedata;
      
  List res=Gibbs_adaptive(gibbs_model, vec_theta, vec_scale, int_iter, 500, 400); //500, 400
  return res;
}

//[[Rcpp::export]]       
List Gibbs_rmdt(vec vec_theta, List dat_list, vec vec_scale, int int_iter)
{
  gibbsModel gibbs_model;
  gibbs_model.dat_list=dat_list;
  gibbs_model.log_post_ratio=&log_post_ratio_rmdt;
      
  List res=Gibbs_adaptive(gibbs_model, vec_theta, vec_scale, int_iter, 100, 500); //500, 400 
  return res;
}


//[[Rcpp::export]]       
List Gibbs_gomp(vec vec_theta, List dat_list, vec vec_scale, int int_iter)
{
  gibbsModel gibbs_model;
  gibbs_model.dat_list=dat_list;
  gibbs_model.log_post_ratio=&log_post_ratio_gomp;
      
  List res=Gibbs_adaptive(gibbs_model, vec_theta, vec_scale, int_iter, 500, 400);   
  return res;
}

//[[Rcpp::export]]       
List Gibbs_wei_altp(vec vec_theta, List dat_list, vec vec_scale, int int_iter)
{
  gibbsModel gibbs_model;
  gibbs_model.dat_list=dat_list;
  gibbs_model.log_post_ratio=&log_post_ratio_wei_altp;
      
  List res=Gibbs_adaptive(gibbs_model, vec_theta, vec_scale, int_iter, 500, 400);   
  return res;
}


vec loglik_recurr_gomp_der(const vec& pars, const List& dat)
{
  vec res(3, fill::zeros);
  
  mat ne_mat=dat["Non_event_mat"];
  mat emat=dat["Event_mat"];
  
  double alpha=exp(pars(0));
  double beta=exp(pars(1))/(1+exp(pars(1))); 
  double kappa=exp(pars(2))/(1+exp(pars(2)));  
  
  vec ett=emat.col(2); 
  vec ewts=emat.col(4);
  
  vec beta_ett=pow_vec(beta, ett);
  
  double pd_11=sum(ewts);
  double pd_12=sum(ewts%(ett+ett%beta_ett*log(kappa)+1/log(beta))*(1-beta));
  double pd_13=sum(ewts%(beta_ett+1/log(kappa))*(1-kappa));  

  vec ne_tau1=ne_mat.col(1);
  vec ne_tau2=ne_mat.col(2);
  vec ne_xt=ne_mat.col(3);
  vec ne_wts=ne_mat.col(4);
  
  vec xt_wts=ne_xt%ne_wts;
  
  vec beta_tau1=pow_vec(beta, ne_tau1);  
  vec beta_tau2=pow_vec(beta, ne_tau2);

  vec kappa_beta_tau1=pow_vec(kappa, beta_tau1);  
  vec kappa_beta_tau2=pow_vec(kappa, beta_tau2);   

  double pd_21=sum(xt_wts%(kappa_beta_tau2-kappa_beta_tau1)*alpha);
  double pd_22=sum(xt_wts%(kappa_beta_tau2%ne_tau2%beta_tau2-kappa_beta_tau1%ne_tau1%beta_tau1)*log(kappa)*alpha*(1-beta));
  double pd_23=sum(xt_wts%(kappa_beta_tau2%beta_tau2-kappa_beta_tau1%beta_tau1)*alpha*(1-kappa));

  res(0)=pd_11-pd_21;
  res(1)=pd_12-pd_22;
  res(2)=pd_13-pd_23;  

  return(res);
}

vec loglik_recurr_wei_altp_der(const vec& pars, const List& dat)
{
  vec res(3, fill::zeros);
  
  mat ne_mat=dat["Non_event_mat"];
  mat emat=dat["Event_mat"];
  
  double alpha=exp(pars(0));
  double mu=pars(1); 
  double sigma=exp(pars(2));  
  
  vec ett=emat.col(2); 
  vec ewts=emat.col(4);
  vec z_ett=(log(ett)-mu)/sigma;
    
  double pd_11=sum(ewts);
  double pd_12=sum(ewts%(1-exp(z_ett))/sigma);
  double pd_13=sum(ewts%(-1+(1-exp(z_ett))%z_ett));  

  vec ne_tau1=ne_mat.col(1);
  vec ne_tau2=ne_mat.col(2);
  vec ne_xt=ne_mat.col(3);
  vec ne_wts=ne_mat.col(4);
  
  vec xt_wts=ne_xt%ne_wts;
  
  vec z_tau1=(log(ne_tau1)-mu)/sigma;  
  vec z_tau2=(log(ne_tau2)-mu)/sigma;  
  
  vec tmp1=psev_vec(z_tau1);
  vec tmp2=psev_vec(z_tau2);

  double pd_21=sum(xt_wts%(psev_vec(z_tau2)-psev_vec(z_tau1))*alpha);
  double pd_22=sum(xt_wts%(tmp2-tmp1)*alpha*(-1)/sigma);
  double pd_23=sum(xt_wts%(tmp2%z_tau2-tmp1%z_tau1)*(-1)*alpha);

  res(0)=pd_11-pd_21;
  res(1)=pd_12-pd_22;
  res(2)=pd_13-pd_23;  
  
  return(res);
}


double loglik_recurr_gomp_joint(const vec& vec_theta, const vec& vec_r, const vec& M_diag, const List& dat_list)
{
   double num_val;
   num_val = loglik_recurr_gomp(vec_theta, dat_list)-0.5 * sum(vec_r % vec_r / M_diag);
   return num_val;
}

double loglik_recurr_wei_altp_joint(const vec& vec_theta, const vec& vec_r, const vec& M_diag, const List& dat_list)
{
   double num_val;
   num_val = loglik_recurr_wei_altp(vec_theta, dat_list)-0.5 * sum(vec_r % vec_r / M_diag);
   return num_val;
}


//[[Rcpp::export]]       
mat Nuts_gomp(vec vec_theta, int int_iter, List dat_list, List nuts_ctrl_list)
{
  nutsModel nuts_model;
  nuts_model.fn_val=&loglik_recurr_gomp;  
  nuts_model.gn_val=&loglik_recurr_gomp_der;
  nuts_model.joint_fn_val=&loglik_recurr_gomp_joint;
    
  vec vec_M_diag=dat_list["vec_M_diag"];
  //read control parameters   
  int int_M_adapt=nuts_ctrl_list["int_M_adapt"];
  double num_delta=nuts_ctrl_list["num_delta"];
  int int_max_treedepth=nuts_ctrl_list["int_max_treedepth"];
  double num_eps=nuts_ctrl_list["num_eps"]; 
   
  //double tmp=gen_loglik_lifedata(vec_theta, dat_list); 
  
  //Rcout<<"Ok here"<<endl;
      
  return Nuts_cpp(nuts_model, dat_list, vec_theta, int_iter, vec_M_diag, int_M_adapt, num_delta, int_max_treedepth, num_eps);
   
  //mat mat1;
  //mat1.zeros();
  
  //return mat1;

  
}


//[[Rcpp::export]]       
mat Nuts_wei_altp(vec vec_theta, int int_iter, List dat_list, List nuts_ctrl_list)
{
  nutsModel nuts_model;
  
  nuts_model.fn_val=&loglik_recurr_wei_altp;  
  nuts_model.gn_val=&loglik_recurr_wei_altp_der;
  nuts_model.joint_fn_val=&loglik_recurr_wei_altp_joint;
    
  vec vec_M_diag=dat_list["vec_M_diag"];
  //read control parameters   
  int int_M_adapt=nuts_ctrl_list["int_M_adapt"];
  double num_delta=nuts_ctrl_list["num_delta"];
  int int_max_treedepth=nuts_ctrl_list["int_max_treedepth"];
  double num_eps=nuts_ctrl_list["num_eps"]; 
   
      
  return Nuts_cpp(nuts_model, dat_list, vec_theta, int_iter, vec_M_diag, int_M_adapt, num_delta, int_max_treedepth, num_eps);
  
}

