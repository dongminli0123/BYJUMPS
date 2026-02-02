simul <- function(seed, sig, time_horizon, h, shift)
{

library(ggplot2)
library(reshape2)
library(Rcpp)
library(mvtnorm)

t1 <- Sys.time()
#config <- tf$ConfigProto(intra_op_parallelism_threads = 4L,
#                         inter_op_parallelism_threads = 4L)
#session = tf$Session(config = config)
#k_set_session(session)

#Sys.setenv(OMP_NUM_THREADS="4")
#parallel::mcaffinity(1:4)
#cl <- makePSOCKcluster(detectCores() - 28)
#registerDoParallel(cl)

#dev.off()

# HiPerGatro environment
Sys.setenv("PKG_CXXFLAGS"="-I/apps/cplex/12.7/cplex/include/ -I/apps/cplex/12.7/concert/include/ -DIL_STD -O2 -DNDEBUG -s")
Sys.setenv("PKG_LIBS" = "-L/apps/cplex/12.7/cplex/lib/x86-64_linux/static_pic/ -lilocplex -lcplex -L/apps/cplex/12.7/concert/lib/x86-64_linux/static_pic/ -lconcert")
sourceCpp(code='
#pragma warning(disable : 4996)
#include <vector>
          #include <ilcplex/ilocplex.h>
          #include <utility>
          #include <algorithm>
          #include <math.h>
          #include <string>
          #include <iostream>
          #include <fstream>
          #include <sstream>
          #include <tuple>  
          #include <Rcpp.h> 
          using namespace Rcpp;
          
          ILOSTLBEGIN
          typedef IloArray<IloNumVarArray> IloNumVarArray2;
          typedef IloArray<IloNumVarArray2> IloNumVarArray3;
          typedef IloArray<IloNumVarArray3> IloNumVarArray4;
          
          struct status
          {
          	int grid_wide;
          	int grid_height;
          	int num_uav;
          
          	int time_horizon;
          
          	double sigma;
          
          
          	std::vector<std::vector<double> >  e;
          	std::vector<std::vector<double> >  v;
          	std::vector<std::vector<double> >  sigma_0;
          	std::vector<std::vector<double> > K;
          
          	double M;
          	double C;
          
          
          	std::vector<int> current_loc_x;
          	std::vector<int> current_loc_y;
          
          	std::vector<std::vector<std::pair<int, int> > >  route;
          
          	//default initializer
          	status() :
          		grid_wide(0),
          		grid_height(0),
          		num_uav(0),
          		time_horizon(0),
          		sigma(0),
          		M(0),
          		//C(0),
          		e(std::vector<std::vector<double> >()),
          		v(std::vector<std::vector<double> >()),
          		sigma_0(std::vector<std::vector<double> >()),
          		K(std::vector<std::vector<double> >()),
          		current_loc_x(std::vector<int>()),
          		current_loc_y(std::vector<int>()),
          		route(std::vector<std::vector<std::pair<int, int> > >())
          	{}
          
          };
          
          
          void optimize_route_grid(status& in_current_status);
          
          // [[Rcpp::export]]
          
          
          NumericMatrix action(int in_grid_wide, int in_grid_height, int in_num_uav, int in_time_horizon, double in_sigma,  double in_C,  NumericMatrix in_K, NumericMatrix sigma_0, NumericMatrix e, NumericMatrix v, NumericMatrix in_current_uav_location)
          {
          
          	status input;
          
          	//basic setup
          	input.grid_wide = in_grid_wide;
          	input.grid_height = in_grid_height;
          	input.num_uav = in_num_uav;
          
          
          	input.time_horizon = in_time_horizon;
          	input.sigma = in_sigma;
          	
          	input.M = 100;
          	input.C = in_C;
          
          	for (int i = 0; i < input.grid_wide; i++)
          	{
          		std::vector<double> temp1;
          		std::vector<double> temp2;
          		std::vector<double> temp3;
          		std::vector<double> temp4;
          
          		for (int j = 0; j < input.grid_height; j++)
          		{
          			temp1.push_back(sigma_0(i,j));
          			temp2.push_back(e(i,j));
          			temp3.push_back(v(i,j));
          			temp4.push_back(in_K(i,j));
          		}
          		input.sigma_0.push_back(temp1);
          		input.e.push_back(temp2);
          		input.v.push_back(temp3);
          		input.K.push_back(temp4);
          
          		temp1.clear();
          		temp2.clear();
          		temp3.clear();
          		temp4.clear();
          	}
          
          
          	//current location of each UAV
          	//current_loc_x[u] correponds to current x coordinate of uav u
          	//current_loc_y[u] correponds to current y coordinate of uav u
          	for (int u = 0; u < input.num_uav; u++)
          	{
          		input.current_loc_x.push_back(in_current_uav_location(u,0));
          		input.current_loc_y.push_back(in_current_uav_location(u,1));
          	}
          
          
          	optimize_route_grid(input);
          
          	NumericMatrix output(input.num_uav, 2);
          
          
          	for (int u = 0; u < input.num_uav; u++)
          	{
          		output(u, 0) = input.route[u][1].first;
          		output(u, 1) = input.route[u][1].second;
          	}
          
          	return output;
          }
          
          
          
          void optimize_route_grid(status& in_current_status)
          {
          	
          	IloEnv env;
          
          
          
          	try {
          		IloModel model(env);
          		IloCplex cplex(model);
          
          
          
          		double M = in_current_status.M;
          		double C = in_current_status.C;
          
          
          
          		IloNumVarArray3 sigma(env, in_current_status.time_horizon);
          		IloNumVarArray3 v(env, in_current_status.time_horizon);
          		IloNumVarArray3 l(env, in_current_status.time_horizon);
          		IloNumVarArray3 l_bar(env, in_current_status.time_horizon);
          
          
          		for (int t = 0; t < in_current_status.time_horizon; t++)
          		{
          			sigma[t] = IloNumVarArray2(env, in_current_status.grid_wide);
          			v[t] = IloNumVarArray2(env, in_current_status.grid_wide);
          			l[t] = IloNumVarArray2(env, in_current_status.grid_wide);
          			l_bar[t] = IloNumVarArray2(env, in_current_status.grid_wide);
          
          			for (int i = 0; i < in_current_status.grid_wide; i++)
          			{
          				sigma[t][i] = IloNumVarArray(env, in_current_status.grid_height, 0, IloInfinity, ILOFLOAT);
          				v[t][i] = IloNumVarArray(env, in_current_status.grid_height, 0, IloInfinity, ILOFLOAT);
          				l[t][i] = IloNumVarArray(env, in_current_status.grid_height, -IloInfinity, IloInfinity, ILOFLOAT);
          				l_bar[t][i] = IloNumVarArray(env, in_current_status.grid_height, -IloInfinity, IloInfinity, ILOFLOAT);
          			}
          		}
          
          
          
          		IloNumVarArray4 x(env, in_current_status.time_horizon);
          
          
          		for (int t = 0; t < in_current_status.time_horizon; t++)
          		{
          			x[t] = IloNumVarArray3(env, in_current_status.num_uav);
          
          			for (int u = 0; u < in_current_status.num_uav; u++)
          			{
          				x[t][u] = IloNumVarArray2(env, in_current_status.grid_wide);
          
          				for (int i = 0; i < in_current_status.grid_wide; i++)
          				{
          					x[t][u][i] = IloNumVarArray(env, in_current_status.grid_height, 0, 1, ILOBOOL);
          
          				}
          			}
          		}
          
          		IloNumVarArray3 y(env, in_current_status.time_horizon);
          		for (int t = 0; t < in_current_status.time_horizon; t++)
          		{
          
          			y[t] = IloNumVarArray2(env, in_current_status.grid_wide);
          
          			for (int i = 0; i < in_current_status.grid_wide; i++)
          			{
          				y[t][i] = IloNumVarArray(env, in_current_status.grid_height, 0, 1, ILOBOOL);
          
          			}
          
          		}
          
          
          		IloExpr Objective(env);
          
          
          		for (int t = 1; t < in_current_status.time_horizon; t++)
          		{
          			for (int i = 0; i < in_current_status.grid_wide; i++)
          			{
          				for (int j = 0; j < in_current_status.grid_height; j++)
          				{
          					Objective += l_bar[t][i][j];
          					
          				}
          			}
          		}
          
          
          		model.add(IloMaximize(env, Objective));
          
          		Objective.end();
          
          
          
          		
          		IloConstraintArray con1(env);
          
          		for (int t = 0; t < in_current_status.time_horizon - 1; t++)
          		{
          			for (int u = 0; u < in_current_status.num_uav; u++)
          			{
          
          				for (int i = 0; i < in_current_status.grid_wide; i++)
          				{
          					for (int j = 0; j < in_current_status.grid_height; j++)
          					{
          						IloExpr sum(env);
          
          						//con13
          						sum += x[t + 1][u][i][j];
          
          						if (i + 1 < in_current_status.grid_wide)
          						{
          							sum += x[t + 1][u][i + 1][j];
          
          							if (j + 1 < in_current_status.grid_height)
          							{
          								sum += x[t + 1][u][i + 1][j + 1];
          							}
          
          							if (j - 1 >= 0)
          							{
          								sum += x[t + 1][u][i + 1][j - 1];
          							}
          
          						}
          
          						if (i - 1 >= 0)
          						{
          							sum += x[t + 1][u][i - 1][j];
          
          							if (j + 1 < in_current_status.grid_height)
          							{
          								sum += x[t + 1][u][i - 1][j + 1];
          							}
          
          							if (j - 1 >= 0)
          							{
          								sum += x[t + 1][u][i - 1][j - 1];
          							}
          
          						}
          
          						if (j + 1 < in_current_status.grid_height)
          						{
          							sum += x[t + 1][u][i][j + 1];
          						}
          
          						if (j - 1 >= 0)
          						{
          							sum += x[t + 1][u][i][j - 1];
          						}
          
          
          
          						con1.add(sum >= x[t][u][i][j]);
          						sum.end();
          					}
          				}
          			}
          		}
          
          		model.add(con1);
          		con1.end();
          
          
          
          		IloConstraintArray con3(env);
          
          		for (int t = 1; t < in_current_status.time_horizon; t++)
          		{
          			for (int u = 0; u < in_current_status.num_uav; u++)
          			{
          
          				IloExpr sum(env);
          
          				for (int i = 0; i < in_current_status.grid_wide; i++)
          				{
          					for (int j = 0; j < in_current_status.grid_height; j++)
          					{
          						sum += x[t][u][i][j];
          					}
          				}
          				con3.add(sum == 1);
          				sum.end();
          			}
          		}
          
          		model.add(con3);
          		con3.end();
          
          
          
          
          		//constraints 15 in the model
          
          		IloConstraintArray con33(env);
          
          		for (int t = 1; t < in_current_status.time_horizon; t++)
          		{
          			for (int i = 0; i < in_current_status.grid_wide; i++)
          			{
          				for (int j = 0; j < in_current_status.grid_height; j++)
          				{
          					IloExpr sum(env);
          
          					for (int u = 0; u < in_current_status.num_uav; u++)
          					{
          
          						sum += x[t][u][i][j];
          
          						con33.add(x[t][u][i][j] <= y[t][i][j]);
          
          					}
          
          					con33.add(sum >= y[t][i][j]);
          					
          					sum.end();
          				}
          			}
          		}
          
          		model.add(con33);
          		con33.end();
          
          
          
          	
          		IloConstraintArray con4(env);
          
          		for (int t = 1; t < in_current_status.time_horizon; t++)
          		{
          
          			for (int i = 0; i < in_current_status.grid_wide; i++)
          			{
          				for (int j = 0; j < in_current_status.grid_height; j++)
          				{
          
          
          					con4.add(sigma[t][i][j] <= in_current_status.K[i][j] * (in_current_status.sigma + sigma[t - 1][i][j]) + M * (1 - y[t][i][j]));
          					con4.add(sigma[t][i][j] <= in_current_status.sigma + sigma[t - 1][i][j]);
          
          
          				}
          			}
          		}
          
          
          
          
          		model.add(con4);
          		con4.end();
          
          
          
          
          
          		IloConstraintArray con44(env);
          
          		for (int t = 1; t < in_current_status.time_horizon; t++)
          		{
          
          			for (int i = 0; i < in_current_status.grid_wide; i++)
          			{
          				for (int j = 0; j < in_current_status.grid_height; j++)
          				{
          
          					con44.add(v[t][i][j] == in_current_status.v[i][j] - in_current_status.sigma_0[i][j] + sigma[t][i][j]);
          					con44.add((l[t][i][j] - in_current_status.e[i][j]) * (l[t][i][j] - in_current_status.e[i][j]) <= C * C * v[t][i][j]);
          					con44.add(l[t][i][j] >= in_current_status.e[i][j]);
          				}
          			}
          		}
          
          		model.add(con44);
          		con44.end();
          
          
          
          
          	
          		IloConstraintArray con5(env);
          
          		for (int t = 1; t < in_current_status.time_horizon; t++)
          		{
          			for (int i = 0; i < in_current_status.grid_wide; i++)
          			{
          				for (int j = 0; j < in_current_status.grid_height; j++)
          				{
          
          					con5.add(l_bar[t][i][j] <= M * y[t][i][j]);
          					con5.add(l_bar[t][i][j] <= l[t][i][j] + M * (1 - y[t][i][j]));
          				}
          			}
          
          		}
          
          		model.add(con5);
          		con5.end();
          
          
          		//initialize sigma and x at time 0
          		IloConstraintArray con6(env);
          
          		//initialize initial location for all UAV
          
          		for (int u = 0; u < in_current_status.num_uav; u++)
          		{
          			con6.add(x[0][u][in_current_status.current_loc_x[u]][in_current_status.current_loc_y[u]] == 1);
          		}
          
          		for (int i = 0; i < in_current_status.grid_wide; i++)
          		{
          			for (int j = 0; j < in_current_status.grid_height; j++)
          			{
          				//initialize sigma
          				con6.add(v[0][i][j] == in_current_status.v[i][j]);
          				con6.add(sigma[0][i][j] == in_current_status.sigma_0[i][j]);
          			}
          		}
          
          		model.add(con6);
          		con6.end();
          
          
          
          		//break symmetry
          		IloConstraintArray con7(env);
          
          		for (int u = 1; u < in_current_status.num_uav; u++)
          		{
          			IloExpr sum(env);
          			for (int i = 0; i < in_current_status.grid_wide; i++)
          			{
          				for (int j = 0; j < in_current_status.grid_height; j++)
          				{
          
          					sum += (x[1][u][i][j] - x[1][u - 1][i][j]) * (i * in_current_status.grid_height + j);
          				}
          			}
          			con7.add(sum >= 0);
          			sum.end();
          		}
          
          		model.add(con7);
          		con7.end();
          
          
          
          	
          		cplex.extract(model);
          
          		cplex.setOut(env.getNullStream());
          		cplex.setWarning(env.getNullStream());
          
          		cplex.setParam(IloCplex::ClockType, 1);
          		cplex.setParam(IloCplex::Threads, 4);
          
          		cplex.setParam(IloCplex::MIPEmphasis, 3);
          		cplex.setParam(IloCplex::EpGap, 5e-2);
          		cplex.setParam(IloCplex::TiLim, 60);
             
          
          
          		double start_time = cplex.getCplexTime();
          
          
          		if (!cplex.solve())
          		{
          			env.error() << "Failed to optimize LP" << std::endl;
          			throw(-1);
          		}
          		//std::cout << "-------------RESULTS---------------" << std::endl;
          		//IloAlgorithm::Status solStatus = cplex.getStatus();
          		//std::cout << "Solution status: " << solStatus << " ";
          
          
          		double finish_time = cplex.getCplexTime();
          
          		double run_time = finish_time - start_time;
          
          		//std::cout<< "CPU time " << run_time << std::endl;
          
          
          		in_current_status.route.clear();
          
          
          		//std::cout << "routes" << std::endl;
          
          		for (int u = 0; u < in_current_status.num_uav; u++)
          		{
          			//std::cout << "UAV #" << u << " ";
          
          			std::vector<std::pair<int, int> > temp;
          
          			for (int t = 0; t < in_current_status.time_horizon; t++)
          			{
          				bool flag = true;
          				for (int i = 0; flag == true && i < in_current_status.grid_wide; i++)
          				{
          					for (int j = 0; j < in_current_status.grid_height; j++)
          					{
          						if (int(round(cplex.getValue(x[t][u][i][j]))) == 1)
          						{
          							//std::cout << "(" << i << "," << j << ") ";
          							temp.push_back(std::make_pair(i, j));
          							flag = false;
          							break;
          						}
          					}
          				}
          
          			}
          			in_current_status.route.push_back(temp);
          			temp.clear();
          
          			//std::cout << std::endl;
 
          		}
          
          
          	}
          	catch (IloException& e) {
          		std::cerr << "Concert exception caught: " << e << std::endl;
          	}
          	catch (...) {
          		std::cerr << "Unknown exception caught" << std::endl;
          	}
          
          	env.end();
          
          
          }
')


# R starts
set.seed(seed)

grid_m = 10
grid_n = 10
d = grid_m * grid_n   #dimension
q = 2   # nuav

## jump model
## spatial structure
distance = matrix(0, d, d)
for (i in 1:d)
{
  y_i=ceiling(i/grid_m)      #y label
  x_i=i - (y_i-1) * grid_m   #x label
  for (j in i:d)
  {
    y_j=ceiling(j/grid_m)
    x_j=j - (y_j-1) * grid_m
    distance[i,j]=sqrt((x_i - x_j)^2 + (y_i - y_j)^2)
    distance[j,i]=sqrt((x_i - x_j)^2 + (y_i - y_j)^2)
  }
  
}
sigma_cov = 1
covar = (sigma_cov^2) * exp(- distance^2/50)
## tapering
r=3
c_taper=pmax(1-distance/r,0)^2 * (1+distance/(2*r))
Cov_taper = covar * c_taper

## prior 
tao = Cov_taper/2  #tao square  x   diagonal value
sigma0 =  Cov_taper/2    #sigma0 square  diagonal value
mu0 = matrix(0, d, 1)  #thet0
delta = matrix(0.5, d, 1) #jump size
p = 0.95
sigma = matrix(0,d,d) + diag(sig,d)  #jump variance =0 equals bayes

t = 1000  #########total time
## m,v,K
thet = matrix(0,d,2^(1))
sigma_n = array(0,dim = c(d,d,t))
m_thet_cut = matrix(0,t,d)   # each row represents m_thet at time t
cov_thet_cut = array(0,dim = c(d,d,t))
v_thet_cut = matrix(0,t,d)   # each row represents v_thet at time t
k_all = matrix(0,t,d)   # each row represents k at time t
## exact lr
l0_ind_cut = matrix(0,t,d)   # each row represents pdf(marginals) of IC at time t
l1_ind_cut = matrix(0,t,d)   # each row represents pdf(marginals) of OC at time t
lr_ind_cut = matrix(0,t,d)    # each row represents likelihood ratio(marginals) at time t
h1 = 1
h0 = 0
mon_lr = rep(0,t)

index = matrix(0,t,q,byrow = T)

n=1   #first obs
mu = mu0
OC_loc = c(45,46,55)
x = rmvnorm(1,mu,tao + sigma0)
layout =  matrix(c(0,grid_m-1,0,grid_m-1), 2,2)     ### initial layout for 2 uavs
index[1,] = layout[,2] * grid_m + layout[,1] + 1   
x_obs = as.matrix(x[index[1,]])
E = matrix(0,q,d)
for (i in 1:q){
  E[i,index[1,i]] = 1
}
tao_obs = tao[index[1,],index[1,]]
inv_tao_obs = solve(tao_obs)
inv_sigman = solve(sigma0 + sigma)
sigma_n[,,1] = solve(t(E) %*% inv_tao_obs %*% E + inv_sigman)
thet[,1] = sigma_n[,,1] %*% (inv_sigman %*% mu0 + t(E) %*% inv_tao_obs %*% x_obs)
thet[,2] = sigma_n[,,1] %*% (inv_sigman %*% (mu0 + delta) + t(E) %*% inv_tao_obs %*% x_obs)
m0 = dmvnorm(t(x_obs), E %*% mu0, E %*% (tao + sigma + sigma0) %*% t(E))
m1 = dmvnorm(t(x_obs), E %*% (mu0 + delta), E %*% (tao + sigma + sigma0) %*% t(E))

a = rep(0,2)
a_h = 10^(-5)
a[1] = p * m0 / (p * m0 + (1 - p) * m1)
a[2] = 1 -  a[1]
a_cut = list(a)

## m,v,k
m_thet_cut[1,] = a_cut[[1]][1] * thet[,1] + a_cut[[1]][2] * thet[,2]
cov_thet_cut[,,n] = a_cut[[1]][1] * thet[,1] %*% t(thet[,1]) + a_cut[[1]][2] * thet[,2] %*% t(thet[,2]) + sigma_n[,,1] - m_thet_cut[n,] %*% t(m_thet_cut[n,])
v_thet_cut[n,] = as.matrix(diag(cov_thet_cut[,,n]))
## exact lr
l1_ind_cut[1,] = a_cut[[1]][1] * dnorm(h1,thet[,1],sqrt(diag(sigma_n[,,n]))) + a_cut[[1]][2] * dnorm(h1,thet[,2],sqrt(diag(sigma_n[,,n])))
l0_ind_cut[1,] = a_cut[[1]][1] * dnorm(h0,thet[,1],sqrt(diag(sigma_n[,,n]))) + a_cut[[1]][2] * dnorm(h0,thet[,2],sqrt(diag(sigma_n[,,n])))
lr_ind_cut[1,]=log(l1_ind_cut[1,]) - log(l0_ind_cut[1,])
## detection stat
mon_lr[1] = sum(head(sort(lr_ind_cut[n,],decreasing=T),q))

c = 0.2/sig     # tao=0.2
K = 1 + 1/(2*c) - sqrt(1/(4*c^2) + 1/c)

thet_cut = thet
n_state = ncol(thet)

for (n in 2:t)
{
  
  #print(n)
  layout1 = layout
  layout = action(grid_m, grid_n, q, time_horizon, sig, 2,  matrix(K,grid_m, grid_n), matrix(diag(sigma_n[,,n-1]), grid_m), matrix(m_thet_cut[n-1,], grid_m), matrix(v_thet_cut[n-1,], grid_m), layout1)
  index[n,] = layout[,2] * grid_m + layout[,1] + 1
  #print(index[n,])
  
  if (n>50){   ######## OC time
     mu[OC_loc] = mu0[OC_loc] + shift
  }
  x = rmvnorm(1,mu,tao + sigma0)
  x_obs = as.matrix(x[index[n,]])
  E = matrix(0,q,d)
  for (i in 1:q){
    E[i,index[n,i]] = 1
  }
  tao_obs = tao[index[n,],index[n,]]
  inv_tao_obs = solve(tao_obs)
  inv_sigman = solve(sigma_n[,,n-1] + sigma)
  sigma_n[,,n] = solve(t(E) %*% inv_tao_obs %*% E + inv_sigman)
  
  thet = matrix(0,d,2*n_state)
  thet[,seq(1,2*n_state-1,by=2)] = sigma_n[,,n] %*% (inv_sigman %*% thet_cut + as.vector(t(E) %*% inv_tao_obs %*% x_obs))
  thet[,seq(2,2*n_state,by=2)] = thet[,seq(1,2*n_state-1,by=2)] + as.vector(sigma_n[,,n] %*% inv_sigman %*% delta)
  
  m0 = matrix(0,n_state,1)
  m1 = matrix(0,n_state,1)
  for (i in 1:n_state){
    v = E %*% (tao + sigma + sigma_n[,,n-1]) %*% t(E)
    m0[i,1] = dmvnorm(t(x_obs), E %*% thet_cut[,i], v)
    m1[i,1] = dmvnorm(t(x_obs), E %*% (thet_cut[,i] + as.vector(delta)), v)
  }
  
  NC_cut = sum( p * a_cut[[n-1]] * m0) + sum( (1 - p) * a_cut[[n-1]] * m1)
  a = rep(0,2*n_state)
  a[seq(1,2*n_state-1,by=2)] = p * a_cut[[n-1]] * m0 / NC_cut
  a[seq(2,2*n_state,by=2)] = (1 - p) * a_cut[[n-1]]  * m1 / NC_cut
  a_cut[[n]] = a[a>a_h]
  a_cut[[n]] = a_cut[[n]] / sum(a_cut[[n]])   ##standardize
  thet_cut = as.matrix(thet[,a>a_h])
  n_state = ncol(thet_cut)
  
  for (j in 1:n_state){
    m_thet_cut[n,] = m_thet_cut[n,] + a_cut[[n]][j] * thet_cut[,j]
    cov_thet_cut[,,n] = cov_thet_cut[,,n] + a_cut[[n]][j] * thet_cut[,j] %*% t(thet_cut[,j])
    l1_ind_cut[n,] = l1_ind_cut[n,] + a_cut[[n]][j] * dnorm(h1,thet_cut[,j],sqrt(diag(sigma_n[,,n])))
    l0_ind_cut[n,] = l0_ind_cut[n,] + a_cut[[n]][j] * dnorm(h0,thet_cut[,j],sqrt(diag(sigma_n[,,n])))
  
  }
  
  cov_thet_cut[,,n] = cov_thet_cut[,,n] + sigma_n[,,n] - (m_thet_cut[n,]) %*% t(m_thet_cut[n,])
  v_thet_cut[n,] = as.matrix(diag(cov_thet_cut[,,n]))
  lr_ind_cut[n,] = log(l1_ind_cut[n,]) - log(l0_ind_cut[n,])
  mon_lr[n] = sum(head(sort(lr_ind_cut[n,],decreasing=T),q))
  
  if( mon_lr[n] > h)
  {
    break
  }

}

#save(mon_lr, file = paste0("OCmon_", seed, "_", sig, "_", time_horizon,".RData"))
return(n - 50)
}


# Based on 500 simulation runs, outputs the magnitude of mean shift (shift), the value of ARL (ARL), and the standard error (std) of BYJUMPS under scenario 3C
sig_3C = c(0.0025, 0.01, 0.03, 0.03)
time_horizon_3C = c(2,5,3,5)
h_3C = c(11.232, 10.255, 9.232, 9.197)    
for (shift in 1:4)
{
  sig = sig_3C[shift]
  time_horizon = time_horizon_3C[shift]
  h = h_3C[shift]
  arl = NULL
  for (seed in 1:500)
  {
    l = simul(seed, sig, time_horizon, h, shift)
    #print(l)
    arl = c(arl, l)
  }
  average = mean(arl)
  std = sd(arl)/sqrt(length(arl))
  print(paste0("3C:shift = ", shift, " ARL = ", average, " std = ", std))
}