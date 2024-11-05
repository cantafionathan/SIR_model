% ADAPTED FROM MODEL IN THIS PAPER
% Turkyilmazoglu M. An extended epidemic model with vaccination: Weak-immune SIRVI. 
% Physica A. 2022;598:127429. doi:10.1016/j.physa.2022.127429

N = 1000000; % total population
I0 = 0.01*N; % number of people sick once the health authority notices an issue
S = N - I0; % total number of susceptible individuals
I = I0; % total number of infected individuals
V = 0; % total number of vaccinated individuals
R = 0; % total number of recovered individuals

% response: max(i)
% treatments: can vary vaccination rate, isolation_factor
% another idea (quarantining sympotomatic infecteds? - have prob that
% infected is asymptomatic or not, then quarantine some portion of those that showing)
% blocking factors???

%probability that an infectious individual is sympomatic
prob_symptomatic = 0.4; % covid
% treatment we can change:
quarantine = 0; % quarantining sympomatic individuals, yes or no

base_num_inter = 20; % number of daily interactions if there was no policies
% treatment we can change:
isolation_factor = 0; % social isolation policy ranging from zero isolation to total isolation
daily_interactions = base_num_inter*(1-isolation_factor)*(1-quarantine*prob_symptomatic); % average number of daily interactions a susceptible individual has
probability_of_spread_unvaccinated = 0.1; % probabilty that an unvaccinated individual gets infected from interacting with an infected individual
probability_of_spread_vaccinated = 0.01*probability_of_spread_unvaccinated; % probability that a vaccinated indidual gets infected from interacting with an infected individual

% beta = transmission rate
beta_unvac = daily_interactions*probability_of_spread_unvaccinated;
beta_vac = daily_interactions*probability_of_spread_vaccinated;

% gamma = recovery rate
gamma = 1/14; % covid % recovery in 2 weeks????
%gamma = 1/31;

% treatment we can change:
% nu = vaccination rate of susceptible individuals
nu = 0; % levels 0, 0.1, 0.2, 0.3

% ds/dt = -beta_unvac*s*i - nu*s
% di/dt beta_unvac*s*i + beta_vac*v*i - gamma*i
% dv/dt = nu*s- beta_vac*v*i
% dr/dt = gamma*i \\ can be solved for from s + i + v + r = 1

% diff equations with forward euler:
k = 0.0001; % look at this
Tf = 90; % days
numsteps = ceil(Tf/k);
k = Tf/numsteps;

s = zeros(1, numsteps);
i = zeros(1, numsteps);
v = zeros(1, numsteps);
r = zeros(1, numsteps);

s(1) = S/N; % prop of population susceptible
i(1) = I/N; % prop of population infected
v(1) = V/N; % prop of vaccinated individuals
r(1) = R/N; % prop of population recovered

for n = 1:(numsteps-1)
    s(n+1) = s(n) + k*(-beta_unvac*s(n)*i(n) - nu*s(n));
    i(n+1) = i(n) + k*(beta_unvac*s(n)*i(n) + beta_vac*v(n)*i(n) - gamma*i(n));
    v(n+1) = v(n) + k*(nu*s(n)- beta_vac*v(n)*i(n));
    r(n+1) = r(n) + k*gamma*i(n);
end

plot(1:numsteps, s, 'b', ...
    1:numsteps, i, 'r', ...
    1:numsteps, v, 'g', ...
    1:numsteps, r, 'k')
legend('Susceptible', 'Infected', 'Vaccinated', 'Recovered')

max(i)