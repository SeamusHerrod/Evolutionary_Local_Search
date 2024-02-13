#include <iostream>
#include <vector>
#include <random>
#include <cmath>
using namespace std;

// the epsilon to be used 
static float EPS = 0.01;
// number of generations to run
static int GEN = 100;
// function to minimize, F 
double F(double x) {
  return 4 + 2*x + 2*sin(20*x) - 4*(x*x);
}

// function for mutating a point 
double mutate (double x) {
  // random value for mutating x 
  // generate random number between 0 and 1
  random_device rd;
  mt19937 gen(rd());
  uniform_real_distribution<> dis(0, 1);
  double r = dis(gen);

  if (r < 0.3) { // x - epsilon w/ prob 0.3 
    return x - EPS;
  }
  else if (r < 0.7) { // copy x w/ prob 0.4 
    return x;
  }
  else { // x + epsilon w/ prob 0.3
    return x + EPS;
  }
}

double selection (const vector<double>& population, const vector<double>& fitness) {
  // compute the total fitness of the population
  double TF = 0.00;
  for (int i = 0; i < fitness.size(); i++) {
    TF += fitness[i];
  }
  //cout << TF <<endl;
  
  // generate random number between 0 and TF
  random_device rd;
  mt19937 gen(rd());
  uniform_real_distribution<> dis(0, TF);
  double r = dis(gen);

  double prev_fitness = 0.0;
  for (int i = 0; i < population.size(); i++) {
    double cur_fitness = prev_fitness + (fitness[i] / TF); // adding fitness to our running total allows us to have the effect of 
                                                           // selecting individuals with probabilities proportional to fitness
    //cout << "individual: " << population[i] << "\tfittness: " << fitness[i] << "\trunning total: " << running_total << "\tr: " << r << endl;
    if (r < cur_fitness) {
      return population[i];
    }
    prev_fitness = cur_fitness;
  }
  // if all individuals in the population are not sufficiently fit, then just return the last element
  return population.back(); 
}

double best_in_gen(const vector<double>& population){
  double cur_best = population[0];
  for (int i = 1; i < population.size(); i++){
    if (F(population[i]) < F(cur_best)) {
      cur_best = population[i];
    }
  }
  return cur_best;
}

void weight_proportional_selection() {

  // initialize population 
  cout << F(0.87) << endl;
  int N = 100;
  vector<double> population(N);
  for (int i = 0; i < N; i++) { // initial population 0, 0.01, 0.02, ... 0.09 
    population[i] = 0.01 * i;
  }

  for (int i = 0; i < GEN; i++) {
    cout << "=================\nBEGINNING GENERATION " << i << "\n=================" << endl;
    // compute fitness of population
    vector<double> fitness(N);
    for (int i = 0; i < N; i++) {
      fitness[i] = F(population[i]);
    }
    // create next population 
    vector<double> next_gen(N);
    for (int i = 0; i < N; i++) {
      double sel = selection(population, fitness); // select a member from population (proportional to its fitness)
      sel = mutate(sel); // mutate based on epsilon
      if (sel < 0.0) {
        sel = 0.0;
      }
      if (sel > 1.0){
        sel = 1.0;
      }
      next_gen[i] = sel;
    }
    population = next_gen;

    cout << "\n\nNEW population " << i << endl;
    for (int i = 0; i < population.size(); i++){
      cout << "individual: " << population[i] << endl;
    }
    double best = best_in_gen(population);
    cout << "best of next generation: \n" << endl;
    cout << "Ind: " << best << "Fitness: " << F(best) << endl;

  }

}

double crossover(double x, double y) {
  // generate random number between 0 and 1
  random_device rd;
  mt19937 gen(rd());
  uniform_real_distribution<> dis(0, 1);
  double r = dis(gen);
  return r*x + (1-r)*y;
}

void weight_proportional_selection_with_crossover() {
  // initialize population 
  cout << F(0.87) << endl;
  int N = 100;
  vector<double> population(N);
  for (int i = 0; i < N; i++) { // initial population 0, 0.01, 0.02, ... 0.09 
    population[i] = 0.01 * i;
  }

  for (int i = 0; i < GEN; i++) {
    cout << "=================\nBEGINNING GENERATION " << i << "\n=================" << endl;
    // compute fitness of population
    vector<double> fitness(N);
    for (int i = 0; i < N; i++) {
      fitness[i] = F(population[i]);
    }
    // create next population 
    vector<double> next_gen(N);
    for (int i = 0; i < N; i++) {
      double selx = selection(population, fitness); // select parent 1 for creating offspring
      double sely = selection(population, fitness); // select a second parent for creating offspring
      double next_sel = crossover(selx, sely); // crossover the parents based on convex combination 
      double baby = mutate(next_sel); // mutate the baby
      // pruning values to lie within the desired range
      if (baby < 0.0) {
        baby = 0.0;
      }
      if (baby > 1.0){
        baby = 1.0;
      }
      // add baby to the next generation 
      next_gen[i] = baby;
    }
    // the next generation takes over
    population = next_gen;

    cout << "\n\nNEW population " << i << endl;
    for (int i = 0; i < population.size(); i++){
      cout << "individual: " << population[i] << endl;
    }
    double best = best_in_gen(population);
    cout << "best of next generation: \n" << endl;
    cout << "Ind: " << best << "Fitness: " << F(best) << endl;

  }
}
int main() {

  weight_proportional_selection();
  cout << "============\nWEIGHT PROPORTIONAL WITH CROSSOVER\n============" << endl;
  weight_proportional_selection_with_crossover();
  return 0;
}
