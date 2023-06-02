
/* To compile optimized:
 *   g++ -O3 -DNDEBUG abm_spaces_tb7.cc -lpthread
 *   Or if you're using clang:
 *   clang++ -O3 -DNDEBUG abm_spaces_tb7.cc -lpthread
 *
 * To compile for debugging:
 *   g++ -g abm_spaces_tb7.cc -lpthread
 *
 * Examples:
 * To run the default version:
 * ./a.out
 * To run it with tracing on 100 times:
 * ./a.out  tracing=1 num_runs=100
 * To run it with tracing on 10 times and test turnaround set to 8 weeks:
 * ./a.out  tracing=1 num_runs=10 tat=8
 * Command line arguments:
 * seed: Random seed (uses time otherwise)
 * cores: Number of CPU cores to use (defaults to number on machine)
 * num_runs: Number of simulations to run (default 1)
 * num_agents: Number of agents (default 20000)
 * num_time: Number of weeks (default 520)
 * tracing: 0 = no tracing, 1 = tracing
 * test_contacts: 0 = no testing of contacts, 1 = testing of contacts
 * mass_x_ray: 1 = implement mass x-ray programme, 0 = no mass x-ray programme (default = 0)
 * tut: Targeted universal testing - 1 = implement tut and 0 = no tut (default = 0)
 * tat: Test turnaround time in weeks (default 1)
 * initial_infection_rate: Proportion of initially exposed agents (default 0.002)
 * report_frequency: How frequently to report (default 52 to allow for anual stats in post-analysis)
 * asymptomatic_rate: Proportion of agents who will be asymptomatic 
 * death_rate: The proportion of agents with active TB who would die if not treated 
 * natural_death_rate: The yearly rate of deaths from causes other than TB
 * uniform_age: Uses default uniform age distribution if = 1, else uses SA age distribution
 * alt_activation: If = 0, samples from default uniform distribution to determine time of TB activation, if = 1, samples from more complex distribution with longer tail
 * alt_non_infectious: if = 0, agents become non-infectious two weeks after starting treatment, if =1, this timing is determined by sampling from a more complex non-uniform sdistribution
 * stay_home_rate: Proportion of agents who stay home when they have symptomatic TB disease  
 * linkage_rate: Proportion of agents who test positive that start treatment
 * testing_rate: Proportion of agents who would test if symptomatic
 * tracing_efficacy: Proportion of contacts traced
 * isolation_rate: Proportion of symptomatic who isolate
 * house_risk, class_risk, block_risk, taxi_risk, and work_risk can be used to set the setting specific risk variables (this is mainly used for model calibration)

 * NOTE: Exposed essentially = latent TB infection. Infectious includes both subclinical and symptomatic/active TB. The two basic pathways are:
 - Exposed > latent > cured
 - exposed > subclinical > symptomatic > death or cure

*/

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <ctime>
#include <future>
#include <iostream>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#define NOT_INFECTED(a) ((a.exposed + a.infectious + a.recovered + a.dead) == 0)

struct Option {
    std::string name;
    std::string description;
    double* value_d;
    int* value_i;
};

struct Parameters {
    int cores = 0;
    int num_runs = 1;
    int num_agents = 20000;
    int num_time = 520;
    int seed = 0;
    int report_frequency = 52;
    int tracing = 0;
    int mass_x_ray = 0;
    int test_contacts = 1;
    int tut = 0; // targeted universal testing
    int reg_tests = 0;
    int reg_tests_freq = 4;
    double reg_tests_sensitivity = 0.7;
    int tat = 1; // test turnaround time
    int output_parameters = 0;
    int sensitivity_tests =  0;
    int uniform_age = 0;
    int alt_activation = 0;
    int alt_non_infectious = 0;
    double asymptomatic_rate = 0;
    double subclinical_rate = 0.1;
    double linkage_rate = 0.8;
    double testing_rate = 0.647;
    double hiv_tb_risk_scaler = 4;
    double asymptomatic_infectiousness_ratio = 0.5; // If symptomatic agents have an infectiousness of 1, asymptomatic agents have infectiousness of 1 * this value
    double death_rate = 0.5;
    double natural_death_rate = 0.01; // fraction of population to die naturally per year
    double birth_rate = 0.011;
    double stay_home_rate = 0.2;
    double x_ray_sensitivity = 0.9;
    double x_ray_coverage = 0.5;
    double tracing_efficacy = 0.8;
    double isolation_rate = 0.5;
    double house_risk = 0.0127;
    double block_risk = 0.00115;
    double class_risk = 0.0076;
    double taxi_risk = 0.022;
    double work_risk = 0.1851;
    double initial_infection_rate = 0.02;
    double sensitivity_range = 0.0;
};

static thread_local std::mt19937 rng;

inline int rand_range_int(int a, int b)
{
    std::uniform_int_distribution<int> dist(a, b);
    return dist(rng);
}

inline double rand_0_1()
{
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    return dist(rng);
}

inline double rand_range(double a, double b)
{
    std::uniform_real_distribution<double> dist(a, b);
    return dist(rng);
}

int age_giver(int uniform_age) {
    if (uniform_age) {
        return rand_range_int(0,90);
    } else {
        std::discrete_distribution<int> distribution  {1027, 1016, 910, 820, 869,
            951, 926, 759, 597, 501, 420,
            357, 288, 217, 148, 96, 98};
        int x = distribution(rng);
        int age = (x * 5) + rand_range_int(0, 4);
        return age;
    }
}

int hiv_giver () {
    std::discrete_distribution<int> distribution {87, 7, 5, 1};
    int hiv = distribution(rng);
    return hiv;
}

int non_infectious_dist () {
  std::discrete_distribution<int> distribution {5, 20, 20, 15,  10, 10, 10, 5, 5}; // These probabilities refer to weeks, starting with 0 to which 2 will be added below
  int week_non_infectious = 2 + distribution(rng);
  return week_non_infectious;  
}

int activation_dist () {
    std::discrete_distribution<int> distribution {690, 110, 70, 35, 35, 25, 10, 15, 10}; // NB: these probabilities do not refer to equal amounts of time, but in seq to 18 months, 6, 6, 6, 12, 12, 12, 24, 24
    int activation = distribution(rng);
    int activation_time = 0;
    if (activation == 0) {
        activation_time = rand_range_int(2, 72);
    } else if (activation == 1) {
        activation_time = rand_range_int(73, 96);
    } else if (activation == 2) {
        activation_time = rand_range_int(97, 120);
  } else if (activation == 3) {
        activation_time = rand_range_int(121, 144);
    } else if (activation == 4) {
        activation_time = rand_range_int(145, 192);
    } else if (activation == 5) {
        activation_time = rand_range_int(193, 240);
    } else if (activation == 6) {
        activation_time = rand_range_int(241, 288);
    } else if (activation == 7) {
        activation_time = rand_range_int(289, 384);
    } else {
        activation_time = rand_range_int(385, 480);
    } 
    return activation_time;
}

struct Agent {
    int id_;
    int gender;
    int exposed = 0;
    int infectious = 0;
    int subclinical = 0;
    int symptomatic = 0;
    int latent = 0;
    int dead = 0;
    int dead_natural = 0;
    int hiv = 0; // 0 = negative, 1 = diagnosed and virally suppressed, 2 = diagnosed an not suppressed, 3 = not diagnosed and not suppressed
    int isolated = 0;
    int tested_positive = 0;
    int times_tested = 0;
    int weeks_infectious = 0;
    int weeks_exposed = 0;
    int weeks_active = 0;
    int time_for_subclinical = -1;
    int time_for_symptoms = -1;
    int time_for_x_ray = -1;
    int time_for_test = -25;
    int time_for_tut = -1;
    int time_for_reg_test = -1;
    int time_for_stay_home = -1;
    int time_for_cure = -1;
    int time_for_end_isolation = -1;
    int time_for_end_stay_home = -1;
    int time_for_non_infectious = -1;
    int time_for_death = -1;
    int time_for_natural_death = -1;
    int will_die = 0;
    int will_stay_home = 0;
    int will_be_subclinical = 0;
    int will_be_symptomatic = 0;
    int stay_home = 0;
    int will_link_to_care = 0;
    int will_test = 0;
    int age;
    int birth = 0;
    int household = 0;
    int block = 0;
    int class_ = -1;
    int workplace = -1;
    int recovered = 0;
    int taxi = 0;
    int times_quarantined = 0;
    int times_activeTB = 0;
    int times_subclinical = 0;
    int times_latent = 0;
    int times_positive = 0;
    int times_x_ray = 0;
    int times_x_ray_ref = 0;
    int cure_date =0;
    int where_infected = 0; // 1=home, 2=block, 3=class, 4=work, 5=taxi
    int birth_week;

    void init_agent(int id, double infection_rate, int uniform_age) {
        id_ = id;
        gender = rand_range_int(0,1);
        age = age_giver(uniform_age);
        hiv = hiv_giver();
        birth_week = rand_range_int(1, 52);
        if (rand_0_1() < infection_rate)
            exposed = 1;
    }
};

struct Simulation {

    std::vector<Agent> agents;
    int current_run;
    int household_size = 4;
    int households_per_block = 40;
    int num_schools = 1;
    unsigned class_size = 40;
    double employment_rate = 0.72;
    double regular_taxi_takers = 0.8;
    unsigned num_workplaces;
    unsigned taxi_capacity = 12;
    int tat;

    int week_first_subclinical = 38;
    int v_week_first_subclinical = 34;
    int weeks_sublinical_to_symptomatic = 26;
    int v_weeks_sublinical_to_symptomatic = 22;
    int symptoms_to_death = 48;
    int v_symptoms_to_death = 24;
    int symptoms_to_test = 3;
    int v_symptoms_to_test = 2;
    int symptoms_to_cure = 36;
    int v_symptoms_to_cure = 12;
    int time_step = 0;

    double asymptomatic_rate;
    double linkage_rate;
    double testing_rate;
    double asymptomatic_infectiousness_ratio;
    double hiv_tb_risk_scaler;
    double subclinical_rate;
    double death_rate;
    double natural_death_rate;
    double birth_rate;
    double stay_home_rate;
    double tracing_efficacy;
    double isolation_rate;
    double x_ray_coverage;
    double x_ray_sensitivity;
    double house_risk;
    double block_risk;
    double class_risk;
    double taxi_risk;
    double work_risk;
    double reg_tests_sensitivity;
    int reg_tests_freq;
    int tracing;
    int test_contacts;
    int mass_x_ray;
    int tut;
    int reg_tests;
    int uniform_age;
    int alt_activation;
    int alt_non_infectious;

    int peak = 0;
    int peak_time = 0;
    int total_infections_at_peak;

    std::vector< std::vector<int> > household_indices;
    std::vector< std::vector<int> > block_indices;
    std::vector< std::vector<int> > class_indices;
    std::vector< std::vector<int> > work_indices;
    std::vector< int > taxi_indices;

    std::ostringstream out;

    int count_exposures() {
        int total = 0;
        for (auto & agent: agents)
            total += agent.exposed;
        return total;
    }

    void make_agents(int num_agents, double infection_rate, int uniform_age)  {
        for (int i = 0; i < num_agents; i++) {
            Agent a;
            a.init_agent(i, infection_rate, uniform_age);
            agents.push_back(a);
        }
        assert(count_exposures() > 0);
    }

    void make_households()  {
        int x = 0;
        int house_num = 0;
        int len = agents.size();
        while (x < len) {
            std::vector<int> indices;
            int size_next_house = rand_range_int(2, household_size * 2 + 2);
            if ( len - x <= (household_size * 2) + 2) {
                size_next_house = agents.size() - x;
            }
            for (int j = x;  j < x + size_next_house; j++) {
                agents[j].household = house_num;
                indices.push_back(j);
            }
            assert(indices.size() > 0);
            household_indices.push_back(indices);
            ++house_num;
            x += size_next_house;
        }
    }

    int check_blocks() {
        int i = 0;
        for (auto & block : block_indices) {
            if (block.size() < 2)
                return 0;
            ++i;
        }
        return 1;
    }

    void make_blocks()  {
        int num_blocks = std::round(agents.size() / households_per_block);
        for (int i = 0; i < num_blocks; i++) {
            std::vector<int> indices;
            block_indices.push_back(indices);
        }
        for (auto & agent: agents) {
            std::vector<int> indices;
            agent.block = agent.household % num_blocks;
            block_indices[agent.block].push_back(agent.id_);
        }
        if (block_indices.back().size() == 0)
            block_indices.pop_back();
        assert(check_blocks());
    }

    int check_classes() {
        for (auto & room: class_indices) {
            if(room.size() < 2) {
                std::cerr << "Too small class" << std::endl;
                return 0;
            }
        }
        for (auto & room : class_indices) {
            int r = agents[room[0]].class_;
            int a = agents[room[0]].age;
            for (auto i: room) {
                if (agents[i].class_ != r) {
                    std::cerr << "Classroom mismatch: "
                              << i << " "
                              << agents[i].class_ << " "
                              << r << std::endl;
                    return 0;
                }
                if (agents[i].age != a) {
                    std::cerr << "Age mismatch: "
                              << i << " "
                              << agents[i].age << " "
                              << a << std::endl;
                    return 0;
                }
            }
        }
        return 1;
    }

    void make_classes()  {
        std::vector<int> indices[12];

        for (auto& agent: agents) {
            if (agent.age >=6 && agent.age < 18) {
                indices[agent.age - 6].push_back(agent.id_);
            }
        }

        for (auto& indice: indices)
            std::shuffle(indice.begin(), indice.end(), rng);

        int room = -1;
        for (unsigned i = 0; i < 12; i++) {
            for (size_t j = 0; j < indices[i].size(); j++) {
                if (j % class_size == 0) {
                    if (indices[i].size() - j > class_size / 2) {
                        std::vector<int> v;
                        class_indices.push_back(v);
                        room++;
                    }
                }
                assert(agents[indices[i][j]].age >= 6 &&
                       agents[indices[i][j]].age < 18);
                agents[indices[i][j]].class_ = room;
                class_indices.back().push_back(indices[i][j]);
            }
        }

        assert(check_classes());
    }

    int check_workplaces() {
        int w = 0;
        for (auto & work: work_indices) {
            if (work.size() < 1) {
                std::cerr << "Wrong workplace size: " << w << " "
                          << work.size() << std::endl;
                return 0;
            }
            for (auto i : work) {
                if (agents[i].workplace != w) {
                    std::cerr << "Wrong workplace: " << agents[i].workplace
                              << " " << w << " " << work_indices.size()
                              << std::endl;
                    return 0;
                }
            }
            ++w;
        }
        return 1;
    }

    void make_workplaces() {
        std::vector<int> indices;

        for (auto& agent: agents) {
            if (agent.age >= 18 && agent.age < 61) {
                if (rand_0_1() < employment_rate) {
                    indices.push_back(agent.id_);
                }
            }
        }
        std::shuffle(indices.begin(), indices.end(), rng);

        int avg_workplace_size = std::round((double) indices.size() / num_workplaces);
        int workplace_size = -1;
        int workplace = -1;
        int j = 0;

        for (auto i: indices) {
            if (j > workplace_size) {
                ++workplace;
                std::vector<int> v;
                work_indices.push_back(v);
                j = 0;
                workplace_size = rand_range_int(2, 2 * avg_workplace_size + 2);
            }
            agents[i].workplace = workplace;
            work_indices.back().push_back(i);
            j++;
        }
        assert(check_workplaces());
    }

    void make_taxis() {
        for (auto& agent: agents) {
            if (agent.workplace > -1) {
                if (rand_0_1() < regular_taxi_takers) {
                    agent.taxi = 1;
                    taxi_indices.push_back(agent.id_);
                }
            }
        }
    }

    void house_transmit() {
        double num_infected;
        for (auto& indices: household_indices) {
            std::vector<int> house;
            num_infected = 0;
            for (auto i: indices) {
                if (agents[i].isolated == 0) {
                    house.push_back(i);
                    if (agents[i].infectious)
                        ++num_infected;
                    if (agents[i].infectious ==1  && agents[i].symptomatic == 0)
                        num_infected = num_infected - asymptomatic_infectiousness_ratio;
                }
            }
            if (num_infected) {
                double this_house_risk = pow(1 - house_risk, num_infected);
                for (auto i : house) {
                    if (NOT_INFECTED(agents[i])) {
                        if (rand_0_1() > this_house_risk) {
                            agents[i].exposed = 1;
                            agents[i].where_infected = 1;
                        }
                    }
                }
            }
        }
    }


    void block_transmit() {
        for (auto& indices: block_indices) {
            std::vector<int> current_block;
            double infections = 0;
            for (auto i: indices) {
                if (agents[i].isolated + agents[i].stay_home == 0) {
                    current_block.push_back(i);
                    infections += agents[i].infectious;
                    if (agents[i].infectious == 1 && agents[i].symptomatic == 0) infections -= asymptomatic_infectiousness_ratio;
                }
            }
            if (current_block.size()) {
                double current_block_risk = pow(1 - block_risk, infections);
                for (auto i: current_block) {
                    if (NOT_INFECTED(agents[i])) {
                        if (rand_0_1() > current_block_risk) {
                            agents[i].exposed = 1;
                            agents[i].where_infected = 2;
                        }
                    }
                }
            }
        }
    }

    void class_transmit() {
        for (auto& room_indices: class_indices) {
            double infected = 0;
            for (auto i: room_indices) {
                infected += agents[i].infectious;
                if (agents[i].infectious == 1 && agents[i].symptomatic == 0) {
                    infected = infected - asymptomatic_infectiousness_ratio;
                }
            }
            double this_class_risk = pow(1 - class_risk,  infected);
            for (auto i: room_indices) {
                if (NOT_INFECTED(agents[i]) &&
                    agents[i].isolated == 0) {
                    if (rand_0_1() > this_class_risk) {
                        agents[i].exposed = 1;
                        agents[i].where_infected = 3;
                    }
                }
            }
        }
    }

    void work_transmit() {
        for (auto &indices : work_indices) {
            double infected = 0;
            for (auto i: indices) {
                if (agents[i].stay_home + agents[i].isolated + agents[i].dead + agents[i].dead_natural == 0)
                    infected += agents[i].infectious;
                if (agents[i].infectious == 1 && agents[i].symptomatic == 0) {
                    infected = infected - asymptomatic_infectiousness_ratio;
                }
            }
            double this_workplace_risk = pow(1 - work_risk, infected);
            for (auto i: indices) {
                if (NOT_INFECTED(agents[i]) &&
                    agents[i].isolated == 0) {
                    if (rand_0_1() > this_workplace_risk) {
                        agents[i].exposed = 1;
                        agents[i].where_infected = 4;
                    }
                }
            }
        }
    }

    void taxi_transmit() {
        std::vector<int> indices;
        for (auto i: taxi_indices) {
            if (agents[i].isolated + agents[i].stay_home + agents[i].dead + agents[i].dead_natural == 0)
                indices.push_back(i);
        }
        shuffle(indices.begin(), indices.end(), rng);
        unsigned current = 0;
        int num_taxis = (double) taxi_indices.size() / taxi_capacity;
        for (int i = 0; i < num_taxis; i++) {
            double infections = 0;
            size_t j = current;
            std::vector<int> taxi;
            for (;j < current + taxi_capacity && j < indices.size(); j++) {
                taxi.push_back(indices[j]);
                infections += agents[indices[j]].infectious;
                if (agents[indices[j]].infectious == 1&& agents[indices[j]].symptomatic == 0) infections -= asymptomatic_infectiousness_ratio;
            }
            if (infections) {
                double current_taxi_risk = pow(1 - taxi_risk, infections);
                for (auto i: taxi) {
                    if (NOT_INFECTED(agents[i])) {
                        if (rand_0_1() > current_taxi_risk) {
                            agents[i].exposed = 1;
                            agents[i].where_infected = 5;
                        }
                    }
                }
            }
            current += taxi_capacity;
        }
    }

    int get_tested(int weeks_infected, int hiv) {
        int test_result = 0;
        int k = rand_range_int(1, 100);
        if (weeks_infected == 1) {
            if (k < 51 && hiv == 0) test_result = 1;
            if (k < 36 && hiv > 0) test_result = 1;
        } else if (weeks_infected == 2) {
            if (k < 76 && hiv == 0) test_result = 1;
            if (k < 54 && hiv > 0) test_result = 1;
        }
        else if (weeks_infected > 2) {
            if (k < 98 && hiv ==0) test_result = 1;
            if (k < 71 && hiv > 0) test_result = 1;
        }
        return test_result;
    }

    int get_reg_tested(double reg_tests_sensitivity) {
        int test_result = 0;
        if (rand_0_1() < reg_tests_sensitivity) {
            test_result  = 1;
}
        return test_result;
    }

    void household_tracing(struct Agent & agent) {
        std::vector<int> household;
        for (auto i: household_indices[agent.household]) {
            if (agents[i].tested_positive == 0 && agents[i].time_for_test < time_step - 12 &&
                agents[i].isolated == 0) {
                assert(agents[i].id_ != agent.id_);
                household.push_back(i);
            }
        }
        shuffle(household.begin(), household.end(), rng);
        int n = std::min(household.size(), (size_t)
                         std::round(tracing_efficacy * household.size()));
        for (int i = 0; i < n; i++) {
            int j = household[i];
            assert(agents[j].household == agent.household);
            if (test_contacts == 1 && agents[j].tested_positive ==0) {
                agents[j].time_for_test = time_step + tat + rand_range_int(1, 2);
            }
            if (rand_0_1() < isolation_rate) {
                agents[j].isolated = 1;
                agents[j].time_for_end_isolation = time_step + 2;
                ++agents[j].times_quarantined;
            } else {
                agents[j].stay_home = 1;
                agents[j].time_for_end_stay_home = time_step + 2;
            }
        }
    }

    void workplace_tracing(struct Agent & agent) {
        if (agent.workplace > -1) {
            std::vector<int> workplace;
            for (auto i: work_indices[agent.workplace]) {
                if (agents[i].tested_positive == 0 &&
                    agents[i].time_for_test < (time_step - 24) &&
                    agents[i].isolated == 0) {
                    assert(agents[i].id_ != agent.id_);
                    assert(agents[i].workplace == agent.workplace);
                    workplace.push_back(i);
                }
            }
            shuffle(workplace.begin(), workplace.end(), rng);
            int n = std::min(workplace.size(), (size_t)
                             std::round(tracing_efficacy * workplace.size()));
            for (int i = 0; i < n; i++) {
                int j = workplace[i];
                assert(agents[j].workplace == agent.workplace);
                if (test_contacts == 1 && agents[j].tested_positive == 0) {
                    agents[j].time_for_test = time_step + tat + rand_range_int(1, 2);}
                if (rand_0_1() < isolation_rate) {
                    agents[j].isolated = 1;
                    agents[j].time_for_end_isolation = time_step + 2;
                    ++agents[j].times_quarantined;
                } else {
                    agents[j].stay_home = 1;
                    agents[j].time_for_end_stay_home = time_step + 2;
                }
            }
        }
    }

    void class_tracing(struct Agent & agent) {
        if (agent.class_ > -1) {
            std::vector<int> room;
            for (auto i: class_indices[agent.class_]) {
                assert(agents[i].class_ == agent.class_);
                assert(agents[i].age == agent.age);
                if (agents[i].tested_positive + agents[i].isolated == 0 &&
                    agents[i].time_for_test < time_step - 24) {
                    assert(agents[i].id_ != agent.id_);
                    room.push_back(i);
                }
            }
            shuffle(room.begin(), room.end(), rng);
            int n = std::min(room.size(), (size_t)
                             std::round(tracing_efficacy * room.size()));
            assert(n <= (int) room.size());
            for (int i = 0; i < n; i++) {
                int j = room[i];
                assert(agents[j].class_ == agent.class_);
                if (test_contacts == 1 && agents[j].tested_positive == 0) {
                    agents[j].time_for_test = time_step + tat + rand_range_int(1, 2);}
                if (rand_0_1() < isolation_rate) {
                    agents[j].isolated = 1;
                    agents[j].time_for_end_isolation = time_step + 2;
                    ++agents[j].times_quarantined;
                } else {
                    agents[j].stay_home = 1;
                    agents[j].time_for_end_stay_home = time_step + 2;
                }
            }
        }
    }

    void give_job (int workplace_num, int id_num, int num_agents) {
        for (auto & agent: agents) {
            if (agent.workplace == -1 &&
                agent.age > 18 && agent.age < 56) {
                agent.workplace = workplace_num;
                for (size_t i = 0; i < work_indices[workplace_num].size(); i++) {
                    if (work_indices[workplace_num][i] == id_num) {
                        work_indices[workplace_num][i] = agent.id_;
                        break;
                    }
                }
                break;
            }
        }
    }

    void school_progress () {
								for (int i = class_indices.size() -1; i >=0; i--) {
												if (agents[class_indices[i][1]].age == 18) {
																for (size_t j =0; j < class_indices[i].size(); j++) {
																				agents[class_indices[i][j]].class_ = -1;
																}
																class_indices.pop_back();
												} else { 
																break;
												}
								}

								std::vector <int> temp_index;
								for (auto& agent: agents) {
												if (agent.age == 6) temp_index.push_back(agent.id_);
								}

								int agent_counter = 0;
								int class_counter = 0;
								int num_classes = std::ceil( (double) temp_index.size() / class_size);
								int modified_class_size = std::ceil( (double) temp_index.size() / num_classes);

								std::vector <std::vector<int>> vec(num_classes);								

								for (size_t k = 0; k < temp_index.size(); k++) {
												agents[temp_index[k]].class_ = class_counter;
												vec[class_counter].push_back(temp_index[k]);
												++agent_counter;
												if (agent_counter == modified_class_size) {
																agent_counter = 0;
																++class_counter;
												}
								}
								class_indices.insert(class_indices.begin(), vec.begin(), vec.end());
    }

    void disease_progress(int num_agents) {
        for (auto& agent: agents) {
            if (agent.dead == 1 || agent.dead_natural ==1)
                continue;
            if (agent.infectious == 1)
                ++agent.weeks_infectious;
            if (agent.symptomatic == 1)
                ++agent.weeks_active;

            if (agent.time_for_subclinical == time_step) {
                agent.subclinical = 1;
                ++ agent.times_subclinical;
                agent.infectious = 1;   
            }
            if (agent.time_for_symptoms == time_step) {
                agent.symptomatic = 1;
                ++ agent.times_activeTB;
                agent.subclinical = 0;
            }
            if (agent.time_for_stay_home == time_step) {
                agent.stay_home = 1;
            }
            if (agent.time_for_end_stay_home == time_step) {
                agent.stay_home = 0;
            }
            if (agent.time_for_end_isolation == time_step) {
                agent.isolated = 0;
            }
            if (agent.time_for_non_infectious == time_step) {
                agent.exposed = 0;
                agent.infectious = 0;
                agent.symptomatic = 0;
                agent.time_for_death = -1;
            }
            if (agent.time_for_x_ray == time_step && agent.tested_positive == 0) {
                if (rand_0_1() < x_ray_coverage && agent.time_for_test < time_step - 2) {
                    ++agent.times_x_ray;
                    if (rand_0_1() < x_ray_sensitivity && (agent.subclinical || agent.symptomatic)) {
                        agent.time_for_test = time_step + 1;
                        ++ agent.times_x_ray_ref;
                    }
                    agent.time_for_x_ray = agent.time_for_x_ray + 52;
                }
            }
            if ((agent.time_for_test == time_step || agent.time_for_tut == time_step || agent.time_for_reg_test == time_step) && agent.tested_positive + agent.dead + agent.dead_natural == 0) {
                int weeks_infected = agent.weeks_exposed + agent.weeks_infectious - tat;
                if (agent.subclinical == 1 || agent.symptomatic == 1) {
                    if (time_step == agent.time_for_reg_test) {agent.tested_positive = get_reg_tested(reg_tests_sensitivity);
                    }
                    else {agent.tested_positive = get_tested(weeks_infected, agent.hiv);
                    }
                }
                ++agent.times_tested;
                if (agent.tested_positive == 1) {
                    ++agent.times_positive;
                    if (agent.will_link_to_care) {
                        if (alt_non_infectious ==0) agent.time_for_non_infectious = time_step + 2;
                        else agent.time_for_non_infectious = time_step + non_infectious_dist(); 
                        agent.time_for_cure = agent.time_for_non_infectious + 14;
                        agent.will_die = 0;
                        agent.time_for_death = -1;
                    }
                }
                if (tracing && agent.tested_positive && agent.will_link_to_care) {
                    agent.isolated = 1;
                    agent.time_for_end_isolation = time_step + 4;
                    ++agent.times_quarantined;
                    household_tracing(agent);
                    // workplace_tracing(agent);
                    // class_tracing(agent);
                }
            }
            if (agent.time_for_cure == time_step) {
                agent.infectious = 0;
                agent.recovered = 1;
                agent.isolated = 0;
                agent.stay_home = 0;
                agent.tested_positive  = 0;
                agent.latent = 0;
                agent.symptomatic = 0;
                agent.cure_date = time_step;
                continue;
            }
            if (agent.time_for_death == time_step) {
                agent.infectious = 0;
                agent.symptomatic = 0;
                agent.dead = 1;
                agent.isolated = 0;
                if (agent.workplace >= 0)
                    give_job(agent.workplace, agent.id_, num_agents);
                agent.workplace = 0;
            }
            if (agent.time_for_natural_death == time_step) {
                agent.infectious = 0;
                agent.symptomatic = 0;
                agent.dead_natural = 1;
                agent.isolated = 0;
                if (agent.workplace >= 0)
                    give_job(agent.workplace, agent.id_, num_agents);
                agent.workplace = 0;
            }
        }
    }

    int check_agent(struct Agent & agent) {
        if (agent.exposed + agent.infectious == 2) {
            std::cout << agent.id_ << " "
                      << agent.exposed << " " << agent.infectious << " "
                      << agent.where_infected << std::endl;
            return 0;
        }
        return 1;
    }

    void exposed_progress() {
        int from_i = week_first_subclinical - v_week_first_subclinical;
        int to_i = week_first_subclinical + v_week_first_subclinical;
        int from_s = weeks_sublinical_to_symptomatic - v_weeks_sublinical_to_symptomatic;
        int to_s = weeks_sublinical_to_symptomatic + v_weeks_sublinical_to_symptomatic;
        int from_t = symptoms_to_test - v_symptoms_to_test;
        int to_t = symptoms_to_test + v_symptoms_to_test;

        for (auto & agent: agents) {
            if (agent.dead + agent.dead_natural == 1)
                continue;
            if (agent.latent == 1)
                continue;
            if (agent.exposed ==1)
                ++agent.weeks_exposed;

            if (agent.infectious == 0 && agent.weeks_exposed == 1) {
                if (agent.hiv < 2) {
                    if (rand_0_1() < subclinical_rate) agent.will_be_subclinical = 1;
                    else agent.will_be_subclinical = 0;
                }
                if (agent.hiv > 1) {
                    if (rand_0_1() < subclinical_rate * hiv_tb_risk_scaler) agent.will_be_subclinical = 1;
                    else agent.will_be_subclinical  = 0;
                }
                if (agent.will_be_subclinical == 0) {
                    agent.time_for_cure = rand_range_int(4, 48);
                    agent.latent = 1;
                    ++agent.times_latent;
                }
                if (agent.will_be_subclinical == 1) {
                    if (rand_0_1() < testing_rate) agent.will_test =1;
                    else agent.will_test =0;
                    if (rand_0_1() < linkage_rate) agent.will_link_to_care =1;
                    else agent.will_link_to_care =0;
                    if (alt_activation == 0) agent.time_for_subclinical = time_step + rand_range_int(from_i, to_i);
                    else agent.time_for_subclinical = time_step + activation_dist();   
                    if (rand_0_1() < 1 - asymptomatic_rate) agent.will_be_symptomatic = 1;
                    else agent.will_be_symptomatic = 0;
                    if (agent.will_be_symptomatic == 1) {
                        agent.time_for_symptoms = agent.time_for_subclinical + rand_range_int(from_s, to_s);
                        if (agent.time_for_test < agent.time_for_symptoms && agent.will_test == 1) agent.time_for_test = agent.time_for_symptoms + rand_range_int(from_t, to_t) + tat;
                        if (rand_0_1() < stay_home_rate) agent.will_stay_home = 1;
                        else agent.will_stay_home = 0;
                        if (agent.will_stay_home == 1) {
                            agent.time_for_stay_home = agent.time_for_symptoms + rand_range_int(1, 5);
                            agent.time_for_end_stay_home += rand_range_int(1, 3);
                        }
                        if (rand_0_1() < death_rate) agent.will_die = 1;
                        else agent.will_die = 0;
                        if (agent.will_die == 0) {
                            agent.time_for_cure = agent.time_for_symptoms + rand_range_int(24, 48);
                        } else {
                            agent.time_for_death = agent.time_for_symptoms + rand_range_int(24, 72);
                        }
                    }
                }
            }

            if (agent.time_for_subclinical  == time_step) {
                agent.exposed = 0;
                agent.infectious = 1;
            }

            assert(check_agent(agent));
        }
    }

    void schedule_natural_deaths (double natural_death_rate, int num_agents, int num_time, int time_step) {
        int num_deaths = std::round(natural_death_rate * num_agents);
        for (int i = 1; i <= num_deaths; i++) {
            int rand_agent = rand_range_int(10, num_agents);
            int k = 0;
            while (k < 10) {
                if (agents[rand_agent - k].time_for_natural_death < 1) {
                    agents[rand_agent - k].time_for_natural_death = time_step + rand_range_int(1, 52);
                    k = 10;
                }
                else ++k;
            }
        }
    }

    void schedule_mass_x_ray () {
        for (auto & agent: agents) {
            if (agent.tested_positive == 0) agent.time_for_x_ray = agent.block + 2;
        }
    }

    void schedule_tut () {
        for (auto & agent: agents) {
            if ((agent.hiv > 0 && agent.hiv < 3 && agent.tested_positive == 0) || (agent.cure_date && agent.cure_date > time_step -48 && agent.tested_positive == 0)) {
                int temp = rand_range_int(1, 52);
                agent.time_for_tut = time_step + temp;
            }
        }
    }

    void schedule_reg_test () {
        for (auto & agent: agents) {
            if (agent.tested_positive == 0 && rand_0_1() > 0.8) {
                int temp = rand_range_int(1, reg_tests_freq);
                agent.time_for_reg_test = time_step + temp;
            }
        }
    }

    void births (double births_per_week, int num_agents, int uniform_age) {
        int births_this_week = floor(births_per_week);
        if (rand_0_1() < births_per_week - births_this_week) ++ births_this_week;
        for (int i=0; i < births_this_week; ++i) {
            Agent a;
            a.init_agent(num_agents + i, 0, uniform_age);
            a.hiv = 0;
            a.exposed = 0;
            a.age=0;
            a.birth = 1;
            int x = rand_range_int(1, agents.size());
            a.household = agents[x].household;
            a.block = agents[x].block;
            agents.push_back(a);
        }
    }

    void agent_ager () {
        for (auto & agent: agents) {
            if (agent.dead == 0) ++ agent.age;
        }
    }

    void retire_agents(int num_agents) {
        for (auto & agent: agents) {
            if (agent.age > 66 && agent.workplace > -1 && agent.dead == 0) {
                give_job(agent.workplace, agent.id_, num_agents);
                agent.workplace = 0;
            }
        }
    }

    void report(int num_time, int num_agents) {
        int total_tests = 0;
        int total_positive_tests = 0;
        int agents_ever_latent = 0;
        int total_latent_cases = 0;
        int total_active_cases = 0;
        int agents_ever_active = 0;
        int total_subclinical_cases = 0;
        int infectious = 0;
        int recovered = 0;
        int dead = 0;
        int dead_natural = 0;
        int births = 0;
        int live_pop = 0;
        int total_infections = 0;
        int infected_initial = 0;
        int infected_home = 0;
        int infected_school = 0;
        int infected_block = 0;
        int infected_work = 0;
        int infected_taxi = 0;
        int quarantines = 0;
        int total_x_ray = 0;
        int agents_ever_x_ray = 0;
        int x_ray_ref =0;
        int active_counter = 0;
        int non_infectious_counter = 0;        
        double mean_subclinical_to_non_infectious = 0;
        double mean_active_to_non_infectious = 0;
        double total_subclinical_to_non_infectious = 0;
        double total_active_to_non_infectious = 0;
        double incidence_100k  = 0;
        double R0 = 1.0 / ( (double) (agents.size() - total_infections_at_peak)
                            / agents.size());
        double percent_detected = 0.0;
        int hiv = 0;
        for (auto & agent: agents) {
            total_tests += agent.times_tested;
            total_positive_tests += agent.times_positive;
            if (agent.times_x_ray) {
                total_x_ray += agent.times_x_ray;
                ++agents_ever_x_ray;
                x_ray_ref += agent.times_x_ray_ref;
            }
            total_latent_cases += agent.times_latent;
            if (agent.times_latent) ++ agents_ever_latent;
            dead += agent.dead;
            dead_natural += agent.dead_natural;
            births += agent.birth;
            recovered += agent.recovered;
            infectious += agent.infectious;
            total_active_cases += agent.times_activeTB;
            if (agent.times_activeTB > 0) ++agents_ever_active;
            total_subclinical_cases  += agent.times_subclinical;
            if (agent.recovered || agent.latent || agent.dead || agent.exposed || agent.infectious)
                ++total_infections;
            infected_initial += (!NOT_INFECTED(agent) && agent.where_infected == 0);
            infected_home += (agent.where_infected == 1);
            infected_block += (agent.where_infected == 2);
            infected_school += (agent.where_infected == 3);
            infected_work += (agent.where_infected == 4);
            infected_taxi += (agent.where_infected == 5);
            quarantines += agent.times_quarantined;
            if (agent.hiv)++hiv;
            if (agent.weeks_infectious > 0) ++non_infectious_counter; 
            if (agent.weeks_active > 0) ++active_counter;
            total_subclinical_to_non_infectious = total_subclinical_to_non_infectious + agent.weeks_infectious;
            total_active_to_non_infectious = total_active_to_non_infectious + agent.weeks_active;   
        }
        mean_subclinical_to_non_infectious = total_subclinical_to_non_infectious  / non_infectious_counter;
        mean_active_to_non_infectious = total_active_to_non_infectious / active_counter;

        if (total_active_cases) {
            percent_detected =  ( (double) total_positive_tests /total_subclinical_cases) * 100.0;
        }

        live_pop = (num_agents + births -dead -dead_natural);

        incidence_100k = (total_active_cases/(num_time/52)) * (100000/num_agents);

        // assert(infected_initial + infected_home + infected_school + infected_block + infected_work + infected_taxi == total_infections);

        out << current_run
            << "," << time_step
            << "," << total_infections
            << "," << agents_ever_latent
            << "," << total_latent_cases
            << "," << agents_ever_active
            << "," << total_active_cases
            << "," << total_subclinical_cases
            << "," << infectious
            << "," << recovered
            << "," << dead
            << "," << dead_natural
            << "," << births
            << "," << live_pop
            << "," << peak
            << "," << peak_time
            << "," << total_infections_at_peak
            << "," << infected_initial
            << "," << hiv
            << "," << infected_home
            << "," << infected_school
            << "," << infected_block
            << "," << infected_work
            << "," << infected_taxi
            << "," << total_tests
            << "," << total_positive_tests
            << "," << percent_detected
            << "," << quarantines
            << "," << R0
            << "," << total_x_ray
            << "," << agents_ever_x_ray
            << "," << x_ray_ref
            << "," << tracing
            << "," << mass_x_ray
            << "," << tut
            << "," << reg_tests
            << "," << tat
            << "," << incidence_100k
            << "," << total_subclinical_to_non_infectious  
            << "," << total_active_to_non_infectious
            << "," << mean_subclinical_to_non_infectious
            << "," << mean_active_to_non_infectious
            << std::endl;
    }

    void run_model(int run, Parameters& parameters) {
        current_run = run;
        tracing = parameters.tracing;
        test_contacts = parameters.test_contacts;
        mass_x_ray = parameters.mass_x_ray;
        tut = parameters.tut;
        reg_tests = parameters.reg_tests; 
        tat = parameters.tat;
        subclinical_rate = parameters.subclinical_rate;
        asymptomatic_rate = parameters.asymptomatic_rate;
        linkage_rate = parameters.linkage_rate;
        testing_rate = parameters.testing_rate;
        asymptomatic_infectiousness_ratio = parameters.asymptomatic_infectiousness_ratio;
        hiv_tb_risk_scaler = parameters.hiv_tb_risk_scaler;
        death_rate = parameters.death_rate;
        natural_death_rate = parameters.natural_death_rate;
        birth_rate = parameters.birth_rate;
        x_ray_sensitivity = parameters.x_ray_sensitivity;
        x_ray_coverage = parameters.x_ray_coverage;
        stay_home_rate = parameters.stay_home_rate;
        tracing_efficacy = parameters.tracing_efficacy;
        isolation_rate = parameters.isolation_rate;
        reg_tests_freq = parameters.        reg_tests_freq;
        reg_tests_sensitivity = parameters.reg_tests_sensitivity;
        uniform_age = parameters.uniform_age;
        alt_activation = parameters.alt_activation;
        alt_non_infectious = parameters.alt_non_infectious;

        house_risk = parameters.house_risk;
        block_risk = parameters.block_risk;
        class_risk = parameters.class_risk;
        taxi_risk = parameters.taxi_risk;
        work_risk = parameters.work_risk;
        num_workplaces = (double) parameters.num_agents / 28.57142857142857142857;

        double births_per_week = (birth_rate * parameters.num_agents) / 52;

        make_agents(parameters.num_agents, parameters.initial_infection_rate, parameters.uniform_age);
        make_households();
        make_blocks();
        make_classes();
        make_workplaces();
        make_taxis();

        if (mass_x_ray) schedule_mass_x_ray();
        if (tut) schedule_tut();
        if (reg_tests) schedule_reg_test();

        for (time_step = 0; time_step < parameters.num_time; time_step++) {

            if (time_step == 0) schedule_natural_deaths(natural_death_rate, parameters.num_agents, parameters.num_time, time_step);

            if (time_step % reg_tests_freq == 0 && time_step > 1) {
                schedule_reg_test();
            }

            if (time_step % 52 == 0 && time_step > 1) {
                schedule_natural_deaths(parameters.natural_death_rate, parameters.num_agents, parameters.num_time, time_step);
                if (tut == 1) schedule_tut();
                agent_ager();
                retire_agents(parameters.num_agents);
                school_progress();
            }

            house_transmit();
            block_transmit();
            class_transmit();
            work_transmit();
            taxi_transmit();

            exposed_progress();
            disease_progress(parameters.num_agents);

            births(births_per_week, parameters.num_agents, uniform_age);

            int infected = 0, infected_total = 0;

            for (auto & agent: agents) {
                if (agent.symptomatic)
                    ++infected;
                if (agent.times_activeTB)
                    ++infected_total;
            }
            if (peak < infected) {
                peak = infected;
                peak_time = time_step;
                total_infections_at_peak = infected_total;
            }
            if ( (time_step + 1) % parameters.report_frequency == 0)
                report(parameters.num_time, parameters.num_agents);
        }
    }
};

double peturb(double mean, double range)
{
    double from = mean * (1.0 - range);
    double to = mean * (1.0 + range);

    return rand_range(from, to);
}

void peturb_parameters(Parameters &p)
{
    double range = p.sensitivity_range;
    p.asymptomatic_rate = peturb(p.asymptomatic_rate, range);
    p.death_rate = peturb(p.death_rate, range);
    p.stay_home_rate = peturb(p.stay_home_rate, range);
    p.tracing_efficacy = peturb(p.tracing_efficacy, range);
    p.isolation_rate = peturb(p.isolation_rate, range);
    p.house_risk = peturb(p.house_risk, range);
    p.block_risk = peturb(p.block_risk, range);
    p.class_risk = peturb(p.class_risk, range);
    p.taxi_risk = peturb(p.taxi_risk, range);
    p.work_risk = peturb(p.work_risk, range);
    p.initial_infection_rate = peturb(p.initial_infection_rate, range);
}

std::string run_one_simulation(int current_run, Parameters parameters)
{
    Simulation s;
    if (parameters.seed > 0) {
        rng.seed(parameters.seed + 11 * current_run);
    } else {
        rng.seed(time(NULL) + 111 * current_run);
    }
    s.run_model(current_run, parameters);
    return s.out.str();
}

void calc_cores(int & cores)
{
    if (cores == 0) {
        cores = std::thread::hardware_concurrency();
        if (cores == 0)
            cores = 1;
    }
}

void print_csv_header()
{
    std::cout << "Run,"
              << "Time_step,"
              << "latent_infections,"
              << "agents_ever_latent,"
              << "total_latent_cases,"
              << "agents_ever_active,"
              << "total_active_cases,"
              << "total_subclinical_cases,"
              << "Infectious,"
              << "Recovered,"
              << "Dead,"
              << "Dead_natural,"
              << "Births,"
              << "Live_pop,"
              << "Peak,"
              << "Peak_time,"
              << "Peak_total,"
              << "Initial,"
              << "HIV,"
              << "Home,"
              << "School,"
              << "Block,"
              << "Work,"
              << "Taxi,"
              << "Tests,"
              << "Positive,"
              << "P_Detected,"
              << "Quarantines,"
              << "R0,"
              << "total_x_ray,"
              << "agents_ever_x_ray,"
              << "x_ray_refs,"
              << "Tracing,"
              << "mass_x_ray,"
              << "tut,"
              << "reg_tests,"
              << "tat,"
              << "incidence_100k,"
              << "total_time_infectious,"
              << "total_time_active,"
              << "mean_time_infectious,"
              << "mean_time_active"
              << std::endl;
}

void run_simulations(Parameters & parameters)
{
    calc_cores(parameters.cores);
    int current_run = 0;
    print_csv_header();

    for (int i = 0; i < parameters.num_runs; i += parameters.cores) {
        std::vector<std::future<std::string> > output;
        int j = 0;
        for (; j < parameters.cores && current_run < parameters.num_runs;
             current_run++, j++) {
            Parameters p = parameters;
            if (parameters.sensitivity_range > 0.0)
                peturb_parameters(p);
            output.push_back(std::async(std::launch::async, run_one_simulation,
                                        current_run, p));
        }
        for (auto & s: output) {
            std::cout << s.get();
        }
    }
}

void print_help(const std::vector<Option>& options, const char *prog_name,
                const char *description)
{
    std::cout << prog_name;
    if (description)
        std::cout << ": " << description;
    std::cout << std::endl;

    std::cout << "Syntax:" << std::endl;
    std::cout << prog_name << " ";
    for (auto & option: options)
        std::cout << "[-" << option.name << "=<value>] ";
    std::cout << std::endl;

    std::cout << "\tOptions:" << std::endl;
    for (auto &option: options)
        std::cout << "\t-" << option.name << ": "
                  << option.description
                  << " (default: " << ( (option.value_d == NULL) ?
                                        *option.value_i : *option.value_d)
                  << ")" << std::endl;
}

void process_options(int argc, char *argv[], std::vector<Option>& options,
                     const char *prog_desc = NULL)
{
    for (int i = 1; i < argc; i++) {
        std::string arg(argv[i]);
        if (arg == "-h" || arg == "--help" || arg == "help") {
            print_help(options, argv[0], prog_desc);
            exit(0);
        }
        size_t p = arg.find('=');
        if (p == std::string::npos || p == arg.size() - 1) {
            print_help(options, argv[0], prog_desc);
            exit(1);
        }
        int start = 0;
        if (arg[0] == '-') {
            if (arg[1] == '-')
                start = 2;
            else
                start = 1;
        }
        std::string name = arg.substr(start, p-start);
        std::string value_s = arg.substr(p + 1);

        int found = 0;
        for (auto & option: options) {
            if (option.name == name) {
                found = 1;
                if (option.value_d) {
                    *option.value_d = std::stod(value_s);
                } else {
                    *option.value_i = std::stoi(value_s);
                }
            }
        }
        if (!found) {
            std::cerr << "Unknown option:" << arg << std::endl;
            std::cerr << "Try: " << std::endl;
            std::cerr << '\t' << argv[0] << " --help" << std::endl;
            exit(1);
        }
    }
}

int main(int argc, char *argv[])
{
    Parameters parameters;
    std::vector<Option> options = {
        {
            "cores",
            "Number of cores to use in multithreaded runs",
            NULL,
            &parameters.cores
        },
        {
            "seed",
            "Random seed",
            NULL,
            &parameters.seed
        },
        {
            "num_runs",
            "Number of simulations",
            NULL,
            &parameters.num_runs
        },
        {
            "num_agents",
            "Number of agents in the simulation",
            NULL,
            &parameters.num_agents
        },
        {
            "num_time",
            "Number of weeks in each simulation",
            NULL,
            &parameters.num_time
        },
        {
            "tracing",
            "Whether or not to do do tracing (0=off, 1=on)",
            NULL,
            &parameters.tracing
        },
        {
            "test_contacts",
            "Whether or not to test contacts (0=off, 1=on)",
            NULL,
            &parameters.test_contacts
        },
        {
            "tut",
            "Targeted universal testing (0=off, 1=on)",
            NULL,
            &parameters.tut
        },
        {
            "mass_x_ray",
            "Whether or not to do mass X-ray screening (0=off, 1=on)",
            NULL,
            &parameters.mass_x_ray
        },
        {
            "tat",
            "Turnaround time for a test",
            NULL,
            &parameters.tat
        },
        {
            "initial_infection_rate",
            "Proportion of agents exposed at beginning of simulation",
            &parameters.initial_infection_rate,
            NULL
        },
        {
            "house_risk",
            "Risk of house exposure",
            &parameters.house_risk,
            NULL
        },
        {
            "block_risk",
            "Risk of exposure in neighborhood",
            &parameters.block_risk,
            NULL
        },
        {
            "class_risk",
            "Risk of classroom exposure",
            &parameters.class_risk,
            NULL
        },
        {
            "taxi_risk",
            "Risk of exposure in taxi",
            &parameters.taxi_risk,
            NULL
        },
        {
            "work_risk",
            "Risk of exposure at work",
            &parameters.work_risk,
            NULL
        },
        {
            "reg_tests_sensitivity",
            "Sensitivity of tests in the high regular testing scenarios",
            &parameters.reg_tests_sensitivity,
            NULL
        },
        {
            "reg_tests_freq",
            "Test frequency in the high regularity testing scenarios",
            NULL,
            &parameters.reg_tests_freq
        },
        {
            "asymptomatic_rate",
            "Asymptomatic infection rate",
            &parameters.asymptomatic_rate,
            NULL
        },
        {
            "death_rate",
            "Infection death rate",
            &parameters.death_rate,
            NULL
        },
        {
            "natural_death_rate",
            "natural death rate",
            &parameters.natural_death_rate,
            NULL
        },
        {
            "uniform_age",
            "Set ages using uniform random",
            NULL,
            &parameters.uniform_age
        },
        {
            "alt_activation",
            "Sample from alternative distribution when determining time to TB activation",
            NULL,
            &parameters.alt_activation
        },
        {
            "alt_non_infectious",
            "Sample from alternative distribution when determining how soon agents on treatment become non-infectious",
            NULL,
            &parameters.alt_non_infectious
        },
        {
            "reg_tests",
            "Sets regular test scenario off or on  - 0 or 1 with 0 being the default",
            NULL,
            &parameters.reg_tests
        },
        {
            "stay_home_rate",
            "Rate at which possibly exposed who stay home",
            &parameters.stay_home_rate,
            NULL
        },
        {
            "linkage_rate",
            "Rate at which agents who test positive start treatment",
            &parameters.linkage_rate,
            NULL
        },
        {
            "testing_rate",
            "likelihood that symptomatic agent will test",
            &parameters.testing_rate,
            NULL
        },
        {
            "tracing_efficacy",
            "Proportion of contacts traced",
            &parameters.tracing_efficacy,
            NULL
        },
        {
            "isolation_rate",
            "Proportion of infected who isolate",
            &parameters.isolation_rate,
            NULL
        },
        {
            "sensitivity_range",
            "Modify parameters on each iteration proportionate to this",
            &parameters.sensitivity_range,
            NULL
        },
        {
            "report_frequency",
            "How often (in simulation weeks) to output the results",
            NULL,
            &parameters.report_frequency
        },
        {
            "output_parameters",
            "Whether to output the parameters",
            NULL,
            &parameters.output_parameters
        }
    };


    process_options(argc, argv, options, "Simulate TB");

    if (parameters.output_parameters) {
        for (auto & option: options) {
            std::cout << option.name << ": "
                      << ( (option.value_d == NULL) ?
                           *option.value_i : *option.value_d)
                      << std::endl;
        }
    }
    run_simulations(parameters);
}
