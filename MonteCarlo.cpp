/*
This file counts as my first completed C++ project, where I wrote a programme which estimates the price of European call options for a given strike price and time to expiry. 
The function requires a csv file containing a column of strike prices and a column of times to expiry, and outputs a csv file containing the price estimate, and confidence interval for each (strike, expiry) pair. 
The theory for this comes from chapter 21 of Hull's: Options, Futures, and Other Derivatives. I relied on a few extra sources, 
which I have included in // comments in within the file where necessary - for example, the inverse standard normal cdf used to compute the quantiles for the confidence interval. 

The sample paths are calculated using Euler-Maruyama simulation where the stock price dynamics follow the Black-Scholes model (under the risk-neutral distribution). 
The risk-free rate, drift, number of sample paths, and number of steps are all global variables which can be altered by the user. There is also another global variable which decides what % confidence interval is wanted;
by default I have set the parameter to obtain a 95% confidence interval. 

NEXT STEPS:
1. The file is very slow if you use a large number of steps for each sample path. I hope to learn about ways to optimise this, either through the method itself, or using tricks within C++ which speed up computations. 
2. I also want to generalise this programme so that I can calculate price estimates for other derivatives such as put options, or even more exotic derivatives with more complex payoffs. 
3. Furthermore, I want to be able to use Monte Carlo simulation for stock prices which don't follow a geometric Brownian motion. For example I want to next update the programme so that the user can decide which model the stock price should follow. 
    I think I could do this with the help of OOP, once I have finished learning about it. 
4. I want to also incorporate some other variance reduction techniques such as importance sampling, which can improve accuracy. 

Overall, I am happy that I have been able to complete this project using the knowledge of C++ I have currently. 
I view this project as an opportunity to revise some of the things I have already learnt, whilst also learn a bit more about some tricks, such as using different libraries. 
I hope to continue improving my C++ knowledge so that I can complete more advanced projects in the future.
*/






#include<iostream>
#include<cmath>
#include<numeric>
#include<math.h>
#include<fstream>
#include<sstream>
#include<vector>
#include<string>
#include<algorithm>
#include<array> 
#include<random> 
#include<chrono>
#include<tuple>

/*
Monte Carlo Pricer for European Call Options. 

Given an array of strike prices and times to expiry, this function, for fixed r, and sigma will compute estimates of the price of the options associtated to those dates, 
and returns the data in a .csv file. 

Strike prices and times to expiry should be read from a csv file. The csv file that this programme returns will have the following columns: 
1. strike price, 
2. time to expiry, t
3. the estimated price of the derivative, 
4. the confidence interval for the price estimate.
*/

// GLOBAL PARAMETERS:
float gRfr=0.04;                    // risk-free-rate ought to be a global variable since it shouldn't change much during the simulation. 
float gIV=0.55;                     // implied volatility should be found, and be constant for this simulation.
float gS0=1.0;                      // spot price of the asset with be constant also. 
const int gN = 5000;                // N = the size of the partition of the time horizon,
const int gn = 100;                 // n is the number of sample paths we wish to simulate. 
const float gAlpha=0.05;            // to generate a 100*(1-gAlpha)% confidence interval for the estimate of the price of the derivative.
const float gPI = 4*std::atan(1);   // pi, used for standard normal cdf inverse. 

// RANDOM NUMBER GENERATOR:
unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator(seed);
std::normal_distribution<float> distribution(0.0, 1.0);


// FUNCTIONS AND FUNCTION SIGNATURES:

    // derivative payoffs:
float callPayoff(const float S, const float K);
    // Euler-Maruyama simulator:
std::vector<std::vector<float>> EL_SamplePaths_GBM(const float T);
    // Monte Carlo Pricer:
std::pair<float, float> MC_Price(const float K, const float T);
    // computes confidence intervals:
std::pair<float, float> confInt(const std::pair<float, float>& mcEst, const float percentile);


    // File Stream:
void samplePaths_to_CSV(const std::vector<std::vector<float>>& paths, const std::string& title, const float T);
std::vector<std::pair<float, float>> read_input_data(const std::string& title);
void price_estimates_to_CSV(const std::vector<std::vector<float>> prices, const std::string& title);

    // Other:

    // standard normal generator:
std::vector<float> generateNormal(int m);
    // inverse erf function: 
float myErfInv2(float x);
    // inverse standard normal cdf:
float invPhi(float x);

int main(){
    // read the csv file with the strikes and expiry dates: 
        // the input file can have any title you want, but the columns must be ['strike price', 'time to expiry'] in that order.
        // The time to expiry column should be in days/252. 
    std::string inputDataTitle;
    std::cout<<"Input tite (remember the '.csv'!): ";
    std::getline(std::cin, inputDataTitle);
    std::vector<std::pair<float, float>> inputData = read_input_data(inputDataTitle);
    // compute the prices and create the new dataset:
    std::vector<std::vector<float>> prices;
    float quantile = invPhi(1-0.5*gAlpha);      // quantile for the confidence intervals
    for (int i=0; i<inputData.size(); i++){
        auto [strike, tau] = inputData[i];
        std::pair<float, float> priceEst = MC_Price(strike, tau);
        std::pair<float, float> confidenceInterval = confInt(priceEst, quantile);
        std::vector<float> row = {
            strike, tau, priceEst.first, priceEst.second, confidenceInterval.first, confidenceInterval.second
        };
        prices.push_back(row);
    }
    // export the dataset to csv:
    price_estimates_to_CSV(prices, "Price_Estimates.csv");

}



// Derivative payoffs: 
float callPayoff(const float S, const float K){
    return (S>K) ? (S-K) : 0.0;
}

// Euler-Maruyama simulator: 
std::vector<std::vector<float>> EL_SamplePaths_GBM(const float T){
    std::vector<std::vector<float>> samplePaths;                                    // to store the sample paths. 
    float dt = T/gN;                                                                // universal time incrememnt across all sample paths
    // for each sample path: 
    for (int i=0; i<gn; i++){
        std::vector<float> path={gS0};                                              // to store the sample path
        std::vector<float> Z = generateNormal(gN+1);                                // generate the standard normal sample
        for (int j=0; j<gN; j++){
            float newS = path[j]*(1+gRfr*dt) + gIV*path[j]*Z[j]*std::sqrt(dt);      // use the Euler-Maruyama simulation
            path.push_back(newS);
        }
        samplePaths.push_back(path);
    }
    return samplePaths;
}

// Monte-Carlo Pricer: 
std::pair<float, float> MC_Price(const float K, const float T){
    /*
    The discount factor is exp(-r*T) where r is gRfr global variable. 

    We also compute the standard deviation of the discounted payoffs, so that we may compute confidence intervals.
    */


    std::vector<std::vector<float>> samplePaths = EL_SamplePaths_GBM(T);                                                // compute the sample paths
    std::vector<float> payoffs;                                                                                         // to store the payoffs
    float sumPayoffs=0;                                                                                                 // cumulative sum of the payoffs
    for (size_t i=0; i<gn; i++){
        float S_T_i = samplePaths[i][gN];
        float payoff = callPayoff(S_T_i, K);
        payoffs.push_back(payoff);
        sumPayoffs+=payoff;
    }
    // estimate the price:
    float meanPayoff = sumPayoffs/gn;
    float disc = exp(-gRfr*T);
    float priceEst = disc*meanPayoff;                                                                                   // Monte Carlo estimator of the Price
    // std.dev of discounted payoffs:
    float sumSq;
    for (size_t k=0; k<gn; k++){
        float x = disc*payoffs[k] - priceEst;
        sumSq+=std::pow(x, 2);
    }
    float stdDev = std::sqrt(sumSq/(gn-1));                                                                            // sample standard deviation of the payoffs
    std::pair<float, float> toReturn;
    toReturn = {priceEst, stdDev};
    return toReturn;
}
    
    // 100(1-percentile)% confidence interval for the Monte Carlo estimate of the price. 
std::pair<float, float> confInt(const std::pair<float, float>& mcEst, const float quantile){
    float lower, upper;
    float mu = mcEst.first;
    float stdDev = mcEst.second; 
    float estStdDev = stdDev*std::sqrt( (float) (gn-1)/ (float) gn);
    lower = mu - (quantile*estStdDev)/std::sqrt(gn);
    upper = mu + (quantile*estStdDev)/std::sqrt(gn);
    return {lower, upper};
}

    
// STANDARD NORMAL SAMPLER:
std::vector<float> generateNormal(int m){
    // generates a sample of size m of the standard normal distribution
    std::vector<float> gens;
    for (int j=0; j<m; j++){
        float z = distribution(generator);
        gens.push_back(z);
    }
    return gens;
}

// INVERSE erf FUNCTION: 
    // (taken from https://stackoverflow.com/questions/27229371/inverse-error-function-in-c at 10:13 on 16/11/25, by nimig18)
float myErfInv2(float x){
    float tt1, tt2, lnx, sgn;
    sgn = (x < 0) ? -1.0f : 1.0f;

    x = (1 - x)*(1 + x);        // x = 1 - x*x;
    lnx = logf(x);

    tt1 = 2/(gPI*0.147) + 0.5f * lnx;
    tt2 = 1/(0.147) * lnx;

    return(sgn*sqrtf(-tt1 + sqrtf(tt1*tt1 - tt2)));
    }

// INVERSE STANDARD NORMAL CDF: 
float invPhi(float x){
    float A = std::sqrt(2);
    return A*myErfInv2(2*x - 1 );
}

// FILE STREAM FUNCTIONS:
    // to export the sample paths to a csv file. 
void samplePaths_to_CSV(const std::vector<std::vector<float>>& paths, std::string& title, const float T){
    std::ofstream simluations(title);
    // add the first row as the partition of [0, T]:
    if (simluations.is_open()){
        for (int j=0; j<gN+1; j++){
            simluations<<j*(T/gN)<<", ";
        }
        simluations<<std::endl;
    }
    // now add the values:
    if (simluations.is_open()){
        for (int i=0; i<gn; i++){
            std::vector<float> path = paths[i];
            for (int j=0; j<path.size(); j++){
                simluations<<path[j]<<", ";
            }
            simluations<<std::endl;
        }
    }
    simluations.close();
}
    // to read the input data 
std::vector<std::pair<float, float>> read_input_data(const std::string& title){
    std::ifstream inputData(title);
    std::vector<std::pair<float, float>> data;
    if (inputData.is_open()){
        std::string line;
        while (std::getline(inputData, line)){
            std::stringstream ss(line);
            std::string K, T;
            std::getline(ss, K, ',');
            std::getline(ss, T, ',');
            std::pair<float, float> datapoint;
            datapoint = {std::stod(K), std::stod(T)};
            data.push_back(datapoint);
        }
    }
    else std::cout<<"File Failed to Open."<<std::endl;
    inputData.close();
    return data;
}

    // to return the Monte Carlo estimates with confidence intervals as a csv. 
void price_estimates_to_CSV(const std::vector<std::vector<float>> prices, const std::string& title){
    std::ofstream outputData(title);
    if (outputData.is_open()){
        outputData<<"strike,"<<"time to expiry,"<<"price,"<<"confidence interval lower,"<< "confidence interval upper"<<std::endl;
        for (int i=0; i<prices.size(); i++){
            std::vector<float> row = prices[i];
            float strike = row[0];
            float tau = row[1];
            float price = row[2];
            float lower = row[4];
            float upper = row[5];
            outputData<<strike<<","<<tau<<","<<price<<","<<lower<<","<<upper<<std::endl;
        }
    }
    outputData.close();

}
