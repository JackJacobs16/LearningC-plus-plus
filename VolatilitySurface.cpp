#include<iostream> 
#include<vector>
#include<tuple> 
#include<algorithm> 
#include<string>
#include<sstream> 
#include<cstring> 
#include<fstream> 
#include<cmath>
#include<math.h>
#include<numeric>
#include<algorithm> 
#include<ctime> 
#include<chrono>

// algorithm to be used to compute difference between dates taken from: https://howardhinnant.github.io/date_algorithms.html#days_from_civil

// Returns number of days since civil 1970-01-01.  Negative values indicate
//    days prior to 1970-01-01.
// Preconditions:  y-m-d represents a date in the civil (Gregorian) calendar
//                 m is in [1, 12]
//                 d is in [1, last_day_of_month(y, m)]
//                 y is "approximately" in
//                   [numeric_limits<Int>::min()/366, numeric_limits<Int>::max()/366]
//                 Exact range of validity is:
//                 [civil_from_days(numeric_limits<Int>::min()),
//                  civil_from_days(numeric_limits<Int>::max()-719468)]
template <class Int>
constexpr
Int
days_from_civil(Int y, unsigned m, unsigned d) noexcept {
    static_assert(std::numeric_limits<unsigned>::digits >= 18,
             "This algorithm has not been ported to a 16 bit unsigned integer");
    static_assert(std::numeric_limits<Int>::digits >= 20,
             "This algorithm has not been ported to a 16 bit signed integer");
    y -= m <= 2;
    const Int era = (y >= 0 ? y : y-399) / 400;
    const unsigned yoe = static_cast<unsigned>(y - era * 400);      // [0, 399]
    const unsigned doy = (153*(m > 2 ? m-3 : m+9) + 2)/5 + d-1;     // [0, 365]
    const unsigned doe = yoe * 365 + yoe/4 - yoe/100 + doy;         // [0, 146096]
    return era * 146097 + static_cast<Int>(doe) - 719468;
}


// Pi: 
const float pi = 4*std::atan(1);

// Standard normal distribution: 
    // standard normal pdf
float phi(float x){
    float A = 1/sqrt(2*pi);
    float B = exp(-0.5*x*x);
    return A*B;
}

    // standard normal cdf
float Phi(float x){
    return 0.5*(1+erf(x/sqrt(2)));
}



class optionEU{
    /*
    Class for a European option. Can be either call or put. The option type must be specified in the constructor - see below. 
    */

    // principle member variables: 
    const std::string expiryDate, underlyingAsset;
    const float strikePrice;
    std::string type;                                        // C for call, P for put. The code will fail if any other input is given. 
    float marketPrice, underlyingPrice, riskFreeRate;

    // derived member variables: 
    std::string id;
    float timeToExp, impliedVol, vega, moneyness, priceBS; 

    // functions for computing the Black-Scholes implied volatility. (Private as the user of this class shouldn't need to use these)
        // NEWTON-RAPHSON: 
            // initial guess:
    float initialGuessNR(){
        float A = std::abs(std::log(moneyness));
        return std::sqrt(2*A/timeToExp);
    }
            // objective function for the Newton-Raphson algorithm: 
    float objectiveNR(float sigma){
        float numerator = getPriceBS()-marketPrice;
        float denominator = getBlackScholesVega();
        return numerator/denominator; 
    }

    public: 
        // primary constructor: 
        optionEU(float strikePrice, std::string underlyingAsset, std::string expiryDate, float marketPrice, std::string type)
            : strikePrice(strikePrice), underlyingAsset(underlyingAsset), expiryDate(expiryDate), marketPrice(marketPrice), type(type){

            // make's sure the correct type is inputted. This shouldn't be an issue if your input csv file is correct. See `inRow` declaration below for correct format for the columns. 
            if (!type.compare("C")&&!type.compare("P")){
                std::cout<<"Error: Incorrect option type. Should be 'C' for call options and 'P' for put options."<<std::endl;
                std::cout<<"Use the function .setType() to correct the error."<<std::endl;
            }
            // all expiry dates should be of the form: "YYYY-MM-DD".
                // all option IDs should be as follows: underlying symbol + expiry Year (YY) + expiry month (MM)  + expiry month day (DD) + type + strike price rounded to nearest whole integer.
                // e.g. a call option on AAPL with strike price 65, expiry date 21st November 2025 should read: AAPL251121C65; 
            std::string id = underlyingAsset;
            std::string expYear, expMon, expMDay;
            expYear = expiryDate.substr(2, 2);
            expMon = expiryDate.substr(5, 2);
            expMDay = expiryDate.substr(8, 2);
            std::string strikeStr = std::to_string((int) strikePrice);
            id+=(expYear+expMon+expMDay)+type+strikeStr;
            this->id = id;

            // we also derive the time to expiry here: 
                // we use the days_from_civil function to compute the time-to-expiry from today's date to the expiry date. 
            std::string expDate = getExpiryDate();
            // get the current date: 
            time_t now = time(nullptr);
            struct tm *currTime = localtime(&now);
            int currY, currM, currD;
            currY = currTime->tm_year+1900;
            currM = currTime->tm_mon+1;
            currD = currTime->tm_mday;
            // get the expiry date: 
            int expY, expM, expD;
            expD = 10*(expDate[8]-'0')+(expDate[9]-'0');
            expM = 10*(expDate[5]-'0')+(expDate[6]-'0');
            expY = 1000*(expDate[0]-'0')+100*(expDate[1]-'0')+10*(expDate[2]-'0')+(expDate[3]-'0');
            int daysToExp = days_from_civil(expY, expM, expD)-days_from_civil(currY, currM, currD); // number of days to expiry
            // return number of days over number of trading days
            this->timeToExp =  (float) daysToExp / (float) 252;
            }
        
        // setters:
        void setMarketPrice(float V_mkt){
            this->marketPrice = V_mkt;
        }
        void setUnderlyingPrice(float S){
            this->underlyingPrice=S;
        }
        void setRiskFreeRate(float r, std::string desc=""){
            if (desc!=""){
                // desc here should describe the risk-free-rate, e.g. 3mo US Tbill, or BofE rate. 
                std::cout<<"Risk free rate defined as: "<<desc<<std::endl;
            }
            this->riskFreeRate=r;
        }
        void setMoneyness(char version){
            // obtain the moneyness of the option given the current price of the underlying asset. 
                // We must stipulate how we want to define moneyness, this can be K/S or K/F where F is the forward price of the underlying asset. 
                // Later I would like to add a way to stipulate moneyness based off of the delta of the option; I have seen volatility surfaces which use this as one axis. 
            // if version="s", we define moneyness as K/S. 
            // if version="f", we define moneyness as K/F.
            float S = getUnderlyingPrice();
            float K = getStrikePrice();
            float moneyness = K/S;
            if (version=='f'){
                float r = getRiskFreeRate();
                float tau = getTimeToExp();
                moneyness*=std::exp(-r*tau);
                this->moneyness=moneyness;
            }
            else if (version=='s'){
                this->moneyness=K/S;
            }
            else std::cout<<"Error: Incorrect version. Should be `s` for spot price or `f` for forward price. "<<std::endl;
        }
        void setBlackScholesVega(float sigma) {
            // compute the current vega of the option: 
                // is the same for both calls and puts thanks to put-call parity. 
            float tau = getTimeToExp();
            float K = getStrikePrice();
            float r = getRiskFreeRate();
            float S = getUnderlyingPrice();
            float A = std::log(S/K);
            float B = (r + (0.5*std::pow(sigma, 2)))*tau;
            float C = sigma*std::sqrt(tau);
            float d1 = (A+B)/C;
            float vega = S*std::sqrt(tau)*phi(d1);
            this->vega = vega;
        }
        void setPriceBS(float sigma) {
            // set's the Black-Scholes price of the option for a given volatility
            float tau = getTimeToExp();
            float K = getStrikePrice();
            float r = getRiskFreeRate();
            float S = getUnderlyingPrice();
            
            float d1 = (std::log(S/K) + tau*(r+0.5*std::pow(sigma, 2)))/(sigma*std::sqrt(tau));
            float d2 = d1 - sigma*std::sqrt(tau);
            float PhiD1 = Phi(d1); 
            float PhiD2 = Phi(d2);

            const char optionType = getType();
            switch(optionType){
                case 'C': 
                    this->priceBS = S*PhiD1 - K*std::exp(-r*tau)*PhiD2;
                    break;
                case 'P':
                    this->priceBS = S*(PhiD1-1) + (1-PhiD2)*K*std::exp(-r*tau);
                    break;
                default: 
                    break;
            }
            
        }
        void setType(std::string type){
            this->type = type;
        }
        
        // getters: 
        std::string getID()const {
            return this->id;
        }
        std::string getExpiryDate() const {
            return this->expiryDate;
        }
        std::string getUnderlyingAsset() const{
            return this->underlyingAsset;
        }
        float getStrikePrice() const {
            return this->strikePrice;
        }
        float getMarketPrice() const{
            return this->marketPrice;
        }
        float getUnderlyingPrice() const {
            return this->underlyingPrice;
        }
        float getTimeToExp() const{
            return this->timeToExp;
        }
        float getRiskFreeRate() const{
            // what we consider to be the risk-free-rate on the underlying asset. 
            return this->riskFreeRate;
        }
        float getMoneyness() const{
            if (this->moneyness) return this->moneyness;
            else std::cout<<"Moneyness has not been calculated yet. Returning 0.";
            return 0;
        }

        const char getType() const{
            // get set the type to be a character here so that we may use a switch statement in the calculation of the price of the option. 
            return this->type[0];
        }

        float getPriceBS() const{
            return this->priceBS;
        }
        float getBlackScholesVega() const{
            return this->vega;
        }
        float getImpliedVolatility() const{
            return this->impliedVol;
        }
        // compute the implied volatility: 
        void impliedVolatility(float epsilon, long unsigned maxIter) {
            /*
            Uses the Newton-Raphson algorithm to estimate the implied volatility of the option. 
            - epsilon is the threshold value for the obective function (once the objective falls below epsilon, return the value). 
            - maxIter is the maximum number of iterations allowed by the algorithm. 
            */ 
            float sigma = initialGuessNR();                                         // setup the initial value of the volatility
            setPriceBS(sigma);                                                      // compute the current Black-Scholes price;
            setBlackScholesVega(sigma);                                             // compute the current vega 
            float currentObj = objectiveNR(sigma);                                  // current objective function
            int N=0;                                                                // number of iterations completed
            while (std::abs(currentObj)>=epsilon&&N<maxIter){                               
                sigma -= currentObj;                                                // Newton-Raphson algorithm step
                setPriceBS(sigma);                                                  // update the Black-Scholes Price
                setBlackScholesVega(sigma);                                         // update the Vega
                float newObjective = objectiveNR(sigma);                            // update the objective functions
                currentObj = newObjective;              
                N++; 
            }
        // update the member variables for the final time: 
        this->impliedVol = sigma;
        setBlackScholesVega(this->impliedVol);
        setPriceBS(this->impliedVol);
        }
};

using inRow = std::tuple<std::string, float, std::string, float>;                                // how our input data should look like (expiry date, strike price, type, market price)
using outRow = std::tuple<std::string, std::string, float, char, float, float, float>;           // how our output data rows will look like (ID, expiry date, strike price, type, market price, implied volatility, implied Price)

int main(){
    /* 
    We use the fstream to read in a csv file of input data with columns outlined by `inRow`. The output is a csv file with columns outlined by `outRow`. The user must input the stock ticker (e.g. AAPL) so that option 
    ID's can be created, as well as the spot price of the underlying asset ($) at the time the option chain data was collected, as well as the constant risk-free-rate. 

    The programme will only output rows to the csv file if the price of the option implied by the estimated implied volatility is not a NaN value. If, for example, the Newton-Raphson algorithm returns + or - inf for the 
    implied volatility, this will lead to the implied price being NaN. In which case, we skip this row and move onto the next. The idea behind using both call and put options in this project is this: 
    For deep OTM call options, the Newton-Raphson algorithm struggles to converge since the Vega of the option falls to 0. However, this problem did not seem to occur when the call option was deep ITM. As a result, I 
    have made the programme in such a way that in most cases, deep OTM calls are ignored, but the equivalent put options, which are deep ITM, are included. By put-call-parity, the implied volatility will be the same for 
    both the call options and put options, as such we do not miss out on any important data. 
    */

    // key parameters: 
    float S, r;                                                             // spot price of underlying asset and the risk-free rate for the Black-Scholes model. 
    std::string ticker;                                                     // stock ticker of the underlying asset, e.g. AAPL or TSLA
    std::cout<<"Input the ticker of the underlying asset: ";
    std::getline(std::cin, ticker); 
    std::cout<<std::endl;
    std::cout<<"Spot price ($): ";
    std::cin>>S;
    std::cout<<std::endl;
    std::cout<<"Risk-free rate (%/100): ";
    std::cin>>r;
    std::cin.ignore(256, '\n');

    // for the Newton-Raphson algorithm: 
    float epsilon = 1e-16;      
    long unsigned maxIter = 10000;

    // read in the csv file: 
    std::vector<inRow> data;
    std::string inputFileTitle; 
    std::cout<<"Input file name: "<<std::endl;
    std::getline(std::cin, inputFileTitle);
    std::ifstream marketData(inputFileTitle);
    if (marketData.is_open()){
        std::string line;
        while (std::getline(marketData, line)){
            std::stringstream ss(line);
            std::string K, V_mkt, T, type;
            std::getline(ss, T, ',');
            std::getline(ss, K, ',');
            std::getline(ss, type, ',');
            std::getline(ss, V_mkt, ',');
            inRow dataRow; 
            float val = std::stof(V_mkt);
            float strike = std::stof(K);
            dataRow = {T, strike, type, val};
            data.push_back(dataRow);
        }
    } else std::cout<<"Error: File did not open."<<std::endl;
    marketData.close();

    std::ofstream outputData("Output_data.csv");
    if (outputData.is_open()){
        outputData<<"ID,"<<"Expiry Date,"<<"Strike Price,"<<"Type,"<<"Market Price,"<<"Implied Volatility,"<<"Implied Price"<<std::endl;
        for (size_t j=0; j<data.size(); j++){
            // compute the implied vol: 
            auto [T, K, type, V_mkt] = data[j];
            float K_copy = K/S;
            float V_mkt_copy = V_mkt/S;
            optionEU option(K_copy, ticker, T, V_mkt_copy, type);
            option.setUnderlyingPrice(1.0);                         // since we have normalised the strike price and market price, we must input the spot price as 1
            option.setRiskFreeRate(r);
            option.setMoneyness('f');
            std::string id = option.getID();
            option.impliedVolatility(epsilon, maxIter);
            float impliedVol = option.getImpliedVolatility();
            float impliedPrice = option.getPriceBS();
            // if we get a NaN value for the implied price, then the volatiliy is not valid, so we ignore this row: 
            if (std::isnan(impliedPrice)){
                continue;
            }
            // output to the file: 
                // remember to multiply the impled price by S so remove the normalisation. 
            outputData<<id<<","<<T<<","<<K<<","<<type<<","<<V_mkt<<","<<impliedVol<<","<<impliedPrice*S<<std::endl;
        }
    } else std::cout<<"Error: File not open."<<std::endl;
    outputData.close();
}