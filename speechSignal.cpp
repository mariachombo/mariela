// Team members: Nicholas, Maria Garcia, Thalia, Raquel
// Date: 10-30-15
// Project: Analysis of speech signals

#include <iostream>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <string>

using namespace std;

double mean(vector<double> &v);
// ***********************************************************
// Precondition: Function recieves a vector of values
// PostCondition: Function returns the average of the values given in the vector
// Summary: Function calculates the average from the values within given vector
// ***********************************************************
double avgMagnitude(vector<double> &v);
// ***********************************************************
// Precondition: Function recieves a vector of values
// PostCondition: Function returns the average magnitude of the values given in the vector
// Summary: Function calculates the average magnitude from the values within given array
// ***********************************************************
double avgPower(vector<double> &v);
// ***********************************************************
// Precondition: Function recieves a vector of values
// PostCondition: Function returns the average power of the values given in vector
// Summary: Function calculates the average power from the values within given array
// ***********************************************************
double stdDev(vector<double> &v);
// ***********************************************************
// Precondition: Function recieves a vector of values
// PostCondition: Function returns the standard deviation
// Summary: Function calculates the standard deviation from the values within given the vector and returns the standard deviation

// ***********************************************************
double var(vector<double> &v);
// ***********************************************************
// Precondition: Function recieves a vector of values
// PostCondition: Function returns the variation
// Summary: Function calculates the variation from the values within given vector
// ***********************************************************
int zeroCross(vector<double> &v1);
// ***********************************************************
// Precondition: Function recieves a vector
// Postcondition: Function finds the number of zero crossings
// Summary: Function is used to find the number of zero crossings and returns the number of zero crossings
// ***********************************************************
void percentDiff(double valA, double valB, ofstream& foutA);
// ***********************************************************
// Precondition: Function recieves two values and an output file reference
// PostCondition: Function writes the percent difference to the ouput file passed in
// Summary: Function takes 2 values and calculates the percent difference and writes to the output file
// ***********************************************************

int main()
{
    ifstream finA, finB;
    ofstream fout;
    vector<double> dataA;
    vector<double> dataB;
    double val, avg1, avg2, avgMag1, avgMag2, avgPow1, avgPow2;
    double standDev1, standDev2, var1, var2;
    int zrCross1, zrCross2;

    finA.open("two_a.txt");
    finB.open("two_b.txt");
    fout.open("comparison.txt");

    if( (finA.fail()) || (finB.fail()) )
    {
        cout << "Error opening input file\n";
        exit(1);
    }
    if(fout.fail())
    {
        cout << "Error opening output file\n";
        exit(1);
    }


    while(finA >> val)
    {
        dataA.push_back(val);
    }
    while(finB >> val)
    {
        dataB.push_back(val);
    }

    fout << "Team Members: Nicholas Rosas, Raquel Figueroa, Thalia Villalobos, Maria Garcia\n\n" << endl;

    avg1 = mean(dataA);
    avg2 = mean(dataB);
    avgMag1 = avgMagnitude(dataA);
    avgMag2 = avgMagnitude(dataB);
    avgPow1 = avgPower(dataA);
    avgPow2 = avgPower(dataB);
    standDev1 = stdDev(dataA);
    standDev2 = stdDev(dataB);
    var1 = var(dataA);
    var2 = var(dataB);
    zrCross1 = zeroCross(dataA);
    zrCross2 = zeroCross(dataB);
//a
    fout << setw(25) << "two_a.txt" << setw(20) << "two_b.txt" << setw(20) << "\% difference" << endl;
    fout << "Mean " << setw(20) << avg1 << setw(20) << avg2;
    percentDiff(avg1, avg2, fout);
    fout << "Avg Magnitude " << setw(11) << avgMag1 << setw(20) << avgMag2;
    percentDiff(avgMag1, avgMag2, fout);
    fout << "Avg Power " << setw(15) << avgPow1 << setw(20) << avgPow2;
    percentDiff(avgPow1, avgPow2, fout);
    fout << "Std Deviation " << setw(11) << standDev1 << setw(20) << standDev2;
    percentDiff(standDev1, standDev2, fout);
    fout << "Variation " << setw(15) << var1 << setw(20) << var2;
    percentDiff(var1, var2, fout);
    fout << "Zero Cross " << setw(14) << zrCross1 << setw(20) << zrCross2;
    percentDiff(zrCross1, zrCross2, fout);

/*
    b. The values for each data file that match closely are the values for standard deviation with only a percent difference of 3.60%.
    c. The values for each data file that match differently are the values for zero crossing with a percent of 49.44%.
    d. Other statistical measures that can be used are mode and or median. Mode is finding the most common value while median finds the middle value.
    e. After analyzing the data, we concluded that the sound recording were from a different person due to the high percent difference in zero crossings and the moderately high percent difference in average magnitude.
*/
    
    finA.close();
    finB.close();
    fout.close();

    return 0;
}
double mean(vector<double> &v)
{
    double sum = 0.0;
    double average;

    for (int i =0; i < v.size(); i ++)
    {
        sum += v[i];
    }
    average = sum / v.size();
    return average;  
} // end of mean()

double avgMagnitude(vector<double> &v)
{
    double sum = 0.0;
    double magnitude;

    for (int i = 0; i < v.size(); i++)
    {
        sum = fabs (sum + (v[i]));
    }
     magnitude = sum / v.size();
    return magnitude;
} // end avgMagnitude()

double avgPower(vector<double> &v)
{
    double aPower;
    double sum=0.0;

    for (int ix=0; ix < v.size(); ix++)
    {
          sum+= pow(v[ix], 2);
    }

    aPower = sum/v.size();

    return aPower;
} // end avgPower()
double var(vector<double> &v)
{
    double sum(0.0), ave(0.0), vSum(0.0), result(0.0);
    
    for (int i = 0; i < v.size(); i++)
    {
        sum += v[i];
    }
    ave = sum / v.size();
    
    for (int i = 0; i < v.size(); i++)
    {
        vSum += pow((v[i] - ave), 2);
    }
    result = vSum / v.size();

    return result;
    
}//End var()
double stdDev(vector<double> &v)
{
    double result(0);
    result = sqrt(var(v));
    return result;
}//End stdDev()
void percentDiff(double valA, double valB, ofstream& foutA)
{
    double result;

    result = ( fabs(valA - valB) / ( (valA + valB) / 2) ) * 100;

    foutA.setf(ios::fixed);
    foutA.setf(ios::showpoint);
    foutA.precision(2);
    foutA << setw(15) << result << " %" << endl;
    foutA.precision(6);

} // end of percent diff

int zeroCross(vector<double> &v1)
{
    int counter = 0;
    for(int ix = 1; ix < v1.size(); ix++)
    {
        if( ((v1[ix - 1] > 0) && (v1[ix] < 0)) || ((v1[ix - 1] < 0) && (v1[ix] > 0)) )
        {
            counter++;
        }
    }
    return counter;

} // end of zeroCross()

