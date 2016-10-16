// Test.cpp : Defines the entry point for the console application.
//


#include <iostream>
#include <vector>
#include <math.h>

using namespace std;

int n;

vector<double> add(vector<double>& vec1, vector<double>& vec2)
{
  vector<double> result;
  result.push_back(vec1.at(0) + vec2.at(0));
  result.push_back(vec1.at(1) + vec2.at(1));
  return result;
}

double f(double l, double s)
{
  return (fabs(s) / sqrt(l));
}

double inc(vector<double> p, vector<double> q)
{
  return ((q.at(1) - p.at(1)) / (q.at(0) - p.at(0)));
}

vector< vector<double> > computeUpperConvexHull(vector<double> * b, int length)
{
  vector< vector<double> > upperConvexHull;
  upperConvexHull.push_back(b[0]);
  if(length > 1)
    upperConvexHull.push_back(b[1]);
  
  for(int i = 2; i < length; i++)
  {
    int l = (int)upperConvexHull.size();
    double delta_1 = inc(upperConvexHull.at(l - 2), upperConvexHull.at(l - 1));
    double delta_2 = inc(upperConvexHull.at(l - 1), b[i]);

    while((delta_2 > delta_1) && (l >= 2))
    {
      upperConvexHull.pop_back();
      l--;
      if(l >= 2)
      {
        delta_2 = inc(upperConvexHull.at(l - 1), b[i]);
        delta_1 = inc(upperConvexHull.at(l - 2), upperConvexHull.at(l - 1));
      }
    }
    upperConvexHull.push_back(b[i]);
  }

  return upperConvexHull;
}

vector<double>* joining(vector< vector<double> >& p, vector< vector<double> >& q)
{
  vector<double>* t = new vector<double>[p.size() + q.size() - 1];
  int i = 1;
  int j = 1;
  t[0] = add(p.at(0), q.at(0));
  while(i+j <= (int)(p.size() + q.size() - 1))
  {
    /** The next run through the loop will define t_{i+j} */
    if(i == (int)p.size())
      j++;
    else if (j == (int)q.size())
      i++;
    /** Now test whether delta_i >= delta'_j */
    else
    {
      double delta_i = inc(p.at(i-1), p.at(i));
      double delta_j = inc(q.at(j-1), q.at(j));
      if(delta_i >= delta_j)
        i++;
      else
        j++;
    }

    t[i+j - 2] = add(p.at(i-1), q.at(j-1));
  }

  return t;
}



double maxSubsequence(double* a, int left, int right)
{
  /** Trivial case: length is 1 and sum is a_i */
  if(left==right)
    return f(1, a[left-1]);
  
  int middle = (int)(right + left - 1)/2;
  int length = middle - left + 1;
  /** CORRECTED*/
  int length2 = right - middle;
  /** CORRECTED*/
  double opt_left = maxSubsequence(a, left, middle);
  double opt_right = maxSubsequence(a, middle+1, right);
  
  /** Now treat the crossing intervals. First treat the left half.
   * Thus, b_k = (k, sum(middle - k + 1, middle)),
   * i.e., the first coordinate of b_k is k and its second coordinate is the sum
   * of the interval that ends with array index middle and has length k.
   */
  vector<double>* b = new vector<double>[length];
  b[0].push_back(1);
  b[0].push_back(a[middle - 1]);
  for(int i = 1; i < length; i++)
  {
    b[i].push_back(b[i-1].at(0) + 1);
    b[i].push_back(b[i-1].at(1) + a[middle - 1 - i]);
  }
  
  /** Now the right half.
   * I.e., the second coordinate of b2_k is the sum of the interval that
   * begins with array index middel + 1 and has length k.
   */
  vector<double> * b2 = new vector<double>[length2];
  b2[0].push_back(1);
  b2[0].push_back(a[middle]);
  for(int i = 1; i < length2; i++)
  {
    b2[i].push_back(b2[i-1].at(0) + 1);
    b2[i].push_back(b2[i-1].at(1) + a[middle + i]);
  }
  
  /** Compute the upper convex hull p_1, ..., p_r of b_1, ..., b_length and
   * the upper convex hull q_1, .., q_s of b2_1, ..., b2_length
   */
  vector< vector<double> > p = computeUpperConvexHull(b, length);
  vector< vector<double> > q = computeUpperConvexHull(b2, length2);
  
  /** Joint the two concave sequences p_1, ..., p_r and q_1, ..., q_s
   * Output is a sequence t_1, ..., t_{r+s-1}
   */
  vector<double>* t = joining(p, q);

  for(int i = 0; i < length; i++)
  {
    b[i].at(1) = -b[i].at(1);
  }
  for(int i = 0; i < length2; i++)
  {
    b2[i].at(1) = -b2[i].at(1);
  }

  /** Do the same for the lower convex hull with output t'_1, ..., t'_{r'+s'-1} */
  vector< vector<double> > p2 = computeUpperConvexHull(b, length);
  vector< vector<double> > q2 = computeUpperConvexHull(b2, length2);
  vector<double>* t2 = joining(p2, q2);
  
  /** Merge the two sequences t and t2*/
  
  /* Evaluate f on all points and determine maximum */
  double opt_crossing = f(t[0].at(0), t[0].at(1));;
  for(int i = 1; i < (int)(p.size() + q.size() - 1); i++)
  {
    double r = f(t[i].at(0), t[i].at(1));
    if(opt_crossing < r)
      opt_crossing = r;
  }
  for(int i = 0; i < (int)(p2.size() + q2.size() - 1); i++)
  {
    double r = f(t2[i].at(0), t2[i].at(1));
    if (opt_crossing < r)
      opt_crossing = r;
  }

  //no memory leaks
  delete[] b;
  delete[] b2;
  delete[] t;
  delete[] t2;

  if(opt_left > opt_right)
    if(opt_left > opt_crossing)
      return opt_left;
    else
      return opt_crossing;
  else
    if(opt_right > opt_crossing)
      return opt_right;
    else
      return opt_crossing;
}

extern "C" {
void multires(double *a , double *erg, int *n) 
     {
     int right;
     right=*n;
     *erg=maxSubsequence(a, 1, right);  
     }
}

