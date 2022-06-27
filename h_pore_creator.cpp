/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *    This program reads an xyz file containing a BN sheet,
 *    creates an n-sided polygon shaped pore at the centre,
 *    and hydrogenizes the dangling bonds at the pore.
 *    % to compile: g++ h_pore_creator.cpp
 *    % to run: ./a.out
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

const float PI = 3.1415;
const double ONLINE = 11111;

double nline_sfinder(vector<double> slope, vector<double> con, double x, double y);
void create_polygon(vector<double> &slope, vector<double> &con, int n, double rad);
void find_neighbours(int i, int size, vector<double> &x, vector<double> &y, vector<double> &z, vector<int> &neigh);

int main()
{
  string line, type, typeName;
  char l;
  double x, y, z, chg, yLo, yHi, zLo, zHi, vol, sigma, ;
  int atoms, typ, i, j, k, id, mol, mn;
  int nAtomsTot = 0, nWallP = 0, nWall = 0, nSheet = 0, nSheet2 = 0;
  double refForce = 6.9477E-11;

  // **** THINGS TO EDIT *******

  mn = 11; // nanotube chirality to calc. radius

  yLo = 0.71;
  zLo = 1.84463;

  yHi = 80.23;
  zHi = 80.549;

  double buffer = 2.0;                      // buffer ~ bond length
  double radius = 0.5 * 3 * 1.42 * mn / PI; // Nanotube radius
  double radiusSquared = (radius + buffer) * (radius + buffer);
  int n_poly = 6; // n sided polygon for pore
  double Ly = yHi - yLo;
  double Lz = zHi - zLo;
  double yMid = (yLo + yHi) * 0.5;
  double zMid = (zLo + zHi) * 0.5;

  sigma = 3.19;

  vector<double> xG;
  vector<double> yG;
  vector<double> zG;
  vector<int> typs;
  vector<double> slope;
  vector<double> con;
  double r_poly, int_val, internal;

  /****************************************************************
  Pore creation section: create pore while maintaining B:N parity
  ****************************************************************/
  int bdel = 0;
  int ndel = 0;
  bool balance = false;
  bool bfailed = false; // variable to check if balancing failed

  while (!balance)
  {

    ifstream RSheet("bnpsheet.xyz"); // input existing sheet
    RSheet >> nSheet;
    getline(RSheet, line);
    getline(RSheet, line);

    // calculate equation of lines describing polygon
    r_poly = sqrt(radiusSquared);
    create_polygon(slope, con, n_poly, r_poly);

    // set value to compare points inside polygon
    internal = nline_sfinder(slope, con, 0, 0);

    bdel = 0;
    ndel = 0;
    for (i = 0; i < nSheet; i++)
    {
      RSheet >> type >> x >> y >> z;

      // if atom outside polygon, then keep else delete

      int_val = nline_sfinder(slope, con, y - yMid, z - zMid);

      if (int_val != internal || int_val == ONLINE)
      {
        xG.push_back(x);
        yG.push_back(y);
        zG.push_back(z);

        if (type == "B")
        {
          typs.push_back(3);
        }
        else if (type == "N")
        {
          typs.push_back(4);
        }
      }
      else
      {
        // if atoms not added to vectors they are automatically deleted
        if (type == "B")
          bdel += 1;
        if (type == "N")
          ndel += 1;
      }
    }

    slope.clear();
    con.clear();
    RSheet.close();

    if (bdel == ndel || bfailed) // conditon to stop while loop, if sheet already balanced or balancing failed
    {
      balance = true;
      cout << "Boron removed = " << bdel << '\t' << "Nitrogen removed = " << ndel << endl;
    }
    else
    {
      cout << "Boron Nitrogen imbalance, increasing buffer to " << buffer + 0.1 << " angstrom" << endl;
      buffer += 0.1;
      radiusSquared = (radius + buffer) * (radius + buffer);

      xG.clear(); // clearing the vectors, to be read in again
      yG.clear();
      zG.clear();
      typs.clear();

      if (buffer > 3.1) // if buffer exceeds 3.1 angstroms and still not balanced, balancing failed
      {
        cout << "failed to balance B and N. Manually set charges zero." << endl;
        cout << "resetting buffer to 1.5" << endl;
        buffer = 1.5;
        radiusSquared = (radius + buffer) * (radius + buffer);
        bfailed = true;
      }
    } // balancing will continue with new buffer size
  }

  /****************************************************************
  Pore hydrogenation: remove dangling bonds by hydrogenation
  ****************************************************************/

  buffer = 4;
  radiusSquared = (r_poly + buffer) * (r_poly + buffer);
  int nSheetNew = xG.size();
  bdel = 0;
  ndel = 0;
  vector<int> neighbours; // store neighbour list transiently for jth particle
  int nH_added = 0;
  int num_neigh, aquadrant, cquadrant;
  double neigh1, neigh2, neigh1y, neigh1z, neigh2y, neigh2z, theta_a, theta_c, theta_calc, theta_act, phi, Dy, Dz;

  // find pore edge atom, if dangling bond then hydrogenise

  for (int j = 0; j < nSheetNew; j++)
  {
    if (((yG[j] - yMid) * (yG[j] - yMid) + (zG[j] - zMid) * (zG[j] - zMid)) < (radiusSquared))
    {
      find_neighbours(j, nSheetNew, xG, yG, zG, neighbours);
      num_neigh = neighbours.size();
      if (num_neigh < 2 && typs[j] != 5)
        typs[j] = 5;
      else if (num_neigh == 2 && typs[j] != 5)
      {
        // estimate hydrogen atom position
        neigh1 = neighbours[0];
        neigh2 = neighbours[1];
        neigh1y = yG[neigh1] - yG[j];
        neigh1z = zG[neigh1] - zG[j];
        neigh2y = yG[neigh2] - yG[j];
        neigh2z = zG[neigh2] - zG[j];
        theta_a = atan2(neigh1z, neigh1y);
        theta_c = atan2(neigh2z, neigh2y);
        if (theta_a >= 0)
          aquadrant = ceil(theta_a / (PI * 0.5));
        else
          aquadrant = ceil((2 * PI + theta_a) / (PI * 0.5));

        if (theta_c >= 0)
          cquadrant = ceil(theta_c / (PI * 0.5));
        else
          cquadrant = ceil((2 * PI + theta_c) / (PI * 0.5));
        if (aquadrant == 0)
          aquadrant = 1;
        if (cquadrant == 0)
          cquadrant = 1;
        theta_calc = acos((neigh1y * neigh2y + neigh1z * neigh2z) / (sqrt(neigh1y * neigh1y + neigh1z * neigh1z) * sqrt(neigh2z * neigh2z + neigh2y * neigh2y)));
        if ((aquadrant * cquadrant) != 4 && (aquadrant * cquadrant) != 3)
        {
          if (theta_a < 0)
            theta_a = 2 * PI + theta_a;
          if (theta_c < 0)
            theta_c = 2 * PI + theta_c;
        }
        phi = min(theta_a, theta_c) + theta_calc / 2;
        Dy = -cos(phi) + yG[j];
        Dz = -sin(phi) + zG[j];
        xG.push_back(xG[j]);
        yG.push_back(Dy);
        zG.push_back(Dz);
        typs.push_back(5);
        nH_added += 1;
      }

      neighbours.clear();
    }
  }

  // print new sheet
  ofstream file("hbn_sheet.xyz", ios::out);
  nSheetNew += nH_added;
  file << nSheetNew << endl;
  file << endl;
  for (int j2 = 0; j2 < nSheetNew; j2++)
  {
    file << j2 + 1 << '\t' << typs[j2] << '\t' << xG[j2] << '\t' << yG[j2] << '\t' << zG[j2] << endl;
  }

  return 0;
}

/****************************************************************
   Find neighbours to check presence of dangling bonds
****************************************************************/
void find_neighbours(int i, int size, vector<double> &x, vector<double> &y, vector<double> &z, vector<int> &neigh)
{
  double distx, disty, distz, dist;
  for (int j = 0; j < size; j++)
  {
    if (j != i)
    {
      distx = x[i] - x[j];
      disty = y[i] - y[j];
      distz = z[i] - z[j];
      dist = sqrt(disty * disty + distz * distz);
      if (dist < 1.44)
      {
        neigh.push_back(j);
      }
    }
  }
}

/****************************************************************
    Generate line equations of n sided polygen
****************************************************************/
void create_polygon(vector<double> &slope, vector<double> &con, int n, double rad)
{
  vector<double> x;
  vector<double> y;
  double theta, nj, x1, x2, y1, y2, m, c;
  for (int i = 0; i < n; i++)
  {
    theta = i * 2 * PI / n;
    x.push_back(rad * cos(theta));
    y.push_back(rad * sin(theta));
  }
  for (int j = 0; j < n; j++)
  {
    if (j != n - 1)
      nj = j + 1;
    else
      nj = 0;
    x1 = x[j];
    x2 = x[nj];
    y1 = y[j];
    y2 = y[nj];
    m = (y2 - y1) / (x2 - x1);
    c = -1 * m * x1 + y1;
    slope.push_back(m);
    con.push_back(c);
  }
}

/****************************************************************
    Find whether atom lies inside, on or outside polygon
****************************************************************/
double nline_sfinder(vector<double> slope, vector<double> con, double x, double y)
{
  int n = slope.size();
  double val, num = 0;
  for (int i = 0; i < n; i++)
  {
    val = slope[i] * x + con[i] - y;
    if (val < 0)
      val = -1;
    else if (val > 0)
      val = 1;
    else
      return ONLINE;
    num += (i + 1) * val;
  }
  return num;
}
