/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *    This program creates a 2-D BN sheet, then rolls it
 *    to form a nanotube.
 *    % to compile: g++ bnntcreator.cpp
 *    % to run: ./a.out m n length
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <string.h>

using namespace std;
int nWallAtoms;
double Pi;

/******************************************************************
	Allocate 2-D arrays
******************************************************************/
double **double2Dmalloc(int nx, int ny)
{
	int i;
	double **ddata;

	ddata = (double **)malloc(nx * sizeof(double *) + nx * ny * sizeof(double));
	ddata[0] = (double *)(ddata + nx);

	for (i = 1; i < nx; i++)
	{
		ddata[i] = ddata[i - 1] + ny;
	}

	return ddata;
}

/******************************************************************
	Add atom to atom array
******************************************************************/
void sheetfiller(double x, double y, double z, double typ, int mol, double **sheet)
{
	sheet[mol - 1][0] = x;
	sheet[mol - 1][1] = y;
	sheet[mol - 1][2] = z;
	sheet[mol - 1][3] = typ;
}

/******************************************************************
	Shift all atoms in atom array
******************************************************************/
void sheetshifter(double shift, int dim, double **sheet)
{
	for (int i = 0; i < nWallAtoms; i++)
		sheet[i][dim] = sheet[i][dim] - shift;
}

/******************************************************************
	Rotate 2-D atom sheet in-plane
******************************************************************/
void rotatesheet(double theta, double **sheet)
{
	double cx, cy;
	for (int i = 0; i < nWallAtoms; i++)
	{
		cx = sheet[i][0];
		cy = sheet[i][1];
		sheet[i][0] = cx * cos(theta) + cy * sin(theta);
		sheet[i][1] = -1 * cx * sin(theta) + cy * cos(theta);
	}
}

/******************************************************************
	Trim excess atoms in sheet not required for nanotube
******************************************************************/
double **trimsheet(double length, double dia, double **sheet)
{

	double miny = 1000;
	for (int ii = 0; ii < nWallAtoms; ii++)
	{
		if (sheet[ii][1] > 0)
		{
			if (sheet[ii][1] < miny)
				miny = sheet[ii][1];
		}
	}
	double t = 1.42 * sqrt(3) / 2;
	int ncell = floor(length / t);
	cout << "n cell = " << ncell << endl;
	if (ncell % 2 == 0)
		ncell -= 1;

	length = ncell * t;
	length = length + miny + 0.3;
	cout << "min y" << miny << endl;
	cout << "act length" << (length - miny) << endl;
	cout << "max edge = " << length << endl;

	double **tsheet = double2Dmalloc(nWallAtoms + 10, 4);
	int count = 0;
	cout << "trim sheet called" << endl;
	for (int i = 0; i < nWallAtoms; i++)
	{
		if (sheet[i][0] <= (Pi * dia + 0.1) && sheet[i][0] >= (0 - 0.1))
		{
			if (sheet[i][1] <= length && sheet[i][1] >= 0)
			{
				count += 1;
				sheetfiller(sheet[i][0], sheet[i][1], sheet[i][2], sheet[i][3], count, tsheet);
			}
		}
	}
	free(sheet);
	nWallAtoms = count;
	sheetshifter(miny, 1, tsheet);
	return tsheet;
}

/******************************************************************
	Roll sheet to form nanotube
******************************************************************/
void rollsheet(double dia, double **sheet)
{
	double phi = 0;
	double minx = 100;

	for (int j = 0; j < nWallAtoms; j++)
	{
		if (sheet[j][0] < minx)
			minx = sheet[j][0];
	}
	sheetshifter(minx, 0, sheet);

	for (int i = 0; i < nWallAtoms; i++)
	{
		phi = sheet[i][0] * 2 / dia;
		sheet[i][0] = 0.5 * dia * sin(phi);
		sheet[i][2] = 0.5 * dia * cos(phi);
	}
}

/******************************************************************
	Remove overlapping atoms due to rolling at the seam
******************************************************************/
double **removeOverlaps(double olim, double **sheet)
{
	double **tcnt = double2Dmalloc(nWallAtoms + 10, 5);
	double **cnt = double2Dmalloc(nWallAtoms + 10, 4);
	double xdist, ydist, zdist, dist;
	int count = 0;
	tcnt[0][4] = 1;
	for (int i = 0; i < nWallAtoms; i++)
	{
		sheetfiller(sheet[i][0], sheet[i][1], sheet[i][2], sheet[i][3], i + 1, tcnt);
		for (int j = i + 1; j < nWallAtoms; j++)
		{
			xdist = sheet[i][0] - sheet[j][0];
			zdist = sheet[i][2] - sheet[j][2];
			ydist = sheet[i][1] - sheet[j][1];
			dist = sqrt(pow(xdist, 2) + pow(zdist, 2) + pow(ydist, 2));
			if (dist > (1.4) && tcnt[j][4] != -1)
				tcnt[j][4] = 1;
			else
				tcnt[j][4] = -1;
		}
	}
	free(sheet);

	for (int s = 0; s < nWallAtoms; s++)
	{
		if (tcnt[s][4] == 1)
		{
			count += 1;
			sheetfiller(tcnt[s][0], tcnt[s][1], tcnt[s][2], sheet[s][3], count, cnt);
		}
	}

	nWallAtoms = count;
	free(tcnt);
	return cnt;
}

/******************************************************************
	Print atom array into file
******************************************************************/
void printsheet(const char *filename, double **sheet)
{
	ofstream printer(filename);
	printer << nWallAtoms << '\n'
			<< endl;
	for (int i = 0; i < nWallAtoms; i++)
		printer << sheet[i][3] << '\t' << sheet[i][0] << '\t' << sheet[i][1] << '\t' << sheet[i][2] << endl;

	printer.close();
}

int main(int argc, char *argv[])
{

	string line, type;
	char l;
	double x, y, z, chg, xlo, xhi, ylo, yhi, zlo, zhi;
	int atoms, typ, i, j, k, id, mol, nWater, bID, aID, nAtomsTot;
	Pi = 3.141592653589793238463;

	// parse command line arguments
	if (argc != 4)
	{
		cout << "Usage = ./executable m n length" << endl;
		return 0;
	}
	double length = strtod(argv[3], NULL);
	double m, n;
	m = strtod(argv[1], NULL);
	n = strtod(argv[2], NULL);
	cout << "input arguments = " << m << " " << n << " " << length << endl;
	mol = 0;

	/******************************************************************
		Create 2-D boron nitride sheet
	******************************************************************/

	double xLo = 0;
	double yLo = 0;
	double zLo = 0;

	double xHi = 60;
	double yHi = 200;
	double zHi = 81;

	double b = 1.42; // bond length Angstroms
	double theta = 120 * Pi / 180.0;
	double c = b * sin(0.5 * theta);
	double d = b * cos(0.5 * theta);

	double D = sqrt(3.0) * b * sqrt(m * m + n * n + m * n) / Pi;
	double ctheta = atan(sqrt(3.0) * n / (2 * m + n));
	cout << "Dia " << D << '\n'
		 << "ctheta " << ctheta << endl;
	double lxright = Pi * D * cos(ctheta) + sqrt(3) * b;
	double lxleft = length * sin(ctheta) + 3 * b;
	double Lx = lxleft + lxright;
	double Ly = length * cos(ctheta) + Pi * D * sin(ctheta) + 4 * b;
	double Lz = zHi - zLo;

	cout << "initial Lx = " << Lx << endl;
	cout << "initial Ly = " << Ly << endl;

	int nX = int(Lx / (2 * c));
	Lx = double(nX) * 2 * c;

	int nY = int(Ly / (2 * d + b));
	Ly = double(nY) * (2 * d + b);

	cout << "new Lx = " << Lx << endl;
	cout << "new Ly = " << Ly << endl;

	xHi = Lx;
	yHi = Ly;

	nWallAtoms = nX * nY * 4;
	nAtomsTot = nWallAtoms;

	double Y;
	typ = 3;
	int nCountWallAtoms = 0;
	double delta = 3.35;

	mol = 1;
	double **sheet = double2Dmalloc(nWallAtoms + 10, 4);
	for (i = 0; i < nWallAtoms; i++)
	{
		sheet[i][0] = 0;
		sheet[i][1] = 0;
		sheet[i][2] = 0;
	}
	cout << "init success" << endl;

	Y = 0;
	for (i = 0; i < nX; i++)
	{
		for (j = 0; j < nY; j++)
		{
			mol++;
			typ = 3;
			x = 0.0 + 2 * c * i;
			z = Y;
			y = 0.0 + (2 * d + 2 * b) * j;
			sheetfiller(x, y, z, typ, mol, sheet);

			mol++;
			typ = 4;
			x = c + 2 * c * i;
			z = Y;
			y = d + (2 * d + 2 * b) * j;
			sheetfiller(x, y, z, typ, mol, sheet);

			mol++;
			typ = 3;
			x = c + 2 * c * i;
			z = Y;
			y = d + b + (2 * d + 2 * b) * j;
			sheetfiller(x, y, z, typ, mol, sheet);

			mol++;
			typ = 4;
			x = 0.0 + 2 * c * i;
			z = Y;
			y = 2 * d + b + (2 * d + 2 * b) * j;
			sheetfiller(x, y, z, typ, mol, sheet);

			nCountWallAtoms += 4;
		}
	}

	printsheet("origsheet", sheet);

	cout << "shifting sheet" << endl;
	sheetshifter(lxleft, 0, sheet);
	printsheet("shiftedsheet", sheet);

	cout << "rotating sheet" << endl;
	rotatesheet(ctheta, sheet);
	printsheet("rotatedsheet", sheet);

	cout << "trimming sheet" << endl;
	double **trimmedsheet = trimsheet(length, D, sheet);
	printsheet("trimmedsheet", trimmedsheet);

	cout << "rolling sheet" << endl;
	rollsheet(D, trimmedsheet);
	printsheet("rolledsheet", trimmedsheet);

	cout << "removing overlaps" << endl;
	double **cnt = removeOverlaps(b, trimmedsheet);
	printsheet("BNNT", cnt);

	return 0;
}
