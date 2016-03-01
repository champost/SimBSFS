/***************************************************************************
© Champak Beeravolu Reddy 2015-now

champak.br@gmail.com

This software is governed by the CeCILL license under French law and
abiding by the rules of distribution of free software.  You can  use,
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability.

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or
data to be ensured and,  more generally, to use and operate it in the
same conditions as regards security.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.

***************************************************************************/

/*
 * main.cpp
 *
 *  Created on: 3 Jul 2015
 *      Author: champost
 */


#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <map>
#include <vector>
#include <algorithm>
#include <numeric>

#include "MersenneTwister.h"
#include "main.h"

using namespace std;

//************ EXTERN **************
int brClass, mutClass, foldBrClass, allBrClasses;
//**********************************

MTRand rMT;
map<vector<int>, int> intVec2BrConfig;
vector<int> npopVec;
map<vector<int>, int> finalTableMap;

int ms_argc;
char **ms_argv;
int ntrees;


double ranMT() { return(rMT()); }

void setMutConfigCount(int *singleConfig) {
        vector<int> singleConfigVec(singleConfig, singleConfig+brClass);
        ++finalTableMap[singleConfigVec];
}


string getMutConfigStr(vector<int> configVec) {
	int size = configVec.size();
	stringstream stst;
	stst << "(";
	for (int j = 0; j < size-1; j++)
		stst << configVec[j] << ",";
	stst << configVec[size-1] << ")";

	return stst.str();
}


int getBrConfigNum(int *brConfVec) {
	vector<int> vec(brConfVec, brConfVec+npopVec.size());
	return intVec2BrConfig[vec];
}


//	conversion from decimal to base-(maxPopSize+1)
void evalBranchConfigs() {
	int quo, rem, maxPopSize, totPopSum, count = 0, sumConfig;
	bool skipConfig;
	maxPopSize = totPopSum = npopVec[0];

	for (size_t i = 1; i < npopVec.size(); i++) {
		totPopSum += npopVec[i];
		if (npopVec[i] > maxPopSize)
			maxPopSize = npopVec[i];
	}
	++maxPopSize;

	for (unsigned long int i = 1; i <= (unsigned long int) pow(maxPopSize,npopVec.size()); i++) {
		quo = i;
		rem = 0;
//		stringstream stst;
//		stst << ")";
		sumConfig = 0;
		skipConfig = false;
		vector<int> vec;

		for (size_t j = 0; j < npopVec.size(); j++) {
			if (quo) {
				rem = quo % (maxPopSize);
				quo /= (maxPopSize);
				if (rem > npopVec[npopVec.size()-1-j]) {
					skipConfig = true;
					break;
				}
//				stst << rem;
				sumConfig += rem;
				vec.push_back(rem);
			}
			else {
//				stst << "0";
				vec.push_back(0);
			}

//			if (j < npopVec.size() - 1)
//				stst << ",";
		}

		if (sumConfig == totPopSum)
			break;

		if (!skipConfig) {
//			stst << "(";
//			string config = stst.str();
//			reverse(config.begin(),config.end());
			reverse(vec.begin(),vec.end());
			intVec2BrConfig[vec] = count;
			++count;
//			printf("%d\t%d\t%s\n", i, count, config.c_str());
//			printf("%d\t%s\n", count, config.c_str());
		}
	}
}


void readPopSizes(int npops) {
	string line;
	ifstream ifs("popconfig.txt",ios::in);
	getline(ifs,line);
	stringstream stst;
	stst << line;
	for (int i = 0; i < npops; i++) {
		int tmp;
		stst >> tmp;
		npopVec.push_back(tmp);
	}
	ifs.close();
}


int main(int argc, char* argv[]) {

//	int nsam = atoi(argv[1]);
	int kmax = atoi(argv[argc-3]), npopSize = atoi(argv[argc-2]);
	char brFold = argv[argc-1][0];
	ntrees = atoi(argv[2]);

	readPopSizes(npopSize);

	brClass = npopVec[0]+1;
	for (size_t i = 1; i < npopVec.size(); i++)
		brClass *= (npopVec[i]+1);
	brClass -= 2;
	allBrClasses = brClass;
	foldBrClass = 0;

	//	fold the branch classes
	if (brFold == 'f') {
		if (brClass % 2)
			brClass = (brClass+1) / 2;
		else
			brClass = brClass / 2;
		foldBrClass = 1;
	}

	if (kmax == 0)
		mutClass = 0;
	else
		mutClass = kmax+2;

//	printf("Total number of (possible) configs: %.0f\n", pow(mutClass, brClass));
//	exit(-1);

	evalBranchConfigs();

	// calling ms
	ms_argc = argc - 3;
	ms_argv = argv;
	main_ms(ms_argc, ms_argv);

	for (map<vector<int>, int>::iterator it = finalTableMap.begin(); it != finalTableMap.end(); it++)
		printf("%s : %.5e\n", getMutConfigStr(it->first).c_str(), (double) it->second/ntrees);

	return 0;
}



