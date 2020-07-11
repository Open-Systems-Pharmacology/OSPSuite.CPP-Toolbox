/*
 * Interface.h
 *
 *  Created on: 06.04.2017
 *      Author: GGKMT
 */

#ifndef SRC_INTERFACE_H_
#define SRC_INTERFACE_H_

void setOutputPtr(int(*)(const char*, ...));
unsigned int initSimulation(const unsigned int, const char*);
void freeSimulation   (const unsigned int);
void exportSimulation (const char*, const char*, const char*, unsigned int, const double*, unsigned int, const double*, bool);
bool checkSimulation(const unsigned int index, const char *xmlFile, unsigned int N_obs, const double *observer, unsigned int N_par, const double *parameter);
bool getDim           (const unsigned int, unsigned int*, unsigned int*, unsigned int*, unsigned int*);
void getTimePoints    (const unsigned int, double*);
void setTimePoints    (const unsigned int, const unsigned int, const double*);
void setAllParameters (const unsigned int, const double*);
void getAllParameters (const unsigned int, double*);
void setParameters    (const unsigned int, const unsigned int, const double*, const double*, const bool);
void getParameterIndex(const unsigned int, const unsigned int, const double*, double* );
void getObserverIndex (const unsigned int, const unsigned int, const double*, double* );
void getStateIndex    (const unsigned int, const unsigned int, const double*, double* );
void setOptions       (const unsigned int, const double, const double, const double, const double, const double, const long int);
void getJacInfo       (const unsigned int, int*, const int**, const int**);
void getJacSparse     (const unsigned int, const double*, double*);
void getJacDense      (const unsigned int, const double*, double*);
void getJacDenseFD    (const unsigned int, const double*, double*);
void getRHS           (const unsigned int, const double*, double*);
bool simulate         (const unsigned int, double*, double*, double*);
void parallelSim      (const unsigned int, const unsigned int, bool*, double*, double*);
void initThreads      (const unsigned int);
void freeThreads      ();

#endif /* SRC_INTERFACE_H_ */
