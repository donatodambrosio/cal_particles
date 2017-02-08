/*
 * parser.h


 *
 *  Created on: 16 gen 2017
 *      Author: giuseppe
 */

#include <OpenCAL/calCommon.h>

#ifndef PARSER_H_
#define PARSER_H_

typedef struct particella{

	CALreal posizione[4];
	CALreal normale [4];
	CALreal velocita [4];
	CALreal accelerazione[4];

	CALreal rho;
	CALreal drhodt;

	CALreal massa;

	CALint imove;

}particellaT;



void parse(particellaT*);

#endif /* PARSER_H_ */
