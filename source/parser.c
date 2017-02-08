/*
 * main.c
 *
 *  Created on: 12 gen 2017
 *      Author: peppenero
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "parser.h"




void parse(particellaT* particelle){

	FILE *particlesFile;

	char line[200];
	char *token;
	int n=0;

	particellaT particella;

	particlesFile = fopen("Particles.dat","r");

	while (fgets(line, 200, particlesFile) != NULL) {
		if(line[0]!='#'){
			int campo=0;

			token = strtok(line," ,");

			while(token != NULL){
				switch (campo) {
				case 0:
					particella.posizione[0]=atof(token);
					campo++;
					break;
				case 1:
					token = strtok(NULL, " ,");
					particella.posizione[1]=atof(token);
					campo++;
					break;
				case 2:
					token = strtok(NULL, " ,");
					particella.posizione[2]=atof(token);
					campo++;
					break;
				case 3:
					token = strtok(NULL, " ,");
					particella.posizione[3]=atof(token);
					campo++;
					break;
				case 4:
					token = strtok(NULL, " ,");
					particella.normale[0]=atof(token);
					campo++;
					break;
				case 5:
					token = strtok(NULL, " ,");
					particella.normale[1]=atof(token);
					campo++;
					break;
				case 6:
					token = strtok(NULL, " ,");
					particella.normale[2]=atof(token);
					campo++;
					break;
				case 7:
					token = strtok(NULL, " ,");
					particella.normale[3]=atof(token);
					campo++;
					break;
				case 8:
					token = strtok(NULL, " ,");
					particella.velocita[0]=atof(token);
					campo++;
					break;
				case 9:
					token = strtok(NULL, " ,");
					particella.velocita[1]=atof(token);
					campo++;
					break;
				case 10:
					token = strtok(NULL, " ,");
					particella.velocita[2]=atof(token);
					campo++;
					break;
				case 11:
					token = strtok(NULL, " ,");
					particella.velocita[3]=atof(token);
					campo++;
					break;
				case 12:
					token = strtok(NULL, " ,");
					particella.accelerazione[0]=atof(token);
					campo++;
					break;
				case 13:
					token = strtok(NULL, " ,");
					particella.accelerazione[1]=atof(token);
					campo++;
					break;
				case 14:
					token = strtok(NULL, " ,");
					particella.accelerazione[2]=atof(token);
					campo++;
					break;
				case 15:
					token = strtok(NULL, " ,");
					particella.accelerazione[3]=atof(token);
					campo++;
					break;
				case 16:
					token = strtok(NULL, " ,");
					particella.rho=atof(token);
					campo++;
					break;
				case 17:
					token = strtok(NULL, " ,");
					particella.drhodt=atof(token);
					campo++;
					break;
				case 18:
					token = strtok(NULL, " ,");
					particella.massa=atof(token);
					campo++;
					break;
				case 19:
					token = strtok(NULL, " ,");

					particella.imove=atoi(token);
					campo=0;
					token=strtok(NULL, " ,");
					particelle[n]=particella;
					n++;
					break;
				}
			}
		}

	}
}

