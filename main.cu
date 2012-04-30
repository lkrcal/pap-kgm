/*
 * main.cpp
 *
 *  Created on: 29.4.2012
 *      Author: jirka
 */

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <cmath>
#include <string.h>
#include <limits.h>
#include <map>
#include "Zasobnik.h"
#include "Graf.h"

using namespace std;

#define CHECK_AMOUNT  100

int POCET_VLAKEN = 32;

// const int NEJMENSI_MOZNY_STUPEN = 2;
int pocetUzlu;
int nejmensiStupen = INT_MAX;
int * zatimNejlepsiKostra = NULL;

int * zasobnik = NULL;
int ZASOBNIK_VELIKOST; // velikost zasobniku v poctu koster
int ZASOBNIK_POCATEK;
int ZASOBNIK_POCET_KOSTER;
int VELIKOST_KOSTRY; // velikost kostry v intech

bool jeZasobnikPrazdny() {
	return ZASOBNIK_POCET_KOSTER == 0;
}

void vytvorZasobnik(int velikost) {
	ZASOBNIK_VELIKOST = velikost;
	zasobnik = new int[velikost * VELIKOST_KOSTRY];
	ZASOBNIK_POCET_KOSTER = 0;
	ZASOBNIK_POCATEK = 0;
}

void roztahniZasobnik() {
	cout << "zvetsuji pole..." << endl;
	ZASOBNIK_VELIKOST *= 2 * VELIKOST_KOSTRY;
	int * novePole = new int [ZASOBNIK_VELIKOST];
	memcpy(novePole, zasobnik + ZASOBNIK_POCATEK,
			ZASOBNIK_POCET_KOSTER * VELIKOST_KOSTRY * sizeof(int));
	delete[] zasobnik;
	zasobnik = novePole;
	ZASOBNIK_POCATEK = 0;
}

void pridejKostruNaZasobnik(int * kostra, int pocatecniUzel, int koncovyUzel) {
	cout << "stack-add: " << (*kostra) << " + [" << pocatecniUzel << ","
			<< koncovyUzel << "]" << endl;
	if (ZASOBNIK_POCET_KOSTER == ZASOBNIK_VELIKOST) {
		roztahniZasobnik();
	}
	int index = (ZASOBNIK_POCATEK + ZASOBNIK_POCET_KOSTER) % ZASOBNIK_VELIKOST;
	index *= VELIKOST_KOSTRY;
	int velikostPoleKostry = 0;
	if (kostra != NULL) {
		while (kostra[velikostPoleKostry] != -1)
			velikostPoleKostry += 2;
	}
	memcpy(zasobnik + index, kostra, velikostPoleKostry * sizeof(int));
	zasobnik[index + velikostPoleKostry] = pocatecniUzel;
	zasobnik[index + velikostPoleKostry + 1] = koncovyUzel;
	if (velikostPoleKostry + 2 < VELIKOST_KOSTRY)
		zasobnik[index + velikostPoleKostry + 2] = -1;
	++ZASOBNIK_POCET_KOSTER;
}

 __device__ int pocetUzluVKostre(int * kostra, int velikost) {
	int pocet = 0;
	for (int i = 0; i < velikost && kostra[i] != -1; i += 2) {
		pocet++;
	}
	return ++pocet;
}

 __device__ int stupenKostry(int * kostra, int velikost) {
	int pole[100];
	for (int i = 0; i < 50; i++) {
		pole[i] = 0;
	}

	for (int i = 0; i < velikost && kostra[i] != -1; i++) {
		++(pole[kostra[i]]);
	}

	int nejStupen = 0;

	for (int i = 0; i < 50; i++) {
		if (pole[i] > nejStupen)
			nejStupen++;
	}
	return nejStupen;
}

 __device__ int stupenKostry(int * kostra, int pocatekDalsiHrany, int konecDalsiHrany,
		int velikost) {
	int pole[100];
		for (int i = 0; i < 50; i++) {
			pole[i] = 0;
		}

		for (int i = 0; i < velikost && kostra[i] != -1; i++) {
			++(pole[kostra[i]]);
		}
		++pole[pocatekDalsiHrany];
		++pole[konecDalsiHrany];


		int nejStupen = 0;

		for (int i = 0; i < 50; i++) {
			if (pole[i] > nejStupen)
				nejStupen++;
		}
		return nejStupen;
}

int * nactiRadek(string s, int nodes) {
	int * radek = new int[nodes];
	istringstream is(s);
	cout << s << ": ";

	for (int i = 0; i < nodes; i++) {
		radek[i] = s[i];
		cout << radek[i] << " ";
	}

	return radek;
}

int ** nactiGraf(char * soubor) {
	string s;
	ifstream in;
	in.open(soubor, ios::binary);
	getline(in, s);
	istringstream is(s);
	is >> pocetUzlu;

	//    in >> pocetUzlu;
	cout << "pocetUzlu: " << pocetUzlu << endl;
	//    in.ignore(INT_MAX, '\n');

	const int ** maticeSousednosti = new const int *[pocetUzlu];

	for (int i = 0; i < pocetUzlu; i++) {
		getline(in, s);
		int * radek = new int[pocetUzlu];
		//        cout << s << ": ";

		for (int j = 0; j < pocetUzlu; j++) {
			if (s.at(j) == '1')
				radek[j] = 1;
			else
				radek[j] = 0;
			cout << radek[j] << " ";
		}

		maticeSousednosti[i] = radek;
		cout << endl;
	}
	in.close();

	//return maticeSousednosti;


	int ** graf;
	    graf = new int * [pocetUzlu];
	    int i, j, k, sousedu;
	    for (i = 0; i < pocetUzlu; i++) {
	        sousedu = 0;
	        for (j = 0; j < pocetUzlu; j++) {
	            if (maticeSousednosti[i][j] == 1)
	                sousedu++;
	        }
	        graf[i] = new int[sousedu + 1];
	        sousedu = 0;
	        for (k = 0; k < pocetUzlu; k++) {
	            if (maticeSousednosti[i][k] == 1)
	                graf[i][sousedu++] = k;
	        }
	        graf[i][sousedu] = -1;
	    }

	    for (int m = 0; m < pocetUzlu; m++) {
	    	delete [] maticeSousednosti[m];
	    }

	    delete [] maticeSousednosti;
	    VELIKOST_KOSTRY = (pocetUzlu - 1) * 2;
	    return graf;
}

void vytiskniKostru() {
	cout << endl;
	cout << "KOSTRA S NEJNIZSIM STUPNEM:" << endl;
	cout << "===========================" << endl;
	cout << (*zatimNejlepsiKostra) << endl;
	cout << endl;
	cout << "NEJVYSSI STUPEN TETO KOSTRY:" << endl;
	cout << "============================" << endl;
	cout << nejmensiStupen << endl;
}

int pocetVychazejicichHran(int ** graf, int uzel) {
	int pocetHran = 0;
	for (int i = 0; i < pocetUzlu; i++) {
		if (graf[uzel][i] == -1)
			return pocetHran;
		else
			pocetHran++;
	}
	return pocetHran;
}

int * sousedniUzly(int ** graf, int uzel, int& pocet) {
	pocet = pocetVychazejicichHran(graf, uzel);
	return graf[uzel];
}

bool obsahujeUzel(int * kostra, int uzel, int velikostKostry) {
	for (int i = 0; i < velikostKostry; i++) {
		if (kostra[i] == -1) {
			return false;
		}
		if (kostra[i] == uzel) {
			return true;
		}
	}
	return false;
}

int zkusExpandovatZUzlu(int ** graf, int * kostra, int uzel,
		int velikostKostry) {
	int pocetSousedu;
	int pocet = 0;
	int * sousedi = sousedniUzly(graf, uzel, pocetSousedu);
	for (int i = 0; i < pocetSousedu; i++) {
		if (!obsahujeUzel(kostra, sousedi[i], velikostKostry)) {
			pocet++;
		}
	}
	return pocet;
}

int posledniPridanyUzelVKostre(int * kostra, int velikostKostry) {
	if (kostra[0] == -1)
		return -1;
	for (int i = 2; i < velikostKostry; i += 2) {
		if (kostra[i] == -1)
			return kostra[i - 1];
	}
	return kostra[velikostKostry];
}

 __device__ int predchoziUzel(int * kostra, int uzel, int velikostKostry) {
	int i = velikostKostry;
	for (int i = 0; i < velikostKostry; i += 2) {
		if (kostra[i] == -1) {
			i -= 1;
			break;
		}
	}
	for (; i > 0; i -= 2) {
		if (kostra[i] == uzel)
			return kostra[i - 1];
	}
	return -1;
}

 __device__ int expandujZUzlu(int ** graf, int * kostra, int uzel, int * zasobnik,
		int writeIndex, int velikostZasobniku, int velikostKostry) {
	int pocet;
	int * sousedi = sousedniUzly(graf, uzel, pocet);
	int i = 0;
	for (; i < pocet; i++) {
		if (!obsahujeUzel(kostra, sousedi[i], velikostKostry)) {
			if (stupenKostry(kostra, uzel, sousedi[i], velikostKostry)
					< nejmensiStupen) { // perspektivni reseni
				int j = 0;
				for (; j < velikostKostry; j++) {
					if (kostra[j] == -1) {
						break;
					}
					zasobnik[((writeIndex + i) % velikostZasobniku)
							* velikostKostry + j] = kostra[j];
				}
				zasobnik[((writeIndex + i) % velikostZasobniku) * velikostKostry
						+ j] = uzel;
				zasobnik[(writeIndex + i % velikostZasobniku) * velikostKostry
						+ j + 1] = sousedi[i];
				if (j + 2 < velikostKostry)
					zasobnik[((writeIndex + i) % velikostZasobniku)
							* velikostKostry + j + 2] = -1;
			}
		}
	}
	return ((writeIndex + i) % velikostZasobniku);
}

__device__ int zkusExpandovat(int ** graf, int * kostra, int velikostKostry) {
	int uzel = posledniPridanyUzelVKostre(kostra, velikostKostry);
	int pocet = zkusExpandovatZUzlu(graf, kostra, uzel, velikostKostry);
	if (pocet == 0) {
		uzel = predchoziUzel(kostra, uzel, velikostKostry);
		while (uzel > -1) {
			cout << "vracim se o hranu zpet..." << endl;
			pocet += zkusExpandovatZUzlu(graf, kostra, uzel, velikostKostry);
			uzel = predchoziUzel(kostra, uzel, velikostKostry);
		}
	}
	return pocet;
}

__device__ void expanduj(int ** graf, int * kostra, int velikostKostry, int * zasobnik, int velikostZasobniku, int writeIndex) {
	int uzel = posledniPridanyUzelVKostre(kostra, velikostKostry);
	int index = expandujZUzlu(graf, kostra, uzel, zasobnik, writeIndex, velikostZasobniku, velikostKostry);
	if (index == writeIndex) {
		uzel = predchoziUzel(kostra, uzel, velikostKostry);
		while (uzel > -1) {
			cout << "vracim se o hranu zpet..." << endl;
			index = expandujZUzlu(graf, kostra, uzel, zasobnik, index, velikostZasobniku, velikostKostry);
			uzel = predchoziUzel(kostra, uzel, velikostKostry);
		}
	}
}

__global__ void kernel(int ** graf, int** zasobnik, int * pocatekZasobniku,
		int * pocetKosterNaZasobniku, int * velikostZasobniku,
		int * reseni_kostry, int * hotovo, int velikostKostry, int checkInterval) {
	int citac = 0;
	bool konec = false;
	int myIndex = threadIdx.x;
	int index;
	int pocetNovychStavu = 0;
	int * kostra;
	__shared__	int reseni_stupne[32];

	__shared__	int cache[32];

	__shared__	int kostry[32];

	kostra = kostry + myIndex;

	reseni_stupne[myIndex] = INT_MAX;

	while (!konec) {
		citac++;
		index = (((*pocatekZasobniku) + myIndex) % (*velikostZasobniku)) * velikostKostry;
		if ((*zasobnik)[index] != -1) {
			// zkopirujeme si kostru do lok. pameti
			for (int i = 0; i < velikostKostry; i++) {
				kostra[i] = (*zasobnik)[index + i];
			}
			if (pocetUzluVKostre(kostra, velikostKostry)
					== pocetUzlu) {
				int stupen = stupenKostry(kostra, velikostKostry);
				if (stupen < reseni_stupne[myIndex]) {
					reseni_stupne[myIndex] = stupen;
					for (int i = 0; i < velikostKostry; i++) {
						reseni_kostry[myIndex * velikostKostry + i] =
								kostra[i];
						if (kostra[i] == -1)
							break;
					}
				}
				if (stupen == 2) {
					(*hotovo) = 1;
				}
				pocetNovychStavu = 0;
			} else {
				pocetNovychStavu = zkusExpandovat(graf,kostra,velikostKostry);
			}
			// odebereme kostru ze zasobniku
			(*zasobnik)[index] = -1;
		} else {
			if (myIndex == 0) {
				(*hotovo) = 1;
			}

		}

		if (citac % checkInterval == 0) {
			if ((*hotovo) == 1)
				konec = true;
		}

		cache[myIndex] = pocetNovychStavu;
		__syncthreads();
		int neaktivnich = 1;

		for (int i = 0; neaktivnich < blockDim.x; i++) {
			if (myIndex < blockDim.x - neaktivnich)
				cache[myIndex] += cache[myIndex + neaktivnich];
			__syncthreads();
			neaktivnich <<= 1;
		}

		if (myIndex == 0) {
			(*pocatekZasobniku) += blockDim.x;
			(*pocatekZasobniku) %= (*velikostZasobniku);
			(*pocetKosterNaZasobniku) = ((*pocetKosterNaZasobniku) < blockDim.x) ? 0 : (*pocetKosterNaZasobniku) - blockDim.x;



//			int staraVelikost = (*velikostZasobniku);
			if ((*pocetKosterNaZasobniku) + cache[0] > (*velikostZasobniku)) {
				cout << "pocet stavu: " << ((*pocetKosterNaZasobniku) + cache[0]) << endl;
			}
//			int * novyZas;
//			cudaMalloc ((void**) &novyZas, (*velikostZasobniku) * velikostKostry * sizeof(int));
//		    cudaMemcpy (novyZas, (*zasobnik), staraVelikost * VELIKOST_KOSTRY * sizeof(int), cudaMemcpyDeviceToDevice);
//			cudaFree(*zasobnik);
//		    (*zasobnik) = novyZas;
		}

		__syncthreads();

		int writeIndex = ((*pocetKosterNaZasobniku) + cache[myIndex])
				% ZASOBNIK_VELIKOST;

		expanduj(graf, kostra, velikostKostry, *zasobnik, *velikostZasobniku, writeIndex);

		__syncthreads();

		if (myIndex == 0) {
			(*pocetKosterNaZasobniku) += cache[0];
		}
		__syncthreads();
	}

	cache[myIndex] = reseni_stupne[myIndex];
	int step = blockDim.x / 2;

	while (step != 0) {
		if (myIndex < step)
			cache[myIndex] =
					cache[myIndex] < cache[myIndex + step] ?
							cache[myIndex] : cache[myIndex + step];
		__syncthreads();
		step /= 2;
	}

	if (myIndex == 0) {
		int nejStupen = cache[0];
		cout << "nejlepsi reseni: stupen " + nejStupen;
		for (int i = 0; i < blockDim.x; i++) {
			if (reseni_stupne[i] == nejStupen) {
				for (int j = 0; j < velikostKostry; j += 2) {
					cout << "[" << reseni_kostry[i * velikostKostry + j] << ", "
							<< reseni_kostry[i * velikostKostry + j] << "]";
					if (j < velikostKostry - 2)
						cout << " ";
				}
			}
		}
	}

}

int main(int argc, char *argv[]) {

	if (argc != 2) {
		cout << "Jako parametr zadejte nazev souboru." << endl;
		return (1);
	}

	int ** graf = nactiGraf(argv[1]);

	vytvorZasobnik(100000);

	int * reseni_kostry = new int[POCET_VLAKEN * VELIKOST_KOSTRY];

	// kopirovani grafu na GPU
	int ** g_graf;
	cudaMalloc ((void**) &g_graf, pocetUzlu * sizeof(int*));
	for (int i = 0; i < pocetUzlu; i++) {
		int size = 0;
		for (int j = 0;; j++) {
			++size;
			if (graf[i][j] == -1)
				break;
		}
		cudaMalloc ((void**) &g_graf[i], size * sizeof(int));
		cudaMemcpy (g_graf, graf, size * sizeof(int), cudaMemcpyHostToDevice);
	}

	int * sousedi = graf[0];
	if (sousedi[0] == -1) {
		cout << "Graf neni souvisly nebo neobsahuje alespon dva uzly" << endl;
		return (1);
	}
	for (int i = 0; sousedi[i] > -1; i++) {
		pridejKostruNaZasobnik(NULL, 0, sousedi[i]);
	}

	// kopirovani zasobniku na GPU
		int * g_zasobnik;
		cudaMalloc ((void**) &g_zasobnik, ZASOBNIK_VELIKOST * VELIKOST_KOSTRY * sizeof(int));
		cudaMemcpy (g_zasobnik, zasobnik, ZASOBNIK_VELIKOST * VELIKOST_KOSTRY * sizeof(int), cudaMemcpyHostToDevice);

		int * g_reseni_kostry;
				cudaMalloc ((void**) &g_reseni_kostry, POCET_VLAKEN * VELIKOST_KOSTRY * sizeof(int));


		int * g_pocatekZasobniku;
		int * g_pocetKoster;
		int * g_velikostZasobniku;
		int * hotovo;
		int * zasobnikPtr;

		cudaMalloc ((void**) &g_pocatekZasobniku, sizeof(int));
		cudaMalloc ((void**) &g_pocetKoster, sizeof(int));
		cudaMalloc ((void**) &g_velikostZasobniku, sizeof(int));
		cudaMalloc ((void**) &hotovo, sizeof(int));
		cudaMalloc ((void**) &zasobnikPtr, sizeof(int*));
		(*g_pocatekZasobniku) = 0;
		(*g_pocetKoster) = ZASOBNIK_POCET_KOSTER;
		(*g_velikostZasobniku) = ZASOBNIK_VELIKOST;
		(*hotovo) = 0;
		zasobnikPtr = zasobnik;

		kernel<<<1, POCET_VLAKEN>>> (g_graf, &zasobnikPtr, g_pocatekZasobniku, g_pocetKoster, g_velikostZasobniku,
				g_reseni_kostry, hotovo, VELIKOST_KOSTRY, 100);

}

