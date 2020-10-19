///////IMPORTANTE Para ejecutar los casos grandes, en mi ordenador, hay que hacerlo de umo en uno.
/////// De todas maneras, para ejecutar algun habrá que cambiar el main.

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <fstream>
#include <vector>
#include <algorithm>
#include <unordered_set>
#include <random>
#include <omp.h>
using namespace std;

vector< vector<double> > matrix;



vector<int> solucion, mundo, usados, aux;
unordered_set<int> sol;
vector<double> distancias_a_la_solucion, distancia_anyadida;
unsigned m;
int n, iteracion;
double scoreTotal;

void input(string file_name){
	ifstream fe(file_name);

  int i, j, cont;
  double distancia;
  vector<double> aux;

  if(fe.is_open()){
    fe >> n >> m;
    cont = 1;

    while (!fe.eof()) {
			for (int l = 0; l < n; l++){
				if (l < cont)
					aux.push_back(0);
				else{
					fe >> i >> j;
					//cout << i << " " << j << endl;
	        fe >> distancia;
	        aux.push_back(distancia);
				}
			}
      matrix.push_back(aux);
			aux.clear();
      cont++;
			if (cont == n)
				fe >> i;
    }

		for (int i = 0; i < n; i++)
			aux.push_back(0);
		matrix.push_back(aux);
  }
  fe.close();
}

void score(){
	double sum = 0;

	for (unsigned i = 0; i < solucion.size(); i++){
		for (unsigned j = 0; j < solucion.size(); j++){
			sum += matrix[solucion[i]][solucion[j]];
		}
	}

	scoreTotal = sum;
}

double distanciaAcumulada(int a){
	double suma = 0;

	for (int i = 0; i < n; i++){
		if (i < a){
			suma += matrix[i][a];
			//cout << matrix.at(i).at(a) << endl;
		}
		else{
			suma += matrix[a][i];
			//cout << matrix.at(a).at(i) << endl;
		}
	}

	return suma;
}

bool remakeScore(int posI, int j){
	double restar, sumar, desfase;
	bool cambio = false;

	restar = distancia_anyadida[posI];
	sumar = distancia_anyadida[j];
	desfase = matrix[j][posI] + matrix[posI][j];

	//cout << "sumar: " << sumar << " restar: " << restar << "  otro: " << desfase << endl << endl;

	if (sumar - desfase > restar){
		cambio = true;
		scoreTotal += sumar - desfase - restar;

		for (int i = 0; i < n; i++){
			if (distancia_anyadida[i] != 0)
				distancia_anyadida[i] += matrix[i][j] + matrix[j][i] - matrix[i][posI] - matrix[posI][i];
				//cout << "sumo al " << i << "  " << matrix[i][j] << " + " << matrix[j][i] <<
				//" y resto al " << i << "  " << matrix[i][posI] << " - " << matrix[posI][i] << endl;
		}
	}

	return cambio;
}

bool intercambio (int i, int j){
	bool cambio = false;
	int pos = -1;

	for (unsigned k = 0; k < m; k++){
		if (solucion[k] == i )
			pos = k;
	}

	if (pos != -1){
		cambio = remakeScore(i, j);

		if (cambio){
			solucion[pos] = j;
		}
	}

	else{
		cout << "Error: pos = -1" << endl;
	}

	return cambio;
}

int buscaPos(int a){
	int pos = 0;

	for (unsigned i = 0; i < mundo.size(); i++){
		if (mundo[i] == a)
			pos = i;
			//cout << mundo.at(i) << " elemento " << a << endl;
	}

	//cout << "la posicion es :" << pos << endl;

	return pos;
}

int buscaMaximo(){
	int max = 0;
	int maxD = distancias_a_la_solucion[0];

	for (unsigned i = 1; i < distancias_a_la_solucion.size(); i++){
		if (maxD < distancias_a_la_solucion[i]){
			max = i;
			maxD = distancias_a_la_solucion[i];
		}
	}

	return max;
}

int buscaMaximoAnyadido(){
	int max = 0;
	int maxD = 0;

	for (unsigned i = 0; i < mundo.size(); i++){
		//cout << endl << "escoger maximo " << mundo[i] << " distancia " << distancia_anyadida[mundo[i]] << endl;
		if (maxD < distancia_anyadida[mundo[i]]){
			max = mundo[i];
			maxD = distancia_anyadida[mundo[i]];
		}
	}

	return max;
}

void actualizaDistancias(int a){
	for (unsigned i = 0; i < distancias_a_la_solucion.size(); i++){
		distancias_a_la_solucion[i] -= matrix[i][a] + matrix[a][i];
	}
	distancias_a_la_solucion[a] = 0;
}

void destruir(){
	matrix.clear();
	solucion.clear();
	mundo.clear();
	usados.clear();
	distancias_a_la_solucion.clear();
	distancia_anyadida.clear();
	iteracion = 0;
}

void destruirLight(){
  solucion.clear();
  sol.clear();
	mundo.clear();
	usados.clear();
	distancias_a_la_solucion.clear();
	distancia_anyadida.clear();
  scoreTotal = 0;
}

void crearMundo(){
  int pos;

  for (int i = 0; i < n; i++){
    mundo.push_back(i);
  }

  for (unsigned i = 0; i < solucion.size(); i++){			//////	Tiene que ser el size()
    pos = buscaPos(solucion[i]);
    mundo.erase(mundo.begin() + pos);
  }
}

void algoritmoGreedy(){
	int max = 0, pos;
	double maxD, posibleMaxD;

	maxD = distanciaAcumulada(0);
	distancias_a_la_solucion.push_back(maxD);

	for (int i = 0; i < n; i++)
		mundo.push_back(i);

	for (unsigned i = 1; i < mundo.size(); i++){
		posibleMaxD = distanciaAcumulada(i);
		distancias_a_la_solucion.push_back(posibleMaxD);
		if (maxD < posibleMaxD){
			max = i;
			maxD = posibleMaxD;
		}
	}

	actualizaDistancias(max);

	solucion.push_back(max);

	pos = buscaPos(max);
	//cout << max << endl << endl;
	mundo.erase(mundo.begin() + pos);

	//cout << solucion.size() << endl;

	while(solucion.size() < m){
		max = buscaMaximo();

		actualizaDistancias(max);

		solucion.push_back(max);

		pos = buscaPos(max);

		mundo.erase(mundo.begin() + pos);
	}
}

void iniciaDistanciaAnyadida(){
	double sum;
	distancia_anyadida.clear();

	for (int k = 0; k < n; k++){
		sum = 0;
		for (unsigned j = 0; j < m; j++){
			//sum += matrix.at(k).at(solucion.at(j)) + matrix.at(solucion.at(j)).at(k);
			if (k < solucion[j])
				sum += matrix[k][solucion[j]];
			else
				sum += matrix[solucion[j]][k];
			//cout << matrix.at(i).at(solucion.at(j)) << " + " << matrix.at(solucion.at(j)).at(i) << endl;
		}
		distancia_anyadida.push_back(sum);
		//cout << endl;
	}
}

void busquedaLocal(int numero_vecinos){
	int max, pos;
	bool cambio;

	iniciaDistanciaAnyadida();

	for(int i = 0; i < numero_vecinos && !mundo.empty(); i++){

		max = buscaMaximoAnyadido();

		aux.clear();

		for (unsigned j = 0; j < m; j++){
			aux.push_back(solucion[j]);
		}

		random_shuffle(aux.begin(), aux.end());

		cambio = intercambio(aux[0], max);

		if (cambio){
			pos = buscaPos(max);
			//cout << "pos " << pos << endl;
			mundo.erase(mundo.begin() + pos);

			mundo.push_back(aux[0]);

			for (unsigned i = 0; i < usados.size(); i++){
				mundo.push_back(usados[i]);
			}

			usados.clear();
		}
		else{
			usados.push_back(max);

			pos = buscaPos(max);
			//cout << "pos " << pos << endl;
			mundo.erase(mundo.begin() + pos);
		}

		iteracion = i;
	}
}

double maxim(vector<double> v){
	double max = 0;

	for (unsigned i = 0; i < v.size(); i++){
		if (v[i] > max)
			max = v[i];
	}

	return max;
}

double min(vector<double> v){
	double min = 100000000;

	for (unsigned i = 0; i < v.size(); i++){
		if (v[i] < min)
			min = v[i];
	}

	return min;
}

int rand(int r) {
  /***GENERADOR DE LOS NUMEROS ALEATORIOS***/
      mt19937 gen(r);
      //uniform_real_distribution<double> ran_u(0.0, 1.0);  //rand for probability
      //uniform_int_distribution<int> ran_i(-1,1);          //rand for neighbour
      uniform_int_distribution<int> ran_l(0, n - 1);         //rand for lattice
      // Para generar el aleatorio ran_u(gen);
  /****************/

  return ran_l(gen);
}

void inicioRandom(){
  int i = rand() % 17;
  while (sol.size() < m) {
    int a = rand(i);
    sol.insert(a);
    i += rand() % 29;
  }

  for (auto it = sol.begin(); it != sol.end(); ++it)
    solucion.push_back(*it);
}


//////////////////////////////////
//////////////////////////////////
//////////////////////////////////
//////Hasta aquí la //////////////
///////práctica 1/////////////////
//////////////////////////////////
//////////////////////////////////
//////////////////////////////////


vector <vector <bool> > matrixSolucionesBin, matrixDeSeleccion, matrixDeHijos;
vector <double> scoreTotalBin;
int cromosomas = 50;
double pAGG = 0.7, pAGE = 1, pGen = 0.001, pBL = 0.1;


void inicioRandomBin (){
	matrixSolucionesBin.clear();
vector <bool> solucionBinP;

	for (int i = 0; i < cromosomas; i++){
		for (int j = 0; j < n; j++)
			solucionBinP.push_back(false);
		matrixSolucionesBin.push_back(solucionBinP);
		solucionBinP.clear();
	}

	unsigned contador = 0;

	for (int i = 0; i < cromosomas; i++){
		int r = rand(contador);
		srand(unsigned(200 * i));
		while (contador < m){
			int a = rand(r);
			if (!matrixSolucionesBin[i][a]){
				matrixSolucionesBin[i][a] = true;
				contador++;
			}
				r += rand() % 29000;
		}
		contador = 0;
	}
}

void inicioScore(){
	vector <double> aux;

	for (int i = 0; i < cromosomas; i++)
		scoreTotalBin.push_back(0);
}

double scoreBinUno(int i, bool sol){
	double sum = 0;

	if (sol){
		for (int j = 0; j < n; j++){
			for (int k = 0; k < n; k++){
				sum += matrixSolucionesBin[i][j] * matrixSolucionesBin[i][k] * (matrix[j][k] + matrix[k][j]);
			}
		}

		return sum / 2;
	}
	else{
		for (int j = 0; j < n; j++){
			for (int k = 0; k < n; k++){
				sum += matrixDeHijos[i][j] * matrixDeHijos[i][k] * (matrix[j][k] + matrix[k][j]);
			}
		}

		return sum / 2;
	}
}

void scoreBin(){
	inicioScore();

	for (int i = 0; i < cromosomas; i++){
		scoreTotalBin[i] = scoreBinUno(i, true);
	}
}

void inicioSeleccionG(){
	vector <bool> aux;

	for (int i = 0; i < cromosomas; i++){
		for (int j = 0; j < n; j++)
			aux.push_back(false);
		matrixDeSeleccion.push_back(aux);
	}
}

void inicioSeleccionE(){
	vector <bool> aux;

	for (int i = 0; i < 2; i++){
		for (int j = 0; j < n; j++)
			aux.push_back(false);
		matrixDeSeleccion.push_back(aux);
	}
}

void inicioHijos(){
	vector <bool> aux;

	for (int i = 0; i < cromosomas; i++){
		for (int j = 0; j < n; j++)
			aux.push_back(false);
		matrixDeHijos.push_back(aux);
	}
}

void seleccionBinaria (int pos){
	srand(unsigned(200 * pos));
	int a = rand(pos) % (29 * (pos + 1)), b = rand(pos) % (23 * (pos + 1));

	int r1 = rand(a) % cromosomas;
	int r2;
	do{
		r2 = rand(b) % cromosomas;
		//cout << "r1: " << r1 << "  r2: " << r2 << endl;
		b = b * 37;
	}while(r1 == r2);

	//cout << "r1: " << r1 << "  r2: " << r2 << endl;
	//cout << scoreTotalBin[r1] << "    " << scoreTotalBin[r2] << endl;

	if (scoreTotalBin[r1] > scoreTotalBin[r2]){
		for (int j = 0; j < n; j++){
			matrixDeSeleccion[pos][j] = matrixSolucionesBin[r1][j];
		}
	}
	else{
		for (int j = 0; j < n; j++){
			matrixDeSeleccion[pos][j] = matrixSolucionesBin[r2][j];
		}
	}
}

void cruceUniforme (int posInA, int posInB, int posOut){
	vector <bool> hijo;

	//cout << "Estoy en la salida: " << posOut << " con los padres: " << posInA << " y " << posInB << endl;

	for (int i = 0; i < n; i++)
		hijo.push_back(false);

	for (int i = 0; i < n; i++){
		if (matrixDeSeleccion[posInA][i] == matrixDeSeleccion[posInB][i]){
			hijo[i] = matrixDeSeleccion[posInA][i];
		}
		else{
			if (rand(posOut) % 2 == 1)
				hijo[i] = true;
			//else
				//hijo.push_back(false);
		}
	}

	///////						///////
	/////// Reparador ///////
	///////						///////

	unsigned unos = 0, pos = -1;

	for (int i = 0; i < n; i++)
		if (hijo[i])
			unos++;

	//cout << "Estoy en la salida: " << posOut << " tengo " << unos << " unos." << endl;

	/*if (unos == m)
		cout << "CANELA" << endl;
	else if (unos < m)
		cout << "HAY MENOS" << endl;
	else
	cout << "HAY MAS" << endl;*/

	for (int i = 0; i < n; i++)
		matrixDeHijos[posOut][i] = hijo[i];

	while (unos != m){
		//cout << "while---->   unos " << unos << " salida " << posOut << endl;
		vector <double> dis;
		double sum = 0;

		if (unos < m){
			if (unos == 0){
				hijo[rand(posInA * 101)] = true;
				unos++;
			}
			else{
				double exceso;
				sum = 0;

				for (int i = 0; i < n; i++){
					for (int j = 0; j < n; j++)
						sum += hijo[i] * hijo[j] * matrix[i][j];
				}

				for (int i = 0; i < n; i++){
					if (hijo[i]){
						dis.push_back(sum);
					}
					else{
						for (int j = 0; j < n; j++)
							exceso += hijo[j] * (matrix[i][j] + matrix[j][i]);
						dis.push_back(sum + exceso);
					}
					exceso = 0;
				}

				int maxD = 0;

				for (int i = 0; i < n; i++){
					if (maxD < dis[i]){
						pos = i;
						maxD = dis[i];
					}
				}

				//cout << "posicion: " << pos << endl;

				hijo[pos] = true;
				matrixDeHijos[posOut][pos] = true;

				unos++;
				//cout << "unos (menores) " << unos << " pos " << pos << " salida " << posOut << endl;
			}
		}
		else{
			for (int j = 0; j < n; j++){
				for (int k = 0; k < n; k++){
					sum += hijo[j] * hijo[k] * matrix[j][k];
				}
				dis.push_back(sum);
				sum = 0;
			}

			int maxD = 0;

			for (int i = 0; i < n; i++){
				if (maxD < dis[i]){
					pos = i;
					maxD = dis[i];
				}
			}

			hijo[pos] = false;
			matrixDeHijos[posOut][pos] = false;
			unos--;
			//cout << "unos(mayores) " << unos << " pos " << pos << " salida " << posOut << endl;
		}
	}

	//cout << unos << endl;
}

void cruceBasadoEnPosicion(int posInA, int posInB, int posOut){
	vector <bool> hijo1, hijo2, recomposicion;
	vector <int> sobrantes;
	int tamanyo;

	//cout << "Estoy en la salida: " << posOut << " con los padres: " << posInA << " y " << posInB << endl;

	for (int i = 0; i < n; i++){
		hijo1.push_back(false);
		hijo2.push_back(false);
	}

	for (int i = 0; i < n; i++){
		if (matrixDeSeleccion[posInA][i] == matrixDeSeleccion[posInB][i]){
			hijo1[i] = matrixDeSeleccion[posInA][i];
			hijo2[i] = matrixDeSeleccion[posInA][i];
		}
		else
			sobrantes.push_back(i);
	}

	tamanyo = sobrantes.size();

	for (int i = 0; i < tamanyo; i++)
		recomposicion.push_back(matrixDeSeleccion[posInA][sobrantes[i]]);

	random_shuffle(recomposicion.begin(), recomposicion.end());

	for (int i = 0; i < tamanyo; i++)
		hijo1[sobrantes[i]] = recomposicion[i];

	random_shuffle(recomposicion.begin(), recomposicion.end());

	for (int i = 0; i < tamanyo; i++)
		hijo2[sobrantes[i]] = recomposicion[i];

	for (int i = 0; i < n; i++){
		matrixDeHijos[posOut][i] = hijo1[i];
		matrixDeHijos[posOut + 1][i] = hijo2[i];
	}
}

void AlgoritmoGeneticoGCU(){
	inicioRandomBin();
	inicioSeleccionG();

	for (int x = 0; x < 50000 / cromosomas; x++){
		unsigned t0, t1;

		//t0 = clock();
		scoreBin();
		//t1 = clock();

		//cout << "Tiempo Score " << (double(t1-t0)/CLOCKS_PER_SEC) << endl;

		/*for (int i = 0; i < cromosomas; i++)
			cout << scoreTotalBin[i] << " ";
		cout << endl << endl;*/


		//if (x % 5 == 0)
			//cout << x << endl;

		/*for (int i = 0; i < cromosomas; i++)
			cout << scoreTotalBin[i] << " ";
		cout << endl;*/
		inicioHijos();

		for (int i = 0; i < cromosomas; i++)
			seleccionBinaria(i);

		/*for (int i = 0; i < cromosomas; i++){
			for (int j = 0; j < n; j++)
				cout << matrixDeSeleccion[i][j] << " ";
			cout << endl;
		}*/

		int numero_de_cruces = cromosomas * pAGG;
		int contador = 0;

		t0 = clock();
		for (int i = 0; i < cromosomas && contador < numero_de_cruces; i = i + 2){
			cruceUniforme(i, i + 1, i);
			contador++;
			if (contador < numero_de_cruces){
				cruceUniforme(i, i + 1, i + 1);
				contador++;
			}
		}
		t1 = clock();
		cout << "Tiempo Cruce " << (double(t1-t0)/CLOCKS_PER_SEC) << endl;

		/*cout << endl << endl;
		for (int i = 0; i < cromosomas; i++){
			for (int j = 0; j < n; j++)
				cout << matrixDeHijos[i][j] << " ";
			cout << endl;
		}*/

		for (int i = 0; i < cromosomas; i++){
			for (int j = 0; j < n; j++)
				if (i < numero_de_cruces)
					matrixSolucionesBin[i][j] = matrixDeHijos[i][j];
				else
					matrixSolucionesBin[i][j] = matrixDeSeleccion[i][j];
		}

		int numero_de_mutaciones = pGen * cromosomas * n;

		for (int i = 0; i < numero_de_mutaciones; i++){
			int a = rand(i) % (317 * (i + 1)), b = rand(i) % (483 * (i + 1));

			int r = (rand(a) + 1) * (rand(a) % cromosomas + 1) - 1;
			int s;

			do{
				s = rand(b);
				b = b * 73;
				//cout << "r " << r << " s " << s << endl;
				//cout << "r % cromosomas " << r % cromosomas << " r / cromosomas " << r / cromosomas << " s / cromosomas " << s / cromosomas << endl;
			}while (matrixSolucionesBin[r % cromosomas][r / cromosomas] == matrixSolucionesBin[r % cromosomas][s]);
			bool swap;
			swap = matrixSolucionesBin[r % cromosomas][r / cromosomas];
			matrixSolucionesBin[r % cromosomas][r / cromosomas] = matrixSolucionesBin[r % cromosomas][s];
			matrixSolucionesBin[r % cromosomas][s] = swap;
		}

		scoreTotalBin.clear();
	}
}

void AlgoritmoGeneticoGCBP(){
	inicioRandomBin();
	inicioSeleccionG();

	for (int x = 0; x < 50000 / cromosomas; x++){
		//unsigned t0, t1;

		//t0 = clock();
		scoreBin();
		//t1 = clock();

		//cout << "Tiempo Score " << (double(t1-t0)/CLOCKS_PER_SEC) << endl;

		//if (x % 5 == 0)
			//cout << x << endl;

		inicioHijos();

		for (int i = 0; i < cromosomas; i++)
			seleccionBinaria(i);

		/*cout << endl << "Selección: " << endl;

		cout << endl << endl;
		for (int i = 0; i < cromosomas; i++){
			for (int j = 0; j < n; j++)
				cout << matrixDeSeleccion[i][j] << " ";
			cout << endl;
		}*/

		int numero_de_cruces = cromosomas * pAGG;
		int contador = 0;

		for (int i = 0; i < cromosomas / 2 && contador < numero_de_cruces; i++){
			cruceBasadoEnPosicion(i, i + 1, 2 * i);
			contador++;
		}

		for (int i = 0; i < cromosomas; i++){
			for (int j = 0; j < n; j++)
				if (i < numero_de_cruces)
					matrixSolucionesBin[i][j] = matrixDeHijos[i][j];
				else
					matrixSolucionesBin[i][j] = matrixDeSeleccion[i][j];
		}

		/*cout << endl << "Nueva Generación: " << endl;

		cout << endl << endl;
		for (int i = 0; i < cromosomas; i++){
			for (int j = 0; j < n; j++)
				cout << matrixSolucionesBin[i][j] << " ";
			cout << endl;
		}*/

		int numero_de_mutaciones = pGen * cromosomas * n;

		for (int i = 0; i < numero_de_mutaciones; i++){
			int a = rand(i) % (317 * (i + 1)), b = rand(i) % (483 * (i + 1));

			int r = (rand(a) + 1) * (rand(a) % cromosomas + 1) - 1;
			int s;

			do{
				s = rand(b);
				b = b * 73;
				//cout << "r " << r << " s " << s << endl;
				//cout << "r % cromosomas " << r % cromosomas << " r / cromosomas " << r / cromosomas << " s / cromosomas " << s / cromosomas << endl;
			}while (matrixSolucionesBin[r % cromosomas][r / cromosomas] == matrixSolucionesBin[r % cromosomas][s]);
			bool swap;
			swap = matrixSolucionesBin[r % cromosomas][r / cromosomas];
			matrixSolucionesBin[r % cromosomas][r / cromosomas] = matrixSolucionesBin[r % cromosomas][s];
			matrixSolucionesBin[r % cromosomas][s] = swap;
		}

		/*cout << endl << "Tras Mutar: " << endl;

		cout << endl << endl;
		for (int i = 0; i < cromosomas; i++){
			for (int j = 0; j < n; j++)
				cout << matrixSolucionesBin[i][j] << " ";
			cout << endl;
		}*/

		scoreTotalBin.clear();
	}
}

void AlgoritmoGeneticoECU (){
	inicioRandomBin();
	inicioSeleccionE();

	for (int x = 0; x < (50000 - cromosomas) / 2; x++){
		//unsigned t0, t1;

		if (x == 0)
			scoreBin();

		//if (x % 5 == 0)
			//cout << x << endl;

		inicioHijos();

		for (int i = 0; i < 2; i++)
			seleccionBinaria(i);

		/*for (int i = 0; i < cromosomas; i++){
			for (int j = 0; j < n; j++)
				cout << matrixDeSeleccion[i][j] << " ";
			cout << endl;
		}*/

		//int numero_de_cruces = cromosomas * pAGE;

		//t0 = clock();
		cruceUniforme(0, 1, 0);
		cruceUniforme(0, 1, 1);
		//t1 = clock();

		//cout << "Tiempo Cruce " << (double(t1-t0)/CLOCKS_PER_SEC) << endl;

		/*cout << endl << endl;
		for (int i = 0; i < cromosomas; i++){
			for (int j = 0; j < n; j++)
				cout << matrixDeHijos[i][j] << " ";
			cout << endl;
		}*/

		int numero_de_mutaciones = pGen * 2 * n;
		//cout << "mutaciones: " << numero_de_mutaciones << endl;

		for (int i = 0; i < numero_de_mutaciones; i++){
			int a = rand(i); // % (1005 * (i + 1)), b = rand(i) % (483 * (i + 1));

			int b = rand(i) % (2 * m), c = 0, d = -1;

			for (int j = 0; j < 2 * n && d == -1; j++){
				if (matrixDeHijos[j / n][j % n])
					c++;
				if (c == b)
					d = j;
			}
			int s;

			do{
				//srand(unsigned(b * 200));
				s = rand(a);
				a*= 97;
				//cout << "d " << d << "    " << matrixDeHijos[d / n][d % n] << " s " << s << " i " << i << endl;
			}while (matrixDeHijos[d / n][d % n] == matrixDeHijos[d / n][s]);
			bool swap;
			swap = matrixDeHijos[d / n][d % n];
			matrixDeHijos[d / n][d % n] = matrixDeHijos[d / n][s];
			matrixDeHijos[d / n][s] = swap;
		}

		double score0, score1, maxG, maxM;
		score0 = scoreBinUno(0, false);
		score1 = scoreBinUno(1, false);

		if (score0 > score1){
			maxG = score0;
			maxM = score1;
		}
		else{
			maxG = score1;
			maxM = score0;
		}

		int minG = 0;
		double minDG = 100000000;

		for (int i = 0; i < cromosomas; i++){
			if (minDG > scoreTotalBin[i]){
				minG = i;
				minDG = scoreTotalBin[i];
			}
		}

		if (maxG < minDG){
			//// Nothing
		}
		else if (maxG > minDG && maxM < minDG){
			scoreTotalBin[minG] = maxG;
			for (int i = 0; i < n; i++){
				if (maxG == score0)
					matrixSolucionesBin[minG][i] = matrixDeHijos[0][i];
				else
					matrixSolucionesBin[minG][i] = matrixDeHijos[1][i];
			}
		}
		else{
			scoreTotalBin[minG] = 0;

			int minM = 0;
			double minDM = 100000000;

			for (int i = 0; i < cromosomas; i++){
				if (minDM > scoreTotalBin[i] && scoreTotalBin[i] != 0){
					minM = i;
					minDM = scoreTotalBin[i];
				}
			}

			scoreTotalBin[minG] = maxG;
			scoreTotalBin[minM] = maxM;
			for (int i = 0; i < n; i++){
				if (maxG == score0){
					matrixSolucionesBin[minG][i] = matrixDeHijos[0][i];
					matrixSolucionesBin[minM][i] = matrixDeHijos[1][i];
				}
				else{
					matrixSolucionesBin[minG][i] = matrixDeHijos[1][i];
					matrixSolucionesBin[minM][i] = matrixDeHijos[0][i];
				}
			}
		}
	}
}

void AlgoritmoGeneticoECBP (){
	inicioRandomBin();
	inicioSeleccionE();

	for (int x = 0; x < (50000 - cromosomas) / 2; x++){
		//unsigned t0, t1;

		if (x == 0)
			scoreBin();

		//if (x % 5 == 0)
			//cout << x << endl;

		inicioHijos();

		for (int i = 0; i < 2; i++)
			seleccionBinaria(i);

		/*for (int i = 0; i < cromosomas; i++){
			for (int j = 0; j < n; j++)
				cout << matrixDeSeleccion[i][j] << " ";
			cout << endl;
		}*/

		//int numero_de_cruces = cromosomas * pAGE;

		cruceBasadoEnPosicion(0, 1, 0);

		/*cout << endl << endl;
		for (int i = 0; i < cromosomas; i++){
			for (int j = 0; j < n; j++)
				cout << matrixDeHijos[i][j] << " ";
			cout << endl;
		}*/

		int numero_de_mutaciones = pGen * 2 * n;
		//cout << "mutaciones: " << numero_de_mutaciones << endl;

		for (int i = 0; i < numero_de_mutaciones; i++){
			int a = rand(i); // % (1005 * (i + 1)), b = rand(i) % (483 * (i + 1));

			int b = rand(i) % (2 * m), c = 0, d = -1;

			for (int j = 0; j < 2 * n && d == -1; j++){
				if (matrixDeHijos[j / n][j % n])
					c++;
				if (c == b)
					d = j;
			}
			int s;

			do{
				//srand(unsigned(b * 200));
				s = rand(a);
				a*= 97;
				//cout << "d " << d << "    " << matrixDeHijos[d / n][d % n] << " s " << s << " i " << i << endl;
			}while (matrixDeHijos[d / n][d % n] == matrixDeHijos[d / n][s]);
			bool swap;
			swap = matrixDeHijos[d / n][d % n];
			matrixDeHijos[d / n][d % n] = matrixDeHijos[d / n][s];
			matrixDeHijos[d / n][s] = swap;
		}

		double score0, score1, maxG, maxM;
		score0 = scoreBinUno(0, false);
		score1 = scoreBinUno(1, false);

		if (score0 > score1){
			maxG = score0;
			maxM = score1;
		}
		else{
			maxG = score1;
			maxM = score0;
		}

		int minG = 0;
		double minDG = 100000000;

		for (int i = 0; i < cromosomas; i++){
			if (minDG > scoreTotalBin[i]){
				minG = i;
				minDG = scoreTotalBin[i];
			}
		}

		if (maxG < minDG){
			//// Nothing
		}
		else if (maxG > minDG && maxM < minDG){
			scoreTotalBin[minG] = maxG;
			for (int i = 0; i < n; i++){
				if (maxG == score0)
					matrixSolucionesBin[minG][i] = matrixDeHijos[0][i];
				else
					matrixSolucionesBin[minG][i] = matrixDeHijos[1][i];
			}
		}
		else{
			scoreTotalBin[minG] = 0;

			int minM = 0;
			double minDM = 100000000;

			for (int i = 0; i < cromosomas; i++){
				if (minDM > scoreTotalBin[i] && scoreTotalBin[i] != 0){
					minM = i;
					minDM = scoreTotalBin[i];
				}
			}

			scoreTotalBin[minG] = maxG;
			scoreTotalBin[minM] = maxM;
			for (int i = 0; i < n; i++){
				if (maxG == score0){
					matrixSolucionesBin[minG][i] = matrixDeHijos[0][i];
					matrixSolucionesBin[minM][i] = matrixDeHijos[1][i];
				}
				else{
					matrixSolucionesBin[minG][i] = matrixDeHijos[1][i];
					matrixSolucionesBin[minM][i] = matrixDeHijos[0][i];
				}
			}
		}
	}
}

void binToInt(int pos){
	vector <int> aux;

	for (int i = 0; i < n; i++){
		if (matrixSolucionesBin[pos][i])
			aux.push_back(i);
	}

	solucion.clear();
	mundo.clear();

	for (unsigned i = 0; i < aux.size(); i++)
		solucion.push_back(aux[i]);

	crearMundo();
}

void intToBin(int pos){
	vector <bool> aux;

	for (int i = 0; i < n; i++)
		aux.push_back(false);

	for (unsigned i = 0; i < m; i++){
			aux[solucion[i]] = true;
	}

	for (int i = 0; i < n; i++)
		matrixSolucionesBin[pos][i] = aux[i];
}

void AlgoritmoMemetico(unsigned generaciones, double proporcion, bool mejor){
	inicioRandomBin();
	inicioSeleccionG();

	for (int x = 0; x < 50000 / (410 * proporcion); x++){
		scoreBin();

		//if (x % 5 == 0)
			//cout << x << endl;

		inicioHijos();

		for (int i = 0; i < cromosomas; i++)
			seleccionBinaria(i);

		int numero_de_cruces = cromosomas * pAGG;
		int contador = 0;

		for (int i = 0; i < cromosomas && contador < numero_de_cruces; i = i + 2){
			cruceUniforme(i, i + 1, i);
			contador++;
			if (contador < numero_de_cruces){
				cruceUniforme(i, i + 1, i + 1);
				contador++;
			}
		}

		for (int i = 0; i < cromosomas; i++){
			for (int j = 0; j < n; j++)
				if (i < numero_de_cruces)
					matrixSolucionesBin[i][j] = matrixDeHijos[i][j];
				else
					matrixSolucionesBin[i][j] = matrixDeSeleccion[i][j];
		}

		int numero_de_mutaciones = pGen * 2 * n;
		//cout << "mutaciones: " << numero_de_mutaciones << endl;

		for (int i = 0; i < numero_de_mutaciones; i++){
			int a = rand(i) % (317 * (i + 1)), b = rand(i) % (483 * (i + 1));

			int r = (rand(a) + 1) * 2 - 1;
			int s;

			do{
				s = rand(b);
				b = b * 73;
				//cout << "r " << r << " s " << s << endl;
			}while (matrixDeHijos[r % cromosomas][r / cromosomas] == matrixDeHijos[r % cromosomas][s]);
			bool swap;
			swap = matrixDeHijos[r % cromosomas][r / cromosomas];
			matrixDeHijos[r % cromosomas][r / cromosomas] = matrixDeHijos[r % cromosomas][s];
			matrixDeHijos[r % cromosomas][s] = swap;
		}



		if (x % cromosomas == cromosomas - 1){
			unsigned numero_de_solucionesBL = cromosomas * proporcion;
			int max = 0;
			vector <int> para_hacerBL;

			for (unsigned i = 0; i < numero_de_solucionesBL; i++){
				if (mejor){
					vector <double> aux;
					double maxD = 0;

					for (int j = 0; j < cromosomas; j++)
						aux.push_back(scoreTotalBin[j]);

					for (int j = 0; j < cromosomas; j++){
						if (maxD < aux[j]){
							max = j;
							maxD = aux[j];
						}
					}
					aux[i] = 0;
					para_hacerBL.push_back(max);
				}
			}

			if(!mejor){
				sol.clear();
				int t = rand(x) % 107;
			  while (sol.size() < numero_de_solucionesBL) {
			    int a = rand(t) % cromosomas;
			    sol.insert(a);
			    t += 79;
			  }

			  for (auto it = sol.begin(); it != sol.end(); ++it)
			    para_hacerBL.push_back(*it);
			}

			for (unsigned i = 0; i < numero_de_solucionesBL; i++){
				binToInt(para_hacerBL[i]);
				scoreTotal = scoreTotalBin[i];
				busquedaLocal(400);
				intToBin(para_hacerBL[i]);
			}
		}
	}
}

int main(int argc, char const *argv[]){

	vector<double> greedy, local, desviacionG, desviacionL, mediaG, mediaL,
		mejorG, mejorL, peorG, peorL, tiempoG, tiempoL, mediaTiempoG, mediaTiempoL;
	vector<double> mejorValor, tiempoGen, scoreGen, vecesGen, desvGen;
	double maxD = 0;

	srand(unsigned(200));

	unsigned t0, t1;

	if (argc >= 2)
		cromosomas = atoi(argv[1]);

	mejorValor.push_back(19485.1875);
	mejorValor.push_back(19701.53711);
	mejorValor.push_back(19547.20703);
	mejorValor.push_back(19596.46875);
	mejorValor.push_back(19602.625);
	mejorValor.push_back(19421.55078);
	mejorValor.push_back(19534.30664);
	mejorValor.push_back(19487.32031);
	mejorValor.push_back(19221.63477);
	mejorValor.push_back(19703.35156);
	mejorValor.push_back(20743);
	mejorValor.push_back(35881);
	mejorValor.push_back(4658);
	mejorValor.push_back(16956);
	mejorValor.push_back(36317);
	mejorValor.push_back(62487);
	mejorValor.push_back(7141);
	mejorValor.push_back(26258);
	mejorValor.push_back(56572);
	mejorValor.push_back(97344);
	mejorValor.push_back(114259);
	mejorValor.push_back(114327);
	mejorValor.push_back(114123);
	mejorValor.push_back(114040);
	mejorValor.push_back(114064);
	mejorValor.push_back(114204);
	mejorValor.push_back(114338);
	mejorValor.push_back(114158);
	mejorValor.push_back(114132);
	mejorValor.push_back(114197);


	/*for (int i = 21; i < 22; i++){
		string a;
		if (i < 11)
			a = "Tablas/GKD-c_" + to_string(i) + "_n500_m50.txt";

		else if (i == 11)
			a = "Tablas/SOM-b_11_n300_m90.txt";

		else if (i == 12)
			a = "Tablas/SOM-b_12_n300_m120.txt";

		else if (i == 13)
			a = "Tablas/SOM-b_13_n400_m40.txt";

		else if (i == 14)
			a = "Tablas/SOM-b_14_n400_m80.txt";

		else if (i == 15)
			a = "Tablas/SOM-b_15_n400_m120.txt";

		else if (i == 16)
			a = "Tablas/SOM-b_16_n400_m160.txt";

		else if (i == 17)
			a = "Tablas/SOM-b_17_n500_m50.txt";

		else if (i == 18)
			a = "Tablas/SOM-b_18_n500_m100.txt";

		else if (i == 19)
			a = "Tablas/SOM-b_19_n500_m150.txt";

		else if (i == 20)
			a = "Tablas/SOM-b_20_n500_m200.txt";

		else if (i > 20)
			a = "Tablas/MDG-a_" + to_string(i) + "_n2000_m200.txt";

		for (int j = 0; j < 1; j++){
			srand(unsigned(j * 100));

			input(a);

			t0 = clock();
			algoritmoGreedy();
			t1 = clock();

			tiempoG.push_back((double(t1-t0)/CLOCKS_PER_SEC));

			score();
			greedy.push_back(scoreTotal);

      destruirLight();
      inicioRandom();
      crearMundo();
      score();

			t0 = clock();
			busquedaLocal(50000);
			t1 = clock();

			tiempoL.push_back((double(t1-t0)/CLOCKS_PER_SEC));

			local.push_back(scoreTotal);

			t0 = clock();
			AlgoritmoGeneticoGCU();
			scoreBin();
			t1 = clock();

			cout << double(t1-t0)/CLOCKS_PER_SEC << endl;

			tiempoGen.push_back((double(t1-t0)/CLOCKS_PER_SEC));

			for (int i = 0; i < cromosomas; i++){
				if (maxD < scoreTotalBin[i]){
					maxD = scoreTotalBin[i];
				}
			}

			scoreGen.push_back(maxD);

			cout << "Genético 1" << endl;

			unsigned cont = 0;
			for (int i = 0; i < cromosomas; i++){
				if (scoreTotalBin[i] == maxD)
					cont++;
			}

			vecesGen.push_back(cont);
			maxD = 0;
			scoreTotalBin.clear();

			t0 = clock();
			AlgoritmoGeneticoGCBP();
			scoreBin();
			t1 = clock();

			tiempoGen.push_back((double(t1-t0)/CLOCKS_PER_SEC));

			for (int i = 0; i < cromosomas; i++){
				if (maxD < scoreTotalBin[i]){
					maxD = scoreTotalBin[i];
				}
			}

			scoreGen.push_back(maxD);

			cout << "Genético 2" << endl;

			cont = 0;
			for (int i = 0; i < cromosomas; i++){
				if (scoreTotalBin[i] == maxD)
					cont++;
			}

			vecesGen.push_back(cont);
			maxD = 0;
			scoreTotalBin.clear();

			t0 = clock();
			AlgoritmoGeneticoECU();
			scoreBin();
			t1 = clock();

			tiempoGen.push_back((double(t1-t0)/CLOCKS_PER_SEC));

			for (int i = 0; i < cromosomas; i++){
				if (maxD < scoreTotalBin[i]){
					maxD = scoreTotalBin[i];
				}
			}

			scoreGen.push_back(maxD);

			cout << "Genético 3" << endl;

			cont = 0;
			for (int i = 0; i < cromosomas; i++){
				if (scoreTotalBin[i] == maxD)
					cont++;
			}

			vecesGen.push_back(cont);
			maxD = 0;
			scoreTotalBin.clear();

			t0 = clock();
			AlgoritmoGeneticoECBP();
			scoreBin();
			t1 = clock();

			tiempoGen.push_back((double(t1-t0)/CLOCKS_PER_SEC));

			for (int i = 0; i < cromosomas; i++){
				if (maxD < scoreTotalBin[i]){
					maxD = scoreTotalBin[i];
				}
			}

			scoreGen.push_back(maxD);

			cout << "Genético 4" << endl;

			cont = 0;
			for (int i = 0; i < cromosomas; i++){
				if (scoreTotalBin[i] == maxD)
					cont++;
			}

			vecesGen.push_back(cont);
			maxD = 0;
			scoreTotalBin.clear();

			if (argc >= 3)
				cromosomas = atoi(argv[2]);

			t0 = clock();
			AlgoritmoMemetico(10, 1, false);
			scoreBin();
			t1 = clock();

			tiempoGen.push_back((double(t1-t0)/CLOCKS_PER_SEC));

			for (int i = 0; i < cromosomas; i++){
				if (maxD < scoreTotalBin[i]){
					maxD = scoreTotalBin[i];
				}
			}

			scoreGen.push_back(maxD);

			cout << "Memético 1" << endl;

			cont = 0;
			for (int i = 0; i < cromosomas; i++){
				if (scoreTotalBin[i] == maxD)
					cont++;
			}

			vecesGen.push_back(cont);
			maxD = 0;
			scoreTotalBin.clear();

			t0 = clock();
			AlgoritmoMemetico(10, 0.1, false);
			scoreBin();
			t1 = clock();

			tiempoGen.push_back((double(t1-t0)/CLOCKS_PER_SEC));

			for (int i = 0; i < cromosomas; i++){
				if (maxD < scoreTotalBin[i]){
					maxD = scoreTotalBin[i];
				}
			}

			scoreGen.push_back(maxD);

			cout << "Memético 2" << endl;

			cont = 0;
			for (int i = 0; i < cromosomas; i++){
				if (scoreTotalBin[i] == maxD)
					cont++;
			}

			vecesGen.push_back(cont);
			maxD = 0;
			scoreTotalBin.clear();

			t0 = clock();
			AlgoritmoMemetico(10, 0.1, true);
			scoreBin();
			t1 = clock();

			tiempoGen.push_back((double(t1-t0)/CLOCKS_PER_SEC));

			for (int i = 0; i < cromosomas; i++){
				if (maxD < scoreTotalBin[i]){
					maxD = scoreTotalBin[i];
				}
			}

			scoreGen.push_back(maxD);

			cout << "Memético 3" << endl;

			cont = 0;
			for (int i = 0; i < cromosomas; i++){
				if (scoreTotalBin[i] == maxD)
					cont++;
			}

			vecesGen.push_back(cont);
			maxD = 0;
			scoreTotalBin.clear();

			destruir(); //////////////////////////////// AL FINAL DEL TODO SIEMPRE
		}

		mejorG.push_back(maxim(greedy));
		mejorL.push_back(maxim(local));
		peorG.push_back(min(greedy));
		peorL.push_back(min(local));

		double desvG = 0, desvL = 0, mG = 0, mL = 0, tG = 0, tL = 0;

		for (int k = 0; k < 1; k++){
			desvG += (mejorValor[i - 1] - greedy[k]) / mejorValor[i - 1];
			desvL += (mejorValor[i - 1] - local[k]) / mejorValor[i - 1];
			mG += greedy[k];
			mL += local[k];
			tG += tiempoG[k];
			tL += tiempoL[k];
		}

		desviacionG.push_back(desvG * 100);
		desviacionL.push_back(desvL * 100);
		mediaG.push_back(mG);
		mediaL.push_back(mL);
		mediaTiempoG.push_back(tG);
		mediaTiempoL.push_back(tL);

		greedy.clear();
		local.clear();
		tiempoG.clear();
		tiempoL.clear();

		for (int j = 0; j < 7; j++){
			desvGen.push_back((mejorValor[i - 1] - scoreGen[7 * (i - 1) + j]) / mejorValor[i - 1]);
		}
	}

	for (int i = 0; i < 1; i++){
		cout << endl << "Caso " << i + 1 << ":" << endl;
		cout << "		-Greedy: " << " 	Desviación Greedy: " << desviacionG[i] << "  Tiempo: " << mediaTiempoG[i] << endl;
		cout << "		-Búsqueda Local: " << " 	Desviación Local: " << desviacionL[i] << " Tiempo: " << mediaTiempoL[i] << endl;
		for (int j = 0; j < 7; j++){
			switch (j) {
				case 0:
					cout << "		-AGG Cruce Uniforme: ";
				break;

				case 1:
					cout << "		-AGG Cruce Basado en Posición: ";
				break;

				case 2:
					cout << "		-AGE Cruce Uniforme: ";
				break;

				case 3:
					cout << "		-AGE Cruce Basado en Posición: ";
				break;

				case 4:
					cout << "		-AM-(10, 1): ";
				break;

				case 5:
					cout << "		-AM-(10, 0.1): ";
				break;

				case 6:
					cout << "		-AM-(10, 0.1mej): ";
				break;
			}
			cout <<  " 	Desviación: " << desvGen[7 * i + j]  << " Tiempo: " << tiempoGen[7 * i + j] << endl;
			cout << "		El mejor se ha repetido: " << vecesGen[7 * i + j] << " veces" << endl;
		}
	}*/

	input("Tablas/MDG-a_21_n2000_m200.txt");

	/*t0 = clock();
	algoritmoGreedy();
	score();
	t1 = clock();

	cout << "Tiempo " << (double(t1-t0)/CLOCKS_PER_SEC) << endl;

	cout << "Desviacion " << (mejorValor[20] - scoreTotal) / mejorValor[20] << endl;

  destruirLight();

	t0 = clock();
	inicioRandom();
  score();
	crearMundo();
  busquedaLocal(50000);
	t1 = clock();

	cout << "Tiempo " << (double(t1-t0)/CLOCKS_PER_SEC) << endl;

	cout << "Desviacion " << (mejorValor[20] - scoreTotal) / mejorValor[20] << endl;*/

	t0 = clock();
	//AlgoritmoGeneticoGCU();
	//AlgoritmoGeneticoGCBP();
	//AlgoritmoGeneticoECU();
	//AlgoritmoGeneticoECBP();
	//AlgoritmoMemetico(10, 1, false);
	//AlgoritmoMemetico(10, 0.1, false);
	AlgoritmoMemetico(10, 0.1, true);
	t1 = clock();

	cout << "Tiempo " << (double(t1-t0)/CLOCKS_PER_SEC) << endl;

	for (int i = 0; i < cromosomas; i++){
		if (maxD < scoreTotalBin[i]){
			maxD = scoreTotalBin[i];
		}
	}

	unsigned cont = 0;
	for (int i = 0; i < cromosomas; i++){
		if (scoreTotalBin[i] == maxD)
			cont++;
	}

	cout << "Desviacion " << (mejorValor[20] - maxD) / mejorValor[20] << endl;

	cout << "Ha salido " << cont << " veces" << endl;

	/*inicioRandomBin();
	t0 = clock();
	scoreBin();
	t1 = clock();

	cout << (double(t1-t0)/CLOCKS_PER_SEC) << endl;*/

	/*cout << endl;
	for (int i = 0; i < cromosomas; i++)
		cout << scoreTotalBin[i] << " ";
	cout << endl;*/

 return 0;
}
