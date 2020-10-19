#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <fstream>
#include <vector>
#include <algorithm>
using namespace std;

vector< vector<double> > matrix;
vector<int> solucion, mundo, usados, aux;
vector<double> distancias_a_la_solucion, distancia_anyadida;
int n, m, iteracion;
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

	for (int i = 0; i < solucion.size(); i++){
		for (int j = 0; j < solucion.size(); j++){
			//if (i != j){
				//if (solucion.at(i) < solucion.at(j)){
					sum += matrix[solucion[i]][solucion[j]];
					//sum += matrix.at(solucion.at(i)).at(solucion.at(j));
				/*}
				else{
					sum += matrix.at(solucion.at(j)).at(solucion.at(i));
				}
			}*/
		}
	}

	scoreTotal = sum;
}

double distanciaAcumulada(int a){
	double suma = 0;

	for (int i = 0; i < n; i++){
		if (i < a){
			suma += matrix[i][a];
			//suma += matrix.at(i).at(a);
			//cout << matrix.at(i).at(a) << endl;
		}
		else{
			suma += matrix[a][i];
			//suma += matrix.at(a).at(i);
			//cout << matrix.at(a).at(i) << endl;
		}
	}

	return suma;
}

bool remakeScore(int posI, int j){
	double restar, sumar;
	bool cambio = false;

	restar = distancia_anyadida[posI];
	//restar = distancia_anyadida.at(posI);
	sumar = distancia_anyadida[j];
	//sumar = distancia_anyadida.at(j);

	//cout << "sumar: " << sumar << " restar: " << restar << endl << endl;

	//for (int i = 0; i < distancia_anyadida.size(); i++)
		//cout << distancia_anyadida.at(i) << " ";
	//cout << endl << endl;

	if (sumar > restar){
		cambio = true;
		scoreTotal += sumar - restar;
	}

	return cambio;
}

bool intercambio (int i, int j){
	bool cambio;
	int pos = -1;

	for (int k = 0; k < m; k++){
		if (/*solucion.at(k)*/ solucion[k] == i )
			pos = k;
	}

	if (pos != -1){
		cambio = remakeScore(i, j);

		if (cambio){
			//solucion.at(pos) = j;
			solucion[pos] = j;
		}
	}

	else{
		cout << "Error: pos = -1" << endl;
	}

	return cambio;
}

int buscaPos(int a){
	int pos;

	for (int i = 0; i < mundo.size(); i++){
		if (/*mundo.at(i)*/ mundo[i] == a)
			pos = i;
			//cout << mundo.at(i) << " elemento " << a << endl;
	}

	//cout << "la posicion es :" << pos << endl;

	return pos;
}

int buscaMaximo(){
	int max = 0;
	int maxD = /*distancias_a_la_solucion.at(0)*/ distancias_a_la_solucion[0];

	for (int i = 1; i < distancias_a_la_solucion.size(); i++){
		if (maxD < /*distancias_a_la_solucion.at(i)*/ distancias_a_la_solucion[i]){
			max = i;
			maxD = /*distancias_a_la_solucion.at(i)*/ distancias_a_la_solucion[i];
		}
	}

	return max;
}

int buscaMaximoAnyadido(){
	int max = 0;
	int maxD = 0;

	for (int i = 0; i < mundo.size(); i++){
		//cout << endl << "escoger maximo " << mundo.at(i) << "distancia" << distancia_anyadida.at(mundo.at(i)) << endl;
		if (maxD < /*distancia_anyadida.at(mundo.at(i))*/ distancia_anyadida[mundo[i]]){
			max = /*mundo.at(i)*/ mundo[i];
			maxD = /*distancia_anyadida.at(mundo.at(i))*/ distancia_anyadida[mundo[i]];
		}
	}

	return max;
}

void actualizaDistancias(int a){
	for (int i = 0; i < distancias_a_la_solucion.size(); i++){
		distancias_a_la_solucion[i] -= matrix[i][a] + matrix[a][i];
		//distancias_a_la_solucion.at(i) -= matrix.at(i).at(a) + matrix.at(a).at(i);
	}
	distancias_a_la_solucion[a] = 0;
	//distancias_a_la_solucion.at(a) = 0;
}

void actualizaDistanciasAnyadidas(int a){
	for (int i = 0; i < distancia_anyadida.size(); i++){
		distancia_anyadida[i] -= matrix[i][a] + matrix[a][i];
		//distancia_anyadida.at(i) -= matrix.at(i).at(a) + matrix.at(a).at(i);
		//cout << "d " << distancia_anyadida.at(i) << " ";
	}
	//cout << endl;
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

void algoritmoGreedy(){
	int max = 0, pos;
	double maxD, posibleMaxD;

	maxD = distanciaAcumulada(0);
	distancias_a_la_solucion.push_back(maxD);

	for (int i = 0; i < n; i++)
		mundo.push_back(i);

	for (int i = 1; i < mundo.size(); i++){
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

void busquedaLocal(){
	int max, pos;
	bool cambio;
	double sum;

	for(int i = 0; i < 50000 && !mundo.empty(); i++){
		distancia_anyadida.clear();

		for (int k = 0; k < n; k++){
			sum = 0;
			for (int j = 0; j < m; j++){
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

		for (int j = 0; j < usados.size(); j++){
			distancia_anyadida[usados[j]] = 0;
			//distancia_anyadida.at(usados.at(j)) = 0;
		}

		/*for (int j = 0; j < mundo.size(); j++)
			cout << mundo.at(j) << " ";
		cout << endl;*/

		max = buscaMaximoAnyadido();

		aux.clear();

		for (int j = 0; j < m; j++){
			aux.push_back(solucion[j]);
			//aux.push_back(solucion.at(j));
		}

		random_shuffle(aux.begin(), aux.end());

		//cout << "cambio el nuevo: " << max << " por el de la solucion: " << aux.at(0) << endl;

		actualizaDistanciasAnyadidas(/*aux.at(0)*/ aux[0]);

		//cout << "max: " << max << endl;

		/*for (int j = 0; j < distancia_anyadida.size(); j++){
			cout << distancia_anyadida.at(j) << " ";
		}
		cout << endl;*/

		cambio = intercambio(/*aux.at(0)*/ aux[0], max);

		if (cambio){
			//actualizaDistancias(max);

			pos = buscaPos(max);
			//cout << "pos " << pos << endl;
			mundo.erase(mundo.begin() + pos);

			mundo.push_back(/*aux.at(0)*/ aux[0]);

			usados.clear();
		}
		else{
			usados.push_back(max);

			pos = buscaPos(max);
			//cout << "pos " << pos << endl;
			mundo.erase(mundo.begin() + pos);
		}

		iteracion = i;

		//if (iteracion % 500 == 0)
			//cout << iteracion << endl;
	}
}

double max(vector<double> v){
	double max = 0;

	for (int i = 0; i < v.size(); i++){
		if (/*v.at(i)*/ v[i] > max)
			max = /*v.at(i)*/ v[i];
	}

	return max;
}

double min(vector<double> v){
	double min = 100000000;

	for (int i = 0; i < v.size(); i++){
		if (/*v.at(i)*/ v[i] < min)
			min = /*v.at(i)*/ v[i];
	}

	return min;
}

int main(){

	vector<double> greedy, local, desviacionG, desviacionL, mediaG, mediaL,
		mejorG, mejorL, peorG, peorL, tiempoG, tiempoL, mediaTiempoG, mediaTiempoL;
	vector<double> mejorValor;

	unsigned t0, t1;

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


	for (int i = 1; i < 31; i++){
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

			t0 = clock();
			busquedaLocal();
			t1 = clock();

			tiempoL.push_back((double(t1-t0)/CLOCKS_PER_SEC));

			local.push_back(scoreTotal);

			destruir();
		}

		/*mejorG.push_back(max_element(greedy, greedy + greedy.size()));
		mejorL.push_back(max_element(local, local + local.size()));
		peorG.push_back(min_element(greedy, greedy + greedy.size()));
		peorL.push_back(min_element(local, local + local.size()));*/

		mejorG.push_back(max(greedy));
		mejorL.push_back(max(local));
		peorG.push_back(min(greedy));
		peorL.push_back(min(local));

		double desvG = 0, desvL = 0, mG = 0, mL = 0, tG = 0, tL = 0;

		for (int k = 0; k < 1; k++){
			desvG += (mejorValor[i - 1] - greedy[k]) / mejorValor[i - 1];
			//desvG += (mejorValor.at(i - 1) - greedy.at(k)) / mejorValor.at(i - 1);
			desvL += (mejorValor[i - 1] - local[k]) / mejorValor[i - 1];
			//desvL += (mejorValor.at(i - 1) - local.at(k)) / mejorValor.at(i - 1);
			mG += greedy[k];
			//mG += greedy.at(k);
			mL += local[k];
			//mL += local.at(k);
			tG += tiempoG[k];
			//tG += tiempoG.at(k);
			tL += tiempoL[k];
			//tL += tiempoL.at(k);
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
	}

	for (int i = 0; i < 30; i++){
		cout << "Greedy: " << endl << "		Caso " << i + 1 << ":" << endl;
		cout << "			Mejor: " << mejorG.at(i) << " Peor: " << peorG.at(i);
		cout << " 	Desviación Greedy: " << desviacionG.at(i) << " Media: " << mediaG.at(i);
		cout << " Media de tiempo: " << mediaTiempoG.at(i) << endl;
		cout << "Búsqueda Local: " << endl << "		Caso " << i + 1 << ":" << endl;
		cout << "			Mejor: " << mejorL.at(i) << " Peor: " << peorL.at(i);
		cout << " 	Desviación Local: " << desviacionL.at(i) << " Media: " << mediaL.at(i);
		cout << " Media de tiempo: " << mediaTiempoL.at(i) << endl;
	}

	/*input("Tablas/MDG-a_21_n2000_m200.txt");

	algoritmoGreedy();
	score();
	busquedaLocal();*/

 return 0;
}
