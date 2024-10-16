#include <iostream>
using namespace std;

int main(){
    int numeros[5] = {1,2,3,4,5}; // "[]" define el numero de elementos
    int suma = 0;

    for (int i = 0; i <=4; i++)
    {
        suma += numeros[i];
    }
    
    cout << "La suma es " << suma << endl;

    int matrix[100][100], filas, columnas;

    cout << "Digite el numero de filas: ";
    cin >> filas;
    cout << "Digite el numero de columns: ";
    cin >> columnas;

    for (int i = 0; i <= filas; i++){
        for (int j = 0; j <= columnas; j++){
            cout << "Digite el numero: [" << i << "][" << j << "]";
            cin >> matrix[i][j];
        }   
    }

    for (int i = 0; i <= filas; i++){
        for (int j = 0; j <= columnas; j++){
            cout << matrix[i][j];
        }  
        cout << '\n'; 
    }
    
    return 0;
}