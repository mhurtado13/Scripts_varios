#include <iostream>

using namespace std;

int encontrar_maximo(int var1, int var2);

int main(){

    int numero1, numero2;
    cout << "Digite dos numeros: ";
    cin >> numero1 >> numero2;

    cout << "El maximo es " << encontrar_maximo(numero1, numero2);
    
    return 0;
}

int encontrar_maximo(int var1, int var2){
    int numMax;

    if(var1 > var2){
        numMax = var1;
    }else{
        numMax = var2;
    }

    return numMax;
}