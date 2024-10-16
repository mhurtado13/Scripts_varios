#include <iostream>

using namespace std;

int main(){
    int numero, dato {5};

    cout << "Digite el numero: "; cin >> numero;

    if(numero == dato){
        cout << "El numero es igual a 5" << endl;
    }else{
        cout << "El numero es diferente a 5" << endl;
    }

    switch(numero){
        case 1: cout << "Es el numero 1" << endl; break;
        case 2: cout << "Es el numero 2" << endl; break;
        default: cout << "No esta en el rango de 1-2" << endl; break;
    }

    int i {0};

    while(i<=10){
        cout << i << endl;
        i++;
    }

    int x {0};

    do
    {
        cout << x << endl;
        x++;
    } while (x<=10);
    
    for (int iterator = 0; iterator <=10; iterator++)
    {
       cout << iterator << endl;
    }
    

    return 0;
}