#include <iostream> // Library for entrance and output of data 

using namespace std; // Avoid doing "std::cout std::endl each time you print something"

int main(){
    int age;
    char sex;
    float height;

    cout << "Digite su edad: "; cin >> age;
    cout << "Sexo: M o F:"; cin >> sex;
    cout << "Height: "; cin >> height;

    cout << "XXXX" << endl;

    float a,b,resultado = 0;
    cout << "Digite el valor de a: "; cin >> a;
    cout << "Digite el valor de b: "; cin >> b;

    resultado = (a/b) + 1;

    cout.precision(2); //Round the output to 2 decimals
    cout << "El resultado es: " << resultado << endl;

    return 0; //Help the OS to detect when the function is finish
}
