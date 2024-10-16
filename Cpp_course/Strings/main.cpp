#include <iostream>
#include <string.h> //Library for strings manipulation

using namespace std;

int main(){
    char palabra[] = "hola";
    int longitud = 0;

    longitud = strlen(palabra);

    cout << "La longitud de la palabra es " << longitud << endl;

    char nombre[] = "Alejandro";
    char nombre2[20];

    //Copiar un string hacia otro 'strcpy'
    strcpy(nombre2, nombre); 

    cout << nombre2 << endl;

    //Comparar dos strings 'strcmp'
    char string[] = "Hola";
    char string2[] = "Hola";

    if(strcmp(string, string2)==0){ // == 0 significa TRUE, == 1 significa FALSE
        cout << "Ambos strings son iguales" << endl;
    }

    //Concatenar una cadena con la otra
    char cad1[] = "Hola ";
    char cad2[] = "mundo";
    char cad3[30];

    strcpy(cad3, cad1);
    strcat(cad3, cad2);  

    cout << cad3 << endl;  
    
    return 0;
}