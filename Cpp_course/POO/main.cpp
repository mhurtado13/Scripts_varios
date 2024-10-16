#include <iostream>
using namespace std;

class Persona{
    private: //Not modifiable but only by the methods
        // Attributes
        int edad;
        string nombre;
    public: //User can access
        Persona(int, string); //Constructor
        // Methods
        void leer();
        void correr();
};

//Constructor - Initialize the class
Persona::Persona(int _edad, string _nombre){
    edad = _edad;
    nombre = _nombre;
}

//Define the methods (no need to declare variables cause they belong to the class - they use the attributes)
void Persona::leer(){
    cout << "Soy " << nombre << " y estoy leyendo un libro" << endl;
}

void Persona::correr(){
    cout << "Soy " << nombre << " y estoy corriend una marathon" << endl;
}

int main(){
    //Create objects
    Persona p1 = Persona(25, "Marcelo");
    Persona p2(19, "Maria");

    //Perform actions/methods
    p1.correr();
    p2.leer();


    return 0;
}