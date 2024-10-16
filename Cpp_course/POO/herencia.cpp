#include <iostream>

using namespace std;

class Persona{
    private:
        string nombre;
        int edad;
    public:
        Persona(string, int);
        void mostrarPersona();
};

class Alumno : public Persona{ //Hereda metodos de la clase Persona (todo lo 'public') 
    private:
        string codigo_alumno;
        float nota_final;
    public:
        Alumno(string, int, string, float); //Primero pones los atributos de la clase 'padre' y luego los de la clase hija
        void mostrarAlumno();
};

Persona::Persona(string _nombre, int _edad){
    nombre = _nombre;
    edad = _edad;
}

Alumno::Alumno(string _nombre, int _edad, string _codigo_alumno, float _nota_final):Persona(_nombre, _edad){
    codigo_alumno = _codigo_alumno;
    nota_final = _nota_final;
}

void Persona::mostrarPersona(){
    cout << "Nombre: " << nombre << endl;
    cout << "Edad: " << edad << endl;
}

void Alumno::mostrarAlumno(){
    mostrarPersona();
    cout << "Codigo Alumno: " << codigo_alumno << endl;
    cout << "Nota Final: " << nota_final << endl;
}


int main(){

    Alumno alumno1("Alejandro", 20, "1231313", 15.6);
    alumno1.mostrarAlumno();
    
    return 0;
}