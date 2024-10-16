#include <iostream>
#include <stdlib.h>

using namespace std;

class Punto{
    private: //Attributes
        int x, y;
    public: //Methods
        Punto();
        void setPunto(int, int); //Setter: set value
        int getPuntoX(); //Getter
        int getPuntoY(); //Getter
};

Punto::Punto(){
}

//Set value to the attributes
void Punto::setPunto(int _x, int _y){
    x = _x;
    y = _y;
}

int Punto::getPuntoX(){
    return x;
}

int Punto::getPuntoY(){
    return y;
}

int main(){
    Punto punto1;
    punto1.setPunto(15,10);

    return 0;
}