#include <iostream>

int main(){
    int number1 = 15; //Decimal
    int number2 = 017; //Octal
    int number3 = 0x0f; //Hexadecimal
    int number4 = 0b00001111; //Binary

    std::cout << "number1: " << number1 << std::endl;
    std::cout << "number2: " << number2 << std::endl;
    std::cout << "number3: " << number3 << std::endl;
    std::cout << "number4: " << number4 << std::endl;

    int elephant_count; //Variable may contain random garbage value
    int lion_count{}; //Initializes to zero
    int dog_count {10}; //Initializes to 10
    int cat_count {15}; //Initializes to 15

    int domesticated_animals {dog_count + cat_count}; //Can use expression as initializer
    
    std::cout << "number1: " << elephant_count << std::endl;
    std::cout << "number2: " << lion_count << std::endl;
    std::cout << "number3: " << dog_count << std::endl;
    std::cout << "number4: " << domesticated_animals << std::endl;

    int value1 {10};
    int value2 {-300};

    std::cout << "value1: " << value1 << std::endl;
    std::cout << "value2: " << value2<< std::endl;
    std::cout << "sizeofvalue1: " << sizeof(value1) << std::endl;
    std::cout << "sizeofvalue2: " << sizeof(value2) << std::endl;

    signed int x {10};
    signed int y {-300};
    //'signed' means you can store positive and negative numbers
    //'unsigned' means you can only store positive numbers

    float x1 {1.12345678};
    double x2 {1.232523532525};
    long double x3 {1.13239140128421};

    std::cout << "sizeof float: " << sizeof(float) << std::endl;
    std::cout << "sizeof double: " << sizeof(double) << std::endl;
    std::cout << "sizeof long double: " << sizeof(long double) << std::endl;

    bool red_light {true};
    bool green_light {false};

    if(red_light == true){
        std::cout << "Stop!" << std::endl;
    }else{
        std::cout << "Go through!" << std::endl;
    }

    char character1 {'a'};
    char character2 {'b'};

    std::cout << character1 << std::endl;
    std::cout << character2 << std::endl;

    //compiler guess the type of variable
    auto var1 {12};
    auto var2 {12.0};


    return 0;
}