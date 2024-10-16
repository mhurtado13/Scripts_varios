
/* Basic workflow
C++ code --> COMPILER --> Executable binary file
*/

#include <iostream> //Library which includes functions like "cout", "std"
#include <string> //Library to include strings

int addNumbers(int first_number, int second_number){
    int sum = first_number + second_number;
    return sum;
}

int main(){
    int first_number = 12;
    int second_number = 6;

    std::cout << "First number: " << first_number << std::endl;
    std::cout << "Second number: " << second_number << std::endl;

    int sum = addNumbers(first_number, second_number);
    std::cout << "Sum: " << sum << std::endl;
    std::cout << "Sum: " << addNumbers(14, 26) << std::endl;

    int age;
    std::string name;

    std::cout << "Please type your name and age: " << std::endl;

    std::cin >> name >> age; // Input from console

    std::cout << "Hello " << name << " you are " << age << " years old!" << std::endl; 
    
    //Data with spaces
    
    /*
    std::string full_name;
    std::cout << "Please type your whole name" << std::endl;

    std::getline(std::cin, full_name); // Take full name input with spaces

    std::cin >> age;

    std::cout << "Hello " << full_name 
    << " you are " << age << " years old!" << std::endl; */

    return 0; //tell the OS that task run succesfully
}

/*
Errors types
- Compile time errors
- Runtime errors
- Warnings
*/ 

