#include "stdio.h"
#include <string>
#include <iostream>

std::string tail(std::string a){
    return a.substr(1, a.length()-1);
}

int Levenshtein(std::string a, std::string b){
    if (a.length()==0){
        return b.length();
    }
    else if (b.length() == 0){
        return a.length();
    }
    if (a[0] == b[0]){
        return Levenshtein(tail(a), tail(b));
    }
    else{
        int temp1 = Levenshtein(tail(a), b);
        int temp2 = Levenshtein(a, tail(b));
        int temp3 = Levenshtein(tail(a), tail(b));
        return 1 + std::min(temp1, std::min(temp2, temp3));
    }

}

int LevDynamic(std::string a, std::string b){
    int len_a = a.length();
    int len_b = b.length();
    int distances[len_a+1][len_b+1];

    for (int i = 0; i<len_a + 1; i++ ){
        distances[i][0] = i;
    }
    for (int i = 0; i<len_b + 1; i++ ){
        distances[0][i] = i;
    }

    for (int i = 0; i<len_a + 1; i++){
        for (int j = 0; j<len_b + 1; j++){
            if (a[i-1] == b[j-1]){
                distances[i][j] = distances[i-1][j-1];
            }
            else {
                int m = distances[i][j-1];
                int n = distances[i-1][j];
                int o = distances[i-1][j-1];

                distances[i][j] = std::min(m, std::min(n, o)) + 1;
            }
        }
    }

    return distances[len_a][len_b];
}

int main(){
    std::string a = "Hello world";
    std::string b = "Hello world";

    int i = LevDynamic(a,b);

    std::cout<< i<< std::endl;;
    return 1;
}
