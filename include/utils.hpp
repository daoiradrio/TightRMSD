#ifndef UTILS_HPP
#define UTILS_HPP

#include <unordered_map>



std::unordered_map<std::string, int> element_numbers({
    {"H", 1},
    {"B", 5},
    {"C", 6},
    {"N", 7},
    {"O", 8},
    {"F", 9},
    {"Cl", 17},
    {"Br", 35},
    {"I", 53}
});



std::unordered_map<std::string, double> valence_radii_single({
    {"C", 0.75},
    {"N", 0.71},
    {"O", 0.66},
    {"H", 0.32},
    {"B", 0.85},
    {"F", 0.64},
    {"Cl", 0.99},
    {"Br", 1.14},
    {"I", 1.33}
});


std::unordered_map<std::string, double> valence_radii_double({
    {"C", 0.67},
    {"N", 0.60},
    {"O", 0.57},
    {"H", NULL},
    {"B", 0.78},
    {"F", NULL},
    {"Cl", NULL},
    {"Br", NULL},
    {"I", NULL}
});


std::unordered_map<std::string, double> valence_radii_triple({
    {"C", 0.60},
    {"N", 0.54},
    {"O", 0.53},
    {"H", NULL},
    {"B", 0.73},
    {"F", NULL},
    {"Cl", NULL},
    {"Br", NULL},
    {"I", NULL}
});



#endif
