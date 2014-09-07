#include <iostream>
#include <random>

#include "op.h"

int main(int argc, char **argv) {
    std::seed_seq sseq{1, 3};
    op::random_device rd;
    std::mt19937_64 rng(rd());

    std::cout << op::randint(0, 2, rng) << "\n";
    std::cout << op::randint(0, 2, rng) << "\n";


    std::cout << std::numeric_limits<float>::digits << "\n";
    std::cout << std::numeric_limits<double>::digits << "\n";

    std::array<char, 5> a { 'a', 'b', 'c', 'd', 'e' };

    for (int i = 0; i < 500; ++i)
    std::cout << *op::random_choice(a.begin(), a.end(), rng);

    unsigned char out[4] = {0};
    for (int i = 0; i < 5; ++i) {
        op::random_sample(a.begin(), a.end(), out, 3, rng);
        std::cout << out << "\n";
    }

    return 0;
}
