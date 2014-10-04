#include <iostream>
#include <fstream>
#include <random>

#include "op.h"

struct A {};
struct B : A {};

void f(A) {}

int main(int argc, char **argv) {
    std::seed_seq sseq{1, 3};
    op::random_device rd;
    std::mt19937_64 rng(rd());

    f(B());

    std::cout << op::randint(0, 2, rng) << "\n";
    std::cout << op::randint(0, 2, rng) << "\n";

    std::cout << std::numeric_limits<float>::digits << "\n";
    std::cout << std::numeric_limits<double>::digits << "\n";

    auto a = op::to_array("abcde");
    std::cout << "array:\n";

    for (auto al : a) {
        std::cout << al << " ";
    }
    std::cout << "\n";

    for (int i = 0; i < 500; ++i)
    std::cout << *op::random_choice(a.begin(), a.end(), rng);

    unsigned char out[4] = {0};
    for (int i = 0; i < 5; ++i) {
        op::random_sample(a.begin(), a.end(), out, 3, rng);
        std::cout << out << "\n";
    }

    op::Image img(1500, 1500);
    
    for (int x = 0; x < img.width()/2; ++x) {
        for (int y = 0; y < img.height(); ++y) {
            img.set_pixel(x, y,
                (x^y)%255, ((x/2)^(y/2))%255, ((x/3)^(y/3))%255);
        }
    }


    for (int x = img.width() / 2; x < img.width(); ++x) {
        for (int y = 0; y < img.height(); ++y) {
            img.set_pixel(x, y,
                op::randint(0, 255, rng),op::randint(0, 255, rng),  op::randint(0, 255, rng)) ;
        }
    }

    std::ofstream file("test.png", std::ios::binary);
    img.write_png(file);

    return 0;
}
