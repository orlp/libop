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

    op::print(op::isqrt(10));
    op::print(op::lcm(-24, 60));

    std::vector<int> primes;
    op::primesbelow(102, std::back_inserter(primes));
    for (int i : op::range(100)) {
        op::print(i, op::isprime(i));
    }

    op::print(op::join(primes.begin(), primes.end()));

    std::cout << op::randint(0, 2, rng) << "\n";
    std::cout << op::randint(0, 2, rng) << "\n";

    std::cout << std::numeric_limits<float>::digits << "\n";
    std::cout << std::numeric_limits<double>::digits << "\n";

    auto a = op::to_array("abcde");
    std::cout << "array:\n";

    for (auto al : a) { std::cout << al << " "; }
    std::cout << "\n";



    for (int i = 0; i < 500; ++i)
    std::cout << *op::random_choice(a.begin(), a.end(), rng);

    unsigned char out[4] = {0};
    for (int i = 0; i < 5; ++i) {
        op::random_sample(a.begin(), a.end(), out, 3, rng);
        std::cout << out << "\n";
    }



    op::Image img(15, 15);
    
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

    op::print(op::startswith("test", "te"), op::startswith("hello", "te"));
    op::print(op::endswith("test", "st"), op::endswith("hello", "he"));
    op::print(op::endswith("", ""), op::endswith("", ""));
    op::print(op::endswith("hah", ""), op::endswith("hah", ""));
    op::print(op::startswith("", ""), op::startswith("", ""));
    op::print(op::startswith("hah", ""), op::startswith("hah", ""));

    for (auto i : op::range(-1, 8ull)) {
        op::print(i);
    }
    
    std::vector<int> vx { 2, 3, 4, 5, 6};
    for (auto i : op::range(1, vx.size())) {
        op::print(i);
    }

    op::print(op::safe_less(-2, 2u));

    op::print();
    op::print();


    struct NonPrintable { };
    auto x = std::make_tuple(2, "test", 3, 8.4, NonPrintable());


    op::tuple_visit(x, 3, op::visit_forward<double&>()) = 42.42;
    auto d = op::tuple_visit(x, 3, op::visit_forward<double>());
    std::cout << d << "\n";

    std::cout << op::tuple_visit(x, 2, op::visit_forward<int>()) << "\n";         // 3
    std::cout << op::tuple_visit(x, 1, op::visit_forward<std::string>()) << "\n"; // test

    // can't implicitly convert from int to string
    try { std::cout << op::tuple_visit(x, 0, op::visit_forward<std::string>()) << "\n"; }
    catch (const std::exception& e) { std::cout << e.what() << "\n"; }

    // out of bounds
    try { op::tuple_visit(x, 5, op::visit_forward<int>()); }
    catch (const std::exception& e) { std::cout << e.what() << "\n"; }

    op::fformat(std::cout, "{0:.*1}\n", 5.1243, 3);

    op::fformat(std::cout, "The answer is {:.*f} {}.\n", 34.2348, 2, "USD");
    
    op::fformat(std::cout, "{:'*<10}\n", 34.2348);
    op::fformat(std::cout, "{:'*10}\n", 34.2348);
    op::fformat(std::cout, "{:'*<10}\n", "test");
    op::fformat(std::cout, "{:'*10}\n", "test");

    for (int i : op::range(10000)) {
        uint64_t n = op::randint<uint64_t>(0, -1, rng);
        auto factors = op::primefactors(n);
        std::string p = op::join(factors.begin(), factors.end());
        op::fprint("{}: {}\n", n, p);
    }

    auto f = op::primefactors(392853920582390413ull);
    op::print(op::join(f.begin(), f.end()));

    return 0;
}
